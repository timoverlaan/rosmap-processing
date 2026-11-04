"""UMAPs of the v2 count files, coloured by the labels of a batch-effects spec.

Reads a spec from config/batch_effects/ and uses its cell filters, cell types,
batch columns, covariates and disease label, so the plots show the labels the
batch-effect tables measure:

1. The cells as in `pipeline batch-effects` (genes per cell, cells without a
   donor, donors with too few cells, up to `max_cells_per_donor` random cells
   per donor), except that cells no cell-type rule matches are kept as "Other".
2. One gene set: log1p(counts / total * target_sum) and PCA without scaling
   over the pooled sample, as in `pipeline batch-effects`, then UMAP on the
   PCs with `k_neighbors` neighbours.
3. One PNG with a panel per label: cell type, dataset (several datasets),
   donor, batch columns, covariates (those of every dataset), disease, extra
   obs columns (`annotate`, e.g. the authors' subclasses), and the counts and
   genes per cell.

Donors, and batches with many levels, get no legend; no panel names a donor.
Outputs (default data/v2/umap/<spec name>/):
    umap_<gene set>.png      the panels
    umap_<gene set>.parquet  per sampled cell: dataset, row in its count file,
                             and the UMAP coordinates
    umap_<gene set>.json     provenance: spec, gene set, sizes, code version
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import h5py
import numpy as np
import pandas as pd

from ..utils.counts import log_normalize
from ..utils.logging import get_logger
from ..utils.obs_summary import obs_column_names, read_obs_column
from ..utils.provenance import provenance_record
from .batch_effects import (
    _as_labels,
    _matrix_node,
    _pca,
    load_spec,
    match_donors,
    prepare_dataset,
    resolve_gene_sets,
    stream_counts,
    subsample,
)

logger = get_logger(__name__)

DEFAULT_GENE_SETS = ("seurat_v3_donor_2000", "mit_seurat_v3_donor_2000")
DEFAULT_OUTPUT_ROOT = Path("data/v2/umap")
OTHER_CELLTYPE = "Other"
# The categorical palette of the dataviz reference instance, in its fixed order,
# and its blue ramp for numbers
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
OTHER_COLOR, MISSING_COLOR, INK = "#a8a7a0", "#e4e3de", "#2b2b29"
# Up to MAX_LEGEND levels each get a colour; up to MAX_FOLD, the largest
# MAX_LEGEND - 1 do and the rest are "Other"; more than that, every level gets a
# colour without a legend (donors, batches) or a text label (annotations).
MAX_LEGEND = len(CATEGORICAL)
MAX_FOLD = 20
MIN_LABEL_SHARE = 0.002   # annotation levels with fewer cells get no text label
UMAP_MIN_DIST = 0.3


@dataclass
class Panel:
    title: str
    values: np.ndarray
    kind: str   # "category", "numeric", "identity" (no legend) or "annotation" (text labels)


def pick_gene_set(spec: Dict[str, Any], name: Optional[str]) -> str:
    """The named gene set, else the first of DEFAULT_GENE_SETS in the spec."""
    if name:
        if name not in spec["gene_sets"]:
            raise ValueError(f"Gene set {name} is not in the spec; it has {sorted(spec['gene_sets'])}")
        return name
    for candidate in DEFAULT_GENE_SETS:
        if candidate in spec["gene_sets"]:
            return candidate
    raise ValueError(f"The spec has none of {DEFAULT_GENE_SETS}; pass --gene-set")


def _numeric_or_none(values: np.ndarray) -> Optional[np.ndarray]:
    """Numbers, if nearly all present values are numbers with many levels."""
    text = pd.Series(values, dtype=object)
    present = text.notna()
    numbers = pd.to_numeric(text.astype(str).str.replace("+", "", regex=False), errors="coerce")
    if present.sum() and numbers[present].notna().mean() > 0.95 and numbers.nunique() > MAX_FOLD:
        return numbers.to_numpy(np.float64)
    return None


def _categorical_panel(title: str, values: np.ndarray, annotation: bool = False) -> Panel:
    n_levels = pd.Series(values, dtype=object).nunique()
    if n_levels <= MAX_FOLD:
        return Panel(title, values, "category")
    return Panel(f"{title} ({n_levels} levels)", values, "annotation" if annotation else "identity")


def build_panels(sample: pd.DataFrame, numeric: Dict[str, bool], batch_columns: List[str],
                 extra: List[str], several: bool) -> List[Panel]:
    """The panels, in a fixed order, from the columns of the pooled sample."""
    panels = [_categorical_panel("cell type", sample["celltype"].to_numpy(), annotation=True)]
    if several:
        panels.append(Panel("dataset", sample["dataset"].to_numpy(), "category"))
    n_donors = sample["donor"].nunique()
    panels.append(Panel(f"donor ({n_donors}, no legend)", sample["donor"].to_numpy(), "identity"))
    for col in batch_columns:
        panel = _categorical_panel(col.replace("batch:", "batch: "), sample[col].to_numpy())
        if panel.kind == "identity":
            panel.title += ", no legend"
        panels.append(panel)
    for key, is_numeric in numeric.items():
        values = sample[key].to_numpy()
        panels.append(Panel(key, values.astype(np.float64), "numeric") if is_numeric
                      else _categorical_panel(key, values))
    if "disease" in sample:
        panels.append(Panel("disease", sample["disease"].to_numpy(), "category"))
    for col in extra:
        values = sample[col].to_numpy()
        as_numbers = _numeric_or_none(values)
        panels.append(Panel(col, as_numbers, "numeric") if as_numbers is not None
                      else _categorical_panel(col, values, annotation=True))
    panels.append(Panel("counts per cell (log10)", np.log10(np.maximum(sample["total_counts"], 1)), "numeric"))
    panels.append(Panel("genes per cell", sample["n_genes"].to_numpy(np.float64), "numeric"))
    return panels


def _level_order(values: pd.Series) -> List[str]:
    """Levels by number if they all are numbers, else by size."""
    counts = values.value_counts()
    numbers = pd.to_numeric(pd.Series(counts.index), errors="coerce")
    if numbers.notna().all():
        return [counts.index[i] for i in np.argsort(numbers.to_numpy(), kind="stable")]
    return list(counts.index)


def category_colors(values: np.ndarray):
    """Per-cell colours and the legend (label, colour) for at most MAX_FOLD
    levels: the largest MAX_LEGEND - 1 in palette order, the rest as "Other",
    when there are more than MAX_LEGEND."""
    s = pd.Series(values, dtype=object)
    levels = _level_order(s.dropna())
    if len(levels) > MAX_LEGEND:
        by_size = list(s.value_counts().index[:MAX_LEGEND - 1])
        kept = [lv for lv in levels if lv in by_size]
    else:
        kept = levels
    palette = iter(CATEGORICAL)
    colors = {lv: OTHER_COLOR if lv == OTHER_CELLTYPE else next(palette) for lv in kept}
    per_cell = s.map(colors).where(s.notna(), MISSING_COLOR).fillna(OTHER_COLOR).to_numpy()
    legend = [(lv, colors[lv]) for lv in kept]
    if len(kept) < len(levels):
        legend.append((f"Other ({len(levels) - len(kept)})", OTHER_COLOR))
    if s.isna().any():
        legend.append(("missing", MISSING_COLOR))
    return per_cell, legend


def identity_colors(values: np.ndarray, seed: int) -> np.ndarray:
    """A random hue per level (no legend), so neighbouring levels differ."""
    import matplotlib.colors as mcolors

    s = pd.Series(values, dtype=object)
    codes, uniques = pd.factorize(s, use_na_sentinel=True)
    rng = np.random.default_rng(seed)
    hsv = np.column_stack([rng.permutation(len(uniques)) / max(len(uniques), 1),
                           rng.uniform(0.45, 0.85, len(uniques)), rng.uniform(0.55, 0.9, len(uniques))])
    table = np.vstack([mcolors.hsv_to_rgb(hsv), mcolors.to_rgb(MISSING_COLOR)])
    return table[codes]


def plot_panels(coords: np.ndarray, panels: List[Panel], path: Path, title: str, seed: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    n = len(coords)
    order = np.random.default_rng(seed).permutation(n)
    size = float(np.clip(3e4 / n, 0.05, 4.0))
    ncols = min(4, len(panels))
    nrows = -(-len(panels) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 5.4 * nrows), squeeze=False)
    fig.patch.set_facecolor("white")
    cmap = mcolors.LinearSegmentedColormap.from_list("ramp", SEQUENTIAL)
    for ax, panel in zip(axes.flat, panels):
        missing = pd.isna(pd.Series(panel.values, dtype=object)).to_numpy()
        # Missing values first, so they lie under the rest
        draw = np.concatenate([order[missing[order]], order[~missing[order]]])
        ax.set_title(f"{panel.title}\n{int((~missing).sum()):,} of {n:,} cells with a value",
                     fontsize=9, color=INK)
        if panel.kind == "numeric":
            values = panel.values.astype(np.float64)
            ok = ~missing
            ax.scatter(*coords[draw[missing[draw]]].T, s=size, c=MISSING_COLOR, linewidths=0, rasterized=True)
            if ok.any():
                lo, hi = np.nanpercentile(values[ok], [1, 99])
                points = ax.scatter(*coords[draw[ok[draw]]].T, s=size, c=values[draw[ok[draw]]], cmap=cmap,
                                    vmin=lo, vmax=hi if hi > lo else lo + 1, linewidths=0, rasterized=True)
                bar = fig.colorbar(points, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03)
                bar.ax.tick_params(labelsize=7, colors=INK)
                bar.outline.set_visible(False)
        elif panel.kind == "category":
            colors, legend = category_colors(panel.values)
            ax.scatter(*coords[draw].T, s=size, c=list(colors[draw]), linewidths=0, rasterized=True)
            handles = [Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=c,
                              markeredgewidth=0, label=str(label)) for label, c in legend]
            ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=min(3, len(handles)),
                      fontsize=7, frameon=False, handletextpad=0.2, columnspacing=0.8, labelcolor=INK)
        else:
            colors = identity_colors(panel.values, seed)
            ax.scatter(*coords[draw].T, s=size, c=colors[draw], linewidths=0, rasterized=True)
            if panel.kind == "annotation":
                s = pd.Series(panel.values, dtype=object)
                share = s.value_counts(normalize=True)
                for level in share.index[share >= MIN_LABEL_SHARE]:
                    x, y = np.median(coords[(s == level).to_numpy()], axis=0)
                    ax.text(x, y, str(level), fontsize=5, color=INK, ha="center", va="center",
                            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    fig.suptitle(title, fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def _read_cell_columns(ds, columns: Sequence[str]) -> Dict[str, np.ndarray]:
    """Genes per cell, and the listed obs columns the file has, for the sampled cells."""
    rows = ds.sample["row"].to_numpy()
    with h5py.File(ds.path, "r") as f:
        n_genes = np.diff(_matrix_node(f, ds.cfg.get("layer"))["indptr"][()])[rows]
        have = set(obs_column_names(f))
        out = {"n_genes": n_genes}
        for col in columns:
            if col in have:
                out[col] = _as_labels(read_obs_column(f, col))[rows]
            else:
                logger.info(f"[{ds.name}] obs has no column {col}; not plotted for this dataset")
    return out


def run_umap(spec_path: Union[str, Path], gene_set: Optional[str] = None, annotate: Sequence[str] = (),
             output: Optional[Union[str, Path]] = None, max_cells_per_donor: Optional[int] = None,
             overwrite: bool = False) -> Path:
    spec = load_spec(spec_path)
    gs = pick_gene_set(spec, gene_set)
    spec["gene_sets"] = {gs: spec["gene_sets"][gs]}
    if max_cells_per_donor:
        spec["max_cells_per_donor"] = max_cells_per_donor
    out = Path(output) if output else DEFAULT_OUTPUT_ROOT / spec["name"]
    png = out / f"umap_{gs}.png"
    if png.exists() and not overwrite:
        raise FileExistsError(f"{png} exists (pass --overwrite to replace it)")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(spec["seed"])
    logger.info(f"UMAP '{spec['name']}', gene set {gs} -> {out}")

    datasets = []
    for name, cfg in spec["datasets"].items():
        # Keep the cells no rule matches, as "Other"
        catch_all = {"column": cfg["celltype"][0]["column"], "match": ".", "label": OTHER_CELLTYPE}
        datasets.append(prepare_dataset(name, {**cfg, "celltype": list(cfg["celltype"]) + [catch_all]}, spec))
    several = len(datasets) > 1
    if several and spec["match_donors"]:
        match_donors(datasets)
    elif several:
        for ds in datasets:
            ds.cells["donor"] = ds.name + ":" + ds.cells["donor"].astype(str)
    for ds in datasets:
        subsample(ds, spec["max_cells_per_donor"], rng)
    positions = resolve_gene_sets(spec, datasets)
    if gs not in positions:
        raise ValueError(f"Gene set {gs} is not usable (see the log)")
    for ds in datasets:
        logger.info(f"[{ds.name}] streaming counts for {len(ds.union)} genes")
        stream_counts(ds, spec["chunk_size"])
        for col, values in _read_cell_columns(ds, annotate).items():
            ds.sample[col] = values
        ds.sample["total_counts"] = ds.totals

    sample = pd.concat([ds.sample.assign(dataset=ds.name) for ds in datasets], ignore_index=True)
    X = np.vstack([log_normalize(ds.counts[:, positions[gs][ds.name]], ds.totals, spec["target_sum"]).toarray()
                   for ds in datasets])
    logger.info(f"PCA of {X.shape[0]} cells x {X.shape[1]} genes")
    pcs = _pca(X, spec["n_pcs"], spec["seed"])
    del X
    import umap

    logger.info(f"UMAP on {pcs.shape[1]} PCs, {spec['k_neighbors']} neighbours")
    coords = umap.UMAP(n_neighbors=spec["k_neighbors"], min_dist=UMAP_MIN_DIST,
                       random_state=spec["seed"]).fit_transform(pcs)

    covariates = [k for k in datasets[0].cfg["covariates"] if all(k in ds.cfg["covariates"] for ds in datasets)]
    numeric = {k: bool(datasets[0].cfg["covariates"][k].get("numeric")) for k in covariates}
    batch_columns = [f"batch:{c}" for c in datasets[0].cfg["batch"]] if not several else []
    extra = [c for c in annotate if c in sample.columns]
    panels = build_panels(sample, numeric, batch_columns, extra, several)
    n_genes = len(positions[gs][datasets[0].name])
    title = (f"{spec['name']}: {len(sample):,} cells from {sample['donor'].nunique():,} donors "
             f"(up to {spec['max_cells_per_donor']} per donor); gene set {gs} ({n_genes} genes), "
             f"{pcs.shape[1]} PCs")
    plot_panels(coords, panels, png, title, spec["seed"])
    logger.info(f"Wrote {png} ({len(panels)} panels)")

    pd.DataFrame({"dataset": sample["dataset"], "row": sample["row"],
                  "umap1": coords[:, 0], "umap2": coords[:, 1]}).to_parquet(out / f"umap_{gs}.parquet")
    record = provenance_record("pipeline umap", {"spec": str(spec_path), "gene_set": gs,
                                                 "annotate": ",".join(annotate)})
    record.update(n_cells=len(sample), n_genes=n_genes, n_pcs=int(pcs.shape[1]),
                  datasets={ds.name: ds.stats for ds in datasets})
    (out / f"umap_{gs}.json").write_text(json.dumps(record, indent=2, default=str))
    return png
