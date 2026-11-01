"""Batch effects: how much of the variation in the v2 data is technical.

Driven by a YAML spec (see config/batch_effects/). Works from the count files,
not the processed outputs, so every gene set is compared on the same cells
with the same normalization:

1. Per dataset, the cell and donor filters of `pipeline scanpy` (genes per
   cell, cells without a donor, donors with too few cells), cells mapped to a
   broad cell type by the spec's rules (unmatched cells are dropped), then up
   to `max_cells_per_donor` random cells per donor.
2. One pass over the count matrix in row chunks: per-cell totals over all
   genes (the size factors), the sampled cells' counts for the genes of all
   gene sets, and pseudobulk sums per donor x cell type over all kept cells.
3. Per gene set: log1p(counts / total * target_sum) and PCA without scaling,
   as in `pipeline scanpy`, over all sampled cells and within each cell type;
   the cell-level metrics of core/batch_metrics.py; pseudobulk log-CPM for the
   donor-level metrics.

One dataset: its technical batch columns, against donor, cell type and the
donor covariates. Several: the dataset label within each broad cell type; with
match_donors, only donors present in every dataset (matched on `match_on`),
plus the paired dataset-vs-donor decomposition.

Most technical batches hold whole donors, so a batch label also carries donor
differences. Every batch, dataset and covariate metric therefore comes with a
null: the same label shuffled between donors (within donors, for the dataset
label of matched donors). `excess` = value - null mean.

Outputs in the spec's output_dir. All are aggregates: rows are gene sets x
cell types x labels, never donors; cell types with fewer than min_donors
donors (in any dataset) are left out.
    cells.tsv          cells and donors per dataset and cell type
    cell_metrics.tsv   PC regression, neighbourhood mixing, kBET-style
                       rejection, silhouettes
    donor_metrics.tsv  pseudobulk partial R², classifier accuracy, paired
                       decomposition
    confounding.tsv    association of batch (or dataset) with the covariates
    summary.tsv        the rows of both metric tables for all cells / the
                       mean over cell types
    provenance.json    the spec, gene set sizes and code version
"""

import json
import re
import time
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from ..core import batch_metrics as bm
from ..core.scanpy_pipeline import MISSING_DONOR_IDS, load_gene_list
from ..utils.constants import (
    DEFAULT_K_NEIGHBORS,
    DEFAULT_PCA_COMPONENTS,
    MIN_CELLS_PER_DONOR,
    MIN_GENES_PER_CELL,
    NORMALIZE_TARGET_SUM,
)
from ..utils.counts import count_non_integer_values, log_normalize, total_counts
from ..utils.logging import get_logger
from ..utils.provenance import provenance_record

try:
    from anndata.io import read_elem, sparse_dataset
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem, sparse_dataset

logger = get_logger(__name__)

DEFAULTS: Dict[str, Any] = {
    "seed": 0,
    "max_cells_per_donor": 300,
    "min_genes": MIN_GENES_PER_CELL,
    "min_cells_per_donor": MIN_CELLS_PER_DONOR,
    "target_sum": NORMALIZE_TARGET_SUM,
    "n_pcs": DEFAULT_PCA_COMPONENTS,
    "k_neighbors": DEFAULT_K_NEIGHBORS,
    "n_permutations": 20,
    "min_donors": 5,
    "min_cells_per_pseudobulk": 10,
    "max_cells_silhouette": 10000,
    "match_donors": False,
    "chunk_size": 50000,
}
REQUIRED = ("output_dir", "datasets", "gene_sets")
DATASET_KEYS = {"counts", "layer", "gene_map", "donor", "match_on", "celltype", "batch", "covariates"}


# ---------------------------------------------------------------- spec


def load_spec(path: Union[str, Path]) -> Dict[str, Any]:
    """The YAML spec with defaults filled in; unknown keys are an error."""
    with open(path) as f:
        raw = yaml.safe_load(f) or {}
    unknown = set(raw) - set(DEFAULTS) - set(REQUIRED) - {"name"}
    if unknown:
        raise ValueError(f"Unknown key(s) in {path}: {sorted(unknown)}")
    missing = [k for k in REQUIRED if not raw.get(k)]
    if missing:
        raise ValueError(f"{path} needs {missing}")
    spec = {**DEFAULTS, "name": Path(path).stem, **raw}
    for name, cfg in spec["datasets"].items():
        bad = set(cfg) - DATASET_KEYS
        if bad:
            raise ValueError(f"Unknown key(s) for dataset {name}: {sorted(bad)}")
        for key in ("counts", "donor", "celltype"):
            if not cfg.get(key):
                raise ValueError(f"Dataset {name} needs '{key}'")
        cfg.setdefault("batch", [])
        cfg.setdefault("covariates", {})
        if spec["match_donors"] and not cfg.get("match_on"):
            raise ValueError(f"match_donors needs 'match_on' for dataset {name}")
    if len(spec["datasets"]) > 1 and any(cfg["batch"] for cfg in spec["datasets"].values()):
        logger.warning("Batch columns are only analysed with a single dataset; ignoring them")
    return spec


def _first_existing(paths) -> Path:
    paths = [paths] if isinstance(paths, str) else list(paths)
    for p in paths:
        if Path(p).exists():
            return Path(p)
    raise FileNotFoundError(f"None of these count files exists: {paths}")


# ---------------------------------------------------------------- obs values


def _label_key(value) -> Optional[str]:
    """An obs value as a label: whole floats lose '.0'; missing markers -> None."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        value = int(value)
    text = str(value).strip()
    return None if text.lower() in MISSING_DONOR_IDS else text


def _as_labels(values, recode: Optional[Dict] = None) -> np.ndarray:
    """Object array of labels (NaN for missing), optionally recoded (values not
    in the recode map become NaN)."""
    codes, uniques = pd.factorize(pd.Series(values), use_na_sentinel=True)
    keys = [_label_key(u) for u in uniques]
    if recode is not None:
        recode = {str(k): v for k, v in recode.items()}
        keys = [recode.get(k) if k is not None else None for k in keys]
    table = np.array([np.nan if k is None else str(k) for k in keys] + [np.nan], dtype=object)
    return table[codes]


def _as_numeric(values) -> np.ndarray:
    """Numbers from obs values; '90+' -> 90 (the ROSMAP age cap)."""
    text = pd.Series(_as_labels(values), dtype=object)
    return pd.to_numeric(text.str.replace("+", "", regex=False), errors="coerce").to_numpy(np.float64)


def _assign_celltypes(columns: Dict[str, np.ndarray], rules: List[Dict], n: int) -> np.ndarray:
    """Broad cell type per cell: the first rule whose regex matches (case-insensitive)."""
    celltype = np.full(n, np.nan, dtype=object)
    for rule in rules:
        codes, uniques = pd.factorize(pd.Series(columns[rule["column"]]), use_na_sentinel=True)
        hits = np.array([bool(re.search(rule["match"], str(u), re.IGNORECASE)) for u in uniques] + [False])
        hit = hits[codes] & pd.isna(celltype)
        celltype[hit] = rule["label"]
    return celltype


def _donor_majority(cells: pd.DataFrame, column: str) -> pd.Series:
    """Most common non-missing value per donor."""
    def mode(s):
        s = s.dropna()
        return s.value_counts().index[0] if len(s) else np.nan
    return cells.groupby("donor", sort=True)[column].agg(mode)


# ---------------------------------------------------------------- datasets


@dataclass
class Dataset:
    name: str
    cfg: Dict[str, Any]
    path: Path
    var_names: List[str]
    cells: pd.DataFrame                  # kept cells; 'row' = position in the file
    stats: Dict[str, int]
    sample: Optional[pd.DataFrame] = None
    union: Optional[np.ndarray] = None   # file columns of all gene sets' genes
    counts: Optional[sp.csr_matrix] = None   # sample x union
    totals: Optional[np.ndarray] = None      # sample, over all genes
    groups: Optional[pd.DataFrame] = None    # donor x cell type
    pseudobulk: Optional[np.ndarray] = None  # groups x union


def _matrix_node(f: h5py.File, layer: Optional[str]):
    if layer is None:
        node = f["X"]
    elif "layers" in f and layer in f["layers"]:
        node = f["layers"][layer]
    else:
        raise ValueError(f"Layer '{layer}' not found")
    if not isinstance(node, h5py.Group) or node.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError("The count matrix must be stored as CSR")
    return node


def prepare_dataset(name: str, cfg: Dict[str, Any], spec: Dict[str, Any]) -> Dataset:
    """Read the needed obs columns and apply the filters; no count data yet."""
    path = _first_existing(cfg["counts"])
    covariates = cfg["covariates"]
    columns = {cfg["donor"], *cfg["batch"], *(r["column"] for r in cfg["celltype"]),
               *(c["column"] for c in covariates.values())}
    if cfg.get("match_on"):
        columns.add(cfg["match_on"])

    logger.info(f"[{name}] reading obs columns and genes per cell from {path}")
    with h5py.File(path, "r") as f:
        var = f["var"]
        var_names = [str(v) for v in read_elem(var[var.attrs["_index"]])]
        node = _matrix_node(f, cfg.get("layer"))
        # Stored entries per row; `pipeline scanpy` counts values > 0, which is
        # the same unless explicit zeros are stored.
        n_genes = np.diff(node["indptr"][()])
        obs = f["obs"]
        missing = sorted(c for c in columns if c not in obs)
        if missing:
            raise ValueError(f"[{name}] obs has no column(s) {missing}")
        raw = {c: read_elem(obs[c]) for c in columns}
    n = len(n_genes)

    stats = {"cells_in_file": n}
    keep = n_genes >= spec["min_genes"]
    stats["cells_min_genes"] = int(keep.sum())
    donors = _as_labels(raw[cfg["donor"]])
    no_donor = pd.isna(donors)
    stats["cells_without_donor"] = int((keep & no_donor).sum())
    keep &= ~no_donor
    per_donor = pd.Series(donors[keep]).value_counts()
    small = set(per_donor.index[per_donor < spec["min_cells_per_donor"]])
    too_small = keep & pd.Series(donors).isin(small).to_numpy()
    stats["donors_too_few_cells"] = len(small)
    stats["cells_of_those_donors"] = int(too_small.sum())
    keep &= ~too_small
    celltype = _assign_celltypes(raw, cfg["celltype"], n)
    stats["cells_without_celltype"] = int((keep & pd.isna(celltype)).sum())
    keep &= ~pd.isna(celltype)
    rows = np.flatnonzero(keep)

    cells = pd.DataFrame({"row": rows, "donor": donors[rows], "celltype": celltype[rows]})
    for col in cfg["batch"]:
        cells[f"batch:{col}"] = _as_labels(raw[col])[rows]
    for key, cov in covariates.items():
        values = _as_numeric(raw[cov["column"]]) if cov.get("numeric") else \
            _as_labels(raw[cov["column"]], cov.get("recode"))
        cells[key] = values[rows]
    if cfg.get("match_on"):
        cells["match"] = _as_labels(raw[cfg["match_on"]])[rows]
    stats["cells_kept"] = len(cells)
    stats["donors_kept"] = cells["donor"].nunique()
    logger.info(f"[{name}] " + ", ".join(f"{k}={v}" for k, v in stats.items()))
    for t, count in cells["celltype"].value_counts().items():
        logger.info(f"[{name}]   {t}: {count} cells")
    return Dataset(name=name, cfg=cfg, path=path, var_names=var_names, cells=cells, stats=stats)


def match_donors(datasets: List[Dataset]) -> None:
    """Keep donors present in every dataset; their donor ID becomes the match key."""
    keys = []
    for ds in datasets:
        per_donor = _donor_majority(ds.cells, "match")
        ds.cells["donor"] = ds.cells["donor"].map(per_donor)
        keys.append(set(per_donor.dropna()))
    common = set.intersection(*keys)
    for ds in datasets:
        before = ds.cells["donor"].nunique()
        ds.cells = ds.cells[ds.cells["donor"].isin(common)].reset_index(drop=True)
        logger.info(f"[{ds.name}] {len(common)} of {before} donors are in every dataset")
    if len(common) < 2:
        raise ValueError("Fewer than 2 donors are present in every dataset")


def subsample(ds: Dataset, max_per_donor: int, rng: np.random.Generator) -> None:
    shuffled = ds.cells.iloc[rng.permutation(len(ds.cells))]
    rank = shuffled.groupby("donor", sort=False).cumcount().to_numpy()
    ds.sample = shuffled[rank < max_per_donor].sort_values("row").reset_index(drop=True)
    logger.info(f"[{ds.name}] sampled {len(ds.sample)} cells (up to {max_per_donor} per donor)")


def stream_counts(ds: Dataset, chunk_size: int) -> None:
    """One pass over the matrix: sampled cells' counts for ds.union, their totals
    over all genes, and pseudobulk sums per donor x cell type."""
    kept = ds.cells["row"].to_numpy()
    sampled = ds.sample["row"].to_numpy()
    group_codes, groups = pd.MultiIndex.from_arrays(
        [ds.cells["donor"], ds.cells["celltype"]]).factorize()
    n_groups = len(groups)
    pseudobulk = np.zeros((n_groups, len(ds.union)))
    library = np.zeros(n_groups)
    totals = np.zeros(len(sampled))
    parts = []
    start_time = time.time()
    with h5py.File(ds.path, "r") as f:
        matrix = sparse_dataset(_matrix_node(f, ds.cfg.get("layer")))
        n_obs = matrix.shape[0]
        for start in range(0, n_obs, chunk_size):
            stop = min(start + chunk_size, n_obs)
            lo, hi = np.searchsorted(kept, [start, stop])
            if lo == hi:
                continue
            block = sp.csr_matrix(matrix[start:stop])[kept[lo:hi] - start]
            bad = count_non_integer_values(block)
            if bad:
                raise ValueError(f"[{ds.name}] rows {start}-{stop}: {bad} values are not counts")
            block_totals = total_counts(block)
            genes = block[:, ds.union]
            indicator = bm._indicator(group_codes[lo:hi], n_groups)
            pseudobulk += (indicator @ genes).toarray()
            library += indicator @ block_totals
            s_lo, s_hi = np.searchsorted(sampled, [start, stop])
            if s_hi > s_lo:
                local = np.searchsorted(kept[lo:hi], sampled[s_lo:s_hi])
                parts.append(genes[local])
                totals[s_lo:s_hi] = block_totals[local]
            if (start // chunk_size) % 10 == 0:
                logger.info(f"[{ds.name}] read {stop}/{n_obs} rows ({time.time() - start_time:.0f} s)")
    ds.counts = sp.vstack(parts).tocsr() if parts else sp.csr_matrix((0, len(ds.union)))
    ds.totals = totals
    ds.pseudobulk = pseudobulk
    ds.groups = pd.DataFrame({
        "donor": groups.get_level_values(0), "celltype": groups.get_level_values(1),
        "n_cells": np.bincount(group_codes, minlength=n_groups), "library_size": library,
    })


# ---------------------------------------------------------------- gene sets


def read_gene_set(path: Union[str, Path]) -> List[str]:
    """Genes of a TXT list, or of an h5ad (its Ensembl IDs if it has them, else
    its var_names). Ensembl version suffixes are dropped."""
    path = Path(path)
    if path.suffix == ".h5ad":
        from ..data.genes import read_var

        var = read_var(path)
        ids = var["gene_id"].astype(str) if "gene_id" in var.columns else None
        entries = ids.tolist() if ids is not None and ids.str.startswith("ENS").mean() > 0.5 \
            else var.index.tolist()
    else:
        entries = load_gene_list(path)
    return [e.split(".")[0] if e.startswith("ENS") else e for e in entries]


def _columns_for(ds: Dataset, entries: List[str]) -> Dict[str, int]:
    """entry -> file column, by Ensembl ID (via the gene map) or by symbol."""
    if sum(e.startswith("ENS") for e in entries) > len(entries) / 2:
        if not ds.cfg.get("gene_map"):
            raise ValueError(f"[{ds.name}] a gene set holds Ensembl IDs; set gene_map")
        from ..data.genes import load_gene_map

        mapping = load_gene_map(ds.cfg["gene_map"])
        lookup: Dict[str, int] = {}
        for col, name in enumerate(ds.var_names):
            gid = mapping.get(name)
            if gid:
                lookup.setdefault(gid.split(".")[0], col)
    else:
        lookup = {}
        for col, name in enumerate(ds.var_names):
            lookup.setdefault(name, col)
    return {e: lookup[e] for e in entries if e in lookup}


def resolve_gene_sets(spec: Dict[str, Any], datasets: List[Dataset]) -> Dict[str, Dict[str, np.ndarray]]:
    """Per gene set and dataset, the file columns of its genes (aligned across
    datasets; genes missing from any dataset are left out). Sets ds.union and
    returns positions into it."""
    columns: Dict[str, Dict[str, np.ndarray]] = {}
    for gs, path in spec["gene_sets"].items():
        if not Path(path).exists():
            logger.warning(f"Gene set {gs}: {path} not found (skipped)")
            continue
        entries = read_gene_set(path)
        found = [_columns_for(ds, entries) for ds in datasets]
        shared = [e for e in dict.fromkeys(entries) if all(e in f for f in found)]
        if len(shared) < len(set(entries)):
            logger.warning(f"Gene set {gs}: {len(set(entries)) - len(shared)} of {len(set(entries))} "
                           f"genes are missing from a dataset; using {len(shared)}")
        if len(shared) < spec["n_pcs"] + 1:
            logger.warning(f"Gene set {gs}: only {len(shared)} genes (skipped)")
            continue
        columns[gs] = {ds.name: np.array([f[e] for e in shared]) for ds, f in zip(datasets, found)}
        logger.info(f"Gene set {gs}: {len(shared)} genes from {path}")
    if not columns:
        return {}
    for ds in datasets:
        ds.union = np.unique(np.concatenate([c[ds.name] for c in columns.values()]))
    return {gs: {ds.name: np.searchsorted(ds.union, by_ds[ds.name]) for ds in datasets}
            for gs, by_ds in columns.items()}


# ---------------------------------------------------------------- metrics


@dataclass
class Label:
    column: str
    numeric: bool = False
    null: Optional[str] = None   # None, "between" or "within" (donors)
    mixing: bool = False         # also neighbourhood mixing, kBET and silhouette


def _pca(X: np.ndarray, n_pcs: int, seed: int) -> np.ndarray:
    from sklearn.decomposition import PCA

    n = min(n_pcs, X.shape[0] - 1, X.shape[1] - 1)
    return PCA(n_components=n, svd_solver="randomized", random_state=seed).fit_transform(X)


def _row(**kw) -> Dict[str, Any]:
    if "null_mean" in kw:
        kw["excess"] = kw["value"] - kw["null_mean"]
        kw["z"] = bm.z_score(kw["value"], kw["null_mean"], kw.get("null_sd", np.nan))
    return kw


def _null(metric, labels, donors, label: Label, spec, rng):
    if label.null is None:
        return np.nan, np.nan
    values = [metric(bm.permute_labels(labels, donors, rng, within=label.null == "within"))
              for _ in range(spec["n_permutations"])]
    return bm.null_summary(values)


def cell_metrics(gs: str, X: np.ndarray, cells: pd.DataFrame, labels: Dict[str, Label],
                 celltypes: List[str], spec, rng) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    scopes = [("all", np.arange(len(cells)))] + \
        [(t, np.flatnonzero(cells["celltype"].to_numpy() == t)) for t in celltypes]
    for scope, idx in scopes:
        pcs = _pca(X[idx], spec["n_pcs"], spec["seed"])
        sub = cells.iloc[idx]
        donors = sub["donor"].to_numpy()
        base = dict(gene_set=gs, celltype=scope, n_cells=len(idx), n_donors=len(set(donors)))
        knn = None
        for name, label in labels.items():
            if scope != "all" and name == "celltype":
                continue
            values = sub[label.column].to_numpy()
            ok = ~pd.isna(values)
            if ok.sum() < 3 or len(set(values[ok])) < 2:
                continue
            v, d, p = values[ok], donors[ok], pcs[ok]
            if label.numeric:
                v = v.astype(np.float64)

            def pcr(lab, p=p, numeric=label.numeric):
                return bm.pc_regression(p, np.asarray(lab, dtype=np.float64) if numeric else lab, numeric)

            null_mean, null_sd = _null(pcr, v, d, label, spec, rng)
            n_levels = np.nan if label.numeric else len(set(v))
            rows.append(_row(**base, label=name, metric="pc_regression", value=pcr(v),
                             null_mean=null_mean, null_sd=null_sd, n_levels=n_levels))
            if name == "celltype":
                rows.append(_row(**base, label=name, metric="silhouette", value=bm.label_asw(
                    pcs, values, rng, spec["max_cells_silhouette"]), n_levels=n_levels))
            if scope == "all" or not label.mixing:
                continue
            if ok.all():
                if knn is None:
                    knn = bm.knn_indices(pcs, spec["k_neighbors"])
                k = knn
            else:  # some cells lack this label: neighbours among the labelled cells
                k = bm.knn_indices(p, spec["k_neighbors"])
            ratio = lambda lab, k=k: bm.mixing_ratio(k, lab)[0]
            null_mean, null_sd = _null(ratio, v, d, label, spec, rng)
            rows.append(_row(**base, label=name, metric="mixing_ratio", value=ratio(v),
                             null_mean=null_mean, null_sd=null_sd, n_levels=n_levels))
            kbet = lambda lab, k=k: bm.knn_chi2_rejection(k, lab)
            null_mean, null_sd = _null(kbet, v, d, label, spec, rng)
            rows.append(_row(**base, label=name, metric="knn_rejection", value=kbet(v),
                             null_mean=null_mean, null_sd=null_sd, n_levels=n_levels))
            rows.append(_row(**base, label=name, metric="batch_silhouette", value=bm.batch_asw(
                p, v, np.zeros(len(v)), rng, spec["max_cells_silhouette"]), n_levels=n_levels))
    return rows


def _log_cpm(pseudobulk: np.ndarray, library: np.ndarray) -> np.ndarray:
    return np.log2(pseudobulk / library[:, None] * 1e6 + 1)


def donor_metrics(gs: str, table: pd.DataFrame, Y: np.ndarray, terms: Dict[str, Label],
                  targets: List[str], celltype: str, spec, rng) -> List[Dict[str, Any]]:
    """Partial R² (and R² alone) of each term in a pseudobulk model, and
    classifier accuracy for each target; nulls permute the term between rows."""
    rows: List[Dict[str, Any]] = []
    base = dict(gene_set=gs, celltype=celltype)
    usable = {}
    for name, label in terms.items():
        values = table[label.column]
        n_levels = values.nunique()
        if n_levels < 2 or (not label.numeric and n_levels > len(table) / 2):
            logger.info(f"{gs} / {celltype}: term {name} has {n_levels} levels for "
                        f"{len(table)} rows; left out of the model")
            continue
        usable[name] = label
    complete = table[[label.column for label in usable.values()]].notna().all(axis=1).to_numpy()
    covs = table.loc[complete, [label.column for label in usable.values()]]
    covs.columns = list(usable)
    Yc = Y[complete]
    numeric = [name for name, label in usable.items() if label.numeric]
    if len(covs) > 3 and usable:
        for name in usable:
            def stats_for(frame):
                Xf, cols = bm.design_matrix(frame, numeric)
                alone, _ = bm.design_matrix(frame[[name]], numeric)
                return (float(np.median(bm.partial_r2(Yc, Xf, cols[name]))),
                        float(np.median(bm.partial_r2(Yc, alone, np.arange(1, alone.shape[1])))))
            partial, alone = stats_for(covs)
            nulls = []
            for _ in range(spec["n_permutations"]):
                shuffled = covs.copy()
                shuffled[name] = rng.permutation(shuffled[name].to_numpy())
                nulls.append(stats_for(shuffled))
            for i, (metric, value) in enumerate([("partial_r2_median", partial), ("r2_alone_median", alone)]):
                null_mean, null_sd = bm.null_summary([n[i] for n in nulls])
                rows.append(_row(**base, label=name, metric=metric, value=value, null_mean=null_mean,
                                 null_sd=null_sd, n_donors=len(covs),
                                 n_levels=np.nan if usable[name].numeric else covs[name].nunique()))
    for name in targets:
        if name not in table.columns:
            continue
        ok = table[name].notna().to_numpy()
        if ok.sum() < 10:
            continue
        accuracy, chance, n_levels = bm.classify(Y[ok], table.loc[ok, name].to_numpy(), rng)
        if np.isfinite(accuracy):
            rows.append(_row(**base, label=name, metric="balanced_accuracy", value=accuracy,
                             null_mean=chance, null_sd=np.nan, n_donors=int(ok.sum()), n_levels=n_levels))
    return rows


def _paired_rows(gs: str, celltype: str, table: pd.DataFrame, Y: np.ndarray, spec) -> List[Dict[str, Any]]:
    """paired_decomposition for each pair of datasets, over the donors in both."""
    rows = []
    for a, b in combinations(list(dict.fromkeys(table["dataset"])), 2):
        rows_a = table.index[table["dataset"] == a]
        rows_b = table.index[table["dataset"] == b]
        donor_a = pd.Series(rows_a, index=table.loc[rows_a, "donor"])
        donor_b = pd.Series(rows_b, index=table.loc[rows_b, "donor"])
        both = sorted(set(donor_a.index) & set(donor_b.index))
        if len(both) < spec["min_donors"]:
            continue
        values = bm.paired_decomposition(Y[donor_a[both].to_numpy()], Y[donor_b[both].to_numpy()])
        rows += [_row(gene_set=gs, celltype=celltype, label=f"{a} vs {b}", metric=metric,
                      value=value, n_donors=len(both)) for metric, value in values.items()]
    return rows


def _mean_rows(rows: List[Dict[str, Any]], weight: str) -> List[Dict[str, Any]]:
    """Per gene set, label and metric: the mean over cell types, weighted by `weight`."""
    df = pd.DataFrame([r for r in rows if r["celltype"] != "all"])
    if df.empty:
        return []
    out = []
    for (gs, label, metric), part in df.groupby(["gene_set", "label", "metric"], sort=False):
        w = part[weight].to_numpy(dtype=np.float64)
        row = dict(gene_set=gs, celltype="mean", label=label, metric=metric,
                   value=float(np.average(part["value"], weights=w)), n_celltypes=len(part))
        if "n_genes" in part:
            row["n_genes"] = part["n_genes"].iloc[0]
        if "null_mean" in part and part["null_mean"].notna().all():
            row["null_mean"] = float(np.average(part["null_mean"], weights=w))
            row["excess"] = row["value"] - row["null_mean"]
        out.append(row)
    return out


# ---------------------------------------------------------------- run


def _celltypes_to_analyse(datasets: List[Dataset], min_donors: int) -> List[str]:
    keep = None
    for ds in datasets:
        donors = ds.sample.groupby("celltype")["donor"].nunique()
        ok = set(donors.index[donors >= min_donors])
        keep = ok if keep is None else keep & ok
    return sorted(keep)


def _cells_table(datasets: List[Dataset], celltypes: List[str]) -> pd.DataFrame:
    rows = []
    for ds in datasets:
        rows.append({"dataset": ds.name, "celltype": "all", **ds.stats,
                     "cells_sampled": len(ds.sample), "donors_sampled": ds.sample["donor"].nunique()})
        for t in celltypes:
            cells, sample = ds.cells[ds.cells["celltype"] == t], ds.sample[ds.sample["celltype"] == t]
            rows.append({"dataset": ds.name, "celltype": t, "cells_kept": len(cells),
                         "donors_kept": cells["donor"].nunique(), "cells_sampled": len(sample),
                         "donors_sampled": sample["donor"].nunique()})
    table = pd.DataFrame(rows)
    counts = table.columns.drop(["dataset", "celltype"])
    table[counts] = table[counts].astype("Int64")
    return table


def _donor_table(ds: Dataset, columns: List[str]) -> pd.DataFrame:
    """One row per donor: the most common value of each column."""
    table = pd.DataFrame({c: _donor_majority(ds.cells, c) for c in columns}) if columns else \
        pd.DataFrame(index=pd.Index(sorted(ds.cells["donor"].unique()), name="donor"))
    table.index.name = "donor"
    return table.reset_index().assign(dataset=ds.name)


def _confounding(donors: pd.DataFrame, groupers: List[str], covariates: Dict[str, Label]) -> pd.DataFrame:
    rows = []
    for g in groupers:
        for name, label in covariates.items():
            pair = donors[[g, label.column]].dropna()
            if pair[g].nunique() < 2 or len(pair) < 3:
                continue
            if label.numeric:
                value, p = bm.eta_squared(pair[label.column].to_numpy(np.float64), pair[g].to_numpy())
                metric = "eta_squared"
            else:
                value, p = bm.cramers_v(pair[g].to_numpy(), pair[label.column].to_numpy())
                metric = "cramers_v"
            rows.append({"grouping": g, "covariate": name, "metric": metric, "value": value,
                         "p_value": p, "n_donors": len(pair), "n_groups": pair[g].nunique()})
    return pd.DataFrame(rows)


def run_batch_effects(spec_path: Union[str, Path], overwrite: bool = False) -> Path:
    spec = load_spec(spec_path)
    out = Path(spec["output_dir"])
    if (out / "summary.tsv").exists() and not overwrite:
        raise FileExistsError(f"{out}/summary.tsv exists (pass --overwrite to replace it)")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(spec["seed"])
    logger.info(f"Batch effects '{spec['name']}' -> {out}")

    datasets = [prepare_dataset(name, cfg, spec) for name, cfg in spec["datasets"].items()]
    single = len(datasets) == 1
    if spec["match_donors"] and not single:
        match_donors(datasets)
    elif not single:
        for ds in datasets:
            ds.cells["donor"] = ds.name + ":" + ds.cells["donor"].astype(str)
    for ds in datasets:
        subsample(ds, spec["max_cells_per_donor"], rng)
    gene_sets = resolve_gene_sets(spec, datasets)
    if not gene_sets:
        raise ValueError("No usable gene sets")
    for ds in datasets:
        logger.info(f"[{ds.name}] streaming counts for {len(ds.union)} genes")
        stream_counts(ds, spec["chunk_size"])

    # Labels: batch columns (one dataset) or the dataset (several), covariates
    # present in every dataset, donor and cell type.
    covariate_names = [k for k in datasets[0].cfg["covariates"]
                       if all(k in ds.cfg["covariates"] for ds in datasets)]
    covariates = {k: Label(k, numeric=bool(datasets[0].cfg["covariates"][k].get("numeric")),
                           null="between") for k in covariate_names}
    labels: Dict[str, Label] = {"celltype": Label("celltype"), "donor": Label("donor")}
    if single:
        groupers = [f"batch:{c}" for c in datasets[0].cfg["batch"]]
        labels.update({g: Label(g, null="between", mixing=True) for g in groupers})
    else:
        groupers = ["dataset"]
        labels["dataset"] = Label("dataset", null="within" if spec["match_donors"] else "between",
                                  mixing=True)
    labels.update(covariates)

    sample = pd.concat([ds.sample.assign(dataset=ds.name) for ds in datasets], ignore_index=True)
    celltypes = _celltypes_to_analyse(datasets, spec["min_donors"])
    logger.info(f"Cell types with >= {spec['min_donors']} donors in every dataset: {celltypes}")
    _cells_table(datasets, celltypes).to_csv(out / "cells.tsv", sep="\t", index=False)

    # One row per dataset x donor. With matched donors every donor is in every
    # dataset, so the dataset cannot be confounded with a donor covariate.
    donor_cols = [c for c in sample.columns if c.startswith("batch:")] + covariate_names
    donors = pd.concat([_donor_table(ds, donor_cols) for ds in datasets], ignore_index=True)
    matched = spec["match_donors"] and not single
    confounding = pd.DataFrame() if matched else _confounding(donors, groupers, covariates)
    confounding.to_csv(out / "confounding.tsv", sep="\t", index=False, float_format="%.4g")

    cell_rows: List[Dict[str, Any]] = []
    donor_rows: List[Dict[str, Any]] = []
    for gs, positions in gene_sets.items():
        t0 = time.time()
        n_genes = len(next(iter(positions.values())))
        X = np.vstack([
            log_normalize(ds.counts[:, positions[ds.name]], ds.totals, spec["target_sum"]).toarray()
            for ds in datasets
        ])
        rows = cell_metrics(gs, X, sample, labels, celltypes, spec, rng)
        cell_rows += [{**r, "n_genes": n_genes} for r in rows]
        del X

        pb = pd.concat([ds.groups.assign(dataset=ds.name) for ds in datasets], ignore_index=True)
        Y_all = np.vstack([_log_cpm(ds.pseudobulk[:, positions[ds.name]], ds.groups["library_size"].to_numpy())
                           for ds in datasets])
        pb = pb.merge(donors, on=["donor", "dataset"], how="left", validate="many_to_one")
        # Matched donors: the paired decomposition replaces the model (whose
        # rows would not be independent), and only the dataset is predicted.
        terms = {} if matched else {g: Label(g) for g in groupers} | covariates
        targets = groupers + ([] if matched else [k for k, c in covariates.items() if not c.numeric])
        for t in celltypes:
            sel = ((pb["celltype"] == t) & (pb["n_cells"] >= spec["min_cells_per_pseudobulk"])).to_numpy()
            table, Y = pb[sel].reset_index(drop=True), Y_all[sel]
            if table["donor"].nunique() < spec["min_donors"]:
                continue
            if matched:
                donor_rows += [{**r, "n_genes": n_genes} for r in _paired_rows(gs, t, table, Y, spec)]
            rows = donor_metrics(gs, table, Y, terms, targets, t, spec, rng)
            donor_rows += [{**r, "n_genes": n_genes} for r in rows]
        logger.info(f"Gene set {gs} ({n_genes} genes) done in {time.time() - t0:.0f} s")

    cell_rows += _mean_rows(cell_rows, "n_cells")
    donor_rows += _mean_rows([r for r in donor_rows if "n_donors" in r], "n_donors")
    cell_df, donor_df = pd.DataFrame(cell_rows), pd.DataFrame(donor_rows)
    cell_df.to_csv(out / "cell_metrics.tsv", sep="\t", index=False, float_format="%.4g")
    donor_df.to_csv(out / "donor_metrics.tsv", sep="\t", index=False, float_format="%.4g")
    summary = pd.concat([
        cell_df[cell_df["celltype"].isin(["all", "mean"])].assign(level="cells"),
        donor_df[donor_df["celltype"] == "mean"].assign(level="donors") if len(donor_df) else None,
    ], ignore_index=True)
    columns = ["gene_set", "level", "celltype", "label", "metric", "value", "null_mean", "excess"]
    summary[[c for c in columns if c in summary.columns]].to_csv(
        out / "summary.tsv", sep="\t", index=False, float_format="%.4g")

    record = provenance_record("pipeline batch-effects", {"spec": str(spec_path)})
    record["spec"] = spec
    record["gene_set_sizes"] = {gs: len(next(iter(p.values()))) for gs, p in gene_sets.items()}
    (out / "provenance.json").write_text(json.dumps(record, indent=2, default=str))
    logger.info(f"Wrote batch-effect tables to {out}")
    return out
