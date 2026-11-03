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
plus the paired dataset-vs-donor decomposition. With several datasets the
cell-level metrics use a balanced part of the sample (balance_datasets): within
each cell type the same number of cells from every dataset, and from each of
their donors, so the larger dataset does not dominate the neighbourhoods.

With `donor_models` (named lists of covariates; `disease` is the disease
label), the donor effect: how much of the variation between donors each
covariate alone, the disease, and each model explains. With one dataset, at
cell level (the share of the PC variance within each cell type, also as a
share of donor's), and at donor level on two random halves of each donor's
cells, whose shared variance is the reproducible part of the donor effect
(cell sampling noise drops out). With matched donors, at donor level between
the datasets: the part of the donor effect both reproduce.

With `exclude_genes` (a TXT list, e.g. the X and Y genes from `genes
on-chromosomes`), those genes are left out of every gene set; `extends` reuses
another spec, so such a variant only names what differs.

A batch that holds whole donors also carries their differences, so every
batch, dataset and covariate metric comes with a null (`null_type`): the same
label shuffled between donors. A label that splits donors is shuffled within
each donor instead: a batch of pooled libraries (which has no donor-level
value, so it stays out of the donor-level metrics), and the dataset of matched
donors (whose donor-level metric is the paired decomposition, plus the
classifier). `excess` = value - null mean.

Outputs in the spec's output_dir. All are aggregates: rows are gene sets x
cell types x labels, never donors; cell types with fewer than min_donors
donors (in any dataset) are left out.
    cells.tsv          cells and donors per dataset and cell type (kept,
                       sampled, and balanced for the cell-level metrics)
    cell_metrics.tsv   PC regression, neighbourhood mixing, kBET-style
                       rejection, silhouettes
    donor_metrics.tsv  pseudobulk partial R², classifier accuracy, paired
                       decomposition
    confounding.tsv    association of batch (or dataset) with the covariates
    disease.tsv        with a `disease` label per dataset: per-gene disease
                       tests (adjusted for `disease_adjust`), classifier
                       accuracy within each dataset, and between datasets
                       the correlation of the disease effects and the
                       accuracy of a classifier fitted on the other dataset
    disease_covariates.tsv  per dataset, how each covariate goes with the
                       disease, and the accuracy of predicting the disease
                       from the covariates alone
    donor_effect.tsv   with donor_models: the donor effect (above)
    summary.tsv        the rows of the metric tables for all cells / the
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
from ..utils.obs_summary import obs_column_names, read_obs_column
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
    "balance_datasets": True,
    "disease_adjust": [],
    "donor_models": {},
    "exclude_genes": None,
    "chunk_size": 50000,
}
REQUIRED = ("output_dir", "datasets", "gene_sets")
CONFOUNDING_COLUMNS = ["grouping", "covariate", "metric", "value", "p_value", "n_donors", "n_groups"]
DATASET_KEYS = {"counts", "layer", "gene_map", "donor", "match_on", "celltype", "batch", "covariates",
                "disease"}


# ---------------------------------------------------------------- spec


def load_spec(path: Union[str, Path]) -> Dict[str, Any]:
    """The YAML spec with defaults filled in; unknown keys are an error. With
    `extends: <spec>` (relative to this one), that spec's keys are the base and
    this spec's keys replace them."""
    with open(path) as f:
        raw = yaml.safe_load(f) or {}
    if raw.get("extends"):
        with open(Path(path).parent / raw["extends"]) as f:
            base = yaml.safe_load(f) or {}
        if base.get("extends"):
            raise ValueError(f"{raw['extends']} extends another spec; only one level is supported")
        raw = {**{k: v for k, v in base.items() if k != "name"}, **raw}
    unknown = set(raw) - set(DEFAULTS) - set(REQUIRED) - {"name", "extends"}
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
    with_disease = [bool(cfg.get("disease")) for cfg in spec["datasets"].values()]
    if any(with_disease) and not all(with_disease):
        raise ValueError("'disease' must be given for every dataset or none")
    for name, cfg in spec["datasets"].items():
        missing = set(spec["disease_adjust"]) - set(cfg["covariates"])
        if missing:
            raise ValueError(f"disease_adjust names covariates that dataset {name} lacks: {sorted(missing)}")
    for model, keys in spec["donor_models"].items():
        if not keys:
            raise ValueError(f"donor_models: model {model} names no covariates")
        for name, cfg in spec["datasets"].items():
            if model in cfg["covariates"] or model in ("disease", "donor"):
                raise ValueError(f"donor_models: model {model} has the name of a label")
            missing = set(keys) - set(cfg["covariates"]) - ({"disease"} if cfg.get("disease") else set())
            if missing:
                raise ValueError(f"donor_models: model {model} names covariates that dataset {name} "
                                 f"lacks: {sorted(missing)}")
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


def _as_numeric(values, recode: Optional[Dict] = None) -> np.ndarray:
    """Numbers from obs values, optionally recoded first (e.g. 'Braak IV' -> 4);
    '90+' -> 90 (the ROSMAP age cap)."""
    text = pd.Series(_as_labels(values, recode), dtype=object)
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
    pseudobulk_half: Optional[np.ndarray] = None  # the same over a random half of the cells


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
    if cfg.get("disease"):
        columns.add(cfg["disease"]["column"])
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
        missing = sorted(set(columns) - set(obs_column_names(f)))
        if missing:
            raise ValueError(f"[{name}] obs has no column(s) {missing}")
        raw = {c: read_obs_column(f, c) for c in columns}
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
        values = _as_numeric(raw[cov["column"]], cov.get("recode")) if cov.get("numeric") else \
            _as_labels(raw[cov["column"]], cov.get("recode"))
        cells[key] = values[rows]
    if cfg.get("disease"):
        cells["disease"] = _as_labels(raw[cfg["disease"]["column"]], cfg["disease"].get("recode"))[rows]
    if cfg.get("match_on"):
        cells["match"] = _as_labels(raw[cfg["match_on"]])[rows]
    stats["cells_kept"] = len(cells)
    stats["donors_kept"] = cells["donor"].nunique()
    logger.info(f"[{name}] " + ", ".join(f"{k}={v}" for k, v in stats.items()))
    for t, count in cells["celltype"].value_counts().items():
        logger.info(f"[{name}]   {t}: {count} cells")
    # Counts only, so a recode map that matches nothing shows up here
    for column in [f"batch:{c}" for c in cfg["batch"]] + list(covariates) + \
            (["disease"] if cfg.get("disease") else []):
        with_value = cells.loc[cells[column].notna(), "donor"].nunique()
        logger.info(f"[{name}]   {column}: {with_value} of {stats['donors_kept']} donors have a value, "
                    f"{cells[column].nunique()} distinct")
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


def balance_datasets(cells: pd.DataFrame, celltypes: List[str], matched: bool, min_donors: int,
                     rng: np.random.Generator) -> np.ndarray:
    """Positions in `cells` (the pooled sample of several datasets) for the
    cell-level metrics: within each cell type, every dataset gets the same
    number of cells, m from each of its donors. With matched donors these are
    the donors with at least m cells in every dataset; otherwise the same
    number of donors per dataset, drawn from those with at least m cells. m is
    the value that keeps the most cells. The larger dataset then does not
    dominate the neighbourhoods, and since every donor has the same weight,
    shuffling the dataset between (or within) donors keeps the split."""
    keep = []
    datasets = cells["dataset"].to_numpy()
    donors = cells["donor"].to_numpy()
    for t in celltypes:
        part = np.flatnonzero(cells["celltype"].to_numpy() == t)
        counts = cells.iloc[part].groupby(["donor", "dataset"]).size().unstack("dataset", fill_value=0)
        best_cells, m = 0, 0
        for candidate in np.unique(counts.to_numpy()):
            enough = counts >= candidate
            n = int(enough.all(axis=1).sum()) if matched else int(enough.sum(axis=0).min())
            if candidate > 0 and n >= min_donors and n * candidate > best_cells:
                best_cells, m = n * candidate, int(candidate)
        if not m:
            logger.warning(f"{t}: cannot balance the datasets with >= {min_donors} donors (left out)")
            continue
        enough = counts >= m
        if matched:
            chosen = {ds: counts.index[enough.all(axis=1)] for ds in counts.columns}
        else:
            n = int(enough.sum(axis=0).min())
            chosen = {ds: rng.choice(counts.index[enough[ds]], size=n, replace=False) for ds in counts.columns}
        for ds, picked in chosen.items():
            rows = part[(datasets[part] == ds) & np.isin(donors[part], picked)]
            for members in pd.Series(rows).groupby(donors[rows]).agg(list):
                keep.append(rng.choice(members, size=m, replace=False))
        logger.info(f"{t}: {m} cells from each of {len(next(iter(chosen.values())))} donors per dataset "
                    f"for the cell-level metrics")
    return np.sort(np.concatenate(keep)) if keep else np.array([], dtype=int)


def stream_counts(ds: Dataset, chunk_size: int, halves_rng: Optional[np.random.Generator] = None) -> None:
    """One pass over the matrix: sampled cells' counts for ds.union, their totals
    over all genes, and pseudobulk sums per donor x cell type (with halves_rng,
    also over a random half of each group's cells)."""
    kept = ds.cells["row"].to_numpy()
    sampled = ds.sample["row"].to_numpy()
    group_codes, groups = pd.MultiIndex.from_arrays(
        [ds.cells["donor"], ds.cells["celltype"]]).factorize()
    n_groups = len(groups)
    pseudobulk = np.zeros((n_groups, len(ds.union)))
    library = np.zeros(n_groups)
    half = halves_rng.random(len(kept)) < 0.5 if halves_rng is not None else None
    pseudobulk_half = np.zeros_like(pseudobulk) if half is not None else None
    library_half = np.zeros(n_groups)
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
            if half is not None:
                in_half = np.flatnonzero(half[lo:hi])
                indicator = bm._indicator(group_codes[lo:hi][in_half], n_groups)
                pseudobulk_half += (indicator @ genes[in_half]).toarray()
                library_half += indicator @ block_totals[in_half]
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
    if half is not None:
        ds.pseudobulk_half = pseudobulk_half
        ds.groups["n_half"] = np.bincount(group_codes[half], minlength=n_groups)
        ds.groups["library_half"] = library_half


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


def _excluded_columns(ds: Dataset, entries: List[str]) -> set:
    """Every file column of the listed genes (Ensembl IDs via the gene map, or symbols)."""
    if sum(e.startswith("ENS") for e in entries) > len(entries) / 2:
        if not ds.cfg.get("gene_map"):
            raise ValueError(f"[{ds.name}] exclude_genes holds Ensembl IDs; set gene_map")
        from ..data.genes import load_gene_map

        mapping = load_gene_map(ds.cfg["gene_map"])
        keys = [mapping.get(name, "").split(".")[0] for name in ds.var_names]
    else:
        keys = list(ds.var_names)
    wanted = set(entries)
    return {col for col, key in enumerate(keys) if key in wanted}


def resolve_gene_sets(spec: Dict[str, Any], datasets: List[Dataset]) -> Dict[str, Dict[str, np.ndarray]]:
    """Per gene set and dataset, the file columns of its genes (aligned across
    datasets; genes missing from any dataset, or on the exclude_genes list in
    any, are left out). Sets ds.union and returns positions into it."""
    excluded = {ds.name: set() for ds in datasets}
    if spec["exclude_genes"]:
        listed = read_gene_set(spec["exclude_genes"])
        excluded = {ds.name: _excluded_columns(ds, listed) for ds in datasets}
        for ds in datasets:
            logger.info(f"[{ds.name}] exclude_genes: {len(excluded[ds.name])} of {len(listed)} listed "
                        f"genes are in the count file")
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
        kept = [e for e in shared if not any(f[e] in excluded[ds.name] for ds, f in zip(datasets, found))]
        if spec["exclude_genes"]:
            logger.info(f"Gene set {gs}: {len(shared) - len(kept)} genes on the exclude_genes list left out")
        shared = kept
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
                 celltypes: List[str], spec, rng,
                 keep_pcs: Optional[Dict[str, np.ndarray]] = None) -> List[Dict[str, Any]]:
    """keep_pcs, if given, receives the PCs of each cell type."""
    rows: List[Dict[str, Any]] = []
    scopes = [("all", np.arange(len(cells)))] + \
        [(t, np.flatnonzero(cells["celltype"].to_numpy() == t)) for t in celltypes]
    for scope, idx in scopes:
        pcs = _pca(X[idx], spec["n_pcs"], spec["seed"])
        if keep_pcs is not None and scope != "all":
            keep_pcs[scope] = pcs
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
                             null_mean=null_mean, null_sd=null_sd, null_type=label.null, n_levels=n_levels))
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
                             null_mean=null_mean, null_sd=null_sd, null_type=label.null, n_levels=n_levels))
            kbet = lambda lab, k=k: bm.knn_chi2_rejection(k, lab)
            null_mean, null_sd = _null(kbet, v, d, label, spec, rng)
            rows.append(_row(**base, label=name, metric="knn_rejection", value=kbet(v),
                             null_mean=null_mean, null_sd=null_sd, null_type=label.null, n_levels=n_levels))
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
        # A term missing for a whole dataset would drop that dataset from the model
        empty = [d for d, v in values.groupby(table["dataset"]) if v.isna().all()]
        if empty:
            logger.info(f"{gs} / {celltype}: term {name} has no values in {empty}; left out of the model")
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
                                 null_sd=null_sd, null_type="between", n_donors=len(covs),
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


def _zscore(Y: np.ndarray) -> np.ndarray:
    sd = Y.std(axis=0)
    return (Y - Y.mean(axis=0)) / np.where(sd > 0, sd, 1.0)


def disease_metrics(gs: str, table: pd.DataFrame, Y: np.ndarray, adjust: Dict[str, Label],
                    celltype: str, spec, rng) -> List[Dict[str, Any]]:
    """Disease signal in the pseudobulk of one cell type (rows: donor x dataset).

    Per dataset: the per-gene disease t statistic, adjusted for `adjust`
    (share of genes at FDR < 0.05, mean t²), and the cross-validated
    classifier accuracy. Per pair of datasets: the correlation of the t
    statistics, and the accuracy of a classifier fitted on one dataset on the
    other, standardized on the training data (as deployed) or each dataset on
    its own (after removing per-gene shifts and scales). Nulls shuffle the
    disease label between donors (of the training dataset for transfer)."""
    rows: List[Dict[str, Any]] = []
    base = dict(gene_set=gs, celltype=celltype)
    numeric = [k for k, label in adjust.items() if label.numeric]
    fits = {}
    for name in dict.fromkeys(table["dataset"]):
        sel = ((table["dataset"] == name) & table["disease"].notna()).to_numpy()
        for label in adjust.values():
            sel &= table[label.column].notna().to_numpy()
        sub = table[sel]
        counts = sub["disease"].value_counts()
        if len(counts) < 2 or counts.min() < 3:
            continue
        covs = sub[["disease"] + [label.column for label in adjust.values()]]
        covs.columns = ["disease"] + list(adjust)
        X, columns = bm.design_matrix(covs, numeric)
        j = int(columns["disease"][0])
        Ys = Y[sel]
        share, mean_t2, t = bm.association(Ys, X, j)
        nulls = []
        for _ in range(spec["n_permutations"]):
            Xp = X.copy()
            Xp[:, j] = rng.permutation(X[:, j])
            nulls.append(bm.association(Ys, Xp, j)[:2])
        extra = dict(n_donors=len(sub), n_positive=int((X[:, j] == 1).sum()), null_type="between")
        for i, (metric, value) in enumerate([("de_share_fdr05", share), ("de_mean_t2", mean_t2)]):
            null_mean, null_sd = bm.null_summary([n[i] for n in nulls])
            rows.append(_row(**base, label=name, metric=metric, value=value,
                             null_mean=null_mean, null_sd=null_sd, **extra))
        accuracy, chance, _ = bm.classify(Ys, sub["disease"].to_numpy(), rng)
        if np.isfinite(accuracy):
            rows.append(_row(**base, label=name, metric="within_accuracy", value=accuracy,
                             null_mean=chance, null_sd=np.nan, **extra))
        fits[name] = dict(Y=Ys, y=sub["disease"].to_numpy(), donors=sub["donor"].to_numpy(),
                          X=X, j=j, t=t)

    for a, b in combinations(fits, 2):
        ta, tb = fits[a]["t"], fits[b]["t"]
        if (np.isfinite(ta) & np.isfinite(tb)).sum() < 3:
            continue

        def corr(t, tb=tb):
            ok = np.isfinite(t) & np.isfinite(tb)
            return float(np.corrcoef(t[ok], tb[ok])[0, 1])

        nulls = []
        for _ in range(spec["n_permutations"]):
            Xp = fits[a]["X"].copy()
            Xp[:, fits[a]["j"]] = rng.permutation(Xp[:, fits[a]["j"]])
            nulls.append(corr(bm.coefficient_t(fits[a]["Y"], Xp, fits[a]["j"])))
        null_mean, null_sd = bm.null_summary(nulls)
        rows.append(_row(**base, label=f"{a} vs {b}", metric="t_correlation", value=corr(ta),
                         null_mean=null_mean, null_sd=null_sd, null_type="between",
                         n_donors=len(fits[a]["y"]) + len(fits[b]["y"])))

    for a, b in combinations(fits, 2):
        for src, dst in ((a, b), (b, a)):
            s, d = fits[src], fits[dst]
            for metric, (Xs, Xd) in (("transfer_accuracy", (s["Y"], d["Y"])),
                                     ("transfer_accuracy_rescaled", (_zscore(s["Y"]), _zscore(d["Y"])))):
                score = lambda y: bm.transfer_accuracy(Xs, y, Xd, d["y"], rng, s["donors"], d["donors"])
                value = score(s["y"])
                if not np.isfinite(value):
                    continue
                null_mean, null_sd = bm.null_summary([score(rng.permutation(s["y"])) for _ in range(5)])
                rows.append(_row(**base, label=f"{src} -> {dst}", metric=metric, value=value,
                                 null_mean=null_mean, null_sd=null_sd, null_type="between",
                                 n_donors=len(d["y"]), n_positive=int((d["X"][:, d["j"]] == 1).sum())))
    return rows


def _effect_design(covs: pd.DataFrame, keys: List[str], numeric: List[str]):
    """Donors (rows of covs) with every covariate in keys, and their design
    matrix (None if the covariates do not vary)."""
    ok = covs[keys].notna().all(axis=1).to_numpy()
    D, _ = bm.design_matrix(covs.loc[ok, keys], [k for k in keys if k in numeric])
    return ok, (D if D.shape[1] > 1 else None)


def donor_effect_cells(gs: str, pcs_by_type: Dict[str, np.ndarray], cells: pd.DataFrame,
                       covs: pd.DataFrame, terms: Dict[str, List[str]], numeric: List[str],
                       spec, rng) -> List[Dict[str, Any]]:
    """Cell level, per cell type: the share of the PC variance that donor
    explains, and per term the share its covariates explain (pc_share), also
    as a share of the donor's (pc_share_of_donor; donor-level covariates
    cannot explain more than donor). covs: one row per donor, indexed by
    donor. Nulls shuffle the covariates between donors."""
    rows: List[Dict[str, Any]] = []
    celltype = cells["celltype"].to_numpy()
    for t, pcs in pcs_by_type.items():
        donors = cells["donor"].to_numpy()[celltype == t]
        base = dict(gene_set=gs, celltype=t)
        rows.append(_row(**base, label="donor", metric="pc_share", value=bm.pc_regression(pcs, donors),
                         n_donors=len(set(donors)), n_cells=len(donors)))
        sub = covs.loc[sorted(set(donors))]
        for name, keys in terms.items():
            ok, D = _effect_design(sub, keys, numeric)
            if D is None or ok.sum() < spec["min_donors"]:
                continue
            kept = sub.index[ok]
            mask = np.isin(donors, kept)
            pos = pd.Index(kept).get_indexer(donors[mask])
            p = pcs[mask]
            donor_share = bm.pc_regression(p, donors[mask])
            value = bm.pc_regression_design(p, D[pos])
            nulls = [bm.pc_regression_design(p, D[rng.permutation(len(D))][pos])
                     for _ in range(spec["n_permutations"])]
            null_mean, null_sd = bm.null_summary(nulls)
            extra = dict(null_type="between", n_donors=int(ok.sum()), n_cells=int(mask.sum()),
                         covariates="+".join(keys))
            rows.append(_row(**base, label=name, metric="pc_share", value=value, null_mean=null_mean,
                             null_sd=null_sd, **extra))
            if donor_share > 0:
                rows.append(_row(**base, label=name, metric="pc_share_of_donor", value=value / donor_share,
                                 null_mean=null_mean / donor_share, null_sd=null_sd / donor_share, **extra))
    return rows


def donor_effect_shared(gs: str, celltype: str, A: np.ndarray, B: np.ndarray, covs: pd.DataFrame,
                        terms: Dict[str, List[str]], numeric: List[str], metric: str,
                        spec, rng) -> List[Dict[str, Any]]:
    """Donor level, the same donors measured twice (rows of A and B, aligned
    with covs): two random halves of their cells (metric
    'split_half_reliability', for the pseudobulk of both halves) or two
    datasets ('cross_dataset_reliability'). The share of the between-donor
    variance the two measurements share, and per term the share of that
    shared variance its covariates explain (reproducible_share). Nulls
    shuffle the covariates between donors."""
    base = dict(gene_set=gs, celltype=celltype)
    reliability = bm.shared_reliability(A, B, spearman_brown=metric == "split_half_reliability")
    rows = [_row(**base, label="donor", metric=metric, value=reliability, n_donors=len(A))]
    for name, keys in terms.items():
        ok, D = _effect_design(covs, keys, numeric)
        if D is None or ok.sum() < spec["min_donors"]:
            continue
        a, b = A[ok], B[ok]
        value = bm.shared_variance_explained(a, b, D)
        nulls = [bm.shared_variance_explained(a, b, D[rng.permutation(len(D))])
                 for _ in range(spec["n_permutations"])]
        null_mean, null_sd = bm.null_summary(nulls)
        rows.append(_row(**base, label=name, metric="reproducible_share", value=value, null_mean=null_mean,
                         null_sd=null_sd, null_type="between", n_donors=int(ok.sum()),
                         covariates="+".join(keys)))
    return rows


def _halves(ds: Dataset, columns: np.ndarray, celltype: str, min_cells: int):
    """log-CPM of the two random halves of each donor's cells of one cell type
    (both with >= min_cells cells), and the donors."""
    g = ds.groups
    other = (g["n_cells"] - g["n_half"]).to_numpy()
    sel = ((g["celltype"] == celltype) & (g["n_half"] >= min_cells)).to_numpy() & (other >= min_cells)
    half, full = ds.pseudobulk_half[sel][:, columns], ds.pseudobulk[sel][:, columns]
    library_half = g.loc[sel, "library_half"].to_numpy()
    A = _log_cpm(half, library_half)
    B = _log_cpm(full - half, g.loc[sel, "library_size"].to_numpy() - library_half)
    return A, B, g.loc[sel, "donor"].to_numpy()


DISEASE_COVARIATE_COLUMNS = ["dataset", "term", "covariates", "metric", "value", "p_value", "null_mean",
                             "excess", "n_donors"]


def disease_covariates(donors: pd.DataFrame, covariates: Dict[str, Label], derived: set,
                       terms: Dict[str, List[str]], numeric: List[str], rng) -> pd.DataFrame:
    """Per dataset, at donor level: how each covariate goes with the disease
    label (Cramér's V, or eta² for a numeric one, with p-values), and the
    cross-validated balanced accuracy of predicting the disease from the
    covariates of each term alone: what a model could get from them as a
    shortcut. Covariates read from the disease column itself (`derived`) are
    left out."""
    rows = []
    for name, part in donors.groupby("dataset", sort=False):
        part = part[part["disease"].notna()]
        if part["disease"].nunique() < 2:
            continue
        for key, label in covariates.items():
            pair = part[["disease", key]].dropna()
            if key in derived or pair[key].nunique() < 2 or len(pair) < 3:
                continue
            if label.numeric:
                value, p = bm.eta_squared(pair[key].to_numpy(np.float64), pair["disease"].to_numpy())
            else:
                value, p = bm.cramers_v(pair["disease"].to_numpy(), pair[key].to_numpy())
            rows.append(dict(dataset=name, term=key, covariates=key,
                             metric="eta_squared" if label.numeric else "cramers_v",
                             value=value, p_value=p, n_donors=len(pair)))
        seen = set()
        for term, keys in terms.items():
            keys = [k for k in keys if k != "disease" and k not in derived]
            sub = part[["disease"] + keys].dropna()
            if not keys or tuple(keys) in seen or len(sub) < 10:
                continue
            seen.add(tuple(keys))
            X, _ = bm.design_matrix(sub[keys], [k for k in keys if k in numeric])
            if X.shape[1] < 2:
                continue
            accuracy, chance, _ = bm.classify(X[:, 1:], sub["disease"].to_numpy(), rng)
            if np.isfinite(accuracy):
                rows.append(dict(dataset=name, term=term, covariates="+".join(keys),
                                 metric="balanced_accuracy", value=accuracy, null_mean=chance,
                                 excess=accuracy - chance, n_donors=len(sub)))
    return pd.DataFrame(rows, columns=DISEASE_COVARIATE_COLUMNS)


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


def _cells_table(datasets: List[Dataset], celltypes: List[str],
                 balanced: Optional[pd.DataFrame] = None) -> pd.DataFrame:
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
    if balanced is not None:
        for i, row in table.iterrows():
            part = balanced[balanced["dataset"] == row["dataset"]]
            if row["celltype"] != "all":
                part = part[part["celltype"] == row["celltype"]]
            table.loc[i, "cells_balanced"] = len(part)
            table.loc[i, "donors_balanced"] = part["donor"].nunique()
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
    return pd.DataFrame(rows, columns=CONFOUNDING_COLUMNS)


def run_batch_effects(spec_path: Union[str, Path], overwrite: bool = False) -> Path:
    spec = load_spec(spec_path)
    out = Path(spec["output_dir"])
    if (out / "summary.tsv").exists() and not overwrite:
        raise FileExistsError(f"{out}/summary.tsv exists (pass --overwrite to replace it)")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(spec["seed"])
    # A separate stream, so the other tables do not depend on the disease metrics
    disease_rng = np.random.default_rng([spec["seed"], 1])
    logger.info(f"Batch effects '{spec['name']}' -> {out}")

    datasets = [prepare_dataset(name, cfg, spec) for name, cfg in spec["datasets"].items()]
    single = len(datasets) == 1
    matched = spec["match_donors"] and not single
    # The donor effect needs each donor measured twice: two random halves of
    # its cells (one dataset) or two datasets (matched donors). Own stream, so
    # the other tables do not depend on it.
    effect = bool(spec["donor_models"]) and (single or matched)
    effect_rng = np.random.default_rng([spec["seed"], 2])
    if spec["donor_models"] and not effect:
        logger.warning("donor_models needs one dataset or matched donors; skipping the donor effect")
    if matched:
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
        stream_counts(ds, spec["chunk_size"], effect_rng if effect and single else None)

    # Labels: batch columns (one dataset) or the dataset (several), covariates
    # present in every dataset, donor and cell type.
    covariate_names = [k for k in datasets[0].cfg["covariates"]
                       if all(k in ds.cfg["covariates"] for ds in datasets)]
    covariates = {k: Label(k, numeric=bool(datasets[0].cfg["covariates"][k].get("numeric")),
                           null="between") for k in covariate_names}
    sample = pd.concat([ds.sample.assign(dataset=ds.name) for ds in datasets], ignore_index=True)
    labels: Dict[str, Label] = {"celltype": Label("celltype"), "donor": Label("donor")}
    if single:
        # A batch that holds whole donors is compared with itself shuffled
        # between donors. One that splits donors (pooled libraries) is
        # compared with itself shuffled within each donor, and has no
        # donor-level value, so it stays out of the donor-level metrics.
        groupers = []
        for col in [f"batch:{c}" for c in datasets[0].cfg["batch"]]:
            labelled = sample.dropna(subset=[col])
            several = float(labelled.groupby("donor")[col].nunique().gt(1).mean()) if len(labelled) else 0.0
            within = several > 0.5
            labels[col] = Label(col, null="within" if within else "between", mixing=True)
            if not within:
                groupers.append(col)
            logger.info(f"{col}: {100 * several:.0f}% of donors have cells in several batches; "
                        f"null shuffles it {'within' if within else 'between'} donors")
    else:
        groupers = ["dataset"]
        labels["dataset"] = Label("dataset", null="within" if spec["match_donors"] else "between",
                                  mixing=True)
    labels.update(covariates)
    has_disease = bool(datasets[0].cfg.get("disease"))
    if has_disease:
        levels = set().union(*(set(ds.cells["disease"].dropna()) for ds in datasets))
        if len(levels) != 2:
            raise ValueError(f"The disease label needs exactly 2 values across the datasets, not {len(levels)}")
    adjust = {k: covariates[k] for k in spec["disease_adjust"]}
    # Donor effect and disease covariates: every covariate alone (except those
    # read from the disease column, which the disease label stands for), the
    # disease, and the spec's joint models
    derived = {k for k in covariate_names for ds in datasets
               if ds.cfg.get("disease") and ds.cfg["covariates"][k]["column"] == ds.cfg["disease"]["column"]}
    effect_terms = {k: [k] for k in covariate_names if k not in derived}
    if has_disease:
        effect_terms["disease"] = ["disease"]
    effect_terms.update(spec["donor_models"])
    numeric = [k for k, c in covariates.items() if c.numeric]
    disease_cov_rng = np.random.default_rng([spec["seed"], 3])

    celltypes = _celltypes_to_analyse(datasets, spec["min_donors"])
    logger.info(f"Cell types with >= {spec['min_donors']} donors in every dataset: {celltypes}")
    # Several datasets: the cell-level metrics use a sample with equal cells
    # per dataset (and per donor) within each cell type.
    cell_idx = balance_datasets(sample, celltypes, matched, spec["min_donors"], rng) \
        if not single and spec["balance_datasets"] else None
    cell_sample = sample if cell_idx is None else sample.iloc[cell_idx].reset_index(drop=True)
    cell_celltypes = [t for t in celltypes if (cell_sample["celltype"] == t).any()]
    _cells_table(datasets, celltypes, None if cell_idx is None else cell_sample).to_csv(
        out / "cells.tsv", sep="\t", index=False)

    # One row per dataset x donor. With matched donors every donor is in every
    # dataset, so the dataset cannot be confounded with a donor covariate.
    donor_cols = [g for g in groupers if g.startswith("batch:")] + covariate_names + \
        (["disease"] if has_disease else [])
    donors = pd.concat([_donor_table(ds, donor_cols) for ds in datasets], ignore_index=True)
    confounding = pd.DataFrame(columns=CONFOUNDING_COLUMNS) if matched else \
        _confounding(donors, groupers, covariates)
    confounding.to_csv(out / "confounding.tsv", sep="\t", index=False, float_format="%.4g")
    if has_disease:
        disease_covariates(donors, covariates, derived, effect_terms, numeric, disease_cov_rng).to_csv(
            out / "disease_covariates.tsv", sep="\t", index=False, float_format="%.4g")
    # Covariates per donor (with matched donors, as the first dataset has them)
    covs = donors.drop_duplicates("donor").set_index("donor")

    cell_rows: List[Dict[str, Any]] = []
    donor_rows: List[Dict[str, Any]] = []
    disease_rows: List[Dict[str, Any]] = []
    effect_rows: List[Dict[str, Any]] = []
    for gs, positions in gene_sets.items():
        t0 = time.time()
        n_genes = len(next(iter(positions.values())))
        X = np.vstack([
            log_normalize(ds.counts[:, positions[ds.name]], ds.totals, spec["target_sum"]).toarray()
            for ds in datasets
        ])
        if cell_idx is not None:
            X = X[cell_idx]
        pcs_by_type = {} if effect and single else None
        rows = cell_metrics(gs, X, cell_sample, labels, cell_celltypes, spec, rng, keep_pcs=pcs_by_type)
        cell_rows += [{**r, "n_genes": n_genes} for r in rows]
        del X
        if pcs_by_type:
            rows = donor_effect_cells(gs, pcs_by_type, cell_sample, covs, effect_terms, numeric, spec,
                                      effect_rng)
            effect_rows += [{**r, "n_genes": n_genes} for r in rows]

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
            if has_disease:
                rows = disease_metrics(gs, table, Y, adjust, t, spec, disease_rng)
                disease_rows += [{**r, "n_genes": n_genes} for r in rows]
            if effect:
                if single:
                    A, B, ids = _halves(datasets[0], positions[datasets[0].name], t,
                                        spec["min_cells_per_pseudobulk"])
                    metric = "split_half_reliability"
                else:  # the donors in both of the first two datasets
                    row_of = {name: pd.Series(part.index, index=part["donor"])
                              for name, part in table.groupby("dataset", sort=False)}
                    row_a, row_b = (row_of.get(ds.name, pd.Series(dtype=int)) for ds in datasets[:2])
                    ids = np.array(sorted(set(row_a.index) & set(row_b.index)), dtype=object)
                    A, B = Y[row_a[ids].to_numpy()], Y[row_b[ids].to_numpy()]
                    metric = "cross_dataset_reliability"
                if len(ids) >= spec["min_donors"]:
                    rows = donor_effect_shared(gs, t, A, B, covs.loc[ids], effect_terms, numeric, metric,
                                               spec, effect_rng)
                    effect_rows += [{**r, "n_genes": n_genes} for r in rows]
        logger.info(f"Gene set {gs} ({n_genes} genes) done in {time.time() - t0:.0f} s")

    cell_rows += _mean_rows(cell_rows, "n_cells")
    donor_rows += _mean_rows([r for r in donor_rows if "n_donors" in r], "n_donors")
    cell_df, donor_df = pd.DataFrame(cell_rows), pd.DataFrame(donor_rows)
    cell_df.to_csv(out / "cell_metrics.tsv", sep="\t", index=False, float_format="%.4g")
    donor_df.to_csv(out / "donor_metrics.tsv", sep="\t", index=False, float_format="%.4g")
    disease_rows += _mean_rows(disease_rows, "n_donors")
    disease_df = pd.DataFrame(disease_rows)
    if has_disease:
        disease_df.to_csv(out / "disease.tsv", sep="\t", index=False, float_format="%.4g")
    effect_rows += _mean_rows(effect_rows, "n_donors")
    effect_df = pd.DataFrame(effect_rows)
    if effect:
        effect_df.to_csv(out / "donor_effect.tsv", sep="\t", index=False, float_format="%.4g")
    summary = pd.concat([
        cell_df[cell_df["celltype"].isin(["all", "mean"])].assign(level="cells"),
        donor_df[donor_df["celltype"] == "mean"].assign(level="donors") if len(donor_df) else None,
        disease_df[disease_df["celltype"] == "mean"].assign(level="disease") if len(disease_df) else None,
        effect_df[effect_df["celltype"] == "mean"].assign(level="donor_effect") if len(effect_df) else None,
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
