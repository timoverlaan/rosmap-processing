"""Pooled pseudobulk over the v2 datasets, with the dataset as a covariate.

Driven by a YAML spec (config/pooled_pseudobulk.yaml) and the global donor
table (pipelines/donor_table.py), which links each dataset's donors to
persons and holds the harmonized covariates.

1. Per dataset, the cell filters of `pipeline batch-effects` (genes per cell,
   cells without a donor, the broad cell type rules), then one pass over the
   count matrix: pseudobulk sums per donor x cell type over the genes every
   dataset has (Ensembl IDs via the gene maps), and the totals over all genes.
2. pseudobulk.h5ad: one row per dataset x donor x cell type (X = summed
   counts; obs = dataset, donor, person, cell type, nuclei, library size and
   the donor table's covariates). Restricted: stays on the cluster.
3. The batch of a sample is its dataset, or, for the datasets in
   batch_split, the dataset and a donor table column (PsychAD: dataset x
   brain bank).
4. Per broad cell type, on the samples with at least min_cells nuclei and
   the genes with a mean CPM of at least min_cpm in every dataset, log2-CPM
   and these ways to remove the batch effect:
     none     as measured
     center   each batch's mean per gene removed
     regress  per gene, a model of batch + the `protect` covariates (sex,
              age, dementia); only the batch terms are removed, so mean
              differences that the covariates explain are kept
     regress_dataset  (only with batch_split) regress with the dataset as
              the batch, for comparison
     anchor   per gene, the batch offsets that best explain the differences
              between the same persons measured in two batches (purely
              technical); batches without such persons are centered
     ruv      RUV-III (as ruv::RUVIII): k unwanted factors alpha from the
              replicate residuals (each person's samples minus their mean)
              and, for the batches without such persons, pseudo-replicates
              (RUV-III-PRPS: per batch, the mean of the samples of one
              dementia x sex x age band group, linked across batches;
              weighted by sqrt(samples)); per sample W from the control
              genes, Y - W alpha. k (ruv_k) and the control genes (all, or
              empirical: the half least associated with the disease in the
              reference model) are chosen per cell type by cross-fitted
              replicate top-1 (ruv_scan.tsv)
   and how each does:
     mixing           share of the PC variance that the dataset, the batch,
                      and within a split dataset its split column explain
                      (pc_r2), and the balanced accuracy of a classifier of
                      each (classifier_accuracy)
     replicates       for persons measured in two datasets: how often their
                      profile in one dataset is the closest (Pearson, over
                      genes, of the deviation from the pooled mean) to their
                      own profile in the other, of all that dataset's samples
                      (top1, chance = 1 / candidates), the median rank as a
                      share of the candidates, and the mean correlation with
                      their own vs other persons' profiles. The anchor
                      offsets are cross-fitted: a person's offsets are
                      estimated without their pairs (n_folds folds); so are
                      the ruv factors (a fold's persons are singletons in M
                      and left out of the pseudo-replicates)
     disease         the reference model: one sample per person (the one
                      with the most nuclei), log2-CPM ~ disease + adjust +
                      batch (de_pooled.tsv); per dataset alone (with its
                      batches); and per correction without the batch term.
                      Share of genes at FDR < 0.05 and the correlation of the
                      t statistics with the reference model's
     offsets          per batch, the size of each method's offsets and how
                      the center, regress and ruv (the batch mean of what it
                      removes) offsets correlate with the anchor (technical)
                      ones
     ruv_scan         per k and control gene set: cross-fitted replicate
                      top-1 and rank share, pc_r2 and disease metrics; the
                      chosen one
     ruv_diagnostic   the within-person differences of the persons in two
                      batches: the shares of their sum of squares that batch
                      offsets, the top k components of the remainder (and the
                      same for a permuted remainder) and the rest explain
     ruv_factors      the top remainder components (pair scores) and the
                      ruv W (per sample) against sample QC (nuclei, library
                      size, UMIs per nucleus): correlations only
5. With umap: true, per cell type a UMAP of the PCs of each of umap_methods,
   coloured by batch, dementia, sex and age, and with the same person's
   samples in two datasets joined by a line (umap/<cell type>.png), and all
   cell types coloured by batch (umap/overview.png). No donor IDs.

Aggregate outputs in output_dir: samples.tsv, genes.tsv, metrics.tsv,
disease.tsv, de_pooled.tsv (per gene), offsets.tsv, ruv_scan.tsv,
ruv_diagnostic.tsv, ruv_factors.tsv, umap/, provenance.json. With the spec
key pseudobulk, an existing pseudobulk.h5ad is read from there (and never
written), so another output_dir can reuse it.

Run: python -m rosmap_processing.pipelines.pooled_pseudobulk config/pooled_pseudobulk.yaml
"""

import argparse
import json
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy import stats

from ..core import batch_metrics as bm
from ..utils.constants import DEFAULT_PCA_COMPONENTS, MIN_GENES_PER_CELL
from ..utils.logging import get_logger
from ..utils.provenance import provenance_record
from .batch_effects import _columns_for, _pca, prepare_dataset, stream_counts
from .umap_plots import CATEGORICAL, INK, MISSING_COLOR, OTHER_COLOR, SEQUENTIAL, UMAP_MIN_DIST

logger = get_logger(__name__)

DEFAULTS: Dict[str, Any] = {
    "seed": 0,
    "min_genes": MIN_GENES_PER_CELL,
    "min_cells": 10,
    "min_cpm": 1.0,
    "min_pairs": 5,
    "n_pcs": DEFAULT_PCA_COMPONENTS,
    "n_folds": 5,
    "chunk_size": 50000,
    "protect": ["sex", "age_capped", "dementia"],
    "disease": "dementia",
    "disease_positive": "yes",
    "disease_adjust": ["sex", "age_capped"],
    "batch_split": {},       # dataset -> donor table column that splits its batch
    "umap": True,
    "umap_methods": ["none", "regress_dataset", "regress"],
    "umap_neighbors": 15,
    "pseudobulk": None,      # an existing pseudobulk.h5ad to read (never written) instead of output_dir's
    # RUV-III (method ruv): numbers of unwanted factors and control gene sets
    # to try (in order of preference), chosen per cell type by cross-fitted
    # replicate top-1; of those within ruv_tolerance of the best, the first
    # control set and the most factors (the replicates cannot see the batches
    # without bridges, which need factors of their own)
    "ruv_k": [1, 2, 4, 8, 12],
    "ruv_controls": ["empirical", "all"],
    "ruv_control_share": 0.5,    # empirical: this share of genes, least associated with the disease
    "ruv_tolerance": 0.005,
    # pseudo-replicates: per batch, the mean of the samples (at least
    # ruv_prps_min) of one group of ruv_prps_by (numeric ones in
    # ruv_prps_bands), linked across batches; weighted by sqrt(samples)
    "ruv_prps_by": ["dementia", "sex", "age_capped"],
    "ruv_prps_bands": [70, 80, 90],
    "ruv_prps_min": 3,
    "ruv_prps_weighted": True,
}
REQUIRED = ("output_dir", "donor_table", "datasets")
DATASET_KEYS = {"counts", "layer", "gene_map", "donor", "celltype"}
METHODS = ("none", "center", "regress", "regress_dataset", "anchor", "ruv")
RUV_CONTROLS = ("empirical", "all")
RUV_REPORTED_FACTORS = 4     # unwanted factors per cell type in ruv_factors.tsv
# Donor table columns copied to the samples (if present)
COVARIATES = ["n_datasets", "cohort", "brain_bank", "sex", "age_death", "age_capped", "race", "hispanic",
              "pmi", "apoe4", "cogdx", "dementia", "cdr", "braak", "cerad", "thal", "adnc", "crossdis_dx"]
NUMERIC = {"n_datasets", "age_death", "age_capped", "pmi", "apoe4", "cogdx", "cdr", "braak", "cerad",
           "thal", "adnc"}
METRIC_COLUMNS = ["celltype", "method", "metric", "label", "value", "chance", "n"]


# ---------------------------------------------------------------- spec


def load_spec(path: Union[str, Path]) -> Dict[str, Any]:
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
        for key in ("counts", "gene_map", "donor", "celltype"):
            if not cfg.get(key):
                raise ValueError(f"Dataset {name} needs '{key}'")
    for key in ("protect", "disease_adjust"):
        bad = [c for c in spec[key] if c not in COVARIATES]
        if bad:
            raise ValueError(f"{key}: {bad} are not donor table covariates ({COVARIATES})")
    for name, column in (spec["batch_split"] or {}).items():
        if name not in spec["datasets"]:
            raise ValueError(f"batch_split: {name} is not one of the datasets")
        if column not in COVARIATES:
            raise ValueError(f"batch_split: {column} is not a donor table covariate ({COVARIATES})")
    bad = [m for m in spec["umap_methods"] if m not in METHODS]
    if bad:
        raise ValueError(f"umap_methods: {bad} are not methods ({METHODS})")
    bad = [c for c in spec["ruv_prps_by"] if c not in COVARIATES]
    if bad:
        raise ValueError(f"ruv_prps_by: {bad} are not donor table covariates ({COVARIATES})")
    bad = [c for c in spec["ruv_controls"] if c not in RUV_CONTROLS]
    if not spec["ruv_controls"] or bad:
        raise ValueError(f"ruv_controls: one or more of {RUV_CONTROLS}, not {spec['ruv_controls']}")
    if not spec["ruv_k"] or any(int(k) < 1 for k in spec["ruv_k"]):
        raise ValueError(f"ruv_k: numbers of factors of at least 1, not {spec['ruv_k']}")
    return spec


# ---------------------------------------------------------------- pseudobulk


def shared_genes(datasets) -> Tuple[List[str], Dict[str, np.ndarray]]:
    """Ensembl IDs that every dataset's file has (via its gene map), sorted,
    and per dataset the file column of each."""
    from ..data.genes import load_gene_map

    found = {}
    for ds in datasets:
        mapping = load_gene_map(ds.cfg["gene_map"])
        found[ds.name] = {gid.split(".")[0] for name in ds.var_names if (gid := mapping.get(name))}
    genes = sorted(set.intersection(*found.values()))
    columns = {}
    for ds in datasets:
        lookup = _columns_for(ds, genes)
        columns[ds.name] = np.array([lookup[g] for g in genes])
        logger.info(f"[{ds.name}] {len(found[ds.name])} genes mapped; {len(genes)} shared by every dataset")
    return genes, columns


def pseudobulk_dataset(ds, columns: np.ndarray, chunk_size: int) -> Tuple[pd.DataFrame, np.ndarray]:
    """Sums per donor x cell type over all kept cells, for the given file
    columns (in their order), and the library sizes over all genes."""
    ds.union = np.unique(columns)
    ds.sample = pd.DataFrame({"row": np.array([], dtype=np.int64)})
    stream_counts(ds, chunk_size)
    order = np.searchsorted(ds.union, columns)
    groups = ds.groups.assign(dataset=ds.name)
    return groups, ds.pseudobulk[:, order]


def _symbols(datasets, columns: Dict[str, np.ndarray]) -> pd.DataFrame:
    return pd.DataFrame({f"{ds.name}_name": np.asarray(ds.var_names, dtype=object)[columns[ds.name]]
                         for ds in datasets})


def read_donor_table(path: Union[str, Path], names: Sequence[str]) -> pd.DataFrame:
    keys = {"person_id": str, **{f"{k}_donor": str for k in names}}
    table = pd.read_csv(path, sep="\t", dtype=keys, low_memory=False)
    missing = [c for c in keys if c not in table.columns]
    if missing:
        raise ValueError(f"{path} has no column(s) {missing}; is it the donor table of these datasets?")
    return table


def attach_persons(groups: pd.DataFrame, table: pd.DataFrame) -> pd.DataFrame:
    """Person and covariates per sample, from each dataset's donor key. Rows
    keep their index; samples of donors not in the table are left out."""
    parts = []
    for name, part in groups.groupby("dataset", sort=False):
        lookup = table.dropna(subset=[f"{name}_donor"]).set_index(f"{name}_donor")["person_id"]
        part = part.assign(person_id=part["donor"].astype(str).map(lookup))
        unmatched = part.loc[part["person_id"].isna(), "donor"].nunique()
        if unmatched:
            logger.warning(f"[{name}] {unmatched} donors are not in the donor table (left out)")
        parts.append(part[part["person_id"].notna()])
    out = pd.concat(parts).sort_index()
    covs = table.set_index("person_id")[[c for c in COVARIATES if c in table.columns]]
    return out.join(covs, on="person_id")


def assign_batches(obs: pd.DataFrame, split: Dict[str, str]) -> np.ndarray:
    """The dataset per sample, or 'dataset:level' for the datasets split by a
    column ('dataset:unknown' where it has no value)."""
    batch = obs["dataset"].astype(str).to_numpy(dtype=object)
    for name, column in (split or {}).items():
        if column not in obs.columns:
            raise ValueError(f"batch_split: the samples have no column {column}; rebuild pseudobulk.h5ad")
        rows = (obs["dataset"] == name).to_numpy()
        level = obs.loc[rows, column].astype(object)
        batch[rows] = (name + ":" + level.where(level.notna(), "unknown").astype(str)).to_numpy()
    return batch


# ---------------------------------------------------------------- corrections


def _fill_unanchored(out: np.ndarray, Y: np.ndarray, ds: np.ndarray, anchored: set) -> np.ndarray:
    """Datasets that could not be corrected otherwise get the mean of the others."""
    rest = [d for d in np.unique(ds) if d not in anchored]
    done = np.isin(ds, list(anchored))
    target = out[done].mean(axis=0) if done.any() else Y.mean(axis=0)
    for d in rest:
        m = ds == d
        out[m] = Y[m] - Y[m].mean(axis=0) + target
    return out


def correct_center(Y: np.ndarray, ds: np.ndarray) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    grand = Y.mean(axis=0)
    offsets = {d: Y[ds == d].mean(axis=0) - grand for d in np.unique(ds)}
    out = Y.copy()
    for d, o in offsets.items():
        out[ds == d] -= o
    return out, offsets


def correct_regress(Y: np.ndarray, ds: np.ndarray, covs: pd.DataFrame, numeric: Sequence[str],
                    min_rows: int = 10) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Remove the dataset terms of Y ~ dataset + covs, fitted on the rows with
    every covariate; covariates that do not vary there are left out."""
    complete = covs.notna().all(axis=1).to_numpy()
    keep = [c for c in covs.columns if covs.loc[complete, c].nunique() > 1]
    levels = [d for d in np.unique(ds) if (complete & (ds == d)).sum() >= min_rows]
    fit = complete & np.isin(ds, levels)
    if len(levels) < 2:
        logger.warning("regress: fewer than 2 datasets with enough complete rows; centering instead")
        return correct_center(Y, ds)
    D = (ds[fit][:, None] == np.array(levels)[None, :]).astype(np.float64)
    if keep:
        C, _ = bm.design_matrix(covs.loc[fit, keep].reset_index(drop=True), [c for c in keep if c in numeric])
        X = np.hstack([D, C[:, 1:]])
    else:
        X = D
    beta, *_ = np.linalg.lstsq(X, Y[fit], rcond=None)
    intercepts = beta[:len(levels)]
    weights = np.array([(ds == d).sum() for d in levels], dtype=np.float64)
    reference = (weights[:, None] * intercepts).sum(axis=0) / weights.sum()
    offsets = {d: intercepts[i] - reference for i, d in enumerate(levels)}
    out = Y.copy()
    for d, o in offsets.items():
        out[ds == d] -= o
    out = _fill_unanchored(out, Y, ds, set(levels))
    for d in set(np.unique(ds)) - set(levels):
        offsets[d] = Y[ds == d].mean(axis=0) - out[ds == d].mean(axis=0)
    return out, offsets


def _replicate_pairs(persons: np.ndarray, ds: np.ndarray, exclude=frozenset()) -> List[Tuple[int, int]]:
    rows: Dict[str, List[int]] = {}
    for i, p in enumerate(persons):
        if p not in exclude:
            rows.setdefault(p, []).append(i)
    return [(i, j) for members in rows.values() if len(members) > 1 for i, j in combinations(members, 2)
            if ds[i] != ds[j]]


def _components(edges: List[Tuple[str, str]], nodes: Sequence[str]) -> List[set]:
    parent = {n: n for n in nodes}

    def find(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    for a, b in edges:
        parent[find(a)] = find(b)
    groups: Dict[str, set] = {}
    for n in nodes:
        groups.setdefault(find(n), set()).add(n)
    return list(groups.values())


def correct_anchor(Y: np.ndarray, ds: np.ndarray, persons: np.ndarray,
                   exclude=frozenset()) -> Tuple[np.ndarray, Dict[str, np.ndarray], set]:
    """Per gene, dataset offsets o that best fit y[p, a] - y[p, b] = o_a - o_b
    over the persons measured in two datasets (least squares, offsets summing
    to zero). Only the largest group of datasets linked by such persons is
    anchored; other datasets are centered on it. Returns the corrected
    values, the offsets and the anchored datasets."""
    pairs = _replicate_pairs(persons, ds, exclude)
    levels = list(np.unique(ds))
    linked = [c for c in _components([(ds[i], ds[j]) for i, j in pairs], levels) if len(c) > 1]
    anchored = max(linked, key=lambda c: sum((ds == d).sum() for d in c)) if linked else set()
    out = Y.copy()
    offsets: Dict[str, np.ndarray] = {}
    if anchored:
        order = sorted(anchored)
        pairs = [(i, j) for i, j in pairs if ds[i] in anchored]
        Dp = np.zeros((len(pairs), len(order)))
        for r, (i, j) in enumerate(pairs):
            Dp[r, order.index(ds[i])] = 1.0
            Dp[r, order.index(ds[j])] = -1.0
        delta = Y[[i for i, _ in pairs]] - Y[[j for _, j in pairs]]
        solved = np.linalg.pinv(Dp) @ delta
        for k, d in enumerate(order):
            offsets[d] = solved[k]
            out[ds == d] -= solved[k]
    out = _fill_unanchored(out, Y, ds, set(anchored))
    for d in set(levels) - set(anchored):
        offsets[d] = Y[ds == d].mean(axis=0) - out[ds == d].mean(axis=0)
    return out, offsets, set(anchored)


def prps_groups(sub: pd.DataFrame, by: Sequence[str], bands: Sequence[float]) -> np.ndarray:
    """Per sample, its group for the pseudo-replicates: the levels of `by`,
    numeric columns cut at `bands`; None where one is missing."""
    labels = np.full(len(sub), "", dtype=object)
    ok = np.ones(len(sub), dtype=bool)
    for c in by:
        if c in NUMERIC:
            x = pd.to_numeric(sub[c], errors="coerce").to_numpy(np.float64)
            ok &= np.isfinite(x)
            part = np.digitize(np.nan_to_num(x), np.asarray(bands, dtype=np.float64)).astype(str)
        else:
            ok &= sub[c].notna().to_numpy()
            part = sub[c].astype(str).to_numpy()
        labels = labels + "|" + part.astype(object)
    labels[~ok] = None
    return labels


def ruv_factors(Y: np.ndarray, batch: np.ndarray, persons: np.ndarray, groups: Optional[np.ndarray] = None,
                exclude=frozenset(), min_prps: int = 3, n_factors: int = 12,
                weighted: bool = True) -> Tuple[Optional[np.ndarray], Dict[str, int]]:
    """RUV-III's unwanted factors (ruv::RUVIII): alpha = the top n_factors
    right singular vectors, times their singular values, of the replicate
    residuals Y0 = residop(Y, M) (each row minus its replicate set's mean).

    M: each person's samples (persons not in exclude, two or more samples;
    the others are singletons, whose residual is 0) and, RUV-III-PRPS style,
    pseudo-replicates: per batch and group (groups, see prps_groups), the mean
    of its samples of persons not in exclude (at least min_prps), linked across
    the batches that have the group. With weighted, a pseudo-replicate's
    residual (from the sample-weighted mean of its set) is weighted by the
    square root of its number of samples, as if each sample had a replicate.

    Returns alpha (None without replicate sets) and the counts of replicate
    sets and pseudo-replicates."""
    blocks = []
    members: Dict[str, List[int]] = {}
    for i, p in enumerate(persons):
        if p not in exclude:
            members.setdefault(p, []).append(i)
    n_sets = 0
    for rows in members.values():
        if len(rows) > 1:
            blocks.append(Y[rows] - Y[rows].mean(axis=0))
            n_sets += 1
    n_prps = n_prps_sets = 0
    if groups is not None:
        cells: Dict[Tuple[str, str], List[int]] = {}
        for i in range(len(Y)):
            if groups[i] is not None and persons[i] not in exclude:
                cells.setdefault((groups[i], batch[i]), []).append(i)
        by_group: Dict[str, List[Tuple[np.ndarray, int]]] = {}
        for (g, _), rows in sorted(cells.items()):
            if len(rows) >= min_prps:
                by_group.setdefault(g, []).append((Y[rows].mean(axis=0), len(rows)))
        for items in by_group.values():
            if len(items) < 2:
                continue
            P = np.vstack([m for m, _ in items])
            n = np.array([c for _, c in items], dtype=np.float64) if weighted else np.ones(len(items))
            residual = P - (n[:, None] * P).sum(axis=0) / n.sum()
            blocks.append(np.sqrt(n)[:, None] * residual)
            n_prps += len(items)
            n_prps_sets += 1
    info = dict(replicate_sets=n_sets, pseudo_replicates=n_prps, pseudo_replicate_sets=n_prps_sets)
    if not blocks:
        return None, info
    Y0 = np.vstack(blocks)
    values, U = np.linalg.eigh(Y0 @ Y0.T)
    order = np.argsort(values)[::-1]
    values, U = values[order], U[:, order]
    rank = int((values > max(values[0], 0) * 1e-10).sum()) if len(values) else 0
    if rank == 0:
        return None, info
    k = min(int(n_factors), rank)
    return U[:, :k].T @ Y0, info


def ruv_apply(Y: np.ndarray, alpha: np.ndarray, k: int,
              controls: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Y - W alpha_k, with W = Y_c alpha_c' (alpha_c alpha_c')^-1 over the
    control genes (all genes if None), as in ruv::RUVIII. Y is centered per
    gene first, so the gene means stay as they are (ruv::RUVIII on uncentered
    Y differs from this only by a constant per gene). Returns the corrected
    values and W (samples x k)."""
    a = alpha[:min(int(k), len(alpha))]
    Yc = Y - Y.mean(axis=0)
    ctl = slice(None) if controls is None else np.asarray(controls)
    W = np.linalg.lstsq(a[:, ctl].T, Yc[:, ctl].T, rcond=None)[0].T
    return Y - W @ a, W


def correct_ruv(Y: np.ndarray, batch: np.ndarray, persons: np.ndarray, k: int,
                groups: Optional[np.ndarray] = None, controls: Optional[np.ndarray] = None,
                exclude=frozenset(), min_prps: int = 3,
                weighted: bool = True) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    """RUV-III with k unwanted factors (ruv_factors, ruv_apply); without
    replicate sets, the batches are centered. Returns the corrected values,
    W and the counts of ruv_factors."""
    alpha, info = ruv_factors(Y, batch, persons, groups, exclude, min_prps, k, weighted)
    if alpha is None:
        logger.warning("ruv: no replicate sets; centering instead")
        return correct_center(Y, batch)[0], np.zeros((len(Y), 0)), info
    out, W = ruv_apply(Y, alpha, k, controls)
    return out, W, info


# ---------------------------------------------------------------- metrics


def _row_standardize(Z: np.ndarray) -> np.ndarray:
    Z = Z - Z.mean(axis=1, keepdims=True)
    sd = np.sqrt((Z ** 2).mean(axis=1, keepdims=True))
    return Z / np.where(sd > 0, sd, 1.0)


def replicate_stats(Yc: np.ndarray, ds: np.ndarray, persons: np.ndarray, queries: Optional[set] = None,
                    center: Optional[np.ndarray] = None) -> Dict[Tuple[str, str], Dict[str, list]]:
    """Per pair of datasets, for each person in both (and in queries): the
    rank of their own profile among the other dataset's samples (both
    directions), the number of candidates, and the correlation with their own
    and (on average) with other persons' profiles."""
    Z = _row_standardize(Yc - (Yc.mean(axis=0) if center is None else center))
    n_genes = Z.shape[1]
    out: Dict[Tuple[str, str], Dict[str, list]] = {}
    for a, b in combinations(list(dict.fromkeys(ds)), 2):
        ia, ib = np.flatnonzero(ds == a), np.flatnonzero(ds == b)
        pa = pd.Series(ia, index=persons[ia])
        pb = pd.Series(ib, index=persons[ib])
        both = sorted(set(pa.index) & set(pb.index) & (queries if queries is not None else set(pa.index)))
        if not both:
            continue
        acc = out.setdefault((a, b), {"rank": [], "n_candidates": [], "own_r": [], "other_r": []})
        for src, dst, src_rows, dst_rows in ((pa, pb, ia, ib), (pb, pa, ib, ia)):
            sims = Z[src[both].to_numpy()] @ Z[dst_rows].T / n_genes
            own_col = pd.Index(dst_rows).get_indexer(dst[both].to_numpy())
            own = sims[np.arange(len(both)), own_col]
            acc["rank"] += list((sims > own[:, None]).sum(axis=1))
            acc["n_candidates"] += [len(dst_rows)] * len(both)
            acc["own_r"] += list(own)
            acc["other_r"] += list((sims.sum(axis=1) - own) / max(len(dst_rows) - 1, 1))
    return out


def _replicate_rows(celltype: str, method: str, acc: Dict[Tuple[str, str], Dict[str, list]],
                    min_pairs: int) -> List[Dict[str, Any]]:
    rows = []
    for (a, b), v in acc.items():
        n = len(v["rank"]) // 2
        if n < min_pairs:
            continue
        rank, cand = np.array(v["rank"], dtype=np.float64), np.array(v["n_candidates"], dtype=np.float64)
        base = dict(celltype=celltype, method=method, label=f"{a} vs {b}", n=n)
        rows += [
            dict(base, metric="replicate_top1", value=float((rank == 0).mean()), chance=float((1 / cand).mean())),
            dict(base, metric="replicate_rank_share", value=float(np.median(rank / np.maximum(cand - 1, 1))),
                 chance=0.5),
            dict(base, metric="replicate_r_own", value=float(np.mean(v["own_r"])), chance=np.nan),
            dict(base, metric="replicate_r_other", value=float(np.mean(v["other_r"])), chance=np.nan),
        ]
    return rows


def _ols(Y: np.ndarray, X: np.ndarray, column: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Coefficient, t and two-sided p of one column of X, per column of Y."""
    n, p = X.shape
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ (X.T @ Y)
    rss = ((Y - X @ beta) ** 2).sum(axis=0)
    df = max(n - p, 1)
    se = np.sqrt(rss / df * xtx_inv[column, column])
    t = np.divide(beta[column], se, out=np.full(Y.shape[1], np.nan), where=se > 0)
    pvalue = np.where(np.isfinite(t), 2 * stats.t.sf(np.abs(np.nan_to_num(t)), df), 1.0)
    return beta[column], t, pvalue


def _disease_design(frame: pd.DataFrame, disease: str, positive: str, adjust: List[str],
                    with_batch: bool) -> np.ndarray:
    blocks = [np.ones((len(frame), 1)), (frame[disease].astype(str) == positive).to_numpy(np.float64)[:, None]]
    if adjust:
        C, _ = bm.design_matrix(frame[adjust].reset_index(drop=True), [c for c in adjust if c in NUMERIC])
        blocks.append(C[:, 1:])
    if with_batch:
        blocks.append(pd.get_dummies(frame["batch"].astype(str), drop_first=True).to_numpy(np.float64))
    return np.hstack(blocks)


def disease_analysis(celltype: str, corrected: Dict[str, np.ndarray], sub: pd.DataFrame, genes: pd.DataFrame,
                     spec) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
    """Rows of disease.tsv and the per-gene reference model (de_pooled.tsv)."""
    disease, positive, adjust = spec["disease"], str(spec["disease_positive"]), list(spec["disease_adjust"])
    needed = [disease] + adjust
    usable = sub[needed].notna().all(axis=1).to_numpy()
    first = sub[usable].sort_values("n_cells", ascending=False).drop_duplicates("person_id").index.to_numpy()
    first = np.sort(first)
    frame = sub.loc[first]
    n_pos = int((frame[disease].astype(str) == positive).sum())
    if n_pos < 3 or len(frame) - n_pos < 3:
        logger.info(f"{celltype}: too few persons with and without {disease} for the disease models")
        return [], pd.DataFrame()

    rows: List[Dict[str, Any]] = []
    X = _disease_design(frame, disease, positive, adjust, with_batch=True)
    beta, t_ref, p = _ols(corrected["none"][first], X, 1)
    fdr = bm.bh_fdr(p)
    de = genes.assign(celltype=celltype, log2fc=beta, t=t_ref, p_value=p, fdr=fdr, n_persons=len(frame),
                      n_positive=n_pos)

    def summary(analysis, t, pvals, n, npos):
        ok = np.isfinite(t) & np.isfinite(t_ref)
        corr = float(np.corrcoef(t[ok], t_ref[ok])[0, 1]) if ok.sum() > 2 else np.nan
        return dict(celltype=celltype, analysis=analysis, share_fdr05=float((bm.bh_fdr(pvals) < 0.05).mean()),
                    n_genes_fdr05=int((bm.bh_fdr(pvals) < 0.05).sum()), n_genes=len(t),
                    t_correlation_with_pooled=corr, n_persons=n, n_positive=npos)

    rows.append(summary("pooled: batch as covariate", t_ref, p, len(frame), n_pos))
    Xm = _disease_design(frame, disease, positive, adjust, with_batch=False)
    for method, Yc in corrected.items():
        _, t, pv = _ols(Yc[first], Xm, 1)
        rows.append(summary(f"pooled: {method}, no batch term", t, pv, len(frame), n_pos))
    for name in dict.fromkeys(sub["dataset"]):
        sel = np.flatnonzero(usable & (sub["dataset"] == name).to_numpy())
        f = sub.iloc[sel]
        npos = int((f[disease].astype(str) == positive).sum())
        if npos < 3 or len(f) - npos < 3:
            continue
        Xd = _disease_design(f, disease, positive, adjust, with_batch=f["batch"].nunique() > 1)
        _, t, pv = _ols(corrected["none"][sel], Xd, 1)
        rows.append(summary(f"{name} alone", t, pv, len(f), npos))
    return rows, de


def _mixing_labels(ds: np.ndarray, batch: np.ndarray,
                   split: Dict[str, str]) -> List[Tuple[str, np.ndarray, Optional[np.ndarray]]]:
    """(label, values, rows or None for all) for the mixing metrics: the
    dataset; with split batches also the batch, and the split column within
    each split dataset."""
    out = [("dataset", ds, None)]
    if (batch != ds).any():
        out.append(("batch", batch, None))
        for name, column in (split or {}).items():
            rows = np.flatnonzero(ds == name)
            if len(np.unique(batch[rows])) > 1:
                out.append((f"{column} within {name}", batch[rows], rows))
    return out


def _umap(pcs: np.ndarray, spec) -> np.ndarray:
    import umap

    n_neighbors = min(int(spec["umap_neighbors"]), len(pcs) - 1)
    return umap.UMAP(n_neighbors=n_neighbors, min_dist=UMAP_MIN_DIST,
                     random_state=spec["seed"]).fit_transform(pcs)


def reliability_rows(celltype: str, method: str, Y: np.ndarray, ds: np.ndarray, persons: np.ndarray,
                     min_pairs: int) -> List[Dict[str, Any]]:
    """cross_dataset_reliability per pair of datasets with at least min_pairs
    persons in both."""
    rows = []
    for a, b in combinations(list(dict.fromkeys(ds)), 2):
        ia, ib = np.flatnonzero(ds == a), np.flatnonzero(ds == b)
        pa, pb = pd.Series(ia, index=persons[ia]), pd.Series(ib, index=persons[ib])
        both = sorted(set(pa.index) & set(pb.index))
        if len(both) >= min_pairs:
            value = bm.shared_reliability(Y[pa[both].to_numpy()], Y[pb[both].to_numpy()])
            rows.append(dict(celltype=celltype, method=method, metric="cross_dataset_reliability",
                             label=f"{a} vs {b}", value=value, chance=0.0, n=len(both)))
    return rows


def score_methods(celltype: str, sub: pd.DataFrame, corrected: Dict[str, np.ndarray], spec,
                  rng: np.random.Generator, Y: Optional[np.ndarray] = None,
                  replicates: Optional[Dict[str, Dict[Tuple[str, str], Dict[str, list]]]] = None):
    """Mixing and replicate metrics of each method's matrix (samples x genes,
    rows as in sub), and each method's PCs. With Y (the uncorrected matrix),
    the anchor method's replicate metrics are cross-fitted. replicates: per
    method, replicate_stats computed elsewhere (e.g. cross-fitted) to use
    instead."""
    ds = sub["dataset"].to_numpy().astype(str)
    batch = sub["batch"].to_numpy().astype(str)
    persons = sub["person_id"].to_numpy().astype(str)
    mixing = _mixing_labels(ds, batch, spec["batch_split"])
    metrics: List[Dict[str, Any]] = []
    all_pcs = {}
    for method, Yc in corrected.items():
        pcs = all_pcs[method] = _pca(Yc - Yc.mean(axis=0), spec["n_pcs"], spec["seed"])
        for label, values, rows in mixing:
            p = pcs if rows is None else pcs[rows]
            base = dict(celltype=celltype, method=method, label=label, n=len(values))
            metrics.append(dict(base, metric="pc_r2", value=bm.pc_regression(p, values), chance=np.nan))
            accuracy, chance, _ = bm.classify(p, values, rng)
            metrics.append(dict(base, metric="classifier_accuracy", value=accuracy, chance=chance))
        if replicates is not None and method in replicates:
            acc = replicates[method]
        elif method != "anchor" or Y is None:
            acc = replicate_stats(Yc, ds, persons)
        else:  # cross-fitted: a person's offsets come from the other persons' pairs
            replicated = sorted({persons[i] for i, _ in _replicate_pairs(persons, batch)})
            folds = np.arange(len(replicated)) % spec["n_folds"]
            rng.shuffle(folds)
            acc = {}
            for k in range(spec["n_folds"]):
                fold = {p for p, f in zip(replicated, folds) if f == k}
                if not fold:
                    continue
                Yk, _, _ = correct_anchor(Y, batch, persons, exclude=frozenset(fold))
                for key, v in replicate_stats(Yk, ds, persons, queries=fold,
                                              center=corrected["anchor"].mean(axis=0)).items():
                    target = acc.setdefault(key, {k2: [] for k2 in v})
                    for k2, values in v.items():
                        target[k2] += values
        metrics += _replicate_rows(celltype, method, acc, spec["min_pairs"])
    return metrics, all_pcs


# ---------------------------------------------------------------- RUV-III: choice and diagnostics


def _merge_replicates(target: Dict, part: Dict) -> None:
    for key, v in part.items():
        acc = target.setdefault(key, {k2: [] for k2 in v})
        for k2, values in v.items():
            acc[k2] += values


def _overall_replicates(acc: Dict, min_pairs: int) -> Tuple[float, float, int]:
    """Top-1 and median rank share over all pairs of datasets with at least
    min_pairs persons (both directions), and the number of queries."""
    rank, cand = [], []
    for v in acc.values():
        if len(v["rank"]) // 2 >= min_pairs:
            rank += v["rank"]
            cand += v["n_candidates"]
    if not rank:
        return np.nan, np.nan, 0
    rank, cand = np.array(rank, dtype=np.float64), np.array(cand, dtype=np.float64)
    return float((rank == 0).mean()), float(np.median(rank / np.maximum(cand - 1, 1))), len(rank)


def _disease_inputs(sub: pd.DataFrame, spec) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """As in disease_analysis: the rows of one sample per person (the most
    nuclei) and the disease design with and without the batch; None with
    fewer than 3 persons with or without the disease."""
    disease, positive, adjust = spec["disease"], str(spec["disease_positive"]), list(spec["disease_adjust"])
    usable = sub[[disease] + adjust].notna().all(axis=1).to_numpy()
    first = np.sort(sub[usable].sort_values("n_cells", ascending=False).drop_duplicates("person_id").index.to_numpy())
    frame = sub.loc[first]
    n_pos = int((frame[disease].astype(str) == positive).sum())
    if n_pos < 3 or len(frame) - n_pos < 3:
        return None
    return (first, _disease_design(frame, disease, positive, adjust, with_batch=True),
            _disease_design(frame, disease, positive, adjust, with_batch=False))


def ruv_control_sets(spec, reference_t: Optional[np.ndarray]) -> Dict[str, Optional[np.ndarray]]:
    """The control gene sets of ruv_controls: all (None) and empirical, the
    ruv_control_share of genes with the smallest |t| for the disease in the
    reference model (batch as covariate); empirical is left out without it."""
    sets: Dict[str, Optional[np.ndarray]] = {}
    for name in spec["ruv_controls"]:
        if name == "all":
            sets[name] = None
        elif reference_t is None:
            logger.info("ruv: no disease model, so no empirical control genes")
        else:
            a = np.abs(reference_t)
            ok = np.isfinite(a)
            sets[name] = np.flatnonzero(ok & (a <= np.quantile(a[ok], float(spec["ruv_control_share"]))))
    if not sets:
        sets["all"] = None
    return sets


def ruv_scan(celltype: str, sub: pd.DataFrame, Y: np.ndarray, spec, rng: np.random.Generator,
             disease_inputs: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None) -> Optional[Dict[str, Any]]:
    """RUV-III with each number of factors (ruv_k) and control gene set
    (ruv_controls). Replicate metrics are cross-fitted: n_folds folds of the
    persons in two batches, a fold's persons singletons in M, left out of the
    pseudo-replicates, and scored as queries. Mixing (pc_r2) and disease
    metrics are of the fit with all persons. The choice: of those within
    ruv_tolerance of the best cross-fitted top-1 over all pairs of datasets,
    the first control set of ruv_controls and the most factors (the
    replicates cannot see the batches without bridges).

    Returns None without replicate sets, else rows (ruv_scan.tsv), choice
    (controls, k), the corrected values Y and W of the choice, its alpha and
    cross-fitted replicate_stats (acc), and the counts of ruv_factors."""
    ds = sub["dataset"].to_numpy().astype(str)
    batch = sub["batch"].to_numpy().astype(str)
    persons = sub["person_id"].to_numpy().astype(str)
    groups = prps_groups(sub, spec["ruv_prps_by"], spec["ruv_prps_bands"])
    ks = sorted({int(k) for k in spec["ruv_k"]})
    prps = dict(min_prps=int(spec["ruv_prps_min"]), weighted=bool(spec["ruv_prps_weighted"]))
    alpha, info = ruv_factors(Y, batch, persons, groups, n_factors=max(ks), **prps)
    if alpha is None:
        logger.warning(f"{celltype}: ruv: no replicate sets or pseudo-replicates; no ruv")
        return None
    reference_t = None
    if disease_inputs is not None:
        first, X_ref, X_plain = disease_inputs
        _, reference_t, _ = _ols(Y[first], X_ref, 1)
    controls = ruv_control_sets(spec, reference_t)
    configs = [(c, k) for c in controls for k in ks[::-1]]   # in order of preference

    replicated = sorted({persons[i] for i, _ in _replicate_pairs(persons, batch)})
    folds = np.arange(len(replicated)) % spec["n_folds"]
    rng.shuffle(folds)
    acc: Dict[Tuple[str, int], Dict] = {cfg: {} for cfg in configs}
    for f in range(spec["n_folds"]):
        fold = frozenset(p for p, x in zip(replicated, folds) if x == f)
        if not fold:
            continue
        alpha_f, _ = ruv_factors(Y, batch, persons, groups, exclude=fold, n_factors=max(ks), **prps)
        if alpha_f is None:
            continue
        for c, k in configs:
            Yk, _ = ruv_apply(Y, alpha_f, k, controls[c])
            _merge_replicates(acc[(c, k)], replicate_stats(Yk, ds, persons, queries=set(fold)))

    mixing = _mixing_labels(ds, batch, spec["batch_split"])
    rows: List[Dict[str, Any]] = []
    score: Dict[Tuple[str, int], float] = {}
    for c, k in configs:
        n_controls = Y.shape[1] if controls[c] is None else len(controls[c])
        base = dict(celltype=celltype, controls=c, n_controls=n_controls, k=k, n_factors=min(k, len(alpha)))
        top1, share, n = _overall_replicates(acc[(c, k)], spec["min_pairs"])
        score[(c, k)] = top1
        rows += [dict(base, metric="replicate_top1", label="all pairs", value=top1, n=n),
                 dict(base, metric="replicate_rank_share", label="all pairs", value=share, n=n)]
        for r in _replicate_rows(celltype, "ruv", acc[(c, k)], spec["min_pairs"]):
            if r["metric"] in ("replicate_top1", "replicate_rank_share"):
                rows.append(dict(base, metric=r["metric"], label=r["label"], value=r["value"], n=r["n"]))
        Yk, _ = ruv_apply(Y, alpha, k, controls[c])
        pcs = _pca(Yk - Yk.mean(axis=0), spec["n_pcs"], spec["seed"])
        for label, values, sel in mixing:
            p = pcs if sel is None else pcs[sel]
            rows.append(dict(base, metric="pc_r2", label=label, value=bm.pc_regression(p, values), n=len(values)))
        if reference_t is not None:
            _, t, pv = _ols(Yk[first], X_plain, 1)
            ok = np.isfinite(t) & np.isfinite(reference_t)
            corr = float(np.corrcoef(t[ok], reference_t[ok])[0, 1]) if ok.sum() > 2 else np.nan
            rows += [dict(base, metric="disease_t_correlation", label="with pooled: batch as covariate",
                          value=corr, n=len(first)),
                     dict(base, metric="disease_share_fdr05", label="no batch term",
                          value=float((bm.bh_fdr(pv) < 0.05).mean()), n=len(first))]
    valid = {cfg: s for cfg, s in score.items() if np.isfinite(s)}
    if valid:
        best = max(valid.values())
        choice = next(cfg for cfg in configs if cfg in valid and valid[cfg] >= best - float(spec["ruv_tolerance"]))
    else:
        choice = (configs[0][0], max(ks))
        logger.warning(f"{celltype}: ruv: no replicates to choose k by; k = {choice[1]}")
    for r in rows:
        r["chosen"] = (r["controls"], r["k"]) == choice
    curve = ", ".join(f"{c} k={k}: {score[(c, k)]:.3f}" for c in controls for k in ks)
    logger.info(f"{celltype}: ruv: {info['replicate_sets']} replicate sets, {info['pseudo_replicates']} "
                f"pseudo-replicates in {info['pseudo_replicate_sets']} sets; cross-fitted top-1 {curve}; "
                f"chose {choice[0]} controls, k = {choice[1]}")
    Y_ruv, W = ruv_apply(Y, alpha, choice[1], controls[choice[0]])
    return dict(rows=rows, choice=choice, Y=Y_ruv, W=W, alpha=alpha[:W.shape[1]], acc=acc[choice], info=info)


def _sample_qc(sub: pd.DataFrame) -> Dict[str, np.ndarray]:
    nuclei = sub["n_cells"].to_numpy(np.float64)
    library = sub["library_size"].to_numpy(np.float64)
    return {"log10 nuclei": np.log10(nuclei), "log10 library size": np.log10(library),
            "log10 UMIs per nucleus": np.log10(library / nuclei)}


def _within(x: np.ndarray, groups: np.ndarray) -> np.ndarray:
    s = pd.Series(np.asarray(x, dtype=np.float64))
    return (s - s.groupby(np.asarray(groups)).transform("mean")).to_numpy()


def _corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def _factor_rows(celltype: str, source: str, scores: np.ndarray, shares: np.ndarray, qc: Dict[str, np.ndarray],
                 groups: np.ndarray) -> List[Dict[str, Any]]:
    """Per factor (column of scores) and QC measure: the Spearman correlation,
    the Pearson correlation within the groups (both centered per group), and
    the share of the factor's variance that the groups explain."""
    rows = []
    for f in range(scores.shape[1]):
        s = scores[:, f]
        within = _within(s, groups)
        total = ((s - s.mean()) ** 2).sum()
        group_r2 = float(1 - (within ** 2).sum() / total) if total > 0 else np.nan
        for name, x in qc.items():
            ok = np.isfinite(x) & np.isfinite(s)
            if ok.sum() < 3:
                continue
            xw = _within(np.where(ok, x, 0.0), groups)
            rows.append(dict(celltype=celltype, source=source, factor=f + 1, variance_share=float(shares[f]),
                             qc=name, r_spearman=_corr(stats.rankdata(s[ok]), stats.rankdata(x[ok])),
                             r_within_batch=_corr(within[ok], xw[ok]), batch_r2=group_r2, n=int(ok.sum())))
    return rows


def bridge_decomposition(celltype: str, Y: np.ndarray, sub: pd.DataFrame, ks: Sequence[int],
                         rng: np.random.Generator, min_pairs: int) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """The within-person differences y_i - y_j of the persons in two batches
    (one row per pair of their samples), as shares of their sum of squares:
    batch_offset (least squares of y_i - y_j = o_a - o_b per gene, as
    anchor), remainder_top_k (the top k components of what is left) and rest;
    remainder_top_k_null: the same with each gene permuted over the pairs (no
    shared structure). Rows of ruv_diagnostic.tsv, and of ruv_factors.tsv for
    the top components: their pair scores against the pairs' differences in
    sample QC, within and between pairs of batches."""
    batch = sub["batch"].to_numpy().astype(str)
    persons = sub["person_id"].to_numpy().astype(str)
    pairs = _replicate_pairs(persons, batch)
    if len(pairs) < min_pairs:
        return [], []
    i = np.array([a for a, _ in pairs])
    j = np.array([b for _, b in pairs])
    D = Y[i] - Y[j]
    levels = sorted(set(batch[i]) | set(batch[j]))
    pos = {b: n for n, b in enumerate(levels)}
    Dp = np.zeros((len(pairs), len(levels)))
    Dp[np.arange(len(pairs)), [pos[b] for b in batch[i]]] = 1.0
    Dp[np.arange(len(pairs)), [pos[b] for b in batch[j]]] = -1.0
    F = Dp @ (np.linalg.pinv(Dp) @ D)
    R = D - F
    total = float((D ** 2).sum())
    values, U = np.linalg.eigh(R @ R.T)
    order = np.argsort(values)[::-1]
    values, U = np.clip(values[order], 0, None), U[:, order]
    Rp = rng.permuted(R, axis=0)
    null = np.sort(np.clip(np.linalg.eigvalsh(Rp @ Rp.T), 0, None))[::-1]
    offset = float((F ** 2).sum()) / total
    rows = []
    for k in ks:
        top = float(values[:k].sum()) / total
        rows.append(dict(celltype=celltype, k=int(k), batch_offset=offset, remainder_top_k=top,
                         rest=1 - offset - top, remainder_top_k_null=float(null[:k].sum()) / total,
                         n_pairs=len(pairs), n_batches=len(levels), n_genes=Y.shape[1]))
    nf = min(RUV_REPORTED_FACTORS, int((values > 0).sum()))
    scores = U[:, :nf] * np.sqrt(values[:nf])
    pair_type = np.array([" vs ".join(sorted((a, b))) for a, b in zip(batch[i], batch[j])])
    qc = {name: x[i] - x[j] for name, x in _sample_qc(sub).items()}
    factors = _factor_rows(celltype, "bridge differences: remainder component", scores, values[:nf] / total, qc,
                           pair_type)
    summary = ", ".join(f"top {r['k']} {r['remainder_top_k']:.3f} (null {r['remainder_top_k_null']:.3f})"
                        for r in rows)
    logger.info(f"{celltype}: within-person differences ({len(pairs)} pairs): batch offsets {offset:.3f} of the "
                f"sum of squares; of the remainder, {summary}")
    return rows, factors


def ruv_factor_rows(celltype: str, sub: pd.DataFrame, Y: np.ndarray, W: np.ndarray,
                    alpha: np.ndarray) -> List[Dict[str, Any]]:
    """ruv_factors.tsv rows for the first RUV_REPORTED_FACTORS columns of W:
    against the samples' QC, within and between batches; variance_share: the
    share of the (centered) sum of squares that factor removes."""
    nf = min(RUV_REPORTED_FACTORS, W.shape[1])
    total = float(((Y - Y.mean(axis=0)) ** 2).sum())
    shares = np.array([(W[:, f] ** 2).sum() * (alpha[f] ** 2).sum() / total for f in range(nf)])
    return _factor_rows(celltype, "ruv W", W[:, :nf], shares, _sample_qc(sub), sub["batch"].to_numpy().astype(str))


def analyse_celltype(celltype: str, obs: pd.DataFrame, counts: np.ndarray, genes: pd.DataFrame, spec,
                     rng: np.random.Generator, extra: Optional[Dict[str, list]] = None):
    """metrics, disease, de_pooled and offsets rows for one broad cell type,
    and (with umap) its samples and their UMAP coordinates per method. extra,
    if given, receives the rows of ruv_scan, ruv_diagnostic and ruv_factors.
    The ruv method has its own random numbers, so the other methods' metrics
    are as without it."""
    sel = ((obs["celltype"] == celltype) & (obs["n_cells"] >= spec["min_cells"])).to_numpy()
    sub = obs[sel].reset_index(drop=True)
    names = list(dict.fromkeys(sub["dataset"]))
    if len(names) < 2:
        logger.info(f"{celltype}: in fewer than 2 datasets; skipped")
        return [], [], pd.DataFrame(), [], None
    C = np.asarray(counts[np.flatnonzero(sel)], dtype=np.float64)
    cpm = C / sub["library_size"].to_numpy()[:, None] * 1e6
    ds = sub["dataset"].to_numpy().astype(str)
    batch = sub["batch"].to_numpy().astype(str)
    persons = sub["person_id"].to_numpy().astype(str)
    keep = np.all([cpm[ds == d].mean(axis=0) >= spec["min_cpm"] for d in names], axis=0)
    Y = np.log2(cpm[:, keep] + 1)
    genes_t = genes[keep].reset_index(drop=True)
    logger.info(f"{celltype}: {len(sub)} samples from {len(names)} datasets in {len(set(batch))} batches, "
                f"{int(keep.sum())} genes (mean CPM >= {spec['min_cpm']} in every dataset)")

    protect = sub[list(spec["protect"])]
    numeric = [c for c in spec["protect"] if c in NUMERIC]
    corrected, offsets = {"none": Y}, {}
    corrected["center"], offsets["center"] = correct_center(Y, batch)
    corrected["regress"], offsets["regress"] = correct_regress(Y, batch, protect, numeric)
    if (batch != ds).any():
        corrected["regress_dataset"], _ = correct_regress(Y, ds, protect, numeric)
    corrected["anchor"], offsets["anchor"], anchored = correct_anchor(Y, batch, persons)
    logger.info(f"{celltype}: anchored batches {sorted(anchored)}")
    ruv_rng = np.random.default_rng([int(spec["seed"]), *map(ord, celltype)])
    ruv = ruv_scan(celltype, sub, Y, spec, ruv_rng, _disease_inputs(sub, spec))

    metrics: List[Dict[str, Any]] = [dict(celltype=celltype, method="any", metric="n_genes", label="",
                                          value=int(keep.sum()), chance=np.nan, n=len(sub))]
    scored, all_pcs = score_methods(celltype, sub, corrected, spec, rng, Y=Y)
    metrics += scored
    if ruv is not None:
        corrected["ruv"] = ruv["Y"]
        removed = Y - ruv["Y"]
        offsets["ruv"] = {b: removed[batch == b].mean(axis=0) for b in np.unique(batch)}
        scored, pcs = score_methods(celltype, sub, {"ruv": ruv["Y"]}, spec, ruv_rng, replicates={"ruv": ruv["acc"]})
        metrics += scored
        all_pcs.update(pcs)
        metrics.append(dict(celltype=celltype, method="ruv", metric="ruv_k", label=ruv["choice"][0],
                            value=ruv["W"].shape[1], chance=np.nan, n=len(sub)))
    metrics += reliability_rows(celltype, "any", Y, ds, persons, spec["min_pairs"])

    disease_rows, de = disease_analysis(celltype, corrected, sub, genes_t, spec)

    offset_rows = []
    dataset_of = dict(zip(batch, ds))
    methods = ("center", "regress", "anchor") + (("ruv",) if "ruv" in offsets else ())
    for b in dict.fromkeys(batch):
        row = dict(celltype=celltype, dataset=dataset_of[b], batch=b, n_samples=int((batch == b).sum()),
                   anchored=b in anchored)
        for method in methods:
            row[f"rms_{method}"] = float(np.sqrt(np.mean(offsets[method][b] ** 2)))
        if b in anchored:
            for method in [m for m in methods if m != "anchor"]:
                row[f"r_{method}_vs_anchor"] = float(np.corrcoef(offsets[method][b], offsets["anchor"][b])[0, 1])
        offset_rows.append(row)

    if extra is not None:
        ks = sorted({int(k) for k in spec["ruv_k"]})
        diagnostic, factors = bridge_decomposition(celltype, Y, sub, ks, ruv_rng, spec["min_pairs"])
        if ruv is not None:
            factors += ruv_factor_rows(celltype, sub, Y, ruv["W"], ruv["alpha"])
            extra.setdefault("ruv_scan", []).extend(ruv["rows"])
        extra.setdefault("ruv_diagnostic", []).extend(diagnostic)
        extra.setdefault("ruv_factors", []).extend(factors)

    embedding = None
    if spec["umap"]:
        coords = {m: _umap(all_pcs[m], spec) for m in spec["umap_methods"] if m in all_pcs}
        embedding = (sub, coords)
    return metrics, disease_rows, de, offset_rows, embedding


# ---------------------------------------------------------------- figures


def method_titles(split: Dict[str, str]) -> Dict[str, str]:
    batch = "dataset" + "".join(f" + {column} in {name}" for name, column in (split or {}).items())
    return {"none": "as measured", "center": f"{batch}: means removed",
            "regress": f"{batch} regressed out", "regress_dataset": "dataset regressed out",
            "anchor": f"{batch}: anchored on replicates",
            "ruv": f"{batch}: RUV-III on replicates"}


def batch_palette(obs: pd.DataFrame, datasets: Sequence[str]) -> Dict[str, str]:
    """A fixed colour per batch, in the order of the datasets in the spec
    (within one, largest batch first); batches beyond the palette are grey."""
    order = []
    for name in datasets:
        part = obs.loc[obs["dataset"] == name, "batch"].astype(str)
        order += list(part.value_counts().index)
    return {b: CATEGORICAL[i] if i < len(CATEGORICAL) else OTHER_COLOR for i, b in enumerate(order)}


def _batch_label(batch: str) -> str:
    return batch.replace(":", ": ")


def _level_colors(values: pd.Series, levels: Sequence[str]) -> Tuple[np.ndarray, List[Tuple[str, str]]]:
    """Colours in palette order for the given levels, grey for others and
    light grey for missing values; the legend counts the samples."""
    s = values.astype(object)
    colors = {lv: CATEGORICAL[i] for i, lv in enumerate(levels)}
    per = s.map(colors).where(s.notna(), MISSING_COLOR).fillna(OTHER_COLOR).to_numpy(dtype=object)
    legend = [(f"{lv} ({int((s == lv).sum()):,})", colors[lv]) for lv in levels]
    if s.notna().any() and (~s.isin(levels) & s.notna()).any():
        legend.append(("other", OTHER_COLOR))
    if s.isna().any():
        legend.append((f"missing ({int(s.isna().sum()):,})", MISSING_COLOR))
    return per, legend


def _legend_handles(items):
    from matplotlib.lines import Line2D

    return [Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=c, markeredgewidth=0,
                   label=label) for label, c in items]


def _clean_axes(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_celltype_umaps(celltype: str, sub: pd.DataFrame, coords: Dict[str, np.ndarray],
                        palette: Dict[str, str], titles: Dict[str, str], mixing: Dict[str, float],
                        path: Path, seed: int) -> None:
    """Rows: the methods; columns: batch, dementia, sex, age, and the same
    person's samples in two datasets joined by a line."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    n = len(sub)
    order = np.random.default_rng(seed).permutation(n)
    size = float(np.clip(1.5e4 / n, 1.0, 8.0))
    batch = sub["batch"].astype(str).to_numpy()
    batch_colors = np.array([palette.get(b, OTHER_COLOR) for b in batch], dtype=object)
    present = [b for b in palette if b in set(batch)]
    batch_legend = [(f"{_batch_label(b)} ({int((batch == b).sum()):,})", palette[b]) for b in present]
    dementia_colors, dementia_legend = _level_colors(sub["dementia"], ["no", "yes"])
    sex_colors, sex_legend = _level_colors(sub["sex"], ["female", "male"])
    age = pd.to_numeric(sub["age_capped"], errors="coerce").to_numpy(np.float64)
    pairs = np.array(_replicate_pairs(sub["person_id"].astype(str).to_numpy(), sub["dataset"].astype(str).to_numpy()))
    linked = np.zeros(n, dtype=bool)
    if len(pairs):
        linked[pairs.ravel()] = True
    cmap = mcolors.LinearSegmentedColormap.from_list("ramp", SEQUENTIAL)
    lo, hi = (np.nanpercentile(age, [1, 99]) if np.isfinite(age).any() else (0, 1))

    columns = ["batch", "dementia", "sex", "age at death (capped at 90)",
               f"same person in two datasets ({len(pairs):,} pairs)"]
    methods = list(coords)
    fig, axes = plt.subplots(len(methods), len(columns), figsize=(3.5 * len(columns), 3.5 * len(methods) + 1.6),
                             squeeze=False)
    fig.patch.set_facecolor("white")
    for r, method in enumerate(methods):
        xy = coords[method]
        last = r == len(methods) - 1
        for c, column in enumerate(columns):
            ax = axes[r, c]
            if c == 0:
                ax.scatter(*xy[order].T, s=size, c=list(batch_colors[order]), linewidths=0, rasterized=True)
                if method in mixing:
                    ax.text(0.02, 0.02, f"batch share of PC variance {mixing[method]:.2f}", transform=ax.transAxes,
                            fontsize=7, color=INK)
                if last:
                    ax.legend(handles=_legend_handles(batch_legend), loc="upper center", bbox_to_anchor=(0.5, -0.02),
                              ncol=1, fontsize=7, frameon=False, handletextpad=0.2, labelcolor=INK)
            elif c in (1, 2):
                colors, legend = (dementia_colors, dementia_legend) if c == 1 else (sex_colors, sex_legend)
                missing = (colors == MISSING_COLOR)
                draw = np.concatenate([order[missing[order]], order[~missing[order]]])
                ax.scatter(*xy[draw].T, s=size, c=list(colors[draw]), linewidths=0, rasterized=True)
                if last:
                    ax.legend(handles=_legend_handles(legend), loc="upper center", bbox_to_anchor=(0.5, -0.02),
                              ncol=1, fontsize=7, frameon=False, handletextpad=0.2, labelcolor=INK)
            elif c == 3:
                ok = np.isfinite(age)
                ax.scatter(*xy[~ok].T, s=size, c=MISSING_COLOR, linewidths=0, rasterized=True)
                draw = order[ok[order]]
                points = ax.scatter(*xy[draw].T, s=size, c=age[draw], cmap=cmap, vmin=lo, vmax=max(hi, lo + 1),
                                    linewidths=0, rasterized=True)
                if last:
                    bar = fig.colorbar(points, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03)
                    bar.ax.tick_params(labelsize=7, colors=INK)
                    bar.outline.set_visible(False)
            else:
                ax.scatter(*xy[~linked].T, s=size, c=MISSING_COLOR, linewidths=0, rasterized=True)
                if len(pairs):
                    ax.add_collection(LineCollection(xy[pairs], colors=INK, linewidths=0.4, alpha=0.35))
                    ax.scatter(*xy[linked].T, s=size, c=list(batch_colors[linked]), linewidths=0, rasterized=True)
            if r == 0:
                ax.set_title(column, fontsize=9, color=INK)
            if c == 0:
                ax.set_ylabel(titles.get(method, method), fontsize=9, color=INK)
            _clean_axes(ax)
    fig.suptitle(f"{celltype}: pseudobulk UMAP, {n:,} samples (dataset x donor), "
                 f"{sub['person_id'].nunique():,} persons", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def plot_overview(embeddings: Dict[str, Tuple[pd.DataFrame, Dict[str, np.ndarray]]], palette: Dict[str, str],
                  titles: Dict[str, str], mixing: Dict[Tuple[str, str], float], path: Path, seed: int) -> None:
    """Columns: the cell types; rows: the methods; coloured by batch."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    celltypes = list(embeddings)
    methods = list(next(iter(embeddings.values()))[1])
    fig, axes = plt.subplots(len(methods), len(celltypes),
                             figsize=(3.0 * len(celltypes), 3.0 * len(methods) + 1.0), squeeze=False)
    fig.patch.set_facecolor("white")
    seen = set()
    for c, celltype in enumerate(celltypes):
        sub, coords = embeddings[celltype]
        batch = sub["batch"].astype(str).to_numpy()
        seen |= set(batch)
        colors = np.array([palette.get(b, OTHER_COLOR) for b in batch], dtype=object)
        order = np.random.default_rng(seed).permutation(len(sub))
        size = float(np.clip(1.5e4 / len(sub), 1.0, 8.0))
        for r, method in enumerate(methods):
            ax = axes[r, c]
            ax.scatter(*coords[method][order].T, s=size * 0.6, c=list(colors[order]), linewidths=0, rasterized=True)
            if (celltype, method) in mixing:
                ax.text(0.02, 0.02, f"batch share {mixing[(celltype, method)]:.2f}", transform=ax.transAxes,
                        fontsize=7, color=INK)
            if r == 0:
                ax.set_title(f"{celltype} ({len(sub):,} samples)", fontsize=9, color=INK)
            if c == 0:
                ax.set_ylabel(titles.get(method, method), fontsize=9, color=INK)
            _clean_axes(ax)
    legend = [(_batch_label(b), col) for b, col in palette.items() if b in seen]
    fig.legend(handles=_legend_handles(legend), loc="lower center", ncol=len(legend), fontsize=8, frameon=False,
               handletextpad=0.2, columnspacing=1.2, labelcolor=INK)
    fig.suptitle("Pooled pseudobulk UMAP per cell type, coloured by batch\n"
                 "batch share: share of the PC variance that the batch explains", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- run


def _sample_table(obs: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows = []
    for (d, b, t), part in obs.groupby(["dataset", "batch", "celltype"], sort=False):
        kept = part[part["n_cells"] >= min_cells]
        rows.append(dict(dataset=d, batch=b, celltype=t, n_samples=len(part), n_kept=len(kept),
                         median_nuclei=float(kept["n_cells"].median()) if len(kept) else np.nan,
                         n_persons_in_other_datasets=int((kept.get("n_datasets", pd.Series(dtype=float)) > 1).sum())))
    return pd.DataFrame(rows)


def run_pooled_pseudobulk(spec_path: Union[str, Path], overwrite: bool = False, rebuild: bool = False) -> Path:
    """overwrite: replace the analysis tables (an existing pseudobulk.h5ad is
    reused); rebuild: also rebuild pseudobulk.h5ad from the count files."""
    import anndata as ad

    spec = load_spec(spec_path)
    out = Path(spec["output_dir"])
    h5ad = Path(spec["pseudobulk"]) if spec["pseudobulk"] else out / "pseudobulk.h5ad"
    if spec["pseudobulk"]:
        if rebuild:
            raise ValueError(f"{spec_path} reads pseudobulk {h5ad} (never written); --rebuild needs a spec "
                             "without 'pseudobulk'")
        if not h5ad.exists():
            raise FileNotFoundError(f"pseudobulk: {h5ad} does not exist")
    if (out / "metrics.tsv").exists() and not (overwrite or rebuild):
        raise FileExistsError(f"{out / 'metrics.tsv'} exists; pass --overwrite to replace the outputs")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(spec["seed"])
    names = list(spec["datasets"])
    table = read_donor_table(spec["donor_table"], names)

    if h5ad.exists() and not rebuild:
        logger.info(f"Reusing {h5ad} (--rebuild to rebuild it)")
        adata = ad.read_h5ad(h5ad)
    else:
        settings = {"min_genes": spec["min_genes"], "min_cells_per_donor": 1}
        datasets = [prepare_dataset(k, {**cfg, "batch": [], "covariates": {}}, settings)
                    for k, cfg in spec["datasets"].items()]
        gene_ids, columns = shared_genes(datasets)
        groups, matrices = [], []
        for ds in datasets:
            g, m = pseudobulk_dataset(ds, columns[ds.name], spec["chunk_size"])
            groups.append(g)
            matrices.append(m)
            ds.cells = None  # free the cell table
        obs = pd.concat(groups, ignore_index=True)
        obs["donor"] = obs["donor"].astype(str)
        matrix = np.vstack(matrices)
        obs = attach_persons(obs, table)
        matched = obs.index.to_numpy()
        obs = obs.reset_index(drop=True)
        var = pd.concat([pd.DataFrame({"gene_id": gene_ids}), _symbols(datasets, columns)], axis=1)
        var.index = var["gene_id"]
        obs.index = (obs["dataset"] + ":" + obs["donor"] + ":" + obs["celltype"]).to_numpy()
        for c in obs.columns:
            if obs[c].dtype == object:
                obs[c] = obs[c].astype("category")
        adata = ad.AnnData(sp.csr_matrix(np.rint(matrix[matched]).astype(np.int32)), obs=obs, var=var)
        adata.uns["provenance_pooled_pseudobulk"] = provenance_record(
            "pooled_pseudobulk", {"spec": str(spec_path), "min_genes": spec["min_genes"]})
        adata.write_h5ad(h5ad, compression="gzip")
        logger.info(f"Wrote {h5ad}: {adata.n_obs} samples x {adata.n_vars} genes")

    obs = adata.obs.copy()
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(object).where(obs[c].notna(), np.nan)
    obs = obs.reset_index(drop=True)
    obs["batch"] = assign_batches(obs, spec["batch_split"])
    for name in spec["batch_split"] or {}:
        persons = obs[obs["dataset"] == name].groupby("batch")["person_id"].nunique()
        logger.info(f"[{name}] persons per batch: {persons.to_dict()}")
    genes = adata.var.reset_index(drop=True)
    counts = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

    _sample_table(obs, spec["min_cells"]).to_csv(out / "samples.tsv", sep="\t", index=False)
    genes.to_csv(out / "genes.tsv", sep="\t", index=False)
    metrics, disease, de, offsets, embeddings = [], [], [], [], {}
    extra: Dict[str, list] = {}
    for celltype in sorted(obs["celltype"].unique()):
        sel = np.flatnonzero((obs["celltype"] == celltype).to_numpy())
        m, d, g, o, e = analyse_celltype(celltype, obs.iloc[sel].reset_index(drop=True), counts[sel].toarray(),
                                         genes, spec, rng, extra=extra)
        metrics += m
        disease += d
        offsets += o
        if len(g):
            de.append(g)
        if e is not None:
            embeddings[celltype] = e
    metrics = pd.DataFrame(metrics, columns=METRIC_COLUMNS)
    metrics.to_csv(out / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(disease).to_csv(out / "disease.tsv", sep="\t", index=False)
    (pd.concat(de, ignore_index=True) if de else pd.DataFrame()).to_csv(out / "de_pooled.tsv", sep="\t",
                                                                          index=False)
    pd.DataFrame(offsets).to_csv(out / "offsets.tsv", sep="\t", index=False)
    for name in ("ruv_scan", "ruv_diagnostic", "ruv_factors"):
        pd.DataFrame(extra.get(name, [])).to_csv(out / f"{name}.tsv", sep="\t", index=False)

    if embeddings:
        figures = out / "umap"
        figures.mkdir(exist_ok=True)
        palette = batch_palette(obs, names)
        titles = method_titles(spec["batch_split"])
        label = "batch" if (obs["batch"] != obs["dataset"]).any() else "dataset"
        share = metrics[(metrics["metric"] == "pc_r2") & (metrics["label"] == label)]
        mixing = {(r.celltype, r.method): float(r.value) for r in share.itertuples()}
        for celltype, (sub, coords) in embeddings.items():
            plot_celltype_umaps(celltype, sub, coords, palette, titles,
                                {m: v for (t, m), v in mixing.items() if t == celltype},
                                figures / f"{celltype}.png", spec["seed"])
        plot_overview(embeddings, palette, titles, mixing, figures / "overview.png", spec["seed"])
        logger.info(f"Wrote the UMAPs to {figures}")
    with open(out / "provenance.json", "w") as f:
        json.dump(provenance_record("pooled_pseudobulk", {"spec": str(spec_path)}) | {"spec": spec}, f,
                  indent=2, default=str)
    logger.info(f"Wrote the aggregate tables to {out}")
    return out


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("spec", help="YAML spec, e.g. config/pooled_pseudobulk.yaml")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace the output tables (reusing an existing pseudobulk.h5ad)")
    parser.add_argument("--rebuild", action="store_true", help="Also rebuild pseudobulk.h5ad from the counts")
    args = parser.parse_args(argv)
    from ..utils.logging import setup_logging

    setup_logging(level="INFO")
    run_pooled_pseudobulk(args.spec, overwrite=args.overwrite, rebuild=args.rebuild)


if __name__ == "__main__":
    main()
