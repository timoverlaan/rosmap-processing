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
3. Per broad cell type, on the samples with at least min_cells nuclei and
   the genes with a mean CPM of at least min_cpm in every dataset, log2-CPM
   and four ways to remove the dataset effect:
     none     as measured
     center   each dataset's mean per gene removed
     regress  per gene, a model of dataset + the `protect` covariates (sex,
              age, dementia); only the dataset terms are removed, so mean
              differences that the covariates explain are kept
     anchor   per gene, the dataset offsets that best explain the differences
              between the same persons measured in two datasets (purely
              technical); datasets without such persons are centered
   and how each does:
     dataset mixing   share of the PC variance that the dataset explains
                      (pc_regression) and the balanced accuracy of a dataset
                      classifier
     replicates       for persons measured in two datasets: how often their
                      profile in one dataset is the closest (Pearson, over
                      genes, of the deviation from the pooled mean) to their
                      own profile in the other, of all that dataset's samples
                      (top1, chance = 1 / candidates), the median rank as a
                      share of the candidates, and the mean correlation with
                      their own vs other persons' profiles. The anchor
                      offsets are cross-fitted: a person's offsets are
                      estimated without their pairs (n_folds folds)
     disease          the reference model: one sample per person (the one
                      with the most nuclei), log2-CPM ~ disease + adjust +
                      dataset (de_pooled.tsv); per dataset alone; and per
                      correction without the dataset term. Share of genes at
                      FDR < 0.05 and the correlation of the t statistics with
                      the reference model's
     offsets          per dataset, the size of each method's offsets and
                      how the center and regress offsets correlate with the
                      anchor (technical) ones

Aggregate outputs in output_dir: samples.tsv, genes.tsv, metrics.tsv,
disease.tsv, de_pooled.tsv (per gene), offsets.tsv, provenance.json.

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
}
REQUIRED = ("output_dir", "donor_table", "datasets")
DATASET_KEYS = {"counts", "layer", "gene_map", "donor", "celltype"}
METHODS = ("none", "center", "regress", "anchor")
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
                    with_dataset: bool) -> np.ndarray:
    blocks = [np.ones((len(frame), 1)), (frame[disease].astype(str) == positive).to_numpy(np.float64)[:, None]]
    if adjust:
        C, _ = bm.design_matrix(frame[adjust].reset_index(drop=True), [c for c in adjust if c in NUMERIC])
        blocks.append(C[:, 1:])
    if with_dataset:
        blocks.append(pd.get_dummies(frame["dataset"].astype(str), drop_first=True).to_numpy(np.float64))
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
    X = _disease_design(frame, disease, positive, adjust, with_dataset=True)
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

    rows.append(summary("pooled: dataset as covariate", t_ref, p, len(frame), n_pos))
    for method, Yc in corrected.items():
        Xm = _disease_design(frame, disease, positive, adjust, with_dataset=False)
        _, t, pv = _ols(Yc[first], Xm, 1)
        rows.append(summary(f"pooled: {method}, no dataset term", t, pv, len(frame), n_pos))
    for name in dict.fromkeys(sub["dataset"]):
        sel = np.flatnonzero(usable & (sub["dataset"] == name).to_numpy())
        f = sub.iloc[sel]
        npos = int((f[disease].astype(str) == positive).sum())
        if npos < 3 or len(f) - npos < 3:
            continue
        _, t, pv = _ols(corrected["none"][sel], _disease_design(f, disease, positive, adjust, False), 1)
        rows.append(summary(f"{name} alone", t, pv, len(f), npos))
    return rows, de


def analyse_celltype(celltype: str, obs: pd.DataFrame, counts: np.ndarray, genes: pd.DataFrame, spec,
                     rng: np.random.Generator):
    """metrics, disease, de_pooled and offsets rows for one broad cell type."""
    sel = ((obs["celltype"] == celltype) & (obs["n_cells"] >= spec["min_cells"])).to_numpy()
    sub = obs[sel].reset_index(drop=True)
    names = list(dict.fromkeys(sub["dataset"]))
    if len(names) < 2:
        logger.info(f"{celltype}: in fewer than 2 datasets; skipped")
        return [], [], pd.DataFrame(), []
    C = np.asarray(counts[np.flatnonzero(sel)], dtype=np.float64)
    cpm = C / sub["library_size"].to_numpy()[:, None] * 1e6
    ds = sub["dataset"].to_numpy().astype(str)
    persons = sub["person_id"].to_numpy().astype(str)
    keep = np.all([cpm[ds == d].mean(axis=0) >= spec["min_cpm"] for d in names], axis=0)
    Y = np.log2(cpm[:, keep] + 1)
    genes_t = genes[keep].reset_index(drop=True)
    logger.info(f"{celltype}: {len(sub)} samples from {len(names)} datasets, {int(keep.sum())} genes "
                f"(mean CPM >= {spec['min_cpm']} in every dataset)")

    protect = sub[list(spec["protect"])]
    corrected, offsets = {"none": Y}, {}
    corrected["center"], offsets["center"] = correct_center(Y, ds)
    corrected["regress"], offsets["regress"] = correct_regress(
        Y, ds, protect, [c for c in spec["protect"] if c in NUMERIC])
    corrected["anchor"], offsets["anchor"], anchored = correct_anchor(Y, ds, persons)
    logger.info(f"{celltype}: anchored datasets {sorted(anchored)}")

    metrics: List[Dict[str, Any]] = [dict(celltype=celltype, method="any", metric="n_genes", label="",
                                          value=int(keep.sum()), chance=np.nan, n=len(sub))]
    for method, Yc in corrected.items():
        pcs = _pca(Yc - Yc.mean(axis=0), spec["n_pcs"], spec["seed"])
        metrics.append(dict(celltype=celltype, method=method, metric="dataset_pc_r2", label="dataset",
                            value=bm.pc_regression(pcs, ds), chance=np.nan, n=len(sub)))
        accuracy, chance, n_levels = bm.classify(pcs, ds, rng)
        metrics.append(dict(celltype=celltype, method=method, metric="dataset_accuracy", label="dataset",
                            value=accuracy, chance=chance, n=len(sub)))
        if method != "anchor":
            acc = replicate_stats(Yc, ds, persons)
        else:  # cross-fitted: a person's offsets come from the other persons' pairs
            replicated = sorted({persons[i] for i, _ in _replicate_pairs(persons, ds)})
            folds = np.arange(len(replicated)) % spec["n_folds"]
            rng.shuffle(folds)
            acc = {}
            for k in range(spec["n_folds"]):
                fold = {p for p, f in zip(replicated, folds) if f == k}
                if not fold:
                    continue
                Yk, _, _ = correct_anchor(Y, ds, persons, exclude=frozenset(fold))
                for key, v in replicate_stats(Yk, ds, persons, queries=fold,
                                              center=corrected["anchor"].mean(axis=0)).items():
                    target = acc.setdefault(key, {k2: [] for k2 in v})
                    for k2, values in v.items():
                        target[k2] += values
        metrics += _replicate_rows(celltype, method, acc, spec["min_pairs"])
    for a, b in combinations(names, 2):
        ia, ib = np.flatnonzero(ds == a), np.flatnonzero(ds == b)
        pa, pb = pd.Series(ia, index=persons[ia]), pd.Series(ib, index=persons[ib])
        both = sorted(set(pa.index) & set(pb.index))
        if len(both) >= spec["min_pairs"]:
            value = bm.shared_reliability(Y[pa[both].to_numpy()], Y[pb[both].to_numpy()])
            metrics.append(dict(celltype=celltype, method="any", metric="cross_dataset_reliability",
                                label=f"{a} vs {b}", value=value, chance=0.0, n=len(both)))

    disease_rows, de = disease_analysis(celltype, corrected, sub, genes_t, spec)

    offset_rows = []
    for d in names:
        row = dict(celltype=celltype, dataset=d, n_samples=int((ds == d).sum()), anchored=d in anchored)
        for method in ("center", "regress", "anchor"):
            row[f"rms_{method}"] = float(np.sqrt(np.mean(offsets[method][d] ** 2)))
        if d in anchored:
            for method in ("center", "regress"):
                row[f"r_{method}_vs_anchor"] = float(np.corrcoef(offsets[method][d], offsets["anchor"][d])[0, 1])
        offset_rows.append(row)
    return metrics, disease_rows, de, offset_rows


# ---------------------------------------------------------------- run


def _sample_table(obs: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    rows = []
    for (d, t), part in obs.groupby(["dataset", "celltype"], sort=False):
        kept = part[part["n_cells"] >= min_cells]
        rows.append(dict(dataset=d, celltype=t, n_samples=len(part), n_kept=len(kept),
                         median_nuclei=float(kept["n_cells"].median()) if len(kept) else np.nan,
                         n_persons_in_other_datasets=int((kept.get("n_datasets", pd.Series(dtype=float)) > 1).sum())))
    return pd.DataFrame(rows)


def run_pooled_pseudobulk(spec_path: Union[str, Path], overwrite: bool = False, rebuild: bool = False) -> Path:
    """overwrite: replace the analysis tables (an existing pseudobulk.h5ad is
    reused); rebuild: also rebuild pseudobulk.h5ad from the count files."""
    import anndata as ad

    spec = load_spec(spec_path)
    out = Path(spec["output_dir"])
    h5ad = out / "pseudobulk.h5ad"
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
    genes = adata.var.reset_index(drop=True)
    counts = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

    _sample_table(obs, spec["min_cells"]).to_csv(out / "samples.tsv", sep="\t", index=False)
    genes.to_csv(out / "genes.tsv", sep="\t", index=False)
    metrics, disease, de, offsets = [], [], [], []
    for celltype in sorted(obs["celltype"].unique()):
        sel = np.flatnonzero((obs["celltype"] == celltype).to_numpy())
        m, d, g, o = analyse_celltype(celltype, obs.iloc[sel].reset_index(drop=True), counts[sel].toarray(),
                                      genes, spec, rng)
        metrics += m
        disease += d
        offsets += o
        if len(g):
            de.append(g)
    pd.DataFrame(metrics, columns=METRIC_COLUMNS).to_csv(out / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(disease).to_csv(out / "disease.tsv", sep="\t", index=False)
    (pd.concat(de, ignore_index=True) if de else pd.DataFrame()).to_csv(out / "de_pooled.tsv", sep="\t",
                                                                          index=False)
    pd.DataFrame(offsets).to_csv(out / "offsets.tsv", sep="\t", index=False)
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
