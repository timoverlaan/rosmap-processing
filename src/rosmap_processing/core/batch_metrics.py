"""Batch-effect metrics on a cell embedding or on donor-level pseudobulk.

Every function takes arrays and labels and returns numbers, so the metrics can
be tested on synthetic data. Labels must not contain missing values (drop
those cells first).

Cell level (an embedding, e.g. PCs of log-normalized expression):
    pc_regression     share of the embedding's variance explained by a label
    mixing_ratio      local label diversity among kNN neighbours (LISI-style)
    knn_chi2_rejection  share of neighbourhoods whose label mix differs from
                      the overall one (kBET-style)
    label_asw / batch_asw  silhouette widths (scIB conventions)

Most technical batches here contain whole donors, so a batch label also
carries the differences between donors. `permute_labels` gives the null for
that: the same label shuffled between donors (or within donors, for a label
that varies inside every donor, such as the dataset of a donor seen twice).

Donor level (pseudobulk, one row per donor or donor x dataset):
    partial_r2        per gene, the variance a term explains beyond the others
    classify          cross-validated balanced accuracy of predicting a label
    cramers_v / eta_squared  association between two donor-level variables
    paired_decomposition  dataset vs donor variance for donors in two datasets
"""

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats


def _codes(labels) -> Tuple[np.ndarray, int]:
    codes, uniques = pd.factorize(pd.Series(np.asarray(labels, dtype=object)), sort=True)
    if (codes < 0).any():
        raise ValueError("labels contain missing values")
    return codes, len(uniques)


def _indicator(codes: np.ndarray, n_levels: int) -> sp.csr_matrix:
    """(n_levels, n_cells) 0/1 matrix: row l marks the cells with code l."""
    n = len(codes)
    return sp.csr_matrix((np.ones(n), (codes, np.arange(n))), shape=(n_levels, n))


def r2_per_column(Y: np.ndarray, labels, numeric: bool = False) -> np.ndarray:
    """R² of each column of Y on a label: group means (one-way ANOVA), or for a
    numeric label a straight-line fit."""
    Y = np.asarray(Y, dtype=np.float64)
    centered = Y - Y.mean(axis=0)
    total = (centered ** 2).sum(axis=0)
    if numeric:
        x = np.asarray(labels, dtype=np.float64)
        x = x - x.mean()
        denom = (x ** 2).sum()
        if denom == 0:
            return np.zeros(Y.shape[1])
        between = (x @ centered) ** 2 / denom
    else:
        codes, n_levels = _codes(labels)
        indicator = _indicator(codes, n_levels)
        sizes = np.asarray(indicator.sum(axis=1)).ravel()
        sums = indicator @ centered
        between = (sums ** 2 / sizes[:, None]).sum(axis=0)
    return np.divide(between, total, out=np.zeros_like(total), where=total > 0)


def pc_regression(pcs: np.ndarray, labels, numeric: bool = False) -> float:
    """Share of the embedding's variance explained by a label: the R² of each
    PC weighted by that PC's variance (scIB's principal component regression)."""
    weights = np.asarray(pcs, dtype=np.float64).var(axis=0)
    r2 = r2_per_column(pcs, labels, numeric=numeric)
    return float((weights * r2).sum() / weights.sum()) if weights.sum() > 0 else float("nan")


def permute_labels(labels, groups, rng: np.random.Generator, within: bool = False) -> np.ndarray:
    """A null version of `labels` that keeps the group (donor) structure.

    within=False: the labels are moved between groups. Group g gets the label
    mix of group pi(g) for a random permutation pi, so a donor-level label
    (one value per donor) is reassigned to other donors with the same group
    sizes per label.
    within=True: the labels are shuffled among the cells of each group, so
    each group keeps its own label mix.
    """
    labels = np.asarray(labels, dtype=object)
    group_codes, n_groups = _codes(groups)
    order = np.argsort(group_codes, kind="stable")
    members = np.split(order, np.cumsum(np.bincount(group_codes, minlength=n_groups))[:-1])
    out = np.empty_like(labels)
    if within:
        for idx in members:
            out[idx] = labels[rng.permutation(idx)]
        return out
    perm = rng.permutation(n_groups)
    for g, source in enumerate(perm):
        pool = labels[members[source]]
        if len(set(pool)) == 1:
            out[members[g]] = pool[0]
        else:
            out[members[g]] = rng.choice(pool, size=len(members[g]), replace=True)
    return out


def knn_indices(X: np.ndarray, k: int) -> np.ndarray:
    """(n_cells, k) indices of each cell's k nearest neighbours (itself excluded)."""
    from sklearn.neighbors import NearestNeighbors

    k = min(k, len(X) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    return nn.kneighbors(X, return_distance=False)[:, 1:]


def _neighbour_counts(knn: np.ndarray, codes: np.ndarray, n_levels: int) -> np.ndarray:
    """(n_cells, n_levels): how many of each cell's neighbours have each label."""
    n_cells = knn.shape[0]
    flat = codes[knn] + n_levels * np.arange(n_cells)[:, None]
    return np.bincount(flat.ravel(), minlength=n_cells * n_levels).reshape(n_cells, n_levels)


def mixing_ratio(knn: np.ndarray, labels) -> Tuple[float, float]:
    """Local label diversity relative to random mixing, and the raw local index.

    For each cell, the inverse Simpson index of the labels among its k
    neighbours (as in LISI, with equal neighbour weights). Divided by the
    inverse of the expected Simpson index if the neighbours were drawn at
    random from all cells, so about 1 means as mixed as random (a little above
    for rare labels), and 1 / (that value) means every neighbourhood holds a
    single label.
    """
    codes, n_levels = _codes(labels)
    k = knn.shape[1]
    counts = _neighbour_counts(knn, codes, n_levels)
    local = 1.0 / ((counts / k) ** 2).sum(axis=1)
    p = np.bincount(codes, minlength=n_levels) / len(codes)
    simpson = (p ** 2).sum()
    # E[sum (c/k)^2] for c ~ Multinomial(k, p)
    expected = 1.0 / (simpson + (1 - simpson) / k)
    return float(local.mean() / expected), float(local.mean())


def knn_chi2_rejection(knn: np.ndarray, labels, alpha: float = 0.05) -> float:
    """Share of neighbourhoods (the cell and its k neighbours) whose label
    counts differ from the overall proportions by a chi-square test, as in
    kBET. 0 = every neighbourhood looks like the whole; 1 = none does."""
    codes, n_levels = _codes(labels)
    if n_levels < 2:
        return float("nan")
    counts = _neighbour_counts(knn, codes, n_levels)
    counts[np.arange(len(codes)), codes] += 1
    size = counts.sum(axis=1, keepdims=True)
    p = np.bincount(codes, minlength=n_levels) / len(codes)
    expected = size * p
    chi2 = ((counts - expected) ** 2 / expected).sum(axis=1)
    pvalues = stats.chi2.sf(chi2, df=n_levels - 1)
    return float((pvalues < alpha).mean())


def _subsample(n: int, max_cells: int, rng: np.random.Generator) -> np.ndarray:
    return np.sort(rng.choice(n, size=max_cells, replace=False)) if n > max_cells else np.arange(n)


def label_asw(X: np.ndarray, labels, rng: np.random.Generator, max_cells: int = 10000) -> float:
    """Silhouette width of a label, scaled to 0-1 as in scIB: (mean + 1) / 2.
    1 = labels form tight separate clusters; 0.5 = no structure."""
    from sklearn.metrics import silhouette_samples

    codes, n_levels = _codes(labels)
    idx = _subsample(len(codes), max_cells, rng)
    if len(np.unique(codes[idx])) < 2:
        return float("nan")
    return float((silhouette_samples(X[idx], codes[idx]).mean() + 1) / 2)


def batch_asw(
    X: np.ndarray, batch, celltypes, rng: np.random.Generator, max_cells: int = 10000
) -> float:
    """Batch silhouette as in scIB: within each cell type, mean of 1 - |s| over
    cells, then averaged over cell types. 1 = batches fully mixed; 0 = every
    batch its own cluster."""
    from sklearn.metrics import silhouette_samples

    batch_codes, _ = _codes(batch)
    type_codes, n_types = _codes(celltypes)
    scores = []
    for t in range(n_types):
        idx = np.flatnonzero(type_codes == t)
        idx = idx[_subsample(len(idx), max_cells, rng)]
        b = batch_codes[idx]
        if len(np.unique(b)) < 2 or len(idx) < 3:
            continue
        scores.append((1 - np.abs(silhouette_samples(X[idx], b))).mean())
    return float(np.mean(scores)) if scores else float("nan")


def design_matrix(covariates: pd.DataFrame, numeric: Sequence[str] = ()) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Intercept + one-hot (drop first) for categorical columns, standardized
    numeric columns. Returns the matrix and, per column, its column indices."""
    blocks, terms, start = [np.ones((len(covariates), 1))], {}, 1
    for name in covariates.columns:
        if name in numeric:
            x = covariates[name].to_numpy(dtype=np.float64)
            sd = x.std()
            block = ((x - x.mean()) / (sd if sd > 0 else 1.0))[:, None]
        else:
            block = pd.get_dummies(covariates[name].astype(str), drop_first=True).to_numpy(dtype=np.float64)
        blocks.append(block)
        terms[name] = np.arange(start, start + block.shape[1])
        start += block.shape[1]
    return np.hstack(blocks), terms


def _rss(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return ((Y - X @ beta) ** 2).sum(axis=0)


def partial_r2(Y: np.ndarray, X: np.ndarray, term_columns: np.ndarray) -> np.ndarray:
    """Per column of Y: (RSS without the term - RSS with it) / RSS without it."""
    Y = np.asarray(Y, dtype=np.float64)
    keep = np.setdiff1d(np.arange(X.shape[1]), term_columns)
    rss_full = _rss(X, Y)
    rss_reduced = _rss(X[:, keep], Y)
    return np.divide(rss_reduced - rss_full, rss_reduced,
                     out=np.zeros_like(rss_full), where=rss_reduced > 0)


def classify(
    X: np.ndarray, labels, rng: np.random.Generator, n_splits: int = 5,
    n_components: int = 20, n_null: int = 5,
) -> Tuple[float, float, int]:
    """Cross-validated balanced accuracy of predicting labels from X (logistic
    regression on PCs), and the same with shuffled labels (chance).

    Levels with fewer than n_splits rows are left out. Returns
    (accuracy, chance, number of levels used)."""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    labels = pd.Series(np.asarray(labels, dtype=object))
    sizes = labels.value_counts()
    keep = labels.isin(sizes.index[sizes >= n_splits]).to_numpy()
    X, y = np.asarray(X)[keep], labels[keep].to_numpy()
    n_levels = len(np.unique(y))
    if n_levels < 2:
        return float("nan"), float("nan"), n_levels
    n_comp = min(n_components, X.shape[1], len(y) - len(y) // n_splits - 1)

    def score(target):
        model = make_pipeline(StandardScaler(), PCA(n_components=n_comp),
                              LogisticRegression(max_iter=2000, class_weight="balanced"))
        folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=int(rng.integers(2**31)))
        return balanced_accuracy_score(target, cross_val_predict(model, X, target, cv=folds))

    accuracy = score(y)
    chance = float(np.mean([score(rng.permutation(y)) for _ in range(n_null)]))
    return float(accuracy), chance, n_levels


def cramers_v(a, b) -> Tuple[float, float]:
    """Cramér's V (bias-corrected) and chi-square p-value of two categorical variables."""
    table = pd.crosstab(pd.Series(np.asarray(a, dtype=object)), pd.Series(np.asarray(b, dtype=object)))
    if min(table.shape) < 2:
        return float("nan"), float("nan")
    chi2, p, _, _ = stats.chi2_contingency(table, correction=False)
    n = table.to_numpy().sum()
    r, c = table.shape
    phi2 = max(0.0, chi2 / n - (r - 1) * (c - 1) / (n - 1))
    r_corr, c_corr = r - (r - 1) ** 2 / (n - 1), c - (c - 1) ** 2 / (n - 1)
    denom = min(r_corr - 1, c_corr - 1)
    return (float(np.sqrt(phi2 / denom)) if denom > 0 else float("nan")), float(p)


def eta_squared(values, groups) -> Tuple[float, float]:
    """Share of a numeric variable's variance between groups, and the
    Kruskal-Wallis p-value."""
    values = np.asarray(values, dtype=np.float64)
    r2 = r2_per_column(values[:, None], groups)[0]
    codes, n_levels = _codes(groups)
    if n_levels < 2:
        return float("nan"), float("nan")
    samples = [values[codes == g] for g in range(n_levels)]
    return float(r2), float(stats.kruskal(*samples).pvalue)


def paired_decomposition(A: np.ndarray, B: np.ndarray) -> Dict[str, float]:
    """For the same donors (rows) measured in two datasets, per gene (column):
    the variance components of a two-way layout without replication, as
    shares of their sum: between datasets, between donors (shared by both
    datasets), and donor x dataset plus noise; and the correlation of donor
    differences between the datasets. Returns medians over genes.

    The components come from the expected mean squares (donor: MS_donor =
    var_residual + 2 var_donor), not the raw sums of squares, whose donor part
    is half the within-dataset variance even when donors share nothing."""
    A, B = np.asarray(A, dtype=np.float64), np.asarray(B, dtype=np.float64)
    n = A.shape[0]
    grand = (A.mean(axis=0) + B.mean(axis=0)) / 2
    total = ((A - grand) ** 2).sum(axis=0) + ((B - grand) ** 2).sum(axis=0)
    ss_dataset = n * ((A.mean(axis=0) - grand) ** 2 + (B.mean(axis=0) - grand) ** 2)
    ss_donor = 2 * (((A + B) / 2 - grand) ** 2).sum(axis=0)
    ms_residual = (total - ss_dataset - ss_donor) / (n - 1)
    var_dataset = np.maximum(ss_dataset - ms_residual, 0) / (2 * n)
    var_donor = np.maximum(ss_donor / (n - 1) - ms_residual, 0) / 2
    var_residual = np.maximum(ms_residual, 0)
    var_total = var_dataset + var_donor + var_residual
    ok = var_total > 0
    a, b = A - A.mean(axis=0), B - B.mean(axis=0)
    denom = np.sqrt((a ** 2).sum(axis=0) * (b ** 2).sum(axis=0))
    corr = np.divide((a * b).sum(axis=0), denom, out=np.full(A.shape[1], np.nan), where=denom > 0)
    return {
        "share_dataset": float(np.median(var_dataset[ok] / var_total[ok])),
        "share_donor": float(np.median(var_donor[ok] / var_total[ok])),
        "share_residual": float(np.median(var_residual[ok] / var_total[ok])),
        "donor_correlation": float(np.nanmedian(corr)),
    }


def null_summary(values: Sequence[float]) -> Tuple[float, float]:
    values = np.asarray([v for v in values if np.isfinite(v)], dtype=np.float64)
    if not len(values):
        return float("nan"), float("nan")
    return float(values.mean()), float(values.std(ddof=1)) if len(values) > 1 else float("nan")


def z_score(observed: float, null_mean: float, null_sd: float) -> Optional[float]:
    if not (np.isfinite(observed) and np.isfinite(null_mean) and np.isfinite(null_sd)) or null_sd == 0:
        return float("nan")
    return float((observed - null_mean) / null_sd)
