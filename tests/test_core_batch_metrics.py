"""Tests for the batch-effect metrics, on synthetic embeddings with known structure."""

import numpy as np
import pandas as pd
import pytest

from rosmap_processing.core import batch_metrics as bm


@pytest.fixture
def rng():
    return np.random.default_rng(0)


def _clusters(rng, n_per=200, shift=5.0, n_dims=10):
    """Two groups; shift moves group 1 along the first axis."""
    X = rng.normal(size=(2 * n_per, n_dims))
    labels = np.repeat(["a", "b"], n_per)
    X[labels == "b", 0] += shift
    return X, labels


def test_r2_per_column_categorical_and_numeric(rng):
    X, labels = _clusters(rng)
    r2 = bm.r2_per_column(X, labels)
    assert r2[0] > 0.8 and np.all(r2[1:] < 0.05)
    x = rng.normal(size=500)
    Y = np.column_stack([2 * x + rng.normal(scale=0.1, size=500), rng.normal(size=500)])
    r2 = bm.r2_per_column(Y, x, numeric=True)
    assert r2[0] > 0.95 and r2[1] < 0.05


def test_pc_regression_separated_vs_random(rng):
    X, labels = _clusters(rng, shift=8)
    assert bm.pc_regression(X, labels) > 0.5
    assert bm.pc_regression(X, rng.permutation(labels)) < 0.05


def test_permute_labels_between_keeps_donor_level_structure(rng):
    donors = np.repeat([f"d{i}" for i in range(10)], 30)
    batch = np.repeat(["b1"] * 4 + ["b2"] * 6, 30)
    null = bm.permute_labels(batch, donors, rng)
    per_donor = pd.Series(null).groupby(donors).nunique()
    assert (per_donor == 1).all()
    assert sorted(pd.Series(null).groupby(donors).first().value_counts()) == [4, 6]


def test_permute_labels_within_keeps_each_donors_mix(rng):
    donors = np.repeat(["d1", "d2"], 50)
    labels = np.array((["x"] * 20 + ["y"] * 30) * 2)
    null = bm.permute_labels(labels, donors, rng, within=True)
    for d in ("d1", "d2"):
        assert sorted(null[donors == d]) == sorted(labels[donors == d])


def test_mixing_ratio_and_knn_rejection(rng):
    X, labels = _clusters(rng, shift=10)
    knn = bm.knn_indices(X, 15)
    separated, _ = bm.mixing_ratio(knn, labels)
    mixed, _ = bm.mixing_ratio(knn, rng.permutation(labels))
    assert separated < 0.6
    assert 0.9 < mixed < 1.1
    assert bm.knn_chi2_rejection(knn, labels) > 0.9
    assert bm.knn_chi2_rejection(knn, rng.permutation(labels)) < 0.2


def test_silhouettes(rng):
    X, labels = _clusters(rng, shift=10)
    assert bm.label_asw(X, labels, rng) > 0.7
    celltypes = np.zeros(len(labels))
    assert bm.batch_asw(X, labels, celltypes, rng) < 0.5
    assert bm.batch_asw(X, rng.permutation(labels), celltypes, rng) > 0.9


def test_partial_r2_finds_the_term(rng):
    n = 200
    covs = pd.DataFrame({"batch": rng.choice(["a", "b", "c"], n), "age": rng.normal(70, 8, n)})
    effect = covs["batch"].map({"a": 0.0, "b": 2.0, "c": -2.0}).to_numpy()
    Y = np.column_stack([effect + rng.normal(size=n), rng.normal(size=n)])
    X, terms = bm.design_matrix(covs, numeric=["age"])
    r2 = bm.partial_r2(Y, X, terms["batch"])
    assert r2[0] > 0.6 and r2[1] < 0.1
    assert bm.partial_r2(Y, X, terms["age"])[0] < 0.1


def test_classify_above_chance(rng):
    labels = np.repeat(["a", "b"], 40)
    X = rng.normal(size=(80, 30))
    X[labels == "b", :5] += 2
    accuracy, chance, n_levels = bm.classify(X, labels, rng, n_null=2)
    assert n_levels == 2 and accuracy > 0.9 and chance < 0.75


def test_cramers_v_and_eta_squared(rng):
    a = np.repeat(["x", "y"], 100)
    v, p = bm.cramers_v(a, a)
    assert v > 0.95 and p < 1e-10
    v, _ = bm.cramers_v(a, rng.permutation(a))
    assert v < 0.2
    values = np.where(a == "x", 0.0, 3.0) + rng.normal(size=200)
    eta2, p = bm.eta_squared(values, a)
    assert eta2 > 0.6 and p < 1e-10


def test_paired_decomposition(rng):
    n, g = 50, 40
    donor = rng.normal(scale=2, size=(n, g))
    A = donor + rng.normal(scale=0.3, size=(n, g))
    B = donor + 1.0 + rng.normal(scale=0.3, size=(n, g))  # same donors, shifted dataset
    out = bm.paired_decomposition(A, B)
    assert out["share_donor"] > 0.7
    assert out["donor_correlation"] > 0.9
    assert 0.02 < out["share_dataset"] < 0.3
    assert abs(sum(out[k] for k in ("share_dataset", "share_donor", "share_residual")) - 1) < 0.1
