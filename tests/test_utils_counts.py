"""Unit tests for the count checks and normalization helper."""

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

from rosmap_processing.utils.counts import (
    assert_integer_counts,
    check_h5ad_counts,
    count_non_integer_values,
    is_integer_counts,
    log_normalize,
    to_int32_counts,
    total_counts,
)


def _counts(n_cells=50, n_genes=30, seed=0):
    rng = np.random.default_rng(seed)
    return sp.csr_matrix(rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32))


def test_integer_counts_detected_sparse_and_dense():
    X = _counts()
    assert is_integer_counts(X)
    assert is_integer_counts(X.toarray())
    assert is_integer_counts(to_int32_counts(X))


def test_non_integer_values_are_counted_anywhere_in_the_matrix():
    X = _counts(n_cells=2000)
    # Corrupt only the last row: a sampled check of the first rows would miss it.
    X = X.tolil()
    X[1999, 0] = 0.5
    X[1999, 1] = -1.0
    X = X.tocsr()
    assert count_non_integer_values(X) == 2
    with pytest.raises(ValueError, match="does not hold raw counts"):
        assert_integer_counts(X)


def test_log_normalized_matrix_is_rejected():
    adata = ad.AnnData(_counts())
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    assert not is_integer_counts(adata.X)


def test_log_normalize_uses_given_totals_on_a_gene_subset():
    X = _counts(n_genes=40)
    totals = total_counts(X)
    expected = ad.AnnData(X.copy())
    sc.pp.normalize_total(expected, target_sum=1e4)
    sc.pp.log1p(expected)

    subset = slice(0, 10)
    got = log_normalize(X[:, subset], totals, 1e4)
    np.testing.assert_allclose(got.toarray(), expected.X[:, subset].toarray(), rtol=1e-5)


def test_log_normalize_keeps_empty_cells_at_zero():
    X = sp.csr_matrix(np.array([[0, 0], [1, 3]], dtype=np.float32))
    got = log_normalize(X, total_counts(X), 1e4).toarray()
    assert np.all(got[0] == 0)
    np.testing.assert_allclose(got[1], np.log1p(np.array([1, 3]) / 4 * 1e4), rtol=1e-6)


def test_check_h5ad_counts(tmp_path):
    counts = _counts()
    adata = ad.AnnData(counts.copy())
    adata.layers["counts"] = counts.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    path = tmp_path / "x.h5ad"
    adata.write_h5ad(path)

    x_stats = check_h5ad_counts(path, key="X")
    layer_stats = check_h5ad_counts(path, key="layers/counts")
    assert not x_stats["is_counts"]
    assert layer_stats["is_counts"]
    assert layer_stats["stored_values"] == counts.nnz
    assert layer_stats["max"] == counts.max()
