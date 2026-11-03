"""Tests for the streamed count copy (data/h5ad_counts.py)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from rosmap_processing.data.h5ad_counts import copy_counts


def _source(tmp_path, values, layer=False):
    X = sp.csr_matrix(values.astype(np.float32))
    obs = pd.DataFrame({"donor": [f"d{i % 2}" for i in range(X.shape[0])]},
                       index=[f"c{i}" for i in range(X.shape[0])])
    var = pd.DataFrame(index=[f"G{j}" for j in range(X.shape[1])])
    adata = ad.AnnData(X=sp.csr_matrix(X.shape, dtype=np.float32) if layer else X, obs=obs, var=var)
    if layer:
        adata.layers["counts"] = X
    path = tmp_path / "src.h5ad"
    adata.write_h5ad(path)
    return path, obs, var


def test_copy_counts_streams_matrix_with_new_obs(tmp_path):
    rng = np.random.default_rng(0)
    values = rng.poisson(0.5, size=(23, 7))
    src, obs, var = _source(tmp_path, values, layer=True)
    obs = obs.assign(sex=["female", "male"] * 11 + ["female"])
    copy_counts(src, tmp_path / "out.h5ad", obs, var, key="layers/counts", uns={"note": "x"}, chunk_rows=5)
    out = ad.read_h5ad(tmp_path / "out.h5ad")
    assert out.X.dtype == np.int32 and sp.isspmatrix_csr(out.X)
    np.testing.assert_array_equal(out.X.toarray(), values)
    assert out.obs["sex"].tolist() == obs["sex"].tolist()
    assert out.var_names.tolist() == var.index.tolist()
    assert out.uns["note"] == "x"


def test_copy_counts_refuses_non_counts_and_leaves_no_output(tmp_path):
    values = np.array([[1.0, 0.0], [0.0, 2.5], [3.0, 0.0]])
    src, obs, var = _source(tmp_path, values)
    with pytest.raises(ValueError, match="not counts"):
        copy_counts(src, tmp_path / "out.h5ad", obs, var, chunk_rows=1)
    assert not (tmp_path / "out.h5ad").exists()
    assert not (tmp_path / "out.h5ad.partial").exists()


def test_copy_counts_checks_obs_length(tmp_path):
    src, obs, var = _source(tmp_path, np.eye(3))
    with pytest.raises(ValueError, match="rows"):
        copy_counts(src, tmp_path / "out.h5ad", obs.iloc[:2], var)
