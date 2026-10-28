"""Tests for `data prepare` and the category-fix guard against dropping counts."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

from rosmap_processing.data.category_fix import fix_categories
from rosmap_processing.data.prepare import prepare_counts, prepare_file


def _seuratdisk_like(normalized: bool, n_cells=40, n_genes=20, seed=0):
    """AnnData shaped like SeuratDisk output: X = Seurat 'data', raw = counts."""
    rng = np.random.default_rng(seed)
    counts = sp.csr_matrix(rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float64))
    obs = pd.DataFrame(
        {
            "projid": np.repeat([101, 102], n_cells // 2).astype(np.int32),
            "cell_type_high_resolution": ["A", "B"] * (n_cells // 2),
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])
    adata = ad.AnnData(counts.copy(), obs=obs, var=var)
    adata.raw = adata.copy()
    if normalized:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata, counts


@pytest.mark.parametrize("normalized", [True, False])
def test_prepare_takes_counts_from_raw(normalized):
    adata, counts = _seuratdisk_like(normalized)
    out = prepare_counts(adata, categorical_columns=["cell_type_high_resolution"])
    assert out.raw is None
    assert out.X.dtype == np.int32
    np.testing.assert_array_equal(out.X.toarray(), counts.toarray())
    assert out.obs["cell_type_high_resolution"].dtype.name == "category"


def test_prepare_fails_without_counts():
    adata, _ = _seuratdisk_like(normalized=True)
    adata = ad.AnnData(adata.X, obs=adata.obs, var=adata.var)  # no .raw
    with pytest.raises(ValueError, match="does not hold raw counts"):
        prepare_counts(adata)


def test_prepare_adds_metadata():
    adata, _ = _seuratdisk_like(normalized=True)
    meta = pd.DataFrame({"projid": [101, 102], "msex": [0, 1]})
    out = prepare_counts(adata, metadata=meta)
    assert list(out.obs.loc[["c0", "c39"], "msex"]) == [0, 1]


def test_prepare_file_never_overwrites(tmp_path):
    adata, _ = _seuratdisk_like(normalized=True)
    src = tmp_path / "converted.h5ad"
    adata.write_h5ad(src)
    with pytest.raises(ValueError, match="in place"):
        prepare_file(src, src)

    dst = tmp_path / "counts.h5ad"
    prepare_file(src, dst)
    assert dst.exists()
    assert "provenance_prepare" in ad.read_h5ad(dst).uns
    with pytest.raises(FileExistsError):
        prepare_file(src, dst)


def test_category_fix_refuses_to_drop_the_only_counts():
    adata, _ = _seuratdisk_like(normalized=True)
    with pytest.raises(ValueError, match="Refusing to remove .raw"):
        fix_categories(adata, columns=[], remove_raw=True)


def test_category_fix_drops_raw_when_x_holds_counts():
    adata, _ = _seuratdisk_like(normalized=False)
    out = fix_categories(adata, columns=[], remove_raw=True)
    assert out.raw is None
    assert out.X.dtype == np.int32
