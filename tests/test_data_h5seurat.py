"""Tests for reading .h5Seurat files directly (synthetic files, same layout as SeuratDisk's)."""

import anndata as ad
import h5py
import numpy as np
import pytest
import scipy.sparse as sp

from rosmap_processing.data.h5seurat import R_NA_INTEGER, read_h5seurat
from rosmap_processing.data.prepare import prepare_file

STR = h5py.string_dtype("utf-8")


def _write_h5seurat(path, counts, cells, features, dims=None):
    """Write `counts` (cells x features) the way SeuratDisk stores a dgCMatrix."""
    csc = sp.csc_matrix(counts.T)  # features x cells, compressed by cell
    with h5py.File(path, "w") as f:
        f.attrs["active.assay"] = "RNA"
        f.create_dataset("cell.names", data=np.array(cells, dtype=object), dtype=STR)
        rna = f.create_group("assays/RNA")
        rna.create_dataset("features", data=np.array(features, dtype=object), dtype=STR)
        for slot, values in [("counts", csc.data), ("data", np.log1p(csc.data))]:
            g = rna.create_group(slot)
            g.attrs["dims"] = np.array(dims or csc.shape, dtype=np.int32)
            g.create_dataset("data", data=values.astype(np.float64))
            g.create_dataset("indices", data=csc.indices.astype(np.int32))
            g.create_dataset("indptr", data=csc.indptr.astype(np.int32))

        md = f.create_group("meta.data")
        md.attrs["colnames"] = np.array(["projid", "cell_type", "batch", "score"], dtype=object)
        md.create_dataset("_index", data=np.array(cells, dtype=object), dtype=STR)
        projid = np.arange(len(cells), dtype=np.int32) % 2 + 101
        projid[-1] = R_NA_INTEGER
        md.create_dataset("projid", data=projid)
        factor = md.create_group("cell_type")
        factor.create_dataset("levels", data=np.array(["Ast", "Mic"], dtype=object), dtype=STR)
        values = np.arange(len(cells), dtype=np.int32) % 2 + 1  # R factors are 1-based
        values[0] = R_NA_INTEGER
        factor.create_dataset("values", data=values)
        md.create_dataset("batch", data=np.array(["b1"] * len(cells), dtype=object), dtype=STR)
        md.create_dataset("score", data=np.linspace(0, 1, len(cells)))


def _example(n_cells=12, n_genes=7, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float64)
    cells = [f"AAAC-{i}" for i in range(n_cells)]
    features = [f"G{i}" for i in range(n_genes)]
    return counts, cells, features


def test_reads_counts_and_metadata(tmp_path):
    counts, cells, features = _example()
    path = tmp_path / "x.h5Seurat"
    _write_h5seurat(path, counts, cells, features)

    adata = read_h5seurat(path)
    assert sp.isspmatrix_csr(adata.X)
    np.testing.assert_array_equal(adata.X.toarray(), counts)
    assert list(adata.obs_names) == cells
    assert list(adata.var_names) == features
    assert list(adata.obs.columns) == ["projid", "cell_type", "batch", "score"]

    cell_type = adata.obs["cell_type"]
    assert cell_type.dtype.name == "category"
    assert list(cell_type.cat.categories) == ["Ast", "Mic"]
    assert cell_type.isna().iloc[0]
    assert list(cell_type.iloc[1:3]) == ["Mic", "Ast"]

    projid = adata.obs["projid"]
    assert projid.isna().iloc[-1]
    assert projid.iloc[0] == 101
    assert projid.astype(str).iloc[0] == "101"
    assert adata.obs["batch"].iloc[0] == "b1"


def test_reads_other_slot(tmp_path):
    counts, cells, features = _example()
    path = tmp_path / "x.h5Seurat"
    _write_h5seurat(path, counts, cells, features)
    np.testing.assert_allclose(read_h5seurat(path, slot="data").X.toarray(), np.log1p(counts))


def test_dims_mismatch_fails(tmp_path):
    counts, cells, features = _example()
    path = tmp_path / "x.h5Seurat"
    _write_h5seurat(path, counts, cells, features, dims=(len(features) + 1, len(cells)))
    with pytest.raises(ValueError, match="dims"):
        read_h5seurat(path)


def test_prepare_file_reads_h5seurat(tmp_path):
    counts, cells, features = _example()
    src = tmp_path / "x.h5Seurat"
    _write_h5seurat(src, counts, cells, features)

    dst = tmp_path / "counts.h5ad"
    prepare_file(src, dst, categorical_columns=["cell_type"])
    out = ad.read_h5ad(dst)
    assert out.X.dtype == np.int32
    np.testing.assert_array_equal(out.X.toarray(), counts)
    assert out.obs["cell_type"].dtype.name == "category"
    assert out.raw is None
