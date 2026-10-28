"""Read the counts and cell metadata of an .h5Seurat file, without R.

SeuratDisk cannot load these files with the SeuratObject in the container
(`GetAssayData(slot = )` is defunct since SeuratObject 5.0), and its h5ad
conversion puts Seurat's `data` slot in X and loses factor labels. The parts we
need are plain HDF5:

- `assays/<assay>/<slot>`: a dgCMatrix (features x cells, compressed by cell)
  as `data`/`indices`/`indptr` with a `dims` attribute. Read as CSR, that is
  the cells x features matrix.
- `assays/<assay>/features` and `cell.names`: the row and column names.
- `meta.data/<column>`: one dataset per column, or a group with `levels` and
  1-based `values` for an R factor; the `colnames` attribute gives the order.
"""

from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..utils.logging import get_logger

logger = get_logger(__name__)

R_NA_INTEGER = np.iinfo(np.int32).min


def _strings(dataset: h5py.Dataset) -> np.ndarray:
    return dataset.asstr()[()]


def _read_column(node, name: str):
    if isinstance(node, h5py.Group):
        if set(node.keys()) != {"levels", "values"}:
            raise ValueError(f"meta.data/{name}: unsupported group with {sorted(node.keys())}")
        levels = _strings(node["levels"])
        values = node["values"][()].astype(np.int64)
        codes = np.where(values == R_NA_INTEGER, -1, values - 1)
        if codes.min(initial=0) < -1 or codes.max(initial=-1) >= len(levels):
            raise ValueError(f"meta.data/{name}: factor codes outside its {len(levels)} levels")
        return pd.Categorical.from_codes(codes, categories=levels)
    if h5py.check_string_dtype(node.dtype) is not None:
        return _strings(node)
    values = node[()]
    if values.dtype.kind == "i" and (values == R_NA_INTEGER).any():
        values = np.where(values == R_NA_INTEGER, np.nan, values)
    return values


def _read_meta_data(group: h5py.Group, cells: np.ndarray) -> pd.DataFrame:
    if "_index" in group and not np.array_equal(_strings(group["_index"]), cells):
        raise ValueError("meta.data/_index does not match cell.names")
    if "colnames" in group.attrs:
        columns = [c.decode() if isinstance(c, bytes) else str(c) for c in np.atleast_1d(group.attrs["colnames"])]
    else:
        columns = sorted(k for k in group.keys() if k != "_index")
    obs = pd.DataFrame(index=pd.Index(cells))
    for name in columns:
        values = _read_column(group[name], name)
        if len(values) != len(cells):
            raise ValueError(f"meta.data/{name} has {len(values)} rows, expected {len(cells)}")
        obs[name] = values
    return obs


def read_h5seurat(path: Path, assay: str = "RNA", slot: str = "counts") -> ad.AnnData:
    """AnnData with X = the assay's `slot` matrix (cells x features), obs = meta.data."""
    with h5py.File(path, "r") as f:
        cells = _strings(f["cell.names"])
        assay_group = f["assays"][assay]
        features = _strings(assay_group["features"])
        matrix = assay_group[slot]
        if not isinstance(matrix, h5py.Group):
            raise ValueError(f"assays/{assay}/{slot} is not a sparse matrix")
        n_features, n_cells = (int(d) for d in matrix.attrs["dims"])
        if (n_features, n_cells) != (len(features), len(cells)):
            raise ValueError(
                f"assays/{assay}/{slot} has dims {n_features} x {n_cells}, but there are "
                f"{len(features)} features and {len(cells)} cells"
            )
        X = sp.csr_matrix(
            (matrix["data"][()], matrix["indices"][()], matrix["indptr"][()]),
            shape=(n_cells, n_features),
        )
        obs = _read_meta_data(f["meta.data"], cells)

    logger.info(f"Read assays/{assay}/{slot}: {n_cells} cells x {n_features} features, {X.nnz} stored values")
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(features)))
