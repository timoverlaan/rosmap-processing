"""Read the HDF5 files written by src/convert_R/rdata_to_h5.R.

The R script stores the CSC arrays of a genes x cells dgCMatrix (`data`,
`indices`, `indptr`), its `shape` as cells x genes, and the `genes` and `cells`
names. Read as CSR, the arrays are the cells x genes matrix, so nothing is
transposed. Cell metadata, if the R file had any, is in CSV files next to it.
"""

from pathlib import Path
from typing import Union

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..utils.logging import get_logger
from .genes import make_unique

logger = get_logger(__name__)


def read_matrix_h5(path: Union[str, Path]) -> ad.AnnData:
    """AnnData with X = the stored matrix (cells x genes), obs and var indexed by name."""
    with h5py.File(path, "r") as f:
        n_cells, n_genes = (int(d) for d in f["shape"][()])
        genes = f["genes"].asstr()[()]
        cells = f["cells"].asstr()[()]
        if (len(cells), len(genes)) != (n_cells, n_genes):
            raise ValueError(f"{path}: shape {n_cells} x {n_genes}, but {len(cells)} cells "
                             f"and {len(genes)} genes")
        indptr = f["indptr"][()].astype(np.int64)
        if len(indptr) != n_cells + 1 or indptr[-1] != f["data"].shape[0]:
            raise ValueError(f"{path}: indptr does not fit {n_cells} cells and the stored values")
        X = sp.csr_matrix((f["data"][()], f["indices"][()], indptr), shape=(n_cells, n_genes))
    logger.info(f"Read {Path(path).name}: {n_cells} cells x {n_genes} genes, {X.nnz} stored values")
    # Duplicate symbols get R's make.unique suffixes (".1", ...). Most
    # DiverseCohorts matrices have them already; gex_mat_Glut3 and _Glut4 keep
    # 20 duplicates bare, with the same genes in the same order.
    unique = make_unique(list(genes), ".")
    renamed = sum(a != b for a, b in zip(genes, unique))
    if renamed:
        logger.info(f"{renamed} duplicate gene names made unique as R's make.unique does")
        genes = unique
    return ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(cells)),
                      var=pd.DataFrame(index=pd.Index(genes)))
