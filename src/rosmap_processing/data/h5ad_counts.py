"""Write a count h5ad from a large one without loading its matrix.

PsychAD's h5ad holds 6.3 million cells; loading the matrix takes more memory
than the conversion needs. `copy_counts` streams a CSR matrix (X, a layer, or
raw/X) into a new file in row chunks, checks that every value is a count,
stores it as int32, and writes the given obs and var next to it.
"""

from pathlib import Path
from typing import Optional, Union

import h5py
import numpy as np
import pandas as pd

from ..utils.counts import count_non_integer_values
from ..utils.logging import get_logger

try:
    from anndata.io import write_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import write_elem

logger = get_logger(__name__)


def csr_node(f: h5py.File, key: str) -> h5py.Group:
    """The CSR group at `key` ("X", "raw/X" or "layers/<name>")."""
    if key not in f:
        raise KeyError(f"{key} not found")
    node = f[key]
    if not isinstance(node, h5py.Group) or node.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError(f"{key} is not stored as a CSR matrix")
    return node


def copy_counts(
    src: Union[str, Path],
    dst: Union[str, Path],
    obs: pd.DataFrame,
    var: pd.DataFrame,
    key: str = "X",
    uns: Optional[dict] = None,
    chunk_rows: int = 200_000,
) -> None:
    """Write dst with X = the counts at src[key] (int32, CSR), and obs, var, uns.

    obs and var must have one row per cell and gene of the matrix, in its
    order. Raises ValueError, and leaves no output, if any stored value is
    negative or not a whole number.
    """
    dst = Path(dst)
    tmp = dst.with_name(dst.name + ".partial")
    with h5py.File(src, "r") as f:
        node = csr_node(f, key)
        n_obs, n_vars = (int(d) for d in node.attrs["shape"])
        if (len(obs), len(var)) != (n_obs, n_vars):
            raise ValueError(f"obs/var have {len(obs)} x {len(var)} rows for a {n_obs} x {n_vars} matrix")
        indptr = node["indptr"][()].astype(np.int64)
        nnz = int(indptr[-1])
        logger.info(f"Copying {key}: {n_obs} cells x {n_vars} genes, {nnz} stored values "
                    f"({node['data'].dtype}) in chunks of {chunk_rows} cells")
        dst.parent.mkdir(parents=True, exist_ok=True)
        try:
            with h5py.File(tmp, "w") as g:
                g.attrs["encoding-type"] = "anndata"
                g.attrs["encoding-version"] = "0.1.0"
                X = g.create_group("X")
                X.attrs["encoding-type"] = "csr_matrix"
                X.attrs["encoding-version"] = "0.1.0"
                X.attrs["shape"] = np.array([n_obs, n_vars], dtype=np.int64)
                chunks = (min(max(nnz, 1), 1 << 20),)
                data = X.create_dataset("data", shape=(nnz,), dtype=np.int32, chunks=chunks, compression="gzip")
                indices = X.create_dataset("indices", shape=(nnz,), dtype=np.int32, chunks=chunks,
                                           compression="gzip")
                X.create_dataset("indptr", data=indptr)
                maximum = 0
                for start in range(0, n_obs, chunk_rows):
                    lo, hi = indptr[start], indptr[min(start + chunk_rows, n_obs)]
                    values = node["data"][lo:hi]
                    bad = count_non_integer_values(values)
                    if bad:
                        raise ValueError(f"{key}: {bad} of the values of cells {start}+ are not counts")
                    if values.size:
                        maximum = max(maximum, values.max().item())
                    if maximum > np.iinfo(np.int32).max:
                        raise ValueError("Counts exceed the int32 range")
                    data[lo:hi] = values.astype(np.int32)
                    indices[lo:hi] = node["indices"][lo:hi]
                    if (start // chunk_rows) % 10 == 0:
                        logger.info(f"  {min(start + chunk_rows, n_obs)}/{n_obs} cells")
                write_elem(g, "obs", obs)
                write_elem(g, "var", var)
                write_elem(g, "uns", uns or {})
                for empty in ("obsm", "varm", "obsp", "varp", "layers"):
                    write_elem(g, empty, {})
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    tmp.replace(dst)
    logger.info(f"Wrote {dst} ({dst.stat().st_size / 1e9:.1f} GB), counts up to {maximum}")
