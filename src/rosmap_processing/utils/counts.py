"""Checks and transforms for raw count matrices.

Every step that expects raw counts calls `assert_integer_counts` on the *whole*
matrix. An earlier version sampled the first 500 cells, which missed that some
per-class files held log-normalized values (see the README, "Data history").
"""

from typing import Optional

import numpy as np
import scipy.sparse as sp

# Values checked per chunk, to bound the temporary arrays on 1e9+ nonzeros.
_CHUNK = 50_000_000


def _values(X) -> np.ndarray:
    """The stored values of X: `.data` for sparse, flattened array for dense."""
    if sp.issparse(X):
        return X.data
    return np.asarray(X).ravel()


def count_non_integer_values(X) -> int:
    """Number of stored values in X that are negative or not whole numbers."""
    values = _values(X)
    if np.issubdtype(values.dtype, np.integer):
        return int(np.count_nonzero(values < 0))
    bad = 0
    for start in range(0, values.size, _CHUNK):
        chunk = values[start:start + _CHUNK]
        bad += int(np.count_nonzero((chunk < 0) | (chunk != np.floor(chunk))))
    return bad


def is_integer_counts(X) -> bool:
    """True if every stored value in X is a non-negative whole number."""
    return count_non_integer_values(X) == 0


def assert_integer_counts(X, name: str = "X") -> None:
    """Raise ValueError unless every stored value in X is a non-negative whole number."""
    bad = count_non_integer_values(X)
    if bad:
        total = _values(X).size
        raise ValueError(
            f"{name} does not hold raw counts: {bad:,} of {total:,} stored values "
            f"are negative or non-integer ({100 * bad / max(total, 1):.1f}%). "
            f"Expected raw UMI counts."
        )


def to_int32_counts(X):
    """Cast a verified count matrix to int32 (sparse stays sparse)."""
    values = _values(X)
    if values.size and values.max() > np.iinfo(np.int32).max:
        raise ValueError("Counts exceed the int32 range")
    if sp.issparse(X):
        X = X.copy()
        X.data = X.data.astype(np.int32)
        return X
    return np.asarray(X).astype(np.int32)


def log_normalize(
    counts,
    size_factors: np.ndarray,
    target_sum: float,
    dtype=np.float32,
) -> sp.csr_matrix:
    """log1p(counts / size_factor * target_sum), row by row.

    `size_factors` are per-cell totals computed *before* any gene subsetting,
    so the result for a gene subset equals the subset of the full-matrix
    normalization (scanpy's normalize_total computes totals over whatever
    genes are present, which is why this does not use it).

    Parameters
    ----------
    counts : sparse or dense (n_cells, n_genes)
    size_factors : (n_cells,) total counts per cell; cells with 0 stay 0
    target_sum : counts per cell after scaling (e.g. 1e4)
    """
    X = sp.csr_matrix(counts, dtype=dtype, copy=True)
    size_factors = np.asarray(size_factors, dtype=np.float64).ravel()
    if size_factors.shape[0] != X.shape[0]:
        raise ValueError(
            f"size_factors has {size_factors.shape[0]} entries for {X.shape[0]} cells"
        )
    scale = np.divide(
        target_sum, size_factors, out=np.zeros_like(size_factors), where=size_factors > 0
    ).astype(dtype)
    X.data *= np.repeat(scale, np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def total_counts(X, dtype: Optional[type] = np.float64) -> np.ndarray:
    """Per-cell total counts as a 1-D array."""
    return np.asarray(X.sum(axis=1), dtype=dtype).ravel()


def check_h5ad_counts(path, key: str = "X") -> dict:
    """Aggregate statistics of a matrix stored in an h5ad, read in chunks.

    `key` is the HDF5 path of the matrix: "X", "raw/X" or "layers/<name>".
    Returns counts only (no values), so the result is safe to log for
    restricted data.
    """
    import h5py

    with h5py.File(path, "r") as f:
        if key not in f:
            raise KeyError(f"{key} not found in {path}")
        node = f[key]
        values = node["data"] if isinstance(node, h5py.Group) else node
        shape = tuple(node.attrs["shape"]) if isinstance(node, h5py.Group) else node.shape
        n_stored = values.shape[0] if values.ndim == 1 else int(np.prod(values.shape))
        stats = {"key": key, "shape": shape, "dtype": str(values.dtype), "stored_values": n_stored,
                 "non_integer": 0, "negative": 0, "min": None, "max": None}
        step = _CHUNK if values.ndim == 1 else max(1, _CHUNK // max(values.shape[1], 1))
        for start in range(0, values.shape[0], step):
            chunk = np.asarray(values[start:start + step]).ravel()
            if not chunk.size:
                continue
            stats["negative"] += int(np.count_nonzero(chunk < 0))
            if not np.issubdtype(chunk.dtype, np.integer):
                stats["non_integer"] += int(np.count_nonzero(chunk != np.floor(chunk)))
            cmin, cmax = chunk.min().item(), chunk.max().item()
            stats["min"] = cmin if stats["min"] is None else min(stats["min"], cmin)
            stats["max"] = cmax if stats["max"] is None else max(stats["max"], cmax)
    stats["is_counts"] = stats["non_integer"] == 0 and stats["negative"] == 0
    return stats
