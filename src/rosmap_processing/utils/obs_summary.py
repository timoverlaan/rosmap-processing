"""Inventory of the obs columns of h5ad files, without any values.

Per column: dtype, number of distinct values, % missing, and how it relates
to the donor: constant within every donor (donor level, e.g. a diagnosis, or a
sequencing batch that holds whole donors) or varying within donors (cell
level, e.g. a cell type). Never reads X and never prints or writes a value,
so the output is safe to share for restricted data.
"""

from pathlib import Path
from typing import List, Optional, Sequence, Union

import h5py
import pandas as pd

from ..utils.logging import get_logger

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem

logger = get_logger(__name__)

DONOR_CANDIDATES = ("Donor ID", "projid", "individualID")


def _matrices(f: h5py.File) -> List[str]:
    """'X', 'layers/<name>', 'raw/X' with their encoding and shape."""
    out = []
    nodes = [("X", f.get("X"))] + [(f"layers/{k}", v) for k, v in f.get("layers", {}).items()]
    if "raw" in f and "X" in f["raw"]:
        nodes.append(("raw/X", f["raw"]["X"]))
    for name, node in nodes:
        if node is None:
            continue
        if isinstance(node, h5py.Group):
            shape = tuple(int(s) for s in node.attrs.get("shape", ()))
            out.append(f"{name}: {node.attrs.get('encoding-type', '?')} {shape}, "
                       f"{node['data'].shape[0]} stored values, {node['data'].dtype}")
        else:
            out.append(f"{name}: dense {node.shape}, {node.dtype}")
    return out


def summarize_obs(path: Union[str, Path], donor_column: Optional[str] = None) -> pd.DataFrame:
    """One row per obs column (see module docstring)."""
    path = Path(path)
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        index_name = obs.attrs["_index"]
        columns = [c for c in obs.keys() if c not in (index_name, "__categories")]
        for line in _matrices(f):
            logger.info(f"{path}: {line}")
        logger.info(f"{path}: uns keys {sorted(f['uns'].keys()) if 'uns' in f else []}")
        if donor_column is None:
            donor_column = next((c for c in DONOR_CANDIDATES if c in columns), None)
        donors = None
        if donor_column is not None:
            donors = pd.Series(read_elem(obs[donor_column])).astype(str).to_numpy()
            n_donors = len(set(donors))
            logger.info(f"{path}: {len(donors)} cells, donor column '{donor_column}', {n_donors} donors")
        rows = []
        for col in columns:
            values = pd.Series(read_elem(obs[col]))
            missing = values.isna() | values.astype(str).str.strip().str.lower().isin(("", "na", "nan", "none"))
            present = values[~missing]
            row = {
                "file": str(path), "column": col, "dtype": str(values.dtype),
                "n_unique": int(present.nunique()),
                "pct_missing": round(100 * float(missing.mean()), 2),
            }
            if donors is not None and col != donor_column:
                per_donor = present.groupby(donors[~missing.to_numpy()]).nunique()
                row["donors_with_values"] = int(len(per_donor))
                row["max_values_per_donor"] = int(per_donor.max()) if len(per_donor) else 0
                row["pct_donors_with_several"] = round(100 * float((per_donor > 1).mean()), 2) \
                    if len(per_donor) else 0.0
                if row["n_unique"] <= 1:
                    row["level"] = "constant"
                elif row["max_values_per_donor"] <= 1:
                    row["level"] = "donor"
                else:
                    row["level"] = "cell"
                row["unique_per_donor_ratio"] = round(row["n_unique"] / max(n_donors, 1), 3)
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_obs_files(
    paths: Sequence[Union[str, Path]], output: Union[str, Path], donor_column: Optional[str] = None
) -> pd.DataFrame:
    table = pd.concat([summarize_obs(p, donor_column) for p in paths], ignore_index=True)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    with pd.option_context("display.max_rows", None, "display.width", 200):
        logger.info(f"obs summary ({len(table)} columns) -> {output}\n"
                    f"{table.drop(columns=['file']).to_string(index=False)}")
    return table
