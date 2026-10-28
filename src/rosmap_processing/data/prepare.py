"""Turn a SeuratDisk-converted h5ad into a raw-count h5ad.

SeuratDisk's `Convert()` writes Seurat's `data` slot to X and the `counts` slot
to `.raw`. Whether `data` holds counts or log-normalized values depends on how
the authors saved each object: in the Mathys 2023 files, Immune_cells,
Inhibitory_neurons and OPCs are log-normalized; in the Green 2024 files, cux2+
and inhibitory are. This step always takes the counts (`.raw` when present),
checks the whole matrix is integer counts, adds clinical metadata, and writes a
new file. It never modifies its input.
"""

from pathlib import Path
from typing import List, Optional

import anndata as ad
import pandas as pd

from ..utils.counts import assert_integer_counts, to_int32_counts
from ..utils.logging import get_logger
from ..utils.provenance import provenance_record
from .category_fix import fix_categories
from .metadata import add_metadata

logger = get_logger(__name__)

# Columns SeuratDisk sometimes leaves in var/raw.var that h5ad reserves.
_RESERVED_VAR_COLUMNS = ["_index", "_categories"]


def promote_raw_counts(adata: ad.AnnData) -> ad.AnnData:
    """Return an AnnData whose X is `.raw.X` (all raw genes), without `.raw`.

    If there is no `.raw`, X is returned unchanged and must already hold counts.
    obs, obsm and uns are kept; var columns of X's genes are carried over.
    """
    if adata.raw is None:
        logger.info("No .raw present; using X as the count matrix")
        return adata

    raw_var = adata.raw.var.drop(columns=_RESERVED_VAR_COLUMNS, errors="ignore").copy()
    var = raw_var
    for col in adata.var.columns:
        if col in _RESERVED_VAR_COLUMNS or col in var.columns:
            continue
        var[col] = adata.var[col].reindex(var.index)

    logger.info(
        f"Using .raw as counts: {adata.raw.shape[1]} genes "
        f"(X had {adata.n_vars}); dropping X (Seurat 'data' slot)"
    )
    return ad.AnnData(
        X=adata.raw.X,
        obs=adata.obs.copy(),
        var=var,
        obsm=dict(adata.obsm),
        uns=dict(adata.uns),
    )


def prepare_counts(
    adata: ad.AnnData,
    categorical_columns: Optional[List[str]] = None,
    metadata: Optional[pd.DataFrame] = None,
) -> ad.AnnData:
    """Counts into X (int32), categorical obs columns, optional clinical metadata."""
    adata = promote_raw_counts(adata)
    assert_integer_counts(adata.X, name="counts matrix")
    adata.X = to_int32_counts(adata.X)
    logger.info(f"Verified raw counts in X: {adata.n_obs} cells x {adata.n_vars} genes (int32)")

    adata = fix_categories(
        adata, columns=categorical_columns or [], remove_raw=False, cast_counts_to_int=False
    )
    if metadata is not None:
        adata = add_metadata(adata, metadata)
    return adata


def prepare_file(
    input_path: Path,
    output_path: Path,
    categorical_columns: Optional[List[str]] = None,
    metadata_path: Optional[Path] = None,
    overwrite: bool = False,
) -> None:
    """Read a converted h5ad, prepare counts + metadata, write to a new file."""
    input_path = Path(input_path)
    output_path = Path(output_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if output_path.resolve() == input_path.resolve():
        raise ValueError("Output must differ from input; this step does not modify files in place")
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path} (pass --overwrite to replace it)")

    logger.info(f"Loading {input_path}")
    adata = ad.read_h5ad(input_path)
    logger.info(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes, raw present: {adata.raw is not None}")

    metadata = pd.read_csv(metadata_path) if metadata_path else None
    adata = prepare_counts(adata, categorical_columns=categorical_columns, metadata=metadata)
    adata.uns["provenance_prepare"] = provenance_record(
        "data prepare",
        {
            "input": input_path,
            "categorical_columns": ",".join(categorical_columns or []),
            "metadata": metadata_path,
        },
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(output_path.name + ".partial")
    adata.write_h5ad(tmp_path, compression="gzip")
    tmp_path.replace(output_path)
    logger.info(f"Wrote {output_path} ({output_path.stat().st_size / 1e6:.1f} MB)")
