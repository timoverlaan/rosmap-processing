#!/usr/bin/env python3
"""
Fix categorical data types in AnnData objects.

This module provides functionality to ensure proper categorical types for columns
and remove raw data layers that may cause issues with downstream processing.
"""

import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import anndata as ad

from ..utils.logging import get_logger

logger = get_logger(__name__)


def _is_integer_counts(adata: ad.AnnData, sample_cells: int = 500, sample_genes: int = 200) -> bool:
    """Return True if X appears to contain integer counts stored as floats (>99% whole numbers)."""
    n_cells = min(sample_cells, adata.shape[0])
    n_genes = min(sample_genes, adata.shape[1])
    sample = adata.X[:n_cells, :n_genes]
    if hasattr(sample, 'toarray'):
        sample = sample.toarray()
    else:
        sample = np.array(sample)
    nonzero = sample[sample > 0].flatten()
    if len(nonzero) == 0:
        return False
    return (np.sum(nonzero == np.floor(nonzero)) / len(nonzero)) > 0.99


def fix_categories(
    adata: ad.AnnData,
    columns: Optional[List[str]] = None,
    remove_raw: bool = True,
    cast_counts_to_int: bool = True,
) -> ad.AnnData:
    """
    Fix categorical data types in AnnData observations.
    
    Converts specified columns to categorical type if they aren't already.
    Optionally removes the raw data layer.
    
    Args:
        adata: Input AnnData object
        columns: List of column names to convert to categorical.
                If None, uses default columns.
        remove_raw: Whether to remove the raw data layer
        
    Returns:
        Modified AnnData object
    """
    logger.info("Fixing categorical data types...")
    
    # Default columns to fix for ROSMAP-MIT data
    if columns is None:
        columns = ["cell_type_high_resolution"]
    
    # Fix categorical types
    fixed_columns = []
    missing_columns = []
    
    for col in columns:
        if col in adata.obs.columns:
            if adata.obs[col].dtype != "category":
                adata.obs[col] = adata.obs[col].astype("category")
                fixed_columns.append(col)
                logger.debug(f"Converted '{col}' to categorical type")
            else:
                logger.debug(f"Column '{col}' already categorical")
        else:
            missing_columns.append(col)
            logger.debug(f"Column '{col}' not found in obs")
    
    if fixed_columns:
        logger.info(f"Converted {len(fixed_columns)} columns to categorical: {fixed_columns}")
    
    if missing_columns:
        logger.info(f"Skipped {len(missing_columns)} missing columns: {missing_columns}")
    
    # Cast float-stored integer counts to int32
    if cast_counts_to_int and np.issubdtype(adata.X.dtype, np.floating):
        if _is_integer_counts(adata):
            logger.info(f"X dtype is {adata.X.dtype} but contains integer counts — casting to int32")
            if hasattr(adata.X, 'data'):
                adata.X.data = adata.X.data.astype(np.int32)
            else:
                adata.X = adata.X.astype(np.int32)
        else:
            logger.info(f"X dtype is {adata.X.dtype} with non-integer values — skipping cast")

    # Remove raw data if requested
    if remove_raw and adata.raw is not None:
        logger.info("Removing raw data layer...")
        del adata.raw
        logger.debug("Raw data removed")
    elif remove_raw:
        logger.debug("No raw data layer to remove")
    
    return adata


def fix_categories_in_file(
    input_path: Path,
    output_path: Optional[Path] = None,
    columns: Optional[List[str]] = None,
    remove_raw: bool = True,
    cast_counts_to_int: bool = True,
    compression: str = "gzip"
) -> None:
    """
    Fix categorical types in an h5ad file.
    
    Args:
        input_path: Path to input h5ad file
        output_path: Path to output h5ad file. If None, overwrites input file.
        columns: List of column names to convert to categorical
        remove_raw: Whether to remove the raw data layer
        compression: Compression method for output file
    """
    logger.info(f"Loading data from {input_path}...")
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    # Load data
    adata = ad.read_h5ad(input_path)
    logger.info(f"Loaded AnnData: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Log current state
    logger.debug(f"obs columns: {list(adata.obs.columns)}")
    logger.debug(f"var columns: {list(adata.var.columns)}")
    if adata.raw is not None:
        logger.debug(f"raw data shape: {adata.raw.shape}")
    
    # Fix categories
    adata = fix_categories(adata, columns=columns, remove_raw=remove_raw, cast_counts_to_int=cast_counts_to_int)
    
    # Determine output path
    if output_path is None:
        output_path = input_path
        logger.info("Overwriting input file...")
    else:
        logger.info(f"Saving to new file: {output_path}")
    
    # Save modified data
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression=compression)
    
    file_size_mb = output_path.stat().st_size / 1e6
    logger.info(f"Saved modified file ({file_size_mb:.1f} MB)")


def main():
    """Command-line interface for fixing categories."""
    import argparse
    from ..utils.logging import setup_logging
    
    parser = argparse.ArgumentParser(
        description="Fix categorical data types in AnnData objects. "
                    "Ensures proper categorical types for specified columns and "
                    "optionally removes raw data layer."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to input h5ad file"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Path to output h5ad file (default: overwrite input)"
    )
    parser.add_argument(
        "--columns",
        type=str,
        nargs="+",
        help="Column names to convert to categorical (default: cell_type_high_resolution)"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep raw data layer (default: remove it)"
    )
    parser.add_argument(
        "--no-cast-int",
        action="store_true",
        help="Skip casting float-stored integer counts to int32 (default: cast if detected)"
    )
    parser.add_argument(
        "--compression",
        type=str,
        default="gzip",
        help="Compression method for output file (default: gzip)"
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(level=args.log_level)
    
    try:
        input_path = Path(args.path)
        output_path = Path(args.output) if args.output else None
        
        fix_categories_in_file(
            input_path=input_path,
            output_path=output_path,
            columns=args.columns,
            remove_raw=not args.keep_raw,
            cast_counts_to_int=not args.no_cast_int,
            compression=args.compression
        )
        
        logger.info("✓ Category fixes complete!")
        
    except Exception as e:
        logger.error(f"Failed to fix categories: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
