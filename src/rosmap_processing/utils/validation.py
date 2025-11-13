#!/usr/bin/env python3
"""
Validate and inspect AnnData h5ad files.

This module provides functionality to check the structure and contents of h5ad files,
including metadata, data matrices, and basic statistics.
"""

import sys
from pathlib import Path
import anndata as ad
import numpy as np
from typing import Union

from ..utils.logging import get_logger

logger = get_logger(__name__)


def check_h5ad_structure(adata: ad.AnnData, verbose: bool = True) -> dict:
    """
    Check and report on the structure of an AnnData object.
    
    Args:
        adata: AnnData object to inspect
        verbose: Whether to log detailed information
        
    Returns:
        Dictionary with structure information
    """
    info = {
        "shape": adata.shape,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "has_raw": adata.raw is not None,
        "layers": list(adata.layers.keys()) if adata.layers else [],
        "obsm_keys": list(adata.obsm.keys()) if adata.obsm else [],
        "varm_keys": list(adata.varm.keys()) if adata.varm else [],
        "uns_keys": list(adata.uns.keys()) if adata.uns else [],
    }
    
    if verbose:
        logger.info("=" * 60)
        logger.info(f"AnnData Structure")
        logger.info("=" * 60)
        logger.info(f"Shape: {info['n_obs']} cells × {info['n_vars']} genes")
        logger.info(f"obs columns ({len(info['obs_columns'])}): {', '.join(info['obs_columns'][:5])}")
        if len(info['obs_columns']) > 5:
            logger.info(f"  ... and {len(info['obs_columns']) - 5} more")
        logger.info(f"var columns ({len(info['var_columns'])}): {', '.join(info['var_columns'][:5])}")
        if len(info['var_columns']) > 5:
            logger.info(f"  ... and {len(info['var_columns']) - 5} more")
        logger.info(f"Has raw data: {info['has_raw']}")
        if info['layers']:
            logger.info(f"Layers: {', '.join(info['layers'])}")
        if info['obsm_keys']:
            logger.info(f"obsm keys: {', '.join(info['obsm_keys'])}")
        if info['varm_keys']:
            logger.info(f"varm keys: {', '.join(info['varm_keys'])}")
        if info['uns_keys']:
            logger.info(f"uns keys: {', '.join(info['uns_keys'])}")
        logger.info("")
    
    return info


def check_obs_metadata(
    adata: ad.AnnData,
    show_counts: bool = True,
    max_categories: int = 20
) -> dict:
    """
    Check and report on observation metadata.
    
    Args:
        adata: AnnData object to inspect
        show_counts: Whether to show value counts for categorical columns
        max_categories: Maximum number of categories to show
        
    Returns:
        Dictionary with metadata statistics
    """
    logger.info("=" * 60)
    logger.info("Observation Metadata")
    logger.info("=" * 60)
    
    metadata = {}
    
    for col in adata.obs.columns:
        col_data = adata.obs[col]
        
        # Collect statistics
        col_info = {
            "dtype": str(col_data.dtype),
            "n_unique": col_data.nunique(),
            "n_missing": col_data.isna().sum()
        }
        
        logger.info(f"\n{col}:")
        logger.info(f"  Type: {col_info['dtype']}")
        logger.info(f"  Unique values: {col_info['n_unique']}")
        
        if col_info['n_missing'] > 0:
            logger.info(f"  Missing: {col_info['n_missing']} ({col_info['n_missing']/len(col_data)*100:.1f}%)")
        
        # Show value counts for categorical/low-cardinality columns
        if show_counts and col_info['n_unique'] <= max_categories:
            counts = col_data.value_counts()
            col_info["value_counts"] = counts.to_dict()
            
            logger.info("  Value counts:")
            for value, count in counts.items():
                logger.info(f"    {value}: {count}")
        elif col_info['n_unique'] > max_categories:
            logger.info(f"  (Too many unique values to display)")
        
        metadata[col] = col_info
    
    logger.info("")
    return metadata


def check_data_matrix(
    adata: ad.AnnData,
    preview_size: int = 20,
    show_raw: bool = True
) -> dict:
    """
    Check and preview data matrices.
    
    Args:
        adata: AnnData object to inspect
        preview_size: Number of rows/columns to preview
        show_raw: Whether to check raw data if available
        
    Returns:
        Dictionary with matrix information
    """
    logger.info("=" * 60)
    logger.info("Data Matrices")
    logger.info("=" * 60)
    
    matrix_info = {}
    
    # Check main matrix
    logger.info(f"\nMain matrix (X):")
    logger.info(f"  Shape: {adata.X.shape}")
    logger.info(f"  Type: {type(adata.X)}")
    logger.info(f"  Dtype: {adata.X.dtype}")
    
    # Get preview data
    n_preview = min(preview_size, adata.X.shape[0], adata.X.shape[1])
    
    if hasattr(adata.X, 'toarray'):
        # Sparse matrix
        preview_data = adata.X[:n_preview, :n_preview].toarray()
        logger.info(f"  Format: Sparse ({adata.X.format if hasattr(adata.X, 'format') else 'unknown'})")
        
        # Calculate sparsity
        if adata.X.shape[0] * adata.X.shape[1] > 0:
            sparsity = 1.0 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
            logger.info(f"  Sparsity: {sparsity*100:.1f}%")
            matrix_info["sparsity"] = sparsity
    else:
        # Dense matrix
        preview_data = adata.X[:n_preview, :n_preview]
        logger.info(f"  Format: Dense")
    
    # Show statistics
    if hasattr(adata.X, 'data'):
        # Sparse matrix
        data_vals = adata.X.data
    else:
        data_vals = adata.X.flatten()
    
    if len(data_vals) > 0:
        logger.info(f"  Min: {np.min(data_vals):.4f}")
        logger.info(f"  Max: {np.max(data_vals):.4f}")
        logger.info(f"  Mean: {np.mean(data_vals):.4f}")
        logger.info(f"  Std: {np.std(data_vals):.4f}")
        
        matrix_info["X"] = {
            "shape": adata.X.shape,
            "min": float(np.min(data_vals)),
            "max": float(np.max(data_vals)),
            "mean": float(np.mean(data_vals)),
            "std": float(np.std(data_vals))
        }
    
    logger.info(f"\n  Preview ({n_preview}×{n_preview}):")
    logger.info(f"{preview_data}")
    
    # Check raw data if available
    if show_raw and adata.raw is not None:
        logger.info(f"\nRaw data (raw.X):")
        logger.info(f"  Shape: {adata.raw.X.shape}")
        logger.info(f"  Type: {type(adata.raw.X)}")
        logger.info(f"  Dtype: {adata.raw.X.dtype}")
        
        raw_preview = adata.raw.X[:n_preview, :n_preview]
        if hasattr(raw_preview, 'toarray'):
            raw_preview = raw_preview.toarray()
        
        logger.info(f"\n  Preview ({n_preview}×{n_preview}):")
        logger.info(f"{raw_preview}")
        
        matrix_info["raw"] = {"shape": adata.raw.X.shape}
    
    logger.info("")
    return matrix_info


def check_h5ad_file(
    file_path: Path,
    show_metadata: bool = True,
    show_matrix: bool = True,
    preview_size: int = 20
) -> dict:
    """
    Check and validate an h5ad file.
    
    Args:
        file_path: Path to h5ad file
        show_metadata: Whether to show metadata details
        show_matrix: Whether to show matrix details
        preview_size: Number of rows/columns to preview
        
    Returns:
        Dictionary with all check results
    """
    logger.info(f"Checking h5ad file: {file_path}")
    logger.info("")
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Load data
    logger.info("Loading AnnData object...")
    adata = ad.read_h5ad(file_path)
    
    results = {}
    
    # Check structure
    results["structure"] = check_h5ad_structure(adata, verbose=True)
    
    # Check metadata
    if show_metadata:
        results["metadata"] = check_obs_metadata(adata, show_counts=True)
    
    # Check data matrix
    if show_matrix:
        results["matrix"] = check_data_matrix(adata, preview_size=preview_size)
    
    return results


def inspect_h5ad(path: Union[str, Path]):
    """Compatibility wrapper used by the CLI.

    Accepts a string or Path and runs the full h5ad checks, logging errors
    and raising exceptions on failure so the CLI can handle the exit.
    """
    file_path = Path(path)
    return check_h5ad_file(file_path)


def main():
    """Command-line interface for checking h5ad files."""
    import argparse
    from ..utils.logging import setup_logging
    
    parser = argparse.ArgumentParser(
        description="Check and validate an h5ad file. "
                    "Displays structure, metadata, and data matrix information."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to h5ad file to check"
    )
    parser.add_argument(
        "--no-metadata",
        action="store_true",
        help="Skip detailed metadata inspection"
    )
    parser.add_argument(
        "--no-matrix",
        action="store_true",
        help="Skip data matrix inspection"
    )
    parser.add_argument(
        "--preview-size",
        type=int,
        default=20,
        help="Number of rows/columns to preview (default: 20)"
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
        file_path = Path(args.path)
        
        check_h5ad_file(
            file_path=file_path,
            show_metadata=not args.no_metadata,
            show_matrix=not args.no_matrix,
            preview_size=args.preview_size
        )
        
        logger.info("✓ File check complete!")
        
    except Exception as e:
        logger.error(f"Failed to check file: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
