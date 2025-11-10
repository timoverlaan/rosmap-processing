#!/usr/bin/env python3
"""
Combine multiple h5ad files into a single AnnData object.

This module provides functionality to concatenate multiple AnnData files with proper
validation and error handling.
"""

import sys
from pathlib import Path
from typing import List, Optional

import anndata as ad
from tqdm import tqdm

from ..utils.logging import get_logger

logger = get_logger(__name__)


def validate_h5ad_files(paths: List[Path]) -> None:
    """
    Validate that all input files exist and are h5ad files.
    
    Args:
        paths: List of file paths to validate
        
    Raises:
        FileNotFoundError: If any file doesn't exist
        ValueError: If any file doesn't have .h5ad extension
    """
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        if path.suffix != ".h5ad":
            raise ValueError(f"File must have .h5ad extension: {path}")


def combine_h5ad_files(
    paths: List[Path],
    join: str = "outer",
    merge: str = "unique",
    label_key: str = "source_file",
    pairwise: bool = True
) -> ad.AnnData:
    """
    Combine multiple h5ad files into a single AnnData object.
    
    Args:
        paths: List of paths to h5ad files to combine
        join: Join strategy for concatenation ('inner' or 'outer')
        merge: Strategy for uns merge ('same', 'unique', 'first', or 'only')
        label_key: Key name for storing source file labels in obs
        pairwise: Whether to combine pairwise obs/var annotations
        
    Returns:
        Combined AnnData object
        
    Raises:
        FileNotFoundError: If any input file doesn't exist
        ValueError: If input files are invalid
    """
    logger.info(f"Combining {len(paths)} h5ad files...")
    
    # Validate all files first
    validate_h5ad_files(paths)
    
    # Load all files
    logger.info("Loading h5ad files...")
    adatas = []
    for path in tqdm(paths, desc="Loading files"):
        try:
            adata = ad.read_h5ad(path)
            logger.debug(f"Loaded {path}: shape={adata.shape}")
            adatas.append(adata)
        except Exception as e:
            logger.error(f"Failed to load {path}: {e}")
            raise
    
    if not adatas:
        raise ValueError("No valid AnnData objects loaded")
    
    # Log shapes before concatenation
    for i, (adata, path) in enumerate(zip(adatas, paths)):
        logger.info(f"  File {i+1}: {path.name} - {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Concatenate all AnnData objects
    logger.info(f"Concatenating AnnData objects (join={join}, merge={merge})...")
    try:
        combined_adata = ad.concat(
            adatas,
            join=join,
            label=label_key,
            keys=[str(p) for p in paths],
            pairwise=pairwise,
            merge=merge
        )
    except Exception as e:
        logger.error(f"Failed to concatenate AnnData objects: {e}")
        raise
    
    logger.info(f"Combined AnnData shape: {combined_adata.shape}")
    logger.info(f"  Total cells: {combined_adata.shape[0]}")
    logger.info(f"  Total genes: {combined_adata.shape[1]}")
    
    # Log basic statistics
    if label_key in combined_adata.obs.columns:
        counts = combined_adata.obs[label_key].value_counts()
        logger.info(f"\nCells per source file:")
        for source, count in counts.items():
            logger.info(f"  {source}: {count} cells")
    
    return combined_adata


def combine_and_save(
    input_paths: List[Path],
    output_path: Path,
    join: str = "outer",
    merge: str = "unique",
    compression: str = "gzip"
) -> None:
    """
    Combine multiple h5ad files and save to output file.
    
    Args:
        input_paths: List of paths to input h5ad files
        output_path: Path to output h5ad file
        join: Join strategy for concatenation
        merge: Strategy for uns merge
        compression: Compression method for output file
    """
    # Combine files
    combined_adata = combine_h5ad_files(input_paths, join=join, merge=merge)
    
    # Save combined file
    logger.info(f"\nSaving combined AnnData to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        combined_adata.write_h5ad(output_path, compression=compression)
        logger.info(f"Successfully saved combined file ({output_path.stat().st_size / 1e6:.1f} MB)")
    except Exception as e:
        logger.error(f"Failed to save combined file: {e}")
        raise


def main():
    """Command-line interface for combining h5ad files."""
    import argparse
    from ..utils.logging import setup_logging
    
    parser = argparse.ArgumentParser(
        description="Combine multiple h5ad files into a single AnnData file."
    )
    parser.add_argument(
        "paths",
        type=str,
        nargs="+",
        help="Paths to input h5ad files"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Path to output h5ad file"
    )
    parser.add_argument(
        "--join",
        type=str,
        default="outer",
        choices=["inner", "outer"],
        help="Join strategy for concatenation (default: outer)"
    )
    parser.add_argument(
        "--merge",
        type=str,
        default="unique",
        choices=["same", "unique", "first", "only"],
        help="Strategy for uns merge (default: unique)"
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
        # Convert paths to Path objects
        input_paths = [Path(p) for p in args.paths]
        output_path = Path(args.output)
        
        # Combine and save
        combine_and_save(
            input_paths=input_paths,
            output_path=output_path,
            join=args.join,
            merge=args.merge,
            compression=args.compression
        )
        
        logger.info("\n✓ Combination complete!")
        
    except Exception as e:
        logger.error(f"Failed to combine files: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
