"""Add clinical metadata to AnnData objects."""

import anndata as ad
import pandas as pd
from pathlib import Path
from typing import Union, Optional
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rosmap_processing.utils.logging import get_logger

logger = get_logger(__name__)


def add_metadata(
    adata: ad.AnnData,
    metadata: pd.DataFrame,
    is_mit: bool = False,
    validate: bool = True,
) -> ad.AnnData:
    """
    Add clinical metadata to AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to add metadata to
    metadata : DataFrame
        Metadata DataFrame with clinical information
    is_mit : bool, default False
        Whether this is MIT data (uses 'individualID' instead of 'projid')
    validate : bool, default True
        Whether to validate the merge operation
        
    Returns
    -------
    AnnData
        AnnData object with added metadata
        
    Raises
    ------
    ValueError
        If join column is not found or validation fails
    """
    logger.info("Adding clinical metadata to AnnData object")
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Metadata shape: {metadata.shape}")
    logger.info(f"Is MIT data: {is_mit}")
    
    # Determine join column based on what's available
    # Both MIT and non-MIT data should use 'projid' to join with rosmap_clinical.csv
    # The rosmap_clinical.csv metadata file has both 'projid' and 'individualID' columns
    
    adata_col = None
    metadata_col = None
    
    # First, try to use 'projid' (preferred for both MIT and non-MIT)
    if "projid" in adata.obs.columns and "projid" in metadata.columns:
        adata_col = "projid"
        metadata_col = "projid"
        logger.info(f"Using 'projid' as join column")
    # Fallback to individualID if projid not available
    elif "individualID" in adata.obs.columns and "individualID" in metadata.columns:
        adata_col = "individualID"
        metadata_col = "individualID"
        logger.info(f"Using 'individualID' as join column")
    else:
        available_adata = list(adata.obs.columns[:10])
        available_metadata = list(metadata.columns[:10])
        raise ValueError(
            f"Cannot find matching join column. "
            f"h5ad columns: {available_adata}... "
            f"metadata columns: {available_metadata}..."
        )
    
    # Validate both columns exist
    if adata_col not in adata.obs.columns:
        raise ValueError(
            f"Column '{adata_col}' not found in h5ad. "
            f"Available: {list(adata.obs.columns[:10])}..."
        )
    
    if metadata_col not in metadata.columns:
        raise ValueError(
            f"Column '{metadata_col}' not found in metadata. "
            f"Available: {list(metadata.columns[:10])}..."
        )
    
    # Store original state
    original_index = adata.obs.index.copy()
    original_obs_cols = set(adata.obs.columns)
    n_obs_before = len(adata.obs)
    
    # Convert columns to compatible types if needed
    # Often projid is int in h5ad but string in metadata
    if adata.obs[adata_col].dtype != metadata[metadata_col].dtype:
        logger.info(
            f"Converting types for merge: "
            f"h5ad[{adata_col}]={adata.obs[adata_col].dtype} → "
            f"metadata[{metadata_col}]={metadata[metadata_col].dtype}"
        )
        # Convert both to string for safer merging
        adata.obs[adata_col] = adata.obs[adata_col].astype(str)
        metadata = metadata.copy()  # Don't modify original
        metadata[metadata_col] = metadata[metadata_col].astype(str)
    
    # If columns have different names, temporarily rename metadata column
    temp_metadata = metadata.copy()
    if adata_col != metadata_col:
        temp_metadata = temp_metadata.rename(columns={metadata_col: adata_col})
        merge_col = adata_col
    else:
        merge_col = adata_col
    
    # Perform merge
    logger.info(f"Merging on column: {merge_col}")
    merged_obs = adata.obs.merge(
        temp_metadata,
        on=merge_col,
        how='left',
        suffixes=('', '_metadata')
    )
    
    # Restore original index (merge resets it)
    merged_obs.index = original_index
    
    # Validation
    if validate:
        n_obs_after = len(merged_obs)
        if n_obs_after != n_obs_before:
            raise ValueError(
                f"Merge changed number of observations: "
                f"{n_obs_before} → {n_obs_after}"
            )
        
        # Check for new NaN values (indicating failed joins)
        new_cols = set(temp_metadata.columns) - set(adata.obs.columns)
        for col in new_cols:
            # Skip the join column
            if col == merge_col:
                continue
            if col in merged_obs.columns:
                n_missing = merged_obs[col].isna().sum()
                if n_missing > 0:
                    logger.warning(
                        f"Column '{col}' has {n_missing} missing values after merge"
                    )
    
    # Update AnnData
    adata.obs = merged_obs
    
    # Clean up data types for h5ad compatibility
    # Convert columns with mixed types to strings
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            # Check if it can be numeric
            try:
                adata.obs[col] = pd.to_numeric(adata.obs[col])
            except (ValueError, TypeError):
                # If not numeric, convert to string
                adata.obs[col] = adata.obs[col].astype(str)
                logger.debug(f"Converted column '{col}' to string for h5ad compatibility")
    
    # Clean up raw data if it exists
    # Remove any columns with reserved names like '_index'
    if adata.raw is not None:
        reserved_names = ['_index', '_categories']
        for col in reserved_names:
            if col in adata.raw.var.columns:
                logger.debug(f"Removing reserved column '{col}' from raw.var")
                adata.raw.var.drop(columns=[col], inplace=True)
    
    # Log summary
    old_cols_set = original_obs_cols
    new_cols = set(merged_obs.columns) - old_cols_set
    logger.info(f"Added {len(new_cols)} new metadata columns")
    if new_cols:
        logger.info(f"New columns: {sorted(new_cols)[:10]}...")  # Show first 10
    
    return adata


def add_metadata_to_file(
    h5ad_path: Union[str, Path],
    metadata_path: Union[str, Path],
    is_mit: bool = False,
    output_path: Optional[Union[str, Path]] = None,
) -> None:
    """
    Add metadata to h5ad file.
    
    Parameters
    ----------
    h5ad_path : str or Path
        Path to h5ad file (will be modified in place if output_path is None)
    metadata_path : str or Path
        Path to metadata CSV file
    is_mit : bool, default False
        Whether this is MIT data
    output_path : str or Path, optional
        Path to save output. If None, modifies file in place
        
    Raises
    ------
    FileNotFoundError
        If input files don't exist
    """
    h5ad_path = Path(h5ad_path)
    metadata_path = Path(metadata_path)
    
    logger.info("="*60)
    logger.info("Adding metadata to h5ad file")
    logger.info("="*60)
    logger.info(f"H5ad file: {h5ad_path}")
    logger.info(f"Metadata file: {metadata_path}")
    
    # Validate inputs
    if not h5ad_path.exists():
        raise FileNotFoundError(f"H5ad file not found: {h5ad_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # Load data
    logger.info("Loading AnnData object...")
    adata = ad.read_h5ad(h5ad_path)
    
    logger.info("Loading metadata CSV...")
    metadata = pd.read_csv(metadata_path)
    
    # Add metadata
    adata = add_metadata(adata, metadata, is_mit=is_mit)
    
    # Save
    if output_path is None:
        output_path = h5ad_path
        logger.info(f"Saving modified AnnData (in place): {output_path}")
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving modified AnnData to: {output_path}")
    
    adata.write_h5ad(output_path, compression='gzip')
    
    logger.info("="*60)
    logger.info("Metadata addition complete!")
    logger.info("="*60)


if __name__ == "__main__":
    import argparse
    from typing import Optional
    
    parser = argparse.ArgumentParser(
        description="Add clinical metadata to an AnnData h5ad file."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to the input h5ad file (will be modified in place unless --output is specified)."
    )
    parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        help="Path to the metadata CSV file."
    )
    parser.add_argument(
        "--MIT",
        action="store_true",
        help="If set, the input data is from MIT (uses 'individualID' for matching)."
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to save output file. If not provided, modifies input file in place."
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
    from rosmap_processing.utils.logging import setup_logging
    setup_logging(level=args.log_level)
    
    try:
        add_metadata_to_file(
            h5ad_path=args.path,
            metadata_path=args.metadata,
            is_mit=args.MIT,
            output_path=args.output,
        )
    except Exception as e:
        logger.error(f"Failed to add metadata: {e}", exc_info=True)
        sys.exit(1)
