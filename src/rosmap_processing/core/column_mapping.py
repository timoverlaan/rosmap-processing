#!/usr/bin/env python3
"""
Column name and type conversion between SeaAD and ROSMAP data formats.

This module provides functionality to convert column names and types between different
single-cell RNA-seq datasets (ROSMAP, ROSMAP-MIT, SeaAD) to enable cross-dataset analysis.
Also derives Wang et al. AD classification labels.
"""

import sys
from pathlib import Path
from typing import Optional

import anndata as ad

from ..utils.logging import get_logger
from ..utils.constants import (
    BRAAK_MAPPING,
    CERAD_MAPPING,
    THAL_MAPPING,
    ROSMAP_MIT_CELL_CLASSES
)

logger = get_logger(__name__)

# Column name constants
CLASS_NAME = "Class"
SUBCLASS_NAME = "Celltype"
SUPERTYPE_NAME = "Subtype"


def convert_rosmap_columns(adata: ad.AnnData) -> ad.AnnData:
    """
    Convert ROSMAP column names to standardized format.
    
    Adds SeaAD-style columns to the AnnData object.
    
    Args:
        adata: AnnData object with ROSMAP data
        
    Returns:
        Modified AnnData object
    """
    logger.info("Converting ROSMAP columns to standard format...")
    
    # Rename columns to standard format
    rename_map = {
        "individualID": "Donor ID",
        "class": CLASS_NAME,
        "subset": SUBCLASS_NAME,
        "state": SUPERTYPE_NAME,
    }
    
    # Only rename columns that exist
    existing_renames = {k: v for k, v in rename_map.items() if k in adata.obs.columns}
    
    if existing_renames:
        adata.obs.rename(columns=existing_renames, inplace=True)
        logger.info(f"Renamed columns: {list(existing_renames.keys())}")
    else:
        logger.warning("No ROSMAP columns found to rename")
    
    logger.info(f"ROSMAP conversion complete. Shape: {adata.shape}")
    return adata


def convert_rosmap_mit_columns(
    adata: ad.AnnData,
    cellclass: str,
    subclass: str
) -> ad.AnnData:
    """
    Convert ROSMAP-MIT column names to standardized format.
    
    This is a special case for MIT data which has different column structure.
    
    Args:
        adata: AnnData object with ROSMAP-MIT data
        cellclass: Cell class label (e.g., 'Neuron', 'Glia')
        subclass: Cell subclass label (e.g., 'Astrocytes', 'Excitatory')
        
    Returns:
        Modified AnnData object
    """
    logger.info(f"Converting ROSMAP-MIT columns (class={cellclass}, subclass={subclass})...")
    
    # Validate cell class
    if cellclass not in ROSMAP_MIT_CELL_CLASSES:
        logger.warning(
            f"Cell class '{cellclass}' not in known classes: {list(ROSMAP_MIT_CELL_CLASSES.keys())}"
        )
    
    # Rename columns
    rename_map = {
        "projid": "Donor ID",
        "cell_type_high_resolution": SUPERTYPE_NAME,
    }
    
    existing_renames = {k: v for k, v in rename_map.items() if k in adata.obs.columns}
    
    if existing_renames:
        adata.obs.rename(columns=existing_renames, inplace=True)
        logger.info(f"Renamed columns: {list(existing_renames.keys())}")
    
    # Add class and subclass columns
    adata.obs[CLASS_NAME] = cellclass
    adata.obs[SUBCLASS_NAME] = subclass
    logger.info(f"Added {CLASS_NAME}='{cellclass}' and {SUBCLASS_NAME}='{subclass}'")
    
    logger.info(f"ROSMAP-MIT conversion complete. Shape: {adata.shape}")
    return adata


def convert_seaad_columns(adata: ad.AnnData) -> ad.AnnData:
    """
    Convert SeaAD column names to ROSMAP-compatible format.
    
    Adds ROSMAP-style numeric columns and derives Wang AD labels.
    
    Args:
        adata: AnnData object with SeaAD data
        
    Returns:
        Modified AnnData object
    """
    logger.info("Converting SeaAD columns to ROSMAP format...")
    
    # Rename columns to standard format
    rename_map = {
        "Class": CLASS_NAME,
        "Subclass": SUBCLASS_NAME,
        "Supertype": SUPERTYPE_NAME,
    }
    
    existing_renames = {k: v for k, v in rename_map.items() if k in adata.obs.columns}
    
    if existing_renames:
        adata.obs.rename(columns=existing_renames, inplace=True)
        logger.info(f"Renamed columns: {list(existing_renames.keys())}")
    
    # Convert categorical columns to numeric using mappings from constants
    numeric_conversions = []
    
    # Braak score
    if "Braak" in adata.obs.columns:
        adata.obs["braaksc"] = adata.obs["Braak"].map(BRAAK_MAPPING)
        numeric_conversions.append("braaksc")
        logger.debug(f"Converted Braak: {adata.obs['Braak'].value_counts().to_dict()}")
    
    # CERAD score (inverted: lower numeric = more pathology in ROSMAP)
    if "CERAD" in adata.obs.columns:
        adata.obs["ceradsc"] = adata.obs["CERAD"].map(CERAD_MAPPING)
        numeric_conversions.append("ceradsc")
        logger.debug(f"Converted CERAD: {adata.obs['CERAD'].value_counts().to_dict()}")
    
    # Thal score
    if "Thal" in adata.obs.columns:
        adata.obs["thalsc"] = adata.obs["Thal"].map(THAL_MAPPING)
        numeric_conversions.append("thalsc")
        logger.debug(f"Converted Thal: {adata.obs['Thal'].value_counts().to_dict()}")
    
    if numeric_conversions:
        logger.info(f"Created numeric columns: {numeric_conversions}")
    
    # Derive Wang et al. AD classification labels
    derive_wang_labels(adata)
    
    logger.info(f"SeaAD conversion complete. Shape: {adata.shape}")
    return adata


def derive_wang_labels(adata: ad.AnnData) -> None:
    """
    Derive Wang et al. AD classification labels based on neuropathology.
    
    Classification criteria:
    - AD: Cognitive dementia + Braak ≥4 + CERAD ≤2
    - Healthy: No dementia + Braak ≤3 + CERAD ≥3
    - Intermediate: Everything else (excluded from analysis)
    
    Args:
        adata: AnnData object (modified in place)
    """
    logger.info("Deriving Wang et al. AD classification labels...")
    
    # Check required columns
    required = ["Cognitive Status", "braaksc", "ceradsc"]
    missing = [col for col in required if col not in adata.obs.columns]
    
    if missing:
        logger.warning(f"Cannot derive Wang labels - missing columns: {missing}")
        return
    
    # Initialize all as "Intermediate"
    adata.obs["Wang"] = "Intermediate"
    adata.obs["Wang_intermediate"] = True
    
    # Identify AD cases
    ad_mask = (
        (adata.obs["Cognitive Status"] == "Dementia") &
        (adata.obs["braaksc"] >= 4) &
        (adata.obs["ceradsc"] <= 2)
    )
    adata.obs.loc[ad_mask, "Wang"] = "AD"
    adata.obs.loc[ad_mask, "Wang_intermediate"] = False
    
    # Identify Healthy cases
    healthy_mask = (
        (adata.obs["Cognitive Status"] != "Dementia") &
        (adata.obs["braaksc"] <= 3) &
        (adata.obs["ceradsc"] >= 3)
    )
    adata.obs.loc[healthy_mask, "Wang"] = "Healthy"
    adata.obs.loc[healthy_mask, "Wang_intermediate"] = False
    
    # Log statistics
    wang_counts = adata.obs["Wang"].value_counts()
    logger.info("Wang label distribution:")
    for label, count in wang_counts.items():
        logger.info(f"  {label}: {count} cells")
    
    # Count unique donors in non-intermediate samples
    if "Donor ID" in adata.obs.columns:
        n_extreme_donors = adata[~adata.obs["Wang_intermediate"]].obs["Donor ID"].nunique()
        logger.info(f"Donors with extreme phenotypes (AD or Healthy): {n_extreme_donors}")


def convert_columns(
    adata: ad.AnnData,
    data_type: str,
    cellclass: Optional[str] = None,
    subclass: Optional[str] = None
) -> ad.AnnData:
    """
    Convert column names and types based on dataset type.
    
    Args:
        adata: Input AnnData object
        data_type: Type of data ('ROSMAP', 'ROSMAP_MIT', or 'SeaAD')
        cellclass: Cell class for ROSMAP-MIT data
        subclass: Cell subclass for ROSMAP-MIT data
        
    Returns:
        Modified AnnData object
        
    Raises:
        ValueError: If invalid data_type or missing required parameters
    """
    if data_type == "ROSMAP":
        return convert_rosmap_columns(adata)
    
    elif data_type == "ROSMAP_MIT":
        if not cellclass or not subclass:
            raise ValueError(
                "For ROSMAP_MIT data, both --cellclass and --subclass must be provided"
            )
        return convert_rosmap_mit_columns(adata, cellclass, subclass)
    
    elif data_type == "SeaAD":
        return convert_seaad_columns(adata)
    
    else:
        raise ValueError(
            f"Invalid data type: {data_type}. Must be 'ROSMAP', 'ROSMAP_MIT', or 'SeaAD'"
        )


def main():
    """Command-line interface for column conversion."""
    import argparse
    from ..utils.logging import setup_logging
    
    parser = argparse.ArgumentParser(
        description="Convert column names and types between SeaAD and ROSMAP data formats. "
                    "Also derives Wang et al. AD classification labels."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to input h5ad file"
    )
    parser.add_argument(
        "--type",
        type=str,
        required=True,
        choices=["ROSMAP", "ROSMAP_MIT", "SeaAD"],
        help="Type of input data"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Path to output h5ad file (if not using --inplace)"
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Modify input file in place"
    )
    parser.add_argument(
        "--cellclass",
        type=str,
        help="Cell class for ROSMAP_MIT data (e.g., 'Neuron', 'Glia')"
    )
    parser.add_argument(
        "--subclass",
        type=str,
        help="Cell subclass for ROSMAP_MIT data (e.g., 'Astrocytes', 'Excitatory')"
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
        # Validate output arguments
        if not args.inplace and not args.output:
            raise ValueError("Either --inplace or --output must be provided")
        
        # Load data
        input_path = Path(args.path)
        logger.info(f"Loading data from {input_path}...")
        adata = ad.read_h5ad(input_path)
        logger.info(f"Loaded AnnData: {adata.shape[0]} cells × {adata.shape[1]} genes")
        
        # Convert columns
        adata = convert_columns(
            adata,
            data_type=args.type,
            cellclass=args.cellclass,
            subclass=args.subclass
        )
        
        # Determine output path
        output_path = input_path if args.inplace else Path(args.output)
        
        # Save result
        logger.info(f"\nSaving converted data to {output_path}...")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_path, compression="gzip")
        
        logger.info("✓ Conversion complete!")
        
    except Exception as e:
        logger.error(f"Failed to convert columns: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
