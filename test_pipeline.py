#!/usr/bin/env python3
"""Quick test for full MIT pipeline."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from rosmap_processing.utils.logging import setup_logging
import anndata as ad
import pandas as pd

logger = setup_logging(level="INFO")

MIT_H5AD = Path("data/raw/ROSMAP_MIT/immune_cells.h5ad")
MIT_METADATA = Path("data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv")
MIT_OUTPUT = Path("test_output/immune_cells_processed.h5ad")

logger.info("Testing full MIT pipeline with subset...")

# Load and take a small subset (1000 cells) for quick testing
logger.info(f"Loading subset from: {MIT_H5AD}")
adata = ad.read_h5ad(MIT_H5AD)
logger.info(f"  Original shape: {adata.shape}")

# Take subset
n_cells = min(1000, adata.shape[0])
adata = adata[:n_cells, :].copy()
logger.info(f"  Subset shape: {adata.shape}")

# Step 1: Fix categories
logger.info("\n  Step 1: Fixing categories...")
if "cell_type_high_resolution" in adata.obs.columns:
    adata.obs["cell_type_high_resolution"] = adata.obs["cell_type_high_resolution"].astype("category")
if adata.raw is not None:
    del adata.raw

# Step 2: Add metadata
logger.info("\n  Step 2: Adding metadata...")
metadata = pd.read_csv(MIT_METADATA)

# MIT data: h5ad has 'projid' (int), metadata has 'individualID' (string)
logger.info(f"    Metadata columns: {list(metadata.columns[:5])}")
logger.info(f"    AnnData obs columns: {list(adata.obs.columns[:5])}")

adata.obs["projid"] = adata.obs["projid"].astype(str)
metadata["individualID"] = metadata["individualID"].astype(str)

# Merge
adata.obs = adata.obs.merge(
    metadata, left_on="projid", right_on="individualID", how="left"
)
adata.obs.index = adata.obs_names
logger.info(f"    Columns after metadata: {len(adata.obs.columns)}")

# Step 3: Convert columns
logger.info("\n  Step 3: Converting columns...")
if "projid" in adata.obs.columns:
    adata.obs.rename(columns={"projid": "Donor ID"}, inplace=True)
if "cell_type_high_resolution" in adata.obs.columns:
    adata.obs.rename(columns={"cell_type_high_resolution": "Subtype"}, inplace=True)
adata.obs["Class"] = "Glia"
adata.obs["Celltype"] = "Immune cells"

# Convert all object columns to string for h5ad compatibility
logger.info("\n  Step 4: Cleaning data types...")
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object':
        adata.obs[col] = adata.obs[col].astype(str)
        logger.info(f"    Converted {col} to string")

# Save
logger.info(f"\n  Saving to: {MIT_OUTPUT}")
MIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

try:
    adata.write_h5ad(MIT_OUTPUT, compression="gzip")
    file_size = MIT_OUTPUT.stat().st_size / 1e6
    logger.info(f"  Output size: {file_size:.1f} MB")
    logger.info("SUCCESS! Full MIT pipeline works!")
except Exception as e:
    logger.error(f"FAILED: {e}", exc_info=True)
    sys.exit(1)
