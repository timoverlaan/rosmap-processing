#!/usr/bin/env python3
"""Quick test for metadata addition."""

import sys
from pathlib import Path
import shutil

sys.path.insert(0, str(Path(__file__).parent / "src"))

from rosmap_processing.utils.logging import setup_logging
from rosmap_processing.data.metadata import add_metadata_to_file

logger = setup_logging(level="DEBUG")

MIT_H5AD = Path("data/raw/ROSMAP_MIT/immune_cells.h5ad")
MIT_METADATA = Path("data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv")
TEST_DIR = Path("test_output")
TEST_DIR.mkdir(exist_ok=True)

logger.info("Testing metadata addition for MIT data...")

# Make a copy
test_file = TEST_DIR / "immune_cells_test_metadata.h5ad"
logger.info(f"Creating test copy: {test_file}")
shutil.copy(MIT_H5AD, test_file)

# Add metadata
logger.info(f"Adding metadata from: {MIT_METADATA}")
try:
    add_metadata_to_file(
        h5ad_path=test_file,
        metadata_path=MIT_METADATA,
        is_mit=True,
        output_path=test_file
    )
    logger.info("SUCCESS!")
except Exception as e:
    logger.error(f"FAILED: {e}", exc_info=True)
    sys.exit(1)

# Verify
import anndata as ad
adata = ad.read_h5ad(test_file)
logger.info(f"Result columns: {len(adata.obs.columns)}")
logger.info(f"Sample columns: {list(adata.obs.columns[:10])}")

test_file.unlink()
logger.info("Test passed!")
