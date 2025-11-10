#!/usr/bin/env python3
"""
Local testing script for refactored modules.

Tests the refactored code on small files before running on cluster.
Tests both ROSMAP-MIT and ROSMAP data.
"""

import sys
from pathlib import Path
import shutil

# Add src to path so we can import the package
sys.path.insert(0, str(Path(__file__).parent / "src"))

from rosmap_processing.utils.logging import setup_logging
from rosmap_processing.utils.validation import check_h5ad_file
from rosmap_processing.data.metadata import add_metadata_to_file
from rosmap_processing.core.column_mapping import convert_columns
from rosmap_processing.data.category_fix import fix_categories_in_file
import anndata as ad

# Setup logging
logger = setup_logging(level="INFO")

# Define test files
TEST_DIR = Path("test_output")
TEST_DIR.mkdir(exist_ok=True)

# ROSMAP-MIT test files (smallest: immune_cells.h5ad ~606 MB)
MIT_H5AD = Path("data/raw/ROSMAP_MIT/immune_cells.h5ad")
MIT_METADATA = Path("data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv")
MIT_OUTPUT = TEST_DIR / "immune_cells_processed.h5ad"

# ROSMAP test files (astrocytes.h5ad ~8GB - we'll use a subset)
ROSMAP_H5AD = Path("data/raw/ROSMAP/astrocytes.h5ad")
ROSMAP_METADATA = Path("data/raw/ROSMAP/rosmap_clinical.csv")
ROSMAP_OUTPUT = TEST_DIR / "astrocytes_processed.h5ad"


def test_check_h5ad():
    """Test 1: Validation module - check h5ad file structure."""
    logger.info("\n" + "="*70)
    logger.info("TEST 1: Validation module (check_h5ad)")
    logger.info("="*70)
    
    try:
        logger.info(f"\nChecking MIT file: {MIT_H5AD}")
        results = check_h5ad_file(
            MIT_H5AD,
            show_metadata=False,  # Skip detailed metadata to keep output short
            show_matrix=False,
            preview_size=5
        )
        logger.info(f"✓ Successfully checked MIT file")
        logger.info(f"  Shape: {results['structure']['n_obs']} cells × {results['structure']['n_vars']} genes")
        return True
    except Exception as e:
        logger.error(f"✗ Failed to check h5ad: {e}", exc_info=True)
        return False


def test_category_fix():
    """Test 2: Category fix module - fix categorical types."""
    logger.info("\n" + "="*70)
    logger.info("TEST 2: Category fix module (fix_categories)")
    logger.info("="*70)
    
    try:
        # Make a copy to test on
        test_file = TEST_DIR / "immune_cells_test.h5ad"
        logger.info(f"Creating test copy: {test_file}")
        shutil.copy(MIT_H5AD, test_file)
        
        # Fix categories
        logger.info("Fixing categories...")
        fix_categories_in_file(
            input_path=test_file,
            output_path=None,  # Overwrite in place
            columns=["cell_type_high_resolution"],
            remove_raw=True
        )
        
        # Verify
        logger.info("Verifying result...")
        adata = ad.read_h5ad(test_file)
        if "cell_type_high_resolution" in adata.obs.columns:
            is_cat = adata.obs["cell_type_high_resolution"].dtype.name == "category"
            logger.info(f"  cell_type_high_resolution is categorical: {is_cat}")
        logger.info(f"  Has raw data: {adata.raw is not None}")
        
        logger.info("✓ Category fix successful")
        test_file.unlink()  # Clean up
        return True
        
    except Exception as e:
        logger.error(f"✗ Failed to fix categories: {e}", exc_info=True)
        return False


def test_add_metadata_mit():
    """Test 3: Add metadata to MIT data."""
    logger.info("\n" + "="*70)
    logger.info("TEST 3: Add metadata module (MIT data)")
    logger.info("="*70)
    
    try:
        # Make a copy to test on
        test_file = TEST_DIR / "immune_cells_with_metadata.h5ad"
        logger.info(f"Creating test copy: {test_file}")
        shutil.copy(MIT_H5AD, test_file)
        
        # Add metadata
        logger.info(f"Adding metadata from: {MIT_METADATA}")
        add_metadata_to_file(
            h5ad_path=test_file,
            metadata_path=MIT_METADATA,
            is_mit=True,
            output_path=test_file
        )
        
        # Verify
        logger.info("Verifying result...")
        adata = ad.read_h5ad(test_file)
        logger.info(f"  Columns after merge: {len(adata.obs.columns)}")
        logger.info(f"  Sample columns: {list(adata.obs.columns)[:5]}")
        
        logger.info("✓ Metadata addition successful (MIT)")
        test_file.unlink()  # Clean up
        return True
        
    except Exception as e:
        logger.error(f"✗ Failed to add MIT metadata: {e}", exc_info=True)
        return False


def test_column_mapping_mit():
    """Test 4: Column mapping for MIT data."""
    logger.info("\n" + "="*70)
    logger.info("TEST 4: Column mapping module (MIT data)")
    logger.info("="*70)
    
    try:
        # Load file
        logger.info(f"Loading: {MIT_H5AD}")
        adata = ad.read_h5ad(MIT_H5AD)
        logger.info(f"  Original columns: {list(adata.obs.columns)[:5]}")
        
        # Convert columns
        logger.info("Converting ROSMAP_MIT columns...")
        adata = convert_columns(
            adata,
            data_type="ROSMAP_MIT",
            cellclass="Glia",
            subclass="Immune cells"
        )
        
        # Verify
        logger.info("Verifying result...")
        expected_cols = ["Class", "Celltype", "Subtype"]
        for col in expected_cols:
            present = col in adata.obs.columns
            logger.info(f"  {col}: {'✓' if present else '✗'}")
            if present and col == "Class":
                logger.info(f"    Value: {adata.obs[col].unique()}")
        
        logger.info("✓ Column mapping successful (MIT)")
        return True
        
    except Exception as e:
        logger.error(f"✗ Failed to map MIT columns: {e}", exc_info=True)
        return False


def test_full_pipeline_mit():
    """Test 5: Full pipeline for MIT data (small subset)."""
    logger.info("\n" + "="*70)
    logger.info("TEST 5: Full pipeline (MIT data - subset)")
    logger.info("="*70)
    
    try:
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
        import pandas as pd
        metadata = pd.read_csv(MIT_METADATA)
        
        # MIT data: h5ad has 'projid' (int), metadata has 'individualID' (string)
        # Convert projid to string for merging
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
        
        # Step 4: Clean data types for h5ad compatibility
        logger.info("\n  Step 4: Cleaning data types...")
        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                adata.obs[col] = adata.obs[col].astype(str)
        
        # Save
        logger.info(f"\n  Saving to: {MIT_OUTPUT}")
        MIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(MIT_OUTPUT, compression="gzip")
        
        file_size = MIT_OUTPUT.stat().st_size / 1e6
        logger.info(f"  Output size: {file_size:.1f} MB")
        logger.info("✓ Full MIT pipeline successful")
        return True
        
    except Exception as e:
        logger.error(f"✗ Failed full MIT pipeline: {e}", exc_info=True)
        return False


def main():
    """Run all tests."""
    logger.info("="*70)
    logger.info("TESTING REFACTORED MODULES LOCALLY")
    logger.info("="*70)
    
    tests = [
        ("Validation (check_h5ad)", test_check_h5ad),
        ("Category fix", test_category_fix),
        ("Add metadata (MIT)", test_add_metadata_mit),
        ("Column mapping (MIT)", test_column_mapping_mit),
        ("Full pipeline (MIT subset)", test_full_pipeline_mit),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            logger.error(f"Test '{name}' crashed: {e}", exc_info=True)
            results.append((name, False))
    
    # Summary
    logger.info("\n" + "="*70)
    logger.info("TEST SUMMARY")
    logger.info("="*70)
    for name, success in results:
        status = "[PASS]" if success else "[FAIL]"
        logger.info(f"{status}: {name}")
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    logger.info(f"\nPassed: {passed}/{total}")
    
    if passed == total:
        logger.info("\nAll tests passed!")
        return 0
    else:
        logger.info(f"\n{total - passed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
