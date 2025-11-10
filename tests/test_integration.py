"""Integration tests using real data files.

These tests require actual data files to be present and may take longer to run.
"""

import pytest
import shutil
from pathlib import Path
from rosmap_processing.data.category_fix import fix_categories_in_file
from rosmap_processing.data.metadata import add_metadata_to_file
from rosmap_processing.core.column_mapping import convert_columns
import anndata as ad


@pytest.mark.integration
def test_category_fix_integration(mit_h5ad_path, temp_output_dir):
    """Test category fix on real data."""
    # Copy file to temp location
    test_file = temp_output_dir / "test.h5ad"
    shutil.copy(mit_h5ad_path, test_file)
    
    # Fix categories
    fix_categories_in_file(
        input_path=test_file,
        output_path=None,
        columns=["cell_type_high_resolution"],
        remove_raw=True
    )
    
    # Verify
    adata = ad.read_h5ad(test_file)
    if "cell_type_high_resolution" in adata.obs.columns:
        assert adata.obs["cell_type_high_resolution"].dtype.name == "category"
    assert adata.raw is None


@pytest.mark.integration
def test_add_metadata_integration_mit(mit_h5ad_path, mit_metadata_path, temp_output_dir):
    """Test adding metadata to MIT data."""
    # Copy file to temp location
    test_file = temp_output_dir / "test.h5ad"
    shutil.copy(mit_h5ad_path, test_file)
    
    # Add metadata
    add_metadata_to_file(
        h5ad_path=test_file,
        metadata_path=mit_metadata_path,
        is_mit=True,
        output_path=test_file
    )
    
    # Verify
    adata = ad.read_h5ad(test_file)
    # Should have more columns after adding metadata
    assert len(adata.obs.columns) > 2
    # Check some expected metadata columns
    expected_cols = ["sex", "ageDeath", "diagnosis"]
    for col in expected_cols:
        if col in adata.obs.columns:
            assert adata.obs[col].notna().any()


@pytest.mark.integration
def test_column_mapping_integration_mit(mit_h5ad_path):
    """Test column mapping on MIT data."""
    # Load data
    adata = ad.read_h5ad(mit_h5ad_path)
    
    # Take small subset for speed
    adata = adata[:1000, :].copy()
    
    # Convert columns
    result = convert_columns(
        adata,
        data_type="ROSMAP_MIT",
        cellclass="Glia",
        subclass="Immune cells"
    )
    
    # Verify standard columns are present
    assert "Class" in result.obs.columns
    assert "Celltype" in result.obs.columns
    assert "Donor ID" in result.obs.columns
    assert "Subtype" in result.obs.columns
    
    # Check values
    assert all(result.obs["Class"] == "Glia")
    assert all(result.obs["Celltype"] == "Immune cells")


@pytest.mark.integration
def test_full_pipeline_mit(mit_h5ad_path, mit_metadata_path, temp_output_dir):
    """Test full processing pipeline on MIT data."""
    import pandas as pd
    
    # Load data (small subset for speed)
    adata = ad.read_h5ad(mit_h5ad_path)
    adata = adata[:500, :].copy()
    
    # Step 1: Fix categories
    if "cell_type_high_resolution" in adata.obs.columns:
        adata.obs["cell_type_high_resolution"] = adata.obs["cell_type_high_resolution"].astype("category")
    if adata.raw is not None:
        del adata.raw
    
    # Step 2: Add metadata
    metadata = pd.read_csv(mit_metadata_path)
    adata.obs["projid"] = adata.obs["projid"].astype(str)
    metadata["individualID"] = metadata["individualID"].astype(str)
    adata.obs = adata.obs.merge(
        metadata, left_on="projid", right_on="individualID", how="left"
    )
    adata.obs.index = adata.obs_names
    
    # Step 3: Convert columns
    if "projid" in adata.obs.columns:
        adata.obs.rename(columns={"projid": "Donor ID"}, inplace=True)
    if "cell_type_high_resolution" in adata.obs.columns:
        adata.obs.rename(columns={"cell_type_high_resolution": "Subtype"}, inplace=True)
    adata.obs["Class"] = "Glia"
    adata.obs["Celltype"] = "Immune cells"
    
    # Step 4: Clean data types
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].astype(str)
    
    # Save
    output_file = temp_output_dir / "processed.h5ad"
    adata.write_h5ad(output_file, compression="gzip")
    
    # Verify
    assert output_file.exists()
    result = ad.read_h5ad(output_file)
    assert result.n_obs == 500
    assert "Class" in result.obs.columns
    assert "Celltype" in result.obs.columns
