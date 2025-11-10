"""Unit tests for the validation module."""

import pytest
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from rosmap_processing.utils.validation import (
    check_h5ad_structure,
    check_obs_metadata,
    check_data_matrix
)
from rosmap_processing.core.combine import validate_h5ad_files


def test_check_h5ad_structure(small_adata):
    """Test checking AnnData structure."""
    result = check_h5ad_structure(small_adata, verbose=False)
    
    assert "shape" in result
    assert "n_obs" in result
    assert "n_vars" in result
    assert "obs_columns" in result
    assert "var_columns" in result
    assert "has_raw" in result
    
    assert result["n_obs"] == small_adata.n_obs
    assert result["n_vars"] == small_adata.n_vars


def test_check_obs_metadata(small_adata):
    """Test checking observation metadata."""
    result = check_obs_metadata(small_adata, show_counts=False)
    
    assert isinstance(result, dict)
    # Check that it has entries for each obs column
    for col in small_adata.obs.columns:
        assert col in result
        assert "dtype" in result[col]
        assert "n_unique" in result[col]


def test_check_data_matrix(small_adata):
    """Test checking data matrix."""
    result = check_data_matrix(small_adata, preview_size=5, show_raw=False)
    
    assert isinstance(result, dict)
    assert "X" in result
    assert result["X"]["shape"] == small_adata.X.shape


def test_validate_h5ad_files_nonexistent():
    """Test validation fails for nonexistent files."""
    with pytest.raises(FileNotFoundError):
        validate_h5ad_files([Path("nonexistent.h5ad")])


def test_validate_h5ad_files_wrong_extension(temp_output_dir):
    """Test validation fails for wrong extension."""
    # Create a dummy file with wrong extension
    test_file = temp_output_dir / "test.txt"
    test_file.write_text("dummy")
    
    with pytest.raises(ValueError, match="must have .h5ad extension"):
        validate_h5ad_files([test_file])


def test_validate_h5ad_files_valid(mit_h5ad_path):
    """Test validation passes for valid file."""
    # Should not raise
    validate_h5ad_files([mit_h5ad_path])
