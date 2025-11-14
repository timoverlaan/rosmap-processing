"""Unit tests for the combine module."""

import pytest
import anndata as ad
import numpy as np
from pathlib import Path
from rosmap_processing.core.combine import validate_h5ad_files, combine_h5ad_files


def test_validate_h5ad_files_valid(mit_h5ad_path):
    """Test validation of valid h5ad files."""
    # Should not raise
    validate_h5ad_files([mit_h5ad_path])


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


def test_combine_h5ad_files_basic(temp_output_dir):
    """Test basic combination of h5ad files."""
    # Create two small test files
    adata1 = ad.AnnData(
        X=np.random.randn(50, 100),
        obs={"cell_id": [f"cell_{i}" for i in range(50)]},
        var={"gene": [f"gene_{i}" for i in range(100)]}
    )
    adata2 = ad.AnnData(
        X=np.random.randn(30, 100),
        obs={"cell_id": [f"cell_{i}" for i in range(50, 80)]},
        var={"gene": [f"gene_{i}" for i in range(100)]}
    )
    
    file1 = temp_output_dir / "test1.h5ad"
    file2 = temp_output_dir / "test2.h5ad"
    
    adata1.write_h5ad(file1)
    adata2.write_h5ad(file2)
    
    # Combine
    result = combine_h5ad_files([file1, file2])
    
    assert result.n_obs == 80  # 50 + 30
    assert result.n_vars == 100
    assert "source_file" in result.obs.columns


def test_combine_h5ad_files_different_genes(temp_output_dir):
    """Test combination with different gene sets."""
    # Create files with overlapping genes
    adata1 = ad.AnnData(
        X=np.random.randn(50, 100),
        obs={"cell_id": [f"cell_{i}" for i in range(50)]},
        var={"gene": [f"gene_{i}" for i in range(100)]}
    )
    adata2 = ad.AnnData(
        X=np.random.randn(30, 80),
        obs={"cell_id": [f"cell_{i}" for i in range(50, 80)]},
        var={"gene": [f"gene_{i}" for i in range(20, 100)]}  # Different genes
    )
    
    file1 = temp_output_dir / "test1.h5ad"
    file2 = temp_output_dir / "test2.h5ad"
    
    adata1.write_h5ad(file1)
    adata2.write_h5ad(file2)
    
    # Combine with outer join (default)
    result = combine_h5ad_files([file1, file2], join="outer")
    
    assert result.n_obs == 80
    assert result.n_vars == 100  # Union of genes


def test_combine_h5ad_files_empty_list():
    """Test that combining empty list raises error."""
    with pytest.raises(ValueError, match="No valid AnnData objects"):
        combine_h5ad_files([])
