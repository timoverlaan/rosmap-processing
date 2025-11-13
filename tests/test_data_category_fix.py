"""Unit tests for the category fix module."""

from rosmap_processing.data.category_fix import fix_categories


def test_fix_categories_converts_to_categorical(small_adata):
    """Test that fix_categories converts specified columns to categorical."""
    # Ensure column is not categorical initially
    if "cell_type_high_resolution" in small_adata.obs.columns:
        small_adata.obs["cell_type_high_resolution"] = small_adata.obs["cell_type_high_resolution"].astype(str)
        assert small_adata.obs["cell_type_high_resolution"].dtype != "category"
        
        # Fix categories
        result = fix_categories(small_adata, columns=["cell_type_high_resolution"], remove_raw=False)
        
        # Check it's now categorical
        assert result.obs["cell_type_high_resolution"].dtype.name == "category"


def test_fix_categories_removes_raw(small_adata):
    """Test that fix_categories removes raw data when requested."""
    # Add some raw data if not present
    if small_adata.raw is None:
        small_adata.raw = small_adata.copy()
    
    assert small_adata.raw is not None
    
    # Fix with remove_raw=True
    result = fix_categories(small_adata, columns=[], remove_raw=True)
    
    assert result.raw is None


def test_fix_categories_keeps_raw(small_adata):
    """Test that fix_categories keeps raw data when requested."""
    # Add some raw data if not present
    if small_adata.raw is None:
        small_adata.raw = small_adata.copy()
    
    # Fix with remove_raw=False
    result = fix_categories(small_adata, columns=[], remove_raw=False)
    
    assert result.raw is not None


def test_fix_categories_missing_column(small_adata):
    """Test that fix_categories handles missing columns gracefully."""
    # This should not raise an error
    result = fix_categories(small_adata, columns=["nonexistent_column"], remove_raw=False)
    
    # Should still return the adata
    assert result is not None


def test_fix_categories_already_categorical(small_adata):
    """Test handling of already-categorical columns."""
    if "cell_type_high_resolution" in small_adata.obs.columns:
        # Make it categorical first
        small_adata.obs["cell_type_high_resolution"] = small_adata.obs["cell_type_high_resolution"].astype("category")
        
        # Should not fail
        result = fix_categories(small_adata, columns=["cell_type_high_resolution"], remove_raw=False)
        
        # Should still be categorical
        assert result.obs["cell_type_high_resolution"].dtype.name == "category"
