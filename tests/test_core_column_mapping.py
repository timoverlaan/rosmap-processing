"""Unit tests for column mapping module."""

import pytest
import anndata as ad
import pandas as pd
from rosmap_processing.core.column_mapping import (
    convert_rosmap_columns,
    convert_rosmap_mit_columns,
    convert_seaad_columns,
    derive_wang_labels
)


def test_convert_rosmap_mit_columns_adds_class_and_subclass(small_adata):
    """Test that MIT conversion adds Class and Celltype columns."""
    result = convert_rosmap_mit_columns(
        small_adata,
        cellclass="Glia",
        subclass="Immune cells"
    )
    
    assert "Class" in result.obs.columns
    assert "Celltype" in result.obs.columns
    assert all(result.obs["Class"] == "Glia")
    assert all(result.obs["Celltype"] == "Immune cells")


def test_convert_rosmap_mit_columns_renames_columns(small_adata):
    """Test that MIT conversion renames columns correctly."""
    # Add required columns if not present
    if "projid" not in small_adata.obs.columns:
        small_adata.obs["projid"] = range(len(small_adata))
    if "cell_type_high_resolution" not in small_adata.obs.columns:
        small_adata.obs["cell_type_high_resolution"] = "Type1"
    
    result = convert_rosmap_mit_columns(
        small_adata,
        cellclass="Glia",
        subclass="Immune cells"
    )
    
    assert "Donor ID" in result.obs.columns
    assert "Subtype" in result.obs.columns


def test_convert_seaad_columns_creates_numeric_scores(small_adata):
    """Test that SeaAD conversion creates numeric score columns."""
    # Add required columns - match length exactly
    n = len(small_adata)
    small_adata.obs["Braak"] = (["Braak 0", "Braak I", "Braak II"] * (n // 3 + 1))[:n]
    small_adata.obs["CERAD"] = (["Absent", "Sparse", "Moderate"] * (n // 3 + 1))[:n]
    small_adata.obs["Thal"] = (["Thal 0", "Thal 1", "Thal 2"] * (n // 3 + 1))[:n]
    small_adata.obs["Cognitive Status"] = "Normal"
    
    result = convert_seaad_columns(small_adata)
    
    assert "braaksc" in result.obs.columns
    assert "ceradsc" in result.obs.columns
    assert "thalsc" in result.obs.columns
    assert pd.api.types.is_numeric_dtype(result.obs["braaksc"])


def test_derive_wang_labels_creates_columns(small_adata):
    """Test that Wang label derivation creates required columns."""
    # Add required columns
    small_adata.obs["Cognitive Status"] = ["Dementia"] * 50 + ["Normal"] * 50
    small_adata.obs["Cognitive Status"] = small_adata.obs["Cognitive Status"][:len(small_adata)]
    small_adata.obs["braaksc"] = [5] * 50 + [2] * 50
    small_adata.obs["braaksc"] = small_adata.obs["braaksc"][:len(small_adata)]
    small_adata.obs["ceradsc"] = [1] * 50 + [4] * 50
    small_adata.obs["ceradsc"] = small_adata.obs["ceradsc"][:len(small_adata)]
    
    derive_wang_labels(small_adata)
    
    assert "Wang" in small_adata.obs.columns
    assert "Wang_intermediate" in small_adata.obs.columns
    assert set(small_adata.obs["Wang"].unique()).issubset({"AD", "Healthy", "Intermediate"})


def test_derive_wang_labels_identifies_ad_cases(small_adata):
    """Test that Wang labels correctly identify AD cases."""
    # Create clear AD case: Dementia + Braak >= 4 + CERAD <= 2
    small_adata.obs["Cognitive Status"] = "Dementia"
    small_adata.obs["braaksc"] = 5
    small_adata.obs["ceradsc"] = 1
    
    derive_wang_labels(small_adata)
    
    assert all(small_adata.obs["Wang"] == "AD")
    assert all(small_adata.obs["Wang_intermediate"] == False)


def test_derive_wang_labels_identifies_healthy_cases(small_adata):
    """Test that Wang labels correctly identify Healthy cases."""
    # Create clear Healthy case: No dementia + Braak <= 3 + CERAD >= 3
    small_adata.obs["Cognitive Status"] = "Normal"
    small_adata.obs["braaksc"] = 2
    small_adata.obs["ceradsc"] = 4
    
    derive_wang_labels(small_adata)
    
    assert all(small_adata.obs["Wang"] == "Healthy")
    assert all(small_adata.obs["Wang_intermediate"] == False)


def test_derive_wang_labels_missing_columns(small_adata):
    """Test that derive_wang_labels handles missing columns gracefully."""
    # Remove required columns if present
    for col in ["Cognitive Status", "braaksc", "ceradsc"]:
        if col in small_adata.obs.columns:
            del small_adata.obs[col]
    
    # Should not raise, just warn
    derive_wang_labels(small_adata)
    
    # Should not have added Wang columns
    assert "Wang" not in small_adata.obs.columns
