"""Unit tests for the constants module."""

from rosmap_processing.utils.constants import (
    MIN_GENES_PER_CELL,
    MIN_CELLS_PER_GENE,
    DEFAULT_HVG_COUNT,
    BRAAK_MAPPING,
    CERAD_MAPPING,
    THAL_MAPPING,
    ROSMAP_MIT_CELL_CLASSES,
    SYNAPSE_IDS_ROSMAP,
    SYNAPSE_IDS_ROSMAP_MIT
)


def test_constants_exist():
    """Test that all constants are defined."""
    assert isinstance(MIN_GENES_PER_CELL, int)
    assert isinstance(MIN_CELLS_PER_GENE, int)
    assert isinstance(DEFAULT_HVG_COUNT, int)


def test_constants_reasonable_values():
    """Test that constants have reasonable values."""
    assert MIN_GENES_PER_CELL > 0
    assert MIN_CELLS_PER_GENE > 0
    assert DEFAULT_HVG_COUNT > 0
    assert DEFAULT_HVG_COUNT <= 10000  # Reasonable upper bound


def test_braak_mapping_structure():
    """Test Braak mapping has correct structure."""
    assert isinstance(BRAAK_MAPPING, dict)
    assert len(BRAAK_MAPPING) > 0
    # Check all values are integers
    for k, v in BRAAK_MAPPING.items():
        assert isinstance(v, int)
        assert 0 <= v <= 6


def test_cerad_mapping_structure():
    """Test CERAD mapping has correct structure."""
    assert isinstance(CERAD_MAPPING, dict)
    assert len(CERAD_MAPPING) > 0
    # Check expected keys
    assert "Absent" in CERAD_MAPPING
    assert "Frequent" in CERAD_MAPPING
    # Values should be 1-4
    for v in CERAD_MAPPING.values():
        assert 1 <= v <= 4


def test_thal_mapping_structure():
    """Test Thal mapping has correct structure."""
    assert isinstance(THAL_MAPPING, dict)
    assert len(THAL_MAPPING) > 0
    # Check all values are integers 0-6
    for k, v in THAL_MAPPING.items():
        assert isinstance(v, int)
        assert 0 <= v <= 6


def test_rosmap_mit_cell_classes():
    """Test ROSMAP MIT cell classes mapping."""
    assert isinstance(ROSMAP_MIT_CELL_CLASSES, dict)
    assert len(ROSMAP_MIT_CELL_CLASSES) > 0
    # Check expected cell types
    expected_keys = ["Astrocytes", "Immune_cells", "Oligodendrocytes"]
    for key in expected_keys:
        assert key in ROSMAP_MIT_CELL_CLASSES


def test_synapse_ids_rosmap():
    """Test ROSMAP Synapse IDs."""
    assert isinstance(SYNAPSE_IDS_ROSMAP, dict)
    assert len(SYNAPSE_IDS_ROSMAP) > 0
    # Check all values are strings starting with 'syn'
    for k, v in SYNAPSE_IDS_ROSMAP.items():
        assert isinstance(v, str)
        assert v.startswith("syn")


def test_synapse_ids_rosmap_mit():
    """Test ROSMAP MIT Synapse IDs."""
    assert isinstance(SYNAPSE_IDS_ROSMAP_MIT, dict)
    assert len(SYNAPSE_IDS_ROSMAP_MIT) > 0
    # Check all values are strings starting with 'syn'
    for k, v in SYNAPSE_IDS_ROSMAP_MIT.items():
        assert isinstance(v, str)
        assert v.startswith("syn")
