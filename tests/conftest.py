"""Test configuration and fixtures for pytest."""

import pytest
from pathlib import Path
import shutil
import tempfile


@pytest.fixture(scope="session")
def test_data_dir():
    """Path to test data directory."""
    return Path("data/raw")


@pytest.fixture(scope="session")
def mit_h5ad_path(test_data_dir):
    """Path to MIT test h5ad file."""
    path = test_data_dir / "ROSMAP_MIT" / "immune_cells.h5ad"
    if not path.exists():
        pytest.skip(f"Test data not found: {path}")
    return path


@pytest.fixture(scope="session")
def mit_metadata_path(test_data_dir):
    """Path to MIT metadata file."""
    path = test_data_dir / "ROSMAP_MIT" / "MIT_ROSMAP_Multiomics_individual_metadata.csv"
    if not path.exists():
        pytest.skip(f"Test metadata not found: {path}")
    return path


@pytest.fixture(scope="session")
def rosmap_h5ad_path(test_data_dir):
    """Path to ROSMAP test h5ad file."""
    path = test_data_dir / "ROSMAP" / "astrocytes.h5ad"
    if not path.exists():
        pytest.skip(f"Test data not found: {path}")
    return path


@pytest.fixture(scope="session")
def rosmap_metadata_path(test_data_dir):
    """Path to ROSMAP metadata file."""
    path = test_data_dir / "ROSMAP" / "rosmap_clinical.csv"
    if not path.exists():
        pytest.skip(f"Test metadata not found: {path}")
    return path


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    # Cleanup
    if temp_dir.exists():
        shutil.rmtree(temp_dir)


@pytest.fixture
def small_adata(mit_h5ad_path):
    """Create a small AnnData object for fast tests (100 cells)."""
    import anndata as ad
    adata = ad.read_h5ad(mit_h5ad_path)
    # Take small subset
    return adata[:100, :].copy()


@pytest.fixture
def medium_adata(mit_h5ad_path):
    """Create a medium AnnData object for integration tests (1000 cells)."""
    import anndata as ad
    adata = ad.read_h5ad(mit_h5ad_path)
    # Take medium subset
    return adata[:1000, :].copy()
