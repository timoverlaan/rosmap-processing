# Repository Restructuring Plan

**Date:** November 10, 2025  
**Repository:** rosmap-processing  
**Status:** PROPOSAL - Not yet implemented

---

## Overview

This document proposes a restructuring of the rosmap-processing repository to improve maintainability, clarity, and ease of onboarding for new contributors.

---

## Proposed New Directory Structure

```
rosmap-processing/
├── README.md                          # Comprehensive project documentation
├── CONTRIBUTING.md                    # Guidelines for contributors
├── CHANGELOG.md                       # Version history and changes
├── LICENSE                            # Project license
├── pyproject.toml                     # Python project metadata
├── pixi.toml                          # Pixi dependency management
├── pixi.lock                          # Locked dependencies
├── .gitignore                         # Git ignore patterns
├── .env.example                       # Example environment variables
│
├── docs/                              # Documentation
│   ├── installation.md                # Setup instructions
│   ├── usage.md                       # Usage examples
│   ├── pipeline.md                    # Pipeline description with diagrams
│   ├── data_sources.md                # Data source documentation
│   ├── troubleshooting.md             # Common issues and solutions
│   └── api_reference.md               # Code API documentation
│
├── src/                               # Source code (Python package)
│   ├── rosmap_processing/             # Main package
│   │   ├── __init__.py
│   │   ├── __version__.py             # Version info
│   │   │
│   │   ├── core/                      # Core processing modules
│   │   │   ├── __init__.py
│   │   │   ├── scanpy_pipeline.py     # Scanpy processing logic
│   │   │   ├── metadata.py            # Metadata handling
│   │   │   ├── combine.py             # Data combination
│   │   │   └── column_mapping.py      # Column name conversion
│   │   │
│   │   ├── data/                      # Data handling
│   │   │   ├── __init__.py
│   │   │   ├── download.py            # Data download from Synapse
│   │   │   ├── loaders.py             # Data loading utilities
│   │   │   └── validators.py          # Input validation
│   │   │
│   │   ├── conversion/                # Format conversion
│   │   │   ├── __init__.py
│   │   │   ├── r_to_h5ad.py          # R-to-h5ad conversion wrapper
│   │   │   └── category_fix.py        # Category fixing utilities
│   │   │
│   │   ├── analysis/                  # Analysis modules
│   │   │   ├── __init__.py
│   │   │   └── deg_analysis.py        # DEG analysis (Mathys 2019)
│   │   │
│   │   ├── utils/                     # Utilities
│   │   │   ├── __init__.py
│   │   │   ├── logging.py             # Logging configuration
│   │   │   ├── config.py              # Configuration management
│   │   │   ├── constants.py           # Constants and magic numbers
│   │   │   └── io.py                  # I/O utilities
│   │   │
│   │   └── cli/                       # Command-line interfaces
│   │       ├── __init__.py
│   │       ├── main.py                # Main CLI entry point
│   │       ├── process.py             # Processing commands
│   │       ├── download.py            # Download commands
│   │       └── convert.py             # Conversion commands
│   │
│   └── r_scripts/                     # R scripts
│       ├── rds_to_h5seurat.R
│       ├── h5seurat_to_h5ad.R
│       └── install_deps.R
│
├── scripts/                           # Operational scripts
│   ├── pipelines/                     # Full pipeline runners
│   │   ├── run_rosmap_pipeline.sh
│   │   ├── run_seaad_pipeline.sh
│   │   └── run_full_pipeline.sh
│   │
│   └── utilities/                     # Utility scripts
│       ├── check_data.sh
│       ├── cleanup.sh
│       └── verify_installation.sh
│
├── slurm/                            # SLURM job submission scripts
│   ├── README.md                     # SLURM usage documentation
│   ├── templates/                    # Job templates
│   │   └── base_job.template.sh
│   ├── jobs/                         # Specific job scripts
│   │   ├── rosmap_mathys_genes.sh
│   │   ├── rosmap_top1k_genes.sh
│   │   ├── seaad_rosmap_genes.sh
│   │   ├── run_all.sh
│   │   └── run_all_seaad.sh
│   └── logs/                         # SLURM output logs (gitignored)
│
├── containers/                       # Container definitions
│   ├── README.md                     # Container documentation
│   ├── pixi/
│   │   ├── pixi.def                  # Singularity definition
│   │   └── CHANGELOG.md              # Container version history
│   └── deseq2/
│       └── deseq2.def
│
├── config/                           # Configuration files
│   ├── default_config.yaml           # Default configuration
│   ├── rosmap_config.yaml            # ROSMAP-specific config
│   └── seaad_config.yaml             # SeaAD-specific config
│
├── tests/                            # Test suite
│   ├── __init__.py
│   ├── conftest.py                   # Pytest configuration
│   ├── fixtures/                     # Test fixtures
│   │   └── sample_data/
│   ├── unit/                         # Unit tests
│   │   ├── test_scanpy_pipeline.py
│   │   ├── test_metadata.py
│   │   ├── test_combine.py
│   │   └── test_validators.py
│   ├── integration/                  # Integration tests
│   │   ├── test_full_pipeline.py
│   │   └── test_conversion.py
│   └── test_data/                    # Small test datasets
│
├── data/                             # Data directory (gitignored)
│   ├── raw/                          # Raw downloaded data
│   │   ├── ROSMAP/
│   │   └── ROSMAP_MIT/
│   ├── processed/                    # Processed data
│   ├── interim/                      # Intermediate files
│   └── metadata/                     # Metadata files
│
├── notebooks/                        # Jupyter notebooks
│   ├── 01_data_exploration.ipynb
│   ├── 02_qa_checks.ipynb
│   └── README.md
│
├── output/                           # Output directory (gitignored)
│   └── .gitkeep
│
└── .github/                          # GitHub-specific files
    ├── workflows/                    # CI/CD workflows
    │   ├── tests.yml
    │   └── build_container.yml
    └── ISSUE_TEMPLATE/
```

---

## Key Structural Changes

### 1. **Convert to Proper Python Package**

**Current:** Loose collection of scripts  
**Proposed:** Structured package with `__init__.py` files

**Benefits:**
- Enables `pip install -e .` for development
- Allows importing between modules
- Better code organization and reusability
- Enables proper testing with pytest

**Implementation:**
```python
# src/rosmap_processing/__init__.py
"""ROSMAP data processing pipeline."""

from .core.scanpy_pipeline import process_h5ad
from .data.download import download_rosmap_data
from .__version__ import __version__

__all__ = ["process_h5ad", "download_rosmap_data", "__version__"]
```

### 2. **Unified CLI with Subcommands**

**Current:** Multiple standalone scripts  
**Proposed:** Single entry point with subcommands

**Example usage:**
```bash
# Instead of:
python src/download_rosmap.py
python src/scanpy_pipeline.py input.h5ad output.h5ad --n-genes 2000

# Use:
rosmap-process download --dataset rosmap-mit
rosmap-process analyze input.h5ad --output output.h5ad --n-genes 2000
rosmap-process convert --format h5seurat-to-h5ad input.h5Seurat
rosmap-process pipeline --config config/rosmap_config.yaml
```

**Implementation:**
```python
# src/rosmap_processing/cli/main.py
import click
from ..core import scanpy_pipeline
from ..data import download

@click.group()
@click.version_option()
def cli():
    """ROSMAP data processing pipeline."""
    pass

@cli.command()
@click.option('--dataset', type=click.Choice(['rosmap', 'rosmap-mit', 'seaad']))
def download(dataset):
    """Download datasets from Synapse or S3."""
    # Implementation

@cli.command()
@click.argument('input_file')
@click.option('--output', required=True)
@click.option('--n-genes', default=2000)
def analyze(input_file, output, n_genes):
    """Run scanpy analysis pipeline."""
    # Implementation

if __name__ == '__main__':
    cli()
```

### 3. **Configuration File System**

**Current:** Hardcoded parameters  
**Proposed:** YAML/TOML configuration files

**Example configuration:**
```yaml
# config/rosmap_config.yaml
dataset:
  name: "ROSMAP"
  type: "single-cell-rnaseq"

processing:
  min_genes: 200
  min_cells: 5
  n_hvgs: 2000
  k_neighbors: 30
  normalization:
    target_sum: 1000000  # CPM
    log_transform: true

paths:
  raw_data: "data/raw/ROSMAP"
  processed: "data/processed"
  metadata: "data/metadata"

synapse:
  use_env_token: true  # Read from SYNAPSE_AUTH_TOKEN env var
  
logging:
  level: "INFO"
  format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
```

### 4. **Centralized Logging**

**Current:** Print statements everywhere  
**Proposed:** Python logging module

**Implementation:**
```python
# src/rosmap_processing/utils/logging.py
import logging
import sys
from pathlib import Path

def setup_logging(level: str = "INFO", log_file: Path = None):
    """Configure logging for the application."""
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    
    # File handler (optional)
    handlers = [console_handler]
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Configure root logger
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        handlers=handlers
    )
    
    return logging.getLogger("rosmap_processing")

# Usage in modules:
# from rosmap_processing.utils.logging import setup_logging
# logger = setup_logging()
# logger.info("Processing started")
# logger.warning("Potential issue detected")
```

### 5. **Shared Constants and Utilities**

**Current:** Magic numbers and duplicated code  
**Proposed:** Centralized constants and utility functions

```python
# src/rosmap_processing/utils/constants.py
"""Constants used throughout the pipeline."""

# Processing thresholds
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 5
DEFAULT_HVG_COUNT = 2000
DEFAULT_K_NEIGHBORS = 30
DEFAULT_PCA_COMPONENTS = 50

# Normalization
CPM_TARGET = 1e6  # Normalize to 1 million reads per cell

# Column names
CLASS_NAME = "Class"
SUBCLASS_NAME = "Celltype"
SUPERTYPE_NAME = "Subtype"
DONOR_ID_COLUMN = "Donor ID"

# Synapse IDs
SYNAPSE_IDS = {
    "rosmap": {
        "astrocytes": "syn53693925",
        "cux2_plus": "syn53694054",
        "cux2_minus": "syn53694068",
        # ... etc
    },
    "rosmap_mit": {
        "astrocytes": "syn52368912",
        # ... etc
    }
}

# File extensions
VALID_H5AD_EXTENSIONS = [".h5ad"]
VALID_H5SEURAT_EXTENSIONS = [".h5Seurat"]
VALID_RDS_EXTENSIONS = [".rds"]
```

```python
# src/rosmap_processing/utils/io.py
"""I/O utilities for loading and saving data."""

import anndata as ad
from pathlib import Path
from typing import Union
import logging

logger = logging.getLogger(__name__)

def load_h5ad(path: Union[str, Path], backed: bool = False) -> ad.AnnData:
    """
    Load an h5ad file with validation and error handling.
    
    Parameters
    ----------
    path : str or Path
        Path to the h5ad file
    backed : bool, default False
        Whether to load in backed mode
        
    Returns
    -------
    AnnData
        Loaded AnnData object
        
    Raises
    ------
    FileNotFoundError
        If the file doesn't exist
    ValueError
        If the file is not a valid h5ad file
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    if path.suffix != ".h5ad":
        raise ValueError(f"Expected .h5ad file, got: {path.suffix}")
    
    logger.info(f"Loading h5ad file: {path}")
    
    try:
        adata = ad.read_h5ad(path, backed=backed)
        logger.info(f"Loaded AnnData object: {adata.shape}")
        return adata
    except Exception as e:
        logger.error(f"Failed to load {path}: {e}")
        raise
```

### 6. **Input Validation Module**

```python
# src/rosmap_processing/data/validators.py
"""Input validation utilities."""

import anndata as ad
from pathlib import Path
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)

class ValidationError(Exception):
    """Raised when validation fails."""
    pass

def validate_adata_columns(
    adata: ad.AnnData,
    required_columns: List[str],
    location: str = "obs"
) -> None:
    """
    Validate that required columns exist in AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        The AnnData object to validate
    required_columns : List[str]
        List of required column names
    location : str
        Where to look for columns: 'obs', 'var', or 'uns'
        
    Raises
    ------
    ValidationError
        If any required columns are missing
    """
    if location == "obs":
        available = adata.obs.columns
    elif location == "var":
        available = adata.var.columns
    else:
        raise ValueError(f"Invalid location: {location}")
    
    missing = [col for col in required_columns if col not in available]
    
    if missing:
        raise ValidationError(
            f"Missing required columns in {location}: {missing}\n"
            f"Available columns: {list(available)}"
        )
    
    logger.info(f"Validation passed: all required {location} columns present")

def validate_file_path(path: Path, must_exist: bool = True) -> Path:
    """
    Validate a file path.
    
    Parameters
    ----------
    path : Path
        Path to validate
    must_exist : bool
        Whether the file must already exist
        
    Returns
    -------
    Path
        Validated path
        
    Raises
    ------
    ValidationError
        If validation fails
    """
    path = Path(path)
    
    if must_exist and not path.exists():
        raise ValidationError(f"File does not exist: {path}")
    
    if not must_exist:
        # Check parent directory exists
        if not path.parent.exists():
            raise ValidationError(f"Parent directory does not exist: {path.parent}")
    
    return path
```

### 7. **Improved Pipeline Scripts**

```bash
# scripts/pipelines/run_rosmap_pipeline.sh
#!/bin/bash
set -euo pipefail  # Exit on error, undefined variables, pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${PROJECT_ROOT}/config/rosmap_config.yaml"
LOG_DIR="${PROJECT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

# Logging
LOG_FILE="${LOG_DIR}/rosmap_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=========================================="
echo "ROSMAP Processing Pipeline"
echo "Started: $(date)"
echo "Config: ${CONFIG_FILE}"
echo "=========================================="

# Load environment
cd "${PROJECT_ROOT}"
eval "$(pixi shell-hook)"

# Function to handle errors
handle_error() {
    echo "ERROR: Pipeline failed at step: $1" >&2
    echo "Check log file: ${LOG_FILE}" >&2
    exit 1
}

# Step 1: Download data
echo -e "\n[Step 1/5] Downloading ROSMAP data..."
rosmap-process download --dataset rosmap --config "${CONFIG_FILE}" || handle_error "download"

# Step 2: Convert formats
echo -e "\n[Step 2/5] Converting data formats..."
rosmap-process convert --format rds-to-h5seurat --input-dir data/raw/ROSMAP || handle_error "convert-rds"
rosmap-process convert --format h5seurat-to-h5ad --input-dir data/raw/ROSMAP || handle_error "convert-h5ad"

# Step 3: Add metadata
echo -e "\n[Step 3/5] Adding metadata..."
rosmap-process metadata add --config "${CONFIG_FILE}" || handle_error "metadata"

# Step 4: Process with scanpy
echo -e "\n[Step 4/5] Running scanpy pipeline..."
rosmap-process analyze --config "${CONFIG_FILE}" || handle_error "analyze"

# Step 5: QA checks
echo -e "\n[Step 5/5] Running QA checks..."
rosmap-process qa --config "${CONFIG_FILE}" || handle_error "qa"

echo -e "\n=========================================="
echo "Pipeline completed successfully!"
echo "Finished: $(date)"
echo "Log file: ${LOG_FILE}"
echo "=========================================="
```

### 8. **Testing Infrastructure**

```python
# tests/conftest.py
"""Pytest configuration and fixtures."""

import pytest
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
import tempfile
import shutil

@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    temp = tempfile.mkdtemp()
    yield Path(temp)
    shutil.rmtree(temp)

@pytest.fixture
def sample_adata():
    """Create a small sample AnnData object for testing."""
    n_obs = 100
    n_vars = 50
    
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
    
    obs = pd.DataFrame({
        'cell_id': [f'cell_{i}' for i in range(n_obs)],
        'Donor ID': np.random.choice(['donor1', 'donor2', 'donor3'], n_obs),
        'Class': np.random.choice(['Neuron', 'Glia'], n_obs),
    }, index=[f'cell_{i}' for i in range(n_obs)])
    
    var = pd.DataFrame({
        'gene_name': [f'gene_{i}' for i in range(n_vars)],
    }, index=[f'gene_{i}' for i in range(n_vars)])
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata

@pytest.fixture
def sample_h5ad_file(sample_adata, temp_dir):
    """Create a temporary h5ad file."""
    path = temp_dir / "test.h5ad"
    sample_adata.write_h5ad(path)
    return path
```

```python
# tests/unit/test_validators.py
"""Tests for validation utilities."""

import pytest
import anndata as ad
from rosmap_processing.data.validators import (
    validate_adata_columns,
    ValidationError
)

def test_validate_adata_columns_success(sample_adata):
    """Test that validation passes with all required columns."""
    required = ['Donor ID', 'Class']
    validate_adata_columns(sample_adata, required, location='obs')
    # Should not raise

def test_validate_adata_columns_missing(sample_adata):
    """Test that validation fails with missing columns."""
    required = ['Donor ID', 'NonexistentColumn']
    
    with pytest.raises(ValidationError) as exc_info:
        validate_adata_columns(sample_adata, required, location='obs')
    
    assert 'NonexistentColumn' in str(exc_info.value)
```

### 9. **Enhanced Documentation**

**README.md structure:**
```markdown
# ROSMAP Processing Pipeline

[![Tests](https://github.com/timoverlaan/rosmap-processing/workflows/tests/badge.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)]()

> Reproducible processing pipeline for ROSMAP Alzheimer's Disease single-cell RNA-seq datasets

## Features

- 🔄 Automated download from Synapse and AWS S3
- 🔬 Scanpy-based single-cell processing
- 📊 Integration of ROSMAP and SeaAD datasets
- 🐳 Containerized for reproducibility
- 🧪 Comprehensive test suite
- 📝 Well-documented API

## Quick Start

### Installation
[detailed steps]

### Basic Usage
[common examples]

### Advanced Usage
[link to docs/]

## Documentation

- [Installation Guide](docs/installation.md)
- [Usage Examples](docs/usage.md)
- [Pipeline Overview](docs/pipeline.md)
- [API Reference](docs/api_reference.md)
- [Troubleshooting](docs/troubleshooting.md)

## Project Structure
[brief overview]

## Citation
[paper citation]

## License
[license info]
```

**Pipeline documentation with diagram:**
```markdown
# docs/pipeline.md

## Pipeline Overview

### ROSMAP Pipeline

```
┌─────────────────┐
│  Download Data  │ ← Synapse API
│  (Synapse)      │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Convert        │ ← RDS → h5Seurat → h5ad
│  Formats        │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Add Metadata   │ ← Join clinical metadata
│                 │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Normalize &    │ ← CPM + log transform
│  Filter         │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Feature        │ ← HVG selection or import
│  Selection      │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  PCA & KNN      │ ← Per-donor or global
│  Graph          │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Save           │ ← Processed h5ad
│  Output         │
└─────────────────┘
```

### Data Flow

1. **Input:** Raw .rds or .h5Seurat files
2. **Intermediate:** .h5ad files with counts
3. **Output:** Processed .h5ad with normalized data, PCA, and KNN graph

[detailed explanation of each step]
```

---

## Migration Strategy

### Phase 1: Structural Changes (Non-Breaking)
1. Create new directory structure
2. Move files to new locations (keep copies in old locations)
3. Add `__init__.py` files
4. Create `pyproject.toml`
5. Add configuration files
6. Update `.gitignore`

### Phase 2: Code Refactoring
1. Extract shared utilities
2. Add constants file
3. Implement logging
4. Add input validation
5. Refactor into modules
6. Add type hints everywhere

### Phase 3: CLI Development
1. Implement CLI with Click
2. Create subcommands
3. Update shell scripts to use CLI
4. Add configuration file support

### Phase 4: Testing & Documentation
1. Write unit tests
2. Write integration tests
3. Add docstrings
4. Write comprehensive docs
5. Add CI/CD workflows

### Phase 5: Cleanup
1. Remove deprecated files
2. Archive old scripts
3. Final testing
4. Release v1.0.0

---

## Implementation Checklist

### Immediate (Week 1)
- [ ] Create new directory structure
- [ ] Add proper `.gitignore` for `data/`, `output/`, `*.sif`, `logs/`, `.env`
- [ ] Create `.env.example` for environment variables
- [ ] Move scripts to appropriate directories
- [ ] Create `pyproject.toml` with basic metadata

### Short-term (Week 2-3)
- [ ] Extract constants to `utils/constants.py`
- [ ] Implement logging in `utils/logging.py`
- [ ] Create validators in `data/validators.py`
- [ ] Add shared I/O utilities
- [ ] Refactor main scripts into modules
- [ ] Add comprehensive README

### Medium-term (Month 1)
- [ ] Implement CLI with Click
- [ ] Create configuration file system
- [ ] Add basic unit tests
- [ ] Write documentation
- [ ] Fix all critical bugs from analysis
- [ ] Add docstrings to all functions

### Long-term (Month 2+)
- [ ] Complete test coverage
- [ ] Add CI/CD workflows
- [ ] Create Jupyter notebook examples
- [ ] Performance optimization
- [ ] Release stable version

---

## Benefits of Restructuring

### For Developers
- **Easier navigation:** Clear separation of concerns
- **Better testing:** Proper module structure enables pytest
- **Reduced duplication:** Shared utilities and constants
- **Faster debugging:** Centralized logging

### For Users
- **Simpler installation:** `pip install -e .`
- **Clearer documentation:** Comprehensive docs in `docs/`
- **Easier customization:** Configuration files
- **Better error messages:** Proper logging and validation

### For Science
- **Reproducibility:** Configuration files track all parameters
- **Transparency:** Well-documented pipeline steps
- **Extensibility:** Modular design allows easy additions
- **Reliability:** Comprehensive testing

---

## Alternative: Minimal Restructuring

If a full restructuring is too much effort, consider these minimal changes:

### Quick Wins (1-2 days)
1. **Add utilities module:**
   - `src/utils.py` with shared functions
   - `src/constants.py` with magic numbers
   
2. **Improve scripts:**
   - Add `set -euo pipefail` to all shell scripts
   - Add error handling and logging
   
3. **Better organization:**
   - `scripts/` directory for operational scripts
   - `slurm/jobs/` for SLURM scripts
   - `slurm/logs/` for outputs
   
4. **Documentation:**
   - Expand README with usage examples
   - Add `docs/` directory with key documentation
   
5. **Configuration:**
   - Create `config/default_config.yaml`
   - Add `.env.example`

This provides 80% of the benefits with 20% of the effort.

---

## Questions to Consider

Before implementing, discuss:

1. **Scope:** Full restructure or minimal changes?
2. **Backward compatibility:** Keep old scripts working?
3. **Timeline:** When should this be completed?
4. **Testing:** Who will test the changes?
5. **Documentation:** Who will write the docs?
6. **Review:** Who will review the changes?

---

## Conclusion

This restructuring will transform the repository from a collection of scripts into a maintainable, documented, and testable software project. While it requires upfront effort, it will:

- Reduce technical debt
- Improve code quality
- Accelerate future development
- Make onboarding new contributors easier
- Increase scientific reproducibility

**Recommended approach:** Start with the minimal restructuring (quick wins), then gradually adopt the full structure as time permits.
