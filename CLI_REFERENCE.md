# ROSMAP Processing CLI Reference

## Overview

The rosmap-processing package now has a unified CLI interface accessible via:
```bash
python -m rosmap_processing <command> <subcommand> [options]
```

## Command Structure

### Data Commands
Manage and prepare data files.

#### Download ROSMAP Data
```bash
python -m rosmap_processing data download
```
Downloads ROSMAP and ROSMAP-MIT data from Synapse.

#### Fix Categorical Columns
```bash
python -m rosmap_processing data category-fix <input.h5ad> [options]
```
Converts specified columns to categorical type and optionally removes raw data layer.

**Options:**
- `--columns COLUMN [COLUMN ...]` - Columns to convert to categorical (default: cell_type_high_resolution)
- `--remove-raw` - Remove raw data layer to save space

**Example:**
```bash
python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    --columns cell_type_high_resolution --remove-raw
```

#### Add Metadata
```bash
python -m rosmap_processing data metadata <input.h5ad> <metadata.csv> [options]
```
Adds metadata from CSV file to h5ad file (modifies in place).

**Options:**
- `--mit` - Use MIT data format (uses different ID column for matching)

**Examples:**
```bash
# For regular ROSMAP data
python -m rosmap_processing data metadata data/raw/ROSMAP/astrocytes.h5ad \
    data/raw/ROSMAP/rosmap_clinical.csv

# For MIT data
python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Immune_cells.h5ad \
    data/raw/ROSMAP/rosmap_clinical.csv --mit
```

---

### Core Processing Commands
Core data processing operations.

#### Combine Files
```bash
python -m rosmap_processing core combine <input1.h5ad> <input2.h5ad> ... --output <output.h5ad>
```
Combines multiple h5ad files into a single file.

**Example:**
```bash
python -m rosmap_processing core combine \
    data/raw/ROSMAP/astrocytes.h5ad \
    data/raw/ROSMAP/cux2+.h5ad \
    data/raw/ROSMAP/microglia.h5ad \
    --output data/raw/ROSMAP/combined.h5ad
```

#### Column Mapping
```bash
python -m rosmap_processing core column-mapping <input.h5ad> --data-type <TYPE> [options]
```
Converts column names between ROSMAP and SeaAD formats.

**Required Options:**
- `--data-type {ROSMAP,ROSMAP_MIT,SeaAD}` - Type of input data

**Optional:**
- `--output <file>` - Output file path (if not using --inplace)
- `--inplace` - Modify file in place
- `--cellclass <class>` - Cell class (required for ROSMAP_MIT)
- `--subclass <subclass>` - Cell subclass (required for ROSMAP_MIT)

**Examples:**
```bash
# Regular ROSMAP data
python -m rosmap_processing core column-mapping data/raw/ROSMAP/astrocytes.h5ad \
    --data-type ROSMAP --inplace

# MIT data (requires cell class/subclass)
python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Immune_cells.h5ad \
    --data-type ROSMAP_MIT --cellclass Glia --subclass Immune --inplace

# SeaAD data
python -m rosmap_processing core column-mapping data/seaAD/sample.h5ad \
    --data-type SeaAD --inplace
```

---

### Pipeline Commands
Complete processing pipelines.

#### Scanpy Pipeline
```bash
python -m rosmap_processing pipeline scanpy <input.h5ad> <output.h5ad> [options]
```
Runs the full scanpy processing pipeline: filtering, normalization, HVG selection, PCA, and KNN graph construction.

**Options:**
- `--layer <name>` - Use specified layer for processing (default: use X)
- `--raw` - Use raw data instead of processed
- `--n-genes <int>` - Number of highly variable genes (default: 2000)
- `--k-neighbors <int>` - Number of neighbors for KNN graph (default: 30)
- `--individual-pca` - Compute PCA separately for each individual
- `--import-genes <file>` - Path to gene list file (.txt or .h5ad) to use instead of HVG selection

**Examples:**
```bash
# Basic usage with defaults (2000 HVGs, 30 neighbors)
python -m rosmap_processing pipeline scanpy \
    data/raw/ROSMAP_MIT/combined.h5ad \
    data/processed/rosmap_processed.h5ad

# With specific parameters
python -m rosmap_processing pipeline scanpy \
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
    data/seaAD/PFC/RNAseq/seaad_2k_k30.h5ad \
    --layer UMIs --n-genes 2000 --k-neighbors 30

# Using imported gene list
python -m rosmap_processing pipeline scanpy \
    data/raw/ROSMAP_MIT/combined.h5ad \
    data/processed/rosmap_mathys_genes.h5ad \
    --n-genes 1000 --k-neighbors 30 \
    --import-genes data/mathys2019_DEGs_genes.txt
```

---

### Utils Commands
Utility functions for inspection and validation.

#### Validate/Inspect File
```bash
python -m rosmap_processing utils validate <input.h5ad>
```
Validates h5ad file and prints detailed information about its contents.

**Example:**
```bash
python -m rosmap_processing utils validate data/raw/ROSMAP/combined.h5ad
```

---

## Full Pipeline Example

Here's how to run the complete processing pipeline:

```bash
# 1. Download data
python -m rosmap_processing data download

# 2. Fix categories and remove raw data
python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    --columns cell_type_high_resolution --remove-raw

# 3. Combine files
python -m rosmap_processing core combine \
    data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    data/raw/ROSMAP_MIT/Immune_cells.h5ad \
    --output data/raw/ROSMAP_MIT/combined.h5ad

# 4. Add metadata
python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/combined.h5ad \
    data/raw/ROSMAP/rosmap_clinical.csv --mit

# 5. Convert column format
python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/combined.h5ad \
    --data-type ROSMAP_MIT --cellclass Mixed --subclass Mixed --inplace

# 6. Run scanpy pipeline
python -m rosmap_processing pipeline scanpy \
    data/raw/ROSMAP_MIT/combined.h5ad \
    data/processed/rosmap_processed.h5ad \
    --n-genes 2000 --k-neighbors 30

# 7. Validate output
python -m rosmap_processing utils validate data/processed/rosmap_processed.h5ad
```

---

## Help

Get help for any command:
```bash
# Main help
python -m rosmap_processing --help

# Command-specific help
python -m rosmap_processing data --help
python -m rosmap_processing core --help
python -m rosmap_processing pipeline --help
python -m rosmap_processing utils --help

# Subcommand-specific help
python -m rosmap_processing data category-fix --help
python -m rosmap_processing pipeline scanpy --help
```

---

## Verbose Logging

Add `-v` or `--verbose` flag for detailed debug logging:
```bash
python -m rosmap_processing -v pipeline scanpy input.h5ad output.h5ad
```

---

## SLURM Usage

The scripts in `slurm/` are already configured to use the new CLI interface. Submit jobs with:

```bash
# Full pipeline
sbatch slurm/run_all.sh

# SeaAD pipeline
sbatch slurm/run_all_seaad.sh

# Specific processing with Mathys genes
sbatch slurm/rosmap_mathys_genes.sh

# Top 1k genes
sbatch slurm/rosmap_top1k_genes.sh

# SeaAD with ROSMAP genes
sbatch slurm/seaad_rosmap_genes.sh
```

All SLURM scripts now use the unified CLI interface internally.
