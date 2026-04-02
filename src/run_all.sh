#!/usr/bin/env bash
# Full pipeline for processing ROSMAP and ROSMAP-MIT data
# Uses unified CLI interface and config system

# Exit on any error
set -e

echo "========================================="
echo "ROSMAP Processing Pipeline"
echo "========================================="
echo ""

# Set up run directories and save config/git info
echo "Setting up run directories..."
RUN_INFO=$(pixi run python src/setup_run.py --config config.yaml --run-name rosmap_processing)
eval "$RUN_INFO"
echo "Output directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo ""

# MIT files
echo "[Step 7/8] Adding metadata to h5ad files..."
echo "  - Adding metadata to combined files..."
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/combined.h5ad data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/combined.h5ad data/raw/ROSMAP/ROSMAP_clinical.csv --mit

# Validate combined files
pixi run python -m rosmap_processing utils validate data/raw/ROSMAP/combined.h5ad
pixi run python -m rosmap_processing utils validate data/raw/ROSMAP_MIT/combined.h5ad
echo "✓ Metadata addition complete"
echo ""

echo "========================================="
echo "Pipeline complete!"
echo "========================================="

# Note: ROSMAP_MIT combined file already has Class/Celltype from individual files
# If you need to process it, use:
# pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/combined.h5ad --data-type ROSMAP_MIT --cellclass Mixed --subclass Mixed --inplace
