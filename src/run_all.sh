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

# First download all the data (both ROSMAP and ROSMAP_MIT)
echo "[Step 1/8] Downloading ROSMAP data from Synapse..."
python -m rosmap_processing data download
echo "✓ Download complete"
echo ""

# Then list all the files
echo "[Step 2/8] Listing downloaded files..."
echo "Files in data/raw/ROSMAP and data/raw/ROSMAP_MIT:"
ls data/raw/ROSMAP/
ls data/raw/ROSMAP_MIT/
echo ""

# The MIT data is in rds format. We first convert to h5Seurat files, then to h5ad. 
# The other ROSMAP data can be converted from h5Seurat to h5ad directly.
echo "[Step 3/8] Converting RDS files to h5Seurat format..."
pixi run Rscript src/convert_R/rds_to_h5Seurat.R \
    data/raw/ROSMAP_MIT/Astrocytes.rds \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set1.rds \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set2.rds \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set3.rds \
    data/raw/ROSMAP_MIT/Immune_cells.rds \
    data/raw/ROSMAP_MIT/Inhibitory_neurons.rds \
    data/raw/ROSMAP_MIT/OPCs.rds \
    data/raw/ROSMAP_MIT/Oligodendrocytes.rds
echo "✓ RDS to h5Seurat conversion complete"
echo ""

# Next, we can convert all the h5Seurat files to h5ad files.
echo "[Step 4/8] Converting h5Seurat files to h5ad format..."
echo "  - Converting MIT data..."
pixi run Rscript src/convert_R/h5Seurat_to_h5ad.R \
    data/raw/ROSMAP_MIT/Astrocytes.h5Seurat \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5Seurat \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5Seurat \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5Seurat \
    data/raw/ROSMAP_MIT/Immune_cells.h5Seurat \
    data/raw/ROSMAP_MIT/Inhibitory_neurons.h5Seurat \
    data/raw/ROSMAP_MIT/OPCs.h5Seurat \
    data/raw/ROSMAP_MIT/Oligodendrocytes.h5Seurat

# And also do this for the other ROSMAP data
echo "  - Converting ROSMAP data..."
pixi run Rscript src/convert_R/h5Seurat_to_h5ad.R \
    data/raw/ROSMAP/astrocytes.h5Seurat \
    data/raw/ROSMAP/cux2+.h5Seurat \
    data/raw/ROSMAP/cux2-.h5Seurat \
    data/raw/ROSMAP/inhibitory.h5Seurat \
    data/raw/ROSMAP/microglia.h5Seurat \
    data/raw/ROSMAP/oligodendroglia.h5Seurat \
    data/raw/ROSMAP/vascular.niche.h5Seurat
echo "✓ h5Seurat to h5ad conversion complete"
echo ""

# Run fix_categories on all the h5ad files, this turns the cell_type column into categorical
# and removes the raw data layer
echo "[Step 5/8] Fixing categorical columns and removing raw data..."
echo "  - Processing MIT files..."
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Astrocytes.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Immune_cells.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/OPCs.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad --columns cell_type_high_resolution --remove-raw

echo "  - Processing ROSMAP files..."
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/astrocytes.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/cux2+.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/cux2-.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/inhibitory.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/microglia.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/oligodendroglia.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing data category-fix data/raw/ROSMAP/vascular.niche.h5ad --columns cell_type_high_resolution --remove-raw
echo "✓ Category fixing complete"
echo ""

# Combine individual h5ad files into two big files, one for ROSMAP and one for ROSMAP_MIT

# Combine ROSMAP_MIT files
echo "[Step 6/8] Combining h5ad files..."
echo "  - Combining MIT files..."
pixi run python -m rosmap_processing core combine \
    data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad \
    data/raw/ROSMAP_MIT/Immune_cells.h5ad \
    data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad \
    data/raw/ROSMAP_MIT/OPCs.h5ad \
    data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad \
    --output data/raw/ROSMAP_MIT/combined.h5ad

pixi run python -m rosmap_processing utils validate data/raw/ROSMAP_MIT/combined.h5ad

# Combine ROSMAP files
echo "  - Combining ROSMAP files..."
pixi run python -m rosmap_processing core combine \
    data/raw/ROSMAP/astrocytes.h5ad \
    data/raw/ROSMAP/cux2+.h5ad \
    data/raw/ROSMAP/cux2-.h5ad \
    data/raw/ROSMAP/inhibitory.h5ad \
    data/raw/ROSMAP/microglia.h5ad \
    data/raw/ROSMAP/oligodendroglia.h5ad \
    data/raw/ROSMAP/vascular.niche.h5ad \
    --output data/raw/ROSMAP/combined.h5ad

pixi run python -m rosmap_processing utils validate data/raw/ROSMAP/combined.h5ad
echo "✓ Combining complete"
echo ""

# Add metadata to all h5ad files using rosmap_clinical.csv

# MIT files
echo "[Step 7/8] Adding metadata to h5ad files..."
echo "  - Adding metadata to MIT files..."
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Astrocytes.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Immune_cells.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/OPCs.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit

# ROSMAP files
echo "  - Adding metadata to ROSMAP files..."
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/astrocytes.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/cux2+.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/cux2-.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/inhibitory.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/microglia.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/oligodendroglia.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/vascular.niche.h5ad data/raw/ROSMAP/rosmap_clinical.csv

# Combined files
echo "  - Adding metadata to combined files..."
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP/combined.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing data metadata data/raw/ROSMAP_MIT/combined.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit

# Validate combined files
pixi run python -m rosmap_processing utils validate data/raw/ROSMAP/combined.h5ad
pixi run python -m rosmap_processing utils validate data/raw/ROSMAP_MIT/combined.h5ad
echo "✓ Metadata addition complete"
echo ""

# Convert column names to standard format (consistent between SeaAD and ROSMAP)

# MIT files with cell class/subclass
echo "[Step 8/8] Converting column formats..."
echo "  - Converting MIT files..."
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Astrocytes.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass Astrocytes --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Excitatory --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Excitatory --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Excitatory --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Immune_cells.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass Immune --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Inhibitory --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/OPCs.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass OPCs --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass Oligodendrocytes --inplace

# ROSMAP files
echo "  - Converting ROSMAP files..."
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/astrocytes.h5ad --data-type ROSMAP --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/cux2+.h5ad --data-type ROSMAP --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/cux2-.h5ad --data-type ROSMAP --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/inhibitory.h5ad --data-type ROSMAP --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/microglia.h5ad --data-type ROSMAP --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/oligodendroglia.h5ad --data-type ROSMAP --inplace
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/vascular.niche.h5ad --data-type ROSMAP --inplace

# Combined ROSMAP file
pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP/combined.h5ad --data-type ROSMAP --inplace
echo "✓ Column format conversion complete"
echo ""

echo "========================================="
echo "Pipeline complete!"
echo "========================================="

# Note: ROSMAP_MIT combined file already has Class/Celltype from individual files
# If you need to process it, use:
# pixi run python -m rosmap_processing core column-mapping data/raw/ROSMAP_MIT/combined.h5ad --data-type ROSMAP_MIT --cellclass Mixed --subclass Mixed --inplace
