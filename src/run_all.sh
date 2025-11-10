#!/usr/bin/env bash
# Full pipeline for processing ROSMAP and ROSMAP-MIT data
# Updated to use modular rosmap_processing package

# First download all the data (both ROSMAP and ROSMAP_MIT)
pixi run python -m rosmap_processing.data.download

# Then list all the files
echo "Files in data/raw/ROSMAP and data/raw/ROSMAP_MIT:"
ls data/raw/ROSMAP/
ls data/raw/ROSMAP_MIT/

# The MIT data is in rds format. We first convert to h5Seurat files, then to h5ad. 
# The other ROSMAP data can be converted from h5Seurat to h5ad directly.
pixi run Rscript src/convert_R/rds_to_h5Seurat.R \
    data/raw/ROSMAP_MIT/Astrocytes.rds \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set1.rds \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set2.rds \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set3.rds \
    data/raw/ROSMAP_MIT/Immune_cells.rds \
    data/raw/ROSMAP_MIT/Inhibitory_neurons.rds \
    data/raw/ROSMAP_MIT/OPCs.rds \
    data/raw/ROSMAP_MIT/Oligodendrocytes.rds

# Next, we can convert all the h5Seurat files to h5ad files.
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
pixi run Rscript src/convert_R/h5Seurat_to_h5ad.R \
    data/raw/ROSMAP/astrocytes.h5Seurat \
    data/raw/ROSMAP/cux2+.h5Seurat \
    data/raw/ROSMAP/cux2-.h5Seurat \
    data/raw/ROSMAP/inhibitory.h5Seurat \
    data/raw/ROSMAP/microglia.h5Seurat \
    data/raw/ROSMAP/oligodendroglia.h5Seurat \
    data/raw/ROSMAP/vascular.niche.h5Seurat

# Run fix_categories on all the h5ad files, this turns the cell_type column into categorical
# and removes the raw data layer
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Astrocytes.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Immune_cells.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/OPCs.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad --columns cell_type_high_resolution --remove-raw

pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/astrocytes.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/cux2+.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/cux2-.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/inhibitory.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/microglia.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/oligodendroglia.h5ad --columns cell_type_high_resolution --remove-raw
pixi run python -m rosmap_processing.data.category_fix data/raw/ROSMAP/vascular.niche.h5ad --columns cell_type_high_resolution --remove-raw

# Combine individual h5ad files into two big files, one for ROSMAP and one for ROSMAP_MIT

# Combine ROSMAP_MIT files
pixi run python -m rosmap_processing.core.combine \
    data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad \
    data/raw/ROSMAP_MIT/Immune_cells.h5ad \
    data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad \
    data/raw/ROSMAP_MIT/OPCs.h5ad \
    data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad \
    --output data/raw/ROSMAP_MIT/combined.h5ad

pixi run python -m rosmap_processing.utils.validation data/raw/ROSMAP_MIT/combined.h5ad

# Combine ROSMAP files
pixi run python -m rosmap_processing.core.combine \
    data/raw/ROSMAP/astrocytes.h5ad \
    data/raw/ROSMAP/cux2+.h5ad \
    data/raw/ROSMAP/cux2-.h5ad \
    data/raw/ROSMAP/inhibitory.h5ad \
    data/raw/ROSMAP/microglia.h5ad \
    data/raw/ROSMAP/oligodendroglia.h5ad \
    data/raw/ROSMAP/vascular.niche.h5ad \
    --output data/raw/ROSMAP/combined.h5ad

pixi run python -m rosmap_processing.utils.validation data/raw/ROSMAP/combined.h5ad

# Add metadata to all h5ad files using rosmap_clinical.csv

# MIT files
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Astrocytes.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Immune_cells.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/OPCs.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit

# ROSMAP files
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/astrocytes.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/cux2+.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/cux2-.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/inhibitory.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/microglia.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/oligodendroglia.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/vascular.niche.h5ad data/raw/ROSMAP/rosmap_clinical.csv

# Combined files
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP/combined.h5ad data/raw/ROSMAP/rosmap_clinical.csv
pixi run python -m rosmap_processing.data.metadata data/raw/ROSMAP_MIT/combined.h5ad data/raw/ROSMAP/rosmap_clinical.csv --mit

# Validate combined files
pixi run python -m rosmap_processing.utils.validation data/raw/ROSMAP/combined.h5ad
pixi run python -m rosmap_processing.utils.validation data/raw/ROSMAP_MIT/combined.h5ad

# Convert column names to standard format (consistent between SeaAD and ROSMAP)

# MIT files with cell class/subclass
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Astrocytes.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass Astrocytes
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Excitatory
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Excitatory
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Excitatory
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Immune_cells.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass Immune
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad --data-type ROSMAP_MIT --cellclass Neuron --subclass Inhibitory
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/OPCs.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass OPCs
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad --data-type ROSMAP_MIT --cellclass Glia --subclass Oligodendrocytes

# ROSMAP files
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/astrocytes.h5ad --data-type ROSMAP
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/cux2+.h5ad --data-type ROSMAP
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/cux2-.h5ad --data-type ROSMAP
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/inhibitory.h5ad --data-type ROSMAP
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/microglia.h5ad --data-type ROSMAP
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/oligodendroglia.h5ad --data-type ROSMAP
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/vascular.niche.h5ad --data-type ROSMAP

# Combined ROSMAP file
pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP/combined.h5ad --data-type ROSMAP

# Note: ROSMAP_MIT combined file already has Class/Celltype from individual files
# If you need to process it, use:
# pixi run python -m rosmap_processing.core.column_mapping data/raw/ROSMAP_MIT/combined.h5ad --data-type ROSMAP_MIT --cellclass Mixed --subclass Mixed
