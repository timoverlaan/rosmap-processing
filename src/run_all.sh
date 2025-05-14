
# First download all the data (both ROSMAP and ROSMAP_MIT)
pixi run python -u src/download_synapse.py

# Then list all the files
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
