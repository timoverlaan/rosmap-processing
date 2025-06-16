
# First download all the data (both ROSMAP and ROSMAP_MIT)
pixi run python -u src/download_rosmap.py

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

# Run fix_categories.py on all the h5ad files, this turns the cell_type column, which is a string, into a categorical column.
# It also just imports and exports again, fixing a backwards compatibility warning that is automatically solved.
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Astrocytes.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Immune_cells.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/OPCs.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad

pixi run python src/fix_categories.py data/raw/ROSMAP/astrocytes.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP/cux2+.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP/cux2-.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP/inhibitory.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP/microglia.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP/oligodendroglia.h5ad
pixi run python src/fix_categories.py data/raw/ROSMAP/vascular.niche.h5ad

# Next, we join the individual h5ad files into two big h5ad files, one for ROSMAP and one for ROSMAP_MIT.
# And check with the check_h5ad.py script that the combined files are valid.
pixi run python src/combine_h5ad.py \
    data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad \
    data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad \
    data/raw/ROSMAP_MIT/Immune_cells.h5ad \
    data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad \
    data/raw/ROSMAP_MIT/OPCs.h5ad \
    data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad \
    --output data/raw/ROSMAP/combined.h5ad

pixi run python src/check_h5ad.py data/raw/ROSMAP/combined.h5ad

pixi run python src/combine_h5ad.py \
    data/raw/ROSMAP/astrocytes.h5ad \
    data/raw/ROSMAP/cux2+.h5ad \
    data/raw/ROSMAP/cux2-.h5ad \
    data/raw/ROSMAP/inhibitory.h5ad \
    data/raw/ROSMAP/microglia.h5ad \
    data/raw/ROSMAP/oligodendroglia.h5ad \
    data/raw/ROSMAP/vascular.niche.h5ad \
    --output data/raw/ROSMAP_MIT/combined.h5ad

pixi run python src/check_h5ad.py data/raw/ROSMAP_MIT/combined.h5ad

# Next, we add the metadata to all the h5ad files (also the combined ones).
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Astrocytes.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Immune_cells.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/OPCs.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT

pixi run python src/add_metadata.py data/raw/ROSMAP/astrocytes.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP/cux2+.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP/cux2-.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP/inhibitory.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP/microglia.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP/oligodendroglia.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP/vascular.niche.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv

# Finally, we combine the metadata for the combined files.
pixi run python src/add_metadata.py data/raw/ROSMAP/combined.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv
pixi run python src/add_metadata.py data/raw/ROSMAP_MIT/combined.h5ad --metadata data/raw/ROSMAP/ROSMAP_clinical.csv --MIT

# We check the combined files again to make sure everything is still valid.
pixi run python src/check_h5ad.py data/raw/ROSMAP/combined.h5ad
pixi run python src/check_h5ad.py data/raw/ROSMAP_MIT/combined.h5ad

# Finally, we rename all the columns, to make the format consistent between SeaAD and ROSMAP(_MIT).
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Astrocytes.h5ad --inplace --type ROSMAP_MIT --cellclass Glia --subclass Astrocytes
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Excitatory_neurons_set1.h5ad --inplace --type ROSMAP_MIT --cellclass Neuron --subclass Excitatory
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Excitatory_neurons_set2.h5ad --inplace --type ROSMAP_MIT --cellclass Neuron --subclass Excitatory
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Excitatory_neurons_set3.h5ad --inplace --type ROSMAP_MIT --cellclass Neuron --subclass Excitatory
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Immune_cells.h5ad --inplace --type ROSMAP_MIT --cellclass Glia --subclass Immune
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Inhibitory_neurons.h5ad --inplace --type ROSMAP_MIT --cellclass Neuron --subclass Inhibitory
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/OPCs.h5ad --inplace --type ROSMAP_MIT --cellclass Glia --subclass OPCs
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/Oligodendrocytes.h5ad --inplace --type ROSMAP_MIT --cellclass Glia --subclass Oligodendrocytes

pixi run python src/match_columns.py data/raw/ROSMAP/astrocytes.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP/cux2+.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP/cux2-.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP/inhibitory.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP/microglia.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP/oligodendroglia.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP/vascular.niche.h5ad --inplace --type ROSMAP

pixi run python src/match_columns.py data/raw/ROSMAP/combined.h5ad --inplace --type ROSMAP
pixi run python src/match_columns.py data/raw/ROSMAP_MIT/combined.h5ad --inplace --type ROSMAP_MIT
