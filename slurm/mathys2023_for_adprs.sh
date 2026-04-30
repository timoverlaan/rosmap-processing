#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy,general
#SBATCH --qos=long
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_mathys2023_for_adprs.out
#SBATCH --error=slurm/out/%j_mathys2023_for_adprs.out
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

# Aggregate Mathys 2023 per-class h5ads -> per-(donor, cell_type, gene) parquet
# for the ad-prs project. Inputs:
#   data/raw/ROSMAP_MIT/{Astrocytes,Excitatory_neurons_set{1,2,3},
#                       Inhibitory_neurons,Immune_cells,Oligodendrocytes,OPCs}.h5ad
#   data/raw/ROSMAP/rosmap_clinical.csv
# Each per-class h5ad needs `projid` and `cell_type_high_resolution` in obs and
# raw counts in X. That's the output of run_all.sh steps 1-4 (download +
# RDS->h5Seurat->h5ad); the later category-fix / metadata / combine /
# column-mapping steps are NOT required by this script. If the per-class h5ads
# don't exist yet, run those steps first (e.g. `sbatch slurm/run_all.sh`).

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
    --bind src/:/opt/app/src/ \
    --bind data/:/opt/app/data/ \
    --env PYTHONPATH=/opt/app/src \
    ./container_pixi_0-2-1.sif pixi run python -m rosmap_processing pipeline mathys2023-for-adprs \
        --input-dir data/raw/ROSMAP_MIT \
        --output-dir data/processed/mathys2023_for_adprs \
        --clinical-csv data/raw/ROSMAP/rosmap_clinical.csv \
        --mit-metadata-csv data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv
