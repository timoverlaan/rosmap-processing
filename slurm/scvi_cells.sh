#!/bin/bash
# scVI integration, stage 1 (CPU): up to max_cells_per_sample nuclei per
# dataset x donor x cell type -> data/v2/scvi/cells.h5ad (restricted).
# Needs the 0.3.0 container. Then slurm/scvi_train.sh.
#   sbatch slurm/scvi_cells.sh [--overwrite]
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_scvi_cells.out
#SBATCH --error=slurm/out/%j_scvi_cells.out

CONTAINER="${CONTAINER:-./container_pixi_0-3-0.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
EXTRA_BINDS="--bind config/:/opt/app/config/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

in_container pixi run -e gpu python -m rosmap_processing pipeline scvi-integration config/scvi_integration.yaml \
    --stage cells "$@"
