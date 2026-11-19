#!/bin/bash
# MrVI integration, scoring (CPU): merges samples_A.h5ad and samples_B.h5ad
# into data/v2/scvi/mrvi/samples.h5ad (restricted) and scores the held-out
# persons as the scVI integration (metrics.tsv, disease.tsv, folds.tsv, umap/).
# Needs both folds of slurm/mrvi_train.sh.
#   sbatch --dependency=afterok:<array job> slurm/mrvi_score.sh
#   sbatch slurm/mrvi_score.sh --smoke     # the smoke run's folds
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_mrvi_score.out
#SBATCH --error=slurm/out/%j_mrvi_score.out

CONTAINER="${CONTAINER:-./container_pixi_0-3-0.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CACHES="--env PYTHONPYCACHEPREFIX=/tmp/pycache --env NUMBA_CACHE_DIR=/tmp/numba --env XDG_CACHE_HOME=/tmp/cache --env MPLCONFIGDIR=/tmp/mpl --env TMPDIR=/tmp"
EXTRA_BINDS="--bind config/:/opt/app/config/ $CACHES --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

in_container bash -c 'pixi run -e gpu python -m rosmap_processing pipeline mrvi-integration "$@"; rc=$?
    df -h / /tmp "$HOME" 2>/dev/null | sed "s/^/df: /"; exit $rc' _ config/mrvi_integration.yaml --stage score "$@"
