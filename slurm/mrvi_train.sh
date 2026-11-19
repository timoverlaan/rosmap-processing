#!/bin/bash
# MrVI integration (GPU): per fold, MrVI on the sampled nuclei of the scVI
# integration with half of the persons in two or more datasets as anchors,
# the per-sample decoded expression (data/v2/scvi/mrvi/samples_<fold>.h5ad,
# restricted) and MrVI's sample distances (distances_<fold>.h5, restricted).
# Needs data/v2/scvi/cells.h5ad (slurm/scvi_cells.sh). Then slurm/mrvi_score.sh.
#   sbatch --array=0-1 slurm/mrvi_train.sh [--overwrite]   # task 0: fold A, task 1: fold B
#   sbatch slurm/mrvi_train.sh --fold A                    # one fold
#   sbatch --time=02:00:00 slurm/mrvi_train.sh --stage all --smoke   # a subset, 2 epochs, both folds, scores
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_mrvi_train.out
#SBATCH --error=slurm/out/%j_mrvi_train.out

CONTAINER="${CONTAINER:-./container_pixi_0-3-0.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Python bytecode and caches (numba, CUDA, matplotlib, ...) to the job's /tmp:
# the writable tmpfs holds only 64 MB (see slurm/scvi_train.sh).
CACHES="--env PYTHONPYCACHEPREFIX=/tmp/pycache --env NUMBA_CACHE_DIR=/tmp/numba --env XDG_CACHE_HOME=/tmp/cache --env CUDA_CACHE_PATH=/tmp/nv --env MPLCONFIGDIR=/tmp/mpl --env TMPDIR=/tmp"
EXTRA_BINDS="--nv --bind config/:/opt/app/config/ $CACHES --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

STAGE=train
if [ "${1:-}" = "--stage" ]; then STAGE="$2"; shift 2; fi
FOLD=()
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    FOLDS=(A B)
    FOLD=(--fold "${FOLDS[$SLURM_ARRAY_TASK_ID]}")
    echo "array task $SLURM_ARRAY_TASK_ID: fold ${FOLDS[$SLURM_ARRAY_TASK_ID]}"
fi
in_container bash -c 'pixi run -e gpu python -m rosmap_processing pipeline mrvi-integration "$@"; rc=$?
    df -h / /tmp "$HOME" 2>/dev/null | sed "s/^/df: /"; exit $rc' _ config/mrvi_integration.yaml --stage "$STAGE" \
    ${FOLD[@]+"${FOLD[@]}"} "$@"
