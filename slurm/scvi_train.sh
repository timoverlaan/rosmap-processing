#!/bin/bash
# scVI integration, stage 2 (GPU): HVGs, scVI and scANVI, per-sample decoded
# expression (data/v2/scvi/samples.h5ad, restricted), then the scores
# (metrics.tsv, disease.tsv, umap/). Needs slurm/scvi_cells.sh first.
#   sbatch slurm/scvi_train.sh [--overwrite]
#   sbatch slurm/scvi_train.sh --stage score     # rescore only (CPU is enough)
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_scvi_train.out
#SBATCH --error=slurm/out/%j_scvi_train.out

CONTAINER="${CONTAINER:-./container_pixi_0-3-0.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Python bytecode and caches (numba, CUDA, matplotlib, ...) to the job's /tmp:
# the writable tmpfs holds only 64 MB, and importing torch, scvi and scanpy
# writes ~90 MB of .pyc files (the gpu env is not precompiled in the image).
CACHES="--env PYTHONPYCACHEPREFIX=/tmp/pycache --env NUMBA_CACHE_DIR=/tmp/numba --env XDG_CACHE_HOME=/tmp/cache --env CUDA_CACHE_PATH=/tmp/nv --env MPLCONFIGDIR=/tmp/mpl --env TMPDIR=/tmp"
EXTRA_BINDS="--nv --bind config/:/opt/app/config/ $CACHES --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

STAGE=train
if [ "${1:-}" = "--stage" ]; then STAGE="$2"; shift 2; fi
in_container bash -c 'pixi run -e gpu python -m rosmap_processing -v pipeline scvi-integration "$@"; rc=$?
    df -h / /tmp "$HOME" 2>/dev/null | sed "s/^/df: /"; exit $rc' _ config/scvi_integration.yaml --stage "$STAGE" "$@"
