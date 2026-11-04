#!/bin/bash
# UMAP of a batch-effects spec's sampled cells, one panel per label, from a YAML
# spec in config/batch_effects/ -> data/v2/umap/<spec name>/umap_<gene set>.png
# Needs what the spec needs (data/v2/counts, data/v2/reference, the gene set).
# Extra arguments go to `pipeline umap`, e.g. --annotate for the authors' subclasses.
# Submit from the repo root:
#   sbatch slurm/umap.sh config/batch_effects/rosmap_mit.yaml --annotate cell_type_high_resolution
#   sbatch slurm/umap.sh config/batch_effects/psychad.yaml --annotate subclass
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_umap.out
#SBATCH --error=slurm/out/%j_umap.out

SPEC="${1:?usage: sbatch slurm/umap.sh config/batch_effects/<name>.yaml [pipeline umap options]}"
shift
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Mount config/ (the spec is read in the container), keep the BLAS and numba
# threads within the allocation, and give numba (umap) a writable cache.
EXTRA_BINDS="--bind config/:/opt/app/config/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS --env NUMBA_NUM_THREADS=$THREADS --env NUMBA_CACHE_DIR=/tmp/numba"
source slurm/common.sh

echo "spec: $SPEC $*"
rp pipeline umap "$SPEC" "$@"
