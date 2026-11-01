#!/bin/bash
# Batch effects per gene set, from a YAML spec in config/batch_effects/
# (aggregate tables only) -> the spec's output_dir, e.g. data/v2/batch_effects/rosmap_mit/
# Needs data/v2/counts, data/v2/reference, data/v2/processed and data/v2/hvg_compare;
# check the spec's column names against slurm/obs_summary.sh first.
# Submit from the repo root:
#   sbatch slurm/batch_effects.sh config/batch_effects/rosmap_mit.yaml
#   sbatch slurm/batch_effects.sh config/batch_effects/rosmap_green.yaml
#   sbatch slurm/batch_effects.sh config/batch_effects/seaad.yaml
#   sbatch slurm/batch_effects.sh config/batch_effects/mit_vs_seaad.yaml
#   sbatch slurm/batch_effects.sh config/batch_effects/mit_vs_green.yaml
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_batch_effects.out
#SBATCH --error=slurm/out/%j_batch_effects.out

SPEC="${1:?usage: sbatch slurm/batch_effects.sh config/batch_effects/<name>.yaml}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Mount config/ (the spec is read in the container) and keep the BLAS threads
# within the allocation.
EXTRA_BINDS="--bind config/:/opt/app/config/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

echo "spec: $SPEC"
rp pipeline batch-effects "$SPEC"
