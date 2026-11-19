#!/bin/bash
# Pooled pseudobulk over the five v2 datasets (config/pooled_pseudobulk.yaml)
# -> data/v2/pooled_pseudobulk/: pseudobulk.h5ad (restricted) and aggregate
# tables (metrics, disease, de_pooled, offsets, samples, genes). Needs the
# donor table (slurm/donor_table.sh). Submit from the repo root:
#   sbatch slurm/pooled_pseudobulk.sh              # first run
#   sbatch slurm/pooled_pseudobulk.sh --overwrite  # redo the analyses on the existing pseudobulk
#   sbatch slurm/pooled_pseudobulk.sh --rebuild    # also re-read the count files
# Another spec as the first argument, e.g. RUV-III on the same pseudobulk:
#   sbatch slurm/pooled_pseudobulk.sh config/pooled_ruv.yaml   # -> data/v2/pooled_ruv/
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_pooled_pseudobulk.out
#SBATCH --error=slurm/out/%j_pooled_pseudobulk.out

THREADS="${SLURM_CPUS_PER_TASK:-8}"
EXTRA_BINDS="--bind config/:/opt/app/config/ --env OMP_NUM_THREADS=$THREADS --env OPENBLAS_NUM_THREADS=$THREADS --env MKL_NUM_THREADS=$THREADS"
source slurm/common.sh

SPEC=config/pooled_pseudobulk.yaml
if [[ "${1:-}" == *.yaml ]]; then SPEC="$1"; shift; fi
rp pipeline pooled-pseudobulk "$SPEC" "$@"
