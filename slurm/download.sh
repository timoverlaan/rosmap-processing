#!/bin/bash
# Download the pinned Synapse files of one or more datasets into
# data/v2/download/<PsychAD|DiverseCohorts|...>/ (MD5-checked; files already
# there are verified, never overwritten). Needs token.txt in the repo root.
# Submit from the repo root:
#   sbatch slurm/download.sh diverse-cohorts psychad
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --output=slurm/out/%j_download.out
#SBATCH --error=slurm/out/%j_download.out

EXTRA_BINDS="--bind token.txt:/opt/app/token.txt"
source slurm/common.sh

[ $# -gt 0 ] || { echo "usage: sbatch slurm/download.sh DATASET [DATASET...]"; exit 1; }
for dataset in "$@"; do
    rp data download --dataset "$dataset" --output-dir "$V2/download"
done
