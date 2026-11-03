#!/bin/bash
# PsychAD: verified download -> data/v2/counts/PsychAD/combined.h5ad.
# Needs the download (sbatch slurm/download.sh psychad) and the ROSMAP
# clinical CSV (sbatch slurm/run_all.sh).
# Submit from the repo root: sbatch slurm/run_all_psychad.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_psychad.out
#SBATCH --error=slurm/out/%j_psychad.out

source slurm/common.sh

in_container sh ./src/run_all_psychad.sh
