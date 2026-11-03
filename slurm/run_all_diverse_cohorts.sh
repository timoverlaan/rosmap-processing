#!/bin/bash
# DiverseCohorts: verified download -> DLPFC count files in data/v2/counts/DiverseCohorts/.
# Needs the download (sbatch slurm/download.sh diverse-cohorts) and the ROSMAP
# clinical CSV (sbatch slurm/run_all.sh).
# Submit from the repo root: sbatch slurm/run_all_diverse_cohorts.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_diverse_cohorts.out
#SBATCH --error=slurm/out/%j_diverse_cohorts.out

source slurm/common.sh

in_container sh ./src/run_all_diverse_cohorts.sh
