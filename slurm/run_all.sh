#!/bin/bash
# ROSMAP + ROSMAP-MIT: download -> raw-count h5ad files in data/v2/counts/.
# Submit from the repo root: sbatch slurm/run_all.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_run_all.out
#SBATCH --error=slurm/out/%j_run_all.out

EXTRA_BINDS="--bind output/:/opt/app/output/ --bind config.yaml:/opt/app/config.yaml --bind token.txt:/opt/app/token.txt"
source slurm/common.sh
mkdir -p output

in_container sh ./src/run_all.sh
