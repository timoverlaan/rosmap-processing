#!/bin/bash
# SeaAD: download + ROSMAP-style columns into data/v2/counts/SeaAD/.
# Submit from the repo root: sbatch slurm/run_all_seaad.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_seaad.out
#SBATCH --error=slurm/out/%j_seaad.out

EXTRA_BINDS="--bind output/:/opt/app/output/ --bind config.yaml:/opt/app/config.yaml"
source slurm/common.sh
mkdir -p output

in_container sh ./src/run_all_seaad.sh
