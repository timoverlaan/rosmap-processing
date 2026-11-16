#!/bin/bash
# Global donor table over the five v2 datasets (config/donor_table.yaml)
# -> data/v2/donors/: donor_table.tsv (restricted, one row per person) and
# aggregate overlap.tsv, summary.tsv and conflicts.tsv. Reads obs columns
# only. Submit from the repo root:
#   sbatch slurm/donor_table.sh [--overwrite]
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_donor_table.out
#SBATCH --error=slurm/out/%j_donor_table.out

EXTRA_BINDS="--bind config/:/opt/app/config/"
source slurm/common.sh

in_container pixi run python -m rosmap_processing.pipelines.donor_table config/donor_table.yaml "$@"
