#!/bin/bash
# Which genes make up the gap between our per-cell UMIs and the authors'
# (scripts/totals_gap.py; needs the spec's qc-scan). Aggregate output only.
#   sbatch slurm/totals_gap.sh config/data_qc/diverse_cohorts.yaml
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --output=slurm/out/%j_totals_gap.out
#SBATCH --error=slurm/out/%j_totals_gap.out

SPEC="${1:?usage: sbatch slurm/totals_gap.sh config/data_qc/<name>.yaml}"
EXTRA_BINDS="--bind config/:/opt/app/config/ --bind scripts/:/opt/app/scripts/"
source slurm/common.sh

in_container pixi run python scripts/totals_gap.py "$SPEC"
