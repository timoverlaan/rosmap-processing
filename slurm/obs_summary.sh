#!/bin/bash
# Step 0 of the batch-effect analysis: the obs columns of the three count files
# (names, types, distinct and missing counts, donor vs cell level; no values)
# -> data/v2/batch_effects/obs_summary_<dataset>.tsv
# Submit from the repo root: sbatch slurm/obs_summary.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_obs_summary.out
#SBATCH --error=slurm/out/%j_obs_summary.out

source slurm/common.sh

DIR="$V2/batch_effects"
skip_if_exists "$DIR/obs_summary_rosmap_mit.tsv" \
    rp utils obs-summary "$MIT_COUNTS" --donor-column projid --output "$DIR/obs_summary_rosmap_mit.tsv"
skip_if_exists "$DIR/obs_summary_rosmap_green.tsv" \
    rp utils obs-summary "$GREEN_COUNTS" --donor-column "$GREEN_DONOR" --output "$DIR/obs_summary_rosmap_green.tsv"
skip_if_exists "$DIR/obs_summary_seaad.tsv" \
    rp utils obs-summary "$SEAAD_H5AD" --donor-column "Donor ID" --output "$DIR/obs_summary_seaad.tsv"
