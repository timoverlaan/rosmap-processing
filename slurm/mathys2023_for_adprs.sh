#!/bin/bash
# Aggregate the Mathys 2023 per-class count files -> per-(donor, cell_type,
# gene) parquet for the ad-prs project, with Ensembl gene IDs.
# Needs data/v2/counts/ROSMAP_MIT (slurm/run_all.sh) and
# data/v2/reference/rosmap_mit_gene_map.tsv (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/mathys2023_for_adprs.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_mathys2023_for_adprs.out
#SBATCH --error=slurm/out/%j_mathys2023_for_adprs.out

source slurm/common.sh

ADPRS_OUT="$OUT/mathys2023_for_adprs"
skip_if_exists "$ADPRS_OUT/README.md" \
    rp pipeline mathys2023-for-adprs \
        --input-dir "$V2/counts/ROSMAP_MIT" \
        --output-dir "$ADPRS_OUT" \
        --gene-map "$REF/rosmap_mit_gene_map.tsv" \
        --clinical-csv data/v2/download/ROSMAP/ROSMAP_clinical.csv \
        --mit-metadata-csv data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv
