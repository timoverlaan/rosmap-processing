#!/bin/bash
# Map ROSMAP-MIT, SeaAD and Green 2024 gene symbols to Ensembl IDs, and list
# the IDs shared by MIT + SeaAD and by all three.
# Needs data/reference/ (bash slurm/download_references.sh) and
# data/v2/counts/{ROSMAP_MIT,ROSMAP}/combined.h5ad (sbatch slurm/run_all.sh).
# Submit from the repo root: sbatch slurm/gene_maps.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --output=slurm/out/%j_gene_maps.out
#SBATCH --error=slurm/out/%j_gene_maps.out

source slurm/common.sh
mkdir -p "$REF"
ENS93=data/reference/Homo_sapiens.GRCh38.93.gtf.gz
GENCODE32=data/reference/gencode.v32.primary_assembly.annotation.gtf.gz

skip_if_exists "$REF/rosmap_mit_gene_map.tsv" \
    rp genes build-map --gtf "$ENS93" --recipe cellranger-GRCh38-3.0.0 \
        --h5ad "$MIT_COUNTS" --output "$REF/rosmap_mit_gene_map.tsv"

skip_if_exists "$REF/seaad_gene_map.tsv" \
    rp genes build-map --gtf "$GENCODE32" --recipe cellranger-GRCh38-2020-A \
        --h5ad "$SEAAD_H5AD" --output "$REF/seaad_gene_map.tsv"

# Green 2024 dropped genes per cell class, so this is name matching against
# the 2020-A list rather than an exact positional match.
skip_if_exists "$REF/rosmap_green_gene_map.tsv" \
    rp genes build-map --gtf "$GENCODE32" --recipe cellranger-GRCh38-2020-A \
        --h5ad "$GREEN_COUNTS" --output "$REF/rosmap_green_gene_map.tsv"

skip_if_exists "$REF/shared_gene_ids.txt" \
    rp genes shared "$REF/rosmap_mit_gene_map.tsv" "$REF/seaad_gene_map.tsv" \
        --output "$REF/shared_gene_ids.txt"

skip_if_exists "$REF/shared_gene_ids_3way.txt" \
    rp genes shared "$REF/rosmap_mit_gene_map.tsv" "$REF/seaad_gene_map.tsv" "$REF/rosmap_green_gene_map.tsv" \
        --output "$REF/shared_gene_ids_3way.txt"
