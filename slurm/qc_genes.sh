#!/bin/bash
# Gene ID lists for the data QC: genes on chromosome MT (mitochondrial UMIs),
# X and Y (sex checks), from the Ensembl 91/93 and GENCODE v32 GTFs
# -> data/v2/reference/{mt,x,y}_gene_ids.txt. Submitted by slurm/submit_qc.sh.
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --output=slurm/out/%j_qc_genes.out
#SBATCH --error=slurm/out/%j_qc_genes.out

source slurm/common.sh

GTFS="data/reference/Homo_sapiens.GRCh38.91.gtf.gz data/reference/Homo_sapiens.GRCh38.93.gtf.gz data/reference/gencode.v32.primary_assembly.annotation.gtf.gz"
for chrom in MT X Y; do
    lower="$(echo "$chrom" | tr '[:upper:]' '[:lower:]')"
    # shellcheck disable=SC2086
    skip_if_exists "$REF/${lower}_gene_ids.txt" \
        rp genes on-chromosomes --gtf $GTFS --chromosomes "$chrom" --output "$REF/${lower}_gene_ids.txt"
done
