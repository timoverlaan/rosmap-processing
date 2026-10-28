#!/bin/bash
# Download the two public gene annotations the 10x references were built from,
# for mapping gene symbols to Ensembl IDs (slurm/gene_maps.sh). Light enough
# for the login node; run from the repo root:
#
#   bash slurm/download_references.sh
#
#   Ensembl 93 GTF     -> Cell Ranger GRCh38-3.0.0 (ROSMAP-MIT / Mathys 2023), ~45 MB
#   GENCODE v32 GTF    -> GRCh38-2020-A (SeaAD, Green 2024 ROSMAP), ~40 MB

set -euo pipefail

DEST=data/reference
mkdir -p "$DEST"

fetch() {
    local url="$1" out="$DEST/$(basename "$1")"
    if [ -e "$out" ]; then
        echo "exists, skipping: $out"
        return
    fi
    curl -fL --retry 3 -o "$out.partial" "$url"
    gzip -t "$out.partial"
    mv "$out.partial" "$out"
    echo "downloaded: $out ($(du -h "$out" | cut -f1))"
}

fetch https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
fetch https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
