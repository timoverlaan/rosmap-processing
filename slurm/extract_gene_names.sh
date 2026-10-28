#!/bin/bash
# Write the gene names of the combined count files (and SeaAD) to
# *_gene_names.txt next to each file.
# Submit from the repo root: sbatch slurm/extract_gene_names.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --output=slurm/out/%j_extract_gene_names.out
#SBATCH --error=slurm/out/%j_extract_gene_names.out

source slurm/common.sh

for h5ad_file in "$V2/counts/ROSMAP/combined.h5ad" "$V2/counts/ROSMAP_MIT/combined.h5ad" "$SEAAD_H5AD"; do
    out_file="${h5ad_file%.h5ad}_gene_names.txt"
    echo "Extracting gene names: $h5ad_file -> $out_file"
    in_container pixi run python -c '
import sys
from rosmap_processing.data.genes import read_var_names
names = read_var_names(sys.argv[1])
with open(sys.argv[2], "w") as f:
    f.write("\n".join(names) + "\n")
print(f"Wrote {len(names)} gene names to {sys.argv[2]}")
' "$h5ad_file" "$out_file"
done
