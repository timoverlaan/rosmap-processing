#!/bin/bash
# Compare HVG selections: seurat_v3 (counts) vs seurat (log-normalized), per
# donor vs pooled, 1k/2k/3k genes.
#   ROSMAP-MIT, within the genes shared with SeaAD, also against the April 2026
#   selections (made before the count/normalization fixes, and not restricted
#   to shared genes) and Lieke's HVG list -> data/v2/hvg_compare/rosmap_mit/
#   Green 2024, on its own genes -> data/v2/hvg_compare/rosmap_green/
# Gene lists and overlap tables only (gene names, no cell or donor data).
# Needs data/v2/counts (slurm/run_all.sh) and data/v2/reference (slurm/gene_maps.sh).
# Submit from the repo root: sbatch slurm/hvg_compare.sh
#SBATCH --account=daic-eemcs-insy-prb
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --output=slurm/out/%j_hvg_compare.out
#SBATCH --error=slurm/out/%j_hvg_compare.out

source slurm/common.sh

OLD=data/processed
COMPARE=(lieke_hvg=data/lieke_hvg_intersection.txt)
for N in 1000 2000 3000; do
    COMPARE+=("april_seurat_v3_pooled_${N}=$OLD/rosmap_mit_hvg${N}_k30.h5ad")
    COMPARE+=("april_seurat_pooled_${N}=$OLD/rosmap_mit_hvg${N}_postnorm_k30.h5ad")
done

skip_if_exists "$V2/hvg_compare/rosmap_mit/hvg_overlap_summary.tsv" \
    rp pipeline hvg-compare "$MIT_COUNTS" \
        --output-dir "$V2/hvg_compare/rosmap_mit" \
        --gene-map "$REF/rosmap_mit_gene_map.tsv" \
        --restrict-genes "$REF/shared_gene_ids.txt" \
        --compare "${COMPARE[@]}"

# Green 2024 on its own genes (as in green_hvg.sh)
skip_if_exists "$V2/hvg_compare/rosmap_green/hvg_overlap_summary.tsv" \
    rp pipeline hvg-compare "$GREEN_COUNTS" \
        --output-dir "$V2/hvg_compare/rosmap_green" \
        --gene-map "$REF/rosmap_green_gene_map.tsv" \
        --donor-column "$GREEN_DONOR"
