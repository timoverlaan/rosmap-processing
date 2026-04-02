#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy,general # Request partition. Default is 'general'
#SBATCH --qos=long         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=12:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=4   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=400GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=END     # Set mail type to 'END' to receive a mail when the job finishes.
#SBATCH --output=slurm/out/%j_lieke_hvg.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm/out/%j_lieke_hvg.out # Set name of error log. %j is the Slurm jobId
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working

# Scenario: intersection of Lieke's HVGs x ROSMAP_MIT x SeaAD (~2195 genes)
# Step 1: ROSMAP_MIT — gene selection is the intersection with Lieke's HVG list
# Step 2: SeaAD    — uses the ROSMAP_MIT output as gene reference (exact same genes)

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
	--bind src/:/opt/app/src/ \
	--bind data/:/opt/app/data/ \
	--env PYTHONPATH=/opt/app/src \
	./container_pixi_0-1-3.sif pixi run python -m rosmap_processing pipeline scanpy \
		data/raw/ROSMAP_MIT/combined.h5ad \
		data/processed/rosmap_mit_lieke_hvg_k30.h5ad \
		--import-genes data/lieke_hvg_intersection.txt \
		--k-neighbors 30

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
	--bind src/:/opt/app/src/ \
	--bind data/:/opt/app/data/ \
	--env PYTHONPATH=/opt/app/src \
	./container_pixi_0-1-3.sif pixi run python -m rosmap_processing pipeline scanpy \
		data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad \
		data/seaAD/PFC/RNAseq/seaad_lieke_hvg_k30.h5ad \
		--layer UMIs \
		--import-genes data/lieke_hvg_intersection.txt \
		--k-neighbors 30

# Extract gene names from output files for comparison
extract_genes() {
    apptainer exec --writable-tmpfs --pwd /opt/app --containall \
        --bind data/:/opt/app/data/ \
        --env H5AD_FILE="$1" \
        --env OUT_FILE="${1%.h5ad}_gene_names.txt" \
        ./container_pixi_0-1-3.sif pixi run python -c '
import anndata, os
adata = anndata.read_h5ad(os.environ["H5AD_FILE"], backed="r")
with open(os.environ["OUT_FILE"], "w") as f:
    f.write("\n".join(adata.var_names) + "\n")
print(f"Wrote {len(adata.var_names)} gene names to {os.environ[\"OUT_FILE\"]}")
'
}

extract_genes data/processed/rosmap_mit_lieke_hvg_k30.h5ad
extract_genes data/seaAD/PFC/RNAseq/seaad_lieke_hvg_k30.h5ad
