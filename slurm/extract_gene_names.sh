#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy,general # Request partition. Default is 'general'
#SBATCH --qos=short        # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=1:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=2   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=16GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=END     # Set mail type to 'END' to receive a mail when the job finishes.
#SBATCH --output=slurm/out/%j_extract_gene_names.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm/out/%j_extract_gene_names.out # Set name of error log. %j is the Slurm jobId
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working

# Each file produces a corresponding *_gene_names.txt next to the input file.
H5AD_FILES=(
    data/raw/ROSMAP/combined.h5ad
    data/raw/ROSMAP_MIT/combined.h5ad
    data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad
)

for h5ad_file in "${H5AD_FILES[@]}"; do
    out_file="${h5ad_file%.h5ad}_gene_names.txt"
    echo "Extracting gene names: $h5ad_file -> $out_file"
    apptainer exec --writable-tmpfs --pwd /opt/app --containall \
        --bind src/:/opt/app/src/ \
        --bind data/:/opt/app/data/ \
        --env PYTHONPATH=/opt/app/src \
        --env H5AD_FILE="$h5ad_file" \
        --env OUT_FILE="$out_file" \
        ./container_pixi_0-1-3.sif pixi run python -c '
import anndata, os
h5ad_file = os.environ["H5AD_FILE"]
out_file = os.environ["OUT_FILE"]
adata = anndata.read_h5ad(h5ad_file, backed="r")
with open(out_file, "w") as f:
    f.write("\n".join(adata.var_names) + "\n")
print(f"Wrote {len(adata.var_names)} gene names to {out_file}")
'
done
