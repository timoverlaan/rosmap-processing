#!/bin/bash
#SBATCH --account=ewi-insy-prb
#SBATCH --partition=insy,general # Request partition. Default is 'general'
#SBATCH --qos=short        # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=1:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=2   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=16GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=END     # Set mail type to 'END' to receive a mail when the job finishes.
#SBATCH --output=slurm/out/%j_inspect_counts.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm/out/%j_inspect_counts.out # Set name of error log. %j is the Slurm jobId
/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working

H5AD_FILES=(
    data/raw/ROSMAP/combined.h5ad
    data/raw/ROSMAP_MIT/combined.h5ad
)

for h5ad_file in "${H5AD_FILES[@]}"; do
    echo ""
    echo "========================================"
    echo "Inspecting: $h5ad_file"
    echo "========================================"
    apptainer exec --writable-tmpfs --pwd /opt/app --containall \
        --bind src/:/opt/app/src/ \
        --bind data/:/opt/app/data/ \
        --env PYTHONPATH=/opt/app/src \
        --env H5AD_FILE="$h5ad_file" \
        ./container_pixi_0-1-3.sif pixi run python -c '
import anndata as ad
import numpy as np
import os

path = os.environ["H5AD_FILE"]
print(f"Loading (backed=r): {path}")
adata = ad.read_h5ad(path, backed="r")
print(f"Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
print(f"X dtype: {adata.X.dtype}")
print(f"has .raw: {adata.raw is not None}")
print(f"layers: {list(adata.layers.keys()) if adata.layers else []}")
print(f"uns keys: {list(adata.uns.keys()) if adata.uns else []}")
print()

# Sample a slice to inspect values without loading full matrix
n_cells = min(500, adata.shape[0])
n_genes = min(200, adata.shape[1])
sample = adata.X[:n_cells, :n_genes]
if hasattr(sample, "toarray"):
    sample = sample.toarray()
else:
    sample = np.array(sample)

nonzero = sample[sample > 0].flatten()
print(f"--- X matrix sample ({n_cells} cells x {n_genes} genes) ---")
print(f"Non-zero values: {len(nonzero):,} ({100*len(nonzero)/(n_cells*n_genes):.1f}% density)")
if len(nonzero) > 0:
    print(f"Min non-zero:  {nonzero.min():.4f}")
    print(f"Max:           {nonzero.max():.4f}")
    print(f"Mean non-zero: {nonzero.mean():.4f}")
    pct_whole = 100 * np.sum(nonzero == np.floor(nonzero)) / len(nonzero)
    print(f"% whole numbers: {pct_whole:.1f}%")
    if pct_whole > 99:
        print("=> Looks like INTEGER COUNT DATA stored as floats")
    elif nonzero.max() < 20:
        print("=> Looks like LOG-NORMALIZED data (max < 20, non-integer)")
    else:
        print("=> Ambiguous - check values above")
    print(f"Sample values (first 20 non-zeros): {nonzero[:20].tolist()}")
print()

# Inspect each layer the same way
for layer_name in (adata.layers.keys() if adata.layers else []):
    layer = adata.layers[layer_name]
    print(f"--- Layer: {layer_name!r} ---")
    print(f"dtype: {layer.dtype}")
    samp = layer[:n_cells, :n_genes]
    if hasattr(samp, "toarray"):
        samp = samp.toarray()
    nz = samp[samp > 0].flatten()
    if len(nz) > 0:
        print(f"Min non-zero: {nz.min():.4f}, Max: {nz.max():.4f}")
        pct_w = 100 * np.sum(nz == np.floor(nz)) / len(nz)
        print(f"% whole numbers: {pct_w:.1f}%")
        print(f"Sample values (first 20 non-zeros): {nz[:20].tolist()}")
    print()
'
done
