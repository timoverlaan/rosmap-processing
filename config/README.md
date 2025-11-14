# Configuration System

The ROSMAP processing pipeline uses a YAML configuration file to control all processing parameters.

## Quick Start

1. Copy the default config to create your own:
   ```bash
   cp config/default_config.yaml config.yaml
   ```

2. Edit `config.yaml` to customize parameters:
   - Processing parameters (number of HVGs, k-neighbors, etc.)
   - Paths for data and output
   - Synapse authentication settings
   - Logging preferences

3. Run the pipeline - it will automatically use `config.yaml`:
   ```bash
   sbatch slurm/run_all.sh
   ```

## Configuration Structure

### Processing Parameters
```yaml
processing:
  min_genes: 200              # Minimum genes per cell
  min_cells: 5                # Minimum cells expressing a gene
  n_hvgs: 2000                # Number of highly variable genes
  k_neighbors: 30             # Number of neighbors for KNN graph
  n_pca_components: 50        # Number of PCA components
  target_sum: 1000000         # CPM normalization target
  log_transform: true         # Apply log1p transformation
  individual_pca: false       # Compute PCA per donor
```

### Paths
```yaml
paths:
  raw_data: "data/raw"
  processed: "data/processed"
  interim: "data/interim"
  metadata: "data/metadata"
  output_base: "output"       # Base directory for run outputs
  logs_base: "logs"            # Base directory for run logs (not used by default)
```

Each run creates a timestamped output directory with logs inside:
- `output/rosmap_processing_job_12345_20250111_143022/`
  - `logs/` - Processing logs
  - Data files and results

SLURM's stdout/stderr logs go to:
- `slurm/out/` - SLURM `.out` and `.err` files (separate from run outputs)

### Synapse Authentication
```yaml
synapse:
  use_env_token: true          # Use SYNAPSE_AUTH_TOKEN environment variable
  token_file: "token.txt"      # Fallback to token file
```

### Logging
```yaml
logging:
  level: "INFO"                # DEBUG, INFO, WARNING, ERROR, CRITICAL
  format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
  log_to_file: true
```

## Run Organization

Each processing run automatically creates:

1. **Timestamped output directory**: `output/{run_name}_job_{job_id}_{timestamp}/`
   - Processed data files
   - Intermediate results
   - Analysis outputs
   - `logs/` subdirectory - Processing logs from the pipeline

2. **SLURM logs** (separate location): `slurm/out/`
   - SLURM stdout (`.out` files)
   - SLURM stderr (`.err` files)

3. **Run metadata** (saved to output directory):
   - `config.yaml`: Copy of configuration used for this run
   - `git_info.txt`: Git commit hash, branch, and status
   - Timestamp and job ID information

This ensures:
- **Reproducibility**: Every run is documented with its exact configuration and code version
- **No overwrites**: Each run gets its own directory
- **Easy comparison**: Compare outputs from different runs side-by-side
- **Traceability**: Know exactly which code version and parameters produced each result

## Example Run Structure

After running `sbatch slurm/run_all.sh` with job ID 12345:

```
output/
  rosmap_processing_job_12345_20250111_143022/
    config.yaml              # Configuration used
    git_info.txt             # Git commit and status
    combined_rosmap.h5ad     # Processed outputs
    logs/
      processing.log         # Processing logs
      ...
slurm/
  out/
    12345_all.out            # SLURM stdout
    12345_all.err            # SLURM stderr
```

## Tips

- Keep `config.yaml` out of git (already in `.gitignore`)
- Use `config/default_config.yaml` as a template
- Create task-specific configs: `config_mathys_genes.yaml`, `config_top1k.yaml`, etc.
- Check `git_info.txt` in output directory to see if run had uncommitted changes

