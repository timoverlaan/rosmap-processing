# ROSMAP Processing Pipeline

> Reproducible processing pipeline for ROSMAP Alzheimer's Disease single-cell RNA-seq datasets

## Overview

This pipeline processes single-cell RNA-seq data from the ROSMAP (Religious Orders Study and Memory and Aging Project) and SeaAD datasets, integrating data from multiple sources and standardizing formats for downstream analysis.

## Features

- 🔄 **Automated data download** from Synapse and AWS S3
- 🔬 **Scanpy-based processing** with customizable parameters
- 📊 **Multi-dataset integration** (ROSMAP, ROSMAP-MIT, SeaAD)
- 🐳 **Containerized workflow** using Singularity for reproducibility
- ⚙️ **Configuration-driven** with YAML configs
- 📝 **Comprehensive logging** for debugging and monitoring
- 🧪 **Modular design** for easy extension

## Quick Start

### Prerequisites

- Access to a SLURM HPC cluster
- [Pixi](https://prefix.dev/) for dependency management
- Synapse account with access to ROSMAP data
- Singularity/Apptainer for containerization

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/timoverlaan/rosmap-processing.git
   cd rosmap-processing
   ```

2. **Install dependencies with Pixi:**
   ```bash
   pixi install
   ```

3. **Set up authentication:**
   
   Create a `token.txt` file with your Synapse token:
   ```bash
   echo "your_synapse_token" > token.txt
   ```
   
   Or set an environment variable:
   ```bash
   export SYNAPSE_AUTH_TOKEN="your_synapse_token"
   ```
   
   **⚠️ Important:** Never commit `token.txt` to git!

4. **Build the container (on cluster):**
   ```bash
   # Submit container build job
   sbatch slurm/build_container.sh
   ```

### Basic Usage

The pipeline provides a unified command-line interface:

```bash
# Download data
pixi run rosmap-process download --dataset rosmap

# Process data with default settings
pixi run rosmap-process analyze input.h5ad --output processed.h5ad

# Run full pipeline with custom config
pixi run rosmap-process pipeline --config config/rosmap_config.yaml
```

For SLURM cluster usage:

```bash
# Submit processing job
sbatch slurm/jobs/rosmap_top1k_genes.sh

# Run complete ROSMAP pipeline
sbatch slurm/jobs/run_all.sh

# Run SeaAD pipeline
sbatch slurm/jobs/run_all_seaad.sh
```

## Project Structure

```
rosmap-processing/
├── src/rosmap_processing/    # Main Python package
│   ├── core/                  # Core processing logic
│   ├── data/                  # Data handling utilities
│   ├── utils/                 # Shared utilities
│   └── cli/                   # Command-line interface
├── scripts/                   # Operational scripts
├── slurm/                     # SLURM job scripts
├── config/                    # Configuration files
├── docs/                      # Documentation
├── tests/                     # Test suite
└── data/                      # Data directory (gitignored)
```

## Configuration

Processing parameters can be customized via YAML configuration files:

```yaml
# config/custom_config.yaml
processing:
  min_genes: 200
  n_hvgs: 2000
  k_neighbors: 30
  
paths:
  raw_data: "data/raw"
  processed: "data/processed"
```

See `config/default_config.yaml` for all available options.

## Documentation

- [Installation Guide](docs/installation.md) - Detailed setup instructions
- [Usage Guide](docs/usage.md) - Common usage patterns
- [Pipeline Overview](docs/pipeline.md) - Pipeline architecture
- [Troubleshooting](docs/troubleshooting.md) - Common issues

## Data Sources

This pipeline processes data from:

- **ROSMAP:** Single-nucleus RNA-seq from human prefrontal cortex
  - Source: [Synapse](https://www.synapse.org/#!Synapse:syn21125841)
  - Citation: [Mathys et al., Nature 2019](https://doi.org/10.1038/s41586-019-1195-2)

- **ROSMAP-MIT:** Extended ROSMAP multiomics dataset
  - Source: [Synapse](https://www.synapse.org/#!Synapse:syn52293417)

- **SeaAD:** Seattle Alzheimer's Disease Brain Cell Atlas
  - Source: [AWS S3](https://registry.opendata.aws/allen-sea-ad-atlas/)
  - Citation: [SEA-AD Consortium](https://portal.brain-map.org/explore/seattle-alzheimers-disease)

## Development

### Setting up development environment

```bash
# Install with dev dependencies
pixi install --all

# Run tests
pixi run pytest

# Format code
pixi run black src/ tests/

# Lint code
pixi run ruff src/ tests/
```

### Running tests

```bash
# Run all tests
pixi run pytest

# Run with coverage
pixi run pytest --cov=rosmap_processing

# Run specific test file
pixi run pytest tests/test_specific.py
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{rosmap_processing,
  author = {Verlaan, Timo},
  title = {ROSMAP Processing Pipeline},
  year = {2025},
  url = {https://github.com/timoverlaan/rosmap-processing}
}
```

And the original ROSMAP data sources as appropriate.

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- ROSMAP study participants and investigators
- Rush Alzheimer's Disease Center
- MIT-ROSMAP consortium
- Allen Institute for Brain Science (SEA-AD)

## Contact

**Timo Verlaan**  
Email: t.verlaan@tudelft.nl  
GitHub: [@timoverlaan](https://github.com/timoverlaan)

## Version History

- **v0.2.0** (2025-11-10): Restructured with modular architecture
- **v0.1.2**: Container updates
- **v0.1.1**: Initial container version
- **v0.1.0**: Initial release
