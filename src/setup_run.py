#!/usr/bin/env python3
"""
Setup script for initializing a processing run.
Creates output directories, saves config and git info.
"""

import sys
import os
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from rosmap_processing.utils.config import load_config

def main():
    """Set up a processing run."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Set up processing run directories")
    parser.add_argument("--config", type=Path, default="config.yaml",
                       help="Path to config file")
    parser.add_argument("--run-name", required=True,
                       help="Name of the run (e.g., rosmap_processing)")
    parser.add_argument("--job-id", 
                       help="SLURM job ID (automatically set from environment)")
    
    args = parser.parse_args()
    
    # Get job ID from environment if not provided
    job_id = args.job_id or os.environ.get("SLURM_JOB_ID")
    
    # Load config
    config = load_config(args.config)
    
    # Set up run directories
    output_dir, log_dir = config.setup_run(
        run_name=args.run_name,
        job_id=job_id,
        config_source=args.config
    )
    
    # Print paths for shell scripts to use
    print(f"OUTPUT_DIR={output_dir}")
    print(f"LOG_DIR={log_dir}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
