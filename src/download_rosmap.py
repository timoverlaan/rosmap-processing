#!/usr/bin/env python3
"""
DEPRECATED: This script has been refactored into the rosmap_processing package.

This wrapper maintains backward compatibility with existing workflows.
The functionality has been moved to rosmap_processing.data.download module.

For new code, please use:
    from rosmap_processing.data.download import download_all_data
"""

import sys
import warnings
from pathlib import Path

warnings.warn(
    "src/download_rosmap.py is deprecated and will be removed in v1.0. "
    "Please use: python -m rosmap_processing.data.download",
    DeprecationWarning,
    stacklevel=2
)

# Import the refactored functionality
try:
    from rosmap_processing.data.download import (
        download_rosmap_data,
        download_rosmap_mit_data
    )
    from rosmap_processing.utils.config import SynapseConfig
    from rosmap_processing.utils.logging import setup_logging
except ImportError as e:
    print(f"Error: Could not import refactored module: {e}")
    print("Please ensure the rosmap_processing package is installed.")
    print("Run: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    # Setup basic logging
    logger = setup_logging(level="INFO")
    
    try:
        # Create synapse config
        synapse_config = SynapseConfig(token_file="token.txt")
        
        # Download both datasets (mimicking original behavior)
        output_base = Path("data/raw")
        
        logger.info("Downloading ROSMAP data...")
        output_dir = output_base / "ROSMAP"
        download_rosmap_data(output_dir, synapse_config)
        logger.info("ROSMAP download complete!")
        
        logger.info("Downloading ROSMAP MIT data...")
        output_dir = output_base / "ROSMAP_MIT"
        download_rosmap_mit_data(output_dir, synapse_config)
        logger.info("ROSMAP MIT download complete!")
            
    except Exception as e:
        logger.error(f"Download failed: {e}", exc_info=True)
        sys.exit(1)
