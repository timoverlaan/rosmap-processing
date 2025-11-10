#!/usr/bin/env python3
"""
DEPRECATED: This script has been refactored into the rosmap_processing package.

This wrapper maintains backward compatibility with existing workflows.
The functionality has been moved to rosmap_processing.data.metadata module.

For new code, please use:
    from rosmap_processing.data.metadata import add_metadata_to_file
"""

import sys
import warnings
from pathlib import Path

warnings.warn(
    "src/add_metadata.py is deprecated and will be removed in v1.0. "
    "Please use: python -m rosmap_processing.data.metadata",
    DeprecationWarning,
    stacklevel=2
)

# Import the refactored functionality
try:
    from rosmap_processing.data.metadata import add_metadata_to_file
    from rosmap_processing.utils.logging import setup_logging
except ImportError as e:
    print(f"Error: Could not import refactored module: {e}")
    print("Please ensure the rosmap_processing package is installed.")
    print("Run: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Script to add metadata to an AnnData object. This should only be necessary for the ROSMAP data, as SeaAD already includes it."
    )
    parser.add_argument("path", type=str, 
                       help="Path to the input h5ad file. Note: this will be overwritten!")
    parser.add_argument("--metadata", type=str, required=True,
                       help="Path to the metadata CSV file. The CSV should have a 'cell_id' column that matches the index of the AnnData object.")
    parser.add_argument("--MIT", action="store_true",
                       help="If set, the input data is from MIT. This will use 'individualID' for matching.")
    parser.add_argument("--log-level", type=str, default="INFO",
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       help="Logging level")
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(level=args.log_level)
    
    try:
        # Call the refactored function
        # Output to same file to maintain backward compatibility
        add_metadata_to_file(
            h5ad_path=Path(args.path),
            metadata_path=Path(args.metadata),
            is_mit=args.MIT,
            output_path=Path(args.path),  # Overwrite original file as per original behavior
            validate=True
        )
        logger.info(f"Metadata added to {args.path}.")
        
    except Exception as e:
        logger.error(f"Failed to add metadata: {e}", exc_info=True)
        sys.exit(1)
