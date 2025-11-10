#!/usr/bin/env python3
"""
DEPRECATED: This script has been refactored into the rosmap_processing package.

This wrapper maintains backward compatibility with existing workflows.
The functionality has been moved to rosmap_processing.core.combine module.

For new code, please use:
    from rosmap_processing.core.combine import combine_and_save
"""

import sys
import warnings
from pathlib import Path

warnings.warn(
    "src/combine_h5ad.py is deprecated and will be removed in v1.0. "
    "Please use: python -m rosmap_processing.core.combine",
    DeprecationWarning,
    stacklevel=2
)

# Import the refactored functionality
try:
    from rosmap_processing.core.combine import combine_and_save
    from rosmap_processing.utils.logging import setup_logging
except ImportError as e:
    print(f"Error: Could not import refactored module: {e}")
    print("Please ensure the rosmap_processing package is installed.")
    print("Run: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Combine all provided h5ad files into a single h5ad file."
    )
    parser.add_argument(
        "paths",
        type=str,
        nargs="+",
        help="Paths to the input h5ad files."
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output h5ad file. This will contain all the data from the input files."
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(level=args.log_level)
    
    try:
        # Convert paths to Path objects
        input_paths = [Path(p) for p in args.paths]
        output_path = Path(args.output)
        
        # Call the refactored function
        combine_and_save(
            input_paths=input_paths,
            output_path=output_path
        )
        
    except Exception as e:
        logger.error(f"Failed to combine files: {e}", exc_info=True)
        sys.exit(1)
