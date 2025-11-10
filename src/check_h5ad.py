#!/usr/bin/env python3
"""
DEPRECATED: This script has been refactored into the rosmap_processing package.

This wrapper maintains backward compatibility with existing workflows.
The functionality has been moved to rosmap_processing.utils.validation module.

For new code, please use:
    from rosmap_processing.utils.validation import check_h5ad_file
"""

import sys
import warnings
from pathlib import Path

warnings.warn(
    "src/check_h5ad.py is deprecated and will be removed in v1.0. "
    "Please use: python -m rosmap_processing.utils.validation",
    DeprecationWarning,
    stacklevel=2
)

# Import the refactored functionality
try:
    from rosmap_processing.utils.validation import check_h5ad_file
    from rosmap_processing.utils.logging import setup_logging
except ImportError as e:
    print(f"Error: Could not import refactored module: {e}")
    print("Please ensure the rosmap_processing package is installed.")
    print("Run: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Check a h5ad file.")
    parser.add_argument(
        "path",
        type=str,
        help="Path to the h5ad file to check."
    )
    
    args = parser.parse_args()
    
    # Setup logging (match original output style with INFO level)
    setup_logging(level="INFO")
    
    try:
        # Call the refactored function
        check_h5ad_file(
            file_path=Path(args.path),
            show_metadata=True,
            show_matrix=True,
            preview_size=20
        )
        
    except Exception as e:
        print(f"Error: Failed to check file: {e}", file=sys.stderr)
        sys.exit(1)
