#!/usr/bin/env python3
"""
DEPRECATED: This script has been refactored into the rosmap_processing package.

This wrapper maintains backward compatibility with existing workflows.
The functionality has been moved to rosmap_processing.data.category_fix module.

For new code, please use:
    from rosmap_processing.data.category_fix import fix_categories_in_file
"""

import sys
import warnings
from pathlib import Path

warnings.warn(
    "src/fix_categories.py is deprecated and will be removed in v1.0. "
    "Please use: python -m rosmap_processing.data.category_fix",
    DeprecationWarning,
    stacklevel=2
)

# Import the refactored functionality
try:
    from rosmap_processing.data.category_fix import fix_categories_in_file
    from rosmap_processing.utils.logging import setup_logging
except ImportError as e:
    print(f"Error: Could not import refactored module: {e}")
    print("Please ensure the rosmap_processing package is installed.")
    print("Run: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    # Setup basic logging
    setup_logging(level="INFO")
    
    # Get path from command line (original behavior)
    if len(sys.argv) != 2:
        print("Usage: python fix_categories.py <path_to_anndata>")
        sys.exit(1)
    
    path = sys.argv[1]
    
    try:
        # Call the refactored function
        # Original behavior: overwrites file in place, removes raw data
        fix_categories_in_file(
            input_path=Path(path),
            output_path=None,  # Overwrite in place
            columns=None,  # Use default columns
            remove_raw=True
        )
        
    except Exception as e:
        print(f"Error: Failed to fix categories: {e}", file=sys.stderr)
        sys.exit(1)
