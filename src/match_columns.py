#!/usr/bin/env python3
"""
DEPRECATED: This script has been refactored into the rosmap_processing package.

This wrapper maintains backward compatibility with existing workflows.
The functionality has been moved to rosmap_processing.core.column_mapping module.

For new code, please use:
    from rosmap_processing.core.column_mapping import convert_columns
"""

import sys
import warnings
from pathlib import Path

warnings.warn(
    "src/match_columns.py is deprecated and will be removed in v1.0. "
    "Please use: python -m rosmap_processing.core.column_mapping",
    DeprecationWarning,
    stacklevel=2
)

# Import the refactored functionality
try:
    from rosmap_processing.core.column_mapping import convert_columns
    from rosmap_processing.utils.logging import setup_logging
    import anndata as ad
except ImportError as e:
    print(f"Error: Could not import refactored module: {e}")
    print("Please ensure the rosmap_processing package is installed.")
    print("Run: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Script converting column names and types between SeaAD and ROSMAP data. "
                    "It also adds the Wang label column that we try to predict"
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to the input h5ad file."
    )
    parser.add_argument(
        "--type",
        type=str,
        required=True,
        help='Type of the input data, either "ROSMAP", "ROSMAP_MIT" or "SeaAD".',
        choices=["ROSMAP", "ROSMAP_MIT", "SeaAD"]
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to the output h5ad file. This will contain the columns in both formats: "
             "ROSMAP format, which is numerically valued, and lowercase. "
             "And the SeaAD format, which contains mainly strings and is capitalized."
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="If set, the input file will be modified in place. "
             "This needs to be set if --output is not provided."
    )
    parser.add_argument(
        "--cellclass",
        type=str,
        help="Cell class for ROSMAP_MIT data (e.g., 'Neuron', 'Glia'). "
             "Required when --type is ROSMAP_MIT."
    )
    parser.add_argument(
        "--subclass",
        type=str,
        help="Cell subclass for ROSMAP_MIT data (e.g., 'Astrocytes', 'Excitatory'). "
             "Required when --type is ROSMAP_MIT."
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
        # Validate output arguments
        if not args.inplace and not args.output:
            raise ValueError("Either --inplace or --output must be provided")
        
        # Load data
        input_path = Path(args.path)
        logger.info(f"Loading data from {input_path}...")
        adata = ad.read_h5ad(input_path)
        
        # Convert columns using refactored function
        adata = convert_columns(
            adata,
            data_type=args.type,
            cellclass=args.cellclass,
            subclass=args.subclass
        )
        
        # Determine output path
        output_path = input_path if args.inplace else Path(args.output)
        
        # Save result
        logger.info(f"Saving converted data to {output_path}...")
        adata.write_h5ad(output_path, compression="gzip")
        logger.info("Done!")
        
    except Exception as e:
        logger.error(f"Failed to convert columns: {e}", exc_info=True)
        sys.exit(1)
