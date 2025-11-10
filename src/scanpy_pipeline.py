"""
DEPRECATED: This script has been refactored.

Use the new modular version instead:
    python -m rosmap_processing.core.scanpy_pipeline

Or for backward compatibility, this wrapper will redirect to the new version.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import and use the new refactored version
from rosmap_processing.core.scanpy_pipeline import process_h5ad
from rosmap_processing.utils.logging import setup_logging
from rosmap_processing.utils.constants import DEFAULT_HVG_COUNT, DEFAULT_K_NEIGHBORS

if __name__ == "__main__":
    import argparse
    
    # Setup logging first
    setup_logging(level="INFO")
    
    parser = argparse.ArgumentParser(description="Performs the basic scanpy processing steps on a given h5ad file.")
    parser.add_argument("path", type=str, help="Path to the h5ad file to process.")
    parser.add_argument("output", type=str, help="Path to save the processed h5ad file.")
    parser.add_argument("--layer", type=str, default=None, help="Layer to select for processing. If None, the main data will be used.")
    parser.add_argument("--raw", action="store_true", help="If set, the raw data will be used instead of the main data.")
    parser.add_argument("--n-genes", type=int, default=DEFAULT_HVG_COUNT, help=f"Number of highly variable genes to select. Default is {DEFAULT_HVG_COUNT}.")
    parser.add_argument("--k-neighbors", type=int, default=DEFAULT_K_NEIGHBORS, help=f"Number of neighbors for the KNN graph. Default is {DEFAULT_K_NEIGHBORS}.")
    parser.add_argument("--individual-pca", action="store_true", help="If set, PCA for the KNN graph will be computed for each individual separately.")
    parser.add_argument("--import-genes", type=str, default=None, help="Path to a TXT or AnnData file to use the genes from. This overrides HVG selection.")
    
    args = parser.parse_args()
    
    try:
        process_h5ad(
            input_path=args.path,
            output_path=args.output,
            layer=args.layer,
            use_raw=args.raw,
            n_hvgs=args.n_genes,
            k_neighbors=args.k_neighbors,
            individual_pca=args.individual_pca,
            import_genes=args.import_genes,
        )
    except Exception as e:
        from rosmap_processing.utils.logging import get_logger
        logger = get_logger(__name__)
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)

