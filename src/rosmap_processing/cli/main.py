"""Main CLI entry point for rosmap-processing."""

import sys
import argparse
from pathlib import Path

# Ensure the package can be imported
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rosmap_processing.utils.logging import setup_logging, get_logger


def main():
    """Main entry point for CLI."""
    parser = argparse.ArgumentParser(
        prog="rosmap-processing",
        description="ROSMAP and SeaAD data processing pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to config file (optional, uses defaults if not provided)"
    )
    
    subparsers = parser.add_subparsers(
        title="commands",
        description="Available processing commands",
        dest="command",
        required=True,
    )
    
    # Data commands
    data_parser = subparsers.add_parser("data", help="Data management commands")
    data_subparsers = data_parser.add_subparsers(dest="data_command", required=True)
    
    # data download
    data_subparsers.add_parser("download", help="Download ROSMAP data from Synapse")
    
    # data category-fix
    category_parser = data_subparsers.add_parser("category-fix", help="Fix categorical columns and remove raw data")
    category_parser.add_argument("input", type=str, help="Path to input h5ad file")
    category_parser.add_argument("--columns", type=str, nargs="+", default=["cell_type_high_resolution"],
                                 help="Columns to convert to categorical (default: cell_type_high_resolution)")
    category_parser.add_argument("--remove-raw", action="store_true", help="Remove raw data layer")
    
    # data metadata
    metadata_parser = data_subparsers.add_parser("metadata", help="Add metadata to h5ad file")
    metadata_parser.add_argument("input", type=str, help="Path to input h5ad file (will be modified in place)")
    metadata_parser.add_argument("metadata", type=str, help="Path to metadata CSV file")
    metadata_parser.add_argument("--mit", action="store_true", help="Use MIT data format (different ID column)")
    
    # Core processing commands
    core_parser = subparsers.add_parser("core", help="Core processing commands")
    core_subparsers = core_parser.add_subparsers(dest="core_command", required=True)
    
    # core combine
    combine_parser = core_subparsers.add_parser("combine", help="Combine multiple h5ad files")
    combine_parser.add_argument("inputs", type=str, nargs="+", help="Input h5ad files to combine")
    combine_parser.add_argument("--output", type=str, required=True, help="Output h5ad file path")
    
    # core column-mapping
    mapping_parser = core_subparsers.add_parser("column-mapping", help="Convert column names between ROSMAP and SeaAD formats")
    mapping_parser.add_argument("input", type=str, help="Path to input h5ad file")
    mapping_parser.add_argument("--data-type", type=str, required=True, choices=["ROSMAP", "ROSMAP_MIT", "SeaAD"],
                                help="Type of input data")
    mapping_parser.add_argument("--output", type=str, help="Path to output h5ad file (optional if --inplace)")
    mapping_parser.add_argument("--inplace", action="store_true", help="Modify file in place")
    mapping_parser.add_argument("--cellclass", type=str, help="Cell class (required for ROSMAP_MIT)")
    mapping_parser.add_argument("--subclass", type=str, help="Cell subclass (required for ROSMAP_MIT)")
    
    # Pipeline commands
    pipeline_parser = subparsers.add_parser("pipeline", help="Processing pipeline commands")
    pipeline_subparsers = pipeline_parser.add_subparsers(dest="pipeline_command", required=True)
    
    # pipeline scanpy
    scanpy_parser = pipeline_subparsers.add_parser("scanpy", help="Run scanpy processing pipeline")
    scanpy_parser.add_argument("input", type=str, help="Path to input h5ad file")
    scanpy_parser.add_argument("output", type=str, help="Path to output h5ad file")
    scanpy_parser.add_argument("--layer", type=str, default=None, help="Layer to use for processing")
    scanpy_parser.add_argument("--raw", action="store_true", help="Use raw data instead of processed")
    scanpy_parser.add_argument("--n-genes", type=int, default=2000, help="Number of highly variable genes (default: 2000)")
    scanpy_parser.add_argument("--k-neighbors", type=int, default=30, help="Number of neighbors for KNN graph (default: 30)")
    scanpy_parser.add_argument("--individual-pca", action="store_true", help="Compute PCA per individual")
    scanpy_parser.add_argument("--import-genes", type=str, default=None, help="Path to gene list file (.txt or .h5ad)")
    scanpy_parser.add_argument("--hvg-after-import", type=int, default=None, help="If set alongside --import-genes, further select N HVGs within the imported gene subset (pre-normalization, seurat_v3 flavor)")
    scanpy_parser.add_argument("--hvg-after-normalize", type=int, default=None, help="If set alongside --import-genes, further select N HVGs after normalization (seurat flavor). Mutually exclusive with --hvg-after-import.")

    # pipeline mathys2023-for-adprs
    adprs_parser = pipeline_subparsers.add_parser(
        "mathys2023-for-adprs",
        help="Aggregate Mathys 2023 per-class h5ads into per-(donor, cell_type, gene) stats for ad-prs",
    )
    adprs_parser.add_argument(
        "--input-dir", type=Path, default=Path("data/raw/ROSMAP_MIT"),
        help="Directory containing per-class h5ad files (default: data/raw/ROSMAP_MIT)",
    )
    adprs_parser.add_argument(
        "--output-dir", type=Path, default=Path("data/processed/mathys2023_for_adprs"),
        help="Directory to write outputs into (default: data/processed/mathys2023_for_adprs)",
    )
    adprs_parser.add_argument(
        "--clinical-csv", type=Path, default=Path("data/raw/ROSMAP/rosmap_clinical.csv"),
        help="ROSMAP clinical CSV (default: data/raw/ROSMAP/rosmap_clinical.csv)",
    )
    adprs_parser.add_argument(
        "--mit-metadata-csv", type=Path,
        default=Path("data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv"),
        help="MIT individual metadata CSV (for individualID -> subject mapping)",
    )
    adprs_parser.add_argument(
        "--keep-intermediate", action="store_true",
        help="Keep per-class parquet chunks under _intermediate/ after merging",
    )

    # Utils commands
    utils_parser = subparsers.add_parser("utils", help="Utility commands")
    utils_subparsers = utils_parser.add_subparsers(dest="utils_command", required=True)
    
    # utils validate
    validate_parser = utils_subparsers.add_parser("validate", help="Validate and inspect h5ad file")
    validate_parser.add_argument("input", type=str, help="Path to h5ad file to validate")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Setup logging
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(level=log_level)
    logger = get_logger(__name__)
    
    # Route to appropriate command
    try:
        if args.command == "data":
            handle_data_command(args)
        elif args.command == "core":
            handle_core_command(args)
        elif args.command == "pipeline":
            handle_pipeline_command(args)
        elif args.command == "utils":
            handle_utils_command(args)
        else:
            parser.print_help()
            sys.exit(1)
    except Exception as e:
        logger.error(f"Command failed: {e}", exc_info=args.verbose)
        sys.exit(1)


def handle_data_command(args):
    """Handle data subcommands."""
    logger = get_logger(__name__)
    
    if args.data_command == "download":
        from rosmap_processing.data.download import download_rosmap_data
        logger.info("Downloading ROSMAP data...")
        download_rosmap_data()
        
    elif args.data_command == "category-fix":
        from rosmap_processing.data.category_fix import fix_categories_in_file
        from pathlib import Path
        logger.info(f"Fixing categories in {args.input}")
        fix_categories_in_file(
            input_path=Path(args.input),
            columns=args.columns,
            remove_raw=args.remove_raw
        )
        
    elif args.data_command == "metadata":
        from rosmap_processing.data.metadata import add_metadata_to_file
        logger.info(f"Adding metadata to {args.input}")
        add_metadata_to_file(
            h5ad_path=args.input,
            metadata_path=args.metadata,
            is_mit=args.mit
        )


def handle_core_command(args):
    """Handle core subcommands."""
    logger = get_logger(__name__)
    
    if args.core_command == "combine":
        from rosmap_processing.core.combine import combine_and_save
        from pathlib import Path
        logger.info(f"Combining {len(args.inputs)} files...")
        combine_and_save(
            input_paths=[Path(p) for p in args.inputs],
            output_path=Path(args.output)
        )
        
    elif args.core_command == "column-mapping":
        from rosmap_processing.core.column_mapping import convert_column_format
        logger.info(f"Converting column format for {args.input}")
        
        if args.data_type == "ROSMAP_MIT" and (not args.cellclass or not args.subclass):
            raise ValueError("--cellclass and --subclass are required for ROSMAP_MIT data")
        
        if not args.inplace and not args.output:
            raise ValueError("Either --inplace or --output must be specified")
        
        convert_column_format(
            input_path=args.input,
            data_type=args.data_type,
            output_path=args.output if not args.inplace else None,
            cellclass=args.cellclass if args.data_type == "ROSMAP_MIT" else None,
            subclass=args.subclass if args.data_type == "ROSMAP_MIT" else None
        )


def handle_pipeline_command(args):
    """Handle pipeline subcommands."""
    logger = get_logger(__name__)
    
    if args.pipeline_command == "scanpy":
        from rosmap_processing.core.scanpy_pipeline import process_h5ad
        logger.info(f"Running scanpy pipeline on {args.input}")
        process_h5ad(
            input_path=args.input,
            output_path=args.output,
            layer=args.layer,
            use_raw=args.raw,
            n_hvgs=args.n_genes,
            k_neighbors=args.k_neighbors,
            individual_pca=args.individual_pca,
            import_genes=args.import_genes,
            hvg_after_import=args.hvg_after_import,
            hvg_after_normalize=args.hvg_after_normalize,
        )

    elif args.pipeline_command == "mathys2023-for-adprs":
        from rosmap_processing.pipelines.mathys2023_for_adprs import run_pipeline
        logger.info(f"Running Mathys 2023 ad-prs aggregation; input={args.input_dir}")
        run_pipeline(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            clinical_csv=args.clinical_csv,
            mit_metadata_csv=args.mit_metadata_csv,
            keep_intermediate=args.keep_intermediate,
        )


def handle_utils_command(args):
    """Handle utils subcommands."""
    logger = get_logger(__name__)
    
    if args.utils_command == "validate":
        from rosmap_processing.utils.validation import inspect_h5ad
        logger.info(f"Validating {args.input}")
        inspect_h5ad(args.input)


if __name__ == "__main__":
    main()
