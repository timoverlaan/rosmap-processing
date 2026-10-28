"""Main CLI entry point for rosmap-processing."""

import sys
import argparse
from pathlib import Path

# Ensure the package can be imported
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rosmap_processing.utils.constants import (
    DEFAULT_PCA_COMPONENTS, DONOR_ID_COLUMN, MIN_CELLS_PER_GENE, MIN_GENES_PER_CELL,
    NORMALIZE_TARGET_SUM,
)
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
    download_parser = data_subparsers.add_parser(
        "download", help="Download pinned Synapse file versions (MD5-checked; present files are verified, never overwritten)"
    )
    download_parser.add_argument("--dataset", choices=["rosmap", "rosmap-mit", "all"], default="all",
                                 help="Which dataset to download (default: all)")
    download_parser.add_argument("--output-dir", type=str, default="data/raw",
                                 help="Base directory; files go to <dir>/ROSMAP and <dir>/ROSMAP_MIT (default: data/raw)")
    download_parser.add_argument("--include-raw-h5ad", action="store_true",
                                 help="Also download the 103 GB unfiltered ROSMAP-MIT h5ad (syn52392369)")

    # data verify
    verify_parser = data_subparsers.add_parser(
        "verify", help="Check that a directory holds the pinned Synapse files (size + MD5)"
    )
    verify_parser.add_argument("--dataset", choices=["rosmap", "rosmap-mit"], required=True)
    verify_parser.add_argument("--dir", type=str, required=True, help="Directory with the downloaded files")
    verify_parser.add_argument("--include-raw-h5ad", action="store_true",
                               help="Also check the unfiltered ROSMAP-MIT h5ad")

    # data prepare
    prepare_parser = data_subparsers.add_parser(
        "prepare",
        help=".h5Seurat (or converted h5ad) -> new h5ad with verified raw counts in X (+ clinical metadata)",
    )
    prepare_parser.add_argument("input", type=str,
                                help=".h5Seurat (RNA counts are read directly) or SeuratDisk-converted h5ad "
                                     "(X = Seurat data, raw = counts)")
    prepare_parser.add_argument("--output", type=str, required=True, help="Output h5ad (must differ from input)")
    prepare_parser.add_argument("--columns", type=str, nargs="*", default=["cell_type_high_resolution"],
                                help="obs columns to make categorical (default: cell_type_high_resolution)")
    prepare_parser.add_argument("--metadata", type=str, default=None, help="Clinical metadata CSV to merge into obs")
    prepare_parser.add_argument("--overwrite", action="store_true", help="Replace an existing output")

    # data category-fix
    category_parser = data_subparsers.add_parser(
        "category-fix", help="Fix categorical columns (refuses to drop .raw unless X holds counts)"
    )
    category_parser.add_argument("input", type=str, help="Path to input h5ad file")
    category_parser.add_argument("--columns", type=str, nargs="+", default=["cell_type_high_resolution"],
                                 help="Columns to convert to categorical (default: cell_type_high_resolution)")
    category_parser.add_argument("--remove-raw", action="store_true", help="Remove raw data layer")

    # data metadata
    metadata_parser = data_subparsers.add_parser("metadata", help="Add metadata to h5ad file")
    metadata_parser.add_argument("input", type=str, help="Path to input h5ad file")
    metadata_parser.add_argument("metadata", type=str, help="Path to metadata CSV file")
    metadata_parser.add_argument("--output", type=str, default=None, help="Output h5ad (default: modify input in place)")
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

    # Gene annotation commands
    genes_parser = subparsers.add_parser("genes", help="Gene ID mapping across references")
    genes_subparsers = genes_parser.add_subparsers(dest="genes_command", required=True)

    from rosmap_processing.data.genes import RECIPES

    build_map_parser = genes_subparsers.add_parser(
        "build-map", help="Map an h5ad's gene symbols to Ensembl IDs via a 10x reference GTF"
    )
    build_map_parser.add_argument("--gtf", type=str, required=True, help="GTF (.gtf or .gtf.gz) the 10x reference was built from")
    build_map_parser.add_argument("--recipe", type=str, required=True, choices=list(RECIPES),
                                  help="10x reference build recipe")
    build_map_parser.add_argument("--h5ad", type=str, required=True, help="h5ad whose var_names to map")
    build_map_parser.add_argument("--output", type=str, required=True, help="Output TSV")

    shared_parser = genes_subparsers.add_parser("shared", help="Ensembl IDs present in every given gene map")
    shared_parser.add_argument("maps", type=str, nargs="+", help="Gene map TSVs (output order follows the first)")
    shared_parser.add_argument("--output", type=str, required=True, help="Output TXT (one Ensembl ID per line)")

    # Pipeline commands
    pipeline_parser = subparsers.add_parser("pipeline", help="Processing pipeline commands")
    pipeline_subparsers = pipeline_parser.add_subparsers(dest="pipeline_command", required=True)

    # pipeline scanpy
    scanpy_parser = pipeline_subparsers.add_parser("scanpy", help="Run scanpy processing pipeline")
    scanpy_parser.add_argument("input", type=str, help="Path to input h5ad file (raw counts)")
    scanpy_parser.add_argument("output", type=str, help="Path to output h5ad file")
    scanpy_parser.add_argument("--layer", type=str, default=None, help="Layer holding the counts (default: X)")
    scanpy_parser.add_argument("--raw", action="store_true", help="Take the counts from .raw")
    scanpy_parser.add_argument("--n-genes", type=int, default=2000, help="Number of highly variable genes (default: 2000)")
    scanpy_parser.add_argument("--k-neighbors", type=int, default=30, help="Number of neighbors for KNN graph (default: 30)")
    scanpy_parser.add_argument("--individual-pca", action="store_true", help="Compute PCA per individual")
    scanpy_parser.add_argument("--import-genes", type=str, default=None,
                               help="TXT gene list (symbols or Ensembl IDs), or a reference .h5ad whose genes and order to reproduce")
    scanpy_parser.add_argument("--hvg-after-import", type=int, default=None, help="With a TXT --import-genes, select N HVGs within it (seurat_v3 on counts)")
    scanpy_parser.add_argument("--hvg-after-normalize", type=int, default=None, help="Select N HVGs on log-normalized values (seurat flavor), within --import-genes if given. Mutually exclusive with --hvg-after-import.")
    scanpy_parser.add_argument("--hvg-batch-key", type=str, default=DONOR_ID_COLUMN,
                               help=f"obs column used as batch_key for HVG selection; 'none' to pool all cells (default: {DONOR_ID_COLUMN})")
    scanpy_parser.add_argument("--gene-map", type=str, default=None, help="Gene map TSV from `genes build-map` (adds var['gene_id'])")
    scanpy_parser.add_argument("--restrict-genes", type=str, default=None,
                               help="TXT list (symbols or Ensembl IDs) all gene choices are made within, e.g. from `genes shared`")
    scanpy_parser.add_argument("--donor-column", type=str, default=None,
                               help=f"obs column holding the donor (default: '{DONOR_ID_COLUMN}', else 'projid'); cells without one are dropped")
    scanpy_parser.add_argument("--allow-missing-genes", action="store_true",
                               help="With a reference .h5ad, zero-fill its genes that are missing here instead of failing")
    scanpy_parser.add_argument("--n-pcs", type=int, default=DEFAULT_PCA_COMPONENTS,
                               help=f"PCA components for the kNN graph (default: {DEFAULT_PCA_COMPONENTS})")
    scanpy_parser.add_argument("--min-genes", type=int, default=MIN_GENES_PER_CELL,
                               help=f"Minimum genes per cell (default: {MIN_GENES_PER_CELL})")
    scanpy_parser.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_GENE,
                               help=f"Minimum cells per gene, when no gene list is imported (default: {MIN_CELLS_PER_GENE})")
    scanpy_parser.add_argument("--target-sum", type=float, default=NORMALIZE_TARGET_SUM,
                               help=f"Counts per cell after normalization (default: {NORMALIZE_TARGET_SUM:g})")
    scanpy_parser.add_argument("--overwrite", action="store_true", help="Replace an existing output")

    # pipeline hvg-compare
    compare_parser = pipeline_subparsers.add_parser(
        "hvg-compare", help="Run HVG variants (flavor x batch handling x size) and tabulate their overlaps"
    )
    compare_parser.add_argument("input", type=str, help="Path to input h5ad file (raw counts)")
    compare_parser.add_argument("--output-dir", type=str, required=True, help="Directory for gene lists and overlap tables")
    compare_parser.add_argument("--n-genes", type=int, nargs="+", default=[1000, 2000, 3000])
    compare_parser.add_argument("--flavors", type=str, nargs="+", default=["seurat_v3", "seurat"],
                                choices=["seurat_v3", "seurat"])
    compare_parser.add_argument("--batch-keys", type=str, nargs="+", default=[DONOR_ID_COLUMN, "none"],
                                help=f"Batch keys to compare; 'none' pools all cells (default: '{DONOR_ID_COLUMN}' none)")
    compare_parser.add_argument("--layer", type=str, default=None, help="Layer holding the counts (default: X)")
    compare_parser.add_argument("--gene-map", type=str, default=None, help="Gene map TSV from `genes build-map`")
    compare_parser.add_argument("--restrict-genes", type=str, default=None, help="TXT list all selections are made within")
    compare_parser.add_argument("--donor-column", type=str, default=None,
                                help=f"obs column holding the donor (default: '{DONOR_ID_COLUMN}', else 'projid')")
    compare_parser.add_argument("--compare", type=str, nargs="*", default=[], metavar="NAME=PATH",
                                help="Extra gene lists to include: TXT, or an h5ad (its genes)")
    compare_parser.add_argument("--min-genes", type=int, default=MIN_GENES_PER_CELL)
    compare_parser.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_GENE)
    compare_parser.add_argument("--target-sum", type=float, default=NORMALIZE_TARGET_SUM)

    # pipeline mathys2023-for-adprs
    from rosmap_processing.pipelines.mathys2023_for_adprs import DEFAULT_INPUT_DIR, DEFAULT_OUTPUT_DIR

    adprs_parser = pipeline_subparsers.add_parser(
        "mathys2023-for-adprs",
        help="Aggregate Mathys 2023 per-class h5ads into per-(donor, cell_type, gene) stats for ad-prs",
    )
    adprs_parser.add_argument(
        "--input-dir", type=Path, default=Path(DEFAULT_INPUT_DIR),
        help=f"Directory containing per-class count h5ad files (default: {DEFAULT_INPUT_DIR})",
    )
    adprs_parser.add_argument(
        "--output-dir", type=Path, default=Path(DEFAULT_OUTPUT_DIR),
        help=f"Directory to write outputs into (default: {DEFAULT_OUTPUT_DIR})",
    )
    adprs_parser.add_argument(
        "--clinical-csv", type=Path, default=Path("data/v2/download/ROSMAP/ROSMAP_clinical.csv"),
        help="ROSMAP clinical CSV (default: data/v2/download/ROSMAP/ROSMAP_clinical.csv)",
    )
    adprs_parser.add_argument(
        "--mit-metadata-csv", type=Path,
        default=Path("data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv"),
        help="MIT individual metadata CSV (for individualID -> subject mapping)",
    )
    adprs_parser.add_argument(
        "--gene-map", type=Path, default=None,
        help="Gene map TSV from `genes build-map`; adds an Ensembl gene_id column",
    )
    adprs_parser.add_argument(
        "--keep-intermediate", action="store_true",
        help="Keep per-class parquet chunks under _intermediate/ after merging",
    )
    adprs_parser.add_argument(
        "--resume", action="store_true",
        help="If the per-donor parquet already exists, skip aggregation and only rebuild donor_metadata.tsv + README.md from it",
    )

    # Utils commands
    utils_parser = subparsers.add_parser("utils", help="Utility commands")
    utils_subparsers = utils_parser.add_subparsers(dest="utils_command", required=True)

    # utils validate
    validate_parser = utils_subparsers.add_parser("validate", help="Validate and inspect h5ad file")
    validate_parser.add_argument("input", type=str, help="Path to h5ad file to validate")

    # utils check-counts
    check_parser = utils_subparsers.add_parser(
        "check-counts", help="Check that a stored matrix holds raw counts (aggregate statistics only)"
    )
    check_parser.add_argument("inputs", type=str, nargs="+", help="h5ad files")
    check_parser.add_argument("--key", type=str, default="X", help="Matrix to check: X, raw/X or layers/<name> (default: X)")

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
        elif args.command == "genes":
            handle_genes_command(args)
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
        from rosmap_processing.data.download import download_dataset
        datasets = ["rosmap", "rosmap-mit"] if args.dataset == "all" else [args.dataset]
        for dataset in datasets:
            download_dataset(dataset, output_base=args.output_dir, include_optional=args.include_raw_h5ad)

    elif args.data_command == "verify":
        from rosmap_processing.data.download import DATASETS, verify_files
        problems = verify_files(DATASETS[args.dataset][1], args.dir, include_optional=args.include_raw_h5ad)
        if problems:
            raise ValueError(f"{len(problems)} file(s) in {args.dir} do not match the pinned Synapse versions")
        logger.info(f"All {args.dataset} files in {args.dir} match the pinned Synapse versions")

    elif args.data_command == "prepare":
        from rosmap_processing.data.prepare import prepare_file
        logger.info(f"Preparing counts from {args.input}")
        prepare_file(
            input_path=Path(args.input),
            output_path=Path(args.output),
            categorical_columns=args.columns,
            metadata_path=Path(args.metadata) if args.metadata else None,
            overwrite=args.overwrite,
        )

    elif args.data_command == "category-fix":
        from rosmap_processing.data.category_fix import fix_categories_in_file
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
            is_mit=args.mit,
            output_path=args.output,
        )


def handle_core_command(args):
    """Handle core subcommands."""
    logger = get_logger(__name__)

    if args.core_command == "combine":
        from rosmap_processing.core.combine import combine_and_save
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


def handle_genes_command(args):
    """Handle genes subcommands."""
    if args.genes_command == "build-map":
        from rosmap_processing.data.genes import build_gene_map_file
        build_gene_map_file(
            gtf_path=args.gtf, recipe=args.recipe, h5ad_path=args.h5ad, output_path=args.output
        )
    elif args.genes_command == "shared":
        from rosmap_processing.data.genes import write_shared_gene_ids
        write_shared_gene_ids(args.maps, args.output)


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
            hvg_batch_key=args.hvg_batch_key,
            gene_map=args.gene_map,
            restrict_genes=args.restrict_genes,
            donor_column=args.donor_column,
            allow_missing_genes=args.allow_missing_genes,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            n_pca_components=args.n_pcs,
            target_sum=args.target_sum,
            overwrite=args.overwrite,
        )

    elif args.pipeline_command == "hvg-compare":
        from rosmap_processing.pipelines.hvg_compare import compare_hvgs
        extra = {}
        for item in args.compare:
            if "=" not in item:
                raise ValueError(f"--compare expects NAME=PATH, got {item!r}")
            name, path = item.split("=", 1)
            extra[name] = path
        compare_hvgs(
            input_path=args.input,
            output_dir=args.output_dir,
            n_genes=args.n_genes,
            flavors=args.flavors,
            batch_keys=args.batch_keys,
            layer=args.layer,
            gene_map=args.gene_map,
            restrict_genes=args.restrict_genes,
            donor_column=args.donor_column,
            extra_lists=extra,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            target_sum=args.target_sum,
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
            resume=args.resume,
            gene_map_path=args.gene_map,
        )


def handle_utils_command(args):
    """Handle utils subcommands."""
    logger = get_logger(__name__)

    if args.utils_command == "validate":
        from rosmap_processing.utils.validation import inspect_h5ad
        logger.info(f"Validating {args.input}")
        inspect_h5ad(args.input)

    elif args.utils_command == "check-counts":
        from rosmap_processing.utils.counts import check_h5ad_counts
        failed = []
        for path in args.inputs:
            stats = check_h5ad_counts(path, key=args.key)
            verdict = "raw counts" if stats["is_counts"] else "NOT raw counts"
            logger.info(
                f"{path} [{stats['key']}]: {verdict}; shape={stats['shape']}, dtype={stats['dtype']}, "
                f"stored={stats['stored_values']:,}, non-integer={stats['non_integer']:,}, "
                f"negative={stats['negative']:,}, min={stats['min']}, max={stats['max']}"
            )
            if not stats["is_counts"]:
                failed.append(path)
        if failed:
            raise ValueError(f"{len(failed)} file(s) do not hold raw counts in {args.key}")


if __name__ == "__main__":
    main()
