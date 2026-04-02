"""Core scanpy processing pipeline for ROSMAP data."""

import scanpy as sc
import anndata as ad
from pathlib import Path
from typing import Optional, List, Set, Union
from tqdm import tqdm

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rosmap_processing.utils.constants import (
    MIN_GENES_PER_CELL,
    MIN_CELLS_PER_GENE,
    DEFAULT_HVG_COUNT,
    DEFAULT_K_NEIGHBORS,
    DEFAULT_PCA_COMPONENTS,
    CPM_TARGET,
    DONOR_ID_COLUMN,
    VALID_H5AD_EXTENSIONS,
    VALID_TXT_EXTENSIONS,
)
from rosmap_processing.utils.logging import get_logger

logger = get_logger(__name__)


def validate_file_path(path: Path, extensions: List[str]) -> None:
    """
    Validate file path and extension.
    
    Parameters
    ----------
    path : Path
        Path to validate
    extensions : List[str]
        List of valid extensions
        
    Raises
    ------
    FileNotFoundError
        If file doesn't exist
    ValueError
        If file has invalid extension
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    if path.suffix not in extensions:
        raise ValueError(
            f"Invalid file extension: {path.suffix}. "
            f"Expected one of: {extensions}"
        )


def detect_hvg_flavor(adata: ad.AnnData, sample_cells: int = 500, sample_genes: int = 200) -> str:
    """
    Detect whether X contains raw counts or log-normalized data and return
    the appropriate HVG flavor for sc.pp.highly_variable_genes.

    Returns 'seurat_v3' for integer count data (stored as floats or ints),
    or 'seurat' for log-normalized continuous data.
    """
    import numpy as np
    n_cells = min(sample_cells, adata.shape[0])
    n_genes = min(sample_genes, adata.shape[1])
    sample = adata.X[:n_cells, :n_genes]
    if hasattr(sample, 'toarray'):
        sample = sample.toarray()
    else:
        sample = np.array(sample)

    nonzero = sample[sample > 0].flatten()
    if len(nonzero) == 0:
        logger.warning("No non-zero values in sample — defaulting to seurat_v3")
        return 'seurat_v3'

    pct_whole = np.sum(nonzero == np.floor(nonzero)) / len(nonzero)
    if pct_whole > 0.99:
        logger.info(f"Data detected as raw counts (%.1f%% whole numbers) — using seurat_v3", pct_whole * 100)
        return 'seurat_v3'
    else:
        logger.info(f"Data detected as log-normalized (%.1f%% whole numbers) — using seurat", pct_whole * 100)
        return 'seurat'


def load_gene_list(gene_file: Path) -> Set[str]:
    """
    Load gene list from file.
    
    Parameters
    ----------
    gene_file : Path
        Path to gene file (TXT or h5ad)
        
    Returns
    -------
    Set[str]
        Set of gene names
        
    Raises
    ------
    ValueError
        If file format is not supported
    """
    logger.info(f"Loading genes from: {gene_file}")
    
    if gene_file.suffix in VALID_TXT_EXTENSIONS:
        with open(gene_file, 'r') as f:
            genes = {line.strip() for line in f if line.strip()}
        logger.info(f"Loaded {len(genes)} genes from TXT file")
        return genes
    
    elif gene_file.suffix in VALID_H5AD_EXTENSIONS:
        adata_genes = ad.read_h5ad(gene_file)
        genes = set(adata_genes.var_names)
        logger.info(f"Loaded {len(genes)} genes from h5ad file")
        return genes
    
    else:
        raise ValueError(
            f"Unsupported gene file format: {gene_file.suffix}. "
            f"Expected TXT or h5ad."
        )


def process_h5ad(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    layer: Optional[str] = None,
    use_raw: bool = False,
    n_hvgs: int = DEFAULT_HVG_COUNT,
    k_neighbors: int = DEFAULT_K_NEIGHBORS,
    individual_pca: bool = False,
    import_genes: Optional[Union[str, Path]] = None,
    hvg_after_import: Optional[int] = None,
    min_genes: int = MIN_GENES_PER_CELL,
    min_cells: int = MIN_CELLS_PER_GENE,
    n_pca_components: int = DEFAULT_PCA_COMPONENTS,
    target_sum: float = CPM_TARGET,
) -> None:
    """
    Process h5ad file through scanpy pipeline.
    
    Parameters
    ----------
    input_path : str or Path
        Path to input h5ad file
    output_path : str or Path
        Path to save processed h5ad file
    layer : str, optional
        Layer to use for processing. If None, uses main X matrix
    use_raw : bool, default False
        Whether to use raw data
    n_hvgs : int, default from constants
        Number of highly variable genes to select
    k_neighbors : int, default from constants
        Number of neighbors for KNN graph
    individual_pca : bool, default False
        Whether to compute PCA per donor
    import_genes : str or Path, optional
        Path to gene list file (TXT or h5ad)
    min_genes : int, default from constants
        Minimum genes per cell
    min_cells : int, default from constants
        Minimum cells per gene
    n_pca_components : int, default from constants
        Number of PCA components
    target_sum : float, default from constants
        Target sum for normalization (CPM)
        
    Raises
    ------
    FileNotFoundError
        If input file doesn't exist
    ValueError
        If layer/raw data is not available
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    
    logger.info("="*60)
    logger.info("Starting scanpy processing pipeline")
    logger.info("="*60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_path}")
    logger.info(f"Layer: {layer}")
    logger.info(f"Use raw: {use_raw}")
    logger.info(f"HVGs: {n_hvgs}")
    logger.info(f"K-neighbors: {k_neighbors}")
    logger.info(f"Individual PCA: {individual_pca}")
    if import_genes:
        logger.info(f"Import genes from: {import_genes}")
    if hvg_after_import is not None:
        logger.info(f"HVG selection after import: {hvg_after_import} genes")
    logger.info("")
    
    # Validate input file
    validate_file_path(input_path, VALID_H5AD_EXTENSIONS)
    
    # Create output directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Load data
    logger.info(f"Loading h5ad file: {input_path}")
    adata = ad.read_h5ad(input_path)
    logger.info(f"Loaded: {adata}")

    # Diagnostic: show a small slice of X to reveal value types
    import numpy as np
    x_slice = adata.X[:5, :5]
    if hasattr(x_slice, 'toarray'):
        x_slice = x_slice.toarray()
    logger.info(f"X[0:5, 0:5] sample values:\n{x_slice}")
    logger.info(f"X dtype: {adata.X.dtype}, min={adata.X.min():.4f}, max={adata.X.max():.4f}")
    if adata.layers:
        logger.info(f"Available layers: {list(adata.layers.keys())}")
    
    # Handle raw data
    if use_raw:
        logger.info("Using raw data")
        if adata.raw is not None:
            adata.X = adata.raw.X.copy()
            adata.var_names = adata.raw.var_names
            adata.obs_names = adata.raw.obs_names
            logger.info("Raw data copied to main matrix")
        else:
            raise ValueError("Raw data is not available in the AnnData object")
    
    # Handle layer selection
    if layer is not None:
        logger.info(f"Using layer: {layer}")
        if layer in adata.layers:
            adata.X = adata.layers[layer].copy()
            del adata.layers[layer]
            if len(adata.layers) == 0:
                del adata.layers
            logger.info(f"Layer '{layer}' moved to main matrix")
        else:
            available_layers = list(adata.layers.keys()) if adata.layers else []
            raise ValueError(
                f"Layer '{layer}' not found. "
                f"Available layers: {available_layers}"
            )
    
    # Filtering
    logger.info("Filtering cells and genes")
    logger.info(f"Shape before filtering: {adata.shape}")
    
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    if not import_genes:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    else:
        logger.info("Skipping gene filtering (importing gene list)")
    
    logger.info(f"Shape after filtering: {adata.shape}")
    
    # Gene selection
    if import_genes:
        import_genes = Path(import_genes)
        select_genes = load_gene_list(import_genes)
        
        # Check overlap
        overlap = set(select_genes) & set(adata.var_names)
        if len(overlap) < len(select_genes):
            missing = len(select_genes) - len(overlap)
            logger.warning(
                f"{missing} genes from import list not found in data. "
                f"Using {len(overlap)} overlapping genes."
            )
        
        logger.info(f"Selecting {len(overlap)} genes from import list")
        adata = adata[:, adata.var_names.isin(overlap)].copy()
        logger.info(f"Shape after gene selection: {adata.shape}")

        if hvg_after_import is not None:
            logger.info(f"Selecting {hvg_after_import} highly variable genes within imported set")
            hvg_flavor = detect_hvg_flavor(adata)
            sc.pp.highly_variable_genes(
                adata,
                flavor=hvg_flavor,
                n_top_genes=hvg_after_import,
                subset=False,
            )

            # Write full ranked HVG list before slicing
            hvg_scores_path = output_path.with_name(output_path.stem + "_hvg_scores.tsv")
            score_cols = [c for c in ['highly_variable', 'highly_variable_rank', 'dispersions_norm', 'dispersions', 'variances_norm', 'variances'] if c in adata.var.columns]
            hvg_df = adata.var[score_cols].copy()
            hvg_df = hvg_df.sort_values('highly_variable_rank')
            hvg_df.to_csv(hvg_scores_path, sep='\t')
            logger.info(f"Wrote full HVG scores ({len(hvg_df)} genes) to {hvg_scores_path}")

            hvg_names_path = output_path.with_name(output_path.stem + "_hvg_names.txt")
            hvg_names_path.write_text('\n'.join(hvg_df.index) + '\n')
            logger.info(f"Wrote full HVG gene names to {hvg_names_path}")

            adata = adata[:, adata.var['highly_variable']].copy()
            logger.info(f"Shape after HVG selection: {adata.shape}")

    else:
        logger.info(f"Selecting {n_hvgs} highly variable genes")
        hvg_flavor = detect_hvg_flavor(adata)
        sc.pp.highly_variable_genes(
            adata,
            flavor=hvg_flavor,
            n_top_genes=n_hvgs,
            subset=True
        )
        n_selected = adata.shape[1]
        logger.info(f"Selected {n_selected} highly variable genes")
    
    # Normalization
    logger.info(f"Normalizing to {target_sum:.0e} total counts (CPM)")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    
    logger.info("Log-transforming data (log1p)")
    sc.pp.log1p(adata)
    
    # PCA (global if not individual)
    if not individual_pca:
        logger.info(f"Computing PCA ({n_pca_components} components)")
        sc.pp.pca(adata, n_comps=n_pca_components)
        logger.info("PCA complete")
    else:
        logger.info("Skipping global PCA (will compute per donor)")
    
    # Build KNN graph per donor
    logger.info("Building KNN graphs per donor")
    
    # Ensure donor ID column exists
    if DONOR_ID_COLUMN not in adata.obs.columns:
        if "projid" in adata.obs.columns:
            logger.info(f"Using 'projid' as {DONOR_ID_COLUMN}")
            adata.obs[DONOR_ID_COLUMN] = adata.obs["projid"]
        else:
            raise ValueError(
                f"Neither '{DONOR_ID_COLUMN}' nor 'projid' found in obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
    
    slices = []
    donors = adata.obs[DONOR_ID_COLUMN].unique()
    skipped_donors = []
    
    for donor in tqdm(donors, desc="Processing donors"):
        donor_slice = adata[adata.obs[DONOR_ID_COLUMN] == donor].copy()
        
        if donor_slice.shape[0] <= 1:
            logger.warning(
                f"Skipping donor {donor}: only {donor_slice.shape[0]} cell(s). "
                f"Cannot construct KNN graph."
            )
            skipped_donors.append(donor)
            continue
        
        # Per-donor PCA if requested
        if individual_pca:
            sc.pp.pca(donor_slice, n_comps=n_pca_components)
        
        # Build KNN graph
        sc.pp.neighbors(donor_slice, n_neighbors=k_neighbors, use_rep="X_pca")
        slices.append(donor_slice)
    
    if skipped_donors:
        logger.warning(
            f"Skipped {len(skipped_donors)} donors with insufficient cells: "
            f"{skipped_donors}"
        )
    
    logger.info(f"Processed {len(slices)} donors successfully")
    
    # Concatenate all donor slices
    logger.info("Concatenating donor data")
    del adata
    adata = ad.concat(
        slices,
        merge="same",
        uns_merge="same",
        label="all",
        index_unique="-",
        pairwise=True,
    )
    logger.info(f"Final shape: {adata.shape}")
    
    # Save
    logger.info(f"Saving processed data to: {output_path}")
    adata.write_h5ad(output_path, compression="gzip")
    
    logger.info("="*60)
    logger.info("Processing complete!")
    logger.info("="*60)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Performs the basic scanpy processing steps on a given h5ad file."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to the h5ad file to process."
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to save the processed h5ad file."
    )
    parser.add_argument(
        "--layer",
        type=str,
        default=None,
        help="Layer to select for processing. If None, the main data will be used."
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="If set, the raw data will be used instead of the main data."
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        default=DEFAULT_HVG_COUNT,
        help=f"Number of highly variable genes to select. Default is {DEFAULT_HVG_COUNT}."
    )
    parser.add_argument(
        "--k-neighbors",
        type=int,
        default=DEFAULT_K_NEIGHBORS,
        help=f"Number of neighbors for the KNN graph. Default is {DEFAULT_K_NEIGHBORS}."
    )
    parser.add_argument(
        "--individual-pca",
        action="store_true",
        help="If set, PCA for the KNN graph will be computed for each individual separately."
    )
    parser.add_argument(
        "--import-genes",
        type=str,
        default=None,
        help="Path to a TXT or AnnData file to use the genes from. This overrides HVG selection."
    )
    parser.add_argument(
        "--hvg-after-import",
        type=int,
        default=None,
        help="If set alongside --import-genes, further select N HVGs within the imported gene subset."
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
    from rosmap_processing.utils.logging import setup_logging
    setup_logging(level=args.log_level)
    
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
            hvg_after_import=args.hvg_after_import,
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)
