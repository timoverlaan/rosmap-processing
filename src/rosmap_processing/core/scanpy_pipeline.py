"""Core scanpy processing pipeline for ROSMAP data.

1. Load and take raw counts (X, a layer, or .raw); the whole matrix must hold
   integer counts.
2. Filter cells (min genes) and, unless a gene list is imported, genes (min cells).
3. Per-cell total counts over all remaining genes. These are the size factors
   for normalization, so the result does not depend on which genes are kept.
4. Optionally attach Ensembl IDs from a gene map (see `data/genes.py`).
5. Choose genes: an imported list (TXT), or a reference h5ad's genes in its
   exact order (missing genes are an error), and/or highly variable genes
   (seurat_v3 on counts, or seurat on log-normalized values), per batch.
6. Subset genes; X = log1p(counts / total * target_sum); counts are kept in
   layers['counts'].
7. PCA (global or per donor) and a kNN graph per donor.
"""

import sys
from collections import Counter
from pathlib import Path
from typing import List, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from tqdm import tqdm

from rosmap_processing.utils.constants import (
    MIN_GENES_PER_CELL,
    MIN_CELLS_PER_GENE,
    DEFAULT_HVG_COUNT,
    DEFAULT_K_NEIGHBORS,
    DEFAULT_PCA_COMPONENTS,
    NORMALIZE_TARGET_SUM,
    DONOR_ID_COLUMN,
    VALID_H5AD_EXTENSIONS,
    VALID_TXT_EXTENSIONS,
)
from rosmap_processing.utils.counts import (
    assert_integer_counts,
    log_normalize,
    to_int32_counts,
    total_counts,
)
from rosmap_processing.utils.logging import get_logger
from rosmap_processing.utils.provenance import provenance_record

logger = get_logger(__name__)


def validate_file_path(path: Path, extensions: List[str]) -> None:
    """
    Validate file path and extension.

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


def load_gene_list(gene_file: Path) -> List[str]:
    """Gene names or Ensembl IDs from a TXT file (one per line), in file order."""
    with open(gene_file, "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    logger.info(f"Loaded {len(genes)} genes from {gene_file}")
    return list(dict.fromkeys(genes))


def _select_counts(adata: ad.AnnData, layer: Optional[str], use_raw: bool) -> ad.AnnData:
    """Put the requested count matrix in X and drop everything else."""
    if use_raw and layer is not None:
        raise ValueError("Use either --raw or --layer, not both")
    if use_raw:
        if adata.raw is None:
            raise ValueError("Raw data is not available in the AnnData object")
        logger.info("Using .raw as counts")
        adata = ad.AnnData(X=adata.raw.X, obs=adata.obs, var=adata.raw.var)
    elif layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")
        logger.info(f"Using layer '{layer}' as counts")
        adata = ad.AnnData(X=adata.layers[layer], obs=adata.obs, var=adata.var)
    else:
        adata = ad.AnnData(X=adata.X, obs=adata.obs, var=adata.var)

    assert_integer_counts(adata.X, name="input counts")
    # float32 holds counts exactly and is what scanpy's HVG code expects.
    adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
    return adata


def _ensure_donor_column(adata: ad.AnnData, donor_column: Optional[str] = None) -> ad.AnnData:
    """Fill obs['Donor ID'] (from `donor_column`, else existing, else projid).

    Cells without a donor are dropped (they cannot go into a per-donor graph);
    only their number is logged.
    """
    if donor_column is not None:
        if donor_column not in adata.obs.columns:
            raise ValueError(f"Donor column '{donor_column}' not in obs")
        logger.info(f"Using '{donor_column}' as {DONOR_ID_COLUMN}")
        adata.obs[DONOR_ID_COLUMN] = adata.obs[donor_column]
    elif DONOR_ID_COLUMN not in adata.obs.columns:
        if "projid" not in adata.obs.columns:
            raise ValueError(f"Neither '{DONOR_ID_COLUMN}' nor 'projid' found in obs")
        logger.info(f"Using 'projid' as {DONOR_ID_COLUMN}")
        adata.obs[DONOR_ID_COLUMN] = adata.obs["projid"]

    missing = adata.obs[DONOR_ID_COLUMN].isna().to_numpy()
    if missing.any():
        logger.warning(f"Dropping {int(missing.sum())} cells without a donor ID")
        adata = adata[~missing].copy()
    return adata


def _attach_gene_ids(adata: ad.AnnData, gene_map: Union[str, Path]) -> None:
    from rosmap_processing.data.genes import load_gene_map

    mapping = load_gene_map(gene_map)
    # "" (not None) for unmapped genes: h5ad string columns cannot hold None.
    adata.var["gene_id"] = [mapping.get(g, "") for g in adata.var_names]
    n = int((adata.var["gene_id"] != "").sum())
    logger.info(f"Attached Ensembl IDs from {gene_map}: {n}/{adata.n_vars} genes mapped")


def _genes_from_list(adata: ad.AnnData, entries: List[str]) -> List[str]:
    """Dataset genes named in a TXT list (symbols, or Ensembl IDs), in dataset order."""
    n_ids = sum(e.startswith("ENS") for e in entries)
    if n_ids > len(entries) / 2:
        if "gene_id" not in adata.var.columns:
            raise ValueError("The gene list holds Ensembl IDs; pass --gene-map for this dataset")
        wanted = set(entries)
        keep = adata.var_names[adata.var["gene_id"].isin(wanted).to_numpy()].tolist()
        n_found = adata.var.loc[keep, "gene_id"].nunique()
        key = "Ensembl IDs"
    else:
        wanted = set(entries)
        keep = [g for g in adata.var_names if g in wanted]
        n_found = len(keep)
        key = "gene symbols"
    if n_found < len(entries):
        logger.warning(
            f"{len(entries) - n_found} of {len(entries)} listed {key} are not in the data; "
            f"using the {n_found} present"
        )
    logger.info(f"Gene list matched on {key}: {len(keep)} genes")
    return keep


def align_to_reference(
    adata: ad.AnnData, ref_var: pd.DataFrame, allow_missing: bool = False
) -> ad.AnnData:
    """Return adata with exactly the reference's genes, in the reference's order.

    Genes are matched on Ensembl ID where both sides have one, otherwise on
    symbol (only if the dataset's gene with that symbol has no conflicting ID).
    var_names become the reference's names; the dataset's own symbol is kept in
    var['symbol_dataset']. Missing genes raise ValueError, or with
    allow_missing are added as all-zero columns.
    """
    def _id(value) -> str:
        return value if isinstance(value, str) else ""

    ref_names = ref_var.index.astype(str).tolist()
    ref_ids = [_id(g) for g in ref_var["gene_id"]] if "gene_id" in ref_var.columns else [""] * len(ref_names)
    data_ids = (
        pd.Series([_id(g) for g in adata.var["gene_id"]], index=adata.var_names)
        if "gene_id" in adata.var.columns
        else pd.Series("", index=adata.var_names)
    )

    by_id = {}
    for name, gid in data_ids.items():
        if gid:
            by_id.setdefault(gid, []).append(name)

    picked: List[Optional[str]] = []
    methods: List[str] = []
    for ref_name, ref_gid in zip(ref_names, ref_ids):
        if ref_gid and len(by_id.get(ref_gid, [])) == 1:
            picked.append(by_id[ref_gid][0])
            methods.append("gene_id")
        elif ref_name in data_ids.index and (
            not ref_gid or not data_ids[ref_name] or data_ids[ref_name] == ref_gid
        ):
            # Same symbol, and no conflicting Ensembl ID on either side.
            picked.append(ref_name)
            methods.append("symbol")
        else:
            picked.append(None)
            methods.append("missing_zero_filled")

    present = [p for p in picked if p is not None]
    if len(set(present)) != len(present):
        raise ValueError("Several reference genes matched the same dataset gene")

    missing = [n for n, p in zip(ref_names, picked) if p is None]
    logger.info(f"Matched {len(ref_names)} reference genes: {dict(Counter(methods))}")
    if missing and not allow_missing:
        raise ValueError(
            f"{len(missing)} reference genes are not in this dataset "
            f"(first: {missing[:20]}). Select genes within the genes both datasets "
            f"share (genes shared), or pass --allow-missing-genes to zero-fill them."
        )
    if missing:
        logger.warning(f"Zero-filling {len(missing)} reference genes missing from this dataset")

    cols = [adata.var_names.get_loc(p) for p in present]
    X = sp.csr_matrix(adata.X)[:, cols]
    if missing:
        X = sp.hstack([X, sp.csr_matrix((adata.n_obs, len(missing)), dtype=X.dtype)], format="csr")
        # Columns are [present..., missing...]; reorder to reference order.
        k_present, k_missing = iter(range(len(present))), iter(range(len(present), X.shape[1]))
        order = [next(k_present) if p is not None else next(k_missing) for p in picked]
        X = X[:, order]

    var = pd.DataFrame(index=pd.Index(ref_names))
    var["gene_id"] = ref_ids
    var["symbol_dataset"] = [p if p is not None else "" for p in picked]
    var["match"] = methods
    return ad.AnnData(X=X, obs=adata.obs, var=var)


def _select_hvgs(
    adata: ad.AnnData,
    n_top_genes: int,
    flavor: str,
    batch_key: Optional[str],
    size_factors: np.ndarray,
    target_sum: float,
    normalized=None,
) -> pd.DataFrame:
    """Run scanpy HVG selection; returns the per-gene result table (all genes).

    flavor 'seurat_v3' works on the counts in X; 'seurat' on log-normalized
    values, computed here with the given size factors unless `normalized`
    (the same matrix, precomputed) is passed.
    """
    # scanpy's batched HVG code needs a categorical batch column; donor IDs
    # derived from projid are integers.
    obs = pd.DataFrame(index=adata.obs_names)
    if batch_key:
        obs[batch_key] = pd.Categorical(adata.obs[batch_key].astype(str))

    if flavor == "seurat_v3":
        X = adata.X
    elif flavor == "seurat":
        X = normalized if normalized is not None else log_normalize(adata.X, size_factors, target_sum)
    else:
        raise ValueError(f"Unsupported HVG flavor: {flavor}")
    work = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=adata.var_names))

    logger.info(
        f"Selecting {n_top_genes} HVGs among {adata.n_vars} genes "
        f"(flavor={flavor}, batch_key={batch_key})"
    )
    sc.pp.highly_variable_genes(
        work, flavor=flavor, n_top_genes=n_top_genes, batch_key=batch_key, subset=False
    )
    table = work.var.copy()

    # The pooled seurat flavor keeps every gene tied at the cutoff, so it can
    # return a few more than n_top_genes; keep exactly the top n.
    n_selected = int(table["highly_variable"].sum())
    if n_selected > n_top_genes:
        keep = rank_hvg_table(table).index[:n_top_genes]
        table["highly_variable"] = table.index.isin(keep)
        logger.info(f"{n_selected} genes tied at the cutoff; kept the top {n_top_genes}")
    return table


def rank_hvg_table(hvg: pd.DataFrame) -> pd.DataFrame:
    """Sort an HVG table the way scanpy chose its genes: selected genes first.

    seurat_v3: by rank (median rank across batches, after the number of
    batches a gene is variable in); seurat: by normalized dispersion (after
    the number of batches).
    """
    keys, ascending = [], []
    if "highly_variable_nbatches" in hvg.columns:
        keys.append("highly_variable_nbatches")
        ascending.append(False)
    if "highly_variable_rank" in hvg.columns:
        keys.append("highly_variable_rank")
        ascending.append(True)
    else:
        keys.append("dispersions_norm")
        ascending.append(False)
    ranked = hvg.sort_values(keys, ascending=ascending, na_position="last", kind="stable")
    return pd.concat([ranked[ranked["highly_variable"]], ranked[~ranked["highly_variable"]]])


def _write_hvg_tables(hvg: pd.DataFrame, output_path: Path) -> None:
    """Write the full HVG table and the ranked gene names next to the output."""
    ranked = rank_hvg_table(hvg)

    scores_path = output_path.with_name(output_path.stem + "_hvg_scores.tsv")
    ranked.to_csv(scores_path, sep="\t")
    names_path = output_path.with_name(output_path.stem + "_hvg_names.txt")
    names_path.write_text("\n".join(ranked.index) + "\n")
    logger.info(f"Wrote HVG scores ({len(ranked)} genes) to {scores_path} and ranked names to {names_path}")


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
    hvg_after_normalize: Optional[int] = None,
    hvg_batch_key: Optional[str] = DONOR_ID_COLUMN,
    gene_map: Optional[Union[str, Path]] = None,
    restrict_genes: Optional[Union[str, Path]] = None,
    donor_column: Optional[str] = None,
    allow_missing_genes: bool = False,
    min_genes: int = MIN_GENES_PER_CELL,
    min_cells: int = MIN_CELLS_PER_GENE,
    n_pca_components: int = DEFAULT_PCA_COMPONENTS,
    target_sum: float = NORMALIZE_TARGET_SUM,
    overwrite: bool = False,
) -> None:
    """
    Process an h5ad file of raw counts into a per-donor kNN graph dataset.

    Parameters
    ----------
    input_path, output_path : str or Path
    layer : str, optional
        Layer holding the counts (default: X)
    use_raw : bool
        Take the counts from .raw
    n_hvgs : int
        Number of HVGs when no gene list is imported (seurat_v3 on counts)
    k_neighbors : int
        Neighbors in the per-donor kNN graph
    individual_pca : bool
        Compute PCA per donor instead of once globally
    import_genes : str or Path, optional
        TXT list of symbols or Ensembl IDs, or a reference h5ad whose genes
        (and their order) the output must reproduce exactly
    hvg_after_import : int, optional
        Select N HVGs (seurat_v3, counts) within the imported TXT list
    hvg_after_normalize : int, optional
        Select N HVGs (seurat, log-normalized) within the imported TXT list,
        or within all genes when nothing is imported
    hvg_batch_key : str or None
        obs column passed to scanpy as batch_key for HVG selection
    gene_map : str or Path, optional
        TSV from `genes build-map`, adds var['gene_id'] for matching
    restrict_genes : str or Path, optional
        TXT list (symbols or Ensembl IDs); all gene choices are made within it,
        e.g. the genes shared with another dataset (`genes shared`)
    donor_column : str, optional
        obs column holding the donor (default: 'Donor ID', else 'projid');
        cells without a donor are dropped
    allow_missing_genes : bool
        With a reference h5ad, zero-fill reference genes missing here
    target_sum : float
        Counts per cell after scaling (totals over all genes)
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if hvg_after_import is not None and hvg_after_normalize is not None:
        raise ValueError("--hvg-after-import and --hvg-after-normalize are mutually exclusive")
    if hvg_batch_key in ("", "none", "None"):
        hvg_batch_key = None

    reference_h5ad = import_genes is not None and Path(import_genes).suffix in VALID_H5AD_EXTENSIONS
    if import_genes is not None:
        validate_file_path(Path(import_genes), VALID_H5AD_EXTENSIONS + VALID_TXT_EXTENSIONS)
        if reference_h5ad and (hvg_after_import or hvg_after_normalize):
            raise ValueError("HVG selection after import needs a TXT gene list, not a reference h5ad")
    if hvg_after_import is not None and import_genes is None:
        raise ValueError("--hvg-after-import requires --import-genes")
    if restrict_genes is not None:
        validate_file_path(Path(restrict_genes), VALID_TXT_EXTENSIONS)
        if reference_h5ad:
            raise ValueError("--restrict-genes cannot be combined with a reference h5ad")

    logger.info("=" * 60)
    logger.info("Starting scanpy processing pipeline")
    logger.info("=" * 60)
    params = dict(
        input=input_path, layer=layer, use_raw=use_raw, n_hvgs=n_hvgs,
        k_neighbors=k_neighbors, individual_pca=individual_pca, import_genes=import_genes,
        hvg_after_import=hvg_after_import, hvg_after_normalize=hvg_after_normalize,
        hvg_batch_key=hvg_batch_key, gene_map=gene_map, restrict_genes=restrict_genes,
        donor_column=donor_column,
        allow_missing_genes=allow_missing_genes,
        min_genes=min_genes, min_cells=min_cells, n_pca_components=n_pca_components,
        target_sum=target_sum,
    )
    for key, value in params.items():
        logger.info(f"{key}: {value}")

    validate_file_path(input_path, VALID_H5AD_EXTENSIONS)
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path} (pass --overwrite to replace it)")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading h5ad file: {input_path}")
    adata = _select_counts(ad.read_h5ad(input_path), layer, use_raw)
    logger.info(f"Counts verified: {adata.n_obs} cells x {adata.n_vars} genes")

    # Filtering
    logger.info(f"Shape before filtering: {adata.shape}")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if import_genes is None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    else:
        logger.info("Skipping gene filtering (importing gene list)")
    logger.info(f"Shape after filtering: {adata.shape}")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"Nothing left after filtering (min_genes={min_genes}, min_cells={min_cells})")

    adata = _ensure_donor_column(adata, donor_column)
    if hvg_batch_key is not None and hvg_batch_key not in adata.obs.columns:
        raise ValueError(f"HVG batch key '{hvg_batch_key}' not in obs")

    # Size factors over all genes, before any gene subsetting.
    adata.obs["total_counts"] = total_counts(adata.X)

    if gene_map is not None:
        _attach_gene_ids(adata, gene_map)

    if restrict_genes is not None:
        keep = _genes_from_list(adata, load_gene_list(Path(restrict_genes)))
        adata = adata[:, keep].copy()
        logger.info(f"Restricted to {adata.n_vars} genes from {restrict_genes}")

    # Gene selection
    hvg_table = None
    if reference_h5ad:
        from rosmap_processing.data.genes import read_var

        ref_var = read_var(import_genes)
        logger.info(f"Aligning to the {len(ref_var)} genes of {import_genes}")
        adata = align_to_reference(adata, ref_var, allow_missing=allow_missing_genes)
    else:
        if import_genes is not None:
            keep = _genes_from_list(adata, load_gene_list(Path(import_genes)))
            adata = adata[:, keep].copy()
            logger.info(f"Shape after gene selection: {adata.shape}")

        n_select, flavor = None, None
        if hvg_after_import is not None:
            n_select, flavor = hvg_after_import, "seurat_v3"
        elif hvg_after_normalize is not None:
            n_select, flavor = hvg_after_normalize, "seurat"
        elif import_genes is None:
            n_select, flavor = n_hvgs, "seurat_v3"

        if n_select is not None:
            hvg_table = _select_hvgs(
                adata, n_select, flavor, hvg_batch_key,
                adata.obs["total_counts"].to_numpy(), target_sum,
            )
            _write_hvg_tables(hvg_table, output_path)
            adata = adata[:, hvg_table["highly_variable"].to_numpy()].copy()
            logger.info(f"Shape after HVG selection: {adata.shape}")

    # Normalization, with totals from all genes
    logger.info(f"Normalizing to {target_sum:g} counts per cell (totals over all genes) and log1p")
    counts = adata.X
    adata.layers["counts"] = to_int32_counts(counts)
    adata.X = log_normalize(counts, adata.obs["total_counts"].to_numpy(), target_sum)
    del counts

    # PCA (global if not individual)
    if not individual_pca:
        logger.info(f"Computing PCA ({n_pca_components} components)")
        sc.pp.pca(adata, n_comps=n_pca_components)
    else:
        logger.info("Skipping global PCA (will compute per donor)")

    adata.uns["provenance_scanpy"] = provenance_record("pipeline scanpy", params)

    # Build KNN graph per donor
    logger.info("Building KNN graphs per donor")
    slices = []
    skipped_donors = 0
    for donor in tqdm(adata.obs[DONOR_ID_COLUMN].unique(), desc="Processing donors"):
        donor_slice = adata[adata.obs[DONOR_ID_COLUMN] == donor].copy()
        if donor_slice.shape[0] <= 1:
            # Donor IDs are restricted information: count only, never log them.
            skipped_donors += 1
            continue
        if individual_pca:
            sc.pp.pca(donor_slice, n_comps=min(n_pca_components, min(donor_slice.shape) - 1))
        sc.pp.neighbors(donor_slice, n_neighbors=k_neighbors, use_rep="X_pca")
        slices.append(donor_slice)

    if skipped_donors:
        logger.warning(f"Skipped {skipped_donors} donors with fewer than 2 cells")
    logger.info(f"Processed {len(slices)} donors successfully")

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

    logger.info(f"Saving processed data to: {output_path}")
    tmp_path = output_path.with_name(output_path.name + ".partial")
    adata.write_h5ad(tmp_path, compression="gzip")
    tmp_path.replace(output_path)

    logger.info("=" * 60)
    logger.info("Processing complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    print("Run via: python -m rosmap_processing pipeline scanpy --help", file=sys.stderr)
    sys.exit(2)
