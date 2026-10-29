"""Compare HVG selections: flavor x batch handling x number of genes.

Loads the counts once, applies the same filtering and gene restriction as
`pipeline scanpy` (without an imported gene list), and for every combination
writes the selected genes in rank order, plus a table of pairwise overlaps.
Extra gene lists (TXT, or the genes of an h5ad) can be added to the
comparison, e.g. earlier outputs. The outputs hold gene names only.

Outputs in the output directory:
    hvg_<flavor>_<donor|pooled>_<n>.txt   selected genes, ranked
    hvg_overlap.tsv                       every pair: sizes, intersection, Jaccard
    hvg_overlap_summary.tsv               per flavor and n: per-donor vs pooled
"""

from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import anndata as ad
import pandas as pd
import scanpy as sc

from ..core.scanpy_pipeline import (
    _attach_gene_ids,
    _ensure_donor_column,
    _genes_from_list,
    _select_counts,
    _select_hvgs,
    load_gene_list,
    rank_hvg_table,
)
from ..utils.constants import (
    DONOR_ID_COLUMN,
    MIN_CELLS_PER_DONOR,
    MIN_CELLS_PER_GENE,
    MIN_GENES_PER_CELL,
    NORMALIZE_TARGET_SUM,
)
from ..utils.counts import log_normalize, total_counts
from ..utils.logging import get_logger

logger = get_logger(__name__)


def variant_name(flavor: str, batch_key: Optional[str], n: int) -> str:
    return f"{flavor}_{'donor' if batch_key else 'pooled'}_{n}"


def _read_list(path: Path) -> List[str]:
    """Genes of a TXT list (one per line) or of an h5ad (its var_names)."""
    if path.suffix == ".h5ad":
        from ..data.genes import read_var_names

        return read_var_names(path)
    return load_gene_list(path)


def overlap_table(lists: Dict[str, Sequence[str]]) -> pd.DataFrame:
    """Pairwise sizes, intersection and Jaccard index of gene lists."""
    rows = []
    for a, b in combinations(lists, 2):
        sa, sb = set(lists[a]), set(lists[b])
        inter = len(sa & sb)
        rows.append({
            "a": a, "b": b, "n_a": len(sa), "n_b": len(sb), "intersection": inter,
            "jaccard": inter / len(sa | sb) if sa | sb else float("nan"),
            "frac_of_a": inter / len(sa) if sa else float("nan"),
            "frac_of_b": inter / len(sb) if sb else float("nan"),
        })
    return pd.DataFrame(rows)


def compare_hvgs(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    n_genes: Sequence[int] = (1000, 2000, 3000),
    flavors: Sequence[str] = ("seurat_v3", "seurat"),
    batch_keys: Sequence[Optional[str]] = (DONOR_ID_COLUMN, None),
    layer: Optional[str] = None,
    gene_map: Optional[Union[str, Path]] = None,
    restrict_genes: Optional[Union[str, Path]] = None,
    donor_column: Optional[str] = None,
    extra_lists: Optional[Dict[str, Union[str, Path]]] = None,
    min_genes: int = MIN_GENES_PER_CELL,
    min_cells: int = MIN_CELLS_PER_GENE,
    min_cells_per_donor: int = MIN_CELLS_PER_DONOR,
    target_sum: float = NORMALIZE_TARGET_SUM,
) -> pd.DataFrame:
    """Run every HVG variant on one dataset and write lists + overlap tables."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    batch_keys = [None if b in (None, "", "none", "None") else b for b in batch_keys]

    logger.info(f"Loading {input_path}")
    adata = _select_counts(ad.read_h5ad(input_path), layer, use_raw=False)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = _ensure_donor_column(adata, donor_column, min_cells_per_donor)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    size_factors = total_counts(adata.X)
    if gene_map is not None:
        _attach_gene_ids(adata, gene_map)
    if restrict_genes is not None:
        adata = adata[:, _genes_from_list(adata, load_gene_list(Path(restrict_genes)))].copy()
    logger.info(f"Candidate genes: {adata.n_vars}, cells: {adata.n_obs}")

    normalized = None
    if "seurat" in flavors:
        normalized = log_normalize(adata.X, size_factors, target_sum)

    lists: Dict[str, List[str]] = {}
    for flavor in flavors:
        for batch_key in batch_keys:
            for n in n_genes:
                name = variant_name(flavor, batch_key, n)
                table = _select_hvgs(
                    adata, n, flavor, batch_key, size_factors, target_sum, normalized=normalized
                )
                ranked = rank_hvg_table(table)
                selected = ranked.index[ranked["highly_variable"].to_numpy()].tolist()
                (output_dir / f"hvg_{name}.txt").write_text("\n".join(selected) + "\n")
                lists[name] = selected
                logger.info(f"{name}: {len(selected)} genes")

    for label, path in (extra_lists or {}).items():
        path = Path(path)
        if not path.exists():
            logger.warning(f"Comparison list {label} not found: {path} (skipped)")
            continue
        lists[label] = _read_list(path)
        logger.info(f"Comparison list {label}: {len(lists[label])} genes from {path}")

    overlaps = overlap_table(lists)
    overlaps.to_csv(output_dir / "hvg_overlap.tsv", sep="\t", index=False, float_format="%.4f")

    summary = []
    if len(batch_keys) == 2 and None in batch_keys:
        donor_key = next(b for b in batch_keys if b)
        for flavor in flavors:
            for n in n_genes:
                a, b = variant_name(flavor, donor_key, n), variant_name(flavor, None, n)
                row = overlaps[(overlaps["a"] == a) & (overlaps["b"] == b)].iloc[0]
                summary.append({"flavor": flavor, "n": n, "intersection": row["intersection"],
                                "jaccard": round(row["jaccard"], 4)})
                logger.info(
                    f"{flavor}, {n} genes: per-donor vs pooled share {row['intersection']} "
                    f"({100 * row['intersection'] / n:.1f}%), Jaccard {row['jaccard']:.3f}"
                )
    summary = pd.DataFrame(summary)
    summary.to_csv(output_dir / "hvg_overlap_summary.tsv", sep="\t", index=False)
    logger.info(f"Wrote {len(lists)} lists and overlap tables to {output_dir}")
    return overlaps
