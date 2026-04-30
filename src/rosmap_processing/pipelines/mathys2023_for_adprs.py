"""Aggregate Mathys 2023 per-class h5ads into per-donor x cell-type expression stats.

Produces inputs for the downstream ad-prs project: a long-format parquet with
per-(donor, cell_type, gene) statistics at both major and subtype granularity,
plus a 1-row-per-donor metadata sidecar and a README. No HVG filter, no donor
filter, no cross-donor pooling — those decisions live downstream.

See the README written by `write_readme()` for the output schema.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from tqdm import tqdm

from ..utils.logging import get_logger

logger = get_logger(__name__)


# Per-class file (basename without .h5ad) -> Yang's 6-class major label.
# Mathys 2023's "Immune" compartment is collapsed to "Mic" — see README sec. on
# cell-type harmonization. Vasc is intentionally absent from this map (Yang/O'Neill
# drop endothelial/pericytes). If the local data ever gains a Vasc h5ad, it will
# be processed at subtype level only as long as it isn't added here.
FILE_TO_MAJOR: Dict[str, str] = {
    "Astrocytes": "Ast",
    "Excitatory_neurons_set1": "Ex",
    "Excitatory_neurons_set2": "Ex",
    "Excitatory_neurons_set3": "Ex",
    "Inhibitory_neurons": "In",
    "Immune_cells": "Mic",
    "Oligodendrocytes": "Oli",
    "OPCs": "Opc",
}

EXPECTED_MAJORS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]

NORMALIZATION_TARGET_SUM = 1e4

DONOR_ID_COL = "projid"
SUBTYPE_COL = "cell_type_high_resolution"

OUTPUT_PARQUET_NAME = "per_donor_celltype_expression.parquet"
DONOR_METADATA_NAME = "donor_metadata.tsv"
README_NAME = "README.md"
INTERMEDIATE_DIR = "_intermediate"


def _resolve_input_files(
    input_dir: Path, file_to_major: Dict[str, str]
) -> Dict[Path, str]:
    """Return {path: major_label} for files that exist on disk.

    Looks for `<basename>.h5ad` case-insensitively (the existing R conversion
    produces both `Astrocytes.h5ad` and `astrocytes.h5ad` depending on which
    branch of run_all.sh ran).
    """
    resolved: Dict[Path, str] = {}
    missing: List[str] = []
    for stem, major in file_to_major.items():
        candidates = [input_dir / f"{stem}.h5ad", input_dir / f"{stem.lower()}.h5ad"]
        match = next((p for p in candidates if p.exists()), None)
        if match is None:
            missing.append(stem)
        else:
            resolved[match] = major
    if missing:
        logger.warning(
            f"{len(missing)} expected per-class h5ad(s) not found and will be "
            f"skipped: {missing}. Run the existing convert/metadata pipeline first."
        )
    if not resolved:
        raise FileNotFoundError(
            f"No per-class h5ads found under {input_dir}. Expected any of "
            f"{list(file_to_major.keys())}."
        )
    return resolved


def _normalize_lognormed(X_raw: sp.csr_matrix) -> sp.csr_matrix:
    """Compute normalize_total(target_sum=1e4) + log1p on a sparse copy.

    Equivalent to scanpy's recipe but bypasses building an intermediate
    AnnData. Operates row-wise, in-place on a copy of X_raw.
    """
    X = X_raw.copy().astype(np.float32, copy=False)
    if not sp.isspmatrix_csr(X):
        X = X.tocsr()
    row_sums = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    scale = (NORMALIZATION_TARGET_SUM / row_sums).astype(np.float32)
    # Multiply each row's data by its scale factor.
    nnz_per_row = np.diff(X.indptr)
    X.data *= np.repeat(scale, nnz_per_row)
    np.log1p(X.data, out=X.data)
    return X


def _aggregate_groups(
    X_raw: sp.csr_matrix,
    X_lognorm: sp.csr_matrix,
    group_codes: np.ndarray,
    n_groups: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Group-by aggregation via sparse indicator matmul.

    Returns
    -------
    mean_lognorm : (n_groups, n_genes) float32 — mean of log-normalized expr
    frac_expressed : (n_groups, n_genes) float32 — fraction of cells with raw count > 0
    n_cells : (n_groups,) int64
    """
    n_cells = X_raw.shape[0]
    indicator = sp.csr_matrix(
        (
            np.ones(n_cells, dtype=np.float32),
            (group_codes.astype(np.int32), np.arange(n_cells, dtype=np.int32)),
        ),
        shape=(n_groups, n_cells),
    )
    n_cells_per_group = np.asarray(indicator.sum(axis=1)).ravel().astype(np.int64)

    # frac_expressed: count nonzero raw values per (group, gene), divide by group size.
    nz = X_raw.copy()
    nz.data = (nz.data > 0).astype(np.float32)
    counts_nz = indicator @ nz
    if sp.issparse(counts_nz):
        counts_nz = counts_nz.toarray()
    frac = counts_nz.astype(np.float32) / n_cells_per_group[:, None].astype(np.float32)

    # mean of log-normalized.
    sum_log = indicator @ X_lognorm
    if sp.issparse(sum_log):
        sum_log = sum_log.toarray()
    mean_log = sum_log.astype(np.float32) / n_cells_per_group[:, None].astype(np.float32)

    return mean_log, frac, n_cells_per_group


def _materialize_long(
    mean_log: np.ndarray,
    frac: np.ndarray,
    n_cells: np.ndarray,
    group_donors: np.ndarray,
    group_cell_types: np.ndarray,
    granularity: str,
    gene_symbols: np.ndarray,
) -> pd.DataFrame:
    """Per-group dense matrices -> long-format DataFrame (one row per gene-group)."""
    n_groups, n_genes = mean_log.shape
    return pd.DataFrame(
        {
            "donor_id": np.repeat(group_donors, n_genes),
            "granularity": granularity,
            "cell_type": np.repeat(group_cell_types, n_genes),
            "gene_symbol": np.tile(gene_symbols, n_groups),
            "mean_expr": mean_log.ravel().astype(np.float32),
            "frac_expressed": frac.ravel().astype(np.float32),
            "n_cells": np.repeat(n_cells.astype(np.int32), n_genes),
        }
    )


def _write_chunk(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df["donor_id"] = df["donor_id"].astype(str)
    df["granularity"] = df["granularity"].astype("category")
    df["cell_type"] = df["cell_type"].astype("category")
    df["gene_symbol"] = df["gene_symbol"].astype("category")
    df.to_parquet(path, compression="zstd", index=False)
    logger.info(
        f"  wrote {path.name} ({path.stat().st_size / 1e6:.1f} MB, {len(df)} rows)"
    )


def process_class_file(
    h5ad_path: Path,
    major_label: str,
    output_dir: Path,
) -> Tuple[Path, Path, pd.DataFrame]:
    """Aggregate one per-class h5ad and write its major + subtype parquet chunks.

    Returns
    -------
    major_chunk_path, subtype_chunk_path, per_donor_cell_counts_df
    """
    logger.info(f"Processing {h5ad_path.name} (major={major_label})")
    adata = ad.read_h5ad(h5ad_path)

    if DONOR_ID_COL not in adata.obs.columns:
        raise ValueError(
            f"{h5ad_path}: missing '{DONOR_ID_COL}' in obs. Did metadata-add run?"
        )
    if SUBTYPE_COL not in adata.obs.columns:
        raise ValueError(
            f"{h5ad_path}: missing '{SUBTYPE_COL}' in obs. Did category-fix run?"
        )

    donor_ids = adata.obs[DONOR_ID_COL].astype(str).values
    subtypes = adata.obs[SUBTYPE_COL].astype(str).values
    gene_symbols = adata.var_names.astype(str).values

    X = adata.X
    if not sp.isspmatrix_csr(X):
        X = sp.csr_matrix(X)
    X_raw = X.astype(np.float32)

    # Sanity: warn if X doesn't look like raw counts.
    if X_raw.nnz:
        sample = X_raw.data[: min(X_raw.nnz, 1000)]
        whole = float(np.mean(sample == np.floor(sample)))
        if whole < 0.99:
            logger.warning(
                f"{h5ad_path.name}: only {whole * 100:.1f}% of nonzero values are "
                f"integers — input may already be normalized. Aggregates assume "
                f"raw counts in X."
            )

    X_lognorm = _normalize_lognormed(X_raw)

    # Major-level: group by donor only (major label is fixed per file).
    donor_codes, donor_uniques = pd.factorize(donor_ids)
    donor_uniques = np.asarray(donor_uniques, dtype=str)
    n_donors = len(donor_uniques)
    logger.info(
        f"  major: {n_donors} donors over {adata.n_obs} cells x {adata.n_vars} genes"
    )
    mean_major, frac_major, n_cells_major = _aggregate_groups(
        X_raw, X_lognorm, donor_codes, n_donors
    )
    major_long = _materialize_long(
        mean_major,
        frac_major,
        n_cells_major,
        group_donors=donor_uniques,
        group_cell_types=np.full(n_donors, major_label),
        granularity="major",
        gene_symbols=gene_symbols,
    )

    # Subtype-level: group by (donor, cell_type_high_resolution).
    pair_index = pd.MultiIndex.from_arrays([donor_ids, subtypes])
    pair_codes, pair_uniques = pd.factorize(pair_index)
    n_pairs = len(pair_uniques)
    pair_donors = np.array([t[0] for t in pair_uniques])
    pair_cts = np.array([t[1] for t in pair_uniques])
    logger.info(f"  subtype: {n_pairs} (donor, subtype) groups")
    mean_sub, frac_sub, n_cells_sub = _aggregate_groups(
        X_raw, X_lognorm, pair_codes, n_pairs
    )
    subtype_long = _materialize_long(
        mean_sub,
        frac_sub,
        n_cells_sub,
        group_donors=pair_donors,
        group_cell_types=pair_cts,
        granularity="subtype",
        gene_symbols=gene_symbols,
    )

    out_dir = output_dir / INTERMEDIATE_DIR
    major_path = out_dir / f"{h5ad_path.stem}_major.parquet"
    subtype_path = out_dir / f"{h5ad_path.stem}_subtype.parquet"
    _write_chunk(major_long, major_path)
    _write_chunk(subtype_long, subtype_path)

    cell_counts = pd.DataFrame(
        {
            "donor_id": donor_uniques,
            f"n_cells_{major_label}": n_cells_major.astype(np.int64),
        }
    )
    return major_path, subtype_path, cell_counts


def merge_chunks(chunk_paths: List[Path], output_path: Path) -> None:
    """Concatenate intermediate parquet chunks into a single file."""
    import pyarrow.parquet as pq

    if not chunk_paths:
        raise ValueError("No chunks to merge")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    first = pq.read_table(str(chunk_paths[0]))
    schema = first.schema
    writer = pq.ParquetWriter(str(output_path), schema, compression="zstd")
    writer.write_table(first)
    for p in chunk_paths[1:]:
        t = pq.read_table(str(p))
        writer.write_table(t.cast(schema))
    writer.close()
    logger.info(
        f"Merged {len(chunk_paths)} chunks -> {output_path} "
        f"({output_path.stat().st_size / 1e6:.1f} MB)"
    )


def build_donor_metadata(
    cell_counts_per_class: List[pd.DataFrame],
    clinical_csv: Path,
    mit_metadata_csv: Optional[Path],
) -> pd.DataFrame:
    """Build per-donor metadata: cell counts + clinical fields, one row per donor."""
    all_donors = set()
    for df in cell_counts_per_class:
        all_donors.update(df["donor_id"].astype(str).tolist())
    donors_df = pd.DataFrame({"donor_id": sorted(all_donors)})

    # Sum across multiple class files that map to the same major label
    # (e.g. Excitatory_neurons_set{1,2,3} all -> n_cells_Ex).
    by_major: Dict[str, pd.DataFrame] = {}
    for df in cell_counts_per_class:
        df_str = df.copy()
        df_str["donor_id"] = df_str["donor_id"].astype(str)
        major_col = next(c for c in df_str.columns if c.startswith("n_cells_"))
        if major_col in by_major:
            by_major[major_col] = (
                pd.concat([by_major[major_col], df_str], ignore_index=True)
                .groupby("donor_id", as_index=False)[major_col]
                .sum()
            )
        else:
            by_major[major_col] = df_str

    for col, df in by_major.items():
        donors_df = donors_df.merge(df, on="donor_id", how="left")

    for m in EXPECTED_MAJORS:
        col = f"n_cells_{m}"
        if col not in donors_df.columns:
            donors_df[col] = 0
        donors_df[col] = donors_df[col].fillna(0).astype(np.int64)

    donors_df["total_cells"] = donors_df[
        [f"n_cells_{m}" for m in EXPECTED_MAJORS]
    ].sum(axis=1)

    clinical = pd.read_csv(clinical_csv, dtype={"projid": str})
    keep_rename = {
        "projid": "donor_id",
        "msex": "sex",
        "age_death": "age_at_death",
        "cts_mmse30_lv": "MMSE",
        "braaksc": "Braak",
        "ceradsc": "CERAD",
        "cogdx": "cogdx",
        "apoe_genotype": "apoe_genotype",
        "individualID": "individualID",
    }
    keep_cols = [c for c in keep_rename if c in clinical.columns]
    clinical = clinical[keep_cols].rename(columns=keep_rename)
    donors_df = donors_df.merge(clinical, on="donor_id", how="left")

    if mit_metadata_csv is not None and mit_metadata_csv.exists():
        mit_meta = pd.read_csv(mit_metadata_csv)
        if "individualID" in mit_meta.columns and "subject" in mit_meta.columns:
            mit_meta = mit_meta[["individualID", "subject"]].drop_duplicates(
                subset=["individualID"]
            )
            donors_df = donors_df.merge(mit_meta, on="individualID", how="left")

    cell_count_cols = ["total_cells"] + [f"n_cells_{m}" for m in EXPECTED_MAJORS]
    other = [c for c in donors_df.columns if c not in cell_count_cols and c != "donor_id"]
    donors_df = donors_df[["donor_id", *other, *cell_count_cols]]
    return donors_df


def write_readme(
    output_dir: Path,
    files_processed: Dict[Path, str],
    n_donors: int,
    n_genes: int,
    total_cells: int,
) -> None:
    files_block = "\n".join(
        f"- {p.name} -> {major}" for p, major in files_processed.items()
    )
    text = f"""# Mathys 2023 per-donor x cell-type aggregates for ad-prs

Output produced by `python -m rosmap_processing pipeline mathys2023-for-adprs`
from this repo (rosmap-processing). Consumer: ad-prs.

## Files

- `{OUTPUT_PARQUET_NAME}` — long-format per-(donor, cell_type, gene) stats
- `{DONOR_METADATA_NAME}` — one row per donor; clinical fields + per-major cell counts
- `{README_NAME}` — this file

## Source

ROSMAP / MIT-lineage snRNA-seq atlas (Mathys et al., Cell 2023). Input was
the per-cell-class h5ads under `data/raw/ROSMAP_MIT/`, after the existing
`run_all.sh` conversion + metadata + category-fix pipeline. Synapse IDs are
in `src/rosmap_processing/utils/constants.py:SYNAPSE_IDS_ROSMAP_MIT`.

Per-class files processed and their major-class harmonization:

{files_block}

The Vasc compartment is **not** included (not present in the local download).

## Normalization recipe

Per cell, raw counts are normalized to a fixed library size of {NORMALIZATION_TARGET_SUM:.0e}
and log1p-transformed:

```python
sc.pp.normalize_total(adata, target_sum={NORMALIZATION_TARGET_SUM:.0e})
sc.pp.log1p(adata)
```

`mean_expr` is the mean of these log-normalized values. `frac_expressed` is
computed from the **raw** counts (fraction of cells in the bucket with
count > 0 for that gene). `n_cells` is the number of cells in the bucket.

## Per-(donor, cell_type, gene) parquet schema

| column | dtype | notes |
|---|---|---|
| donor_id | str (categorical) | `projid` from the AnnData |
| granularity | category | "major" or "subtype" |
| cell_type | category | major: Ex/In/Ast/Mic/Oli/Opc; subtype: verbatim `cell_type_high_resolution` |
| gene_symbol | category | verbatim `var_names` from the AnnData (HGNC symbol from the 10x reference; ~33,538 features in the GRCh38 / GENCODE v32 build that 10x uses) |
| mean_expr | float32 | mean log-normalized expression |
| frac_expressed | float32 | fraction with raw count > 0 |
| n_cells | int32 | cells in (donor, cell_type) bucket |

Empty buckets are not emitted. To get per-major counts including zeros, use
`{DONOR_METADATA_NAME}`.

## Cell-type harmonization

Major labels follow Yang 2023's 6-class scheme:

- Astrocytes -> Ast
- Excitatory_neurons_set{{1,2,3}} -> Ex (three input files merged into one major bucket)
- Inhibitory_neurons -> In
- Immune_cells -> Mic — **note**: Mathys 2023 grouped microglia + macrophages + T cells into one "Immune" compartment. We pass through the *full* Immune compartment as Mic. If you need microglia-only, filter on the subtype-level rows where `cell_type` starts with `Mic.`.
- Oligodendrocytes -> Oli
- OPCs -> Opc

Subtype labels are the verbatim `cell_type_high_resolution` strings from the
atlas (e.g. `Ast.2`, `Mic.1`, `Exc.RELN-LAMP5`, etc.). No subtype -> major
remap is encoded in this output beyond the file partition.

## Gene IDs

The atlas h5ads only carry HGNC symbols in `var_names` (the `gene_id` column
was dropped during the RDS -> h5Seurat -> h5ad conversion). We pass symbols
through verbatim. **Resolve to Ensembl in the consumer (ad-prs).** The
underlying 10x reference is GRCh38 / GENCODE v32; ad-prs should use a
matching GTF for the +/-30 kb windowing step.

## Donor selection

All donors present in the input are emitted, with no pathology filter.
Mathys 2019 reference-control flagging (`is_yang2023_reference_control`,
`in_mathys2019`, `pathology_group_mathys2019`) is **not done here** — the
deidentified -> ROSMAP ID mapping needed to flag the original 24 controls
was deferred. Apply downstream using a Mathys 2019 donor list.

## Sanity numbers from this run

- donors: {n_donors}
- genes: {n_genes}
- total cells: {total_cells:,}
"""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / README_NAME).write_text(text)
    logger.info(f"Wrote {output_dir / README_NAME}")


def run_pipeline(
    input_dir: Path,
    output_dir: Path,
    clinical_csv: Path,
    mit_metadata_csv: Optional[Path] = None,
    file_to_major: Optional[Dict[str, str]] = None,
    keep_intermediate: bool = False,
) -> None:
    """End-to-end aggregation: per-class h5ads -> long parquet + donor metadata + README."""
    file_to_major = file_to_major or FILE_TO_MAJOR
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    clinical_csv = Path(clinical_csv)
    mit_metadata_csv = Path(mit_metadata_csv) if mit_metadata_csv else None

    logger.info("=" * 60)
    logger.info("Mathys 2023 -> ad-prs aggregation")
    logger.info("=" * 60)
    logger.info(f"Input dir: {input_dir}")
    logger.info(f"Output dir: {output_dir}")
    logger.info(f"Clinical CSV: {clinical_csv}")
    logger.info(f"MIT metadata CSV: {mit_metadata_csv}")

    output_dir.mkdir(parents=True, exist_ok=True)

    files = _resolve_input_files(input_dir, file_to_major)
    logger.info(f"Will process {len(files)} per-class h5ad(s)")

    major_chunks: List[Path] = []
    subtype_chunks: List[Path] = []
    cell_counts: List[pd.DataFrame] = []

    n_genes_seen: Optional[int] = None
    for path, major in tqdm(files.items(), desc="Per-class aggregation"):
        major_path, subtype_path, cc = process_class_file(path, major, output_dir)
        major_chunks.append(major_path)
        subtype_chunks.append(subtype_path)
        cell_counts.append(cc)

        n_genes_chunk = pd.read_parquet(major_path, columns=["gene_symbol"])[
            "gene_symbol"
        ].nunique()
        if n_genes_seen is None:
            n_genes_seen = n_genes_chunk
        elif n_genes_chunk != n_genes_seen:
            logger.warning(
                f"{path.name}: gene count {n_genes_chunk} differs from prior "
                f"{n_genes_seen}. The merged parquet will be a union; "
                f"per-(donor, cell_type) gene rows may be inconsistent across "
                f"the two granularities. Verify input."
            )

    out_parquet = output_dir / OUTPUT_PARQUET_NAME
    merge_chunks(major_chunks + subtype_chunks, out_parquet)

    donor_meta = build_donor_metadata(cell_counts, clinical_csv, mit_metadata_csv)
    out_meta = output_dir / DONOR_METADATA_NAME
    donor_meta.to_csv(out_meta, sep="\t", index=False)
    logger.info(f"Wrote {out_meta} ({len(donor_meta)} donors)")

    n_genes = int(n_genes_seen or 0)
    total_cells = int(donor_meta["total_cells"].sum())
    write_readme(
        output_dir=output_dir,
        files_processed=files,
        n_donors=len(donor_meta),
        n_genes=n_genes,
        total_cells=total_cells,
    )

    if not keep_intermediate:
        intermediate = output_dir / INTERMEDIATE_DIR
        for p in major_chunks + subtype_chunks:
            p.unlink(missing_ok=True)
        if intermediate.exists() and not any(intermediate.iterdir()):
            intermediate.rmdir()

    summary = {
        "n_donors": int(len(donor_meta)),
        "n_genes": n_genes,
        "total_cells": total_cells,
        "files_processed": [p.name for p in files],
        "output_parquet": str(out_parquet),
        "donor_metadata": str(out_meta),
    }
    logger.info("Summary: " + json.dumps(summary))
    logger.info("=" * 60)
    logger.info("Done.")
    logger.info("=" * 60)


def main() -> None:
    import argparse

    from ..utils.logging import setup_logging

    parser = argparse.ArgumentParser(
        description=(
            "Aggregate per-class Mathys 2023 h5ads into per-(donor, cell_type, "
            "gene) statistics for ad-prs."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/raw/ROSMAP_MIT"),
        help="Directory containing per-class h5ad files (default: data/raw/ROSMAP_MIT)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/mathys2023_for_adprs"),
        help="Directory to write outputs into (default: data/processed/mathys2023_for_adprs)",
    )
    parser.add_argument(
        "--clinical-csv",
        type=Path,
        default=Path("data/raw/ROSMAP/rosmap_clinical.csv"),
        help="ROSMAP clinical CSV (default: data/raw/ROSMAP/rosmap_clinical.csv)",
    )
    parser.add_argument(
        "--mit-metadata-csv",
        type=Path,
        default=Path("data/raw/ROSMAP_MIT/MIT_ROSMAP_Multiomics_individual_metadata.csv"),
        help="MIT individual metadata CSV (for individualID -> subject mapping)",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep per-class parquet chunks under _intermediate/ after merging",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    args = parser.parse_args()

    setup_logging(level=args.log_level)
    run_pipeline(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        clinical_csv=args.clinical_csv,
        mit_metadata_csv=args.mit_metadata_csv,
        keep_intermediate=args.keep_intermediate,
    )


if __name__ == "__main__":
    main()
