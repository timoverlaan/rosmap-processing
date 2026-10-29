# ROSMAP Processing CLI Reference

All commands run as `python -m rosmap_processing <group> <command>` (inside the
container: `pixi run python -m rosmap_processing ...`). Add `-v` for debug
logging; every command has `--help`.

Commands that write files never modify their input, and refuse to replace an
existing output unless `--overwrite` is given.

## data

### `data download`
```bash
python -m rosmap_processing data download [--dataset all|rosmap|rosmap-mit] \
    [--output-dir data/raw] [--include-raw-h5ad]
```
Downloads the pinned versions in `SYNAPSE_FILES` (`utils/constants.py`) into
`<output-dir>/ROSMAP` and `<output-dir>/ROSMAP_MIT`, and checks each file's size
and MD5. A file that is already present is kept if it matches and is an error
if it does not (it is never overwritten). The 103 GB unfiltered ROSMAP-MIT h5ad
is only fetched with `--include-raw-h5ad`. Needs `SYNAPSE_AUTH_TOKEN` or
`token.txt` only when something has to be downloaded.

### `data verify`
```bash
python -m rosmap_processing data verify --dataset rosmap|rosmap-mit --dir <directory>
```
Checks size and MD5 of the pinned files in a directory; exits non-zero on any
missing or mismatched file. Needs no Synapse login.

### `data prepare`
```bash
python -m rosmap_processing data prepare <input.h5Seurat|converted.h5ad> --output <counts.h5ad> \
    [--columns cell_type_high_resolution ...] [--metadata data/v2/download/ROSMAP/ROSMAP_clinical.csv]
```
Takes an `.h5Seurat` file (the RNA `counts` slot and `meta.data` are read
directly, R factors become categoricals) or a SeuratDisk-converted h5ad
(X = Seurat `data` slot, `.raw` = counts),
puts the counts in X as int32, checks that every stored value is a
non-negative integer (fails otherwise), makes the given obs columns
categorical, and merges clinical metadata on `projid` (or `individualID`).

### `data category-fix`
```bash
python -m rosmap_processing data category-fix <input.h5ad> [--columns ...] [--remove-raw]
```
Makes obs columns categorical, in place. `--remove-raw` is refused when X does
not hold counts, because `.raw` is then the only copy of them.

### `data metadata`
```bash
python -m rosmap_processing data metadata <input.h5ad> <metadata.csv> [--output <out.h5ad>]
```
Merges a clinical CSV into obs. Without `--output` the input is modified in place.

## core

### `core combine`
```bash
python -m rosmap_processing core combine <a.h5ad> <b.h5ad> ... --output <combined.h5ad>
```
Concatenates cells (outer join on genes; a gene missing from one input is zero
in its cells).

### `core column-mapping`
```bash
python -m rosmap_processing core column-mapping <input.h5ad> --data-type SeaAD --output <out.h5ad>
```
Adds ROSMAP-style columns to SeaAD (`braaksc`, `ceradsc`, `thalsc`, Wang labels).
Also supports `ROSMAP` and `ROSMAP_MIT` renames (not used by the pipeline).

## genes

### `genes build-map`
```bash
python -m rosmap_processing genes build-map --gtf <annotation.gtf.gz> \
    --recipe cellranger-GRCh38-3.0.0|cellranger-GRCh38-2020-A --h5ad <data.h5ad> --output <map.tsv>
```
Rebuilds the 10x reference's gene list from the GTF it was built from and maps
each gene symbol of the h5ad to an Ensembl ID. Output TSV columns: `var_name`,
`gene_id`, `ref_name`, `method` (`position` when the rebuilt list matches the
data exactly, else `name`, `unique_symbol`, `unmatched` or `ambiguous`).

| Dataset | Reference | GTF |
|---|---|---|
| ROSMAP-MIT (Mathys 2023), 33,538 genes | `cellranger-GRCh38-3.0.0` | Ensembl 93 |
| SeaAD, Green 2024 ROSMAP, 36,601 genes | `cellranger-GRCh38-2020-A` | GENCODE v32 |

### `genes shared`
```bash
python -m rosmap_processing genes shared <map_a.tsv> <map_b.tsv> [<map_c.tsv> ...] --output <shared_ids.txt>
```
Ensembl IDs mapped in every given map, in the first map's order.

## pipeline

### `pipeline scanpy`
```bash
python -m rosmap_processing pipeline scanpy <counts.h5ad> <output.h5ad> [options]
```
1. Takes the counts from X (or `--layer`, `--raw`) and checks they are integers.
2. Filters cells (`--min-genes`, default 200), then drops cells without a
   donor (missing, or a placeholder such as `NA`) and donors with fewer than
   `--min-cells-per-donor` cells (default 300), then, unless a gene list is
   imported, filters genes (`--min-cells`, default 5).
3. Computes per-cell totals over all genes, used for normalization.
4. Chooses genes (below), then sets X = log1p(counts / total × `--target-sum`)
   (default 1e4) and keeps the counts in `layers['counts']`.
5. PCA (`--n-pcs`, default 50; `--individual-pca` for per donor) and a kNN graph
   per donor (`--k-neighbors`, default 30).

Gene choice:

| Options | Genes |
|---|---|
| (none) | `--n-genes` HVGs, seurat_v3 on counts |
| `--hvg-after-normalize N` | N HVGs, seurat flavor on log-normalized values |
| `--import-genes list.txt` | the listed genes (symbols or Ensembl IDs) |
| `--import-genes list.txt --hvg-after-import N` | N HVGs (seurat_v3) within the list |
| `--import-genes list.txt --hvg-after-normalize N` | N HVGs (seurat) within the list |
| `--import-genes reference.h5ad` | exactly the reference's genes, in its order; fails on missing genes unless `--allow-missing-genes` (zero-fills) |

Other options: `--hvg-batch-key` (default `Donor ID`, `none` to pool),
`--gene-map map.tsv` (adds `var['gene_id']`, used for matching),
`--restrict-genes list.txt` (all choices are made within this list, e.g. the
output of `genes shared`), `--donor-column COL` (obs column with the donor,
copied to `Donor ID`; default `Donor ID`, else `projid`; cells without a donor
are dropped), `--overwrite`.

Writes `<output>_hvg_scores.tsv` and `<output>_hvg_names.txt` when HVGs are
selected, and records parameters and the git commit in `uns['provenance_scanpy']`.

### `pipeline hvg-compare`
```bash
python -m rosmap_processing pipeline hvg-compare <counts.h5ad> --output-dir <dir> \
    [--n-genes 1000 2000 3000] [--flavors seurat_v3 seurat] [--batch-keys "Donor ID" none] \
    [--gene-map map.tsv] [--restrict-genes list.txt] [--compare NAME=PATH ...]
```
Loads the counts once and runs every combination of flavor, batch handling and
size, with the same filtering and normalization as `pipeline scanpy`. Writes
`hvg_<flavor>_<donor|pooled>_<n>.txt` (selected genes, ranked),
`hvg_overlap.tsv` (every pair of lists: sizes, intersection, Jaccard) and
`hvg_overlap_summary.tsv` (per-donor vs pooled per flavor and size). `--compare`
adds other lists (TXT, or an h5ad's genes). Outputs hold gene names only.

### `pipeline mathys2023-for-adprs`
```bash
python -m rosmap_processing pipeline mathys2023-for-adprs \
    [--input-dir data/v2/counts/ROSMAP_MIT] [--output-dir data/v2/processed/mathys2023_for_adprs] \
    [--gene-map data/v2/reference/rosmap_mit_gene_map.tsv] [--clinical-csv ...] [--mit-metadata-csv ...]
```
Per-(donor, cell type, gene) mean log-normalized expression and fraction
expressed, for the ad-prs project. Fails if any input does not hold raw counts.
See the README it writes into the output directory.

## utils

### `utils check-counts`
```bash
python -m rosmap_processing utils check-counts <a.h5ad> ... [--key X|raw/X|layers/<name>]
```
Reads the stored matrix in chunks and reports only aggregate statistics (shape,
dtype, number of non-integer and negative values, min, max). Exits non-zero if
any file does not hold raw counts. Safe for restricted data.

### `utils validate`
```bash
python -m rosmap_processing utils validate <file.h5ad>
```
Prints structure, obs value counts and a 20×20 preview of X. **The preview shows
data values: do not use it on restricted data where logs are shared.**
