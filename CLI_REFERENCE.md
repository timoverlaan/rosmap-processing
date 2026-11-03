# ROSMAP Processing CLI Reference

All commands run as `python -m rosmap_processing <group> <command>` (inside the
container: `pixi run python -m rosmap_processing ...`). Add `-v` for debug
logging; every command has `--help`.

Commands that write files never modify their input, and refuse to replace an
existing output unless `--overwrite` is given.

## data

### `data download`
```bash
python -m rosmap_processing data download [--dataset all|rosmap|rosmap-mit|psychad|diverse-cohorts] \
    [--output-dir data/raw] [--include-raw-h5ad]
```
Downloads the pinned versions in `SYNAPSE_FILES` (`utils/constants.py`) into
`<output-dir>/ROSMAP`, `ROSMAP_MIT`, `PsychAD` and `DiverseCohorts`, and checks
each file's size and MD5. A file that is already present is kept if it matches
and is an error if it does not (it is never overwritten). The unfiltered h5ad
files (ROSMAP-MIT 103 GB, PsychAD before QC 106 GB) are only fetched with
`--include-raw-h5ad`. Needs `SYNAPSE_AUTH_TOKEN` or `token.txt` only when
something has to be downloaded. On DAIC: `sbatch slurm/download.sh <dataset>...`.

### `data verify`
```bash
python -m rosmap_processing data verify --dataset rosmap|rosmap-mit|psychad|diverse-cohorts --dir <directory>
```
Checks size and MD5 of the pinned files in a directory; exits non-zero on any
missing or mismatched file. Needs no Synapse login.

### `data prepare`
```bash
python -m rosmap_processing data prepare <input.h5Seurat|converted.h5ad|matrix.h5> --output <counts.h5ad> \
    [--columns cell_type_high_resolution ...] [--metadata data/v2/download/ROSMAP/ROSMAP_clinical.csv]
```
Takes an `.h5Seurat` file (the RNA `counts` slot and `meta.data` are read
directly, R factors become categoricals), a SeuratDisk-converted h5ad
(X = Seurat `data` slot, `.raw` = counts), or the `.h5` that
`src/convert_R/rdata_to_h5.R` writes from the sparse matrix in an `.rds`/`.rda`,
puts the counts in X as int32, checks that every stored value is a
non-negative integer (fails otherwise), makes the given obs columns
categorical, and merges clinical metadata on `projid` (or `individualID`).

### `data diverse-cohorts`
```bash
python -m rosmap_processing data diverse-cohorts <gex_mat_X.h5> --output <X.h5ad> \
    [--download-dir data/v2/download/DiverseCohorts] [--clinical <ROSMAP_clinical.csv>] \
    [--region "dorsolateral prefrontal cortex"|all]
```
One DiverseCohorts cell class (the `.h5` that `rdata_to_h5.R` writes from a
`gex_mat_*.rda`) -> a count file. The libraries pool three donors, so cells are
kept only if the 2024 cell table has them (and they are from `--region`); the
table gives donor, specimen, region, pool (`batch`), cell classes and Rush's
cogdx/Braak/CERAD. Adds the individual metadata, projid (Rush mapping, else the
clinical CSV) and the ROSMAP clinical columns joined on individualID (ROS and
MAP donors only). Run for every class by `slurm/run_all_diverse_cohorts.sh`.

### `data psychad`
```bash
python -m rosmap_processing data psychad <merged_final_clean.h5ad> --output <combined.h5ad> \
    [--download-dir data/v2/download/PsychAD] [--clinical <ROSMAP_clinical.csv>]
```
PsychAD's counts (`layers['counts']`; X there is normalized) -> X as int32,
streamed in chunks of cells (6.3 million cells, 24.6 billion values) and checked
to be counts. obs gains `individualID`, the NPS-AD individual metadata, the
tiered diagnoses, `hashing_set` and `channel` (hashing set + pool replicate from
the cell name), and the ROSMAP clinical columns for Rush donors. Genes keep
their Ensembl IDs. Run by `slurm/run_all_psychad.sh`.

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
| DiverseCohorts, 27,984 genes (name matching) | `cellranger-GRCh38-2020-A` | GENCODE v32 |
| PsychAD, 34,890 genes named by Ensembl ID (kept; symbols added) | `cellranger-GRCh38-2020-A` | GENCODE v32 |

### `genes shared`
```bash
python -m rosmap_processing genes shared <map_a.tsv> <map_b.tsv> [<map_c.tsv> ...] --output <shared_ids.txt>
```
Ensembl IDs mapped in every given map, in the first map's order.

### `genes on-chromosomes`
```bash
python -m rosmap_processing genes on-chromosomes --gtf <a.gtf.gz> [<b.gtf.gz> ...] \
    --chromosomes X Y --output <ids.txt>
```
Ensembl IDs (without versions; a chrY PAR copy counts as its X gene) of every
gene on the given chromosomes in any of the GTFs (`X` and `chrX` alike), e.g.
to leave the sex-chromosome genes out of an analysis.

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

### `pipeline batch-effects`
```bash
python -m rosmap_processing pipeline batch-effects config/batch_effects/<name>.yaml [--overwrite]
```
How much variation is technical, per gene set, from the count files (so every
gene set is compared on the same cells and normalization). The YAML spec lists
the datasets (count file, layer, gene map, donor column, broad cell-type rules,
technical batch columns, donor covariates) and the gene sets (processed h5ad
files or TXT lists); see `config/batch_effects/` and the module docstring of
`pipelines/batch_effects.py`. Applies the cell and donor filters of `pipeline
scanpy`, samples up to 300 cells per donor, and reads the matrix once in row
chunks (pseudobulk over all kept cells). With one dataset it scores its batch
columns; with several, the dataset label (`match_donors: true` keeps only
donors present in every dataset and adds a paired dataset-vs-donor
decomposition: variance components and the correlation of each donor's
deviation between the datasets). With several datasets the cell-level metrics
use a balanced sample (`balance_datasets`, default on): within each cell type
the same number of cells per dataset and per donor, so the larger dataset does
not dominate the neighbourhoods (`cells_balanced` in `cells.tsv`).

Cell level, overall and within each broad cell type (PCA as in `pipeline
scanpy`): PC regression, kNN neighbourhood mixing relative to random,
kBET-style rejection, silhouettes. Donor level (pseudobulk log-CPM per donor ×
cell type): partial R² per term, cross-validated classifier accuracy, and the
association of batch with each covariate. Batch, dataset and covariate metrics
come with a null in which the label is shuffled between donors (`excess` =
value − null mean), because batches hold whole donors. Writes `cells.tsv`,
`cell_metrics.tsv`, `donor_metrics.tsv`, `confounding.tsv`, `summary.tsv` and
`provenance.json`: aggregate rows only (gene set × cell type × label, never a
donor); cell types with fewer than `min_donors` donors are left out.

With a two-level `disease` label for every dataset (e.g. dementia from
`cogdx`), `disease.tsv` measures the signal a correction must keep, per gene
set and cell type on the pseudobulk: per dataset, the share of genes with a
disease effect at FDR < 0.05 and the mean t² (adjusted for the covariates in
`disease_adjust`) and the cross-validated classifier accuracy; per pair of
datasets, the correlation of the per-gene disease t statistics and the
balanced accuracy of a classifier fitted on one dataset and scored on the
other (`transfer_accuracy`: standardized on the training data;
`transfer_accuracy_rescaled`: each dataset standardized on its own). With the
same donors in both datasets, transfer is scored on donors left out of the
fit. Nulls shuffle the disease label between donors. `disease_covariates.tsv`
gives, per dataset, the association of each covariate with the disease
(Cramér's V or η²) and the accuracy of predicting the disease from covariates
alone (a shortcut a model could take).

`donor_models` (named covariate lists; `disease` is the disease label) adds
`donor_effect.tsv`: how much of the variation between donors each covariate,
the disease and each model explains. With one dataset: at cell level, the
share of the PC variance within a cell type (`pc_share`) and that share
relative to donor's (`pc_share_of_donor`); at donor level, each donor's cells
are split in two random halves, whose shared variance is the reproducible
donor effect (`split_half_reliability`: the share of the between-donor
pseudobulk variance that is not cell sampling noise), and
`reproducible_share` is the part of it the covariates explain. With matched
donors, the same between the two datasets (`cross_dataset_reliability`). Nulls
shuffle the covariates between donors. Numeric covariates take a `recode` map
too (e.g. `Braak IV` → 4).

`exclude_genes` (a TXT list of Ensembl IDs or symbols) leaves those genes out
of every gene set. `extends: <spec>.yaml` (relative to the spec) takes another
spec as the base, so a variant names only what differs: the `*_no_xy.yaml`
specs rerun the others without the X and Y genes (`genes on-chromosomes`, made
by `slurm/gene_maps.sh`) into their own output directories.

Specs for the newer datasets: `diverse_cohorts.yaml` and `psychad.yaml`
(within each), `mit_vs_diverse_cohorts.yaml` and `mit_vs_psychad.yaml` (all
donors, as `mit_vs_seaad.yaml`), and the shared-donor specs
`{mit,green}_vs_{diverse_cohorts,psychad}_shared.yaml` and
`diverse_cohorts_vs_psychad_shared.yaml` (as `mit_vs_green.yaml`, matched on
individualID). PsychAD shares only 6 donors with MIT and 4 with Green; the
Green one stops with "The disease label needs exactly 2 values".

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

### `utils obs-summary`
```bash
python -m rosmap_processing utils obs-summary <a.h5ad> ... --output summary.tsv [--donor-column COL]
```
One row per obs column: dtype, number of distinct values, % missing, and
whether it is constant within donors (`donor`), varies within them (`cell`) or
is `constant`, with the number of distinct values per donor. Also logs the
stored matrices and uns keys. Never prints or writes a value, so it is safe
for restricted data.

### `utils validate`
```bash
python -m rosmap_processing utils validate <file.h5ad>
```
Prints structure, obs value counts and a 20×20 preview of X. **The preview shows
data values: do not use it on restricted data where logs are shared.**
