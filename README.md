# ROSMAP Processing Pipeline

Processing of single-nucleus RNA-seq data from ROSMAP (two atlases) and SeaAD
into raw-count h5ad files, per-donor kNN graph datasets, and per-donor
aggregates for the ad-prs project.

> **Restricted data.** ROSMAP data are access-controlled. Do not copy clinical
> metadata, obs tables or expression values out of the cluster, into logs that
> are shared, or into issues and chats. `utils check-counts` reports aggregate
> statistics only; `utils validate` prints data values.

## Data sources

| Name in this repo | Dataset | Cells / donors | Genes (10x reference) | Source |
|---|---|---|---|---|
| `ROSMAP_MIT` | Mathys et al., *Cell* 2023 (MIT) | 2,342,020 / 427 | 33,538 (Cell Ranger GRCh38-3.0.0, Ensembl 93) | [Synapse syn52293417](https://www.synapse.org/#!Synapse:syn52293417) |
| `ROSMAP` | Green et al., *Nature* 2024 (De Jager / Habib labs) | 1,649,672 / 437 | 18,745–20,326 per cell class (GRCh38-2020-A, filtered per class) | Synapse syn53366818 (see `SYNAPSE_IDS_ROSMAP`) |
| `SeaAD` | SEA-AD, prefrontal cortex (A9) | 1,395,601 / 83 | 36,601 (GRCh38-2020-A, GENCODE v32) | [AWS open data](https://registry.opendata.aws/allen-sea-ad-atlas/) |
| `DiverseCohorts` | AMP-AD DiverseCohorts (Rush), 10x multiome gene expression, DLPFC only | 1,812,921 / 412 | 27,984 (Ensembl 91 names; all mapped) | Synapse, study AMP-AD_DiverseCohorts (see `SYNAPSE_IDS_DIVERSE_COHORTS`) |
| `PsychAD` | PsychAD (Lee et al. 2024, Mount Sinai), hashed donor pools, DLPFC | 6,320,459 / 1,494 | 34,890 (named by Ensembl ID) | Synapse, study MSSM_PsychAD (see `SYNAPSE_IDS_PSYCHAD`) |

In the Green 2024 objects each cell class keeps only genes detected in at least
15 of its nuclei ([`subset.cell.type.R`](https://github.com/naomihabiblab/BEYOND_DLPFC/blob/487325e9624e53d34dca516da58767c386677d60/2.%20Cell-type%20analysis/utils/subset.cell.type.R#L57-L69)),
so combining the classes fills those genes with zeros, which is a close
approximation.

`data prepare` reads the RNA `counts` slot and `meta.data` of each `.h5Seurat`
file in Python (`data/h5seurat.py`). SeuratDisk is not used for this: with the
container's SeuratObject its `LoadH5Seurat` fails (`GetAssayData(slot = )` is
defunct), and its h5ad conversion puts Seurat's `data` slot in X and loses
factor labels. The ROSMAP-MIT `.h5Seurat` files were made from the RDS files
with SeuratDisk in April 2026 (`rds_to_h5Seurat.R`, which only runs for missing
files).

Every Synapse file is pinned to a version and MD5 in `SYNAPSE_FILES`
(`utils/constants.py`); `data download` and `data verify` check them. All of
these entities got a new version in July 2026 (2026-07-08/09/14), but each new
version has exactly the previous version's MD5 and size: no file changed.

## Layout on the cluster (v2)

```
data/raw/ROSMAP_MIT/          Mathys 2023 RDS + metadata (MD5-verified in place) and the .h5Seurat made from them
data/raw/ROSMAP/              old Green copies, NOT used (rewritten in 2025, see "Data history")
data/v2/download/ROSMAP/      Green 2024 objects + ROSMAP metadata, fresh pinned downloads
data/v2/download/PsychAD/     PsychAD h5ad (after the authors' QC) + donor metadata
data/v2/download/DiverseCohorts/  DiverseCohorts multiome expression matrices (.rda), cell tables, donor metadata
data/v2/intermediate/DiverseCohorts/  the .rda matrices as HDF5 (src/convert_R/rdata_to_h5.R)
data/reference/               Ensembl 93 and GENCODE v32 GTFs
data/v2/counts/<DATASET>/     raw counts in X (int32) + clinical metadata, per class and combined.h5ad
                              (RNA counts read straight from the .h5Seurat files, see below)
data/v2/counts/SeaAD/         SeaAD with ROSMAP-style columns (counts in layer 'UMIs')
data/v2/reference/            gene maps (symbol -> Ensembl ID), shared_gene_ids.txt (MIT+SeaAD), shared_gene_ids_3way.txt (+Green)
data/v2/processed/            processed datasets (below), mathys2023_for_adprs/
data/v2/hvg_compare/          HVG variant lists and overlaps, per dataset
data/v2/batch_effects/        obs column summaries and batch-effect tables (aggregates only)
```

Processed datasets:

| Files | Genes |
|---|---|
| `rosmap_mit_<name>.h5ad` + `seaad_<name>.h5ad` (pairs) | chosen on ROSMAP-MIT among genes SeaAD also has; SeaAD gets exactly the same genes in the same order |
| `rosmap_mit_hvg2000_3way_k30.h5ad` + `seaad_...` + `rosmap_green_...` (triplet) | chosen on ROSMAP-MIT among genes all three have; all three identical |
| `rosmap_green_<name>.h5ad` (standalone) | Green 2024's own HVGs |

Genes are matched by Ensembl ID (`var['gene_id']`); `var_names` are the
ROSMAP-MIT symbols in pairs and triplets. Each processed file has
`layers['counts']`, X = log1p(counts / total × 1e4) with totals over all genes,
`obs['total_counts']`, `obs['Donor ID']` (projid; for Green 2024 its
individualID), a per-donor kNN graph, and `uns['provenance_*']` with parameters
and the git commit. Cells without a donor ID (in Green 2024, the text `NA`) and
donors with fewer than 300 cells are dropped (counted in the log): the per-donor
seurat_v3 HVG fit fails on the smallest ROSMAP-MIT donors.

Green 2024 and Mathys 2023 are both ROSMAP DLPFC atlases and very likely share
many participants: do not use one to validate a model trained on the other.
SeaAD is an independent cohort. DiverseCohorts is all Rush (ROS, MAP, MARS,
CLINCORE) and PsychAD has 152 Rush donors; their Rush donors carry the ROSMAP
individualID. DiverseCohorts shares 100 donors with MIT and 109 with Green;
PsychAD only 6 and 4, and the two share 87 with each other.

## Running on DAIC

Submit from the repo root, with a valid Kerberos ticket (`klist`; `kinit -R`).
The scripts use account `daic-eemcs-insy-prb`, partition `all`, and the
container `container_pixi_0-2-1.sif` (code is mounted from `src/`, so code
changes need no rebuild). Each step skips outputs that already exist, so a
failed job can be resubmitted.

```bash
bash slurm/download_references.sh                 # login node: the two GTFs (~85 MB)
sbatch slurm/smoke_test.sh                        # unit tests + R check in the container (minutes)
sbatch slurm/run_all.sh                           # downloads Green (~112 GB), verifies MIT, -> data/v2/counts (long)
sbatch slurm/run_all_seaad.sh                     # SeaAD -> data/v2/counts/SeaAD (optional, see below)
sbatch slurm/inspect_counts.sh                    # check the count files (aggregate output only)
sbatch slurm/gene_maps.sh                         # gene maps + shared Ensembl IDs
sbatch slurm/hvg_compare.sh                       # HVG variants (per donor vs pooled, ...) and their overlaps
sbatch slurm/regular_hvg.sh                       # hvg{1000,2000,3000}_k30 pairs
sbatch slurm/regular_hvg_postnorm.sh              # hvg{...}_postnorm_k30 pairs
sbatch slurm/lieke_hvg.sh                         # lieke_hvg_k30 pair
sbatch slurm/lieke_allgenes_hvg.sh                # lieke_allgenes_hvg2k_k30 pair
sbatch slurm/lieke_allgenes_hvg_postnorm.sh       # lieke_allgenes_hvg2k_postnorm_k30 pair
sbatch slurm/rosmap_mathys_genes.sh               # mathys2019degs_k30 pair
sbatch slurm/aligned_hvg.sh                       # hvg2000_3way_k30 triplet (MIT, SeaAD, Green)
sbatch slurm/green_hvg.sh                         # Green 2024 standalone: hvg{1000,2000,3000}[_postnorm]_k30
sbatch slurm/mathys2023_for_adprs.sh              # ad-prs aggregates
sbatch slurm/obs_summary.sh                       # obs column names/types of the count files (no values)
sbatch slurm/batch_effects.sh config/batch_effects/rosmap_mit.yaml    # batch effects within MIT
sbatch slurm/batch_effects.sh config/batch_effects/rosmap_green.yaml  # within Green 2024
sbatch slurm/batch_effects.sh config/batch_effects/seaad.yaml         # within SeaAD (library, chemistry)
sbatch slurm/batch_effects.sh config/batch_effects/mit_vs_seaad.yaml  # dataset effect MIT vs SeaAD
sbatch slurm/batch_effects.sh config/batch_effects/mit_vs_green.yaml  # MIT vs Green, shared donors
sbatch slurm/download.sh diverse-cohorts psychad  # pinned DiverseCohorts and PsychAD files (~175 GB)
sbatch slurm/run_all_diverse_cohorts.sh           # -> data/v2/counts/DiverseCohorts (DLPFC; about 1.5 h)
sbatch slurm/run_all_psychad.sh                   # -> data/v2/counts/PsychAD (about 1 h)
sbatch slurm/gene_maps.sh                         # again, for their gene maps
sbatch slurm/batch_effects.sh config/batch_effects/<spec>.yaml  # diverse_cohorts*, psychad*, mit_vs_*, green_vs_*
```

The batch-effect specs name obs columns (batch, cell type, covariates); check
them against the obs summary before submitting. Their outputs are aggregate
tables (see `pipeline batch-effects` in CLI_REFERENCE.md).

HVGs are selected per donor by default (`--hvg-batch-key "Donor ID"`): genes
must vary between cells within donors, so donor-level differences (sex genes,
batch) do not drive the selection, but neither do genuine between-donor
(e.g. disease) differences. For a pooled variant of any pair, add
`--hvg-batch-key none` with a new name, e.g.
`process_pair "hvg2000_pooled_k30" --n-genes 2000 --hvg-batch-key none`.

Chain them with `--dependency=afterok:<jobid>`. Without `run_all_seaad.sh`, the
SeaAD jobs read the original download in `data/seaAD/` (set `SEAAD_H5AD` to
override).

## Data history: the April 2026 outputs

The files in `data/raw/*/`*.h5ad, `data/processed/` and `data/seaAD/.../seaad_*`
from April 2026 should not be used:

- SeuratDisk writes Seurat's `data` slot to X and the counts to `.raw`. The old
  pipeline then deleted `.raw`. For ROSMAP-MIT Immune_cells, Inhibitory_neurons
  and OPCs (504,090 cells), and Green cux2+ and inhibitory, `data` was
  log-normalized, so X held log-normalized values that were normalized again.
  The ad-prs `mean_expr` for In, Mic and Opc is affected the same way.
- Library sizes were computed after HVG subsetting, with a 1e6 target.
- The ROSMAP-MIT and SeaAD outputs had different gene orders, and genes were
  matched by symbol across different 10x references.
- From May to November 2025 the h5Seurat conversion script saved the loaded
  object back over the downloaded file (`SaveH5Seurat(..., overwrite = TRUE)`,
  commits b011398 to 68bd6c0). Six of the seven Green `.h5Seurat` files in
  `data/raw/ROSMAP/` therefore match no Synapse version (only `cux2-`, downloaded
  again afterwards, does). v2 downloads them fresh into `data/v2/download/ROSMAP/`.
- The ad-prs `donor_metadata.tsv` joined clinical data on projid as text, but the
  clinical CSV pads projid to 8 characters, so donors whose projid starts with a
  zero got no clinical fields.

The v2 pipeline fixes all of these; `data prepare`, `pipeline scanpy` and the
ad-prs aggregation now fail on input that is not raw counts, and every input is
checked against its pinned Synapse MD5.

## Development

```bash
pixi install
pixi run pytest -m "not integration"
```

The seurat_v3 test is skipped when scikit-misc's loess does not work (the
Windows build); the container's Linux build is the one that matters.

See [CLI_REFERENCE.md](CLI_REFERENCE.md) for all commands and options.
