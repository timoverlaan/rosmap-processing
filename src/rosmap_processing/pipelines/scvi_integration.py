"""scVI and scANVI on the cells of the pooled pseudobulk datasets, scored
like the pseudobulk corrections (pipelines/pooled_pseudobulk.py).

Spec (config/scvi_integration.yaml): pooled_spec (the pooled pseudobulk spec,
for the datasets, cell filters, cell types, donor table, batch_split, the
sample and gene filters, protect and the disease models), output_dir, and the
settings in DEFAULTS.

Stages:
  cells  (CPU) per dataset, the cell filters and broad cell types of the
         pooled pseudobulk, then up to max_cells_per_sample nuclei per
         dataset x donor x cell type, at random; their counts for the genes
         every dataset has -> cells.h5ad (restricted: donor IDs), written
         one row chunk at a time (all genes of all cells do not fit in memory)
  train  (GPU) seurat_v3 HVGs over the batches, as scanpy's batch_key
         version but reading one batch at a time; scVI with the batch (the
         dataset; PsychAD x brain bank) as batch_key; scANVI from it with the
         broad cell types as labels. Per model, each cell's decoded
         expression averaged over all batches, and its latent position,
         averaged per sample (dataset x donor x cell type) -> samples.h5ad
         (restricted). Then the scoring below.
  score  (CPU or GPU) per cell type, on the samples with at least min_cells
         sampled nuclei and the genes with a mean CPM of at least min_cpm in
         every dataset: the same nuclei's raw pseudobulk as measured (none)
         and with the regress correction, against scVI and scANVI. Mixing,
         replicates, reliability between datasets and the disease models as
         in the pooled pseudobulk (metrics.tsv, disease.tsv), and UMAPs of the
         samples (umap/overview.png, umap/<cell type>.png).

The raw pseudobulk here sums at most max_cells_per_sample nuclei, so its
numbers are not those of pipeline pooled-pseudobulk, which uses all nuclei.

Run: python -m rosmap_processing pipeline scvi-integration config/scvi_integration.yaml --stage cells
"""

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from ..core import batch_metrics as bm
from ..utils.counts import count_non_integer_values, total_counts
from ..utils.logging import get_logger
from ..utils.provenance import provenance_record
from . import pooled_pseudobulk as pp
from .batch_effects import _matrix_node, prepare_dataset

try:
    from anndata.io import sparse_dataset, write_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import sparse_dataset, write_elem

logger = get_logger(__name__)

DEFAULTS: Dict[str, Any] = {
    "seed": 0,
    "max_cells_per_sample": 100,
    "n_hvg": 3000,
    "n_latent": 30,
    "n_layers": 2,
    "n_hidden": 128,
    "max_epochs": 30,
    "scanvi_epochs": 10,
    "batch_size": 1024,
    "decode_chunk": 20000,
    "read_chunk": 200000,
    "umap_methods": ["none", "regress", "scvi", "scanvi"],
}
REQUIRED = ("pooled_spec", "output_dir")
MODELS = ("scvi", "scanvi")
STAGES = ("cells", "train", "score")


def load_spec(path: Union[str, Path]) -> Dict[str, Any]:
    with open(path) as f:
        raw = yaml.safe_load(f) or {}
    unknown = set(raw) - set(DEFAULTS) - set(REQUIRED) - {"name"}
    if unknown:
        raise ValueError(f"Unknown key(s) in {path}: {sorted(unknown)}")
    missing = [k for k in REQUIRED if not raw.get(k)]
    if missing:
        raise ValueError(f"{path} needs {missing}")
    spec = {**DEFAULTS, "name": Path(path).stem, **raw}
    spec["pooled"] = pp.load_spec(spec["pooled_spec"])
    bad = [m for m in spec["umap_methods"] if m not in ("none", "regress") + MODELS]
    if bad:
        raise ValueError(f"umap_methods: {bad} are not methods")
    return spec


# ---------------------------------------------------------------- cells


def sample_cells(cells: pd.DataFrame, max_per_sample: int, rng: np.random.Generator) -> pd.DataFrame:
    """Up to max_per_sample random cells per donor x cell type, in file order."""
    shuffled = cells.iloc[rng.permutation(len(cells))]
    rank = shuffled.groupby(["donor", "celltype"], sort=False).cumcount().to_numpy()
    return shuffled[rank < max_per_sample].sort_values("row").reset_index(drop=True)


class CSRWriter:
    """X of a new h5ad file, appended in row blocks (float32 data, int32
    indices, int64 indptr), so the whole matrix is never in memory."""

    def __init__(self, f: h5py.File, n_vars: int, chunk: int = 1 << 20):
        self.group = f.create_group("X")
        self.group.attrs["encoding-type"] = "csr_matrix"
        self.group.attrs["encoding-version"] = "0.1.0"
        opts = dict(chunks=(chunk,), maxshape=(None,), compression="gzip", compression_opts=1)
        self.data = self.group.create_dataset("data", shape=(0,), dtype=np.float32, **opts)
        self.indices = self.group.create_dataset("indices", shape=(0,), dtype=np.int32, **opts)
        self.indptr = self.group.create_dataset("indptr", data=np.zeros(1, dtype=np.int64), **opts)
        self.n_obs, self.n_vars = 0, n_vars
        self.group.attrs["shape"] = (0, n_vars)

    def append(self, block) -> None:
        block = sp.csr_matrix(block, dtype=np.float32)
        if block.shape[1] != self.n_vars:
            raise ValueError(f"block has {block.shape[1]} columns, not {self.n_vars}")
        block.sort_indices()
        nnz, rows = self.data.shape[0], self.indptr.shape[0]
        self.data.resize((nnz + block.nnz,))
        self.data[nnz:] = block.data
        self.indices.resize((nnz + block.nnz,))
        self.indices[nnz:] = block.indices.astype(np.int32)
        self.indptr.resize((rows + block.shape[0],))
        self.indptr[rows:] = block.indptr[1:].astype(np.int64) + nnz
        self.n_obs += block.shape[0]
        self.group.attrs["shape"] = (self.n_obs, self.n_vars)


def write_h5ad_rest(f: h5py.File, obs: pd.DataFrame, var: pd.DataFrame, uns: Dict[str, Any]) -> None:
    """Everything but X, to a file whose X a CSRWriter wrote."""
    if len(obs) != f["X"].attrs["shape"][0] or len(var) != f["X"].attrs["shape"][1]:
        raise ValueError("obs and var do not match X")
    f.attrs["encoding-type"] = "anndata"
    f.attrs["encoding-version"] = "0.1.0"
    write_elem(f, "obs", obs)
    write_elem(f, "var", var)
    write_elem(f, "uns", uns)
    for key in ("layers", "obsm", "obsp", "varm", "varp"):
        write_elem(f, key, {})


def stream_sample(ds, rows: np.ndarray, columns: np.ndarray, chunk_size: int, writer: CSRWriter) -> np.ndarray:
    """The counts of the given file rows (sorted) for the given file columns
    (in their order), appended to writer one row chunk at a time; returns the
    rows' totals over all genes."""
    totals = np.zeros(len(rows))
    start_time = time.time()
    with h5py.File(ds.path, "r") as f:
        matrix = sparse_dataset(_matrix_node(f, ds.cfg.get("layer")))
        n_obs = matrix.shape[0]
        for start in range(0, n_obs, chunk_size):
            stop = min(start + chunk_size, n_obs)
            lo, hi = np.searchsorted(rows, [start, stop])
            if lo == hi:
                continue
            block = sp.csr_matrix(matrix[start:stop])[rows[lo:hi] - start]
            bad = count_non_integer_values(block)
            if bad:
                raise ValueError(f"[{ds.name}] rows {start}-{stop}: {bad} values are not counts")
            totals[lo:hi] = total_counts(block)
            writer.append(block[:, columns])
            if (start // chunk_size) % 20 == 0:
                logger.info(f"[{ds.name}] read {stop}/{n_obs} rows ({time.time() - start_time:.0f} s)")
    return totals


def build_cells(spec: Dict[str, Any], overwrite: bool = False) -> Path:
    pooled = spec["pooled"]
    out = Path(spec["output_dir"])
    path = out / "cells.h5ad"
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass --overwrite to replace it")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(spec["seed"])
    names = list(pooled["datasets"])
    table = pp.read_donor_table(pooled["donor_table"], names)
    settings = {"min_genes": pooled["min_genes"], "min_cells_per_donor": 1}
    datasets = [prepare_dataset(k, {**cfg, "batch": [], "covariates": {}}, settings)
                for k, cfg in pooled["datasets"].items()]
    gene_ids, columns = pp.shared_genes(datasets)
    var = pd.concat([pd.DataFrame({"gene_id": gene_ids}), pp._symbols(datasets, columns)], axis=1)
    var.index = var["gene_id"].to_numpy()
    partial = path.with_suffix(".h5ad.partial")
    parts = []
    with h5py.File(partial, "w") as f:
        writer = CSRWriter(f, len(gene_ids))
        for ds in datasets:
            sample = sample_cells(ds.cells, spec["max_cells_per_sample"], rng)
            nuclei = ds.cells.groupby(["donor", "celltype"], sort=False, observed=True).size()
            part = pd.DataFrame({"dataset": ds.name, "donor": sample["donor"].astype(str).to_numpy(),
                                 "celltype": sample["celltype"].to_numpy()})
            part["sample_nuclei"] = nuclei.reindex(
                pd.MultiIndex.from_arrays([sample["donor"], sample["celltype"]])).to_numpy()
            part = pp.attach_persons(part, table)   # cells of donors in the donor table; index = row
            sample = sample.iloc[part.index.to_numpy()]
            logger.info(f"[{ds.name}] sampled {len(sample)} of {len(ds.cells)} cells (up to "
                        f"{spec['max_cells_per_sample']} per donor x cell type, donors in the donor table)")
            part = part.reset_index(drop=True)
            part["library_size"] = stream_sample(ds, sample["row"].to_numpy(), columns[ds.name],
                                                 pooled["chunk_size"], writer)
            parts.append(part)
            ds.cells = None
        obs = pd.concat(parts, ignore_index=True)
        obs["batch"] = pp.assign_batches(obs, pooled["batch_split"])
        obs.index = obs.index.astype(str)
        for c in obs.columns:
            if obs[c].dtype == object:
                obs[c] = obs[c].astype("category")
        uns = {"provenance_scvi_cells": provenance_record(
            "scvi_integration cells", {"spec": spec["name"], "max_cells_per_sample": spec["max_cells_per_sample"]})}
        write_h5ad_rest(f, obs, var, uns)
        nnz = writer.data.shape[0]
    partial.replace(path)
    logger.info(f"Wrote {path}: {len(obs)} cells x {len(var)} genes ({nnz} stored counts), "
                f"{obs.groupby(['dataset', 'donor', 'celltype'], observed=True).ngroups} samples")
    return path


# ---------------------------------------------------------------- train


def _read_rows(X, rows: np.ndarray, chunk: int):
    """Rows (sorted) of a backed CSR matrix, read in contiguous chunks."""
    parts = []
    for start in range(int(rows[0]), int(rows[-1]) + 1, chunk):
        lo, hi = np.searchsorted(rows, [start, start + chunk])
        if hi > lo:
            first, stop = int(rows[lo]), int(rows[hi - 1]) + 1
            parts.append(sp.csr_matrix(X[first:stop])[rows[lo:hi] - first])
    return sp.vstack(parts, format="csr")


def combine_hvg_ranks(norm_gene_vars: np.ndarray, index: pd.Index, n_top: int) -> pd.DataFrame:
    """seurat_v3 HVGs from each batch's normalized variances (batches x genes),
    as scanpy's highly_variable_genes(flavor="seurat_v3", batch_key=...)."""
    ranked = np.argsort(np.argsort(-norm_gene_vars, axis=1), axis=1).astype(np.float32)
    nbatches = np.sum((ranked < n_top).astype(int), axis=0)
    ranked[ranked >= n_top] = np.nan
    median = np.ma.median(np.ma.masked_invalid(ranked), axis=0).filled(np.nan)
    df = pd.DataFrame({"highly_variable_nbatches": nbatches, "highly_variable_rank": median,
                       "variances_norm": np.mean(norm_gene_vars, axis=0)}, index=index)
    order = df.sort_values(["highly_variable_rank", "highly_variable_nbatches"],
                           ascending=[True, False], na_position="last").index
    df["highly_variable"] = False
    df.loc[order[:int(n_top)], "highly_variable"] = True
    return df


def batch_hvgs(X, var: pd.DataFrame, batch: np.ndarray, n_top: int, chunk: int) -> pd.DataFrame:
    """seurat_v3 HVGs over the batches of a backed count matrix, reading one
    batch at a time (scanpy's batch_key version needs the whole matrix)."""
    import anndata as ad
    import scanpy as sc

    norm = []
    for b in np.unique(batch):   # scanpy's order
        counts = _read_rows(X, np.flatnonzero(batch == b), chunk)
        stats = sc.pp.highly_variable_genes(ad.AnnData(counts, var=pd.DataFrame(index=var.index)),
                                            flavor="seurat_v3", n_top_genes=n_top, inplace=False,
                                            check_values=False)
        norm.append(stats["variances_norm"].to_numpy())
        logger.info(f"HVG statistics for batch {b}: {counts.shape[0]} cells")
        del counts
    return combine_hvg_ranks(np.vstack(norm), var.index, n_top)


def read_hvg_cells(path: Path, n_top: int, chunk: int):
    """cells.h5ad restricted to its seurat_v3 HVGs over the batches, without
    loading all of its genes at once."""
    import anndata as ad

    cells = ad.read_h5ad(path, backed="r")
    try:
        logger.info(f"{cells.n_obs} cells x {cells.n_vars} genes")
        batch = cells.obs["batch"].astype(str).to_numpy()
        hvg = batch_hvgs(cells.X, cells.var, batch, n_top, chunk)
        keep = np.flatnonzero(hvg["highly_variable"].to_numpy())
        parts = [sp.csr_matrix(cells.X[start:min(start + chunk, cells.n_obs)])[:, keep]
                 for start in range(0, cells.n_obs, chunk)]
        var = cells.var.iloc[keep].join(hvg.iloc[keep].drop(columns="highly_variable"))
        adata = ad.AnnData(sp.vstack(parts, format="csr"), obs=cells.obs.copy(), var=var)
    finally:
        cells.file.close()
    logger.info(f"{adata.n_vars} HVGs (seurat_v3 over {len(set(batch))} batches)")
    return adata


def _sample_codes(obs: pd.DataFrame):
    keys = ["dataset", "donor", "celltype"]
    codes, groups = pd.MultiIndex.from_frame(obs[keys].astype(str)).factorize()
    return codes, groups


def decode_samples(model, adata, codes: np.ndarray, n_groups: int, chunk: int):
    """Per sample, the mean over its cells of the decoded expression (averaged
    over all batches; proportions over the model's genes) and of the latent
    mean."""
    batches = list(adata.obs["batch"].cat.categories)
    expr = np.zeros((n_groups, adata.n_vars))
    latent = None
    for start in range(0, adata.n_obs, chunk):
        idx = np.arange(start, min(start + chunk, adata.n_obs))
        indicator = bm._indicator(codes[idx], n_groups)
        values = model.get_normalized_expression(adata, indices=idx, transform_batch=batches,
                                                 library_size=1.0, return_numpy=True)
        expr += indicator @ np.asarray(values, dtype=np.float64)
        z = model.get_latent_representation(adata, indices=idx)
        latent = np.zeros((n_groups, z.shape[1])) if latent is None else latent
        latent += indicator @ z
        if (start // chunk) % 10 == 0:
            logger.info(f"decoded {idx[-1] + 1}/{adata.n_obs} cells")
    size = np.bincount(codes, minlength=n_groups)[:, None]
    return expr / size, latent / size


def train_models(spec: Dict[str, Any], overwrite: bool = False) -> Path:
    import anndata as ad
    import scvi
    import torch

    out = Path(spec["output_dir"])
    path = out / "samples.h5ad"
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass --overwrite to replace it")
    scvi.settings.seed = spec["seed"]
    torch.set_float32_matmul_precision("high")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    logger.info(f"scvi-tools {scvi.__version__}, torch {torch.__version__}, accelerator {accelerator}"
                + (f" ({torch.cuda.get_device_name(0)})" if accelerator == "gpu" else ""))

    adata = read_hvg_cells(out / "cells.h5ad", spec["n_hvg"], spec["read_chunk"])
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["celltype"] = adata.obs["celltype"].astype("category")

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    models, histories = {}, []
    model = scvi.model.SCVI(adata, n_latent=spec["n_latent"], n_layers=spec["n_layers"],
                            n_hidden=spec["n_hidden"], gene_likelihood="nb")
    train_kwargs = dict(batch_size=spec["batch_size"], accelerator=accelerator, devices=1, early_stopping=True)
    model.train(max_epochs=spec["max_epochs"], **train_kwargs)
    model.save(str(out / "scvi_model"), overwrite=True)
    models["scvi"] = model
    scanvi = scvi.model.SCANVI.from_scvi_model(model, unlabeled_category="Unknown", labels_key="celltype")
    scanvi.train(max_epochs=spec["scanvi_epochs"], n_samples_per_label=100, **train_kwargs)
    scanvi.save(str(out / "scanvi_model"), overwrite=True)
    models["scanvi"] = scanvi
    for name, m in models.items():
        h = pd.concat([v.rename(columns=lambda c: k) for k, v in m.history.items()
                       if k.startswith(("elbo", "reconstruction_loss", "kl_local"))], axis=1)
        histories.append(h.assign(model=name).reset_index(names="epoch"))
        logger.info(f"{name}: {len(h)} epochs")
    pd.concat(histories, ignore_index=True).to_csv(out / "history.tsv", sep="\t", index=False)

    codes, groups = _sample_codes(adata.obs)
    n_groups = len(groups)
    first = pd.Series(np.arange(adata.n_obs)).groupby(codes).first().to_numpy()
    obs = adata.obs.iloc[first].reset_index(drop=True)
    obs["n_cells"] = np.bincount(codes, minlength=n_groups)
    obs["library_size"] = np.bincount(codes, weights=adata.obs["library_size"].to_numpy(np.float64),
                                      minlength=n_groups)
    counts = bm._indicator(codes, n_groups) @ adata.X
    samples = ad.AnnData(sp.csr_matrix(np.rint(counts.toarray() if sp.issparse(counts) else counts)
                                       .astype(np.int32)), obs=obs, var=adata.var.copy())
    for name, m in models.items():
        expr, latent = decode_samples(m, adata, codes, n_groups, spec["decode_chunk"])
        samples.layers[name] = expr.astype(np.float32)
        samples.obsm[f"latent_{name}"] = latent.astype(np.float32)
    samples.obs.index = (obs["dataset"].astype(str) + ":" + obs["donor"].astype(str) + ":"
                         + obs["celltype"].astype(str)).to_numpy()
    samples.uns["provenance_scvi_train"] = provenance_record(
        "scvi_integration train", {"spec": spec["name"], "scvi_tools": scvi.__version__})
    samples.write_h5ad(path, compression="gzip")
    logger.info(f"Wrote {path}: {samples.n_obs} samples x {samples.n_vars} genes")
    return path


# ---------------------------------------------------------------- score


def score_samples(spec: Dict[str, Any]) -> Path:
    """metrics.tsv, disease.tsv and the UMAPs from samples.h5ad."""
    import anndata as ad

    pooled = spec["pooled"]
    out = Path(spec["output_dir"])
    samples = ad.read_h5ad(out / "samples.h5ad")
    obs = samples.obs.copy()
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(object).where(obs[c].notna(), np.nan)
    obs = obs.reset_index(drop=True)
    counts = samples.X.toarray() if sp.issparse(samples.X) else np.asarray(samples.X)
    genes = samples.var.reset_index(drop=True)
    models = [m for m in MODELS if m in samples.layers]
    numeric = [c for c in pooled["protect"] if c in pp.NUMERIC]
    rng = np.random.default_rng(spec["seed"])
    umap_spec = {**pooled, "umap_methods": spec["umap_methods"]}

    metrics, disease, embeddings = [], [], {}
    for celltype in sorted(obs["celltype"].unique()):
        sel = ((obs["celltype"] == celltype) & (obs["n_cells"] >= pooled["min_cells"])).to_numpy()
        sub = obs[sel].reset_index(drop=True)
        ds = sub["dataset"].to_numpy().astype(str)
        if len(set(ds)) < 2:
            continue
        cpm = counts[sel] / sub["library_size"].to_numpy()[:, None] * 1e6
        keep = np.all([cpm[ds == d].mean(axis=0) >= pooled["min_cpm"] for d in dict.fromkeys(ds)], axis=0)
        Y = np.log2(cpm[:, keep] + 1)
        batch = sub["batch"].to_numpy().astype(str)
        corrected = {"none": Y, "regress": pp.correct_regress(Y, batch, sub[list(pooled["protect"])], numeric)[0]}
        for m in models:
            corrected[m] = np.log2(np.asarray(samples.layers[m])[sel][:, keep] * 1e6 + 1)
        logger.info(f"{celltype}: {len(sub)} samples, {int(keep.sum())} of {len(keep)} HVGs "
                    f"(mean CPM >= {pooled['min_cpm']} in every dataset)")
        persons = sub["person_id"].to_numpy().astype(str)
        rows, pcs = pp.score_methods(celltype, sub, corrected, pooled, rng)
        metrics.append(dict(celltype=celltype, method="any", metric="n_genes", label="", value=int(keep.sum()),
                            chance=np.nan, n=len(sub)))
        metrics += rows
        for method, M in corrected.items():
            metrics += pp.reliability_rows(celltype, method, M, ds, persons, pooled["min_pairs"])
        d, _ = pp.disease_analysis(celltype, corrected, sub, genes[keep].reset_index(drop=True), pooled)
        disease += d
        coords = {m: pp._umap(pcs[m], umap_spec) for m in spec["umap_methods"] if m in pcs}
        embeddings[celltype] = (sub, coords)

    metrics = pd.DataFrame(metrics, columns=pp.METRIC_COLUMNS)
    metrics.to_csv(out / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(disease).to_csv(out / "disease.tsv", sep="\t", index=False)

    figures = out / "umap"
    figures.mkdir(exist_ok=True)
    palette = pp.batch_palette(obs, list(pooled["datasets"]))
    titles = {**pp.method_titles(pooled["batch_split"]),
              "none": f"as measured ({spec['max_cells_per_sample']} nuclei max)",
              "scvi": "scVI (batch key: batch)", "scanvi": "scANVI (+ broad cell types)"}
    label = "batch" if (obs["batch"] != obs["dataset"]).any() else "dataset"
    share = metrics[(metrics["metric"] == "pc_r2") & (metrics["label"] == label)]
    mixing = {(r.celltype, r.method): float(r.value) for r in share.itertuples()}
    for celltype, (sub, coords) in embeddings.items():
        pp.plot_celltype_umaps(celltype, sub, coords, palette, titles,
                               {m: v for (t, m), v in mixing.items() if t == celltype},
                               figures / f"{celltype}.png", spec["seed"])
    if embeddings:
        pp.plot_overview(embeddings, palette, titles, mixing, figures / "overview.png", spec["seed"])
    with open(out / "provenance.json", "w") as f:
        json.dump(provenance_record("scvi_integration score", {"spec": spec["name"]})
                  | {"spec": {k: v for k, v in spec.items() if k != "pooled"}}, f, indent=2, default=str)
    logger.info(f"Wrote metrics.tsv, disease.tsv and umap/ to {out}")
    return out


def run_scvi_integration(spec_path: Union[str, Path], stage: str, overwrite: bool = False) -> Path:
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES}")
    spec = load_spec(spec_path)
    if stage == "cells":
        return build_cells(spec, overwrite)
    if stage == "train":
        train_models(spec, overwrite)
    return score_samples(spec)


def main(argv: Optional[List[str]] = None) -> None:
    import argparse

    from ..utils.logging import setup_logging

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("spec")
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(level="INFO")
    run_scvi_integration(args.spec, args.stage, args.overwrite)


if __name__ == "__main__":
    main()
