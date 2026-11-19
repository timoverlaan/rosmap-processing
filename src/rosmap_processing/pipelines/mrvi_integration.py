"""MrVI (scvi-tools, torch) on the sampled nuclei of the scVI integration
(pipelines/scvi_integration.py), with the persons measured in two or more
datasets as anchors, cross-fitted so that a person's replicate scores never
come from a model that was told which of their samples are the same person.

Spec (config/mrvi_integration.yaml): scvi_spec (the scVI integration spec,
for cells.h5ad, the HVGs and, through its pooled_spec, the scoring settings),
output_dir, and the settings in DEFAULTS.

Design:
  sample_key  the person, so one person's nuclei in two datasets share a
              sample and anchor the datasets to each other; batch_key the
              batch (dataset; PsychAD x brain bank), a nuisance that only the
              decoder sees.
  folds       the bridge persons (nuclei in two or more datasets) are split
              in halves A and B, balanced within each combination of
              datasets. Model A keeps the person as the sample for half A
              (anchors) and gives half B's nuclei one sample per dataset x
              donor (held out); model B the other way round. Each bridge
              person is scored from the model that held them out.

Stages:
  train  (GPU) per fold (--fold A or B; both without): the seurat_v3 HVGs of
         the scVI run (computed once, cached in hvg_<n>.tsv), MrVI, then
         - per sample (dataset x donor x cell type), the mean over its nuclei
           of the decoded expression at library size 1, averaged over all
           batches (as the scVI decoding), and of u and z
           -> samples_<fold>.h5ad (restricted)
         - per cell type, MrVI's local sample-sample distances (l2 between
           the counterfactual z of every sample) averaged over up to
           distance_cells of its nuclei -> distances_<fold>.h5 (restricted)
         - history_<fold>.tsv, model_<fold>/
  score  (CPU) samples.h5ad with layers mrvi_A, mrvi_B and mrvi (a bridge
         person's samples from the model that held them out, everyone else
         the mean of the two models). Per cell type, on the samples and genes
         of the scVI scoring, none and regress against:
           mrvi_A, mrvi_B  one model's decoded expression; its replicate
                           metrics and cross-dataset reliability only over the
                           persons it held out
           mrvi            mixing and disease: the mean of the two models;
                           replicates: both models' held-out persons pooled
                           (every bridge person once, as for none and
                           regress); reliability: the mean of the two
           mrvi_anchors    replicate metrics over the persons a model was
                           given as anchors (circular; only to show by how
                           much)
           mrvi_distance   replicate metrics from MrVI's sample distances:
                           the rank of a held-out person's own sample in the
                           other dataset among all of that dataset's samples
         -> metrics.tsv, disease.tsv, folds.tsv, history.tsv, umap/.
  all    train both folds, then score.

--smoke: smoke_bridge_persons bridge persons and smoke_other_persons others
(all their sampled nuclei), smoke_epochs epochs, in output_dir/smoke.

Run: python -m rosmap_processing pipeline mrvi-integration config/mrvi_integration.yaml --stage train --fold A
"""

import gc
import json
import time
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from ..utils.logging import get_logger
from ..utils.provenance import provenance_record
from . import pooled_pseudobulk as pp
from . import scvi_integration as svi

logger = get_logger(__name__)

DEFAULTS: Dict[str, Any] = {
    "seed": 0,
    "n_latent": 30,
    "n_latent_u": 10,
    "labels_key": None,          # e.g. celltype: one component of the u prior per label
    "max_epochs": 60,
    "early_stopping_patience": 10,
    "batch_size": 1024,
    "lr": 2e-3,
    "n_epochs_kl_warmup": 20,
    "decode_batch": 4096,
    "distance_cells": 2000,      # nuclei per cell type for the sample distances
    "distance_batch": 8,         # nuclei at a time (x samples^2 distances on the GPU)
    "umap_methods": ["none", "regress", "mrvi_A"],
    "smoke_bridge_persons": 30,
    "smoke_other_persons": 20,
    "smoke_epochs": 2,
}
REQUIRED = ("scvi_spec", "output_dir")
FOLDS = ("A", "B")
STAGES = ("train", "score", "all")
SCORED = ("none", "regress", "mrvi_A", "mrvi_B")


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
    spec["scvi"] = svi.load_spec(spec["scvi_spec"])
    spec["pooled"] = spec["scvi"]["pooled"]
    bad = [m for m in spec["umap_methods"] if m not in SCORED]
    if bad:
        raise ValueError(f"umap_methods: {bad} are not among {SCORED}")
    return spec


def run_dir(spec: Dict[str, Any], smoke: bool) -> Path:
    out = Path(spec["output_dir"])
    return out / "smoke" if smoke else out


# ---------------------------------------------------------------- folds


def _dataset_sets(obs: pd.DataFrame) -> pd.Series:
    """Per person, the datasets they have rows in, joined by '+'."""
    pairs = obs[["person_id", "dataset"]].astype(str).drop_duplicates()
    return pairs.groupby("person_id")["dataset"].agg(lambda s: "+".join(sorted(s)))


def assign_folds(obs: pd.DataFrame, seed: int) -> pd.Series:
    """Half A or B for every person with rows in two or more datasets (index
    person_id): at random, alternating within each combination of datasets,
    so each combination is split as evenly as it can be."""
    sets = _dataset_sets(obs)
    bridges = sets[sets.str.contains("+", regex=False)]
    rng = np.random.default_rng(seed)
    fold: Dict[str, str] = {}
    sizes = {"A": 0, "B": 0}
    for combo in sorted(set(bridges)):
        members = sorted(bridges.index[(bridges == combo).to_numpy()])
        first = "A" if sizes["A"] <= sizes["B"] else "B"
        second = "B" if first == "A" else "A"
        for i, k in enumerate(rng.permutation(len(members))):
            f = first if i % 2 == 0 else second
            fold[members[k]] = f
            sizes[f] += 1
    return pd.Series(fold, name="fold", dtype=object).sort_index()


def with_folds(obs: pd.DataFrame, seed: int) -> pd.DataFrame:
    folds = assign_folds(obs, seed)
    return obs.assign(fold=obs["person_id"].astype(str).map(folds).fillna("").astype(str).to_numpy())


def _fold(obs: pd.DataFrame) -> np.ndarray:
    return obs["fold"].astype(object).where(obs["fold"].notna(), "").astype(str).to_numpy()


def mrvi_sample_ids(obs: pd.DataFrame, model: str) -> np.ndarray:
    """MrVI's sample for each row of obs in the given model: the person
    ('p:<person>'), except for the bridge persons of the other half, held
    out, who get one sample per dataset x donor ('d:<dataset>:<donor>')."""
    fold = _fold(obs)
    held = (fold != "") & (fold != model)
    person = ("p:" + obs["person_id"].astype(str)).to_numpy()
    own = ("d:" + obs["dataset"].astype(str) + ":" + obs["donor"].astype(str)).to_numpy()
    return np.where(held, own, person)


def scored_persons(obs: pd.DataFrame, model: str, anchors: bool = False) -> np.ndarray:
    """The person of each row for the replicate metrics of a model: the
    persons it held out keep their person (anchors=True: the persons it
    anchored), the others get a label per dataset, so they form no pairs."""
    fold = _fold(obs)
    other = (fold == model) if not anchors else ((fold != "") & (fold != model))
    person = obs["person_id"].astype(str).to_numpy()
    return np.where(other, ("~" + obs["dataset"].astype(str) + ":" + obs["person_id"].astype(str)).to_numpy(),
                    person)


def fold_table(obs: pd.DataFrame) -> pd.DataFrame:
    """Persons per combination of datasets and half (no IDs)."""
    sets = _dataset_sets(obs)
    first = obs.assign(person_id=obs["person_id"].astype(str)).drop_duplicates("person_id")
    fold = pd.Series(_fold(first), index=first["person_id"].to_numpy())
    frame = pd.DataFrame({"datasets": sets, "fold": fold.reindex(sets.index).fillna("").to_numpy()})
    table = frame.groupby(["datasets", "fold"]).size().unstack(fill_value=0)
    return table.rename(columns={"": "not split"}).reset_index()


def smoke_rows(obs: pd.DataFrame, n_bridge: int, n_other: int, seed: int) -> np.ndarray:
    """Rows (sorted) of n_bridge random bridge persons and n_other others."""
    sets = _dataset_sets(obs)
    bridge = sets.str.contains("+", regex=False)
    rng = np.random.default_rng(seed)
    pick = []
    for pool, n in ((sets.index[bridge.to_numpy()], n_bridge), (sets.index[~bridge.to_numpy()], n_other)):
        pick += list(rng.choice(np.asarray(pool), size=min(n, len(pool)), replace=False))
    return np.flatnonzero(obs["person_id"].astype(str).isin(set(pick)).to_numpy())


# ---------------------------------------------------------------- train


def load_cells(spec: Dict[str, Any], smoke: bool):
    """cells.h5ad on its seurat_v3 HVGs over the batches (computed as in the
    scVI run, cached in output_dir/hvg_<n>.tsv): all nuclei, or with smoke
    those of a subset of persons. obs gets the fold of each bridge person."""
    import anndata as ad

    sv = spec["scvi"]
    path = Path(sv["output_dir"]) / "cells.h5ad"
    cache = Path(spec["output_dir"]) / f"hvg_{sv['n_hvg']}.tsv"
    chunk = sv["read_chunk"]
    cells = ad.read_h5ad(path, backed="r")
    try:
        logger.info(f"{path}: {cells.n_obs} cells x {cells.n_vars} genes")
        obs = cells.obs.copy()
        if cache.exists():
            hvg = pd.read_csv(cache, sep="\t", index_col=0)
            hvg.index = hvg.index.astype(str)
            logger.info(f"{len(hvg)} HVGs from {cache}")
        else:
            stats = svi.batch_hvgs(cells.X, cells.var, obs["batch"].astype(str).to_numpy(), sv["n_hvg"], chunk)
            hvg = stats[stats["highly_variable"].to_numpy()].drop(columns="highly_variable")
            cache.parent.mkdir(parents=True, exist_ok=True)
            partial = cache.with_name(cache.name + ".partial")
            hvg.to_csv(partial, sep="\t")
            partial.replace(cache)
            logger.info(f"{len(hvg)} HVGs (seurat_v3 over {obs['batch'].nunique()} batches) -> {cache}")
        keep = cells.var_names.get_indexer(hvg.index)
        if (keep < 0).any():
            raise ValueError(f"{cache} has {int((keep < 0).sum())} genes that {path} does not; delete it")
        keep = np.sort(keep)
        if smoke:
            rows = smoke_rows(obs, spec["smoke_bridge_persons"], spec["smoke_other_persons"], spec["seed"])
            X = svi._read_rows(cells.X, rows, chunk)[:, keep]
            obs = obs.iloc[rows]
            logger.info(f"smoke: {len(rows)} cells of {obs['person_id'].nunique()} persons")
        else:
            X = sp.vstack([sp.csr_matrix(cells.X[start:min(start + chunk, cells.n_obs)])[:, keep]
                           for start in range(0, cells.n_obs, chunk)], format="csr")
        var = cells.var.iloc[keep].join(hvg)
    finally:
        cells.file.close()
    obs = with_folds(obs, spec["seed"])
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].cat.remove_unused_categories()   # no empty batches to average over
    obs["batch"] = obs["batch"].astype("category")
    obs["celltype"] = obs["celltype"].astype("category")
    folds = obs.drop_duplicates("person_id")["fold"].value_counts()
    logger.info(f"bridge persons: {int(folds.get('A', 0))} in half A, {int(folds.get('B', 0))} in half B; "
                f"{int(folds.get('', 0))} persons in one dataset")
    return ad.AnnData(X, obs=obs, var=var)


def _gpu_peak(torch) -> str:
    if not torch.cuda.is_available():
        return "no GPU"
    return f"peak GPU memory {torch.cuda.max_memory_allocated() / 2**30:.1f} GiB allocated, " \
           f"{torch.cuda.max_memory_reserved() / 2**30:.1f} GiB reserved"


def decode_samples(model, adata, codes: np.ndarray, n_groups: int, batch_size: int):
    """Per sample, the mean over its cells of the decoded expression (the
    cell's own MrVI sample; library size 1, averaged over all batches) and
    of the posterior means of u and z."""
    import torch
    from scvi import REGISTRY_KEYS as K

    module = model.module
    module.eval()
    device = model.device
    n_batch = model.summary_stats.n_batch
    expr = torch.zeros((n_groups, adata.n_vars), dtype=torch.float64, device=device)
    lat_u = lat_z = None
    codes_t = torch.as_tensor(np.asarray(codes), dtype=torch.long, device=device)
    loader = model._make_data_loader(adata=adata, batch_size=batch_size, shuffle=False)
    done = 0
    start_time = time.time()
    with torch.inference_mode():
        for step, tensors in enumerate(loader):
            idx = tensors[K.INDICES_KEY].flatten().long().to(device)
            inf = module.inference(tensors[K.X_KEY], tensors[K.SAMPLE_KEY], use_mean=True)
            z = inf["z"]
            ones = torch.ones((z.shape[0], 1), device=z.device)
            h = torch.zeros((z.shape[0], adata.n_vars), device=z.device)
            for b in range(n_batch):
                h += module.px(z, torch.full((z.shape[0],), b, device=z.device), size_factor=ones).mean
            group = codes_t[idx]
            expr.index_add_(0, group, (h / n_batch).double())
            if lat_u is None:
                lat_u = torch.zeros((n_groups, inf["u"].shape[1]), dtype=torch.float64, device=device)
                lat_z = torch.zeros((n_groups, z.shape[1]), dtype=torch.float64, device=device)
            lat_u.index_add_(0, group, inf["u"].double())
            lat_z.index_add_(0, group, z.double())
            done += z.shape[0]
            if step % 100 == 0:
                logger.info(f"decoded {done}/{adata.n_obs} cells ({time.time() - start_time:.0f} s)")
    size = np.bincount(codes, minlength=n_groups)[:, None].astype(np.float64)
    return expr.cpu().numpy() / size, lat_u.cpu().numpy() / size, lat_z.cpu().numpy() / size


def sample_distances(model, adata, n_cells: int, batch_size: int,
                     rng: np.random.Generator) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Per cell type, MrVI's local sample distances averaged over up to
    n_cells random nuclei of it: for each nucleus, the l2 distances between
    the posterior-mean z it would have in each of the model's samples.
    Returns the sample names (in matrix order) and the matrices."""
    import torch

    from scvi import REGISTRY_KEYS as K

    module = model.module
    module.eval()
    device = model.device
    names = np.asarray(model.sample_order).astype(str)
    n_sample = len(names)
    all_samples = torch.arange(n_sample, device=device)
    out = {}
    celltype = adata.obs["celltype"].astype(str).to_numpy()
    with torch.inference_mode():
        for t in sorted(set(celltype)):
            rows = np.flatnonzero(celltype == t)
            rows = np.sort(rng.choice(rows, size=min(n_cells, len(rows)), replace=False))
            total = torch.zeros((n_sample, n_sample), dtype=torch.float64, device=device)
            start_time = time.time()
            for tensors in model._make_data_loader(adata=adata, indices=rows, batch_size=batch_size, shuffle=False):
                u = module.inference(tensors[K.X_KEY], tensors[K.SAMPLE_KEY], use_mean=True)["u"]
                n = u.shape[0]
                z_base, eps = module.qz(u.repeat_interleave(n_sample, dim=0), all_samples.repeat(n))
                if eps.shape[-1] == 2 * module.n_latent:   # use_map=False: the mean of q(eps)
                    eps = eps[..., :module.n_latent]
                z = (z_base + eps).reshape(n, n_sample, -1)
                total += torch.cdist(z, z).sum(dim=0).double()
            out[t] = (total / len(rows)).float().cpu().numpy()
            logger.info(f"{t}: sample distances over {len(rows)} nuclei, {n_sample} samples "
                        f"({time.time() - start_time:.0f} s)")
    return names, out


def write_distances(path: Path, names: np.ndarray, distances: Dict[str, np.ndarray]) -> None:
    partial = path.with_name(path.name + ".partial")
    with h5py.File(partial, "w") as f:
        f.create_dataset("samples", data=np.asarray(names, dtype=object), dtype=h5py.string_dtype())
        group = f.create_group("celltypes")
        for t, D in distances.items():
            group.create_dataset(t, data=D.astype(np.float32), compression="gzip", compression_opts=1)
    partial.replace(path)


def read_distances(path: Path) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    with h5py.File(path, "r") as f:
        names = f["samples"].asstr()[:]
        return np.asarray(names), {t: f["celltypes"][t][:] for t in f["celltypes"]}


def train_fold(spec: Dict[str, Any], adata, fold: str, out: Path, smoke: bool) -> Path:
    import scvi
    import torch
    from scvi.external import TorchMRVI

    path = out / f"samples_{fold}.h5ad"
    scvi.settings.seed = spec["seed"] + FOLDS.index(fold)
    scvi.settings.logging_dir = str(out / f"scvi_log_{fold}")   # not the container's small tmpfs; per fold
    torch.set_float32_matmul_precision("high")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    if accelerator == "gpu":
        torch.cuda.reset_peak_memory_stats()
    logger.info(f"[{fold}] scvi-tools {scvi.__version__}, torch {torch.__version__}, accelerator {accelerator}"
                + (f" ({torch.cuda.get_device_name(0)})" if accelerator == "gpu" else ""))

    adata.obs["mrvi_sample"] = pd.Categorical(mrvi_sample_ids(adata.obs, fold))
    fold_of = _fold(adata.obs)
    held = (fold_of != "") & (fold_of != fold)
    logger.info(f"[{fold}] {adata.obs['mrvi_sample'].nunique()} MrVI samples: "
                f"{adata.obs.loc[fold_of == fold, 'person_id'].nunique()} anchor persons, "
                f"{adata.obs.loc[held, 'person_id'].nunique()} held-out persons as "
                f"{adata.obs.loc[held, 'mrvi_sample'].nunique()} dataset x donor samples; "
                f"{adata.obs['batch'].nunique()} batches")
    TorchMRVI.setup_anndata(adata, sample_key="mrvi_sample", batch_key="batch", labels_key=spec["labels_key"])
    model = TorchMRVI(adata, n_latent=spec["n_latent"], n_latent_u=spec["n_latent_u"])
    epochs = spec["smoke_epochs"] if smoke else spec["max_epochs"]
    start = time.time()
    model.train(max_epochs=epochs, accelerator=accelerator, devices=1, batch_size=spec["batch_size"],
                train_size=0.9, early_stopping=True, early_stopping_patience=spec["early_stopping_patience"],
                plan_kwargs={"lr": spec["lr"], "n_epochs_kl_warmup": spec["n_epochs_kl_warmup"]})
    elapsed = time.time() - start
    history = pd.concat([v.rename(columns=lambda c: k) for k, v in model.history.items()
                         if k.startswith(("elbo", "reconstruction_loss", "kl_local"))], axis=1)
    n_epochs = len(history)
    logger.info(f"[{fold}] trained {n_epochs} epochs in {elapsed / 60:.1f} min "
                f"({elapsed / max(n_epochs, 1):.0f} s/epoch); {_gpu_peak(torch)}")
    history.assign(model=fold).reset_index(names="epoch").to_csv(out / f"history_{fold}.tsv", sep="\t", index=False)
    try:
        model.save(str(out / f"model_{fold}"), overwrite=True)
    except Exception as e:   # the decoded samples below are what the scoring needs
        logger.warning(f"[{fold}] could not save the model: {type(e).__name__}: {e}")

    codes, groups = svi._sample_codes(adata.obs)
    samples = svi.sample_sums(adata, codes, len(groups))
    expr, lat_u, lat_z = decode_samples(model, adata, codes, len(groups), spec["decode_batch"])
    samples.layers["mrvi"] = expr.astype(np.float32)
    samples.obsm["latent_u"] = lat_u.astype(np.float32)
    samples.obsm["latent_z"] = lat_z.astype(np.float32)
    samples.obs = samples.obs.drop(columns=[c for c in samples.obs.columns
                                            if c.startswith("_") or c == "mrvi_sample"])
    settings = {k: spec[k] for k in DEFAULTS if not k.startswith(("smoke", "umap"))}
    samples.uns["mrvi"] = {"fold": fold, "epochs": n_epochs, "minutes": round(elapsed / 60, 1), "smoke": smoke,
                           **{k: ("" if v is None else v) for k, v in settings.items()}}
    samples.uns["provenance_mrvi_train"] = provenance_record(
        "mrvi_integration train", {"spec": spec["name"], "fold": fold, "scvi_tools": scvi.__version__})
    partial = path.with_name(path.name + ".partial")
    samples.write_h5ad(partial, compression="gzip")
    partial.replace(path)
    logger.info(f"[{fold}] wrote {path}: {samples.n_obs} samples x {samples.n_vars} genes; {_gpu_peak(torch)}")

    names, distances = sample_distances(model, adata, spec["distance_cells"], spec["distance_batch"],
                                        np.random.default_rng(spec["seed"]))
    write_distances(out / f"distances_{fold}.h5", names, distances)
    logger.info(f"[{fold}] wrote distances_{fold}.h5; {_gpu_peak(torch)}")
    del model
    gc.collect()
    if accelerator == "gpu":
        torch.cuda.empty_cache()
    return path


def train(spec: Dict[str, Any], folds, smoke: bool, overwrite: bool) -> None:
    out = run_dir(spec, smoke)
    for fold in folds:
        if (out / f"samples_{fold}.h5ad").exists() and not overwrite:
            raise FileExistsError(f"{out / f'samples_{fold}.h5ad'} exists; pass --overwrite to replace it")
    out.mkdir(parents=True, exist_ok=True)
    adata = load_cells(spec, smoke)
    for fold in folds:
        train_fold(spec, adata, fold, out, smoke)


# ---------------------------------------------------------------- score


def merge_folds(out: Path) -> Path:
    """samples.h5ad from samples_A.h5ad and samples_B.h5ad: layers mrvi_A,
    mrvi_B and mrvi (each bridge person from the model that held them out,
    the others the mean of the two)."""
    import anndata as ad

    path = out / "samples.h5ad"
    parts = {k: ad.read_h5ad(out / f"samples_{k}.h5ad") for k in FOLDS}
    a, b = parts["A"], parts["B"]
    if not (a.obs_names.equals(b.obs_names) and a.var_names.equals(b.var_names)):
        raise ValueError("samples_A.h5ad and samples_B.h5ad do not have the same samples and genes")
    if not (_fold(a.obs) == _fold(b.obs)).all():
        raise ValueError("samples_A.h5ad and samples_B.h5ad do not have the same folds")
    samples = ad.AnnData(a.X, obs=a.obs.copy(), var=a.var.copy())
    for k, part in parts.items():
        samples.layers[f"mrvi_{k}"] = np.asarray(part.layers["mrvi"])
        for key, value in part.obsm.items():
            samples.obsm[f"{key}_{k}"] = value
        samples.uns[f"mrvi_{k}"] = part.uns.get("mrvi", {})
    fold = _fold(samples.obs)[:, None]
    A, B = samples.layers["mrvi_A"], samples.layers["mrvi_B"]
    samples.layers["mrvi"] = np.where(fold == "B", A, np.where(fold == "A", B, (A + B) / 2)).astype(np.float32)
    samples.uns["provenance_mrvi_merge"] = provenance_record("mrvi_integration merge", {})
    samples.write_h5ad(path, compression="gzip")
    logger.info(f"Wrote {path}: {samples.n_obs} samples, layers {sorted(samples.layers)}")
    return path


def distance_stats(D: np.ndarray, names: np.ndarray, ds: np.ndarray, ids: np.ndarray,
                   persons: np.ndarray) -> Dict[Tuple[str, str], Dict[str, list]]:
    """As pp.replicate_stats, from a sample-sample distance matrix (names:
    its samples; ids: the sample of each row): per pair of datasets, for
    each person in both, the number of the other dataset's samples closer
    than their own (both directions), the number of candidates, and the
    distance to their own and (on average) to the other candidates."""
    pos = pd.Index(names).get_indexer(ids)
    if (pos < 0).any():
        raise ValueError(f"{int((pos < 0).sum())} rows have no sample in the distance matrix")
    out: Dict[Tuple[str, str], Dict[str, list]] = {}
    for a, b in combinations(list(dict.fromkeys(ds)), 2):
        ia, ib = np.flatnonzero(ds == a), np.flatnonzero(ds == b)
        pa, pb = pd.Series(pos[ia], index=persons[ia]), pd.Series(pos[ib], index=persons[ib])
        both = sorted(set(pa.index) & set(pb.index))
        if not both:
            continue
        acc = out.setdefault((a, b), {"rank": [], "n_candidates": [], "own_r": [], "other_r": []})
        for src, dst, dst_rows in ((pa, pb, ib), (pb, pa, ia)):
            candidates = np.unique(pos[dst_rows])
            query, own_pos = src[both].to_numpy(), dst[both].to_numpy()
            d = D[np.ix_(query, candidates)]
            own = D[query, own_pos]
            acc["rank"] += list((d < own[:, None]).sum(axis=1))
            acc["n_candidates"] += [len(candidates)] * len(both)
            acc["own_r"] += list(own)
            acc["other_r"] += list((d.sum(axis=1) - own) / max(len(candidates) - 1, 1))
    return out


def _merge_acc(target: Dict, acc: Dict) -> None:
    for key, v in acc.items():
        into = target.setdefault(key, {k: [] for k in v})
        for k, values in v.items():
            into[k] += values


def _fold_means(rows: List[Dict[str, Any]], method: str) -> List[Dict[str, Any]]:
    """The mean over the folds of the mixing and reliability rows (n: the
    samples for mixing, the persons of both folds for reliability)."""
    frame = pd.DataFrame(rows)
    frame = frame[frame["metric"].isin(["pc_r2", "classifier_accuracy", "cross_dataset_reliability"])]
    out = []
    for (celltype, metric, label), part in frame.groupby(["celltype", "metric", "label"], sort=False):
        n = int(part["n"].sum()) if metric == "cross_dataset_reliability" else int(part["n"].iloc[0])
        out.append(dict(celltype=celltype, method=method, metric=metric, label=label,
                        value=float(part["value"].mean()), chance=float(part["chance"].mean()), n=n))
    return out


def _disease_mean(rows: List[Dict[str, Any]], celltype: str) -> List[Dict[str, Any]]:
    folds = [r for r in rows if r["celltype"] == celltype
             and r["analysis"] in {f"pooled: mrvi_{k}, no batch term" for k in FOLDS}]
    if len(folds) != len(FOLDS):
        return []
    mean = dict(folds[0], analysis="pooled: mrvi, no batch term (mean of the two models)")
    for key in ("share_fdr05", "n_genes_fdr05", "t_correlation_with_pooled"):
        mean[key] = float(np.mean([r[key] for r in folds]))
    return [mean]


def score(spec: Dict[str, Any], out: Path) -> Path:
    """metrics.tsv, disease.tsv, folds.tsv, history.tsv and the UMAPs."""
    import anndata as ad

    pooled = spec["pooled"]
    merge_folds(out)
    samples = ad.read_h5ad(out / "samples.h5ad")
    obs = samples.obs.copy()
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(object).where(obs[c].notna(), np.nan)
    obs = obs.reset_index(drop=True)
    obs["fold"] = _fold(obs)
    counts = samples.X.toarray() if sp.issparse(samples.X) else np.asarray(samples.X)
    genes = samples.var.reset_index(drop=True)
    layers = {f"mrvi_{k}": np.asarray(samples.layers[f"mrvi_{k}"]) for k in FOLDS}
    distances = {k: read_distances(out / f"distances_{k}.h5") for k in FOLDS
                 if (out / f"distances_{k}.h5").exists()}
    numeric = [c for c in pooled["protect"] if c in pp.NUMERIC]
    rng = np.random.default_rng(spec["seed"])
    umap_spec = {**pooled, "umap_methods": spec["umap_methods"]}
    min_pairs = pooled["min_pairs"]

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
        persons = sub["person_id"].to_numpy().astype(str)
        corrected = {"none": Y, "regress": pp.correct_regress(Y, batch, sub[list(pooled["protect"])], numeric)[0]}
        models = {m: np.log2(L[sel][:, keep] * 1e6 + 1) for m, L in layers.items()}
        logger.info(f"{celltype}: {len(sub)} samples, {int(keep.sum())} of {len(keep)} HVGs "
                    f"(mean CPM >= {pooled['min_cpm']} in every dataset); held out per model: "
                    + ", ".join(f"{k} {len(set(persons[(sub['fold'] != '').to_numpy() & (sub['fold'] != k).to_numpy()]))}"
                                for k in FOLDS))
        metrics.append(dict(celltype=celltype, method="any", metric="n_genes", label="", value=int(keep.sum()),
                            chance=np.nan, n=len(sub)))
        rows, pcs = pp.score_methods(celltype, sub, corrected, pooled, rng)
        metrics += rows
        for method, M in corrected.items():
            metrics += pp.reliability_rows(celltype, method, M, ds, persons, min_pairs)

        fold_rows: List[Dict[str, Any]] = []
        held_out, anchored, by_distance = {}, {}, {}
        for k in FOLDS:
            name, M = f"mrvi_{k}", models[f"mrvi_{k}"]
            scored = scored_persons(sub, k)
            rows_k, pcs_k = pp.score_methods(celltype, sub.assign(person_id=scored), {name: M}, pooled, rng)
            pcs[name] = pcs_k[name]
            fold_rows += rows_k + pp.reliability_rows(celltype, name, M, ds, scored, min_pairs)
            _merge_acc(held_out, pp.replicate_stats(M, ds, scored))
            _merge_acc(anchored, pp.replicate_stats(M, ds, scored_persons(sub, k, anchors=True)))
            if k in distances:
                names, matrices = distances[k]
                if celltype in matrices:
                    _merge_acc(by_distance, distance_stats(matrices[celltype], names, ds,
                                                           mrvi_sample_ids(sub, k), scored))
        metrics += fold_rows + _fold_means(fold_rows, "mrvi")
        metrics += pp._replicate_rows(celltype, "mrvi", held_out, min_pairs)
        metrics += pp._replicate_rows(celltype, "mrvi_anchors", anchored, min_pairs)
        metrics += [dict(r, metric=r["metric"].replace("replicate_r_", "replicate_distance_"))
                    for r in pp._replicate_rows(celltype, "mrvi_distance", by_distance, min_pairs)]

        d, _ = pp.disease_analysis(celltype, {**corrected, **models}, sub, genes[keep].reset_index(drop=True), pooled)
        disease += d + _disease_mean(d, celltype)
        coords = {m: pp._umap(pcs[m], umap_spec) for m in spec["umap_methods"] if m in pcs}
        embeddings[celltype] = (sub, coords)

    metrics = pd.DataFrame(metrics, columns=pp.METRIC_COLUMNS)
    metrics.to_csv(out / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(disease).to_csv(out / "disease.tsv", sep="\t", index=False)
    fold_table(obs).to_csv(out / "folds.tsv", sep="\t", index=False)
    histories = [pd.read_csv(out / f"history_{k}.tsv", sep="\t") for k in FOLDS if (out / f"history_{k}.tsv").exists()]
    if histories:
        pd.concat(histories, ignore_index=True).to_csv(out / "history.tsv", sep="\t", index=False)

    figures = out / "umap"
    figures.mkdir(exist_ok=True)
    palette = pp.batch_palette(obs, list(pooled["datasets"]))
    titles = {**pp.method_titles(pooled["batch_split"]),
              "none": f"as measured ({spec['scvi']['max_cells_per_sample']} nuclei max)",
              "mrvi_A": "MrVI, model A (half A anchored)", "mrvi_B": "MrVI, model B (half B anchored)"}
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
        json.dump(provenance_record("mrvi_integration score", {"spec": spec["name"]})
                  | {"spec": {k: v for k, v in spec.items() if k not in ("pooled", "scvi")}}, f, indent=2, default=str)
    logger.info(f"Wrote metrics.tsv, disease.tsv, folds.tsv and umap/ to {out}")
    return out


def run_mrvi_integration(spec_path: Union[str, Path], stage: str, fold: Optional[str] = None,
                         smoke: bool = False, overwrite: bool = False) -> Path:
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES}")
    if fold is not None and fold not in FOLDS:
        raise ValueError(f"fold must be one of {FOLDS}")
    spec = load_spec(spec_path)
    out = run_dir(spec, smoke)
    if stage in ("train", "all"):
        train(spec, [fold] if fold else list(FOLDS), smoke, overwrite)
    if stage in ("score", "all"):
        return score(spec, out)
    return out


def main(argv: Optional[List[str]] = None) -> None:
    import argparse

    from ..utils.logging import setup_logging

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("spec")
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--fold", choices=FOLDS)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(level="INFO")
    run_mrvi_integration(args.spec, args.stage, args.fold, args.smoke, args.overwrite)


if __name__ == "__main__":
    main()
