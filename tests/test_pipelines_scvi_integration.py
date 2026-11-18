"""The cells and score stages of the scVI integration on the synthetic
datasets of the pooled pseudobulk tests (training needs scvi-tools, which is
only in the gpu environment; here a stand-in for its decoded expression)."""

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from rosmap_processing.pipelines.scvi_integration import build_cells, load_spec, sample_cells, score_samples
from .test_pipelines_pooled_pseudobulk import N_GENES, PERSONS, _setup


def _spec(tmp_path, max_cells=12):
    pooled = _setup(tmp_path, split=True)
    path = tmp_path / "scvi.yaml"
    path.write_text(yaml.safe_dump({"pooled_spec": str(pooled), "output_dir": str(tmp_path / "scvi"),
                                    "max_cells_per_sample": max_cells, "umap_methods": ["none", "scvi"]}))
    return load_spec(path)


def test_sample_cells_caps_each_sample():
    cells = pd.DataFrame({"row": np.arange(50), "donor": ["d1"] * 30 + ["d2"] * 20,
                          "celltype": (["Exc", "Ast"] * 25)})
    out = sample_cells(cells, 8, np.random.default_rng(0))
    assert out.groupby(["donor", "celltype"]).size().max() == 8
    assert (np.diff(out["row"]) > 0).all()   # file order


def test_cells_and_score_stages(tmp_path):
    spec = _spec(tmp_path)
    cells = ad.read_h5ad(build_cells(spec))
    n_samples = 2 * sum(len(p) for p in PERSONS.values())
    assert cells.n_obs == 12 * n_samples and cells.n_vars == N_GENES - 1
    assert set(cells.obs["batch"]) == {"a", "b:Rush", "b:Sinai", "c"}
    assert cells.obs["sample_nuclei"].min() == 30 and cells.obs["library_size"].min() > 0

    # stand-in for train_models: per-sample sums and, as the "scvi" layer, proportions
    keys = cells.obs[["dataset", "donor", "celltype"]].astype(str)
    codes, groups = pd.MultiIndex.from_frame(keys).factorize()
    indicator = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))))
    summed = np.asarray((indicator @ cells.X).todense())
    first = pd.Series(np.arange(cells.n_obs)).groupby(codes).first().to_numpy()
    obs = cells.obs.iloc[first].reset_index(drop=True)
    obs["n_cells"] = np.bincount(codes)
    obs["library_size"] = np.bincount(codes, weights=cells.obs["library_size"].to_numpy())
    samples = ad.AnnData(sp.csr_matrix(summed.astype(np.int32)), obs=obs, var=cells.var.copy())
    samples.layers["scvi"] = (summed / summed.sum(axis=1, keepdims=True)).astype(np.float32)
    samples.obs.index = samples.obs.index.astype(str)
    samples.write_h5ad(tmp_path / "scvi" / "samples.h5ad")

    out = score_samples(spec)
    m = pd.read_csv(out / "metrics.tsv", sep="\t")
    assert {"none", "regress", "scvi"} <= set(m["method"])
    pc = m[(m["metric"] == "pc_r2") & (m["label"] == "dataset")].set_index(["celltype", "method"])["value"]
    assert pc[("Exc", "none")] > 0.3 and pc[("Exc", "regress")] < 0.05
    rel = m[m["metric"] == "cross_dataset_reliability"]
    assert set(rel["method"]) >= {"none", "regress", "scvi"}
    disease = pd.read_csv(out / "disease.tsv", sep="\t")
    assert "pooled: scvi, no batch term" in set(disease["analysis"])
    assert (out / "umap" / "overview.png").exists() and (out / "umap" / "Exc.png").exists()
    for name in ("metrics.tsv", "disease.tsv"):
        text = (out / name).read_text()
        assert "P15" not in text and "a3" not in text
