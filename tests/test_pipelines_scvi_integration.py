"""The cells and score stages of the scVI integration on the synthetic
datasets of the pooled pseudobulk tests (training needs scvi-tools, which is
only in the gpu environment; here a stand-in for its decoded expression)."""

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from rosmap_processing.pipelines.scvi_integration import (CSRWriter, build_cells, load_spec, read_hvg_cells,
                                                          sample_cells, score_samples, write_h5ad_rest)
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


def test_csr_writer_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    blocks = [sp.random(n, 9, density=0.3, format="csr", random_state=rng, data_rvs=lambda k: rng.integers(1, 9, k))
              for n in (5, 0, 7, 3)]
    obs = pd.DataFrame({"dataset": pd.Categorical(["a"] * 8 + ["b"] * 7)}, index=[str(i) for i in range(15)])
    var = pd.DataFrame({"gene_id": [f"G{i}" for i in range(9)]}, index=[f"G{i}" for i in range(9)])
    path = tmp_path / "x.h5ad"
    with h5py.File(path, "w") as f:
        writer = CSRWriter(f, 9, chunk=4)
        for b in blocks:
            writer.append(b)
        write_h5ad_rest(f, obs, var, {"note": "test"})
    back = ad.read_h5ad(path)
    assert back.shape == (15, 9) and back.X.dtype == np.float32
    assert (back.X != sp.vstack(blocks)).nnz == 0
    assert list(back.obs["dataset"]) == list(obs["dataset"]) and back.uns["note"] == "test"
    backed = ad.read_h5ad(path, backed="r")
    assert (sp.csr_matrix(backed.X[6:12]) != sp.vstack(blocks).tocsr()[6:12]).nnz == 0
    backed.file.close()


def test_hvgs_per_batch_match_scanpy(tmp_path):
    try:
        import skmisc.loess  # noqa: F401  (seurat_v3; broken on Windows)
    except Exception as e:
        pytest.skip(f"skmisc.loess: {e}")
    import scanpy as sc

    spec = _spec(tmp_path)
    path = build_cells(spec)
    hvg = read_hvg_cells(path, n_top=15, chunk=7)
    cells = ad.read_h5ad(path)
    expected = sc.pp.highly_variable_genes(cells, n_top_genes=15, flavor="seurat_v3", batch_key="batch",
                                           inplace=False)
    assert list(hvg.var_names) == list(cells.var_names[expected["highly_variable"].to_numpy()])
    assert hvg.n_obs == cells.n_obs
    assert (hvg.X != cells[:, hvg.var_names].X).nnz == 0
    np.testing.assert_allclose(hvg.var["highly_variable_rank"],
                               expected.loc[expected["highly_variable"], "highly_variable_rank"])


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
