"""Tests for the HVG comparison pipeline (seurat flavor only: no loess needed)."""

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from rosmap_processing.pipelines.hvg_compare import compare_hvgs, overlap_table


def test_overlap_table():
    t = overlap_table({"a": ["G1", "G2", "G3"], "b": ["G2", "G3", "G4", "G5"]}).iloc[0]
    assert (t["n_a"], t["n_b"], t["intersection"]) == (3, 4, 2)
    assert t["jaccard"] == 2 / 5


def test_compare_hvgs_writes_lists_and_overlaps(tmp_path):
    rng = np.random.default_rng(0)
    n_donors, per_donor, n_genes = 4, 50, 80
    rates = rng.gamma(0.5, 2.0, size=n_genes)
    # A donor-level effect on some genes, so per-donor and pooled selections differ
    donor_effect = rng.uniform(0.2, 5.0, size=(n_donors, n_genes))
    X = np.vstack([
        rng.poisson(rates * donor_effect[d] * rng.uniform(0.5, 2, size=(per_donor, 1)))
        for d in range(n_donors)
    ])
    obs = pd.DataFrame({"projid": np.repeat(np.arange(n_donors) + 100, per_donor)},
                       index=[f"c{i}" for i in range(n_donors * per_donor)])
    src = tmp_path / "in.h5ad"
    ad.AnnData(sp.csr_matrix(X.astype(np.float32)), obs=obs,
               var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])).write_h5ad(src)
    extra = tmp_path / "list.txt"
    extra.write_text("G1\nG2\nG3\n")

    out = tmp_path / "cmp"
    overlaps = compare_hvgs(src, out, n_genes=[10, 20], flavors=["seurat"],
                            extra_lists={"mine": extra}, min_genes=1, min_cells=1,
                            min_cells_per_donor=10)

    for name in ["seurat_donor_10", "seurat_pooled_10", "seurat_donor_20", "seurat_pooled_20"]:
        genes = (out / f"hvg_{name}.txt").read_text().split()
        assert len(genes) == int(name.rsplit("_", 1)[1])
    assert set(overlaps["a"]) | set(overlaps["b"]) >= {"mine", "seurat_donor_10"}
    summary = pd.read_csv(out / "hvg_overlap_summary.tsv", sep="\t")
    assert list(summary["n"]) == [10, 20]
    assert summary["intersection"].between(0, 20).all()
