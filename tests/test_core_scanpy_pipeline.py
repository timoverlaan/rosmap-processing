"""End-to-end tests of the scanpy pipeline on small synthetic data."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from rosmap_processing.core.scanpy_pipeline import process_h5ad

N_DONORS, CELLS_PER_DONOR, N_GENES = 3, 40, 60


def _dataset(seed=0, gene_names=None, counts=True):
    rng = np.random.default_rng(seed)
    n_cells = N_DONORS * CELLS_PER_DONOR
    # Gene-specific rates so some genes are clearly more variable than others.
    rates = rng.gamma(0.5, 2.0, size=N_GENES)
    X = rng.poisson(rates * rng.uniform(0.2, 3.0, size=(n_cells, 1)))
    obs = pd.DataFrame(
        {"Donor ID": np.repeat([f"D{i}" for i in range(N_DONORS)], CELLS_PER_DONOR)},
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names or [f"G{i}" for i in range(N_GENES)])
    adata = ad.AnnData(sp.csr_matrix(X.astype(np.float32)), obs=obs, var=var)
    if not counts:
        adata.X = adata.X / 3.0
    return adata


def _gene_map(path, names, ids):
    pd.DataFrame({"var_name": names, "gene_id": ids}).to_csv(path, sep="\t", index=False)


COMMON = dict(k_neighbors=5, n_pca_components=2, min_genes=1)


def test_normalization_uses_totals_over_all_genes(tmp_path):
    src = tmp_path / "in.h5ad"
    adata = _dataset()
    adata.write_h5ad(src)
    genes = tmp_path / "genes.txt"
    genes.write_text("\n".join(["G1", "G5", "G9"]) + "\n")

    out = tmp_path / "out.h5ad"
    process_h5ad(src, out, import_genes=genes, min_cells=0, **COMMON)
    res = ad.read_h5ad(out)

    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    cells = [c.rsplit("-", 1)[0] for c in res.obs_names]
    expected = np.log1p(
        adata[cells, ["G1", "G5", "G9"]].X.toarray() / totals[adata.obs_names.get_indexer(cells), None] * 1e4
    )
    np.testing.assert_allclose(res.X.toarray(), expected, rtol=1e-5)
    np.testing.assert_array_equal(res.layers["counts"].toarray(), adata[cells, ["G1", "G5", "G9"]].X.toarray())
    assert "provenance_scanpy" in res.uns
    assert "connectivities" in res.obsp


def test_rejects_normalized_input(tmp_path):
    src = tmp_path / "in.h5ad"
    _dataset(counts=False).write_h5ad(src)
    with pytest.raises(ValueError, match="does not hold raw counts"):
        process_h5ad(src, tmp_path / "out.h5ad", **COMMON)


def _loess_works() -> bool:
    """seurat_v3 needs scikit-misc's loess; some Windows builds import but fail."""
    try:
        from skmisc.loess import loess

        x = np.linspace(0, 1, 50)
        loess(x, x ** 2).fit()
        return True
    except Exception:
        return False


def test_hvg_per_donor_writes_ranked_tables(tmp_path):
    if not _loess_works():
        pytest.skip("scikit-misc loess unavailable (needed for seurat_v3)")
    src = tmp_path / "in.h5ad"
    _dataset().write_h5ad(src)
    out = tmp_path / "out.h5ad"
    process_h5ad(src, out, n_hvgs=10, min_cells=1, **COMMON)
    assert ad.read_h5ad(out).n_vars == 10
    assert (tmp_path / "out_hvg_scores.tsv").exists()
    assert len((tmp_path / "out_hvg_names.txt").read_text().split()) > 10


def test_postnorm_hvg_within_restricted_genes(tmp_path):
    src = tmp_path / "in.h5ad"
    adata = _dataset()
    # Integer donor IDs, like projid in ROSMAP-MIT (scanpy needs a categorical batch key)
    adata.obs["Donor ID"] = np.repeat(np.arange(N_DONORS), CELLS_PER_DONOR)
    adata.write_h5ad(src)
    restrict = tmp_path / "shared.txt"
    allowed = [f"G{i}" for i in range(0, N_GENES, 2)]
    restrict.write_text("\n".join(allowed) + "\n")
    out = tmp_path / "out.h5ad"
    process_h5ad(src, out, hvg_after_normalize=8, restrict_genes=restrict, min_cells=1, **COMMON)
    res = ad.read_h5ad(out)
    assert res.n_vars == 8
    assert set(res.var_names) <= set(allowed)


def test_second_dataset_aligned_to_reference_by_gene_id(tmp_path):
    # Reference: genes G0..G59 with IDs E0..E59.
    ref_src = tmp_path / "ref_in.h5ad"
    _dataset(seed=0).write_h5ad(ref_src)
    ref_map = tmp_path / "ref_map.tsv"
    _gene_map(ref_map, [f"G{i}" for i in range(N_GENES)], [f"E{i}" for i in range(N_GENES)])
    genes = tmp_path / "genes.txt"
    genes.write_text("\n".join(["G3", "G7", "G11", "G20"]) + "\n")
    ref_out = tmp_path / "ref_out.h5ad"
    process_h5ad(ref_src, ref_out, import_genes=genes, gene_map=ref_map, **COMMON)

    # Second dataset: reversed gene order, G7 renamed (same Ensembl ID).
    names = [f"G{i}" for i in reversed(range(N_GENES))]
    names = ["G7_NEW" if n == "G7" else n for n in names]
    other = _dataset(seed=1, gene_names=names)
    other_src = tmp_path / "other_in.h5ad"
    other.write_h5ad(other_src)
    other_map = tmp_path / "other_map.tsv"
    ids = [f"E{i}" for i in reversed(range(N_GENES))]
    _gene_map(other_map, names, ids)

    other_out = tmp_path / "other_out.h5ad"
    process_h5ad(other_src, other_out, import_genes=ref_out, gene_map=other_map, **COMMON)

    ref = ad.read_h5ad(ref_out)
    res = ad.read_h5ad(other_out)
    assert list(res.var_names) == list(ref.var_names)
    assert list(res.var["gene_id"]) == list(ref.var["gene_id"])
    assert res.var.loc["G7", "symbol_dataset"] == "G7_NEW"
    cells = [c.rsplit("-", 1)[0] for c in res.obs_names]
    np.testing.assert_array_equal(
        res.layers["counts"][:, list(res.var_names).index("G7")].toarray().ravel(),
        other[cells, "G7_NEW"].X.toarray().ravel(),
    )


def test_missing_reference_genes_fail_or_zero_fill(tmp_path):
    ref_src = tmp_path / "ref_in.h5ad"
    _dataset().write_h5ad(ref_src)
    genes = tmp_path / "genes.txt"
    genes.write_text("G1\nG2\nG3\n")
    ref_out = tmp_path / "ref_out.h5ad"
    process_h5ad(ref_src, ref_out, import_genes=genes, **COMMON)

    names = [f"G{i}" for i in range(N_GENES) if i != 2]
    other = _dataset(seed=2)[:, names].copy()
    other_src = tmp_path / "other_in.h5ad"
    other.write_h5ad(other_src)

    with pytest.raises(ValueError, match="not in this dataset"):
        process_h5ad(other_src, tmp_path / "fail.h5ad", import_genes=ref_out, **COMMON)

    filled = tmp_path / "filled.h5ad"
    process_h5ad(other_src, filled, import_genes=ref_out, allow_missing_genes=True, **COMMON)
    res = ad.read_h5ad(filled)
    assert list(res.var_names) == ["G1", "G2", "G3"]
    assert res.layers["counts"][:, 1].nnz == 0
    assert res.var.loc["G2", "match"] == "missing_zero_filled"


def test_donor_column_and_cells_without_donor(tmp_path):
    adata = _dataset()
    adata.obs = adata.obs.drop(columns=["Donor ID"])
    donors = np.repeat(["R1", "R2", "R3"], CELLS_PER_DONOR).astype(object)
    donors[:5] = None  # e.g. cells whose donor could not be identified
    adata.obs["individualID"] = donors
    src = tmp_path / "in.h5ad"
    adata.write_h5ad(src)
    genes = tmp_path / "genes.txt"
    genes.write_text("G1\nG2\nG3\n")

    out = tmp_path / "out.h5ad"
    process_h5ad(src, out, import_genes=genes, donor_column="individualID", **COMMON)
    res = ad.read_h5ad(out)
    assert res.n_obs == adata.n_obs - 5
    assert set(res.obs["Donor ID"].astype(str)) == {"R1", "R2", "R3"}


def test_refuses_to_overwrite(tmp_path):
    src = tmp_path / "in.h5ad"
    _dataset().write_h5ad(src)
    out = tmp_path / "out.h5ad"
    out.write_text("existing")
    with pytest.raises(FileExistsError):
        process_h5ad(src, out, **COMMON)
