"""Tests of `pipeline umap` on the synthetic count files of the batch-effects tests."""

import json

import numpy as np
import pandas as pd
import pytest

from rosmap_processing.pipelines.umap_plots import (
    CATEGORICAL,
    MISSING_COLOR,
    OTHER_COLOR,
    build_panels,
    category_colors,
    pick_gene_set,
    run_umap,
)

from .test_pipelines_batch_effects import COVARIATES, DISEASE, _counts, _gene_files, _write_spec

ONLY_EXC = [{"column": "source_file", "match": "Excitatory", "label": "Exc"}]


def test_umap_of_one_dataset(tmp_path):
    rng = np.random.default_rng(0)
    donors = [f"R{100 + i}" for i in range(12)]
    batch_of = {d: ["B1", "B2", "B3"][i // 4] for i, d in enumerate(donors)}
    names = _counts(tmp_path / "counts.h5ad", rng, donors, batch_of)
    _gene_files(tmp_path, list(names))
    spec = _write_spec(tmp_path, "one", {"one": {
        "counts": str(tmp_path / "counts.h5ad"), "donor": "projid", "celltype": ONLY_EXC,
        "batch": ["batch"], "covariates": COVARIATES, "disease": DISEASE,
    }}, {"seurat_v3_donor_2000": str(tmp_path / "with_batch_genes.txt")})

    png = run_umap(spec, annotate=["source_file", "not_a_column"], output=tmp_path / "umap")
    assert png == tmp_path / "umap" / "umap_seurat_v3_donor_2000.png"
    assert png.stat().st_size > 10_000
    # Astrocytes match no rule but are kept (as "Other"): every sampled cell
    cells = pd.read_parquet(tmp_path / "umap" / "umap_seurat_v3_donor_2000.parquet")
    assert len(cells) == 12 * 80
    assert list(cells.columns) == ["dataset", "row", "umap1", "umap2"]
    assert np.isfinite(cells[["umap1", "umap2"]].to_numpy()).all()
    record = json.loads((tmp_path / "umap" / "umap_seurat_v3_donor_2000.json").read_text())
    assert record["n_cells"] == 12 * 80
    assert record["params"]["gene_set"] == "seurat_v3_donor_2000"
    assert not any(d in json.dumps(record) for d in donors)

    with pytest.raises(FileExistsError):
        run_umap(spec, output=tmp_path / "umap")


def test_umap_of_two_datasets(tmp_path):
    rng = np.random.default_rng(1)
    first, second = [f"A{i}" for i in range(6)], [f"B{i}" for i in range(6)]
    names = _counts(tmp_path / "a.h5ad", rng, first, {d: "B1" for d in first})
    _counts(tmp_path / "b.h5ad", rng, second, {d: "B1" for d in second}, dataset_shift=True)
    _gene_files(tmp_path, list(names))
    cfg = {"donor": "projid", "celltype": ONLY_EXC, "covariates": COVARIATES}
    spec = _write_spec(tmp_path, "two", {"a": {"counts": str(tmp_path / "a.h5ad"), **cfg},
                                         "b": {"counts": str(tmp_path / "b.h5ad"), **cfg}},
                       {"mit_seurat_v3_donor_2000": str(tmp_path / "with_batch_genes.txt")})

    run_umap(spec, output=tmp_path / "umap", max_cells_per_donor=40)
    cells = pd.read_parquet(tmp_path / "umap" / "umap_mit_seurat_v3_donor_2000.parquet")
    assert cells["dataset"].value_counts().to_dict() == {"a": 6 * 40, "b": 6 * 40}


def test_pick_gene_set():
    spec = {"gene_sets": {"lieke_hvg": "x", "seurat_v3_donor_2000": "y"}}
    assert pick_gene_set(spec, None) == "seurat_v3_donor_2000"
    assert pick_gene_set(spec, "lieke_hvg") == "lieke_hvg"
    with pytest.raises(ValueError):
        pick_gene_set(spec, "missing")
    with pytest.raises(ValueError):
        pick_gene_set({"gene_sets": {"lieke_hvg": "x"}}, None)


def test_category_colors_fold_small_levels():
    values = np.array([f"L{i}" for i in range(12) for _ in range(20 - i)] + [np.nan], dtype=object)
    colors, legend = category_colors(values)
    labels = [label for label, _ in legend]
    # The 7 largest in palette order, then the other 5 and the missing values
    assert labels == [f"L{i}" for i in range(7)] + ["Other (5)", "missing"]
    assert [c for _, c in legend[:7]] == CATEGORICAL[:7]
    assert colors[-1] == MISSING_COLOR
    assert set(colors[values == "L11"]) == {OTHER_COLOR}


def test_category_colors_order_and_other_celltype():
    _, legend = category_colors(np.array(["4", "1", "1", "2"] * 3, dtype=object))
    assert legend == [("1", CATEGORICAL[0]), ("2", CATEGORICAL[1]), ("4", CATEGORICAL[2])]
    # Text levels by size; the "Other" cell type is grey and takes no palette slot
    _, legend = category_colors(np.array(["Exc"] * 5 + ["Other"] * 3 + ["Ast"] * 2, dtype=object))
    assert legend == [("Exc", CATEGORICAL[0]), ("Other", OTHER_COLOR), ("Ast", CATEGORICAL[1])]


def test_panels_hide_donors_and_large_batches():
    n = 300
    sample = pd.DataFrame({
        "celltype": np.array(["Exc", "Ast", "Other"] * 100, dtype=object),
        "donor": np.array([f"D{i % 30}" for i in range(n)], dtype=object),
        "batch:pool": np.array([f"P{i % 25}" for i in range(n)], dtype=object),
        "sex": np.array(["0", "1"] * 150, dtype=object),
        "age": np.linspace(60, 90, n),
        "disease": np.array(["yes", "no", np.nan] * 100, dtype=object),
        "subclass": np.array([f"S{i % 40}" for i in range(n)], dtype=object),
        "pct_mito": np.array([str(x) for x in np.linspace(0, 5, n)], dtype=object),
        "total_counts": np.full(n, 1000.0),
        "n_genes": np.full(n, 500),
    })
    panels = build_panels(sample, {"sex": False, "age": True}, ["batch:pool"], ["subclass", "pct_mito"],
                          several=False)
    kinds = {p.title: p.kind for p in panels}
    assert kinds == {
        "cell type": "category", "donor (30, no legend)": "identity",
        "batch: pool (25 levels), no legend": "identity", "sex": "category", "age": "numeric",
        "disease": "category", "subclass (40 levels)": "annotation", "pct_mito": "numeric",
        "counts per cell (log10)": "numeric", "genes per cell": "numeric",
    }
