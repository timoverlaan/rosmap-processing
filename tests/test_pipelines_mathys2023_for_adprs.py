"""Unit tests for the Mathys 2023 -> ad-prs aggregation pipeline."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

from rosmap_processing.pipelines.mathys2023_for_adprs import (
    EXPECTED_MAJORS,
    _aggregate_groups,
    _normalize_lognormed,
    build_donor_metadata,
    consolidate_chunks,
    merge_chunks,
    process_class_file,
    run_pipeline,
)


def _make_adata(donors, subtypes, X):
    """Build a tiny AnnData with the obs/var columns the aggregator expects."""
    obs = pd.DataFrame(
        {"projid": donors, "cell_type_high_resolution": subtypes},
        index=[f"cell_{i}" for i in range(len(donors))],
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(X.shape[1])])
    return ad.AnnData(X=sp.csr_matrix(X.astype(np.float32)), obs=obs, var=var)


def test_aggregate_groups_matches_pandas_groupby():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 60, 8
    X = (rng.poisson(0.5, size=(n_cells, n_genes))).astype(np.float32)
    X[X < 0] = 0
    X_sp = sp.csr_matrix(X)
    X_log = X_sp.copy()
    X_log.data = np.log1p(X_log.data)

    group_codes = rng.integers(0, 4, size=n_cells)
    n_groups = 4

    mean_log, frac, n_cells_per_group = _aggregate_groups(
        X_sp, X_log, group_codes, n_groups
    )

    # Reference via pandas groupby on dense
    df = pd.DataFrame(X.copy())
    df["group"] = group_codes
    expected_n = df.groupby("group").size().reindex(range(n_groups), fill_value=0).values
    np.testing.assert_array_equal(n_cells_per_group, expected_n)

    # frac = mean of (X > 0) per group
    expected_frac = (
        (df.drop(columns="group") > 0).groupby(df["group"]).mean().reindex(range(n_groups)).values
    )
    np.testing.assert_allclose(frac, expected_frac, rtol=1e-6, atol=1e-7)

    # mean_log = mean of log1p(X) per group
    df_log = pd.DataFrame(np.log1p(X))
    df_log["group"] = group_codes
    expected_mean_log = (
        df_log.drop(columns="group").groupby(df_log["group"]).mean().reindex(range(n_groups)).values
    )
    np.testing.assert_allclose(mean_log, expected_mean_log, rtol=1e-6, atol=1e-7)


def test_aggregate_groups_handles_singleton_group():
    # One cell in its own group — mean and frac equal that cell's values.
    X = np.array([[2.0, 0.0, 5.0], [3.0, 1.0, 0.0]], dtype=np.float32)
    X_sp = sp.csr_matrix(X)
    X_log = X_sp.copy()
    X_log.data = np.log1p(X_log.data)

    mean_log, frac, n = _aggregate_groups(X_sp, X_log, np.array([0, 1]), 2)
    assert n.tolist() == [1, 1]
    np.testing.assert_allclose(mean_log[0], np.log1p(X[0]))
    np.testing.assert_allclose(frac[1], (X[1] > 0).astype(float))


def test_normalize_lognormed_matches_scanpy():
    rng = np.random.default_rng(1)
    X = rng.poisson(1.0, size=(40, 10)).astype(np.float32)
    X_sp = sp.csr_matrix(X)

    ours = _normalize_lognormed(X_sp).toarray()

    obs = pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])])
    var = pd.DataFrame(index=[f"g{i}" for i in range(X.shape[1])])
    adata = ad.AnnData(X=sp.csr_matrix(X.copy()), obs=obs, var=var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    theirs = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)

    np.testing.assert_allclose(ours, theirs, rtol=1e-5, atol=1e-6)


def test_normalize_lognormed_handles_zero_row():
    # A cell with all-zero counts must not produce NaN/inf.
    X = np.array([[0, 0, 0], [1, 2, 3]], dtype=np.float32)
    X_sp = sp.csr_matrix(X)
    out = _normalize_lognormed(X_sp).toarray()
    assert np.isfinite(out).all()
    np.testing.assert_array_equal(out[0], np.zeros(3, dtype=np.float32))


def test_process_class_file_writes_parquet_with_expected_schema(tmp_path):
    # Two donors x two subtypes, three genes. Donor "A" has cells in subtype "x",
    # donor "B" has cells in both "x" and "y".
    donors = ["A", "A", "B", "B", "B"]
    subtypes = ["x", "x", "x", "y", "y"]
    X = np.array(
        [
            [1, 0, 2],
            [3, 0, 0],
            [0, 1, 0],
            [4, 0, 1],
            [2, 2, 0],
        ],
        dtype=np.float32,
    )
    adata = _make_adata(donors, subtypes, X)
    in_path = tmp_path / "Astrocytes.h5ad"
    adata.write_h5ad(in_path)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    major_p, sub_p, cc = process_class_file(in_path, "Ast", out_dir)

    assert major_p.exists() and sub_p.exists()
    df_major = pd.read_parquet(major_p)
    df_sub = pd.read_parquet(sub_p)

    expected_cols = {
        "donor_id",
        "granularity",
        "cell_type",
        "gene_symbol",
        "mean_expr",
        "frac_expressed",
        "n_cells",
    }
    assert set(df_major.columns) == expected_cols
    assert set(df_sub.columns) == expected_cols

    # Major: 2 donors x 3 genes = 6 rows.
    assert len(df_major) == 6
    assert set(df_major["cell_type"].astype(str).unique()) == {"Ast"}
    assert set(df_major["granularity"].astype(str).unique()) == {"major"}

    # Subtype: (A,x), (B,x), (B,y) -> 3 groups x 3 genes = 9 rows.
    assert len(df_sub) == 9
    sub_groups = set(
        zip(df_sub["donor_id"].astype(str), df_sub["cell_type"].astype(str))
    )
    assert sub_groups == {("A", "x"), ("B", "x"), ("B", "y")}

    # Cell counts roll-up: 2 cells for A, 3 for B
    counts = dict(zip(cc["donor_id"], cc["n_cells_Ast"]))
    assert counts == {"A": 2, "B": 3}

    # Spot-check: donor A subtype x, gene GENE0
    # cells [1,3] -> raw counts [1,3], frac_expressed = 2/2 = 1.0
    # n_cells = 2
    row = df_sub[
        (df_sub["donor_id"].astype(str) == "A")
        & (df_sub["cell_type"].astype(str) == "x")
        & (df_sub["gene_symbol"].astype(str) == "GENE0")
    ].iloc[0]
    assert row["n_cells"] == 2
    np.testing.assert_allclose(row["frac_expressed"], 1.0, atol=1e-6)


def test_build_donor_metadata_zero_fills_missing_majors_and_merges_clinical(tmp_path):
    # Cell counts for two majors only (Ast, Mic); other 4 must be zero-filled.
    cc_ast = pd.DataFrame({"donor_id": ["1", "2"], "n_cells_Ast": [10, 5]})
    cc_mic = pd.DataFrame({"donor_id": ["1", "3"], "n_cells_Mic": [4, 7]})

    # Two class files contribute to the same major -> sum.
    cc_ex_a = pd.DataFrame({"donor_id": ["1", "2"], "n_cells_Ex": [100, 50]})
    cc_ex_b = pd.DataFrame({"donor_id": ["1", "3"], "n_cells_Ex": [25, 80]})

    clinical_csv = tmp_path / "rosmap_clinical.csv"
    pd.DataFrame(
        {
            "projid": ["1", "2", "3"],
            "msex": [1, 0, 1],
            "age_death": ["89", "92", "90+"],
            "cts_mmse30_lv": [25, 28, 12],
            "braaksc": [4, 2, 6],
            "ceradsc": [2, 4, 1],
            "cogdx": [4, 1, 4],
            "apoe_genotype": [33, 23, 34],
            "individualID": ["R1", "R2", "R3"],
        }
    ).to_csv(clinical_csv, index=False)

    df = build_donor_metadata([cc_ast, cc_mic, cc_ex_a, cc_ex_b], clinical_csv, None)

    assert sorted(df["donor_id"].tolist()) == ["1", "2", "3"]
    by_donor = df.set_index("donor_id")

    # Ex sums correctly across the two class files.
    assert by_donor.loc["1", "n_cells_Ex"] == 125
    assert by_donor.loc["3", "n_cells_Ex"] == 80
    # Donor 2 had no Ex_b entry — uses Ex_a only.
    assert by_donor.loc["2", "n_cells_Ex"] == 50

    # Missing majors -> zero, not NaN.
    for m in EXPECTED_MAJORS:
        assert by_donor[f"n_cells_{m}"].notna().all()
    assert by_donor.loc["3", "n_cells_Ast"] == 0
    assert by_donor.loc["2", "n_cells_Mic"] == 0
    assert by_donor.loc["1", "n_cells_In"] == 0

    # total_cells = sum across the 6 majors.
    assert by_donor.loc["1", "total_cells"] == 10 + 4 + 125
    assert by_donor.loc["3", "total_cells"] == 0 + 7 + 80

    # Clinical fields renamed and merged.
    assert by_donor.loc["1", "Braak"] == 4
    assert by_donor.loc["2", "CERAD"] == 4
    assert by_donor.loc["1", "sex"] == 1
    assert by_donor.loc["3", "individualID"] == "R3"


def test_normalize_projid():
    from rosmap_processing.pipelines.mathys2023_for_adprs import normalize_projid

    got = normalize_projid(["00123456", "123456", "123456.0", 123456, "10000001"])
    assert list(got) == ["123456", "123456", "123456", "123456", "10000001"]


def test_build_donor_metadata_joins_zero_padded_clinical_projid(tmp_path):
    # The clinical CSV pads projid to 8 characters; the AnnData has integers.
    cc = pd.DataFrame({"donor_id": ["123456", "20000001"], "n_cells_Ast": [3, 4]})
    clinical_csv = tmp_path / "ROSMAP_clinical.csv"
    pd.DataFrame({"projid": ["00123456", "20000001"], "braaksc": [3, 5]}).to_csv(clinical_csv, index=False)

    df = build_donor_metadata([cc], clinical_csv, None).set_index("donor_id")
    assert df.loc["123456", "Braak"] == 3
    assert df.loc["20000001", "Braak"] == 5


def _write_chunk_for_test(path, rows):
    """Write a small parquet chunk with the pipeline's schema from a list of dicts."""
    df = pd.DataFrame(rows)
    df["donor_id"] = df["donor_id"].astype(str)
    df["granularity"] = df["granularity"].astype("category")
    df["cell_type"] = df["cell_type"].astype("category")
    df["gene_symbol"] = df["gene_symbol"].astype("category")
    df["mean_expr"] = df["mean_expr"].astype(np.float32)
    df["frac_expressed"] = df["frac_expressed"].astype(np.float32)
    df["n_cells"] = df["n_cells"].astype(np.int32)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, compression="zstd", index=False)


def test_consolidate_chunks_single_chunk_passthrough(tmp_path):
    in_path = tmp_path / "in.parquet"
    _write_chunk_for_test(
        in_path,
        [
            {"donor_id": "A", "granularity": "major", "cell_type": "Ex",
             "gene_symbol": "G1", "mean_expr": 1.5, "frac_expressed": 0.5, "n_cells": 10},
        ],
    )
    out_path = tmp_path / "out.parquet"
    consolidate_chunks([in_path], out_path)
    df = pd.read_parquet(out_path)
    assert len(df) == 1
    assert df.iloc[0]["mean_expr"] == pytest.approx(1.5)
    assert df.iloc[0]["n_cells"] == 10


def test_consolidate_chunks_weighted_merge_for_duplicate_keys(tmp_path):
    # Donor A appears in both chunks for major Ex / gene G1 with different
    # n and means; B appears only in chunk2; G2 only in chunk1.
    chunk1 = tmp_path / "set1.parquet"
    chunk2 = tmp_path / "set2.parquet"
    _write_chunk_for_test(
        chunk1,
        [
            {"donor_id": "A", "granularity": "major", "cell_type": "Ex",
             "gene_symbol": "G1", "mean_expr": 2.0, "frac_expressed": 0.4, "n_cells": 100},
            {"donor_id": "A", "granularity": "major", "cell_type": "Ex",
             "gene_symbol": "G2", "mean_expr": 3.0, "frac_expressed": 0.6, "n_cells": 100},
        ],
    )
    _write_chunk_for_test(
        chunk2,
        [
            {"donor_id": "A", "granularity": "major", "cell_type": "Ex",
             "gene_symbol": "G1", "mean_expr": 4.0, "frac_expressed": 0.8, "n_cells": 50},
            {"donor_id": "B", "granularity": "major", "cell_type": "Ex",
             "gene_symbol": "G1", "mean_expr": 1.0, "frac_expressed": 0.2, "n_cells": 25},
        ],
    )

    out_path = tmp_path / "out.parquet"
    consolidate_chunks([chunk1, chunk2], out_path)
    df = pd.read_parquet(out_path)

    # 3 unique key tuples after consolidation: (A,G1), (A,G2), (B,G1)
    assert len(df) == 3
    keys = set(zip(df["donor_id"].astype(str), df["gene_symbol"].astype(str)))
    assert keys == {("A", "G1"), ("A", "G2"), ("B", "G1")}

    a_g1 = df[(df["donor_id"].astype(str) == "A") & (df["gene_symbol"].astype(str) == "G1")].iloc[0]
    # n_cells = 100 + 50 = 150
    assert a_g1["n_cells"] == 150
    # mean = (100*2 + 50*4) / 150 = 400/150 = 2.6667
    assert a_g1["mean_expr"] == pytest.approx(400 / 150, rel=1e-5)
    # frac = (100*0.4 + 50*0.8) / 150 = 80/150
    assert a_g1["frac_expressed"] == pytest.approx(80 / 150, rel=1e-5)

    a_g2 = df[(df["donor_id"].astype(str) == "A") & (df["gene_symbol"].astype(str) == "G2")].iloc[0]
    assert a_g2["n_cells"] == 100
    assert a_g2["mean_expr"] == pytest.approx(3.0)

    b_g1 = df[df["donor_id"].astype(str) == "B"].iloc[0]
    assert b_g1["n_cells"] == 25
    assert b_g1["mean_expr"] == pytest.approx(1.0)


def test_run_pipeline_validates_clinical_csv_up_front(tmp_path):
    """A bad --clinical-csv must fail before the (expensive) aggregation runs."""
    rng = np.random.default_rng(11)
    input_dir = tmp_path / "raw"
    input_dir.mkdir()
    n_cells, n_genes = 6, 3
    donors = rng.choice(["d1", "d2"], size=n_cells).tolist()
    subtypes = ["Ast.1"] * n_cells
    X = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)
    _make_adata(donors, subtypes, X).write_h5ad(input_dir / "Astrocytes.h5ad")

    bad_clinical = tmp_path / "does_not_exist.csv"
    with pytest.raises(FileNotFoundError, match="case-sensitive"):
        run_pipeline(
            input_dir=input_dir,
            output_dir=tmp_path / "out",
            clinical_csv=bad_clinical,
            mit_metadata_csv=None,
            keep_intermediate=False,
        )

    # No aggregation output should exist.
    assert not (tmp_path / "out" / "per_donor_celltype_expression.parquet").exists()


def test_run_pipeline_resume_rebuilds_metadata_from_existing_parquet(tmp_path):
    """--resume must skip aggregation and rebuild only metadata + README."""
    rng = np.random.default_rng(13)
    input_dir = tmp_path / "raw"
    input_dir.mkdir()
    out_dir = tmp_path / "out"
    n_cells, n_genes = 10, 4

    donors_ast = rng.choice(["d1", "d2", "d3"], size=n_cells).tolist()
    _make_adata(
        donors_ast,
        ["Ast.1"] * n_cells,
        rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32),
    ).write_h5ad(input_dir / "Astrocytes.h5ad")

    donors_mic = rng.choice(["d1", "d2"], size=n_cells).tolist()
    _make_adata(
        donors_mic,
        ["Mic.1"] * n_cells,
        rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32),
    ).write_h5ad(input_dir / "Immune_cells.h5ad")

    clinical_csv = tmp_path / "clinical.csv"
    pd.DataFrame(
        {
            "projid": ["d1", "d2", "d3"],
            "msex": [1, 0, 1],
            "age_death": ["80", "85", "90"],
            "cts_mmse30_lv": [27, 24, 22],
            "braaksc": [2, 4, 6],
            "ceradsc": [3, 2, 1],
            "cogdx": [1, 4, 4],
        }
    ).to_csv(clinical_csv, index=False)

    # First, a regular run that produces the parquet + metadata + README.
    run_pipeline(
        input_dir=input_dir,
        output_dir=out_dir,
        clinical_csv=clinical_csv,
        mit_metadata_csv=None,
        keep_intermediate=False,
    )
    parquet = out_dir / "per_donor_celltype_expression.parquet"
    meta = out_dir / "donor_metadata.tsv"
    assert parquet.exists() and meta.exists()
    parquet_mtime_before = parquet.stat().st_mtime_ns
    meta_first = pd.read_csv(meta, sep="\t")

    # Wipe metadata + README to simulate the post-aggregation crash, leave the parquet.
    meta.unlink()
    (out_dir / "README.md").unlink()

    # Resume: must rebuild metadata + README using counts from the existing parquet
    # without recomputing the parquet itself.
    run_pipeline(
        input_dir=input_dir,
        output_dir=out_dir,
        clinical_csv=clinical_csv,
        mit_metadata_csv=None,
        keep_intermediate=False,
        resume=True,
    )
    assert meta.exists() and (out_dir / "README.md").exists()
    # Parquet wasn't touched.
    assert parquet.stat().st_mtime_ns == parquet_mtime_before

    # Metadata content matches the original run on the cell-count + total columns.
    meta_resume = pd.read_csv(meta, sep="\t")
    count_cols = [c for c in meta_first.columns if c.startswith("n_cells_") or c == "total_cells"]
    pd.testing.assert_frame_equal(
        meta_first[["donor_id"] + count_cols].sort_values("donor_id").reset_index(drop=True),
        meta_resume[["donor_id"] + count_cols].sort_values("donor_id").reset_index(drop=True),
    )


def test_run_pipeline_consolidates_excitatory_sets(tmp_path):
    """End-to-end: two Ex source files sharing a donor must produce one merged row."""
    rng = np.random.default_rng(7)
    input_dir = tmp_path / "raw"
    input_dir.mkdir()
    out_dir = tmp_path / "out"
    n_genes = 4

    # Both files include donor "shared"; only set2 has donor "only_set2".
    for stem, donors in [
        ("Excitatory_neurons_set1", ["shared", "shared", "shared"]),
        ("Excitatory_neurons_set2", ["shared", "only_set2", "only_set2"]),
    ]:
        subtypes = ["Exc.A"] * len(donors)
        X = rng.poisson(2.0, size=(len(donors), n_genes)).astype(np.float32)
        adata = _make_adata(donors, subtypes, X)
        adata.write_h5ad(input_dir / f"{stem}.h5ad")

    clinical_csv = tmp_path / "clinical.csv"
    pd.DataFrame(
        {
            "projid": ["shared", "only_set2"],
            "msex": [1, 0],
            "age_death": ["80", "85"],
            "cts_mmse30_lv": [27, 24],
            "braaksc": [2, 4],
            "ceradsc": [3, 2],
            "cogdx": [1, 4],
        }
    ).to_csv(clinical_csv, index=False)

    run_pipeline(
        input_dir=input_dir,
        output_dir=out_dir,
        clinical_csv=clinical_csv,
        mit_metadata_csv=None,
        keep_intermediate=False,
    )

    parquet_df = pd.read_parquet(out_dir / "per_donor_celltype_expression.parquet")

    # No duplicate (donor, granularity, cell_type, gene_symbol) tuples post-merge.
    key_cols = ["donor_id", "granularity", "cell_type", "gene_symbol"]
    for col in key_cols:
        parquet_df[col] = parquet_df[col].astype(str)
    dup_count = parquet_df.duplicated(subset=key_cols).sum()
    assert dup_count == 0, f"found {dup_count} duplicate key tuples"

    # The 'shared' donor has cells in both files; n_cells at major Ex = 3 + 1 = 4.
    shared_major = parquet_df[
        (parquet_df["donor_id"] == "shared")
        & (parquet_df["granularity"] == "major")
        & (parquet_df["cell_type"] == "Ex")
    ]
    assert len(shared_major) == n_genes
    assert (shared_major["n_cells"] == 4).all()

    # The 'only_set2' donor sits in only one file; n_cells at major Ex = 2.
    only_major = parquet_df[
        (parquet_df["donor_id"] == "only_set2")
        & (parquet_df["granularity"] == "major")
        & (parquet_df["cell_type"] == "Ex")
    ]
    assert len(only_major) == n_genes
    assert (only_major["n_cells"] == 2).all()

    # donor_metadata Ex count for 'shared' is also 4 (build_donor_metadata sums per-major).
    donors_df = pd.read_csv(out_dir / "donor_metadata.tsv", sep="\t")
    by_donor = donors_df.set_index("donor_id")
    assert by_donor.loc["shared", "n_cells_Ex"] == 4
    assert by_donor.loc["only_set2", "n_cells_Ex"] == 2


def test_run_pipeline_end_to_end(tmp_path):
    rng = np.random.default_rng(42)
    input_dir = tmp_path / "raw"
    input_dir.mkdir()
    out_dir = tmp_path / "out"

    # Tiny synthetic atlas: two class files (Astrocytes + Immune_cells), shared donors.
    n_cells, n_genes = 30, 5
    for stem, subtypes in [("Astrocytes", ["Ast.1", "Ast.2"]), ("Immune_cells", ["Mic.1"])]:
        donors = rng.choice(["d1", "d2", "d3"], size=n_cells).tolist()
        ct = rng.choice(subtypes, size=n_cells).tolist()
        X = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)
        adata = _make_adata(donors, ct, X)
        adata.write_h5ad(input_dir / f"{stem}.h5ad")

    clinical_csv = tmp_path / "clinical.csv"
    pd.DataFrame(
        {
            "projid": ["d1", "d2", "d3"],
            "msex": [1, 0, 1],
            "age_death": ["80", "85", "90"],
            "cts_mmse30_lv": [27, 24, 22],
            "braaksc": [2, 4, 6],
            "ceradsc": [3, 2, 1],
            "cogdx": [1, 4, 4],
            "individualID": ["R1", "R2", "R3"],
        }
    ).to_csv(clinical_csv, index=False)

    run_pipeline(
        input_dir=input_dir,
        output_dir=out_dir,
        clinical_csv=clinical_csv,
        mit_metadata_csv=None,
        keep_intermediate=False,
    )

    # Outputs exist
    assert (out_dir / "per_donor_celltype_expression.parquet").exists()
    assert (out_dir / "donor_metadata.tsv").exists()
    assert (out_dir / "README.md").exists()
    # Intermediate cleaned up
    assert not (out_dir / "_intermediate").exists()

    # Parquet loads cleanly and has expected granularities + cell_types.
    parquet_df = pd.read_parquet(out_dir / "per_donor_celltype_expression.parquet")
    assert set(parquet_df["granularity"].astype(str).unique()) == {"major", "subtype"}
    majors_emitted = set(
        parquet_df.loc[parquet_df["granularity"] == "major", "cell_type"]
        .astype(str)
        .unique()
    )
    assert majors_emitted == {"Ast", "Mic"}
    sub_emitted = set(
        parquet_df.loc[parquet_df["granularity"] == "subtype", "cell_type"]
        .astype(str)
        .unique()
    )
    assert sub_emitted.issubset({"Ast.1", "Ast.2", "Mic.1"})

    # Donor metadata: 3 donors, all 6 majors present (zeros for absent ones).
    donors_df = pd.read_csv(out_dir / "donor_metadata.tsv", sep="\t")
    assert len(donors_df) == 3
    for m in EXPECTED_MAJORS:
        assert f"n_cells_{m}" in donors_df.columns
    # Only Ast and Mic have any cells; the rest must be zero.
    assert donors_df["n_cells_In"].sum() == 0
    assert donors_df["n_cells_Ex"].sum() == 0
    assert donors_df["n_cells_Oli"].sum() == 0
    assert donors_df["n_cells_Opc"].sum() == 0
    # total_cells should equal Ast + Mic for each donor.
    np.testing.assert_array_equal(
        donors_df["total_cells"].values,
        (donors_df["n_cells_Ast"] + donors_df["n_cells_Mic"]).values,
    )
