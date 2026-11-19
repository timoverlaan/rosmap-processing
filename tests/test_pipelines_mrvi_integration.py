"""The fold assignment, MrVI sample IDs and scoring of the MrVI integration,
on synthetic obs and on the synthetic datasets of the pooled pseudobulk tests
(training needs scvi-tools, which is only in the gpu environment; here a
stand-in for its decoded expression and sample distances)."""

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from rosmap_processing.pipelines import scvi_integration as svi
from rosmap_processing.pipelines.mrvi_integration import (FOLDS, assign_folds, distance_stats, fold_table,
                                                          load_cells, load_spec, mrvi_sample_ids, score,
                                                          scored_persons, smoke_rows, with_folds, write_distances)
from .test_pipelines_pooled_pseudobulk import PERSONS
from .test_pipelines_scvi_integration import _spec as _scvi_spec


def _obs():
    """Rows per dataset x person: X1-X9 in a+b, Y1-Y3 in a+b+c, Z1-Z4 in b+c,
    and persons in one dataset."""
    rows = []
    for p in [f"X{i}" for i in range(9)]:
        rows += [("a", p), ("b", p)]
    for p in [f"Y{i}" for i in range(3)]:
        rows += [("a", p), ("b", p), ("c", p)]
    for p in [f"Z{i}" for i in range(4)]:
        rows += [("b", p), ("c", p)]
    rows += [("a", f"S{i}") for i in range(5)] + [("c", f"T{i}") for i in range(5)]
    obs = pd.DataFrame(rows * 2, columns=["dataset", "person_id"])   # two cell types
    obs["celltype"] = ["Exc"] * len(rows) + ["Ast"] * len(rows)
    obs["donor"] = obs["dataset"] + "_" + obs["person_id"].str.lower()
    return obs


def test_folds_split_each_combination_of_datasets():
    obs = _obs()
    folds = assign_folds(obs, seed=0)
    assert set(folds.index) == {f"X{i}" for i in range(9)} | {f"Y{i}" for i in range(3)} | {f"Z{i}" for i in range(4)}
    assert set(folds) == {"A", "B"}
    for prefix, n in (("X", 9), ("Y", 3), ("Z", 4)):
        counts = folds[folds.index.str.startswith(prefix)].value_counts()
        assert counts.sum() == n and abs(counts.get("A", 0) - counts.get("B", 0)) <= 1
    assert abs((folds == "A").sum() - (folds == "B").sum()) <= 1
    assert folds.equals(assign_folds(obs.sample(frac=1, random_state=1), seed=0))   # row order does not matter
    assert not all(assign_folds(obs, seed=s).equals(folds) for s in range(1, 6))


def test_sample_ids_anchor_one_half_and_hold_out_the_other():
    obs = with_folds(_obs(), seed=0)
    fold = obs["fold"].to_numpy()
    for model, other in (("A", "B"), ("B", "A")):
        ids = mrvi_sample_ids(obs, model)
        per_person = pd.Series(ids).groupby(obs["person_id"].to_numpy()).nunique()
        n_datasets = obs.groupby("person_id")["dataset"].nunique()
        anchored = obs.loc[fold == model, "person_id"].unique()
        held = obs.loc[fold == other, "person_id"].unique()
        assert (per_person[anchored] == 1).all()                        # one sample over the datasets
        assert (per_person[held] == n_datasets[held]).all()             # one per dataset
        assert (per_person[obs.loc[fold == "", "person_id"].unique()] == 1).all()
        assert all(i.startswith("d:") for i in ids[fold == other]) and all(i.startswith("p:") for i in ids[fold != other])
        # the replicate metrics of this model only pair the persons it held out
        scored = scored_persons(obs, model)
        pairs = {p for p in set(scored) if obs.loc[scored == p, "dataset"].nunique() > 1}
        assert pairs == set(held)
        circular = scored_persons(obs, model, anchors=True)
        assert {p for p in set(circular) if obs.loc[circular == p, "dataset"].nunique() > 1} == set(anchored)


def test_smoke_rows_take_all_rows_of_the_chosen_persons():
    obs = _obs()
    rows = smoke_rows(obs, n_bridge=4, n_other=2, seed=0)
    chosen = obs.iloc[rows]
    assert chosen["person_id"].nunique() == 6 and (np.diff(rows) > 0).all()
    assert len(chosen) == obs["person_id"].isin(set(chosen["person_id"])).sum()
    sets = chosen.groupby("person_id")["dataset"].nunique()
    assert (sets > 1).sum() == 4


def test_distance_stats_ranks_the_own_sample():
    # samples: persons p0-p3 in a and b (held out: one sample each), q0-q2 only in b
    names = np.array([f"a{i}" for i in range(4)] + [f"b{i}" for i in range(4)] + ["q0", "q1", "q2"])
    rng = np.random.default_rng(0)
    points = rng.normal(size=(4, 5)) * 10
    coords = np.vstack([points, points + 0.1, rng.normal(size=(3, 5)) * 10])
    coords[4 + 3] += 100   # p3's sample in b is far from everything
    D = np.linalg.norm(coords[:, None] - coords[None], axis=2)
    ds = np.array(["a"] * 4 + ["b"] * 7)
    ids = names.copy()
    persons = np.array([f"p{i}" for i in range(4)] * 2 + ["q0", "q1", "q2"])
    acc = distance_stats(D, names, ds, ids, persons)
    v = acc[("a", "b")]
    assert v["n_candidates"] == [7] * 4 + [4] * 4
    assert v["rank"][:3] == [0, 0, 0] and v["rank"][3] == 6      # a -> b: p3's own sample is the farthest
    assert v["rank"][4:7] == [0, 0, 0] and v["rank"][7] == (D[7, :4] < D[7, 3]).sum()   # b -> a
    assert np.allclose(v["own_r"][:3], np.sqrt(5 * 0.01))


def _mrvi_spec(tmp_path):
    scvi = _scvi_spec(tmp_path)
    path = tmp_path / "mrvi.yaml"
    path.write_text(yaml.safe_dump({"scvi_spec": str(tmp_path / "scvi.yaml"),
                                    "output_dir": str(tmp_path / "scvi" / "mrvi"),
                                    "umap_methods": ["none", "mrvi_A"]}))
    return scvi, load_spec(path)


def test_load_cells_reads_the_cached_hvgs(tmp_path):
    scvi, spec = _mrvi_spec(tmp_path)
    cells = ad.read_h5ad(svi.build_cells(scvi))
    genes = list(cells.var_names[::3])
    cache = tmp_path / "scvi" / "mrvi" / "hvg_3000.tsv"
    cache.parent.mkdir(parents=True)
    pd.DataFrame({"highly_variable_rank": np.arange(len(genes), dtype=float)}, index=genes).to_csv(cache, sep="\t")
    adata = load_cells(spec, smoke=False)
    assert list(adata.var_names) == genes and adata.n_obs == cells.n_obs
    assert (adata.X != cells[:, genes].X).nnz == 0
    assert set(adata.obs["fold"]) == {"A", "B", ""}
    assert adata.obs.loc[adata.obs["fold"] != "", "person_id"].nunique() == 15   # P15-P29 are in a and b
    small = load_cells({**spec, "smoke_bridge_persons": 4, "smoke_other_persons": 3}, smoke=True)
    assert small.obs["person_id"].nunique() == 7 and small.n_obs < cells.n_obs
    assert (small.X != cells[small.obs_names, genes].X).nnz == 0
    assert set(small.obs["batch"].cat.categories) == set(small.obs["batch"])   # no empty batches


def test_score_stage(tmp_path):
    scvi, spec = _mrvi_spec(tmp_path)
    cells = ad.read_h5ad(svi.build_cells(scvi))
    cells.obs = with_folds(cells.obs, seed=0)
    codes, groups = svi._sample_codes(cells.obs)
    out = tmp_path / "scvi" / "mrvi"
    out.mkdir(parents=True)
    rng = np.random.default_rng(1)
    for k in FOLDS:   # stand-in for train_fold: proportions with a model-specific offset per gene
        samples = svi.sample_sums(cells, codes, len(groups))
        props = np.asarray(samples.X.todense(), dtype=np.float64) + 1
        props *= rng.lognormal(0, 0.1, size=props.shape[1])
        samples.layers["mrvi"] = (props / props.sum(axis=1, keepdims=True)).astype(np.float32)
        samples.write_h5ad(out / f"samples_{k}.h5ad")
        # distances: one point per MrVI sample, a held-out person's samples near their person's point
        ids = mrvi_sample_ids(samples.obs, k)
        persons = samples.obs["person_id"].astype(str).to_numpy()
        names = np.unique(ids)
        base = {p: rng.normal(size=4) * 10 for p in set(persons)}
        point = {i: base[p] + rng.normal(size=4) * 0.1 for i, p in zip(ids, persons)}
        coords = np.array([point[n] for n in names])
        D = np.linalg.norm(coords[:, None] - coords[None], axis=2)
        write_distances(out / f"distances_{k}.h5", names, {"Exc": D, "Ast": D})
        pd.DataFrame({"epoch": [0, 1], "elbo_train": [2.0, 1.0], "model": k}).to_csv(out / f"history_{k}.tsv",
                                                                                     sep="\t", index=False)
    score(spec, out)

    merged = ad.read_h5ad(out / "samples.h5ad")
    fold = merged.obs["fold"].astype(str).to_numpy()
    for f, source in (("A", "mrvi_B"), ("B", "mrvi_A")):   # each bridge person from the model that held them out
        assert np.array_equal(merged.layers["mrvi"][fold == f], merged.layers[source][fold == f])
    m = pd.read_csv(out / "metrics.tsv", sep="\t")
    assert {"none", "regress", "mrvi_A", "mrvi_B", "mrvi", "mrvi_anchors", "mrvi_distance"} <= set(m["method"])
    top1 = m[(m["metric"] == "replicate_top1") & (m["label"] == "a vs b")].set_index(["celltype", "method"])
    held = {k: int(((fold != "") & (fold != k) & (merged.obs["celltype"] == "Exc").to_numpy()).sum() // 2)
            for k in FOLDS}
    assert top1.loc[("Exc", "none"), "n"] == top1.loc[("Exc", "mrvi"), "n"] == 15 == sum(held.values())
    for k in FOLDS:
        assert top1.loc[("Exc", f"mrvi_{k}"), "n"] == held[k]
    assert top1.loc[("Exc", "mrvi_distance"), "value"] == 1.0
    assert "replicate_distance_own" in set(m.loc[m["method"] == "mrvi_distance", "metric"])
    pc = m[(m["metric"] == "pc_r2") & (m["label"] == "dataset")].set_index(["celltype", "method"])["value"]
    assert np.isclose(pc[("Exc", "mrvi")], (pc[("Exc", "mrvi_A")] + pc[("Exc", "mrvi_B")]) / 2)
    rel = m[(m["metric"] == "cross_dataset_reliability") & (m["celltype"] == "Exc")].set_index("method")
    assert rel.loc["mrvi", "n"] == 15
    disease = pd.read_csv(out / "disease.tsv", sep="\t")
    assert "pooled: mrvi, no batch term (mean of the two models)" in set(disease["analysis"])
    folds = pd.read_csv(out / "folds.tsv", sep="\t").set_index("datasets")
    assert folds.loc["a+b", ["A", "B"]].sum() == 15 and folds.loc["a+b", "not split"] == 0
    assert (out / "umap" / "overview.png").exists() and len(pd.read_csv(out / "history.tsv", sep="\t")) == 4
    persons = {p for ps in PERSONS.values() for p in ps}
    for name in ("metrics.tsv", "disease.tsv", "folds.tsv"):
        text = (out / name).read_text()
        assert not any(p in text for p in persons if len(p) > 2) and "a3" not in text


def test_fold_table_counts_persons_per_combination_of_datasets():
    table = fold_table(with_folds(_obs(), seed=0)).set_index("datasets")
    assert set(table.index) == {"a", "c", "a+b", "a+b+c", "b+c"}
    assert table.loc["a+b", ["A", "B"]].sum() == 9 and table.loc["a", "not split"] == 5
    assert table.loc["a+b", "not split"] == 0 and table.loc["a", ["A", "B"]].sum() == 0
