"""Tests of the pooled pseudobulk on synthetic count files: three datasets
with known per-gene dataset offsets, persons measured twice, and a disease
effect."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from rosmap_processing.pipelines.pooled_pseudobulk import (
    correct_anchor,
    correct_center,
    correct_regress,
    correct_ruv,
    prps_groups,
    ruv_apply,
    ruv_factors,
    run_pooled_pseudobulk,
)

N_GENES = 60
DISEASE_GENES = list(range(0, 8))     # x2 in persons with dementia
BANK_GENES = list(range(50, 56))       # x3 in the b persons of brain bank "Sinai"
OFFSET_GENES = {"a": [], "b": list(range(20, 30)), "c": list(range(30, 40))}   # x4 in b, x3 in c
# a: P0-P29; b: P15-P44 (P15-P29 in both); c: its own S0-S14
PERSONS = {"a": [f"P{i}" for i in range(0, 30)], "b": [f"P{i}" for i in range(15, 45)],
           "c": [f"S{i}" for i in range(15)]}
SINAI = frozenset(f"P{i}" for i in range(30, 45))   # only in b


def _person_effects(rng):
    everyone = sorted({p for ps in PERSONS.values() for p in ps})
    return {p: rng.lognormal(0, 0.3, size=N_GENES) for p in everyone}


def _dementia(person):
    return "yes" if int(person[1:]) % 2 == 0 else "no"


def _write_counts(path, rng, name, persons, effects, rates, cells=30, bank=frozenset()):
    blocks, obs = [], []
    for k, person in enumerate(persons):
        for celltype in ("Exc", "Ast"):
            r = np.tile(rates * effects[person], (cells, 1))
            if celltype == "Ast":
                r[:, 40:50] *= 5
            if _dementia(person) == "yes":
                r[:, DISEASE_GENES] *= 2
            r[:, OFFSET_GENES[name]] *= {"a": 1, "b": 4, "c": 3}[name]
            if person in bank:
                r[:, BANK_GENES] *= 3
            blocks.append(rng.poisson(r * rng.uniform(0.5, 2.0, size=(cells, 1))))
            obs.append(pd.DataFrame({"donor_key": f"{name}{k}", "class": celltype},
                                    index=[f"{name}{k}_{celltype}_{i}" for i in range(cells)]))
    obs = pd.concat(obs)
    obs["donor_key"] = obs["donor_key"].astype("category")
    obs["class"] = obs["class"].astype("category")
    var = pd.DataFrame(index=[f"{name}_G{i}" for i in range(N_GENES)])
    ad.AnnData(sp.csr_matrix(np.vstack(blocks).astype(np.float32)), obs=obs, var=var).write_h5ad(path)
    # gene map: the same Ensembl IDs under each dataset's own names; c lacks the last gene
    ids = [f"ENSG{i:011d}" for i in range(N_GENES)]
    if name == "c":
        ids[-1] = None
    pd.DataFrame({"var_name": var.index, "gene_id": ids, "ref_name": var.index, "method": "symbol"}) \
        .to_csv(path.with_suffix(".map.tsv"), sep="\t", index=False)


def _setup(tmp_path, split=False):
    rng = np.random.default_rng(3)
    effects = _person_effects(rng)
    rates = rng.gamma(2.0, 2.0, size=N_GENES) + 1
    rows = {}
    for name, persons in PERSONS.items():
        _write_counts(tmp_path / f"{name}.h5ad", rng, name, persons, effects, rates,
                      bank=SINAI if split and name == "b" else frozenset())
        for k, p in enumerate(persons):
            rows.setdefault(p, {"person_id": p})[f"{name}_donor"] = f"{name}{k}"
    table = pd.DataFrame(list(rows.values()))
    table["sex"] = ["female" if i % 3 else "male" for i in range(len(table))]
    table["age_capped"] = 70 + np.arange(len(table)) % 20
    table["dementia"] = table["person_id"].map(_dementia)
    table["n_datasets"] = table[[f"{k}_donor" for k in PERSONS]].notna().sum(axis=1)
    table["brain_bank"] = np.where(table["person_id"].isin(SINAI), "Sinai", "Rush")
    table.to_csv(tmp_path / "donor_table.tsv", sep="\t", index=False)
    celltype = [{"column": "class", "match": "^Exc$", "label": "Exc"},
                {"column": "class", "match": "^Ast$", "label": "Ast"}]
    spec = {"output_dir": str(tmp_path / "out"), "donor_table": str(tmp_path / "donor_table.tsv"),
            "min_genes": 1, "min_cells": 5, "min_cpm": 1.0, "n_pcs": 10, "chunk_size": 500, "n_folds": 3,
            "ruv_prps_min": 2,
            "umap": split, "batch_split": {"b": "brain_bank"} if split else {},
            "umap_methods": ["none", "regress", "anchor", "ruv"],
            "datasets": {name: {"counts": str(tmp_path / f"{name}.h5ad"),
                                "gene_map": str(tmp_path / f"{name}.map.tsv"),
                                "donor": "donor_key", "celltype": celltype} for name in PERSONS}}
    path = tmp_path / "spec.yaml"
    path.write_text(yaml.safe_dump(spec))
    return path


def test_pooled_pseudobulk_end_to_end(tmp_path):
    out = run_pooled_pseudobulk(_setup(tmp_path))
    adata = ad.read_h5ad(out / "pseudobulk.h5ad")
    assert adata.n_vars == N_GENES - 1                      # c lacks one gene
    assert adata.n_obs == 2 * sum(len(p) for p in PERSONS.values())
    assert adata.obs["n_cells"].min() == 30
    assert set(adata.obs.columns) >= {"dataset", "donor", "person_id", "celltype", "library_size", "dementia"}
    # the counts are the sums over the donor's cells
    a = ad.read_h5ad(tmp_path / "a.h5ad")
    cells = (a.obs["donor_key"] == "a3").to_numpy() & (a.obs["class"] == "Exc").to_numpy()
    row = adata.obs_names.get_loc("a:a3:Exc")
    assert np.allclose(adata.X[row].toarray().ravel(), np.asarray(a.X[cells].sum(axis=0)).ravel()[:N_GENES - 1])

    m = pd.read_csv(out / "metrics.tsv", sep="\t")
    pc = m[(m["metric"] == "pc_r2") & (m["label"] == "dataset")].set_index(["celltype", "method"])["value"]
    for t in ("Exc", "Ast"):
        assert pc[(t, "none")] > 0.3
        assert pc[(t, "center")] < 0.05 and pc[(t, "anchor")] < 0.05 and pc[(t, "regress")] < 0.05
    top1 = m[(m["metric"] == "replicate_top1") & (m["label"] == "a vs b")].set_index(["celltype", "method"])
    assert (top1["n"] == 15).all()
    for t in ("Exc", "Ast"):
        assert top1.loc[(t, "anchor"), "value"] >= 0.8 and top1.loc[(t, "center"), "value"] >= 0.8
        assert top1.loc[(t, "anchor"), "value"] > top1.loc[(t, "anchor"), "chance"] * 5
    # c shares no persons, so no replicate rows involve it
    replicate_labels = m.loc[m["metric"].str.startswith("replicate"), "label"].astype(str)
    assert set(replicate_labels) == {"a vs b"}

    offsets = pd.read_csv(out / "offsets.tsv", sep="\t").set_index(["celltype", "dataset"])
    assert offsets.loc[("Exc", "b"), "anchored"] and not offsets.loc[("Exc", "c"), "anchored"]
    assert offsets.loc[("Exc", "b"), "r_center_vs_anchor"] > 0.9

    de = pd.read_csv(out / "de_pooled.tsv", sep="\t")
    exc = de[de["celltype"] == "Exc"].set_index("gene_id")
    hits = exc.index[exc["fdr"] < 0.05]
    disease_ids = {f"ENSG{i:011d}" for i in DISEASE_GENES}
    assert disease_ids <= set(hits)
    assert (exc.loc[sorted(disease_ids), "log2fc"] > 0.5).all()
    disease = pd.read_csv(out / "disease.tsv", sep="\t")
    pooled = disease[(disease["celltype"] == "Exc") & (disease["analysis"] == "pooled: batch as covariate")]
    assert pooled["n_persons"].item() == 45 + 15   # one sample per person
    assert set(disease["analysis"]) >= {"a alone", "pooled: anchor, no batch term"}

    # ruv: the bridged offset and, through the pseudo-replicates, c's are removed
    for t in ("Exc", "Ast"):
        assert pc[(t, "ruv")] < 0.05
        assert top1.loc[(t, "ruv"), "value"] >= 0.8
    assert "pooled: ruv, no batch term" in set(disease["analysis"])
    assert offsets.loc[("Exc", "b"), "r_ruv_vs_anchor"] > 0.9
    scan = pd.read_csv(out / "ruv_scan.tsv", sep="\t")
    assert set(scan["k"]) == {1, 2, 4, 8, 12} and set(scan["controls"]) == {"empirical", "all"}
    assert set(scan["metric"]) >= {"replicate_top1", "replicate_rank_share", "pc_r2", "disease_t_correlation"}
    chosen = scan[scan["chosen"]].drop_duplicates(["celltype", "controls", "k"])
    assert sorted(chosen["celltype"]) == ["Ast", "Exc"]          # one choice per cell type
    k_rows = m[m["metric"] == "ruv_k"].set_index("celltype")
    assert (k_rows.loc[chosen["celltype"], "value"].to_numpy() <= chosen["k"].to_numpy()).all()
    diagnostic = pd.read_csv(out / "ruv_diagnostic.tsv", sep="\t")
    assert np.allclose(diagnostic[["batch_offset", "remainder_top_k", "rest"]].sum(axis=1), 1.0)
    assert (diagnostic["batch_offset"] > 0.3).all()            # b's x4 genes dominate the within-person gaps
    assert (diagnostic["n_pairs"] == 15).all()
    factors = pd.read_csv(out / "ruv_factors.tsv", sep="\t")
    assert set(factors["source"]) == {"ruv W", "bridge differences: remainder component"}
    assert factors["r_spearman"].abs().max() <= 1.0

    # aggregate tables hold no person or donor IDs
    for name in ("samples.tsv", "metrics.tsv", "disease.tsv", "offsets.tsv", "de_pooled.tsv", "genes.tsv",
                 "ruv_scan.tsv", "ruv_diagnostic.tsv", "ruv_factors.tsv"):
        text = (out / name).read_text()
        assert "P15" not in text and "a3" not in text and "S1" not in text

    with pytest.raises(FileExistsError):
        run_pooled_pseudobulk(tmp_path / "spec.yaml")
    run_pooled_pseudobulk(tmp_path / "spec.yaml", overwrite=True)   # reuses pseudobulk.h5ad


def test_batch_split_by_brain_bank(tmp_path):
    """b's donors come from two brain banks, one with extra offsets; with b
    split by brain bank, regress removes them and regress_dataset does not."""
    path = _setup(tmp_path, split=True)
    spec = yaml.safe_load(path.read_text())
    path.write_text(yaml.safe_dump({**spec, "ruv_method_controls": "all"}))
    out = run_pooled_pseudobulk(path)
    scan = pd.read_csv(out / "ruv_scan.tsv", sep="\t")
    assert set(scan.loc[scan["chosen"], "controls"]) == {"all"}      # ruv_method_controls
    assert set(scan["controls"]) == {"all", "empirical"}             # both still in the scan
    with pytest.raises(ValueError):
        path.write_text(yaml.safe_dump({**spec, "ruv_method_controls": "housekeeping"}))
        run_pooled_pseudobulk(path, overwrite=True)
    samples = pd.read_csv(out / "samples.tsv", sep="	")
    assert set(samples["batch"]) == {"a", "b:Rush", "b:Sinai", "c"}
    m = pd.read_csv(out / "metrics.tsv", sep="	")
    within = m[(m["metric"] == "pc_r2") & (m["label"] == "brain_bank within b")].set_index(["celltype", "method"])
    for t in ("Exc", "Ast"):
        assert within.loc[(t, "none"), "value"] > 0.1
        assert within.loc[(t, "regress_dataset"), "value"] > 0.1
        assert within.loc[(t, "regress"), "value"] < 0.05
        assert within.loc[(t, "center"), "value"] < 0.05
        # Sinai has no bridges; its pseudo-replicates link it to Rush
        assert within.loc[(t, "ruv"), "value"] < 0.05
    offsets = pd.read_csv(out / "offsets.tsv", sep="	").set_index(["celltype", "batch"])
    assert offsets.loc[("Exc", "b:Rush"), "anchored"] and not offsets.loc[("Exc", "b:Sinai"), "anchored"]
    disease = pd.read_csv(out / "disease.tsv", sep="	")
    assert "pooled: regress_dataset, no batch term" in set(disease["analysis"])
    assert (out / "umap" / "Exc.png").stat().st_size > 0 and (out / "umap" / "overview.png").stat().st_size > 0


def test_anchor_recovers_technical_offsets_despite_composition():
    """b has more diseased persons, which shifts its mean; the anchor offsets
    (from the same persons in both datasets) are the technical part only."""
    rng = np.random.default_rng(0)
    n_genes, shift = 50, np.zeros(50)
    shift[:10] = 1.0                          # technical offset of b
    disease_effect = np.zeros(n_genes)
    disease_effect[10:20] = 1.5
    persons_a = [f"P{i}" for i in range(40)]
    persons_b = [f"P{i}" for i in range(20, 60)]
    sick = {f"P{i}" for i in range(60) if i >= 30}   # b: 30 of 40 sick; a: 10 of 40
    base = {p: rng.normal(0, 0.3, n_genes) + (disease_effect if p in sick else 0)
            for p in sorted(set(persons_a + persons_b))}   # sorted: set order varies between runs
    Y = np.vstack([base[p] + rng.normal(0, 0.05, n_genes) for p in persons_a] +
                  [base[p] + shift + rng.normal(0, 0.05, n_genes) for p in persons_b])
    ds = np.array(["a"] * 40 + ["b"] * 40)
    persons = np.array(persons_a + persons_b)
    _, offsets, anchored = correct_anchor(Y, ds, persons)
    assert anchored == {"a", "b"}
    assert np.allclose(offsets["b"] - offsets["a"], shift, atol=0.1)
    _, center = correct_center(Y, ds)
    # centering also removes the composition difference in the disease genes
    assert (center["b"] - center["a"])[10:20].mean() > 0.5
    covs = pd.DataFrame({"dementia": ["yes" if p in sick else "no" for p in persons]})
    _, regress = correct_regress(Y, ds, covs, [])
    assert np.allclose(regress["b"] - regress["a"], shift, atol=0.15)


def test_ruv_matches_ruviii_formula():
    """Without pseudo-replicates, ruv_factors + ruv_apply are ruv::RUVIII
    (up to a constant per gene: Y is centered first)."""
    rng = np.random.default_rng(1)
    n, g, k = 30, 40, 3
    Y = rng.normal(size=(n, g)) + rng.normal(size=g)
    persons = np.array([f"P{i // 2}" if i < 20 else f"S{i}" for i in range(n)])   # 10 pairs, 10 singletons
    batch = np.array(["a" if i % 2 == 0 else "b" for i in range(n)])
    ctl = np.arange(0, g, 2)
    levels = list(dict.fromkeys(persons))
    M = (persons[:, None] == np.array(levels)[None, :]).astype(np.float64)
    Y0 = Y - M @ np.linalg.solve(M.T @ M, M.T @ Y)                       # residop(Y, M)
    values, U = np.linalg.eigh(Y0 @ Y0.T)
    alpha = U[:, np.argsort(values)[::-1][:k]].T @ Y
    ac = alpha[:, ctl]
    W = Y[:, ctl] @ ac.T @ np.linalg.inv(ac @ ac.T)
    expected = Y - W @ alpha

    a, info = ruv_factors(Y, batch, persons, n_factors=k)
    assert info["replicate_sets"] == 10 and info["pseudo_replicates"] == 0
    out, W2 = ruv_apply(Y, a, k, ctl)
    assert W2.shape == (n, k)
    assert np.allclose(out - out.mean(axis=0), expected - expected.mean(axis=0))
    assert np.allclose(out.mean(axis=0), Y.mean(axis=0))                 # gene means kept


def _bridged(rng, n_genes=200, quality=0.5):
    """a: P0-P59, b: P30-P89; b has an offset on genes 0-19, every sample a
    technical factor q (e.g. RNA quality) with its own loading, persons with
    odd numbers a disease effect on genes 100-119."""
    shift = np.zeros(n_genes)
    shift[:20] = 1.0
    loading = rng.normal(0, quality, n_genes)
    disease = np.zeros(n_genes)
    disease[100:120] = 1.0
    persons = np.array([f"P{i}" for i in range(60)] + [f"P{i}" for i in range(30, 90)])
    ds = np.array(["a"] * 60 + ["b"] * 60)
    sick = np.array([int(p[1:]) % 2 == 1 for p in persons])
    base = {p: rng.normal(0, 0.3, n_genes) for p in sorted(set(persons))}
    q = rng.normal(0, 1, len(persons))
    Y = np.vstack([base[p] + disease * s + shift * (d == "b") + qi * loading + rng.normal(0, 0.05, n_genes)
                   for p, d, s, qi in zip(persons, ds, sick, q)])
    return Y, ds, persons, sick, q


def test_ruv_removes_a_per_sample_factor_that_anchor_leaves():
    rng = np.random.default_rng(0)
    Y, ds, persons, sick, q = _bridged(rng)
    out, W, info = correct_ruv(Y, ds, persons, k=2)
    assert info["replicate_sets"] == 30
    assert max(abs(np.corrcoef(W[:, f], q)[0, 1]) for f in range(2)) > 0.9
    anchored, _, _ = correct_anchor(Y, ds, persons)
    rows_a, rows_b = np.arange(30, 60), np.arange(60, 90)                  # P30-P59 in a and in b

    def gap(X):
        return float(((X[rows_a] - X[rows_b]) ** 2).mean())

    assert gap(out) < gap(anchored) / 5
    effect = out[sick].mean(axis=0) - out[~sick].mean(axis=0)
    assert abs(effect[100:120].mean() - 1.0) < 0.15                      # the disease effect is kept
    assert abs(np.delete(effect, np.arange(100, 120)).mean()) < 0.1


def test_ruv_held_out_persons_do_not_shape_the_factors():
    """Cross-fitting: a fold's persons are singletons in M and left out of the
    pseudo-replicates, so their values do not change the factors."""
    rng = np.random.default_rng(2)
    Y, ds, persons, sick, _ = _bridged(rng)
    groups = np.where(sick, "|yes", "|no").astype(object)
    fold = frozenset(f"P{i}" for i in range(30, 40))
    a1, info = ruv_factors(Y, ds, persons, groups, exclude=fold, min_prps=3, n_factors=4)
    assert info["replicate_sets"] == 20 and info["pseudo_replicates"] == 4
    Y2 = Y.copy()
    held = np.isin(persons, list(fold))
    Y2[held] += rng.normal(0, 1, size=(held.sum(), Y.shape[1]))
    a2, _ = ruv_factors(Y2, ds, persons, groups, exclude=fold, min_prps=3, n_factors=4)

    def projector(a):
        return a.T @ np.linalg.solve(a @ a.T, a)

    assert np.allclose(projector(a1), projector(a2))
    a3, _ = ruv_factors(Y2, ds, persons, groups, min_prps=3, n_factors=4)
    assert not np.allclose(projector(a1), projector(a3))


def test_ruv_pseudo_replicates_remove_the_offset_of_a_batch_without_bridges():
    rng = np.random.default_rng(4)
    n_genes = 150
    Y, ds, persons, sick, _ = _bridged(rng, n_genes=n_genes, quality=0.0)
    offset_c = np.zeros(n_genes)
    offset_c[40:60] = 1.5                                                  # c: its own persons S0-S39
    persons_c = np.array([f"S{i}" for i in range(40)])
    sick_c = np.arange(40) % 2 == 1
    Yc = rng.normal(0, 0.3, size=(40, n_genes)) + offset_c + rng.normal(0, 0.05, size=(40, n_genes))
    Yc[sick_c, 100:120] += 1.0
    Y = np.vstack([Y, Yc])
    ds = np.concatenate([ds, ["c"] * 40])
    persons = np.concatenate([persons, persons_c])
    sub = pd.DataFrame({"dementia": np.where(np.concatenate([sick, sick_c]), "yes", "no"),
                        "sex": ["female" if i % 3 else "male" for i in range(len(ds))],
                        "age_capped": 75 + (np.arange(len(ds)) % 4) * 4})
    groups = prps_groups(sub, ["dementia", "sex", "age_capped"], [70, 80, 90])
    assert len(set(groups)) == 2 * 2 * 2

    def c_gap(X):
        return float((X[ds == "c", 40:60].mean() - X[ds != "c", 40:60].mean()))

    with_prps, _, info = correct_ruv(Y, ds, persons, k=2, groups=groups, min_prps=3)
    assert info["pseudo_replicate_sets"] == 8
    without, _, _ = correct_ruv(Y, ds, persons, k=2)
    assert abs(c_gap(without)) > 1.0
    assert abs(c_gap(with_prps)) < 0.15
    assert abs((with_prps[ds == "b", :20].mean() - with_prps[ds == "a", :20].mean())) < 0.1


def test_spec_pseudobulk_reads_an_existing_one(tmp_path):
    """pseudobulk: another output_dir reads the pseudobulk.h5ad and never writes it."""
    path = _setup(tmp_path)
    first = run_pooled_pseudobulk(path)
    h5ad = first / "pseudobulk.h5ad"
    stamp = h5ad.stat().st_mtime_ns
    spec = yaml.safe_load(path.read_text())
    spec.update(output_dir=str(tmp_path / "again"), pseudobulk=str(h5ad))
    other = tmp_path / "again.yaml"
    other.write_text(yaml.safe_dump(spec))
    out = run_pooled_pseudobulk(other)
    assert not (out / "pseudobulk.h5ad").exists() and h5ad.stat().st_mtime_ns == stamp
    a = pd.read_csv(first / "metrics.tsv", sep="\t")
    b = pd.read_csv(out / "metrics.tsv", sep="\t")
    pd.testing.assert_frame_equal(a, b)
    with pytest.raises(ValueError):
        run_pooled_pseudobulk(other, rebuild=True)
