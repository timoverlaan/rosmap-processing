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
            "umap": split, "batch_split": {"b": "brain_bank"} if split else {},
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

    # aggregate tables hold no person or donor IDs
    for name in ("samples.tsv", "metrics.tsv", "disease.tsv", "offsets.tsv", "de_pooled.tsv", "genes.tsv"):
        text = (out / name).read_text()
        assert "P15" not in text and "a3" not in text and "S1" not in text

    with pytest.raises(FileExistsError):
        run_pooled_pseudobulk(tmp_path / "spec.yaml")
    run_pooled_pseudobulk(tmp_path / "spec.yaml", overwrite=True)   # reuses pseudobulk.h5ad


def test_batch_split_by_brain_bank(tmp_path):
    """b's donors come from two brain banks, one with extra offsets; with b
    split by brain bank, regress removes them and regress_dataset does not."""
    out = run_pooled_pseudobulk(_setup(tmp_path, split=True))
    samples = pd.read_csv(out / "samples.tsv", sep="	")
    assert set(samples["batch"]) == {"a", "b:Rush", "b:Sinai", "c"}
    m = pd.read_csv(out / "metrics.tsv", sep="	")
    within = m[(m["metric"] == "pc_r2") & (m["label"] == "brain_bank within b")].set_index(["celltype", "method"])
    for t in ("Exc", "Ast"):
        assert within.loc[(t, "none"), "value"] > 0.1
        assert within.loc[(t, "regress_dataset"), "value"] > 0.1
        assert within.loc[(t, "regress"), "value"] < 0.05
        assert within.loc[(t, "center"), "value"] < 0.05
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
