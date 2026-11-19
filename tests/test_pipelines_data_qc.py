"""End-to-end tests of `pipeline qc-scan` and `qc-check` on synthetic count
files: two datasets of the same donors, one donor with the wrong recorded sex,
the authors' per-cell totals, and strong donor effects."""

import shutil

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from rosmap_processing.pipelines import data_qc

N_DONORS = 12
SWAPPED = 3            # recorded sex is wrong for this donor
NAMED = ["XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "MT-CO1", "MT-ND1", "SLC17A7", "GAD1", "AQP4", "SNAP25",
         "PTPRC"]
GENES = NAMED + [f"G{i}" for i in range(150)]


def _counts(path, seed, donor_effects, authors=False):
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 1.0, size=len(GENES))
    blocks, obs = [], []
    for d in range(N_DONORS):
        male = d % 2 == 0
        for celltype, n in (("Excitatory", 40), ("Astrocytes", 40), ("Immune", 6)):
            r = np.tile(base * donor_effects[d], (n, 1))
            r[:, 0] = 0.01 if male else 5.0               # XIST
            r[:, 1:5] = 3.0 if male else 0.0              # Y genes
            r[:, 5:7] = 2.0                               # mitochondrial
            r[:, 7] = 8.0 if celltype == "Excitatory" else 0.05
            r[:, 10] = 6.0 if celltype == "Excitatory" else 0.2
            r[:, 9] = 8.0 if celltype == "Astrocytes" else 0.05
            r[:, 11] = 6.0 if celltype == "Immune" else 0.0
            depth = rng.uniform(0.5, 2.0, size=(n, 1))
            blocks.append(rng.poisson(r * depth))
            recorded_male = male != (d == SWAPPED)
            obs.append(pd.DataFrame({
                "projid": 100 + d, "individualID": f"R{d:03d}",
                "source_file": f"x/{celltype}.h5ad", "fine": f"{celltype}_sub",
                "msex": float(recorded_male), "age_death": str(70 + d), "cogdx": 4.0 if d % 3 == 0 else 1.0,
                "braaksc": 5.0 if d % 3 == 0 else 1.0, "ceradsc": 1.0 if d % 3 == 0 else 4.0,
                "braak_other": 5.0 if d % 3 == 0 else (1.0 if d != 5 else 2.0),
            }, index=[f"d{d}_{celltype}_{i}" for i in range(n)]))
    counts = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    obs = pd.concat(obs)
    if authors:
        obs["nCount_RNA"] = np.asarray(counts.sum(axis=1)).ravel()
        obs["nFeature_RNA"] = np.diff(counts.indptr)
        obs["mt_percent"] = 100 * np.asarray(counts[:, 5:7].sum(axis=1)).ravel() / obs["nCount_RNA"]
    ad.AnnData(counts, obs=obs, var=pd.DataFrame(index=GENES)).write_h5ad(path)


@pytest.fixture
def setup(tmp_path):
    rng = np.random.default_rng(0)
    donor_effects = rng.lognormal(0, 0.6, size=(N_DONORS, len(GENES)))
    donor_effects[:, :12] = 1.0
    _counts(tmp_path / "a.h5ad", 1, donor_effects, authors=True)
    _counts(tmp_path / "b.h5ad", 2, donor_effects)
    pd.DataFrame({"var_name": GENES, "gene_id": [f"ENSG{i:011d}" for i in range(len(GENES))],
                  "ref_name": GENES, "method": "symbol"}).to_csv(tmp_path / "map.tsv", sep="\t", index=False)
    (tmp_path / "mt.txt").write_text("ENSG00000000005\nENSG00000000006\n")
    (tmp_path / "x.txt").write_text("ENSG00000000000\n")
    (tmp_path / "y.txt").write_text("\n".join(f"ENSG{i:011d}" for i in range(1, 5)) + "\n")
    (tmp_path / "genes.txt").write_text("\n".join(GENES) + "\n")
    specs = []
    for name in ("a", "b"):
        be_spec = {
            "output_dir": str(tmp_path / f"be_{name}"), "gene_sets": {"all": str(tmp_path / "genes.txt")},
            "datasets": {name: {
                "counts": str(tmp_path / f"{name}.h5ad"), "gene_map": str(tmp_path / "map.tsv"),
                "donor": "projid" if name == "a" else "individualID",
                "celltype": [{"column": "source_file", "match": "Excitatory", "label": "Exc"},
                             {"column": "source_file", "match": "Astrocytes", "label": "Ast"}],
                "covariates": {"sex": {"column": "msex"}, "age": {"column": "age_death", "numeric": True},
                               "braak": {"column": "braaksc", "numeric": True},
                               "cerad": {"column": "ceradsc", "numeric": True,
                                         "recode": {"1": 3, "2": 2, "3": 1, "4": 0}}},
                "disease": {"column": "cogdx", "recode": {"1": "no", "4": "yes"}},
            }},
        }
        (tmp_path / f"be_{name}.yaml").write_text(yaml.safe_dump(be_spec))
        qc = {"output_dir": str(tmp_path / f"qc_{name}"), "batch_effects": f"be_{name}.yaml",
              "match_on": "individualID", "author_label": "fine", "min_cells": 10, "min_cells_identity": 10,
              "min_donors_per_group": 3, "identity_genes": 100, "chunk_size": 200,
              "covariates": {"cogdx": {"column": "cogdx"}},
              "mt_genes": str(tmp_path / "mt.txt"), "x_genes": str(tmp_path / "x.txt"),
              "y_genes": str(tmp_path / "y.txt")}
        if name == "a":
            qc["author_qc"] = {"umis": "nCount_RNA", "genes": "nFeature_RNA", "mito_pct": "mt_percent"}
            qc["metadata_sources"] = {"braak": [{"column": "braaksc", "numeric": True},
                                                {"column": "braak_other", "numeric": True}]}
        (tmp_path / f"qc_{name}.yaml").write_text(yaml.safe_dump(qc))
        specs.append(str(tmp_path / f"qc_{name}.yaml"))
    for s in specs:
        data_qc.run_scan(s)
    return tmp_path, specs


def test_totals_match_the_authors(setup):
    tmp_path, specs = setup
    data_qc.run_check("totals", specs[:1])
    t = pd.read_csv(tmp_path / "qc_a" / "totals.tsv", sep="\t")
    overall = t[t["celltype"] == "all"].set_index("field")
    assert set(overall.index) == {"umis", "genes", "mito_pct"}
    assert (overall["pct_exact"] == 100).all()
    assert overall.loc["umis", "n_compared"] == N_DONORS * 86


def test_sex_check_finds_the_swapped_donor(setup):
    tmp_path, specs = setup
    data_qc.run_check("sex", specs[:1])
    summary = pd.read_csv(tmp_path / "qc_a" / "sex_summary.tsv", sep="\t").iloc[0]
    assert summary["n_discordant"] == 1 and summary["n_ambiguous"] == 0
    follow_up = pd.read_csv(tmp_path / "qc_a" / "restricted" / "sex_follow_up.tsv", sep="\t")
    assert list(follow_up["donor"].astype(str)) == [str(100 + SWAPPED)]
    cell = pd.read_csv(tmp_path / "qc_a" / "sex_cell_level.tsv", sep="\t")
    female = cell[cell["measure"] == "share of cells with Y-gene UMIs"].iloc[0]
    assert female["max"] == 0


def test_cells_markers_and_de(setup):
    tmp_path, specs = setup
    for check in ("cells", "markers", "sex-de", "ad-de"):
        data_qc.run_check(check, specs[:1])
    cells = pd.read_csv(tmp_path / "qc_a" / "qc_cells.tsv", sep="\t").set_index("group")
    assert cells.loc["all", "n_cells"] == N_DONORS * 86
    assert cells.loc["Other: Immune_sub", "n_cells"] == N_DONORS * 6
    spec_table = pd.read_csv(tmp_path / "qc_a" / "marker_specificity.tsv", sep="\t")
    assert spec_table.set_index("gene").loc["SLC17A7", "own_is_highest"]
    assert spec_table.set_index("gene").loc["AQP4", "own_is_highest"]
    other = pd.read_csv(tmp_path / "qc_a" / "unmatched_labels.tsv", sep="\t").iloc[0]
    assert other["top_marker_set"] == "Immune"
    de = pd.read_csv(tmp_path / "qc_a" / "sex_de.tsv", sep="\t")
    exc = de[de["celltype"] == "Exc"].set_index("symbol")
    assert exc.loc["XIST", "t"] < 0 and exc.loc["RPS4Y1", "t"] > 0
    summary = pd.read_csv(tmp_path / "qc_a" / "ad_de_summary.tsv", sep="\t")
    assert set(summary["contrast"]) == {"pathology", "disease"}


def test_metadata_and_identity(setup):
    tmp_path, specs = setup
    data_qc.run_check("metadata", specs[:1])
    agreement = pd.read_csv(tmp_path / "qc_a" / "metadata_agreement.tsv", sep="\t").iloc[0]
    assert agreement["n_both"] == N_DONORS and agreement["n_agree"] == N_DONORS - 1
    data_qc.run_check("sex-de", specs)
    data_qc.run_check("sex", specs)
    out = data_qc.run_check("identity", specs, output=tmp_path / "identity")
    meta = pd.read_csv(out / "identity_metadata.tsv", sep="\t").set_index("covariate")
    assert meta.loc["sex_recorded", "n_agree"] == N_DONORS
    assert meta.loc["sex_expression", "n_agree"] == N_DONORS
    # The swapped donor's record is wrong in both datasets, so neither backs it
    assert meta.loc["sex_expression_off_in_a", "n_both"] == 1
    assert meta.loc["sex_expression_off_in_a", "n_agree"] == 0
    expr = pd.read_csv(out / "identity_expression.tsv", sep="\t").set_index("celltype")
    assert expr.loc["combined", "pct_own_best_a_to_b"] == 100
    data_qc.run_check("sex-de-compare", specs, output=tmp_path / "compare")
    compare = pd.read_csv(tmp_path / "compare" / "sex_de_compare.tsv", sep="\t")
    assert compare[compare["genes"] == "all"]["pearson_t"].min() > 0


GENO_DONORS, GENO_PERSONS, GENO_EQTLS = 60, 100, 150


def _genotype_setup(tmp_path):
    """60 donors, 100 genotyped persons, 150 strong eQTLs. Donor 10's sample is
    person 80's; donors 20 and 21 are swapped; donor 30 has no genotyped
    label. Writes the counts, the scan and eqtls.tsv; returns the dosages
    (persons x eQTLs)."""
    rng = np.random.default_rng(7)
    n_donors, n_persons, k = GENO_DONORS, GENO_PERSONS, GENO_EQTLS
    dosage = rng.binomial(2, rng.uniform(0.25, 0.5, size=k), size=(n_persons, k)).astype(float)
    true = np.arange(n_donors)
    true[10], true[20], true[21] = 80, 21, 20
    genes = NAMED + [f"E{i}" for i in range(k)] + [f"G{i}" for i in range(100)]
    base = rng.gamma(2.0, 1.0, size=len(genes))
    base[len(NAMED):len(NAMED) + k] = 2.0
    blocks, obs = [], []
    for d in range(n_donors):
        for celltype in ("Excitatory", "Astrocytes"):
            r = np.tile(base * rng.lognormal(0, 0.1, size=len(genes)), (25, 1))
            r[:, len(NAMED):len(NAMED) + k] *= np.exp(1.2 * (dosage[true[d]] - 1))
            r[:, 7] = 8.0 if celltype == "Excitatory" else 0.05
            r[:, 9] = 8.0 if celltype == "Astrocytes" else 0.05
            blocks.append(rng.poisson(r * rng.uniform(0.5, 2.0, size=(25, 1))))
            obs.append(pd.DataFrame({"individualID": f"R{d:03d}" if d != 30 else "X030",
                                     "source_file": f"x/{celltype}.h5ad", "fine": f"{celltype}_sub",
                                     "msex": float(d % 2)}, index=[f"d{d}_{celltype}_{i}" for i in range(25)]))
    ad.AnnData(sp.csr_matrix(np.vstack(blocks).astype(np.float32)), obs=pd.concat(obs),
               var=pd.DataFrame(index=genes)).write_h5ad(tmp_path / "g.h5ad")
    gene_ids = [f"ENSG{i:011d}" for i in range(len(genes))]
    pd.DataFrame({"var_name": genes, "gene_id": gene_ids, "ref_name": genes, "method": "symbol"}).to_csv(
        tmp_path / "map.tsv", sep="\t", index=False)
    (tmp_path / "genes.txt").write_text("\n".join(genes) + "\n")
    for name in ("mt", "x", "y"):
        (tmp_path / f"{name}.txt").write_text("\n")
    be_spec = {"output_dir": str(tmp_path / "be"), "gene_sets": {"all": str(tmp_path / "genes.txt")},
               "datasets": {"g": {"counts": str(tmp_path / "g.h5ad"), "gene_map": str(tmp_path / "map.tsv"),
                                  "donor": "individualID",
                                  "celltype": [{"column": "source_file", "match": "Excitatory", "label": "Exc"},
                                               {"column": "source_file", "match": "Astrocytes", "label": "Ast"}],
                                  "covariates": {"sex": {"column": "msex"}}}}}
    (tmp_path / "be.yaml").write_text(yaml.safe_dump(be_spec))
    qc = {"output_dir": str(tmp_path / "qc"), "batch_effects": "be.yaml", "match_on": "individualID",
          "author_label": "fine", "min_cells": 10, "chunk_size": 500,
          "mt_genes": str(tmp_path / "mt.txt"), "x_genes": str(tmp_path / "x.txt"), "y_genes": str(tmp_path / "y.txt")}
    (tmp_path / "qc.yaml").write_text(yaml.safe_dump(qc))
    data_qc.run_scan(tmp_path / "qc.yaml")
    # A quarter of the eQTLs give the alleles on the other strand
    pd.DataFrame({"gene_id": gene_ids[len(NAMED):len(NAMED) + k], "symbol": [f"E{i}" for i in range(k)],
                  "tissue": "t", "rank": 1, "rsid": [f"rs{i}" for i in range(k)],
                  "variant_b38": [f"chr1_{5000 + i}_A_G_b38" for i in range(k)], "chrom": "1",
                  "pos_b37": [1000 + i for i in range(k)],
                  "ref": ["T" if i % 4 == 0 else "A" for i in range(k)],
                  "alt": ["C" if i % 4 == 0 else "G" for i in range(k)],
                  "maf": 0.3, "pvalue": 1e-20, "nes": 1.0}).to_csv(tmp_path / "eqtls.tsv", sep="\t", index=False)
    return dosage


def _write_dosages(path, dosage, snp, samples):
    geno = pd.DataFrame(dosage.T, columns=samples)
    geno.insert(0, "snp", snp)
    geno.insert(1, "a1", "A")
    geno.insert(2, "a2", "G")
    geno.to_csv(path, sep="\t", index=False)


def _assert_swaps_found(out, id_col):
    n_donors = GENO_DONORS
    s = pd.read_csv(out / "genotype_summary.tsv", sep="\t").iloc[0]
    assert s["n_eqtl_pairs"] == GENO_EQTLS and s["n_genotyped_label"] == n_donors - 1
    assert s["n_other_person"] == 3 and s["n_other_person_in_dataset"] == 2 and s["n_mutual_swaps"] == 1
    assert s["n_own_genotype"] == n_donors - 4 and s["n_identified"] == 1
    calls = pd.read_csv(out / "restricted" / "genotype_calls_qc.tsv", sep="\t", dtype=str).set_index("donor")
    assert calls.loc["R020", "best_is_donor"] == "R021"
    return calls[["call", f"best_{id_col}"]]


def _rosmap_genotypes(tmp_path, dosage, duplicate=None):
    """ROSMAP-style dosages (samples ROS/MAP + projid; half the variants by
    rsID, half by GRCh37 chr:pos) and an ID map with msex (projid p: p % 2).
    duplicate: a person genotyped a second time, as projid 100 + that person."""
    n_persons, k = GENO_PERSONS, GENO_EQTLS
    samples = [f"{'ROS' if p % 2 else 'MAP'}{p:08d}" for p in range(n_persons)]
    if duplicate is not None:
        dosage = np.vstack([dosage, dosage[duplicate]])
        samples.append(f"ROS{100 + duplicate:08d}")
    _write_dosages(tmp_path / "geno.tsv.gz", dosage, [f"rs{i}" if i % 2 else f"chr1:{1000 + i}" for i in range(k)],
                   samples)
    pd.DataFrame({"individualID": [f"R{p:03d}" for p in range(n_persons)],
                  "projid": [f"{p:08d}" for p in range(n_persons)],
                  "msex": [str(p % 2) for p in range(n_persons)]}).to_csv(tmp_path / "ids.csv", index=False)


def test_genotype_check_finds_swaps(tmp_path, monkeypatch):
    """ROSMAP-style: samples named ROS/MAP + projid, labels mapped to projids."""
    monkeypatch.setattr(data_qc, "GENOTYPE_MIN_VARIANTS", 50)
    _rosmap_genotypes(tmp_path, _genotype_setup(tmp_path))
    out = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "genotype",
                            genotypes=tmp_path / "geno.tsv.gz", eqtls=tmp_path / "eqtls.tsv",
                            id_map=tmp_path / "ids.csv")
    calls = _assert_swaps_found(out, "projid")
    assert calls.loc["R010"].tolist() == ["other person", "80"]
    assert calls.loc["X030"].tolist() == ["identified", "30"]


def test_genotype_check_needs_enough_variants(tmp_path):
    """Fewer kept variants than GENOTYPE_MIN_VARIANTS (150 < 200): no calls."""
    _rosmap_genotypes(tmp_path, _genotype_setup(tmp_path))
    out = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "genotype",
                            genotypes=tmp_path / "geno.tsv.gz", eqtls=tmp_path / "eqtls.tsv",
                            id_map=tmp_path / "ids.csv")
    s = pd.read_csv(out / "genotype_summary.tsv", sep="\t").iloc[0]
    assert not s["calls_made"] and 0 < s["n_kept_variants"] <= GENO_EQTLS
    assert s["n_no_data"] == GENO_DONORS


def _write_sex_by_donor(qc_dir):
    """The sex check's per-donor file: expression sex is the true person's
    (persons: p % 2 = male), recorded is the label's."""
    true = np.arange(GENO_DONORS)
    true[10], true[20], true[21] = 80, 21, 20
    donor = [f"R{d:03d}" if d != 30 else "X030" for d in range(GENO_DONORS)]
    pd.DataFrame({"donor": donor, "recorded": ["M" if d % 2 else "F" for d in range(GENO_DONORS)],
                  "expression": ["M" if t % 2 else "F" for t in true]}).to_csv(
        qc_dir / "restricted" / "sex_by_donor.tsv", sep="\t", index=False)


def _genotype_actions(monkeypatch, *argv):
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "scripts" / "genotype_actions.py"
    spec = importlib.util.spec_from_file_location("genotype_actions", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr("sys.argv", ["genotype_actions.py", *map(str, argv)])
    module.main()


def test_genotype_actions(tmp_path, monkeypatch):
    """scripts/genotype_actions.py: the mutual swap is confirmed by sex and
    relabelled; donor 10 (same sex as its person 80) and X030 (no genotyped
    label) are unconfirmed and dropped; the rest are kept. With a second
    dataset holding the same samples, the swaps are contradicted."""
    monkeypatch.setattr(data_qc, "GENOTYPE_MIN_VARIANTS", 50)
    _rosmap_genotypes(tmp_path, _genotype_setup(tmp_path))
    _write_sex_by_donor(tmp_path / "qc")
    kwargs = dict(genotypes=tmp_path / "geno.tsv.gz", eqtls=tmp_path / "eqtls.tsv", id_map=tmp_path / "ids.csv")
    one = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "one", **kwargs)
    _genotype_actions(monkeypatch, one, "--id-map", tmp_path / "ids.csv")
    actions = pd.read_csv(one / "restricted" / "genotype_actions.tsv", sep="\t", dtype=str).set_index("donor")
    assert actions.loc["R020", ["action", "relabel_to", "relabel_to_individualID"]].tolist() == \
        ["relabel", "21", "R021"]
    assert actions.loc["R021", ["action", "relabel_to"]].tolist() == ["relabel", "20"]
    assert "mutual swap" in actions.loc["R020", "motivation"]
    assert actions.loc["R010", ["action", "reason"]].tolist() == ["drop", "unconfirmed"]
    assert actions.loc["X030", "action"] == "drop"
    assert (actions.drop(index=["R010", "R020", "R021", "X030"])["action"] == "keep").all()
    summary = pd.read_csv(one / "genotype_actions_summary.tsv", sep="\t")
    assert summary["n_donors"].sum() == GENO_DONORS

    shutil.copytree(tmp_path / "qc", tmp_path / "qc2")
    qc2 = yaml.safe_load((tmp_path / "qc.yaml").read_text())
    qc2["output_dir"] = str(tmp_path / "qc2")
    (tmp_path / "qc2.yaml").write_text(yaml.safe_dump(qc2))
    two = data_qc.run_check("genotype", [tmp_path / "qc.yaml", tmp_path / "qc2.yaml"], output=tmp_path / "two",
                            **kwargs)
    _genotype_actions(monkeypatch, two, "--id-map", tmp_path / "ids.csv", "--output", tmp_path / "actions")
    actions = pd.read_csv(tmp_path / "actions" / "restricted" / "genotype_actions.tsv", sep="\t", dtype=str)
    swapped = actions[actions["donor"].isin(["R020", "R021"])]
    assert len(swapped) == 4 and (swapped["reason"] == "contradicted").all()


def test_genotype_check_validation_tables(tmp_path, monkeypatch):
    """The tables for checking the calls: matched-person sex, a second dataset
    with the same samples (cross detail), a duplicated genotype sample
    (kinship) and expected z / weak flags."""
    monkeypatch.setattr(data_qc, "GENOTYPE_MIN_VARIANTS", 50)
    _rosmap_genotypes(tmp_path, _genotype_setup(tmp_path), duplicate=90)
    _write_sex_by_donor(tmp_path / "qc")
    # A second dataset with the same nuclei
    shutil.copytree(tmp_path / "qc", tmp_path / "qc2")
    qc2 = yaml.safe_load((tmp_path / "qc.yaml").read_text())
    qc2["output_dir"] = str(tmp_path / "qc2")
    (tmp_path / "qc2.yaml").write_text(yaml.safe_dump(qc2))

    out = data_qc.run_check("genotype", [tmp_path / "qc.yaml", tmp_path / "qc2.yaml"],
                            output=tmp_path / "genotype", genotypes=tmp_path / "geno.tsv.gz",
                            eqtls=tmp_path / "eqtls.tsv", id_map=tmp_path / "ids.csv")
    s = pd.read_csv(out / "genotype_summary.tsv", sep="\t").set_index("dataset")
    assert s.loc["qc", "calls_made"] and s.loc["qc", "n_other_person"] == 3
    sex = pd.read_csv(out / "genotype_matched_sex.tsv", sep="\t").set_index(["dataset", "call"])
    assert sex.loc[("qc", "other person"), ["n", "n_matched_sex_known", "n_matched_sex_agrees"]].tolist() == [3, 3, 3]
    assert sex.loc[("qc", "other person"), "n_label_sex_agrees"] == 1  # donor 10's person 80 has its sex
    assert sex.loc[("qc", "identified"), "n_matched_sex_agrees"] == 1
    own = sex.loc[("qc", "own genotype")]
    assert own["n_label_sex_agrees"] == own["n_label_sex_known"] == own["n"]
    kin = pd.read_csv(out / "genotype_kinship.tsv", sep="\t").iloc[0]
    assert kin["n_samples"] == GENO_PERSONS + 1 and kin["n_pairs_duplicate"] == 1
    pairs = pd.read_csv(out / "restricted" / "genotype_kinship_pairs.tsv", sep="\t", dtype=str)
    assert {"MAP00000090", "ROS00000190"} <= set(pairs["sample_a"]) | set(pairs["sample_b"])
    flagged = pd.read_csv(out / "genotype_kinship_flagged.tsv", sep="\t").set_index(["dataset", "call"])
    assert flagged.loc[("qc", "other person"), "n_best_has_duplicate"] == 0
    cross = pd.read_csv(out / "genotype_cross.tsv", sep="\t").iloc[0]
    assert cross["n_same_genotype"] == cross["n_persons_identified_in_both"] == GENO_DONORS
    assert cross["n_same_genotype_not_label"] == 3  # the same swapped samples in both
    detail = pd.read_csv(out / "genotype_cross_detail.tsv", sep="\t").set_index(["dataset_a", "call_a", "call_b"])
    row = detail.loc[("qc", "other person", "other person")]
    assert row["n"] == row["n_b_same_best"] == row["n_b_z_at_best_a_ge6"] == 3
    calls = pd.read_csv(out / "restricted" / "genotype_calls_qc.tsv", sep="\t").set_index("donor")
    mine = calls[calls["call"] == "own genotype"]
    assert mine["expected_z"].notna().all() and mine["strength"].between(0.4, 2.5).all()
    assert {"next_z", "weak", "best_sex", "best_max_kinship", "label_best_kinship"} <= set(calls.columns)


def test_genotype_check_with_sample_ids_and_grch38(tmp_path, monkeypatch):
    """PsychAD-style: the samples are named by the donors' match key and the
    variants by GRCh38 chr:pos (half by rsID); no ID map."""
    monkeypatch.setattr(data_qc, "GENOTYPE_MIN_VARIANTS", 50)
    dosage = _genotype_setup(tmp_path)
    n_persons, k = GENO_PERSONS, GENO_EQTLS
    _write_dosages(tmp_path / "geno.tsv.gz", dosage, [f"rs{i}" if i % 2 else f"chr1:{5000 + i}" for i in range(k)],
                   [f"R{p:03d}" for p in range(n_persons)])
    kwargs = dict(genotypes=tmp_path / "geno.tsv.gz", eqtls=tmp_path / "eqtls.tsv", sample_ids="match")
    out = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "genotype", genotype_build=38,
                            **kwargs)
    calls = _assert_swaps_found(out, "sample")
    assert calls.loc["R010"].tolist() == ["other person", "R080"]
    assert calls.loc["X030"].tolist() == ["identified", "R030"]
    # GRCh37 positions find only the rsID-named half
    out = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "genotype37", **kwargs)
    assert pd.read_csv(out / "genotype_summary.tsv", sep="\t").iloc[0]["n_eqtl_pairs"] == k // 2


def test_genotype_actions_with_wgs(tmp_path, monkeypatch):
    """Array and WGS combined. The array has persons 40 and 41's DNA swapped
    (WGS is right): their donors are kept. Donor 10's sample is person 80, of
    the same sex: unconfirmed by either source alone, confirmed by both
    matching 80. Donors 20 and 21 are a mutual swap."""
    monkeypatch.setattr(data_qc, "GENOTYPE_MIN_VARIANTS", 50)
    dosage = _genotype_setup(tmp_path)
    _write_sex_by_donor(tmp_path / "qc")
    array = dosage.copy()
    array[[40, 41]] = array[[41, 40]]
    _rosmap_genotypes(tmp_path, array)
    arr = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "array",
                            genotypes=tmp_path / "geno.tsv.gz", eqtls=tmp_path / "eqtls.tsv",
                            id_map=tmp_path / "ids.csv")
    _write_dosages(tmp_path / "wgs.tsv.gz", dosage, [f"chr1:{1000 + i}" for i in range(GENO_EQTLS)],
                   [f"R{p:03d}" for p in range(GENO_PERSONS)])
    wgs = data_qc.run_check("genotype", [tmp_path / "qc.yaml"], output=tmp_path / "wgs",
                            genotypes=tmp_path / "wgs.tsv.gz", eqtls=tmp_path / "eqtls.tsv", sample_ids="match")
    status = tmp_path / "avw" / "restricted"
    status.mkdir(parents=True)
    pd.DataFrame({"projid": [str(p) for p in range(GENO_PERSONS)],
                  "status": ["matches another person's WGS" if p in (40, 41) else "same person"
                             for p in range(GENO_PERSONS)]}).to_csv(status / "array_vs_wgs_by_person.tsv",
                                                                     sep="\t", index=False)
    _genotype_actions(monkeypatch, arr, "--wgs", f"wgs={wgs}", "--array-vs-wgs", tmp_path / "avw",
                      "--id-map", tmp_path / "ids.csv", "--output", tmp_path / "actions")
    a = pd.read_csv(tmp_path / "actions" / "restricted" / "genotype_actions.tsv", sep="\t", dtype=str,
                    keep_default_na=False).set_index("donor")
    for d in ("R040", "R041"):
        assert a.loc[d, ["action", "reason", "identity"]].tolist() == ["keep", "own genotype, array entry wrong",
                                                                       "verified"]
    assert a.loc["R010", ["action", "reason", "relabel_to_individualID"]].tolist() == \
        ["relabel", "confirmed by two genotype sources", "R080"]
    assert a.loc["R020", ["action", "relabel_to_individualID"]].tolist() == ["relabel", "R021"]
    assert a.loc["R021", ["action", "relabel_to_individualID"]].tolist() == ["relabel", "R020"]
    assert "mutual swap" in a.loc["R020", "motivation"]
    assert a.loc["X030", "action"] == "drop"
    rest = a.drop(index=["R010", "R020", "R021", "X030", "R040", "R041"])
    assert (rest["action"] == "keep").all() and (rest["identity"] == "verified").all()
    assert {"array_call", "wgs_call", "decided_by"} <= set(a.columns)
