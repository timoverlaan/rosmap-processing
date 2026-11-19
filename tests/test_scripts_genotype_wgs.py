"""scripts/genotype_wgs_check.py on synthetic array and WGS genotypes."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _script():
    spec = importlib.util.spec_from_file_location("genotype_wgs_check", ROOT / "scripts" / "genotype_wgs_check.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write(path, snp, a1, a2, dosage, samples):
    table = pd.DataFrame(dosage, columns=samples)
    table.insert(0, "snp", snp)
    table.insert(1, "a1", a1)
    table.insert(2, "a2", a2)
    table.to_csv(path, sep="\t", index=False)


def test_genotype_wgs_check(tmp_path, monkeypatch):
    """40 people on the array, 35 with WGS. WGS sample SM-5 holds the DNA of
    person 38 (who has no WGS of their own); half the WGS variants count the
    other allele."""
    rng = np.random.default_rng(3)
    n_people, n_var = 40, 300
    dna = rng.binomial(2, rng.uniform(0.2, 0.5, size=n_var)[:, None], size=(n_var, n_people)).astype(float)
    _write(tmp_path / "array.tsv.gz", [f"rs{i}" for i in range(n_var)], "G", "A", dna,
           [f"{'ROS' if p % 2 else 'MAP'}{p:08d}" for p in range(n_people)])
    wgs = dna[:, :35].copy()
    wgs[:, 5] = dna[:, 38]
    flip = np.arange(n_var) % 2 == 1
    wgs[flip] = 2 - wgs[flip]
    _write(tmp_path / "wgs.tsv.gz", [f"chr1:{1000 + i}" for i in range(n_var)],
           np.where(flip, "A", "G"), np.where(flip, "G", "A"), wgs, [f"SM-{p}" for p in range(35)])
    pd.DataFrame({"gene_id": [f"ENSG{i}" for i in range(n_var)], "rsid": [f"rs{i}" for i in range(n_var)],
                  "chrom": "1", "pos_b37": [1000 + i for i in range(n_var)],
                  "variant_b38": [f"chr1_{5000 + i}_A_G_b38" for i in range(n_var)], "ref": "A", "alt": "G"}
                 ).to_csv(tmp_path / "eqtls.tsv", sep="\t", index=False)
    pd.DataFrame({"individualID": [f"R{p:03d}" for p in range(n_people)],
                  "projid": [f"{p:08d}" for p in range(n_people)],
                  "msex": [p % 2 for p in range(n_people)]}).to_csv(tmp_path / "clinical.csv", index=False)
    pd.DataFrame({"specimenID": [f"SM-{p}" for p in range(35)], "individualID": [f"R{p:03d}" for p in range(35)],
                  "assay": "wholeGenomeSeq"}).to_csv(tmp_path / "biospecimen.csv", index=False)
    rate = lambda people: [0.01 if p % 2 else 0.3 for p in people]  # noqa: E731
    pd.DataFrame({"IID": [f"ROS{p:08d}" for p in range(n_people)], "het_rate": rate(range(n_people))}).to_csv(
        tmp_path / "array_chrx.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": [f"SM-{p}" for p in range(35)],
                  "het_rate": rate([38 if p == 5 else p for p in range(35)])}).to_csv(
        tmp_path / "wgs_chrx.tsv", sep="\t", index=False)
    pd.DataFrame({"dataset": ["rosmap_mit"], "donor": ["D1"], "action": ["drop"], "reason": ["unconfirmed"],
                  "label": ["5"], "best": ["38"], "expression_sex": ["F"]}).to_csv(
        tmp_path / "actions.tsv", sep="\t", index=False)

    monkeypatch.setattr("sys.argv", ["genotype_wgs_check.py", "--array", str(tmp_path / "array.tsv.gz"),
                                     "--array-chrx", str(tmp_path / "array_chrx.tsv"),
                                     "--wgs", str(tmp_path / "wgs.tsv.gz"), "--wgs-chrx", str(tmp_path / "wgs_chrx.tsv"),
                                     "--eqtls", str(tmp_path / "eqtls.tsv"), "--clinical", str(tmp_path / "clinical.csv"),
                                     "--biospecimen", str(tmp_path / "biospecimen.csv"),
                                     "--actions", str(tmp_path / "actions.tsv"), "--output", str(tmp_path / "out")])
    _script().main()
    people = pd.read_csv(tmp_path / "out" / "restricted" / "array_vs_wgs_by_person.tsv", sep="\t",
                         dtype={"projid": str, "best_wgs_projid": str}).set_index("projid")
    assert (people.drop(index=["5", "38"]).loc[[str(p) for p in range(35) if p != 5], "status"] == "same person").all()
    assert people.loc["0", "self_r"] > 0.95
    assert people.loc["5", "status"] == "matches no WGS"
    assert people.loc["38", "status"] == "no WGS under this projid, matches a WGS sample"
    assert people.loc["38", "best_wgs_projid"] == "5"
    assert (people.loc[[str(p) for p in range(36, 40) if p != 38], "status"] == "no WGS").all()
    assert people.loc["5", "wgs_dna_sex"] == "F" and people.loc["5", "recorded_sex"] == "M"  # 38's DNA in 5's WGS
    sex = pd.read_csv(tmp_path / "out" / "dna_sex_summary.tsv", sep="\t").set_index("dna_sex")
    assert sex.loc["array_dna_sex", "n_differ"] == 0 and sex.loc["wgs_dna_sex", "n_differ"] == 1
    renamed = pd.read_csv(tmp_path / "out" / "restricted" / "wgs.by_individual.tsv.gz", sep="\t", nrows=1)
    assert list(renamed.columns[3:6]) == ["R000", "R001", "R002"]
    flagged = pd.read_csv(tmp_path / "out" / "restricted" / "flagged_donors_wgs.tsv", sep="\t").iloc[0]
    assert flagged["label_array_vs_wgs"] == "matches no WGS"
    assert flagged["best_array_vs_wgs"] == "no WGS under this projid, matches a WGS sample"
