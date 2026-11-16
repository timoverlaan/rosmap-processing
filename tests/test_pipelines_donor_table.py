"""Tests of the global donor table on synthetic count files."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from rosmap_processing.pipelines.donor_table import load_spec, run_donor_table


def _write(path, obs_by_donor, cells=3):
    """obs_by_donor: donor -> dict of obs values, repeated over its cells."""
    rows = []
    for donor, values in obs_by_donor.items():
        rows += [{**values, "_cell": f"{donor}_{i}"} for i in range(cells)]
    obs = pd.DataFrame(rows).set_index("_cell")
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = obs[c].astype("category")
    adata = ad.AnnData(sp.csr_matrix(np.ones((len(obs), 3), dtype=np.float32)), obs=obs,
                       var=pd.DataFrame(index=["G0", "G1", "G2"]))
    adata.write_h5ad(path)


ROSMAP_VARS = {
    "projid": {"column": "projid"},
    "sex": {"column": "msex", "recode": {"0": "female", "1": "male"}},
    "age_death": {"column": "age_death", "numeric": True},
    "apoe_genotype": {"column": "apoe_genotype"},
    "dementia": {"column": "cogdx", "recode": {"1": "no", "2": "no", "3": "no", "4": "yes", "5": "yes"}},
    "braak": {"column": "braaksc"},
}


def _setup(tmp_path):
    # dataset a: persons R1-R4; dataset b: R3-R6 (R3, R4 shared), R4's age differs
    a = {f"P{i}": {"projid": 100 + i, "individualID": f"R{i}", "msex": i % 2, "age_death": "90+" if i == 1 else
                   str(80 + i), "apoe_genotype": 34.0 if i == 2 else 33.0, "cogdx": 4 if i % 2 else 1,
                   "braaksc": float(i)} for i in range(1, 5)}
    b = {f"R{i}": {"projid": f"{100 + i:08d}", "individualID": f"R{i}", "sex": "male" if i % 2 else "female",
                   "ageDeath": str(80 + i + (5 if i == 4 else 0)), "cdr": 1.0 if i == 6 else 0.5,
                   "cogdx": np.nan if i > 4 else (4 if i % 2 else 1), "Braak": "Stage II"}
         for i in range(3, 7)}
    # dataset c: its own donors, race as checkboxes (D2 has two boxes checked)
    c = {f"D{i}": {"Donor ID": f"D{i}", "White": "Checked" if i != 3 else "Unchecked",
                   "Asian": "Checked" if i >= 2 else "Unchecked", "Sex": "Female",
                   "Status": "Reference" if i == 1 else "Dementia"} for i in range(1, 4)}
    _write(tmp_path / "a.h5ad", a)
    _write(tmp_path / "b.h5ad", b)
    _write(tmp_path / "c.h5ad", c)
    spec = {
        "output_dir": str(tmp_path / "out"),
        "datasets": {
            "a": {"counts": str(tmp_path / "a.h5ad"), "donor": "projid", "person": "individualID",
                  "variables": ROSMAP_VARS},
            "b": {"counts": str(tmp_path / "b.h5ad"), "donor": "individualID", "person": "individualID",
                  "variables": {
                      "projid": {"column": "projid"}, "sex": {"column": "sex"},
                      "age_death": {"column": "ageDeath", "numeric": True}, "cdr": {"column": "cdr"},
                      "dementia": [{"column": "cogdx", "recode": {"1": "no", "4": "yes"}},
                                   {"column": "cdr", "recode": {"0": "no", "0.5": "no", "1": "yes"}}],
                      "braak": {"column": "Braak", "recode": {"Stage II": 2}}}},
            "c": {"counts": str(tmp_path / "c.h5ad"), "donor": "Donor ID", "person": "Donor ID",
                  "variables": {
                      "sex": {"column": "Sex", "recode": {"Female": "female"}},
                      "race": {"checked": "Checked", "checkboxes": {"White": "White", "Asian": "Asian"}},
                      "dementia": {"column": "Status", "recode": {"Dementia": "yes", "Reference": "no"}}}},
        },
    }
    path = tmp_path / "spec.yaml"
    path.write_text(yaml.safe_dump(spec))
    return path


def test_donor_table_merges_persons_across_datasets(tmp_path):
    out = run_donor_table(_setup(tmp_path)).parent
    table = pd.read_csv(out / "donor_table.tsv", sep="\t").set_index("person_id")
    assert len(table) == 9   # R1-R6 and D1-D3
    assert table.loc["R3", "datasets"] == "a+b" and table.loc["R5", "datasets"] == "b"
    assert float(table.loc["R3", "a_donor"]) == 103 and table.loc["R3", "b_donor"] == "R3"
    assert table.loc["R3", "a_n_nuclei"] == 3
    # padded and float projids match
    assert table.loc["R3", "projid"] == 103 and table.loc["R5", "projid"] == 105
    # the first dataset with a value wins
    assert table.loc["R4", "age_death"] == 84 and table.loc["R4", "age_death_from"] == "a"
    assert table.loc["R1", "age_death"] == 90 and table.loc["R1", "age_capped"] == 90
    assert table.loc["R2", "apoe4"] == 1 and table.loc["R1", "apoe4"] == 0
    # second source fills in: CDR for donors without cogdx
    assert table.loc["R6", "dementia"] == "yes" and table.loc["R5", "dementia"] == "no"
    assert table.loc["R5", "braak"] == 2
    # checkboxes: one box -> that race, two -> Multiple
    assert table.loc["D1", "race"] == "White" and table.loc["D3", "race"] == "Asian"
    assert table.loc["D2", "race"] == "Multiple"
    assert table.loc["D1", "dementia"] == "no"

    overlap = pd.read_csv(out / "overlap.tsv", sep="\t")
    regions = overlap[overlap["kind"] == "region"].set_index("datasets")["n_persons"]
    assert regions["a"] == 2 and regions["b"] == 2 and regions["a+b"] == 2 and regions["c"] == 3
    assert overlap.query("kind == 'n_datasets' and datasets == '2'")["n_persons"].item() == 2

    conflicts = pd.read_csv(out / "conflicts.tsv", sep="\t").set_index("variable")
    assert conflicts.loc["age_death", "n_compared"] == 2 and conflicts.loc["age_death", "n_differ"] == 1
    assert conflicts.loc["projid", "n_differ"] == 0 and conflicts.loc["sex", "n_differ"] == 0

    summary = pd.read_csv(out / "summary.tsv", sep="\t")
    assert set(summary["metric"]) <= {"n", "n_with_value", "count", "q1", "median", "q3"}
    count = summary.query("dataset == 'all' and variable == 'dementia' and level == 'yes'")["value"].item()
    assert count == (table["dementia"] == "yes").sum()
    # aggregate outputs hold no person IDs
    for name in ("overlap.tsv", "summary.tsv", "conflicts.tsv"):
        text = (out / name).read_text()
        assert not any(p in text for p in ("R1", "R4", "P3", "D2"))


def test_unrecoded_labels_are_an_error(tmp_path):
    path = _setup(tmp_path)
    spec = yaml.safe_load(path.read_text())
    spec["datasets"]["b"]["variables"]["sex"] = {"column": "ageDeath"}   # ages are not sexes
    path.write_text(yaml.safe_dump(spec))
    with pytest.raises(ValueError, match="sex: values"):
        run_donor_table(path)


def test_spec_rejects_unknown_variables(tmp_path):
    path = _setup(tmp_path)
    spec = yaml.safe_load(path.read_text())
    spec["datasets"]["a"]["variables"]["height"] = {"column": "msex"}
    path.write_text(yaml.safe_dump(spec))
    with pytest.raises(ValueError, match="unknown variable height"):
        load_spec(path)


def test_existing_output_needs_overwrite(tmp_path):
    path = _setup(tmp_path)
    run_donor_table(path)
    with pytest.raises(FileExistsError):
        run_donor_table(path)
    run_donor_table(path, overwrite=True)
