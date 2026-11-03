"""Tests for the DiverseCohorts annotation (data/diverse_cohorts.py), on made-up tables."""

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from rosmap_processing.data import diverse_cohorts as dc


@pytest.fixture
def download_dir(tmp_path):
    d = tmp_path / "dc"
    d.mkdir()
    # Pool VR001 holds donors R1 (DLPFC) and R2 (STG); VR002 holds R3 (DLPFC)
    pd.DataFrame({
        "Unnamed: 0": range(5),
        "individualID": ["R1", "R1", "R2", "R3", "R3"],
        "specimenID": ["R1_MF", "R1_MF", "R2_ST", "R3_MF", "R3_MF"],
        "region": [dc.DLPFC, dc.DLPFC, "superior temporal gyrus", dc.DLPFC, dc.DLPFC],
        "batch": ["VR001", "VR001", "VR001", "VR002", "VR002"],
        "barcode": ["AAAC-1", "AAAG-1", "AAAT-1", "AAAC-1", "CCCA-1"],
        "nCount_RNA": [10, 20, 30, 40, 50],
        "broadcellclass": ["Astro"] * 5,
        "subcellclass": ["Astro_A"] * 5,
        "cogdx_rush": [1, 1, 4, 4, 4],
        "braak_rush": [2, 2, 5, 4, 4],
        "cerad_rush": [4, 4, 1, 2, 2],
        "sex": ["female", "female", "male", "male", "male"],
    }).to_csv(d / dc.CELL_TABLE, index=False)
    individual = pd.DataFrame({"individualID": ["R1", "R2", "R3"], "sex": ["female", "male", "male"]})
    for col in dc.DONOR_COLUMNS:
        individual[col] = ["90+", "85", "77"] if col == "ageDeath" else ["x", "y", "z"]
    individual.to_csv(d / dc.INDIVIDUAL_METADATA, index=False)
    pd.DataFrame({"individualID": ["R1", "R3"], "projid": [12345, 20000001]}).to_csv(d / dc.PROJID_MAPPING, index=False)
    # R3 is not in the clinical CSV (a MARS donor, say); R2 has no Rush mapping
    clinical = pd.DataFrame({"projid": ["00012345", "00000777"], "individualID": ["R1", "R2"],
                             "cogdx": [1.0, 5.0], "msex": [0, 1], "race": [1, 2]})
    clinical.to_csv(tmp_path / "clinical.csv", index=False)
    return d


def _matrix(path, cells, n_genes=3):
    counts = np.arange(len(cells) * n_genes).reshape(len(cells), n_genes) % 4
    m = sp.csc_matrix(counts.T.astype(np.float64))
    with h5py.File(path, "w") as f:
        f["data"], f["indices"], f["indptr"] = m.data, m.indices, m.indptr
        f["shape"] = np.array([len(cells), n_genes])
        f["genes"] = np.array([f"G{j}" for j in range(n_genes)], dtype="S")
        f["cells"] = np.array(cells, dtype="S")
    return counts


def test_annotate_keeps_region_cells_with_donor_metadata(download_dir, tmp_path):
    # VR001_TTTT-1 is not in the cell table (e.g. a doublet)
    cells = ["VR001_AAAC-1", "VR001_AAAT-1", "VR002_AAAC-1", "VR001_TTTT-1", "VR001_AAAG-1"]
    counts = _matrix(tmp_path / "m.h5", cells)
    dc.prepare_file(tmp_path / "m.h5", tmp_path / "out.h5ad", download_dir, clinical=tmp_path / "clinical.csv")
    out = ad.read_h5ad(tmp_path / "out.h5ad")
    assert out.obs_names.tolist() == ["VR001_AAAC-1", "VR002_AAAC-1", "VR001_AAAG-1"]
    np.testing.assert_array_equal(out.X.toarray(), counts[[0, 2, 4]])
    assert out.obs["individualID"].tolist() == ["R1", "R3", "R1"]
    # projid without zero padding; the clinical columns joined on individualID (R3 has none)
    assert out.obs["projid"].tolist() == ["12345", "20000001", "12345"]
    np.testing.assert_array_equal(out.obs["cogdx"].to_numpy(), [1.0, np.nan, 1.0])
    assert "rosmap_race" in out.obs and out.obs["ageDeath"].tolist() == ["90+", "77", "90+"]


def test_projid_from_clinical_when_not_mapped(download_dir, tmp_path):
    donors = dc.donor_table(download_dir, tmp_path / "clinical.csv").set_index("individualID")
    assert donors.loc["R2", "projid"] == "777" and donors.loc["R3", "projid"] == "20000001"


def test_conflicting_projids_are_an_error(download_dir, tmp_path):
    clinical = pd.read_csv(tmp_path / "clinical.csv", dtype={"projid": str})
    clinical.loc[0, "projid"] = "00099999"
    clinical.to_csv(tmp_path / "clinical.csv", index=False)
    with pytest.raises(ValueError, match="different projids"):
        dc.donor_table(download_dir, tmp_path / "clinical.csv")


def test_annotate_all_regions(download_dir, tmp_path):
    cells = ["VR001_AAAC-1", "VR001_AAAT-1"]
    _matrix(tmp_path / "m.h5", cells)
    dc.prepare_file(tmp_path / "m.h5", tmp_path / "out.h5ad", download_dir, region=None)
    out = ad.read_h5ad(tmp_path / "out.h5ad")
    assert out.n_obs == 2 and "cogdx" not in out.obs


def test_cell_table_must_be_unique(download_dir):
    table = pd.read_csv(download_dir / dc.CELL_TABLE)
    pd.concat([table, table.iloc[:1]]).to_csv(download_dir / dc.CELL_TABLE, index=False)
    with pytest.raises(ValueError, match="twice"):
        dc.read_cell_table(download_dir / dc.CELL_TABLE)
