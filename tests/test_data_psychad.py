"""Tests for the PsychAD conversion (data/psychad.py), on a made-up file."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from rosmap_processing.data import psychad


@pytest.fixture
def source(tmp_path):
    d = tmp_path / "psychad"
    d.mkdir()
    donors = ["AMPAD_MSSM_0000000001", "AMPAD_MSSM_0000000001", "R1000001", "R1000001", "AMPAD_HBCC_0000000009"]
    cells = [f"{donor}-{rep}-ACGT{i}-{rep - 1}" for i, (donor, rep) in enumerate(zip(donors, [1, 2, 1, 2, 1]))]
    counts = np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0], [0, 0, 1], [2, 2, 2]])
    adata = ad.AnnData(
        X=sp.csr_matrix(np.log1p(counts).astype(np.float32)),
        obs=pd.DataFrame({psychad.DONOR_COLUMN: pd.Categorical(donors), "class": ["EN", "IN", "EN", "Oligo", "EN"]},
                         index=cells),
        var=pd.DataFrame(index=["ENSG01", "ENSG02", "ENSG03"]),
    )
    adata.layers["counts"] = sp.csr_matrix(counts.astype(np.float32))
    adata.write_h5ad(d / psychad.H5AD)
    ids = ["AMPAD_MSSM_0000000001", "R1000001", "AMPAD_HBCC_0000000009"]
    individual = pd.DataFrame({"individualID": ids})
    for col in psychad.INDIVIDUAL_COLUMNS:
        individual[col] = ["90+", "77", "81"] if col == "ageDeath" else ["a", "b", "c"]
    individual.to_csv(d / psychad.INDIVIDUAL_METADATA, index=False)
    dx = pd.DataFrame({"IndividualID": ids})
    for col in psychad.DIAGNOSIS_COLUMNS:
        dx[col] = ["AD", "CTRL", np.nan]
    dx.to_csv(d / psychad.DIAGNOSES, index=False)
    # the HBCC donor has no hashing set
    pd.DataFrame({"individualID": ids[:2], "Set_name": ["NPSAD_001_A", "NPSAD_002_A"]}).to_csv(
        d / psychad.HASHING, index=False)
    pd.DataFrame({"projid": [5], "individualID": ["R1000001"], "cogdx": [4], "race": [2]}).to_csv(
        tmp_path / "clinical.csv", index=False)
    return d, counts


def test_prepare_file_counts_and_donor_columns(source, tmp_path):
    d, counts = source
    psychad.prepare_file(d / psychad.H5AD, tmp_path / "out.h5ad", d, clinical=tmp_path / "clinical.csv",
                         chunk_rows=2)
    out = ad.read_h5ad(tmp_path / "out.h5ad")
    assert out.X.dtype == np.int32
    np.testing.assert_array_equal(out.X.toarray(), counts)
    obs = out.obs
    assert obs["individualID"].tolist() == obs[psychad.DONOR_COLUMN].astype(str).tolist()
    assert obs["ageDeath"].tolist() == ["90+", "90+", "77", "77", "81"]
    assert obs["Tier1_crossDis_dx"].tolist()[:3] == ["AD", "AD", "CTRL"] and pd.isna(obs["Tier1_crossDis_dx"].iloc[4])
    # channel = hashing set + pool replicate; none without a hashing set
    assert obs["channel"].tolist()[:4] == ["NPSAD_001_A_1", "NPSAD_001_A_2", "NPSAD_002_A_1", "NPSAD_002_A_2"]
    assert pd.isna(obs["channel"].iloc[4])
    # ROSMAP clinical columns only for the Rush donor
    assert obs["cogdx"].tolist()[2:4] == [4.0, 4.0] and obs["cogdx"].isna().sum() == 3
    assert "rosmap_race" in obs


def test_prepare_file_never_overwrites(source, tmp_path):
    d, _ = source
    (tmp_path / "out.h5ad").write_text("x")
    with pytest.raises(FileExistsError):
        psychad.prepare_file(d / psychad.H5AD, tmp_path / "out.h5ad", d)
