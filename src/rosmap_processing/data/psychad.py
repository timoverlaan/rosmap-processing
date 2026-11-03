"""PsychAD (MSSM_PsychAD, Lee et al. 2024) -> a ROSMAP-style count file.

`merged_final_clean.h5ad` holds the authors' QC'd nuclei (6.3 million x 34,890
genes named by Ensembl ID) with normalized values in X and the counts in
layers['counts']. `prepare_file` writes one count file with those counts in X
(int32, streamed in chunks: data/h5ad_counts.py) and the file's obs plus:

- individualID (the file's SubID_export_synapse) and the donor metadata: the
  NPS-AD individual metadata, the paper's tiered diagnoses, the hashing set;
- channel: the hashing set and the pool replicate from the cell name
  (`<donor>-<replicate>-<barcode>-<index>`). A set of about six donors runs in
  two channels and nearly every donor is in both, so channel is a batch
  within donors, the hashing set one between donors;
- the ROSMAP clinical columns for the Rush donors, joined on individualID,
  which they share with ROSMAP (the ROS and MAP participants, 65 of 152).
"""

from pathlib import Path
from typing import Optional, Union

import h5py
import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from .h5ad_counts import copy_counts

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem

logger = get_logger(__name__)

H5AD = "merged_final_clean.h5ad"
INDIVIDUAL_METADATA = "NPS-AD_individual_metadata.csv"
DIAGNOSES = "NPS-AD_neuropathology_and_diagnosis_metadata.csv"
HASHING = "NPS-AD_nucleus_hashing_metadata.csv"
DONOR_COLUMN = "SubID_export_synapse"
COUNTS_KEY = "layers/counts"

INDIVIDUAL_COLUMNS = ["cohort", "dataContributionGroup", "sex", "race", "isHispanic", "ageDeath", "PMI",
                      "apoeGenotype", "apoe4Status", "CDR", "amyCerad", "amyAny", "amyThal", "Braak", "bScore",
                      "PD", "FTD", "MS", "SCZ", "BD", "PTSD", "ADHD", "OCD"]
DIAGNOSIS_COLUMNS = ["Diagnosis", "Tier1_crossDis_dx", "Tier2_AD_dx", "Tier3_NPS_dx", "dxPC1"]
# ROSMAP clinical columns that would clash (race is coded differently there)
CLINICAL_RENAME = {"race": "rosmap_race"}


def donor_table(download_dir: Union[str, Path], clinical: Optional[Union[str, Path]] = None) -> pd.DataFrame:
    """One row per individualID: individual metadata, diagnoses, hashing set
    and (Rush ROS and MAP donors) the ROSMAP clinical columns."""
    d = Path(download_dir)
    donors = pd.read_csv(d / INDIVIDUAL_METADATA, usecols=["individualID"] + INDIVIDUAL_COLUMNS)
    dx = pd.read_csv(d / DIAGNOSES, usecols=["IndividualID"] + DIAGNOSIS_COLUMNS)
    donors = donors.merge(dx.rename(columns={"IndividualID": "individualID"}), on="individualID", how="left")
    hashing = pd.read_csv(d / HASHING, usecols=["individualID", "Set_name"])
    donors = donors.merge(hashing.rename(columns={"Set_name": "hashing_set"}), on="individualID", how="left")
    if clinical is not None:
        table = pd.read_csv(clinical, low_memory=False).rename(columns=CLINICAL_RENAME)
        donors = donors.merge(table, on="individualID", how="left")
    if donors["individualID"].duplicated().any():
        raise ValueError("A donor appears twice in the PsychAD metadata")
    return donors


def annotate_obs(obs: pd.DataFrame, download_dir: Union[str, Path],
                 clinical: Optional[Union[str, Path]] = None) -> pd.DataFrame:
    """obs plus individualID, the donor columns and channel (one row per cell, same order)."""
    donor = obs[DONOR_COLUMN].astype(str).to_numpy()
    codes, uniques = pd.factorize(donor)
    table = donor_table(download_dir, clinical).set_index("individualID").reindex(uniques)
    missing = int(table["sex"].isna().sum())
    if missing:
        logger.warning(f"{missing} of {len(uniques)} donors are not in the individual metadata")
    out = obs.copy()
    out["individualID"] = pd.Categorical(donor)
    for col in table.columns:
        values = table[col].to_numpy()[codes]
        if table[col].dtype == object:
            # text, or values with gaps from the join: categories as text
            text = pd.Series(values, dtype=object)
            values = pd.Categorical(text.where(text.isna(), text.astype(str)))
        out[col] = values
    replicate = pd.Series(obs.index.astype(str)).str.rsplit("-", n=3, expand=True)[1].to_numpy()
    channel = np.where(pd.isna(out["hashing_set"]), None,
                       out["hashing_set"].astype(str).to_numpy() + "_" + replicate.astype(str))
    out["channel"] = pd.Categorical(channel)
    logger.info(f"{len(uniques)} donors; {int(table['hashing_set'].notna().sum())} with a hashing set, "
                f"{out['channel'].nunique()} channels; "
                f"{int(table['Tier1_crossDis_dx'].notna().sum())} with a Tier 1 diagnosis")
    if clinical is not None:
        logger.info(f"{int(table['cogdx'].notna().sum()) if 'cogdx' in table else 0} donors have ROSMAP "
                    f"clinical data")
    return out


def prepare_file(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    download_dir: Union[str, Path],
    clinical: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    chunk_rows: int = 100_000,
) -> None:
    """merged_final_clean.h5ad -> count h5ad (int32 counts in X, annotated obs)."""
    from ..utils.provenance import provenance_record

    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path} (pass --overwrite to replace it)")
    with h5py.File(input_path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
    logger.info(f"Read obs and var of {input_path}: {len(obs)} cells, {len(var)} genes")
    obs = annotate_obs(obs, download_dir, clinical)
    uns = {"provenance_prepare": provenance_record(
        "data psychad", {"input": input_path, "counts": COUNTS_KEY, "download_dir": download_dir,
                         "clinical": clinical})}
    copy_counts(input_path, output_path, obs, var, key=COUNTS_KEY, uns=uns, chunk_rows=chunk_rows)
