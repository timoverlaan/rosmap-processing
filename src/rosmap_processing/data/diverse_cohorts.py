"""AMP-AD DiverseCohorts 10x multiome, gene expression -> ROSMAP-style count files.

The 2025 release has one `gex_mat_<class>.rda` per cell class: a dgCMatrix
(genes x cells) whose cells are named `<pool>_<barcode>`, converted to .h5 by
src/convert_R/rdata_to_h5.R. Each library pools three specimens (donor x
region) of different donors, so which donor and region a cell is from is only
in the cell table `rosmap_diversity_snmultiome_clusterdata_20241107.csv`
(pool as `batch`, `barcode`, `individualID`, `specimenID`, `region`, cell
classes, and Rush's cogdx/Braak/CERAD). The 2023 matrices have no such table
for their cells and are not used.

`annotate` keeps the cells that are in the cell table and from the chosen
region, and adds donor metadata: the DiverseCohorts individual metadata (age,
APOE, race, cohort, ... for every donor), projid (Rush mapping) and, where the
donor is a ROS or MAP participant, the ROSMAP clinical columns. A Rush donor's
individualID is the ROSMAP one, so the clinical CSV is joined on it; most
DiverseCohorts donors are from MARS and CLINCORE, which it does not cover.
Rush's cogdx, Braak and CERAD are in the cell table for every donor.
"""

from pathlib import Path
from typing import Optional, Union

import anndata as ad
import numpy as np
import pandas as pd

from ..pipelines.mathys2023_for_adprs import normalize_projid
from ..utils.logging import get_logger

logger = get_logger(__name__)

CELL_TABLE = "rosmap_diversity_snmultiome_clusterdata_20241107.csv"
INDIVIDUAL_METADATA = "AMP-AD_DiverseCohorts_individual_metadata.csv"
PROJID_MAPPING = "AMP-AD_DiverseCohorts_Rush_projid_mapping.csv"
DLPFC = "dorsolateral prefrontal cortex"

CELL_COLUMNS = ["individualID", "specimenID", "region", "batch", "broadcellclass", "subcellclass",
                "cogdx_rush", "braak_rush", "cerad_rush", "sex"]
# From the individual metadata; 'sex' is already in the cell table
DONOR_COLUMNS = ["dataContributionGroup", "cohort", "race", "isHispanic", "ageDeath", "PMI",
                 "apoeGenotype", "amyCerad", "Braak", "reag", "ADoutcome"]
# ROSMAP clinical columns that would clash: race is coded differently there;
# its projid is kept to check the Rush mapping
CLINICAL_RENAME = {"race": "rosmap_race", "projid": "rosmap_projid"}


def read_cell_table(path: Union[str, Path]) -> pd.DataFrame:
    """The cell table, indexed by `<pool>_<barcode>` (the matrices' cell names)."""
    cells = pd.read_csv(path, usecols=CELL_COLUMNS + ["barcode"], low_memory=False)
    cells.index = pd.Index(cells["batch"].astype(str) + "_" + cells["barcode"].astype(str))
    if not cells.index.is_unique:
        raise ValueError(f"{path}: {int(cells.index.duplicated().sum())} cells appear twice")
    return cells.drop(columns="barcode")


def donor_table(download_dir: Union[str, Path], clinical: Optional[Union[str, Path]] = None) -> pd.DataFrame:
    """One row per DiverseCohorts individualID: individual metadata, projid
    (Rush mapping, else the clinical CSV's; as text without leading zeros) and
    the ROSMAP clinical columns (joined on individualID)."""
    download_dir = Path(download_dir)
    donors = pd.read_csv(download_dir / INDIVIDUAL_METADATA, usecols=["individualID"] + DONOR_COLUMNS)
    mapping = pd.read_csv(download_dir / PROJID_MAPPING)
    if mapping["individualID"].duplicated().any():
        raise ValueError(f"{PROJID_MAPPING}: an individualID has several projids")
    mapping["projid"] = normalize_projid(mapping["projid"])
    donors = donors.merge(mapping[["individualID", "projid"]], on="individualID", how="left")
    if clinical is not None:
        table = pd.read_csv(clinical, low_memory=False).rename(columns=CLINICAL_RENAME)
        table["rosmap_projid"] = normalize_projid(table["rosmap_projid"])
        if table["individualID"].duplicated().any():
            raise ValueError(f"{clinical}: duplicate individualID")
        donors = donors.merge(table, on="individualID", how="left")
        both = donors["projid"].notna() & donors["rosmap_projid"].notna()
        differ = int((donors.loc[both, "projid"] != donors.loc[both, "rosmap_projid"]).sum())
        if differ:
            raise ValueError(f"The Rush mapping and {clinical} give different projids for {differ} donors")
        donors["projid"] = donors["projid"].fillna(donors.pop("rosmap_projid"))
    return donors


def annotate(
    adata: ad.AnnData,
    download_dir: Union[str, Path],
    clinical: Optional[Union[str, Path]] = None,
    region: Optional[str] = DLPFC,
    cells: Optional[pd.DataFrame] = None,
) -> ad.AnnData:
    """The cells of `adata` that are in the cell table (and `region`, unless
    None), with cell and donor metadata in obs. Logs counts only."""
    if cells is None:
        cells = read_cell_table(Path(download_dir) / CELL_TABLE)
    in_table = adata.obs_names.isin(cells.index)
    logger.info(f"{int(in_table.sum())} of {adata.n_obs} cells are in the cell table")
    keep = in_table.copy()
    if region is not None:
        cell_region = cells["region"].reindex(adata.obs_names).to_numpy()
        keep &= cell_region == region
        logger.info(f"{int(keep.sum())} of them are from the {region}")
    adata = adata[keep].copy()

    obs = cells.loc[adata.obs_names].copy()
    donors = donor_table(download_dir, clinical)
    missing = set(obs["individualID"]) - set(donors["individualID"])
    if missing:
        logger.warning(f"{len(missing)} donors are not in the individual metadata")
    obs = obs.reset_index(names="cell").merge(donors, on="individualID", how="left").set_index("cell")
    obs.index.name = None
    if len(obs) != adata.n_obs:
        raise ValueError("The donor metadata has duplicate individualIDs")
    n_donors = obs["individualID"].nunique()
    with_projid = obs.loc[obs["projid"].notna(), "individualID"].nunique()
    logger.info(f"{n_donors} donors; {with_projid} with a projid")
    if clinical is not None:
        with_clinical = obs.loc[obs["cogdx"].notna(), "individualID"].nunique() if "cogdx" in obs else 0
        logger.info(f"{with_clinical} donors have ROSMAP clinical data")
        agree = (obs["cogdx"].astype(float) == obs["cogdx_rush"].astype(float)) | obs["cogdx"].isna()
        logger.info(f"cogdx and cogdx_rush differ for {obs.loc[~agree, 'individualID'].nunique()} donors")
    # Text columns (including ages like "90+") as categoricals, which h5ad
    # stores with their missing values
    adata.obs = obs.astype({c: "category" for c in obs.columns if obs[c].dtype == object})
    return adata


def prepare_file(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    download_dir: Union[str, Path],
    clinical: Optional[Union[str, Path]] = None,
    region: Optional[str] = DLPFC,
    overwrite: bool = False,
) -> None:
    """One class: the .h5 from rdata_to_h5.R -> annotated count h5ad (int32 counts in X)."""
    from ..utils.provenance import provenance_record
    from .prepare import prepare_counts
    from .rmatrix import read_matrix_h5

    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path} (pass --overwrite to replace it)")
    adata = annotate(read_matrix_h5(input_path), download_dir, clinical, region)
    adata = prepare_counts(adata, categorical_columns=[])
    adata.uns["provenance_prepare"] = provenance_record(
        "data diverse-cohorts",
        {"input": input_path, "download_dir": download_dir, "clinical": clinical, "region": region},
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_name(output_path.name + ".partial")
    adata.write_h5ad(tmp, compression="gzip")
    tmp.replace(output_path)
    logger.info(f"Wrote {output_path}: {adata.n_obs} cells x {adata.n_vars} genes")
