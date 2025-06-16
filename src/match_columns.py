import anndata as ad
import pandas as pd
import numpy as np

import argparse


parser = argparse.ArgumentParser(description="Script converting column names and types between SeaAD and ROSMAP data. \
                                  It also adds the Wang label column that we try to predict")
parser.add_argument(
    "path",
    type=str,
    required=True,
    help="Path to the input h5ad file.",
)
parser.add_argument(
    "--type",
    type=str,
    required=True,
    help="Type of the input data, either \"ROSMAP\", \"ROSMAP_MIT\" or \"SeaAD\".",
    choices=["ROSMAP", "ROSMAP_MIT", "SeaAD"],
)
parser.add_argument(
    "--output",
    type=str,
    help="Path to the output h5ad file. This will contain the columns in both formats: \
        ROSMAP format, which is numerically valued, and lowercase. \
        And the SeaAD format, which contains mainly strings and is capitalized.",
)
parser.add_argument(
    "--inplace",
    action="store_true",
    help="If set, the input file will be modified in place. This needs to be set if --output is not provided.",
)
args = parser.parse_args()


CLASS_NAME = "Class"
SUBCLASS_NAME = "Celltype" 
SUPERTYPE_NAME = "Subtype"  # the "state" or cluster


def convert_rosmap(adata: ad.AnnData) -> None:
    """
    Adds the missing SeaAD-style columns (in-place) to the AnnData object.

    - mit: bool, True if MIT data.
    """

    # ROSMAP (without metadata) has these:
    # 'nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT', 'batch', 'individualID', 'DoubletFinder.score', 'subset', 'class', 'state'

    adata.obs.rename(columns={
        "individualID": "Donor ID",
        # "Sex": "Sex",
        # "Overall AD neuropathological Change": "ADNC",
        # "Thal": "Thal",
        # "CERAD score": "CERAD",
        # "Overall CAA score": "CAA",
        # "LATE": "LATE",
        "class": CLASS_NAME,
        "subset": SUBCLASS_NAME,
        "state": SUPERTYPE_NAME,
    }, inplace=True)

    # TODO: also do the metadata columns



def convert_rosmap_mit(adata: ad.AnnData, cellclass: str, subclass: str) -> None:
    """
    Adds the missing SeaAD-style columns (in-place) to the AnnData object.
    This is a special case for the MIT data, which has different column names.
    """

    adata.obs.rename(columns={
        "projid": "Donor ID",  # TODO: this might be a problem
        "cell_type_high_resolution": SUPERTYPE_NAME,
    }, inplace=True)

    adata.obs[CLASS_NAME] = cellclass
    adata.obs[SUBCLASS_NAME] = subclass

    # TODO: also do the metadata columns


def convert_seaad(adata: ad.AnnData) -> None:
    """
    Adds the missing ROSMAP-style columns (in-place) to the AnnData object.
    """

    adata.obs.rename(columns={
        "Donor ID": "Donor ID",
        "Sex": "Sex",
        "Overall AD neuropathological Change": "ADNC",
        "Thal": "Thal",
        "CERAD score": "CERAD",
        "Overall CAA score": "CAA",
        "LATE": "LATE",
        "Class": CLASS_NAME,
        "Subclass": SUBCLASS_NAME,
        "Supertype": SUPERTYPE_NAME,
    }, inplace=True)
    
    adata.obs["braaksc"] = adata.obs["Braak"].map({
        "Braak 0": 0,
        "Braak I": 1,
        "Braak II": 2,
        "Braak III": 3,
        "Braak IV": 4,
        "Braak V": 5,
        "Braak VI": 6,
    })

    # CERAD Score in SeaAD is the "amount of neuritic plaques", so we need to invert the mapping to correspond to the ROSMAP ceradsc column
    adata.obs["ceradsc"] = adata.obs["CERAD"].map({  # ['Moderate', 'Absent', 'Reference', 'Frequent', 'Sparse']
        "Absent": 4,
        "Sparse": 3,
        "Moderate": 2,
        "Frequent": 1,
        "Reference": 4,
    })  

    adata.obs["thalsc"] = adata.obs["Thal"].map({
        "Thal 0": 0,
        "Thal 1": 1,
        "Thal 2": 2,
        "Thal 3": 3,
        "Thal 4": 4,
        "Thal 5": 5,
        "Thal 6": 6,
    })

    # Derive the Wang labels
    adata.obs["Wang"] = None
    adata.obs["Wang_intermediate"] = True  # These will be excluded

    # SeaAD doesn't have the cogdx column, so we derive it from the cognitive status
    adata.obs.loc[(adata.obs["Cognitive Status"] == "Dementia") & (adata.obs["braaksc"] >= 4) & (adata.obs["ceradsc"] <= 2), "Wang"] = "AD"
    adata.obs.loc[(adata.obs["Cognitive Status"] == "Dementia") & (adata.obs["braaksc"] >= 4) & (adata.obs["ceradsc"] <= 2), "Wang_intermediate"] = False

    adata.obs.loc[(adata.obs["Cognitive Status"] != "Dementia") & (adata.obs["braaksc"] <= 3) & (adata.obs["ceradsc"] >= 3), "Wang"] = "Healthy"
    adata.obs.loc[(adata.obs["Cognitive Status"] != "Dementia") & (adata.obs["braaksc"] <= 3) & (adata.obs["ceradsc"] >= 3), "Wang_intermediate"] = False



if __name__ == "__main__":

    print(f"Loading data from {args.path}...")
    adata = ad.read_h5ad(args.path)

    if not args.inplace and not args.output:
        raise ValueError("Either --inplace or --output must be provided.")

    if args.type == "ROSMAP":
        print("Converting ROSMAP data to SeaAD format...")
        convert_rosmap(adata)

    elif args.type == "ROSMAP_MIT":
    
        # Then cellclass and subclass need to be added as an argument
        if not args.cellclass or not args.cellsubclass:
            raise ValueError("For ROSMAP_MIT data, --cellclass and --cellsubclass must be provided.")

        print("Converting ROSMAP MIT data to SeaAD format...")
        convert_rosmap_mit(adata, cellclass=args.cellclass, cellsubclass=args.cellsubclass)

    else:
        print("Converting SeaAD data to ROSMAP format...")
        convert_seaad(adata)

    print("Finished converting metadata columns.\n")
    
    n_Wang_samples = adata[~adata.obs["Wang_intermediate"]].obs["Donor ID"].nunique()
    print(f"Number of samples with Wang labels (so remaining samples after selecting extremes): {n_Wang_samples}")

    print(f"\nSaving data to {args.output}...")
    adata.write_h5ad(args.output, compression="gzip")
    print("Done!")