import anndata as ad
import pandas as pd
import numpy as np

import argparse


parser = argparse.ArgumentParser(description="Script to add metadata to an AnnData object.")
parser.add_argument(
    "path",
    type=str,
    required=True,
    help="Path to the input h5ad file. Note: this will be overwritten!",
)
parser.add_argument(
    "--metadata",
    type=str,
    required=True,
    help="Path to the metadata CSV file. The CSV should have a 'cell_id' column that matches the index of the AnnData object.",
)
parser.add_argument(
    "--MIT",
    action="store_true",
    help="If set, the input data is from MIT. This will use 'individualID' for matching.",
)
args = parser.parse_args()


def add_metadata(adata: ad.AnnData, metadata: pd.DataFrame, mit: bool) -> None:

    # if ROSMAP (non-MIT), we have only the projid. This is also in the metadata csv.
    # for MIT, we have to use the individualID in the h5ad, which should also be present in the metadata csv.

    join_col = "individualID" if mit else "projid"
    if join_col not in adata.obs.columns:
        raise ValueError(f"Column '{join_col}' not found in AnnData object, which is required for matching with metadata.")
    
    for col in metadata.columns:
        adata.obs[col] = np.nan
        

        # TODO: correctly join the dataframes here




if __name__ == "__main__":
    # Read the AnnData object
    adata = ad.read_h5ad(args.path)

    # Read the metadata CSV file
    metadata = pd.read_csv(args.metadata)

    # Add metadata to the AnnData object
    add_metadata(adata, metadata, mit=args.MIT)

    # Write the modified AnnData object back to the file
    print("Writing modified AnnData object...")
    adata.write_h5ad(args.path, compression='gzip')
    print(f"Metadata added to {args.path}.")