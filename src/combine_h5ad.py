import anndata as ad
import pandas as pd
import numpy as np

import argparse

from tqdm import tqdm


parser = argparse.ArgumentParser(description="Combine all provided h5ad files into a single h5ad file.")
parser.add_argument(
    "paths",
    type=str,
    nargs="+",
    required=True,
    help="Paths to the input h5ad files.",
)
parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="Path to the output h5ad file. This will contain all the data from the input files.",
)
args = parser.parse_args()


def combine_h5ad(paths: list[str]) -> ad.AnnData:
    """
    Combines multiple h5ad files into a single AnnData object.
    """
    adatas = [ad.read_h5ad(path) for path in tqdm(paths, desc="Loading h5ad files")]
    
    # Concatenate all AnnData objects
    print("Combining AnnData objects...")
    # TODO: check if the joining is correct here.
    combined_adata = ad.concat(adatas, join='outer', label='celltype_object', keys=paths, pairwise=True, merge="unique")
    print(f"Combined AnnData object shape: {combined_adata.shape}")

    return combined_adata


if __name__ == "__main__":
    print("Starting h5ad combining pipeline...")
    print("Arguments:")
    print(f"  Paths:")
    for path in args.paths:
        print(f"    - {path}")
    print(f"  Output: {args.output}")
    print("\n")

    combined_adata = combine_h5ad(args.paths)

    print("Saving combined AnnData object to:", args.output)
    combined_adata.write_h5ad(args.output, compression='gzip')
    print("Done.")