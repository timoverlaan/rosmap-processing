import anndata as ad
import argparse


parser = argparse.ArgumentParser(description="Check a h5ad file.")
parser.add_argument(
    "path",
    type=str,
    help="Path to the h5ad file to check.",
)
args = parser.parse_args()


if __name__ == "__main__":

    print("Importing h5ad file: ", args.path)

    adata = ad.read_h5ad(args.path)

    print(adata)
    print(adata.var_names)
    print(adata.obs_names)
    print()

    print("   =========================    \n")

    for obs_col in adata.obs.columns:
        print(obs_col)
        print(adata.obs[obs_col].value_counts())
        print()

    print("   =========================    \n")

    # and preview the data
    print("X:")
    print(adata.X.shape)
    print(adata.X[:20, :20].toarray())

    if adata.raw is not None:
        print("raw.X:")
        print(adata.raw.X.shape)
        print(adata.raw.X[:20, :20].toarray())

    