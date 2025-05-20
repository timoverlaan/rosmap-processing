import anndata as ad
import sys


if __name__ == "__main__":

    # get path from command line
    if len(sys.argv) != 2:
        print("Usage: python fix_categories.py <path_to_anndata>")
        sys.exit(1)

    path = sys.argv[1]
    print(f"Fixing categories in {path}")

    adata = ad.read_h5ad(path)

    COLUMN_NAME = "cell_type_high_resolution"
    if COLUMN_NAME in adata.obs.columns:
        if adata.obs[COLUMN_NAME].dtype != "category":
            adata.obs[COLUMN_NAME] = adata.obs[COLUMN_NAME].astype("category")
    else:
        print(f"Column {COLUMN_NAME} not found in adata.obs")  # (just ignore, and re-save, for non-MIT data)

    print(adata)

    print(adata.obs)

    print(adata.var)
    
    del adata.raw

    adata.write_h5ad(path)
