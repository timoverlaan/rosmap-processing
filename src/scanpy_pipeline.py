import scanpy as sc
import anndata as ad
import argparse 

from tqdm import tqdm


parser = argparse.ArgumentParser(description="Performs the basic scanpy processing steps on a given h5ad file.")
parser.add_argument("path", type=str, help="Path to the h5ad file to process.")
parser.add_argument("output", type=str, help="Path to save the processed h5ad file.")
parser.add_argument("--layer", type=str, default=None, help="Layer to select for processing. If None, the main data will be used.")
parser.add_argument("--raw", action="store_true", help="If set, the raw data will be used instead of the main data.")
parser.add_argument("--n-genes", type=int, default=2000, help="Number of highly variable genes to select. Default is 2000.")
parser.add_argument("--k-neighbors", type=int, default=30, help="Number of neighbors for the KNN graph. Default is 30.")
parser.add_argument("--individual-pca", action="store_true", help="If set, PCA for the KNN graph will be computed for each individual separately.")
parser.add_argument("--import-genes", type=str, default=None, help="Path to a TXT or AnnData file to use the genes from. This overrides HVG selection.")
# parser.add_argument("--export-overlap", action="store_true", help="If set, and the genes from --import-genes are not a subset of the input data, \
#                     only the overlap is used, and we also export an h5ad for the refence dataset with only the overlapping genes.")
args = parser.parse_args()


if __name__ == "__main__":

    print("Starting scanpy processing pipeline...")
    print("Arguments:")
    print(f"  Path: {args.path}")
    print(f"  Output: {args.output}")
    print(f"  Layer: {args.layer}")
    print(f"  Use raw data: {args.raw}")
    print(f"  Number of highly variable genes: {args.n_genes}")
    print(f"  Number of neighbors for KNN graph: {args.k_neighbors}")
    print(f"  Individual PCA: {args.individual_pca}")
    if args.import_genes:
        print(f"  Importing genes from: {args.import_genes}")
    print("\n")

    print("Importing h5ad file: ", args.path)
    adata = ad.read_h5ad(args.path)

    print(adata)
    print()

    if args.raw:
        print("Using raw data.")
        if adata.raw is not None:
            # TODO: this is not tested yet, but it should work
            adata.X = adata.raw.X.copy()
            adata.var_names = adata.raw.var_names
            adata.obs_names = adata.raw.obs_names
        else:
            raise ValueError("Raw data is not available in the AnnData object.")

    if args.layer is not None:
        print(f"Using layer: {args.layer}")
        if args.layer in adata.layers:
            del adata.X
            adata.X = adata.layers[args.layer].copy()
            del adata.layers[args.layer]
            del adata.layers
        else:
            raise ValueError(f"Layer '{args.layer}' not found in the AnnData object.")
    
    print("Processing data...")
    
    print(f"Filtering cells and genes with minimum thresholds. Shape before filtering: {adata.shape}")
    sc.pp.filter_cells(adata, min_genes=200)
    if not args.import_genes:
        sc.pp.filter_genes(adata, min_cells=5)
    else:
        print("Not filtering genes, as we are importing genes from another file.")
    print(f"Shape after filtering: {adata.shape}")

    if args.import_genes:
        print(f"Importing genes from: {args.import_genes}")
        if args.import_genes.endswith('.txt'):
            with open(args.import_genes, 'r') as f:
                select_genes = [line.strip() for line in f if line.strip()]
            select_genes = set(select_genes)
        elif args.import_genes.endswith('.h5ad'):
            adata_import_genes = ad.read_h5ad(args.import_genes)
            select_genes = adata_import_genes.var_names

        else:
            raise ValueError("The import file must be a TXT or AnnData file.")
        
        if not set(select_genes).issubset(adata.var_names):
            print(f"Warning: The genes in the import file are not a subset of the input data. Using only the overlapping genes.")
        
        print(f"Importing {len(select_genes)} genes: ")
        print(", ".join(list(select_genes)))

        adata_overlap = adata[:, adata.var_names.isin(select_genes)].copy()
        del adata
        adata = adata_overlap

        print(f"Shape after importing genes: {adata.shape}")
        
    else:
        print(f"Selecting highly variable genes...")
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=args.n_genes, subset=True)
        print(f"Number of highly variable genes selected: {adata.var['highly_variable'].sum()}")

    print("Normalizing total (CPM) and log-transforming data...")
    sc.pp.normalize_total(adata, target_sum=1e6)  # Normalize to 1 million reads per cell
    sc.pp.log1p(adata)
    print("Done.")

    if not args.individual_pca:
        print("Calculating PCA...")
        sc.pp.pca(adata, n_comps=50)
        print("Done.")
    else:
        print("Skipping full PCA, will compute PCA for each donor separately.")

    print("Building KNN graph...")
    slices = []

    if not "Donor ID" in adata.obs.columns:
	adata.obs["Donor ID"] = adata.obs["projid"]
    for donor in tqdm(adata.obs["Donor ID"].unique(), desc="Building donor KNN graphs"):
        donor_slice = adata[adata.obs["Donor ID"] == donor].copy()
        if args.individual_pca:  # if individual PCA is requested
            sc.pp.pca(donor_slice, n_comps=50)
        sc.pp.neighbors(donor_slice, n_neighbors=args.k_neighbors, use_rep="X_pca")
        slices.append(donor_slice)

    del adata
    adata = ad.concat(
        slices,
        merge="same",
        uns_merge="same",
        label="all",
        index_unique="-",
        pairwise=True,
    )

    # Save the processed data
    print(f"\nSaving processed data to: {args.output}")
    adata.write_h5ad(args.output, compression="gzip")
    print("Processing complete.\n\n")
