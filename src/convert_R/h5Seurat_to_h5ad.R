# Load required library
library(SeuratDisk)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(args) == 0) {
  stop("Usage: Rscript h5Seurat_to_h5ad.R file1.h5Seurat file2.h5Seurat ...")
}

# Iterate over each provided .h5Seurat file
for (h5Seurat_path in args) {
  # Derive the base name without the .h5Seurat extension
  base_name <- sub("\\.h5Seurat$", "", h5Seurat_path)
  # Define the output .h5ad file path
  h5ad_path <- paste0(base_name, ".h5ad")

  # Display conversion message
  cat("Converting:", h5Seurat_path, "→", h5ad_path, "\n")

  # Perform the conversion
  Convert(h5Seurat_path, dest = "h5ad", overwrite = TRUE)

  # free up memory
  gc()
}
