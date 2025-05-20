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

  # Load the .h5Seurat file
  seurat_obj <- LoadH5Seurat(h5Seurat_path)

  # Convert all factor columns in metadata to character to preserve labels
  # Without this, the categories are saved as integers in the .h5ad file, 
  #    and we can't find the celltype labels!
  seurat_obj@meta.data[] <- lapply(seurat_obj@meta.data, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # Save changes back to .h5Seurat (overwrite is necessary)
  SaveH5Seurat(seurat_obj, filename = h5Seurat_path, overwrite = TRUE)


  # Perform the conversion
  Convert(h5Seurat_path, dest = "h5ad", assay = "RNA", overwrite = TRUE)

  # free up memory
  gc()
}
