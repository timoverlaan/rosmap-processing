# Load required library
library(SeuratDisk)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(args) == 0) {
  stop("Usage: Rscript h5Seurat_to_h5ad.R file1.h5Seurat file2.h5Seurat ...")
}

# Use current directory for temp files instead of /tmp (which may be full)
# This ensures we have enough space for large h5Seurat files
temp_dir <- file.path(getwd(), ".temp_conversion")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = temp_dir)

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
  
  # Save directly to h5ad using SaveH5Seurat -> Convert workflow
  # We use a temporary file to avoid corrupting the original
  temp_h5seurat <- tempfile(fileext = ".h5Seurat")
  SaveH5Seurat(seurat_obj, filename = temp_h5seurat, overwrite = TRUE)
  
  # Perform the conversion from temp file
  Convert(temp_h5seurat, dest = "h5ad", assay = "RNA", overwrite = TRUE)
  
  # Move the resulting h5ad file to the correct location
  temp_h5ad <- sub("\\.h5Seurat$", ".h5ad", temp_h5seurat)
  file.rename(temp_h5ad, h5ad_path)
  
  # Clean up temp file
  file.remove(temp_h5seurat)

  # free up memory
  gc()
}

# Clean up temp directory at the end
cat("Cleaning up temporary directory...\n")
unlink(temp_dir, recursive = TRUE)
cat("✓ Conversion complete\n")
