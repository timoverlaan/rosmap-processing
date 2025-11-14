# Load required library
library(SeuratDisk)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(args) == 0) {
  stop("Usage: Rscript h5Seurat_to_h5ad.R file1.h5Seurat file2.h5Seurat ...")
}

## Iterate over each provided .h5Seurat file and convert using a local temp dir
for (h5Seurat_path in args) {
  base_name <- sub("\\.h5Seurat$", "", h5Seurat_path)
  h5ad_path <- paste0(base_name, ".h5ad")

  cat("Converting:", h5Seurat_path, "→", h5ad_path, "\n")

  # create a temporary directory next to the source file to avoid /tmp
  src_dir <- dirname(h5Seurat_path)
  local_temp_dir <- file.path(src_dir, ".temp_conversion")
  dir.create(local_temp_dir, showWarnings = FALSE, recursive = TRUE)

  # Make sure temp files are cleaned up on error or exit for this file
  temp_h5seurat <- tempfile(pattern = "conv_", tmpdir = local_temp_dir, fileext = ".h5Seurat")
  temp_h5ad <- sub("\\.h5Seurat$", ".h5ad", temp_h5seurat)
  cleanup_file <- function() {
    if (file.exists(temp_h5seurat)) try(file.remove(temp_h5seurat), silent = TRUE)
    if (file.exists(temp_h5ad)) try(file.remove(temp_h5ad), silent = TRUE)
  }

  tryCatch({
    # Load the .h5Seurat file (read-only)
    seurat_obj <- LoadH5Seurat(h5Seurat_path)

    # Convert all factor columns in metadata to character to preserve labels
    seurat_obj@meta.data[] <- lapply(seurat_obj@meta.data, function(x) {
      if (is.factor(x)) as.character(x) else x
    })

    # Save to a local temp h5Seurat and convert from that temp file
    SaveH5Seurat(seurat_obj, filename = temp_h5seurat, overwrite = TRUE)
    Convert(temp_h5seurat, dest = "h5ad", assay = "RNA", overwrite = TRUE)

    # Move the result into place (overwrite if present)
    if (!file.rename(temp_h5ad, h5ad_path)) {
      stop(sprintf("Failed to move converted file %s -> %s", temp_h5ad, h5ad_path))
    }

    cat("✓ Successfully converted:", h5ad_path, "\n")
  }, error = function(e) {
    cat("✗ Failed to convert:", h5Seurat_path, "\n")
    cat("  Error:", conditionMessage(e), "\n")
    cleanup_file()
    # rethrow so calling shells can catch the error if needed
    stop(e)
  }, finally = {
    # attempt best-effort cleanup of temporary files for this file
    cleanup_file()
  })

  # free up memory
  gc()
}

cat("All conversions attempted.\n")
