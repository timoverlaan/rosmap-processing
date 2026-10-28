# Convert .h5Seurat files to .h5ad with SeuratDisk.
#
# Usage: Rscript h5Seurat_to_h5ad.R [--outdir=DIR] [--overwrite] file1.h5Seurat file2.h5Seurat ...
#
# SeuratDisk writes the RNA `data` slot to X and the `counts` slot to .raw.
# Depending on how the object was saved, `data` is either a copy of the counts
# or log-normalized values, so the h5ad's X is NOT guaranteed to hold counts.
# Always run `rosmap_processing data prepare` on the result, which moves the
# counts from .raw into X and checks them.
#
# Without --outdir the h5ad is written next to its .h5Seurat. Existing outputs
# are skipped unless --overwrite is given, so a failed run can be resumed.

library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)
outdir <- NULL
overwrite <- FALSE
files <- character(0)
for (a in args) {
  if (startsWith(a, "--outdir=")) {
    outdir <- sub("^--outdir=", "", a)
  } else if (a == "--overwrite") {
    overwrite <- TRUE
  } else {
    files <- c(files, a)
  }
}

if (length(files) == 0) {
  stop("Usage: Rscript h5Seurat_to_h5ad.R [--outdir=DIR] [--overwrite] file1.h5Seurat ...")
}
if (!is.null(outdir)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Iterate over each provided .h5Seurat file and convert using a local temp dir
for (h5Seurat_path in files) {
  base_name <- sub("\\.h5Seurat$", "", basename(h5Seurat_path))
  target_dir <- if (is.null(outdir)) dirname(h5Seurat_path) else outdir
  h5ad_path <- file.path(target_dir, paste0(base_name, ".h5ad"))

  if (file.exists(h5ad_path) && !overwrite) {
    cat("Skipping (output exists):", h5ad_path, "\n")
    next
  }

  cat("Converting:", h5Seurat_path, "→", h5ad_path, "\n")

  # create a temporary directory next to the output to avoid /tmp
  local_temp_dir <- file.path(target_dir, ".temp_conversion")
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
  rm(seurat_obj)
  gc()
}

cat("All conversions attempted.\n")
