library(Seurat)
library(SeuratDisk)

# Usage: [pixi run] Rscript rds_to_h5Seurat.R [--overwrite] file1.rds file2.rds ...
# Existing .h5Seurat outputs are skipped unless --overwrite is given.

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
files <- args[args != "--overwrite"]

if (length(files) == 0) {
  stop("Usage: [pixi run] Rscript rds_to_h5seurat.R [--overwrite] file1.rds file2.rds ...")
}

for (rds_path in files) {
  base_name <- sub("\\.rds$", "", rds_path)
  h5seurat_path <- paste0(base_name, ".h5Seurat")

  if (file.exists(h5seurat_path) && !overwrite) {
    cat("Skipping (output exists):", h5seurat_path, "\n")
    next
  }

  cat("Converting:", rds_path, "→", h5seurat_path, "\n")

  seurat_obj <- readRDS(rds_path)

  # first we have to update the Seurat object to the latest version
  seurat_obj <- UpdateSeuratObject(seurat_obj)

  SaveH5Seurat(seurat_obj, filename = h5seurat_path, overwrite = TRUE)

  # free up memory
  rm(seurat_obj)
  gc()
}
