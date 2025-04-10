library(Seurat)
library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: [pixi run] Rscript rds_to_h5seurat.R file1.rds file2.rds ...")
}

for (rds_path in args) {
  base_name <- sub("\\.rds$", "", rds_path)
  h5seurat_path <- paste0(base_name, ".h5Seurat")

  cat("Converting:", rds_path, "→", h5seurat_path, "\n")

  seurat_obj <- readRDS(rds_path)

  # first we have to update the Seurat object to the latest version
  seurat_obj <- UpdateSeuratObject(seurat_obj)

  SaveH5Seurat(seurat_obj, filename = h5seurat_path, overwrite = TRUE)

  # free up memory
  rm(seurat_obj)
  gc()
}
