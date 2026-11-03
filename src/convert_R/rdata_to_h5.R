suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
})

# Write the count matrix of an .rds or .rda file to HDF5, for `data
# from-matrix` to read in Python: the sparse arrays of the genes x cells matrix
# (its CSC arrays are the CSR arrays of cells x genes), the gene names and the
# cell names. The matrix is the file's only sparse matrix, or the RNA counts of
# its only Seurat object. Data frames with a row per cell are written next to
# it as <output>.<name>.csv.gz. Values are written as they are; `data
# from-matrix` checks that they are counts.
#
# Usage: [pixi run] Rscript rdata_to_h5.R [--overwrite] input.(rds|rda) output.h5

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- args[args != "--overwrite"]
if (length(args) != 2) stop("Usage: Rscript rdata_to_h5.R [--overwrite] input.(rds|rda) output.h5")
input <- args[1]
output <- args[2]
if (file.exists(output) && !overwrite) {
  cat("Skipping (output exists):", output, "\n")
  quit(save = "no")
}

if (grepl("\\.rds$", input, ignore.case = TRUE)) {
  objects <- list(readRDS(input))
  names(objects) <- "object"
} else {
  env <- new.env()
  loaded <- load(input, envir = env)
  objects <- mget(loaded, envir = env)
}

counts_of <- function(obj) {
  if (inherits(obj, "Seurat")) {
    m <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"),
                  error = function(e) obj@assays$RNA@counts)
    return(m)
  }
  if (inherits(obj, "sparseMatrix")) return(obj)
  NULL
}
matrices <- Filter(Negate(is.null), lapply(objects, counts_of))
if (length(matrices) != 1) {
  stop(sprintf("%s: expected one sparse matrix or Seurat object, found %d among %s",
               input, length(matrices), paste(names(objects), collapse = ", ")))
}
m <- as(matrices[[1]], "CsparseMatrix")
if (is.null(rownames(m)) || is.null(colnames(m))) stop(input, ": the matrix has no gene or cell names")
cat(sprintf("%s: %s, %d genes x %d cells, %.0f stored values\n", basename(input), names(matrices),
            nrow(m), ncol(m), length(m@x)))

tmp <- paste0(output, ".partial")
if (file.exists(tmp)) file.remove(tmp)
h5 <- H5File$new(tmp, mode = "w")
write_vector <- function(name, x) {
  h5$create_dataset(name, robj = x, chunk_dims = min(length(x), 1048576L), gzip_level = 1)
}
write_vector("data", m@x)
write_vector("indices", m@i)
write_vector("indptr", m@p)
h5[["shape"]] <- as.integer(c(ncol(m), nrow(m)))  # cells x genes
write_vector("genes", rownames(m))
write_vector("cells", colnames(m))
h5$close_all()
file.rename(tmp, output)

# Cell metadata: data frames whose rows are the matrix's cells
for (name in names(objects)) {
  obj <- objects[[name]]
  meta <- if (inherits(obj, "Seurat")) obj@meta.data else if (is.data.frame(obj)) obj else NULL
  if (!is.null(meta) && nrow(meta) == ncol(m)) {
    path <- sprintf("%s.%s.csv.gz", sub("\\.h5$", "", output), name)
    write.csv(cbind(cell = rownames(meta), meta), gzfile(path), row.names = FALSE)
    cat(sprintf("  cell metadata %s: %d columns -> %s\n", name, ncol(meta), path))
  }
}
cat("Wrote", output, "\n")
