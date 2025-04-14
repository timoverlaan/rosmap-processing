
install.packages("BiocManager")
BiocManager::install("rhdf5")  # SeuratDisk needs this

install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
