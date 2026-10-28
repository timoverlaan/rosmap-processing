"""Constants used throughout the ROSMAP processing pipeline."""

from typing import NamedTuple

# Processing thresholds
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 5
DEFAULT_HVG_COUNT = 2000
DEFAULT_K_NEIGHBORS = 30
DEFAULT_PCA_COMPONENTS = 50

# Normalization: counts per cell after scaling, before log1p. The per-cell
# totals are computed over all genes, before any gene subsetting.
NORMALIZE_TARGET_SUM = 1e4

# Column names - standardized across ROSMAP and SeaAD
CLASS_NAME = "Class"
SUBCLASS_NAME = "Celltype"
SUPERTYPE_NAME = "Subtype"
DONOR_ID_COLUMN = "Donor ID"

# Synapse IDs for ROSMAP data
SYNAPSE_IDS_ROSMAP = {
    "astrocytes": "syn53693925",
    "cux2_plus": "syn53694054",
    "cux2_minus": "syn53694068",
    "inhibitory": "syn53693978",
    "microglia": "syn53693904",
    "oligodendroglia": "syn53693961",
    "vascular_niche": "syn53693877",
    "cell_annotations": "syn55219673",
    "cell_annotations_n437": "syn53694215",
    "clinical": "syn3191087",
    "clinical_codebook": "syn3191090",
    "scrna_metadata": "syn21073536",
    "specimen_metadata": "syn21323366",
    "demultiplex": "syn34572333",
}

# Synapse IDs for ROSMAP MIT data
SYNAPSE_IDS_ROSMAP_MIT = {
    "astrocytes": "syn52368912",
    "excitatory_neurons_set1": "syn52368925",
    "excitatory_neurons_set2": "syn52368950",
    "excitatory_neurons_set3": "syn52368932",
    "immune_cells": "syn52368905",
    "inhibitory_neurons": "syn52368921",
    "opcs": "syn52368910",
    "oligodendrocytes": "syn52368918",
    "raw_data": "syn52392369",  # Unfiltered counts (103 GB, only with --include-raw-h5ad)
    "sc_metadata": "syn52368902",
    "specimen_metadata": "syn52430345",
    "individual_metadata": "syn52430346",
}


class SynapseFile(NamedTuple):
    """A pinned Synapse file: local file name, version, MD5 and size in bytes."""

    name: str
    version: int
    md5: str
    size: int


# Every downloaded file is pinned to a version and checked against its MD5.
# Checked on 2026-10-28 against the public Synapse metadata: the re-versioning
# in July 2026 (all of these got a new version on 2026-07-08/09/14) did not
# change any file; each new version has the previous version's MD5. Earlier
# local copies of six Green h5Seurat files differ from every version because
# a 2025 script version rewrote them in place (see README, "Data history").
SYNAPSE_FILES = {
    # ROSMAP, Green et al. 2024
    "syn53693925": SynapseFile("astrocytes.h5Seurat", 2, "2dff03136c38395d2e89efdaff648389", 12485424164),
    "syn53694054": SynapseFile("cux2+.h5Seurat", 2, "5411682bf6b40286819fa6300c51ad56", 32711871130),
    "syn53694068": SynapseFile("cux2-.h5Seurat", 2, "ebb757ad0e5a2536d8d6466a7c50793e", 23307977746),
    "syn53693978": SynapseFile("inhibitory.h5Seurat", 2, "b9e7a39013c8ca09db75d74f2896cce8", 17522992653),
    "syn53693904": SynapseFile("microglia.h5Seurat", 2, "9399ef699a50937ac1ec9efd9a1d5d71", 4037731221),
    "syn53693961": SynapseFile("oligodendroglia.h5Seurat", 4, "ffe0ba4ab150eac4fe1bf65d2d784951", 19559590497),
    "syn53693877": SynapseFile("vascular.niche.h5Seurat", 4, "6de8aa49c7a903064bf541394f3fb1e1", 699748427),
    "syn55219673": SynapseFile("cell-annotation.full-atlas.csv", 2, "738cf4bf18c59c9a3e15d3bcad5b3e97", 207369610),
    "syn53694215": SynapseFile("cell-annotation.n437.csv", 4, "cfc7274a003c500ed0d25257c5bb8268", 213417913),
    "syn3191087": SynapseFile("ROSMAP_clinical.csv", 12, "c91cc4985f7bb4908f356c674c57ae52", 335210),
    "syn3191090": SynapseFile("ROSMAP_clinical_codebook.pdf", 3, "8962216bde62ce5dd9e148cb74b9600a", 36942),
    "syn21073536": SynapseFile("ROSMAP_assay_scrnaSeq_metadata.csv", 10, "be8c68d5b05d101122608f0fee022a3d", 326675),
    "syn21323366": SynapseFile("ROSMAP_biospecimen_metadata.csv", 20, "8004dce267cc8163a74952c296e02703", 1203320),
    "syn34572333": SynapseFile("ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv", 2, "47fcfd3124d20d3d0ad6add92360f50b", 967621293),
    # ROSMAP-MIT, Mathys et al. 2023
    "syn52368912": SynapseFile("Astrocytes.rds", 3, "1b973527aecf0a3d0b8548f723633d56", 1702950062),
    "syn52368925": SynapseFile("Excitatory_neurons_set1.rds", 3, "3518280f854c43b3cb5488973335e3a4", 6066311071),
    "syn52368950": SynapseFile("Excitatory_neurons_set2.rds", 3, "5f9a1acd501fdeeb4ca34bca5d6c5a09", 10476855962),
    "syn52368932": SynapseFile("Excitatory_neurons_set3.rds", 3, "413fbfa5c302369c9efe5380ec3e27b4", 6226421585),
    "syn52368905": SynapseFile("Immune_cells.rds", 3, "ad08f5bf317f49ee7bcc040c3c318374", 645762382),
    "syn52368921": SynapseFile("Inhibitory_neurons.rds", 3, "bfab3e416935772982fecac816f6ffb0", 5594328392),
    "syn52368910": SynapseFile("OPCs.rds", 3, "31506333ef0e17a600e03f8aecda9a57", 1222819608),
    "syn52368918": SynapseFile("Oligodendrocytes.rds", 3, "eeacbce6e85b6cb1bbc60a1fb2ea2e0d", 5420754177),
    "syn52392369": SynapseFile("PFC427_raw_data.h5ad", 2, "f2525f6a001c42c15d232fed53308269", 103066397330),
    # The Synapse entity is named without ".csv"; the downloaded file has it.
    "syn52368902": SynapseFile("MIT_ROSMAP_Multiomics_assay_snRNAseq_metadata.csv", 5, "a029c558c796b47c6c52e94ae2bd74a4", 532555),
    "syn52430345": SynapseFile("MIT_ROSMAP_Multiomics_biospecimen_metadata.csv", 2, "7fbb7a0e29f639b3f91ba7aed1295c81", 541180),
    "syn52430346": SynapseFile("MIT_ROSMAP_Multiomics_individual_metadata.csv", 3, "7a825d858470ea50400d1596f6e54087", 27192),
}

# Large files that are only downloaded on request.
SYNAPSE_OPTIONAL = {"syn52392369"}

# File extensions
VALID_H5AD_EXTENSIONS = [".h5ad"]
VALID_H5SEURAT_EXTENSIONS = [".h5Seurat", ".h5seurat"]
VALID_RDS_EXTENSIONS = [".rds", ".RDS"]
VALID_TXT_EXTENSIONS = [".txt"]

# SeaAD AWS S3 paths
SEAAD_S3_BUCKET = "sea-ad-single-cell-profiling"
SEAAD_S3_PATHS = {
    "pfc_rnaseq": "PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad",
    "pfc_metadata": "PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv",
}

# Cell type mappings
ROSMAP_MIT_CELL_CLASSES = {
    "Astrocytes": {"class": "Glia", "subclass": "Astrocytes"},
    "Excitatory_neurons_set1": {"class": "Neuron", "subclass": "Excitatory"},
    "Excitatory_neurons_set2": {"class": "Neuron", "subclass": "Excitatory"},
    "Excitatory_neurons_set3": {"class": "Neuron", "subclass": "Excitatory"},
    "Immune_cells": {"class": "Glia", "subclass": "Immune"},
    "Inhibitory_neurons": {"class": "Neuron", "subclass": "Inhibitory"},
    "OPCs": {"class": "Glia", "subclass": "OPCs"},
    "Oligodendrocytes": {"class": "Glia", "subclass": "Oligodendrocytes"},
}

# Braak stage mapping
BRAAK_MAPPING = {
    "Braak 0": 0,
    "Braak I": 1,
    "Braak II": 2,
    "Braak III": 3,
    "Braak IV": 4,
    "Braak V": 5,
    "Braak VI": 6,
}

# CERAD score mapping (inverted for SeaAD compatibility)
CERAD_MAPPING = {
    "Absent": 4,
    "Sparse": 3,
    "Moderate": 2,
    "Frequent": 1,
    "Reference": 4,
}

# Thal phase mapping
THAL_MAPPING = {
    "Thal 0": 0,
    "Thal 1": 1,
    "Thal 2": 2,
    "Thal 3": 3,
    "Thal 4": 4,
    "Thal 5": 5,
    "Thal 6": 6,
}

# Define __all__ to control what gets exported with "from constants import *"
__all__ = [
    # Processing thresholds
    "MIN_GENES_PER_CELL",
    "MIN_CELLS_PER_GENE",
    "DEFAULT_HVG_COUNT",
    "DEFAULT_K_NEIGHBORS",
    "DEFAULT_PCA_COMPONENTS",
    # Normalization
    "NORMALIZE_TARGET_SUM",
    # Column names
    "CLASS_NAME",
    "SUBCLASS_NAME",
    "SUPERTYPE_NAME",
    "DONOR_ID_COLUMN",
    # Synapse IDs
    "SYNAPSE_IDS_ROSMAP",
    "SYNAPSE_IDS_ROSMAP_MIT",
    "SynapseFile",
    "SYNAPSE_FILES",
    "SYNAPSE_OPTIONAL",
    # File extensions
    "VALID_H5AD_EXTENSIONS",
    "VALID_H5SEURAT_EXTENSIONS",
    "VALID_RDS_EXTENSIONS",
    "VALID_TXT_EXTENSIONS",
    # SeaAD S3
    "SEAAD_S3_BUCKET",
    "SEAAD_S3_PATHS",
    # Cell type mappings
    "ROSMAP_MIT_CELL_CLASSES",
    # Pathology mappings
    "BRAAK_MAPPING",
    "CERAD_MAPPING",
    "THAL_MAPPING",
]
