"""Constants used throughout the ROSMAP processing pipeline."""

# Processing thresholds
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 5
DEFAULT_HVG_COUNT = 2000
DEFAULT_K_NEIGHBORS = 30
DEFAULT_PCA_COMPONENTS = 50

# Normalization
CPM_TARGET = 1e6  # Normalize to 1 million reads per cell (CPM)

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
    "raw_data": "syn52392369",  # Unfiltered counts
    "sc_metadata": "syn52368902",
    "specimen_metadata": "syn52430345",
    "individual_metadata": "syn52430346",
}

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
