"""Constants used throughout the ROSMAP processing pipeline."""

from typing import NamedTuple

# Processing thresholds
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 5
# Donors with fewer cells are dropped. The per-donor seurat_v3 HVG fit (loess)
# fails on the smallest, shallowest ROSMAP-MIT donors (172, 219 and 298 cells).
MIN_CELLS_PER_DONOR = 300
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

# Synapse IDs for PsychAD (MSSM_PsychAD in the AD Knowledge Portal, Lee et al.
# 2024): the h5ad after the authors' QC, and the donor metadata
SYNAPSE_IDS_PSYCHAD = {
    "h5ad": "syn68559919",
    "h5ad_raw_no_qc": "syn68554605",  # Before QC (106 GB, only with --include-raw-h5ad)
    "individual_metadata": "syn55251012",
    "biospecimen_metadata": "syn55251032",
    "neuropathology_diagnosis": "syn74913889",
    "neuropathology_dictionary": "syn57373364",
    "hashing_metadata": "syn53446426",
    "assay_metadata": "syn57341886",
    "donor_id_mapping": "syn62808264",
    "deidentified_donors": "syn68498949",
}

# Synapse IDs for AMP-AD_DiverseCohorts: the gene expression part of the 10x
# multiome (filtered counts per cell class, 2023, and the gene expression
# matrices per subclass, 2025), the cell metadata of both releases, and the
# donor metadata. The libraries pool several donors; only the cell metadata
# says which donor a cell is from. The ATAC files are not used.
SYNAPSE_IDS_DIVERSE_COHORTS = {
    "counts_astro": "syn52138701",
    "counts_endo": "syn52138703",
    "counts_epend": "syn52138711",
    "counts_gaba": "syn52138732",
    "counts_glut": "syn52138744",
    "counts_micro": "syn52138745",
    "counts_olig": "syn55446335",
    "counts_opc": "syn52138751",
    "counts_other_immune": "syn52138752",
    "gex_astro": "syn71722418",
    "gex_endo": "syn71722427",
    "gex_epend": "syn71722437",
    "gex_gad1": "syn71722501",
    "gex_gad2": "syn71722506",
    "gex_glut1": "syn71722509",
    "gex_glut2": "syn71722518",
    "gex_glut3": "syn71722534",
    "gex_glut4": "syn71722579",
    "gex_micro": "syn71722580",
    "gex_olig1": "syn71722585",
    "gex_olig2": "syn71722588",
    "gex_opc": "syn71722591",
    "cell_metadata_2023": "syn52136317",
    "cell_clusters_2024": "syn72367414",
    "individual_metadata": "syn51757646",
    "biospecimen_metadata": "syn51757645",
    "assay_multiome_metadata": "syn53649030",
    "rush_projid_mapping": "syn55128291",
    "extended_neuropath": "syn61089081",
    "data_dictionary": "syn51757664",
    "data_dictionary_neuropath": "syn61089558",
    "outcomes": "syn51757663",
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
    # PsychAD, Lee et al. 2024 (pinned 2026-11-03)
    "syn68559919": SynapseFile("merged_final_clean.h5ad", 2, "8e45971c37cc65f04c5b3f73770d16f5", 113394136209),
    "syn68554605": SynapseFile("merged_raw_no_qc.h5ad", 2, "2b8609f60568c54639e9f3a9ce6eaf86", 106118643495),
    "syn55251012": SynapseFile("NPS-AD_individual_metadata.csv", 9, "6f8b814195a3eed3818ed7a394b8fb9c", 464473),
    "syn55251032": SynapseFile("NPS-AD_biospecimen_metadata.csv", 4, "a0770eb197260e5216dc180a8ef00ad0", 304014),
    "syn74913889": SynapseFile("NPS-AD_neuropathology_and_diagnosis_metadata.csv", 3, "8ed554cb2d5882c8ffd1cc94cefbf78e", 147581),
    "syn57373364": SynapseFile("NPS-AD_neuropathology_diagnosis_data_dictionary.pdf", 2, "0c78f4b66c0d7d30242cb5c594381d38", 56842),
    "syn53446426": SynapseFile("NPS-AD_nucleus_hashing_metadata.csv", 6, "a4b2d769f1b02fd2250ae2375477fd3c", 826279),
    "syn57341886": SynapseFile("NPS-AD_snRNAseq_assay_metadata.csv", 4, "9855fd26900b5124fdb4f03ea16df7a1", 324059),
    "syn62808264": SynapseFile("NPS-AD_DonorID_mapping.csv", 3, "8db47327c70fe0cd5a04d867790270ee", 28892),
    "syn68498949": SynapseFile("de-identify_AllDonors_simple.csv", 4, "b296e9e86d1bc6fc4a0e94935f1500e5", 46246),
    # AMP-AD_DiverseCohorts (pinned 2026-11-03)
    "syn52138701": SynapseFile("counts_Astro_filtered.rds", 2, "1540fb5bc49c2fcff1b58fa5d25f6e5a", 1136644278),
    "syn52138703": SynapseFile("counts_Endo_filtered.rds", 2, "b808b94573407a4d7d6f3659a32daf73", 144873133),
    "syn52138711": SynapseFile("counts_Epend_filtered.rds", 2, "1023427aef218e675108374cf9330b96", 61379706),
    "syn52138732": SynapseFile("counts_Gaba_filtered.rds", 2, "977cb2c94ebe15a024300089609d6a8a", 2985174589),
    "syn52138744": SynapseFile("counts_Glut_filtered.rds", 2, "91d8451eab67884f8c16796b7122ed8d", 4500865811),
    "syn52138745": SynapseFile("counts_Micro_filtered.rds", 2, "c3a75277adc80c940a93c4f663d7f2c8", 361854107),
    "syn55446335": SynapseFile("counts_Olig_filtered.rds", 2, "9a156f1b86e97d350cbee61e2f069ca2", 2279492674),
    "syn52138751": SynapseFile("counts_Opc_filtered.rds", 2, "4531dd38960894856698b0109cc25bd7", 370599733),
    "syn52138752": SynapseFile("counts_OtherImmune_filtered.rds", 2, "cb240a5662cf287b32a95222d417dd04", 10823490),
    "syn71722418": SynapseFile("gex_mat_Astro.rda", 2, "efc6b246d926b501b5faf58b8e2293c8", 4092103846),
    "syn71722427": SynapseFile("gex_mat_Endo.rda", 2, "6220a052be91dae98fd2ad2372e2f6d2", 993699602),
    "syn71722437": SynapseFile("gex_mat_Epend.rda", 2, "0e7b6308346e6bce4a0dad9ca0ad6a4e", 433205871),
    "syn71722501": SynapseFile("gex_mat_Gad1.rda", 2, "c747c5a7b219526ba2b156abc39832dc", 5060497765),
    "syn71722506": SynapseFile("gex_mat_Gad2.rda", 2, "ac51d052a428799d3131f7feb5d9d629", 4316680072),
    "syn71722509": SynapseFile("gex_mat_Glut1.rda", 2, "663444bbb6e482665c421e31fe694c35", 4951859862),
    "syn71722518": SynapseFile("gex_mat_Glut2.rda", 2, "3daee380d939c381d02e2918f8b47cc9", 4357006421),
    "syn71722534": SynapseFile("gex_mat_Glut3.rda", 2, "4865e44cae8db578d93a4129d9d746aa", 3958081346),
    "syn71722579": SynapseFile("gex_mat_Glut4.rda", 2, "cdb66186adb7f76816b268a1270e0865", 4198494011),
    "syn71722580": SynapseFile("gex_mat_Micro.rda", 2, "13025a47a7fa8b2165fe6d1d75933ecf", 1536175278),
    "syn71722585": SynapseFile("gex_mat_Olig1.rda", 2, "2981f8e5a15d295974421877c25a8ad2", 3239210039),
    "syn71722588": SynapseFile("gex_mat_Olig2.rda", 2, "abadc0ee5607dd0fee68078221251ef9", 2911968239),
    "syn71722591": SynapseFile("gex_mat_Opc.rda", 2, "9632a427ebba2f5920c038646dc44ddb", 1643163991),
    "syn52136317": SynapseFile("rosmap_diversity_snmultiome_celltypes_metadata_20230410.csv", 2, "cd7f3defbd2d55ccead51f3ab2d91b4d", 351486305),
    "syn72367414": SynapseFile("rosmap_diversity_snmultiome_clusterdata_20241107.csv", 2, "72c8c7a8f285cee96dc2b03b374a4138", 735715383),
    "syn51757646": SynapseFile("AMP-AD_DiverseCohorts_individual_metadata.csv", 22, "82e8473a617dc84312201e26be0855c3", 359693),
    "syn51757645": SynapseFile("AMP-AD_DiverseCohorts_biospecimen_metadata.csv", 15, "8d5c3cf33864c24bdc38586e63d67a08", 651374),
    "syn53649030": SynapseFile("AMP-AD_DiverseCohorts_assay_multiome_metadata.csv", 6, "5022b2010dcba1735a80775aafd61aa8", 78640),
    "syn55128291": SynapseFile("AMP-AD_DiverseCohorts_Rush_projid_mapping.csv", 3, "b7c464d80757433e3bb72bdc151b51d6", 6482),
    "syn61089081": SynapseFile("Harmonized_Extended_Neuropath_Diverse_Cohorts.csv", 3, "8544d084b5d02e34edc6fb8cb25cf776", 91462),
    "syn51757664": SynapseFile("DiverseCohorts_DataDictionary.pdf", 2, "31c795101132552349f0eff21880f3d0", 75905),
    "syn61089558": SynapseFile("DiverseCohorts_DataDictionary_ExtNeuropath.pdf", 3, "aa54632b12e0f5371ce85ffa9354e2eb", 99611),
    "syn51757663": SynapseFile("DiverseCohortsOutcomes.pdf", 2, "bca81fd587465007f143bce5be369579", 52960),
}

# Large files that are only downloaded on request.
SYNAPSE_OPTIONAL = {"syn52392369", "syn68554605"}

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
    "SYNAPSE_IDS_PSYCHAD",
    "SYNAPSE_IDS_DIVERSE_COHORTS",
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
