import synapseclient


# For the ROSMAP data
syn_id_ast = "syn53693925"  # astrocytes.h5Seurat
syn_id_cux_2plus = "syn53694054"  # cux2+.h5Seurat
syn_id_cux_2min = "syn53694068"  # cux2-.h5Seurat
syn_id_inh = "syn53693978"  #  inhibitory.h5Seurat
syn_id_micro = "syn53693904"  # microglia.h5Seurat
syn_id_oligo = "syn53693961"  # oligodendrocytes.h5Seurat
syn_id_vasc = "syn53693877"  # vascular.niche.h5Seurat

syn_id_cell_annotations = "syn55219673"  # cell-annotation.full.atlas.csv
syn_id_cell_annocations_n437 = "syn53694215"  # cell-annotation.n437.csv -> This is the one used in the paper

# And for the ROSMAP MIT data
syn_id_mit_ast = "syn52368912"  # Astrocytes.rds
syn_id_mit_exc1 = "syn52368925"  # Excitatory_neurons_set1.rds
syn_id_mit_exc2 = "syn52368950"  # Excitatory_neurons_set2.rds
syn_id_mit_exc3 = "syn52368932"  # Excitatory_neurons_set3.rds
syn_id_mit_immune = "syn52368905"  # Immune_cells.rds
syn_id_mit_inh = "syn52368921"  # Inhibitory_neurons.rds
syn_id_mit_opc = "syn52368910"  # OPCs.rds
syn_id_mit_oligo = "syn52368918"  # Oligodendrocytes.rds

syn_id_mit_raw = "syn52392369"  # PFC427_raw_data.h5ad --> Unfiltered counts (not using this one for now)


if __name__ == "__main__":

    with open("token.txt", "r") as f:  # NOTE: Make sure to add input the token to this file when reproducing
        token = f.read().strip()
    syn = synapseclient.login(silent=True, authToken=token)


    # first download the regular rosmap data
    for syn_id in [
        syn_id_cell_annotations,
        syn_id_cell_annocations_n437,
        syn_id_ast,
        syn_id_cux_2plus,
        syn_id_cux_2min,
        syn_id_inh,
        syn_id_micro,
        syn_id_oligo,
        syn_id_vasc,
    ]:
        syn.get(syn_id_ast, downloadLocation="data/raw/ROSMAP/", ifcollision="keep.local")

    # then download the mit data
    for syn_id in [
        syn_id_mit_ast,
        syn_id_mit_exc1,
        syn_id_mit_exc2,
        syn_id_mit_exc3,
        syn_id_mit_immune,
        syn_id_mit_inh,
        syn_id_mit_opc,
        syn_id_mit_oligo,
    ]:
        syn.get(syn_id, downloadLocation="data/raw/ROSMAP_MIT/", ifcollision="keep.local")
