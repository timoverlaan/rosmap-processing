"""
Print the gene set of each h5ad file: reads only `var` (the gene table) and the
HVG / provenance entries of `uns`. Never reads obs, X or layers, so the output
holds no cell or donor data.

Run on the cluster, in the container, and save the output locally:

    ssh daic 'cd <repo> && apptainer exec --writable-tmpfs --containall --pwd /opt/app \\
        --bind data/:/opt/app/data/ ./container_pixi_0-2-1.sif \\
        pixi run python - $(ls data/v2/processed/*.h5ad)' \\
        < scripts/extract_gene_sets.py > gene_sets.tsv

Output, tab-separated:
    FILE  path  n_vars  var_columns  uns_keys  hvg_json  provenance_json
    GENE  path  var_name  gene_id     (gene_id empty if var has no ID column)
"""

import json
import sys

import h5py
import numpy as np


def _decode(values):
    return [v.decode() if isinstance(v, bytes) else str(v) for v in values]


def read_column(group, name):
    obj = group[name]
    if isinstance(obj, h5py.Group):  # categorical
        categories = _decode(obj["categories"][()])
        return [categories[c] if c >= 0 else "" for c in obj["codes"][()]]
    return _decode(obj[()])


def read_tree(obj):
    if isinstance(obj, h5py.Dataset):
        value = obj[()]
        if isinstance(value, bytes):
            return value.decode()
        if isinstance(value, np.ndarray):
            return _decode(value.ravel()[:20])
        return value.item() if hasattr(value, "item") else value
    return {key: read_tree(obj[key]) for key in obj.keys()}


def main(paths):
    for path in paths:
        with h5py.File(path, "r") as f:
            var = f["var"]
            index_name = var.attrs["_index"]
            names = read_column(var, index_name)
            columns = [c for c in var.keys() if c != index_name]
            id_column = next((c for c in ("gene_id", "gene_ids") if c in columns), None)
            gene_ids = read_column(var, id_column) if id_column else [""] * len(names)

            uns_keys = list(f["uns"].keys()) if "uns" in f else []
            hvg = read_tree(f["uns/hvg"]) if "hvg" in uns_keys else {}
            provenance = {}
            if "provenance_scanpy" in uns_keys:
                record = read_tree(f["uns/provenance_scanpy"])
                provenance = {key: record.get(key) for key in ("params", "git_commit", "created")}

            print("\t".join([
                "FILE", path, str(len(names)), ",".join(columns), ",".join(uns_keys),
                json.dumps(hvg, default=str), json.dumps(provenance, default=str),
            ]))
            for name, gene_id in zip(names, gene_ids):
                print(f"GENE\t{path}\t{name}\t{gene_id}")


if __name__ == "__main__":
    main(sys.argv[1:])
