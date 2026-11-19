"""One row per donor of a processed h5ad, to ship with it: the donor-level
obs columns (the most common value over the donor's cells, with a count of
donors where the cells disagree), the number of cells, the pools/batches, and
the genotype check's verdict for the donor (scripts/genotype_actions.py).

    python scripts/export_donor_table.py PROCESSED.h5ad OUTPUT.tsv --dataset rosmap_green \
        --actions data/v2/qc/genotype_actions/restricted/genotype_actions.tsv \
        [--donor-column "Donor ID"] [--batch-column batch]

OUTPUT has donor IDs and clinical values (restricted). Prints counts only.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from rosmap_processing.utils.obs_summary import obs_column_names, read_obs_column

# Per-cell columns, not donor attributes
CELL_LEVEL = {"_index", "DoubletFinder.score", "nCount_RNA", "nCount_SCT", "nFeature_RNA", "nFeature_SCT", "n_genes",
              "total_counts", "cell.type", "class", "subset", "state", "source_file", "all"}
GENOTYPE_COLUMNS = ["call", "action", "reason", "identity", "relabel_to_individualID", "motivation", "own_z", "best_z",
                    "expected_z", "strength", "expression_sex", "label_sex"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--dataset", required=True, help="the dataset's name in genotype_actions.tsv")
    parser.add_argument("--actions", type=Path, required=True)
    parser.add_argument("--donor-column", default="Donor ID")
    parser.add_argument("--batch-column", default="batch")
    args = parser.parse_args()

    with h5py.File(args.h5ad, "r") as f:
        names = obs_column_names(f)
        donor = pd.Series(np.asarray(read_obs_column(f, args.donor_column)), dtype=object).astype(str)
        columns = {c: pd.Series(np.asarray(read_obs_column(f, c)), dtype=object) for c in names
                   if c not in CELL_LEVEL and c != args.donor_column}
    groups = donor.groupby(donor).groups
    rows, mixed = [], {}
    for d, idx in groups.items():
        row = {"donor": d, "n_cells": len(idx)}
        for c, values in columns.items():
            v = values.iloc[idx].dropna()
            if c == args.batch_column:
                counts = v.astype(str).value_counts()
                row["n_batches"] = len(counts)
                row["batches"] = ";".join(f"{b}:{n}" for b, n in counts.items())
                continue
            if v.empty:
                row[c] = np.nan
                continue
            counts = v.astype(str).value_counts()
            row[c] = v[v.astype(str) == counts.index[0]].iloc[0]
            if len(counts) > 1:
                mixed[c] = mixed.get(c, 0) + 1
        rows.append(row)
    table = pd.DataFrame(rows)

    actions = pd.read_csv(args.actions, sep="\t", dtype=str, keep_default_na=False)
    actions = actions[actions["dataset"] == args.dataset].drop_duplicates("donor").set_index("donor")
    per_source = [c for c in actions.columns if c.endswith("_call") or c.endswith("_reason")]
    for c in GENOTYPE_COLUMNS + ["decided_by"] + per_source:
        table[f"genotype_{c}"] = table["donor"].map(actions[c]) if c in actions else np.nan
    table["genotype_checked"] = table["donor"].isin(actions.index)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)

    print(f"{args.h5ad.name}: {len(donor)} cells, {len(table)} donors; {len(table.columns)} columns")
    print(f"donor-level columns whose cells disagree within a donor (most common value kept): {mixed or 'none'}")
    print(f"donors in the genotype check: {int(table['genotype_checked'].sum())} of {len(table)}")
    print(table.groupby(["genotype_action", "genotype_identity"], dropna=False).size().rename("n_donors")
          .reset_index().to_string(index=False))
    print(f"cells per donor: median {table['n_cells'].median():.0f}, min {table['n_cells'].min()}, "
          f"max {table['n_cells'].max()}; batches per donor: median {table['n_batches'].median():.0f}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
