"""
Compare gene list overlaps across multiple .txt files (one gene per line).

Usage:
    python scripts/compare_gene_lists.py <file1> <file2> [file3 ...]

Each file is labelled by its path. Add more files (e.g. outputs of
slurm/extract_gene_names.sh) simply by appending them as arguments.

Output: a pairwise overlap table printed to stdout.
"""

import sys
from pathlib import Path


def load_genes(path: str) -> set:
    return {
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def short_label(path: str) -> str:
    return Path(path).name


def main(files: list[str]) -> None:
    gene_sets = {f: load_genes(f) for f in files}
    labels = [short_label(f) for f in files]
    sizes = {f: len(gene_sets[f]) for f in files}

    col_w = max(len(l) for l in labels) + 2

    # Header
    header = f"{'':>{col_w}}" + "".join(f"{l:>{col_w}}" for l in labels)
    print(header)
    print("-" * len(header))

    for f_a in files:
        row = f"{short_label(f_a):>{col_w}}"
        for f_b in files:
            if f_a == f_b:
                row += f"{'—':>{col_w}}"
            else:
                overlap = gene_sets[f_a] & gene_sets[f_b]
                n = len(overlap)
                pct_a = 100 * n / sizes[f_a] if sizes[f_a] else 0
                cell = f"{n} ({pct_a:.0f}%)"
                row += f"{cell:>{col_w}}"
        print(row)

    print()
    print("Cell format: overlap_count (% of row set)")
    print()
    print("Set sizes:")
    for f in files:
        print(f"  {short_label(f)}: {sizes[f]}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python scripts/compare_gene_lists.py <file1> <file2> [file3 ...]")
        sys.exit(1)
    main(sys.argv[1:])
