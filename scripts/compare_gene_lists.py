"""
Compare gene list overlaps across multiple .txt files (one gene per line).

Usage:
    python scripts/compare_gene_lists.py <file1> <file2> [file3 ...]
    python scripts/compare_gene_lists.py --heatmap overlap.png <file1> <file2> [file3 ...]

Each file is labelled by its path. Add more files (e.g. outputs of
slurm/extract_gene_names.sh) simply by appending them as arguments.

Output: a pairwise overlap table printed to stdout, and optionally a heatmap image.
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
    name = Path(path).stem  # drop extension
    for suffix in ("_gene_names", "_genes"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def build_matrix(files: list[str], gene_sets: dict, sizes: dict):
    """Returns (pct_matrix, count_matrix) as lists-of-lists, row = set A, col = set B."""
    n = len(files)
    pct = [[0.0] * n for _ in range(n)]
    counts = [[0] * n for _ in range(n)]
    for i, f_a in enumerate(files):
        for j, f_b in enumerate(files):
            if i == j:
                pct[i][j] = float("nan")
                counts[i][j] = sizes[f_a]
            else:
                overlap = len(gene_sets[f_a] & gene_sets[f_b])
                counts[i][j] = overlap
                pct[i][j] = 100 * overlap / sizes[f_a] if sizes[f_a] else 0.0
    return pct, counts


def print_table(files: list[str], pct: list, counts: list) -> None:
    labels = [short_label(f) for f in files]
    sizes = [counts[i][i] for i in range(len(files))]
    col_w = max(len(l) for l in labels) + 2

    header = f"{'':>{col_w}}" + "".join(f"{l:>{col_w}}" for l in labels)
    print(header)
    print("-" * len(header))

    for i, f_a in enumerate(files):
        row = f"{short_label(f_a):>{col_w}}"
        for j in range(len(files)):
            if i == j:
                row += f"{'—':>{col_w}}"
            else:
                cell = f"{counts[i][j]} ({pct[i][j]:.0f}%)"
                row += f"{cell:>{col_w}}"
        print(row)

    print()
    print("Cell format: overlap_count (% of row set)")
    print()
    print("Set sizes:")
    for i, f in enumerate(files):
        print(f"  {short_label(f)}: {sizes[i]}")


def save_heatmap(files: list[str], pct: list, counts: list, out_path: str) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    labels = [short_label(f) for f in files]
    sizes = [counts[i][i] for i in range(len(files))]
    n = len(files)

    pct_arr = np.array(pct, dtype=float)
    count_arr = np.array(counts, dtype=int)

    fig, ax = plt.subplots(figsize=(max(6, n * 1.4), max(5, n * 1.2)))

    # Mask diagonal for colour scaling
    masked = np.ma.masked_invalid(pct_arr)
    im = ax.imshow(masked, cmap="YlOrRd", vmin=0, vmax=100, aspect="auto")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("% of row set", fontsize=10)

    # Annotate cells
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, f"{sizes[i]:,}", ha="center", va="center",
                        fontsize=9, color="grey", style="italic")
            else:
                p = pct_arr[i, j]
                c = count_arr[i, j]
                color = "white" if p > 60 else "black"
                ax.text(j, i, f"{c:,}\n({p:.0f}%)", ha="center", va="center",
                        fontsize=8.5, color=color, linespacing=1.4)

    # Axis labels with set sizes
    row_labels = [f"{l}\n(n={s:,})" for l, s in zip(labels, sizes)]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=10)
    ax.set_yticklabels(row_labels, fontsize=10)

    ax.set_xlabel("column set", fontsize=10)
    ax.set_ylabel("row set  (% = overlap / row set)", fontsize=10)
    ax.set_title("Pairwise gene list overlap", fontsize=13, pad=12)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Heatmap saved to {out_path}")


def main(args: list[str]) -> None:
    heatmap_out = None
    files = []

    i = 0
    while i < len(args):
        if args[i] == "--heatmap":
            heatmap_out = args[i + 1]
            i += 2
        else:
            files.append(args[i])
            i += 1

    if len(files) < 2:
        print("Usage: python scripts/compare_gene_lists.py [--heatmap out.png] <file1> <file2> [file3 ...]")
        sys.exit(1)

    gene_sets = {f: load_genes(f) for f in files}
    sizes = {f: len(gene_sets[f]) for f in files}
    pct, counts = build_matrix(files, gene_sets, sizes)

    print_table(files, pct, counts)

    if heatmap_out:
        save_heatmap(files, pct, counts, heatmap_out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python scripts/compare_gene_lists.py [--heatmap out.png] <file1> <file2> [file3 ...]")
        sys.exit(1)
    main(sys.argv[1:])
