"""
Compare gene list overlaps across multiple .txt or .h5ad files (one gene per line / var_names).

Usage:
    # Simple (no groups):
    python scripts/compare_gene_lists.py [--heatmap out.png] file1 file2 [file3 ...]

    # With groups (files appear in group order; ungrouped files appended last):
    python scripts/compare_gene_lists.py --heatmap out.png \\
        --group "Lieke ref"   data/Lieke/hvg.txt data/Lieke/allgenes.txt \\
        --group "Lieke HVGs"  data/processed/rosmap_mit_lieke_hvg_k30_gene_names.txt \\
        --group "Regular 2k"  data/processed/rosmap_mit_hvg2000_k30.h5ad

Cells show: overlap count and % of the *row* set that is covered.
Subset annotations in each cell:
    ⊆  row set is a strict subset of col set   (100% of row is in col, sizes differ)
    ⊇  col set is a strict subset of row set
    =  sets are identical
"""

import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _infer_top_n(stem: str) -> int | None:
    """
    Parse the intended top-N from a _hvg_names.txt filename, e.g.:
      rosmap_mit_hvg1000_postnorm_k30  -> 1000
      rosmap_mit_lieke_allgenes_hvg2k  -> 2000
    Returns None if no recognisable pattern found.
    """
    import re
    m = re.search(r'hvg(\d+k?)', stem)
    if not m:
        return None
    raw = m.group(1)
    if raw.endswith('k'):
        return int(raw[:-1]) * 1000
    return int(raw)


def _load_top_n_from_scores(scores_path: Path, n: int) -> set:
    """Read a _hvg_scores.tsv and return the top-N gene names by rank/dispersion."""
    import pandas as pd
    df = pd.read_csv(scores_path, sep='\t', index_col=0)
    if 'highly_variable_rank' in df.columns:
        df = df.sort_values('highly_variable_rank', ascending=True)
    elif 'dispersions_norm' in df.columns:
        df = df.sort_values('dispersions_norm', ascending=False)
    else:
        raise ValueError(f"Cannot determine sort column in {scores_path}")
    return set(df.index[:n])


def load_genes(path: str) -> set:
    p = Path(path)
    if p.suffix == ".h5ad":
        import anndata as ad
        adata = ad.read_h5ad(p, backed="r")
        return set(adata.var_names)

    # For _hvg_names.txt files: check if we need to slice to top-N
    if p.stem.endswith("_hvg_names"):
        base_stem = p.stem[: -len("_hvg_names")]
        n = _infer_top_n(base_stem)
        if n is not None:
            scores_path = p.with_name(base_stem + "_hvg_scores.tsv")
            if scores_path.exists():
                all_genes = [
                    line.strip()
                    for line in p.read_text().splitlines()
                    if line.strip() and not line.startswith("#")
                ]
                if len(all_genes) > n:
                    return _load_top_n_from_scores(scores_path, n)

    return {
        line.strip()
        for line in p.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def short_label(path: str) -> str:
    name = Path(path).stem
    for suffix in ("_gene_names", "_hvg_names", "_genes"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    # Clean up common trailing tokens for readability
    for token in ("_k30", "_k15"):
        if name.endswith(token):
            name = name[: -len(token)]
    return name


# ---------------------------------------------------------------------------
# Matrix
# ---------------------------------------------------------------------------

def build_matrix(files, gene_sets, sizes):
    """Returns (pct_matrix, count_matrix) as lists-of-lists, row=A, col=B."""
    n = len(files)
    pct    = [[0.0] * n for _ in range(n)]
    counts = [[0]   * n for _ in range(n)]
    for i, f_a in enumerate(files):
        for j, f_b in enumerate(files):
            if i == j:
                pct[i][j]    = float("nan")
                counts[i][j] = sizes[f_a]
            else:
                overlap      = len(gene_sets[f_a] & gene_sets[f_b])
                counts[i][j] = overlap
                pct[i][j]    = 100 * overlap / sizes[f_a] if sizes[f_a] else 0.0
    return pct, counts


def subset_symbol(i, j, counts):
    """Return subset/superset/equal indicator string for cell (i,j)."""
    if i == j:
        return ""
    size_i = counts[i][i]
    size_j = counts[j][j]
    overlap = counts[i][j]
    row_in_col = (overlap == size_i)  # row set is subset of col set
    col_in_row = (overlap == size_j)  # col set is subset of row set
    if row_in_col and col_in_row:
        return "(=)"
    if row_in_col:
        return "(row<col)"
    if col_in_row:
        return "(row>col)"
    return ""


def subset_symbol_unicode(i, j, counts):
    """Unicode version for heatmap rendering (matplotlib handles it fine)."""
    if i == j:
        return ""
    size_i = counts[i][i]
    size_j = counts[j][j]
    overlap = counts[i][j]
    row_in_col = (overlap == size_i)
    col_in_row = (overlap == size_j)
    if row_in_col and col_in_row:
        return "="
    if row_in_col:
        return "\u2286"   # ⊆
    if col_in_row:
        return "\u2287"   # ⊇
    return ""


# ---------------------------------------------------------------------------
# Text output
# ---------------------------------------------------------------------------

def print_table(files, pct, counts):
    labels = [short_label(f) for f in files]
    sizes  = [counts[i][i] for i in range(len(files))]
    col_w  = max(len(l) for l in labels) + 2

    header = f"{'':>{col_w}}" + "".join(f"{l:>{col_w}}" for l in labels)
    print(header)
    print("-" * len(header))

    for i, f_a in enumerate(files):
        row = f"{short_label(f_a):>{col_w}}"
        for j in range(len(files)):
            if i == j:
                row += f"{'—':>{col_w}}"
            else:
                sym  = subset_symbol(i, j, counts)
                cell = f"{counts[i][j]} ({pct[i][j]:.0f}%){sym}"
                row += f"{cell:>{col_w}}"
        print(row)

    print()
    print("Cell format: overlap_count (% of row set) [(row<col) = row subset of col, (row>col) = superset, (=) = identical]")
    print()
    print("Set sizes:")
    for i, f in enumerate(files):
        print(f"  {short_label(f)}: {sizes[i]:,}")


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------

# Palette for group colour bars — cycles if more than len(PALETTE) groups
PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52",
    "#8172B2", "#937860", "#DA8BC3", "#8C8C8C",
]


def save_heatmap(files, pct, counts, out_path, groups=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import to_rgba

    labels = [short_label(f) for f in files]
    sizes  = [counts[i][i] for i in range(len(files))]
    n      = len(files)

    pct_arr   = np.array(pct,    dtype=float)
    count_arr = np.array(counts, dtype=int)

    # --- figure layout ---
    cell_size  = 1.5          # inches per cell
    label_pad  = 2.0          # left margin for y-labels
    group_bar  = 0.25 if groups else 0.0  # width of coloured group bar
    fig_w = label_pad + group_bar + n * cell_size + 1.5   # +1.5 for colorbar
    fig_h = group_bar + n * cell_size + 1.5               # +1.5 for x-labels

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # --- main heatmap ---
    masked = np.ma.masked_invalid(pct_arr)
    im = ax.imshow(masked, cmap="YlOrRd", vmin=0, vmax=100, aspect="auto")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("% of row set", fontsize=10)

    # --- cell annotations ---
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, f"{sizes[i]:,}", ha="center", va="center",
                        fontsize=9, color="grey", style="italic")
            else:
                p   = pct_arr[i, j]
                c   = count_arr[i, j]
                sym = subset_symbol_unicode(i, j, counts)
                fg  = "white" if p > 60 else "black"
                txt = f"{c:,}\n({p:.0f}%)"
                if sym:
                    txt += f"\n{sym}"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=8, color=fg, linespacing=1.3)

    # --- group separators (thick lines between groups) ---
    if groups:
        boundary = 0
        for g_idx, (g_label, g_files) in enumerate(groups):
            boundary += len(g_files)
            if boundary < n:
                ax.axhline(boundary - 0.5, color="black", linewidth=2.5)
                ax.axvline(boundary - 0.5, color="black", linewidth=2.5)

    # --- axis tick labels with set sizes ---
    row_labels = [f"{l}\n(n={s:,})" for l, s in zip(labels, sizes)]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_yticklabels(row_labels, fontsize=9)

    # --- coloured group bars on top and left ---
    if groups:
        ax_pos = ax.get_position()   # in figure fraction — computed later via tight_layout trick

        # We draw group spans as coloured rectangles just outside the axes,
        # using ax.transData + ax.get_xlim()/get_ylim() for positioning.
        # Easier: use secondary broken-axis trick with inset_axes, OR just
        # draw patches in data coordinates slightly outside the heatmap.
        pos = 0
        legend_handles = []
        for g_idx, (g_label, g_files) in enumerate(groups):
            sz    = len(g_files)
            color = PALETTE[g_idx % len(PALETTE)]
            start = pos - 0.5
            end   = pos + sz - 0.5

            # top bar (above column axis)
            ax.annotate(
                "", xy=(end, -0.5), xytext=(start, -0.5),
                xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-", color=color, lw=4),
                annotation_clip=False,
            )
            ax.text(
                (start + end) / 2, -0.52, g_label,
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8.5,
                color=color, fontweight="bold",
                clip_on=False,
            )

            # left bar (beside row axis)
            ax.annotate(
                "", xy=(-0.5, start), xytext=(-0.5, end),
                xycoords=("axes fraction", "data"),
                textcoords=("axes fraction", "data"),
                arrowprops=dict(arrowstyle="-", color=color, lw=4),
                annotation_clip=False,
            )
            ax.text(
                -0.52, (start + end) / 2, g_label,
                transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8.5,
                color=color, fontweight="bold", rotation=90,
                clip_on=False,
            )

            legend_handles.append(mpatches.Patch(color=color, label=g_label))
            pos += sz

        ax.legend(handles=legend_handles, loc="upper left",
                  bbox_to_anchor=(1.18, 1), fontsize=9, title="Groups",
                  title_fontsize=9, framealpha=0.8)

    ax.set_xlabel("column set", fontsize=10)
    ax.set_ylabel("row set  (% = overlap / row set)", fontsize=10)
    ax.set_title("Pairwise gene list overlap", fontsize=13, pad=14)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Heatmap saved to {out_path}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv):
    """
    Returns (files, groups, heatmap_out).
    files  : ordered list of all file paths (grouped first, then ungrouped)
    groups : list of (label, [files]) — only set when --group is used
    """
    heatmap_out = None
    groups      = []   # list of (label, [files])
    ungrouped   = []

    i = 0
    while i < len(argv):
        if argv[i] == "--heatmap":
            heatmap_out = argv[i + 1]
            i += 2
        elif argv[i] == "--group":
            label      = argv[i + 1]
            i += 2
            group_files = []
            while i < len(argv) and not argv[i].startswith("--"):
                group_files.append(argv[i])
                i += 1
            groups.append((label, group_files))
        else:
            ungrouped.append(argv[i])
            i += 1

    # Build ordered file list
    files = []
    for _, gf in groups:
        files.extend(gf)
    if ungrouped:
        if groups:
            groups.append(("Other", ungrouped))
        files.extend(ungrouped)

    return files, (groups if groups else None), heatmap_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args):
    files, groups, heatmap_out = parse_args(args)

    if len(files) < 2:
        print(__doc__)
        sys.exit(1)

    gene_sets = {f: load_genes(f) for f in files}
    sizes     = {f: len(gene_sets[f]) for f in files}
    pct, counts = build_matrix(files, gene_sets, sizes)

    print_table(files, pct, counts)

    if heatmap_out:
        save_heatmap(files, pct, counts, heatmap_out, groups=groups)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
