"""
Heatmap of the overlaps between the gene sets of all exported h5ad files,
labelled by version, HVG selection and dataset.

Input is the output of scripts/extract_gene_sets.py (gene tables only). Genes
are matched by Ensembl ID; files that store only symbols (the April 2026 and
older outputs) are mapped with their dataset's v2 gene map
(data/v2/reference/{rosmap_mit,seaad}_gene_map.tsv).

Usage:
    python scripts/hvg_overlap_heatmap.py gene_sets.tsv \\
        --mit-gene-map rosmap_mit_gene_map.tsv --seaad-gene-map seaad_gene_map.tsv \\
        --out out/hvg_overlap/hvg_overlap_all

Writes <out>.png, <out>.pdf, <out>_files.tsv (one row per file) and
<out>_pairs.tsv (every ordered pair: intersection, Jaccard, fraction of a).
"""

import argparse
import json
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
SURFACE, DIAGONAL = "#fcfcfb", "#f0efec"
BLUES = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DATASET_COLORS = {"ROSMAP-MIT": "#eb6834", "SeaAD": "#1baf7a", "Green 2024": "#6250d6"}
DATASET_ORDER = list(DATASET_COLORS)

VERSIONS = ["v2 (Oct 2026)", "v1 (Apr 2026)", "Jun 2025"]
VERSION_TAG = {"v2 (Oct 2026)": "v2", "v1 (Apr 2026)": "v1", "Jun 2025": "2025"}

# (pattern on the file name without dataset prefix, selection label), in plot order
SELECTIONS = [
    (r"^hvg(\d+)_k30$", "seurat_v3 on counts"),
    (r"^hvg(\d+)_postnorm_k30$", "seurat on log-normalized"),
    (r"^hvg(2000)_3way_k30$", "seurat_v3, genes all 3 datasets have"),
    (r"^lieke_allgenes_hvg(2)k_k30$", "seurat_v3 within Lieke's all-genes"),
    (r"^lieke_allgenes_hvg(2)k_postnorm_k30$", "seurat within Lieke's all-genes"),
    (r"^lieke_hvg_k30$", "Lieke's HVG list (fixed)"),
    (r"^mathys2019degs_k30$", "Mathys 2019 DEGs (fixed)"),
    (r"^(1000)_k30$", "1000 HVGs, method not recorded"),
]
SELECTION_ORDER = [label for _, label in SELECTIONS]


def describe(path: str) -> dict:
    """Version, dataset, HVG selection and where the genes were chosen, from the path."""
    stem = Path(path).stem
    if path.startswith("data/v2/processed/"):
        version = VERSIONS[0]
    elif stem == "seaad1000_k30":
        version = VERSIONS[2]
    else:
        version = VERSIONS[1]

    for prefix, dataset in (("rosmap_mit_", "ROSMAP-MIT"), ("rosmap_green_", "Green 2024"),
                            ("seaad_", "SeaAD"), ("seaad", "SeaAD")):
        if stem.startswith(prefix):
            rest = stem[len(prefix):]
            break
    else:
        raise ValueError(f"unknown dataset: {path}")

    for pattern, label in SELECTIONS:
        m = re.match(pattern, rest)
        if m:
            selection = label
            n_target = int(m.group(1)) * (1000 if rest.startswith("lieke_allgenes") else 1) if m.groups() else 0
            break
    else:
        raise ValueError(f"unknown selection: {path}")

    # SeaAD files, and Green's 3-way file, copy the genes of the ROSMAP-MIT file
    if version == VERSIONS[2]:
        chosen_on = "unknown"
    elif dataset == "SeaAD" or (dataset == "Green 2024" and "3way" in rest):
        chosen_on = "ROSMAP-MIT"
    else:
        chosen_on = dataset

    label = f"{VERSION_TAG[version]} · {stem}"
    if chosen_on not in (dataset, "unknown"):
        label += "  ← MIT genes"
    return dict(path=path, stem=stem, version=version, dataset=dataset, selection=selection,
                n_target=n_target, chosen_on=chosen_on, label=label)


def load(gene_sets_tsv: Path, gene_maps: dict):
    files, genes = [], {}
    with open(gene_sets_tsv, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "FILE":
                provenance = json.loads(parts[6]) if len(parts) > 6 and parts[6] else {}
                files.append((parts[1], provenance))
                genes[parts[1]] = []
            elif parts[0] == "GENE":
                genes[parts[1]].append((parts[2], parts[3]))

    symbol_to_id = {
        dataset: pd.read_csv(path, sep="\t").set_index("var_name")["gene_id"].to_dict()
        for dataset, path in gene_maps.items()
    }
    rows, sets = [], []
    for path, provenance in files:
        info = describe(path)
        keys, unmapped = set(), 0
        for name, gene_id in genes[path]:
            # Older files have no ID column, or symbols in it
            if gene_id.startswith("ENSG"):
                gene_id = gene_id.split(".")[0]
            else:
                gene_id = symbol_to_id.get(info["dataset"], {}).get(name, "")
            if isinstance(gene_id, str) and gene_id.startswith("ENSG"):
                keys.add(gene_id)
            else:
                keys.add(f"symbol:{name}")  # kept, so set sizes stay true; matches nothing
                unmapped += 1
        params = provenance.get("params", {})
        info.update(n_genes=len(keys), unmapped=unmapped,
                    hvg_batch_key=params.get("hvg_batch_key", ""),
                    git_commit=(provenance.get("git_commit") or "")[:7])
        rows.append(info)
        sets.append(keys)

    table = pd.DataFrame(rows)
    order = sorted(range(len(rows)), key=lambda i: (
        VERSIONS.index(rows[i]["version"]), SELECTION_ORDER.index(rows[i]["selection"]),
        rows[i]["n_target"], DATASET_ORDER.index(rows[i]["dataset"])))
    return table.iloc[order].reset_index(drop=True), [sets[i] for i in order]


def overlaps(sets):
    n = len(sets)
    inter = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i, n):
            inter[i, j] = inter[j, i] = len(sets[i] & sets[j])
    return inter


def groups(values):
    """(start, end, value) runs of equal consecutive values."""
    runs, start = [], 0
    for k in range(1, len(values) + 1):
        if k == len(values) or values[k] != values[start]:
            runs.append((start, k - 1, values[start]))
            start = k
    return runs


def plot(table, inter, out_base: Path):
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    n = len(table)
    sizes = np.diag(inter)
    pct = 100 * inter / sizes[:, None]
    identical = (inter == sizes[:, None]) & (inter == sizes[None, :])

    cell = 0.30
    M = n * cell
    L_br, L_lab, S, G, R = 2.0, 2.9, 0.14, 0.03, 1.6
    T_title, T_br, T_lab, B = 0.9, 0.45, 2.9, 2.4
    W = L_br + L_lab + S + G + M + R
    H = T_title + T_br + T_lab + S + G + M + B
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)

    def axes(x, y, w, h):
        ax = fig.add_axes([x / W, y / H, w / W, h / H])
        ax.set_facecolor(SURFACE)
        return ax

    mx, my = L_br + L_lab + S + G, B

    # Heatmap
    ax = axes(mx, my, M, M)
    cmap = LinearSegmentedColormap.from_list("blues", BLUES)
    cmap.set_bad(DIAGONAL)
    shown = np.where(np.eye(n, dtype=bool), np.nan, pct)
    ax.imshow(np.ma.masked_invalid(shown), cmap=cmap, vmin=0, vmax=100, aspect="auto",
              interpolation="nearest")
    for k in range(n + 1):
        ax.axhline(k - 0.5, color=SURFACE, lw=0.6)
        ax.axvline(k - 0.5, color=SURFACE, lw=0.6)
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, f"{sizes[i]:,}", ha="center", va="center", fontsize=4.6,
                        color=INK2, style="italic")
                continue
            p = pct[i, j]
            text = "=" if identical[i, j] else str(int(np.floor(p)))
            ax.text(j, i, text, ha="center", va="center", fontsize=5.6 if text != "=" else 7,
                    color="white" if p >= 55 else INK, fontweight="bold" if text == "=" else None)

    versions, selections = table["version"].tolist(), table["selection"].tolist()
    keys = list(zip(versions, selections))
    for start, end, _ in groups(keys)[:-1]:
        ax.axhline(end + 0.5, color=INK2, lw=0.7)
        ax.axvline(end + 0.5, color=INK2, lw=0.7)
    for start, end, _ in groups(versions)[:-1]:
        ax.axhline(end + 0.5, color=INK, lw=2.2)
        ax.axvline(end + 0.5, color=INK, lw=2.2)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1)

    # Dataset strips, carrying the row and column labels
    colors = np.array([matplotlib.colors.to_rgb(DATASET_COLORS[d]) for d in table["dataset"]])
    left = axes(mx - G - S, my, S, M)
    left.imshow(colors[:, None, :], aspect="auto", interpolation="nearest")
    left.set_xticks([])
    left.set_yticks(range(n))
    left.set_yticklabels(table["label"], fontsize=6.6, color=INK)
    top = axes(mx, my + M + G, M, S)
    top.imshow(colors[None, :, :], aspect="auto", interpolation="nearest")
    top.set_yticks([])
    top.set_xticks(range(n))
    top.xaxis.tick_top()
    top.set_xticklabels(table["label"], fontsize=6.6, color=INK, rotation=90)
    for strip in (left, top):
        strip.tick_params(length=0, pad=3)
        for spine in strip.spines.values():
            spine.set_visible(False)

    # Version and selection brackets (left), version brackets (top)
    br = axes(0.15, my, L_br - 0.15, M)
    br.set_xlim(0, 1)
    br.set_ylim(n - 0.5, -0.5)
    br.axis("off")
    for start, end, version in groups(versions):
        br.plot([0.09, 0.09], [start - 0.4, end + 0.4], color=INK, lw=2.4, solid_capstyle="butt")
        br.text(0.04, (start + end) / 2, version, rotation=90, ha="center", va="center",
                fontsize=9.5, fontweight="bold", color=INK)
    for start, end, (_, selection) in groups(keys):
        br.plot([0.17, 0.17], [start - 0.38, end + 0.38], color=INK2, lw=1.4, solid_capstyle="butt")
        br.text(0.2, (start + end) / 2, selection, ha="left", va="center", fontsize=6.9, color=INK2)

    tb = axes(mx, my + M + G + S + T_lab, M, T_br)
    tb.set_xlim(-0.5, n - 0.5)
    tb.set_ylim(0, 1)
    tb.axis("off")
    for start, end, version in groups(versions):
        tb.plot([start - 0.4, end + 0.4], [0.25, 0.25], color=INK, lw=2.4, solid_capstyle="butt")
        tb.text((start + end) / 2, 0.45, version, ha="center", va="bottom", fontsize=9.5,
                fontweight="bold", color=INK)

    # Colour scale
    cax = axes(mx + M + 0.35, my + M * 0.55, 0.16, M * 0.4)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 100)), cax=cax)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, colors=INK2, length=0)
    cb.set_label("% of the row file's genes\nthat are also in the column file", fontsize=8.5, color=INK2)

    # Title, legend and notes
    fig.text(0.15 / W, 1 - 0.35 / H, f"Gene-set overlap between all {n} exported h5ad files",
             fontsize=17, color=INK, va="top")
    fig.text(0.15 / W, 1 - 0.75 / H,
             "Cell: % of the row file's genes that are also in the column file (floored; genes "
             "matched by Ensembl ID).  '=': identical gene sets.  Diagonal: number of genes.  "
             "Thick lines separate versions, thin lines HVG selections.",
             fontsize=9, color=INK2, va="top")
    fig.legend(handles=[Patch(color=c, label=d) for d, c in DATASET_COLORS.items()],
               loc="upper left", bbox_to_anchor=(mx / W, (B - 0.25) / H), ncol=3, frameon=False,
               fontsize=8.5, title="Dataset (strip colour; also in the file name)",
               title_fontsize=8.5, alignment="left", labelcolor=INK2)
    older = table[table["version"] != VERSIONS[0]]
    notes = [
        "v2 (Oct 2026, commit {commit}): HVGs chosen per donor (batch key '{key}'). Pairs: chosen on "
        "ROSMAP-MIT among the genes SeaAD also has; the SeaAD file copies them (← MIT genes). 3-way: "
        "chosen on ROSMAP-MIT among the genes all three datasets have; SeaAD and Green copy them. "
        "Green standalone: Green's own genes.".format(
            commit=",".join(sorted(set(table.loc[table["version"] == VERSIONS[0], "git_commit"]))),
            key=",".join(sorted(set(table.loc[table["version"] == VERSIONS[0], "hvg_batch_key"])))),
        "v1 (Apr 2026, superseded): HVGs chosen on all cells pooled; SeaAD files copy the ROSMAP-MIT "
        "genes matched by symbol, so they have fewer. v1 and Jun 2025 files store symbols only; they "
        f"were mapped to Ensembl IDs with their dataset's v2 gene map ({int(older['unmapped'].sum())} of "
        f"{int(older['n_genes'].sum()):,} symbols unmapped; those match nothing).",
        "Jun 2025: data/seaAD/PFC/RNAseq/seaad1000_k30.h5ad, an older SeaAD output; how its genes were "
        "chosen is not recorded in the file.",
    ]
    y = B - 1.0
    for note in notes:
        lines = textwrap.wrap(note, width=int((W - mx - 0.3) * 15))
        fig.text(mx / W, y / H, "\n".join(lines), fontsize=8, color=INK2, va="top", linespacing=1.4)
        y -= 0.17 * len(lines) + 0.12

    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=200, facecolor=SURFACE)
    fig.savefig(out_base.with_suffix(".pdf"), facecolor=SURFACE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gene_sets", type=Path)
    parser.add_argument("--mit-gene-map", type=Path, required=True)
    parser.add_argument("--seaad-gene-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    table, sets = load(args.gene_sets, {"ROSMAP-MIT": args.mit_gene_map, "SeaAD": args.seaad_gene_map})
    inter = overlaps(sets)
    plot(table, inter, args.out)

    sizes = np.diag(inter)
    table.drop(columns="label").to_csv(f"{args.out}_files.tsv", sep="\t", index=False)
    pairs = [
        dict(a=table.at[i, "label"], b=table.at[j, "label"], n_a=sizes[i], n_b=sizes[j],
             intersection=inter[i, j], jaccard=round(inter[i, j] / (sizes[i] + sizes[j] - inter[i, j]), 4),
             frac_of_a=round(inter[i, j] / sizes[i], 4))
        for i in range(len(table)) for j in range(len(table)) if i != j
    ]
    pd.DataFrame(pairs).to_csv(f"{args.out}_pairs.tsv", sep="\t", index=False)
    print(table[["version", "stem", "selection", "n_genes", "unmapped", "chosen_on"]].to_string())


if __name__ == "__main__":
    main()
