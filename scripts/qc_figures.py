"""Figures of the data QC, from the aggregate tables of `pipeline qc-check`.

    python scripts/qc_figures.py QC_DIR OUT_DIR

QC_DIR holds data/v2/qc/ as written on the cluster (one folder per dataset,
plus identity/ and sex-de-compare/); only the top-level TSVs are read, never
restricted/. Writes PNGs to OUT_DIR.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DATASETS = ["rosmap_mit", "rosmap_green", "seaad", "diverse_cohorts", "psychad"]
NAMES = {"rosmap_mit": "ROSMAP-MIT", "rosmap_green": "Green 2024", "seaad": "SEA-AD",
         "diverse_cohorts": "DiverseCohorts", "psychad": "PsychAD"}
# Categorical slots 1-5 of the reference palette, fixed per dataset
COLOURS = {"rosmap_mit": "#2a78d6", "rosmap_green": "#eb6834", "seaad": "#1baf7a",
           "diverse_cohorts": "#eda100", "psychad": "#e87ba4"}
CELLTYPES = ["Exc", "Inh", "Ast", "Oli", "OPC", "Mic", "Vasc", "Other"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 9, "axes.titlesize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False,
})


def read(qc: Path, dataset: str, name: str) -> pd.DataFrame:
    path = qc / dataset / f"{name}.tsv"
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def legend(fig, datasets, **kw):
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOURS[d]) for d in datasets]
    fig.legend(handles, [NAMES[d] for d in datasets], loc="upper center", ncol=len(datasets),
               bbox_to_anchor=(0.5, 1.0), handlelength=1.2, columnspacing=1.4, **kw)


def qc_distributions(qc: Path, out: Path) -> None:
    """Box per dataset and cell type: 5-95% whiskers, 25-75% box, median."""
    tables = {d: read(qc, d, "qc_cells") for d in DATASETS}
    tables = {d: t[t["level"] == "celltype"] for d, t in tables.items() if len(t)}
    metrics = [("umis", "UMIs per nucleus", True), ("genes", "Genes per nucleus", True),
               ("mito_pct", "Mitochondrial UMIs (%)", False)]
    fig, axes = plt.subplots(len(metrics), 1, figsize=(9, 8.2), sharex=True)
    width = 0.8 / len(tables)
    for ax, (metric, title, log) in zip(axes, metrics):
        for k, (d, t) in enumerate(tables.items()):
            for i, ct in enumerate(CELLTYPES):
                row = t[t["group"] == ct]
                if not len(row) or row["n_cells"].iloc[0] < 100:
                    continue
                r = row.iloc[0]
                x = i - 0.4 + width * (k + 0.5)
                ax.plot([x, x], [r[f"{metric}_p05"], r[f"{metric}_p95"]], color=COLOURS[d], lw=1)
                ax.add_patch(plt.Rectangle((x - width * 0.38, r[f"{metric}_p25"]), width * 0.76,
                                           r[f"{metric}_p75"] - r[f"{metric}_p25"], color=COLOURS[d], lw=0))
                ax.plot([x - width * 0.38, x + width * 0.38], [r[f"{metric}_p50"]] * 2, color=SURFACE, lw=1.4)
        if log:
            ax.set_yscale("log")
        ax.set_ylabel(title)
        ax.grid(axis="x", visible=False)
    axes[-1].set_xticks(range(len(CELLTYPES)), CELLTYPES)
    axes[-1].set_xlim(-0.55, len(CELLTYPES) - 0.45)
    legend(fig, list(tables))
    fig.text(0.01, 0.005, "Boxes: 25th-75th percentile, line: median, whiskers: 5th-95th. Cell types with "
             "fewer than 100 nuclei are left out. All nuclei in the released files.", color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(out / "qc_distributions.png", dpi=200)
    plt.close(fig)


def composition(qc: Path, out: Path) -> None:
    tables = {d: read(qc, d, "composition") for d in DATASETS}
    tables = {d: t for d, t in tables.items() if len(t)}
    shades = ["#2a78d6", "#7aaee8", "#1baf7a", "#eda100", "#e87ba4", "#6250d6", "#008300", "#b5b4ad"]
    fig, ax = plt.subplots(figsize=(9, 2.9))
    for k, (d, t) in enumerate(tables.items()):
        left = 0.0
        for ct, colour in zip(CELLTYPES, shades):
            row = t[t["celltype"] == ct]
            pct = float(row["pct_cells"].iloc[0]) if len(row) else 0.0
            if pct <= 0:
                continue
            ax.barh(k, pct, left=left, color=colour, edgecolor=SURFACE, linewidth=1.5, height=0.7)
            if pct >= 4:
                ax.text(left + pct / 2, k, f"{ct}\n{pct:.0f}%", ha="center", va="center", fontsize=7,
                        color="#ffffff" if ct in ("Exc", "Ast", "Oli", "OPC", "Mic") else INK)
            left += pct
    ax.set_yticks(range(len(tables)), [NAMES[d] for d in tables])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of nuclei (%)")
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in shades]
    ax.legend(handles, CELLTYPES, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=8, handlelength=1.0)
    fig.tight_layout()
    fig.savefig(out / "composition.png", dpi=200)
    plt.close(fig)


def sex_check(qc: Path, out: Path) -> None:
    datasets = [d for d in DATASETS if len(read(qc, d, "sex_histograms"))]
    fig, axes = plt.subplots(2, len(datasets), figsize=(2.1 * len(datasets) + 0.6, 4.6), squeeze=False)
    for j, d in enumerate(datasets):
        h = read(qc, d, "sex_histograms")
        s = read(qc, d, "sex_summary").iloc[0]
        ax = axes[0, j]
        for sex, colour in (("F", "#e34948"), ("M", "#2a78d6"), ("missing", "#b5b4ad")):
            part = h[h["recorded"] == sex]
            if part["n_donors"].sum():
                ax.bar(part["bin_left"], part["n_donors"], width=0.5, align="edge", color=colour, alpha=0.85,
                       label={"F": "recorded female", "M": "recorded male", "missing": "no sex"}[sex])
        ax.set_title(f"{NAMES[d]}\n{int(s['n_discordant'])} discordant, {int(s['n_ambiguous'])} ambiguous",
                     fontsize=8.5)
        ax.set_xlabel("Y genes − XIST (log2 CPM)", fontsize=8)
        if j == 0:
            ax.set_ylabel("Donors")
        # XIST has no Y copy, so XIST UMIs in a male donor's nuclei come from
        # other (female) donors: ambient RNA or doublets in pooled libraries.
        # (Y-gene UMIs in female nuclei also come from X-Y paralog mapping.)
        ax = axes[1, j]
        part = h[h["recorded"] == "xist_in_male"]
        ax.bar(100 * part["bin_left"], part["n_donors"], width=1.0, align="edge", color="#2a78d6")
        cell = read(qc, d, "sex_cell_level")
        row = cell[cell["measure"] == "share of cells with XIST UMIs"].iloc[0]
        ax.set_title(f"median {100 * row['value_p50']:.1f}% of nuclei", fontsize=8.5)
        ax.set_xlim(0, 60)
        ax.set_xlabel("Nuclei with XIST UMIs (%)", fontsize=8)
        if j == 0:
            ax.set_ylabel("Male donors")
    axes[0, 0].legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "sex_check.png", dpi=200)
    plt.close(fig)


def markers(qc: Path, out: Path) -> None:
    datasets = [d for d in DATASETS if len(read(qc, d, "markers"))]
    fig, axes = plt.subplots(len(datasets), 1, figsize=(10.5, 1.55 * len(datasets) + 0.8), squeeze=False)
    for ax, d in zip(axes[:, 0], datasets):
        m = read(qc, d, "markers")
        m = m[m["label"].isin(CELLTYPES)]
        genes = list(dict.fromkeys(m["gene"]))
        rows = [ct for ct in CELLTYPES if ct in set(m["label"])]
        for i, ct in enumerate(rows):
            part = m[m["label"] == ct].set_index("gene").reindex(genes)
            ax.scatter(range(len(genes)), [i] * len(genes), s=60 * part["frac_detected"].fillna(0),
                       c=part["mean_log_norm"].fillna(0), cmap="Blues", vmin=0, vmax=3, edgecolors="none")
        ax.set_yticks(range(len(rows)), rows, fontsize=7)
        ax.set_ylim(len(rows) - 0.4, -0.6)
        ax.set_xticks(range(len(genes)), genes if ax is axes[-1, 0] else [], rotation=90, fontsize=7)
        ax.set_xlim(-0.6, len(genes) - 0.4)
        ax.set_ylabel(NAMES[d], fontsize=8)
        ax.grid(visible=False)
    fig.text(0.01, 0.005, "Dot size: share of nuclei detecting the gene; colour: mean log-normalized "
             "expression (0-3).", color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / "markers.png", dpi=200)
    plt.close(fig)


def identity(qc: Path, out: Path) -> None:
    t = pd.read_csv(qc / "identity" / "identity_expression.tsv", sep="\t")
    pairs = list(dict.fromkeys(zip(t["dataset_a"], t["dataset_b"])))
    fig, ax = plt.subplots(figsize=(9, 3.2))
    width = 0.8 / len(pairs)
    cts = [c for c in CELLTYPES[:-1] + ["combined"] if c in set(t["celltype"])]
    shades = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
    for k, (a, b) in enumerate(pairs):
        part = t[(t["dataset_a"] == a) & (t["dataset_b"] == b)].set_index("celltype").reindex(cts)
        x = np.arange(len(cts)) - 0.4 + width * (k + 0.5)
        best = (part["pct_own_best_a_to_b"] + part["pct_own_best_b_to_a"]) / 2
        n = part["n_donors"].dropna()
        ax.bar(x, best, width=width * 0.9, color=shades[k],
               label=f"{NAMES[a]} – {NAMES[b]} ({int(n.max()) if len(n) else 0} donors)")
    ax.set_xticks(range(len(cts)), cts)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Donors whose own profile\nmatches best (%)")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.3), ncol=3, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out / "identity.png", dpi=200)
    plt.close(fig)


def de_controls(qc: Path, out: Path) -> None:
    rows = [read(qc, d, "ad_de_mathys2019") for d in DATASETS]
    m = pd.concat([r for r in rows if len(r)], ignore_index=True) if any(len(r) for r in rows) else pd.DataFrame()
    known = pd.read_csv(qc / "sex-de-compare" / "sex_de_known_genes.tsv", sep="\t") \
        if (qc / "sex-de-compare" / "sex_de_known_genes.tsv").exists() else pd.DataFrame()
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    ax = axes[0]
    if len(known):
        genes = ["XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "KDM5C", "KDM6A", "DDX3X",
                 "ZFX"]
        datasets = [d for d in DATASETS if d in set(known["dataset"])]
        width = 0.8 / len(datasets)
        for k, d in enumerate(datasets):
            part = known[known["dataset"] == d]
            med = part.groupby("gene")["rank_abs_t"].median().reindex(genes)
            x = np.arange(len(genes)) - 0.4 + width * (k + 0.5)
            ax.bar(x, med, width=width * 0.9, color=COLOURS[d], label=NAMES[d])
        ax.set_yscale("log")
        ax.set_xticks(range(len(genes)), genes, rotation=60, fontsize=7)
        ax.set_ylabel("Rank by |t| (median over cell types)")
        ax.set_title("Sex DE: known sex genes rank first", fontsize=9)
        ax.grid(axis="x", visible=False)
        ax.legend(fontsize=7, loc="upper left")
    ax = axes[1]
    if len(m):
        m = m[m["contrast"] == "pathology"]
        cts = [c for c in CELLTYPES if c in set(m["celltype"])]
        datasets = [d for d in DATASETS if d in set(m["dataset"])]
        width = 0.8 / len(datasets)
        for k, d in enumerate(datasets):
            part = m[m["dataset"] == d].set_index("celltype").reindex(cts)
            x = np.arange(len(cts)) - 0.4 + width * (k + 0.5)
            ax.bar(x, part["pct_same_direction"], width=width * 0.9, color=COLOURS[d], label=NAMES[d])
        ax.axhline(50, color=INK2, lw=1, ls="--")
        ax.text(-0.45, 47, "chance", ha="left", va="top", fontsize=7, color=INK2)
        ax.set_xticks(range(len(cts)), cts)
        ax.set_ylim(0, 100)
        ax.set_ylabel("Mathys 2019 DEGs in\nthe same direction (%)")
        ax.set_title("AD pathology vs none: Mathys 2019 DEGs", fontsize=9)
        ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(out / "de_controls.png", dpi=200)
    plt.close(fig)


def main() -> None:
    qc, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for fn in (qc_distributions, composition, sex_check, markers, identity, de_controls):
        try:
            fn(qc, out)
            print(f"{fn.__name__}: done")
        except FileNotFoundError as e:
            print(f"{fn.__name__}: skipped ({e})")


if __name__ == "__main__":
    main()
