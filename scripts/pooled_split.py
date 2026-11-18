"""What splits one dataset's pooled pseudobulk samples into two groups
(aggregate output: counts, medians, gene names; no donor IDs).

    python scripts/pooled_split.py config/pooled_pseudobulk.yaml \
        [--dataset diverse_cohorts] [--pool-column batch] [--class-column broadcellclass]

Needs `pipeline pooled-pseudobulk` (its pseudobulk.h5ad). Per broad cell
type, the dataset's samples (log2-CPM of the genes the pipeline keeps, with
the sex means removed so that sex does not drive the split) are split in two
by k-means on their top PCs; group B is the smaller one. Then:
  separation   silhouette of the split, next to the same split of the other
               datasets
  consistency  are the same donors in group B in every cell type (adjusted
               Rand index between cell types; donors by the number of cell
               types in which they are in B). A donor's group: B if it is in
               B in most of its cell types
  donors       per group: cohort, sex, race, dementia, in other datasets, age,
               PMI, APOE4, library pool number (VR<n>) and whether pools are
               pure, nuclei and genes per nucleus, cell class shares
  samples      per cell type and group: counts per nucleus, ambient RNA
               (neuronal genes in glia, glial genes in neurons)
  replicates   after the regress correction, how often a person in group A or
               B is closest to their own sample in another dataset
  genes        the genes that differ most between the groups, per cell type
               and in most cell types
  figure       UMAP as measured and after the regress correction (as in the
               pipeline), the dataset's samples coloured by group

Outputs in <output_dir>/split_<dataset>/: report.txt (what is printed) and
umap.png.
"""

import argparse
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from rosmap_processing.pipelines import pooled_pseudobulk as pp
from rosmap_processing.pipelines.batch_effects import _matrix_node, _pca
from rosmap_processing.utils.obs_summary import read_obs_column

NEURONAL = ["SNAP25", "SYT1", "NRGN", "STMN2", "GAP43", "RBFOX3", "SLC17A7", "GAD1"]
GLIAL = ["PLP1", "MBP", "MOBP", "ST18", "GFAP", "AQP4"]
NEURONS = {"Exc", "Inh"}
N_SPLIT_PCS = 20
LINES = []


def say(*parts) -> None:
    text = " ".join(str(p) for p in parts)
    print(text, flush=True)
    LINES.append(text)


def auc(b: np.ndarray, a: np.ndarray) -> float:
    """P(a value in B > one in A), ties counted half."""
    b, a = b[np.isfinite(b)], a[np.isfinite(a)]
    if not len(a) or not len(b):
        return np.nan
    u = stats.mannwhitneyu(b, a, alternative="two-sided").statistic
    return float(u / (len(a) * len(b)))


def numeric_row(name, values_b, values_a) -> dict:
    vb, va = np.asarray(values_b, dtype=np.float64), np.asarray(values_a, dtype=np.float64)
    return dict(variable=name, median_A=np.nanmedian(va) if np.isfinite(va).any() else np.nan,
                median_B=np.nanmedian(vb) if np.isfinite(vb).any() else np.nan, auc_B_gt_A=auc(vb, va),
                n_A=int(np.isfinite(va).sum()), n_B=int(np.isfinite(vb).sum()))


def remove_sex_means(Y: np.ndarray, sex: pd.Series) -> np.ndarray:
    out = Y.copy()
    s = sex.astype(object).where(sex.notna(), "missing").to_numpy()
    for level in np.unique(s):
        rows = s == level
        out[rows] -= Y[rows].mean(axis=0) - Y.mean(axis=0)
    return out


def split_two(Y: np.ndarray, seed: int):
    """k-means (k = 2) labels on the top PCs (1 = the smaller group) and the
    silhouette."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    pcs = _pca(Y - Y.mean(axis=0), N_SPLIT_PCS, seed)
    labels = KMeans(2, n_init=10, random_state=seed).fit(pcs).labels_
    if (labels == 1).sum() > (labels == 0).sum():
        labels = 1 - labels
    return labels, float(silhouette_score(pcs, labels))


def pool_number(pool: str) -> float:
    m = re.search(r"\d+", str(pool))
    return float(m.group()) if m else np.nan


def cell_level(path: str, donor: str, pool: str, klass: str, layer):
    """Per donor of the count file: nuclei, median genes per nucleus, the
    majority pool and the number of pools, and the share of each class."""
    with h5py.File(path, "r") as f:
        donors = pd.Series(np.asarray(read_obs_column(f, donor)).astype(str))
        pools = pd.Series(np.asarray(read_obs_column(f, pool)).astype(str))
        classes = pd.Series(np.asarray(read_obs_column(f, klass)).astype(str))
        n_genes = np.diff(_matrix_node(f, layer)["indptr"][()])
    cells = pd.DataFrame({"donor": donors, "pool": pools, "klass": classes, "n_genes": n_genes})
    g = cells.groupby("donor")
    per = pd.DataFrame({"nuclei": g.size(), "genes_per_nucleus": g["n_genes"].median(),
                        "pool": g["pool"].agg(lambda s: s.value_counts().index[0]),
                        "n_pools": g["pool"].nunique()})
    shares = pd.crosstab(cells["donor"], cells["klass"], normalize="index")
    shares.columns = [f"share_{c}" for c in shares.columns]
    return per.join(shares)


def main() -> None:
    import anndata as ad

    parser = argparse.ArgumentParser()
    parser.add_argument("spec")
    parser.add_argument("--dataset", default="diverse_cohorts")
    parser.add_argument("--pool-column", default="batch")
    parser.add_argument("--class-column", default="broadcellclass")
    args = parser.parse_args()
    spec = pp.load_spec(args.spec)
    name = args.dataset
    seed = spec["seed"]
    out = Path(spec["output_dir"]) / f"split_{name}"
    out.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(Path(spec["output_dir"]) / "pseudobulk.h5ad")
    obs = adata.obs.copy()
    for c in obs.columns:
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            obs[c] = obs[c].astype(object).where(obs[c].notna(), np.nan)
    obs = obs.reset_index(drop=True)
    obs["batch"] = pp.assign_batches(obs, spec["batch_split"])
    var = adata.var.reset_index(drop=True)
    symbol_col = next((f"{d}_name" for d in ("rosmap_mit", name) if f"{d}_name" in var.columns), "gene_id")
    symbols = var[symbol_col].astype(str).to_numpy()
    counts = adata.X.tocsr()
    numeric = [c for c in spec["protect"] if c in pp.NUMERIC]

    say(f"== {name}: what splits its pooled pseudobulk samples in two (k-means on {N_SPLIT_PCS} PCs, "
        f"sex means removed; B = the smaller group)")
    labels, figure_data, sample_rows, gene_tables, replicate_rows, sil_rows = {}, {}, [], {}, [], []
    for celltype in sorted(obs["celltype"].unique()):
        sel = ((obs["celltype"] == celltype) & (obs["n_cells"] >= spec["min_cells"])).to_numpy()
        sub = obs[sel].reset_index(drop=True)
        if name not in set(sub["dataset"]):
            continue
        C = counts[np.flatnonzero(sel)].toarray().astype(np.float64)
        cpm = C / sub["library_size"].to_numpy()[:, None] * 1e6
        ds = sub["dataset"].astype(str).to_numpy()
        keep = np.all([cpm[ds == d].mean(axis=0) >= spec["min_cpm"] for d in dict.fromkeys(ds)], axis=0)
        Y = np.log2(cpm[:, keep] + 1)
        sym = symbols[keep]

        # the split, and the same split of the other datasets for reference
        for d in dict.fromkeys(ds):
            rows = ds == d
            lab, sil = split_two(remove_sex_means(Y[rows], sub.loc[rows, "sex"]), seed)
            sil_rows.append(dict(celltype=celltype, dataset=d, silhouette=sil, share_B=float(lab.mean()),
                                 n=int(rows.sum())))
            if d == name:
                own, own_lab = np.flatnonzero(rows), lab
        donors = sub.loc[own, "donor"].astype(str).to_numpy()
        labels[celltype] = pd.Series(own_lab, index=donors)

        # per sample: counts per nucleus and ambient RNA
        part = sub.iloc[own]
        lab = own_lab
        per_nucleus = (part["library_size"] / part["n_cells"]).to_numpy()
        markers = GLIAL if celltype in NEURONS else NEURONAL
        idx = [i for i, s in enumerate(symbols) if s in markers]
        ambient = np.log2(cpm[own][:, idx].mean(axis=1) + 1) if idx else np.full(len(own), np.nan)
        mito = [i for i, s in enumerate(symbols) if s.upper().startswith("MT-")]
        for row in (numeric_row("nuclei in sample", part["n_cells"].to_numpy()[lab == 1],
                                part["n_cells"].to_numpy()[lab == 0]),
                    numeric_row("counts per nucleus", per_nucleus[lab == 1], per_nucleus[lab == 0]),
                    numeric_row(f"ambient: log2 CPM of {'glial' if celltype in NEURONS else 'neuronal'} genes",
                                ambient[lab == 1], ambient[lab == 0])):
            sample_rows.append(dict(celltype=celltype, **row))
        if mito:
            share = C[own][:, mito].sum(axis=1) / part["library_size"].to_numpy()
            sample_rows.append(dict(celltype=celltype, **numeric_row("mitochondrial share", share[lab == 1],
                                                                     share[lab == 0])))

        # genes that differ between the groups (sex means removed)
        Yd = remove_sex_means(Y[own], part["sex"])
        t, _ = stats.ttest_ind(Yd[lab == 1], Yd[lab == 0], equal_var=False)
        gene_tables[celltype] = pd.DataFrame({"gene": sym, "t": t,
                                              "log2fc": Yd[lab == 1].mean(axis=0) - Yd[lab == 0].mean(axis=0)})

        # replicates after the regress correction, by group
        batch = sub["batch"].astype(str).to_numpy()
        persons = sub["person_id"].astype(str).to_numpy()
        Yc, _ = pp.correct_regress(Y, batch, sub[list(spec["protect"])], numeric)
        for g, label in ((0, "A"), (1, "B")):
            queries = set(persons[own][lab == g])
            for (a, b), v in pp.replicate_stats(Yc, ds, persons, queries=queries).items():
                if name not in (a, b):
                    continue
                rank = np.array(v["rank"])
                replicate_rows.append(dict(celltype=celltype, group=label, other=b if a == name else a,
                                           n_persons=len(rank) // 2, top1=float((rank == 0).mean())))

        # UMAP as in the pipeline
        coords = {}
        for method, M in (("none", Y), ("regress", Yc)):
            coords[method] = pp._umap(_pca(M - M.mean(axis=0), spec["n_pcs"], seed), spec)
        group = np.full(len(sub), "other", dtype=object)
        group[own] = np.where(lab == 1, "B", "A")
        figure_data[celltype] = (group, coords)

    say("\n-- separation: silhouette of the k = 2 split, per dataset (share in B)")
    sil = pd.DataFrame(sil_rows)
    say(sil.pivot_table(index="dataset", columns="celltype", values="silhouette").round(2).to_string())
    say(sil.pivot_table(index="dataset", columns="celltype", values="share_B").round(2).to_string())

    say("\n-- consistency: adjusted Rand index of the split between cell types (shared donors)")
    from sklearn.metrics import adjusted_rand_score

    cts = list(labels)
    ari = pd.DataFrame(index=cts, columns=cts, dtype=float)
    for a in cts:
        for b in cts:
            shared = labels[a].index.intersection(labels[b].index)
            ari.loc[a, b] = adjusted_rand_score(labels[a][shared], labels[b][shared]) if len(shared) > 2 else np.nan
    say(ari.round(2).to_string())
    table = pd.DataFrame(labels)
    n_b, n_ct = table.sum(axis=1), table.notna().sum(axis=1)
    say("donors by the number of cell types in which they are in B:",
        n_b.astype(int).value_counts().sort_index().to_dict())
    donor_group = pd.Series(np.where(n_b > n_ct / 2, "B", "A"), index=table.index)
    say("donor groups:", donor_group.value_counts().to_dict())

    say("\n-- donors per group")
    first = obs[obs["dataset"] == name].drop_duplicates("donor").set_index("donor")
    first.index = first.index.astype(str)
    first = first.loc[donor_group.index]
    first["group"] = donor_group
    first["in_other_datasets"] = np.where(pd.to_numeric(first["n_datasets"], errors="coerce") > 1, "yes", "no")
    for col in ("cohort", "sex", "race", "dementia", "in_other_datasets", "hispanic"):
        if col in first.columns:
            say(pd.crosstab(first[col].fillna("missing"), first["group"]).to_string(), "\n")
    cfg = spec["datasets"][name]
    cells = cell_level(cfg["counts"], cfg["donor"], args.pool_column, args.class_column, cfg.get("layer"))
    cells = cells.reindex(donor_group.index)
    cells["pool_number"] = cells["pool"].map(pool_number)
    rows = []
    for col in ("age_death", "age_capped", "pmi", "apoe4", "braak", "cerad"):
        if col in first.columns:
            v = pd.to_numeric(first[col], errors="coerce").to_numpy()
            rows.append(numeric_row(col, v[donor_group == "B"], v[donor_group == "A"]))
    for col in ["pool_number", "n_pools", "nuclei", "genes_per_nucleus"] + \
               [c for c in cells.columns if c.startswith("share_")]:
        v = cells[col].to_numpy(np.float64)
        rows.append(numeric_row(col, v[donor_group.to_numpy() == "B"], v[donor_group.to_numpy() == "A"]))
    say(pd.DataFrame(rows).round(3).to_string(index=False))

    # does a pool-number range separate the groups, and are pools pure?
    pn = cells["pool_number"].to_numpy()
    isb = (donor_group == "B").to_numpy()
    ok = np.isfinite(pn)
    if ok.sum() > 2 and 0 < isb[ok].sum() < ok.sum():
        accs = []
        for t in np.unique(pn[ok]):
            pred = pn[ok] > t
            accs.append(((pred & isb[ok]).sum() / isb[ok].sum() + (~pred & ~isb[ok]).sum() / (~isb[ok]).sum()) / 2)
        say(f"\nbest single pool-number threshold: balanced accuracy {max(max(accs), 1 - min(accs)):.2f} "
            f"(0.5 = no better than chance)")
        bins = pd.cut(pd.Series(pn[ok]), bins=10)
        say("share in B by pool-number range (tenths of the range):")
        say(pd.DataFrame({"bin": bins, "B": isb[ok]}).groupby("bin", observed=True)["B"].agg(["mean", "size"])
            .round(2).to_string())
    pools = pd.DataFrame({"pool": cells["pool"].to_numpy(), "B": isb}).groupby("pool")["B"]
    multi = pools.agg(["size", "mean"])
    multi = multi[multi["size"] >= 2]
    pure = float(((multi["mean"] == 0) | (multi["mean"] == 1)).mean()) if len(multi) else np.nan
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(200):
        perm = pd.DataFrame({"pool": cells["pool"].to_numpy(), "B": rng.permutation(isb)}).groupby("pool")["B"]
        m = perm.agg(["size", "mean"])
        m = m[m["size"] >= 2]
        null.append(((m["mean"] == 0) | (m["mean"] == 1)).mean())
    say(f"pools with 2+ donors: {len(multi)}; share with all donors in one group {pure:.2f} "
        f"(shuffled groups: {np.mean(null):.2f})")

    say("\n-- samples per cell type and group")
    say(pd.DataFrame(sample_rows).round(3).to_string(index=False))

    say("\n-- replicates after the regress correction: top-1 retrieval by group")
    rep = pd.DataFrame(replicate_rows)
    if len(rep):
        say(rep.pivot_table(index=["other", "group"], columns="celltype", values="top1").round(2).to_string())
        say(rep.pivot_table(index=["other", "group"], columns="celltype", values="n_persons").to_string())

    say("\n-- genes that differ most between the groups (B vs A, sex means removed)")
    for celltype, g in gene_tables.items():
        g = g.dropna().sort_values("t")
        up = ", ".join(f"{r.gene} {r.log2fc:+.2f}" for r in g.tail(12)[::-1].itertuples())
        down = ", ".join(f"{r.gene} {r.log2fc:+.2f}" for r in g.head(12).itertuples())
        say(f"{celltype} higher in B: {up}")
        say(f"{celltype} lower in B: {down}")
    top = {}
    for celltype, g in gene_tables.items():
        g = g.dropna()
        for direction, part in (("up", g.nlargest(100, "t")), ("down", g.nsmallest(100, "t"))):
            for gene in part["gene"]:
                top.setdefault((gene, direction), set()).add(celltype)
    common = sorted(((len(v), k) for k, v in top.items() if len(v) >= 4), reverse=True)
    say("genes in the top 100 of 4+ cell types:",
        "; ".join(f"{gene} ({direction}, {n})" for n, (gene, direction) in common[:60]))

    plot(figure_data, name, out / "umap.png", seed)
    (out / "report.txt").write_text("\n".join(LINES) + "\n")
    say(f"\nWrote {out / 'report.txt'} and {out / 'umap.png'}")


def plot(figure_data, name, path, seed) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from rosmap_processing.pipelines.umap_plots import CATEGORICAL, INK, MISSING_COLOR

    colors = {"other": MISSING_COLOR, "A": CATEGORICAL[0], "B": CATEGORICAL[1]}
    cts = list(figure_data)
    fig, axes = plt.subplots(2, len(cts), figsize=(3.0 * len(cts), 6.8), squeeze=False)
    fig.patch.set_facecolor("white")
    for c, celltype in enumerate(cts):
        group, coords = figure_data[celltype]
        order = np.random.default_rng(seed).permutation(len(group))
        order = np.concatenate([order[group[order] == "other"], order[group[order] != "other"]])
        for r, method in enumerate(("none", "regress")):
            ax = axes[r, c]
            xy = coords[method]
            ax.scatter(*xy[order].T, s=2.5, c=[colors[g] for g in group[order]], linewidths=0, rasterized=True)
            if r == 0:
                ax.set_title(f"{celltype}: A {int((group == 'A').sum())}, B {int((group == 'B').sum())}",
                             fontsize=9, color=INK)
            if c == 0:
                ax.set_ylabel("as measured" if method == "none" else "batch regressed out", fontsize=9, color=INK)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=col, markeredgewidth=0,
                      label=lab) for lab, col in ((f"{name} group A", colors["A"]), (f"{name} group B", colors["B"]),
                                                  ("other datasets", colors["other"]))]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False, labelcolor=INK)
    fig.suptitle(f"{name}: two groups of samples (k-means within the dataset), on the pooled pseudobulk UMAP",
                 fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
