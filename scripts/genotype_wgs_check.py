"""Check the array genotypes against WGS, and give the WGS a donor's label.

1. Maps the WGS sample IDs to individuals (specimenID or individualID in the
   biospecimen metadata) and to projids (clinical CSV), and writes the WGS
   dosages renamed by individualID, for `qc-check genotype --sample-ids match`.
2. Array vs WGS: for every pair of an array sample (ROS/MAP + projid) and a
   WGS sample, the correlation of standardized dosages over the shared eQTL
   variants (palindromic SNPs left out): about 0.9 or more for the same
   person, 0.5 for first-degree relatives, 0 for others. Per person in both:
   does the array sample match its own WGS best?
3. DNA sex from chrX heterozygosity (array and WGS) vs the recorded sex.
4. For the donors that genotype_actions.tsv relabels, drops or flags: what the
   WGS says about the label's and the matched person's genotype entries.

    python scripts/genotype_wgs_check.py --array ARRAY.tsv.gz --array-chrx ARRAY_CHRX.tsv \
        --wgs WGS.tsv.gz --wgs-chrx WGS_CHRX.tsv --eqtls EQTLS --clinical CLINICAL.csv \
        --biospecimen BIOSPECIMEN.csv [...] --actions ACTIONS.tsv --output DIR

Writes DIR/restricted/ (per-person tables, the renamed WGS dosages) and DIR/
aggregate tables; prints aggregates only.
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from rosmap_processing.pipelines.data_qc import _match_eqtls, _projid, read_genotypes

COMPLEMENT = str.maketrans("ACGT", "TGCA")
SAME, RELATIVE = 0.8, 0.35  # array-WGS correlation: the same person; a first-degree relative or closer
MALE_HET, FEMALE_HET = 0.05, 0.15


def dna_sex(rate) -> np.ndarray:
    rate = np.asarray(rate, dtype=float)
    return np.where(rate < MALE_HET, "M", np.where(rate >= FEMALE_HET, "F", np.where(np.isfinite(rate), "?", "")))


def wgs_individuals(samples, biospecimen_paths):
    """WGS sample ID -> individualID (by specimenID, else by individualID itself)."""
    spec_to_ind, individuals = {}, set()
    for path in biospecimen_paths:
        table = pd.read_csv(path, dtype=str)
        print(f"{Path(path).name}: columns {', '.join(table.columns[:12])}{' ...' if len(table.columns) > 12 else ''}")
        if {"specimenID", "individualID"} <= set(table.columns):
            for s, i in zip(table["specimenID"], table["individualID"]):
                if isinstance(s, str) and isinstance(i, str):
                    spec_to_ind.setdefault(s.strip(), i.strip())
            individuals |= set(table["individualID"].dropna().str.strip())
    out = []
    for s in samples:
        out.append(spec_to_ind.get(s) or (s if s in individuals else ""))
    out = np.asarray(out, dtype=object)
    print(f"WGS samples: {len(samples)}; mapped to an individual: {int((out != '').sum())} "
          f"(by specimenID {sum(s in spec_to_ind for s in samples)}, by individualID "
          f"{sum(s not in spec_to_ind and s in individuals for s in samples)}); distinct individuals "
          f"{len(set(out) - {''})}")
    return out


def aligned(eqtl, array_variants, wgs_variants, wgs_build=37):
    """Rows of the shared, non-palindromic eQTL variants in both, and +1/-1 for
    the WGS row's counted allele being the array's a1 or a2."""
    a_rows = _match_eqtls(eqtl, array_variants, 37)
    w_rows = _match_eqtls(eqtl, wgs_variants, wgs_build)
    keep, sign = [], []
    seen = set()
    for a, w in zip(a_rows, w_rows):
        if a < 0 or w < 0 or (a, w) in seen:
            continue
        seen.add((a, w))
        a1, a2 = array_variants.at[a, "a1"], array_variants.at[a, "a2"]
        w1 = wgs_variants.at[w, "a1"]
        if {a1, a2} in ({"A", "T"}, {"C", "G"}):
            continue
        if w1 in (a1, a1.translate(COMPLEMENT)):
            sign.append(1)
        elif w1 in (a2, a2.translate(COMPLEMENT)):
            sign.append(-1)
        else:
            continue
        keep.append((a, w))
    return np.asarray(keep, dtype=int).reshape(-1, 2), np.asarray(sign)


def standardize(X):
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = (X - np.nanmean(X, axis=1, keepdims=True)) / np.nanstd(X, axis=1, keepdims=True)
    return Z


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--array", type=Path, required=True)
    parser.add_argument("--array-chrx", type=Path)
    parser.add_argument("--wgs", type=Path, required=True)
    parser.add_argument("--wgs-chrx", type=Path)
    parser.add_argument("--wgs-build", type=int, choices=(37, 38), default=37, help="of the WGS chrN:pos IDs")
    parser.add_argument("--eqtls", type=Path, required=True)
    parser.add_argument("--clinical", type=Path, required=True)
    parser.add_argument("--biospecimen", type=Path, nargs="+", required=True)
    parser.add_argument("--actions", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    restricted = args.output / "restricted"
    restricted.mkdir(parents=True, exist_ok=True)

    clinical = pd.read_csv(args.clinical, dtype=str)
    ind_to_pid = dict(zip(clinical["individualID"], _projid(clinical["projid"])))
    pid_sex = dict(zip(_projid(clinical["projid"]),
                       pd.to_numeric(clinical.get("msex"), errors="coerce").map({1.0: "M", 0.0: "F"})))

    # 1. WGS sample -> individual -> projid; WGS dosages renamed by individualID
    wv, ws, wd = read_genotypes(args.wgs)
    w_ind = wgs_individuals(ws, args.biospecimen)
    w_pid = np.array([ind_to_pid.get(i, "") if i else "" for i in w_ind], dtype=object)
    first = (w_ind != "") & ~pd.Series(w_ind).duplicated().to_numpy()
    print(f"WGS individuals with a ROSMAP projid: {len(set(w_pid) - {''})}; individuals with several WGS samples: "
          f"{int(pd.Series(w_ind[w_ind != '']).duplicated().sum())} (first kept)")
    renamed = restricted / (args.wgs.name.replace(".tsv.gz", "") + ".by_individual.tsv.gz")
    with gzip.open(renamed, "wt") as out:
        out.write("\t".join(["snp", "a1", "a2"] + list(w_ind[first])) + "\n")
        for k in range(len(wv)):
            values = ["" if np.isnan(v) else str(int(v)) for v in wd[k, first]]
            out.write("\t".join([wv.at[k, "snp"], wv.at[k, "a1"], wv.at[k, "a2"]] + values) + "\n")
    print(f"wrote {renamed}: {len(wv)} sites x {int(first.sum())} individuals")

    # 2. Array vs WGS
    av, as_, ad = read_genotypes(args.array)
    a_pid = _projid(as_)
    eqtl = pd.read_csv(args.eqtls, sep="\t", dtype={"chrom": str, "rsid": str})
    rows, sign = aligned(eqtl, av, wv, args.wgs_build)
    A = standardize(ad[rows[:, 0]])
    W = standardize(wd[rows[:, 1]]) * sign[:, None]
    ma, mw = np.isfinite(A), np.isfinite(W)
    n_common = ma.T.astype(float) @ mw.astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        R = (np.nan_to_num(A).T @ np.nan_to_num(W)) / n_common  # array samples x WGS samples
    print(f"array vs WGS: {len(rows)} shared non-palindromic eQTL variants; {len(as_)} array x {len(ws)} WGS samples")
    w_of_pid = {}
    for j, p in enumerate(w_pid):
        if p and p not in w_of_pid:
            w_of_pid[p] = j
    people = []
    for i, p in enumerate(a_pid):
        j = w_of_pid.get(p, -1)
        order = np.argsort(-np.nan_to_num(R[i], nan=-9))
        best = int(order[0])
        other = int(next(k for k in order if k != j and (not w_pid[k] or w_pid[k] != p)))
        self_r = float(R[i, j]) if j >= 0 else np.nan
        if j < 0:
            status = "no WGS" if R[i, best] < SAME else "no WGS under this projid, matches a WGS sample"
        elif best == j and self_r >= SAME:
            status = "same person"
        elif self_r >= SAME:
            status = "same person, and a duplicate"
        elif R[i, best] >= SAME:
            status = "matches another person's WGS"
        elif self_r >= RELATIVE:
            status = "relative-level match"
        else:
            status = "matches no WGS"
        people.append({"projid": p, "wgs_sample": ws[j] if j >= 0 else "", "self_r": self_r,
                       "best_wgs_projid": w_pid[best], "best_r": float(R[i, best]),
                       "best_other_r": float(R[i, other]), "best_other_projid": w_pid[other],
                       "n_variants": int(n_common[i, j]) if j >= 0 else 0, "status": status})
    people = pd.DataFrame(people)
    if args.array_chrx and args.array_chrx.exists():
        ax = pd.read_csv(args.array_chrx, sep="\t", dtype={"IID": str})
        people["array_dna_sex"] = pd.Series(dna_sex(ax["het_rate"]), index=_projid(ax["IID"])).reindex(
            people["projid"]).to_numpy()
    if args.wgs_chrx and args.wgs_chrx.exists():
        wx = pd.read_csv(args.wgs_chrx, sep="\t", dtype={"sample": str})
        wsex = dict(zip(wx["sample"], dna_sex(wx["het_rate"])))
        people["wgs_dna_sex"] = [wsex.get(s, "") if s else "" for s in people["wgs_sample"]]
    people["recorded_sex"] = people["projid"].map(pid_sex).fillna("")
    people.to_csv(restricted / "array_vs_wgs_by_person.tsv", sep="\t", index=False, float_format="%.4f")

    summary = people["status"].value_counts().rename_axis("status").reset_index(name="n_array_samples")
    print(summary.to_string(index=False))
    same = people["status"] == "same person"
    print(f"same person: median r {people.loc[same, 'self_r'].median():.3f} (min {people.loc[same, 'self_r'].min():.3f}); "
          f"best other person's WGS among them: median {people.loc[same, 'best_other_r'].median():.3f}, "
          f"p99 {people.loc[same, 'best_other_r'].quantile(0.99):.3f}, max {people.loc[same, 'best_other_r'].max():.3f}")
    sex_rows = []
    for col in ("array_dna_sex", "wgs_dna_sex"):
        if col in people:
            k = people[people[col].isin(["M", "F"]) & people["recorded_sex"].isin(["M", "F"])]
            sex_rows.append({"dna_sex": col, "n_known": len(k), "n_agree": int((k[col] == k["recorded_sex"]).sum()),
                             "n_differ": int((k[col] != k["recorded_sex"]).sum()),
                             "n_unclear": int((people[col] == "?").sum())})
    if {"array_dna_sex", "wgs_dna_sex"} <= set(people.columns):
        k = people[people["array_dna_sex"].isin(["M", "F"]) & people["wgs_dna_sex"].isin(["M", "F"])]
        sex_rows.append({"dna_sex": "array vs wgs", "n_known": len(k),
                         "n_agree": int((k["array_dna_sex"] == k["wgs_dna_sex"]).sum()),
                         "n_differ": int((k["array_dna_sex"] != k["wgs_dna_sex"]).sum()), "n_unclear": 0})
    sex = pd.DataFrame(sex_rows)
    print(sex.to_string(index=False))
    summary.to_csv(args.output / "array_vs_wgs_summary.tsv", sep="\t", index=False)
    sex.to_csv(args.output / "dna_sex_summary.tsv", sep="\t", index=False)

    # 4. The flagged donors
    if args.actions and args.actions.exists():
        actions = pd.read_csv(args.actions, sep="\t", dtype=str, keep_default_na=False)
        actions = actions[actions["dataset"].isin(["rosmap_mit", "rosmap_green", "diverse_cohorts"])]
        flagged = actions[(actions["action"] != "keep") | actions["reason"].isin(["identified weak"])]
        by_pid = people.set_index("projid")
        out = []
        for r in flagged.to_dict("records"):
            row = {k: r[k] for k in ("dataset", "donor", "action", "reason", "label", "best", "expression_sex")}
            for who in ("label", "best"):
                p = _projid(r[who])[0] if r[who] else ""
                info = by_pid.loc[p] if p in by_pid.index else None
                row[f"{who}_array_vs_wgs"] = info["status"] if info is not None else "not on the array"
                row[f"{who}_array_wgs_r"] = info["self_r"] if info is not None else np.nan
                for col in ("array_dna_sex", "wgs_dna_sex", "recorded_sex"):
                    row[f"{who}_{col}"] = info[col] if info is not None and col in info else (
                        pid_sex.get(p, "") if col == "recorded_sex" else "")
            out.append(row)
        flagged_out = pd.DataFrame(out)
        flagged_out.to_csv(restricted / "flagged_donors_wgs.tsv", sep="\t", index=False, float_format="%.3f")
        agg = flagged_out.groupby(["dataset", "reason", "label_array_vs_wgs", "best_array_vs_wgs"]).size() \
            .rename("n_donors").reset_index() if len(flagged_out) else pd.DataFrame()
        agg.to_csv(args.output / "flagged_donors_wgs_summary.tsv", sep="\t", index=False)
        with pd.option_context("display.width", 250, "display.max_rows", 100):
            print(agg.to_string(index=False))
            cols = [c for c in flagged_out.columns if c.endswith(("dna_sex", "recorded_sex"))] + ["expression_sex"]
            sexes = flagged_out.groupby(["dataset", "reason"] + cols).size().rename("n").reset_index()
            print(sexes.to_string(index=False))
            sexes.to_csv(args.output / "flagged_donors_sex_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
