"""Per-donor action from the genotype check (`pipeline qc-check genotype`):
keep, relabel or drop, with the motivation and the evidence.

Reads the restricted outputs of one or more genotype-check folders (the calls
per dataset, and the cross-dataset rows of a multi-dataset run) and writes
OUTPUT/restricted/genotype_actions.tsv (one row per donor; donor IDs, so
restricted) and OUTPUT/genotype_actions_summary.tsv (counts per dataset,
action and reason; no IDs). Prints only the summary.

    python scripts/genotype_actions.py data/v2/qc/genotype_v2 data/v2/qc/genotype_psychad \
        --id-map data/v2/download/ROSMAP/ROSMAP_clinical.csv --output data/v2/qc/genotype_actions

Rules (a sample "matches" a genotype when the check calls it so):
- own genotype -> keep (a weak match or a sex record that disagrees is noted).
- unclear -> keep if the label is still the best match (identity probable),
  else drop (identity ambiguous).
- no signal -> keep (unverified), except drop when the donor has plenty of
  nuclei (>= MISMATCH_NUCLEI) and its own genotype scores below MISMATCH_Z,
  or its expression sex contradicts the record.
- other person / identified (the sample matches another genotyped person X):
  - X is another donor of the dataset whose sample matches X too -> drop
    (a second sample of X);
  - identified at under GENOTYPE_WEAK of the expected z -> keep, unverified
    (a relative of X or a mixture, not a clean swap; evidence that the sample
    is not X supports the label), unless the sex contradicts the record;
  - evidence against an RNA swap -> drop: the expression sex fits the label's
    record but not X's, or another dataset's sample of the labelled person
    also matches X (then the label's genotype or ID is wrong, or X is a
    relative);
  - evidence for an RNA swap -> relabel to X: the expression sex contradicts
    the label's record and fits X's, or another dataset's sample of the
    labelled person matches the label (or, without a genotype, not X);
  - otherwise -> drop (unconfirmed: an RNA swap and a genotype filed under the
    wrong ID look alike).
- not genotyped / no data -> keep (unchecked), except drop when the expression
  sex contradicts the record.
A relabel whose target is kept under its own label elsewhere in the dataset,
or that two donors share, becomes a drop.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

GENOTYPE_WEAK = 0.65   # as in data_qc: a confident call below this share of the expected z is weak
NO_MATCH_Z = 3.0       # another dataset's sample "does not match" a genotype below this z
MATCH_Z = 6.0          # ... and matches it at this z or more
MISMATCH_NUCLEI = 2000  # with this many nuclei a sample of the labelled person scores z ~10
MISMATCH_Z = 2.0
CONFIDENT = ("own genotype", "other person", "identified")


def _projid(value) -> str:
    m = re.search(r"(\d+)(?:\.0+)?$", str(value))
    return str(int(m.group(1))) if m else ""


def _known(x) -> bool:
    return x in ("M", "F")


def _f(x, digits=1) -> str:
    return "n/a" if x is None or not np.isfinite(x) else f"{x:.{digits}f}"


def read_folder(folder: Path):
    """(calls per dataset, cross rows keyed by (dataset, donor), label column kind)."""
    restricted = folder / "restricted"
    calls, kinds = {}, set()
    for path in sorted(restricted.glob("genotype_calls_*.tsv")):
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        kind = "projid" if "label_projid" in table else "sample"
        kinds.add(kind)
        table = table.rename(columns={f"label_{kind}": "label", f"best_{kind}": "best"})
        for c in ("own_z", "best_z", "next_z", "expected_z", "strength", "n_nuclei"):
            table[c] = pd.to_numeric(table[c], errors="coerce")
        table["weak"] = table["weak"].eq("True")
        table["kind"] = kind
        calls[path.stem.replace("genotype_calls_", "")] = table
    cross = defaultdict(list)
    path = restricted / "genotype_cross_by_donor.tsv"
    if path.exists():
        rows = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        for c in ("own_z_b", "z_b_at_best_a", "best_z_a"):
            rows[c] = pd.to_numeric(rows[c], errors="coerce")
        for r in rows.to_dict("records"):
            cross[(r["dataset_a"], r["donor_a"])].append(r)
    return calls, cross


def _evidence(r: dict, cross_rows: list):
    """Evidence about a donor whose sample matches another person X: texts for
    an RNA swap, against it from sex, against it from another dataset, and
    other observations (another dataset's z against X between NO_MATCH_Z and
    MATCH_Z)."""
    ex, lab, bs = r["expression_sex"], r["label_sex"], r["best_sex"]
    pro, con_sex, con_cross, other = [], [], [], []
    if _known(ex) and _known(lab) and _known(bs):
        if ex != lab and ex == bs:
            pro.append(f"its expression sex ({ex}) contradicts the label's record ({lab}) and fits the matched "
                       f"person's")
        elif ex == lab and ex != bs:
            con_sex.append(f"its expression sex ({ex}) fits the label's record, not the matched person's ({bs})")
    for c in cross_rows:
        z_x, own_b = c["z_b_at_best_a"], c["own_z_b"]
        where = f"{c['dataset_b']} has a sample of the labelled person ({c['donor_b']}, called {c['call_b']})"
        if not np.isfinite(z_x):
            continue
        if z_x >= MATCH_Z:
            con_cross.append(f"{where} that also matches the matched person's genotype (z {_f(z_x)})")
        elif c["call_b"] == "own genotype" and z_x < NO_MATCH_Z:
            pro.append(f"{where} that matches the label's genotype (z {_f(own_b)}) and not the matched "
                       f"person's (z {_f(z_x)})")
        elif c["call_b"] in ("not genotyped", "identified") and z_x < NO_MATCH_Z and r["call"] == "identified":
            pro.append(f"{where} that does not match the matched person's genotype (z {_f(z_x)})")
        elif z_x >= NO_MATCH_Z:
            other.append(f"{where} that scores z {_f(z_x)} against the matched person's genotype "
                         f"(suggestive, below {MATCH_Z:.0f})")
    return pro, con_sex, con_cross, other


def decide(r: dict, partner: dict, cross_rows: list, target: str):
    """(action, reason, identity, motivation) of one donor."""
    call = r["call"]
    ex, lab = r["expression_sex"], r["label_sex"]
    sex_off = _known(ex) and _known(lab) and ex != lab
    n = int(r["n_nuclei"]) if np.isfinite(r["n_nuclei"]) else 0
    own_z, best_z, exp_z, strength = r["own_z"], r["best_z"], r["expected_z"], r["strength"]
    label_best = r["best"] == r["label"] and r["label"] != ""
    sex_text = f"its expression sex ({ex}) contradicts the recorded sex ({lab})"
    match_text = (f"The sample matches the genotype of {target} at z {_f(best_z)} ({_f(100 * strength, 0)}% of "
                  f"the z {_f(exp_z)} expected at its {n} nuclei)")

    if call == "own genotype":
        notes = []
        if r["weak"]:
            confirm = [c["dataset_b"] for c in cross_rows if c["call_b"] == "own genotype"]
            notes.append(f"weak match: z {_f(own_z)}, {_f(100 * strength, 0)}% of the expected {_f(exp_z)}"
                         + (f"; {', '.join(confirm)} confirms the person" if confirm else ""))
        if sex_off:
            notes.append(f"matches its own genotype, but {sex_text}: the sex record is probably wrong")
        reason = "own genotype, sex record" if sex_off else "own genotype, weak" if r["weak"] else "own genotype"
        return "keep", reason, "verified", "; ".join(notes)

    if call == "unclear":
        if label_best and not sex_off:
            return ("keep", "unclear, label best", "probable",
                    f"its label is the best match (z {_f(own_z)}) but leads the next person by less than 3 "
                    f"(next z {_f(r['next_z'])})")
        why = (f"{sex_text}, and its label matches only weakly (z {_f(own_z)})" if label_best else
               f"another person's genotype ({target}, z {_f(best_z)}) scores above its own label (z {_f(own_z)}), "
               f"but not decisively")
        return "drop", "unclear, ambiguous", "ambiguous", f"Identity ambiguous: {why}."

    if call == "no signal":
        if np.isfinite(own_z) and own_z < MISMATCH_Z and n >= MISMATCH_NUCLEI:
            return ("drop", "no signal despite nuclei", "probably wrong",
                    f"With {n} nuclei a sample of the labelled person would score about z {_f(exp_z)}, but its "
                    f"own genotype scores z {_f(own_z)} and no one else reaches z 6 (best z {_f(best_z)}): the "
                    f"sample is probably someone without a genotype here.")
        if sex_off:
            return ("drop", "no signal, sex off", "probably wrong",
                    f"Too little expression signal to check the identity (best z {_f(best_z)}, {n} nuclei), and "
                    f"{sex_text}.")
        return "keep", "no signal", "unverified", f"too little signal to check (best z {_f(best_z)}, {n} nuclei)"

    if call in ("other person", "identified"):
        if partner is not None and partner["call"] == "own genotype":
            return ("drop", "second sample of a donor", "wrong",
                    f"{match_text}; {target} is also donor {partner['donor']} of this dataset, whose sample "
                    f"matches their genotype too. This is a second sample of that person; the labelled "
                    f"person's own sample is missing.")
        pro, con_sex, con_cross, other = _evidence(r, cross_rows)
        if call == "identified" and r["weak"]:
            # A weak match to X with no genotype for the label: a relative of X or
            # a mixture; evidence that the sample is not X supports the label
            if sex_off:
                return ("drop", "identified weak, sex off", "probably wrong",
                        f"The label has no genotype; the sample matches {target} only weakly (z {_f(best_z)}, "
                        f"{_f(100 * strength, 0)}% of expected), and {sex_text}.")
            notes = [f"label has no genotype; the sample matches {target} at about half strength "
                     f"(z {_f(best_z)}, {_f(100 * strength, 0)}% of expected): a relative of {target} or a "
                     f"mixture, not a clean swap"]
            notes += [f"{t}, which argues against a swap to {target}" for t in con_sex]
            notes += [f"{t}, which points to a relative rather than a mixture" for t in con_cross + other]
            if not con_sex and _known(ex) and ex == lab:
                notes.append("its expression sex fits the label")
            return "keep", "identified weak", "unverified", "; ".join(notes)
        if con_sex or con_cross:
            return ("drop", "contradicted", "uncertain",
                    f"{match_text}, but {'; '.join(con_sex + con_cross)}. The error may be in the genotype file, "
                    f"the ID map or a sex record rather than the RNA (or {target} is a relative), so neither "
                    f"identity can be trusted.")
        if r["weak"]:
            return ("drop", "weak match", "uncertain",
                    f"{match_text}; the match is weak, so it is no basis for a relabel, and the label does not "
                    f"match.")
        chain = ""
        if partner is not None:
            chain = (f" The donor labelled as {target} in this dataset ({partner['donor']}) holds a sample of "
                     f"someone else (called {partner['call']}), so the labels are shifted.")
        if pro:
            label_part = (f" and its label's genotype at z {_f(own_z)}" if np.isfinite(own_z) else
                          " (its label has no genotype)")
            return ("relabel", "confirmed swap", "wrong",
                    f"{match_text}{label_part}. RNA swap confirmed: {'; '.join(pro)}.{chain}")
        unlabelled = "its label has no genotype" if call == "identified" else f"its label scores z {_f(own_z)}"
        return ("drop", "unconfirmed", "wrong or uncertain",
                f"{match_text}, and {unlabelled}. No independent check (no informative sex difference, no "
                f"sample of the person in another dataset): an RNA swap and a genotype filed under the wrong ID "
                f"look alike, so it is not relabelled.{chain}")

    # not genotyped, no data
    if sex_off:
        return ("drop", f"{call}, sex off", "probably wrong",
                f"The label has no genotype to check the identity, and {sex_text}.")
    return "keep", call, "unchecked", "no genotype to check" if call == "not genotyped" else "no expression data"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folders", type=Path, nargs="+", help="genotype-check output folders")
    parser.add_argument("--id-map", type=Path, help="CSV with individualID and projid (for projid labels)")
    parser.add_argument("--output", type=Path, help="default: the first folder")
    args = parser.parse_args()

    to_iid = {}
    if args.id_map:
        ids = pd.read_csv(args.id_map, dtype=str, usecols=["individualID", "projid"]).dropna()
        for iid, pid in zip(ids["individualID"], ids["projid"]):
            to_iid.setdefault(_projid(pid), iid)

    rows = []
    for folder in args.folders:
        calls, cross = read_folder(folder)
        for name, table in calls.items():
            by_label = {lab: r for r in table.to_dict("records") if (lab := r["label"])}
            by_donor = {r["donor"]: r for r in table.to_dict("records")}
            for r in table.to_dict("records"):
                partner = by_donor.get(r["best_is_donor"]) if r["best_is_donor"] else None
                if partner is None and r["best"] in by_label and r["best"] != r["label"]:
                    partner = by_label[r["best"]]
                best_iid = to_iid.get(r["best"], "") if r["kind"] == "projid" else r["best"]
                target = (f"{best_iid or 'an unmapped individualID'} (projid {r['best']})" if r["kind"] == "projid"
                          else r["best"])
                action, reason, identity, motivation = decide(r, partner, cross.get((name, r["donor"]), []), target)
                rows.append({"dataset": name, "donor": r["donor"], "match": r["match"], "label": r["label"],
                             "call": r["call"], "action": action, "reason": reason, "identity": identity,
                             "relabel_to": r["best"] if action == "relabel" else "",
                             "relabel_to_individualID": best_iid if action == "relabel" else "",
                             "motivation": motivation, "own_z": r["own_z"], "best_z": r["best_z"],
                             "next_z": r["next_z"], "expected_z": r["expected_z"], "strength": r["strength"],
                             "n_nuclei": r["n_nuclei"], "expression_sex": r["expression_sex"],
                             "label_sex": r["label_sex"], "best_sex": r["best_sex"], "best": r["best"],
                             "best_is_donor": r["best_is_donor"]})
    out = pd.DataFrame(rows)

    # A relabel onto a person kept under their own label, or shared by two donors, is dropped
    for name, table in out.groupby("dataset"):
        kept = set(table.loc[table["action"] == "keep", "label"]) - {""}
        targets = table.loc[table["action"] == "relabel", "relabel_to"]
        for i, t in targets.items():
            why = ("is kept under their own label in this dataset" if t in kept else
                   "is the relabel target of another donor too" if (targets == t).sum() > 1 else "")
            if why:
                out.loc[i, ["action", "reason", "relabel_to", "relabel_to_individualID"]] = \
                    ["drop", "relabel conflict", "", ""]
                out.loc[i, "motivation"] += f" Not relabelled: the target {why}."

    output = args.output or args.folders[0]
    (output / "restricted").mkdir(parents=True, exist_ok=True)
    out.to_csv(output / "restricted" / "genotype_actions.tsv", sep="\t", index=False, float_format="%.3g")
    summary = out.groupby(["dataset", "action", "reason"]).size().rename("n_donors").reset_index()
    summary.to_csv(output / "genotype_actions_summary.tsv", sep="\t", index=False)
    with pd.option_context("display.max_rows", 200, "display.width", 200):
        print(summary.to_string(index=False))
        print(out.pivot_table(index="dataset", columns="action", values="donor", aggfunc="size", fill_value=0))
    print(f"wrote {output / 'restricted' / 'genotype_actions.tsv'} ({len(out)} donors) and "
          f"{output / 'genotype_actions_summary.tsv'}")


if __name__ == "__main__":
    main()
