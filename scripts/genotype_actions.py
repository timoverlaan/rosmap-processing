"""Per-donor action from the genotype check (`pipeline qc-check genotype`):
keep, relabel or drop, with the motivation and the evidence.

Reads the restricted outputs of one or more genotype-check folders (the calls
per dataset, and the cross-dataset rows of a multi-dataset run) and writes
OUTPUT/restricted/genotype_actions.tsv (one row per donor; donor IDs, so
restricted) and OUTPUT/genotype_actions_summary.tsv (counts per dataset,
action and reason; no IDs). Prints only the summary.

    python scripts/genotype_actions.py data/v2/qc/genotype_v2 data/v2/qc/genotype_psychad \
        --wgs amp1_wgs=data/v2/qc/genotype_wgs/check_amp1 divco_wgs=data/v2/qc/genotype_wgs/check_divco \
        --array-vs-wgs data/v2/qc/genotype_wgs/amp1 data/v2/qc/genotype_wgs/divco \
        --id-map data/v2/download/ROSMAP/ROSMAP_clinical.csv --output data/v2/qc/genotype_actions

Each genotype source (the positional folders: array or a cohort's own
genotypes; --wgs: the same check against WGS) gives a decision per donor by
the rules below; `combine` then merges a donor's decisions: a full-strength
own-genotype match in any source keeps the label (if another source finds a
swap, the WGS wins only where the array entry under the label is shown wrong
by WGS, else the donor is dropped); a relabel needs its sources to agree on
the person, and a swap unconfirmed in one source is confirmed when another
source that genotypes the label matches the same person; otherwise the most
specific drop, else the best-supported keep. Sources where the label has no
genotype and nothing matches are ignored. A mutual swap (two donors each
holding the other's sample) counts as evidence for an RNA swap.

Rules per source (a sample "matches" a genotype when the check calls it so):
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
            if partner["call"] == "other person" and partner["best"] == r["label"] and r["label"]:
                pro.append(f"the donor labelled as {target} ({partner['donor']}) holds a sample of this donor's "
                           f"labelled person: a mutual swap")
            else:
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


IDENTITY_RANK = {"verified": 0, "probable": 1, "unverified": 2, "unchecked": 3}
DROP_RANK = ["contradicted", "genotype sources disagree", "no signal despite nuclei", "second sample of a donor",
             "unconfirmed", "unclear, ambiguous", "weak match"]
LABELLED = ("own genotype", "other person", "unclear", "no signal")  # the label has a genotype in the source


def source_decisions(folder: Path, source: str, to_iid: dict):
    """{(dataset, donor): decision of this source} for one genotype-check folder."""
    out = {}
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
            out[(name, r["donor"])] = {"source": source, "r": r, "action": action, "reason": reason,
                                       "identity": identity, "motivation": motivation, "best_iid": best_iid,
                                       "target": target}
    return out


def combine(decisions: list, array_status: dict):
    """One donor's final (action, reason, identity, motivation, deciding decision)
    from the decisions of its genotype sources (array first, then WGS).

    - A full-strength own-genotype match in any source verifies the label,
      unless another source finds a swap or a second sample: then the WGS wins
      when the array genotype under the label is someone else's (array vs
      WGS), else the sources disagree and the donor is dropped.
    - A relabel needs every relabelling source to name the same person; a
      swap that one source leaves unconfirmed is confirmed when another
      source, in which the label is genotyped, matches the same person.
    - Otherwise the most specific drop wins, else the best-supported keep."""
    if len(decisions) == 1:
        d = decisions[0]
        return d["action"], d["reason"], d["identity"], d["motivation"], d
    tag = lambda d: f"[{d['source']}] "  # noqa: E731

    def text(ds):
        return " ".join(tag(d) + (d["motivation"] or f"{d['r']['call']}.") for d in ds)

    own = [d for d in decisions if d["reason"] in ("own genotype", "own genotype, weak", "own genotype, sex record")]
    full_own = [d for d in own if d["reason"] != "own genotype, weak"]
    swaps = [d for d in decisions if d["action"] == "relabel" or d["reason"] == "second sample of a donor"]
    if full_own or (own and not swaps):
        chosen = (full_own or own)[0]
        if swaps and len(full_own) >= 2 and len(swaps) == 1:
            return ("keep", chosen["reason"], "verified",
                    f"{', '.join(d['source'] for d in full_own)} match its own genotype; only "
                    f"{tag(swaps[0])}disagrees: {swaps[0]['motivation']}", chosen)
        if swaps:
            wgs_own = [d for d in full_own if d["source"] != "array"]
            array_swaps = [d for d in swaps if d["source"] == "array"]
            status = array_status.get(_projid(array_swaps[0]["r"]["label"])) if array_swaps else None
            if wgs_own and len(array_swaps) == len(swaps) and status and status != "same person":
                return ("keep", "own genotype, array entry wrong", "verified",
                        f"{tag(wgs_own[0])}matches its own genotype (z {_f(wgs_own[0]['r']['own_z'])}). "
                        f"{tag(array_swaps[0])}{array_swaps[0]['motivation']} But the array genotype under this "
                        f"label is not this person's own (array vs WGS: {status}): the array entry is wrong, "
                        f"not the RNA.", wgs_own[0])
            return ("drop", "genotype sources disagree", "uncertain",
                    f"The sources disagree on who this is: {text(own + swaps)}", chosen)
        notes = [tag(d) + d["motivation"] for d in decisions if d is not chosen and d["action"] != "keep"]
        reason = chosen["reason"]
        if any(d["source"] == "array" and d["action"] != "keep" for d in decisions):
            status = array_status.get(_projid(decisions[0]["r"]["label"]))
            if status and status not in ("same person", "no WGS"):
                notes.append(f"the array genotype under this label is not this person's (array vs WGS: {status}), "
                             f"so the array's mismatch is an error on the array side")
                reason = "own genotype, array entry wrong"
        others = [d["source"] for d in decisions if d is not chosen and d["reason"].startswith("own genotype")]
        motivation = "; ".join(filter(None, [
            chosen["motivation"],
            f"confirmed by {', '.join(others)}" if others else "",
            ("other sources: " + " ".join(notes)) if notes else ""]))
        return "keep", reason, "verified", motivation, chosen

    relabels = [d for d in decisions if d["action"] == "relabel"]
    if relabels:
        targets = {d["best_iid"] or d["r"]["best"] for d in relabels}
        others = [d for d in decisions if d["reason"] in ("contradicted", "second sample of a donor")]
        if len(targets) == 1 and not others:
            d = relabels[0]
            agree = [x["source"] for x in decisions if x is not d and x["r"]["call"] in ("other person", "identified")
                     and (x["best_iid"] or x["r"]["best"]) in targets]
            return ("relabel", d["reason"], "wrong",
                    d["motivation"] + (f" {', '.join(agree)} matches the same person." if agree else ""), d)
        return ("drop", "genotype sources disagree", "uncertain",
                f"The sources disagree on who this is: {text(relabels + others)}", relabels[0])

    # A swap left unconfirmed in one source is confirmed by another source that
    # genotypes the label and matches the same person
    matched = [d for d in decisions if d["r"]["call"] in ("other person", "identified") and not d["r"]["weak"]]
    labelled = [d for d in matched if d["r"]["call"] == "other person"]
    seconds = [d for d in decisions if d["reason"] == "second sample of a donor"]
    if not seconds and len(matched) >= 2 and labelled and len({d["best_iid"] or d["r"]["best"] for d in matched}) == 1 \
            and not any(d["reason"] == "contradicted" for d in decisions):
        d = labelled[0]
        return ("relabel", "confirmed by two genotype sources", "wrong",
                f"{text(matched)} Two genotype sources match the sample to the same person, and the label's own "
                f"genotype does not match it, so the RNA is that person's.", d)

    drops = [d for d in decisions if d["action"] == "drop"]
    if drops:
        drops.sort(key=lambda d: DROP_RANK.index(d["reason"]) if d["reason"] in DROP_RANK else len(DROP_RANK))
        d = drops[0]
        return "drop", d["reason"], d["identity"], text(drops), d
    keeps = sorted(decisions, key=lambda d: IDENTITY_RANK.get(d["identity"], 9))
    d = keeps[0]
    return "keep", d["reason"], d["identity"], text([x for x in keeps if x["motivation"]]), d


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folders", type=Path, nargs="+", help="genotype-check output folders (array or own genotypes)")
    parser.add_argument("--wgs", nargs="*", default=[], metavar="NAME=FOLDER",
                        help="genotype-check folders run against WGS (labels by individualID), as name=folder")
    parser.add_argument("--array-vs-wgs", type=Path, nargs="*", default=[],
                        help="scripts/genotype_wgs_check.py output folders (array vs WGS per person)")
    parser.add_argument("--id-map", type=Path, help="CSV with individualID and projid (for projid labels)")
    parser.add_argument("--output", type=Path, help="default: the first folder")
    args = parser.parse_args()

    to_iid = {}
    if args.id_map:
        ids = pd.read_csv(args.id_map, dtype=str, usecols=["individualID", "projid"]).dropna()
        for iid, pid in zip(ids["individualID"], ids["projid"]):
            to_iid.setdefault(_projid(pid), iid)
    array_status = {}
    for folder in args.array_vs_wgs:
        path = folder / "restricted" / "array_vs_wgs_by_person.tsv"
        for r in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).to_dict("records"):
            if r["status"] != "no WGS" or r["projid"] not in array_status:
                array_status[r["projid"]] = r["status"]

    sources = [("array", f) for f in args.folders] + [tuple(s.split("=", 1)) for s in args.wgs]
    per_donor = defaultdict(list)
    for source, folder in sources:
        for key, d in source_decisions(Path(folder), source, to_iid).items():
            per_donor[key].append(d)

    rows = []
    for (name, donor), decisions in per_donor.items():
        # A source in which the label has no genotype and nothing matches adds nothing
        useful = [d for d in decisions if d["r"]["call"] not in ("not genotyped", "no data")] or decisions[:1]
        action, reason, identity, motivation, d = combine(useful, array_status)
        r = d["r"]
        row = {"dataset": name, "donor": donor, "match": r["match"], "label": decisions[0]["r"]["label"],
               "call": r["call"], "action": action, "reason": reason, "identity": identity,
               "relabel_to": r["best"] if action == "relabel" else "",
               "relabel_to_individualID": d["best_iid"] if action == "relabel" else "",
               "decided_by": d["source"], "motivation": motivation, "own_z": r["own_z"], "best_z": r["best_z"],
               "next_z": r["next_z"], "expected_z": r["expected_z"], "strength": r["strength"],
               "n_nuclei": r["n_nuclei"], "expression_sex": r["expression_sex"], "label_sex": r["label_sex"],
               "best_sex": r["best_sex"], "best": r["best"], "best_is_donor": r["best_is_donor"]}
        for x in decisions:
            row[f"{x['source']}_call"] = x["r"]["call"]
            row[f"{x['source']}_reason"] = x["reason"]
        rows.append(row)
    out = pd.DataFrame(rows)
    out["relabel_to_individualID"] = out["relabel_to_individualID"].fillna("")

    # A relabel onto a person kept under their own label, or shared by two donors,
    # is dropped (people compared by individualID where known)
    person = out["relabel_to_individualID"].where(out["relabel_to_individualID"] != "", out["relabel_to"])
    for name, table in out.groupby("dataset"):
        kept = set(table.loc[table["action"] == "keep", "match"].fillna("")) - {""}
        kept |= set(table.loc[table["action"] == "keep", "label"].fillna("")) - {""}
        targets = person[table.index][table["action"] == "relabel"]
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
