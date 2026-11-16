"""Global donor table: one row per person across the v2 datasets.

Driven by a YAML spec (config/donor_table.yaml). For every dataset, the
spec names the count file, the dataset's own donor key (`donor`, as in its
processed files), the key that identifies a person across datasets
(`person`: the ROSMAP individualID for the Rush donors), and how each
harmonized variable (VARIABLES) is read from its obs columns. A variable is
one source or a list of sources, the first with a value winning per donor:

    {column: msex, recode: {"0": female, "1": male}}   # a label, recoded
    {column: age_death, numeric: true}                 # a number ('90+' -> 90)
    {value: Rush}                                      # the same for every donor
    {checkboxes: {White: "Race (choice=White)", ...}, checked: Checked}

A recode map turns values it does not list into missing values. Donor-level
values are the most common value over the donor's cells.

Where one person is in several datasets, the first dataset in the spec that
has a value gives it; conflicts.tsv counts how often the datasets disagree.

Outputs in the spec's output_dir:
    donor_table.tsv   one row per person (restricted: stays on the cluster)
    overlap.tsv       persons per combination of datasets (exclusive Venn
                      regions) and per number of datasets
    summary.tsv       per dataset and over all persons: values per level,
                      medians and quartiles (aggregates only)
    conflicts.tsv     per variable and pair of datasets: persons in both with
                      a value in both, and how many differ
    provenance.json   the spec and code version

Run: python -m rosmap_processing.pipelines.donor_table config/donor_table.yaml
"""

import argparse
import json
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import h5py
import numpy as np
import pandas as pd
import yaml

from ..utils.logging import get_logger
from ..utils.obs_summary import obs_column_names, read_obs_column
from ..utils.provenance import provenance_record
from .batch_effects import _as_labels, _as_numeric, _first_existing, _label_key

logger = get_logger(__name__)

RACES = ["White", "Black or African American", "Asian", "American Indian or Alaska Native",
         "Native Hawaiian or Pacific Islander", "Other", "Multiple"]
# name -> allowed labels (a list), or (min, max) for numbers, or None (any label)
VARIABLES: Dict[str, Any] = {
    "projid": None,
    "cohort": None,
    "brain_bank": None,
    "sex": ["female", "male"],
    "age_death": (0, 120),
    "race": RACES,
    "hispanic": ["no", "yes"],
    "educ": (0, 30),
    "pmi": (0, 500),
    "apoe_genotype": ["22", "23", "24", "33", "34", "44"],
    "cogdx": (1, 6),
    "dementia": ["no", "yes"],
    "cdr": (0, 5),
    "braak": (0, 6),
    "cerad": (0, 3),
    "thal": (0, 5),
    "adnc": (0, 3),
    "crossdis_dx": None,
}
NUMERIC = {k for k, v in VARIABLES.items() if isinstance(v, tuple)}
DERIVED = ["age_capped", "apoe4"]   # min(age, 90); number of e4 alleles
AGE_CAP = 90
# Differences up to this much are not a conflict (ages capped at 90 first)
TOLERANCE = {"age_death": 1.0, "age_capped": 1.0, "pmi": 0.5, "educ": 0.5}
DATASET_KEYS = {"counts", "donor", "person", "variables"}
SUMMARY_COLUMNS = ["dataset", "variable", "level", "metric", "value"]


# ---------------------------------------------------------------- spec


def load_spec(path: Union[str, Path]) -> Dict[str, Any]:
    with open(path) as f:
        spec = yaml.safe_load(f) or {}
    unknown = set(spec) - {"output_dir", "datasets", "name"}
    if unknown:
        raise ValueError(f"Unknown key(s) in {path}: {sorted(unknown)}")
    if not spec.get("output_dir") or not spec.get("datasets"):
        raise ValueError(f"{path} needs output_dir and datasets")
    spec.setdefault("name", Path(path).stem)
    for name, cfg in spec["datasets"].items():
        bad = set(cfg) - DATASET_KEYS
        if bad:
            raise ValueError(f"Unknown key(s) for dataset {name}: {sorted(bad)}")
        for key in ("counts", "donor", "person"):
            if not cfg.get(key):
                raise ValueError(f"Dataset {name} needs '{key}'")
        cfg.setdefault("variables", {})
        for var, sources in cfg["variables"].items():
            if var not in VARIABLES:
                raise ValueError(f"Dataset {name}: unknown variable {var} (known: {list(VARIABLES)})")
            cfg["variables"][var] = sources if isinstance(sources, list) else [sources]
            for src in cfg["variables"][var]:
                kinds = [k for k in ("column", "value", "checkboxes") if k in src]
                if len(kinds) != 1:
                    raise ValueError(f"Dataset {name}, {var}: a source needs one of column, value, checkboxes")
    return spec


def _source_columns(cfg: Dict[str, Any]) -> List[str]:
    cols = {cfg["donor"], cfg["person"]}
    for sources in cfg["variables"].values():
        for src in sources:
            if "column" in src:
                cols.add(src["column"])
            if "checkboxes" in src:
                cols.update(src["checkboxes"].values())
    return sorted(cols)


# ---------------------------------------------------------------- per donor


def _majority(codes: np.ndarray, n: int, values: np.ndarray) -> Tuple[np.ndarray, int]:
    """Most common non-missing value per donor code (NaN without one), and the
    number of donors with more than one distinct value."""
    numeric = values.dtype.kind == "f"
    out = np.full(n, np.nan, dtype=np.float64 if numeric else object)
    ok = ~pd.isna(values)
    if not ok.any():
        return out, 0
    counts = pd.DataFrame({"d": codes[ok], "v": values[ok]}).groupby(["d", "v"], sort=False).size()
    several = int((counts.groupby(level=0).size() > 1).sum())
    best = counts.groupby(level=0).idxmax()
    out[np.array([d for d, _ in best])] = [v for _, v in best]
    return out, several


def _checkboxes(raw: Dict[str, Any], src: Dict[str, Any], n: int) -> np.ndarray:
    """One label from several checkbox columns: the checked one, 'Multiple'
    if several are, missing if none is."""
    checked = str(src.get("checked", "Checked"))
    hits = np.zeros((n, len(src["checkboxes"])), dtype=bool)
    for j, column in enumerate(src["checkboxes"].values()):
        hits[:, j] = pd.Series(_as_labels(raw[column])).eq(checked).to_numpy()
    names = np.array(list(src["checkboxes"]), dtype=object)
    out = np.full(n, np.nan, dtype=object)
    one = hits.sum(axis=1) == 1
    out[one] = names[hits[one].argmax(axis=1)]
    out[hits.sum(axis=1) > 1] = src.get("several", "Multiple")
    return out


def _cell_values(var: str, src: Dict[str, Any], raw: Dict[str, Any], n: int) -> np.ndarray:
    if "value" in src:
        return np.full(n, float(src["value"]) if var in NUMERIC else str(src["value"]),
                       dtype=np.float64 if var in NUMERIC else object)
    if "checkboxes" in src:
        return _checkboxes(raw, src, n)
    if var in NUMERIC or src.get("numeric"):
        return _as_numeric(raw[src["column"]], src.get("recode"))
    return _as_labels(raw[src["column"]], src.get("recode"))


def _check_values(name: str, var: str, values: np.ndarray) -> np.ndarray:
    """Labels outside the allowed set are an error (naming the labels, which
    are level names, not data); numbers outside the range become missing."""
    allowed = VARIABLES[var]
    if isinstance(allowed, list):
        bad = sorted(set(values[~pd.isna(values)]) - set(allowed))
        if bad:
            raise ValueError(f"[{name}] {var}: values {bad} are not among {allowed}; recode them")
    elif isinstance(allowed, tuple):
        values = values.astype(np.float64)
        out = ~np.isnan(values) & ((values < allowed[0]) | (values > allowed[1]))
        if out.any():
            logger.warning(f"[{name}] {var}: {int(out.sum())} donors outside {allowed} set to missing")
            values = np.where(out, np.nan, values)
    return values


def dataset_donors(name: str, cfg: Dict[str, Any]) -> pd.DataFrame:
    """One row per donor of one dataset: person, donor key, nuclei and the
    harmonized variables."""
    path = _first_existing(cfg["counts"])
    columns = _source_columns(cfg)
    logger.info(f"[{name}] reading {len(columns)} obs columns from {path}")
    with h5py.File(path, "r") as f:
        missing = sorted(set(columns) - set(obs_column_names(f)))
        if missing:
            raise ValueError(f"[{name}] obs has no column(s) {missing}")
        raw = {c: read_obs_column(f, c) for c in columns}
    donors = _as_labels(raw[cfg["donor"]])
    has_donor = ~pd.isna(donors)
    codes, keys = pd.factorize(pd.Series(donors[has_donor]))
    n = len(keys)
    logger.info(f"[{name}] {len(donors)} cells, {int((~has_donor).sum())} without a donor, {n} donors")
    raw = {c: v[has_donor] for c, v in raw.items()}  # arrays and Categoricals alike
    n_cells = int(has_donor.sum())

    table = pd.DataFrame({"donor": np.asarray(keys, dtype=object),
                          "n_nuclei": np.bincount(codes, minlength=n)})
    person, several = _majority(codes, n, _as_labels(raw[cfg["person"]]))
    if several:
        logger.warning(f"[{name}] {several} donors have cells with more than one {cfg['person']}")
    table["person_id"] = person
    for var, sources in cfg["variables"].items():
        merged = np.full(n, np.nan, dtype=np.float64 if var in NUMERIC else object)
        for src in sources:
            values, several = _majority(codes, n, _cell_values(var, src, raw, n_cells))
            if several:
                logger.warning(f"[{name}] {var} from {src.get('column', src)}: {several} donors have "
                               f"more than one value (the most common is used)")
            fill = pd.isna(merged) & ~pd.isna(values)
            merged[fill] = values[fill]
        merged = _check_values(name, var, merged)
        table[var] = merged
        logger.info(f"[{name}]   {var}: {int((~pd.isna(merged)).sum())} of {n} donors have a value")

    no_person = table["person_id"].isna()
    if no_person.any():
        logger.warning(f"[{name}] {int(no_person.sum())} donors without a {cfg['person']} are left out")
        table = table[~no_person]
    dup = table["person_id"].duplicated(keep=False)
    if dup.any():
        logger.warning(f"[{name}] {int(table.loc[dup, 'person_id'].nunique())} persons have several donor "
                       f"keys; keeping the one with the most nuclei")
        table = table.sort_values("n_nuclei", ascending=False).drop_duplicates("person_id")
    return table.reset_index(drop=True)


# ---------------------------------------------------------------- merge


def _normalize(var: str, values: pd.Series) -> pd.Series:
    if var == "projid":  # '00123456', 123456.0 and '123456' alike
        return values.map(lambda v: _label_key(int(v)) if _label_key(v) is not None
                          and str(_label_key(v)).isdigit() else _label_key(v))
    return values


def _derived(table: pd.DataFrame) -> None:
    if "age_death" in table:
        table["age_capped"] = np.minimum(table["age_death"].astype(np.float64), AGE_CAP)
    if "apoe_genotype" in table:
        table["apoe4"] = table["apoe_genotype"].map(lambda g: float(str(g).count("4")) if pd.notna(g) else np.nan)


def _differ(var: str, a: pd.Series, b: pd.Series) -> np.ndarray:
    if var in NUMERIC or var in DERIVED:
        return np.abs(a.astype(np.float64) - b.astype(np.float64)).to_numpy() > TOLERANCE.get(var, 1e-9)
    return (a.astype(str) != b.astype(str)).to_numpy()


def merge_datasets(per_dataset: Dict[str, pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The donor table (first dataset with a value wins) and the conflicts."""
    names = list(per_dataset)
    for t in per_dataset.values():
        for var in VARIABLES:
            if var in t:
                t[var] = _normalize(var, t[var])
        _derived(t)
    persons = list(dict.fromkeys(p for t in per_dataset.values() for p in t["person_id"]))
    table = pd.DataFrame(index=pd.Index(persons, name="person_id"))
    indexed = {k: t.set_index("person_id") for k, t in per_dataset.items()}
    for k in names:
        table[f"in_{k}"] = table.index.isin(indexed[k].index)
    table["n_datasets"] = table[[f"in_{k}" for k in names]].sum(axis=1)
    table["datasets"] = table[[f"in_{k}" for k in names]].apply(
        lambda row: "+".join(k for k in names if row[f"in_{k}"]), axis=1)
    for k in names:
        table[f"{k}_donor"] = indexed[k]["donor"].reindex(table.index)
        table[f"{k}_n_nuclei"] = indexed[k]["n_nuclei"].reindex(table.index).astype("Int64")

    variables = [v for v in list(VARIABLES) + DERIVED if any(v in t for t in indexed.values())]
    conflicts = []
    for var in variables:
        merged = pd.Series(np.nan, index=table.index, dtype=np.float64 if var in NUMERIC | set(DERIVED)
                           else object)
        source = pd.Series(np.nan, index=table.index, dtype=object)
        for k in names:
            if var not in indexed[k]:
                continue
            values = indexed[k][var].reindex(table.index)
            fill = merged.isna() & values.notna()
            merged[fill] = values[fill]
            source[fill] = k
        table[var] = merged
        if var not in DERIVED:
            table[f"{var}_from"] = source
        for a, b in combinations(names, 2):
            if var not in indexed[a] or var not in indexed[b]:
                continue
            both = indexed[a].index.intersection(indexed[b].index)
            va, vb = indexed[a].loc[both, var], indexed[b].loc[both, var]
            ok = va.notna().to_numpy() & vb.notna().to_numpy()
            if not len(both):
                continue
            conflicts.append({"variable": var, "dataset_a": a, "dataset_b": b, "n_persons_in_both": len(both),
                              "n_compared": int(ok.sum()), "n_differ": int(_differ(var, va[ok], vb[ok]).sum())})
    return table.reset_index(), pd.DataFrame(conflicts, columns=[
        "variable", "dataset_a", "dataset_b", "n_persons_in_both", "n_compared", "n_differ"])


# ---------------------------------------------------------------- summaries


def overlap_table(table: pd.DataFrame, names: List[str]) -> pd.DataFrame:
    rows = []
    member = table[[f"in_{k}" for k in names]].to_numpy()
    for r in range(1, len(names) + 1):
        for combo in combinations(range(len(names)), r):
            want = np.zeros(len(names), dtype=bool)
            want[list(combo)] = True
            rows.append({"kind": "region", "datasets": "+".join(names[i] for i in combo),
                         "n_persons": int((member == want).all(axis=1).sum())})
    for k in range(1, len(names) + 1):
        rows.append({"kind": "n_datasets", "datasets": str(k), "n_persons": int((member.sum(axis=1) == k).sum())})
    rows.append({"kind": "total", "datasets": "any", "n_persons": len(table)})
    return pd.DataFrame(rows)


def summary_table(table: pd.DataFrame, names: List[str]) -> pd.DataFrame:
    """Per dataset (its persons, with the merged values) and over all persons:
    counts per level, or median and quartiles. No minima or maxima, which are
    single persons' values."""
    rows = []
    groups = [("all", np.ones(len(table), dtype=bool))] + [(k, table[f"in_{k}"].to_numpy()) for k in names]
    variables = [v for v in list(VARIABLES) + DERIVED if v in table and v != "projid"]
    for group, mask in groups:
        sub = table[mask]
        rows.append({"dataset": group, "variable": "persons", "level": "", "metric": "n", "value": len(sub)})
        for var in variables:
            values = sub[var]
            rows.append({"dataset": group, "variable": var, "level": "", "metric": "n_with_value",
                         "value": int(values.notna().sum())})
            if var in NUMERIC or var in DERIVED:
                v = values.dropna().astype(np.float64)
                if var in ("cogdx", "braak", "cerad", "thal", "adnc", "apoe4", "cdr"):  # ordered stages
                    for level, count in v.value_counts().sort_index().items():
                        rows.append({"dataset": group, "variable": var, "level": f"{level:g}",
                                     "metric": "count", "value": int(count)})
                if len(v) >= 5:
                    for metric, q in (("q1", 0.25), ("median", 0.5), ("q3", 0.75)):
                        rows.append({"dataset": group, "variable": var, "level": "", "metric": metric,
                                     "value": float(v.quantile(q))})
            else:
                for level, count in values.value_counts().items():
                    rows.append({"dataset": group, "variable": var, "level": str(level),
                                 "metric": "count", "value": int(count)})
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


# ---------------------------------------------------------------- run


def run_donor_table(spec_path: Union[str, Path], overwrite: bool = False) -> Path:
    spec = load_spec(spec_path)
    out = Path(spec["output_dir"])
    target = out / "donor_table.tsv"
    if target.exists() and not overwrite:
        raise FileExistsError(f"{target} exists; pass --overwrite to replace it")
    out.mkdir(parents=True, exist_ok=True)
    names = list(spec["datasets"])
    per_dataset = {k: dataset_donors(k, cfg) for k, cfg in spec["datasets"].items()}
    table, conflicts = merge_datasets(per_dataset)
    table.to_csv(target, sep="\t", index=False)
    overlap = overlap_table(table, names)
    overlap.to_csv(out / "overlap.tsv", sep="\t", index=False)
    summary_table(table, names).to_csv(out / "summary.tsv", sep="\t", index=False)
    conflicts.to_csv(out / "conflicts.tsv", sep="\t", index=False)
    with open(out / "provenance.json", "w") as f:
        json.dump(provenance_record("donor_table", {"spec": str(spec_path)}) | {"spec": spec}, f, indent=2,
                  default=str)
    shared = int((table["n_datasets"] > 1).sum())
    logger.info(f"Wrote {target}: {len(table)} persons, {shared} in more than one dataset")
    for _, row in conflicts[conflicts["n_differ"] > 0].iterrows():
        logger.info(f"  conflict {row['variable']}: {row['dataset_a']} vs {row['dataset_b']} differ for "
                    f"{row['n_differ']} of {row['n_compared']} persons")
    return target


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("spec", help="YAML spec, e.g. config/donor_table.yaml")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing donor table")
    args = parser.parse_args(argv)
    from ..utils.logging import setup_logging

    setup_logging(level="INFO")
    run_donor_table(args.spec, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
