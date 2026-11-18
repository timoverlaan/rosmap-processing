"""Data QC: how well the v2 count files were processed (aggregate tables only).

Driven by a YAML spec per dataset (config/data_qc/), which names the dataset
block of a batch-effects spec (count file, gene map, donor, cell-type rules,
covariates) and adds what the QC needs: the authors' per-cell QC values, the
authors' label for the cells no cell-type rule matches, the donor match key,
and for some datasets the metadata sources to compare.

`pipeline qc-scan SPEC`: one pass over the whole count matrix (every cell, no
filters) -> SPEC's output_dir/restricted/ (per cell and per donor, so it stays
on the cluster and is only read by the checks):
    scan.h5     per cell: UMIs, genes detected, mitochondrial, XIST and Y-gene
                UMIs, donor, cell type, label, and the authors' values; sums
                per donor x cell type over all genes (pseudobulk); per label,
                marker genes' mean log-expression and detection; per cell
                type, the cells detecting each gene
    donors.tsv  one row per donor: covariates (as in the batch-effects spec),
                recorded sex, match key

`pipeline qc-check CHECK SPEC [SPEC ...]`: aggregate tables (rows are cell
types, genes, bins or labels, never donors or cells) -> --output, by default
the spec's output_dir (one spec) or data/v2/qc/<check> (several):
    totals          our per-cell UMIs, genes and mitochondrial % vs the
                    authors' (exact matches, ratios)
    cells           QC distributions per cell type and label (quantiles,
                    histograms), donor-level outliers
    sex             sex from XIST and Y genes per donor vs the recorded sex;
                    per donor, the share of cells with the other sex's
                    transcripts (pooled-library contamination)
    markers         canonical markers per cell type and label; composition
    sex-de          per-gene sex effect per cell type (pseudobulk, + age)
    ad-de           per-gene AD effect per cell type (pathology, and the
                    disease label), against the Mathys 2019 DEGs
    metadata        agreement between metadata sources in the file (e.g.
                    Rush, AMP-AD and ROSMAP Braak stages); pools listed in
                    the spec (e.g. without ATAC fragments)
    sex-de-compare  (several specs) concordance of the sex effects
    identity        (several specs) donors in two datasets: agreement of
                    their metadata, and whether a donor's expression matches
                    its own profile in the other dataset best
    genotype        (one or more specs) whether a donor's expression at
                    strong cis-eQTLs matches the genotype of the person it is
                    labelled as, scored against every genotyped sample (so a
                    swapped sample is also identified); -> data/v2/qc/genotype

Individual-level lists for follow-up (donor IDs of sex or identity
mismatches) go to <output>/restricted/ and are never logged.
"""

import json
import time
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy import stats

from ..core import batch_metrics as bm
from ..utils.counts import count_non_integer_values, log_normalize, total_counts
from ..utils.logging import get_logger
from ..utils.obs_summary import obs_column_names, read_obs_column
from ..utils.provenance import provenance_record
from . import batch_effects as be

try:
    from anndata.io import read_elem, sparse_dataset
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem, sparse_dataset

logger = get_logger(__name__)

OTHER = "Other"
XIST = "XIST"
# Y genes outside the pseudoautosomal regions that are broadly expressed
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY", "NLGN4Y", "TXLNGY"]
# Known sex-biased genes: XIST, X genes that escape X inactivation, Y genes
SEX_GENES = [XIST, "TSIX", "JPX", "KDM5C", "KDM6A", "DDX3X", "EIF1AX", "ZFX", "USP9X", "TXLNG",
             "NLGN4X", "PRKX"] + Y_GENES
MARKERS = {
    "Exc": ["SLC17A7", "SATB2", "CAMK2A", "NRGN"],
    "Inh": ["GAD1", "GAD2", "SLC32A1", "LHX6"],
    "Neuron": ["SNAP25", "SYT1", "RBFOX3"],
    "Ast": ["AQP4", "GFAP", "SLC1A2", "ALDH1L1"],
    "Oli": ["MBP", "PLP1", "MOBP", "ST18"],
    "OPC": ["PDGFRA", "VCAN", "CSPG4"],
    "Mic": ["CSF1R", "P2RY12", "CX3CR1", "C3"],
    "Vasc": ["CLDN5", "FLT1", "PDGFRB", "DCN"],
    "Immune": ["PTPRC", "CD3E", "CD247", "MS4A1", "F13A1", "MRC1"],
}
MATHYS2019_CELLTYPES = {"Ex": "Exc", "In": "Inh", "Ast": "Ast", "Oli": "Oli", "Opc": "OPC", "Mic": "Mic"}
QUANTILES = (0.05, 0.25, 0.5, 0.75, 0.95)
CHECKS = ("totals", "cells", "sex", "markers", "sex-de", "ad-de", "metadata", "sex-de-compare", "identity",
          "genotype")
MULTI_CHECKS = ("sex-de-compare", "identity")

QC_DEFAULTS: Dict[str, Any] = {
    "dataset": None,            # the dataset of the batch-effects spec (default: its only one)
    "celltype": None,           # cell-type rules replacing the batch-effects ones
    "match_on": None,           # obs column that identifies a donor across datasets
    "covariates": {},           # more covariates (batch-effects format), e.g. cogdx
    "author_label": None,       # the authors' label for cells no rule matches
    "author_qc": {},            # umis / genes / mito_pct / mito_fraction: obs column
    "author_qc_table": None,    # the same from a CSV: {path, index: [cols], sep, columns}
    "sex_recode": None,         # raw sex value -> M / F (default: male/female, 1/0)
    "metadata_sources": {},     # key: [sources] (see check_metadata)
    "pool_column": None,
    "pools": {},                # name: [pool values], counted by check_metadata
    "chunk_size": 50000,
    "target_sum": 1e4,
    "min_cells": 20,            # per donor x cell type, for the pseudobulk checks
    "min_cells_identity": 30,
    "identity_genes": 2000,
    "min_donors_per_group": 5,
    "mt_genes": "data/v2/reference/mt_gene_ids.txt",
    "x_genes": "data/v2/reference/x_gene_ids.txt",
    "y_genes": "data/v2/reference/y_gene_ids.txt",
}
QC_REQUIRED = ("output_dir", "batch_effects")


# ---------------------------------------------------------------- spec


def load_qc_spec(path: Union[str, Path]) -> Dict[str, Any]:
    """The QC spec with defaults, plus 'cfg': the batch-effects dataset block
    (paths relative to the repo root; `batch_effects` relative to this file)."""
    with open(path) as f:
        raw = yaml.safe_load(f) or {}
    unknown = set(raw) - set(QC_DEFAULTS) - set(QC_REQUIRED) - {"name"}
    if unknown:
        raise ValueError(f"Unknown key(s) in {path}: {sorted(unknown)}")
    missing = [k for k in QC_REQUIRED if not raw.get(k)]
    if missing:
        raise ValueError(f"{path} needs {missing}")
    be_spec = be.load_spec(Path(path).parent / raw["batch_effects"])
    name = raw.get("dataset") or next(iter(be_spec["datasets"]))
    if name not in be_spec["datasets"]:
        raise ValueError(f"{raw['batch_effects']} has no dataset {name}")
    if not raw.get("dataset") and len(be_spec["datasets"]) > 1:
        raise ValueError(f"{raw['batch_effects']} has several datasets; name one with 'dataset'")
    cfg = dict(be_spec["datasets"][name])
    if raw.get("celltype"):
        cfg["celltype"] = raw["celltype"]
    if raw.get("match_on"):
        cfg["match_on"] = raw["match_on"]
    cfg["covariates"] = {**cfg.get("covariates", {}), **(raw.get("covariates") or {})}
    spec = {**QC_DEFAULTS, "name": Path(path).stem, **raw, "dataset": name, "cfg": cfg}
    if "sex" not in cfg["covariates"]:
        raise ValueError(f"{path}: the dataset needs a 'sex' covariate")
    for key in spec["author_qc"]:
        if key not in ("umis", "genes", "mito_pct", "mito_fraction"):
            raise ValueError(f"{path}: author_qc key {key!r} is not umis, genes, mito_pct or mito_fraction")
    return spec


def _restricted(directory: Path) -> Path:
    path = Path(directory) / "restricted"
    path.mkdir(parents=True, exist_ok=True)
    readme = path / "README.txt"
    if not readme.exists():
        readme.write_text("Per-cell and per-donor data (donor IDs). Input of the QC checks and for "
                          "follow-up on the cluster only; not for sharing.\n")
    return path


def _write_json(path: Path, step: str, params: Dict[str, Any], extra: Optional[Dict] = None) -> None:
    record = provenance_record(step, params)
    if extra:
        record.update(extra)
    path.write_text(json.dumps(record, indent=2, default=str) + "\n")


# ---------------------------------------------------------------- labels


def normalize_sex(values, recode: Optional[Dict] = None) -> np.ndarray:
    """'M' / 'F' / NaN. Without a recode map: male/female (any case) and the
    ROSMAP msex codes 1 (male) / 0 (female)."""
    labels = be._as_labels(values, recode)
    default = {"m": "M", "male": "M", "man": "M", "1": "M", "f": "F", "female": "F", "woman": "F", "0": "F"}
    out = np.array([default.get(str(v).strip().lower(), np.nan) if not pd.isna(v) else np.nan
                    for v in labels], dtype=object)
    unmapped = int((~pd.isna(labels) & pd.isna(out)).sum())
    if unmapped:
        logger.warning(f"{unmapped} sex values are neither male nor female (left missing)")
    return out


def _read_ids(path: Optional[Union[str, Path]]) -> set:
    if not path or not Path(path).exists():
        if path:
            logger.warning(f"{path} not found")
        return set()
    return {line.strip().split(".")[0] for line in Path(path).read_text().splitlines() if line.strip()}


def gene_table(var_names: Sequence[str], gene_map: Optional[Union[str, Path]], spec: Dict[str, Any]) -> pd.DataFrame:
    """Per file column: var_name, gene_id (Ensembl, if mapped), symbol (the
    reference's name, else var_name) and chrom: MT, X (with the PAR genes), Y
    (Y only), auto(somal), or unmapped."""
    genes = pd.DataFrame({"var_name": [str(v) for v in var_names]})
    if gene_map:
        gm = pd.read_csv(gene_map, sep="\t", dtype=str).drop_duplicates("var_name")
        genes = genes.merge(gm[["var_name", "gene_id", "ref_name"]], on="var_name", how="left")
    else:
        genes["gene_id"], genes["ref_name"] = np.nan, np.nan
    genes["symbol"] = genes["ref_name"].fillna(genes["var_name"]).astype(str)
    mt, x, y = _read_ids(spec["mt_genes"]), _read_ids(spec["x_genes"]), _read_ids(spec["y_genes"])
    ids = genes["gene_id"]
    chrom = np.where(ids.isna(), "unmapped", "auto").astype(object)
    chrom[ids.isin(y).to_numpy()] = "Y"
    chrom[ids.isin(x).to_numpy()] = "X"
    is_mt = ids.isin(mt).to_numpy() if mt else genes["symbol"].str.upper().str.startswith("MT-").to_numpy()
    chrom[is_mt] = "MT"
    genes["chrom"] = chrom
    return genes.drop(columns="ref_name")


def _first_column(genes: pd.DataFrame, symbol: str) -> Optional[int]:
    hits = np.flatnonzero(genes["symbol"].str.upper().to_numpy() == symbol.upper())
    return int(hits[0]) if len(hits) else None


# ---------------------------------------------------------------- scan


def _author_values(spec: Dict[str, Any], raw: Dict[str, Any], obs_names) -> Dict[str, np.ndarray]:
    """The authors' per-cell QC values: umis, genes, mito_pct (as % of UMIs)."""
    out = {}
    for key, column in spec["author_qc"].items():
        out[key] = pd.to_numeric(pd.Series(np.asarray(raw[column])), errors="coerce").to_numpy(np.float64)
    table = spec["author_qc_table"]
    if table:
        index, columns = list(table["index"]), dict(table["columns"])
        logger.info(f"reading the authors' per-cell QC from {table['path']}")
        df = pd.read_csv(table["path"], usecols=index + list(columns.values()), low_memory=False)
        key = df[index[0]].astype(str)
        for c in index[1:]:
            key = key + table.get("sep", "_") + df[c].astype(str)
        df.index = pd.Index(key)
        duplicated = int(df.index.duplicated().sum())
        if duplicated:
            logger.warning(f"{duplicated} rows of the QC table have a duplicate cell name (first kept)")
            df = df[~df.index.duplicated()]
        aligned = df.reindex(pd.Index([str(n) for n in obs_names]))
        logger.info(f"{int(aligned.notna().any(axis=1).sum())} of {len(aligned)} cells are in the QC table")
        for k, column in columns.items():
            out[k] = pd.to_numeric(aligned[column], errors="coerce").to_numpy(np.float64)
    if "mito_fraction" in out:
        out["mito_pct"] = 100 * out.pop("mito_fraction")
    return out


def _donor_table(cfg: Dict[str, Any], spec: Dict[str, Any], raw: Dict[str, Any], donors: np.ndarray,
                 donor_ids: Sequence[str]) -> pd.DataFrame:
    """One row per donor (in donor-code order): n_cells, covariates (the most
    common value), recorded sex, match key."""
    cells = pd.DataFrame({"donor": donors})
    for key, cov in cfg["covariates"].items():
        cells[key] = be._as_numeric(raw[cov["column"]], cov.get("recode")) if cov.get("numeric") else \
            be._as_labels(raw[cov["column"]], cov.get("recode"))
    cells["sex_recorded"] = normalize_sex(raw[cfg["covariates"]["sex"]["column"]], spec["sex_recode"])
    if cfg.get("disease"):
        cells["disease"] = be._as_labels(raw[cfg["disease"]["column"]], cfg["disease"].get("recode"))
    if cfg.get("match_on"):
        cells["match"] = be._as_labels(raw[cfg["match_on"]])
    cells = cells[cells["donor"].notna()]
    table = pd.DataFrame(index=pd.Index(donor_ids, name="donor"))
    table["n_cells"] = cells.groupby("donor").size()
    for column in cells.columns.drop("donor"):
        table[column] = be._donor_majority(cells, column)
    return table.reset_index()


def run_scan(spec_path: Union[str, Path], overwrite: bool = False) -> Path:
    spec = load_qc_spec(spec_path)
    cfg = spec["cfg"]
    restricted = _restricted(Path(spec["output_dir"]))
    scan_path = restricted / "scan.h5"
    if scan_path.exists() and not overwrite:
        logger.info(f"{scan_path} exists (use --overwrite to redo)")
        return scan_path
    path = be._first_existing(cfg["counts"])
    columns = {cfg["donor"], *(r["column"] for r in cfg["celltype"]),
               *(c["column"] for c in cfg["covariates"].values()), *spec["author_qc"].values()}
    for extra in (cfg.get("match_on"), spec["author_label"], (cfg.get("disease") or {}).get("column")):
        if extra:
            columns.add(extra)

    logger.info(f"[{spec['name']}] reading obs from {path}")
    with h5py.File(path, "r") as f:
        var = f["var"]
        var_names = [str(v) for v in read_elem(var[var.attrs["_index"]])]
        n_obs = int(be._matrix_node(f, cfg.get("layer")).attrs["shape"][0])
        missing = sorted(columns - set(obs_column_names(f)))
        if missing:
            raise ValueError(f"obs has no column(s) {missing}")
        raw = {c: read_obs_column(f, c) for c in columns}
        obs_names = read_elem(f["obs"][f["obs"].attrs["_index"]]) if spec["author_qc_table"] else None
    author = _author_values(spec, raw, obs_names)

    donors = be._as_labels(raw[cfg["donor"]])
    donor_codes, donor_ids = pd.factorize(pd.Series(donors), use_na_sentinel=True)
    donor_ids = [str(d) for d in donor_ids]
    celltype = be._assign_celltypes(raw, cfg["celltype"], n_obs)
    celltype = np.where(pd.isna(celltype), OTHER, celltype).astype(object)
    ct_levels = list(dict.fromkeys(r["label"] for r in cfg["celltype"])) + [OTHER]
    ct_codes = pd.Categorical(celltype, categories=ct_levels).codes.astype(np.int64)
    if spec["author_label"]:
        author_label = pd.Series(be._as_labels(raw[spec["author_label"]])).fillna("unlabelled").astype(str)
        label = np.where(celltype == OTHER, ("Other: " + author_label).to_numpy(), celltype)
    else:
        label = celltype
    label_codes, label_levels = pd.factorize(pd.Series(label), sort=True)
    group_codes = np.full(n_obs, -1, dtype=np.int64)
    valid = donor_codes >= 0
    codes, groups = pd.MultiIndex.from_arrays([donor_codes[valid], ct_codes[valid]]).factorize()
    group_codes[valid] = codes
    group_donor = np.asarray(groups.get_level_values(0), dtype=np.int64)
    group_ct = np.asarray(groups.get_level_values(1), dtype=np.int64)
    logger.info(f"[{spec['name']}] {n_obs} cells, {len(donor_ids)} donors, {int((~valid).sum())} cells "
                f"without a donor, {len(groups)} donor x cell type groups")
    for t in ct_levels:
        logger.info(f"[{spec['name']}]   {t}: {int((celltype == t).sum())} cells")

    genes = gene_table(var_names, cfg.get("gene_map"), spec)
    mt_cols = np.flatnonzero(genes["chrom"].to_numpy() == "MT")
    xist_col = _first_column(genes, XIST)
    y_cols = [c for c in (_first_column(genes, g) for g in Y_GENES) if c is not None]
    marker_symbols = list(dict.fromkeys(g for gs in MARKERS.values() for g in gs))
    marker_cols = {g: _first_column(genes, g) for g in marker_symbols}
    marker_symbols = [g for g in marker_symbols if marker_cols[g] is not None]
    marker_idx = [marker_cols[g] for g in marker_symbols]
    logger.info(f"[{spec['name']}] {len(var_names)} genes ({int(genes['gene_id'].notna().sum())} mapped), "
                f"{len(mt_cols)} mitochondrial, XIST {'found' if xist_col is not None else 'MISSING'}, "
                f"{len(y_cols)} of {len(Y_GENES)} Y genes, {len(marker_symbols)} of "
                f"{sum(len(g) for g in MARKERS.values())} markers")

    n_genes, n_groups, n_ct, n_lab = len(var_names), len(groups), len(ct_levels), len(label_levels)
    pseudobulk = np.zeros((n_groups, n_genes))
    library = np.zeros(n_groups)
    detected = np.zeros((n_ct, n_genes))
    marker_sum = np.zeros((n_lab, len(marker_idx)))
    marker_det = np.zeros((n_lab, len(marker_idx)))
    umis = np.zeros(n_obs)
    ngenes = np.zeros(n_obs, dtype=np.int32)
    mito = np.zeros(n_obs)
    xist = np.full(n_obs, np.nan)
    ycount = np.zeros(n_obs)
    start_time = time.time()
    chunk = int(spec["chunk_size"])
    with h5py.File(path, "r") as f:
        matrix = sparse_dataset(be._matrix_node(f, cfg.get("layer")))
        for start in range(0, n_obs, chunk):
            stop = min(start + chunk, n_obs)
            block = sp.csr_matrix(matrix[start:stop])
            bad = count_non_integer_values(block)
            if bad:
                raise ValueError(f"rows {start}-{stop}: {bad} values are not counts")
            tot = total_counts(block)
            positive = block.copy()
            positive.data = (positive.data > 0).astype(np.float64)
            positive.eliminate_zeros()
            umis[start:stop] = tot
            ngenes[start:stop] = np.diff(positive.indptr)
            if len(mt_cols):
                mito[start:stop] = np.asarray(block[:, mt_cols].sum(axis=1)).ravel()
            if xist_col is not None:
                xist[start:stop] = block[:, xist_col].toarray().ravel()
            if y_cols:
                ycount[start:stop] = np.asarray(block[:, y_cols].sum(axis=1)).ravel()
            g = group_codes[start:stop]
            rows = np.flatnonzero(g >= 0)
            if len(rows):
                local, inverse = np.unique(g[rows], return_inverse=True)
                indicator = bm._indicator(inverse, len(local))
                pseudobulk[local] += (indicator @ block[rows]).toarray()
                library[local] += indicator @ tot[rows]
            detected += (bm._indicator(ct_codes[start:stop], n_ct) @ positive).toarray()
            if marker_idx:
                m = block[:, marker_idx]
                indicator = bm._indicator(label_codes[start:stop], n_lab)
                marker_sum += (indicator @ log_normalize(m, tot, spec["target_sum"], dtype=np.float64)).toarray()
                m_pos = m.copy()
                m_pos.data = (m_pos.data > 0).astype(np.float64)
                marker_det += (indicator @ m_pos).toarray()
            if (start // chunk) % 10 == 0:
                logger.info(f"[{spec['name']}] read {stop}/{n_obs} rows ({time.time() - start_time:.0f} s)")
    logger.info(f"[{spec['name']}] read all {n_obs} rows in {time.time() - start_time:.0f} s")

    donor_table = _donor_table(cfg, spec, raw, donors, donor_ids)
    tmp = scan_path.with_suffix(".h5.tmp")
    str_type = h5py.string_dtype()
    gz = dict(compression="gzip", compression_opts=1)
    with h5py.File(tmp, "w") as h:
        c = h.create_group("cells")
        for key, values in {"umis": umis, "genes": ngenes, "mito": mito, "xist": xist, "y": ycount,
                            "donor": donor_codes.astype(np.int32), "celltype": ct_codes.astype(np.int16),
                            "label": label_codes.astype(np.int32)}.items():
            c.create_dataset(key, data=values, **gz)
        for key, values in author.items():
            c.create_dataset(f"author_{key}", data=values, **gz)
        h.create_dataset("donors", data=np.array(donor_ids, dtype=object), dtype=str_type)
        h.create_dataset("celltypes", data=np.array(ct_levels, dtype=object), dtype=str_type)
        h.create_dataset("labels", data=np.array([str(x) for x in label_levels], dtype=object), dtype=str_type)
        g = h.create_group("groups")
        g.create_dataset("donor", data=group_donor)
        g.create_dataset("celltype", data=group_ct)
        g.create_dataset("n_cells", data=np.bincount(group_codes[valid], minlength=n_groups))
        g.create_dataset("library", data=library)
        h.create_dataset("pseudobulk", data=pseudobulk, chunks=(max(1, min(64, n_groups)), n_genes), **gz)
        gg = h.create_group("genes")
        for column in ("var_name", "gene_id", "symbol", "chrom"):
            gg.create_dataset(column, data=genes[column].fillna("").astype(str).to_numpy(dtype=object),
                              dtype=str_type)
        h.create_dataset("detected", data=detected)
        h.create_dataset("ct_n_cells", data=np.bincount(ct_codes, minlength=n_ct))
        mk = h.create_group("markers")
        mk.create_dataset("genes", data=np.array(marker_symbols, dtype=object), dtype=str_type)
        mk.create_dataset("sum", data=marker_sum)
        mk.create_dataset("detected", data=marker_det)
        mk.create_dataset("label_n_cells", data=np.bincount(label_codes, minlength=n_lab))
        h.attrs["dataset"] = spec["name"]
        h.attrs["counts"] = str(path)
        h.attrs["provenance"] = json.dumps(provenance_record("qc-scan", {"spec": str(spec_path)}), default=str)
    tmp.replace(scan_path)
    donor_table.to_csv(restricted / "donors.tsv", sep="\t", index=False)
    logger.info(f"[{spec['name']}] wrote {scan_path} and {restricted / 'donors.tsv'}")
    return scan_path


# ---------------------------------------------------------------- loading


class Scan:
    """A scan.h5 opened for the checks (cells and small tables in memory; the
    pseudobulk is read by row on demand)."""

    def __init__(self, spec: Dict[str, Any]):
        self.spec = spec
        self.name = spec["name"]
        restricted = Path(spec["output_dir"]) / "restricted"
        self.path = restricted / "scan.h5"
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} missing: run `pipeline qc-scan` first")
        with h5py.File(self.path, "r") as h:
            self.cells = {k: h["cells"][k][()] for k in h["cells"]}
            self.celltypes = list(h["celltypes"].asstr()[()])
            self.labels = list(h["labels"].asstr()[()])
            self.groups = pd.DataFrame({k: h["groups"][k][()] for k in h["groups"]})
            self.genes = pd.DataFrame({k: h["genes"][k].asstr()[()] for k in h["genes"]})
            self.detected = h["detected"][()]
            self.ct_n_cells = h["ct_n_cells"][()]
            self.marker_genes = list(h["markers"]["genes"].asstr()[()])
            self.marker_sum = h["markers"]["sum"][()]
            self.marker_det = h["markers"]["detected"][()]
            self.label_n_cells = h["markers"]["label_n_cells"][()]
            self.n_donor_ids = len(h["donors"])
        self.genes["gene_id"] = self.genes["gene_id"].replace("", np.nan)
        self.groups["celltype_name"] = [self.celltypes[c] for c in self.groups["celltype"]]
        self.donors = pd.read_csv(restricted / "donors.tsv", sep="\t", dtype={"donor": str, "match": str})
        if len(self.donors) != self.n_donor_ids:
            raise ValueError(f"{restricted / 'donors.tsv'} does not match {self.path}")

    def author(self, key: str) -> Optional[np.ndarray]:
        return self.cells.get(f"author_{key}")

    def pseudobulk_rows(self, rows: np.ndarray, columns: Optional[np.ndarray] = None) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        order = np.argsort(rows)
        with h5py.File(self.path, "r") as h:
            data = h["pseudobulk"][rows[order]] if len(rows) else np.zeros((0, len(self.genes)))
        out = np.empty_like(data)
        out[order] = data
        return out if columns is None else out[:, columns]

    def pseudobulk_columns(self, columns: Sequence[int]) -> np.ndarray:
        with h5py.File(self.path, "r") as h:
            ds = h["pseudobulk"]
            out = np.zeros((ds.shape[0], len(columns)))
            step = ds.chunks[0] if ds.chunks else 1024
            cols = np.asarray(columns)
            for start in range(0, ds.shape[0], step * 16):
                out[start:start + step * 16] = ds[start:start + step * 16][:, cols]
        return out


def _quantile_row(values: np.ndarray, prefix: str) -> Dict[str, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return {f"{prefix}_p{int(q * 100):02d}": np.nan for q in QUANTILES}
    qs = np.quantile(values, QUANTILES)
    return {f"{prefix}_p{int(q * 100):02d}": float(v) for q, v in zip(QUANTILES, qs)}


def _otsu(values: np.ndarray) -> float:
    """Threshold between two groups that maximizes the between-group variance."""
    v = np.sort(values[np.isfinite(values)])
    if len(v) < 2 or v[0] == v[-1]:
        return float(np.nan)
    n = len(v)
    csum = np.cumsum(v)
    k = np.arange(1, n)
    w0, w1 = k / n, 1 - k / n
    m0 = csum[:-1] / k
    m1 = (csum[-1] - csum[:-1]) / (n - k)
    best = int(np.argmax(w0 * w1 * (m0 - m1) ** 2))
    return float((v[best] + v[best + 1]) / 2)


def _logcpm(counts: np.ndarray, library: np.ndarray) -> np.ndarray:
    return np.log2((counts + 0.5) / (library[:, None] + 1.0) * 1e6)


def _ols_term(Y: np.ndarray, term: np.ndarray, adjust: pd.DataFrame, numeric: Sequence[str]):
    """Per column of Y: coefficient, t and two-sided p of `term` (0/1 or
    numeric) in Y ~ 1 + term + adjust."""
    blocks = [np.ones((len(term), 1)), np.asarray(term, dtype=np.float64)[:, None]]
    if adjust.shape[1]:
        blocks.append(bm.design_matrix(adjust, numeric)[0][:, 1:])
    X = np.hstack(blocks)
    n, p = X.shape
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ (X.T @ Y)
    df = max(n - np.linalg.matrix_rank(X), 1)
    rss = ((Y - X @ beta) ** 2).sum(axis=0)
    se = np.sqrt(rss / df * xtx_inv[1, 1])
    t = np.divide(beta[1], se, out=np.full(Y.shape[1], np.nan), where=se > 0)
    pval = np.where(np.isfinite(t), 2 * stats.t.sf(np.abs(np.nan_to_num(t)), df), np.nan)
    return beta[1], t, pval


def _fdr(p: np.ndarray) -> np.ndarray:
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    out[ok] = bm.bh_fdr(p[ok])
    return out


def _celltype_matrix(scan: Scan, celltype: str, min_cells: int, donor_mask: Optional[np.ndarray] = None):
    """Pseudobulk of one cell type: (donor codes, counts, library) for the
    donors with at least min_cells cells (and in donor_mask)."""
    g = scan.groups
    keep = (g["celltype_name"] == celltype) & (g["n_cells"] >= min_cells)
    if donor_mask is not None:
        keep &= donor_mask[g["donor"].to_numpy()]
    rows = np.flatnonzero(keep.to_numpy())
    return g["donor"].to_numpy()[rows], scan.pseudobulk_rows(rows), g["library"].to_numpy()[rows]


def _analysed_celltypes(scan: Scan) -> List[str]:
    return [t for t in scan.celltypes if t != OTHER]


# ---------------------------------------------------------------- checks: one dataset


def check_totals(scan: Scan) -> Dict[str, pd.DataFrame]:
    c = scan.cells
    with np.errstate(divide="ignore", invalid="ignore"):
        ours = {"umis": c["umis"], "genes": c["genes"].astype(np.float64),
                "mito_pct": np.where(c["umis"] > 0, 100 * c["mito"] / c["umis"], np.nan)}
    tolerance = {"umis": 0.5, "genes": 0.5, "mito_pct": 0.01}
    rows = []
    for key in ("umis", "genes", "mito_pct"):
        theirs = scan.author(key)
        if theirs is None:
            continue
        for t in ["all"] + scan.celltypes:
            mask = np.ones(len(theirs), bool) if t == "all" else c["celltype"] == scan.celltypes.index(t)
            ok = mask & np.isfinite(theirs) & np.isfinite(ours[key])
            a, b = ours[key][ok], theirs[ok]
            d = a - b
            row = {"dataset": scan.name, "field": key, "celltype": t, "n_cells": int(mask.sum()),
                   "n_compared": int(ok.sum())}
            if len(d):
                row.update({
                    "pct_exact": 100 * float(np.mean(np.abs(d) <= tolerance[key])),
                    "pct_within_1pct": 100 * float(np.mean(np.abs(d) <= np.maximum(tolerance[key], 0.01 * np.abs(b)))),
                    "median_ours": float(np.median(a)), "median_authors": float(np.median(b)),
                    "median_ratio": float(np.median(a[b > 0] / b[b > 0])) if (b > 0).any() else np.nan,
                    "pearson_r": float(np.corrcoef(a, b)[0, 1]) if len(d) > 2 and a.std() > 0 and b.std() > 0 else np.nan,
                    "diff_p01": float(np.quantile(d, 0.01)), "diff_p50": float(np.median(d)),
                    "diff_p99": float(np.quantile(d, 0.99)), "max_abs_diff": float(np.abs(d).max()),
                })
            rows.append(row)
    if not rows:
        logger.info(f"[{scan.name}] the authors give no per-cell UMIs, genes or mitochondrial %")
        rows.append({"dataset": scan.name, "field": "none", "celltype": "all", "n_cells": len(c["umis"]),
                     "n_compared": 0})
    table = pd.DataFrame(rows)
    _log_table(scan.name, "totals", table[table["celltype"] == "all"])
    return {"totals": table}


def _log_table(name: str, title: str, table: pd.DataFrame) -> None:
    with pd.option_context("display.max_rows", 200, "display.width", 250, "display.max_columns", 40):
        logger.info(f"[{name}] {title}\n{table.to_string(index=False)}")


# Smallest spread used for the donor outlier rule (> 3 MAD from the median
# donor), so a metric that is nearly constant (e.g. mitochondrial % in nuclei
# preps with almost none) does not flag donors over trivial differences
MAD_FLOOR = {"median_log10_umis": 0.05, "median_genes": 50.0, "median_mito_pct": 0.5, "neuron_share": 0.02}

HISTOGRAMS = {
    "log10_umis": np.arange(1.5, 5.75, 0.05),
    "log10_genes": np.arange(1.5, 4.55, 0.05),
    "mito_pct": np.arange(0, 20.25, 0.25),
}


def check_cells(scan: Scan, min_genes: int = 200) -> Dict[str, pd.DataFrame]:
    c = scan.cells
    with np.errstate(divide="ignore", invalid="ignore"):
        mito_pct = np.where(c["umis"] > 0, 100 * c["mito"] / c["umis"], np.nan)
        metrics = {"umis": c["umis"], "genes": c["genes"].astype(np.float64), "mito_pct": mito_pct}
        hist_values = {"log10_umis": np.log10(np.maximum(c["umis"], 1)),
                       "log10_genes": np.log10(np.maximum(c["genes"], 1)), "mito_pct": mito_pct}

    def summary(mask, kind, name):
        row = {"dataset": scan.name, "level": kind, "group": name, "n_cells": int(mask.sum()),
               "n_donors": int(len(np.unique(c["donor"][mask & (c["donor"] >= 0)])))}
        for key, values in metrics.items():
            row.update(_quantile_row(values[mask], key))
        row["pct_below_min_genes"] = 100 * float(np.mean(c["genes"][mask] < min_genes)) if mask.any() else np.nan
        for cut in (5, 10):
            row[f"pct_mito_above_{cut}"] = 100 * float(np.nanmean(mito_pct[mask] > cut)) if mask.any() else np.nan
        return row

    rows = [summary(np.ones(len(c["umis"]), bool), "all", "all")]
    rows += [summary(c["celltype"] == i, "celltype", t) for i, t in enumerate(scan.celltypes)]
    rows += [summary(c["label"] == i, "label", t) for i, t in enumerate(scan.labels) if t.startswith(OTHER)]
    table = pd.DataFrame(rows)
    _log_table(scan.name, "per-cell QC", table[table["level"] != "label"][
        ["group", "n_cells", "n_donors", "umis_p50", "genes_p50", "mito_pct_p50", "pct_below_min_genes",
         "pct_mito_above_5"]])

    hist = []
    for i, t in enumerate(["all"] + scan.celltypes):
        mask = np.ones(len(c["umis"]), bool) if t == "all" else c["celltype"] == i - 1
        for key, edges in HISTOGRAMS.items():
            v = np.clip(hist_values[key][mask], edges[0], edges[-1] - 1e-9)
            counts, _ = np.histogram(v[np.isfinite(v)], bins=edges)
            hist += [{"dataset": scan.name, "celltype": t, "metric": key, "bin_left": float(lo),
                      "bin_right": float(hi), "n_cells": int(n)} for lo, hi, n in zip(edges[:-1], edges[1:], counts)]

    # Donors: cells, medians, neuron share; outliers > 3 MAD from the median donor
    has = c["donor"] >= 0
    frame = pd.DataFrame({"donor": c["donor"][has], "log10_umis": hist_values["log10_umis"][has],
                          "genes": metrics["genes"][has], "mito_pct": mito_pct[has],
                          "neuron": np.isin(c["celltype"][has], [scan.celltypes.index(t) for t in ("Exc", "Inh")
                                                                 if t in scan.celltypes])})
    per_donor = frame.groupby("donor").agg(n_cells=("donor", "size"), median_log10_umis=("log10_umis", "median"),
                                           median_genes=("genes", "median"), median_mito_pct=("mito_pct", "median"),
                                           neuron_share=("neuron", "mean"))
    donor_rows = []
    for metric in per_donor.columns:
        v = per_donor[metric].to_numpy(np.float64)
        v = v[np.isfinite(v)]
        med = float(np.median(v))
        mad = max(1.4826 * float(np.median(np.abs(v - med))), MAD_FLOOR.get(metric, 0.0))
        row = {"dataset": scan.name, "metric": metric, "n_donors": len(v), **_quantile_row(v, "value"),
               "mad": mad, "n_low_outliers": int(np.sum(v < med - 3 * mad)) if mad > 0 else 0,
               "n_high_outliers": int(np.sum(v > med + 3 * mad)) if mad > 0 else 0}
        if metric == "n_cells":
            for cut in (50, 300, 1000):
                row[f"n_below_{cut}"] = int(np.sum(v < cut))
        donor_rows.append(row)
    any_outlier = np.zeros(len(per_donor), bool)
    for metric in ("median_log10_umis", "median_genes", "median_mito_pct", "neuron_share"):
        v = per_donor[metric].to_numpy(np.float64)
        med = np.nanmedian(v)
        mad = max(1.4826 * np.nanmedian(np.abs(v - med)), MAD_FLOOR.get(metric, 0.0))
        if mad > 0:
            any_outlier |= np.abs(v - med) > 3 * mad
    donor_rows.append({"dataset": scan.name, "metric": "any_outlier", "n_donors": len(per_donor),
                       "n_high_outliers": int(any_outlier.sum())})
    donors = pd.DataFrame(donor_rows)
    _log_table(scan.name, "donors", donors)
    return {"qc_cells": table, "qc_histograms": pd.DataFrame(hist), "qc_donors": donors}


def check_sex(scan: Scan) -> Dict[str, pd.DataFrame]:
    genes = scan.genes
    xist_col = _first_column(genes, XIST)
    y_cols = [c for c in (_first_column(genes, g) for g in Y_GENES) if c is not None]
    if xist_col is None or not y_cols:
        raise ValueError(f"[{scan.name}] XIST or the Y genes are missing")
    g = scan.groups
    n = scan.n_donor_ids
    sums = scan.pseudobulk_columns([xist_col] + y_cols)
    donor = g["donor"].to_numpy()
    xist = np.bincount(donor, weights=sums[:, 0], minlength=n)
    ysum = np.bincount(donor, weights=sums[:, 1:].sum(axis=1), minlength=n)
    library = np.bincount(donor, weights=g["library"].to_numpy(), minlength=n)
    with np.errstate(divide="ignore", invalid="ignore"):
        x_score = np.log2(1 + xist / library * 1e6)
        y_score = np.log2(1 + ysum / library * 1e6)
    tx, ty = _otsu(x_score), _otsu(y_score)
    called = np.where((y_score > ty) & (x_score <= tx), "M",
                      np.where((x_score > tx) & (y_score <= ty), "F", "ambiguous")).astype(object)
    recorded = scan.donors["sex_recorded"].fillna("missing").to_numpy(dtype=object)
    logger.info(f"[{scan.name}] thresholds: log2(1 + XIST CPM) {tx:.2f}, log2(1 + Y CPM) {ty:.2f}")
    calls = pd.crosstab(pd.Series(recorded, name="recorded"), pd.Series(called, name="expression"))
    calls = calls.reset_index().melt(id_vars="recorded", var_name="expression", value_name="n_donors")
    calls.insert(0, "dataset", scan.name)
    discordant = (recorded != "missing") & (called != "ambiguous") & (recorded != called)
    summary = {"dataset": scan.name, "n_donors": n, "n_recorded": int((recorded != "missing").sum()),
               "n_male_expr": int((called == "M").sum()), "n_female_expr": int((called == "F").sum()),
               "n_ambiguous": int((called == "ambiguous").sum()), "n_discordant": int(discordant.sum()),
               "xist_threshold": tx, "y_threshold": ty,
               "xist_score_gap": float(np.median(x_score[called == "F"]) - np.median(x_score[called == "M"]))
               if (called == "F").any() and (called == "M").any() else np.nan,
               "y_score_gap": float(np.median(y_score[called == "M"]) - np.median(y_score[called == "F"]))
               if (called == "F").any() and (called == "M").any() else np.nan}

    # Cell level: share of a donor's cells with the other sex's transcripts
    c = scan.cells
    has = c["donor"] >= 0
    n_cells = np.bincount(c["donor"][has], minlength=n)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_y = np.bincount(c["donor"][has], weights=(c["y"][has] > 0), minlength=n) / n_cells
        frac_xist = np.bincount(c["donor"][has], weights=(np.nan_to_num(c["xist"][has]) > 0), minlength=n) / n_cells
    cell_rows = []
    for sex, values, what in (("F", frac_y, "share of cells with Y-gene UMIs"),
                              ("M", frac_xist, "share of cells with XIST UMIs"),
                              ("M", frac_y, "share of cells with Y-gene UMIs (own sex)"),
                              ("F", frac_xist, "share of cells with XIST UMIs (own sex)")):
        v = values[(called == sex) & (n_cells >= 50)]
        row = {"dataset": scan.name, "donors": f"expression {sex}", "measure": what, "n_donors": len(v),
               **_quantile_row(v, "value"), "max": float(v.max()) if len(v) else np.nan}
        for cut in (0.02, 0.05, 0.1, 0.2):
            row[f"n_above_{cut}"] = int(np.sum(v > cut))
        cell_rows.append(row)

    score = y_score - x_score
    edges = np.arange(-16, 16.5, 0.5)
    hist = []
    for sex in ("M", "F", "missing"):
        v = np.clip(score[(recorded == sex) & np.isfinite(score)], edges[0], edges[-1] - 1e-9)
        counts, _ = np.histogram(v, bins=edges)
        hist += [{"dataset": scan.name, "recorded": sex, "bin_left": float(lo), "bin_right": float(hi),
                  "n_donors": int(k)} for lo, hi, k in zip(edges[:-1], edges[1:], counts)]
    edges = np.arange(0, 1.0001, 0.01)
    for sex, values, what in (("F", frac_y, "y_in_female"), ("M", frac_xist, "xist_in_male")):
        v = values[(called == sex) & (n_cells >= 50)]
        counts, _ = np.histogram(np.clip(v, 0, 1 - 1e-9), bins=edges)
        hist += [{"dataset": scan.name, "recorded": what, "bin_left": float(lo), "bin_right": float(hi),
                  "n_donors": int(k)} for lo, hi, k in zip(edges[:-1], edges[1:], counts)]

    flagged = discordant | (called == "ambiguous")
    by_donor = pd.DataFrame({"donor": scan.donors["donor"], "recorded": recorded, "expression": called,
                             "log2_xist_cpm": x_score, "log2_y_cpm": y_score, "share_cells_y": frac_y,
                             "share_cells_xist": frac_xist, "n_cells": n_cells})
    restricted = _restricted(Path(scan.spec["output_dir"]))
    by_donor[flagged].to_csv(restricted / "sex_follow_up.tsv", sep="\t", index=False)
    # Every donor's call, for the identity check (same donors in other datasets)
    by_donor[["donor", "recorded", "expression"]].to_csv(restricted / "sex_by_donor.tsv", sep="\t", index=False)
    logger.info(f"[{scan.name}] {int(flagged.sum())} donors to follow up -> {restricted / 'sex_follow_up.tsv'}")
    out = {"sex_summary": pd.DataFrame([summary]), "sex_calls": calls, "sex_cell_level": pd.DataFrame(cell_rows),
           "sex_histograms": pd.DataFrame(hist)}
    _log_table(scan.name, "sex", out["sex_summary"])
    _log_table(scan.name, "sex calls", calls)
    _log_table(scan.name, "cell-level sex signal", out["sex_cell_level"])
    return out


def check_markers(scan: Scan) -> Dict[str, pd.DataFrame]:
    rows = []
    marker_of = {g: k for k, gs in MARKERS.items() for g in gs}
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = scan.marker_sum / scan.label_n_cells[:, None]
        frac = scan.marker_det / scan.label_n_cells[:, None]
    for i, label in enumerate(scan.labels):
        celltype = OTHER if label.startswith(OTHER) else label
        for j, gene in enumerate(scan.marker_genes):
            rows.append({"dataset": scan.name, "label": label, "celltype": celltype,
                         "n_cells": int(scan.label_n_cells[i]), "gene": gene, "marker_of": marker_of[gene],
                         "mean_log_norm": float(mean[i, j]), "frac_detected": float(frac[i, j])})
    markers = pd.DataFrame(rows)

    spec_rows = []
    harmonized = [t for t in scan.celltypes if t != OTHER and t in scan.labels]
    for t in harmonized:
        own = [g for g in MARKERS.get(t, []) if g in scan.marker_genes]
        for gene in own:
            j = scan.marker_genes.index(gene)
            values = {u: mean[scan.labels.index(u), j] for u in harmonized}
            others = {u: v for u, v in values.items() if u != t}
            top_other = max(others, key=others.get) if others else None
            spec_rows.append({"dataset": scan.name, "celltype": t, "gene": gene, "mean_own": float(values[t]),
                              "top_other": top_other, "mean_top_other": float(others[top_other]) if top_other else np.nan,
                              "own_is_highest": bool(not others or values[t] > max(others.values())),
                              "frac_own": float(frac[scan.labels.index(t), j])})
    specificity = pd.DataFrame(spec_rows)
    if len(specificity):
        _log_table(scan.name, "marker specificity",
                   specificity.groupby("celltype").agg(n_markers=("gene", "size"),
                                                       n_highest_in_own=("own_is_highest", "sum")).reset_index())

    # Composition: pooled and per donor
    g = scan.groups
    per_donor = g.pivot_table(index="donor", columns="celltype_name", values="n_cells", aggfunc="sum",
                              fill_value=0)
    share = per_donor.div(per_donor.sum(axis=1), axis=0)
    comp = []
    total = scan.ct_n_cells.sum()
    for i, t in enumerate(scan.celltypes):
        v = share[t].to_numpy() if t in share else np.zeros(len(share))
        comp.append({"dataset": scan.name, "celltype": t, "n_cells": int(scan.ct_n_cells[i]),
                     "pct_cells": 100 * float(scan.ct_n_cells[i] / total), "n_donors": int(np.sum(v > 0)),
                     **{k: 100 * x for k, x in _quantile_row(v, "donor_pct").items()}})
    composition = pd.DataFrame(comp)
    _log_table(scan.name, "composition", composition)

    # Cells no rule matched: their authors' label, size, and the marker set they express most
    other_rows = []
    for i, label in enumerate(scan.labels):
        if not label.startswith(OTHER):
            continue
        set_means = {k: float(np.nanmean([mean[i, scan.marker_genes.index(g)] for g in gs if g in scan.marker_genes]))
                     for k, gs in MARKERS.items() if any(g in scan.marker_genes for g in gs)}
        top = max(set_means, key=set_means.get)
        other_rows.append({"dataset": scan.name, "label": label, "n_cells": int(scan.label_n_cells[i]),
                           "top_marker_set": top, "top_marker_set_mean": set_means[top]})
    other = pd.DataFrame(other_rows)
    if len(other):
        _log_table(scan.name, "cells without a harmonized type", other)
    return {"markers": markers, "marker_specificity": specificity, "composition": composition,
            "unmatched_labels": other}


def _expressed(counts: np.ndarray, library: np.ndarray, min_cpm: float = 1.0, min_share: float = 0.2) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cpm = counts / library[:, None] * 1e6
    return (cpm >= min_cpm).mean(axis=0) >= min_share


def check_sex_de(scan: Scan) -> Dict[str, pd.DataFrame]:
    d = scan.donors
    sex = d["sex_recorded"].to_numpy(dtype=object)
    age = pd.to_numeric(d.get("age"), errors="coerce").to_numpy(np.float64) if "age" in d else \
        np.full(len(d), np.nan)
    usable = np.isin(sex, ["M", "F"])
    use_age = np.isfinite(age).mean() > 0.5
    if use_age:
        usable &= np.isfinite(age)
    rows, summary = [], []
    for t in _analysed_celltypes(scan):
        donors, counts, library = _celltype_matrix(scan, t, scan.spec["min_cells"], usable)
        n_m, n_f = int(np.sum(sex[donors] == "M")), int(np.sum(sex[donors] == "F"))
        if min(n_m, n_f) < scan.spec["min_donors_per_group"]:
            logger.warning(f"[{scan.name}] {t}: {n_m} male and {n_f} female donors, left out")
            continue
        keep = _expressed(counts, library)
        Y = _logcpm(counts[:, keep], library)
        adjust = pd.DataFrame({"age": age[donors]}) if use_age else pd.DataFrame(index=range(len(donors)))
        logfc, t_stat, p = _ols_term(Y, (sex[donors] == "M").astype(float), adjust, ["age"])
        genes = scan.genes[keep]
        part = pd.DataFrame({"dataset": scan.name, "celltype": t, "gene_id": genes["gene_id"].to_numpy(),
                             "symbol": genes["symbol"].to_numpy(), "chrom": genes["chrom"].to_numpy(),
                             "mean_log2cpm": Y.mean(axis=0), "log2fc_male": logfc, "t": t_stat, "p": p,
                             "fdr": _fdr(p), "n_donors": len(donors), "n_male": n_m, "n_female": n_f})
        rows.append(part)
        sig = part["fdr"] < 0.05
        top = part.reindex(part["t"].abs().sort_values(ascending=False).index).head(10)
        summary.append({"dataset": scan.name, "celltype": t, "n_donors": len(donors), "n_male": n_m,
                        "n_female": n_f, "n_genes": len(part), "n_fdr05": int(sig.sum()),
                        **{f"n_fdr05_{ch}": int((sig & (part["chrom"] == ch)).sum()) for ch in ("X", "Y", "auto")},
                        "top10": ",".join(top["symbol"])})
    out = {"sex_de": pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(),
           "sex_de_summary": pd.DataFrame(summary)}
    _log_table(scan.name, "sex DE", out["sex_de_summary"])
    return out


def pathology_label(braak: np.ndarray, cerad: np.ndarray) -> np.ndarray:
    """1 = AD pathology (Braak >= III and CERAD moderate or frequent), 0 = none
    or little (Braak <= II and CERAD none or sparse), NaN otherwise. CERAD on
    the 0 (none) - 3 (frequent) scale."""
    out = np.full(len(braak), np.nan)
    out[(braak >= 3) & (cerad >= 2)] = 1
    out[(braak <= 2) & (cerad <= 1)] = 0
    return out


def read_mathys2019(path: Union[str, Path]) -> Dict[str, pd.DataFrame]:
    """Mathys 2019 supplementary table 2 (as data/mathys2019_DEGs.csv): per
    cell type, the genes with their fold change, mixed-model z and whether
    they are DEGs in both models."""
    table = pd.read_csv(path, header=[0, 1], index_col=0)
    out = {}
    for theirs, ours in MATHYS2019_CELLTYPES.items():
        part = table[theirs]
        deg = part["DEGs.Ind.Mix.models"].astype(str).str.lower() == "true"
        frame = pd.DataFrame({"symbol": [str(s).upper() for s in part.index], "fc": part["IndModel.FC"],
                              "z": part["MixedModel.z"], "deg": deg.to_numpy()}).dropna(subset=["fc"])
        out[ours] = frame[frame["deg"]].drop_duplicates("symbol")
    return out


def check_ad_de(scan: Scan, mathys2019: Optional[Union[str, Path]] = None) -> Dict[str, pd.DataFrame]:
    d = scan.donors
    num = lambda col: pd.to_numeric(d[col], errors="coerce").to_numpy(np.float64) if col in d else \
        np.full(len(d), np.nan)
    contrasts = {"pathology": pathology_label(num("braak"), num("cerad"))}
    if "disease" in d:
        contrasts["disease"] = np.where(d["disease"] == "yes", 1.0, np.where(d["disease"] == "no", 0.0, np.nan))
    sex = d["sex_recorded"].to_numpy(dtype=object)
    age = num("age")
    reference = read_mathys2019(mathys2019) if mathys2019 else {}
    rows, summary, versus = [], [], []
    for contrast, label in contrasts.items():
        usable = np.isfinite(label) & np.isin(sex, ["M", "F"]) & np.isfinite(age)
        for t in _analysed_celltypes(scan):
            donors, counts, library = _celltype_matrix(scan, t, scan.spec["min_cells"], usable)
            n_ad, n_ctrl = int(np.sum(label[donors] == 1)), int(np.sum(label[donors] == 0))
            if min(n_ad, n_ctrl) < scan.spec["min_donors_per_group"]:
                logger.warning(f"[{scan.name}] {contrast} {t}: {n_ad} vs {n_ctrl} donors, left out")
                continue
            keep = _expressed(counts, library)
            Y = _logcpm(counts[:, keep], library)
            adjust = pd.DataFrame({"sex": sex[donors], "age": age[donors]})
            logfc, t_stat, p = _ols_term(Y, label[donors], adjust, ["age"])
            genes = scan.genes[keep]
            part = pd.DataFrame({"dataset": scan.name, "contrast": contrast, "celltype": t,
                                 "gene_id": genes["gene_id"].to_numpy(), "symbol": genes["symbol"].to_numpy(),
                                 "mean_log2cpm": Y.mean(axis=0), "log2fc": logfc, "t": t_stat, "p": p,
                                 "fdr": _fdr(p), "n_ad": n_ad, "n_control": n_ctrl})
            rows.append(part)
            summary.append({"dataset": scan.name, "contrast": contrast, "celltype": t, "n_ad": n_ad,
                            "n_control": n_ctrl, "n_genes": len(part), "n_fdr05": int((part["fdr"] < 0.05).sum()),
                            "n_p05": int((part["p"] < 0.05).sum())})
            if t in reference:
                ref = reference[t]
                ours = part.assign(symbol=part["symbol"].str.upper()).drop_duplicates("symbol").set_index("symbol")
                both = ref[ref["symbol"].isin(ours.index)]
                if len(both):
                    our_t = ours.loc[both["symbol"], "t"].to_numpy()
                    our_p = ours.loc[both["symbol"], "p"].to_numpy()
                    agree = np.sign(our_t) == np.sign(both["fc"].to_numpy())
                    nominal = our_p < 0.05
                    versus.append({
                        "dataset": scan.name, "contrast": contrast, "celltype": t, "n_degs": len(ref),
                        "n_tested": len(both), "n_up": int((both["fc"] > 0).sum()),
                        "pct_same_direction": 100 * float(agree.mean()),
                        "binom_p": float(stats.binomtest(int(agree.sum()), len(agree), 0.5).pvalue),
                        "n_nominal": int(nominal.sum()),
                        "pct_same_direction_nominal": 100 * float(agree[nominal].mean()) if nominal.any() else np.nan,
                        "spearman_t_vs_z": float(stats.spearmanr(our_t, both["z"].to_numpy())[0])
                        if len(both) > 2 else np.nan,
                    })
    out = {"ad_de": pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(),
           "ad_de_summary": pd.DataFrame(summary), "ad_de_mathys2019": pd.DataFrame(versus)}
    _log_table(scan.name, "AD DE", out["ad_de_summary"])
    if len(out["ad_de_mathys2019"]):
        _log_table(scan.name, "AD DE vs Mathys 2019", out["ad_de_mathys2019"])
    return out


def _source_values(raw_values, source: Dict[str, Any]) -> np.ndarray:
    if source.get("sex"):
        return normalize_sex(raw_values, source.get("recode"))
    if source.get("numeric"):
        return be._as_numeric(raw_values, source.get("recode"))
    return be._as_labels(raw_values, source.get("recode"))


def check_metadata(spec: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    """metadata_sources: {key: [source, ...]}, a source being an obs column or
    {column, recode, numeric, tolerance, sex}. Every pair of sources of a key
    is compared per donor (the most common value of each)."""
    cfg = spec["cfg"]
    sources = {k: [s if isinstance(s, dict) else {"column": s} for s in v]
               for k, v in (spec["metadata_sources"] or {}).items()}
    columns = {cfg["donor"]} | {s["column"] for v in sources.values() for s in v}
    if spec["pool_column"]:
        columns.add(spec["pool_column"])
    path = be._first_existing(cfg["counts"])
    with h5py.File(path, "r") as f:
        missing = sorted(columns - set(obs_column_names(f)))
        if missing:
            raise ValueError(f"obs has no column(s) {missing}")
        raw = {c: read_obs_column(f, c) for c in columns}
    donors = be._as_labels(raw[cfg["donor"]])
    agreement, confusion = [], []
    for key, srcs in sources.items():
        values = {}
        for s in srcs:
            frame = pd.DataFrame({"donor": donors, "v": _source_values(raw[s["column"]], s)})
            values[s["column"]] = be._donor_majority(frame[frame["donor"].notna()], "v")
        for a, b in combinations(srcs, 2):
            va, vb = values[a["column"]], values[b["column"]]
            both = va.notna() & vb.notna()
            numeric = a.get("numeric") and b.get("numeric")
            x, y = va[both], vb[both]
            if numeric:
                tol = float(a.get("tolerance", b.get("tolerance", 0)))
                agree = (np.abs(x.astype(float) - y.astype(float)) <= tol + 1e-9)
            else:
                agree = x.astype(str) == y.astype(str)
            agreement.append({"dataset": spec["name"], "key": key, "source_a": a["column"], "source_b": b["column"],
                              "n_donors_a": int(va.notna().sum()), "n_donors_b": int(vb.notna().sum()),
                              "n_both": int(both.sum()), "n_agree": int(agree.sum()),
                              "pct_agree": 100 * float(agree.mean()) if both.any() else np.nan,
                              "max_abs_diff": float(np.abs(x.astype(float) - y.astype(float)).max())
                              if numeric and both.any() else np.nan})
            if not numeric or len(set(x) | set(y)) <= 12:
                pairs = pd.DataFrame({"value_a": x.astype(str), "value_b": y.astype(str)}).value_counts()
                confusion += [{"dataset": spec["name"], "key": key, "source_a": a["column"], "source_b": b["column"],
                               "value_a": va_, "value_b": vb_, "n_donors": int(n)}
                              for (va_, vb_), n in pairs.items()]
    pools = []
    if spec["pools"] and spec["pool_column"]:
        pool = be._as_labels(raw[spec["pool_column"]])
        for name, members in spec["pools"].items():
            hit = pd.Series(pool).isin([str(m) for m in members]).to_numpy()
            pools.append({"dataset": spec["name"], "pools": name, "n_listed": len(members),
                          "n_in_file": int(len(set(pool[hit]))), "n_cells": int(hit.sum()),
                          "n_donors": int(len(set(donors[hit & pd.notna(donors)]))),
                          "n_pools_in_file": int(len(set(pool[pd.notna(pool)])))})
    out = {"metadata_agreement": pd.DataFrame(agreement), "metadata_confusion": pd.DataFrame(confusion),
           "metadata_pools": pd.DataFrame(pools)}
    _log_table(spec["name"], "metadata agreement", out["metadata_agreement"])
    if pools:
        _log_table(spec["name"], "pools", out["metadata_pools"])
    return out


# ---------------------------------------------------------------- checks: several datasets


def check_sex_de_compare(specs: List[Dict[str, Any]]) -> Dict[str, pd.DataFrame]:
    tables = {}
    for spec in specs:
        path = Path(spec["output_dir"]) / "sex_de.tsv"
        if not path.exists():
            raise FileNotFoundError(f"{path} missing: run `pipeline qc-check sex-de` first")
        tables[spec["name"]] = pd.read_csv(path, sep="\t")
    rows, known = [], []
    for (na, a), (nb, b) in combinations(tables.items(), 2):
        for t in sorted(set(a["celltype"]) & set(b["celltype"])):
            x = a[(a["celltype"] == t) & a["gene_id"].notna()].drop_duplicates("gene_id").set_index("gene_id")
            y = b[(b["celltype"] == t) & b["gene_id"].notna()].drop_duplicates("gene_id").set_index("gene_id")
            shared = x.index.intersection(y.index)
            for chrom in ("all", "X", "Y", "auto"):
                ids = shared if chrom == "all" else shared[x.loc[shared, "chrom"] == chrom]
                if len(ids) < 3:
                    continue
                tx, ty = x.loc[ids, "t"].to_numpy(), y.loc[ids, "t"].to_numpy()
                ok = np.isfinite(tx) & np.isfinite(ty)
                top_a = set(x.loc[ids].assign(a=np.abs(tx)).nlargest(50, "a").index)
                top_b = set(y.loc[ids].assign(a=np.abs(ty)).nlargest(50, "a").index)
                both_sig = (x.loc[ids, "fdr"].to_numpy() < 0.05) & (y.loc[ids, "fdr"].to_numpy() < 0.05)
                rows.append({"dataset_a": na, "dataset_b": nb, "celltype": t, "genes": chrom, "n_genes": int(ok.sum()),
                             "pearson_t": float(np.corrcoef(tx[ok], ty[ok])[0, 1]),
                             "spearman_t": float(stats.spearmanr(tx[ok], ty[ok])[0]),
                             "top50_overlap": len(top_a & top_b), "n_fdr05_both": int(both_sig.sum()),
                             "pct_same_sign_fdr05_both": 100 * float(np.mean(np.sign(tx[both_sig]) ==
                                                                             np.sign(ty[both_sig])))
                             if both_sig.any() else np.nan})
    for name, table in tables.items():
        for gene in SEX_GENES:
            hits = table[table["symbol"].str.upper() == gene]
            for _, r in hits.iterrows():
                part = table[table["celltype"] == r["celltype"]]
                rank = int((part["t"].abs() > abs(r["t"])).sum()) + 1
                known.append({"dataset": name, "celltype": r["celltype"], "gene": gene, "chrom": r["chrom"],
                              "log2fc_male": r["log2fc_male"], "t": r["t"], "fdr": r["fdr"],
                              "rank_abs_t": rank, "n_genes": len(part)})
    out = {"sex_de_compare": pd.DataFrame(rows), "sex_de_known_genes": pd.DataFrame(known)}
    _log_table("all", "sex DE concordance", out["sex_de_compare"][out["sex_de_compare"]["genes"] != "Y"])
    return out


def _metadata_pair(a: pd.DataFrame, b: pd.DataFrame) -> List[Dict[str, Any]]:
    """Agreement of the shared donors' covariates (a and b indexed by match key)."""
    rows = []
    checks = {"sex_recorded": "exact", "age": "age", "braak": "exact", "cerad": "exact", "apoe4": "exact",
              "cogdx": "exact"}
    for column, kind in checks.items():
        if column not in a or column not in b:
            continue
        x, y = a[column], b[column]
        both = x.notna() & y.notna()
        if kind == "age":
            xv = np.minimum(pd.to_numeric(x[both], errors="coerce"), 90)
            yv = np.minimum(pd.to_numeric(y[both], errors="coerce"), 90)
            agree = np.abs(xv - yv) <= 1
        else:
            agree = x[both].astype(str).str.replace(r"\.0$", "", regex=True) == \
                y[both].astype(str).str.replace(r"\.0$", "", regex=True)
        rows.append({"covariate": column, "n_both": int(both.sum()), "n_agree": int(agree.sum()),
                     "pct_agree": 100 * float(agree.mean()) if both.any() else np.nan})
    return rows


def _identity_pair(sa: Scan, sb: Scan, keys: List[str], spec: Dict[str, Any], restricted: Path):
    """Per shared cell type: correlation of the donors' centered pseudobulk
    profiles between the datasets, and the rank of each donor's own profile."""
    da = sa.donors.reset_index().rename(columns={"index": "code"}).dropna(subset=["match"]).drop_duplicates("match")
    db = sb.donors.reset_index().rename(columns={"index": "code"}).dropna(subset=["match"]).drop_duplicates("match")
    code_a = da.set_index("match").loc[keys, "code"].to_numpy()
    code_b = db.set_index("match").loc[keys, "code"].to_numpy()
    ga, gb = sa.genes, sb.genes
    ok_a = ga["gene_id"].notna() & ga["chrom"].eq("auto")
    ok_b = gb["gene_id"].notna() & gb["chrom"].eq("auto")
    ids = pd.Index(ga.loc[ok_a, "gene_id"]).drop_duplicates().intersection(
        pd.Index(gb.loc[ok_b, "gene_id"]).drop_duplicates())
    col_a = pd.Series(np.arange(len(ga)), index=ga["gene_id"])[~ga["gene_id"].duplicated().to_numpy()].loc[ids].to_numpy()
    col_b = pd.Series(np.arange(len(gb)), index=gb["gene_id"])[~gb["gene_id"].duplicated().to_numpy()].loc[ids].to_numpy()
    n = len(keys)
    rows, stacked = [], []
    for t in [t for t in _analysed_celltypes(sa) if t in _analysed_celltypes(sb)]:
        rows_by_donor = []
        for scan in (sa, sb):
            g = scan.groups
            sel = ((g["celltype_name"] == t) & (g["n_cells"] >= spec["min_cells_identity"])).to_numpy()
            rows_by_donor.append(pd.Series(np.flatnonzero(sel), index=g["donor"].to_numpy()[sel]))
        present = np.isin(code_a, rows_by_donor[0].index) & np.isin(code_b, rows_by_donor[1].index)
        if present.sum() < spec["min_donors_per_group"]:
            continue
        profiles = []
        for scan, by_donor, codes, cols in ((sa, rows_by_donor[0], code_a, col_a), (sb, rows_by_donor[1], code_b, col_b)):
            r = by_donor.loc[codes[present]].to_numpy()
            profiles.append((scan.pseudobulk_rows(r, cols), scan.groups["library"].to_numpy()[r]))
        with np.errstate(divide="ignore", invalid="ignore"):
            cpm = [c / lib[:, None] * 1e6 for c, lib in profiles]
        expressed = (cpm[0].mean(axis=0) >= 10) & (cpm[1].mean(axis=0) >= 10)
        Z = []
        for (c, lib) in profiles:
            y = _logcpm(c[:, expressed], lib)
            Z.append(y - y.mean(axis=0))
        var_rank = stats.rankdata(Z[0].var(axis=0)) + stats.rankdata(Z[1].var(axis=0))
        top = np.argsort(-var_rank)[:spec["identity_genes"]]
        A, B = Z[0][:, top], Z[1][:, top]
        A = (A - A.mean(axis=1, keepdims=True)) / A.std(axis=1, keepdims=True)
        B = (B - B.mean(axis=1, keepdims=True)) / B.std(axis=1, keepdims=True)
        corr = A @ B.T / A.shape[1]
        full = np.full((n, n), np.nan)
        idx = np.flatnonzero(present)
        full[np.ix_(idx, idx)] = corr
        stacked.append(np.arctanh(np.clip(full, -0.999999, 0.999999)))
        rows.append(_rank_row(corr, t, len(top)))
    if not stacked:
        return rows, pd.Series(dtype=float)
    # Mean Fisher z over the cell types each pair has, back to a correlation
    combined = np.tanh(np.nanmean(np.stack(stacked), axis=0))
    observed = np.isfinite(np.diag(combined))
    sub = combined[np.ix_(observed, observed)]
    sub = np.where(np.isfinite(sub), sub, -np.inf)
    rows.append(_rank_row(sub, "combined", np.nan))
    ranks = 1 + (sub > np.diag(sub)[:, None]).sum(axis=1)
    ranks_ba = 1 + (sub > np.diag(sub)[None, :]).sum(axis=0)
    keys_obs = np.array(keys)[observed]
    pd.DataFrame({"match": keys_obs, "rank_of_own_profile": ranks})[ranks > 1].to_csv(
        restricted / f"identity_follow_up_{sa.name}_{sb.name}.tsv", sep="\t", index=False)
    return rows, pd.Series(np.maximum(ranks, ranks_ba), index=keys_obs)


def _rank_row(corr: np.ndarray, celltype: str, n_genes) -> Dict[str, Any]:
    n = corr.shape[0]
    own = np.diag(corr)
    rank_ab = 1 + (corr > own[:, None]).sum(axis=1)
    rank_ba = 1 + (corr > own[None, :]).sum(axis=0)
    finite = np.where(np.isfinite(corr), corr, np.nan)
    off = finite[~np.eye(n, dtype=bool)]
    # A donor without any own signal: its own correlation is not above the
    # 95th percentile of its correlations with the other donors
    others = np.where(np.eye(n, dtype=bool), np.nan, finite)
    p95 = np.nanpercentile(others, 95, axis=1) if n > 2 else np.full(n, np.nan)
    return {"celltype": celltype, "n_donors": n, "n_genes": n_genes,
            "pct_own_best_a_to_b": 100 * float(np.mean(rank_ab == 1)),
            "pct_own_best_b_to_a": 100 * float(np.mean(rank_ba == 1)),
            "pct_own_top5": 100 * float(np.mean(np.concatenate([rank_ab, rank_ba]) <= 5)),
            "n_rank_above_10_a_to_b": int(np.sum(rank_ab > 10)),
            "n_no_own_signal_a_to_b": int(np.sum(~(own > p95))),
            "median_rank": float(np.median(np.concatenate([rank_ab, rank_ba]))),
            "pct_chance": 100.0 / n,
            "mean_own_corr": float(np.nanmean(np.where(np.isfinite(own), own, np.nan))),
            "mean_other_corr": float(np.nanmean(off[np.isfinite(off)])) if np.isfinite(off).any() else np.nan}


def _sex_expression_pair(sa: Scan, sb: Scan, keys: List[str], restricted: Path) -> List[Dict[str, Any]]:
    """With both datasets' sex check: does a shared donor's expression sex
    agree between the datasets, and where one dataset's expression contradicts
    the recorded sex, does the other dataset's expression side with the
    record (then that one dataset's sample is the odd one out)?"""
    calls = []
    for scan in (sa, sb):
        path = Path(scan.spec["output_dir"]) / "restricted" / "sex_by_donor.tsv"
        if not path.exists():
            return [], []
        table = pd.read_csv(path, sep="\t", dtype={"donor": str})
        table["match"] = scan.donors["match"].to_numpy()
        calls.append(table.dropna(subset=["match"]).drop_duplicates("match").set_index("match").loc[keys])
    a, b = calls
    sure = a["expression"].isin(["M", "F"]) & b["expression"].isin(["M", "F"])
    rows = [{"covariate": "sex_expression", "n_both": int(sure.sum()),
             "n_agree": int((a["expression"] == b["expression"])[sure].sum())}]
    odd = []
    for name, x, y in (("a", a, b), ("b", b, a)):
        off = x["expression"].isin(["M", "F"]) & (x["expression"] != x["recorded"])
        backed = off & (y["expression"] == y["recorded"])
        rows.append({"covariate": f"sex_expression_off_in_{name}", "n_both": int(off.sum()),
                     "n_agree": int(backed.sum())})
        odd += [(k, name) for k in x.index[off]]
    for r in rows:
        r["pct_agree"] = 100 * r["n_agree"] / r["n_both"] if r["n_both"] else np.nan
    if odd:
        pd.DataFrame(odd, columns=["match", "expression_off_in"]).to_csv(
            restricted / f"sex_follow_up_{sa.name}_{sb.name}.tsv", sep="\t", index=False)
    return rows, [k for k, _ in odd]


def check_identity(specs: List[Dict[str, Any]], output: Path) -> Dict[str, pd.DataFrame]:
    scans = [Scan(s) for s in specs]
    restricted = _restricted(output)
    meta_rows, expr_rows = [], []
    for sa, sb in combinations(scans, 2):
        if "match" not in sa.donors or "match" not in sb.donors:
            continue
        a = sa.donors.dropna(subset=["match"]).drop_duplicates("match").set_index("match")
        b = sb.donors.dropna(subset=["match"]).drop_duplicates("match").set_index("match")
        keys = sorted(set(a.index) & set(b.index))
        logger.info(f"{sa.name} vs {sb.name}: {len(keys)} donors in both")
        if len(keys) < specs[0]["min_donors_per_group"]:
            continue
        sex_rows, sex_off = _sex_expression_pair(sa, sb, keys, restricted)
        rows, ranks = _identity_pair(sa, sb, keys, specs[0], restricted)
        expr_rows += [{"dataset_a": sa.name, "dataset_b": sb.name, **r} for r in rows]
        # Donors whose expression sex contradicts the record in one dataset:
        # does their expression still match their own profile in the other?
        if sex_off:
            off = ranks.reindex(sex_off).dropna()
            for name, hit in (("own_profile_best_of_sex_off_donors", off == 1),
                              ("own_profile_rank_above_10_of_sex_off_donors", off > 10)):
                sex_rows.append({"covariate": name, "n_both": len(off), "n_agree": int(hit.sum()),
                                 "pct_agree": 100 * float(hit.mean()) if len(off) else np.nan})
        for row in _metadata_pair(a.loc[keys], b.loc[keys]) + sex_rows:
            meta_rows.append({"dataset_a": sa.name, "dataset_b": sb.name, "n_shared": len(keys), **row})
    out = {"identity_metadata": pd.DataFrame(meta_rows), "identity_expression": pd.DataFrame(expr_rows)}
    _log_table("all", "shared donors: metadata", out["identity_metadata"])
    _log_table("all", "shared donors: expression", out["identity_expression"])
    return out


# ---------------------------------------------------------------- checks: genotype

GENOTYPE_Z = 6.0          # a confident identity: robust z of the score among all genotyped people
GENOTYPE_MARGIN = 3.0     # ... at least this far above the next person
GENOTYPE_PCS = 10         # expression PCs removed per layer (at most one per 10 donors)
GENOTYPE_MAX_P = 1e-3     # features kept: in-sample eQTL p at most this
GENOTYPE_MIN_CPM = 5.0    # eQTL genes expressed at least this much in a layer
GENOTYPE_MIN_DONORS = 30
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _projid(values) -> np.ndarray:
    """projid as digits without leading zeros: '00123456', 123456.0 and
    'ROS00123456' (a genotype sample) alike; '' when there is none."""
    s = pd.Series(values, dtype=object).astype(str).str.extract(r"(\d+)(?:\.0+)?$")[0]
    return s.map(lambda x: str(int(x)) if isinstance(x, str) else "").to_numpy(dtype=object)


def read_genotypes(path: Union[str, Path]):
    """A dosage table (snp, a1, a2, then one column per genotyped sample) ->
    (variants, sample IDs, dosage matrix variants x samples)."""
    table = pd.read_csv(path, sep="\t", dtype={"snp": str, "a1": str, "a2": str})
    return (table[["snp", "a1", "a2"]].reset_index(drop=True), np.asarray(table.columns[3:], dtype=str),
            table.iloc[:, 3:].to_numpy(dtype=np.float64))


def _match_eqtls(eqtls: pd.DataFrame, variants: pd.DataFrame) -> np.ndarray:
    """Per eQTL, the genotype row of its variant (by rsID or GRCh37 chr:pos,
    with the same alleles on either strand), else -1."""
    rows: Dict[str, List[int]] = {}
    for i, snp in enumerate(variants["snp"]):
        rows.setdefault(snp, []).append(i)
    a1, a2 = variants["a1"].to_numpy(), variants["a2"].to_numpy()
    out = []
    for e in eqtls.itertuples():
        alleles = {str(e.ref), str(e.alt)}
        flipped = {a.translate(_COMPLEMENT) for a in alleles}
        candidates = list(rows.get(str(e.rsid), []))
        if pd.notna(e.pos_b37):
            candidates += rows.get(f"chr{e.chrom}:{int(e.pos_b37)}", [])
        out.append(next((i for i in candidates if {a1[i], a2[i]} in (alleles, flipped)), -1))
    return np.asarray(out, dtype=np.int64)


def _genotype_layers(scan: Scan, gene_cols: np.ndarray, spec: Dict[str, Any]):
    """Per layer (all nuclei of a donor, and each cell type): the donors in it
    and the eQTL genes' log-CPM with the top expression PCs regressed out,
    standardized (NaN for genes expressed too little in the layer)."""
    g = scan.groups
    n = scan.n_donor_ids
    total = np.zeros((n, len(scan.genes)))
    total_cells = np.zeros(n)
    layers = []
    for t in scan.celltypes:
        rows = np.flatnonzero((g["celltype_name"] == t).to_numpy())
        if not len(rows):
            continue
        counts = scan.pseudobulk_rows(rows)
        donors = g["donor"].to_numpy()[rows]
        np.add.at(total, donors, counts)
        n_cells = g["n_cells"].to_numpy()[rows]
        np.add.at(total_cells, donors, n_cells)
        if t != OTHER:
            keep = n_cells >= spec["min_cells"]
            layers.append((t, donors[keep], counts[keep]))
    keep = total_cells >= spec["min_cells"]
    layers.insert(0, ("all", np.flatnonzero(keep), total[keep]))
    del total
    out = []
    for name, donors, counts in layers:
        if len(donors) < GENOTYPE_MIN_DONORS:
            continue
        library = counts.sum(axis=1)
        cpm = counts / np.maximum(library, 1)[:, None] * 1e6
        mean_cpm = cpm.mean(axis=0)
        y = np.log2(cpm[:, mean_cpm >= 10] + 1)
        sd = y.std(axis=0)
        y = (y[:, sd > 0] - y[:, sd > 0].mean(axis=0)) / sd[sd > 0]
        k = min(GENOTYPE_PCS, len(donors) // 10, y.shape[1])
        X = np.hstack([np.ones((len(donors), 1)), np.linalg.svd(y, full_matrices=False)[0][:, :k]])
        feature = np.log2(cpm[:, gene_cols] + 1)
        resid = feature - X @ np.linalg.lstsq(X, feature, rcond=None)[0]
        sd = resid.std(axis=0)
        ok = (mean_cpm[gene_cols] >= GENOTYPE_MIN_CPM) & (sd > 0)
        z = np.full(resid.shape, np.nan)
        z[:, ok] = (resid[:, ok] - resid[:, ok].mean(axis=0)) / sd[ok]
        out.append((name, donors, z, k))
    return out


def _genotype_scores(E: np.ndarray, Zf: np.ndarray, own: np.ndarray):
    """Score of each donor's expression against each genotyped sample.

    E: donors x features (standardized expression, NaN = no data); Zf:
    samples x features (standardized dosage of the feature's variant); own:
    each donor's labelled sample (-1 for none). The eQTL effect r of each
    feature is estimated over the labelled donors (left out for the donor
    scored); features with p > GENOTYPE_MAX_P get no weight. Score =
    sum_f r_f * E_df * Zf_sf, as a robust z over all samples."""
    has = own >= 0
    Zown = np.full(E.shape, np.nan)
    Zown[has] = Zf[own[has]]
    valid = np.isfinite(E) & np.isfinite(Zown)
    n = valid.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        e = np.where(valid, E, np.nan)
        e = (e - np.nanmean(e, axis=0)) / np.nanstd(e, axis=0)
        z = np.where(valid, Zown, np.nan)
        z = (z - np.nanmean(z, axis=0)) / np.nanstd(z, axis=0)
        r = np.nanmean(e * z, axis=0)
        t = r * np.sqrt((n - 2) / np.maximum(1 - r ** 2, 1e-12))
    p = np.where(n >= GENOTYPE_MIN_DONORS, 2 * stats.t.sf(np.abs(np.nan_to_num(t)), np.maximum(n - 2, 1)), 1.0)
    kept = np.isfinite(r) & (p <= GENOTYPE_MAX_P)
    W = np.where(valid, (n * r - np.nan_to_num(e * z)) / np.maximum(n - 1, 1), r)
    W[:, ~kept] = 0
    S = (np.nan_to_num(W) * np.nan_to_num(E)) @ np.nan_to_num(Zf).T
    med = np.median(S, axis=1, keepdims=True)
    mad = 1.4826 * np.median(np.abs(S - med), axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        Zs = np.where(mad > 0, (S - med) / mad, np.nan)
    return Zs, kept, r, n


def _genotype_calls(Zs: np.ndarray, own: np.ndarray) -> pd.DataFrame:
    rows = []
    for d in range(Zs.shape[0]):
        z = Zs[d]
        if not np.isfinite(z).any():
            rows.append({"call": "no data", "own_z": np.nan, "best": -1, "best_z": np.nan, "next_z": np.nan})
            continue
        order = np.argsort(-np.nan_to_num(z, nan=-np.inf))
        best, best_z, next_z = int(order[0]), float(z[order[0]]), float(z[order[1]])
        if own[d] >= 0:
            own_z = float(z[own[d]])
            other_z = float(z[order[1]] if best == own[d] else best_z)
            if best == own[d] and own_z >= GENOTYPE_Z and own_z - other_z >= GENOTYPE_MARGIN:
                call = "own genotype"
            elif best != own[d] and best_z >= GENOTYPE_Z and best_z - own_z >= GENOTYPE_MARGIN \
                    and best_z - next_z >= GENOTYPE_MARGIN:
                call = "other person"
            elif max(own_z, best_z) < GENOTYPE_Z:
                call = "no signal"
            else:
                call = "unclear"
        else:
            own_z = np.nan
            call = ("identified" if best_z >= GENOTYPE_Z and best_z - next_z >= GENOTYPE_MARGIN
                    else "not genotyped")
        rows.append({"call": call, "own_z": own_z, "best": best, "best_z": best_z, "next_z": next_z})
    return pd.DataFrame(rows)


def check_genotype(specs: List[Dict[str, Any]], output: Path, genotypes: Union[str, Path],
                   eqtls: Union[str, Path], id_map: Union[str, Path]) -> Dict[str, pd.DataFrame]:
    """Does each donor's expression match the genotype of the person it is
    labelled as? Over strong cis-eQTLs (e.g. GTEx cortex, scripts/gtex_eqtls.py),
    each donor is scored against every genotyped sample, as in MixupMapper
    (Westra et al. 2011). id_map: CSV with individualID and projid (the match
    key -> the projid in the genotype sample IDs)."""
    variants, samples, dosage = read_genotypes(genotypes)
    eqtl = pd.read_csv(eqtls, sep="\t", dtype={"chrom": str, "rsid": str, "gene_id": str})
    eqtl["row"] = _match_eqtls(eqtl, variants)
    eqtl = eqtl[eqtl["row"] >= 0].reset_index(drop=True)
    sample_pid = _projid(samples)
    first = ~pd.Series(sample_pid).duplicated().to_numpy() & (sample_pid != "")
    samples, sample_pid, dosage = samples[first], sample_pid[first], dosage[:, first]
    with np.errstate(invalid="ignore", divide="ignore"):
        Zall = ((dosage - dosage.mean(axis=1, keepdims=True)) / dosage.std(axis=1, keepdims=True)).T
    logger.info(f"genotypes: {len(samples)} samples, {len(variants)} variants; {len(eqtl)} eQTL pairs found "
                f"({eqtl['gene_id'].nunique()} genes)")
    ids = pd.read_csv(id_map, dtype=str, usecols=["individualID", "projid"]).dropna()
    to_pid = pd.Series(_projid(ids["projid"]), index=ids["individualID"].to_numpy())
    to_pid = to_pid[~to_pid.index.duplicated()]
    sample_index = pd.Series(np.arange(len(samples)), index=sample_pid)
    restricted = _restricted(output)

    summary, features, sex_rows, hist, by_person = [], [], [], [], {}
    for spec in specs:
        scan = Scan(spec)
        donors = scan.donors
        match = donors["match"] if "match" in donors else pd.Series([np.nan] * len(donors))
        label_pid = match.map(to_pid).fillna("").to_numpy(dtype=object)
        own = sample_index.reindex(label_pid).fillna(-1).to_numpy(dtype=np.int64)
        gid = pd.Series(np.arange(len(scan.genes)), index=scan.genes["gene_id"].fillna(""))
        gid = gid[~gid.index.duplicated()]
        pairs = eqtl[eqtl["gene_id"].isin(gid.index)].reset_index(drop=True)
        gene_cols = gid.loc[pairs["gene_id"]].to_numpy()
        blocks, snp_rows, layer_of = [], [], []
        for name, layer_donors, z, k in _genotype_layers(scan, gene_cols, spec):
            ok = np.flatnonzero(np.isfinite(z).any(axis=0))
            block = np.full((len(donors), len(ok)), np.nan)
            block[layer_donors] = z[:, ok]
            blocks.append(block)
            snp_rows.append(pairs["row"].to_numpy()[ok])
            layer_of += [(name, k)] * len(ok)
        E = np.hstack(blocks)
        Zf = Zall[:, np.concatenate(snp_rows)]
        Zs, kept, r, n_f = _genotype_scores(E, Zf, own)
        calls = _genotype_calls(Zs, own)
        calls.insert(0, "donor", donors["donor"].to_numpy())
        calls.insert(1, "match", match.to_numpy())
        calls.insert(2, "label_projid", label_pid)
        calls["best_projid"] = np.where(calls["best"] >= 0, sample_pid[np.maximum(calls["best"], 0)], "")
        # A best match that is another donor of this dataset: a swap or a duplicate label
        in_data = pd.Series(np.arange(len(donors)), index=label_pid)
        in_data = in_data[(in_data.index != "") & ~in_data.index.duplicated()]
        other = calls["call"].isin(["other person", "identified"])
        calls["best_is_donor"] = ""
        hit = other & calls["best_projid"].isin(in_data.index)
        calls.loc[hit, "best_is_donor"] = donors["donor"].to_numpy()[in_data.loc[calls.loc[hit, "best_projid"]].to_numpy()]
        mutual = 0
        for i in np.flatnonzero(hit.to_numpy()):
            j = in_data.loc[calls.at[i, "best_projid"]]
            mutual += int(calls.at[j, "best_projid"] == label_pid[i] and calls.at[j, "call"] == "other person")
        confident = calls["call"].isin(["own genotype", "other person", "identified"])
        same_twice = int(calls.loc[confident, "best_projid"].duplicated(keep=False).sum())
        calls.drop(columns="best").to_csv(restricted / f"genotype_calls_{scan.name}.tsv", sep="\t", index=False)
        for m, pid, c in zip(match, calls["best_projid"], calls["call"]):
            if isinstance(m, str) and c in ("own genotype", "other person", "identified"):
                by_person.setdefault(m, {})[scan.name] = pid

        labelled = own >= 0
        matched = calls["call"] == "own genotype"
        other_z = np.where(calls["best"] == own, calls["next_z"], calls["best_z"])
        counts = calls["call"].value_counts()
        row = {"dataset": scan.name, "n_donors": len(donors), "n_genotyped_label": int(labelled.sum()),
               "n_features": int(E.shape[1]), "n_features_kept": int(kept.sum()),
               "n_eqtl_pairs": len(pairs), "n_genotype_samples": len(samples),
               "median_own_z": float(np.nanmedian(calls.loc[labelled, "own_z"])) if labelled.any() else np.nan}
        for c in ("own genotype", "other person", "unclear", "no signal", "no data", "identified", "not genotyped"):
            row[f"n_{c.replace(' ', '_')}"] = int(counts.get(c, 0))
        row["n_other_person_in_dataset"] = int((hit & (calls["call"] == "other person")).sum())
        row["n_mutual_swaps"] = mutual // 2
        row["n_same_person_twice"] = same_twice
        row["p99_next_z_of_own"] = float(np.nanpercentile(other_z[matched], 99)) if matched.any() else np.nan
        row["max_next_z_of_own"] = float(np.nanmax(other_z[matched])) if matched.any() else np.nan
        summary.append(row)

        layer_names = np.array([name for name, _ in layer_of])
        for name in dict.fromkeys(layer_names):
            sel = layer_names == name
            features.append({"dataset": scan.name, "layer": name, "n_pcs": dict(layer_of)[name],
                             "n_features": int(sel.sum()), "n_kept": int((kept & sel).sum()),
                             "median_donors": float(np.median(n_f[sel])),
                             "median_abs_r_kept": float(np.median(np.abs(r[kept & sel]))) if (kept & sel).any()
                             else np.nan})
        edges = np.arange(-10, 60.5, 1.0)
        for measure, values in (("own_z", calls.loc[labelled, "own_z"].to_numpy()),
                                ("best_other_z", other_z[labelled.astype(bool)])):
            v = np.clip(values[np.isfinite(values)], edges[0], edges[-1] - 1e-9)
            h, _ = np.histogram(v, bins=edges)
            hist += [{"dataset": scan.name, "measure": measure, "bin_left": float(lo), "bin_right": float(hi),
                      "n_donors": int(k)} for lo, hi, k in zip(edges[:-1], edges[1:], h)]

        sex_path = Path(spec["output_dir"]) / "restricted" / "sex_by_donor.tsv"
        if sex_path.exists():
            sex = pd.read_csv(sex_path, sep="\t", dtype={"donor": str})
            off = (sex["expression"].isin(["M", "F"]) & sex["recorded"].isin(["M", "F"])
                   & (sex["expression"] != sex["recorded"])).to_numpy()
            amb = (sex["expression"] == "ambiguous").to_numpy()
            for what, mask in (("expression sex differs from the record", off), ("expression sex ambiguous", amb)):
                vc = calls.loc[mask, "call"].value_counts()
                sex_rows.append({"dataset": scan.name, "donors": what, "n": int(mask.sum()),
                                 **{f"n_{c.replace(' ', '_')}": int(vc.get(c, 0)) for c in
                                    ("own genotype", "other person", "unclear", "no signal", "identified",
                                     "not genotyped")}})
        logger.info(f"[{scan.name}] genotype calls: {counts.to_dict()}")

    cross = []
    for a, b in combinations([s["name"] for s in specs], 2):
        both = [(p[a], p[b]) for p in by_person.values() if a in p and b in p]
        cross.append({"dataset_a": a, "dataset_b": b, "n_persons_identified_in_both": len(both),
                      "n_same_genotype": sum(x == y for x, y in both)})
    out = {"genotype_summary": pd.DataFrame(summary), "genotype_features": pd.DataFrame(features),
           "genotype_sex": pd.DataFrame(sex_rows), "genotype_cross": pd.DataFrame(cross),
           "genotype_histograms": pd.DataFrame(hist)}
    for name in ("genotype_summary", "genotype_features", "genotype_sex", "genotype_cross"):
        _log_table("all", name, out[name])
    return out


# ---------------------------------------------------------------- entry point


def run_check(check: str, spec_paths: Sequence[Union[str, Path]], output: Optional[Union[str, Path]] = None,
              mathys2019: Optional[Union[str, Path]] = None, genotypes: Optional[Union[str, Path]] = None,
              eqtls: Optional[Union[str, Path]] = None, id_map: Optional[Union[str, Path]] = None) -> Path:
    if check not in CHECKS:
        raise ValueError(f"Unknown check {check!r}; one of {CHECKS}")
    specs = [load_qc_spec(p) for p in spec_paths]
    if check == "genotype":
        if not (genotypes and eqtls and id_map):
            raise ValueError("genotype needs --genotypes, --eqtls and --id-map")
        output = Path(output or Path("data/v2/qc") / check)
        output.mkdir(parents=True, exist_ok=True)
        _write_tables(check_genotype(specs, output, genotypes, eqtls, id_map), output, check, spec_paths)
        return output
    if check in MULTI_CHECKS:
        if len(specs) < 2:
            raise ValueError(f"{check} needs two or more specs")
        output = Path(output or Path("data/v2/qc") / check)
        output.mkdir(parents=True, exist_ok=True)
        tables = check_sex_de_compare(specs) if check == "sex-de-compare" else check_identity(specs, output)
        _write_tables(tables, output, check, spec_paths)
        return output
    written = None
    for spec, spec_path in zip(specs, spec_paths):
        out = Path(output) if output and len(specs) == 1 else Path(spec["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        if check == "metadata":
            tables = check_metadata(spec)
        else:
            scan = Scan(spec)
            if check == "ad-de":
                tables = check_ad_de(scan, mathys2019)
            else:
                tables = {"totals": check_totals, "cells": check_cells, "sex": check_sex,
                          "markers": check_markers, "sex-de": check_sex_de}[check](scan)
        _write_tables(tables, out, check, [spec_path])
        written = out
    return written


def _write_tables(tables: Dict[str, pd.DataFrame], output: Path, check: str, spec_paths) -> None:
    for name, table in tables.items():
        path = output / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False, float_format="%.6g")
        logger.info(f"wrote {path} ({len(table)} rows)")
    _write_json(output / f"qc_{check}.json", f"qc-check {check}", {"specs": [str(p) for p in spec_paths]},
                {"tables": sorted(tables)})
