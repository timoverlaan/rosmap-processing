"""Map dataset gene symbols to Ensembl gene IDs.

The datasets carry only gene symbols, from different 10x references:

- ROSMAP-MIT (Mathys 2023): 33,538 features = Cell Ranger GRCh38-3.0.0
  (Ensembl 93 GTF, filtered on gene_biotype).
- SeaAD and the Green 2024 ROSMAP atlas: 36,601 features = GRCh38-2020-A
  (GENCODE v32 primary assembly, filtered on gene and transcript type, without
  readthrough and chrY PAR transcripts).

Symbols change between releases, so matching across datasets by symbol loses
genes. This module rebuilds each reference's feature list from the public GTF
with 10x's published filters, applies the same name de-duplication the data
went through (R `make.unique` or anndata's `-1` suffixes), and maps each
dataset gene to its Ensembl ID. When the rebuilt list equals the dataset's genes
in order, the map is exact by position; otherwise genes are matched by name and
the unmatched ones are reported.
"""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

from ..utils.logging import get_logger

logger = get_logger(__name__)

_IG_TR = [
    "IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene",
    "IG_LV_gene", "IG_V_gene", "IG_V_pseudogene",
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene",
]

# 10x's documented build steps for each reference.
RECIPES: Dict[str, Dict] = {
    "cellranger-GRCh38-3.0.0": {
        "gtf": "Homo_sapiens.GRCh38.93.gtf.gz (Ensembl release 93)",
        "type_attr": "gene_biotype",
        "biotypes": frozenset(["protein_coding", "lincRNA", "antisense", *_IG_TR]),
        "transcript_filter": False,
        "n_genes": 33538,
    },
    "cellranger-GRCh38-2020-A": {
        "gtf": "gencode.v32.primary_assembly.annotation.gtf.gz (GENCODE release 32)",
        "type_attr": "gene_type",
        "biotypes": frozenset(["protein_coding", "lncRNA", *_IG_TR]),
        "transcript_filter": True,
        "n_genes": 36601,
    },
}

_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
_VERSION_RE = re.compile(r"^(ENS[GTE]\d+)\.\d+")
_EXCLUDED_TAGS = {"readthrough_transcript", "PAR"}


def _open_text(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_attributes(field: str) -> Tuple[Dict[str, str], List[str]]:
    attrs: Dict[str, str] = {}
    tags: List[str] = []
    for key, value in _ATTR_RE.findall(field):
        if key == "tag":
            tags.append(value)
        else:
            attrs.setdefault(key, value)
    return attrs, tags


def strip_version(gene_id: str) -> str:
    """ENSG00000223972.5 -> ENSG00000223972 (a _PAR_Y suffix is kept)."""
    return _VERSION_RE.sub(r"\1", gene_id)


def read_reference_genes(gtf_path: Union[str, Path], recipe: str) -> pd.DataFrame:
    """Rebuild a 10x reference's gene list (gene_id, gene_name) in GTF order."""
    if recipe not in RECIPES:
        raise ValueError(f"Unknown recipe {recipe!r}; choose from {list(RECIPES)}")
    spec = RECIPES[recipe]
    type_attr = spec["type_attr"]
    biotypes = spec["biotypes"]

    genes: List[Tuple[str, str]] = []
    genes_with_ok_transcript = set()
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature == "gene":
                attrs, _ = _parse_attributes(parts[8])
                if attrs.get(type_attr) in biotypes:
                    gene_id = strip_version(attrs["gene_id"])
                    genes.append((gene_id, attrs.get("gene_name", gene_id)))
            elif feature == "transcript" and spec["transcript_filter"]:
                attrs, tags = _parse_attributes(parts[8])
                if (
                    attrs.get("gene_type") in biotypes
                    and attrs.get("transcript_type") in biotypes
                    and not _EXCLUDED_TAGS.intersection(tags)
                ):
                    genes_with_ok_transcript.add(strip_version(attrs["gene_id"]))

    if spec["transcript_filter"]:
        genes = [g for g in genes if g[0] in genes_with_ok_transcript]

    df = pd.DataFrame(genes, columns=["gene_id", "gene_name"])
    logger.info(
        f"{recipe}: rebuilt {len(df)} genes from {Path(gtf_path).name} "
        f"(10x reference has {spec['n_genes']})"
    )
    if len(df) != spec["n_genes"]:
        logger.warning(
            f"{recipe}: gene count differs from the 10x reference; mapping falls back "
            f"to name matching where positions disagree"
        )
    return df


def make_unique(names: Sequence[str], sep: str = ".") -> List[str]:
    """De-duplicate names like R's make.unique / anndata's make_index_unique.

    Later occurrences of a name get `sep` + the first counter (from 1) that does
    not collide with any existing name.
    """
    used = set(names)
    counters: Dict[str, int] = {}
    seen = set()
    out = []
    for name in names:
        if name not in seen:
            seen.add(name)
            out.append(name)
            continue
        i = counters.get(name, 0)
        while True:
            i += 1
            candidate = f"{name}{sep}{i}"
            if candidate not in used:
                break
        counters[name] = i
        used.add(candidate)
        out.append(candidate)
    return out


def build_gene_map(var_names: Sequence[str], reference: pd.DataFrame) -> pd.DataFrame:
    """Map each dataset gene to an Ensembl ID using a rebuilt reference list.

    Returns a DataFrame with columns var_name, gene_id (None if unmatched),
    ref_name and method ('position', 'name', 'unique_symbol' or 'unmatched'),
    in dataset order.
    """
    var_names = [str(v) for v in var_names]
    ref_ids = reference["gene_id"].tolist()
    ref_names = reference["gene_name"].tolist()

    candidates = {sep: make_unique(ref_names, sep) for sep in (".", "-")}

    for sep, unique_names in candidates.items():
        if unique_names == var_names:
            logger.info(f"Exact positional match with the reference (duplicate suffix '{sep}')")
            return pd.DataFrame({
                "var_name": var_names,
                "gene_id": ref_ids,
                "ref_name": ref_names,
                "method": "position",
            })

    # Name matching: use the de-duplication convention that matches most genes.
    sep = max(candidates, key=lambda s: len(set(candidates[s]) & set(var_names)))
    by_unique_name = dict(zip(candidates[sep], zip(ref_ids, ref_names)))
    name_counts = pd.Series(ref_names).value_counts()
    by_symbol = {
        name: (gid, name)
        for gid, name in zip(ref_ids, ref_names)
        if name_counts[name] == 1
    }

    rows = []
    for v in var_names:
        if v in by_unique_name:
            gid, ref_name = by_unique_name[v]
            rows.append((v, gid, ref_name, "name"))
        elif v in by_symbol:
            gid, ref_name = by_symbol[v]
            rows.append((v, gid, ref_name, "unique_symbol"))
        else:
            rows.append((v, None, None, "unmatched"))
    gene_map = pd.DataFrame(rows, columns=["var_name", "gene_id", "ref_name", "method"])

    dup = gene_map["gene_id"].notna() & gene_map["gene_id"].duplicated(keep=False)
    if dup.any():
        logger.warning(f"{int(dup.sum())} genes share an Ensembl ID; these are marked ambiguous")
        gene_map.loc[dup, "method"] = "ambiguous"
        gene_map.loc[dup, "gene_id"] = None

    counts = gene_map["method"].value_counts().to_dict()
    logger.info(f"Name matching (duplicate suffix '{sep}'): {counts}")
    return gene_map


def read_var(h5ad_path: Union[str, Path]) -> pd.DataFrame:
    """The var table of an h5ad file without loading X or obs."""
    import h5py

    try:
        from anndata.io import read_elem
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem

    with h5py.File(h5ad_path, "r") as f:
        var = read_elem(f["var"])
    var.index = var.index.astype(str)
    return var


def read_var_names(h5ad_path: Union[str, Path]) -> List[str]:
    """Gene names of an h5ad file without loading X or obs."""
    return read_var(h5ad_path).index.tolist()


def build_gene_map_file(
    gtf_path: Union[str, Path],
    recipe: str,
    h5ad_path: Union[str, Path],
    output_path: Union[str, Path],
) -> pd.DataFrame:
    """Write the gene map TSV for one dataset and log a summary."""
    reference = read_reference_genes(gtf_path, recipe)
    var_names = read_var_names(h5ad_path)
    logger.info(f"{Path(h5ad_path).name}: {len(var_names)} genes")
    gene_map = build_gene_map(var_names, reference)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gene_map.to_csv(output_path, sep="\t", index=False)
    n_mapped = int(gene_map["gene_id"].notna().sum())
    logger.info(
        f"Wrote {output_path}: {n_mapped}/{len(gene_map)} genes mapped to Ensembl IDs"
    )
    return gene_map


def load_gene_map(path: Union[str, Path]) -> Dict[str, str]:
    """var_name -> Ensembl ID, for mapped genes only."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["gene_id"].notna()]
    return dict(zip(df["var_name"], df["gene_id"]))


def shared_gene_ids(map_a: Union[str, Path], *other_maps: Union[str, Path]) -> List[str]:
    """Ensembl IDs mapped in every gene map, in map_a's order."""
    ids = list(load_gene_map(map_a).values())
    for other in other_maps:
        other_ids = set(load_gene_map(other).values())
        ids = [gid for gid in ids if gid in other_ids]
    return ids


def write_shared_gene_ids(
    maps: Sequence[Union[str, Path]], output_path: Union[str, Path]
) -> List[str]:
    if len(maps) < 2:
        raise ValueError("Need at least two gene maps")
    shared = shared_gene_ids(*maps)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(shared) + "\n")
    names = ", ".join(Path(m).name for m in maps)
    logger.info(f"Wrote {output_path}: {len(shared)} Ensembl IDs shared by {names}")
    return shared
