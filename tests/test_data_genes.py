"""Tests for rebuilding 10x reference gene lists and mapping symbols to Ensembl IDs."""

import gzip

import pandas as pd
import pytest

from rosmap_processing.data.genes import (
    build_gene_map,
    make_unique,
    read_reference_genes,
    shared_gene_ids,
)


def _gtf_line(feature, gene_id, name, gene_type, type_key, tx_type=None, tags=()):
    attrs = f'gene_id "{gene_id}"; gene_name "{name}"; {type_key} "{gene_type}";'
    if feature == "transcript":
        attrs += f' transcript_id "T{gene_id}"; transcript_type "{tx_type or gene_type}";'
        attrs += "".join(f' tag "{t}";' for t in tags)
    return "\t".join(["chr1", "TEST", feature, "1", "100", ".", "+", ".", attrs]) + "\n"


def _write(path, lines):
    with gzip.open(path, "wt") as f:
        f.write("#!genome-build GRCh38\n")
        f.writelines(lines)


def test_make_unique_matches_r_semantics():
    assert make_unique(["A", "A", "A"]) == ["A", "A.1", "A.2"]
    # An existing name is never reused as a suffix
    assert make_unique(["A", "A", "A.1"]) == ["A", "A.2", "A.1"]
    assert make_unique(["A", "B", "A"], sep="-") == ["A", "B", "A-1"]


def test_cellranger_3_recipe_filters_on_gene_biotype(tmp_path):
    gtf = tmp_path / "ens93.gtf.gz"
    _write(gtf, [
        _gtf_line("gene", "ENSG01", "AAA", "protein_coding", "gene_biotype"),
        _gtf_line("gene", "ENSG02", "BBB", "processed_pseudogene", "gene_biotype"),
        _gtf_line("gene", "ENSG03", "AAA", "lincRNA", "gene_biotype"),
    ])
    ref = read_reference_genes(gtf, "cellranger-GRCh38-3.0.0")
    assert ref["gene_id"].tolist() == ["ENSG01", "ENSG03"]


def test_2020a_recipe_drops_readthrough_par_and_versions(tmp_path):
    gtf = tmp_path / "gencode32.gtf.gz"
    t = "gene_type"
    _write(gtf, [
        _gtf_line("gene", "ENSG01.5", "AAA", "protein_coding", t),
        _gtf_line("transcript", "ENSG01.5", "AAA", "protein_coding", t),
        # readthrough-only gene: excluded
        _gtf_line("gene", "ENSG02.1", "AAA-BBB", "protein_coding", t),
        _gtf_line("transcript", "ENSG02.1", "AAA-BBB", "protein_coding", t, tags=["readthrough_transcript"]),
        # chrY PAR copy: excluded; X copy kept
        _gtf_line("gene", "ENSG03.2", "PARX", "protein_coding", t),
        _gtf_line("transcript", "ENSG03.2", "PARX", "protein_coding", t),
        _gtf_line("gene", "ENSG03.2_PAR_Y", "PARX", "protein_coding", t),
        _gtf_line("transcript", "ENSG03.2_PAR_Y", "PARX", "protein_coding", t, tags=["PAR"]),
        # lncRNA kept, pseudogene not
        _gtf_line("gene", "ENSG04.1", "LNC", "lncRNA", t),
        _gtf_line("transcript", "ENSG04.1", "LNC", "lncRNA", t),
        _gtf_line("gene", "ENSG05.1", "PSG", "unprocessed_pseudogene", t),
        _gtf_line("transcript", "ENSG05.1", "PSG", "unprocessed_pseudogene", t),
    ])
    ref = read_reference_genes(gtf, "cellranger-GRCh38-2020-A")
    assert ref["gene_id"].tolist() == ["ENSG01", "ENSG03", "ENSG04"]


def test_build_gene_map_positional_with_duplicates():
    ref = pd.DataFrame({"gene_id": ["E1", "E2", "E3"], "gene_name": ["A", "B", "A"]})
    gm = build_gene_map(["A", "B", "A.1"], ref)
    assert gm["method"].eq("position").all()
    assert gm["gene_id"].tolist() == ["E1", "E2", "E3"]


def test_build_gene_map_by_name_reports_unmatched():
    ref = pd.DataFrame({"gene_id": ["E1", "E2", "E3"], "gene_name": ["A", "B", "C"]})
    gm = build_gene_map(["C", "A", "NEWNAME"], ref)
    assert gm["gene_id"].tolist()[:2] == ["E3", "E1"]
    assert gm.loc[2, "method"] == "unmatched"


def test_shared_gene_ids(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    pd.DataFrame({"var_name": ["A", "B", "C"], "gene_id": ["E1", "E2", None]}).to_csv(a, sep="\t", index=False)
    pd.DataFrame({"var_name": ["B2", "A", "D"], "gene_id": ["E2", "E1", "E9"]}).to_csv(b, sep="\t", index=False)
    assert shared_gene_ids(a, b) == ["E1", "E2"]

    c = tmp_path / "c.tsv"
    pd.DataFrame({"var_name": ["X"], "gene_id": ["E2"]}).to_csv(c, sep="\t", index=False)
    assert shared_gene_ids(a, b, c) == ["E2"]
