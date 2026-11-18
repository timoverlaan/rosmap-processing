"""Strong cis-eQTLs in GTEx v8 cortex, for the genotype check of the data QC
(`pipeline qc-check genotype`): public summary statistics from the GTEx API.

For the GENES strongest eGenes (lowest nominal p of the lead variant, over
Brain_Frontal_Cortex_BA9 and Brain_Cortex), all conditionally independent
eQTLs with MAF >= 0.05 and p <= 1e-6, with their GRCh37 position (from GTEx's
variant table) so they can be found in GRCh37 genotype files by rsID or by
chr:pos.

    python scripts/gtex_eqtls.py data/v2/reference/gtex_cortex_eqtls.tsv [--genes 2000]
"""

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://gtexportal.org/api/v2"
TISSUES = ("Brain_Frontal_Cortex_BA9", "Brain_Cortex")
BATCH = 50


def get(endpoint: str, params) -> list:
    """All pages of a GTEx API v2 endpoint (params: list of (key, value))."""
    out, page = [], 0
    while True:
        query = urllib.parse.urlencode(list(params) + [("datasetId", "gtex_v8"), ("itemsPerPage", 100000),
                                                       ("page", page)])
        request = urllib.request.Request(f"{API}/{endpoint}?{query}", headers={"User-Agent": "eqtl-list/1.0"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    body = json.load(response)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        out += body["data"]
        page += 1
        if page >= body["paging_info"]["numberOfPages"]:
            return out


def get_batched(endpoint: str, key: str, values: list, extra=()) -> list:
    """`get` over values in batches; a batch the server fails on is split in
    halves (a lone value it fails on is skipped and counted)."""
    out, failed = [], []

    def run(batch):
        try:
            out.extend(get(endpoint, [(key, v) for v in batch] + list(extra)))
        except urllib.error.HTTPError:
            if len(batch) == 1:
                failed.append(batch[0])
            else:
                run(batch[:len(batch) // 2])
                run(batch[len(batch) // 2:])

    for i in range(0, len(values), BATCH):
        run(values[i:i + BATCH])
    if failed:
        print(f"{endpoint}: the server failed on {len(failed)} of {len(values)} {key}s; skipped")
    return out


def add_grch37(table: pd.DataFrame) -> pd.Series:
    """GRCh37 position of each rsID from Ensembl's GRCh37 REST server (the
    mapping on the same chromosome; none if it is not unique there)."""
    rsids = sorted(set(table["rsid"].dropna()) - {""})
    positions = {}
    for i in range(0, len(rsids), 200):
        body = json.dumps({"ids": rsids[i:i + 200]}).encode()
        request = urllib.request.Request("https://grch37.rest.ensembl.org/variation/homo_sapiens", data=body,
                                         headers={"Content-Type": "application/json", "Accept": "application/json",
                                                  "User-Agent": "eqtl-list/1.0"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    found = json.load(response)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        for rsid, record in found.items():
            positions[rsid] = {(m["seq_region_name"], m["start"]) for m in record.get("mappings", [])
                               if m.get("assembly_name") == "GRCh37"}
    out = []
    for rsid, chrom in zip(table["rsid"], table["chrom"]):
        here = [p for c, p in positions.get(rsid, ()) if c == chrom]
        out.append(here[0] if len(here) == 1 else pd.NA)
    return pd.array(out, dtype="Int64")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path)
    parser.add_argument("--genes", type=int, default=2000)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--max-p", type=float, default=1e-6)
    args = parser.parse_args()

    egenes = pd.concat([pd.DataFrame(get("association/egene", [("tissueSiteDetailId", t)])) for t in TISSUES])
    best = egenes.groupby("gencodeId")["pValue"].min().sort_values()
    genes = list(best.index[:args.genes])
    print(f"{egenes['gencodeId'].nunique()} eGenes in {', '.join(TISSUES)}; keeping the {len(genes)} strongest")

    eqtls = []
    for tissue in TISSUES:
        eqtls += get_batched("association/independentEqtl", "gencodeId", genes, [("tissueSiteDetailId", tissue)])
    eqtls = pd.DataFrame(eqtls)
    eqtls = eqtls[(eqtls["maf"] >= args.min_maf) & (eqtls["pValue"] <= args.max_p)]
    print(f"{len(eqtls)} independent eQTLs (MAF >= {args.min_maf}, p <= {args.max_p})")

    parts = eqtls["variantId"].str.split("_", expand=True)
    table = pd.DataFrame({
        "gene_id": eqtls["gencodeId"].str.split(".").str[0], "symbol": eqtls["geneSymbol"],
        "tissue": eqtls["tissueSiteDetailId"], "rank": eqtls["rank"], "rsid": eqtls["snpId"],
        "variant_b38": eqtls["variantId"], "chrom": parts[0].str.replace("chr", "", regex=False),
        "ref": parts[2], "alt": parts[3], "maf": eqtls["maf"], "pvalue": eqtls["pValue"], "nes": eqtls["nes"],
    })
    table = table[table["chrom"].isin([str(c) for c in range(1, 23)])]
    table.insert(7, "pos_b37", add_grch37(table))
    table = table.sort_values("pvalue").drop_duplicates(["gene_id", "variant_b38"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output}: {len(table)} gene-variant pairs, {table['gene_id'].nunique()} genes, "
          f"{table['variant_b38'].nunique()} variants ({table['pos_b37'].notna().sum()} pairs with a GRCh37 position)")


if __name__ == "__main__":
    main()
