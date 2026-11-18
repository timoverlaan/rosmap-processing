"""Dosages of an imputed genotype set at the variants of an eQTL table, for
`pipeline qc-check genotype`.

The input is PLINK 1 dosage files (format=1, no header: SNP A1 A2 then one
dosage per sample in .fam order), one per chromosome, as for ROSMAP's
Affymetrix 6 imputation (GRCh37; SNP IDs are rsIDs or chrN:pos). A variant is
kept when its ID is an eQTL's rsID or its GRCh37 chrN:pos (scripts/gtex_eqtls.py).
The output is a gzipped TSV: snp, a1, a2, then one column per sample (the
.fam IIDs). It holds genotypes: keep it with the restricted data. Only counts
are printed.

    python scripts/extract_dosages.py EQTLS DOSAGE_PATTERN FAM OUTPUT [--processes 8]

DOSAGE_PATTERN contains {chr}, e.g.
.../amp-ad_rosmap_rush-broad_affymetrixgenechip6_imputed_chr{chr}.dosage.gz
"""

import argparse
import gzip
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd


def scan_chromosome(job):
    chrom, path, keys, n_samples = job
    kept, bad = [], 0
    with gzip.open(path, "rb") as f:
        for line in f:
            snp = line.split(None, 1)[0]
            if snp in keys:
                fields = line.split()
                if len(fields) != 3 + n_samples:
                    bad += 1
                    continue
                kept.append(b"\t".join(fields))
    return chrom, kept, bad


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("eqtls", type=Path)
    parser.add_argument("dosage_pattern")
    parser.add_argument("fam", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()

    eqtls = pd.read_csv(args.eqtls, sep="\t", dtype={"chrom": str, "rsid": str})
    iids = pd.read_csv(args.fam, sep=r"\s+", header=None, dtype=str)[1].tolist()
    jobs = []
    for chrom in [str(c) for c in range(1, 23)]:
        here = eqtls[eqtls["chrom"] == chrom]
        keys = {s.encode() for s in here["rsid"].dropna()}
        keys |= {f"chr{chrom}:{int(p)}".encode() for p in here["pos_b37"].dropna()}
        if keys:
            jobs.append((chrom, args.dosage_pattern.format(chr=chrom), keys, len(iids)))
    results = {}
    with ProcessPoolExecutor(args.processes) as pool:
        for chrom, kept, bad in pool.map(scan_chromosome, jobs):
            results[chrom] = kept
            print(f"chr{chrom}: {len(kept)} variants kept" + (f", {bad} malformed rows skipped" if bad else ""),
                  flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wb") as out:
        out.write(("\t".join(["snp", "a1", "a2"] + iids) + "\n").encode())
        for chrom, *_ in jobs:
            for row in results[chrom]:
                out.write(row + b"\n")
    total = sum(len(v) for v in results.values())
    print(f"wrote {args.output}: {total} variants x {len(iids)} samples "
          f"(the eQTL table has {eqtls['variant_b38'].nunique()} variants)")


if __name__ == "__main__":
    main()
