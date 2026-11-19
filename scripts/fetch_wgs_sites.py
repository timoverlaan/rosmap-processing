"""Genotypes of WGS samples at the eQTL variants, and chrX heterozygosity
(DNA sex), from a Synapse folder of per-chromosome tabix-indexed joint-call
VCFs. Only the records near the wanted sites are transferred, through byte
ranges (rosmap_processing.data.remote_vcf), not the files themselves.

    python scripts/fetch_wgs_sites.py SYNAPSE_FOLDER EQTLS OUTPUT --build 37|38 \
        [--threads 4] [--max-gb 15] [--chromosomes 21,22]

OUTPUT (restricted; gzipped TSV): snp (chrN:pos in that build), a1 (the
counted allele, an ALT), a2 (REF), then one column per VCF sample: the count
of a1 (empty when missing), like scripts/extract_dosages.py writes.
OUTPUT stem + .chrx.tsv: per sample, the heterozygous share of common
biallelic chrX SNPs outside the pseudo-autosomal regions (about 0 for men).
Needs token.txt (Synapse) in the working directory. Prints counts and the
bytes transferred only: no sample IDs, no URLs.
"""

import argparse
import gzip
import json
import logging
import re
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

logging.disable(logging.CRITICAL)  # synapseclient logs presigned URLs at debug level

from rosmap_processing.data.download import _login, _token  # noqa: E402
from rosmap_processing.data.remote_vcf import (BgzfLines, HttpRanges, allele_counts, fetch_records,  # noqa: E402
                                               heterozygous, read_tbi, vcf_samples)

COMPLEMENT = str.maketrans("ACGT", "TGCA")
VCF_NAME = re.compile(r"[._](?:chr)?([0-9]{1,2}|X|Y)\.(?:recalibrated_variants\.)?vcf\.gz$")
CHRX_FROM, CHRX_TO = 5_000_000, 150_000_000  # inside the non-PAR region in GRCh37 and GRCh38


class Budget:
    def __init__(self, max_bytes):
        self.max_bytes, self.used, self.lock = max_bytes, 0, threading.Lock()

    def add(self, n):
        with self.lock:
            self.used += n
            if self.used > self.max_bytes:
                raise RuntimeError(f"transfer budget of {self.max_bytes / 1e9:.1f} GB exceeded")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", help="Synapse folder with the per-chromosome VCFs and their .tbi")
    parser.add_argument("eqtls", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--build", type=int, choices=(37, 38), required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-gb", type=float, default=15.0, help="stop when this much has been transferred")
    parser.add_argument("--chromosomes", help="comma-separated subset (default: 1-22 and X)")
    parser.add_argument("--chrx-windows", type=int, default=60)
    parser.add_argument("--chrx-records", type=int, default=60, help="records read per chrX window")
    args = parser.parse_args()

    syn = _login(_token())
    syn_lock = threading.Lock()
    files = {c["name"]: c["id"] for c in syn.getChildren(args.folder, includeTypes=["file"])}
    vcfs = {}
    for name, entity in files.items():
        m = VCF_NAME.search(name)
        if m and "annotated" not in name.lower() and name + ".tbi" in files:
            vcfs[m.group(1)] = (entity, files[name + ".tbi"])
    wanted = args.chromosomes.split(",") if args.chromosomes else [str(c) for c in range(1, 23)] + ["X"]
    missing = [c for c in wanted if c not in vcfs]
    print(f"{args.folder}: indexed VCFs for {len(vcfs)} chromosomes; missing of those wanted: {missing or 'none'}")

    eqtl = pd.read_csv(args.eqtls, sep="\t", dtype={"chrom": str})
    if args.build == 37:
        pos = pd.to_numeric(eqtl["pos_b37"], errors="coerce")
    else:
        pos = pd.to_numeric(eqtl["variant_b38"].str.split("_").str[1], errors="coerce")
    eqtl = eqtl.assign(pos=pos).dropna(subset=["pos"]).drop_duplicates(["chrom", "pos"])
    alleles = {(c, int(p)): {str(r), str(a)} for c, p, r, a in zip(eqtl["chrom"], eqtl["pos"], eqtl["ref"], eqtl["alt"])}
    budget = Budget(args.max_gb * 1e9)
    tmp = Path(tempfile.mkdtemp())

    def open_vcf(chrom):
        entity, tbi_entity = vcfs[chrom]
        with syn_lock:
            tbi_path = syn.get(tbi_entity, downloadLocation=str(tmp / chrom), ifcollision="overwrite.local").path
            info = syn.get(entity, downloadFile=False)
        linear = read_tbi(Path(tbi_path).read_bytes())
        seq = chrom if chrom in linear else f"chr{chrom}"

        def url():
            body = {"requestedFiles": [{"fileHandleId": info.dataFileHandleId, "associateObjectId": info.id,
                                        "associateObjectType": "FileEntity"}],
                    "includePreSignedURLs": True, "includeFileHandles": False}
            with syn_lock:
                result = syn.restPOST("/fileHandle/batch", body=json.dumps(body), endpoint=syn.fileHandleEndpoint)
            return result["requestedFiles"][0]["preSignedURL"]

        ranges = HttpRanges(url, int(info._file_handle["contentSize"]))
        return ranges, linear, seq

    def metered(ranges):
        read = ranges.read

        def counted(start, length):
            data = read(start, length)
            budget.add(len(data))
            return data
        ranges.read = counted
        return ranges

    def eqtl_sites(chrom):
        ranges, linear, seq = open_vcf(chrom)
        ranges = metered(ranges)
        samples = vcf_samples(ranges)
        targets = sorted(p for c, p in alleles if c == chrom)
        rows, n_filtered = [], 0
        for p, fields in fetch_records(ranges, linear, seq, targets):
            if (chrom, p) in {(r[0], r[1]) for r in rows[-3:]}:
                continue
            ref, alts = fields[3].decode(), fields[4].decode().split(",")
            want = alleles[(chrom, p)]
            flipped = {a.translate(COMPLEMENT) for a in want}
            for k, alt in enumerate(alts, start=1):
                if len(ref) == 1 and len(alt) == 1 and {ref, alt} in (want, flipped):
                    n_filtered += fields[6] not in (b"PASS", b".")
                    rows.append((chrom, p, alt, ref, allele_counts(fields, k)))
                    break
        print(f"chr{chrom}: {len(rows)} of {len(targets)} eQTL sites found ({n_filtered} not PASS); "
              f"{ranges.bytes_read / 1e6:.0f} MB read")
        return samples, rows

    def chrx_het():
        ranges, linear, seq = open_vcf("X")
        ranges = metered(ranges)
        samples = vcf_samples(ranges)
        index = linear[seq]
        n_sites = np.zeros(len(samples))
        n_het = np.zeros(len(samples))
        used = 0
        stream = BgzfLines(ranges)
        for start in np.linspace(CHRX_FROM, CHRX_TO, args.chrx_windows).astype(int):
            window = (start - 1) >> 14
            if window >= len(index):
                continue
            stream.seek(index[window])
            taken = 0
            for line in stream.lines():
                if line.startswith(b"#"):
                    continue
                fields = line.split(b"\t")
                if int(fields[1]) < start:
                    continue
                taken += 1
                if taken > args.chrx_records:
                    break
                ref, alt = fields[3], fields[4]
                if len(ref) != 1 or len(alt) != 1 or fields[6] not in (b"PASS", b"."):
                    continue
                counts = allele_counts(fields, 1)
                freq = np.nanmean(counts) / 2
                if not 0.05 <= freq <= 0.95:
                    continue
                het = heterozygous(fields)
                ok = np.isfinite(het)
                n_sites += ok
                n_het += np.where(ok, het, 0)
                used += 1
        print(f"chrX: {used} common biallelic SNPs used for heterozygosity; {ranges.bytes_read / 1e6:.0f} MB read")
        return samples, n_sites, n_het

    autosomes = [c for c in wanted if c != "X" and c in vcfs]
    with ThreadPoolExecutor(args.threads) as pool:
        x_future = pool.submit(chrx_het) if "X" in wanted and "X" in vcfs else None
        results = list(pool.map(eqtl_sites, autosomes))
        x_result = x_future.result() if x_future else None

    samples = results[0][0]
    for s, _ in results[1:]:
        if s != samples:
            raise ValueError("the per-chromosome VCFs have different samples")
    rows = [r for _, rs in results for r in rs]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt") as out:
        out.write("\t".join(["snp", "a1", "a2"] + samples) + "\n")
        for chrom, p, a1, a2, d in rows:
            values = ["" if np.isnan(v) else str(int(v)) for v in d]
            out.write("\t".join([f"chr{chrom}:{p}", a1, a2] + values) + "\n")
    missing_share = np.mean([np.isnan(d).mean() for *_, d in rows]) if rows else np.nan
    print(f"wrote {args.output}: {len(rows)} sites x {len(samples)} samples (mean missing share "
          f"{missing_share:.3f})")
    if x_result:
        xs, n_sites, n_het = x_result
        if xs != samples:
            raise ValueError("chrX VCF has different samples")
        with np.errstate(invalid="ignore", divide="ignore"):
            rate = n_het / n_sites
        x_path = args.output.with_name(args.output.name.replace(".tsv.gz", "") + ".chrx.tsv")
        pd.DataFrame({"sample": samples, "n_sites": n_sites.astype(int), "n_het": n_het.astype(int),
                      "het_rate": rate}).to_csv(x_path, sep="\t", index=False, float_format="%.4f")
        r = rate[np.isfinite(rate)]
        print(f"wrote {x_path}: chrX heterozygosity per sample; median sites {np.median(n_sites):.0f}; "
              f"het rate < 0.03: {int((r < 0.03).sum())}, 0.03-0.10: {int(((r >= 0.03) & (r < 0.10)).sum())}, "
              f">= 0.10: {int((r >= 0.10).sum())}")
    print(f"transferred {budget.used / 1e9:.2f} GB in all")


if __name__ == "__main__":
    main()
