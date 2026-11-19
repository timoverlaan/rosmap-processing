"""Per sample, the heterozygous share of common chrX SNPs (PLINK chromosome
23, outside the pseudo-autosomal regions, which PLINK codes as 25) from a
PLINK .bed/.bim/.fam: about 0 for men, about 0.3 for women. A DNA-side sex for
the genotype check (does a genotype sample's DNA fit its recorded sex?).

    python scripts/plink_chrx_het.py PREFIX OUTPUT.tsv [--min-maf 0.05]

OUTPUT (restricted): IID, n_sites, n_het, het_rate, fam_sex (the .fam's sex
code: 1 male, 2 female, 0 unknown). Prints counts only.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("prefix", help="PLINK prefix (PREFIX.bed, .bim, .fam)")
    parser.add_argument("output", type=Path)
    parser.add_argument("--min-maf", type=float, default=0.05)
    args = parser.parse_args()

    fam = pd.read_csv(f"{args.prefix}.fam", sep=r"\s+", header=None, dtype=str)
    bim = pd.read_csv(f"{args.prefix}.bim", sep=r"\s+", header=None, usecols=[0], dtype=str)[0]
    n = len(fam)
    row_bytes = (n + 3) // 4
    rows = np.flatnonzero(bim.isin(["23", "X"]).to_numpy())
    with open(f"{args.prefix}.bed", "rb") as f:
        if f.read(3) != b"\x6c\x1b\x01":
            raise ValueError("not a SNP-major PLINK .bed")
        n_sites, n_het = np.zeros(n), np.zeros(n)
        used = 0
        for i in rows:
            f.seek(3 + int(i) * row_bytes)
            packed = np.frombuffer(f.read(row_bytes), dtype=np.uint8)
            codes = np.stack([(packed >> s) & 3 for s in (0, 2, 4, 6)], axis=1).ravel()[:n]
            called = codes != 1  # 00 hom A1, 10 het, 11 hom A2, 01 missing
            het = codes == 2
            a2 = np.where(codes == 3, 2, np.where(het, 1, 0))
            freq = a2[called].sum() / max(2 * called.sum(), 1)
            if not args.min_maf <= freq <= 1 - args.min_maf:
                continue
            n_sites += called
            n_het += het
            used += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = n_het / n_sites
    out = pd.DataFrame({"IID": fam[1], "n_sites": n_sites.astype(int), "n_het": n_het.astype(int),
                        "het_rate": rate, "fam_sex": fam[4]})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False, float_format="%.4f")
    male, female = rate < 0.05, rate >= 0.15
    print(f"{n} samples, {len(rows)} chrX SNPs, {used} common ones used; het rate < 0.05: {int(male.sum())}, "
          f"0.05-0.15: {int((~male & ~female & np.isfinite(rate)).sum())}, >= 0.15: {int(female.sum())}")
    fam_sex = fam[4].to_numpy()
    print(f"vs the .fam sex: DNA male & fam 1: {int((male & (fam_sex == '1')).sum())}, DNA female & fam 2: "
          f"{int((female & (fam_sex == '2')).sum())}, DNA male & fam 2: {int((male & (fam_sex == '2')).sum())}, "
          f"DNA female & fam 1: {int((female & (fam_sex == '1')).sum())}, fam 0: {int((fam_sex == '0').sum())}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
