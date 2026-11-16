"""Which genes make up the difference between our per-cell UMIs and the
authors' (aggregate output: gene names and counts of cells, no cell or donor).

    python scripts/totals_gap.py config/data_qc/<name>.yaml [--cells 20000]

Needs the spec's qc-scan. On a random sample of cells with the authors'
totals, a gene can only be part of the gap if it never has more UMIs in a
cell than that cell's gap (ours - theirs). The candidates' UMIs are then
summed: if the sum equals the gap in (nearly) every cell, the authors' totals
leave out exactly these genes (up to genes too rare to tell apart).
"""

import argparse

import h5py
import numpy as np
import scipy.sparse as sp

from rosmap_processing.pipelines import batch_effects as be
from rosmap_processing.pipelines.data_qc import Scan, load_qc_spec

try:
    from anndata.io import sparse_dataset
except ImportError:  # anndata < 0.11
    from anndata.experimental import sparse_dataset


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("spec")
    parser.add_argument("--cells", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    spec = load_qc_spec(args.spec)
    scan = Scan(spec)
    ours, theirs = scan.cells["umis"], scan.author("umis")
    if theirs is None:
        raise SystemExit("the spec has no author UMIs")
    ok = np.flatnonzero(np.isfinite(theirs))
    rng = np.random.default_rng(args.seed)
    rows = np.sort(rng.choice(ok, size=min(args.cells, len(ok)), replace=False))
    gap = ours[rows] - theirs[rows]
    print(f"{len(rows)} sampled cells; gap (ours - authors) quantiles 1/50/99%: "
          f"{np.quantile(gap, [0.01, 0.5, 0.99])}; negative in {int((gap < 0).sum())} cells")

    cfg = spec["cfg"]
    with h5py.File(be._first_existing(cfg["counts"]), "r") as f:
        matrix = sparse_dataset(be._matrix_node(f, cfg.get("layer")))
        X = sp.vstack([sp.csr_matrix(matrix[r:r + 1]) for r in rows]).tocsr()
    per_row = np.diff(X.indptr)
    too_many = X.data > np.repeat(gap, per_row)
    violates = np.zeros(X.shape[1], bool)
    violates[np.unique(X.indices[too_many])] = True
    detected = np.zeros(X.shape[1], bool)
    detected[np.unique(X.indices)] = True
    candidates = np.flatnonzero(detected & ~violates)
    total = np.asarray(X[:, candidates].sum(axis=1)).ravel()
    print(f"{len(candidates)} of {int(detected.sum())} detected genes never exceed a cell's gap")
    print(f"sum over them equals the gap in {100 * np.mean(total == gap):.2f}% of cells, "
          f"is below it in {100 * np.mean(total < gap):.2f}%, above in {100 * np.mean(total > gap):.2f}%")

    genes = scan.genes
    umis = np.asarray(X[:, candidates].sum(axis=0)).ravel()
    order = np.argsort(-umis)
    share = umis[order].cumsum() / max(umis.sum(), 1)
    n80 = int(np.searchsorted(share, 0.8)) + 1
    print(f"{n80} candidates carry 80% of their UMIs; the top 40:")
    for i in order[:40]:
        g = candidates[i]
        print(f"  {genes['symbol'].iloc[g]:<18} {genes['gene_id'].iloc[g]:<17} {genes['chrom'].iloc[g]:<6} "
              f"{int(umis[i])} UMIs in the sample")
    symbols = genes["symbol"].iloc[candidates].astype(str)
    patterns = {"MT-": symbols.str.startswith("MT-"), "RP[SL]": symbols.str.match(r"^RP[SL]\d"),
                "A?.../ENSG (unnamed)": symbols.str.match(r"^(A[CLP]\d{6}|ENSG\d+|LINC|Z\d{5})"),
                "other": None}
    rest = np.ones(len(symbols), bool)
    for name, hit in patterns.items():
        hit = rest if hit is None else hit.to_numpy() & rest
        print(f"  {name}: {int(hit.sum())} candidate genes, {int(umis[hit].sum())} UMIs")
        rest &= ~hit
    if scan.author("genes") is not None:
        gene_gap = scan.cells["genes"][rows] - scan.author("genes")[rows]
        n_cand = np.diff((X[:, candidates] > 0).tocsr().indptr)
        print(f"genes detected among the candidates equal the genes gap in {100 * np.mean(n_cand == gene_gap):.2f}% "
              f"of cells")


if __name__ == "__main__":
    main()
