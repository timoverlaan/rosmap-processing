"""Size, MD5 and S3 multipart ETags of a file, read once.

An S3 object uploaded in parts has the ETag md5(md5(part 1) ... md5(part n))-n,
so a local copy can be checked against an object version's listed ETag when
the part size is known; the usual ones are tried.

    python scripts/s3_etag.py FILE [--expected-etag ETAG] [--expected-size BYTES]
"""

import argparse
import hashlib
import os
import sys

MIB = 1024 * 1024


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("--part-size-mib", type=int, nargs="+", default=[5, 8, 16, 32, 64, 100, 128])
    parser.add_argument("--expected-etag", default=None)
    parser.add_argument("--expected-size", type=int, default=None)
    args = parser.parse_args()

    size = os.path.getsize(args.path)
    print(f"{args.path}: {size} bytes")
    full = hashlib.md5()
    parts = {p: [] for p in args.part_size_mib}
    current = {p: hashlib.md5() for p in args.part_size_mib}
    filled = {p: 0 for p in args.part_size_mib}
    with open(args.path, "rb") as f:
        while True:
            block = f.read(MIB)
            if not block:
                break
            full.update(block)
            for p in args.part_size_mib:
                current[p].update(block)
                filled[p] += len(block)
                if filled[p] == p * MIB:
                    parts[p].append(current[p].digest())
                    current[p], filled[p] = hashlib.md5(), 0
    print(f"MD5 {full.hexdigest()}")
    etags = {}
    for p in args.part_size_mib:
        digests = parts[p] + ([current[p].digest()] if filled[p] else [])
        etags[p] = f"{hashlib.md5(b''.join(digests)).hexdigest()}-{len(digests)}"
        print(f"ETag with {p} MiB parts: {etags[p]}")
    ok = True
    if args.expected_size is not None:
        ok &= size == args.expected_size
        print(f"size {'matches' if size == args.expected_size else 'DIFFERS FROM'} {args.expected_size}")
    if args.expected_etag:
        expected = args.expected_etag.strip('"')
        hits = [p for p, e in etags.items() if e == expected]
        ok &= bool(hits)
        print(f"ETag {expected}: " + (f"matches with {hits[0]} MiB parts" if hits else "NO MATCH"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
