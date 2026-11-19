"""remote_vcf: sites of a BGZF VCF through its tabix linear index."""

import gzip
import struct
import zlib

import numpy as np
import pytest

from rosmap_processing.data import remote_vcf

SAMPLES = ["s1", "s2", "s3"]


def _block(data: bytes) -> bytes:
    deflate = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = deflate.compress(data) + deflate.flush()
    bsize = 18 + len(payload) + 8 - 1
    header = (b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6) + b"BC"
              + struct.pack("<HH", 2, bsize))
    return header + payload + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))


def _write_indexed_vcf(path, records, block_size=150):
    """A BGZF VCF (small blocks, so records span blocks) and its tabix
    linear index, as htslib lays them out (empty windows take the previous
    window's offset)."""
    header = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
              + "\t".join(SAMPLES) + "\n")
    text = header.encode()
    starts = []
    for chrom, pos, ref, alt, gts in records:
        starts.append((pos, len(text)))
        text += f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:DP\t".encode()
        text += "\t".join(f"{g}:9" for g in gts).encode() + b"\n"
    blob, block_starts = b"", []
    for i in range(0, len(text), block_size):
        block_starts.append((i, len(blob)))
        blob += _block(text[i:i + block_size])
    blob += _block(b"")  # EOF marker
    path.write_bytes(blob)

    def voffset(u):
        ustart, cstart = max(b for b in block_starts if b[0] <= u)
        return (cstart << 16) | (u - ustart)

    n_win = (max(p for p, _ in starts) - 1 >> remote_vcf.TBI_SHIFT) + 1
    linear = [None] * n_win
    for pos, u in starts:
        w = (pos - 1) >> remote_vcf.TBI_SHIFT
        if linear[w] is None:
            linear[w] = voffset(u)
    last = voffset(starts[0][1])
    for w in range(n_win):
        linear[w] = last = linear[w] if linear[w] is not None else last
    name = b"1\x00"
    tbi = b"TBI\x01" + struct.pack("<8i", 1, 2, 1, 2, 0, ord("#"), 0, len(name)) + name
    tbi += struct.pack("<i", 0) + struct.pack("<i", n_win) + struct.pack(f"<{n_win}Q", *linear)
    return gzip.compress(tbi)


RECORDS = [("1", 100, "A", "G", ["0/0", "0/1", "1/1"]),
           ("1", 5000, "C", "T", ["0|1", "./.", "1|1"]),
           ("1", 17000, "G", "A", ["0/0", "0/0", "0/1"]),
           ("1", 17000, "G", "GA", ["0/1", "0/0", "0/0"]),
           ("1", 40000, "T", "C,G", ["1/2", "2/2", "0/1"]),
           ("1", 90000, "A", "C", ["1", "0", "./."]),
           ("1", 90001, "C", "T", ["0/1", "0/1", "0/1"])]


@pytest.mark.parametrize("jump", [0, 1 << 20])
def test_fetch_records(tmp_path, jump):
    tbi = _write_indexed_vcf(tmp_path / "x.vcf.gz", RECORDS)
    linear = remote_vcf.read_tbi(tbi)
    ranges = remote_vcf.FileRanges(tmp_path / "x.vcf.gz")
    assert remote_vcf.vcf_samples(ranges) == SAMPLES
    got = list(remote_vcf.fetch_records(ranges, linear, "1", [90001, 17000, 100, 6000, 40000, 200000], jump=jump))
    assert [(p, f[4].decode()) for p, f in got] == [(100, "G"), (17000, "A"), (17000, "GA"), (40000, "C,G"),
                                                     (90001, "T")]
    assert list(remote_vcf.fetch_records(ranges, linear, "2", [100])) == []


def test_allele_counts_and_heterozygous(tmp_path):
    tbi = _write_indexed_vcf(tmp_path / "x.vcf.gz", RECORDS)
    ranges = remote_vcf.FileRanges(tmp_path / "x.vcf.gz")
    records = dict((p, f) for p, f in remote_vcf.fetch_records(ranges, remote_vcf.read_tbi(tbi), "1",
                                                                [5000, 40000, 90000]))
    np.testing.assert_array_equal(remote_vcf.allele_counts(records[5000], 1), [1, np.nan, 2])
    np.testing.assert_array_equal(remote_vcf.allele_counts(records[40000], 2), [1, 2, 0])
    np.testing.assert_array_equal(remote_vcf.heterozygous(records[90000]), [0, 0, np.nan])  # haploid calls
    np.testing.assert_array_equal(remote_vcf.heterozygous(records[40000]), [1, 0, 1])
