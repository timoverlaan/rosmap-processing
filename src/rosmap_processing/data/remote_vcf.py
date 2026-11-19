"""Selected sites of a bgzipped, tabix-indexed VCF read through byte ranges
(a local file, or an HTTPS URL such as a Synapse presigned URL), without
downloading the whole file.

The tabix linear index gives, per 16 kb window of a sequence, the virtual
offset of the first record there; from that offset only the BGZF blocks up to
the wanted positions are read. For a joint-called WGS VCF of ~2,000 samples a
site costs roughly the records of its window (tens of kB compressed each).
No htslib needed: BGZF blocks are plain deflate members with a size field.
"""

import gzip
import http.client
import struct
import time
import urllib.parse
import zlib
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

TBI_SHIFT = 14  # tabix linear-index windows of 16 kb


def read_tbi(data: bytes) -> Dict[str, List[int]]:
    """A .tbi index -> {sequence name: linear index (virtual offsets per window)}."""
    raw = gzip.decompress(data)
    if raw[:4] != b"TBI\x01":
        raise ValueError("not a tabix index")
    n_ref, _fmt, _seq, _beg, _end, _meta, _skip, l_nm = struct.unpack_from("<8i", raw, 4)
    names = [n.decode() for n in raw[36:36 + l_nm].split(b"\x00")[:n_ref]]
    pos = 36 + l_nm
    linear = {}
    for name in names:
        (n_bin,) = struct.unpack_from("<i", raw, pos)
        pos += 4
        for _ in range(n_bin):
            _bin, n_chunk = struct.unpack_from("<Ii", raw, pos)
            pos += 8 + 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", raw, pos)
        pos += 4
        linear[name] = list(struct.unpack_from(f"<{n_intv}Q", raw, pos))
        pos += 8 * n_intv
    return linear


def _block_size(buf: bytes, i: int) -> int:
    if buf[i:i + 4] != b"\x1f\x8b\x08\x04":
        raise ValueError("not a BGZF block")
    (xlen,) = struct.unpack_from("<H", buf, i + 10)
    j, end = i + 12, i + 12 + xlen
    while j < end:
        (slen,) = struct.unpack_from("<H", buf, j + 2)
        if buf[j] == 66 and buf[j + 1] == 67:
            return struct.unpack_from("<H", buf, j + 4)[0] + 1
        j += 4 + slen
    raise ValueError("BGZF block without a size field")


def _inflate(block: bytes) -> bytes:
    (xlen,) = struct.unpack_from("<H", block, 10)
    return zlib.decompress(block[12 + xlen:-8], -15)


class FileRanges:
    """Byte ranges of a local file."""

    def __init__(self, path):
        self.handle = open(path, "rb")
        self.handle.seek(0, 2)
        self.size = self.handle.tell()
        self.bytes_read = 0

    def read(self, start: int, length: int) -> bytes:
        self.handle.seek(start)
        data = self.handle.read(length)
        self.bytes_read += len(data)
        return data


class HttpRanges:
    """Byte ranges of a URL over one kept-alive connection (a range costs a
    round trip, not a new TLS handshake). url() gives a fresh URL (presigned
    URLs expire); it is called again after a failed request. A server that
    ignores the Range header is an error and the connection is dropped unread,
    so a whole file is never streamed by accident."""

    def __init__(self, url: Callable[[], str], size: int, retries: int = 6, timeout: float = 180):
        self.url, self.size, self.retries, self.timeout = url, size, retries, timeout
        self._conn, self.bytes_read, self.requests = None, 0, 0
        self._open(url())

    def _open(self, current: str) -> None:
        parts = urllib.parse.urlsplit(current)
        self._target = parts.path + (f"?{parts.query}" if parts.query else "")
        connection = http.client.HTTPSConnection if parts.scheme == "https" else http.client.HTTPConnection
        self._close()
        self._conn = connection(parts.netloc, timeout=self.timeout)

    def _close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def read(self, start: int, length: int) -> bytes:
        end = min(start + length, self.size) - 1
        if end < start:
            return b""
        for attempt in range(self.retries):
            try:
                self._conn.request("GET", self._target, headers={"Range": f"bytes={start}-{end}"})
                response = self._conn.getresponse()
                if response.status != 206:
                    self._close()
                    if response.status == 200:
                        raise RuntimeError("the server ignored the byte range (HTTP 200)")
                    raise ConnectionError(f"HTTP {response.status}")
                data = response.read()
                self.bytes_read += len(data)
                self.requests += 1
                return data
            except (http.client.HTTPException, ConnectionError, TimeoutError, OSError):
                if attempt == self.retries - 1:
                    raise
                time.sleep(2 ** attempt)
                self._open(self.url())
        raise RuntimeError("unreachable")


class BgzfLines:
    """Lines of a BGZF file from a virtual offset, read block by block."""

    def __init__(self, ranges, chunk: int = 1 << 18):
        self.ranges, self.chunk = ranges, chunk
        self.coff, self.buf, self.text = 0, b"", b""

    def seek(self, voffset: int) -> None:
        self.coff, self.buf, self.text = voffset >> 16, b"", b""
        if self._next_block():
            self.text = self.text[voffset & 0xFFFF:]

    def _more(self) -> bool:
        start = self.coff + len(self.buf)
        if start >= self.ranges.size:
            return False
        self.buf += self.ranges.read(start, self.chunk)
        return True

    def _next_block(self) -> bool:
        while len(self.buf) < 18:
            if not self._more():
                return False
        size = _block_size(self.buf, 0)
        while len(self.buf) < size:
            if not self._more():
                raise ValueError("truncated BGZF block")
        block, self.buf = self.buf[:size], self.buf[size:]
        self.coff += size
        data = _inflate(block)
        self.text += data
        return True

    def lines(self) -> Iterator[bytes]:
        while True:
            i = self.text.find(b"\n")
            while i < 0:
                if not self._next_block():
                    if self.text:
                        line, self.text = self.text, b""
                        yield line
                    return
                i = self.text.find(b"\n")
            line, self.text = self.text[:i], self.text[i + 1:]
            yield line


def vcf_samples(ranges) -> List[str]:
    """The sample IDs of a BGZF VCF (its #CHROM line)."""
    stream = BgzfLines(ranges)
    stream.seek(0)
    for line in stream.lines():
        if line.startswith(b"#CHROM"):
            return [s.decode() for s in line.split(b"\t")[9:]]
        if not line.startswith(b"#"):
            break
    raise ValueError("no #CHROM line in the VCF header")


def _record_after(ranges, coff: int, probe: int = 1 << 16) -> Optional[Tuple[int, int, str, int]]:
    """The first record that starts in a BGZF block at or after compressed
    offset `coff`: (block offset, offset of the line in the block, chrom,
    position); None at the end of the file. Block starts are found by their
    header and checked by inflating (the size in the trailer must agree)."""
    buf, base = b"", coff
    while True:
        while len(buf) < 18:
            more = ranges.read(base + len(buf), probe)
            if not more:
                return None
            buf += more
        k = buf.find(b"\x1f\x8b\x08\x04")
        if k < 0 or len(buf) - k < 18:
            keep = max(len(buf) - 17, 0) if k < 0 else k
            base, buf = base + keep, buf[keep:]
            more = ranges.read(base + len(buf), probe)
            if not more:
                return None
            buf += more
            continue
        try:
            size = _block_size(buf, k)
        except ValueError:
            base, buf = base + k + 1, buf[k + 1:]
            continue
        while len(buf) - k < size:
            more = ranges.read(base + len(buf), probe)
            if not more:
                return None
            buf += more
        block = buf[k:k + size]
        try:
            data = _inflate(block)
            valid = len(data) == struct.unpack_from("<I", block, size - 4)[0]
        except (zlib.error, struct.error):
            valid = False
        if not valid:
            base, buf = base + k + 1, buf[k + 1:]
            continue
        # A line that starts in this block (after its first newline)
        nl = data.find(b"\n")
        while nl >= 0 and nl + 1 < len(data):
            head = data[nl + 1:nl + 200]
            if not head.startswith(b"#") and head.count(b"\t") >= 2:
                chrom, pos = head.split(b"\t", 2)[:2]
                return base + k, nl + 1, chrom.decode(), int(pos)
            if head.count(b"\t") < 2:
                break  # the line's start continues in the next block: try the next block
            nl = data.find(b"\n", nl + 1)
        base, buf = base + k + size, buf[k + size:]


def _start_before(ranges, index: List[int], chrom: str, p: int, close: int = 1 << 17) -> int:
    """A virtual offset of a record start before position p, near it: from the
    linear index's window offset, bisected over the compressed blocks of the
    window (the next window's offset bounds it) until within `close` bytes."""
    window = (p - 1) >> TBI_SHIFT
    best = index[window]
    lo = best >> 16
    later = [v for v in index[window + 1:window + 64] if v > best]
    hi = (later[0] >> 16) if later else ranges.size
    while hi - lo > close:
        mid = (lo + hi) // 2
        found = _record_after(ranges, mid)
        if found is None or found[0] >= hi or found[2] != chrom or found[3] >= p:
            hi = mid
        else:
            best, lo = (found[0] << 16) | found[1], found[0]
    return best


def fetch_records(ranges, linear: Dict[str, List[int]], chrom: str, positions: Iterable[int],
                  jump: int = 1 << 18, chunk: int = 1 << 16, bisect: bool = True,
                  close: int = 1 << 17) -> Iterator[Tuple[int, List[bytes]]]:
    """All records of `chrom` at the given 1-based positions, in order:
    (position, the record's tab-separated fields). Streams on between nearby
    positions; seeks when the next one's window starts more than `jump`
    compressed bytes ahead, to a record start found by bisecting the window
    (`bisect`) or to the window start. Reads `chunk` bytes per range request."""
    index = linear.get(chrom)
    if not index:
        return
    stream = BgzfLines(ranges, chunk=chunk)
    lines: Optional[Iterator[bytes]] = None
    pending: Optional[Tuple[int, bytes]] = None  # a line read past the previous position
    for p in sorted(set(int(x) for x in positions)):
        window = (p - 1) >> TBI_SHIFT
        if window >= len(index):
            continue
        if pending is not None and pending[0] > p:
            continue  # the stream is past p: no record there
        target = index[window]
        if lines is None or (target >> 16) > stream.coff + jump:
            stream.seek(_start_before(ranges, index, chrom, p, close) if bisect else target)
            lines, pending = stream.lines(), None
        while True:
            if pending is not None:
                pos, line = pending
                pending = None
            else:
                line = next(lines, None)
                if line is None:
                    lines = None
                    break
                if line.startswith(b"#"):
                    continue
                head = line.split(b"\t", 2)
                if head[0].decode() != chrom:
                    pending = (1 << 62, line)
                    break
                pos = int(head[1])
            if pos < p:
                continue
            if pos == p:
                yield pos, line.split(b"\t")
                continue
            pending = (pos, line)
            break


def _genotypes(fields: Sequence[bytes]) -> List[List[bytes]]:
    gt = fields[8].split(b":").index(b"GT")
    return [sample.split(b":", gt + 1)[gt].replace(b"|", b"/").split(b"/") for sample in fields[9:]]


def allele_counts(fields: Sequence[bytes], allele: int) -> np.ndarray:
    """Per sample, the number of copies of allele `allele` (0 = REF, 1.. = ALT)
    from the GT field; NaN when missing."""
    code = str(allele).encode()
    return np.array([np.nan if b"." in calls else sum(c == code for c in calls) for calls in _genotypes(fields)])


def heterozygous(fields: Sequence[bytes]) -> np.ndarray:
    """Per sample, 1 for a heterozygous diploid call, 0 for a homozygous or
    haploid one, NaN when missing."""
    return np.array([np.nan if b"." in calls else float(len(calls) == 2 and calls[0] != calls[1])
                     for calls in _genotypes(fields)])
