"""Helpers for bumbl index files (``.bumbl.bi``).

Schemas: see bumbl_index.md for description
"""

from dataclasses import dataclass
import hashlib
import os
import struct
from typing import Protocol
import urllib.error
import urllib.request

import numpy as np
import mumemto.utils as utils 

_U64 = struct.Struct("<Q")
_U64x2 = struct.Struct("<QQ")

_HEADER_MAGIC = b"bumblbi"
_HEADER_SIZE = 4
_CHECKSUM_SAMPLE_SIZE = 1000

# Empty-bin genomic span sentinels (per-bin MIN/MAX in index files)
_BIN_SPAN_EMPTY_MIN = np.uint64(np.iinfo(np.uint64).max)
_BIN_SPAN_EMPTY_MAX = np.uint64(0)


def _bin_span_for_row_interval(mums, seq_idx: int, lo: int, hi: int):
    """Genomic [min_start, max_end) on column seq_idx for half-open rows [lo, hi)."""
    if lo >= hi:
        return _BIN_SPAN_EMPTY_MIN, _BIN_SPAN_EMPTY_MAX
    col = mums.starts[lo:hi, seq_idx].astype(np.int64, copy=False)
    ends = col + mums.lengths[lo:hi].astype(np.int64, copy=False)
    span_min = np.uint64(col[0])
    span_max = np.uint64(ends.max())
    return span_min, span_max


def _bin_span_for_bin_pairs(mums, seq_idx: int, pairs):
    """
    Genomic min start / max end for union of half-open row ranges ``pairs`` list of (s, e).

    Each range is contiguous rows with non-decreasing starts on seq_idx; min per range
    is O(1) from the first row; max end is one numpy reduction per range.
    """
    if not pairs:
        return _BIN_SPAN_EMPTY_MIN, _BIN_SPAN_EMPTY_MAX
    starts_col = mums.starts[:, seq_idx]
    lens = mums.lengths
    range_mins = []
    range_maxs = []
    for s, e in pairs:
        s, e = int(s), int(e)
        if s >= e:
            continue
        range_mins.append(starts_col[s])
        range_maxs.append(
            np.uint64((starts_col[s:e].astype(np.int64, copy=False) + lens[s:e].astype(np.int64, copy=False)).max())
        )
    if not range_mins:
        return _BIN_SPAN_EMPTY_MIN, _BIN_SPAN_EMPTY_MAX
    span_min = np.uint64(min(int(x) for x in range_mins))
    span_max = max(range_maxs)
    return span_min, span_max


def _parse_bins_arg(bins):
    if type(bins) == int:
        return bins, bins
    if isinstance(bins, (tuple, list)) and len(bins) == 2:
        return int(bins[0]), int(bins[1])
    raise ValueError("Invalid type for bins. Must be int or (start, end).")


def _build_header_word(index_format):
    if not (0 <= index_format <= 0xFF):
        raise ValueError("index format must fit in one byte")
    return np.uint64(int.from_bytes(_HEADER_MAGIC + bytes([index_format]), "little"))


def _parse_header_word(header_word):
    raw = int(header_word).to_bytes(8, "little")
    if raw[:7] != _HEADER_MAGIC:
        raise AssertionError("invalid bumbl index header")
    return raw[7]


def bumbl_lengths_checksum(lengths):
    """
    Compute the checksum stored in `.bumbl.bi` files.

    Checksum is SHA-256 of the first 1000 lengths values (as uint64),
    truncated to 8 bytes (64 bits). The truncated bytes are interpreted as a
    big-endian uint64.
    """
    sample = np.asarray(lengths[:_CHECKSUM_SAMPLE_SIZE], dtype=np.uint64)
    digest = hashlib.sha256(sample.tobytes(order="C")).digest()
    checksum = np.uint64(int.from_bytes(digest[:8], "big"))
    return checksum

def checksum_from_bumbl(bumbl_path):
    """
    Compute the `(num_seqs, checksum)` pair for a `.bumbl` file without loading the whole file.

    Reads the first min(1000, n_mums) lengths (uint32) after the bumbl header.
    """
    with open(bumbl_path, "rb") as fin:
        # flags (uint16), then n_seqs (uint64), n_mums (uint64)
        fin.read(2)  # flags
        n_seqs = _U64.unpack(fin.read(8))[0]
        n_mums = _U64.unpack(fin.read(8))[0]
        sample_n = min(_CHECKSUM_SAMPLE_SIZE, int(n_mums))
        lengths32 = np.fromfile(fin, dtype=np.uint32, count=sample_n)
    return n_seqs, bumbl_lengths_checksum(lengths32.astype(np.uint64, copy=False))


def verify_bumbl_sorted_column(bumbl_path, seq_idx, chunk_rows=65536):
    """
    Stream-check that ``starts[:, seq_idx]`` is non-decreasing (sorted for interval queries).

    Returns True if sorted; raises ValueError if ``seq_idx`` is out of range; returns False
    if any adjacent pair violates monotonicity.
    """
    
    prev = None
    for _lengths, starts_col, _strands_col in utils.parse_bumbl_generator(
        bumbl_path,
        seq_idx=int(seq_idx),
        verbose=False,
        chunksize=int(chunk_rows),
        return_chunk=True,
        return_blocks=False,
    ):
        if starts_col.size == 0:
            continue
        if prev is not None and starts_col[0] < prev:
            return False
        if np.any(np.diff(starts_col) < 0):
            return False
        prev = int(starts_col[-1])
    return True


def verify_bumbl_index(bumbl_path, index_path):
    """
    Verify that `index_path` matches `bumbl_path` by comparing the stored
    checksum to a checksum computed from the `.bumbl` file.

    Returns True if it matches; otherwise raises AssertionError.
    """
    with open(index_path, "rb") as fin:
        fin.seek(0)
        _parse_header_word(_U64.unpack(fin.read(8))[0])  # validates magic
        stored_checksum = np.uint64(_U64.unpack(fin.read(8))[0])
        fin.read(8)  # bin_width
        stored_num_seqs = np.uint64(_U64.unpack(fin.read(8))[0])
    expected_num_seqs, expected_checksum = checksum_from_bumbl(bumbl_path)
    if stored_num_seqs != expected_num_seqs:
        raise AssertionError("bumbl/index num_seqs mismatch")
    if stored_checksum != expected_checksum:
        raise AssertionError("bumbl/index checksum mismatch")
    return True

def rle(col):
    n = len(col)
    if n == 0:
        return []
    change = np.r_[True, col[1:] != col[:-1]]
    starts = np.flatnonzero(change)
    vals = col[starts]
    ends = np.r_[starts[1:], n]
    return zip(vals, starts, ends)
    
def generate_ranges(mums, bin_width=1_000_000, verbose=False):
    """Find all ranges of MUMs for each bin in each sequence."""
    rmaps = []
    seq_iter = range(mums.num_seqs)
    if verbose:
        import importlib

        tqdm = importlib.import_module("tqdm.auto").tqdm
        seq_iter = tqdm(seq_iter, desc="index:sequences", unit="seq")
    for i in seq_iter:
        bins = mums.starts[:, i] // bin_width
        rmap = [[] for _ in range(bins.max() + 1)]
        runs = rle(bins)
        for b, s, e in runs:
            rmap[b].append((s,e))
        rmaps.append(rmap)
    return rmaps

def build_bumbl_multiindex(mums, bin_width=1_000_000, verbose=False):
    """Compute bumbl index given a set of MUMs as a MUMdata object."""
    range_maps = generate_ranges(mums, bin_width, verbose=verbose)
    INT_SIZE = 8 # measured in uint64 words
    
    # HEADER
    checksum = bumbl_lengths_checksum(mums.lengths)
    index = [_build_header_word(1), checksum]
    index.append(bin_width)
    index.append(mums.num_seqs)
    
    # DOC_OFFSETS table
    doc_offsets = [0] * (mums.num_seqs + 1)  # offset to each doc
    doc_index = []
    num_bins = [len(r) for r in range_maps]
    for seq_idx in range(mums.num_seqs):
        doc_offsets[seq_idx] = len(doc_index) * INT_SIZE # store the current offset of doc_index as the start of doc seq_idx
        doc_index.append(num_bins[seq_idx]) # number of bins for doc seq_idx
        rmap = range_maps[seq_idx]
        for b in range(len(rmap)):
            smin, smax = _bin_span_for_bin_pairs(mums, seq_idx, rmap[b])
            doc_index.append(smin)
            doc_index.append(smax)
        bin_offsets = []
        flattened_ranges = []
        for pairs in rmap:
            n = len(pairs)
            bin_offsets.append(n * INT_SIZE * 2) # pair of int64s. Offset relative to start of range array
            for pair in pairs:
                flattened_ranges.extend(pair)
        bin_offsets = np.array(bin_offsets)
        bin_offsets = np.insert(np.cumsum(bin_offsets), 0, 0)
        doc_index.extend(bin_offsets)
        doc_index.extend(flattened_ranges)
    doc_offsets[-1] = len(doc_index) * INT_SIZE
    # four header uint64s + doc_offsets
    doc_offsets = np.array(doc_offsets)
    doc_offsets = doc_offsets + (_HEADER_SIZE * INT_SIZE) + (INT_SIZE * len(doc_offsets))
    index.extend(doc_offsets)
    index.extend(doc_index)
    index = np.array(index, dtype=np.uint64)
    return index
    
def build_bumbl_singleindex(mums, seq_idx, bin_width = 1_000_000):
    """Compute bumbl index given a set of MUMs as a MUMdata object."""
    edges = np.arange(0, mums.starts[-1, seq_idx] + bin_width, bin_width)
    offsets = np.searchsorted(
        mums.starts[:, seq_idx], edges, side="left"
    ).astype(np.uint64)
    num_bins = int(offsets.size) - 1
    span_flat = []
    for b in range(num_bins):
        lo, hi = int(offsets[b]), int(offsets[b + 1])
        smin, smax = _bin_span_for_row_interval(mums, seq_idx, lo, hi)
        span_flat.extend((smin, smax))
    checksum = bumbl_lengths_checksum(mums.lengths)
    header = np.array(
        [_build_header_word(0), checksum, bin_width, mums.num_seqs, np.uint64(seq_idx)],
        dtype=np.uint64,
    )
    return np.concatenate([header, offsets, np.array(span_flat, dtype=np.uint64)])

def write_bumbl_index(path, index):
    """Write bumbl index to file"""
    assert index.dtype == np.uint64
    with open(path, 'wb') as f:
        f.write(index.tobytes())


####### Functions to parse a BumblBi index #######


@dataclass(frozen=True, slots=True)
class BumblMultiIndex:
    """In-memory view of one document (one seq_idx) inside a multi-index."""

    seq_idx: int
    bin_width: np.uint64
    num_seqs: np.uint64
    bin_num: np.uint64
    max_bin: int  # largest valid bin index: int(bin_num) - 1
    bin_spans: np.ndarray  # uint64, shape (bin_num, 2): [min_start, max_end] per bin
    boundaries: np.ndarray  # uint64, shape (bin_num+1,), byte offsets into ranges section
    ranges: np.ndarray  # uint64, shape (k,2), (mum_start, mum_end) pairs

    def get_bins(self, bins):
        bin_start, bin_end = _parse_bins_arg(bins)
        bin_num = int(self.bin_num)
        if not (0 <= bin_start <= bin_end < bin_num):
            raise AssertionError("invalid bin")

        off0 = int(self.boundaries[bin_start])
        off1 = int(self.boundaries[bin_end + 1])
        if off1 < off0:
            raise AssertionError("corrupt index (bin offsets)")
        if (off0 % 16) != 0 or (off1 % 16) != 0:
            raise AssertionError("corrupt index (bin offsets not multiple of 2*u64)")
        i0 = off0 // 16
        i1 = off1 // 16
        return self.ranges[i0:i1]

    def bin_is_empty(self, b: int) -> bool:
        """True if bin ``b`` has no indexed row ranges (zero byte span in ranges section)."""
        n = int(self.bin_num)
        bb = int(b)
        if not (0 <= bb < n):
            raise AssertionError("invalid bin")
        return int(self.boundaries[bb]) == int(self.boundaries[bb + 1])

    def contains_left_bound(self, b: int, pos: int) -> bool:
        """
        Whether bin ``b`` can bound a query *start* at ``pos`` using stored spans.

        True iff the bin is non-empty and ``pos >= span_min`` (min MUM start in the bin).
        ``pos > span_max`` is still True (all MUMs lie left of ``pos``).
        """
        if self.bin_is_empty(b):
            return False
        return int(pos) >= int(self.bin_spans[int(b), 0])

    def contains_right_bound(self, b: int, pos: int) -> bool:
        """
        Whether bin ``b`` can bound a query *end* at ``pos`` using stored spans.

        True iff the bin is non-empty and ``pos <= span_max`` (max MUM end in the bin).
        ``pos < span_min`` is still True (all MUMs lie right of ``pos``).
        """
        if self.bin_is_empty(b):
            return False
        return int(pos) <= int(self.bin_spans[int(b), 1])

    def closest_nonzero_bin_left(self, bin_idx: int) -> int | None:
        n = int(self.bin_num)
        b = int(bin_idx)
        if not (0 <= b < n):
            raise AssertionError("invalid bin")
        for i in range(b, -1, -1):
            if not self.bin_is_empty(i):
                return i
        return None

    def closest_nonzero_bin_right(self, bin_idx: int) -> int | None:
        n = int(self.bin_num)
        b = int(bin_idx)
        if not (0 <= b < n):
            raise AssertionError("invalid bin")
        for i in range(b, n):
            if not self.bin_is_empty(i):
                return i
        return None

    def coord_to_bin(self, coord) -> int:
        bw = int(self.bin_width)
        return int(coord) // bw


@dataclass(frozen=True, slots=True)
class BumblSingleIndex:
    """In-memory view of a single-index."""

    bin_width: np.uint64
    num_seqs: np.uint64
    seq_idx: np.uint64
    offsets: np.ndarray  # uint64, shape (bin_num+1,)
    bin_spans: np.ndarray  # uint64, shape (bin_num, 2): [min_start, max_end] per bin
    max_bin: int  # largest valid bin index: offsets.size - 2

    def get_bins(self, bins):
        bin_start, bin_end = _parse_bins_arg(bins)
        offsets = self.offsets
        bin_num = int(offsets.size - 1)
        if not (0 <= bin_start <= bin_end < bin_num):
            raise AssertionError("invalid bin range.")
        starts = offsets[bin_start : bin_end + 1]
        ends = offsets[bin_start + 1 : bin_end + 2]
        return np.stack((starts, ends), axis=1)

    def bin_is_empty(self, b: int) -> bool:
        """True if bin ``b`` has no MUM rows (half-open row interval is empty)."""
        offsets = self.offsets
        n = int(offsets.size - 1)
        bb = int(b)
        if not (0 <= bb < n):
            raise AssertionError("invalid bin")
        return int(offsets[bb]) == int(offsets[bb + 1])

    def contains_left_bound(self, b: int, pos: int) -> bool:
        """
        Whether bin ``b`` can bound a query *start* at ``pos`` using stored spans.

        True iff the bin is non-empty and ``pos >= span_min``.
        ``pos > span_max`` is still True (all MUMs lie left of ``pos``).
        """
        if self.bin_is_empty(b):
            return False
        return int(pos) >= int(self.bin_spans[int(b), 0])

    def contains_right_bound(self, b: int, pos: int) -> bool:
        """
        Whether bin ``b`` can bound a query *end* at ``pos`` using stored spans.

        True iff the bin is non-empty and ``pos <= span_max``.
        ``pos < span_min`` is still True (all MUMs lie right of ``pos``).
        """
        if self.bin_is_empty(b):
            return False
        return int(pos) <= int(self.bin_spans[int(b), 1])

    def closest_nonzero_bin_left(self, bin_idx: int) -> int | None:
        offsets = self.offsets
        n = int(offsets.size - 1)
        b = int(bin_idx)
        if not (0 <= b < n):
            raise AssertionError("invalid bin")
        for i in range(b, -1, -1):
            if not self.bin_is_empty(i):
                return i
        return None

    def closest_nonzero_bin_right(self, bin_idx: int) -> int | None:
        offsets = self.offsets
        n = int(offsets.size - 1)
        b = int(bin_idx)
        if not (0 <= b < n):
            raise AssertionError("invalid bin")
        for i in range(b, n):
            if not self.bin_is_empty(i):
                return i
        return None

    def coord_to_bin(self, coord) -> int:
        bw = int(self.bin_width)
        return int(coord) // bw


def _url_size(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    try:
        with urllib.request.urlopen(req) as resp:
            cr = resp.headers.get("Content-Range")
            if cr is None:
                raise RuntimeError("Missing Content-Range header")
    except urllib.error.HTTPError as e:
        raise RuntimeError(
            f"Range request failed ({e.code}). URL may not support byte ranges."
        ) from e
    # format: "bytes 0-0/12345"
    try:
        return int(cr.split("/")[-1])
    except Exception as e:
        raise RuntimeError(f"Could not parse Content-Range: {cr!r}") from e


class _RandomAccessReader(Protocol):
    def read_at(self, offset: int, nbytes: int) -> bytes: ...
    def size(self) -> int: ...


class _FileReader:
    def __init__(self, fin):
        self._fin = fin
        self._size = int(os.fstat(fin.fileno()).st_size)

    def read_at(self, offset: int, nbytes: int) -> bytes:
        self._fin.seek(int(offset))
        data = self._fin.read(int(nbytes))
        if len(data) != int(nbytes):
            raise AssertionError("short read")
        return data

    def size(self) -> int:
        return self._size


class _UrlReader:
    def __init__(self, url: str):
        self._url = url
        self._size = None

    def read_at(self, offset: int, nbytes: int) -> bytes:
        start = int(offset)
        n = int(nbytes)
        req = urllib.request.Request(
            self._url, headers={"Range": f"bytes={start}-{start + n - 1}"}
        )
        try:
            with urllib.request.urlopen(req) as resp:
                data = resp.read()
        except urllib.error.HTTPError as e:
            raise RuntimeError(
                f"Range request failed ({e.code}). URL may not support byte ranges."
            ) from e
        if len(data) != n:
            raise RuntimeError(f"Short read: wanted {n} bytes @ {start}, got {len(data)}")
        return data

    def size(self) -> int:
        if self._size is None:
            self._size = _url_size(self._url)
        return int(self._size)


def _open_reader(src):
    if isinstance(src, (_FileReader, _UrlReader)):
        return src
    is_url = isinstance(src, str) and (
        src.startswith("http://") or src.startswith("https://")
    )
    if is_url:
        return _UrlReader(src)
    return src


def _parse_multi_index_reader(reader: _RandomAccessReader, seq_idx):
    header = reader.read_at(0, 8 * _HEADER_SIZE)
    _format = _parse_header_word(np.frombuffer(header[:8], dtype=np.uint64, count=1)[0])
    assert _format == 1, "Invalid FORMAT, multi-index format is type 1."
    bin_width = np.frombuffer(header[16:24], dtype=np.uint64, count=1)[0]
    num_seqs = np.frombuffer(header[24:32], dtype=np.uint64, count=1)[0]
    if not (0 <= int(seq_idx) < int(num_seqs)):
        raise AssertionError("invalid seq_idx")

    off = (_HEADER_SIZE * 8) + 8 * int(seq_idx)
    doc_offs = reader.read_at(off, 16)
    doc_start, doc_end = np.frombuffer(doc_offs, dtype=np.uint64, count=2)
    if not (int(doc_start) < int(doc_end)):
        raise AssertionError("corrupt index (doc offsets)")

    doc_raw = reader.read_at(int(doc_start), int(doc_end - doc_start))
    if len(doc_raw) < 8:
        raise AssertionError("corrupt index (doc too small)")

    bin_num = np.uint64(_U64.unpack_from(doc_raw, 0)[0])
    bn = int(bin_num)
    span_bytes = bn * 16
    boundaries_bytes = (bn + 1) * 8
    header_bytes = 8 + span_bytes + boundaries_bytes
    if len(doc_raw) < header_bytes:
        raise AssertionError("corrupt index (doc truncated)")

    if bn == 0:
        bin_spans = np.empty((0, 2), dtype=np.uint64)
    else:
        bin_spans = np.frombuffer(
            doc_raw, dtype=np.uint64, count=2 * bn, offset=8
        ).reshape(bn, 2, order="C")
    boundaries = np.frombuffer(
        doc_raw, dtype=np.uint64, count=bn + 1, offset=8 + span_bytes
    )
    ranges_raw = memoryview(doc_raw)[header_bytes:]
    if (len(ranges_raw) % 16) != 0:
        raise AssertionError("corrupt index (ranges not multiple of 2*u64)")
    ranges = (
        np.frombuffer(ranges_raw, dtype=np.uint64).reshape(-1, 2)
        if len(ranges_raw)
        else np.empty((0, 2), dtype=np.uint64)
    )
    if boundaries.size and int(boundaries[-1]) != len(ranges_raw):
        raise AssertionError("corrupt index (boundary/range size mismatch)")

    return BumblMultiIndex(
        seq_idx=int(seq_idx),
        bin_width=np.uint64(bin_width),
        num_seqs=np.uint64(num_seqs),
        bin_num=bin_num,
        max_bin=int(bin_num) - 1,
        bin_spans=bin_spans.copy(),
        boundaries=boundaries,
        ranges=ranges,
    )


def _parse_single_index_reader(reader: _RandomAccessReader, seq_idx):
    header = reader.read_at(0, 8 * 5)
    _format = _parse_header_word(np.frombuffer(header[:8], dtype=np.uint64, count=1)[0])
    assert _format == 0, "Invalid FORMAT, single-index format is type 0."
    bin_width = np.frombuffer(header[16:24], dtype=np.uint64, count=1)[0]
    num_seqs = np.frombuffer(header[24:32], dtype=np.uint64, count=1)[0]
    stored_seq_idx = np.frombuffer(header[32:40], dtype=np.uint64, count=1)[0]
    if seq_idx is not None and int(seq_idx) != int(stored_seq_idx):
        raise AssertionError("seq_idx mismatch for single-index")

    offsets_off = 8 * 5
    tail = reader.read_at(offsets_off, reader.size() - offsets_off)
    tail_len = len(tail)
    if tail_len % 8 != 0:
        raise AssertionError("corrupt index (tail not multiple of 8)")
    n_u64 = tail_len // 8
    if n_u64 < 1 or (n_u64 - 1) % 3 != 0:
        raise AssertionError("corrupt index (single-index tail length)")
    bin_num = (n_u64 - 1) // 3
    offsets = np.frombuffer(tail, dtype=np.uint64, count=bin_num + 1, offset=0).copy()
    if bin_num == 0:
        bin_spans = np.empty((0, 2), dtype=np.uint64)
    else:
        bin_spans = np.frombuffer(
            tail, dtype=np.uint64, count=2 * bin_num, offset=8 * (bin_num + 1)
        ).reshape(bin_num, 2, order="C").copy()

    return BumblSingleIndex(
        bin_width=np.uint64(bin_width),
        num_seqs=np.uint64(num_seqs),
        seq_idx=np.uint64(stored_seq_idx),
        offsets=offsets,
        bin_spans=bin_spans,
        max_bin=bin_num - 1,
    )


def parse_index(src, seq_idx=None):
    reader = _open_reader(src)
    if isinstance(reader, (_UrlReader, _FileReader)):
        r = reader
    else:
        fin = open(reader, "rb")
        try:
            r = _FileReader(fin)
        except Exception:
            fin.close()
            raise

    try:
        index_format = _parse_header_word(
            np.frombuffer(r.read_at(0, 8), dtype=np.uint64, count=1)[0]
        )
        if index_format == 0:
            return _parse_single_index_reader(r, seq_idx=seq_idx)
        if index_format == 1:
            if seq_idx is None:
                raise ValueError("seq_idx is required for multi-index")
            return _parse_multi_index_reader(r, seq_idx=seq_idx)
        raise AssertionError("unknown index format")
    finally:
        if not isinstance(reader, (_UrlReader, _FileReader)):
            fin.close()

# Version of BUMBL range helpers below (bump when layout, dtypes, or HTTP behavior changes).
_BUMBL_RANGE_HELPERS_VERSION = 2


def parse_bumbl_range(mumfile, mum_ranges):
    """
    Load MUM rows for half-open index intervals mum_ranges[i] = [start, end).

    mum_ranges: (N, 2) array; each row is [start, end) in MUM index order.
    Negative indices count from n_mums (same rules as single-slice indexing).
    Strands are read from the bumbl packed strand block (same layout as mumemto.utils.parse_bumbl_generator).
    """
    length_size = 4
    start_size = 8

    length_chunks = []
    start_chunks = []
    strand_chunks = []

    reader = _open_reader(mumfile)
    if isinstance(reader, (_UrlReader, _FileReader)):
        r = reader
    else:
        fin = open(reader, "rb")
        try:
            r = _FileReader(fin)
        except Exception:
            fin.close()
            raise

    try:
        header = r.read_at(0, 2 + 8 + 8)
        np.frombuffer(header[:2], dtype=np.uint16, count=1)  # flags (unused)
        n_seqs = int(np.frombuffer(header[2:10], dtype=np.uint64, count=1)[0])
        n_mums = int(np.frombuffer(header[10:18], dtype=np.uint64, count=1)[0])

        lengths_pos = 2 + 8 + 8
        offsets_pos = lengths_pos + (n_mums * length_size)
        strands_pos = offsets_pos + (n_mums * n_seqs * start_size)

        for mum_start, mum_end in mum_ranges:
            mum_start = int(mum_start)
            mum_end = int(mum_end)
            if mum_start < 0:
                mum_start = n_mums + mum_start
            if mum_end < 0:
                mum_end = n_mums + mum_end
            mum_start = max(0, min(mum_start, n_mums))
            mum_end = max(mum_start, min(mum_end, n_mums))
            n_sel = mum_end - mum_start
            if n_sel == 0:
                continue

            lengths_bytes = r.read_at(lengths_pos + mum_start * length_size, n_sel * length_size)
            length_chunks.append(np.frombuffer(lengths_bytes, dtype=np.uint32, count=n_sel).copy())

            starts_off = offsets_pos + mum_start * n_seqs * start_size
            starts_nbytes = n_sel * n_seqs * start_size
            starts_bytes = r.read_at(starts_off, starts_nbytes)
            start_chunks.append(
                np.frombuffer(starts_bytes, dtype=np.int64, count=n_sel * n_seqs)
                .reshape((n_sel, n_seqs))
                .copy()
            )

            bit0 = mum_start * n_seqs
            n_bits = n_sel * n_seqs
            byte0 = bit0 // 8
            byte1 = (bit0 + n_bits + 7) // 8
            packed = r.read_at(strands_pos + byte0, byte1 - byte0)
            bits = np.unpackbits(np.frombuffer(packed, dtype=np.uint8))
            off = bit0 % 8
            strand_chunks.append(bits[off : off + n_bits].reshape((n_sel, n_seqs)).copy())
    finally:
        if not isinstance(reader, (_UrlReader, _FileReader)):
            fin.close()

    if not length_chunks:
        lengths = np.array([], dtype=np.uint32)
        starts = np.empty((0, n_seqs), dtype=np.int64)
        strands = np.empty((0, n_seqs), dtype=bool)
    else:
        lengths = np.concatenate(length_chunks)
        starts = np.vstack(start_chunks)
        strands = np.vstack(strand_chunks)

    return utils.MUMdata.from_arrays(lengths, starts, strands)


def find_chr(starts, lengths):
    offsets = np.cumsum(lengths)
    contig_idx = np.searchsorted(offsets, starts, side="right")
    left_start = np.hstack((0, offsets[:-1]))
    rel_offsets = starts - left_start[contig_idx]
    return contig_idx, rel_offsets


def convert_local_to_global_coords(coords, names, lengths):
    """Convert contig-local `contig:start-end` to global offsets."""
    coords = coords.split(":")
    contig = coords[0]
    start, end = int(coords[1].split("-")[0]), int(coords[1].split("-")[1])
    assert contig in names, f"sequence {contig} not found in indicated FASTA file"
    contig_idx = names.index(contig)
    assert start <= end < lengths[contig_idx], f'Region {start}-{end} is invalid for contig {contig} with length {lengths[contig_idx]}'
    offset = sum(lengths[: contig_idx])
    return offset + start, offset + end


def convert_global_to_local_coords(start, end, names, lengths):
    contig, rel_offsets = find_chr((start, end), lengths)
    assert contig[0] == contig[1], (
        f"start and end coords are in different contigs: {names[contig[0]]} and {names[contig[1]]}"
    )
    return names[contig[0]], rel_offsets
