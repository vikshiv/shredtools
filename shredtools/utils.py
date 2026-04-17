"""Helpers for bumbl index files (``.bumbl.bi``).

Schemas: see bumbl_index.md for description
"""

from dataclasses import dataclass
import hashlib
import struct

import numpy as np
import mumemto.utils as utils 

_U64 = struct.Struct("<Q")
_U64x2 = struct.Struct("<QQ")

_HEADER_MAGIC = b"bumblbi"
_HEADER_SIZE = 4
_CHECKSUM_SAMPLE_SIZE = 1000


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
        bin_offsets = []
        flattened_ranges = []
        for pairs in range_maps[seq_idx]:
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
    offsets = np.searchsorted(mums.starts[:, seq_idx], edges, side="left")[:-1].astype(np.uint64)
    # write header word, checksum, bin width, num_seqs, seq_idx, then offset array
    checksum = bumbl_lengths_checksum(mums.lengths)
    return np.insert(offsets, 0, [_build_header_word(0), checksum, bin_width, mums.num_seqs, np.uint64(seq_idx)])

def write_bumbl_index(path, index):
    """Write bumbl index to file"""
    assert index.dtype == np.uint64
    with open(path, 'wb') as f:
        f.write(index.tobytes())


def get_bin_multi_file(index_path, seq_idx, bins):
    """
    Fast selective read from *.bumbl.bi. Must be multi-index.
    Returns np.ndarray shape (n, 2), dtype=np.uint64 (zero-copy view of read bytes).
    """
    bin_start, bin_end = _parse_bins_arg(bins)
    with open(index_path, "rb") as fin:
        # header
        fin.seek(0)
        _format = _parse_header_word(_U64.unpack(fin.read(8))[0])
        assert _format == 1, "Invalid FORMAT, multi-index format is type 1."
        fin.read(8)  # checksum (uint64)
        _bin_width = _U64.unpack(fin.read(8))[0]
        num_seqs   = _U64.unpack(fin.read(8))[0]
        if not (0 <= seq_idx < num_seqs):
            raise AssertionError("invalid seq_idx")

        # doc_offsets[seq_idx], doc_offsets[seq_idx+1]
        fin.seek((_HEADER_SIZE * 8) + 8 * seq_idx)
        doc_start, doc_end = _U64x2.unpack(fin.read(16))
        if not (doc_start < doc_end):
            raise AssertionError("corrupt index (doc offsets)")

        # bin_num
        fin.seek(doc_start)
        bin_num = _U64.unpack(fin.read(8))[0]
        ranges = []
        for bin_idx in range(bin_start, bin_end + 1):
            if not (0 <= bin_idx < bin_num):
                raise AssertionError("invalid bin")

            # off0, off1 (byte offsets into ranges area)
            fin.seek(doc_start + 8 * (1 + bin_idx))
            off0, off1 = _U64x2.unpack(fin.read(16))
            if off1 < off0:
                raise AssertionError("corrupt index (bin offsets)")

            ranges_base = doc_start + 8 * (bin_num + 2)  # 1 + (bin_num+1)
            byte_start = ranges_base + off0
            byte_end   = ranges_base + off1

            fin.seek(byte_start)
            raw = fin.read(byte_end - byte_start)
            if (len(raw) % 16) != 0:
                raise AssertionError("corrupt index (ranges not multiple of 2*u64)")

            ranges.append(np.frombuffer(raw, dtype=np.uint64).reshape(-1, 2))
        ranges = np.concatenate(ranges, axis=0) if ranges else np.empty((0, 2), dtype=np.uint64)
    return ranges

def get_bin_single_file(index_path, bins):
    """
    Fast selective read from *.bumbl.bi. Must be single-index.
    Returns np.ndarray shape (n, 2), dtype=np.uint64 (zero-copy view of read bytes).
    """
    bin_start, bin_end = _parse_bins_arg(bins)
    with open(index_path, "rb") as fin:
        # header
        fin.seek(0)
        _format = _parse_header_word(_U64.unpack(fin.read(8))[0])
        assert _format == 0, "Invalid FORMAT, single-index format is type 0."
        fin.seek(8 * 4, 1)  # checksum, bin_width, num_seqs, seq_idx

        # Remaining bytes are a uint64 offsets array. For bin i, the MUM row range is:
        # [offsets[i], offsets[i+1]).
        offsets = np.fromfile(fin, dtype=np.uint64)
        assert offsets.size >= 2, "index is empty."
        bin_num = offsets.size - 1
        assert 0 <= bin_start <= bin_end < bin_num, "invalid bin range."
        starts = offsets[bin_start : bin_end + 1]
        ends = offsets[bin_start + 1 : bin_end + 2]
        return np.stack((starts, ends), axis=1)
    
    

####### Functions to parse a BumblBi index #######



### Convention: only store the relevant document index in memory for a multi-index

@dataclass(frozen=True, slots=True)
class BumblMultiIndex:
    """In-memory view of one document (one seq_idx) inside a multi-index."""

    seq_idx: int
    bin_width: np.uint64
    num_seqs: np.uint64
    bin_num: np.uint64
    boundaries: np.ndarray  # uint64, shape (bin_num+1,), byte offsets into ranges section
    ranges: np.ndarray  # uint64, shape (k,2), (mum_start, mum_end) pairs


@dataclass(frozen=True, slots=True)
class BumblSingleIndex:
    """In-memory view of a single-index."""

    bin_width: np.uint64
    num_seqs: np.uint64
    seq_idx: np.uint64
    offsets: np.ndarray  # uint64, shape (bin_num+1,)


def parse_multi_index(index_path, seq_idx):
    """
    Parse a multi-index (*.bumbl.bi, FORMAT=1) into an in-memory object
    containing only the requested document's slice.
    """
    with open(index_path, "rb") as fin:
        fin.seek(0)
        _format = _parse_header_word(_U64.unpack(fin.read(8))[0])
        assert _format == 1, "Invalid FORMAT, multi-index format is type 1."
        fin.read(8)  # checksum
        bin_width = np.uint64(_U64.unpack(fin.read(8))[0])
        num_seqs = np.uint64(_U64.unpack(fin.read(8))[0])
        if not (0 <= int(seq_idx) < int(num_seqs)):
            raise AssertionError("invalid seq_idx")

        # doc_offsets[seq_idx], doc_offsets[seq_idx+1]
        fin.seek((_HEADER_SIZE * 8) + 8 * int(seq_idx))
        doc_start, doc_end = _U64x2.unpack(fin.read(16))
        if not (doc_start < doc_end):
            raise AssertionError("corrupt index (doc offsets)")

        fin.seek(doc_start)
        doc_raw = fin.read(doc_end - doc_start)

    if len(doc_raw) < 8:
        raise AssertionError("corrupt index (doc too small)")

    bin_num = np.uint64(_U64.unpack_from(doc_raw, 0)[0])
    boundaries_bytes = (int(bin_num) + 1) * 8
    header_bytes = 8 + boundaries_bytes
    if len(doc_raw) < header_bytes:
        raise AssertionError("corrupt index (doc truncated)")

    boundaries = np.frombuffer(doc_raw, dtype=np.uint64, count=int(bin_num) + 1, offset=8)
    ranges_raw = memoryview(doc_raw)[header_bytes:]
    if (len(ranges_raw) % 16) != 0:
        raise AssertionError("corrupt index (ranges not multiple of 2*u64)")
    ranges = (
        np.frombuffer(ranges_raw, dtype=np.uint64).reshape(-1, 2)
        if len(ranges_raw)
        else np.empty((0, 2), dtype=np.uint64)
    )

    # Optional integrity check: last boundary should match ranges section length.
    if boundaries.size and int(boundaries[-1]) != len(ranges_raw):
        raise AssertionError("corrupt index (boundary/range size mismatch)")

    return BumblMultiIndex(
        seq_idx=int(seq_idx),
        bin_width=bin_width,
        num_seqs=num_seqs,
        bin_num=bin_num,
        boundaries=boundaries,
        ranges=ranges,
    )


def parse_single_index(index_path, seq_idx):
    """
    Parse a single-index (*.bumbl.bi, FORMAT=0) into an in-memory object.

    If `seq_idx` is provided, it is validated against the value stored in the file.
    """
    with open(index_path, "rb") as fin:
        fin.seek(0)
        _format = _parse_header_word(_U64.unpack(fin.read(8))[0])
        assert _format == 0, "Invalid FORMAT, single-index format is type 0."
        fin.read(8)  # checksum
        bin_width = np.uint64(_U64.unpack(fin.read(8))[0])
        num_seqs = np.uint64(_U64.unpack(fin.read(8))[0])
        stored_seq_idx = np.uint64(_U64.unpack(fin.read(8))[0])
        offsets = np.fromfile(fin, dtype=np.uint64)

    if offsets.size < 2:
        raise AssertionError("index is empty.")
    if seq_idx is not None and int(seq_idx) != int(stored_seq_idx):
        raise AssertionError("seq_idx mismatch for single-index")

    return BumblSingleIndex(
        bin_width=bin_width,
        num_seqs=num_seqs,
        seq_idx=stored_seq_idx,
        offsets=offsets,
    )

def get_bin_multi(index, seq_idx, bins):
    """
    Query ranges from an in-memory multi-index doc slice.

    Returns np.ndarray shape (n, 2), dtype=np.uint64 of (mum_start, mum_end) pairs.
    """
    if not isinstance(index, BumblMultiIndex):
        raise TypeError("index must be a BumblMultiIndex (use parse_multi_index)")
    if int(seq_idx) != int(index.seq_idx):
        raise AssertionError("seq_idx mismatch for this parsed multi-index doc")

    bin_start, bin_end = _parse_bins_arg(bins)
    bin_num = int(index.bin_num)
    if not (0 <= bin_start <= bin_end < bin_num):
        raise AssertionError("invalid bin")

    out = []
    for bin_idx in range(bin_start, bin_end + 1):
        off0 = int(index.boundaries[bin_idx])
        off1 = int(index.boundaries[bin_idx + 1])
        if off1 < off0:
            raise AssertionError("corrupt index (bin offsets)")
        if (off0 % 16) != 0 or (off1 % 16) != 0:
            raise AssertionError("corrupt index (bin offsets not multiple of 2*u64)")
        i0 = off0 // 16
        i1 = off1 // 16
        out.append(index.ranges[i0:i1])

    return np.concatenate(out, axis=0) if out else np.empty((0, 2), dtype=np.uint64)


def get_bin_single(index, bins):
    """
    Query row ranges from an in-memory single-index.

    Returns np.ndarray shape (n, 2), dtype=np.uint64 of (mum_start, mum_end) row ranges.
    """
    if not isinstance(index, BumblSingleIndex):
        raise TypeError("index must be a BumblSingleIndex (use parse_single_index)")

    bin_start, bin_end = _parse_bins_arg(bins)
    offsets = index.offsets
    bin_num = int(offsets.size - 1)
    if not (0 <= bin_start <= bin_end < bin_num):
        raise AssertionError("invalid bin range.")
    starts = offsets[bin_start : bin_end + 1]
    ends = offsets[bin_start + 1 : bin_end + 2]
    return np.stack((starts, ends), axis=1)

    