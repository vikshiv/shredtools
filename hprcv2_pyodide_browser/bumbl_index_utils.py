# Vendored from shredtools/utils.py (index + range I/O for Pyodide / static hosting).
# HPRCv2 uses a multi-index only (FORMAT=1); single-index parsing is omitted.
#
# To avoid a slow Pyodide micropip install, this file also includes a tiny subset
# of `mumemto.utils` needed by this browser:
# - multilengths parsing (`get_sequence_lengths`, `get_contig_names`)
# - a minimal `MUMdata` container (row indexing)

from dataclasses import dataclass
from array import array
from bisect import bisect_right
import os
import struct
import sys
from typing import Protocol
import urllib.error
import urllib.request


def get_sequence_lengths(lengths_file, multilengths=False):
    # Copied from https://github.com/vikshiv/mumemto/blob/main/mumemto/utils.py
    def get_lengths(lengths_file):
        return [int(l.split()[1]) for l in open(lengths_file, "r").read().splitlines()]

    def get_multilengths(lengths_file):
        offset = []
        cur_offset = []
        for l in open(lengths_file, "r").readlines():
            l = l.strip().split()
            if l[1] == "*":
                if cur_offset:
                    offset.append(cur_offset)
                cur_offset = []
                continue
            cur_offset.append(int(l[2]))
        offset.append(cur_offset)
        return offset

    simple = True
    try:
        with open(lengths_file, "r") as f:
            first_line = f.readline().strip().split()
            if len(first_line) > 1 and first_line[1] == "*":
                simple = False
    except FileNotFoundError:
        raise FileNotFoundError(f"File {lengths_file} not found.")
    if simple and multilengths:
        raise ValueError("Multi-FASTA lengths not available in ", lengths_file)
    if not simple:
        offsets = get_multilengths(lengths_file)
        return offsets if multilengths else [sum(o) for o in offsets]
    else:
        return get_lengths(lengths_file)


def get_contig_names(lengths_file):
    """
    Copied from https://github.com/vikshiv/mumemto/blob/main/mumemto/utils.py

    Get contig names from a multilengths file.
    Returns a list of lists where each inner list contains the contig names for one sequence.

    Args:
        lengths_file: Path to the lengths file in multilengths format

    Returns:
        List of lists where names[i] is the list of contig names for sequence i
    """
    names = []
    cur_name = []
    first_line = True
    for l in open(lengths_file, "r").readlines():
        l = l.strip().split()
        if first_line and l[1] != "*":
            raise ValueError("Lengths file must be formatted as multilengths.")
        first_line = False
        if l[1] == "*":
            if cur_name:
                names.append(cur_name)
            cur_name = []
            continue
        cur_name.append(l[1])
    names.append(cur_name)
    return names


@dataclass(frozen=True, slots=True)
class _MUMRow:
    length: int
    starts: list[int]
    strands: list[bool]


@dataclass(slots=True)
class MUMdata:
    """
    Minimal MUM container for the browser (pure Python, no numpy).

    Storage is row-major:
    - lengths[i] is uint32 length
    - starts[(i*n_seqs)+j] is int64 start for seq j
    - strands[(i*n_seqs)+j] is 0/1
    """

    n_seqs: int
    lengths: array  # 'I'
    starts: array  # 'q'
    strands: bytearray  # 0/1 per entry, length = n_mums*n_seqs

    @property
    def num_mums(self) -> int:
        return len(self.lengths)

    @property
    def num_seqs(self) -> int:
        return int(self.n_seqs)

    @classmethod
    def empty(cls, n_seqs: int):
        return cls(int(n_seqs), array("I"), array("q"), bytearray())

    def start(self, i: int, j: int) -> int:
        return int(self.starts[i * self.n_seqs + j])

    def strand(self, i: int, j: int) -> bool:
        return bool(self.strands[i * self.n_seqs + j])

    def starts_col(self, j: int) -> list[int]:
        n = self.num_mums
        ns = self.n_seqs
        base = int(j)
        return [int(self.starts[i * ns + base]) for i in range(n)]

    def __getitem__(self, idx):
        if isinstance(idx, int):
            i = int(idx)
            ns = self.n_seqs
            s0 = i * ns
            starts = [int(self.starts[s0 + j]) for j in range(ns)]
            strands = [bool(self.strands[s0 + j]) for j in range(ns)]
            return _MUMRow(length=int(self.lengths[i]), starts=starts, strands=strands)
        raise TypeError("Only integer row indexing is supported in this minimal MUMdata.")

_U64 = struct.Struct("<Q")
_HEADER_MAGIC = b"bumblbi"
_HEADER_SIZE = 4


def _array_from_bytes(typecode: str, b: bytes) -> array:
    a = array(typecode)
    a.frombytes(b)
    # Ensure little-endian interpretation.
    if sys.byteorder != "little":
        a.byteswap()
    return a


def _parse_bins_arg(bins):
    if type(bins) == int:
        return bins, bins
    if isinstance(bins, (tuple, list)) and len(bins) == 2:
        return int(bins[0]), int(bins[1])
    raise ValueError("Invalid type for bins. Must be int or (start, end).")


def _parse_header_word(header_word):
    raw = int(header_word).to_bytes(8, "little")
    if raw[:7] != _HEADER_MAGIC:
        raise AssertionError("invalid bumbl index header")
    return raw[7]


@dataclass(frozen=True, slots=True)
class BumblMultiIndex:
    """In-memory view of one document (one seq_idx) inside a multi-index."""

    src: str
    seq_idx: int
    bin_width: int
    num_seqs: int
    bin_num: int
    boundaries: array  # 'Q' byte offsets
    range_set_base: int  # absolute byte offset of the ranges section start

    async def get_bins(self, bins):
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
        if off1 == off0:
            return []

        # Fetch only the needed range bytes for these bins.
        reader = _open_reader(self.src)
        if not isinstance(reader, _UrlReader):
            # Keep the browser app URL-only; local file support can be added if needed.
            raise RuntimeError("BumblMultiIndex.get_bins requires a URL src in this app")
        raw = await reader.read_at(self.range_set_base + off0, off1 - off0)
        u64s = _array_from_bytes("Q", raw)
        if (len(u64s) % 2) != 0:
            raise AssertionError("corrupt index (range pairs)")
        out = []
        for i in range(0, len(u64s), 2):
            out.append((int(u64s[i]), int(u64s[i + 1])))
        return out

    def closest_nonzero_bin_left(self, bin_idx: int):
        n = int(self.bin_num)
        b = int(bin_idx)
        if not (0 <= b < n):
            raise AssertionError("invalid bin")
        boundaries = self.boundaries
        for i in range(b, -1, -1):
            if boundaries[i] != boundaries[i + 1]:
                return i
        return None

    def closest_nonzero_bin_right(self, bin_idx: int):
        n = int(self.bin_num)
        b = int(bin_idx)
        if not (0 <= b < n):
            raise AssertionError("invalid bin")
        boundaries = self.boundaries
        for i in range(b, n):
            if boundaries[i] != boundaries[i + 1]:
                return i
        return None

    def coord_to_bin(self, coord) -> int:
        bw = int(self.bin_width)
        max_bin = int(self.bin_num) - 1
        b = int(coord) // bw
        return max(0, min(b, max_bin))


def _url_size(url):
    raise RuntimeError("sync _url_size is not supported in Pyodide; use _url_size_async")


class _RandomAccessReader(Protocol):
    def read_at(self, offset: int, nbytes: int) -> bytes: ...
    def size(self) -> int: ...


class _AsyncRandomAccessReader(Protocol):
    async def read_at(self, offset: int, nbytes: int) -> bytes: ...
    async def size(self) -> int: ...


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

    async def read_at(self, offset: int, nbytes: int) -> bytes:
        # In Pyodide, urllib does not support https (no SSL); use browser fetch via pyfetch.
        # pyodide.http exists at runtime inside Pyodide (not in CPython).
        from pyodide.http import pyfetch  # type: ignore

        start = int(offset)
        n = int(nbytes)
        end = start + n - 1
        resp = await pyfetch(self._url, headers={"Range": f"bytes={start}-{end}"})
        if resp.status not in (200, 206):
            raise RuntimeError(f"Range request failed ({resp.status}). URL may not support byte ranges.")
        data = await resp.bytes()
        if len(data) != n:
            raise RuntimeError(f"Short read: wanted {n} bytes @ {start}, got {len(data)}")
        return data

    async def size(self) -> int:
        if self._size is None:
            self._size = await _url_size_async(self._url)
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


async def _url_size_async(url):
    from pyodide.http import pyfetch  # type: ignore

    resp = await pyfetch(url, headers={"Range": "bytes=0-0"})
    if resp.status not in (200, 206):
        raise RuntimeError(f"Range request failed ({resp.status}). URL may not support byte ranges.")
    cr = resp.headers.get("Content-Range")
    if cr is None:
        raise RuntimeError("Missing Content-Range header")
    try:
        return int(cr.split("/")[-1])
    except Exception as e:
        raise RuntimeError(f"Could not parse Content-Range: {cr!r}") from e


async def _parse_multi_index_reader(reader: _AsyncRandomAccessReader, seq_idx):
    header = await reader.read_at(0, 8 * _HEADER_SIZE)
    _format = _parse_header_word(_U64.unpack_from(header, 0)[0])
    assert _format == 1, "Invalid FORMAT, multi-index format is type 1."
    bin_width = int(_U64.unpack_from(header, 16)[0])
    num_seqs = int(_U64.unpack_from(header, 24)[0])
    if not (0 <= int(seq_idx) < int(num_seqs)):
        raise AssertionError("invalid seq_idx")

    off = (_HEADER_SIZE * 8) + 8 * int(seq_idx)
    doc_offs = await reader.read_at(off, 16)
    doc_start = int(_U64.unpack_from(doc_offs, 0)[0])
    doc_end = int(_U64.unpack_from(doc_offs, 8)[0])
    if not (doc_start < doc_end):
        raise AssertionError("corrupt index (doc offsets)")

    # Read only doc header + boundaries; do not fetch the full ranges table.
    doc_head = await reader.read_at(doc_start, 8)
    if len(doc_head) != 8:
        raise AssertionError("corrupt index (doc too small)")
    bin_num = int(_U64.unpack_from(doc_head, 0)[0])
    boundaries_bytes = (bin_num + 1) * 8
    boundaries_raw = await reader.read_at(doc_start + 8, boundaries_bytes)
    if len(boundaries_raw) != boundaries_bytes:
        raise AssertionError("corrupt index (doc truncated)")
    boundaries = _array_from_bytes("Q", boundaries_raw)
    range_set_base = doc_start + 8 + boundaries_bytes

    return BumblMultiIndex(
        src=getattr(reader, "_url", ""),  # best-effort; replaced in parse_index wrapper
        seq_idx=int(seq_idx),
        bin_width=int(bin_width),
        num_seqs=int(num_seqs),
        bin_num=int(bin_num),
        boundaries=boundaries,
        range_set_base=int(range_set_base),
    )


async def parse_index(src, seq_idx=None):
    reader = _open_reader(src)
    if isinstance(reader, _UrlReader):
        r = reader
        fin = None
    elif isinstance(reader, _FileReader):
        r = reader
        fin = None
    else:
        fin = open(reader, "rb")
        r = _FileReader(fin)

    try:
        index_format = _parse_header_word(_U64.unpack(await r.read_at(0, 8))[0])
        if index_format != 1:
            raise ValueError(
                "Expected multi-index (FORMAT=1) for this app; "
                f"index has FORMAT={index_format}."
            )
        if seq_idx is None:
            raise ValueError("seq_idx is required for multi-index")
        idx = await _parse_multi_index_reader(r, seq_idx=seq_idx)
        # Ensure the returned index knows its URL for on-demand bin range fetches.
        return BumblMultiIndex(
            src=str(src),
            seq_idx=idx.seq_idx,
            bin_width=idx.bin_width,
            num_seqs=idx.num_seqs,
            bin_num=idx.bin_num,
            boundaries=idx.boundaries,
            range_set_base=idx.range_set_base,
        )
    finally:
        if fin is not None:
            fin.close()


async def parse_bumbl_range(mumfile, mum_ranges):
    length_size = 4
    start_size = 8

    lengths_out = array("I")
    starts_out = array("q")
    strands_out = bytearray()

    reader = _open_reader(mumfile)
    if isinstance(reader, _UrlReader):
        r = reader
        fin = None
    elif isinstance(reader, _FileReader):
        r = reader
        fin = None
    else:
        fin = open(reader, "rb")
        r = _FileReader(fin)

    try:
        header = await r.read_at(0, 2 + 8 + 8)
        # flags (uint16) unused
        n_seqs = int(_U64.unpack_from(header, 2)[0])
        n_mums = int(_U64.unpack_from(header, 10)[0])

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

            lengths_bytes = await r.read_at(
                lengths_pos + mum_start * length_size, n_sel * length_size
            )
            lengths_chunk = _array_from_bytes("I", lengths_bytes)
            lengths_out.extend(lengths_chunk)

            starts_off = offsets_pos + mum_start * n_seqs * start_size
            starts_nbytes = n_sel * n_seqs * start_size
            starts_bytes = await r.read_at(starts_off, starts_nbytes)
            starts_chunk = _array_from_bytes("q", starts_bytes)
            starts_out.extend(starts_chunk)

            bit0 = mum_start * n_seqs
            n_bits = n_sel * n_seqs
            byte0 = bit0 // 8
            byte1 = (bit0 + n_bits + 7) // 8
            packed = await r.read_at(strands_pos + byte0, byte1 - byte0)
# numpy's unpackbits default is MSB-first; replicate that mapping.
            off = bit0 % 8
            base_bit = off
            for k in range(n_bits):
                bit_index = base_bit + k
                b = packed[bit_index // 8]
                strands_out.append((b >> (7 - (bit_index % 8))) & 1)
    finally:
        if fin is not None:
            fin.close()

    return MUMdata(int(n_seqs), lengths_out, starts_out, strands_out)


def find_chr(starts, lengths):
    # Pure Python replacement for the previous numpy-based implementation.
    offsets = []
    s = 0
    for L in lengths:
        s += int(L)
        offsets.append(s)
    contig_idx = [bisect_right(offsets, int(x)) for x in starts]
    left_start = [0] + offsets[:-1]
    rel_offsets = [int(starts[i]) - left_start[contig_idx[i]] for i in range(len(starts))]
    return contig_idx, rel_offsets


def convert_local_to_global_coords(coords, names, lengths):
    """Convert contig-local `contig:start-end` to global offsets."""
    coords = coords.split(":")
    contig = coords[0]
    start, end = int(coords[1].split("-")[0]), int(coords[1].split("-")[1])
    if contig not in names:
        raise ValueError(f"sequence {contig} not found in indicated FASTA file")
    offset = sum(lengths[: names.index(contig)])
    return offset + start, offset + end


def convert_global_to_local_coords(start, end, names, lengths):
    contig, rel_offsets = find_chr((start, end), lengths)
    if contig[0] != contig[1]:
        # Match shredtools/utils.py exception type/message so callers can consistently
        # detect and handle split-contig projections.
        raise AssertionError(
            f"start and end coords are in different contigs: {names[contig[0]]} and {names[contig[1]]}"
        )
    return names[contig[0]], rel_offsets


async def get_mum_ranges_flanks(index, coords):
    idx = index
    s, e = coords
    bin_start = idx.coord_to_bin(s)
    bin_end = idx.coord_to_bin(e)

    left_bin = idx.closest_nonzero_bin_left(bin_start)
    right_bin = idx.closest_nonzero_bin_right(bin_end)
    if left_bin is None or right_bin is None:
        return [], (None, None), (bin_start, bin_end)

    left_bin = int(left_bin)
    right_bin = int(right_bin)

    left_ranges = await idx.get_bins(left_bin)
    if right_bin == left_bin:
        ranges = left_ranges
    else:
        right_ranges = await idx.get_bins(right_bin)
        ranges = left_ranges + right_ranges

    return ranges, (left_bin, right_bin), (bin_start, bin_end)
