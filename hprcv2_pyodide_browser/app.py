# HPRCv2 region browser — same data path as mod_scripts/index_extract.py (no CLI, no files out).
import json
from bisect import bisect_left
import re
from bisect import bisect_right

import bumbl_index_utils as sutils

BUMBL_BI_URL = "https://genome-idx.s3.amazonaws.com/mumemto/hprcv2_enhanced_merged.bumbl.bi"
BUMBL_URL = "https://genome-idx.s3.amazonaws.com/mumemto/hprcv2_enhanced_merged.bumbl"

# Caches (modest memory): lengths metadata and indices reused across runs.
_LENGTHS_META = {}  # lengths_path -> (seq_lengths_multi, contig_names, num_seqs)
_INDEX_BY_SEQ = {}  # seq_idx -> parsed index for BUMBL_BI_URL


def _get_lengths_meta(lengths_path: str):
    meta = _LENGTHS_META.get(lengths_path)
    if meta is not None:
        return meta
    seq_lengths_multi = sutils.get_sequence_lengths(lengths_path, multilengths=True)
    contig_names = sutils.get_contig_names(lengths_path)
    meta = (seq_lengths_multi, contig_names, len(seq_lengths_multi))
    _LENGTHS_META[lengths_path] = meta
    return meta


async def _get_index(seq_idx: int):
    idx = _INDEX_BY_SEQ.get(int(seq_idx))
    if idx is not None:
        return idx, True
    idx = await sutils.parse_index(BUMBL_BI_URL, seq_idx=int(seq_idx))
    _INDEX_BY_SEQ[int(seq_idx)] = idx
    return idx, False


async def warm_index(seq_idx: int) -> bool:
    """Preload the index for seq_idx; used on genome dropdown change."""
    idx = await sutils.parse_index(BUMBL_BI_URL, seq_idx=int(seq_idx))
    _INDEX_BY_SEQ[int(seq_idx)] = idx
    return True


def find_target_region(coll_mums, coords, seq_idx, sequences, right_key=None):
    """From mod_scripts/index_extract.py (no stderr prints in browser)."""
    n = int(coll_mums.num_mums)
    if n <= 0:
        raise ValueError("No MUMs loaded for this query (empty range slice).")
    starts_col = coll_mums.starts_col(seq_idx)
    left_mum_idx = bisect_left(starts_col, coords[0]) - 1
    if right_key is None:
        right_key = [starts_col[i] + int(coll_mums.lengths[i]) for i in range(n)]
    right_mum_idx = bisect_left(right_key, coords[1])
    # Clamp to valid row indices; bisect can yield -1 or n depending on coords.
    left_mum_idx = max(0, min(int(left_mum_idx), n - 1))
    right_mum_idx = max(0, min(int(right_mum_idx), n - 1))
    mum_bounds = (left_mum_idx, right_mum_idx)
    left_mum, right_mum = coll_mums[mum_bounds[0]], coll_mums[mum_bounds[1]]
    left_offset, right_offset = 0, 0
    if coords[0] < left_mum.starts[seq_idx] + left_mum.length:
        left_offset = coords[0] - left_mum.starts[seq_idx]
    if coords[1] > right_mum.starts[seq_idx]:
        right_offset = coords[1] - right_mum.starts[seq_idx]
    other_coords = [
        (
            coll_mums.start(mum_bounds[0], i) + left_offset,
            coll_mums.start(mum_bounds[1], i) + right_offset,
        )
        for i in sequences
    ]
    return mum_bounds, other_coords


def compute_margins(coll_mums, coords, seq_idx, right_key=None):
    """Match the left/right margin prints from mod_scripts/index_extract.py."""
    n = int(coll_mums.num_mums)
    if n <= 0:
        raise ValueError("No MUMs loaded for this query (empty range slice).")
    starts_col = coll_mums.starts_col(seq_idx)
    left_mum_idx = bisect_left(starts_col, coords[0]) - 1
    if right_key is None:
        right_key = [starts_col[i] + int(coll_mums.lengths[i]) for i in range(n)]
    right_mum_idx = bisect_left(right_key, coords[1])
    left_mum_idx = max(0, min(int(left_mum_idx), n - 1))
    right_mum_idx = max(0, min(int(right_mum_idx), n - 1))
    left_mum, right_mum = coll_mums[left_mum_idx], coll_mums[right_mum_idx]

    left_bound = left_mum.starts[seq_idx]
    right_bound = right_mum.starts[seq_idx]

    left_offset, right_offset = 0, 0
    if coords[0] < left_mum.starts[seq_idx] + left_mum.length:
        left_offset = coords[0] - left_mum.starts[seq_idx]
    if coords[1] > right_mum.starts[seq_idx]:
        right_offset = coords[1] - right_mum.starts[seq_idx]

    left_margin = coords[0] - left_bound - left_offset
    right_margin = right_bound + right_offset - coords[1]
    return int(left_margin), int(right_margin)


def _bracket_ok(mums, coords, seq_idx, right_key):
    """
    Return (left_ok, right_ok) indicating whether this MUM slice contains
    a MUM row before coords[0] and a MUM row whose end is after coords[1].
    """
    n = int(mums.num_mums)
    if n <= 0:
        return False, False
    starts_col = mums.starts_col(seq_idx)
    li = bisect_left(starts_col, coords[0]) - 1
    ri = bisect_left(right_key, coords[1])
    return (li >= 0), (ri < n)


async def _get_mums_expanding(idx, coords, seq_idx, max_steps: int = 8):
    """
    Fetch MUMs for bins around coords, expanding outward until bracketing holds.

    Starts with the same flanking-bin logic (<=2 bins), then widens left/right
    to additional non-empty bins if needed.
    """
    s, e = coords
    bin_start = idx.coord_to_bin(s)
    bin_end = idx.coord_to_bin(e)

    left_bin = idx.closest_nonzero_bin_left(bin_start)
    right_bin = idx.closest_nonzero_bin_right(bin_end)
    if left_bin is None or right_bin is None:
        return None, [], (bin_start, bin_end)
    left_bin = int(left_bin)
    right_bin = int(right_bin)

    steps = 0
    while True:
        ranges = await idx.get_bins(left_bin)
        if right_bin != left_bin:
            ranges = ranges + (await idx.get_bins(right_bin))
        if not ranges:
            return None, [], (bin_start, bin_end)

        mums = await sutils.parse_bumbl_range(BUMBL_URL, ranges)
        starts_col = mums.starts_col(seq_idx)
        right_key = [starts_col[i] + int(mums.lengths[i]) for i in range(int(mums.num_mums))]
        left_ok, right_ok = _bracket_ok(mums, coords, seq_idx, right_key)
        if left_ok and right_ok:
            return (mums, right_key), ranges, (bin_start, bin_end)

        if steps >= int(max_steps):
            # Give up with a clear message rather than crashing.
            raise ValueError(
                f"Could not bracket region from flanking bins after {max_steps} expansions "
                f"(bin_start={bin_start}, bin_end={bin_end}, left_bin={left_bin}, right_bin={right_bin})."
            )

        # Expand outward only on the side(s) that are missing.
        if not left_ok:
            lb2 = idx.closest_nonzero_bin_left(left_bin - 1)
            if lb2 is None:
                left_ok = True  # cannot expand further; stop trying left
            else:
                left_bin = int(lb2)
        if not right_ok:
            rb2 = idx.closest_nonzero_bin_right(right_bin + 1)
            if rb2 is None:
                right_ok = True
            else:
                right_bin = int(rb2)
        steps += 1


def format_bed(
    contig_names, seq_lengths_multi, other_coords, sequences, path_placeholder="."
):
    """Same rows as index_extract extract_bed; last column is a placeholder in the browser."""
    lines = []
    for i, seq in enumerate(sequences):
        name, rel_offsets = sutils.convert_global_to_local_coords(
            other_coords[i][0],
            other_coords[i][1],
            contig_names[int(seq)],
            seq_lengths_multi[int(seq)],
        )
        lines.append(
            f"{name}\t{int(rel_offsets[0])}\t{int(rel_offsets[1])}\t{path_placeholder}\n"
        )
    return "".join(lines)


def describe_ui(lengths_path: str) -> str:
    """JSON for genome / contig dropdowns. Call after writing lengths to MEMFS."""
    seq_lengths_multi, contig_names, n = _get_lengths_meta(lengths_path)
    genomes = []
    for i in range(n):
        cnames = contig_names[i]
        label = cnames[0] if cnames else f"seq_{i}"
        genomes.append(
            {
                "seq_idx": i,
                "label": label,
                "contigs": cnames,
            }
        )
    return json.dumps({"n_seqs": n, "genomes": genomes})


async def run(
    lengths_path: str,
    seq_idx: int,
    range_str: str,
    sequences=None,
) -> str:
    """
    Returns BED text or raises with a clear error message.
    `range_str` is chr:start-end (same as -r in index_extract).
    """
    seq_lengths_multi, contig_names, num_seqs = _get_lengths_meta(lengths_path)
    if not (0 <= seq_idx < num_seqs):
        raise ValueError(f"seq_idx {seq_idx} invalid (N = {num_seqs})")
    if sequences is not None and any(s >= num_seqs for s in sequences):
        raise ValueError(f"Invalid sequence index (N = {num_seqs})")
    if sequences is None:
        sequences = list(range(num_seqs))

    coords = sutils.convert_local_to_global_coords(
        range_str, contig_names[seq_idx], seq_lengths_multi[seq_idx]
    )
    idx, _cached = await _get_index(seq_idx)
    got, ranges, requested_bins = await _get_mums_expanding(idx, coords, seq_idx)
    if got is None:
        raise ValueError(
            f"No MUM ranges for bins {requested_bins} (index could not find flanking bins)"
        )
    mums, right_key = got
    _mum_bounds, other_coords = find_target_region(mums, coords, seq_idx, sequences, right_key=right_key)
    return format_bed(
        contig_names, seq_lengths_multi, other_coords, sequences, path_placeholder="."
    )


async def run_with_bounds(
    lengths_path: str,
    seq_idx: int,
    range_str: str,
    sequences=None,
) -> str:
    """
    Like `run`, but returns JSON: { bed: str, bounds: { contig, start, end } }.
    Bounds are the extracted interval for the selected genome (seq_idx), in contig-local coords.
    """
    seq_lengths_multi, contig_names, num_seqs = _get_lengths_meta(lengths_path)
    if not (0 <= seq_idx < num_seqs):
        raise ValueError(f"seq_idx {seq_idx} invalid (N = {num_seqs})")
    # For the browser app we always compute all sequences; UI-side filtering is handled in JS.
    sequences = list(range(num_seqs))

    coords = sutils.convert_local_to_global_coords(
        range_str, contig_names[seq_idx], seq_lengths_multi[seq_idx]
    )
    idx, index_cached = await _get_index(seq_idx)
    got, ranges, requested_bins = await _get_mums_expanding(idx, coords, seq_idx)
    if got is None:
        raise ValueError(
            f"No MUM ranges for bins {requested_bins} (index could not find flanking bins)"
        )
    # Each element of `ranges` is a half-open [mum_start, mum_end) slice.
    mum_chunks = int(len(ranges))
    mums_sliced = int(sum(int(b) - int(a) for a, b in ranges))

    mums, right_key = got
    _mum_bounds, other_coords = find_target_region(mums, coords, seq_idx, sequences, right_key=right_key)
    rows = []
    unavailable = []
    _span_re = re.compile(r"start and end coords are in different contigs:\s+(.+?)\s+and\s+(.+)$")
    for i, seq in enumerate(sequences):
        try:
            name, rel_offsets = sutils.convert_global_to_local_coords(
                other_coords[i][0],
                other_coords[i][1],
                contig_names[int(seq)],
                seq_lengths_multi[int(seq)],
            )
        except AssertionError as e:
            msg = e.args[0] if e.args else str(e)
            m = _span_re.search(str(msg))
            c1, c2 = (m.group(1), m.group(2)) if m else ("", "")
            label = contig_names[int(seq)][0] if contig_names[int(seq)] else f"seq_{int(seq)}"
            unavailable.append(
                {
                    "seq_idx": int(seq),
                    "label": label,
                    "contig_a": str(c1),
                    "contig_b": str(c2),
                    "reason": str(msg),
                }
            )
            continue
        rows.append(
            {
                "seq_idx": int(seq),
                "contig": name,
                "start": int(rel_offsets[0]),
                "end": int(rel_offsets[1]),
            }
        )

    # other_coords is aligned with `sequences`; when sequences is default, index==seq_idx.
    try:
        self_i = sequences.index(seq_idx)
    except ValueError:
        self_i = 0
    b0, b1 = other_coords[self_i]
    try:
        contig, rel = sutils.convert_global_to_local_coords(
            b0, b1, contig_names[seq_idx], seq_lengths_multi[seq_idx]
        )
    except AssertionError as e:
        # This happens when the *requested region itself* spans contigs on the selected genome
        # (e.g. end runs off the contig). Return a clean error instead of a Pyodide traceback.
        msg = e.args[0] if e.args else str(e)
        raise ValueError(f"Selected-genome region spans multiple contigs: {msg}") from None
    left_margin, right_margin = compute_margins(mums, coords, seq_idx, right_key=right_key)

    return json.dumps(
        {
            "rows": rows,
            "bounds": {"contig": contig, "start": int(rel[0]), "end": int(rel[1])},
            "margins": {"left": left_margin, "right": right_margin},
            "unavailable": unavailable,
            "mum_slices": {"chunks": mum_chunks, "mums": mums_sliced},
        }
    )
