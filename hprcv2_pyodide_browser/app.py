# HPRCv2 region browser — same data path as mod_scripts/index_extract.py (no CLI, no files out).
import json
from bisect import bisect_left

import bumbl_index_utils as sutils

BUMBL_BI_URL = "https://genome-idx.s3.amazonaws.com/mumemto/hprcv2_enhanced_merged.bumbl.bi"
BUMBL_URL = "https://genome-idx.s3.amazonaws.com/mumemto/hprcv2_enhanced_merged.bumbl"

async def warm_index(seq_idx: int) -> bool:
    """Preload the index for seq_idx; used on genome dropdown change."""
    await sutils.parse_index(BUMBL_BI_URL, seq_idx=int(seq_idx))
    return True


def find_target_region(coll_mums, coords, seq_idx, sequences):
    """From mod_scripts/index_extract.py (no stderr prints in browser)."""
    starts_col = coll_mums.starts_col(seq_idx)
    left_mum_idx = bisect_left(starts_col, coords[0]) - 1
    right_key = [starts_col[i] + int(coll_mums.lengths[i]) for i in range(coll_mums.num_mums)]
    right_mum_idx = bisect_left(right_key, coords[1])
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


def compute_margins(coll_mums, coords, seq_idx):
    """Match the left/right margin prints from mod_scripts/index_extract.py."""
    starts_col = coll_mums.starts_col(seq_idx)
    left_mum_idx = bisect_left(starts_col, coords[0]) - 1
    right_key = [starts_col[i] + int(coll_mums.lengths[i]) for i in range(coll_mums.num_mums)]
    right_mum_idx = bisect_left(right_key, coords[1])
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
    seq_lengths_multi = sutils.get_sequence_lengths(lengths_path, multilengths=True)
    contig_names = sutils.get_contig_names(lengths_path)
    n = len(seq_lengths_multi)
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
    seq_lengths_multi = sutils.get_sequence_lengths(lengths_path, multilengths=True)
    contig_names = sutils.get_contig_names(lengths_path)
    num_seqs = len(seq_lengths_multi)
    if not (0 <= seq_idx < num_seqs):
        raise ValueError(f"seq_idx {seq_idx} invalid (N = {num_seqs})")
    if sequences is not None and any(s >= num_seqs for s in sequences):
        raise ValueError(f"Invalid sequence index (N = {num_seqs})")
    if sequences is None:
        sequences = list(range(num_seqs))

    coords = sutils.convert_local_to_global_coords(
        range_str, contig_names[seq_idx], seq_lengths_multi[seq_idx]
    )
    idx = await sutils.parse_index(BUMBL_BI_URL, seq_idx=seq_idx)
    ranges, _snapped_bins, requested_bins = await sutils.get_mum_ranges_flanks(idx, coords)
    if len(ranges) == 0:
        raise ValueError(
            f"No MUM ranges for bins {requested_bins} (index could not find flanking bins)"
        )
    mums = await sutils.parse_bumbl_range(BUMBL_URL, ranges)
    _mum_bounds, other_coords = find_target_region(mums, coords, seq_idx, sequences)
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
    seq_lengths_multi = sutils.get_sequence_lengths(lengths_path, multilengths=True)
    contig_names = sutils.get_contig_names(lengths_path)
    num_seqs = len(seq_lengths_multi)
    if not (0 <= seq_idx < num_seqs):
        raise ValueError(f"seq_idx {seq_idx} invalid (N = {num_seqs})")
    # For the browser app we always compute all sequences; UI-side filtering is handled in JS.
    sequences = list(range(num_seqs))

    coords = sutils.convert_local_to_global_coords(
        range_str, contig_names[seq_idx], seq_lengths_multi[seq_idx]
    )
    idx = await sutils.parse_index(BUMBL_BI_URL, seq_idx=seq_idx)
    ranges, _snapped_bins, requested_bins = await sutils.get_mum_ranges_flanks(idx, coords)
    if len(ranges) == 0:
        raise ValueError(
            f"No MUM ranges for bins {requested_bins} (index could not find flanking bins)"
        )

    mums = await sutils.parse_bumbl_range(BUMBL_URL, ranges)
    _mum_bounds, other_coords = find_target_region(mums, coords, seq_idx, sequences)
    rows = []
    for i, seq in enumerate(sequences):
        name, rel_offsets = sutils.convert_global_to_local_coords(
            other_coords[i][0],
            other_coords[i][1],
            contig_names[int(seq)],
            seq_lengths_multi[int(seq)],
        )
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
    contig, rel = sutils.convert_global_to_local_coords(
        b0, b1, contig_names[seq_idx], seq_lengths_multi[seq_idx]
    )
    left_margin, right_margin = compute_margins(mums, coords, seq_idx)

    return json.dumps(
        {
            "rows": rows,
            "bounds": {"contig": contig, "start": int(rel[0]), "end": int(rel[1])},
            "margins": {"left": left_margin, "right": right_margin},
        }
    )
