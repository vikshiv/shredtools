# HPRCv2 region browser — same extract path as shredtools/extract_from_mums.py (no CLI, no files out).
import base64
import json
import re
from bisect import bisect_right

import bumbl_index_utils as sutils

# Shown under the page title after Pyodide loads (`index.html` renders Markdown).
# Split on the first blank line: text *after* it stays visible; the *first* paragraph is
# inside a collapsed **How does it work?** toggle. Leave INTRO empty to hide the block.
# Markdown is HTML in the browser (trusted — only you should edit this string).
INTRO = """
Shredtools uses **multi-MUMs**, collinear exact match markers along a pangenome, to identify and extract homologous regions to a query region of interest. Sometimes for a query region, the nearest flanking multi-MUM markers are some distance away from the requested interval. We report this distance on each side as "bounds". Bounds of 0 on each side indicate that the interval falls directly on multi-MUMs and the exact position of the corresponding sequence in each assembly can be found. This means the extracted region is *likely* a homologous sequence to your region of interest.

We use two different datasets from the Human Pangenome Reference Consortium (HPRC). Release 1 contains 92 assemblies and release 2 contains 476. For release 2, we also include an additional multi-MUM index with improved coverage, which we recommend as the default index.

All indexes are hosted thanks to the AWS Open Data Sponsorship Program and are freely available to query and for download and offline use.
"""
_S3_MUMEMTO = "https://genome-idx.s3.amazonaws.com/mumemto"

# Each preset: .bumbl and .bumbl.bi on S3. Contig/length metadata lives in ``pangenome_lengths.json``
# (built with ``lengths_to_json.py``) and is loaded at startup via ``load_lengths_bundle_path``.
# Order preserved for dropdown: enhanced default, then HPRCr2 merged, then HPRCr1.
PANGENOMES: dict[str, dict[str, str]] = {
    "hprcv2_enhanced": {
        "label": "HPRCr2 (enhanced)",
        "bumbl": f"{_S3_MUMEMTO}/hprcv2_enhanced_merged.bumbl",
        "bi": f"{_S3_MUMEMTO}/hprcv2_enhanced_merged.bumbl.bi",
    },
    "hprcv2_merged": {
        "label": "HPRCr2",
        "bumbl": f"{_S3_MUMEMTO}/hprcv2_merged.bumbl",
        "bi": f"{_S3_MUMEMTO}/hprcv2_merged.bumbl.bi",
    },
    "hprcv1": {
        "label": "HPRCr1",
        "bumbl": f"{_S3_MUMEMTO}/hprcv1.bumbl",
        "bi": f"{_S3_MUMEMTO}/hprcv1.bumbl.bi",
    },
}

ACTIVE_PANGENOME = "hprcv2_enhanced"

# Loaded once from ``pangenome_lengths.json`` (keys match ``PANGENOMES``).
_LENGTHS_BUNDLE: dict | None = None
# Index cache: "pangenome_key:seq_idx" → parsed multi-index document.
_INDEX_BY_SEQ: dict = {}

# Last successful extract (for synteny plot without re-fetching S3).
_LAST_EXTRACT: dict | None = None


def load_lengths_bundle_path(path: str) -> None:
    """
    Load merged lengths JSON into memory. Call from the browser after writing the file to MEMFS.

    Each top-level value must include ``seq_lengths_multi`` and ``contig_names`` (see ``lengths_to_json.py``).
    """
    global _LENGTHS_BUNDLE
    with open(path, "r", encoding="utf-8") as f:
        bundle = json.load(f)
    for k in PANGENOMES:
        if k not in bundle:
            raise ValueError(f"Lengths bundle missing preset key {k!r}")
        ent = bundle[k]
        sm = ent.get("seq_lengths_multi")
        cn = ent.get("contig_names")
        if not isinstance(sm, list) or not isinstance(cn, list):
            raise ValueError(f"Invalid lengths entry for {k!r}")
        if len(sm) != len(cn):
            raise ValueError(f"seq_lengths_multi / contig_names length mismatch for {k!r}")
    _LENGTHS_BUNDLE = bundle


def pangenome_options_json() -> str:
    """Dropdown options: key and UI label."""
    order = ["hprcv2_enhanced", "hprcv2_merged", "hprcv1"]
    out = []
    for k in order:
        if k in PANGENOMES:
            p = PANGENOMES[k]
            out.append({"key": k, "label": p["label"]})
    return json.dumps(out)


def set_active_pangenome(key: str) -> None:
    """Switch active S3 bumbl/bi pair; clears index cache only (lengths stay in the loaded bundle)."""
    global ACTIVE_PANGENOME
    if key not in PANGENOMES:
        raise ValueError(f"Unknown pangenome {key!r}")
    ACTIVE_PANGENOME = key
    _INDEX_BY_SEQ.clear()


def _active_bumbl_bi() -> str:
    return PANGENOMES[ACTIVE_PANGENOME]["bi"]


def _active_bumbl() -> str:
    return PANGENOMES[ACTIVE_PANGENOME]["bumbl"]


def _get_lengths_meta():
    if _LENGTHS_BUNDLE is None:
        raise RuntimeError("Lengths bundle not loaded; call load_lengths_bundle_path first.")
    ent = _LENGTHS_BUNDLE[ACTIVE_PANGENOME]
    seq_lengths_multi = ent["seq_lengths_multi"]
    contig_names = ent["contig_names"]
    n = len(seq_lengths_multi)
    return seq_lengths_multi, contig_names, n


def _index_cache_key(seq_idx: int) -> str:
    return f"{ACTIVE_PANGENOME}:{int(seq_idx)}"


async def _get_index(seq_idx: int):
    ck = _index_cache_key(seq_idx)
    idx = _INDEX_BY_SEQ.get(ck)
    if idx is not None:
        return idx, True
    idx = await sutils.parse_index(_active_bumbl_bi(), seq_idx=int(seq_idx))
    _INDEX_BY_SEQ[ck] = idx
    return idx, False


async def warm_index(seq_idx: int) -> bool:
    """Preload the index for seq_idx; used on genome dropdown change."""
    ck = _index_cache_key(seq_idx)
    idx = await sutils.parse_index(_active_bumbl_bi(), seq_idx=int(seq_idx))
    _INDEX_BY_SEQ[ck] = idx
    return True


def find_target_region(coll_mums, coords, seq_idx, sequences, right_key=None):
    """Aligned with ``shredtools/extract_from_mums.find_target_region`` (bisect_right + bounds)."""
    n = int(coll_mums.num_mums)
    if n <= 0:
        raise ValueError("No MUMs loaded for this query (empty range slice).")
    starts_col = coll_mums.starts_col(seq_idx)
    left_mum_idx = bisect_right(starts_col, coords[0]) - 1
    if right_key is None:
        right_key = [starts_col[i] + int(coll_mums.lengths[i]) for i in range(n)]
    right_mum_idx = bisect_right(right_key, coords[1])
    left_mum_idx = max(0, min(int(left_mum_idx), n - 1))
    right_mum_idx = max(0, min(int(right_mum_idx), n - 1))
    mum_bounds = (left_mum_idx, right_mum_idx)
    left_mum, right_mum = coll_mums[mum_bounds[0]], coll_mums[mum_bounds[1]]
    if not (
        coords[0] >= left_mum.starts[seq_idx]
        and coords[1] < right_mum.starts[seq_idx] + right_mum.length
    ):
        raise ValueError(
            "Loaded MUM slice does not bound the requested coordinates (try a wider region)."
        )
    left_offset, right_offset = 0, 0
    if coords[0] < left_mum.starts[seq_idx] + left_mum.length:
        left_offset = coords[0] - left_mum.starts[seq_idx]
    if coords[1] >= right_mum.starts[seq_idx]:
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
    """Match ``extract_from_mums`` stderr margins (exclusive outer bounds)."""
    n = int(coll_mums.num_mums)
    if n <= 0:
        raise ValueError("No MUMs loaded for this query (empty range slice).")
    starts_col = coll_mums.starts_col(seq_idx)
    left_mum_idx = bisect_right(starts_col, coords[0]) - 1
    if right_key is None:
        right_key = [starts_col[i] + int(coll_mums.lengths[i]) for i in range(n)]
    right_mum_idx = bisect_right(right_key, coords[1])
    left_mum_idx = max(0, min(int(left_mum_idx), n - 1))
    right_mum_idx = max(0, min(int(right_mum_idx), n - 1))
    left_mum, right_mum = coll_mums[left_mum_idx], coll_mums[right_mum_idx]

    left_bound = left_mum.starts[seq_idx] + left_mum.length - 1
    right_bound = right_mum.starts[seq_idx]

    left_offset, right_offset = 0, 0
    if coords[0] < left_mum.starts[seq_idx] + left_mum.length:
        left_offset = coords[0] - left_mum.starts[seq_idx]
        left_margin = 0
    else:
        left_margin = coords[0] - left_bound
    if coords[1] >= right_mum.starts[seq_idx]:
        right_offset = coords[1] - right_mum.starts[seq_idx]
        right_margin = 0
    else:
        right_margin = right_bound - coords[1]
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
    li = bisect_right(starts_col, coords[0]) - 1
    ri = bisect_right(right_key, coords[1])
    return (li >= 0), (ri < n)


async def _get_mums_expanding(idx, coords, seq_idx, max_steps: int = 8):
    """
    Fetch MUMs for bins around coords: first use index flanks + span bounds (extract path),
    sort by query column, then optionally widen bins if the slice still fails to bracket.
    """
    s, e = int(coords[0]), int(coords[1])
    bin_start = idx.coord_to_bin(s)
    bin_end = idx.coord_to_bin(e)

    got = await sutils.get_mum_ranges_flanks(idx, (s, e))
    if got is None:
        return None, [], (bin_start, bin_end)

    ranges, (left_bin, right_bin), _ = got
    steps = 0
    while True:
        if not ranges:
            return None, [], (bin_start, bin_end)

        mums = await sutils.parse_bumbl_range(_active_bumbl(), ranges)
        mums = sutils.sort_mums_by_seq_column(mums, seq_idx)
        starts_col = mums.starts_col(seq_idx)
        right_key = [starts_col[i] + int(mums.lengths[i]) for i in range(int(mums.num_mums))]
        left_ok, right_ok = _bracket_ok(mums, (s, e), seq_idx, right_key)
        if left_ok and right_ok:
            return (mums, right_key), ranges, (bin_start, bin_end)

        if steps >= int(max_steps):
            raise ValueError(
                f"Could not bracket region from flanking bins after {max_steps} expansions "
                f"(bin_start={bin_start}, bin_end={bin_end}, left_bin={left_bin}, right_bin={right_bin})."
            )

        if not left_ok:
            lb2 = idx.closest_nonzero_bin_left(left_bin - 1)
            if lb2 is None:
                left_ok = True
            else:
                left_bin = int(lb2)
                more = await idx.get_bins(left_bin)
                ranges = more + ranges
        if not right_ok:
            rb2 = idx.closest_nonzero_bin_right(right_bin + 1)
            if rb2 is None:
                right_ok = True
            else:
                right_bin = int(rb2)
                more = await idx.get_bins(right_bin)
                ranges = ranges + more
        steps += 1


def format_bed(
    contig_names, seq_lengths_multi, other_coords, sequences, path_placeholder="."
):
    """Same rows as extract_from_mums.extract_bed; last column is a placeholder in the browser."""
    lines = []
    for i, seq in enumerate(sequences):
        try:
            name, rel_offsets = sutils.convert_global_to_local_coords(
                other_coords[i][0],
                other_coords[i][1],
                contig_names[int(seq)],
                seq_lengths_multi[int(seq)],
            )
        except AssertionError:
            continue
        lines.append(
            f"{name}\t{int(rel_offsets[0])}\t{int(rel_offsets[1])}\t{path_placeholder}\n"
        )
    return "".join(lines)


def describe_ui() -> str:
    """JSON for genome / contig dropdowns (uses the loaded lengths bundle for ``ACTIVE_PANGENOME``)."""
    seq_lengths_multi, contig_names, n = _get_lengths_meta()
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


def _genome_labels() -> list[str]:
    seq_lengths_multi, contig_names, n = _get_lengths_meta()
    del seq_lengths_multi
    labels = []
    for i in range(n):
        cnames = contig_names[i]
        labels.append(cnames[0] if cnames else f"seq_{i}")
    return labels


def _seq_lengths_totals() -> list[int]:
    seq_lengths_multi, _cn, _n = _get_lengths_meta()
    return [sum(x) for x in seq_lengths_multi]


def plot_extract_png(seq_indices: list[int]) -> str:
    """
    Render extract synteny for the last ``run_with_bounds`` result.
    Returns JSON: ``{png_b64, n_mums, n_rows}`` or ``{error}``.
    Matplotlib must already be loaded in Pyodide (``loadPackage`` from JS).
    """
    global _LAST_EXTRACT
    if _LAST_EXTRACT is None:
        return json.dumps({"error": "No extract cached. Run a query first."})

    ctx = _LAST_EXTRACT
    seq_idx = int(ctx["seq_idx"])
    coords = ctx["coords"]
    mums = ctx["mums"]
    mum_bounds = ctx["mum_bounds"]
    other_coords = ctx["other_coords"]
    num_seqs = len(other_coords)

    if not seq_indices:
        return json.dumps({"error": "No sequences selected for plot."})

    seen = set()
    plot_seqs = []
    for s in seq_indices:
        si = int(s)
        if si < 0 or si >= num_seqs or si in seen:
            continue
        seen.add(si)
        plot_seqs.append(si)
    if seq_idx not in seen:
        plot_seqs.insert(0, seq_idx)
        seen.add(seq_idx)
    if not plot_seqs:
        return json.dumps({"error": "No valid sequence indices for plot."})

    labels_all = _genome_labels()
    genome_labels = [labels_all[s] for s in plot_seqs]
    seq_lengths = _seq_lengths_totals()

    try:
        import synteny_plot
    except ImportError as e:
        return json.dumps({"error": f"Plot module not available: {e}"})

    try:
        png = synteny_plot.plot_extract(
            coords,
            mums,
            mum_bounds,
            other_coords,
            seq_idx,
            plot_seqs,
            seq_lengths,
            genome_labels=genome_labels,
        )
    except Exception as e:
        return json.dumps({"error": str(e)})

    n_mums = int(mum_bounds[1]) - int(mum_bounds[0]) + 1
    return json.dumps(
        {
            "png_b64": base64.b64encode(png).decode("ascii"),
            "n_mums": n_mums,
            "n_rows": len(plot_seqs),
        }
    )


async def run(
    seq_idx: int,
    range_str: str,
    sequences=None,
) -> str:
    """
    Returns BED text or raises with a clear error message.
    `range_str` is chr:start-end (same as shredtools extract -r).
    """
    seq_lengths_multi, contig_names, num_seqs = _get_lengths_meta()
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
            f"No bounding MUMs found for region {range_str!r} (bins {requested_bins})."
        )
    mums, right_key = got
    _mum_bounds, other_coords = find_target_region(mums, coords, seq_idx, sequences, right_key=right_key)
    return format_bed(
        contig_names, seq_lengths_multi, other_coords, sequences, path_placeholder="."
    )


async def run_with_bounds(
    seq_idx: int,
    range_str: str,
    sequences=None,
) -> str:
    """
    Like `run`, but returns JSON: { bed: str, bounds: { contig, start, end } }.
    Bounds are the extracted interval for the selected genome (seq_idx), in contig-local coords.
    """
    global _LAST_EXTRACT
    _LAST_EXTRACT = None

    seq_lengths_multi, contig_names, num_seqs = _get_lengths_meta()
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
            f"No bounding MUMs found for region {range_str!r} (bins {requested_bins})."
        )
    # Each element of `ranges` is a half-open [mum_start, mum_end) slice.
    mum_chunks = int(len(ranges))
    mums_sliced = int(sum(int(b) - int(a) for a, b in ranges))

    mums, right_key = got
    mum_bounds, other_coords = find_target_region(mums, coords, seq_idx, sequences, right_key=right_key)
    _LAST_EXTRACT = {
        "seq_idx": int(seq_idx),
        "coords": coords,
        "mums": mums,
        "mum_bounds": mum_bounds,
        "other_coords": other_coords,
        "sequences": sequences,
    }
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
