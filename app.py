# HPRCv2 region browser — same extract path as shredtools/extract_from_mums.py (no CLI, no files out).
import base64
import json
import re
from bisect import bisect_right

import bumbl_index_utils as sutils

_S3_MUMEMTO = "https://genome-idx.s3.amazonaws.com/mumemto"
_S3_HPRC_ASSEMBLIES = "https://human-pangenomics.s3.us-west-2.amazonaws.com/"

# Each preset: .bumbl, .bumbl.bi, and .lengths on S3. Contig/length metadata for the browser
# also lives in ``data/pangenome_lengths.json`` (built with ``lengths_to_json.py``).
# Order preserved for dropdown: enhanced default, then HPRCr2 merged, then HPRCr1.
PANGENOMES: dict[str, dict[str, str]] = {
    "hprcv2_enhanced": {
        "label": "HPRCr2 (enhanced)",
        "bumbl": f"{_S3_MUMEMTO}/hprcv2_enhanced_merged.bumbl",
        "bi": f"{_S3_MUMEMTO}/hprcv2_enhanced_merged.bumbl.bi",
        "lengths": f"{_S3_MUMEMTO}/hprcv2_enhanced_merged.lengths",
    },
    "hprcv2_merged": {
        "label": "HPRCr2",
        "bumbl": f"{_S3_MUMEMTO}/hprcv2_merged.bumbl",
        "bi": f"{_S3_MUMEMTO}/hprcv2_merged.bumbl.bi",
        "lengths": f"{_S3_MUMEMTO}/hprcv2_merged.lengths",
    },
    "hprcv1": {
        "label": "HPRCr1",
        "bumbl": f"{_S3_MUMEMTO}/hprcv1.bumbl",
        "bi": f"{_S3_MUMEMTO}/hprcv1.bumbl.bi",
        "lengths": f"{_S3_MUMEMTO}/hprcv1.lengths",
    },
}

# Built-in HPRC presets only (custom remote datasets are registered at runtime).
_BUILTIN_PANGENOME_KEYS = ("hprcv2_enhanced", "hprcv2_merged", "hprcv1")


def _bumbl_stem(bumbl_url: str) -> str:
    return bumbl_url.rsplit("/", 1)[-1].removesuffix(".bumbl")


INTRO_VISIBLE_SEP = "\n\n---\n\n"


def _build_intro() -> str:
    """Markdown for the page intro (see comment on ``INTRO``)."""
    order = ["hprcv2_enhanced", "hprcv2_merged", "hprcv1"]
    collapsed = [
        "Shredtools uses **multi-MUMs**, collinear exact match markers along a pangenome, "
        "to identify and extract homologous regions to a query region of interest. Sometimes "
        "for a query region, the nearest flanking multi-MUM markers are some distance away "
        'from the requested interval. We report this distance on each side as "bounds". '
        "Bounds of 0 on each side indicate that the interval falls directly on multi-MUMs and "
        "the exact position of the corresponding sequence in each assembly can be found. This "
        "means the extracted region is *likely* a homologous sequence to your region of interest.",
        f"**Index files** — hosted at [`s3://genome-idx/mumemto/`]({_S3_MUMEMTO}/). "
        "Each pangenome preset uses a `.bumbl`, `.bumbl.bi`, and `.lengths` trio:",
    ]
    for key in order:
        p = PANGENOMES[key]
        stem = _bumbl_stem(p["bumbl"])
        collapsed.append(f"- **{p['label']}** (`{stem}`):")
        collapsed.append(f"  - [`{stem}.bumbl`]({p['bumbl']})")
        collapsed.append(f"  - [`{stem}.bumbl.bi`]({p['bi']})")
        collapsed.append(f"  - [`{stem}.lengths`]({p['lengths']})")

    collapsed.append("")
    collapsed.append(
        f"FASTA sequences are pulled directly from HPRC hosted assemblies "
        f"([`s3://human-pangenomics/`]({_S3_HPRC_ASSEMBLIES}))."
    )

    visible = [
        "We use two different datasets from the Human Pangenome Reference Consortium (HPRC). "
        "Release 1 contains 92 assemblies and release 2 contains 476. For release 2, we also "
        "include an additional multi-MUM index with improved coverage, which we recommend as "
        "the default index.",
        "All indexes are hosted thanks to the AWS Open Data Sponsorship Program and are freely "
        "available to query and for download and offline use.",
    ]
    return "\n".join(collapsed) + INTRO_VISIBLE_SEP + "\n\n".join(visible)


# Shown under the page title after Pyodide loads (`index.html` renders Markdown).
# Split on ``INTRO_VISIBLE_SEP``: text before it is inside **How does it work?**;
# text after stays visible. Leave INTRO empty to hide the block.
# Markdown is HTML in the browser (trusted — only you should edit this string).
INTRO = _build_intro()

ACTIVE_PANGENOME = "hprcv2_enhanced"

# Loaded once from ``data/pangenome_lengths.json`` (keys match ``PANGENOMES``).
_LENGTHS_BUNDLE: dict | None = None
# Index cache: "pangenome_key:seq_idx" → parsed multi-index document.
_INDEX_BY_SEQ: dict = {}

# Last successful extract (BED from flanks; plot reloads bins left_bin..right_bin).
_LAST_EXTRACT: dict | None = None


def load_lengths_bundle_path(path: str) -> None:
    """
    Load merged lengths JSON into memory. Call from the browser after writing the file to MEMFS.

    Each top-level value must include ``seq_lengths_multi`` and ``contig_names`` (see ``lengths_to_json.py``).
    Only built-in HPRC keys are required; custom datasets are added via ``register_custom_pangenome``.
    """
    global _LENGTHS_BUNDLE
    with open(path, "r", encoding="utf-8") as f:
        bundle = json.load(f)
    for k in _BUILTIN_PANGENOME_KEYS:
        if k not in bundle:
            raise ValueError(f"Lengths bundle missing preset key {k!r}")
        ent = bundle[k]
        sm = ent.get("seq_lengths_multi")
        cn = ent.get("contig_names")
        if not isinstance(sm, list) or not isinstance(cn, list):
            raise ValueError(f"Invalid lengths entry for {k!r}")
        if len(sm) != len(cn):
            raise ValueError(f"seq_lengths_multi / contig_names length mismatch for {k!r}")
    # Preserve any already-registered custom entries if this is called again.
    if _LENGTHS_BUNDLE is not None:
        for k, ent in _LENGTHS_BUNDLE.items():
            if k not in _BUILTIN_PANGENOME_KEYS:
                bundle[k] = ent
    _LENGTHS_BUNDLE = bundle


def register_custom_pangenome(
    key: str, label: str, bumbl: str, bi: str, lengths_text: str
) -> None:
    """
    Register a remote custom dataset from multilengths text + bumbl/bi URLs.

    Contig metadata is parsed from ``lengths_text`` (same format as Mumemto ``.lengths``).
    """
    global _LENGTHS_BUNDLE, ACTIVE_PANGENOME, _LAST_EXTRACT
    if _LENGTHS_BUNDLE is None:
        raise RuntimeError("Lengths bundle not loaded; call load_lengths_bundle_path first.")
    k = str(key or "").strip()
    if not k:
        raise ValueError("Custom pangenome key is required.")
    if k in _BUILTIN_PANGENOME_KEYS:
        raise ValueError(f"Cannot overwrite built-in pangenome {k!r}")
    lab = str(label or "").strip() or k
    bumbl_url = str(bumbl or "").strip()
    bi_url = str(bi or "").strip()
    if not bumbl_url or not bi_url:
        raise ValueError("bumbl and bi URLs are required.")
    text = str(lengths_text or "")
    sm = sutils.get_sequence_lengths_from_text(text, multilengths=True, label=k)
    cn = sutils.get_contig_names_from_text(text, label=k)
    asm = sutils.get_assembly_names_from_text(text, label=k)
    if len(sm) != len(cn):
        raise ValueError(f"seq_lengths_multi / contig_names length mismatch for {k!r}")
    if len(asm) != len(sm):
        raise ValueError(f"assembly_labels / seq_lengths_multi length mismatch for {k!r}")
    if not sm:
        raise ValueError(f"No sequences found in lengths for {k!r}")
    PANGENOMES[k] = {
        "label": lab,
        "bumbl": bumbl_url,
        "bi": bi_url,
        "lengths": "",
    }
    _LENGTHS_BUNDLE[k] = {
        "seq_lengths_multi": sm,
        "contig_names": cn,
        "assembly_labels": asm,
    }
    # Drop cached indexes for this key (re-register / update).
    stale = [ck for ck in _INDEX_BY_SEQ if ck.startswith(f"{k}:")]
    for ck in stale:
        del _INDEX_BY_SEQ[ck]
    if ACTIVE_PANGENOME == k:
        _LAST_EXTRACT = None


def unregister_custom_pangenome(key: str) -> None:
    """Remove a previously registered custom dataset."""
    global ACTIVE_PANGENOME, _LAST_EXTRACT
    k = str(key or "").strip()
    if k in _BUILTIN_PANGENOME_KEYS:
        raise ValueError(f"Cannot unregister built-in pangenome {k!r}")
    if k not in PANGENOMES:
        return
    del PANGENOMES[k]
    if _LENGTHS_BUNDLE is not None and k in _LENGTHS_BUNDLE:
        del _LENGTHS_BUNDLE[k]
    stale = [ck for ck in _INDEX_BY_SEQ if ck.startswith(f"{k}:")]
    for ck in stale:
        del _INDEX_BY_SEQ[ck]
    if ACTIVE_PANGENOME == k:
        ACTIVE_PANGENOME = "hprcv2_enhanced"
        _LAST_EXTRACT = None
        _INDEX_BY_SEQ.clear()


def pangenome_options_json() -> str:
    """Dropdown options for built-in HPRC presets: key and UI label."""
    out = []
    for k in _BUILTIN_PANGENOME_KEYS:
        if k in PANGENOMES:
            p = PANGENOMES[k]
            out.append({"key": k, "label": p["label"]})
    return json.dumps(out)


def set_active_pangenome(key: str) -> None:
    """Switch active bumbl/bi pair; clears index cache and last extract (lengths stay in the bundle)."""
    global ACTIVE_PANGENOME, _LAST_EXTRACT
    if key not in PANGENOMES:
        raise ValueError(f"Unknown pangenome {key!r}")
    ACTIVE_PANGENOME = key
    _INDEX_BY_SEQ.clear()
    _LAST_EXTRACT = None


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
    assembly_labels = ent.get("assembly_labels")
    if not isinstance(assembly_labels, list) or len(assembly_labels) != n:
        assembly_labels = None
    return seq_lengths_multi, contig_names, n, assembly_labels


def _genome_label_for(i: int, contig_names, assembly_labels) -> str:
    """Prefer cleaned assembly labels (custom datasets); else first contig (HPRC PanSN)."""
    if assembly_labels is not None and i < len(assembly_labels):
        lab = str(assembly_labels[i] or "").strip()
        if lab:
            return lab
    cnames = contig_names[i] if i < len(contig_names) else None
    return cnames[0] if cnames else f"seq_{i}"


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

    Returns
    -------
    (mums, right_key), ranges, (bin_start, bin_end), (left_bin, right_bin)
        On success. ``(left_bin, right_bin)`` are the final expanded flank bins (inclusive).
    None, [], (bin_start, bin_end), None
        When no usable flanks.
    """
    s, e = int(coords[0]), int(coords[1])
    bin_start = idx.coord_to_bin(s)
    bin_end = idx.coord_to_bin(e)

    got = await sutils.get_mum_ranges_flanks(idx, (s, e))
    if got is None:
        return None, [], (bin_start, bin_end), None

    ranges, (left_bin, right_bin), _ = got
    steps = 0
    while True:
        if not ranges:
            return None, [], (bin_start, bin_end), None

        mums = await sutils.parse_bumbl_range(_active_bumbl(), ranges)
        mums = sutils.sort_mums_by_seq_column(mums, seq_idx)
        starts_col = mums.starts_col(seq_idx)
        right_key = [starts_col[i] + int(mums.lengths[i]) for i in range(int(mums.num_mums))]
        left_ok, right_ok = _bracket_ok(mums, (s, e), seq_idx, right_key)
        if left_ok and right_ok:
            return (mums, right_key), ranges, (bin_start, bin_end), (int(left_bin), int(right_bin))

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


def describe_ui() -> str:
    """JSON for genome / contig dropdowns (uses the loaded lengths bundle for ``ACTIVE_PANGENOME``)."""
    seq_lengths_multi, contig_names, n, assembly_labels = _get_lengths_meta()
    genomes = []
    for i in range(n):
        cnames = contig_names[i]
        genomes.append(
            {
                "seq_idx": i,
                "label": _genome_label_for(i, contig_names, assembly_labels),
                "contigs": cnames,
                "contig_lengths": [int(x) for x in seq_lengths_multi[i]],
            }
        )
    return json.dumps({"n_seqs": n, "genomes": genomes})


def _genome_labels() -> list[str]:
    _sm, contig_names, n, assembly_labels = _get_lengths_meta()
    return [_genome_label_for(i, contig_names, assembly_labels) for i in range(n)]


def _seq_lengths_totals() -> list[int]:
    seq_lengths_multi, _cn, _n, _asm = _get_lengths_meta()
    return [sum(x) for x in seq_lengths_multi]


def _chunk_to_plot_mums(chunk, seq_idx, bound_lo, bound_hi, seq_list, other_coords):
    """
    Project MUM rows in ``[bound_lo, bound_hi]`` on ``seq_idx`` onto ``seq_list`` columns,
    subtracting each sequence's ``other_coords`` window start (window-relative coords).
    Returns ``(plot_mums, n_kept)``.
    """
    from array import array

    out_lengths = array("I")
    out_starts = array("q")
    out_strands = bytearray()
    n_kept = 0
    bound_lo = int(bound_lo)
    bound_hi = int(bound_hi)
    for i in range(int(chunk.num_mums)):
        st = chunk.start(i, seq_idx)
        if st < 0 or st < bound_lo or st > bound_hi:
            continue
        out_lengths.append(int(chunk.lengths[i]))
        for src_seq in seq_list:
            x = chunk.start(i, src_seq)
            if x != -1:
                x = x - int(other_coords[src_seq][0])
            out_starts.append(int(x))
            out_strands.append(1 if chunk.strand(i, src_seq) else 0)
        n_kept += 1
    return sutils.MUMdata(len(seq_list), out_lengths, out_starts, out_strands), n_kept


async def plot_extract_png(seq_indices: list[int], dark: bool = False) -> str:
    """
    Render extract synteny for the last ``run_with_bounds`` result.

    Loads MUM bins ``left_bin..right_bin`` one at a time (subset-style byte ranges),
    converts each chunk to polygons, then discards the chunk so peak memory stays ~one bin.

    Returns JSON: ``{png_b64, n_mums, n_rows}`` or ``{error}``.
    Matplotlib must already be loaded in Pyodide (``loadPackage`` from JS).
    """
    global _LAST_EXTRACT
    if _LAST_EXTRACT is None:
        return json.dumps({"error": "No extract cached. Run a query first."})

    ctx = _LAST_EXTRACT
    seq_idx = int(ctx["seq_idx"])
    coords = ctx["coords"]
    other_coords = ctx["other_coords"]
    left_bin = int(ctx["left_bin"])
    right_bin = int(ctx["right_bin"])
    bound_lo, bound_hi = ctx["bound_starts"]
    bound_lo, bound_hi = int(bound_lo), int(bound_hi)
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
        import viz_mums
    except ImportError as e:
        return json.dumps({"error": f"Plot module not available: {e}"})

    try:
        idx, _ = await _get_index(seq_idx)
        polygons = []
        colors = []
        n_mums = 0
        centering = [0] * len(plot_seqs)

        for b in range(left_bin, right_bin + 1):
            if idx.bin_is_empty(b):
                continue
            ranges = await idx.get_bins(b)
            if not ranges:
                continue
            chunk = await sutils.parse_bumbl_range(_active_bumbl(), ranges)
            plot_mums, n_kept = _chunk_to_plot_mums(
                chunk, seq_idx, bound_lo, bound_hi, plot_seqs, other_coords
            )
            del chunk
            if n_kept == 0:
                continue
            poly, cols = viz_mums.get_mum_polygons(
                plot_mums, centering, inv_color="green"
            )
            del plot_mums
            polygons.extend(poly)
            colors.extend(cols)
            n_mums += n_kept

        png = synteny_plot.plot_extract_polygons(
            coords,
            other_coords,
            seq_idx,
            plot_seqs,
            seq_lengths,
            polygons,
            colors,
            genome_labels=genome_labels,
            dark=bool(dark),
        )
    except Exception as e:
        return json.dumps({"error": str(e)})

    return json.dumps(
        {
            "png_b64": base64.b64encode(png).decode("ascii"),
            "n_mums": int(n_mums),
            "n_rows": len(plot_seqs),
        }
    )


def _run_error_message(exc: BaseException) -> str:
    """User-facing message for expected extract failures (no Pyodide traceback)."""
    msg = exc.args[0] if exc.args else str(exc)
    s = str(msg)
    if "No bounding MUMs" in s:
        return "No flanking multi-MUMs were found for this interval."
    if "Selected-genome region spans multiple contigs" in s:
        return "No flanking multi-MUM between interval and contig end."
    if "Loaded MUM slice does not bound" in s:
        return (
            "MUM markers near this interval do not fully bracket the request. "
            "Try a wider region."
        )
    return s


async def run_with_bounds(
    seq_idx: int,
    range_str: str,
    sequences=None,
) -> str:
    """
    Like `run`, but returns JSON: { bed: str, bounds: { contig, start, end } }.
    Bounds are the extracted interval for the selected genome (seq_idx), in contig-local coords.
    On failure returns ``{"error": "<message>"}`` instead of raising.
    """
    global _LAST_EXTRACT
    _LAST_EXTRACT = None

    try:
        seq_lengths_multi, contig_names, num_seqs, assembly_labels = _get_lengths_meta()
        if not (0 <= seq_idx < num_seqs):
            raise ValueError(f"seq_idx {seq_idx} invalid (N = {num_seqs})")
        # For the browser app we always compute all sequences; UI-side filtering is handled in JS.
        sequences = list(range(num_seqs))

        coords = sutils.convert_local_to_global_coords(
            range_str, contig_names[seq_idx], seq_lengths_multi[seq_idx]
        )
        idx, _ = await _get_index(seq_idx)
        got, ranges, requested_bins, flank_bins = await _get_mums_expanding(
            idx, coords, seq_idx
        )
        if got is None:
            raise ValueError(
                f"No bounding MUMs found for region {range_str!r} (bins {requested_bins})."
            )
        # Each element of `ranges` is a half-open [mum_start, mum_end) slice.
        mum_chunks = int(len(ranges))
        mums_sliced = int(sum(int(b) - int(a) for a, b in ranges))

        mums, right_key = got
        mum_bounds, other_coords = find_target_region(
            mums, coords, seq_idx, sequences, right_key=right_key
        )
        left_bin, right_bin = flank_bins
        bound_starts = (
            int(mums.start(mum_bounds[0], seq_idx)),
            int(mums.start(mum_bounds[1], seq_idx)),
        )
        _LAST_EXTRACT = {
            "seq_idx": int(seq_idx),
            "coords": coords,
            "mums": mums,
            "mum_bounds": mum_bounds,
            "other_coords": other_coords,
            "sequences": sequences,
            "left_bin": int(left_bin),
            "right_bin": int(right_bin),
            "bound_starts": bound_starts,
        }
        rows = []
        unavailable = []
        _span_re = re.compile(
            r"start and end coords are in different contigs:\s+(.+?)\s+and\s+(.+)$"
        )
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
                label = _genome_label_for(int(seq), contig_names, assembly_labels)
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
            # MUM-bounded interval spans contigs on the selected genome.
            msg = e.args[0] if e.args else str(e)
            raise ValueError(
                f"Selected-genome region spans multiple contigs: {msg}"
            ) from None
        left_margin, right_margin = compute_margins(
            mums, coords, seq_idx, right_key=right_key
        )

        return json.dumps(
            {
                "rows": rows,
                "bounds": {"contig": contig, "start": int(rel[0]), "end": int(rel[1])},
                "margins": {"left": left_margin, "right": right_margin},
                "unavailable": unavailable,
                "mum_slices": {"chunks": mum_chunks, "mums": mums_sliced},
            }
        )
    except (ValueError, AssertionError) as e:
        return json.dumps({"error": _run_error_message(e)})
