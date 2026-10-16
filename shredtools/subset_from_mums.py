#!/usr/bin/env python3

import argparse
import os
import sys

import numpy as np

import mumemto.utils as mutils

from shredtools import utils as sutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Subset multi-MUM rows overlapping a region (or BED) on one assembly."
    )
    parser.add_argument("mum_file", type=str, help="Path or URL to input MUM file (.mums or .bumbl).")

    region = parser.add_mutually_exclusive_group(required=True)
    region.add_argument("--range", "-r", default=None,
                        help="Single region coordinates. Format -> contig:start-end.")
    region.add_argument("--bed", "-B", default=None,
                        help="BED file of regions in the target genome (many intervals).")

    parser.add_argument("--seq-idx", "-s", type=int, required=True,
                        help="Index of the genome the region/BED refers to (filelist order).")
    parser.add_argument("--mode", "-m", choices=["keep", "exclude", "both"], default="keep",
                        help="keep = only overlapping MUMs (default); exclude = drop overlapping "
                             "MUMs; both = write both.")
    parser.add_argument("--format", "-f", choices=["mums", "bumbl", "auto"], default=None,
                        help="Output format. Default infers from the output extension, else 'mums'. "
                             "'auto' follows the input file's format.")
    parser.add_argument("--output", "-o", default=None,
                        help="Output path (single mode) or prefix -> <prefix>.<mode>.<ext> "
                             "(--mode both). Required for --mode both and for bumbl output to a "
                             "non-file destination; text single-mode goes to stdout if omitted.")
    parser.add_argument("--lengths", "-l", dest="lens", default=None,
                        help="Path or URL to (multi)lengths file (default: <mum_file>.lengths).")
    parser.add_argument("--bumblbi", "-b", dest="bi", default=None,
                        help="Path or URL to bumbl index (default: <mum_file>.bi). Enables the "
                             "indexed fast path for --mode keep with --range.")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print progress messages to stderr (implied when writing to a file).")

    args = parser.parse_args(args)

    # Sanity checks: file existence, lengths file, bumbl index, output requirements.
    if not sutils.is_url(args.mum_file) and not os.path.exists(args.mum_file):
        print(f"MUM file {args.mum_file} not found", file=sys.stderr)
        raise SystemExit(1)
    if args.bed is not None and not os.path.exists(args.bed):
        print(f"BED file {args.bed} not found", file=sys.stderr)
        raise SystemExit(1)

    if args.lens is None:
        args.lens = os.path.splitext(args.mum_file)[0] + ".lengths"
    if not sutils.is_url(args.lens) and not os.path.exists(args.lens):
        print(f"Lengths file {args.lens} not found, and no lengths file provided", file=sys.stderr)
        raise SystemExit(1)

    # Bumbl index is optional: used only to accelerate keep+range on a .bumbl input.
    if args.bi is None and args.mum_file.endswith(".bumbl"):
        cand = args.mum_file + ".bi"
        if sutils.is_url(cand) or os.path.exists(cand):
            args.bi = cand
    if args.bi is not None and not sutils.is_url(args.bi) and not os.path.exists(args.bi):
        print(f"Bumbl index {args.bi} not found", file=sys.stderr)
        raise SystemExit(1)

    if args.mode == "both" and args.output is None:
        print("--output prefix is required when --mode both", file=sys.stderr)
        raise SystemExit(1)
    return args


def _resolve_output_format(args) -> str:
    """Return 'mums' or 'bumbl' for the output."""
    if args.format in ("mums", "bumbl"):
        return args.format
    if args.format == "auto":
        return "bumbl" if args.mum_file.endswith(".bumbl") else "mums"
    # Not specified: infer from the output name, else default to text.
    out = args.output
    if out and out != "-":
        if out.endswith(".bumbl"):
            return "bumbl"
        if out.endswith(".mums"):
            return "mums"
    return "mums"


def build_contig_offsets(contig_names, contig_lengths):
    """Return {contig_name: cumulative_offset} for one sequence."""
    offsets, cum = {}, 0
    for name, length in zip(contig_names, contig_lengths):
        offsets[name] = cum
        cum += length
    return offsets


def load_bed(bed_file, offsets, contig_lengths):
    """Load BED file and convert to merged global [start, end) intervals. Returns (starts, ends)."""
    ivs = []
    for lineno, line in enumerate(open(bed_file), 1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split()
        chrom, start, end = f[0], int(f[1]), int(f[2])
        if chrom not in offsets:
            sys.exit(f"BED contig '{chrom}' not found in sequence's contigs: "
                     f"{sorted(offsets)[:10]}...")
        length = contig_lengths[chrom]
        if not (0 <= start < end <= length):
            sys.exit(f"Invalid BED interval on line {lineno}: {chrom}:{start}-{end} "
                     f"is invalid for contig {chrom} with length {length} "
                     f"(need 0 <= start < end <= length)")
        base = offsets[chrom]
        ivs.append((base + start, base + end))
    return _merge_intervals(ivs)


def _merge_intervals(ivs):
    """Sort and merge overlapping intervals. Returns (starts, ends) of merged intervals."""
    ivs = sorted(ivs)
    merged = []
    for s, e in ivs:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    m_starts = np.array([s for s, _ in merged], dtype=np.int64)
    m_ends = np.array([e for _, e in merged], dtype=np.int64)
    return m_starts, m_ends


def overlap_mask(pos, lengths, m_starts, m_ends):
    """Return boolean mask of which MUMs overlap the merged intervals."""
    n = len(pos)
    hit = np.zeros(n, dtype=bool)
    if len(m_starts) == 0:
        return hit
    ends = pos + lengths.astype(np.int64)
    # index of the last merged interval starting at/before pos
    idx = np.searchsorted(m_starts, pos, side="right") - 1
    valid = idx >= 0
    left = np.zeros(n, dtype=bool)
    left[valid] = m_ends[idx[valid]] > pos[valid]          # pos falls inside that interval
    nxt = idx + 1
    has_next = nxt < len(m_starts)
    right = np.zeros(n, dtype=bool)
    right[has_next] = m_starts[nxt[has_next]] < ends[has_next]  # MUM extends into the next interval
    hit = left | right
    hit[pos < 0] = False
    return hit


def _empty_mumdata(num_seqs):
    """Return an empty MUMdata object with the given number of sequences."""
    return mutils.MUMdata.from_arrays(
        np.array([], dtype=np.uint32),
        np.empty((0, num_seqs), dtype=np.int64),
        np.empty((0, num_seqs), dtype=bool),
    )


def _write(mumdata, path, fmt):
    """Write MUMdata to path in the given format."""
    if fmt == "bumbl":
        mumdata.write_bums(path, blocks=None)
    else:
        mumdata.write_mums(path, blocks=None)


def main(args=None):
    args = parse_arguments(args)

    out_fmt = _resolve_output_format(args)
    ext = out_fmt
    single_dest = args.output if (args.output and args.output != "-") else "/dev/stdout"
    log_progress = args.verbose or (args.output is not None and args.output != "-")

    if out_fmt == "bumbl" and args.mode != "both" and single_dest == "/dev/stdout":
        print("bumbl output requires --output (cannot stream binary to stdout)", file=sys.stderr)
        raise SystemExit(1)

    with sutils.local_file(args.lens) as lens_path:
        contig_names = mutils.get_contig_names(lens_path)
        contig_lengths = mutils.get_sequence_lengths(lens_path, multilengths=True)
        num_seqs = len(contig_names)
        if not (0 <= args.seq_idx < num_seqs):
            print(f"Sequence index {args.seq_idx} is invalid (N = {num_seqs})", file=sys.stderr)
            raise SystemExit(1)

        # Build the merged global [start, end) intervals to match against.
        coords = None
        if args.range is not None:
            try:
                coords = sutils.convert_local_to_global_coords(
                    args.range, contig_names[args.seq_idx], contig_lengths[args.seq_idx]
                )
            except AssertionError as e:
                print(f"Invalid region {args.range}: {e}", file=sys.stderr)
                raise SystemExit(1)
            region_start, region_end_incl = coords
            m_starts = np.array([region_start], dtype=np.int64)
            m_ends = np.array([region_end_incl + 1], dtype=np.int64)
            region_desc = args.range
        else:
            names_i, lengths_i = contig_names[args.seq_idx], contig_lengths[args.seq_idx]
            offsets = build_contig_offsets(names_i, lengths_i)
            m_starts, m_ends = load_bed(args.bed, offsets, dict(zip(names_i, lengths_i)))
            region_desc = f"{len(m_starts)} BED intervals ({os.path.basename(args.bed)})"
            if log_progress:
                print(f"loaded {len(m_starts)} merged BED intervals for seq {args.seq_idx}",
                      file=sys.stderr)

        # Indexed fast path: keep-only, single range, .bumbl input, index available.
        use_index = (
            args.mode == "keep"
            and args.range is not None
            and args.mum_file.endswith(".bumbl")
            and args.bi is not None
        )

        if use_index:
            idx = sutils.parse_index(args.bi, seq_idx=args.seq_idx)
            ranges = sutils.get_mum_ranges_region(idx, coords)
            if ranges is None:
                if log_progress:
                    print(f"No indexed bins found for region {region_desc}.", file=sys.stderr)
                subset = _empty_mumdata(num_seqs)
            else:
                mums = sutils.parse_bumbl_range(args.mum_file, ranges)
                if log_progress:
                    print(f"loaded {len(mums)} candidate MUM rows from index", file=sys.stderr)
                pos = mums.starts[:, args.seq_idx].astype(np.int64)
                hit = overlap_mask(pos, mums.lengths, m_starts, m_ends)
                subset = mums.slice(hit)
                subset.sort(args.seq_idx)
            if log_progress:
                print(f"kept {len(subset)} MUM rows overlapping {region_desc}", file=sys.stderr)
            _write(subset, single_dest, out_fmt)
            if log_progress and single_dest != "/dev/stdout":
                print(f"wrote {len(subset)} MUM rows to {single_dest}", file=sys.stderr)
            return

        # Full-file path: required for exclude/both, BED input, .mums input, or no index.
        if sutils.is_url(args.mum_file):
            print("Reading the full MUM file is not supported for URL inputs; use --mode keep "
                  "--range with a .bumbl index instead.", file=sys.stderr)
            raise SystemExit(1)

        md = mutils.MUMdata(args.mum_file, sort=False)
        if md.starts.shape[1] != num_seqs:
            print(f"[warn] MUM file has {md.starts.shape[1]} sequences but lengths file has "
                  f"{num_seqs}; check that they match.", file=sys.stderr)

        pos = md.starts[:, args.seq_idx].astype(np.int64)
        hit = overlap_mask(pos, md.lengths, m_starts, m_ends)

        if args.mode in ("keep", "both"):
            sub = md.slice(hit)
            dest = f"{args.output}.keep.{ext}" if args.mode == "both" else single_dest
            _write(sub, dest, out_fmt)
            if log_progress and dest != "/dev/stdout":
                print(f"wrote {len(sub)} MUM rows overlapping {region_desc} to {dest}",
                      file=sys.stderr)
        if args.mode in ("exclude", "both"):
            sub = md.slice(~hit)
            dest = f"{args.output}.exclude.{ext}" if args.mode == "both" else single_dest
            _write(sub, dest, out_fmt)
            if log_progress and dest != "/dev/stdout":
                print(f"wrote {len(sub)} MUM rows not overlapping {region_desc} to {dest}",
                      file=sys.stderr)

        if log_progress:
            print(f"kept {int(hit.sum())} / excluded {int((~hit).sum())} MUM rows "
                  f"(of {len(hit)} total)", file=sys.stderr)


if __name__ == "__main__":
    main()
