#!/usr/bin/env python3

import argparse
import sys
import os

import numpy as np
import pysam
from pysam import samtools as pysam_samtools
from multiprocessing import Pool

import mumemto.utils as mutils
from tqdm.auto import tqdm

from mumemto import mum as run_mumemto
from shredtools import utils as sutils

_CTX: dict = {}


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description=(
            "Enhance multi-MUM collection by finding local MUMs in gaps "
            "between collinear global MUMs (optionally also block flanks)."
        )
    )
    parser.add_argument("mum_file", type=str, help="Path to the input MUM file.")
    parser.add_argument(
        "min_gap_length",
        type=int,
        help="Minimum gap length (bp). If any genome has a gap >= this, the gap is enhanced.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="output",
        help="Output file (.mums or .bumbl).",
    )
    parser.add_argument(
        "--sequences",
        "-s",
        type=int,
        nargs="*",
        default=None,
        help="Sequence indices to include. Default: all.",
    )
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Number of parallel workers (one gap per task).",
    )
    parser.add_argument(
        "--lengths",
        "-l",
        dest="lens",
        help="Multilengths file. Default: <mum_file>.lengths",
    )
    parser.add_argument(
        "--flanks",
        action="store_true",
        help=(
            "Also search terminal flanks (before the first block and after "
            "the last) and gaps between consecutive collinear blocks. "
            "Off by default."
        ),
    )
    parser.add_argument(
        "--outlier-factor",
        dest="outlier_factor",
        type=float,
        default=None,
        help=(
            "Outlier threshold factor for too-small gaps: mark if gap < "
            "(median - factor*MAD). Off by default (no outlier detection); "
            "suggested value: 5."
        ),
    )
    parser.add_argument(
        "--no-partial",
        dest="no_partial",
        action="store_true",
        help=(
            "Do not add partial MUMs when some assemblies lack a valid gap "
            "sequence. Only valid when --outlier-factor is not set."
        ),
    )
    parser.add_argument(
        "--min-match-len",
        dest="min_match_len",
        type=int,
        default=10,
        help="Minimum MUM length for gap-local mumemto (default: 10).",
    )

    args = parser.parse_args(args)

    if not os.path.exists(args.mum_file):
        print(f"MUM file {args.mum_file} not found", file=sys.stderr)
        sys.exit(1)

    if args.lens is None:
        args.lens = os.path.splitext(args.mum_file)[0] + ".lengths"
        if not os.path.exists(args.lens):
            print(
                f"Lengths file {args.lens} not found, and no lengths file provided",
                file=sys.stderr,
            )
            sys.exit(1)

    if not args.output.endswith(".mums") and not args.output.endswith(".bumbl"):
        args.output += ".mums"

    if args.threads < 1:
        print("--threads must be >= 1", file=sys.stderr)
        sys.exit(1)

    if args.outlier_factor is not None and args.outlier_factor < 0:
        print("--outlier-factor must be >= 0", file=sys.stderr)
        sys.exit(1)

    if args.no_partial and args.outlier_factor is not None:
        print(
            "--no-partial conflicts with --outlier-factor "
            "(outlier mode intentionally produces partial MUMs)",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.min_match_len < 1:
        print("--min-match-len must be >= 1", file=sys.stderr)
        sys.exit(1)

    return args


def _load_lengths_layout(lengths_path):
    seq_paths = mutils.get_seq_paths(lengths_path)
    seq_lengths_multi = mutils.get_sequence_lengths(lengths_path, multilengths=True)
    contig_names = mutils.get_contig_names(lengths_path)
    total_lens = [sum(int(x) for x in row) for row in seq_lengths_multi]
    return seq_lengths_multi, contig_names, total_lens, seq_paths


def ensure_fai_indexes(seq_paths):
    seen = set()
    for path in seq_paths:
        if path in seen:
            continue
        seen.add(path)
        if not os.path.isfile(path):
            print(f"FASTA not found: {path}", file=sys.stderr)
            raise SystemExit(1)
        fai = path + ".fai"
        if not os.path.isfile(fai) or os.path.getmtime(fai) < os.path.getmtime(path):
            try:
                pysam_samtools.faidx(path)
            except Exception as e:
                print(f"pysam.samtools.faidx failed for {path}: {e}", file=sys.stderr)
                raise SystemExit(1) from e


def _subset_sequences(
    args,
    mums,
    seq_lengths_multi,
    contig_names,
    total_lens,
    seq_paths,
):
    num_seqs = len(total_lens)
    if args.sequences is not None and any(s >= num_seqs or s < 0 for s in args.sequences):
        print(f"Invalid sequence index in {args.sequences} (N = {num_seqs})", file=sys.stderr)
        raise SystemExit(1)
    if args.sequences is None:
        args.sequences = list(range(num_seqs))
        return mums, seq_lengths_multi, contig_names, total_lens, seq_paths
    seq_lengths_multi = [seq_lengths_multi[i] for i in args.sequences]
    contig_names = [contig_names[i] for i in args.sequences]
    total_lens = [total_lens[i] for i in args.sequences]
    seq_paths = [seq_paths[i] for i in args.sequences]
    mums = mums[:, args.sequences]
    return mums, seq_lengths_multi, contig_names, total_lens, seq_paths


def _pair_gap_coords(mums, l, r):
    """Orientation-aware [left, right) coords and lengths between MUM indices l and r."""
    left_starts = mums.starts[l]
    right_starts = mums.starts[r]
    left_lens = mums.lengths[l]
    right_lens = mums.lengths[r]

    flipped = right_starts < left_starts
    lens = np.where(flipped, right_lens, left_lens)
    gap_len = np.abs(left_starts - right_starts) - lens
    left_coords = np.where(flipped, right_starts + right_lens, left_starts + left_lens)
    right_coords = np.where(flipped, left_starts, right_starts)
    return np.column_stack((left_coords, right_coords)), gap_len


def _block_coverage_endpoints(mums, s, e):
    """Per-genome leftmost start and rightmost end of MUMs s..e (orientation-aware)."""
    first_starts = mums.starts[s]
    last_starts = mums.starts[e]
    first_lens = mums.lengths[s]
    last_lens = mums.lengths[e]
    flipped = last_starts < first_starts
    leftmost = np.where(flipped, last_starts, first_starts).astype(np.int64, copy=False)
    rightmost = np.where(
        flipped, first_starts + first_lens, last_starts + last_lens
    ).astype(np.int64, copy=False)
    return leftmost, rightmost


def compute_gap_coords_and_lengths(mums, include_flanks=False, total_lens=None):
    """
    Build per-interval [left, right) coords and lengths for collinear gaps.

    Internal gaps are consecutive MUM pairs (i, i+1) within each block.

    With include_flanks=True, also append (once each):
      - left terminal:  [0, start of first MUM in the first block)
      - inter-block:    gaps between last MUM of block_k and first of block_{k+1}
      - right terminal: [end of last MUM in the last block, genome_len)
    Intermediate blocks do not get flanks back to 0 / out to genome end.
    """
    if include_flanks and total_lens is None:
        raise ValueError("total_lens is required when include_flanks=True")

    num_seqs = int(mums.starts.shape[1]) if len(mums) else 0
    genome_lens = (
        np.asarray(total_lens, dtype=np.int64) if total_lens is not None else None
    )

    gaps = []
    coords = []
    gap_lengths = []

    blocks = list(mums.blocks)
    for s, e in blocks:
        for i in range(s, e):
            pair_coords, gap_len = _pair_gap_coords(mums, i, i + 1)
            coords.append(pair_coords)
            gap_lengths.append(gap_len)
            gaps.append((i, i + 1))

    if include_flanks and blocks:
        s0, e0 = blocks[0]
        sN, eN = blocks[-1]

        # Left terminal: [0, leftmost start of first block)
        leftmost, _ = _block_coverage_endpoints(mums, s0, e0)
        left_lo = np.zeros(num_seqs, dtype=np.int64)
        left_hi = leftmost
        coords.append(np.column_stack((left_lo, left_hi)))
        gap_lengths.append(np.maximum(left_hi - left_lo, 0))
        gaps.append(("left", s0, e0))

        # Inter-block gaps between consecutive collinear blocks
        for bi in range(len(blocks) - 1):
            _s_a, e_a = blocks[bi]
            s_b, _e_b = blocks[bi + 1]
            pair_coords, gap_len = _pair_gap_coords(mums, e_a, s_b)
            coords.append(pair_coords)
            gap_lengths.append(gap_len)
            gaps.append(("inter", e_a, s_b))

        # Right terminal: [rightmost end of last block, genome_len)
        _, rightmost = _block_coverage_endpoints(mums, sN, eN)
        right_lo = rightmost
        right_hi = genome_lens
        coords.append(np.column_stack((right_lo, right_hi)))
        gap_lengths.append(np.maximum(right_hi - right_lo, 0))
        gaps.append(("right", sN, eN))

    if len(gaps) == 0:
        print("No collinear MUM gaps found.", file=sys.stderr)
        raise SystemExit(1)

    gap_lengths = np.asarray(gap_lengths)
    coords = np.asarray(coords)
    return gaps, coords, gap_lengths


def compute_outlier_mask(gap_lengths, outlier_factor=None):
    """
    Mark genomes that are *too small* for a gap vs per-gap median (med - outlier_factor*MAD).
    Returns None when outlier_factor is None (no outlier detection).
    """
    if outlier_factor is None:
        return None
    if gap_lengths.size == 0:
        return np.zeros_like(gap_lengths, dtype=bool)
    med = np.median(gap_lengths, axis=1, keepdims=True)
    mad = np.median(np.abs(gap_lengths - med), axis=1, keepdims=True)
    mad = np.where(mad == 0, 1, mad)
    return gap_lengths < (med - float(outlier_factor) * mad)


def build_gap_work(gap_ids, coords, outlier_mask, num_seqs):
    work = []
    for gap_id in gap_ids:
        row = coords[gap_id]
        if outlier_mask is None:
            keep = np.ones(num_seqs, dtype=bool)
        else:
            keep = ~outlier_mask[gap_id]
        work.append((int(gap_id), row, keep))
    return work


def _worker_init(
    seq_paths,
    seq_lengths_multi,
    contig_names_per_seq,
    total_lens_per_seq,
    min_match_len,
    num_seqs,
    no_partial,
):
    global _CTX
    fastas = [pysam.FastaFile(p) for p in seq_paths]
    _CTX = {
        "fastas": fastas,
        "seq_lengths_multi": [list(x) for x in seq_lengths_multi],
        "contig_names_per_seq": [list(x) for x in contig_names_per_seq],
        "total_lens_per_seq": list(total_lens_per_seq),
        "min_match_len": int(min_match_len),
        "num_seqs": int(num_seqs),
        "no_partial": bool(no_partial),
    }


def _fetch_gap_intervals(coords_row, keep_mask):
    """Return list length N: None or (global_g0, subseq) per assembly."""
    ctx = _CTX
    fastas = ctx["fastas"]
    total_lens = ctx["total_lens_per_seq"]
    seq_lengths_multi = ctx["seq_lengths_multi"]
    contig_names = ctx["contig_names_per_seq"]
    n = int(ctx["num_seqs"])
    gap_seq = [None] * n
    coords_row = np.asarray(coords_row, dtype=np.int64)

    for j in range(n):
        if not bool(keep_mask[j]):
            continue
        g0, g1 = int(coords_row[j, 0]), int(coords_row[j, 1])
        if g1 <= g0 or g0 < 0 or g1 < 0 or g1 > int(total_lens[j]):
            continue
        try:
            cname, rel = sutils.convert_global_to_local_coords(
                g0,
                g1 - 1,
                contig_names[j],
                seq_lengths_multi[j],
            )
        except AssertionError:
            continue
        lo, hi_incl = int(rel[0]), int(rel[1])
        subseq = fastas[j].fetch(cname, lo, hi_incl + 1)
        if len(subseq) == 0:
            continue
        gap_seq[j] = (g0, subseq)
    return gap_seq


def enhance_gap(gap_seqs, min_match_len, no_partial=False):
    input_seqs = []
    present_idx = []
    present_offsets = []
    for i, payload in enumerate(gap_seqs):
        if payload is None:
            continue
        start, s = payload
        present_idx.append(i)
        present_offsets.append(start)
        input_seqs.append([s])

    num_seqs_total = len(gap_seqs)
    if not input_seqs:
        return None
    # Some assemblies had invalid/empty gap sequence → would yield partial MUMs.
    if no_partial and len(present_idx) < num_seqs_total:
        return None

    out_mums = run_mumemto(input_seqs, min_match_len=int(min_match_len))
    num_mums = len(out_mums)
    lengths = np.empty(num_mums, dtype=np.uint32)
    starts = np.full((num_mums, num_seqs_total), -1, dtype=np.int64)
    strands = np.zeros((num_mums, num_seqs_total), dtype=bool)
    present_offsets_a = np.asarray(present_offsets, dtype=np.int64)

    for j in range(num_mums):
        length, offsets, match_strands = out_mums[j]
        lengths[j] = length
        starts[j, present_idx] = offsets + present_offsets_a
        strands[j, present_idx] = match_strands

    mums = mutils.MUMdata.from_arrays(lengths, starts, strands)
    mums.sort()
    strict_mums = mums[:, present_idx]
    blocks = mutils.find_coll_blocks(strict_mums)
    coll_mums = [i for s, e in blocks for i in range(s, e + 1)]
    return mums[coll_mums, :]


def _worker_one(task):
    _gap_id, coords_row, keep_mask = task
    ctx = _CTX
    payloads = _fetch_gap_intervals(coords_row, keep_mask)
    return enhance_gap(payloads, ctx["min_match_len"], no_partial=ctx["no_partial"])


def main(args=None):
    args = parse_arguments(args)

    seq_lengths_multi, contig_names, total_lens, seq_paths = _load_lengths_layout(args.lens)

    mums = mutils.MUMdata(args.mum_file)

    mums, seq_lengths_multi, contig_names, total_lens, seq_paths = _subset_sequences(
        args, mums, seq_lengths_multi, contig_names, total_lens, seq_paths
    )

    num_seqs = len(total_lens)
    ensure_fai_indexes(seq_paths)

    if mums.blocks is None:
        mums.blocks = mutils.find_coll_blocks(mums)
        print(f"found blocks: {len(mums.blocks)}", file=sys.stderr)
    else:
        print(f"using provided blocks: {len(mums.blocks)}", file=sys.stderr)

    if len(mums) == 0 or mums.blocks is None or len(mums.blocks) == 0:
        print("No collinear MUMs found.", file=sys.stderr)
        raise SystemExit(1)

    gaps, coords, gap_lengths = compute_gap_coords_and_lengths(
        mums,
        include_flanks=bool(args.flanks),
        total_lens=total_lens if args.flanks else None,
    )
    gap_ids = np.where((gap_lengths >= int(args.min_gap_length)).any(axis=1))[0].tolist()

    kind = "gaps + flanks" if args.flanks else "gap"
    print(f"Total MUM {kind} count: {len(gaps)}", file=sys.stderr)
    print(f"Gaps to enhance: {len(gap_ids)}", file=sys.stderr)

    outlier_mask = compute_outlier_mask(gap_lengths, outlier_factor=args.outlier_factor)
    work = build_gap_work(gap_ids, coords, outlier_mask, num_seqs)

    initargs = (
        seq_paths,
        seq_lengths_multi,
        contig_names,
        total_lens,
        int(args.min_match_len),
        num_seqs,
        bool(args.no_partial),
    )

    threads = min(int(args.threads), max(1, len(work)))
    with Pool(processes=threads, initializer=_worker_init, initargs=initargs) as pool:
        enhanced_mums = list(
            tqdm(
                pool.imap_unordered(_worker_one, work),
                total=len(work),
                desc="enhance gaps",
            )
        )

    combined = mums.copy()
    for em in enhanced_mums:
        if em is None or len(em) == 0:
            continue
        combined = combined + em

    combined.sort(ref_col=0, copy=False)
    if args.output.endswith(".bumbl"):
        combined.write_bumbl(args.output)
    elif args.output.endswith(".mums"):
        combined.write_mums(args.output)
    print(f"wrote enhanced MUMs to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
