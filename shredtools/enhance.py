#!/usr/bin/env python3

import argparse
import sys, os
import subprocess
import numpy as np
from multiprocessing import Pool
from collections import defaultdict

import mumemto.utils as mutils
from tqdm.auto import tqdm

from mumemto import mum as run_mumemto

_seq_paths = None
_slices_by_seq = None
_gap_seqs = None
_output_dir = None

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Enhance multi-MUM collection by finding local MUMs in gaps between collinear global MUMs."
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
        help="Output directory.",
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
        help="Number of threads to use.",
    )
    parser.add_argument(
        "--lengths",
        "-l",
        dest="lens",
        help="Lengths file. Default: <mum_file>.lengths",
    )
    parser.add_argument(
        "--outlier-factor",
        dest="outlier_factor",
        type=float,
        default=None,
        help="Outlier threshold factor for too-small gaps: mark if gap < (median - factor*MAD). Default: 5.",
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

    if not args.output.endswith('.mums') and not args.output.endswith('.bumbl'):
        args.output += '.mums'

    if args.threads < 1:
        print("--threads must be >= 1", file=sys.stderr)
        sys.exit(1)

    if args.outlier_factor is not None and args.outlier_factor < 0:
        print("--outlier-factor must be >= 0", file=sys.stderr)
        sys.exit(1)

    return args


def _subset_sequences(args, mums, seq_lengths, seq_paths):
    NUM_SEQS = len(seq_lengths)
    if args.sequences is not None and any([s >= NUM_SEQS or s < 0 for s in args.sequences]):
        print(f"Invalid sequence index in {args.sequences} (N = {NUM_SEQS})", file=sys.stderr)
        raise SystemExit(1)
    if args.sequences is None:
        args.sequences = list(range(NUM_SEQS))
        return mums, seq_lengths, seq_paths
    seq_lengths = [seq_lengths[i] for i in args.sequences]
    seq_paths = [seq_paths[i] for i in args.sequences]
    mums = mums[:, args.sequences]
    return mums, seq_lengths, seq_paths


def compute_gap_coords_and_lengths(mums):
    gaps = []
    for s, e in mums.blocks:
        for i in range(s, e):
            gaps.append((i, i + 1))
    if len(gaps) == 0:
        print("No collinear MUM gaps found.", file=sys.stderr)
        raise SystemExit(1)
    
    coords = []
    gap_lengths = []
    for l, r in gaps:
        left_starts = mums.starts[l]
        right_starts = mums.starts[r]
        left_lens = mums.lengths[l]
        right_lens = mums.lengths[r]

        flipped = right_starts < left_starts

        lens = np.where(flipped, right_lens, left_lens)
        gap_lengths.append(np.abs(left_starts - right_starts) - lens)

        left_coords = np.where(flipped, right_starts + right_lens, left_starts + left_lens)
        right_coords = np.where(flipped, left_starts, right_starts)
        coords.append(np.column_stack((left_coords, right_coords)))

    gap_lengths = np.asarray(gap_lengths)
    coords = np.asarray(coords)
    return gaps, coords, gap_lengths


def select_gaps(gap_lengths, min_gap_length):
    if gap_lengths.size == 0:
        return []
    return np.where((gap_lengths >= int(min_gap_length)).any(axis=1))[0].tolist()


def compute_outlier_mask(gap_lengths, outlier_factor=5.0):
    """
    Mark genomes that are *too small* for a gap vs per-gap median (med - outlier_factor*MAD).
    """
    if outlier_factor is None:
        return None
    if gap_lengths.size == 0:
        return np.zeros_like(gap_lengths, dtype=bool)
    med = np.median(gap_lengths, axis=1, keepdims=True)
    mad = np.median(np.abs(gap_lengths - med), axis=1, keepdims=True)
    mad = np.where(mad == 0, 1, mad)
    return gap_lengths < (med - float(outlier_factor) * mad)


def parse_single_fasta(path):
    with open(path, "r") as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))

def filter_gap_slices(coords, gap_ids, outlier_mask, seq_lengths):
    slices_by_seq = [{} for _ in range(len(seq_lengths))]
    for gap_id in gap_ids:
        c = coords[gap_id]
        for seq_local_idx in range(c.shape[0]):
            if outlier_mask is not None and outlier_mask[gap_id, seq_local_idx]:
                continue
            start, end = c[seq_local_idx]
            if end <= start:
                continue
            if start < 0 or end < 0:
                continue
            if end > seq_lengths[seq_local_idx]:
                continue
            slices_by_seq[seq_local_idx][gap_id] = (start, end)
    return slices_by_seq

def _init_collect_gap_sequences(seq_paths, slices_by_seq):
    global _seq_paths, _slices_by_seq
    _seq_paths = seq_paths
    _slices_by_seq = slices_by_seq


def _collect_gap_slices_for_one_seq(seq_local_idx):
    slices = _slices_by_seq[seq_local_idx]
    if not slices:
        return []
    fasta_path = _seq_paths[seq_local_idx]
    seq = parse_single_fasta(fasta_path)
    out = []
    for gap_id, (start, end) in slices.items():
        out.append((gap_id, seq_local_idx, start, seq[start:end]))
    return out

def collect_gap_sequences_from_fastas(seq_paths, slices_by_seq, gap_ids, threads):
    gap_seqs = {gap_id: [None for _ in range(len(seq_paths))] for gap_id in gap_ids}
    threads = min(int(threads), max(1, len(seq_paths)))
    with Pool(
        processes=threads,
        initializer=_init_collect_gap_sequences,
        initargs=(seq_paths, slices_by_seq),
    ) as pool:
        for items in tqdm(
            pool.imap_unordered(_collect_gap_slices_for_one_seq, range(len(seq_paths))),
            total=len(seq_paths),
            desc="load gap sequences",
        ):
            for gap_id, seq_local_idx, start, subseq in items:
                gap_seqs[gap_id][seq_local_idx] = (start, subseq)
    return gap_seqs

def enhance_gap(gap_seqs):
    missing = np.zeros(len(gap_seqs), dtype=bool)
    input_seqs = []
    present_idx = []
    present_offsets = []
    for i, payload in enumerate(gap_seqs):
        if payload is None:
            missing[i] = True
            continue
        start, s = payload
        present_idx.append(i)
        present_offsets.append(start)
        input_seqs.append([s])

    out_mums = run_mumemto(input_seqs, min_match_len=10)
    num_mums = len(out_mums)
    num_seqs_total = len(gap_seqs)
    lengths = np.empty(num_mums, dtype=np.uint32)
    starts = np.full((num_mums, num_seqs_total), -1, dtype=np.int64)
    strands = np.zeros((num_mums, num_seqs_total), dtype=bool)
    present_offsets = np.asarray(present_offsets, dtype=np.int64)

    for j in range(num_mums):
        length, offsets, match_strands = out_mums[j]
        lengths[j] = length
        starts[j, present_idx] = offsets + present_offsets
        strands[j, present_idx] = match_strands
        
    mums = mutils.MUMdata.from_arrays(lengths, starts, strands)
    mums.sort()
    strict_mums = mums[:, present_idx]
    blocks = mutils.find_coll_blocks(strict_mums)
    coll_mums = [i for s, e in blocks for i in range(s, e+1)]   
    return mums[coll_mums, :]

def main(args=None):
    args = parse_arguments(args)

    seq_lengths = mutils.get_sequence_lengths(args.lens)
    seq_paths = mutils.get_seq_paths(args.lens)
    
    mums = mutils.MUMdata(args.mum_file)
    
    mums, seq_lengths, seq_paths = _subset_sequences(args, mums, seq_lengths, seq_paths)

    if mums.blocks is None:
        mums.blocks = mutils.find_coll_blocks(mums)
        print(f"found blocks: {len(mums.blocks)}", file=sys.stderr)
    else:
        print(f"using provided blocks: {len(mums.blocks)}", file=sys.stderr)
    
    if len(mums) == 0 or mums.blocks is None or len(mums.blocks) == 0:
        print("No collinear MUMs found.", file=sys.stderr)
        raise SystemExit(1)

    gaps, coords, gap_lengths = compute_gap_coords_and_lengths(mums)
    gap_ids = np.where((gap_lengths >= int(args.min_gap_length)).any(axis=1))[0].tolist()

    print(f"Total MUM gap count: {len(gaps)}", file=sys.stderr)
    print(f"Gaps to enhance: {len(gap_ids)}", file=sys.stderr)


    outlier_mask = compute_outlier_mask(gap_lengths, outlier_factor=args.outlier_factor)
    slices_by_seq = filter_gap_slices(coords, gap_ids, outlier_mask, seq_lengths)

    gap_seqs = collect_gap_sequences_from_fastas(seq_paths, slices_by_seq, gap_ids, args.threads)

    print(f"loaded gap sequences for {len(gap_seqs)} gaps", file=sys.stderr)

    threads = min(int(args.threads), max(1, len(gap_ids)))
    with Pool(processes=threads) as pool:
        enhanced_mums = list(
            tqdm(
                pool.imap_unordered(enhance_gap, (gap_seqs[gap_id] for gap_id in gap_ids)),
                total=len(gap_ids),
                desc="Running mumemto on gaps",
            )
        )

    combined = mums.copy()
    for em in enhanced_mums:
        if em is None or len(em) == 0:
            continue
        combined = combined + em

    combined.sort(ref_col=0, copy=False)
    if args.output.endswith('.bumbl'):
        combined.write_bumbl(args.output)
    elif args.output.endswith('.mums'):
        combined.write_mums(args.output)
    print(f"wrote enhanced MUMs to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
