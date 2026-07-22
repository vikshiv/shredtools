#!/usr/bin/env python3

import argparse
import os
import sys

import mumemto.utils as mutils
import mumemto.viz_mums
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from mumemto.utils import MUMdata, get_seq_paths

from shredtools import utils as sutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Shred a MUM file into smaller fragments and optionally visualize."
    )
    parser.add_argument("mum_file", type=str, help="Path to the input MUM file.")
    parser.add_argument(
        "--shred-size",
        type=int,
        default=10_000_000,
        help="Target size of each shred in base pairs (default: 10Mbp).",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="If set, generate a visualization of the shreds.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="shreds",
        help="Output directory for shredded MUMs.",
    )
    parser.add_argument(
        "--sequences",
        "-s",
        type=int,
        nargs="*",
        default=None,
        help="Sequence indices to include in output. Default: all.",
    )
    parser.add_argument(
        "--lengths",
        "-l",
        dest="lens",
        help="Lengths file (multilengths format).",
    )
    parser.add_argument(
        "--per-sequence-bed",
        action="store_true",
        help="Write one BED per sequence (sequence{i}.bed). Default: one BED per shred (shred{j}.bed).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress messages to stderr.",
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

    return args


def _global_contig_idx(pos, seq_lengths):
    contig_idx, _ = sutils.find_chr(np.asarray([pos], dtype=np.int64), seq_lengths)
    return int(contig_idx[0])


def _mum_end_at_contig_boundary(end_pos, seq_lengths):
    """True if exclusive global end_pos sits on a contig or genome boundary."""
    cum = np.cumsum(seq_lengths)
    return end_pos in cum


def _contig_break_at_mum(row, next_row, mums, seq_lengths_multi, num_seqs):
    """Return True if any present sequence ends a contig at this MUM."""
    for s in range(num_seqs):
        if mums.starts[row, s] < 0:
            continue
        end_pos = int(mums.starts[row, s]) + int(mums.lengths[row])
        if _mum_end_at_contig_boundary(end_pos, seq_lengths_multi[s]):
            return True
        if next_row is not None and mums.starts[next_row, s] >= 0:
            end_contig = _global_contig_idx(end_pos, seq_lengths_multi[s])
            next_contig = _global_contig_idx(
                int(mums.starts[next_row, s]), seq_lengths_multi[s]
            )
            if end_contig != next_contig:
                return True
    return False


def _safe_convert(global_start, global_end_excl, contig_names, seq_lengths_multi, fasta_path):
    """Convert half-open global [start, end) to contig-local BED coords."""
    if global_end_excl <= global_start:
        return None
    try:
        return sutils.convert_global_to_local_coords(
            global_start, global_end_excl - 1, contig_names, seq_lengths_multi
        )
    except AssertionError as e:
        reason = e.args[0] if e.args else str(e)
        print(
            f"Skipping BED line for {os.path.basename(fasta_path)}: {reason}",
            file=sys.stderr,
        )
        return None


def coll_row_indices(mums):
    return np.array(
        [i for s, e in mums.blocks for i in range(s, e + 1)],
        dtype=np.intp,
    )


def shred_mums(mums, shred_size, seq_lengths_multi, verbose=False):
    if mums.blocks is None:
        mums.blocks = mutils.find_coll_blocks(mums, verbose=verbose)
        if verbose:
            print(f"found blocks: {len(mums.blocks)}", file=sys.stderr)
    elif verbose:
        print(f"using provided blocks: {len(mums.blocks)}", file=sys.stderr)

    coll_rows = coll_row_indices(mums)
    if len(coll_rows) == 0:
        print("No collinear MUMs found.", file=sys.stderr)
        sys.exit(1)

    num_seqs = mums.strands.shape[1]
    shred_starts = np.zeros((num_seqs, 1), dtype=np.int64)
    genome_ends = np.array(
        [sum(seq_lengths_multi[s]) for s in range(num_seqs)], dtype=np.int64
    )

    for k, row in enumerate(coll_rows):
        ends = mums.starts[row, :].astype(np.int64) + int(mums.lengths[row])
        valid = mums.starts[row, :] >= 0
        next_row = coll_rows[k + 1] if k + 1 < len(coll_rows) else None

        cut_size = False
        if np.all(mums.strands[row, :]) and valid.any():
            span_mean = (ends[valid] - shred_starts[valid, -1]).mean()
            cut_size = span_mean > shred_size

        cut_contig = _contig_break_at_mum(
            row, next_row, mums, seq_lengths_multi, num_seqs
        )

        if cut_size or cut_contig:
            shred_starts = np.column_stack((shred_starts, ends.reshape(-1, 1)))

    if not np.all(shred_starts[:, -1] == genome_ends):
        shred_starts = np.column_stack((shred_starts, genome_ends.reshape(-1, 1)))
    return mums, shred_starts


def _write_bed_line(handle, global_start, global_end, contig_names, seq_lengths_multi, fasta_path):
    converted = _safe_convert(
        global_start, global_end, contig_names, seq_lengths_multi, fasta_path
    )
    if converted is None:
        return
    name, rel_offsets = converted
    handle.write(f"{name}\t{rel_offsets[0]}\t{rel_offsets[1] + 1}\t{fasta_path}\n")


def shred_bed(args, shredded_mums, contig_names, seq_lengths_multi):
    paths = get_seq_paths(args.lens)
    num_shreds = shredded_mums.shape[1] - 1

    if args.per_sequence_bed:
        for seq_i in args.sequences:
            bed_path = os.path.join(args.output, f"sequence{seq_i}.bed")
            with open(bed_path, "w") as bed_file:
                for j in range(num_shreds):
                    _write_bed_line(
                        bed_file,
                        int(shredded_mums[seq_i, j]),
                        int(shredded_mums[seq_i, j + 1]),
                        contig_names[seq_i],
                        seq_lengths_multi[seq_i],
                        paths[seq_i],
                    )
    else:
        for j in range(num_shreds):
            bed_path = os.path.join(args.output, f"shred{j}.bed")
            with open(bed_path, "w") as bed_file:
                for seq_i in args.sequences:
                    _write_bed_line(
                        bed_file,
                        int(shredded_mums[seq_i, j]),
                        int(shredded_mums[seq_i, j + 1]),
                        contig_names[seq_i],
                        seq_lengths_multi[seq_i],
                        paths[seq_i],
                    )


def plot(genome_lengths, polygons, colors, centering, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    for idx, g in enumerate(genome_lengths):
        ax.plot(
            [centering[idx] + 0, centering[idx] + g],
            [idx, idx],
            alpha=0.2,
            linewidth=0.75,
            c="black",
        )

    ax.add_collection(
        PolyCollection(
            polygons, linewidths=0, alpha=0.8, edgecolors=colors, facecolors=colors
        )
    )

    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis="y", which="both", length=0)
    if genomes:
        ax.set_yticklabels(genomes)
    else:
        ax.yaxis.set_ticklabels([])

    ax.set_xlabel("genomic position")
    ax.set_ylabel("sequences")
    ax.set_ylim(-0.25, len(genome_lengths) - 1 + 0.25)
    ax.set_xlim(0, max_length)
    fig.set_tight_layout(True)
    ax.invert_yaxis()

    if size:
        fig.set_size_inches(*size)
    return fig, ax


def plot_shreds(args, mums, shred_starts):
    num_seqs = mums.strands.shape[1]
    seq_lengths = [sum(x) for x in mutils.get_sequence_lengths(args.lens, multilengths=True)]
    centering = [0] * num_seqs

    poly, colors = mumemto.viz_mums.get_block_polygons(
        mums.blocks, mums, centering, inv_color="green"
    )
    fig, ax = plot(seq_lengths, poly, colors, centering, size=(10, 15))

    for i in range(shred_starts.shape[1]):
        ax.plot(shred_starts[:, i], np.arange(num_seqs), color="red", linewidth=2)

    fig.savefig(os.path.join(args.output, "shred_synteny.pdf"))


def main(args=None):
    args = parse_arguments(args)
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    seq_lengths_multi = mutils.get_sequence_lengths(args.lens, multilengths=True)
    num_seqs = len(seq_lengths_multi)
    if args.sequences is not None and any(s >= num_seqs for s in args.sequences):
        print(
            f"Sequence index {max(args.sequences)} is invalid (N = {num_seqs})",
            file=sys.stderr,
        )
        sys.exit(1)
    if args.sequences is None:
        args.sequences = list(range(num_seqs))

    mums = MUMdata(args.mum_file)
    mums, shredded_mums = shred_mums(
        mums, args.shred_size, seq_lengths_multi, verbose=args.verbose
    )
    contig_names = mutils.get_contig_names(args.lens)

    shred_bed(args, shredded_mums, contig_names, seq_lengths_multi)

    if args.plot:
        plot_shreds(args, mums, shredded_mums)


if __name__ == "__main__":
    main()
