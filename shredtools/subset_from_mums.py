#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import mumemto.utils as mutils

from shredtools import utils as sutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Subset multi-MUM rows overlapping a query region on one assembly."    )
    parser.add_argument(
        "mum_file",
        type=str,
        help="Path or URL to input .bumbl file.",
    )
    parser.add_argument(
        "--seq-idx",
        "-s",
        type=int,
        help="Index of sequence with target region",
        required=True,
    )
    parser.add_argument(
        "--range",
        "-r",
        type=str,
        help="Region coordinates. Format -> contig:start-end",
        required=True,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output path (.mums or .bumbl). Writes .mums to stdout if omitted.",
    )
    parser.add_argument(
        "--format",
        "-f",
        type=str,
        choices=("mums", "bumbl"),
        default=None,
        help="Output format when writing to stdout or when -o has no extension.",
    )
    parser.add_argument(
        "--lengths",
        "-l",
        dest="lens",
        help="Path or URL to lengths file (default: <stem>.lengths)",
    )
    parser.add_argument(
        "--bumblbi",
        "-b",
        dest="bi",
        help="Path or URL to bumbl index (default: <mum_file>.bi)",
    )

    args = parser.parse_args(args)
    args.bi, args.lens = sutils.resolve_bumbl_sidecars(args.mum_file, args.bi, args.lens)

    return args


def _resolve_output_format(args) -> str:
    if args.format is not None:
        return args.format
    if args.output is None or args.output == "-":
        return "mums"
    if args.output.endswith(".bumbl"):
        return "bumbl"
    if args.output.endswith(".mums"):
        return "mums"
    return "mums"


def _resolve_output_dest(args) -> str:
    if args.output is None or args.output == "-":
        return "/dev/stdout"
    return args.output


def filter_mums_in_region(mums, seq_idx, region_start, region_end_excl):
    """Return MUMdata rows overlapping [region_start, region_end_excl) on seq_idx."""
    starts = mums.starts[:, seq_idx]
    present = starts >= 0
    overlap = present & (starts < region_end_excl) & (starts + mums.lengths > region_start)
    return mums[overlap]


def _empty_mumdata(num_seqs: int) -> mutils.MUMdata:
    return mutils.MUMdata.from_arrays(
        np.array([], dtype=np.uint32),
        np.empty((0, num_seqs), dtype=np.int64),
        np.empty((0, num_seqs), dtype=bool),
    )


def main(args=None):
    args = parse_arguments(args)

    with sutils.local_file(args.lens) as lens_path:
        seq_lengths_multi = mutils.get_sequence_lengths(lens_path, multilengths=True)
        contig_names = mutils.get_contig_names(lens_path)
        num_seqs = len(seq_lengths_multi)
        if not (0 <= args.seq_idx < num_seqs):
            print(
                f"Sequence index {args.seq_idx} is invalid (N = {num_seqs})",
                file=sys.stderr,
            )
            raise SystemExit(1)

        coords = sutils.convert_local_to_global_coords(
            args.range,
            contig_names[args.seq_idx],
            seq_lengths_multi[args.seq_idx],
        )
        region_start, region_end_incl = coords
        region_end_excl = region_end_incl + 1

        idx = sutils.parse_index(args.bi, seq_idx=args.seq_idx)
        ranges = sutils.get_mum_ranges_region(idx, coords)
        if ranges is None:
            print(
                f"No indexed bins found for region {args.range}.",
                file=sys.stderr,
            )
            subset = _empty_mumdata(num_seqs)
        else:
            mums = sutils.parse_bumbl_range(args.mum_file, ranges)
            print(f"loaded {len(mums)} candidate MUM rows from index", file=sys.stderr)
            subset = filter_mums_in_region(
                mums, args.seq_idx, region_start, region_end_excl
            )
            subset.sort(args.seq_idx)

        print(f"kept {len(subset)} MUM rows overlapping {args.range}", file=sys.stderr)

        out_fmt = _resolve_output_format(args)
        dest = _resolve_output_dest(args)

        if out_fmt == "mums":
            subset.write_mums(dest, blocks=None)
        else:
            subset.write_bums(dest, blocks=None)

        if dest != "/dev/stdout":
            print(f"wrote {len(subset)} MUM rows to {dest}", file=sys.stderr)


if __name__ == "__main__":
    main()
