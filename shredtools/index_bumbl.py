#!/usr/bin/env python3

"""
Build a ``.bumbl.bi`` index for a ``.bumbl`` file.

Verifies that MUM rows are sorted (if needed), then writes ``<bumbl_path>.bi`` by default.
"""

import argparse
import os
import sys

from mumemto.utils import MUMdata

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import utils as shredutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Verify bumbl row order and write a .bumbl.bi index (single or multi)."
    )
    parser.add_argument("bumbl_file", type=str, help="Path to the input .bumbl file.")
    parser.add_argument(
        "--seq-idx",
        "-s",
        type=int,
        default=0,
        help="Sequence column in starts used for sortedness check and single-index axis (default: 0).",
    )
    fmt = parser.add_mutually_exclusive_group()
    fmt.add_argument(
        "--multi",
        action="store_true",
        help="Write multi-index format (default).",
    )
    fmt.add_argument(
        "--single",
        action="store_true",
        help="Write single-index format.",
    )
    parser.add_argument(
        "--bin-width",
        "-w",
        type=int,
        default=1_000_000,
        help="Genomic bin width in bp (default: 1_000_000).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output index path (default: <bumbl_file>.bi).",
    )
    parser.add_argument(
        "--no-verify-sorted",
        action="store_true",
        help="Skip streaming check that offsets are sorted.",
    )
    parser.add_argument(
        "--no-verify-checksum",
        action="store_true",
        help="Skip post-write verify_bumbl_index checksum check.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print verbose status messages and tqdm progress (if available) to stderr.",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Only verify an existing index checksum, then exit.",
    )
    args = parser.parse_args(args)

    if not os.path.exists(args.bumbl_file):
        print(f"Bumbl file {args.bumbl_file} not found", file=sys.stderr)
        sys.exit(1)

    if not args.bumbl_file.endswith(".bumbl"):
        print("Warning: input does not end with .bumbl", file=sys.stderr)

    if args.output is None:
        args.output = args.bumbl_file + ".bi"

    # Default: multi-index; use --single for single-index
    args.multi_index = not args.single

    return args


def main(args=None):
    args = parse_arguments(args)

    def log(msg):
        if args.verbose:
            print(msg, file=sys.stderr)

    log(f"Loading bumbl: {args.bumbl_file}")
    mums = MUMdata(args.bumbl_file, sort=False)
    if mums.num_mums == 0:
        print("Cannot build index: no MUMs in file (n_mums=0).", file=sys.stderr)
        sys.exit(1)

    if not (0 <= args.seq_idx < mums.num_seqs):
        print(
            f"seq_idx {args.seq_idx} invalid for num_seqs={mums.num_seqs}",
            file=sys.stderr,
        )
        sys.exit(1)

    if not args.no_verify_sorted and mums.num_mums > 1:
        log(f"Checking sortedness: starts[:, {args.seq_idx}]")
        col = mums.starts[:, args.seq_idx]
        if not (col[1:] >= col[:-1]).all():
            print(
                f"Sortedness check failed: starts[:, {args.seq_idx}] is not non-decreasing.",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.multi_index:
        log("Building multi-index")
        index = shredutils.build_bumbl_multiindex(
            mums, bin_width=args.bin_width, verbose=args.verbose
        )
        fmt_name = "multi"
    else:
        log("Building single-index")
        index = shredutils.build_bumbl_singleindex(
            mums, args.seq_idx, bin_width=args.bin_width
        )
        fmt_name = "single"

    log(f"Writing index: {args.output}")
    shredutils.write_bumbl_index(args.output, index)

    print(
        f"Wrote {fmt_name} index: {args.output} "
        f"(n_mums={mums.num_mums}, n_seqs={mums.num_seqs}, "
        f"bin_width={args.bin_width}, seq_idx={args.seq_idx})",
        file=sys.stderr,
    )

    if not args.no_verify_checksum:
        try:
            log("Verifying checksum")
            shredutils.verify_bumbl_index(args.bumbl_file, args.output)
        except AssertionError as e:
            print(f"Index verification failed: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
