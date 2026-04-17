#!/usr/bin/env python3

"""
Sort a ``.bumbl`` file by start position and write the result.

Notes:
Collinear blocks are dropped upon sorting.
Sorting is only for strict multi-MUMs, no partial MUMs or MEMs.
"""

import argparse
import os
import sys
import tempfile

from mumemto.utils import MUMdata

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import utils as shredutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Sort a .bumbl file by start position and write output."
    )
    parser.add_argument("bumbl_file", type=str, help="Path to the input .bumbl file.")
    parser.add_argument(
        "--seq-idx",
        "-s",
        type=int,
        default=0,
        help="Sequence column to sort by (default: 0).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output .bumbl path (default: <input>.sorted.bumbl).",
    )
    parser.add_argument(
        "--in-place",
        "-i",
        action="store_true",
        help="Rewrite the input file after writing to a temporary file.",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Only verify sortedness of starts[:, seq_idx], then exit.",
    )
    args = parser.parse_args(args)

    if not os.path.exists(args.bumbl_file):
        print(f"Bumbl file {args.bumbl_file} not found", file=sys.stderr)
        sys.exit(1)

    if args.in_place and args.output is not None:
        print("Cannot use both --in-place and --output.", file=sys.stderr)
        sys.exit(1)

    if args.output is None and not args.in_place:
        base, ext = os.path.splitext(args.bumbl_file)
        args.output = f"{base}.sorted{ext or '.bumbl'}"

    return args


def main(args=None):
    args = parse_arguments(args)

    if args.verify:
        ok = shredutils.verify_bumbl_sorted_column(args.bumbl_file, args.seq_idx)
        if not ok:
            print(
                f"Bumbl file {args.bumbl_file} is not sorted by position in sequence {args.seq_idx}.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Bumbl file {args.bumbl_file} is sorted!", file=sys.stderr)
        sys.exit(0)

    mums = MUMdata(args.bumbl_file, sort=False)
    if mums.num_mums == 0 or mums.num_seqs == 0:
        print("Bumbl file empty.", file=sys.stderr)
        sys.exit(0)

    if not (0 <= args.seq_idx < mums.num_seqs):
        print(
            f"seq_idx {args.seq_idx} invalid for num_seqs={mums.num_seqs}",
            file=sys.stderr,
        )
        sys.exit(1)

    mums.sort(ref_col=args.seq_idx, copy=False)
    mums.blocks = None

    dest = args.bumbl_file if args.in_place else args.output
    out_dir = os.path.dirname(os.path.abspath(dest)) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=".bumbl.sort.", suffix=".tmp", dir=out_dir)
    os.close(fd)
    try:
        mums.write_bums(tmp_path, blocks=None)
        os.replace(tmp_path, dest)
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise

    print(
        f"Wrote sorted bumbl: {dest} (n_mums={mums.num_mums}, seq_idx={args.seq_idx})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
