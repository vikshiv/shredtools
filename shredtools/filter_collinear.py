#!/usr/bin/env python3

"""
Filter MUMs that are not part of collinear blocks (for indexing with shredtools)
"""

import argparse
import os
import sys
import tempfile

import mumemto.utils as mutils

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import utils as shredutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Filter out non-collinear MUMs"
    )
    parser.add_argument("bumbl_file", type=str, help="Path to the input .bumbl file.")
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output .bumbl path (default: <input>.coll.bumbl).",
    )
    parser.add_argument(
        "--in-place",
        "-i",
        action="store_true",
        help="Rewrite the input file after writing to a temporary file.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose mode",
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

    mums = mutils.MUMdata(args.bumbl_file, sort=True, verbose=args.verbose)
    if mums.num_mums == 0 or mums.num_seqs == 0:
        print("Bumbl file empty.", file=sys.stderr)
        sys.exit(0)

    if mums.blocks is None:
        mums.blocks = mutils.find_coll_blocks(mums, verbose=args.verbose)
        print(f"found blocks: {len(mums.blocks)}", file=sys.stderr)
    else:
        print(f"using provided blocks: {len(mums.blocks)}", file=sys.stderr)

    coll_mums = [i for s, e in mums.blocks for i in range(s, e+1)]
    print(f"Collinear MUMs: {len(coll_mums)} / {len(mums)} ({len(coll_mums) * 100 / len(mums):.3f}%)", file=sys.stderr)
    
    new_blocks = []
    last = -1
    for s, e in mums.blocks:
        cur = e - s + 1
        new_blocks.append((last+1, last+1+cur))
        last = last+1+cur

    mums = mums[coll_mums, :]

    dest = args.bumbl_file if args.in_place else args.output
    out_dir = os.path.dirname(os.path.abspath(dest)) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=".bumbl.coll.", suffix=".tmp", dir=out_dir)
    os.close(fd)
    try:
        mums.write_bums(tmp_path, blocks=new_blocks)
        os.replace(tmp_path, dest)
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise

    print(
        f"Wrote filtered bumbl: {dest} (n_mums={mums.num_mums}, num_blocks={len(new_blocks)})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
