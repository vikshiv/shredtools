#!/usr/bin/env python3

"""
Filter a .bumbl file:

- default: keep only MUMs in collinear blocks (coordinate-system MUMs)
- --partial: drop partial MUMs (any start == -1)
"""

import argparse
import os
import sys
import tempfile

import mumemto.utils as mutils


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description=(
            "Filter a .bumbl file. Default: keep only collinear-block MUMs. "
            "With --partial: drop partial MUMs (starts containing -1)."
        )
    )
    parser.add_argument("bumbl_file", type=str, help="Path to the input .bumbl file.")
    parser.add_argument(
        "--partial",
        action="store_true",
        help="Drop partial MUMs (any assembly start == -1). Skips collinear-block filtering.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help=(
            "Output .bumbl path (default: <input>.coll.bumbl, or "
            "<input>.strict.bumbl with --partial)."
        ),
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
        suffix = ".strict" if args.partial else ".coll"
        args.output = f"{base}{suffix}{ext or '.bumbl'}"

    return args


def filter_partial_mums(mums):
    """Drop MUMs with any absent start (-1). Mutates and returns mums."""
    n_before = len(mums)
    mums.filter_pmums()
    n_after = len(mums)
    removed = n_before - n_after
    pct = (n_after * 100 / n_before) if n_before else 0.0
    print(
        f"Strict MUMs: {n_after} / {n_before} ({pct:.3f}%; removed {removed} partial)",
        file=sys.stderr,
    )
    mums.blocks = None
    return mums


def filter_collinear_mums(mums, verbose=False):
    """Keep only MUMs in collinear blocks. Returns (mums, new_blocks)."""
    if mums.blocks is None:
        mums.blocks = mutils.find_coll_blocks(mums, verbose=verbose)
        print(f"found blocks: {len(mums.blocks)}", file=sys.stderr)
    else:
        print(f"using provided blocks: {len(mums.blocks)}", file=sys.stderr)

    coll_mums = [i for s, e in mums.blocks for i in range(s, e + 1)]
    print(
        f"Collinear MUMs: {len(coll_mums)} / {len(mums)} "
        f"({len(coll_mums) * 100 / len(mums):.3f}%)",
        file=sys.stderr,
    )

    new_blocks = []
    last = -1
    for s, e in mums.blocks:
        cur = e - s + 1
        new_blocks.append((last + 1, last + cur))
        last = last + cur

    mums = mums[coll_mums, :]
    return mums, new_blocks


def _write_bumbl(mums, dest, blocks=None):
    out_dir = os.path.dirname(os.path.abspath(dest)) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=".bumbl.filter.", suffix=".tmp", dir=out_dir)
    os.close(fd)
    try:
        mums.write_bums(tmp_path, blocks=blocks)
        os.replace(tmp_path, dest)
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def main(args=None):
    args = parse_arguments(args)

    # Collinear mode sorts; --partial only needs a row scan for -1 starts.
    mums = mutils.MUMdata(
        args.bumbl_file,
        sort=not args.partial,
        verbose=args.verbose,
    )
    if mums.num_mums == 0 or mums.num_seqs == 0:
        print("Bumbl file empty.", file=sys.stderr)
        sys.exit(0)

    dest = args.bumbl_file if args.in_place else args.output
    blocks = None

    if args.partial:
        mums = filter_partial_mums(mums)
    else:
        mums, blocks = filter_collinear_mums(mums, verbose=args.verbose)

    _write_bumbl(mums, dest, blocks=blocks)

    if args.partial:
        print(
            f"Wrote strict bumbl: {dest} (n_mums={mums.num_mums})",
            file=sys.stderr,
        )
    else:
        print(
            f"Wrote filtered bumbl: {dest} "
            f"(n_mums={mums.num_mums}, num_blocks={len(blocks)})",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
