#!/usr/bin/env python3
"""Unified entrypoint for shredtools subcommands."""

from __future__ import annotations

import importlib
import sys

COMMANDS: dict[str, tuple[str, str, str]] = {
    "enhance": (
        "shredtools.enhance",
        "main",
        "Enhance multi-MUM collection by finding local MUMs in gaps between collinear global MUMs.",
    ),
    "extract": (
        "shredtools.extract_from_mums",
        "main",
        "Extract syntenic regions from a query interval (BED).",
    ),
    "fasta": (
        "shredtools.extract_fastas",
        "main",
        "Extract FASTA regions from a BED file produced by extract_from_mums.",
    ),
    "filter": (
        "shredtools.filter_collinear",
        "main",
        "Filter out non-collinear MUMs.",
    ),
    "index": (
        "shredtools.index_bumbl",
        "main",
        "Verify bumbl row order and write a .bumbl.bi index (single or multi).",
    ),
    "shred": (
        "shredtools.shred_from_mums",
        "main",
        "Shred a MUM file into smaller fragments and optionally visualize.",
    ),
    "sort": (
        "shredtools.sort_bumbl",
        "main",
        "Sort a .bumbl file by start position and write output.",
    ),
    "stats": (
        "shredtools.stats",
        "main",
        "Report header metadata and associated files for a .mums or .bumbl file.",
    ),
    "subset": (
        "shredtools.subset_from_mums",
        "main",
        "Subset multi-MUM rows overlapping a query region on one assembly.",
    ),
}

__version__ = "0.1.0"

_COMMAND_HELP = "\n".join(
    f"  {name:8}  {desc}" for name, (_, _, desc) in sorted(COMMANDS.items())
)


def _print_top_level_help() -> None:
    print("usage: shredtools <command> [options]")
    print()
    print("Manipulate multi-MUMs and pangenome collections.")
    print()
    print("commands:")
    print(_COMMAND_HELP)
    print()
    print("Run `shredtools <command> -h` for command-specific help.")
    print()
    print("Version:", __version__)


def main(argv: list[str] | None = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    if not argv or argv[0] in ("-h", "--help"):
        _print_top_level_help()
        return

    if argv[0] in ("--version", "-V", "version"):
        print(f"shredtools version {__version__}")
        return

    cmd = argv[0]
    rest = argv[1:]

    if cmd in ("-h", "--help"):
        _print_top_level_help()
        return

    spec = COMMANDS.get(cmd)
    if spec is None:
        print(f"Unknown command: {cmd}", file=sys.stderr)
        _print_top_level_help()
        raise SystemExit(1)

    module_name, attr, _ = spec
    module = importlib.import_module(module_name)
    entry = getattr(module, attr)

    if cmd == "fasta":
        entry(argv=rest)
    else:
        entry(rest)


if __name__ == "__main__":
    main()
