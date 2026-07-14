#!/usr/bin/env python3
"""Unified entrypoint for shredtools subcommands."""

from __future__ import annotations

import importlib
import os
import signal
import sys

# When stdout is closed early (e.g. `| less` then quit), exit cleanly instead
# of raising BrokenPipeError. Python ignores SIGPIPE by default, which turns
# closed-pipe writes into exceptions; restoring SIG_DFL matches typical Unix CLI
# behavior (no traceback).
try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

COMMANDS: dict[str, tuple[str, str, str]] = {
    "enhance": (
        "shredtools.enhance",
        "main",
        "Enhance multi-MUM collection by finding local MUMs in gaps between collinear global MUMs.",
    ),
    "extract": (
        "shredtools.extract_from_mums",
        "main",
        "Extract syntenic regions given a query interval from any assembly in the collection using multi-MUMs.",
    ),
    "fasta": (
        "shredtools.extract_fastas",
        "main",
        "Extract FASTA regions from a BED file produced by shredtools extract.",
    ),
    "filter": (
        "shredtools.filter_collinear",
        "main",
        "Filter out non-collinear MUMs.",
    ),
    "index": (
        "shredtools.index_bumbl",
        "main",
        "Write a .bumbl.bi index for a .bumbl file (single or multi), optionally verifying row order.",
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
        "Subset to only multi-MUM rows contained in a query region from a specific assembly.",
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


def _exit_on_broken_pipe() -> None:
    """Avoid flush errors during interpreter shutdown after stdout closes."""
    try:
        sys.stdout.flush()
    except BrokenPipeError:
        pass
    try:
        fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(fd, sys.stdout.fileno())
    except OSError:
        pass
    raise SystemExit(0)


def main(argv: list[str] | None = None) -> None:
    try:
        _run(argv)
    except BrokenPipeError:
        _exit_on_broken_pipe()


def _run(argv: list[str] | None = None) -> None:
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
