#!/usr/bin/env python3
"""Unified entrypoint for shredtools subcommands."""

from __future__ import annotations

import importlib
import sys

COMMANDS: dict[str, tuple[str, str]] = {
    "extract": ("shredtools.extract_from_mums", "main"),
    "fasta": ("shredtools.extract_fastas", "main"),
    "filter": ("shredtools.filter_collinear", "main"),
    "index": ("shredtools.index_bumbl", "main"),
    "shred": ("shredtools.shred_from_mums", "main"),
    "sort": ("shredtools.sort_bumbl", "main"),
    "enhance": ("shredtools.enhance", "main"),
    "stats": ("shredtools.stats", "main"),
}

_COMMAND_HELP = "\n".join(f"  {name:8}  {mod.rsplit('.', 1)[-1]}" for name, (mod, _) in sorted(COMMANDS.items()))


def _print_top_level_help() -> None:
    print("usage: shredtools <command> [options]")
    print()
    print("Manipulate multi-MUMs and pangenome collections.")
    print()
    print("commands:")
    print(_COMMAND_HELP)
    print()
    print("Run `shredtools <command> -h` for command-specific help.")


def main(argv: list[str] | None = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    if not argv or argv[0] in ("-h", "--help"):
        _print_top_level_help()
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

    module_name, attr = spec
    module = importlib.import_module(module_name)
    entry = getattr(module, attr)

    if cmd == "fasta":
        entry(argv=rest)
    else:
        entry(rest)


if __name__ == "__main__":
    main()
