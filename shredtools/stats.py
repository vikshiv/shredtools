#!/usr/bin/env python3
"""Print metadata and associated sidecar files for a MUM / BUMBL file."""

from __future__ import annotations

import argparse
import json
import os
import sys

from shredtools.utils import MumHeaderInfo, read_mum_header


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Report header metadata and associated files for a .mums or .bumbl file."
    )
    parser.add_argument("mum_file", type=str, help="Path to input .mums or .bumbl file.")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a single JSON object instead of human-readable text.",
    )
    parsed = parser.parse_args(args)

    if not os.path.exists(parsed.mum_file):
        print(f"MUM file {parsed.mum_file} not found", file=sys.stderr)
        sys.exit(1)

    return parsed


def _format_from_path(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".bumbl", ".mums"):
        return ext.lstrip(".")
    return ext.lstrip(".") or "unknown"


def _probe_sidecars(mum_file: str) -> dict[str, dict]:
    prefix = os.path.splitext(mum_file)[0]
    candidates = {
        "lengths": f"{prefix}.lengths",
        "athresh": f"{prefix}.athresh",
        "thresh": f"{prefix}.thresh",
        "filelist": f"{prefix}.filelist",
        "bumbl_index": f"{mum_file}.bi",
        "sorted_bumbl": f"{prefix}.sorted.bumbl",
        "coll_bumbl": f"{prefix}.coll.bumbl",
        "enhanced_bumbl": f"{prefix}.enhanced.bumbl",
    }
    out = {}
    for name, path in candidates.items():
        out[name] = {
            "path": os.path.abspath(path),
            "present": os.path.exists(path),
        }
    return out


def _build_report(mum_file: str, header: MumHeaderInfo) -> dict:
    return {
        "file": header.path,
        "format": _format_from_path(mum_file),
        "num_mums": header.n_mums,
        "num_seqs": header.n_seqs,
        "flags": header.flags,
        "collinear_blocks_stored": header.blocks_stored,
        "num_collinear_blocks": header.num_blocks,
        "associated_files": _probe_sidecars(mum_file),
    }


def _print_human(report: dict) -> None:
    print(f"file: {report['file']}")
    print(f"format: {report['format']}")
    print(f"num_mums: {report['num_mums']}")
    print(f"num_seqs: {report['num_seqs']}")
    blocks = "yes" if report["collinear_blocks_stored"] else "no"
    print(f"collinear_blocks_stored: {blocks}")
    if report["collinear_blocks_stored"]:
        print(f"num_collinear_blocks: {report['num_collinear_blocks']}")
    print("associated_files:")
    for name, info in report["associated_files"].items():
        status = "present" if info["present"] else "missing"
        print(f"  {name}: {info['path']} ({status})")


def main(args=None):
    parsed = parse_arguments(args)
    header = read_mum_header(parsed.mum_file)
    report = _build_report(parsed.mum_file, header)

    if parsed.json:
        print(json.dumps(report, indent=2))
    else:
        _print_human(report)


if __name__ == "__main__":
    main()
