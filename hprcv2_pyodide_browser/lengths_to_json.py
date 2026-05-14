#!/usr/bin/env python3
"""
Convert multilengths ``.lengths`` files to JSON for tooling / the Pyodide browser.

**Single file** — one object with ``seq_lengths_multi`` and ``contig_names``::

    conda run -n shredtools python lengths_to_json.py -i path/to/index.lengths -o one.json

**Merged bundle** (e.g. ``pangenome_lengths.json`` for the browser) — top-level keys, each value is that same shape::

    conda run -n shredtools python lengths_to_json.py -o pangenome_lengths.json \\
      --preset hprcv2_enhanced=~/vast/shredtools/index/hprcv2_enhanced_merged.lengths \\
      --preset hprcv2_merged=~/vast/shredtools/index/hprcv2_merged.lengths \\
      --preset hprcv1=~/vast/shredtools/index/hprcv1.lengths

Use ``mumemto.utils`` when available; otherwise uses ``bumbl_index_utils`` in this directory
(same multilengths parsing rules).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _load_mutils():
    try:
        import mumemto.utils as mutils  # type: ignore
    except ImportError:
        here = Path(__file__).resolve().parent
        if str(here) not in sys.path:
            sys.path.insert(0, str(here))
        import bumbl_index_utils as mutils

    return mutils


def lengths_path_to_record(path: str | Path) -> dict:
    """Parse one multilengths file into a JSON-serializable dict."""
    mutils = _load_mutils()
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(p)
    seq_lengths_multi = mutils.get_sequence_lengths(str(p), multilengths=True)
    contig_names = mutils.get_contig_names(str(p))
    return {
        "seq_lengths_multi": [[int(x) for x in row] for row in seq_lengths_multi],
        "contig_names": [[str(c) for c in row] for row in contig_names],
    }


def parse_preset_arg(spec: str) -> tuple[str, Path]:
    """Parse ``KEY=PATH`` from ``--preset``."""
    if "=" not in spec:
        raise ValueError(f"--preset must be KEY=PATH, got {spec!r}")
    key, _, path = spec.partition("=")
    key = key.strip()
    if not key:
        raise ValueError(f"Missing preset key in {spec!r}")
    return key, Path(path.strip())


def build_merged_bundle(presets: dict[str, Path]) -> dict[str, dict]:
    """``{preset_key: lengths_path_to_record(path), ...}``."""
    return {k: lengths_path_to_record(p) for k, p in presets.items()}


def write_json(data: dict, output: Path) -> str:
    """Write compact JSON; returns the serialized string (for size reporting)."""
    output.parent.mkdir(parents=True, exist_ok=True)
    blob = json.dumps(data, separators=(",", ":"))
    output.write_text(blob, encoding="utf-8")
    return blob


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output ``.json`` path.",
    )
    p.add_argument(
        "-i",
        "--input",
        type=Path,
        default=None,
        help="One multilengths file → single-record JSON (mutually exclusive with --preset).",
    )
    p.add_argument(
        "--preset",
        action="append",
        default=None,
        metavar="KEY=PATH",
        help="Repeat for merged bundle JSON (mutually exclusive with -i).",
    )
    args = p.parse_args()

    if args.preset and args.input is not None:
        p.error("Use either -i (single record) or one or more --preset (merged), not both.")

    if args.preset:
        presets: dict[str, Path] = {}
        for raw in args.preset:
            key, path = parse_preset_arg(raw)
            if key in presets:
                p.error(f"Duplicate preset key {key!r}")
            presets[key] = path
        for key, path in presets.items():
            if not path.is_file():
                raise SystemExit(f"Missing lengths file for preset {key!r}: {path}")
        data = build_merged_bundle(presets)
    else:
        if args.input is None:
            p.error("Provide -i for a single file, or one or more --preset KEY=PATH for merged output.")
        data = lengths_path_to_record(args.input)

    blob = write_json(data, args.output)
    print(f"Wrote {args.output} ({len(blob.encode('utf-8'))} bytes)")


if __name__ == "__main__":
    main()
