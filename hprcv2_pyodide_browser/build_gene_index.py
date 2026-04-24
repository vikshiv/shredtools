#!/usr/bin/env python3
"""
Build a compact gene/annotation lookup JSON for the Pyodide browser UI.

Supports inputs:
- BED3+ (contig, start, end, [name])
- GTF (gene_id/gene_name in attributes; feature "gene" preferred)
- GFF3 (ID/Name in attributes; feature "gene" preferred)

Output JSON shape (default):
{
  "<build>": {
    "<gene>": [{"contig": "...", "start": 0, "end": 123, "label": "<gene>"}],
    ...
  }
}

Coordinates:
- GTF/GFF3 are 1-based inclusive on disk.
- By default we emit 0-based, half-open intervals (BED-style) for easier use in apps.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


@dataclass(frozen=True, slots=True)
class Interval:
    contig: str
    start: int
    end: int
    label: str


def _open_text(path: str) -> io.TextIOBase:
    if path == "-":
        return io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8", newline="")
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _infer_format(path: str, fmt: Optional[str]) -> str:
    if fmt:
        return fmt.lower()
    p = path.lower()
    if p.endswith((".bed", ".bed.gz")):
        return "bed"
    if p.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    if p.endswith((".gff", ".gff.gz", ".gff3", ".gff3.gz")):
        return "gff3"
    raise SystemExit("Could not infer format; pass --format {bed,gtf,gff3}.")


def _parse_gtf_attrs(s: str) -> Dict[str, str]:
    # Example: gene_id "ENSG..."; gene_name "BRCA1"; transcript_id "..."
    out: Dict[str, str] = {}
    for part in s.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        k, v = part.split(" ", 1)
        v = v.strip().strip('"')
        out[k] = v
    return out


def _parse_gff3_attrs(s: str) -> Dict[str, str]:
    # Example: ID=gene:ENSG...;Name=BRCA1;...
    out: Dict[str, str] = {}
    for part in s.strip().split(";"):
        if not part:
            continue
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def _iter_bed(fin: io.TextIOBase) -> Iterator[Tuple[str, int, int, str]]:
    for line in fin:
        if not line or line.startswith(("#", "track", "browser")):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        contig = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            continue
        name = parts[3] if len(parts) >= 4 and parts[3] else f"{contig}:{start}-{end}"
        if start < 0 or end < 0 or end <= start:
            continue
        yield contig, start, end, name


def _iter_gtf_or_gff(
    fin: io.TextIOBase, fmt: str, feature_preference: Tuple[str, ...]
) -> Iterator[Tuple[str, int, int, str, str]]:
    # Yields: contig, start1, end1, feature, name
    # start/end are still 1-based inclusive here.
    for line in fin:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        contig, _source, feature, start_s, end_s, _score, _strand, _frame, attrs_s = parts[:9]
        try:
            start1 = int(start_s)
            end1 = int(end_s)
        except ValueError:
            continue
        if start1 <= 0 or end1 <= 0 or end1 < start1:
            continue

        if fmt == "gtf":
            attrs = _parse_gtf_attrs(attrs_s)
            name = attrs.get("gene_name") or attrs.get("gene_id") or ""
        else:
            attrs = _parse_gff3_attrs(attrs_s)
            name = attrs.get("Name") or attrs.get("gene_name") or attrs.get("ID") or ""

        if not name:
            continue
        yield contig, start1, end1, feature, name


def _pick_best_features(
    records: Iterable[Tuple[str, int, int, str, str]],
    feature_preference: Tuple[str, ...],
) -> Iterable[Tuple[str, int, int, str]]:
    """
    Prefer exact feature types (e.g. "gene") if present; otherwise fall back to whatever exists.
    Returns: contig, start1, end1, name
    """
    # Two-pass would be expensive; do one pass collecting counts of preferred features.
    kept: List[Tuple[str, int, int, str]] = []
    saw_preferred = False
    preferred = set(feature_preference)
    for contig, s, e, feature, name in records:
        is_pref = feature in preferred
        if is_pref:
            saw_preferred = True
        kept.append((contig, s, e, name, feature))

    if not kept:
        return []
    if not saw_preferred:
        return [(c, s, e, n) for (c, s, e, n, _f) in kept]
    return [(c, s, e, n) for (c, s, e, n, f) in kept if f in preferred]


def _aggregate_intervals(
    it: Iterable[Tuple[str, int, int, str]],
    *,
    input_coords: str,
    output_coords: str,
) -> Dict[str, Dict[str, Tuple[int, int]]]:
    """
    Returns nested dict: gene -> contig -> (min_start, max_end) in output_coords system.
    """
    if input_coords not in ("bed", "gtf"):
        raise ValueError("input_coords must be bed or gtf")
    if output_coords not in ("bed", "gtf"):
        raise ValueError("output_coords must be bed or gtf")

    def to_out(s: int, e: int) -> Tuple[int, int]:
        # input is either:
        # - bed: 0-based, half-open [s,e)
        # - gtf: 1-based inclusive [s,e]
        if input_coords == "gtf":
            # convert to bed first
            s0 = s - 1
            e0 = e
        else:
            s0, e0 = s, e
        if output_coords == "bed":
            return s0, e0
        # output gtf: 1-based inclusive
        return s0 + 1, e0

    agg: Dict[str, Dict[str, Tuple[int, int]]] = {}
    for contig, s, e, name in it:
        s2, e2 = to_out(s, e)
        g = agg.setdefault(name, {})
        if contig not in g:
            g[contig] = (s2, e2)
        else:
            s0, e0 = g[contig]
            g[contig] = (min(s0, s2), max(e0, e2))
    return agg


def _emit_json(
    agg: Dict[str, Dict[str, Tuple[int, int]]],
    *,
    build: str,
    out_path: str,
    label_mode: str,
) -> None:
    out: Dict[str, Dict[str, List[Dict[str, object]]]] = {build: {}}
    for gene, contigs in agg.items():
        hits = []
        for contig, (start, end) in contigs.items():
            if label_mode == "gene":
                label = gene
            elif label_mode == "gene_contig":
                label = f"{gene} {contig}:{start}-{end}"
            else:
                raise ValueError("invalid label_mode")
            hits.append({"contig": contig, "start": int(start), "end": int(end), "label": label})
        # stable output
        hits.sort(key=lambda h: (str(h["contig"]), int(h["start"]), int(h["end"])))
        out[build][gene] = hits

    if out_path == "-":
        json.dump(out, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return
    tmp = out_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, out_path)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build gene/annotation JSON index for the Pyodide browser."
    )
    ap.add_argument("input", help="Input annotation file (BED/GTF/GFF3). Use '-' for stdin.")
    ap.add_argument(
        "--format",
        choices=["bed", "gtf", "gff3"],
        default=None,
        help="Input format (auto-infer from filename if omitted).",
    )
    ap.add_argument(
        "--build",
        default="GRCh38",
        help="Build label key in output JSON (e.g. GRCh38, CHM13v2.0).",
    )
    ap.add_argument(
        "--output",
        "-o",
        default="annotations.genes.json",
        help="Output JSON path (default: annotations.genes.json). Use '-' for stdout.",
    )
    ap.add_argument(
        "--output-coords",
        choices=["bed", "gtf"],
        default="bed",
        help="Emit coordinates as 'bed' (0-based half-open) or 'gtf' (1-based inclusive).",
    )
    ap.add_argument(
        "--label",
        choices=["gene", "gene_contig"],
        default="gene",
        help="Label field to store per hit.",
    )
    ap.add_argument(
        "--prefer-features",
        default="gene",
        help="Comma-separated feature types to prefer (default: gene).",
    )
    args = ap.parse_args(argv)

    fmt = _infer_format(args.input, args.format)
    prefer = tuple([x.strip() for x in str(args.prefer_features).split(",") if x.strip()])
    if not prefer:
        prefer = ("gene",)

    with _open_text(args.input) as fin:
        if fmt == "bed":
            records = _iter_bed(fin)
            agg = _aggregate_intervals(records, input_coords="bed", output_coords=args.output_coords)
        else:
            raw = _iter_gtf_or_gff(fin, fmt=fmt, feature_preference=prefer)
            picked = _pick_best_features(raw, feature_preference=prefer)
            # picked yields (contig, start1, end1, name) with gtf-style coords
            agg = _aggregate_intervals(picked, input_coords="gtf", output_coords=args.output_coords)

    _emit_json(agg, build=str(args.build), out_path=str(args.output), label_mode=str(args.label))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

