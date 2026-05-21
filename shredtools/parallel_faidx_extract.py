#!/usr/bin/env python3
"""Extract FASTA regions in parallel with ``samtools faidx``.

Input lines (whitespace-separated): contig  start  end  path/to/reference.fa

Each reference must have a ``.fai`` index. All extracts go into the output directory
``--out-dir``, one FASTA per line, named ``<fasta_basename_stem>_extract<suffix>`` (for
example ``100000_CHM13.pri.fa`` → ``100000_CHM13.pri_extract.fa``). If several lines use
the same FASTA path, a short ``_<contig>_<start>_<end>`` fragment is inserted before the
extension so names stay unique. Use ``-v`` for a tqdm progress bar; successful runs do not print each path.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import NamedTuple, Sequence

from tqdm.auto import tqdm


class RegionLine(NamedTuple):
    contig: str
    start: int
    end: int
    fasta: str


def parse_regions_file(path: str) -> list[RegionLine]:
    rows: list[RegionLine] = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 3)
            if len(parts) < 4:
                print(f"{path}:{lineno}: need at least 4 fields (contig start end fasta)", file=sys.stderr)
                raise SystemExit(1)
            contig, s, e, fasta = parts[0], parts[1], parts[2], parts[3]
            try:
                start, end = int(s), int(e)
            except ValueError as ex:
                print(f"{path}:{lineno}: invalid start/end: {ex}", file=sys.stderr)
                raise SystemExit(1) from ex
            rows.append(RegionLine(contig=contig, start=start, end=end, fasta=fasta))
    if not rows:
        print(f"{path}: no data rows", file=sys.stderr)
        raise SystemExit(1)
    return rows


def _sanitize_filename(s: str, max_len: int = 200) -> str:
    out = re.sub(r"[^\w.\-]+", "_", s, flags=re.ASCII)
    out = out.strip("._") or "region"
    return out[:max_len]


def region_to_samtools(contig: str, start: int, end: int, bed: bool) -> str:
    """Return region string for samtools faidx (1-based inclusive)."""
    if bed:
        # BED-style: 0-based start, half-open end
        lo, hi = start + 1, end
    else:
        lo, hi = start, end
    if lo < 1 or hi < lo:
        raise ValueError(f"invalid interval after conversion: {contig}:{lo}-{hi}")
    return f"{contig}:{lo}-{hi}"


def assign_output_filenames(rows: list[RegionLine]) -> list[str]:
    """One output basename per row: ``{stem}_extract{ext}``, disambiguated if needed."""
    base_names: list[str] = []
    for r in rows:
        base = os.path.basename(r.fasta)
        root, ext = os.path.splitext(base)
        if not ext:
            ext = ".fa"
        base_names.append(f"{root}_extract{ext}")

    dup_groups: defaultdict[str, list[int]] = defaultdict(list)
    for i, bn in enumerate(base_names):
        dup_groups[bn].append(i)

    out: list[str] = list(base_names)
    for bn, idxs in dup_groups.items():
        if len(idxs) <= 1:
            continue
        for i in idxs:
            r = rows[i]
            root, ext = os.path.splitext(os.path.basename(r.fasta))
            if not ext:
                ext = ".fa"
            slug = _sanitize_filename(f"{r.contig}_{r.start}_{r.end}")
            out[i] = f"{root}_extract_{slug}{ext}"
    return out


def _one_job(args: tuple[RegionLine, str, str, str, bool]) -> tuple[str, str | None]:
    row, out_dir, out_basename, samtools, bed = args
    fai = row.fasta + ".fai"
    if not os.path.isfile(row.fasta):
        return (row.fasta, f"missing FASTA: {row.fasta}")
    if not os.path.isfile(fai):
        return (row.fasta, f"missing FAI: {fai}")

    os.makedirs(out_dir, exist_ok=True)

    try:
        region = region_to_samtools(row.contig, row.start, row.end, bed)
    except ValueError as e:
        return (row.fasta, str(e))

    out_fa = os.path.join(out_dir, out_basename)

    cmd = [samtools, "faidx", row.fasta, region]
    try:
        with open(out_fa, "wb") as out:
            subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        err = e.stderr.decode("utf-8", errors="replace") if e.stderr else str(e)
        return (row.fasta, f"{' '.join(cmd)} failed: {err.strip()}")

    return (out_fa, None)


def main(argv: Sequence[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("regions_file", help="Whitespace-separated: contig start end path/to.fa")
    p.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        required=True,
        help="Output directory; one *_extract*.fa per input line is written here",
    )
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1, help="Parallel workers (default: CPU count)")
    p.add_argument("--samtools", default="samtools", help="samtools executable (default: samtools)")
    p.add_argument(
        "--bed",
        action="store_true",
        help="Interpret start/end as BED (0-based start, exclusive end); convert for samtools",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Show a tqdm progress bar on stderr")
    args = p.parse_args(list(argv) if argv is not None else None)

    rows = parse_regions_file(args.regions_file)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_names = assign_output_filenames(rows)

    work = [(r, out_dir, out_names[i], args.samtools, args.bed) for i, r in enumerate(rows)]
    failures: list[tuple[str, str]] = []

    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futures = {ex.submit(_one_job, w): w for w in work}
        done_iter = as_completed(futures)
        if args.verbose:
            done_iter = tqdm(done_iter, total=len(futures), desc="faidx", unit="region")
        for fut in done_iter:
            path_or_fa, err = fut.result()
            if err:
                failures.append((path_or_fa, err))

    if failures:
        print("\nErrors:", file=sys.stderr)
        for target, msg in failures:
            print(f"  {target}: {msg}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
