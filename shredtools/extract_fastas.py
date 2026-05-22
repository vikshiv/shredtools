#!/usr/bin/env python3
"""Extract FASTA regions from a BED file produced by ``extract_from_mums``.

Input is tab-separated BED (4 columns): contig, start, end, path/to/reference.fa.
Coordinates are 0-based start and exclusive end (standard BED).

Without ``--agc``, each row is extracted with ``pysam`` (requires ``.fai`` next to each
reference FASTA). With ``--agc``, all regions are fetched in one ``agc getctg`` call and
optionally split into per-row FASTAs (see ``--no-split``).
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

import pysam
from tqdm.auto import tqdm


class RegionLine(NamedTuple):
    contig: str
    start: int
    end: int
    fasta: str


def parse_bed_file(path: str) -> list[RegionLine]:
    rows: list[RegionLine] = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                print(f"{path}:{lineno}: need 4 tab-separated fields (contig start end fasta)", file=sys.stderr)
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


def region_to_interval(contig: str, start: int, end: int) -> str:
    """Return 1-based inclusive region string for agc getctg."""
    lo, hi = start + 1, end
    if lo < 1 or hi < lo:
        raise ValueError(f"invalid interval after conversion: {contig}:{lo}-{hi}")
    return f"{contig}:{lo}-{hi}"


def fasta_header(contig: str, start: int, end: int) -> str:
    """FASTA header matching samtools faidx style (1-based inclusive coords)."""
    lo, hi = start + 1, end
    return f">{contig}:{lo}-{hi}"


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


def validate_fai_paths(rows: list[RegionLine]) -> list[str]:
    """Return error messages for missing FASTA or FAI files."""
    errors: list[str] = []
    seen: set[str] = set()
    for r in rows:
        if r.fasta in seen:
            continue
        seen.add(r.fasta)
        if not os.path.isfile(r.fasta):
            errors.append(f"missing FASTA: {r.fasta}")
        elif not os.path.isfile(r.fasta + ".fai"):
            errors.append(f"missing FAI: {r.fasta}.fai")
    return errors


def list_agc_contigs(agc: str, archive: str) -> set[str]:
    result = subprocess.run(
        [agc, "listctg", archive],
        check=True,
        capture_output=True,
        text=True,
    )
    contigs: set[str] = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            contigs.add(parts[-1])
        elif len(parts) == 1:
            contigs.add(parts[0])
    return contigs


def validate_agc_contigs(rows: list[RegionLine], agc_contigs: set[str]) -> list[str]:
    missing = sorted({r.contig for r in rows if r.contig not in agc_contigs})
    if not missing:
        return []
    return [f"contig not in AGC archive: {c}" for c in missing]


def split_multifasta(combined_path: str, out_names: list[str], out_dir: str) -> None:
    """Split a multi-record FASTA into one file per output name."""
    if len(out_names) == 0:
        return
    idx = 0
    out_fh = None
    with open(combined_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if out_fh is not None:
                    out_fh.close()
                if idx >= len(out_names):
                    print(
                        f"{combined_path}: more FASTA records than BED rows ({len(out_names)} expected)",
                        file=sys.stderr,
                    )
                    raise SystemExit(1)
                out_path = os.path.join(out_dir, out_names[idx])
                out_fh = open(out_path, "w", encoding="utf-8")
                idx += 1
            if out_fh is not None:
                out_fh.write(line)
    if out_fh is not None:
        out_fh.close()
    if idx != len(out_names):
        print(
            f"{combined_path}: expected {len(out_names)} FASTA records, got {idx}",
            file=sys.stderr,
        )
        raise SystemExit(1)


def _extract_one_fasta(
    args: tuple[str, str, list[tuple[RegionLine, str]]],
) -> list[tuple[str, str | None]]:
    """Extract all BED rows that share one reference FASTA (one open handle per file)."""
    fasta, out_dir, items = args
    results: list[tuple[str, str | None]] = []
    try:
        ff = pysam.FastaFile(fasta)
    except (OSError, ValueError) as e:
        return [(fasta, str(e))] * len(items)

    try:
        for row, out_basename in items:
            out_fa = os.path.join(out_dir, out_basename)
            try:
                if row.start < 0 or row.end <= row.start:
                    raise ValueError(
                        f"invalid BED interval: {row.contig}:{row.start}-{row.end}"
                    )
                seq = ff.fetch(row.contig, row.start, row.end)
                with open(out_fa, "w", encoding="utf-8") as out:
                    out.write(f"{fasta_header(row.contig, row.start, row.end)}\n{seq}\n")
                results.append((out_fa, None))
            except (KeyError, ValueError) as e:
                results.append((row.fasta, f"{row.contig}:{row.start}-{row.end}: {e}"))
    finally:
        ff.close()

    return results


def run_pysam_extract(
    rows: list[RegionLine],
    out_dir: str,
    out_names: list[str],
    threads: int,
    verbose: bool,
) -> None:
    fai_errors = validate_fai_paths(rows)
    if fai_errors:
        print("File errors:", file=sys.stderr)
        for msg in fai_errors:
            print(f"  {msg}", file=sys.stderr)
        raise SystemExit(1)

    by_fasta: defaultdict[str, list[tuple[RegionLine, str]]] = defaultdict(list)
    for i, r in enumerate(rows):
        by_fasta[r.fasta].append((r, out_names[i]))

    work = [(fasta, out_dir, items) for fasta, items in by_fasta.items()]
    failures: list[tuple[str, str]] = []

    with ProcessPoolExecutor(max_workers=max(1, threads)) as ex:
        futures = {ex.submit(_extract_one_fasta, w): w for w in work}
        done_iter = as_completed(futures)
        if verbose:
            done_iter = tqdm(done_iter, total=len(futures), desc="extract", unit="file")
        for fut in done_iter:
            for path_or_fa, err in fut.result():
                if err:
                    failures.append((path_or_fa, err))

    if failures:
        print("\nErrors:", file=sys.stderr)
        for target, msg in failures:
            print(f"  {target}: {msg}", file=sys.stderr)
        raise SystemExit(1)


def run_agc_extract(
    rows: list[RegionLine],
    out_dir: str,
    out_names: list[str],
    archive: str,
    agc: str,
    threads: int,
    no_split: bool,
) -> None:
    try:
        agc_contigs = list_agc_contigs(agc, archive)
    except subprocess.CalledProcessError as e:
        err = e.stderr.decode("utf-8", errors="replace") if e.stderr else str(e)
        print(f"agc listctg failed: {err.strip()}", file=sys.stderr)
        raise SystemExit(1) from e

    contig_errors = validate_agc_contigs(rows, agc_contigs)
    if contig_errors:
        print("Preflight errors:", file=sys.stderr)
        for msg in contig_errors:
            print(f"  {msg}", file=sys.stderr)
        raise SystemExit(1)

    regions: list[str] = []
    for r in rows:
        try:
            regions.append(region_to_interval(r.contig, r.start, r.end))
        except ValueError as e:
            print(f"invalid region {r.contig}:{r.start}-{r.end}: {e}", file=sys.stderr)
            raise SystemExit(1) from e

    combined_name = "combined_extract.fa"
    combined_path = os.path.join(out_dir, combined_name)
    cmd = [agc, "getctg", "-t", str(threads), "-o", combined_path, archive, *regions]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        err = e.stderr or e.stdout or str(e)
        print(f"{' '.join(cmd[:6])} ... failed: {err.strip()}", file=sys.stderr)
        raise SystemExit(1) from e

    if not no_split:
        split_multifasta(combined_path, out_names, out_dir)
        os.remove(combined_path)


def main(argv: Sequence[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bed_file", help="BED from extract_from_mums (contig, start, end, fasta path)")
    p.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        required=True,
        help="Output directory; one *_extract*.fa per BED row (unless --no-split with --agc)",
    )
    p.add_argument("--agc", metavar="PATH", help="AGC archive; extract via agc getctg instead of pysam")
    p.add_argument(
        "-t",
        "--threads",
        type=int,
        default=os.cpu_count() or 1,
        help="Threads: agc getctg -t, or parallel reference FASTAs via pysam (default: CPU count)",
    )
    p.add_argument(
        "--no-split",
        action="store_true",
        help="With --agc: keep single combined FASTA instead of splitting per BED row",
    )
    p.add_argument("--agc-bin", default="agc", help="agc executable (default: agc)")
    p.add_argument("-v", "--verbose", action="store_true", help="Show a tqdm progress bar on stderr")
    args = p.parse_args(list(argv) if argv is not None else None)

    if args.no_split and not args.agc:
        print("--no-split requires --agc", file=sys.stderr)
        raise SystemExit(1)

    rows = parse_bed_file(args.bed_file)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_names = assign_output_filenames(rows)
    threads = max(1, args.threads)

    if args.agc:
        if not os.path.isfile(args.agc):
            print(f"AGC archive not found: {args.agc}", file=sys.stderr)
            raise SystemExit(1)
        run_agc_extract(
            rows, out_dir, out_names, args.agc, args.agc_bin, threads, args.no_split
        )
    else:
        run_pysam_extract(rows, out_dir, out_names, threads, args.verbose)


if __name__ == "__main__":
    main()
