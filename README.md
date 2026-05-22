# **shredtools**: pangenome coordinates from multi-MUMs

Shredtools is a toolkit for querying homologous regions in a [mumemto](https://github.com/vikshiv/mumemto) pangenome. It reads binary `.bumbl` multi-MUM files, builds interval indexes (`.bumbl.bi`) for fast lookup, and extracts corresponding coordinates across the collection as BED (and optionally FASTA).

Shredtools expects a mumemto run that produced a `.bumbl` file (use `mumemto -b` or `mumemto convert`) and the matching `.lengths` sidecar from the same output prefix.

---

## Installation

Requires [mumemto](https://github.com/vikshiv/mumemto) at runtime (not installed automatically by pip). Python 3.9+.

```bash
git clone https://github.com/vikshiv/shredtools.git
cd shredtools
pip install -e .
shredtools -h
```

You can also run `python -m shredtools`.

**Pip dependencies:** `numpy`, `pysam`, `tqdm`

**Optional:** `matplotlib` (for `extract --plot` / `--plot-full`); `.fai` indexes beside reference FASTAs, or the `agc` tool (`fasta --agc`) for sequence extraction.

---

## Quick start

After running mumemto on a collection (e.g. `mumemto assemblies/*.fa -o collection -b`), build an index and extract a region:

```bash
# Prepare .bumbl for indexing (sort required; filter optional)
shredtools sort collection.bumbl -s 0
shredtools filter collection.sorted.bumbl -o collection.coll.bumbl   # optional

# Build .bumbl.bi (default: multi-index)
shredtools index collection.coll.bumbl -s 0 -v

# Extract homologous BED, then fetch FASTA
shredtools extract collection.coll.bumbl -s 0 -r chr1:1000000-2000000 \
  -o regions/prefix -l collection.lengths
shredtools fasta regions/prefix.bed -o fasta_out/
```

---

## Getting started

### Index a `.bumbl` file

`extract` on a `.bumbl` input requires a `.bumbl.bi` index alongside it. The index maps genomic bins to row ranges in the `.bumbl` file so only overlapping MUMs are loaded (including HTTP range requests for remote files).

**Prerequisites**

1. **Sorted rows.** MUM rows must be non-decreasing on the reference column you index (`starts[:, seq_idx]`). Use `shredtools sort`; this drops embedded collinear-block metadata from the file.
2. **Optional filter.** `shredtools filter` keeps only MUMs that belong to collinear blocks, which can yield a cleaner index for synteny-style queries.

```bash
shredtools sort collection.bumbl -s 0
# → collection.sorted.bumbl (default output name)

shredtools filter collection.sorted.bumbl -o collection.coll.bumbl
# optional; use -o explicitly (default output name differs)
```

**Build the index**

```bash
shredtools index collection.coll.bumbl -s 0
# → collection.coll.bumbl.bi (default)
```

| Flag | Default | Description |
|------|---------|-------------|
| `-s` / `--seq-idx` | `0` | Column used for sortedness check; reference axis for single-index |
| `--multi` | on (default) | Multi-index: per-sequence bins (used by `extract`) |
| `--single` | off | Single-index over one reference column |
| `-w` / `--bin-width` | `1000000` | Genomic bin width (bp) |
| `-o` / `--output` | `<bumbl>.bi` | Output path |
| `--no-verify-sorted` | off | Skip non-decreasing check on `starts[:, seq_idx]` |
| `--no-verify-checksum` | off | Skip post-write checksum verification |
| `-v` / `--verbose` | off | Progress on stderr |

**Examples**

```bash
# Default multi-index (recommended for extract)
shredtools index collection.coll.bumbl -s 0 -v

# Finer bins for dense genomes
shredtools index collection.coll.bumbl -s 0 -w 500000 -o collection.coll.bumbl.bi

# Single-reference index only
shredtools index collection.coll.bumbl -s 0 --single -o collection.coll.single.bi
```

> [!TIP]
> Run `shredtools stats collection.bumbl` to list expected sidecar files (`.lengths`, `.bi`, sorted/filtered variants).

> [!NOTE]
> Verify sort order without rewriting: `shredtools sort collection.bumbl -s 0 --verify`

---

### `.bumbl.bi` index (overview)

A `.bumbl.bi` file is a binary sidecar keyed to a specific `.bumbl` **row order**. It stores a small header (magic `bumblbi`, format byte, checksum over the first 1000 MUM lengths, bin width, sequence count) followed by bin structures that map genomic intervals to half-open MUM row ranges `[mum_start, mum_end)`.

- **Multi-index (default, format byte 1):** per-sequence genomic bins. `extract` uses this for any reference column (`-s`).
- **Single-index (format byte 0):** one reference column only; requires rows sorted on that column.

For the full on-disk layout, query algorithms, and diagrams, see **[bumbl_index/bumbl_index.md](bumbl_index/bumbl_index.md)**.

---

### Extract homologous regions

Given a query interval on one pangenome sequence, `extract` finds the bounding multi-MUMs and reports the homologous interval on each selected genome as BED.

**Required inputs**

| File | Resolution |
|------|------------|
| `collection.bumbl` | Positional argument |
| `collection.bumbl.bi` | Default: `<mum_file>.bi`, or `-b` |
| `collection.lengths` | Default: `<stem>.lengths`, or `-l` (mumemto multilengths format) |

**Region format.** `--range` is **contig-local** on the sequence named by `-s`: `contig:start-end` (e.g. `chr1:1000000-2000000`). Shredtools converts this to global linear coordinates using the `.lengths` manifest, then queries the index.

```bash
shredtools extract collection.coll.bumbl \
  -s 0 \
  -r chr1:1000000-2000000 \
  -o regions/prefix \
  -l collection.lengths
```

| Flag | Required | Description |
|------|----------|-------------|
| `-s` / `--seq-idx` | yes | Pangenome column that defines the query coordinate system |
| `-r` / `--range` | yes | `contig:start-end` on that sequence |
| `-o` / `--output` | no | Output prefix → `<prefix>.bed`; writes to stdout if omitted |
| `-l` / `--lengths` | no | Multilengths file (default: `<stem>.lengths`) |
| `-b` / `--bumblbi` | no | Index path (default: `<mum_file>.bi`) |
| `-x` / `--sequences` | no | Subset of sequence indices (default: all) |
| `--plot` | no | Write `<prefix>_extract_synteny.pdf` (requires `-o`, needs matplotlib) |
| `--plot-full` | no | Write `<prefix>_full_synteny.pdf` (requires `-o`) |

**More examples**

```bash
# BED to stdout; only genomes 0, 1, 2
shredtools extract collection.coll.bumbl -s 0 -r chr1:1-50000 -x 0 1 2

# With synteny plot
shredtools extract collection.coll.bumbl -s 0 -r chr1:1000000-2000000 \
  -o regions/prefix -l collection.lengths --plot

# Remote .bumbl and index (HTTP range requests)
shredtools extract https://example.org/collection.bumbl -s 0 -r chr1:1-1000 \
  -b https://example.org/collection.bumbl.bi -l collection.lengths
```

When the query interval does not align exactly to MUM boundaries, `extract` prints `left margin` and/or `right margin` on stderr (distance in bp from the query edge to the enclosing MUM).

**BED output format**

Four tab-separated columns (0-based start, exclusive end):

```
contig    start    end    /path/to/sequence.fa
```

Paths come from the mumemto filelist embedded in the `.lengths` file.

**Fetch sequences**

```bash
shredtools fasta regions/prefix.bed -o fasta_out/
# or from an AGC archive:
shredtools fasta regions/prefix.bed -o fasta_out/ --agc archive.agc -t 8
```

---

## Other commands

| Command | Description |
|---------|-------------|
| `shredtools shred` | Vertically slice a pangenome into homologous strips (BED/FASTA) |
| `shredtools sort` | Sort a `.bumbl` by start on a reference column |
| `shredtools filter` | Keep only collinear-block MUMs in a `.bumbl` |
| `shredtools index` | Build a `.bumbl.bi` index for interval queries |
| `shredtools extract` | Extract homologous regions from a query interval (BED) |
| `shredtools fasta` | Fetch FASTA sequences for an extract BED |
| `shredtools enhance` | Fill gaps between collinear MUMs with local mumemto |
| `shredtools stats` | Print MUM file metadata and associated sidecar paths |

Run `shredtools <command> -h` for full options.

---

## Getting help

If you run into issues or have questions, open a GitHub issue on this repository or reach out to vshivak1 [at] jhu.edu.
