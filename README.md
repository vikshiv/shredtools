# **shredtools**: pangenome coordinates from multi-MUMs

<img src="logo/shredtools_logo.png" alt="Shredtools logo" width="130" align="right"/>

Shredtools is a toolkit for querying a pangenome using multi-MUMs computed by [Mumemto](https://github.com/vikshiv/mumemto). It takes in a set of multi-MUMs (`.bumbl)`, builds an index (`.bumbl.bi`) for fast lookup, and enables operations like extracting syntenic regions across the pangenome, slicing a set of assemblies into smaller subunits (shreds), and refining the multiple alignment by finding more exact matches recursively.

Shredtools expects a set of multi-MUMs in a `.bumbl` file (use `mumemto -b` or `mumemto convert`) and a corresponding `.lengths`. See [Mumemto](https://github.com/vikshiv/mumemto) for details on file formats.

**Browser:** try the `shredtools` extract web app at [https://vikshiv.github.io/shredtools/](https://vikshiv.github.io/shredtools/)

---

## Installation

Requires [mumemto](https://github.com/vikshiv/mumemto) at runtime (not installed automatically by pip, can be installed via conda, pip, or from source). Python 3.9+.

```bash
git clone https://github.com/vikshiv/shredtools.git
cd shredtools
pip install -e .
shredtools -h
```

You can also run `python -m shredtools`.

**Pip dependencies:** `numpy`, `pysam`, `tqdm`

**Optional:** `matplotlib` (for `extract --plot` / `--plot-full`); `.fai` indexes beside reference FASTAs, or the `agc` tool (`fasta --agc`) for sequence extraction (can be installed via bioconda).

---

## Quick start

After running mumemto on a collection of assemblies (e.g. `mumemto assemblies/*.fa -o pangenome -b`), build an index and extract a region:

```bash
# Prepare .bumbl for indexing (sort and compute collinear blocks to filter out non-coordinate system MUMs)
shredtools filter pangenome.bumbl -o pangenome_sorted.bumbl

# Build .bumbl.bi (default: multi-index)
shredtools index --multi pangenome_sorted.bumbl -v

# Extract syntenic regions, then fetch FASTA
shredtools extract pangenome_sorted.bumbl -s 0 -r chr1:1000000-2000000 -o regions/prefix -l pangenome.lengths
shredtools fasta regions/prefix.bed -o fasta_out/
```

> [!TIP]  
> You can overwrite the original `bumbl` file with `-i`

---

## Getting started

### Index a `.bumbl` file

`extract` on a `.bumbl` input requires a `.bumbl.bi` index alongside it. The index maps genomic bins to row ranges in the `.bumbl` file so only subsets of MUMs are loaded at a time.

**Prerequisites**

The input `bumbl` file should be pre-processed before indexing. `shredtools filter` can both sort and filter for collinear multi-MUMs in a single command. For completeness, the following prerequisites are required prior to indexing:

1. **Sorted rows.** MUM rows must be sorted by position in one of the sequences. Use `shredtools sort`.
2. **Optional filter.** `shredtools filter` keeps only MUMs that belong to collinear blocks, which form the coordinate system and filters out potentially spurious matches.

```bash
# achieves both pre-processing steps in one command, in place
shredtools filter -i pangenome.bumbl

# split into two steps:
shredtools sort pangenome.bumbl -s 0 -o pangenome.sorted.bumbl
shredtools filter -i pangenome.sorted.bumbl
# optional; use -o explicitly
```

**Build the index**
There are two index types. A single-index is based on one reference assembly and can handle queries with respect to regions in the reference coordinates. A multi-index is a slightly larger index that can accept queries of any region from any assembly.

```bash
shredtools index pangenome.coll.bumbl --multi # multi-index
shredtools index pangenome.coll.bumbl -s 0 --single # single-index
```


| Flag                   | Default             | Description                                                           |
| ---------------------- | ------------------- | --------------------------------------------------------------------- |
| `-s` / `--seq-idx`     | `0`                 | Reference assembly for single-index                                   |
| `--multi` / `--single` | `--multi` (default) | Multi-index: per-sequence bins (used by `extract`)                    |
| `-w` / `--bin-width`   | `1000000`           | Genomic bin width (bp) (large = smaller index, more memory per query) |
| `-o` / `--output`      | `<bumbl>.bi`        | Output path                                                           |
| `-v` / `--verbose`     | off                 | Progress on stderr                                                    |


> [!WARNING]
> A `.bumbl.bi` file is a binary index corresponding to an associated `.bumbl` **row order**. If the order of the `.bumbl` file changes (such as running `shredtools sort`), then the index is invalid and needs to be re-generated. 

For the full index specifications, see **[bumbl_index/bumbl_index.md](bumbl_index/bumbl_index.md)**.

---

### Extract homologous regions

Given a query interval on one pangenome sequence, `extract` finds the bounding multi-MUMs and reports the syntenic interval on each selected genome as BED. Input files can be remote URLs (see the [Index Zone](https://benlangmead.github.io/aws-indexes/mumemto) for pre-built indexes over HPRC2).

**Required inputs**


| File                 | Description                                                      |
| -------------------- | ---------------------------------------------------------------- |
| `pangenome.bumbl`    | Positional argument, can be a remote file url                    |
| `pangenome.bumbl.bi` | Default: `<mum_file>.bi`, or `-b`                                |
| `pangenome.lengths`  | Default: `<stem>.lengths`, or `-l` (mumemto multilengths format) |


**Region format.** `--range` is **contig-local** on the sequence named by `-s`: `contig:start-end` (e.g. `chr1:1000000-2000000`).

```bash
shredtools extract pangenome.coll.bumbl \
  -s 0 \
  -r chr1:1000000-2000000 \
  -o regions/prefix \
  -l pangenome.lengths
```


| Flag               | Required | Description                                                            |
| ------------------ | -------- | ---------------------------------------------------------------------- |
| `-s` / `--seq-idx` | yes      | Pangenome column that defines the query coordinate system              |
| `-r` / `--range`   | yes      | `contig:start-end` on that sequence                                    |
| `-o` / `--output`  | no       | Output prefix → `<prefix>.bed`; writes to stdout if omitted            |
| `--plot`           | no       | Write `<prefix>_extract_synteny.pdf` (requires `-o`, needs matplotlib) |


> [!TIP]
> A `lengths` file or `bumbl.bi` index can be passed in with `-l` or `-b` respectively, if not automatically detected with the same input prefix.

**More examples**

```bash
# BED to stdout; only genomes 0, 1, 2
shredtools extract pangenome.coll.bumbl -s 0 -r chr1:1-50000 -x 0 1 2

# With synteny plot
shredtools extract pangenome.coll.bumbl -s 0 -r chr1:1000000-2000000 \
  -o regions/prefix -l pangenome.lengths --plot

# Remote .bumbl and index (HTTP range requests)
shredtools extract https://url/to/pangenome.bumbl -s 0 -r chr1:1-1000 \
  -b https://url/to/pangenome.bumbl.bi -l pangenome.lengths
```

When the query interval does not align exactly to MUM boundaries, `extract` prints `left margin` and/or `right margin` on stderr (distance in bp from the query edge to the enclosing MUM). This serves as a guide for the bounding region around the extracted regions.

**BED output format**

Four tab-separated columns (0-based start, exclusive end):

```
contig    start    end    /path/to/sequence.fa
```

Paths come from the `.lengths` file.

**Fetch sequences**

You can extract the FASTA sequences of extracted regions using `shredtools fasta`. If the pangenome is compressed with an AGC archive (see [AGC](https://github.com/refresh-bio/agc)), you can provide it via the `--agc` option to fetch sequences directly from the archive.

```bash
shredtools fasta regions/prefix.bed -o fasta_out/
# or from an AGC archive:
shredtools fasta regions/prefix.bed -o fasta_out/ --agc archive.agc -t 8
```

---

## Other commands


| Command              | Description                                             |
| -------------------- | ------------------------------------------------------- |
| `shredtools shred`   | Vertically slice a pangenome into syntenic strips (BED) |
| `shredtools sort`    | Sort a `.bumbl` by position in an assembly              |
| `shredtools filter`  | Keep only collinear-block MUMs in a `.bumbl`            |
| `shredtools index`   | Build a `.bumbl.bi` index for interval queries          |
| `shredtools extract` | Extract syntenic regions from a query interval (BED)    |
| `shredtools fasta`   | Fetch FASTA sequences for an extract/shred BED          |
| `shredtools enhance` | Fill gaps between collinear MUMs with local mumemto     |
| `shredtools stats`   | Print MUM file metadata and associated sidecar paths    |


Run `shredtools <command> -h` for full options.

---

## Getting help

If you run into issues or have questions, open a GitHub issue on this repository or reach out to vs632 [at] cam.ac.uk.
