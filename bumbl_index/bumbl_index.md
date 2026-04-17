# BUMBL index (`.bumbl.bi`)

A `.bumbl.bi` file is a **binary index** over a `.bumbl` MUM file. It stores information for accessing MUMs in a specific genomic interval, without reading the whole file into memory.

All multi-byte integers are **little-endian `uint64`** unless noted otherwise.

## Index file header

Every `.bumbl.bi` file starts with a fixed **4 × `uint64`** header:

```
┌────────────────────────────────────────────────────────────────────────┐
│  FORMAT        (uint64)  Bytes: b"bumblbi" + FORMAT (1 byte)           │
│  CHECKSUM      (uint64)  SHA-256(lengths[:min(1000,n_mums)]) → u64     │
│  BIN_WIDTH     (uint64)  Genomic bin width, e.g. 1_000_000 bp          │
│  NUM_SEQS      (uint64)  Number of sequences in the corresponding      │
│                          `.bumbl` file                                 │
└────────────────────────────────────────────────────────────────────────┘
```

### `FORMAT`

`FORMAT` is a single `uint64` whose **raw 8 bytes** (little-endian) are:

- first 7 bytes: ASCII `b"bumblbi"`
- last 1 byte: `FORMAT` (0 = single-index, 1 = multi-index)

### `CHECKSUM`

`CHECKSUM` is a hash of the first 1000 `.bumbl` length values, used to sanity-check that the index corresponds to the same `.bumbl` **MUM row order** as at index-build time (the index depends on that order).


| Variant          | When to use                                                          | Coordinate system                                                    |
| ---------------- | -------------------------------------------------------------------- | -------------------------------------------------------------------- |
| **Single-index** | Querying regions in a *single* reference sequence                    | Linear positions on that reference                                   |
| **Multi-index**  | Querying regions across any sequence; interval queries use `seq_idx` | Per-sequence linear coordinates (same as your `.lengths` / manifest) |


---

## Structure of the `.bumbl` file

The on-disk layout is what `MUMdata.parse_bums` and `parse_bumbl_generator` in `mumemto/utils.py` read. Row **i** is one MUM; columns **j** are sequences.

### Layout

```
HEADER
┌────────────────────────────────────────────────────────────────────────┐
│  FLAGS (uint16)               Packed bit flags (see FLAGS table)       │
│  N_SEQS (uint64)              Number of aligned sequences              │
│  N_MUMS (uint64)              Number of MUM rows                       │
└────────────────────────────────────────────────────────────────────────┘
DATA (row-major MUM order)
┌────────────────────────────────────────────────────────────────────────┐
│  LENGTHS                      N_MUMS × uint32 — length per MUM         │
│  STARTS                       N_MUMS × N_SEQS × int64 (-1 = absent)    │
│  STRANDS                      ⌈N_MUMS × N_SEQS / 8⌉ bytes, bit-packed  │
└────────────────────────────────────────────────────────────────────────┘
OPTIONAL (if FLAGS.coll_blocks)
┌────────────────────────────────────────────────────────────────────────┐
│  NUM_BLOCKS (uint64)          If FLAGS.coll_blocks                     │
│  BLOCKS                       NUM_BLOCKS × 2 × uint32 (see mumemto)    │
└────────────────────────────────────────────────────────────────────────┘
```

### Semantics

- **Lengths** and **starts** are stored in **row-major MUM order**: all lengths first, then starts `(mum 0, seq 0)…(mum 0, seq N-1)`, then `(mum 1, seq 0)…`, etc.
- **Strands** are **packed bits** in the same linear order as `strands.reshape(n_mums * n_seqs)`.
- A start of **-1** means the MUM does not occur on that sequence.
- **MUMdata** loads into `lengths`, `starts` (`num_mums × num_seqs`), `strands` (`bool`), and optional block metadata when `coll_blocks` is set.

The `.bumbl.bi` index only stores **indices into this row order** (half-open ranges `[mum_start, mum_end)`), not raw coordinates. The checksum in the index is used to verify the row order in bumbl file is intact.

---

## Single-index format (`FORMAT = 0`)

**Requirement:** MUM rows must be sorted by the reference column you index on, conventionally the first sequence, strictly increasing.

### Layout

```
HEADER
┌────────────────────────────────────────────────────────────────────────┐
│  FORMAT (uint64)        b"bumblbi" + 0                                  │
│  CHECKSUM (uint64)            See "Index file header"                  │
│  BIN_WIDTH (uint64)           Genomic bin width, e.g. 1_000_000 bp     │
│  NUM_SEQS (uint64)            Number of sequences in the `.bumbl`      │
│  SEQ_IDX (uint64)             Which `.bumbl` sequence column is indexed│
└────────────────────────────────────────────────────────────────────────┘
DATA
┌────────────────────────────────────────────────────────────────────────┐
│  OFFSET[0] … OFFSET[NUM_BINS]       First MUM row index per bin        │
└────────────────────────────────────────────────────────────────────────┘
```

- **BIN_WIDTH**: size of each genomic bin (smaller **BIN_WIDTH** ⇒ larger index).
- **NUM_BINS**: `len(OFFSET) - 1` (not stored explicitly).
- **OFFSET[i]**: **MUM row index** (0-based) of the **first** MUM whose start in the indexed column (`SEQ_IDX`) is **≥ `i × BIN_WIDTH`**.

To query an interval `[s, e)` on the indexed sequence coordinate system:

1. `bin_lo = s // BIN_WIDTH`, `bin_hi = (e − 1) // BIN_WIDTH`.
2. MUM rows to consider lie in `[OFFSET[bin_lo], OFFSET[bin_hi + 1])` in the `.bumbl` row order.

---

## Multi-index format (`FORMAT = 1`)

Use when the `.bumbl` has many sequences and each query supplies **seq_idx** and an interval in **that sequence’s** linear coordinates (consistent with your `.lengths` / manifest).

### Layout

```
HEADER
┌────────────────────────────────────────────────────────────────────────┐
│  FORMAT (uint64)        b"bumblbi" + 1                                 │
│  CHECKSUM (uint64)            See "Index file header"                  │
│  BIN_WIDTH (uint64)           Bin width in per-sequence coordinates    │
│  NUM_SEQS (uint64)            Number of sequences (documents)          │
└────────────────────────────────────────────────────────────────────────┘
DOC TABLE  (NUM_SEQS + 1) × uint64
┌────────────────────────────────────────────────────────────────────────┐
│  DOC_OFF[0] … DOC_OFF[NUM_SEQS]   Byte offsets into this index file    │
│  Sequence s spans [DOC_OFF[s], DOC_OFF[s + 1])                         │
└────────────────────────────────────────────────────────────────────────┘
```

For one sequence, `doc_start = DOC_OFF[seq_idx]`:

```
PER-SEQUENCE HEADER
┌────────────────────────────────────────────────────────────────────────┐
│  BIN_NUM (uint64)             Bin count along this sequence’s axis     │
└────────────────────────────────────────────────────────────────────────┘
PER-SEQUENCE DATA
┌────────────────────────────────────────────────────────────────────────┐
│  BOUNDARY[0] … BOUNDARY[BIN_NUM]   (BIN_NUM + 1) × uint64              │
│  range_set_base = doc_start + 8 × (BIN_NUM + 2)                        │
└────────────────────────────────────────────────────────────────────────┘
SET OF RANGES (per bin; byte offsets relative to range_set_base)
┌────────────────────────────────────────────────────────────────────────┐
│  Per bin b: byte range [BOUNDARY[b], BOUNDARY[b + 1]) relative to      │
│  range_set_base — concatenated (MUM_START, MUM_END) pairs, 16 B each   │
└────────────────────────────────────────────────────────────────────────┘
```

For bin index `b`: `off0 = BOUNDARY[b]`, `off1 = BOUNDARY[b + 1]` (seek to `doc_start + 8 × (1 + b)`, read two `uint64`s). The byte range `[range_set_base + off0, range_set_base + off1)` is **one bin’s** set of ranges (a contiguous run of `(mum_start, mum_end)` pairs); concatenate across bins for the full query.

### Query steps

1. Map `[s, e)` to bins: `bin_start = s // BIN_WIDTH`, `bin_end = (e − 1) // BIN_WIDTH`.
2. For each bin from `bin_start` through `bin_end` (inclusive), decode that bin’s set of ranges and concatenate the `(mum_start, mum_end)` pairs.
3. Resolve coordinates in `.bumbl`, e.g. `starts[:, 0][mum_start:mum_end]` (sort if needed).

**Checks** (reference reader): `DOC_OFF[seq_idx] < DOC_OFF[seq_idx + 1]`, `off1 ≥ off0`, each bin’s range slice length divisible by 16.

### Visual (one sequence)

Linear axis (same bins as the query uses):

```
        0              BIN_WIDTH           2 × BIN_WIDTH           3 × BIN_WIDTH
        |------------------|------------------|------------------|------------------|
              bin 0              bin 1              bin 2              bin 3
```

Contiguous **set of ranges** in file order (each cell is one `(mum_start, mum_end)` pair, 16 bytes; vertical bars separate bins). A bin can hold several pairs when MUM rows are split into non-contiguous runs in `.bumbl` row order:

```
  byte offset from range_set_base
                0        16       32        48       64       80            96
              ┌───────────┬───────────────────┬───────────────────────────────────┬───────────┐
              │ (m0, m1)  │ (m2, m3) (m4, m5) │ (m6, m7) (m8, m9) (m10, m11)      │ (m12, m13)│
              └───────────┴───────────────────┴───────────────────────────────────┴───────────┘
  bin         :      0    │         1         │                2                  │     3
              ↑           ↑                   ↑                                   ↑
              BIN[0]     BIN[1]              BIN[2]                              BIN[3]
```

`BOUNDARY[b]` and `BOUNDARY[b + 1]` delimit the byte span for bin `b` inside the set of ranges (not the number of pairs—decode pairs until that span ends).