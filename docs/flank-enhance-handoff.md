# Handoff: enhance flanks (branch `enhance-flanks`)

## Goal

Extend `shredtools enhance` so it also searches **left/right flanks** of each collinear block — sequence before the first MUM and after the last MUM — not only internal gaps between consecutive MUMs.

Today flanks are never considered. On cichlid `plus.enhanced.mums` that leaves ~1.7 kb left and ~11 kb right empty on genome 0 even after enhance; purple fills in `archive/viz_enhanced_mums.py` only appear in interior gaps.

## Repo / branch

- Repo: `/homes/shivakumar/software/shredtools`
- Branch: `enhance-flanks` (from `dev` @ viz commit)
- Core file: [`shredtools/enhance.py`](../shredtools/enhance.py)
- Viz (optional check): `archive/viz_enhanced_mums.py`
- Example data: `/homes/shivakumar/hps-vshiv/cichlid/analysis/plus.enhanced.mums` + `plus.lengths`

## Current behavior (do not break)

1. Load MUMs + multilengths; find/use collinear `mums.blocks`.
2. `compute_gap_coords_and_lengths`: for each block, only pairs `(i, i+1)` for `i in [s, e)`.
3. Keep gaps where any genome has length ≥ `min_gap_length`.
4. Optional `--outlier-factor`: drop genomes with gap `< med - factor*MAD` (MAD 0 → 1). Outlier assemblies stay `-1` → partial MUMs.
5. Fetch gap FASTA per kept genome → local `mumemto` → keep collinear locals → map to global coords → merge + sort.

## Proposed flank design

Treat each block’s flanks as two extra “gaps”:

| flank | per-genome interval |
|-------|---------------------|
| Left | `[0, start_of_first_MUM_in_block)` (respect strand/orientation the same way internal gaps do) |
| Right | `[end_of_last_MUM_in_block, genome_length)` |

Reuse the same pipeline as internal gaps: outlier mask, fetch, `enhance_gap`, merge. Prefer extending `compute_gap_coords_and_lengths` (or a sibling that returns flank rows in the same `coords` / `gap_lengths` shape) so `build_gap_work` / workers stay unchanged.

### Design choices (pick sensible defaults)

- **CLI**: e.g. `--flanks` / `--no-flanks` (default on or off — prefer **opt-in `--flanks`** first so existing runs stay identical).
- **Orientation**: internal gaps already handle flipped MUMs (`right_starts < left_starts`). Flanks must use the same left/right endpoint logic for reverse-strand runs.
- **Multi-contig**: flanks are in **concat-global** coords (same as MUM starts); `_fetch_gap_intervals` already converts via multilengths. A flank that spans contig boundaries may need the same handling as long internal gaps (fetch may fail / AssertionError → that genome becomes partial).
- **Outliers**: apply the same `--outlier-factor` to flank lengths (short flank = possible truncation/assembly difference).
- **Empty / tiny flanks**: still gated by `min_gap_length` (user can pass `0`).
- **Multiple blocks**: emit left+right flank per block (not only the global first/last MUM).

## Acceptance checks

1. Without `--flanks`, output matches previous enhance on a fixture / `plus.enhanced.mums` with large `min_gap_length` (no new work) or bit-identical when no gaps enhanced.
2. With `--flanks` and `min_gap_length 0` on cichlid plus: new MUMs appear with genome-0 starts `< first_MUM` and/or `>= last_MUM_end` (or partials with other genomes present in those regions).
3. Plot with `archive/viz_enhanced_mums.py` shows purple (or new blue if strict) extending into previous white flanks when matches exist.
4. Tests: unit-test flank coord construction (synthetic 2–3 genomes, one block); optional CLI smoke like existing `test_shred_enhance_fasta.py`.

## Out of scope

- Changing MAD outlier semantics (separate issue: MAD≈0 makes factor ≈ bp-below-median).
- Wiring viz into the main CLI.
- Enhancing *between* separate collinear blocks (only flanks of each block).
