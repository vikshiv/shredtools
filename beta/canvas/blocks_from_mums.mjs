/**
 * Build Canvas left/right/strands from denseloaded MUMs.
 * Port of browser_extract/py/query.py find_coll_blocks_simple + blocks_to_endpoint_arrays,
 * with mumsToSingletonBlocks fallback (mumemto browser/extract_zoom.mjs).
 */

/**
 * @typedef {{ nSeqs: number, numMums: number, lengths: number[], starts: number[], strands: number[] }} MumArrays
 */

function mumStart(mums, i, g) {
  return mums.starts[i * mums.nSeqs + g];
}

function mumStrand(mums, i, g) {
  return !!mums.strands[i * mums.nSeqs + g];
}

function mumLen(mums, i) {
  return mums.lengths[i];
}

/** True if every genome column is present on row i. */
function allPresent(mums, i) {
  for (let g = 0; g < mums.nSeqs; g++) {
    if (mumStart(mums, i, g) < 0) return false;
  }
  return true;
}

/**
 * Lightweight collinear blocks (strict: all seqs present).
 * @returns {Array<[number, number]>} inclusive MUM index pairs
 */
export function findCollBlocksSimple(mums, maxBreak = 1000) {
  const n = mums.numMums;
  if (n === 0) return [];
  const ns = mums.nSeqs;
  const blocks = [];
  let i = 0;
  while (i < n) {
    if (!allPresent(mums, i)) {
      i += 1;
      continue;
    }
    let j = i;
    while (j + 1 < n) {
      if (!allPresent(mums, j + 1)) break;
      let ok = true;
      for (let k = 0; k < ns; k++) {
        const prevEnd = mumStart(mums, j, k) + mumLen(mums, j);
        const nxt = mumStart(mums, j + 1, k);
        const gap = nxt - prevEnd;
        if (gap < 0 || gap > maxBreak) {
          ok = false;
          break;
        }
        if (mumStrand(mums, j, k) !== mumStrand(mums, j + 1, k)) {
          ok = false;
          break;
        }
      }
      if (!ok) break;
      j += 1;
    }
    if (j > i) blocks.push([i, j]);
    else if (mumLen(mums, i) >= maxBreak) blocks.push([i, i]);
    i = j + 1;
  }
  return blocks;
}

/** Convert (L,R) blocks → left/right/strands lists-of-lists (-1 = absent). */
export function blocksToEndpointArrays(mums, blocks) {
  const left = [];
  const right = [];
  const strands = [];
  for (const [L, R] of blocks) {
    const ls = [];
    const rs = [];
    const st = [];
    for (let k = 0; k < mums.nSeqs; k++) {
      let a = mumStart(mums, L, k);
      let b = mumStart(mums, R, k);
      if (a < 0 || b < 0) {
        ls.push(-1);
        rs.push(-1);
        st.push(true);
        continue;
      }
      let end = b + mumLen(mums, R);
      if (a > end) [a, end] = [end, a];
      ls.push(a);
      rs.push(end);
      st.push(mumStrand(mums, L, k));
    }
    left.push(ls);
    right.push(rs);
    strands.push(st);
  }
  return { left, right, strands, nBlocks: left.length };
}

/**
 * One ribbon per MUM with ≥2 genomes overlapping windows.
 * Coords may be absolute (then windows shift) or already window-relative
 * (pass windows as [0, span] per genome).
 */
export function mumsToSingletonBlocks(mums, windows) {
  const nSeqs = mums.nSeqs;
  const left = [];
  const right = [];
  const strands = [];
  const shifts = windows.map((w) => (w ? w[0] : 0));

  for (let i = 0; i < mums.numMums; i++) {
    const nl = new Array(nSeqs).fill(null);
    const nr = new Array(nSeqs).fill(null);
    const st = new Array(nSeqs).fill(true);
    let present = 0;
    for (let g = 0; g < nSeqs; g++) {
      const x = mumStart(mums, i, g);
      if (x < 0) continue;
      const L = mumLen(mums, i);
      const w = windows[g];
      if (!w) continue;
      const a = x;
      const b = x + L;
      if (a < w[1] && b > w[0]) {
        const ca = Math.max(a, w[0]);
        const cb = Math.min(b, w[1]);
        if (cb > ca) {
          nl[g] = ca - shifts[g];
          nr[g] = cb - shifts[g];
          st[g] = mumStrand(mums, i, g);
          present++;
        }
      }
    }
    if (present < 2) continue;
    left.push(nl);
    right.push(nr);
    strands.push(st);
  }
  const span = Math.max(1, ...windows.map((w) => (w ? w[1] - w[0] : 0)));
  return { left, right, strands, dataXMax: span, nBlocks: left.length };
}

function singletonWindows(mums, windowSpans) {
  const spans = windowSpans || [];
  const windows = [];
  for (let g = 0; g < mums.nSeqs; g++) {
    const span =
      spans[g] != null
        ? spans[g]
        : Math.max(
            1,
            ...Array.from({ length: mums.numMums }, (_, i) => {
              const s = mumStart(mums, i, g);
              return s < 0 ? 0 : s + mumLen(mums, i);
            })
          );
    windows.push([0, span]);
  }
  return windows;
}

/**
 * Build ribbons from denseloaded MUMs.
 * Expects window-relative MUM starts (shredtools denseload).
 *
 * Default for extract: **singleton** (matches mumemto densify + matplotlib
 * per-MUM ribbons). Collinear merges of all-present MUMs across hundreds of
 * genomes paint huge blocks and exaggerate any strand flip into false
 * whole-window "inversions".
 *
 * @param {MumArrays} mums
 * @param {{ maxBreak?: number, minCollinearBlocks?: number, windowSpans?: number[], preferCollinear?: boolean }} [opts]
 */
export function blocksFromMums(mums, opts = {}) {
  const windows = singletonWindows(mums, opts.windowSpans);
  if (!opts.preferCollinear) {
    const single = mumsToSingletonBlocks(mums, windows);
    return { ...single, source: "singleton" };
  }
  const maxBreak = opts.maxBreak ?? 1000;
  const minCollinear = opts.minCollinearBlocks ?? 1;
  const coll = findCollBlocksSimple(mums, maxBreak);
  if (coll.length >= minCollinear) {
    const ep = blocksToEndpointArrays(mums, coll);
    return { ...ep, source: "collinear" };
  }
  const single = mumsToSingletonBlocks(mums, windows);
  return { ...single, source: "singleton" };
}
