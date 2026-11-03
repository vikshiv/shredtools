/**
 * Nested extract navigation (mumemto extract_zoom.mjs), adapted for null absents.
 */

function bisectRight(arr, x) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] <= x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

function absent(a, b) {
  return a == null || b == null || a === -1 || b === -1;
}

/**
 * Derive per-genome windows from ribbons overlapping [start,end] on seqIdx.
 * Coords are in the same space as left/right (e.g. home-justified extract).
 */
export function extractWindowsFromBlocks(left, right, seqIdx, start, end, nSeqs) {
  const s = Math.min(start, end);
  const e = Math.max(start, end);
  if (!(e > s)) {
    return {
      ok: false,
      error: "empty selection",
      windows: [],
      hitBlocks: [],
      query: { seqIdx, start: s, end: e },
    };
  }

  const ids = [];
  const starts = [];
  const ends = [];
  for (let i = 0; i < left.length; i++) {
    let a = left[i][seqIdx];
    let b = right[i][seqIdx];
    if (absent(a, b)) continue;
    if (a > b) [a, b] = [b, a];
    ids.push(i);
    starts.push(a);
    ends.push(b);
  }
  if (!ids.length) {
    return {
      ok: false,
      error: "no blocks on genome",
      windows: [],
      hitBlocks: [],
      query: { seqIdx, start: s, end: e },
    };
  }

  const order = Array.from({ length: ids.length }, (_, i) => i);
  order.sort((i, j) => starts[i] - starts[j] || ends[i] - ends[j]);
  const S = order.map((i) => starts[i]);
  const E = order.map((i) => ends[i]);

  let leftPos = bisectRight(S, s) - 1;
  if (leftPos < 0) leftPos = 0;
  let rightPos = bisectRight(E, e);
  if (rightPos >= order.length) rightPos = order.length - 1;
  if (rightPos < leftPos) rightPos = leftPos;

  const hit = [];
  for (let k = 0; k < order.length; k++) {
    if (S[k] < e && E[k] > s) hit.push(ids[order[k]]);
  }
  if (hit.length) {
    const hitSet = new Set(hit);
    let first = -1;
    let last = -1;
    for (let k = 0; k < order.length; k++) {
      if (!hitSet.has(ids[order[k]])) continue;
      if (first < 0) first = k;
      last = k;
    }
    if (first >= 0) {
      leftPos = first;
      rightPos = last;
    }
  }

  const iL = order[leftPos];
  const iR = order[rightPos];
  const blockL = ids[iL];
  const blockR = ids[iR];

  let leftOffset = 0;
  let rightOffset = 0;
  if (s < ends[iL]) leftOffset = s - starts[iL];
  if (e >= starts[iR]) rightOffset = e - starts[iR];
  leftOffset = Math.max(0, Math.min(leftOffset, ends[iL] - starts[iL]));
  rightOffset = Math.max(0, Math.min(rightOffset, ends[iR] - starts[iR]));

  const windows = new Array(nSeqs).fill(null);
  for (let g = 0; g < nSeqs; g++) {
    let la = left[blockL][g];
    let lb = right[blockL][g];
    let ra = left[blockR][g];
    let rb = right[blockR][g];
    if (absent(la, lb) || absent(ra, rb)) {
      windows[g] = null;
      continue;
    }
    if (la > lb) [la, lb] = [lb, la];
    if (ra > rb) [ra, rb] = [rb, ra];
    const w0 = la + leftOffset;
    const w1 = ra + rightOffset;
    windows[g] = w0 <= w1 ? [w0, w1] : [w1, w0];
  }

  let q0 = Math.min(s, e);
  let q1 = Math.max(s, e);
  const gLo = Math.min(...starts);
  const gHi = Math.max(...ends);
  q0 = Math.min(gHi, Math.max(gLo, q0));
  q1 = Math.min(gHi, Math.max(gLo, q1));
  if (!(q1 > q0)) {
    return {
      ok: false,
      error: "selection outside genome blocks",
      windows: [],
      hitBlocks: [],
      query: { seqIdx, start: s, end: e },
    };
  }
  windows[seqIdx] = [q0, q1];

  return {
    ok: true,
    windows,
    hitBlocks: hit.length ? hit : ids.slice(Math.min(iL, iR), Math.max(iL, iR) + 1),
    blockL,
    blockR,
    query: { seqIdx, start: q0, end: q1 },
  };
}

/**
 * Justify blocks so each genome's window start → 0. Keep flanks when context=true.
 */
export function projectBlocksToExtract(
  left,
  right,
  strands,
  windows,
  genomeOrder,
  { context = true } = {}
) {
  const nSeqs = windows.length;
  const shifts = windows.map((w) => (w ? w[0] : 0));
  const span = Math.max(1, ...windows.map((w) => (w ? w[1] - w[0] : 0)));

  const outL = [];
  const outR = [];
  const outS = [];
  let xMin = 0;
  let xMax = span;

  for (let i = 0; i < left.length; i++) {
    if (!context) {
      let keep = false;
      for (let g = 0; g < nSeqs; g++) {
        const w = windows[g];
        if (!w) continue;
        let a = left[i][g];
        let b = right[i][g];
        if (absent(a, b)) continue;
        if (a > b) [a, b] = [b, a];
        if (a < w[1] && b > w[0]) {
          keep = true;
          break;
        }
      }
      if (!keep) continue;
    }

    const nl = new Array(nSeqs).fill(null);
    const nr = new Array(nSeqs).fill(null);
    for (let g = 0; g < nSeqs; g++) {
      const w = windows[g];
      let a = left[i][g];
      let b = right[i][g];
      if (!w || absent(a, b)) continue;
      if (a > b) [a, b] = [b, a];
      if (context) {
        nl[g] = a - shifts[g];
        nr[g] = b - shifts[g];
      } else {
        const ca = Math.max(a, w[0]);
        const cb = Math.min(b, w[1]);
        if (cb <= ca) continue;
        nl[g] = ca - shifts[g];
        nr[g] = cb - shifts[g];
      }
      xMin = Math.min(xMin, nl[g], nr[g]);
      xMax = Math.max(xMax, nl[g], nr[g]);
    }
    let present = 0;
    for (let g = 0; g < nSeqs; g++) if (nl[g] != null) present++;
    if (present < 2) continue;
    outL.push(nl);
    outR.push(nr);
    outS.push(strands[i].slice());
  }

  const centering = genomeOrder.map(() => 0);
  const regionEnds = windows.map((w) => (w ? w[1] - w[0] : null));
  const pad = Math.max(span * 0.05, 1000);
  return {
    left: outL,
    right: outR,
    strands: outS,
    centering,
    dataXMin: xMin - pad,
    dataXMax: xMax + pad,
    regionX0: 0,
    regionX1: span,
    regionEnds,
    nBlocks: outL.length,
  };
}
