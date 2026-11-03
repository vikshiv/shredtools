/** Decode overview.json and build ribbon polygons in JS (fast path for overview). */

function b64ToBytes(b64) {
  const bin = atob(b64);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

export function decodeOverview(raw) {
  const nSeqs = raw.n_seqs;
  const nBlocks = raw.n_blocks;
  const leftBuf = b64ToBytes(raw.left_i64);
  const rightBuf = b64ToBytes(raw.right_i64);
  const strandBuf = b64ToBytes(raw.strands_u8);
  // Copy into aligned buffers for BigInt64Array
  const leftAligned = new ArrayBuffer(leftBuf.byteLength);
  const rightAligned = new ArrayBuffer(rightBuf.byteLength);
  new Uint8Array(leftAligned).set(leftBuf);
  new Uint8Array(rightAligned).set(rightBuf);
  const leftView = new BigInt64Array(leftAligned);
  const rightView = new BigInt64Array(rightAligned);

  const left = new Array(nBlocks);
  const right = new Array(nBlocks);
  const strands = new Array(nBlocks);
  for (let i = 0; i < nBlocks; i++) {
    const ls = new Array(nSeqs);
    const rs = new Array(nSeqs);
    const st = new Array(nSeqs);
    const base = i * nSeqs;
    for (let j = 0; j < nSeqs; j++) {
      ls[j] = Number(leftView[base + j]);
      rs[j] = Number(rightView[base + j]);
      st[j] = !!strandBuf[base + j];
    }
    left[i] = ls;
    right[i] = rs;
    strands[i] = st;
  }
  return {
    n_seqs: nSeqs,
    n_blocks: nBlocks,
    n_mums: raw.n_mums || 0,
    seq_lengths: raw.seq_lengths,
    labels: raw.labels,
    left,
    right,
    strands,
  };
}

function pointsToPoly(points) {
  const starts = points.map((p) => p[0]);
  const ends = points.map((p) => p[1]).reverse();
  return starts.concat(ends);
}

/** Reject zero-area / empty ribbons. Intentional inversion bowties are allowed. */
function polyOk(poly) {
  if (!poly || poly.length < 3) return false;
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  for (const [x, y] of poly) {
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
  }
  // Must span at least one genome step and have positive width
  if (maxY - minY < 0.5) return false;
  if (maxX - minX <= 0) return false;
  return true;
}

/**
 * @param {boolean} twist  Swap start/end on genomes after the first so an
 *   inverted pair becomes a bowtie (left↔right cross) instead of a straight ribbon.
 *   Bundle-of-MUMs inversions look twisted naturally; a single collinear block
 *   needs this geometric cross or it draws as a flat red bar.
 */
function pushPoly(polygons, colors, points, color, twist = false) {
  if (points.length < 2) return;
  let pts = points;
  if (twist) {
    pts = points.map((p, i) => (i === 0 ? p : [p[1], p[0]]));
  }
  const poly = pointsToPoly(pts);
  if (!polyOk(poly)) return;
  // Skip ribbons that flare wildly between adjacent genomes (artifact)
  let maxEdge = 0;
  let minEdge = Infinity;
  for (let i = 0; i < pts.length; i++) {
    const w = Math.abs(pts[i][1][0] - pts[i][0][0]);
    maxEdge = Math.max(maxEdge, w);
    minEdge = Math.min(minEdge, w);
  }
  // Only drop grossly inconsistent ribbons (real collinear blocks are similar width)
  if (minEdge > 0 && maxEdge / minEdge > 80) return;
  polygons.push(poly);
  colors.push(color);
}

/**
 * Build ribbons for display rows genomeOrder.
 *
 * Inversion is **pairwise by current adjacency**: the edge between display rows
 * i and i+1 is red+twisted iff those two genomes disagree on strand for the
 * block. Reordering genomes therefore changes which edges are inverted — this
 * is not baked into geometry and must not seed from a fixed reference genome.
 */
export function buildBlockPolygons(
  left,
  right,
  strands,
  genomeOrder,
  centering,
  {
    color = "#2F6FED",
    invColor = "#C62828",
    maxBlocks = 50000,
    /** Overview collinear blocks need a geometric bowtie; extract MUMs match viz_mums (color only). */
    twistInversions = true,
  } = {}
) {
  const polygons = [];
  const colors = [];
  const nDisp = genomeOrder.length;
  const step = left.length > maxBlocks ? Math.ceil(left.length / maxBlocks) : 1;

  for (let bi = 0; bi < left.length; bi += step) {
    let points = [];
    let prevStrand = null;
    for (let idx = 0; idx < nDisp; idx++) {
      const g = genomeOrder[idx];
      let a = left[bi][g];
      let b = right[bi][g];
      const strand = !!strands[bi][g];
      // Absent: null (extract flanks) or classic -1 (overview)
      if (a == null || b == null || a === -1 || b === -1) {
        if (points.length >= 2) pushPoly(polygons, colors, points, color);
        points = [];
        prevStrand = null;
        continue;
      }
      if (a > b) [a, b] = [b, a];
      const c = centering[idx] || 0;
      points.push([
        [c + a, idx],
        [c + b, idx],
      ]);

      if (prevStrand !== null && points.length >= 2) {
        if (prevStrand !== strand) {
          // Orientation flips between neighbors → red edge (optional bowtie twist)
          if (points.length > 2) {
            pushPoly(polygons, colors, points.slice(0, -1), color);
          }
          pushPoly(polygons, colors, points.slice(-2), invColor, twistInversions);
          points = [points[points.length - 1]];
        }
      }
      prevStrand = strand;
    }
    if (points.length >= 2) pushPoly(polygons, colors, points, color);
  }
  return { polygons, colors };
}

export function contigOffsets(order, lengths, spacer = 0) {
  const left = new Array(lengths.length).fill(0);
  let pos = 0;
  for (let i = 0; i < order.length; i++) {
    const cid = order[i];
    left[cid] = pos;
    pos += lengths[cid];
    if (i < order.length - 1) pos += spacer;
  }
  return left;
}

export function globalToContig(pos, contigLengths) {
  if (pos < 0) return [-1, -1];
  let acc = 0;
  for (let cid = 0; cid < contigLengths.length; cid++) {
    const L = contigLengths[cid];
    if (pos < acc + L) return [cid, pos - acc];
    acc += L;
  }
  const last = contigLengths.length - 1;
  return [last, Math.max(0, contigLengths[last] - 1)];
}

/**
 * Per-contig covered bp from concat-global block endpoints.
 * Assigns each present block's span to the contig containing its midpoint
 * (mirrors mumemto weighting MUMs by contig_id).
 *
 * @returns {number[][]} covered[g][cid]
 */
export function contigCoverageFromBlocks(left, right, multilengths) {
  const nSeqs = multilengths.length;
  const covered = multilengths.map((ml) => Array(ml.length).fill(0));
  for (let bi = 0; bi < left.length; bi++) {
    for (let g = 0; g < nSeqs; g++) {
      const a = left[bi][g];
      const b = right[bi][g];
      if (a < 0 || b < 0) continue;
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      const span = Math.max(0, hi - lo);
      if (span <= 0) continue;
      const ml = multilengths[g];
      if (!ml || !ml.length) continue;
      const [cid] = globalToContig(Math.floor((lo + hi) / 2), ml);
      if (cid < 0 || cid >= ml.length) continue;
      covered[g][cid] += span;
    }
  }
  return covered;
}

/**
 * Drop a trailing suffix of empty/sparse contigs per genome (mumemto
 * ``trim_trailing_empty_contigs``). Mid-genome sparse contigs are kept.
 * Concat coords of the kept prefix are unchanged.
 *
 * @returns {{
 *   multilengths: number[][],
 *   contigNames: string[][],
 *   left: number[][],
 *   right: number[][],
 *   strands: boolean[][],
 *   keepN: number[],
 *   dropped: number,
 * }}
 */
export function trimTrailingEmptyContigs(
  left,
  right,
  strands,
  multilengths,
  contigNames,
  { minFrac = 1e-4 } = {}
) {
  const nSeqs = multilengths.length;
  const nBlocks = left.length;
  const coveredAll = contigCoverageFromBlocks(left, right, multilengths);
  const keepN = new Array(nSeqs);
  const newMls = new Array(nSeqs);
  const newNames = new Array(nSeqs);
  let dropped = 0;

  for (let g = 0; g < nSeqs; g++) {
    const ml = multilengths[g] || [];
    const names = (contigNames && contigNames[g]) || [];
    const n = ml.length;
    if (n === 0) {
      keepN[g] = 0;
      newMls[g] = [];
      newNames[g] = [];
      continue;
    }
    const covered = coveredAll[g];
    const total = covered.reduce((a, b) => a + b, 0);
    // No hits at all → leave contig list unchanged (mumemto docstring)
    if (total <= 0) {
      keepN[g] = n;
      newMls[g] = ml.slice();
      newNames[g] = names.slice(0, n);
      continue;
    }
    const thresh = Math.floor(minFrac * total);
    let last = n - 1;
    while (last > 0) {
      const L = Math.max(0, Math.floor(Number(ml[last]) || 0));
      const frac = L > 0 ? covered[last] / L : 0;
      if (covered[last] <= thresh || frac < 0.02) last -= 1;
      else break;
    }
    const kn = last + 1;
    keepN[g] = kn;
    dropped += n - kn;
    newMls[g] = ml.slice(0, kn);
    newNames[g] = names.slice(0, kn);
  }

  // Nullify endpoints past the kept concat length (suffix drop only)
  const newLeft = Array.from({ length: nBlocks }, () => Array(nSeqs).fill(-1));
  const newRight = Array.from({ length: nBlocks }, () => Array(nSeqs).fill(-1));
  const newStrands = Array.from({ length: nBlocks }, () => Array(nSeqs).fill(true));
  for (let g = 0; g < nSeqs; g++) {
    const newLen = newMls[g].reduce((a, b) => a + b, 0);
    for (let bi = 0; bi < nBlocks; bi++) {
      const a = left[bi][g];
      const b = right[bi][g];
      newStrands[bi][g] = strands[bi][g];
      if (a < 0 || b < 0) continue;
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      if (lo >= newLen || hi > newLen) continue;
      newLeft[bi][g] = a;
      newRight[bi][g] = b;
    }
  }

  return {
    multilengths: newMls,
    contigNames: newNames,
    left: newLeft,
    right: newRight,
    strands: newStrands,
    keepN,
    dropped,
  };
}

export function applyLayoutToBlocks(left, right, strands, multilengths, order, orient, spacer = 0) {
  const nSeqs = multilengths.length;
  const nBlocks = left.length;
  const newLeft = Array.from({ length: nBlocks }, () => Array(nSeqs).fill(-1));
  const newRight = Array.from({ length: nBlocks }, () => Array(nSeqs).fill(-1));
  const newStrands = Array.from({ length: nBlocks }, () => Array(nSeqs).fill(true));
  const offsets = multilengths.map((ml, g) => contigOffsets(order[g], ml, spacer));

  for (let g = 0; g < nSeqs; g++) {
    const ml = multilengths[g];
    for (let bi = 0; bi < nBlocks; bi++) {
      const a = left[bi][g];
      const b = right[bi][g];
      if (a < 0 || b < 0) continue;
      const [cidA, locA] = globalToContig(a, ml);
      const [cidB] = globalToContig(Math.max(0, b - 1), ml);
      if (cidA !== cidB) {
        newLeft[bi][g] = a;
        newRight[bi][g] = b;
        newStrands[bi][g] = strands[bi][g];
        continue;
      }
      const L = ml[cidA];
      const o = orient[g][cidA];
      if (o < 0) {
        const span = Math.max(1, b - a);
        let locA2 = L - locA - 1;
        let locB2 = locA2 - span + 1;
        if (locB2 > locA2) [locA2, locB2] = [locB2, locA2];
        newLeft[bi][g] = offsets[g][cidA] + Math.min(locA2, locB2);
        newRight[bi][g] = offsets[g][cidA] + Math.max(locA2, locB2) + 1;
        newStrands[bi][g] = !strands[bi][g];
      } else {
        const [_, locB] = globalToContig(Math.max(0, b - 1), ml);
        newLeft[bi][g] = offsets[g][cidA] + Math.min(locA, locB);
        newRight[bi][g] = offsets[g][cidA] + Math.max(locA, locB) + 1;
        newStrands[bi][g] = strands[bi][g];
      }
    }
  }
  return { left: newLeft, right: newRight, strands: newStrands };
}

export function displayLengths(multilengths, order, spacer = 0) {
  return multilengths.map((ml, g) => {
    const ordered = order[g].map((c) => ml[c]);
    return ordered.reduce((a, b) => a + b, 0) + spacer * Math.max(0, ordered.length - 1);
  });
}

/**
 * Per-genome contig segments in display coordinates (after order/orient/spacer).
 * ``center`` is added to every segment start/end for that genome.
 *
 * @returns {Array<Array<{cid, name, start, end, length}>>} indexed by genome id
 */
export function buildContigSegments(
  multilengths,
  contigNames,
  order,
  spacer = 0,
  centeringByGenome = null
) {
  const out = [];
  for (let g = 0; g < multilengths.length; g++) {
    const ml = multilengths[g];
    const names = contigNames[g] || [];
    const ord = order[g];
    const c0 = centeringByGenome ? centeringByGenome[g] || 0 : 0;
    let pos = c0;
    const segs = [];
    for (let i = 0; i < ord.length; i++) {
      const cid = ord[i];
      const L = intOr(ml[cid]);
      segs.push({
        cid,
        name: names[cid] || `ctg${cid}`,
        start: pos,
        end: pos + L,
        length: L,
        local0: 0,
        isContigStart: true,
        isContigEnd: true,
      });
      pos += L;
      if (i < ord.length - 1) pos += spacer;
    }
    out.push(segs);
  }
  return out;
}

/**
 * Clip display-space contig segments to an extract window [w0, w1], remapped to
 * window-relative [0, w1-w0]. Only edges that are true contig ends keep markers.
 */
export function clipSegmentsToWindow(segmentsByGenome, windows) {
  const out = [];
  for (let g = 0; g < windows.length; g++) {
    const w = windows[g];
    const segs = segmentsByGenome[g] || [];
    if (!w) {
      out.push([]);
      continue;
    }
    const [w0, w1] = w;
    const clipped = [];
    for (const seg of segs) {
      const a = Math.max(seg.start, w0);
      const b = Math.min(seg.end, w1);
      if (!(b > a)) continue;
      clipped.push({
        cid: seg.cid,
        name: seg.name,
        start: a - w0,
        end: b - w0,
        length: b - a,
        local0: a - seg.start,
        isContigStart: a === seg.start,
        isContigEnd: b === seg.end,
      });
    }
    out.push(clipped);
  }
  return out;
}

/**
 * Justify contig segments to extract space (window start → 0) without clipping,
 * so flanks remain for pan context. Contig boundary flags stay true.
 */
export function shiftSegmentsToExtract(segmentsByGenome, windows) {
  const out = [];
  for (let g = 0; g < windows.length; g++) {
    const w = windows[g];
    const segs = segmentsByGenome[g] || [];
    if (!w) {
      out.push([]);
      continue;
    }
    const w0 = w[0];
    out.push(
      segs.map((seg) => ({
        cid: seg.cid,
        name: seg.name,
        start: seg.start - w0,
        end: seg.end - w0,
        length: seg.length,
        local0: seg.local0 != null ? seg.local0 : 0,
        isContigStart: seg.isContigStart !== false,
        isContigEnd: seg.isContigEnd !== false,
      }))
    );
  }
  return out;
}

function intOr(x) {
  return Math.max(0, Math.floor(Number(x) || 0));
}

export function identityOrder(multilengths) {
  return multilengths.map((ml) => [...Array(ml.length).keys()]);
}

export function identityOrient(multilengths) {
  return multilengths.map((ml) => Array(ml.length).fill(1));
}

/** Median-based untangle vs genome 0 (JS port of untangle_bridge.untangle_from_blocks). */
export function untangleFromBlocks(left, right, strands, multilengths, ref = 0) {
  const nSeqs = multilengths.length;
  const order = identityOrder(multilengths);
  const orient = identityOrient(multilengths);
  for (let g = 0; g < nSeqs; g++) {
    if (g === ref) continue;
    const nCtg = multilengths[g].length;
    const scores = Array.from({ length: nCtg }, () => []);
    const votes = Array(nCtg).fill(0);
    for (let bi = 0; bi < left.length; bi++) {
      const aRef = left[bi][ref];
      const bRef = right[bi][ref];
      const aG = left[bi][g];
      const bG = right[bi][g];
      if (aRef < 0 || bRef < 0 || aG < 0 || bG < 0) continue;
      const midRef = (aRef + bRef) / 2;
      const [cidG] = globalToContig(Math.floor((aG + bG) / 2), multilengths[g]);
      if (cidG < 0 || cidG >= nCtg) continue;
      scores[cidG].push(midRef);
      votes[cidG] += strands[bi][g] === strands[bi][ref] ? 1 : -1;
    }
    const medians = scores.map((arr, cid) => {
      if (!arr.length) return [Infinity, cid];
      const s = arr.slice().sort((a, b) => a - b);
      const m = s[Math.floor(s.length / 2)];
      return [m, cid];
    });
    medians.sort((a, b) => a[0] - b[0]);
    const newOrder = medians.filter(([m]) => m !== Infinity).map(([, cid]) => cid);
    const seen = new Set(newOrder);
    for (let cid = 0; cid < nCtg; cid++) if (!seen.has(cid)) newOrder.push(cid);
    order[g] = newOrder;
    orient[g] = votes.map((v) => (v < 0 ? -1 : 1));
  }
  return { order, orient };
}
