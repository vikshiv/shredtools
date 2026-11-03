/**
 * View culling + optional stride for collinear blocks.
 *
 * Overview must still draw essentially all blocks. Tiny genomic blocks are
 * kept and given a minimum screen width later in the canvas — we must NOT
 * drop them when zoomed out (that left the plot empty until zoom-in).
 */

function absent(a, b) {
  // Overview uses -1; extract-justified flanks use null (coords may be < 0)
  return a == null || b == null || a === -1 || b === -1;
}

function overlapsView(left, right, genomeOrder, x0, x1) {
  for (const g of genomeOrder) {
    let a = left[g];
    let b = right[g];
    if (absent(a, b)) continue;
    if (a > b) [a, b] = [b, a];
    if (a < x1 && b > x0) return true;
  }
  return false;
}

function widthStats(left, right, genomeOrder) {
  let minW = Infinity;
  let maxW = 0;
  let n = 0;
  for (const g of genomeOrder) {
    let a = left[g];
    let b = right[g];
    if (absent(a, b)) continue;
    const w = Math.abs(b - a);
    minW = Math.min(minW, w);
    maxW = Math.max(maxW, w);
    n++;
  }
  return { minW: n ? minW : 0, maxW, n };
}

/**
 * @returns {{left, right, strands, nIn, nOut, stride}}
 */
export function mergeBlocksForLod(
  left,
  right,
  strands,
  genomeOrder,
  { x0, x1, maxDraw = 20000, maxWidthRatio = 40 } = {}
) {
  const nIn = left.length;
  if (nIn === 0) return { left: [], right: [], strands: [], nIn: 0, nOut: 0, stride: 1 };

  const visible = [];
  for (let i = 0; i < nIn; i++) {
    if (!overlapsView(left[i], right[i], genomeOrder, x0, x1)) continue;
    const { minW, maxW, n } = widthStats(left[i], right[i], genomeOrder);
    if (n < 2) continue;
    if (minW > 0 && maxW / minW > maxWidthRatio) continue;
    visible.push(i);
  }

  // Prefer drawing everything; stride only if truly huge
  const stride = Math.max(1, Math.ceil(visible.length / maxDraw));
  const outL = [];
  const outR = [];
  const outS = [];
  for (let k = 0; k < visible.length; k += stride) {
    const i = visible[k];
    outL.push(left[i]);
    outR.push(right[i]);
    outS.push(strands[i]);
  }
  return { left: outL, right: outR, strands: outS, nIn, nOut: outL.length, stride };
}

/**
 * Expand each block's interval to at least ``minPx`` on screen so overview
 * ribbons stay visible without inventing merged mega-blocks.
 * Returns new left/right arrays (strands reused).
 */
export function ensureMinScreenWidth(left, right, genomeOrder, pxPerBp, minPx = 1.0) {
  const minBp = minPx / Math.max(pxPerBp, 1e-12);
  const outL = new Array(left.length);
  const outR = new Array(right.length);
  for (let i = 0; i < left.length; i++) {
    const ls = left[i].slice();
    const rs = right[i].slice();
    for (const g of genomeOrder) {
      let a = ls[g];
      let b = rs[g];
      if (absent(a, b)) continue;
      if (a > b) [a, b] = [b, a];
      if (b - a < minBp) {
        const mid = (a + b) / 2;
        a = mid - minBp / 2;
        b = mid + minBp / 2;
      }
      ls[g] = a;
      rs[g] = b;
    }
    outL[i] = ls;
    outR[i] = rs;
  }
  return { left: outL, right: outR };
}
