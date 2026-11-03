/** Native pan/zoom canvas for mumemto-style synteny ribbons. */

import { buildBlockPolygons } from "./overview_geom.mjs?v=8";
import { mergeBlocksForLod } from "./lod.mjs?v=8";

const AXIS_H_CONCAT = 40;
const AXIS_H_SPACED = 48;
const PAD_L = 6;
const PAD_R = 10;
const PAD_T = 6;

/** Compact label when step is unknown (LOD / tips). */
function formatBp(x) {
  const v = Math.abs(x);
  if (v >= 1e6) return (x / 1e6).toFixed(v >= 1e7 ? 0 : 1) + " Mb";
  if (v >= 1e3) return (x / 1e3).toFixed(v >= 1e4 ? 0 : 1) + " kb";
  return Math.round(x) + " bp";
}

/**
 * Tick labels that resolve neighboring ticks. At ~57 Mb with a 5 kb step,
 * plain formatBp collapses every label to "57 Mb".
 */
function formatBpTick(x, step) {
  const s = Math.max(1, Number(step) || 1);
  const v = Math.abs(x);
  if (s >= 1e6) {
    return (x / 1e6).toFixed(s >= 5e6 ? 0 : 1) + " Mb";
  }
  if (v >= 1e6) {
    const dec = Math.min(6, Math.max(1, Math.ceil(-Math.log10(s / 1e6))));
    return (x / 1e6).toFixed(dec) + " Mb";
  }
  if (s >= 1e3 || v >= 1e3) {
    const unit = 1e3;
    const dec = Math.min(3, Math.max(0, Math.ceil(-Math.log10(s / unit))));
    return (x / unit).toFixed(dec) + " kb";
  }
  return Math.round(x) + " bp";
}

function niceTicks(x0, x1, target = 6) {
  const span = Math.max(1, x1 - x0);
  const raw = span / target;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / pow;
  let step;
  if (n < 1.5) step = pow;
  else if (n < 3.5) step = 2 * pow;
  else if (n < 7.5) step = 5 * pow;
  else step = 10 * pow;
  const start = Math.ceil(x0 / step) * step;
  const ticks = [];
  for (let t = start; t <= x1 + 1e-6; t += step) ticks.push(t);
  return { ticks, step };
}

export class NativeSyntenyView {
  constructor(canvas, wrap, { onViewChange, onExtractSelect, onSelectDrag } = {}) {
    this.canvas = canvas;
    this.wrap = wrap;
    this.ctx = canvas.getContext("2d", { alpha: false });
    this.onViewChange = onViewChange;
    this.onExtractSelect = onExtractSelect;
    this.onSelectDrag = onSelectDrag;

    this.left = [];
    this.right = [];
    this.strands = [];
    this.genomeOrder = [];
    this.seqLengths = [];
    this.centering = [];
    /** @type {Array<Array<{cid,name,start,end,length,local0?,isContigStart?,isContigEnd?}>>} */
    this.segmentsByGenome = [];
    this.mode = "concat"; // concat | spaced
    /** @type {"overview"|"extract"} */
    this.viewKind = "overview";
    /** Query-genome axis for extract: contig-local labels (shredtools extract style). */
    this.queryAxis = null; // { seqIdx, contigName, local0 } | null
    /** Requested extract region; yellow markers while set. endsByGenome = jagged right. */
    this.extractRegion = null; // { x0, x1, endsByGenome?: (number|null)[] } | null
    this.rowH = 16;
    this.minRowH = 6;
    this.maxRowH = 36;

    this.viewX0 = 0;
    this.viewX1 = 1;
    this.dataXMin = 0;
    this.dataXMax = 1;

    this._drag = null;
    this._raf = 0;
    this._bind();
  }

  get axisH() {
    return this.mode === "spaced" ? AXIS_H_SPACED : AXIS_H_CONCAT;
  }

  _bind() {
    const w = this.wrap;
    w.addEventListener("wheel", (e) => this._onWheel(e), { passive: false });
    w.addEventListener("mousedown", (e) => this._onDown(e));
    window.addEventListener("mousemove", (e) => this._onMove(e));
    window.addEventListener("mouseup", () => this._onUp());
    window.addEventListener("resize", () => this.scheduleDraw());
  }

  setData({
    left,
    right,
    strands,
    genomeOrder,
    seqLengths,
    centering,
    dataXMin,
    dataXMax,
    segmentsByGenome,
    mode,
    viewKind,
    queryAxis,
    extractRegion,
    focusRegion,
  }) {
    this.left = left;
    this.right = right;
    this.strands = strands;
    this.genomeOrder = genomeOrder.slice();
    this.seqLengths = seqLengths.slice();
    this.centering = centering ? centering.slice() : genomeOrder.map(() => 0);
    this.segmentsByGenome = segmentsByGenome || [];
    this.mode = mode || "concat";
    this.viewKind = viewKind || "overview";
    this.queryAxis = queryAxis || null;
    this.extractRegion = extractRegion || null;
    this.dataXMin = dataXMin != null ? dataXMin : 0;
    this.dataXMax = Math.max(this.dataXMin + 1, dataXMax || 1);
    if (focusRegion && focusRegion.x1 > focusRegion.x0) {
      this.viewX0 = focusRegion.x0;
      this.viewX1 = focusRegion.x1;
    } else if (!(this.viewX1 > this.viewX0)) {
      this.viewX0 = this.dataXMin;
      this.viewX1 = this.dataXMax;
    }
    this.scheduleDraw();
  }

  resetView() {
    if (this.extractRegion && this.extractRegion.x1 > this.extractRegion.x0) {
      this.viewX0 = this.extractRegion.x0;
      this.viewX1 = this.extractRegion.x1;
    } else {
      this.viewX0 = this.dataXMin;
      this.viewX1 = this.dataXMax;
    }
    this.rowH = 16;
    this.scheduleDraw();
    this._emit();
  }

  setView(x0, x1) {
    const spanAvail = Math.max(1, this.dataXMax - this.dataXMin);
    const minSpan = Math.max(50, spanAvail / 1e6);
    let a = Math.min(x0, x1);
    let b = Math.max(x0, x1);
    if (b - a < minSpan) b = a + minSpan;
    this.viewX0 = a;
    this.viewX1 = b;
    this.scheduleDraw();
    this._emit();
  }

  /** Focus the yellow extract region without clearing it. */
  focusExtractRegion() {
    if (!this.extractRegion) return;
    this.setView(this.extractRegion.x0, this.extractRegion.x1);
  }

  /**
   * Change display order without resetting the view or reloading geometry.
   * Next draw() rebuilds ribbons via buildBlockPolygons(genomeOrder): each
   * adjacent pair that disagrees on strand gets a twisted red edge.
   */
  setGenomeOrder(genomeOrder) {
    this.genomeOrder = genomeOrder.slice();
    // Centering is per display row; keep length in sync (values stay 0 here)
    this.centering = this.genomeOrder.map((_, i) => this.centering[i] || 0);
    this.scheduleDraw();
    this._emit();
  }

  _emit() {
    if (this.onViewChange) {
      this.onViewChange({
        x0: this.viewX0,
        x1: this.viewX1,
        rowH: this.rowH,
        nGenomes: this.genomeOrder.length,
        mode: this.mode,
      });
    }
  }

  scheduleDraw() {
    if (this._raf) return;
    this._raf = requestAnimationFrame(() => {
      this._raf = 0;
      this.draw();
    });
  }

  _layout() {
    const dpr = window.devicePixelRatio || 1;
    const cssW = Math.max(320, this.wrap.clientWidth);
    const n = Math.max(1, this.genomeOrder.length);
    const plotH = PAD_T + n * this.rowH;
    const cssH = plotH + this.axisH;
    this.canvas.style.width = cssW + "px";
    this.canvas.style.height = cssH + "px";
    this.canvas.width = Math.floor(cssW * dpr);
    this.canvas.height = Math.floor(cssH * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.cssW = cssW;
    this.cssH = cssH;
    this.plotW = cssW - PAD_L - PAD_R;
    this.plotH = plotH;
    return { cssW, cssH, plotH };
  }

  xToPx(x) {
    return PAD_L + ((x - this.viewX0) / (this.viewX1 - this.viewX0)) * this.plotW;
  }

  pxToX(px) {
    return this.viewX0 + ((px - PAD_L) / this.plotW) * (this.viewX1 - this.viewX0);
  }

  /** Map display x → contig-local label using a reference genome's segments. */
  displayToContigLocal(x, genomeId) {
    const segs = this.segmentsByGenome[genomeId] || [];
    for (const s of segs) {
      if (x >= s.start && x < s.end) {
        return { name: s.name, local: x - s.start, inContig: true, seg: s };
      }
    }
    // in spacer / outside
    for (let i = 0; i < segs.length - 1; i++) {
      if (x >= segs[i].end && x < segs[i + 1].start) {
        return { name: null, local: null, inContig: false, between: [segs[i].name, segs[i + 1].name] };
      }
    }
    return { name: null, local: null, inContig: false };
  }

  draw() {
    this._layout();
    const ctx = this.ctx;
    const { cssW, cssH, plotH } = this;
    const n = this.genomeOrder.length;

    ctx.fillStyle = "#f3efe6";
    ctx.fillRect(0, 0, cssW, cssH);

    this._drawTracksAndContigs(ctx, cssW, n);

    const pxPerBp = this.plotW / Math.max(1, this.viewX1 - this.viewX0);
    const lod = mergeBlocksForLod(this.left, this.right, this.strands, this.genomeOrder, {
      x0: this.viewX0,
      x1: this.viewX1,
      maxDraw: 25000,
    });
    // Match mumemto viz defaults for collinear blocks: alpha=0.8, linewidth=0,
    // color=#00A2FF — real genomic widths, no min-px padding, no stroke.
    const draw = buildBlockPolygons(
      lod.left,
      lod.right,
      lod.strands,
      this.genomeOrder,
      this.centering,
      {
        color: "#00A2FF",
        invColor: "#FF0000",
        maxBlocks: 50000,
        // Extract: color-only like viz_mums — twisting collinear blocks paints
        // huge false X-ribbons across the whole window.
        twistInversions: this.viewKind !== "extract",
      }
    );

    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, cssW, plotH);
    ctx.clip();

    ctx.globalAlpha = 0.8;
    for (let i = 0; i < draw.polygons.length; i++) {
      const poly = draw.polygons[i];
      if (!poly || poly.length < 3) continue;
      ctx.beginPath();
      for (let k = 0; k < poly.length; k++) {
        const [x, yi] = poly[k];
        const px = this.xToPx(x);
        const py = PAD_T + yi * this.rowH + this.rowH / 2;
        if (k === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.closePath();
      ctx.fillStyle = draw.colors[i] || "#00A2FF";
      ctx.fill();
    }
    ctx.globalAlpha = 1;
    ctx.restore();

    // Contig cascades are overview-layout (index-aligned). In extract they zig-zag
    // across unrelated contigs and look like false inversions — skip them.
    if (this.viewKind !== "extract") {
      this._drawContigBoundaryRibbons(ctx, n, plotH);
    }
    // Yellow request bounds (persist while panning; cleared on reset / new zoom)
    this._drawExtractRegionMarkers(ctx, n, plotH);
    this._drawSelection(ctx);

    this._drawAxis(ctx, cssW, plotH);
    this._lodInfo = {
      nIn: lod.nIn,
      nOut: lod.nOut,
      stride: lod.stride || 1,
      pxPerBp,
      nPolys: draw.polygons.length,
      minPx: 0,
    };
  }

  _drawTracksAndContigs(ctx, cssW, n) {
    for (let i = 0; i < n; i++) {
      const y = PAD_T + i * this.rowH;
      ctx.fillStyle = i % 2 === 0 ? "rgba(0,0,0,0.035)" : "transparent";
      ctx.fillRect(0, y, cssW, this.rowH);

      const g = this.genomeOrder[i];
      const segs = this.segmentsByGenome[g] || [];
      const mid = y + this.rowH / 2;

      if (this.mode === "spaced") {
        // Shade spacer gaps
        for (let s = 0; s < segs.length - 1; s++) {
          const gap0 = segs[s].end;
          const gap1 = segs[s + 1].start;
          if (gap1 <= gap0) continue;
          const xA = this.xToPx(gap0);
          const xB = this.xToPx(gap1);
          if (xB < PAD_L || xA > cssW - PAD_R) continue;
          ctx.fillStyle = "rgba(180, 160, 120, 0.18)";
          ctx.fillRect(Math.max(PAD_L, xA), y, Math.min(cssW - PAD_R, xB) - Math.max(PAD_L, xA), this.rowH);
        }
        // Contig backbones as separate segments
        ctx.strokeStyle = "#6a6358";
        ctx.lineWidth = 1.25;
        ctx.lineCap = "butt";
        for (const seg of segs) {
          const xA = this.xToPx(seg.start);
          const xB = this.xToPx(seg.end);
          if (xB < PAD_L || xA > cssW - PAD_R) continue;
          ctx.beginPath();
          ctx.moveTo(Math.max(PAD_L, xA), mid);
          ctx.lineTo(Math.min(cssW - PAD_R, xB), mid);
          ctx.stroke();
        }
      } else {
        // Concat: continuous backbone (prefer segment extent — works for justified extract)
        let xA;
        let xB;
        if (segs.length) {
          xA = this.xToPx(segs[0].start);
          xB = this.xToPx(segs[segs.length - 1].end);
        } else {
          const c = this.centering[i] || 0;
          const len = this.seqLengths[g] || 0;
          xA = this.xToPx(c);
          xB = this.xToPx(c + len);
        }
        ctx.strokeStyle = "#cfc7b8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(Math.max(PAD_L, xA), mid);
        ctx.lineTo(Math.min(cssW - PAD_R, xB), mid);
        ctx.stroke();
      }
    }
  }

  /**
   * Yellow request markers: left is justified (straight at x0), right follows
   * each genome's extract window end (jagged cascade along boundary MUMs).
   */
  _drawExtractRegionMarkers(ctx, n, plotH) {
    const reg = this.extractRegion;
    if (!reg || !(reg.x1 > reg.x0)) return;
    const lw = Math.max(2.75, Math.min(4, this.rowH * 0.28));
    const ends = reg.endsByGenome || null;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, this.cssW, plotH);
    ctx.clip();

    // Left: all genomes share justified start → vertical cascade
    this._strokeCascade(ctx, n, () => reg.x0, "#E6B800", lw);
    // Right: per-genome window end (w1 - w0) → jagged along the boundary
    this._strokeCascade(
      ctx,
      n,
      (_segs, displayRow) => {
        const g = this.genomeOrder[displayRow];
        if (ends && ends[g] != null) return ends[g];
        return reg.x1;
      },
      "#E6B800",
      lw
    );

    // Per-row wash outside the request (respects jagged right edge)
    ctx.fillStyle = "rgba(230, 184, 0, 0.06)";
    for (let i = 0; i < n; i++) {
      const g = this.genomeOrder[i];
      const y = PAD_T + i * this.rowH;
      const endX = ends && ends[g] != null ? ends[g] : reg.x1;
      const xA = this.xToPx(reg.x0);
      const xB = this.xToPx(endX);
      if (xA > PAD_L) {
        ctx.fillRect(PAD_L, y, Math.max(0, xA - PAD_L), this.rowH);
      }
      if (xB < this.cssW - PAD_R) {
        ctx.fillRect(xB, y, Math.max(0, this.cssW - PAD_R - xB), this.rowH);
      }
    }
    ctx.restore();
  }

  /**
   * Shred-style cascading boundary: a thick 1D polyline through each genome row
   * at corresponding contig edges (not per-row ticks). Concat = red, spaced = black.
   * Extract: only true contig ends (isContigStart / isContigEnd) participate.
   */
  _drawContigBoundaryRibbons(ctx, n, plotH) {
    let maxSegs = 0;
    for (let i = 0; i < n; i++) {
      const g = this.genomeOrder[i];
      maxSegs = Math.max(maxSegs, (this.segmentsByGenome[g] || []).length);
    }
    if (maxSegs < 1) return;

    const color = this.mode === "concat" ? "#C62828" : "#1c1b19";
    const lw = Math.max(2.25, Math.min(3.5, this.rowH * 0.2));

    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, this.cssW, plotH);
    ctx.clip();

    // One cascade per contig-start column (aligns by display order index)
    for (let si = 0; si < maxSegs; si++) {
      this._strokeCascade(ctx, n, (segs) => {
        const seg = segs[si];
        if (!seg || !seg.isContigStart) return null;
        return seg.start;
      }, color, lw);
    }
    // Contig ends: needed for spaced (gap sides) and for the final genome end in concat.
    // Skip ends that coincide with the next contig's start (concat internal duplicates).
    for (let si = 0; si < maxSegs; si++) {
      this._strokeCascade(ctx, n, (segs) => {
        const seg = segs[si];
        if (!seg || !seg.isContigEnd) return null;
        const next = segs[si + 1];
        if (next && next.isContigStart && Math.abs(next.start - seg.end) < 0.5) return null;
        return seg.end;
      }, color, lw);
    }

    ctx.restore();
    ctx.lineWidth = 1;
  }

  /** Polyline through genome rows at getX(segs, displayRow); breaks where null. */
  _strokeCascade(ctx, n, getX, color, lw) {
    const pts = [];
    for (let i = 0; i < n; i++) {
      const g = this.genomeOrder[i];
      const x = getX(this.segmentsByGenome[g] || [], i);
      if (x == null || !Number.isFinite(x)) {
        pts.push(null);
        continue;
      }
      pts.push({
        px: this.xToPx(x),
        py: PAD_T + i * this.rowH + this.rowH / 2,
      });
    }
    const live = pts.filter((p) => p);
    if (live.length < 2) return;

    ctx.strokeStyle = color;
    ctx.lineWidth = lw;
    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    ctx.globalAlpha = 0.92;

    let started = false;
    ctx.beginPath();
    for (let i = 0; i < pts.length; i++) {
      const p = pts[i];
      if (!p) {
        if (started) {
          ctx.stroke();
          ctx.beginPath();
          started = false;
        }
        continue;
      }
      if (!started) {
        ctx.moveTo(p.px, p.py);
        started = true;
      } else {
        ctx.lineTo(p.px, p.py);
      }
    }
    if (started) ctx.stroke();
    ctx.globalAlpha = 1;
  }

  _drawAxis(ctx, cssW, plotH) {
    const y0 = plotH;
    const axisH = this.axisH;
    ctx.fillStyle = "#fffdf8";
    ctx.fillRect(0, y0, cssW, axisH);
    ctx.strokeStyle = "#d9d2c5";
    ctx.beginPath();
    ctx.moveTo(0, y0 + 0.5);
    ctx.lineTo(cssW, y0 + 0.5);
    ctx.stroke();

    ctx.font = '11px "IBM Plex Mono", ui-monospace, monospace';
    ctx.textBaseline = "top";

    if (this.viewKind === "extract" && this.queryAxis) {
      this._drawExtractAxis(ctx, cssW, y0);
      return;
    }

    const refG = this.genomeOrder[0];
    const segs = this.segmentsByGenome[refG] || [];

    if (this.mode === "spaced") {
      this._drawSpacedAxis(ctx, cssW, y0, segs);
    } else {
      this._drawConcatAxis(ctx, cssW, y0, segs);
    }
  }

  /**
   * Extract: single query contig-local axis (shredtools style).
   * Do not walk full-assembly shifted segments — those paint every chromosome's
   * ticks into the strip (e.g. repeating "13 Mb") and collide with the scale bar.
   */
  _drawExtractAxis(ctx, cssW, y0) {
    const q = this.queryAxis;
    const local0 = q.local0 || 0;
    ctx.textAlign = "center";
    const { ticks, step } = niceTicks(this.viewX0 + local0, this.viewX1 + local0, 7);
    ctx.fillStyle = "#5c574e";
    ctx.font = '10px "IBM Plex Mono", ui-monospace, monospace';
    for (const loc of ticks) {
      const x = loc - local0;
      const px = this.xToPx(x);
      if (px < PAD_L - 2 || px > cssW - PAD_R + 2) continue;
      ctx.strokeStyle = "#b9b1a2";
      ctx.beginPath();
      ctx.moveTo(px, y0);
      ctx.lineTo(px, y0 + 5);
      ctx.stroke();
      ctx.fillText(formatBpTick(loc, step), px, y0 + 7);
    }
    if (q.contigName) {
      ctx.fillStyle = "#2a4f9e";
      ctx.font = '9px "IBM Plex Sans", sans-serif';
      ctx.textAlign = "left";
      ctx.fillText(shortName(q.contigName), PAD_L, y0 + 22);
      ctx.textAlign = "center";
    }
    this._drawScaleBar(ctx, cssW, y0, step);
  }

  _drawConcatAxis(ctx, cssW, y0, segs) {
    // Continuous global (concat) coordinates
    const { ticks, step } = niceTicks(this.viewX0, this.viewX1, 7);
    ctx.fillStyle = "#5c574e";
    ctx.textAlign = "center";
    for (const t of ticks) {
      const px = this.xToPx(t);
      if (px < PAD_L - 2 || px > cssW - PAD_R + 2) continue;
      ctx.strokeStyle = "#b9b1a2";
      ctx.beginPath();
      ctx.moveTo(px, y0);
      ctx.lineTo(px, y0 + 6);
      ctx.stroke();
      ctx.fillText(formatBp(t), px, y0 + 8);
    }

    // Contig names only (boundaries are the thick row lines)
    for (let i = 0; i < segs.length; i++) {
      const mid = (segs[i].start + segs[i].end) / 2;
      const midPx = this.xToPx(mid);
      const wPx = (segs[i].end - segs[i].start) * (this.plotW / Math.max(1, this.viewX1 - this.viewX0));
      if (wPx > 48 && midPx > PAD_L && midPx < cssW - PAD_R) {
        ctx.fillStyle = "#8a3b3b";
        ctx.font = '9px "IBM Plex Mono", ui-monospace, monospace';
        ctx.fillText(shortName(segs[i].name), midPx, y0 + 20);
        ctx.font = '11px "IBM Plex Mono", ui-monospace, monospace';
      }
    }

    this._drawScaleBar(ctx, cssW, y0, step);
  }

  _drawSpacedAxis(ctx, cssW, y0, segs) {
    // Per-contig local coordinates on the reference genome; blank in spacers
    ctx.fillStyle = "#5c574e";
    ctx.textAlign = "center";

    for (const seg of segs) {
      const xA = this.xToPx(seg.start);
      const xB = this.xToPx(seg.end);
      if (xB < PAD_L || xA > cssW - PAD_R) continue;

      const clipA = Math.max(PAD_L, xA);
      const clipB = Math.min(cssW - PAD_R, xB);
      ctx.fillStyle = "rgba(47, 111, 237, 0.06)";
      ctx.fillRect(clipA, y0 + 1, Math.max(0, clipB - clipA), this.axisH - 2);

      const visLocal0 = Math.max(0, this.viewX0 - seg.start);
      const visLocal1 = Math.min(seg.length, this.viewX1 - seg.start);
      if (visLocal1 <= visLocal0) continue;
      const { ticks, step } = niceTicks(visLocal0, visLocal1, Math.max(2, Math.floor((clipB - clipA) / 80)));
      ctx.fillStyle = "#5c574e";
      ctx.font = '10px "IBM Plex Mono", ui-monospace, monospace';
      for (const loc of ticks) {
        if (loc < 0 || loc > seg.length) continue;
        const px = this.xToPx(seg.start + loc);
        if (px < PAD_L - 2 || px > cssW - PAD_R + 2) continue;
        ctx.strokeStyle = "#b9b1a2";
        ctx.beginPath();
        ctx.moveTo(px, y0);
        ctx.lineTo(px, y0 + 5);
        ctx.stroke();
        ctx.fillText(formatBp(loc), px, y0 + 7);
      }

      const midPx = (clipA + clipB) / 2;
      if (clipB - clipA > 36) {
        ctx.fillStyle = "#2a4f9e";
        ctx.font = '9px "IBM Plex Sans", sans-serif';
        ctx.fillText(shortName(seg.name), midPx, y0 + 22);
      }
    }

    // shade spacer regions on axis
    for (let i = 0; i < segs.length - 1; i++) {
      const g0 = segs[i].end;
      const g1 = segs[i + 1].start;
      const xA = this.xToPx(g0);
      const xB = this.xToPx(g1);
      if (xB <= xA) continue;
      const clipA = Math.max(PAD_L, xA);
      const clipB = Math.min(cssW - PAD_R, xB);
      if (clipB <= clipA) continue;
      ctx.fillStyle = "rgba(180, 160, 120, 0.22)";
      ctx.fillRect(clipA, y0 + 1, clipB - clipA, this.axisH - 2);
      if (clipB - clipA > 28) {
        ctx.fillStyle = "#8a7d66";
        ctx.font = '8px "IBM Plex Sans", sans-serif';
        ctx.fillText("gap", (clipA + clipB) / 2, y0 + 14);
      }
    }

    const { step } = niceTicks(0, Math.max(1, (segs[0] && segs[0].length) || 1e6), 5);
    this._drawScaleBar(ctx, cssW, y0, step);
  }

  _drawScaleBar(ctx, cssW, y0, step) {
    const barBp = step;
    const barPx = barBp * (this.plotW / Math.max(1, this.viewX1 - this.viewX0));
    // Sit in the bottom of the axis band (below tick labels / contig name).
    if (barPx > 28 && barPx < this.plotW * 0.35) {
      const bx = cssW - PAD_R - barPx;
      const by = y0 + this.axisH - 3;
      ctx.strokeStyle = "#1c1b19";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(bx, by);
      ctx.lineTo(bx + barPx, by);
      ctx.stroke();
      ctx.lineWidth = 1;
      ctx.textAlign = "right";
      ctx.textBaseline = "bottom";
      ctx.fillStyle = "#1c1b19";
      ctx.font = '10px "IBM Plex Mono", ui-monospace, monospace';
      ctx.fillText(formatBpTick(barBp, barBp), bx + barPx, by - 3);
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
    }
  }

  _local(e) {
    const rect = this.canvas.getBoundingClientRect();
    return { x: e.clientX - rect.left, y: e.clientY - rect.top };
  }

  /** Keep extract-select X inside the plot so pxToX cannot invent off-genome coords. */
  _clampPlotX(px) {
    const lo = PAD_L;
    const hi = Math.max(lo, (this.cssW || 0) - PAD_R);
    if (!Number.isFinite(px)) return lo;
    return Math.min(hi, Math.max(lo, px));
  }

  _rowFromY(y) {
    const r = Math.floor((y - PAD_T) / this.rowH);
    if (r < 0 || r >= this.genomeOrder.length) return -1;
    return r;
  }

  _onWheel(e) {
    // Plain wheel: let #syntenyCanvasWrap scroll vertically (hundreds of genomes).
    // Shift+wheel: row height only.
    if (!e.shiftKey) return;
    e.preventDefault();
    const factor = e.deltaY < 0 ? 1.08 : 1 / 1.08;
    this.rowH = Math.min(this.maxRowH, Math.max(this.minRowH, this.rowH * factor));
    this.scheduleDraw();
    this._emit();
  }

  _onDown(e) {
    if (e.button !== 0) return;
    const p = this._local(e);
    const row = this._rowFromY(p.y);

    // Alt/Meta+drag = pan. Plain drag = region select (zoom / nested extract).
    if (e.altKey || e.metaKey) {
      this.wrap.classList.add("panning");
      this._drag = {
        kind: "pan",
        x: e.clientX,
        viewX0: this.viewX0,
        viewX1: this.viewX1,
      };
      return;
    }

    if (row < 0) return;
    this.wrap.style.cursor = "col-resize";
    const x = this._clampPlotX(p.x);
    this._drag = {
      kind: "select",
      row,
      x0: x,
      x1: x,
      clientX: e.clientX,
      clientY: e.clientY,
    };
    this._emitSelectDrag(e);
    this.scheduleDraw();
  }

  _onMove(e) {
    if (!this._drag) return;
    if (this._drag.kind === "pan") {
      const dx = e.clientX - this._drag.x;
      const span = this._drag.viewX1 - this._drag.viewX0;
      const dGen = (-dx / Math.max(1, this.plotW)) * span;
      this.viewX0 = this._drag.viewX0 + dGen;
      this.viewX1 = this._drag.viewX1 + dGen;
      this.scheduleDraw();
      return;
    }
    if (this._drag.kind === "select") {
      const p = this._local(e);
      this._drag.x1 = this._clampPlotX(p.x);
      this._drag.clientX = e.clientX;
      this._drag.clientY = e.clientY;
      this._emitSelectDrag(e);
      this.scheduleDraw();
    }
  }

  _emitSelectDrag(e) {
    if (!this.onSelectDrag || !this._drag || this._drag.kind !== "select") return;
    const a = Math.min(this._drag.x0, this._drag.x1);
    const b = Math.max(this._drag.x0, this._drag.x1);
    const genomeIdx = this.genomeOrder[this._drag.row];
    this.onSelectDrag({
      displayRow: this._drag.row,
      genomeIdx,
      xStart: this.pxToX(a),
      xEnd: this.pxToX(b),
      clientX: e?.clientX ?? this._drag.clientX,
      clientY: e?.clientY ?? this._drag.clientY,
      active: true,
    });
  }

  _onUp(e) {
    if (!this._drag) return;
    const drag = this._drag;
    this._drag = null;
    this.wrap.classList.remove("panning");
    this.wrap.style.cursor = "crosshair";

    if (drag.kind === "pan") {
      this._emit();
      return;
    }
    if (drag.kind === "select") {
      if (this.onSelectDrag) this.onSelectDrag({ active: false });
      const a = Math.min(drag.x0, drag.x1);
      const b = Math.max(drag.x0, drag.x1);
      this.scheduleDraw();
      if (b - a < 6) return;
      const xStart = this.pxToX(a);
      const xEnd = this.pxToX(b);
      const genomeIdx = this.genomeOrder[drag.row];
      if (this.onExtractSelect) {
        this.onExtractSelect({
          displayRow: drag.row,
          genomeIdx,
          xStart,
          xEnd,
        });
      }
    }
  }

  _drawSelection(ctx) {
    if (!this._drag || this._drag.kind !== "select") return;
    const a = Math.min(this._drag.x0, this._drag.x1);
    const b = Math.max(this._drag.x0, this._drag.x1);
    const y = PAD_T + this._drag.row * this.rowH;
    ctx.fillStyle = "rgba(47, 111, 237, 0.22)";
    ctx.strokeStyle = "rgba(47, 111, 237, 0.9)";
    ctx.lineWidth = 1;
    ctx.fillRect(a, y, b - a, this.rowH);
    ctx.strokeRect(a + 0.5, y + 0.5, b - a - 1, this.rowH - 1);
  }

  getLodStats() {
    return this._lodInfo || null;
  }
}

function shortName(name) {
  if (!name) return "";
  // Prefer last accession-like token; keep short
  const s = String(name);
  if (s.length <= 14) return s;
  return s.slice(0, 12) + "…";
}

export { PAD_T, formatBp };
