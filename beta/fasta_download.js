/**
 * Client-side multi-FASTA download from indexed HPRC Release 2 assemblies.
 * Fetches subsequences via HTTP Range on .fa.gz + .fai + .gzi (BGZF).
 */

const ASSEMBLY_INDEX_URL =
  "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv";

const S3_HTTPS_PREFIX = "https://human-pangenomics.s3.us-west-2.amazonaws.com/";

const HPRC_R2_PRESETS = new Set(["hprcv2_enhanced", "hprcv2_merged"]);

const SCRATCH_FILE = "shredtools_extract.fa";
const SCRATCH_FILE_GZ = "shredtools_extract.fa.gz";

/** @type {Map<string, {fa: string, fai: string, gzi: string}> | null} */
let assemblyIndex = null;

/** @type {Map<string, string>} */
const faiCache = new Map();

/** @type {Map<string, Array<{compressed: number, uncompressed: number}>>} */
const gziCache = new Map();

/** @type {{ fetchBytes: (url: string, start: number|null, end: number|null) => Promise<Uint8Array>, fetchText: (url: string) => Promise<string> } | null} */
let httpClient = null;

/** @type {((blockBytes: Uint8Array) => Promise<Uint8Array>) | null} */
let bgzfDecompressor = null;

const MAX_INFLIGHT_FETCHES = 6;
const FETCH_RETRIES = 3;
let inflightFetches = 0;
const fetchWaiters = [];

/** @type {AbortSignal | null} */
let activeBuildAbortSignal = null;

function checkBuildAborted() {
  if (activeBuildAbortSignal?.aborted) {
    throw new DOMException("FASTA build cancelled", "AbortError");
  }
}

export function isFastaBuildCancelled(err) {
  return err?.name === "AbortError";
}

function acquireFetchSlot() {
  if (inflightFetches < MAX_INFLIGHT_FETCHES) {
    inflightFetches++;
    return Promise.resolve();
  }
  return new Promise((resolve) => fetchWaiters.push(resolve));
}

function releaseFetchSlot() {
  inflightFetches--;
  const next = fetchWaiters.shift();
  if (next) {
    inflightFetches++;
    next();
  }
}

async function withFetchSlot(fn) {
  await acquireFetchSlot();
  try {
    return await fn();
  } finally {
    releaseFetchSlot();
  }
}

function isRetryableFetchError(err) {
  const msg = String(err?.message ?? err ?? "");
  return (
    err instanceof TypeError ||
    msg.includes("Failed to fetch") ||
    msg.includes("NetworkError") ||
    msg.includes("Load failed")
  );
}

async function withRetries(label, fn) {
  let lastErr;
  for (let attempt = 0; attempt < FETCH_RETRIES; attempt++) {
    try {
      checkBuildAborted();
      return await fn();
    } catch (err) {
      if (err?.name === "AbortError") throw err;
      lastErr = err;
      if (!isRetryableFetchError(err) || attempt + 1 >= FETCH_RETRIES) break;
      await new Promise((r) => setTimeout(r, 250 * (attempt + 1)));
    }
  }
  const detail = lastErr?.message ?? String(lastErr);
  throw new Error(`${label}: ${detail}`);
}

export function setHttpClient(client) {
  httpClient = client;
}

export function setBgzfDecompressor(fn) {
  bgzfDecompressor = fn;
}

export function isFastaDownloadSupported(pangenomeKey) {
  return HPRC_R2_PRESETS.has(String(pangenomeKey ?? ""));
}

function s3ToHttps(url) {
  const s = String(url ?? "").trim();
  if (s.startsWith("s3://human-pangenomics/")) {
    return S3_HTTPS_PREFIX + s.slice("s3://human-pangenomics/".length);
  }
  if (s.startsWith("https://")) return s;
  throw new Error(`Unsupported assembly URL: ${s}`);
}

function parseContigKey(contig) {
  const parts = String(contig ?? "").split("#");
  if (parts.length < 2) return null;
  return { sampleId: parts[0], hap: parts[1] };
}

function indexKey(sampleId, hap) {
  return `${sampleId}:${hap}`;
}

function parseCsvLine(line) {
  const out = [];
  let cur = "";
  let inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (inQuotes) {
      if (ch === '"') {
        if (line[i + 1] === '"') {
          cur += '"';
          i++;
        } else {
          inQuotes = false;
        }
      } else {
        cur += ch;
      }
    } else if (ch === '"') {
      inQuotes = true;
    } else if (ch === ",") {
      out.push(cur);
      cur = "";
    } else {
      cur += ch;
    }
  }
  out.push(cur);
  return out;
}

async function defaultFetchBytes(url, start, end) {
  checkBuildAborted();
  const headers = start != null && end != null ? { Range: `bytes=${start}-${end}` } : undefined;
  const resp = await fetch(url, {
    headers,
    signal: activeBuildAbortSignal ?? undefined,
  });
  if (resp.status !== 200 && resp.status !== 206) {
    throw new Error(`HTTP ${resp.status}`);
  }
  return new Uint8Array(await resp.arrayBuffer());
}

async function defaultFetchText(url) {
  const bytes = await defaultFetchBytes(url, null, null);
  return new TextDecoder().decode(bytes);
}

async function httpFetchBytes(url, start, end) {
  const label = `fetch ${url}`;
  return withRetries(label, () =>
    withFetchSlot(() => {
      if (httpClient) return httpClient.fetchBytes(url, start, end);
      return defaultFetchBytes(url, start, end);
    })
  );
}

async function httpFetchText(url) {
  const label = `fetch ${url}`;
  return withRetries(label, () =>
    withFetchSlot(() => {
      if (httpClient) return httpClient.fetchText(url);
      return defaultFetchText(url);
    })
  );
}

export async function loadAssemblyIndex() {
  if (assemblyIndex) return assemblyIndex;
  const text = await httpFetchText(ASSEMBLY_INDEX_URL);
  const lines = text.split(/\r?\n/).filter((l) => l.trim());
  const map = new Map();
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (i === 0 && line.startsWith("sample_id,")) continue;
    const cols = parseCsvLine(line);
    if (cols.length < 13) continue;
    const [sampleId, hap, , , , , , , , , assemblyFai, assemblyGzi, assembly] = cols;
    if (!sampleId || sampleId === "sample_id") continue;
    map.set(indexKey(sampleId, hap), {
      fa: s3ToHttps(assembly),
      fai: s3ToHttps(assemblyFai),
      gzi: s3ToHttps(assemblyGzi),
    });
  }
  assemblyIndex = map;
  return map;
}

export async function fetchRange(url, start, end) {
  return httpFetchBytes(url, start, end);
}

function parseFaiEntry(contig, faiText) {
  const lines = faiText.split(/\r?\n/);
  for (const line of lines) {
    if (!line) continue;
    const parts = line.split("\t");
    if (parts[0] === contig) {
      return {
        length: Number(parts[1]),
        offset: Number(parts[2]),
        linebases: Number(parts[3]),
        linelength: Number(parts[4]),
      };
    }
  }
  return null;
}

async function getFaiText(faiUrl) {
  if (faiCache.has(faiUrl)) return faiCache.get(faiUrl);
  const text = await httpFetchText(faiUrl);
  faiCache.set(faiUrl, text);
  return text;
}

function parseGzi(buffer) {
  const view = new DataView(buffer.buffer, buffer.byteOffset, buffer.byteLength);
  const n = Number(view.getBigUint64(0, true));
  const entries = [{ compressed: 0, uncompressed: 0 }];
  let off = 8;
  for (let i = 0; i < n; i++) {
    const compressed = Number(view.getBigUint64(off, true));
    const uncompressed = Number(view.getBigUint64(off + 8, true));
    entries.push({ compressed, uncompressed });
    off += 16;
  }
  return entries;
}

async function getGziEntries(gziUrl) {
  if (gziCache.has(gziUrl)) return gziCache.get(gziUrl);
  const buf = await httpFetchBytes(gziUrl, null, null);
  const entries = parseGzi(buf);
  gziCache.set(gziUrl, entries);
  return entries;
}

function findBgzfBlock(gzi, offset) {
  let lo = 0;
  let hi = gzi.length - 1;
  while (lo < hi) {
    const mid = Math.ceil((lo + hi) / 2);
    if (gzi[mid].uncompressed <= offset) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

function sequenceByteRange(fai, start, end) {
  if (end <= start) return { startByte: fai.offset, endByte: fai.offset };
  const startByte = fai.offset + start + Math.floor(start / fai.linebases);
  const endByte = fai.offset + end + Math.floor((end - 1) / fai.linebases);
  return { startByte, endByte };
}

function bgzfBlockSize(buf, offset) {
  const xlen = buf[offset + 10] | (buf[offset + 11] << 8);
  let pos = offset + 12;
  const endExtra = offset + 12 + xlen;
  while (pos + 4 <= endExtra) {
    const subId = buf[pos] | (buf[pos + 1] << 8);
    const subLen = buf[pos + 2] | (buf[pos + 3] << 8);
    if (subId === 0x4342) {
      return buf[pos + 4] | (buf[pos + 5] << 8);
    }
    pos += 4 + subLen;
  }
  throw new Error("BC subfield not found in BGZF block");
}

async function decompressGzipBlock(blockBytes) {
  if (bgzfDecompressor) {
    return bgzfDecompressor(blockBytes);
  }
  const xlen = blockBytes[10] | (blockBytes[11] << 8);
  const deflateStart = 12 + xlen;
  const deflateEnd = blockBytes.length - 8;
  if (deflateEnd <= deflateStart) {
    throw new Error("Invalid BGZF block");
  }
  const deflateData = blockBytes.subarray(deflateStart, deflateEnd);
  const stream = new Blob([deflateData]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function fetchBgzfBlock(faUrl, compressedOff) {
  const headerEnd = compressedOff + 31;
  let chunk = await fetchRange(faUrl, compressedOff, headerEnd);
  const blockSize = bgzfBlockSize(chunk, 0);
  if (chunk.length < blockSize) {
    chunk = await fetchRange(faUrl, compressedOff, compressedOff + blockSize - 1);
  }
  if (chunk.length < blockSize) {
    throw new Error(`Short BGZF read: wanted ${blockSize} bytes, got ${chunk.length}`);
  }
  return chunk.subarray(0, blockSize);
}

async function fetchUncompressedRange(faUrl, gziUrl, startByte, endByte) {
  const gzi = await getGziEntries(gziUrl);
  const startBlockIdx = findBgzfBlock(gzi, startByte);
  const uncompressedBase = gzi[startBlockIdx].uncompressed;
  const parts = [];
  let blockIdx = startBlockIdx;
  let uncompressedEnd = uncompressedBase;

  while (uncompressedEnd < endByte) {
    const block = await fetchBgzfBlock(faUrl, gzi[blockIdx].compressed);
    const inflated = await decompressGzipBlock(block);
    parts.push(inflated);
    blockIdx++;
    if (blockIdx < gzi.length) {
      uncompressedEnd = gzi[blockIdx].uncompressed;
    } else {
      uncompressedEnd += inflated.length;
    }
  }

  const total = parts.reduce((n, p) => n + p.length, 0);
  const merged = new Uint8Array(total);
  let pos = 0;
  for (const p of parts) {
    merged.set(p, pos);
    pos += p.length;
  }

  const relStart = startByte - uncompressedBase;
  const relEnd = endByte - uncompressedBase;
  return merged.subarray(relStart, relEnd);
}

function bytesToSequence(bytes) {
  let out = "";
  for (let i = 0; i < bytes.length; i++) {
    const c = bytes[i];
    if (c === 10 || c === 13) continue;
    out += String.fromCharCode(c);
  }
  return out;
}

export async function fetchSubsequence(urls, contig, start, end) {
  const faiText = await getFaiText(urls.fai);
  const fai = parseFaiEntry(contig, faiText);
  if (!fai) {
    throw new Error(`Contig ${contig} not found in FAI`);
  }
  if (end > fai.length) {
    throw new Error(`Interval end ${end} exceeds contig length ${fai.length}`);
  }
  const { startByte, endByte } = sequenceByteRange(fai, start, end);
  const raw = await fetchUncompressedRange(urls.fa, urls.gzi, startByte, endByte);
  return bytesToSequence(raw);
}

function formatFastaRecord(contig, start, end, seq) {
  return `>${contig}:${start}-${end}\n${seq}\n`;
}

async function createGzipWriterBackend(scratchName) {
  const encoder = new TextEncoder();
  const cs = new CompressionStream("gzip");
  const compressWriter = cs.writable.getWriter();

  if (navigator.storage?.getDirectory) {
    const root = await navigator.storage.getDirectory();
    try {
      await root.removeEntry(scratchName);
    } catch {
      /* missing scratch file */
    }
    const handle = await root.getFileHandle(scratchName, { create: true });
    const opfsWritable = await handle.createWritable();
    const pipeDone = cs.readable.pipeTo(opfsWritable);
    return {
      mode: "opfs",
      async append(text) {
        await compressWriter.write(encoder.encode(text));
      },
      async finish() {
        await compressWriter.close();
        await pipeDone;
        const ab = await (await handle.getFile()).arrayBuffer();
        return new Blob([ab], { type: "application/gzip" });
      },
      async cleanup() {
        try {
          await root.removeEntry(scratchName);
        } catch {
          /* already removed */
        }
      },
    };
  }

  const chunks = [];
  const pump = (async () => {
    const reader = cs.readable.getReader();
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      chunks.push(value);
    }
  })();

  return {
    mode: "memory",
    async append(text) {
      await compressWriter.write(encoder.encode(text));
    },
    async finish() {
      await compressWriter.close();
      await pump;
      return new Blob(chunks, { type: "application/gzip" });
    },
    async cleanup() {},
  };
}

async function createOutputWriter(opts = {}) {
  const gzip = !!opts.gzip;
  const scratchName = gzip ? SCRATCH_FILE_GZ : SCRATCH_FILE;

  if (gzip) {
    return createGzipWriterBackend(scratchName);
  }

  if (navigator.storage?.getDirectory) {
    const root = await navigator.storage.getDirectory();
    try {
      await root.removeEntry(scratchName);
    } catch {
      /* missing scratch file */
    }
    const handle = await root.getFileHandle(scratchName, { create: true });
    const writable = await handle.createWritable();
    return {
      mode: "opfs",
      async append(text) {
        await writable.write(text);
      },
      async finish() {
        await writable.close();
        const file = await handle.getFile();
        const ab = await file.arrayBuffer();
        return new Blob([ab], { type: "text/plain" });
      },
      async cleanup() {
        try {
          await root.removeEntry(scratchName);
        } catch {
          /* already removed */
        }
      },
    };
  }

  const chunks = [];
  return {
    mode: "memory",
    async append(text) {
      chunks.push(text);
    },
    async finish() {
      return new Blob(chunks, { type: "text/plain" });
    },
    async cleanup() {},
  };
}

async function runPool(items, worker, concurrency, onItemDone, signal) {
  const results = new Array(items.length);
  let next = 0;
  async function runOne() {
    while (!signal?.aborted) {
      const i = next++;
      if (i >= items.length) return;
      try {
        results[i] = await worker(items[i], i);
        if (onItemDone) onItemDone(results[i], i);
      } catch (err) {
        if (err?.name === "AbortError") return;
        throw err;
      }
    }
  }
  const n = Math.min(concurrency, items.length);
  await Promise.all(Array.from({ length: n }, runOne));
  if (signal?.aborted) {
    throw new DOMException("FASTA build cancelled", "AbortError");
  }
  return results;
}

/**
 * Build a multi-FASTA file from BED rows.
 * @param {Array<{contig: string, start: number, end: number}>} rows
 * @param {{ concurrency?: number, gzip?: boolean, signal?: AbortSignal, onProgress?: (info: {done: number, total: number, label: string, phase: string}) => void }} opts
 */
export async function buildMultifasta(rows, opts = {}) {
  const onProgress = opts.onProgress ?? (() => {});
  const concurrency = opts.concurrency ?? 6;
  const gzip = !!opts.gzip;
  const signal = opts.signal ?? null;
  activeBuildAbortSignal = signal;
  let writer = null;
  try {
    checkBuildAborted();
    const index = await loadAssemblyIndex();
    onProgress({ done: 0, total: rows.length, label: "assembly index", phase: "start" });
    writer = await createOutputWriter({ gzip });
    const skipped = [];
    let written = 0;
    let fetched = 0;

    const results = await runPool(
      rows,
      async (row) => {
        checkBuildAborted();
        const key = parseContigKey(row.contig);
        if (!key) {
          return { ok: false, label: row.contig, reason: "invalid contig name" };
        }
        const entry = index.get(indexKey(key.sampleId, key.hap));
        if (!entry) {
          return {
            ok: false,
            label: `${key.sampleId}_hap${key.hap}`,
            reason: "not in Release 2 assembly index",
          };
        }
        try {
          const seq = await fetchSubsequence(entry, row.contig, row.start, row.end);
          return {
            ok: true,
            label: `${key.sampleId}_hap${key.hap}`,
            text: formatFastaRecord(row.contig, row.start, row.end, seq),
          };
        } catch (err) {
          if (err?.name === "AbortError") throw err;
          return {
            ok: false,
            label: `${key.sampleId}_hap${key.hap}`,
            reason: err?.message ?? String(err),
          };
        }
      },
      concurrency,
      (result) => {
        fetched++;
        onProgress({
          done: fetched,
          total: rows.length,
          label: result.label,
          phase: "fetch",
        });
      },
      signal
    );

    checkBuildAborted();
    onProgress({ done: 0, total: rows.length, label: "", phase: "write" });
    for (let i = 0; i < results.length; i++) {
      checkBuildAborted();
      const r = results[i];
      if (!r) continue;
      if (r.ok) {
        await writer.append(r.text);
        written++;
      } else {
        skipped.push({ label: r.label, reason: r.reason });
      }
      onProgress({
        done: i + 1,
        total: rows.length,
        label: r.label,
        phase: "write",
      });
    }

    checkBuildAborted();
    const outputMode = writer.mode;
    const file = await writer.finish();
    await writer.cleanup();
    writer = null;
    return { file, written, skipped, mode: outputMode, gzip };
  } finally {
    activeBuildAbortSignal = null;
    if (signal?.aborted && writer) {
      await writer.cleanup();
    }
  }
}

export function triggerDownload(file, filename) {
  const url = URL.createObjectURL(file);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.rel = "noopener";
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 60_000);
}

export function fastaDownloadFilename(bounds, gzip = false) {
  const ext = gzip ? ".fa.gz" : ".fa";
  if (bounds?.contig != null && bounds.start != null && bounds.end != null) {
    const contig = String(bounds.contig).replace(/[^\w.#-]+/g, "_");
    return `extract_${contig}_${bounds.start}-${bounds.end}${ext}`;
  }
  return `extract${ext}`;
}
