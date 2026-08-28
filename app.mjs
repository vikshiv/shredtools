import {
  isFastaDownloadSupported,
  isFastaBuildCancelled,
  buildMultifasta,
  triggerDownload,
  fastaDownloadFilename,
  setBgzfDecompressor,
} from "./fasta_download.js";

const $ = (id) => document.getElementById(id);
    const THEME_KEY = "shredtools-theme";
    const CUSTOM_DATASETS_KEY = "shredtools.customDatasets";
    const HUMAN_TAB_ID = "human";
    const CICHLID_TAB_ID = "cichlid";
    const ADD_TAB_ID = "__add__";

    const HUMAN_EMOJI = "🧍";
    const CICHLID_EMOJI = "🐟";
    const CICHLID_BUMBL_URL =
      "https://ftp.ebi.ac.uk/pub/databases/metagenomics/research-team/shivakumar/cichlids.bumbl";

    const ROUTE_BY_TAB = {
      [HUMAN_TAB_ID]: "hprc",
      [CICHLID_TAB_ID]: "cichlid",
    };
    const TAB_BY_ROUTE = {
      hprc: HUMAN_TAB_ID,
      cichlid: CICHLID_TAB_ID,
    };

    function pathnameLastSegment() {
      let p = location.pathname.replace(/\/+$/, "");
      if (p.endsWith("/index.html")) p = p.slice(0, -"/index.html".length);
      const parts = p.split("/").filter(Boolean);
      return (parts[parts.length - 1] || "").toLowerCase();
    }

    function isRoutedPath() {
      const last = pathnameLastSegment();
      return last === "hprc" || last === "cichlid";
    }

    function appBasePath() {
      let p = location.pathname.replace(/\/+$/, "");
      if (p.endsWith("/index.html")) p = p.slice(0, -"/index.html".length);
      const parts = p.split("/").filter(Boolean);
      const last = (parts[parts.length - 1] || "").toLowerCase();
      if (last === "hprc" || last === "cichlid") parts.pop();
      return parts.length ? "/" + parts.join("/") : "";
    }

    function urlForRoute(route) {
      const base = appBasePath();
      return `${base}/${route}`.replace(/\/{2,}/g, "/");
    }

    function routeFromLocation() {
      const forced = String(window.__SHREDTOOLS_ROUTE__ || "").trim().toLowerCase();
      if (forced && TAB_BY_ROUTE[forced]) return forced;
      const last = pathnameLastSegment();
      if (TAB_BY_ROUTE[last]) return last;
      return ROUTE_BY_TAB[HUMAN_TAB_ID];
    }

    function tabIdFromRoute(route) {
      return TAB_BY_ROUTE[route] || HUMAN_TAB_ID;
    }

    function appAssetBase() {
      return new URL(".", import.meta.url).href;
    }

    function updateDocumentTitle(tab) {
      if (!tab) return;
      let suffix;
      if (tab.kind === "human") suffix = "HPRC";
      else if (tab.id === CICHLID_TAB_ID) suffix = "Cichlid";
      else suffix = tab.label;
      const pageTitle = `Shredtools extract - ${suffix}`;
      document.title = `Shredtools - ${suffix}`;
      const h1 = $("appTitle");
      if (h1) h1.textContent = pageTitle;
    }

    function updateBrowserRoute(tabId, { replace = false } = {}) {
      const route = ROUTE_BY_TAB[tabId];
      if (!route) return;
      const targetPath = urlForRoute(route);
      const currentPath = location.pathname.replace(/\/+$/, "") || "/";
      const targetNorm = targetPath.replace(/\/+$/, "") || "/";
      if (currentPath === targetNorm) return;
      const target = targetPath + location.search + location.hash;
      const state = { tabId };
      if (replace) history.replaceState(state, "", target);
      else history.pushState(state, "", target);
    }

    function getTheme() {
      return localStorage.getItem(THEME_KEY) || "dark";
    }
    document.documentElement.dataset.theme = getTheme();
    const t0 = performance.now();
    let statusEpoch = 0;
    const status = (t, sub = "") => {
      $("status").textContent = t;
      $("substatus").textContent = sub;
    };
    const beginStatus = () => ++statusEpoch;
    const statusAt = (epoch, t, sub = "") => {
      if (epoch !== statusEpoch) return;
      status(t, sub);
    };
    const secs = (tStart = t0) => ((performance.now() - tStart) / 1000).toFixed(1) + "s";

    function renderIntroMarkdown(md, wrap) {
      const restEl = $("appIntroRest");
      const firstEl = $("appIntroFirst");
      const detailsEl = $("appIntroDetails");
      const s = String(md || "").trim();
      if (!s) {
        wrap.style.display = "none";
        restEl.innerHTML = "";
        firstEl.innerHTML = "";
        restEl.hidden = true;
        detailsEl.hidden = true;
        return;
      }
      let first;
      let rest;
      const sep = "\n\n---\n\n";
      const sepIdx = s.indexOf(sep);
      if (sepIdx === -1) {
        const m = /\n\s*\n/.exec(s);
        if (!m) {
          first = s;
          rest = "";
        } else {
          first = s.slice(0, m.index).trim();
          rest = s.slice(m.index).replace(/^\s*\n+/, "").trim();
        }
      } else {
        first = s.slice(0, sepIdx).trim();
        rest = s.slice(sepIdx + sep.length).trim();
      }

      const parse =
        typeof marked !== "undefined" && typeof marked.parse === "function"
          ? marked.parse.bind(marked)
          : null;

      function setHtml(el, src) {
        if (!src) {
          el.innerHTML = "";
          return;
        }
        if (parse) {
          const out = parse(src, { breaks: true });
          if (out instanceof Promise) {
            void out.then((html) => {
              el.innerHTML = html;
            });
          } else {
            el.innerHTML = out;
          }
        } else {
          el.textContent = src;
        }
      }

      wrap.style.display = "";
      restEl.hidden = !rest;
      setHtml(restEl, rest);
      detailsEl.hidden = !first;
      setHtml(firstEl, first);
      wrap.dataset.hasIntro = "1";
    }

    function prettyGenome(label) {
      // Standard: GENOME#HAP_NUM#CONTIG -> GENOME_hapHAP_NUM
      const s = String(label ?? "");
      const parts = s.split("#");
      if (parts.length >= 2) return `${parts[0]}_hap${parts[1]}`;
      return s;
    }

    const ASSEMBLY_EXTS = [
      ".fasta.gz",
      ".fastq.gz",
      ".fa.gz",
      ".fna.gz",
      ".fq.gz",
      ".fasta.bz2",
      ".fa.bz2",
      ".fna.bz2",
      ".fasta.xz",
      ".fa.xz",
      ".fna.xz",
      ".fasta",
      ".fastq",
      ".fna",
      ".fa",
      ".fq",
    ];

    function cleanAssemblyLabel(name) {
      let baseName = String(name || "").trim();
      if (!baseName) return "";
      baseName = baseName.split(/[/\\]/).pop();
      const lower = baseName.toLowerCase();
      for (const ext of ASSEMBLY_EXTS) {
        if (lower.endsWith(ext)) return baseName.slice(0, -ext.length);
      }
      return baseName;
    }

    /** Assembly labels from multilengths ``*`` header lines (column 0). */
    function parseAssemblyLabelsFromLengths(text) {
      const labels = [];
      for (const raw of String(text || "").split(/\r?\n/)) {
        const line = raw.trim();
        if (!line) continue;
        const parts = line.split(/\s+/);
        if (parts.length >= 2 && parts[1] === "*") {
          labels.push(cleanAssemblyLabel(parts[0]));
        }
      }
      return labels;
    }

    function applyAssemblyLabelsToManifest(mani, labels) {
      if (!mani?.genomes || !labels?.length) return mani;
      for (const g of mani.genomes) {
        const lab = labels[g.seq_idx];
        if (lab) g.label = lab;
      }
      return mani;
    }

    /** Assembly display id for a seq_idx (custom datasets); null for HPRC / unknown. */
    function assemblyNameForSeq(seqIdx) {
      if (isHumanTab()) return null;
      const i = Number(seqIdx);
      if (!Number.isFinite(i)) return null;
      const tab = activeTab();
      if (tab?.assemblyLabels?.[i]) return tab.assemblyLabels[i];
      const g = getGenome(i);
      if (g?.label) return String(g.label);
      return null;
    }

    function genomeDisplayForRow(r) {
      const raw = String(r?.contig ?? "");
      const asm = assemblyNameForSeq(r?.seq_idx);
      if (asm) return asm;
      // HPRC / PanSN contigs: GENOME#HAP#CONTIG
      return prettyGenome(raw || r?.label || `seq_${r?.seq_idx ?? ""}`);
    }

    function contigDisplayForRow(r) {
      return prettyContig(String(r?.contig ?? ""));
    }

    function bedLineForRow(r) {
      const contig = contigDisplayForRow(r);
      const start = r?.start ?? "";
      const end = r?.end ?? "";
      if (isHumanTab()) {
        // HPRC: keep full contig id in col1; name field unused.
        return `${r.contig}\t${start}\t${end}\t.\n`;
      }
      const asm = assemblyNameForSeq(r?.seq_idx) || ".";
      return `${contig}\t${start}\t${end}\t${asm}\n`;
    }

    function prettyContig(label) {
      // Standard: GENOME#HAP_NUM#CONTIG -> CONTIG
      const s = String(label ?? "");
      const parts = s.split("#");
      if (parts.length >= 3) return parts.slice(2).join("#");
      return s;
    }

    function regionLength(r) {
      const d = Number(r?.end) - Number(r?.start);
      return Number.isFinite(d) ? Math.max(0, d) : 0;
    }

    status("Loading Pyodide…", `elapsed ${secs()}`);
    const pyodide = await loadPyodide();
    // No micropip installs needed: we vendor the tiny mumemto subset we need.
    const base = appAssetBase();
    // Bust CDN/browser caches of MEMFS sources (GitHub Pages often serves stale .py).
    const _assetQ = `v=${Date.now()}`;
    status("Loading app code…", `fetching local scripts (elapsed ${secs()})`);
    const uRes = await fetch(new URL(`bumbl_index_utils.py?${_assetQ}`, base), { cache: "no-store" });
    const aRes = await fetch(new URL(`app.py?${_assetQ}`, base), { cache: "no-store" });
    if (!uRes.ok) throw new Error("Could not load bumbl_index_utils.py (serve this folder over HTTP, not file://).");
    if (!aRes.ok) throw new Error("Could not load app.py.");
    const [uText, aText] = [await uRes.text(), await aRes.text()];
    pyodide.FS.writeFile("/bumbl_index_utils.py", uText);
    pyodide.FS.writeFile("/app.py", aText);
    status("Initializing Python…", `importing app (elapsed ${secs()})`);
    await pyodide.runPythonAsync(`
import sys
if "/" not in sys.path:
    sys.path.insert(0, "/")
# Ensure a fresh import after writing MEMFS sources (avoids stale module cache).
for _mod in ("app", "bumbl_index_utils"):
    sys.modules.pop(_mod, None)
    `);
    await pyodide.runPythonAsync("import app");

    function bytesFromB64(b64) {
      const bin = atob(b64);
      const out = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
      return out;
    }

    setBgzfDecompressor(async (blockBytes) => {
      pyodide.FS.writeFile("/bgzf_block.bin", blockBytes);
      const b64 = await pyodide.runPythonAsync(`
import base64, zlib
with open("/bgzf_block.bin", "rb") as f:
    block = f.read()
xlen = block[10] | (block[11] << 8)
raw = block[12 + xlen : len(block) - 8]
data = zlib.decompressobj(-15).decompress(raw)
base64.b64encode(data).decode()
      `);
      return bytesFromB64(b64);
    });

    const jRes = await fetch(new URL("pangenome_lengths.json", base));
    if (!jRes.ok) {
      throw new Error(
        `Missing pangenome_lengths.json (HTTP ${jRes.status}). Build it with lengths_to_json.py (see docstring in that file).`
      );
    }
    const jBuf = await jRes.arrayBuffer();
    pyodide.FS.writeFile("/pangenome_lengths.json", new Uint8Array(jBuf));
    await pyodide.runPythonAsync("import app; app.load_lengths_bundle_path('/pangenome_lengths.json')");
    const intro = String(
      await pyodide.runPythonAsync("import app as _a; getattr(_a, 'INTRO', '') or ''")
    ).trim();
    renderIntroMarkdown(intro, $("appIntroWrap"));

    let manifest = null;
    let selectedSeqs = null; // Set<number>; null means not initialized yet
    let lastRows = null;
    let lastMeta = null;
    let lastUnavailable = null;
    let lastFastaBlob = null;
    let lastFastaFilename = "extract.fa";
    let lastFastaGzip = false;
    let geneIndex = null; // { [gene: string]: Array<{contig,start,end,label}> }
    let pangenomeList = null;
    /** Lengths: ``pangenome_lengths.json`` → MEMFS → ``app.load_lengths_bundle_path``. */
    let mplReady = null;
    let plotModulesReady = null;

    /** @type {{ id: string, kind: "human"|"builtin"|"custom", label: string, emoji?: string, bumblUrl?: string, lengthsUrl?: string, biUrl?: string, annotationUrl?: string, pyKey?: string, annotationJson?: object|null, loadError?: string|null, registered?: boolean }} */
    let datasetTabs = [
      { id: HUMAN_TAB_ID, kind: "human", label: "Human" },
      {
        id: CICHLID_TAB_ID,
        kind: "builtin",
        label: "Cichlid",
        emoji: CICHLID_EMOJI,
        bumblUrl: CICHLID_BUMBL_URL,
        lengthsUrl: CICHLID_BUMBL_URL.replace(/\.bumbl$/i, ".lengths"),
        biUrl: `${CICHLID_BUMBL_URL}.bi`,
        pyKey: "builtin_cichlid",
        annotationJson: null,
        loadError: null,
        registered: false,
      },
    ];
    let activeTabId = HUMAN_TAB_ID;
    let lengthsUrlTouched = false;
    let switchingTab = false;
    let emojiSelectReady = null;
    /** @type {Array<{e: string, n: string, g: string, q: string}>} */
    let emojiCatalog = [];
    let emojiHighlight = -1;

    function setSelectedEmoji(emoji) {
      const val = normalizeEmojiIcon(emoji || "");
      $("addEmoji").value = val;
      $("addEmojiPreview").textContent = val;
      if (!val) {
        $("addEmojiSearch").value = "";
      } else if (!$("addEmojiSearch").value.trim()) {
        $("addEmojiSearch").value = val;
      }
    }

    function hideEmojiList() {
      const list = $("addEmojiList");
      list.hidden = true;
      list.innerHTML = "";
      emojiHighlight = -1;
    }

    function showEmojiResults(query) {
      const list = $("addEmojiList");
      const q = String(query || "").trim().toLowerCase();
      list.innerHTML = "";
      emojiHighlight = -1;
      if (!emojiCatalog.length) {
        list.hidden = true;
        return;
      }
      // Empty query: show a short starter set from each group; otherwise filter.
      let hits;
      if (!q) {
        hits = emojiCatalog.slice(0, 80);
      } else {
        hits = [];
        for (const item of emojiCatalog) {
          if (item.e === q || item.q.includes(q) || item.e.includes(query.trim())) {
            hits.push(item);
            if (hits.length >= 120) break;
          }
        }
      }
      if (!hits.length) {
        const empty = document.createElement("div");
        empty.className = "emoji-picker-empty";
        empty.textContent = "No matching emoji.";
        list.appendChild(empty);
        list.hidden = false;
        return;
      }
      const frag = document.createDocumentFragment();
      hits.forEach((item, i) => {
        const btn = document.createElement("button");
        btn.type = "button";
        btn.className = "emoji-picker-item";
        btn.role = "option";
        btn.dataset.index = String(i);
        btn.dataset.emoji = item.e;
        btn.innerHTML =
          `<span class="e">${item.e}</span><span class="n">${item.n} · ${item.g}</span>`;
        btn.addEventListener("mousedown", (ev) => {
          // mousedown before blur so selection sticks
          ev.preventDefault();
          setSelectedEmoji(item.e);
          $("addEmojiSearch").value = `${item.e} ${item.n}`;
          hideEmojiList();
        });
        frag.appendChild(btn);
      });
      list.appendChild(frag);
      list.hidden = false;
    }

    function moveEmojiHighlight(delta) {
      const list = $("addEmojiList");
      if (list.hidden) return;
      const items = Array.from(list.querySelectorAll(".emoji-picker-item"));
      if (!items.length) return;
      emojiHighlight = (emojiHighlight + delta + items.length) % items.length;
      items.forEach((el, i) => el.setAttribute("aria-selected", i === emojiHighlight ? "true" : "false"));
      items[emojiHighlight].scrollIntoView({ block: "nearest" });
    }

    async function ensureEmojiSelect() {
      if (emojiSelectReady) return emojiSelectReady;
      emojiSelectReady = (async () => {
        const hint = $("addEmojiHint");
        try {
          const res = await fetch(new URL("emoji-by-group.json", base));
          if (!res.ok) throw new Error(`HTTP ${res.status}`);
          const groups = await res.json();
          const catalog = [];
          for (const g of groups) {
            const gName = String(g.name || "");
            for (const item of g.emojis || []) {
              const e = String(item.e || "");
              const n = String(item.n || "");
              if (!e) continue;
              catalog.push({
                e,
                n,
                g: gName,
                q: `${n} ${gName}`.toLowerCase(),
              });
            }
          }
          emojiCatalog = catalog;
          if (hint) hint.textContent = "Type to search by name, or pick from the list.";
          wireEmojiPicker();
        } catch (e) {
          console.warn("[shredtools] emoji list load failed", e);
          if (hint) hint.textContent = "Emoji list unavailable; leave blank to use the name.";
        }
      })();
      return emojiSelectReady;
    }

    let emojiPickerWired = false;
    function wireEmojiPicker() {
      if (emojiPickerWired) return;
      emojiPickerWired = true;
      const search = $("addEmojiSearch");
      const clearBtn = $("addEmojiClear");
      search.addEventListener("focus", () => {
        showEmojiResults(search.value);
      });
      search.addEventListener("input", () => {
        // Typing clears a prior selection until the user picks again.
        if ($("addEmoji").value && !search.value.includes($("addEmoji").value)) {
          $("addEmoji").value = "";
          $("addEmojiPreview").textContent = "";
        }
        showEmojiResults(search.value);
      });
      search.addEventListener("keydown", (ev) => {
        if (ev.key === "ArrowDown") {
          ev.preventDefault();
          if ($("addEmojiList").hidden) showEmojiResults(search.value);
          moveEmojiHighlight(1);
        } else if (ev.key === "ArrowUp") {
          ev.preventDefault();
          moveEmojiHighlight(-1);
        } else if (ev.key === "Enter") {
          const list = $("addEmojiList");
          if (!list.hidden && emojiHighlight >= 0) {
            ev.preventDefault();
            const item = list.querySelector(`.emoji-picker-item[data-index="${emojiHighlight}"]`);
            if (item) {
              setSelectedEmoji(item.dataset.emoji);
              search.value = `${item.dataset.emoji} ${item.querySelector(".n")?.textContent?.split(" · ")[0] || ""}`.trim();
              hideEmojiList();
            }
          }
        } else if (ev.key === "Escape") {
          hideEmojiList();
        }
      });
      search.addEventListener("blur", () => {
        // Delay so option mousedown can run first.
        setTimeout(() => hideEmojiList(), 120);
      });
      clearBtn.addEventListener("click", () => {
        setSelectedEmoji("");
        hideEmojiList();
        $("addEmojiSearch").focus();
      });
    }

    function activeTab() {
      return datasetTabs.find((t) => t.id === activeTabId) || datasetTabs[0];
    }

    function isHumanTab() {
      return activeTabId === HUMAN_TAB_ID;
    }

    function defaultLengthsUrl(bumblUrl) {
      const s = String(bumblUrl || "").trim();
      if (/\.bumbl$/i.test(s)) return s.replace(/\.bumbl$/i, ".lengths");
      return s ? `${s}.lengths` : "";
    }

    function deriveBiUrl(bumblUrl) {
      return `${String(bumblUrl || "").trim()}.bi`;
    }

    function shortTabLabel(name) {
      const s = String(name || "").trim() || "?";
      if (s.length <= 6) return s;
      return s.slice(0, 5) + "…";
    }

    /** Keep at most one visible emoji / short grapheme cluster for the sidebar. */
    function normalizeEmojiIcon(raw) {
      const s = String(raw || "").trim();
      if (!s) return "";
      const chars = Array.from(s);
      // Allow a single emoji that may use a ZWJ sequence (take first ~8 code points max).
      if (chars.length <= 8) return chars.join("");
      return chars.slice(0, 8).join("");
    }

    function loadStoredCustomDatasets() {
      try {
        const raw = localStorage.getItem(CUSTOM_DATASETS_KEY);
        if (!raw) return [];
        const parsed = JSON.parse(raw);
        if (!Array.isArray(parsed)) return [];
        return parsed
          .filter((e) => e && typeof e === "object" && e.bumblUrl && e.name)
          .map((e) => ({
            id: String(e.id || `custom_${Math.random().toString(36).slice(2, 10)}`),
            name: String(e.name),
            bumblUrl: String(e.bumblUrl),
            lengthsUrl: String(e.lengthsUrl || defaultLengthsUrl(e.bumblUrl)),
            annotationUrl: e.annotationUrl ? String(e.annotationUrl) : "",
            emoji: normalizeEmojiIcon(e.emoji || ""),
          }));
      } catch {
        return [];
      }
    }

    function saveStoredCustomDatasets() {
      const payload = datasetTabs
        .filter((t) => t.kind === "custom")
        .map((t) => ({
          id: t.id,
          name: t.label,
          bumblUrl: t.bumblUrl,
          lengthsUrl: t.lengthsUrl,
          annotationUrl: t.annotationUrl || "",
          emoji: t.emoji || "",
        }));
      try {
        localStorage.setItem(CUSTOM_DATASETS_KEY, JSON.stringify(payload));
      } catch (e) {
        console.warn("[shredtools] localStorage save failed", e);
      }
    }

    function clearQueryResults() {
      lastRows = null;
      lastMeta = null;
      lastUnavailable = null;
      $("out").value = "";
      $("bounds").textContent = "";
      renderIntervalsTable([]);
      renderUnavailableTable([]);
      renderLengthsView([]);
      showPlotButton(false);
      showFastaDownloadButton(false);
      clearFastaSave();
      clearSyntenyPlot();
    }

    function setMainView(mode) {
      const add = $("viewAddDataset");
      const query = $("viewQuery");
      if (mode === "add") {
        add.classList.add("active");
        query.classList.add("hidden-view");
      } else {
        add.classList.remove("active");
        query.classList.remove("hidden-view");
      }
    }

    function renderSidebar() {
      const host = $("sidebarTabs");
      host.innerHTML = "";
      for (const tab of datasetTabs) {
        const btn = document.createElement("button");
        btn.type = "button";
        const useEmoji = tab.kind === "human" || !!(tab.emoji && String(tab.emoji).trim());
        btn.className =
          "sidebar-btn" +
          (useEmoji ? " sidebar-emoji" : " sidebar-name");
        btn.dataset.tabId = tab.id;
        btn.setAttribute("aria-selected", tab.id === activeTabId && activeTabId !== ADD_TAB_ID ? "true" : "false");
        btn.title = tab.label;
        btn.setAttribute("aria-label", tab.label);
        if (tab.kind === "human") {
          btn.textContent = HUMAN_EMOJI;
        } else if (tab.emoji) {
          btn.textContent = tab.emoji;
        } else {
          btn.textContent = shortTabLabel(tab.label);
        }
        btn.addEventListener("click", () => {
          void switchToTab(tab.id);
        });
        host.appendChild(btn);
      }
      const addBtn = $("sidebarAdd");
      addBtn.setAttribute("aria-selected", activeTabId === ADD_TAB_ID ? "true" : "false");
    }

    function setHumanUiChrome(visible) {
      $("pangenomeRow").style.display = visible ? "" : "none";
      const clearEl = $("headerClear");
      if (!visible) {
        $("appIntroWrap").style.display = "none";
        // Let query controls sit beside the floated logo (avoids a logo-height gap).
        clearEl.classList.add("logo-clear-skip");
      } else {
        clearEl.classList.remove("logo-clear-skip");
        if ($("appIntroWrap").dataset.hasIntro === "1") {
          $("appIntroWrap").style.display = "";
        }
      }
    }

    function setCustomRemoveVisible(visible) {
      $("customDatasetRemoveRow").classList.toggle("visible", !!visible);
    }

    function fillBuildOptionsFromKeys(keys, selected) {
      const sel = $("build");
      sel.innerHTML = "";
      for (const k of keys) {
        const o = document.createElement("option");
        o.value = k;
        o.textContent = k;
        sel.appendChild(o);
      }
      if (selected && keys.includes(selected)) sel.value = selected;
      else if (keys.length) sel.value = keys[0];
    }

    function resetHumanBuildOptions() {
      fillBuildOptionsFromKeys(["CHM13v2.0", "GRCh38"], "CHM13v2.0");
    }

    function findGenomeForContig(contigWanted) {
      const want = String(contigWanted || "");
      const hasChr = want.startsWith("chr");
      const alt = hasChr ? want.slice(3) : `chr${want}`;
      if (!manifest?.genomes) return null;
      for (const g of manifest.genomes) {
        for (const v of g.contigs || []) {
          const pc = prettyContig(v);
          if (v === want || v === alt || pc === want || pc === alt) {
            return { seqIdx: g.seq_idx, contig: v };
          }
        }
      }
      return null;
    }

    async function registerCustomInPython(tab) {
      const lengthsRes = await fetch(tab.lengthsUrl);
      if (!lengthsRes.ok) {
        throw new Error(`Lengths fetch failed (HTTP ${lengthsRes.status}). Check URL and CORS.`);
      }
      const lengthsText = await lengthsRes.text();
      tab.assemblyLabels = parseAssemblyLabelsFromLengths(lengthsText);
      const biUrl = tab.biUrl || deriveBiUrl(tab.bumblUrl);
      tab.biUrl = biUrl;
      const memPath = `/custom_lengths_${tab.pyKey}.txt`;
      pyodide.FS.writeFile(memPath, lengthsText);
      await pyodide.runPythonAsync(`
import app
with open(${JSON.stringify(memPath)}, "r", encoding="utf-8") as _f:
    _lengths_text = _f.read()
app.register_custom_pangenome(
    ${JSON.stringify(tab.pyKey)},
    ${JSON.stringify(tab.label)},
    ${JSON.stringify(tab.bumblUrl)},
    ${JSON.stringify(biUrl)},
    _lengths_text,
)
      `);

      try {
        if (tab.annotationUrl) {
          const aRes = await fetch(tab.annotationUrl);
          if (!aRes.ok) {
            throw new Error(`Annotation fetch failed (HTTP ${aRes.status}).`);
          }
          const parsed = await aRes.json();
          if (!parsed || typeof parsed !== "object" || Array.isArray(parsed)) {
            throw new Error("Annotation JSON must be an object of build → genes.");
          }
          tab.annotationJson = parsed;
        } else {
          tab.annotationJson = null;
        }
      } catch (e) {
        try {
          await pyodide.runPythonAsync(
            `import app; app.unregister_custom_pangenome(${JSON.stringify(tab.pyKey)})`
          );
        } catch (_) {
          /* ignore */
        }
        tab.registered = false;
        throw e;
      }
      tab.registered = true;
      tab.loadError = null;
    }

    async function getActivePangenomeKey() {
      return String(await pyodide.runPythonAsync("import app; app.ACTIVE_PANGENOME"));
    }

    async function setActivePangenomeKey(key) {
      await pyodide.runPythonAsync(`import app; app.set_active_pangenome(${JSON.stringify(key)})`);
    }

    /** Invalidate query UI so Run cannot use a stale manifest after a failed switch. */
    function markQueryUiUnloaded(hint) {
      manifest = null;
      selectedSeqs = null;
      $("genome").innerHTML = "";
      $("contig").innerHTML = "";
      $("genome").disabled = true;
      $("contig").disabled = true;
      $("start").disabled = true;
      $("end").disabled = true;
      $("region").disabled = true;
      $("run").disabled = true;
      $("lucky").disabled = true;
      if (hint) $("lengthsHint").textContent = hint;
      setGeneUiVisible(false);
      showPlotButton(false);
      showFastaDownloadButton(false);
    }

    async function switchToTab(tabId, { fromPopstate = false } = {}) {
      if (switchingTab) return;
      if (tabId === ADD_TAB_ID) {
        activeTabId = ADD_TAB_ID;
        renderSidebar();
        setMainView("add");
        $("addDatasetError").textContent = "";
        return;
      }
      const tab = datasetTabs.find((t) => t.id === tabId);
      if (!tab) return;
      switchingTab = true;
      const prevPyKey = await getActivePangenomeKey();
      try {
        activeTabId = tab.id;
        renderSidebar();
        setMainView("query");
        clearQueryResults();
        setCustomRemoveVisible(tab.kind === "custom");
        updateDocumentTitle(tab);
        if (!fromPopstate && (tab.kind === "human" || tab.kind === "builtin")) {
          updateBrowserRoute(tab.id);
        }

        if (tab.kind === "human") {
          setHumanUiChrome(true);
          resetHumanBuildOptions();
          $("pangenome").disabled = true;
          const key = $("pangenome").value || "hprcv2_enhanced";
          await loadPangenomeAndInitUI(key);
          $("pangenome").disabled = false;
          status("Ready.", "Human (HPRC)");
          return;
        }

        setHumanUiChrome(false);
        if (!tab.registered || !tab.assemblyLabels?.length) {
          status("Loading dataset…", tab.label);
          try {
            await registerCustomInPython(tab);
          } catch (e) {
            tab.loadError = String(e?.message ?? e);
            tab.registered = false;
            markQueryUiUnloaded("Dataset failed to load.");
            status("Could not load dataset.", tab.loadError);
            return;
          }
        }
        if (tab.loadError && !tab.registered) {
          markQueryUiUnloaded("Dataset failed to load.");
          status("Could not load dataset.", tab.loadError);
          return;
        }
        if (tab.annotationJson) {
          const builds = Object.keys(tab.annotationJson);
          fillBuildOptionsFromKeys(builds, builds[0]);
        } else {
          setGeneUiVisible(false);
        }
        await setActivePangenomeKey(tab.pyKey);
        try {
          await initLengthsUI();
          if (tab.annotationJson) {
            setGeneUiVisible(true);
            $("geneQuery").value = "";
            await loadGeneIndexFromAnnotation($("build").value);
          } else {
            setGeneUiVisible(false);
          }
        } catch (e) {
          // Roll back Python so ACTIVE_PANGENOME cannot disagree with JS manifest.
          try {
            if (prevPyKey && prevPyKey !== tab.pyKey) await setActivePangenomeKey(prevPyKey);
          } catch (_) {
            /* ignore */
          }
          markQueryUiUnloaded("Dataset switch failed (see status).");
          throw e;
        }
        showFastaDownloadButton(false);
        status("Ready.", tab.label);
      } catch (e) {
        status("Error", String(e?.message ?? e));
      } finally {
        switchingTab = false;
      }
    }

    async function addCustomDatasetFromForm({ name, bumblUrl, lengthsUrl, annotationUrl, emoji }) {
      const id = `custom_${Date.now().toString(36)}_${Math.random().toString(36).slice(2, 7)}`;
      const pyKey = id;
      const tab = {
        id,
        kind: "custom",
        label: name,
        bumblUrl,
        lengthsUrl,
        biUrl: deriveBiUrl(bumblUrl),
        annotationUrl: annotationUrl || "",
        emoji: normalizeEmojiIcon(emoji || ""),
        pyKey,
        annotationJson: null,
        loadError: null,
        registered: false,
      };
      status("Loading custom dataset…", name);
      await registerCustomInPython(tab);
      datasetTabs.push(tab);
      saveStoredCustomDatasets();
      renderSidebar();
      await switchToTab(tab.id);
    }

    async function removeActiveCustomDataset() {
      const tab = activeTab();
      if (!tab || tab.kind !== "custom") return;
      if (!window.confirm(`Remove dataset “${tab.label}”? This only deletes the saved links.`)) return;
      try {
        if (tab.registered) {
          await pyodide.runPythonAsync(`import app; app.unregister_custom_pangenome(${JSON.stringify(tab.pyKey)})`);
        }
      } catch (e) {
        console.warn(e);
      }
      datasetTabs = datasetTabs.filter((t) => t.id !== tab.id);
      saveStoredCustomDatasets();
      renderSidebar();
      await switchToTab(HUMAN_TAB_ID);
    }

    async function restoreCustomDatasetsOnBoot() {
      const stored = loadStoredCustomDatasets();
      for (const e of stored) {
        datasetTabs.push({
          id: e.id,
          kind: "custom",
          label: e.name,
          bumblUrl: e.bumblUrl,
          lengthsUrl: e.lengthsUrl,
          biUrl: deriveBiUrl(e.bumblUrl),
          annotationUrl: e.annotationUrl || "",
          emoji: e.emoji || "",
          pyKey: e.id,
          annotationJson: null,
          loadError: null,
          registered: false,
        });
      }
      // Eagerly register custom datasets; built-in remote tabs load on first select.
      for (const tab of datasetTabs.filter((t) => t.kind === "custom")) {
        try {
          await registerCustomInPython(tab);
        } catch (err) {
          tab.loadError = String(err?.message ?? err);
          tab.registered = false;
          console.warn("[shredtools] restore failed for", tab.label, err);
        }
      }
      renderSidebar();
    }

    const LENGTH_HIST_NUM_BINS = 40;

    async function ensurePlotModules() {
      if (!plotModulesReady) {
        plotModulesReady = (async () => {
          const plotBase = appAssetBase();
          const plotQ = `v=${Date.now()}`;
          const [vRes, sRes] = await Promise.all([
            fetch(new URL(`viz_mums.py?${plotQ}`, plotBase), { cache: "no-store" }),
            fetch(new URL(`synteny_plot.py?${plotQ}`, plotBase), { cache: "no-store" }),
          ]);
          if (!vRes.ok) throw new Error("Could not load viz_mums.py.");
          if (!sRes.ok) throw new Error("Could not load synteny_plot.py.");
          const [vText, sText] = await Promise.all([vRes.text(), sRes.text()]);
          pyodide.FS.writeFile("/viz_mums.py", vText);
          pyodide.FS.writeFile("/synteny_plot.py", sText);
          await pyodide.runPythonAsync(`
import sys
for _mod in ("viz_mums", "synteny_plot"):
    sys.modules.pop(_mod, None)
import synteny_plot
          `);
        })();
      }
      return plotModulesReady;
    }

    function ensureMpl() {
      if (!mplReady) {
        mplReady = pyodide.loadPackage(["matplotlib", "numpy"]);
      }
      return mplReady;
    }

    function showPlotButton(show) {
      const btn = $("plotBtn");
      btn.classList.toggle("visible", !!show);
      btn.disabled = !show;
    }

    function showFastaDownloadButton(show) {
      const btn = $("fastaDownloadBtn");
      const gzipWrap = $("fastaGzipWrap");
      const gzipCb = $("fastaGzip");
      const supported = isHumanTab() && isFastaDownloadSupported($("pangenome").value);
      const visible = !!show && supported;
      btn.classList.toggle("visible", visible);
      btn.disabled = !visible;
      gzipWrap.classList.toggle("visible", visible);
      gzipCb.disabled = !visible;
      btn.title = supported
        ? ""
        : "FASTA download is available for HPRCr2 pangenomes only.";
    }

    function refreshFastaDownloadButton() {
      showFastaDownloadButton(!!lastRows?.length);
    }

    function clearFastaSave() {
      lastFastaBlob = null;
      lastFastaGzip = false;
      setFastaProgress(false);
      const saveBtn = $("fastaSaveBtn");
      saveBtn.textContent = "Save FASTA";
      saveBtn.classList.remove("visible");
      saveBtn.disabled = true;
    }

    function showFastaSaveButton(show) {
      const btn = $("fastaSaveBtn");
      btn.textContent = lastFastaGzip ? "Save FASTA.gz" : "Save FASTA";
      btn.classList.toggle("visible", !!show);
      btn.disabled = !show || !lastFastaBlob;
    }

    let fastaBuildAbort = null;

    function showFastaCancelButton(show) {
      const btn = $("fastaCancelBtn");
      btn.classList.toggle("visible", !!show);
      btn.hidden = !show;
      btn.disabled = !show;
    }

    function setFastaProgress(on, done = 0, total = 0, label = "", phase = "fetch") {
      const wrap = $("fastaProgress");
      const bar = $("fastaProgressBar");
      const text = $("fastaProgressText");
      const pct = total > 0 ? Math.round((done / total) * 100) : 0;
      wrap.classList.toggle("visible", !!on);
      wrap.hidden = !on;
      bar.style.width = `${pct}%`;
      bar.setAttribute("aria-valuenow", String(pct));
      if (!on) {
        text.textContent = "";
        return;
      }
      const phaseLabel =
        phase === "write" ? "Writing" : phase === "start" ? "Loading" : "Fetching";
      const detail = label ? ` ${label}` : "";
      text.textContent = `${phaseLabel} ${done}/${total} (${pct}%)${detail}`;
    }

    function showPlotDownload(show) {
      const btn = $("plotDownloadBtn");
      btn.classList.toggle("visible", !!show);
      btn.disabled = !show;
    }

    function syntenyDownloadFilename() {
      const b = lastMeta?.bounds;
      if (b?.contig != null && b.start != null && b.end != null) {
        const contig = String(b.contig).replace(/[^\w.-]+/g, "_");
        return `synteny_${contig}_${b.start}-${b.end}.png`;
      }
      return "synteny_extract.png";
    }

    function setPlotLoading(on, text) {
      const wrap = $("plotLoading");
      wrap.classList.toggle("visible", !!on);
      wrap.hidden = !on;
      wrap.setAttribute("aria-busy", String(!!on));
      if (text) $("plotLoadingText").textContent = text;
      if (on) {
        $("syntenyImg").classList.remove("visible");
        $("syntenyHint").hidden = true;
      } else {
        $("syntenyHint").hidden = false;
      }
    }

    function clearSyntenyPlot() {
      setPlotLoading(false);
      showPlotDownload(false);
      const img = $("syntenyImg");
      img.removeAttribute("src");
      img.classList.remove("visible");
      $("syntenyHint").textContent =
        "Run a query, then click Plot to render MUM extract synteny (matplotlib loads on first plot).";
      $("syntenyHint").hidden = false;
    }

    function getSelectedSeqIndices() {
      if (!manifest) return [];
      if (selectedSeqs) return [...selectedSeqs].sort((a, b) => a - b);
      return manifest.genomes.map((g) => g.seq_idx);
    }

    function getFilteredRows() {
      if (!lastRows || !manifest) return [];
      const keep = selectedSeqs ? selectedSeqs : new Set(manifest.genomes.map((g) => g.seq_idx));
      return lastRows.filter((r) => keep.has(r.seq_idx));
    }

    function fmtBp(n) {
      return Number.isFinite(n) ? Math.round(n).toLocaleString() : "";
    }

    function themeColor(name) {
      return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    }

    function updateThemeToggleLabel() {
      const dark = document.documentElement.dataset.theme === "dark";
      const btn = $("themeToggle");
      btn.setAttribute("aria-pressed", String(dark));
      btn.setAttribute("aria-label", dark ? "Dark mode on" : "Dark mode off");
      btn.title = dark ? "Dark mode on (click for light)" : "Light mode (click for dark)";
    }

    function applyTheme(theme) {
      document.documentElement.dataset.theme = theme;
      localStorage.setItem(THEME_KEY, theme);
      updateThemeToggleLabel();
      if ($("tabLengths").getAttribute("aria-selected") === "true") {
        renderLengthsView(getFilteredRows());
      }
    }

    updateThemeToggleLabel();
    $("themeToggle").addEventListener("click", () => {
      applyTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark");
    });

    /** Integer tick positions from 0 to maxY inclusive (for count axis). */
    function countAxisTicks(maxY, desiredSteps = 5) {
      if (maxY <= 0) return [0];
      const step = Math.max(1, Math.ceil(maxY / desiredSteps));
      const out = [];
      for (let v = 0; v <= maxY; v += step) out.push(v);
      if (out[out.length - 1] !== maxY) out.push(maxY);
      return out;
    }

    function renderLengthsView(rows) {
      const canvas = $("lengthsHistCanvas");
      const hint = $("lengthsHistHint");
      if (!canvas || !hint) return;
      const ctx = canvas.getContext("2d");
      const rr = Array.isArray(rows) ? rows : [];

      const panel = $("panelLengths");
      const rect = canvas.getBoundingClientRect();
      const cssW = Math.max(320, (panel && panel.offsetWidth) || rect.width || canvas.clientWidth || 800);
      const cssH = 240;
      const dpr = window.devicePixelRatio || 1;
      canvas.style.height = `${cssH}px`;
      canvas.width = Math.floor(cssW * dpr);
      canvas.height = Math.floor(cssH * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

      if (rr.length === 0) {
        hint.textContent = "No output. Run a query or select assemblies in the sequence list.";
        ctx.fillStyle = themeColor("--hist-bg");
        ctx.fillRect(0, 0, cssW, cssH);
        return;
      }

      const lengths = rr.map(regionLength);
      const n = lengths.length;
      const minL = Math.min(...lengths);
      const maxL = Math.max(...lengths);
      hint.textContent = `Region lengths (bp) for ${n.toLocaleString()} assembly row(s): min ${fmtBp(minL)}, max ${fmtBp(maxL)}.`;

      let numBins = Math.min(LENGTH_HIST_NUM_BINS, Math.max(1, n));
      const counts = new Array(numBins).fill(0);
      const lo = minL;
      const hi = maxL;
      let step;
      if (hi === lo) {
        numBins = 1;
        step = 1;
        counts[0] = n;
      } else {
        step = (hi - lo) / numBins;
        for (const L of lengths) {
          let i = Math.floor((L - lo) / step);
          if (i < 0) i = 0;
          if (i >= numBins) i = numBins - 1;
          counts[i]++;
        }
      }

      const maxCount = Math.max(1, ...counts);
      const yTicks = countAxisTicks(maxCount, 5);
      const margin = { L: 58, R: 12, T: 20, B: 36 };
      const pw = cssW - margin.L - margin.R;
      const ph = cssH - margin.T - margin.B;
      const plotLeft = margin.L;
      const plotRight = margin.L + pw;
      const plotTop = margin.T;
      const plotBot = margin.T + ph;

      ctx.fillStyle = themeColor("--hist-bg");
      ctx.fillRect(0, 0, cssW, cssH);

      ctx.strokeStyle = themeColor("--hist-grid");
      ctx.lineWidth = 1;
      for (const tick of yTicks) {
        const y = plotBot - (tick / maxCount) * ph;
        ctx.beginPath();
        ctx.moveTo(plotLeft, y);
        ctx.lineTo(plotRight, y);
        ctx.stroke();
      }

      const bw = pw / numBins;
      ctx.fillStyle = themeColor("--accent-bar");
      for (let i = 0; i < numBins; i++) {
        const h = (counts[i] / maxCount) * ph;
        const x = plotLeft + i * bw + 1;
        const y = plotBot - h;
        const wbar = Math.max(0, bw - 2);
        if (h > 0) ctx.fillRect(x, y, wbar, h);
      }

      ctx.strokeStyle = themeColor("--hist-axis");
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(plotLeft, plotTop);
      ctx.lineTo(plotLeft, plotBot);
      ctx.lineTo(plotRight, plotBot);
      ctx.stroke();

      ctx.fillStyle = themeColor("--hist-label");
      ctx.font = "12px system-ui,sans-serif";
      ctx.textBaseline = "middle";
      ctx.textAlign = "right";
      for (const tick of yTicks) {
        const y = plotBot - (tick / maxCount) * ph;
        ctx.fillText(String(tick.toLocaleString()), plotLeft - 6, y);
      }

      ctx.save();
      ctx.translate(14, plotTop + ph / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("Count", 0, 0);
      ctx.restore();

      ctx.textBaseline = "alphabetic";
      ctx.textAlign = "center";
      ctx.font = "12px system-ui,sans-serif";
      ctx.fillStyle = themeColor("--hist-label");
      const xTickY = cssH - 20;
      ctx.fillText(fmtBp(lo), plotLeft, xTickY);
      ctx.fillText(fmtBp(hi), plotRight, xTickY);
      if (hi !== lo) {
        const mid = (lo + hi) / 2;
        ctx.fillText(fmtBp(mid), plotLeft + pw / 2, xTickY);
      }
      ctx.font = "11px system-ui,sans-serif";
      ctx.fillStyle = themeColor("--hist-label-sub");
      ctx.fillText("Length (bp)", plotLeft + pw / 2, cssH - 6);
    }

    function setActiveTab(which) {
      const isTable = which === "table";
      const isBed = which === "bed";
      const isUnavail = which === "unavailable";
      const isLengths = which === "lengths";
      const isSynteny = which === "synteny";
      $("tabTable").setAttribute("aria-selected", String(isTable));
      $("tabBed").setAttribute("aria-selected", String(isBed));
      $("tabUnavailable").setAttribute("aria-selected", String(isUnavail));
      $("tabLengths").setAttribute("aria-selected", String(isLengths));
      $("tabSynteny").setAttribute("aria-selected", String(isSynteny));
      $("panelTable").classList.toggle("active", isTable);
      $("panelBed").classList.toggle("active", isBed);
      $("panelUnavailable").classList.toggle("active", isUnavail);
      $("panelLengths").classList.toggle("active", isLengths);
      $("panelSynteny").classList.toggle("active", isSynteny);
      if (isLengths) renderLengthsView(getFilteredRows());
    }

    $("tabTable").addEventListener("click", () => setActiveTab("table"));
    $("tabBed").addEventListener("click", () => setActiveTab("bed"));
    $("tabUnavailable").addEventListener("click", () => setActiveTab("unavailable"));
    $("tabLengths").addEventListener("click", () => setActiveTab("lengths"));
    $("tabSynteny").addEventListener("click", () => setActiveTab("synteny"));

    function renderIntervalsTable(rows) {
      const body = $("bedTableBody");
      body.innerHTML = "";
      const rr = Array.isArray(rows) ? rows : [];
      if (rr.length === 0) {
        const tr = document.createElement("tr");
        const td = document.createElement("td");
        td.colSpan = 5;
        td.className = "hint";
        td.textContent = "No output.";
        tr.appendChild(td);
        body.appendChild(tr);
        return;
      }
      for (const r of rr) {
        const genome = genomeDisplayForRow(r);
        const contig = contigDisplayForRow(r);
        const tr = document.createElement("tr");
        const tdG = document.createElement("td");
        tdG.textContent = genome;
        const tdC = document.createElement("td");
        tdC.className = "col-contig";
        tdC.textContent = contig;
        tdC.title = contig;
        const tdS = document.createElement("td");
        tdS.className = "mono";
        tdS.textContent = String(r?.start ?? "");
        const tdE = document.createElement("td");
        tdE.className = "mono";
        tdE.textContent = String(r?.end ?? "");
        const tdL = document.createElement("td");
        tdL.className = "mono";
        tdL.textContent = fmtBp(regionLength(r));
        tr.appendChild(tdG);
        tr.appendChild(tdC);
        tr.appendChild(tdS);
        tr.appendChild(tdE);
        tr.appendChild(tdL);
        body.appendChild(tr);
      }
    }

    function renderUnavailableTable(unavailable) {
      const body = $("unavailTableBody");
      body.innerHTML = "";
      const rows = Array.isArray(unavailable) ? unavailable : [];
      if (rows.length === 0) {
        const tr = document.createElement("tr");
        const td = document.createElement("td");
        td.colSpan = 3;
        td.className = "hint";
        td.textContent = "None.";
        tr.appendChild(td);
        body.appendChild(tr);
        return;
      }
      for (const r of rows) {
        const tr = document.createElement("tr");
        const asm = assemblyNameForSeq(r?.seq_idx);
        const genome = asm || prettyGenome(r?.label ?? `seq_${r?.seq_idx ?? ""}`);
        const c1 = prettyContig(String(r?.contig_a ?? ""));
        const c2 = prettyContig(String(r?.contig_b ?? ""));
        for (const v of [genome, c1, c2]) {
          const td = document.createElement("td");
          td.textContent = v;
          tr.appendChild(td);
        }
        body.appendChild(tr);
      }
    }

    const buildToSeqIdx = {
      "CHM13v2.0": 0,
      "GRCh38": 11,
    };
    const seqIdxToBuild = {
      0: "CHM13v2.0",
      11: "GRCh38",
    };
    const buildToGeneJson = {
      "CHM13v2.0": "annotations.genes.chm13.json",
      "GRCh38": "annotations.genes.grch38.json",
    };
    let syncingBuildGenome = false;
    let suppressGeneClear = false;
    let suppressRegionSync = false;

    function clearGeneIfUserEdit() {
      if (suppressGeneClear) return;
      if ($("geneQuery").value) $("geneQuery").value = "";
    }

    function clearIntervalFieldsIfGeneEdit() {
      if (suppressGeneClear) return;
      suppressRegionSync = true;
      try {
        $("region").value = "";
        $("start").value = "";
        $("end").value = "";
      } finally {
        queueMicrotask(() => {
          suppressRegionSync = false;
          refreshRegionValidation();
        });
      }
    }

    function setGeneStatus(msg) {
      $("geneStatus").textContent = msg || "";
    }

    async function loadGeneIndexForBuild(build) {
      const buildKey = String(build || "");
      const rel = buildToGeneJson[buildKey];
      if (!rel) {
        geneIndex = null;
        setGeneStatus("Unknown build.");
        return;
      }
      try {
        setGeneStatus("Loading genes…");
        const res = await fetch(new URL(rel, base));
        if (!res.ok) throw new Error(`fetch failed (${res.status})`);
        const parsed = await res.json();
        const idx = parsed?.[buildKey] ?? null;
        if (!idx || typeof idx !== "object") throw new Error("missing build key in JSON");
        geneIndex = idx;
        const n = Object.keys(geneIndex).length;
        setGeneStatus(`Loaded ${n.toLocaleString()} genes.`);
        fillGeneDatalist();
      } catch (e) {
        geneIndex = null;
        setGeneStatus(`Gene load error: ${String(e?.message ?? e)}`);
      }
    }

    async function loadGeneIndexFromAnnotation(build) {
      const tab = activeTab();
      const buildKey = String(build || "");
      const root = tab?.annotationJson;
      if (!root || typeof root !== "object") {
        geneIndex = null;
        setGeneStatus("No annotation loaded.");
        return;
      }
      const idx = root[buildKey];
      if (!idx || typeof idx !== "object") {
        geneIndex = null;
        setGeneStatus("Unknown build in annotation.");
        return;
      }
      geneIndex = idx;
      const n = Object.keys(geneIndex).length;
      setGeneStatus(`Loaded ${n.toLocaleString()} genes.`);
      fillGeneDatalist();
    }

    function fillGeneDatalist() {
      const dl = $("geneList");
      dl.innerHTML = "";
      if (!geneIndex) return;
      const genes = Object.keys(geneIndex);
      genes.sort();
      const frag = document.createDocumentFragment();
      for (const g of genes) {
        const opt = document.createElement("option");
        opt.value = g;
        frag.appendChild(opt);
      }
      dl.appendChild(frag);
    }

    function ensureContigValue(contig) {
      const sel = $("contig");
      const want = String(contig || "");
      const hasChr = want.startsWith("chr");
      const alt = hasChr ? want.slice(3) : `chr${want}`;
      for (const o of Array.from(sel.options)) {
        const v = String(o.value || "");
        const pc = prettyContig(v);
        if (v === want || v === alt || pc === want || pc === alt) return v;
      }
      return null;
    }

    function setGeneUiVisible(visible) {
      const show = !!visible;
      $("buildRow").style.display = show ? "" : "none";
      $("geneRow").style.display = show ? "" : "none";
      $("build").disabled = !show;
      $("geneQuery").disabled = !show;
      if (!show) {
        $("geneQuery").value = "";
        $("geneList").innerHTML = "";
        geneIndex = null;
        setGeneStatus("");
      }
    }

    function parseRegion(text) {
      const s = String(text || "").trim();
      // Accept UCSC-like contig:start-end with optional commas/underscores/spaces.
      const m = /^([^:]+)\s*:\s*([0-9][0-9,_]*)\s*-\s*([0-9][0-9,_]*)\s*$/.exec(s);
      if (!m) return { ok: false, error: "Expected CONTIG:START-END" };
      const contig = m[1].trim();
      const start = Number(String(m[2]).replace(/[,_]/g, ""));
      const end = Number(String(m[3]).replace(/[,_]/g, ""));
      if (!Number.isFinite(start) || !Number.isFinite(end)) return { ok: false, error: "Start/end must be numbers" };
      if (start < 0 || end < 0 || start >= end) {
        return { ok: false, error: "Need start < end (half-open: end is exclusive, same as shredtools -r)" };
      }
      return { ok: true, contig, start: Math.trunc(start), end: Math.trunc(end) };
    }

    function getGenome(seqIdx) {
      if (!manifest?.genomes) return null;
      return manifest.genomes.find((x) => x.seq_idx === seqIdx) ?? null;
    }

    function contigLength(seqIdx, contigValue) {
      const g = getGenome(seqIdx);
      if (!g?.contigs || !g?.contig_lengths) return null;
      const i = g.contigs.indexOf(contigValue);
      if (i < 0) return null;
      return g.contig_lengths[i];
    }

    function validateHalfOpenInterval(contig, start, end, L) {
      const c = String(contig || "");
      const s = Math.trunc(Number(start));
      const e = Math.trunc(Number(end));
      if (!Number.isFinite(L) || L < 0) {
        return { ok: false, error: "Contig length unknown." };
      }
      if (!(0 <= s && s < e && e <= L)) {
        return {
          ok: false,
          error:
            `Region ${s}-${e} is invalid for contig ${prettyContig(c)} with length ${fmtBp(L)} ` +
            "(expected half-open [start, end) in contig coordinates)",
        };
      }
      return { ok: true, start: s, end: e };
    }

    function resolveIntervalInputs() {
      if (!manifest) return { ok: false, error: "Still loading…" };
      const seqIdx = parseInt($("genome").value, 10);
      if (!Number.isFinite(seqIdx)) return { ok: false, error: "No genome selected." };

      const regionText = String($("region").value || "").trim();
      if (regionText) {
        const parsed = parseRegion(regionText);
        if (!parsed.ok) return { ok: false, error: parsed.error };
        const contigWanted = ensureContigValue(parsed.contig);
        if (!contigWanted) {
          return { ok: false, error: `Invalid contig for selected genome: ${parsed.contig}` };
        }
        return {
          ok: true,
          seqIdx,
          contig: contigWanted,
          start: parsed.start,
          end: parsed.end,
          source: "region",
        };
      }

      const contig = $("contig").value;
      const start = Number($("start").value);
      const end = Number($("end").value);
      if (!contig) return { ok: false, error: "No contig selected." };
      if (!Number.isFinite(start) || !Number.isFinite(end)) {
        return { ok: false, error: "Start/end must be numbers." };
      }
      if (start >= end) {
        return { ok: false, error: "Need start < end (half-open: end is exclusive)." };
      }
      return {
        ok: true,
        seqIdx,
        contig,
        start: Math.trunc(start),
        end: Math.trunc(end),
        source: "fields",
      };
    }

    function validateResolvedInterval() {
      const resolved = resolveIntervalInputs();
      if (!resolved.ok) return resolved;
      const L = contigLength(resolved.seqIdx, resolved.contig);
      if (L == null) {
        return { ok: false, error: "Contig not found for selected genome." };
      }
      const bounds = validateHalfOpenInterval(
        resolved.contig,
        resolved.start,
        resolved.end,
        L
      );
      if (!bounds.ok) return bounds;
      return {
        ok: true,
        seqIdx: resolved.seqIdx,
        contig: resolved.contig,
        start: bounds.start,
        end: bounds.end,
        contigLen: L,
      };
    }

    function setInputInvalid(el, invalid) {
      el.classList.toggle("input-invalid", !!invalid);
    }

    function updateContigLengthLabel(contigValue) {
      const el = $("contigLen");
      if (!manifest) {
        el.textContent = "";
        return;
      }
      const seqIdx = parseInt($("genome").value, 10);
      const contig = contigValue ?? $("contig").value;
      const L = contigLength(seqIdx, contig);
      el.textContent = L != null ? `${fmtBp(L)} bp` : "";
    }

    function refreshRegionValidation() {
      if (!manifest || $("start").disabled) return;
      const resolved = resolveIntervalInputs();
      const v = resolved.ok ? validateResolvedInterval() : resolved;
      const invalid = !v.ok;
      setInputInvalid($("start"), invalid);
      setInputInvalid($("end"), invalid);
      setInputInvalid($("region"), invalid && String($("region").value || "").trim() !== "");
      $("run").disabled = invalid;
      if (resolved.ok) {
        updateContigLengthLabel(resolved.contig);
      } else {
        updateContigLengthLabel();
      }
    }

    function updateRegionFromFields() {
      if (suppressRegionSync) return;
      const contigRaw = String($("contig").value || "");
      const contig = prettyContig(contigRaw);
      const start = String($("start").value || "").trim();
      const end = String($("end").value || "").trim();
      if (!contig || !start || !end) return;
      const sNum = Number(start);
      const eNum = Number(end);
      if (!Number.isFinite(sNum) || !Number.isFinite(eNum) || sNum >= eNum) return;
      suppressRegionSync = true;
      try {
        $("region").value = `${contig}:${Math.trunc(sNum)}-${Math.trunc(eNum)}`;
      } finally {
        queueMicrotask(() => {
          suppressRegionSync = false;
        });
      }
    }

    function maybeUpdateFieldsFromRegionLive() {
      if (suppressRegionSync) return;
      const txt = String($("region").value || "").trim();
      if (!txt) return;
      const parsed = parseRegion(txt);
      if (!parsed.ok) return; // allow partial typing without fighting the user
      const contigWanted = ensureContigValue(parsed.contig);
      if (!contigWanted) return;
      suppressRegionSync = true;
      try {
        $("contig").value = contigWanted;
        $("start").value = String(parsed.start);
        $("end").value = String(parsed.end);
      } finally {
        queueMicrotask(() => {
          suppressRegionSync = false;
        });
      }
    }

    function applyGene(gene) {
      const g = String(gene || "").trim();
      if (!g) return;
      if (!geneIndex) {
        status("Genes not loaded yet.");
        return;
      }
      const hits = geneIndex[g];
      if (!hits || hits.length === 0) {
        status("Gene not found.", g);
        return;
      }
      const hit = hits[0];
      let contigWanted = ensureContigValue(hit.contig);
      if (!contigWanted) {
        const found = findGenomeForContig(hit.contig);
        if (found) {
          syncingBuildGenome = true;
          try {
            $("genome").value = String(found.seqIdx);
            syncContigOptions();
          } finally {
            syncingBuildGenome = false;
          }
          contigWanted = found.contig;
        }
      }
      if (!contigWanted) {
        status("Contig not available for selected genome.", String(hit.contig));
        return;
      }
      suppressGeneClear = true;
      suppressRegionSync = true;
      try {
        $("contig").value = contigWanted;
        $("start").value = String(hit.start);
        $("end").value = String(hit.end);
      } finally {
        // let the current call stack finish (change/input events) before re-enabling clears
        queueMicrotask(() => {
          suppressGeneClear = false;
          suppressRegionSync = false;
          // Now that sync is re-enabled, update Region to match the gene-selected interval.
          updateRegionFromFields();
          refreshRegionValidation();
        });
      }
      status("Interval set from gene.", `${g} → ${contigWanted}:${hit.start}-${hit.end}`);
    }

    async function initLengthsUI() {
      status("Parsing lengths…", `building dropdowns (elapsed ${secs()})`);
      const jsonStr = String(await pyodide.runPythonAsync(`import app; app.describe_ui()`));
      manifest = JSON.parse(jsonStr);
      const tab = activeTab();
      if (tab?.kind === "custom" && tab.assemblyLabels?.length) {
        applyAssemblyLabelsToManifest(manifest, tab.assemblyLabels);
      }
      fillGenomeOptions();
      syncContigOptions();
      initSeqPicker();
      $("start").disabled = false;
      $("end").disabled = false;
      $("region").disabled = false;
      $("run").disabled = false;
      $("lucky").disabled = false;
      $("lengthsHint").textContent = "";
      refreshRegionValidation();
      status("Lengths loaded. Set interval and click Run.", `elapsed ${secs()}`);
      // Gene UI visibility + loading is driven by the selected genome (seq_idx).
      $("genome").dispatchEvent(new Event("change"));
    }

    async function loadPangenomeAndInitUI(key) {
      const opt = pangenomeList.find((o) => o.key === key);
      if (!opt) throw new Error(`Unknown pangenome: ${key}`);
      status("Switching pangenome…", opt.label);
      const prevPyKey = await getActivePangenomeKey();
      await setActivePangenomeKey(key);
      try {
        await initLengthsUI();
      } catch (e) {
        try {
          if (prevPyKey && prevPyKey !== key) await setActivePangenomeKey(prevPyKey);
        } catch (_) {
          /* ignore */
        }
        markQueryUiUnloaded("Pangenome switch failed (see status).");
        throw e;
      }
    }

    function updateSeqSummary() {
      const n = manifest?.genomes?.length ?? 0;
      const selN = selectedSeqs ? selectedSeqs.size : 0;
      const all = selN === n;
      $("seqSummary").textContent = all ? `Sequences: all (${n})` : `Sequences: ${selN}/${n}`;
    }

    function applySeqFilterToOutput() {
      if (!lastRows) return;
      const keep = selectedSeqs ? selectedSeqs : new Set(manifest.genomes.map((g) => g.seq_idx));
      const rows = lastRows.filter((r) => keep.has(r.seq_idx));
      const bedText = rows.map((r) => bedLineForRow(r)).join("");
      $("out").value = bedText;
      renderIntervalsTable(rows);
      if ($("tabLengths").getAttribute("aria-selected") === "true") {
        renderLengthsView(rows);
      }
    }

    function renderSeqList() {
      const q = String($("seqFilter").value || "").toLowerCase().trim();
      const list = $("seqList");
      list.innerHTML = "";
      const frag = document.createDocumentFragment();
      for (const g of manifest.genomes) {
        const label = prettyGenome(g.label);
        if (q && !label.toLowerCase().includes(q)) continue;
        const id = `seq_${g.seq_idx}`;
        const row = document.createElement("div");
        row.style.display = "flex";
        row.style.alignItems = "center";
        row.style.gap = "0.4rem";
        row.style.margin = "0.15rem 0";
        const cb = document.createElement("input");
        cb.type = "checkbox";
        cb.id = id;
        cb.checked = selectedSeqs.has(g.seq_idx);
        cb.addEventListener("change", () => {
          if (cb.checked) selectedSeqs.add(g.seq_idx);
          else selectedSeqs.delete(g.seq_idx);
          updateSeqSummary();
          applySeqFilterToOutput();
        });
        const lab = document.createElement("label");
        lab.htmlFor = id;
        lab.textContent = label;
        row.appendChild(cb);
        row.appendChild(lab);
        frag.appendChild(row);
      }
      list.appendChild(frag);
    }

    function initSeqPicker() {
      const n = manifest.genomes.length;
      selectedSeqs = new Set();
      for (let i = 0; i < n; i++) selectedSeqs.add(manifest.genomes[i].seq_idx);
      updateSeqSummary();
      renderSeqList();

      $("seqAll").onclick = () => {
        selectedSeqs = new Set(manifest.genomes.map((g) => g.seq_idx));
        updateSeqSummary();
        renderSeqList();
        applySeqFilterToOutput();
      };
      $("seqNone").onclick = () => {
        selectedSeqs = new Set();
        updateSeqSummary();
        renderSeqList();
        applySeqFilterToOutput();
      };
      $("seqFilter").oninput = () => renderSeqList();
    }

    function fillGenomeOptions() {
      const sel = $("genome");
      sel.innerHTML = "";
      for (const g of manifest.genomes) {
        const o = document.createElement("option");
        o.value = String(g.seq_idx);
        o.textContent = `${prettyGenome(g.label)} (seq ${g.seq_idx})`;
        sel.appendChild(o);
      }
      sel.disabled = false;
    }

    function syncContigOptions() {
      const gi = parseInt($("genome").value, 10);
      const g = manifest.genomes.find((x) => x.seq_idx === gi);
      const sel = $("contig");
      sel.innerHTML = "";
      for (const c of g.contigs) {
        const o = document.createElement("option");
        o.value = c;
        o.textContent = prettyContig(c);
        sel.appendChild(o);
      }
      sel.disabled = false;
      refreshRegionValidation();
    }

    $("pangenome").addEventListener("change", async () => {
      const key = $("pangenome").value;
      $("pangenome").disabled = true;
      $("genome").disabled = true;
      $("contig").disabled = true;
      $("run").disabled = true;
      $("lucky").disabled = true;
      $("start").disabled = true;
      $("end").disabled = true;
      $("region").disabled = true;
      lastRows = null;
      lastMeta = null;
      lastUnavailable = null;
      $("out").value = "";
      $("bounds").textContent = "";
      renderIntervalsTable([]);
      renderUnavailableTable([]);
      renderLengthsView([]);
      showPlotButton(false);
      showFastaDownloadButton(false);
      clearFastaSave();
      clearSyntenyPlot();
      try {
        await loadPangenomeAndInitUI(key);
        status("Pangenome switched.", `elapsed ${secs()}`);
      } catch (e) {
        $("lengthsHint").textContent = "Pangenome switch failed (see status).";
        status("Error", String(e?.message ?? e));
      } finally {
        $("pangenome").disabled = false;
      }
    });

    $("genome").addEventListener("change", () => {
      if (!manifest) return;
      syncContigOptions();
      refreshRegionValidation();
      if (selectedSeqs) renderSeqList();
      const seqIdx = parseInt($("genome").value, 10);
      const tab = activeTab();
      if (tab?.kind === "custom") {
        if (tab.annotationJson && !syncingBuildGenome) {
          setGeneUiVisible(true);
          // Keep current build selection; do not force-hide.
        } else if (!tab.annotationJson) {
          setGeneUiVisible(false);
        }
      } else {
        const build = seqIdxToBuild[seqIdx];
        const eligible = build != null;
        setGeneUiVisible(eligible);
        if (eligible && !syncingBuildGenome) {
          syncingBuildGenome = true;
          try {
            $("build").value = build;
          } finally {
            syncingBuildGenome = false;
          }
          // Load genes for this build (async) and clear any previous gene selection.
          $("geneQuery").value = "";
          loadGeneIndexForBuild(build);
        }
      }
      // Warm the bumbl index for this genome so subsequent Runs are faster.
      (async () => {
        const epoch = beginStatus();
        const w0 = performance.now();
        statusAt(epoch, "Loading genome index…", `seq ${seqIdx} (elapsed ${secs(w0)})`);
        try {
          await pyodide.runPythonAsync(`
import app
await app.warm_index(${seqIdx})
          `);
          statusAt(epoch, "Ready.", `index cached for seq ${seqIdx}`);
        } catch (e) {
          statusAt(epoch, "Error", String(e?.message ?? e));
        }
      })();
    });

    $("build").addEventListener("change", async () => {
      if (!manifest) return;
      const build = $("build").value;
      const tab = activeTab();
      if (tab?.kind === "custom") {
        $("geneQuery").value = "";
        await loadGeneIndexFromAnnotation(build);
        return;
      }
      const seqIdx = buildToSeqIdx[build];
      if (!syncingBuildGenome && Number.isFinite(seqIdx)) {
        syncingBuildGenome = true;
        $("genome").value = String(seqIdx);
        $("genome").dispatchEvent(new Event("change"));
        syncingBuildGenome = false;
      }
      $("geneQuery").value = "";
      await loadGeneIndexForBuild(build);
    });

    $("geneQuery").addEventListener("input", () => {
      if (!manifest) return;
      clearIntervalFieldsIfGeneEdit();
    });

    $("geneQuery").addEventListener("change", () => {
      if (!manifest) return;
      applyGene($("geneQuery").value);
    });

    $("region").addEventListener("input", () => {
      clearGeneIfUserEdit();
      maybeUpdateFieldsFromRegionLive();
      refreshRegionValidation();
    });
    $("contig").addEventListener("change", () => {
      clearGeneIfUserEdit();
      updateRegionFromFields();
      refreshRegionValidation();
    });
    $("start").addEventListener("input", () => {
      clearGeneIfUserEdit();
      updateRegionFromFields();
      refreshRegionValidation();
    });
    $("end").addEventListener("input", () => {
      clearGeneIfUserEdit();
      updateRegionFromFields();
      refreshRegionValidation();
    });

    function randomInt(maxExclusive) {
      return Math.floor(Math.random() * maxExclusive);
    }

    function pickWeightedIndex(weights) {
      let total = 0;
      for (const w of weights) total += w;
      if (!(total > 0)) return 0;
      let r = Math.random() * total;
      for (let i = 0; i < weights.length; i++) {
        r -= weights[i];
        if (r <= 0) return i;
      }
      return weights.length - 1;
    }

    function pickLuckyInterval() {
      const genomes = manifest?.genomes;
      if (!genomes?.length) return null;
      const g = genomes[randomInt(genomes.length)];
      const lengths = (g.contig_lengths || []).map((x) => Number(x));
      if (!g.contigs?.length || lengths.length !== g.contigs.length) return null;
      // Length-weighted so short scaffolds are not over-sampled.
      const minUseful = 1000;
      const weights = lengths.map((L) => (Number.isFinite(L) && L >= minUseful ? L : 0));
      const ci =
        weights.some((w) => w > 0) ? pickWeightedIndex(weights) : randomInt(g.contigs.length);
      const contig = g.contigs[ci];
      const L = lengths[ci];
      if (!(Number.isFinite(L) && L >= 1)) return null;
      const maxWin = Math.min(100_000, L);
      const minWin = Math.min(10_000, maxWin);
      const win = minWin + randomInt(maxWin - minWin + 1);
      const start = randomInt(L - win + 1);
      return {
        seqIdx: g.seq_idx,
        contig,
        start,
        end: start + win,
        label: prettyGenome(g.label),
      };
    }

    function applyLuckyInterval(pick) {
      if (!pick) return false;
      suppressGeneClear = true;
      suppressRegionSync = true;
      try {
        if ($("geneQuery").value) $("geneQuery").value = "";
        $("genome").value = String(pick.seqIdx);
        $("genome").dispatchEvent(new Event("change"));
        $("contig").value = pick.contig;
        $("start").value = String(pick.start);
        $("end").value = String(pick.end);
        $("region").value = `${prettyContig(pick.contig)}:${pick.start}-${pick.end}`;
      } finally {
        suppressGeneClear = false;
        suppressRegionSync = false;
      }
      refreshRegionValidation();
      return true;
    }

    $("lucky").addEventListener("click", () => {
      if (!manifest) {
        status("Still loading…");
        return;
      }
      const pick = pickLuckyInterval();
      if (!pick || !applyLuckyInterval(pick)) {
        status("Could not pick a random region.");
        return;
      }
      status(
        "Feeling lucky…",
        `${pick.label} ${prettyContig(pick.contig)}:${pick.start}-${pick.end}`
      );
      $("run").click();
    });

    $("run").addEventListener("click", async () => {
      if (!manifest) {
        status("Still loading…");
        return;
      }
      const v = validateResolvedInterval();
      if (!v.ok) {
        refreshRegionValidation();
        status("Invalid interval.", v.error);
        return;
      }
      const { seqIdx, contig, start, end } = v;
      const rangeStr = `${contig}:${start}-${end}`;
      const epoch = beginStatus();
      const q0 = performance.now();
      statusAt(epoch, "Querying…", `HTTP Range reads to S3 (query ${secs(q0)})`);
      $("out").value = "";
      $("bounds").textContent = "";
      showPlotButton(false);
      showFastaDownloadButton(false);
      clearFastaSave();
      clearSyntenyPlot();
      lastRows = null;
      renderIntervalsTable([]);
      renderUnavailableTable([]);
      renderLengthsView([]);
      try {
        const resultJson = String(
          await pyodide.runPythonAsync(`
import app
await app.run_with_bounds(${seqIdx}, ${JSON.stringify(rangeStr)}, None)
          `)
        );
        const result = JSON.parse(resultJson);
        if (result.error) {
          statusAt(epoch, "Could not run query.", result.error);
          $("out").value = "";
          showPlotButton(false);
          showFastaDownloadButton(false);
          clearFastaSave();
          clearSyntenyPlot();
          renderIntervalsTable([]);
          renderUnavailableTable([]);
          renderLengthsView([]);
          lastRows = null;
          lastMeta = null;
          lastUnavailable = null;
          return;
        }
        lastRows = result.rows ?? [];
        lastMeta = result;
        lastUnavailable = result.unavailable ?? [];
        applySeqFilterToOutput();
        renderUnavailableTable(lastUnavailable);
        let mumSliceSummary = "";
        if (result.mum_slices) {
          mumSliceSummary = `Loaded ${result.mum_slices.mums} MUMs, ${result.mum_slices.chunks} chunks`;
          console.log(`[shredtools] mum slices: ${mumSliceSummary}`);
        }
        const b = [];
        if (result.bounds) {
          b.push(`Bounds (selected genome): ${result.bounds.contig}:${result.bounds.start}-${result.bounds.end}`);
        }
        if (result.margins) {
          b.push(`Margins: left ${result.margins.left}, right ${result.margins.right}`);
        }
        $("bounds").textContent = b.join(" | ");
        statusAt(epoch, "Done.", `query ${secs(q0)}${mumSliceSummary ? ` | ${mumSliceSummary}` : ""}`);
        showPlotButton(true);
        showFastaDownloadButton(true);
        clearFastaSave();
        clearSyntenyPlot();
        setActiveTab("table");
      } catch (err) {
        console.error(err);
        const msg = err?.message ?? String(err);
        statusAt(epoch, "Error", msg);
        $("out").value = "";
        showPlotButton(false);
        showFastaDownloadButton(false);
        clearFastaSave();
        clearSyntenyPlot();
        renderIntervalsTable([]);
        renderUnavailableTable([]);
        lastRows = null;
        lastMeta = null;
        lastUnavailable = null;
        renderLengthsView([]);
      }
    });

    $("fastaSaveBtn").addEventListener("click", () => {
      if (!lastFastaBlob) return;
      triggerDownload(lastFastaBlob, lastFastaFilename);
      status("FASTA saved.", lastFastaFilename);
    });

    $("fastaCancelBtn").addEventListener("click", () => {
      if (!fastaBuildAbort) return;
      fastaBuildAbort.abort();
      status("Cancelling FASTA build…", "");
    });

    $("fastaDownloadBtn").addEventListener("click", async () => {
      const rows = getFilteredRows();
      if (!rows.length) {
        status("Nothing to download.", "Select at least one assembly with results.");
        return;
      }
      const btn = $("fastaDownloadBtn");
      const runBtn = $("run");
      const gzipCb = $("fastaGzip");
      btn.disabled = true;
      runBtn.disabled = true;
      gzipCb.disabled = true;
      const useGzip = gzipCb.checked;
      const q0 = performance.now();
      fastaBuildAbort = new AbortController();
      showFastaCancelButton(true);
      setFastaProgress(true, 0, rows.length, "", "start");
      status("Building FASTA…", `0/${rows.length}`);
      try {
        const { file, written, skipped, mode, gzip } = await buildMultifasta(rows, {
          gzip: useGzip,
          signal: fastaBuildAbort.signal,
          onProgress({ done, total, label, phase }) {
            setFastaProgress(true, done, total, label, phase);
            const phaseLabel = phase === "write" ? "Writing" : "Fetching";
            status("Building FASTA…", `${phaseLabel} ${done}/${total}${label ? `: ${label}` : ""}`);
          },
        });
        setFastaProgress(false);
        if (!written) {
          const detail = skipped.map((s) => `${s.label}: ${s.reason}`).join("; ");
          status("FASTA build failed.", detail || "No sequences fetched.");
          return;
        }
        clearFastaSave();
        lastFastaGzip = !!gzip;
        const mime = lastFastaGzip ? "application/gzip" : "text/plain";
        lastFastaBlob = file instanceof Blob ? file : new Blob([file], { type: mime });
        lastFastaFilename = fastaDownloadFilename(lastMeta?.bounds, lastFastaGzip);
        showFastaSaveButton(true);
        const fmt = lastFastaGzip ? "FASTA.gz" : "FASTA";
        let sub = `Built ${written} sequence${written === 1 ? "" : "s"} (${mode}${lastFastaGzip ? ", gzip" : ""}) in ${secs(q0)} — click Save ${fmt}`;
        if (skipped.length) {
          sub += ` | skipped ${skipped.length}: ${skipped.map((s) => s.label).join(", ")}`;
        }
        status("FASTA ready.", sub);
      } catch (err) {
        console.error(err);
        setFastaProgress(false);
        if (isFastaBuildCancelled(err)) {
          status("FASTA build cancelled.", "");
        } else {
          status("FASTA download error", err?.message ?? String(err));
        }
      } finally {
        showFastaCancelButton(false);
        fastaBuildAbort = null;
        btn.disabled = false;
        runBtn.disabled = false;
        refreshFastaDownloadButton();
      }
    });

    $("plotDownloadBtn").addEventListener("click", () => {
      const src = $("syntenyImg").src;
      if (!src) return;
      const a = document.createElement("a");
      a.href = src;
      a.download = syntenyDownloadFilename();
      a.rel = "noopener";
      document.body.appendChild(a);
      a.click();
      a.remove();
    });

    $("plotBtn").addEventListener("click", async () => {
      const q0 = performance.now();
      const btn = $("plotBtn");
      btn.disabled = true;
      setActiveTab("synteny");
      setPlotLoading(true, "Loading plot modules…");
      status("Plotting…", "loading plot modules if needed");
      try {
        await ensurePlotModules();
        setPlotLoading(true, "Loading matplotlib…");
        status("Plotting…", "loading matplotlib if needed");
        await ensureMpl();
        setPlotLoading(true, "Fetching MUM bins for plot…");
        status("Plotting…", `matplotlib ready (${secs(q0)})`);
        const seqs = getSelectedSeqIndices();
        const dark = document.documentElement.dataset.theme === "dark";
        const plotJson = String(
          await pyodide.runPythonAsync(`
import inspect
import app
_result = app.plot_extract_png(${JSON.stringify(seqs)}, ${dark ? "True" : "False"})
if inspect.isawaitable(_result):
    _result = await _result
_result
          `)
        );
        const plotResult = JSON.parse(plotJson);
        if (plotResult.error) {
          setPlotLoading(false);
          status("Plot error", plotResult.error);
          return;
        }
        setPlotLoading(false);
        $("syntenyImg").src = "data:image/png;base64," + plotResult.png_b64;
        $("syntenyImg").classList.add("visible");
        showPlotDownload(true);
        $("syntenyHint").textContent =
          `Extract synteny: ${plotResult.n_mums} MUMs across ${plotResult.n_rows} assemblies.`;
        $("syntenyHint").hidden = false;
        status("Plot done.", `render ${secs(q0)}`);
      } catch (err) {
        setPlotLoading(false);
        status("Plot error", err?.message ?? String(err));
      } finally {
        btn.disabled = false;
      }
    });

    // Load default pangenome (HPRCr2 enhanced): refresh UI from the already-loaded lengths bundle.
    (async () => {
      try {
        pangenomeList = JSON.parse(
          String(await pyodide.runPythonAsync("import app; app.pangenome_options_json()"))
        );
        const ps = $("pangenome");
        ps.innerHTML = "";
        for (const o of pangenomeList) {
          const op = document.createElement("option");
          op.value = o.key;
          op.textContent = o.label;
          ps.appendChild(op);
        }
        ps.value = "hprcv2_enhanced";
        ps.disabled = true;
        renderSidebar();
        await restoreCustomDatasetsOnBoot();
        void ensureEmojiSelect();
        if (!isRoutedPath()) {
          updateBrowserRoute(HUMAN_TAB_ID, { replace: true });
        }
        const initialTabId = tabIdFromRoute(routeFromLocation());
        await switchToTab(initialTabId, { fromPopstate: true });
      } catch (e) {
        $("lengthsHint").textContent = "Could not load pangenome_lengths.json or init UI.";
        status("Error", String(e?.message ?? e));
      } finally {
        $("pangenome").disabled = false;
      }
    })();

    $("sidebarAdd").addEventListener("click", () => {
      void ensureEmojiSelect();
      void switchToTab(ADD_TAB_ID);
    });

    $("addBumbl").addEventListener("input", () => {
      if (lengthsUrlTouched) return;
      $("addLengths").value = defaultLengthsUrl($("addBumbl").value);
    });
    $("addLengths").addEventListener("input", () => {
      lengthsUrlTouched = true;
    });

    $("addDatasetCancel").addEventListener("click", () => {
      void switchToTab(HUMAN_TAB_ID);
    });

    $("addDatasetForm").addEventListener("submit", async (ev) => {
      ev.preventDefault();
      const name = String($("addName").value || "").trim();
      let emoji = String($("addEmoji").value || "").trim();
      if (!emoji) {
        // Allow typing/pasting an emoji directly without clicking a result.
        const typed = String($("addEmojiSearch").value || "").trim();
        const firstToken = typed.split(/\s+/)[0] || "";
        if (firstToken && emojiCatalog.some((x) => x.e === firstToken)) {
          emoji = firstToken;
        } else if (firstToken && /\p{Extended_Pictographic}/u.test(firstToken)) {
          emoji = firstToken;
        }
      }
      const bumblUrl = String($("addBumbl").value || "").trim();
      const lengthsUrl = String($("addLengths").value || "").trim() || defaultLengthsUrl(bumblUrl);
      const annotationUrl = String($("addAnnotation").value || "").trim();
      const errEl = $("addDatasetError");
      errEl.textContent = "";
      if (!name) {
        errEl.textContent = "Name is required.";
        return;
      }
      if (!bumblUrl) {
        errEl.textContent = "Bumbl URL is required.";
        return;
      }
      if (!lengthsUrl) {
        errEl.textContent = "Lengths URL is required.";
        return;
      }
      const submitBtn = $("addDatasetSubmit");
      submitBtn.disabled = true;
      try {
        await addCustomDatasetFromForm({ name, bumblUrl, lengthsUrl, annotationUrl, emoji });
        $("addDatasetForm").reset();
        setSelectedEmoji("");
        lengthsUrlTouched = false;
        errEl.textContent = "";
      } catch (e) {
        errEl.textContent = String(e?.message ?? e);
        status("Could not load dataset.", String(e?.message ?? e));
      } finally {
        submitBtn.disabled = false;
      }
    });

    $("removeCustomDataset").addEventListener("click", () => {
      void removeActiveCustomDataset();
    });

    window.addEventListener("popstate", () => {
      if (activeTabId === ADD_TAB_ID) return;
      const tabId = tabIdFromRoute(routeFromLocation());
      if (tabId !== activeTabId) void switchToTab(tabId, { fromPopstate: true });
    });
