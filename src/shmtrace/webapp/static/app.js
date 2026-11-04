/* shmtrace desktop app: set up a project, run it, read the results, export. */
"use strict";

const S = {                       // everything the page knows
  options: null, reference: null, contaminant: null, germlineVcf: null,
  samples: [],                    // {name, path, group, control}
  sheet: null,                    // the loaded sample sheet (provenance for the Sanger QC)
  job: null, result: null, table: "samples", sort: null, sortDir: 1,
};
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmtInt = (n) => Number(n).toLocaleString("en-GB");
const plural = (n, word) => `${fmtInt(n)} ${word}${n === 1 ? "" : "s"}`;

/* ── tabs ─────────────────────────────────────────────────── */
function showTab(name) {
  for (const t of ["setup", "results", "reads", "tables", "calls", "export", "switch", "repair", "power", "screen"]) {
    $("page-" + t).classList.toggle("active", t === name);
    $("tab-" + t).setAttribute("aria-selected", String(t === name));
  }
  if (name === "results") drawFigures();
  if (name === "reads") rdInit();
  if (name === "calls") rvInit();
  if (name === "repair") rpInit();
  if (name === "power") pwInit();
}
function enableResultTabs(on) {
  for (const t of ["results", "tables", "calls", "export", "repair", "screen"]) $("tab-" + t).disabled = !on;
  $("tab-reads").disabled = !(on && S.result && S.result.has_alignments);
}

/* ── start up ─────────────────────────────────────────────── */
async function init() {
  const r = await fetch("/api/options");
  S.options = await r.json();
  $("version").textContent = "v" + S.options.version;
  fill($("region"), S.options.presets, "generic_amplicon");
  $("region").onchange = () => $("region-hint").textContent = S.options.preset_help[$("region").value] || "";
  $("region").onchange();
  buildReferencePicker();
  fill($("platform"), S.options.platforms, "ont");
  fill($("mode"), S.options.modes, "recommended");
  $("mode").onchange = () => $("mode-hint").textContent = S.options.mode_help[$("mode").value] || "";
  $("mode").onchange();
  $("platform").onchange = () => {
    if ($("platform").value === "sanger" && $("mode").value === "recommended") {
      $("mode-hint").textContent = "Sanger clones: every high-quality mismatch counts as a mutation.";
    } else $("mode").onchange();
  };
  dropZone($("drop-ref"), false, (f) => upload(f).then(([x]) => {
    $("ref-pick").value = "";
    S.reference = x; $("ref-name").textContent = "Reference: " + x.name;
  }));
  dropZone($("drop-vcf"), false, (f) => upload(f).then(([x]) => {
    if (!x) return;
    S.germlineVcf = x; $("vcf-name").textContent = "Germline VCF: " + x.name;
  }));
  dropZone($("drop-contam"), false, (f) => upload(f).then(([x]) => {
    S.contaminant = x; $("contam-name").textContent = "Contaminant: " + x.name;
  }));
  dropZone($("drop-samples"), true, (f) => upload(f).then(addSamples), "pick-folder");
  dropZone($("drop-model"), false, (f) => upload(f).then(([x]) => {
    S.model = x; $("model-name").textContent = "Targeting model: " + x.name;
  }));
  if (new URLSearchParams(location.search).has("desktop")) $("save-local").style.display = "";
  swInit();
}

function fill(sel, items, def) {
  sel.innerHTML = items.map(x => `<option${x === def ? " selected" : ""}>${esc(x)}</option>`).join("");
}

function buildReferencePicker() {
  const refs = S.options.references || [];
  if (!refs.length) return;
  $("ref-pick").innerHTML = '<option value="">Use my own file…</option>'
    + refs.map(r => `<option value="${esc(r.id)}">${esc(r.name)} (${r.length} bp)</option>`).join("");
  $("ref-pick").onchange = () => {
    const r = refs.find(x => x.id === $("ref-pick").value);
    S.reference = r ? { name: r.name, path: r.path } : null;
    $("ref-name").textContent = r ? "Reference: " + r.name : "No reference chosen.";
    // a contaminant added for the previous built-in reference goes with it; one you chose stays
    if (S.contaminant && S.contaminant.auto) {
      S.contaminant = null; $("contam-name").textContent = "";
    }
    const c = (S.options.contaminants || []).find(x => r && r.id.startsWith(x.for));
    if (c && !S.contaminant) {
      S.contaminant = { name: c.name, path: c.path, auto: true };
      $("contam-name").textContent = "Contaminant: " + c.name + " (added for you)";
    }
  };
}

/* ── files ────────────────────────────────────────────────── */
function dropZone(el, multiple, onFiles, folderBtn) {
  const input = document.createElement("input");
  input.type = "file"; input.multiple = multiple; input.style.display = "none";
  input.dataset.role = el.id + "-files";
  el.after(input);
  el.onclick = () => input.click();
  input.onchange = () => { if (input.files.length) onFiles([...input.files]); input.value = ""; };
  if (folderBtn) {
    // a separate input, because webkitdirectory turns a file picker into a
    // folder picker and the two cannot be the same control
    const dirInput = document.createElement("input");
    dirInput.type = "file"; dirInput.multiple = true; dirInput.style.display = "none";
    dirInput.setAttribute("webkitdirectory", "");
    dirInput.dataset.role = el.id + "-folder";
    el.after(dirInput);
    $(folderBtn).onclick = () => dirInput.click();
    dirInput.onchange = () => { if (dirInput.files.length) onFiles([...dirInput.files]); dirInput.value = ""; };
  }
  el.ondragover = (e) => { e.preventDefault(); el.classList.add("over"); };
  el.ondragleave = () => el.classList.remove("over");
  el.ondrop = async (e) => {
    e.preventDefault(); el.classList.remove("over");
    const files = await filesFromDrop(e.dataTransfer);
    if (files.length) onFiles(files);
  };
}

/** Files from a drop, each tagged with its path inside a dropped folder. */
async function filesFromDrop(dt) {
  const out = [];
  const items = [...(dt.items || [])].map(i => i.webkitGetAsEntry && i.webkitGetAsEntry()).filter(Boolean);
  if (!items.length) return [...dt.files];
  const walk = (entry, prefix) => new Promise((res) => {
    if (entry.isFile) entry.file((f) => { f._rel = prefix + entry.name; out.push(f); res(); });
    else if (entry.isDirectory) {
      const rd = entry.createReader();
      const all = [];
      const readMore = () => rd.readEntries(async (es) => {
        if (es.length) { all.push(...es); readMore(); }       // readEntries pages
        else { for (const e of all) await walk(e, prefix + entry.name + "/"); res(); }
      });
      readMore();
    } else res();
  });
  for (const it of items) await walk(it, "");
  return out;
}

async function upload(files) {
  const fd = new FormData();
  for (const f of files) {
    fd.append("files", f, f.name);
    fd.append("paths", f._rel || f.webkitRelativePath || f.name);
  }
  const r = await fetch("/api/files", { method: "POST", body: fd });
  if (!r.ok) { error(await r.text()); return []; }
  return (await r.json()).files;
}

/* A library is a sample, or a control of some kind. All but a baseline are taken to be free of
   mutations; a baseline (the cells before the experiment) is the control only when nothing
   mutation-free is given (project.py). "Control" alone is the plain yes of a sample sheet. */
const ROLES = [["no", "sample"], ["yes", "control"], ["aid_null", "AID-null"],
               ["empty_vector", "empty vector"], ["plasmid", "plasmid"],
               ["germline_tissue", "germline tissue"], ["baseline", "baseline"]];
function roleOf(s) { return !s.control ? "no" : (s.control_kind || "yes"); }
function setRole(i, v) {
  const s = S.samples[i];
  s.control = v !== "no";
  s.control_kind = (v === "no" || v === "yes") ? "" : v;
  renderSamples();
}

/* ── samples ──────────────────────────────────────────────── */
// "sample_R1.fastq.gz", "sample.R1.fq", "sample_1.fastq": the 1 is the mate number
const R1R2 = /([._-])(R)?1(?=(\.|_|-|$))/i;
const strip = (n) => n.replace(/\.(gz|bz2)$/i, "")
  .replace(/\.(fastq|fq|fasta|fa|fna|fas|fsa|bam|sam|cram|ab1|abi|seq|gb|gbk)$/i, "");

function addSamples(files) {
  // A dropped folder is one sample: a folder of Sanger traces is one library,
  // and the engine reads a directory as one set of reads.
  const folders = new Map();
  const loose = [];
  for (const f of files) {
    const rel = f.rel || "";
    const top = rel.includes("/") ? rel.split("/")[0] : null;
    if (top) {
      if (!folders.has(top)) folders.set(top, { name: top, path: f.folder || f.path, n: 0 });
      folders.get(top).n++;
    } else loose.push(f);
  }
  for (const d of folders.values()) {
    // the server put the folder's files under one directory; use that directory
    const dir = d.path.replace(/\/[^/]*$/, "");
    S.samples.push({ name: d.name, path: dir, group: "", control: false, note: `${d.n} files` });
  }
  // pair Illumina R1/R2 by name; the engine takes "R1+R2"
  const used = new Set();
  for (const f of loose) {
    if (used.has(f)) continue;
    const base = strip(f.name);
    const m = base.match(R1R2);
    if (m) {
      const mateName = base.slice(0, m.index) + m[1] + (m[2] || "") + "2"
        + base.slice(m.index + m[0].length);
      const mate = loose.find(g => g !== f && !used.has(g) && strip(g.name) === mateName);
      if (mate) {
        used.add(f); used.add(mate);
        const label = (base.slice(0, m.index) + base.slice(m.index + m[0].length))
          .replace(/[._-]+$/, "") || base;
        S.samples.push({ name: label, path: f.path + "+" + mate.path,
                         group: "", control: false, note: "paired reads" });
        continue;
      }
    }
    used.add(f);
    S.samples.push({ name: strip(f.name), path: f.path, group: "", control: false });
  }
  autoGroup();
}

function clearSamples() { S.samples = []; S.sheet = null; $("sheet-msg").innerHTML = ""; renderSamples(); }

// A sample sheet names every file's sample, group and control status, and where
// each clone came from (order, well, label), which the Sanger QC reports with
// every flag. See shmtrace.samplesheet.
async function loadSheet(given) {
  if (given) $("sheet-path").value = given;
  const path = $("sheet-path").value.trim();
  if (!path) return;
  const r = await fetch("/api/samplesheet", {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ path }),
  });
  const d = await r.json();
  if (!r.ok) {
    const probs = (d.detail && d.detail.problems) || [d.detail || "could not read the sheet"];
    $("sheet-msg").innerHTML = `<div class="msg err" style="margin-top:8px">The sheet has problems:
      <ul>${probs.map(p => `<li>${esc(p)}</li>`).join("")}</ul></div>`;
    return;
  }
  S.samples = d.samples.map(s => ({ name: s.name, path: s.path, group: s.group, control: s.control,
                                    control_kind: s.control_kind || "", note: `${s.path.length} file(s) from the sheet` }));
  S.sheet = d.sheet;
  const ex = d.excluded.length ? `; ${d.excluded.length} row(s) excluded by the sheet` : "";
  const pooledFiles = Object.keys(d.pooled || {});
  $("demux-opts").hidden = !pooledFiles.length;
  const pooledMsg = pooledFiles.length ? `<div class="msg" style="margin-top:6px">Pooled run: ${pooledFiles.length} file(s) shared by
    ${Object.values(d.pooled).reduce((n, x) => n + Object.keys(x).length, 0)} samples, told apart by barcode. The reads are
    split when you run; give the primers next to the barcodes below if they are not the reference ends.</div>` : "";
  $("sheet-msg").innerHTML = `<div class="hint" style="margin-top:6px">${d.samples.length} samples
    from ${esc(path)}${ex}</div>` + pooledMsg + d.warnings.map(w => `<div class="msg" style="margin-top:6px">${esc(w)}</div>`).join("");
  renderSamples();
}

function autoGroup() {
  // Guess the group from the file name: drop the parts every sample shares
  // (the run ID, the chemistry, "sub500") and the bare numbers, and what is
  // left is what distinguishes this sample. "DHJGKB_5_AID_KO_1_C57BL6_sub500"
  // next to its siblings becomes "AID KO".
  const split = (s) => s.split(/[_.\-\s]+/).filter(Boolean);
  const all = S.samples.map(s => split(s.name));
  const common = new Set(all.length > 1
    ? all[0].filter(tok => all.every(toks => toks.includes(tok)))
    : []);
  S.samples.forEach((s, i) => {
    if (s.group) return;
    let toks = all[i].filter(t => !common.has(t) && !/^\d+$/.test(t));
    toks = toks.filter(t => /[A-Za-z]/.test(t));
    if (!toks.length) toks = all[i].filter(t => /[A-Za-z]/.test(t));
    // a trailing replicate number ("Wildtype2", "rep3") is not part of the group
    const label = toks.map(t => t.replace(/\d+$/, "")).filter(Boolean).join(" ");
    s.group = label || s.name;
  });
  renderSamples();
}

function markControls() {
  // whole words only: "ko" must not match "Knockin" or a mouse called "Koala"
  const re = /(^|[^a-z])(ko|knockout|empty|empty_?vector|ev|aicda|aid_?ko|null|control|ctrl|germline|untransfected|mock|uninfected)([^a-z]|$)/i;
  for (const s of S.samples) if (re.test(s.name)) s.control = true;
  renderSamples();
}

function renderSamples() {
  const groups = [...new Set(S.samples.map(s => s.group).filter(Boolean))];
  $("samples").innerHTML = S.samples.map((s, i) => `
    <div class="sample">
      <span class="nm" title="${esc(s.path)}">${esc(s.name)}${s.note ? ` <span class="pill ctrl">${esc(s.note)}</span>` : ""}</span>
      <input list="group-list" value="${esc(s.group)}" oninput="S.samples[${i}].group=this.value"
             placeholder="group">
      <select style="width:auto" title="Is this library a control, and of what kind?"
              onchange="setRole(${i}, this.value)">
        ${ROLES.map(([v, t]) => `<option value="${v}" ${roleOf(s) === v ? "selected" : ""}>${t}</option>`).join("")}
      </select>
      <button class="x" title="Remove" onclick="S.samples.splice(${i},1);renderSamples()">✕</button>
    </div>`).join("")
    + `<datalist id="group-list">${groups.map(g => `<option>${esc(g)}</option>`).join("")}</datalist>`;
  const n = S.samples.length, c = S.samples.filter(s => s.control).length;
  $("sample-count").textContent = `${n} sample${n === 1 ? "" : "s"}${c ? `, ${c} control${c === 1 ? "" : "s"}` : ""}`;
}

/* ── run ──────────────────────────────────────────────────── */
function error(msg) {
  $("run-error").innerHTML = `<div class="msg err">${esc(msg)}</div>`;
}

async function runProject() {
  $("run-error").innerHTML = "";
  if (!S.reference) return error("Choose a reference FASTA first.");
  if (!S.samples.length) return error("Add at least one sample.");
  const names = S.samples.map(s => s.name);
  if (new Set(names).size !== names.length) return error("Two samples have the same name.");
  const over = {};
  if ($("window").value.trim()) {
    const m = $("window").value.match(/(\d+)\s*-\s*(\d+)/);
    if (!m) return error("The analysis window should look like 100-630.");
    over.analysis_window = [+m[1], +m[2]];
  }
  if ($("maxp").value.trim()) over.max_pvalue = +$("maxp").value;
  over.jobs = +$("jobs").value;
  if ($("call-indels").checked) over.exclude_indels = false;
  if ($("short-reads").checked) over.rescue_short_reads = true;
  if ($("umi").value !== "off") {
    over.umi = $("umi").value;
    over.umi_ends = $("umi-ends").value;
    over.umi_length = +$("umi-length").value || 18;
    over.umi_min_reads = +$("umi-min").value || 3;
    over.umi_offset = +$("umi-offset").value || 0;
    const L = $("umi-left").value.trim().toUpperCase(), R = $("umi-right").value.trim().toUpperCase();
    if (/[^ACGTRYSWKMBDHVN]/.test(L + R)) return error("UMI anchors may only contain DNA letters and IUPAC codes.");
    if (L) over.umi_left_anchor = L;
    if (R) over.umi_right_anchor = R;
  }
  const saturation = $("saturation").checked;

  const project = {
    name: $("proj-name").value || "SHMtrace project",
    region: $("region").value, platform: $("platform").value, mode: $("mode").value,
    reference: S.reference.path,
    contaminants: S.contaminant ? [S.contaminant.path] : [],
    depth: $("depth").value ? +$("depth").value : null,
    saturation,
    browse: $("keep-alignment").checked,
    targeting_model: S.model ? S.model.path : null,
    overrides: over,
    samples: S.samples.map(s => ({ name: s.name, path: s.path, group: s.group || s.name, control: s.control,
                                   ...(s.control && s.control_kind ? { control_kind: s.control_kind } : {}) })),
    germline_vcf: S.germlineVcf ? S.germlineVcf.path : null,
    sample_sheet: S.sheet || null,
    demux: { left_anchor: $("demux-left").value.trim(), right_anchor: $("demux-right").value.trim() },
    // the bundled sample data, still under its own name: the server does not keep the run (main.py)
    ...(S.example && $("proj-name").value === S.example ? { example: true } : {}),
  };
  $("run-btn").disabled = true;
  $("run-status").style.display = "";
  const r = await fetch("/api/run", {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ project }),
  });
  if (!r.ok) { $("run-btn").disabled = false; return error((await r.json()).detail || "could not start"); }
  S.job = (await r.json()).job;
  poll();
}

async function poll() {
  const r = await fetch("/api/job/" + S.job);
  const j = await r.json();
  $("run-bar").style.width = Math.round(100 * j.progress) + "%";
  $("run-msg").textContent = j.message;
  if (j.state === "done") return loadResult();
  if (j.state === "error") { $("run-btn").disabled = false; return error(j.error); }
  setTimeout(poll, 700);
}

async function loadResult() {
  S.result = await (await fetch("/api/result/" + S.job)).json();
  $("run-btn").disabled = false;
  enableResultTabs(true);
  buildHeadline();
  buildTablePicker();
  buildExport();
  showTab("results");
  if (typeof pjAfterLoad === "function") pjAfterLoad();
}

/* ── results ──────────────────────────────────────────────── */
const MS_BLUE = "#1F6FB4", MS_RED = "#C0392B";
const PALETTE = [MS_BLUE, MS_RED, "#3C7D4F", "#B07A1E", "#6A4C93", "#2A7F8F",
  "#8C564B", "#C2185B", "#5D6D7E", "#7F8C1E"];
// blue for the test groups, red for the control, as the paper's figures do
function colourOf(g) {
  const r = S.result;
  const ctrl = new Set((r.groups || []).filter(x => x.control).map(x => x.group));
  const spare = PALETTE.slice(2);
  const seen = {};
  let used = false, k = 0;
  for (const name of r.group_order) {
    if (ctrl.has(name)) seen[name] = MS_RED;
    else if (!used) { seen[name] = MS_BLUE; used = true; }
    else seen[name] = spare[k++ % spare.length];
  }
  return seen[g] || MS_BLUE;
}
// a sample with nothing called has no logarithm: it is shown at the density one
// single mutation would have given it, the smallest difference it could resolve
function valueOf(s, key) {
  const v = s[key];
  if (v === 0 && key === "unique_density_x1e5" && s.one_event_density_x1e5)
    return { v: s.one_event_density_x1e5, censored: true };
  return { v, censored: false };
}
const geoMean = (xs) => Math.pow(10, xs.reduce((a, v) => a + Math.log10(v), 0) / xs.length);
function pFor(r, group, key, log) {
  const row = (r.groups || []).find(x => x.group === group);
  if (!row || !row.compared_with) return null;
  const want = (log ? `${key}_p_welch_log_vs_` : `${key}_p_welch_vs_`);
  for (const k of Object.keys(row)) if (k.startsWith(want) && row[k] != null)
    return { p: row[k], ctrl: row.compared_with };
  for (const k of Object.keys(row)) if (k.startsWith(`${key}_p_welch_vs_`) && row[k] != null)
    return { p: row[k], ctrl: row.compared_with };
  return null;
}
const dark = () => matchMedia("(prefers-color-scheme: dark)").matches;
const layoutBase = () => ({
  paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
  font: { color: dark() ? "#e8eaee" : "#15181d", size: 12 },
  margin: { l: 60, r: 16, t: 10, b: 50 }, showlegend: false,
  xaxis: { gridcolor: dark() ? "#2b313b" : "#e9ecf1", zeroline: false },
  yaxis: { gridcolor: dark() ? "#2b313b" : "#e9ecf1", zeroline: false },
});
const CONFIG = {
  displaylogo: false, responsive: true, displayModeBar: "hover",
  modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
  toImageButtonOptions: { scale: 3, filename: "shmtrace" },
};
// a legend above the plot needs room, and the hover toolbar must not sit on it
const withLegend = (L) => ({
  ...L, showlegend: true, margin: { ...L.margin, t: 46 },
  legend: { orientation: "h", y: 1.16, x: 0, yanchor: "bottom" },
  modebar: { orientation: "v" },
});

function buildHeadline() {
  const r = S.result, p = r.params;
  const ctrls = r.samples.filter(s => s.control);
  const key = "unique_density_x1e5";
  const geo = (g) => {
    const v = r.samples.filter(s => s.group === g && s[key] != null)
      .map(s => valueOf(s, key).v).filter(x => x > 0);
    return v.length ? geoMean(v) : null;
  };
  // the verdict: every non-control group against the control, strongest first
  const ctrlGroup = (r.groups || []).find(x => x.control);
  const lines = (r.groups || []).filter(x => !x.control && x.compared_with).map(x => {
    const st = pFor(r, x.group, key, true);
    const a = geo(x.group), b = geo(x.compared_with);
    const fold = (a != null && b) ? a / b : null;
    return { g: x.group, ctrl: x.compared_with, a, b, fold, p: st ? st.p : null, n: x.n };
  }).sort((u, v) => (v.fold || 0) - (u.fold || 0));

  const verdict = lines.map(l => {
    const strong = l.p != null && l.p < 0.05 && l.fold && l.fold >= 3;
    const word = strong ? "clearly above" : (l.fold && l.fold >= 2 ? "above"
                 : (l.fold != null && l.fold <= 0.5 ? "below" : "not separated from"));
    // 73-fold, 2.4-fold, 0.43-fold: whole numbers only where the decimals no longer matter
    const fx = (f) => f >= 10 ? f.toFixed(0) : (f >= 1 ? f.toFixed(1) : f.toFixed(2));
    return `<div class="verdict ${strong ? "yes" : "maybe"}">
      <b>${esc(l.g)}</b> is ${word} <b>${esc(l.ctrl)}</b>:
      ${l.a == null ? "—" : l.a.toFixed(1)} vs ${l.b == null ? "—" : l.b.toFixed(2)} per 10⁵ bp
      ${l.fold ? `· <b>${fx(l.fold)}-fold</b>` : ""}
      ${l.p != null ? `· P = ${l.p < 1e-4 ? l.p.toExponential(1) : l.p.toFixed(4)}
        <span class="hint">(Welch on log₁₀ densities, n = ${l.n} vs ${ctrlGroup ? ctrlGroup.n : "?"})</span>` : ""}
    </div>`;
  }).join("");

  // the names the set-up guide uses; the internal ones stay in a tooltip for anyone matching files
  const regionName = (typeof REGION_TITLES !== "undefined" && REGION_TITLES[p.region]) || p.region;
  const platformName = (typeof PLATFORM_TITLES !== "undefined" && PLATFORM_TITLES[p.platform]) || p.platform;
  const modeHelp = (S.options.mode_help || {})[p.mode] || "";
  $("headline").innerHTML = `
    <div class="pj-head"><h2>${esc(r.name)}</h2>
      <span class="hint" id="pj-saved"></span>
      <button class="btn small ghost" onclick="pjSave()"
              title="The run, its review marks and the reads Check calls draws, as one .shmtrace file">Save project file</button></div>
    ${verdict || '<p class="hint">No control group, so there is nothing to compare against.</p>'}
    <dl class="kv">
      <dt>Region</dt><dd><span title="preset ${esc(p.region)}">${esc(regionName)}</span> · reference ${esc(p.reference_name)} ·
        positions ${p.analysis_window[0]}–${p.analysis_window[1]} of ${fmtInt(p.reference_length)} bp</dd>
      <dt>Data</dt><dd>${esc(platformName)} ·
        <span title="${esc(`${modeHelp} Mode ${p.mode}, error model ${p.error_model}.`.trim())}">${esc(p.mode)} calling</span>${p.depth ? " · " + fmtInt(p.depth) + (p.umi && p.umi !== "off" ? " molecules per sample (UMI consensus)" : " reads per sample") : ""}</dd>
      <dt>Samples</dt><dd>${plural(r.samples.length, "sample")} in ${plural(r.group_order.length, "group")},
        ${ctrls.length ? `${ctrls.length} of them ${ctrls.length === 1 ? "a control" : "controls"}` : "no controls"}</dd>
    </dl>
    ${(r.notes || []).map(n => `<div class="msg${n.startsWith("UMI consensus:") ? "" : " err"}" style="margin-top:10px">${esc(n)}</div>`).join("")}`;
}

// each chart's numbers as a tidy CSV (plotdata.py), to draw it another way
const FIG_DATA = {
  "f-density": "density_by_group.csv", "f-freq": "mutation_frequency_by_group.csv",
  "f-matrix": "substitution_matrix.csv", "f-spectrum": "spectrum_by_sample.csv",
  "f-load": "mutations_per_read.csv", "f-map": "positional_by_sample.csv",
  "f-mutated": "mutated_reads_by_group.csv", "f-hot": "hotspot_targeting.csv",
  "f-motifs": "motif_targeting.csv",
};
function figShell(id, title, note) {
  const data = FIG_DATA[id];
  return `<div class="fig"><div class="cap"><h4>${esc(title)}</h4>
      <span class="hint">${esc(note || "")}</span>
      ${data ? `<button class="btn small ghost fig-data" onclick="dl('plot_data/${data}')"
        title="The numbers behind this chart, one row per value, to graph however you like (${data})">Data</button>` : ""}</div>
    <div class="plot" id="${id}"></div><p class="hint fig-note" id="${id}-note" hidden></p></div>`;
}
// a chart drawn from a handful of mutations shows noise as a pattern: below this many
// calls a group (or, for a per-sample rate, a sample) is left out and named underneath
const MIN_CALLS = 20, MIN_SAMPLE_CALLS = 10;
function figNote(id, text) {
  const el = $(id + "-note");
  if (el) { el.textContent = text || ""; el.hidden = !text; }
}
// draws ahead of a layout change (a tab shown, a window narrowed) and would keep the old width
const plotResizer = typeof ResizeObserver !== "undefined"
  ? new ResizeObserver(es => es.forEach(e => e.target.data && Plotly.Plots.resize(e.target))) : null;

function drawFigures() {
  if (!S.result || $("figs").dataset.job === S.job) return;
  $("figs").dataset.job = S.job;
  $("figs").innerHTML =
    figShell("f-density", "Mutation density",
      "unique mutations per 10⁵ bp; one dot per sample, bar is the geometric mean ± geometric SD")
    + figShell("f-freq", "Mutation frequency", "mutations per 10³ bases sequenced")
    + figShell("f-matrix", "From → to substitutions", "unique mutations, and % of the group's total")
    + figShell("f-spectrum", "Substitution spectrum", "group mean ± SD, every sample shown")
    + figShell("f-load", "Mutations per read", "distribution per group")
    + figShell("f-map", "Position along the amplicon", "% of reads mutated at each position")
    + figShell("f-mutated", "Reads carrying a mutation", "% of reads with at least one")
    + figShell("f-hot", "AID hotspot targeting", "% of G:C mutations in WRC/GYW")
    + figShell("f-motifs", "Motif targeting",
      "% of G:C (AID) or A:T (Pol η) mutations in each motif; group mean ± SD, every sample shown; grey: chance");
  plotDots("f-density", "unique_density_x1e5", "Unique mutations / 10⁵ bp");
  plotDots("f-freq", "mutation_frequency_x1e3", "Mutations / 10³ bp");
  plotMatrix();
  plotSpectrum();
  plotLoad();
  plotMap();
  plotDots("f-mutated", "mutated_reads_pct", "% of reads mutated", { pct: true });
  plotDots("f-hot", "pct_GC_in_WRC_GYW", "% of G:C mutations at WRC/GYW",
           { pct: true, countKey: "gc_mutations", countWord: "G:C mutations" });
  plotMotifs();
  if (plotResizer) document.querySelectorAll("#figs .plot").forEach(el => plotResizer.observe(el));
}

// pct: a percentage, drawn on a linear 0–100 axis. countKey: the calls a sample's value rests on;
// a sample with fewer than MIN_SAMPLE_CALLS is left out and named under the chart.
function plotDots(el, key, title, opts = {}) {
  const r = S.result, traces = [], shapes = [], annotations = [];
  const per = {}, left = [];
  r.group_order.forEach(g => {
    per[g] = r.samples.filter(s => s.group === g && s[key] != null).filter(s => {
      if (!opts.countKey || (s[opts.countKey] || 0) >= MIN_SAMPLE_CALLS) return true;
      left.push(`${s.sample} (${s[opts.countKey] || 0})`);
      return false;
    }).map(s => ({ name: s.sample, ...valueOf(s, key) }));
  });
  figNote(el, left.length ? `Not shown, fewer than ${MIN_SAMPLE_CALLS} ${opts.countWord}: ${left.join(", ")}.` : "");
  const all = Object.values(per).flat().map(x => x.v);
  if (!all.length) {
    if (left.length) $(el).style.minHeight = "0";      // keep the card: the note says why it is empty
    else $(el).parentElement.style.display = "none";
    return;
  }
  const pos = all.filter(v => v > 0);
  // the test follows the data's spread, whatever the axis: a percentage spanning a 30-fold range is
  // still compared on the log scale, as the groups table does, and drawn on a linear 0–100 axis
  const logTest = pos.length > 0 && (Math.max(...pos) / Math.min(...pos) > 8 || all.some(v => v <= 0));
  const log = logTest && !opts.pct;
  const floor = pos.length ? Math.min(...pos) / 3 : 1e-3;
  const centre = {};

  r.group_order.forEach(g => {
    const rows = per[g];
    if (!rows.length) return;
    const vals = rows.map(x => Math.max(x.v, log ? floor : -Infinity));
    let m, lo, hi;
    if (log) {
      const ls = vals.map(Math.log10);
      m = ls.reduce((a, b) => a + b, 0) / ls.length;
      const sd = ls.length > 1 ? Math.sqrt(ls.reduce((a, v) => a + (v - m) ** 2, 0) / (ls.length - 1)) : 0;
      [m, lo, hi] = [10 ** m, 10 ** (m - sd), 10 ** (m + sd)];
    } else {
      m = vals.reduce((a, b) => a + b, 0) / vals.length;
      const sd = vals.length > 1 ? Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / (vals.length - 1)) : 0;
      [lo, hi] = [m - sd, m + sd];
    }
    // the fold beside the P is on the scale of the test (zeros at the floor, as a log axis draws them)
    centre[g] = log || !logTest ? m : geoMean(vals.map(v => Math.max(v, floor)));
    traces.push({                       // the mean, drawn as a wide horizontal tick
      type: "scatter", mode: "markers", x: [g], y: [m], showlegend: false,
      marker: { symbol: "line-ew", size: 44, color: colourOf(g),
                line: { color: colourOf(g), width: 3 } },
      error_y: { type: "data", symmetric: false, array: [hi - m], arrayminus: [m - lo],
                 color: colourOf(g), thickness: 1.6, width: 9 },
      hovertemplate: `${g}: ${log ? "geometric mean" : "mean"} %{y:.3g}<extra></extra>`,
    });
    traces.push({
      type: "scatter", mode: "markers", x: rows.map(() => g), y: vals,
      text: rows.map(x => x.name + (x.censored ? " (no mutations called)" : "")),
      hovertemplate: "%{text}: %{y:.3g}<extra></extra>", showlegend: false,
      marker: {
        size: 11, color: rows.map(x => x.censored ? (dark() ? "#1b1f26" : "#ffffff") : colourOf(g)),
        line: { color: colourOf(g), width: 2 },
      },
      jitter: 0.5, pointpos: 0,
    });
  });

  // the group difference, written where the eye already is
  const idx = (g) => r.group_order.indexOf(g);
  let step = 0;
  const top = Math.max(...all.map(v => Math.max(v, log ? floor : v)));
  r.group_order.forEach(g => {
    const st = pFor(r, g, key, logTest);
    if (!st || idx(st.ctrl) < 0 || !(per[g] || []).length || !(per[st.ctrl] || []).length) return;
    const a2 = Math.min(idx(g), idx(st.ctrl)), b2 = Math.max(idx(g), idx(st.ctrl));
    const y = log ? top * Math.pow(2.0, 1 + step) : top * (1.12 + 0.14 * step);
    const yl = log ? y / 1.3 : y - top * 0.04;
    // the bracket is a trace, not a shape: a trace reads the log axis correctly
    traces.push({
      type: "scatter", mode: "lines", showlegend: false, hoverinfo: "skip", cliponaxis: false,
      x: [r.group_order[a2], r.group_order[a2], r.group_order[b2], r.group_order[b2]],
      y: [yl, y, y, yl], line: { color: dark() ? "#9aa3b2" : "#555", width: 1 },
    });
    const fold = centre[st.ctrl] ? centre[g] / centre[st.ctrl] : null;
    annotations.push({
      x: (a2 + b2) / 2, y: log ? Math.log10(y * 1.2) : y + top * 0.03, xref: "x", yref: "y",
      text: `P = ${st.p < 1e-4 ? st.p.toExponential(1) : st.p.toFixed(4)}` +
            (fold && fold >= 1.5 ? ` (~${fold.toFixed(0)}×)` : ""),
      showarrow: false, font: { size: 12 },
    });
    step++;
  });

  const L = { ...layoutBase(), shapes, annotations, margin: { l: 64, r: 16, t: 34, b: 44 } };
  L.yaxis = { ...L.yaxis, title, type: log ? "log" : "linear" };
  L.xaxis = { ...L.xaxis, type: "category" };
  // a label at each power of ten only: Plotly's default also labels the 2s and 5s between them
  if (log) L.yaxis = { ...L.yaxis, range: [Math.log10(floor / 2), Math.log10(top * Math.pow(2.0, step + 1.6))],
                       dtick: 1, tickformat: "~g" };
  if (opts.pct) L.yaxis.range = [0, Math.max(100, step ? top * (1.12 + 0.14 * step) * 1.08 : 0)];
  Plotly.newPlot(el, traces, L, CONFIG);
}

function plotSpectrum() {
  const r = S.result;
  const subs = [];
  for (const a of "ATCG") for (const b of "ATCG") if (a !== b) subs.push(a + ">" + b);
  const counts = {};
  r.group_order.forEach(g => counts[g] = Object.values(r.spectrum_group[g] || {}).reduce((a, b) => a + b, 0));
  // 1 of 4 calls is "25%": a group with a handful of calls has no spectrum
  let shown = r.group_order.filter(g => counts[g] >= MIN_CALLS);
  const dropped = r.group_order.filter(g => counts[g] < MIN_CALLS);
  if (!shown.length) shown = r.group_order;

  const traces = [];
  shown.forEach(g => {
    const rows = r.samples.filter(s => s.group === g)
      .map(s => r.spectrum_sample[s.sample] || {})
      .filter(c => Object.values(c).reduce((a, b) => a + b, 0) > 0);
    const pcts = subs.map(k => rows.map(c => 100 * (c[k] || 0) /
      Object.values(c).reduce((a, b) => a + b, 0)));
    const means = pcts.map(v => v.length ? v.reduce((a, b) => a + b, 0) / v.length : 0);
    const sds = pcts.map((v, i) => v.length > 1
      ? Math.sqrt(v.reduce((a, x) => a + (x - means[i]) ** 2, 0) / (v.length - 1)) : 0);
    traces.push({
      type: "bar", name: g, x: subs.map(k => k.replace(">", "→")), y: means,
      marker: { color: colourOf(g) },
      error_y: { type: "data", array: sds, color: colourOf(g), thickness: 1.2, width: 4 },
      hovertemplate: `${g} %{x}: %{y:.1f}%<extra></extra>`,
    });
    // every sample shown as its own open point, as the paper does
    const xs = [], ys = [], tx = [];
    subs.forEach((k, i) => pcts[i].forEach((v, j) => {
      xs.push(k.replace(">", "→")); ys.push(v); tx.push(rows[j] ? g : g);
    }));
    traces.push({
      type: "scatter", mode: "markers", x: xs, y: ys, showlegend: false,
      marker: { size: 6, color: "rgba(0,0,0,0)", line: { color: colourOf(g), width: 1.4 } },
      hovertemplate: `${g} %{x}: %{y:.1f}%<extra></extra>`,
    });
  });
  const L = withLegend({ ...layoutBase(), barmode: "group" });
  L.yaxis.title = "% of unique mutations";
  figNote("f-spectrum", dropped.length && dropped.length < r.group_order.length
    ? "Not shown, too few calls for a spectrum: " + dropped.map(g => `${g} (n = ${counts[g]})`).join(", ") + "." : "");
  Plotly.newPlot("f-spectrum", traces, L, CONFIG);
}

function plotMatrix() {
  const r = S.result, bb = ["A", "T", "C", "G"];
  const rows = bb.slice().reverse();   // drawn bottom-up, so A ends up on top
  const total = (g) => Object.values(r.spectrum_group[g] || {}).reduce((a, b) => a + b, 0);
  // the same rule as the spectrum: two calls make a "100%" cell, not a pattern
  let shown = r.group_order.filter(g => total(g) >= MIN_CALLS);
  const dropped = r.group_order.filter(g => total(g) < MIN_CALLS);
  if (!shown.length) shown = r.group_order;
  figNote("f-matrix", dropped.length && shown.length < r.group_order.length
    ? "Not shown, too few calls for a matrix: " + dropped.map(g => `${g} (n = ${total(g)})`).join(", ") + "." : "");
  const n = shown.length;
  // light: white to blue; dark: the panel colour to blue, so an empty cell is not a white box
  const scale = dark() ? [[0, "#232a35"], [1, "#3d8fe0"]] : [[0, "#ffffff"], [0.5, "#6baed6"], [1, "#08306b"]];
  const traces = [], annotations = [];
  shown.forEach((g, gi) => {
    const c = r.spectrum_group[g] || {};
    const tot = total(g);
    const pct = (a, b) => a === b ? null : 100 * (c[`${a}>${b}`] || 0) / (tot || 1);
    const z = rows.map(a => bb.map(b => pct(a, b)));
    const zmax = Math.max(20, ...z.flat().map(v => v || 0));
    traces.push({
      type: "heatmap", z, x: bb, y: rows, xaxis: `x${gi + 1}`, yaxis: `y${gi + 1}`,
      colorscale: scale, showscale: false, zmin: 0, zmax,
      xgap: 2, ygap: 2,
      hovertemplate: `${g}  %{y}→%{x}: %{z:.1f}%<extra></extra>`,
    });
    rows.forEach(a => bb.forEach(b => {
      const v = a === b ? null : (c[`${a}>${b}`] || 0);
      annotations.push({
        xref: `x${gi + 1}`, yref: `y${gi + 1}`, x: b, y: a, showarrow: false,
        text: v === null ? "—"
          : `<b>${v}</b><br><span style="font-size:9px">${(pct(a, b)).toFixed(1)}%</span>`,
        font: { size: 12, color: dark() || (pct(a, b) || 0) > 0.55 * zmax ? "#fff" : "#222" },
      });
    }));
    annotations.push({
      xref: `x${gi + 1} domain`, yref: "paper", x: 0.5, y: 1.18, showarrow: false,
      text: `<b>${esc(g)}</b>  (n = ${tot})`, font: { size: 13 },
    });
  });
  const base = layoutBase();
  const L = { ...base, annotations, height: 330, showlegend: false,
    grid: { rows: 1, columns: n, pattern: "independent", xgap: 0.22 },
    margin: { l: 40, r: 14, t: 72, b: 40 } };
  for (let i = 1; i <= n; i++) {
    L[`xaxis${i}`] = { ...base.xaxis, side: "top", type: "category",
      title: i === 1 ? "To" : "", showgrid: false, ticks: "" };
    L[`yaxis${i}`] = { ...base.yaxis, type: "category",
      title: i === 1 ? "From" : "", showgrid: false, ticks: "" };
  }
  Plotly.newPlot("f-matrix", traces, L, CONFIG);
}

const LOAD_CLASSES = [["0", 0, 0], ["1–3", 1, 3], ["4–10", 4, 10], ["11–20", 11, 20], ["20+", 21, 1e9]];
const LOAD_SHADES = ["#F4F7FB", "#D6E6F4", "#A8CBE8", "#5B9BD5", "#1F4E79"];
// on a dark page the unmutated share is the darkest slice, so most of a control's donut is not a white ring
const LOAD_SHADES_DARK = ["#2b313b", "#29476a", "#33679b", "#4f8fd0", "#9cc3ee"];

function plotLoad() {
  const r = S.result, byGroup = {};
  const gmap = Object.fromEntries(r.samples.map(s => [s.sample, s.group]));
  for (const [s, hist] of Object.entries(r.mutation_load)) {
    const g = gmap[s]; byGroup[g] = byGroup[g] || {};
    for (const [k, n] of Object.entries(hist)) byGroup[g][+k] = (byGroup[g][+k] || 0) + n;
  }
  const groups = r.group_order.filter(g => byGroup[g]);
  const traces = groups.map((g, i) => {
    const vals = LOAD_CLASSES.map(([, a, b]) => Object.entries(byGroup[g])
      .reduce((acc, [k, n]) => acc + (+k >= a && +k <= b ? n : 0), 0));
    const tot = vals.reduce((a, b) => a + b, 0);
    return {
      type: "pie", hole: 0.52, values: vals, labels: LOAD_CLASSES.map(c => c[0]),
      name: g, sort: false, direction: "clockwise", rotation: 0,
      domain: { row: 0, column: i },
      marker: { colors: dark() ? LOAD_SHADES_DARK : LOAD_SHADES, line: { color: dark() ? "#1b1f26" : "#333", width: 1 } },
      textinfo: "percent", texttemplate: "%{percent:.0%}",
      textposition: "inside", insidetextorientation: "horizontal",
      hovertemplate: `${g}: %{label} mutations — %{value} reads (%{percent})<extra></extra>`,
      title: { text: `<b>${esc(g)}</b><br><span style="font-size:11px">n = ${tot.toLocaleString()}</span>`, position: "bottom center" },
    };
  });
  const L = withLegend({ ...layoutBase(), grid: { rows: 1, columns: groups.length } });
  L.legend = { ...L.legend, title: { text: "Mutations per read  " } };
  Plotly.newPlot("f-load", traces, L, CONFIG);
}

// every motif family side by side: what share of the calls at its bases it holds in each group,
// against what chance gives (its share of those bases in the window). A sample is left out of a
// family when it has fewer than MIN_SAMPLE_CALLS calls at the family's bases.
function plotMotifs() {
  const r = S.result, fams = r.motif_families || [], el = "f-motifs";
  if (!fams.length) { $(el).parentElement.style.display = "none"; return; }
  const labels = fams.map(f => f.label), traces = [], left = new Set();
  const groups = r.group_order, w = 0.8 / Math.max(groups.length, 1);
  groups.forEach((g, gi) => {
    const mean = [], sd = [], xs = [], ys = [], text = [];
    fams.forEach((f, fi) => {
      const nKey = f.bases === "GC" ? "gc_mutations" : "at_mutations", key = `pct_${f.bases}_in_${f.key}`;
      const vals = r.samples.filter(s => s.group === g && s[key] != null).filter(s => {
        if ((s[nKey] || 0) >= MIN_SAMPLE_CALLS) return true;
        left.add(s.sample);
        return false;
      });
      const v = vals.map(s => s[key]);
      const m = v.length ? v.reduce((a, b) => a + b, 0) / v.length : null;
      mean.push(m);
      sd.push(v.length > 1 ? Math.sqrt(v.reduce((a, x) => a + (x - m) ** 2, 0) / (v.length - 1)) : 0);
      vals.forEach(s => { xs.push(fi - 0.4 + w * (gi + 0.5)); ys.push(s[key]); text.push(`${s.sample}, ${f.label}`); });
    });
    const x = fams.map((_f, fi) => fi - 0.4 + w * (gi + 0.5));
    traces.push({ type: "bar", name: g, x, y: mean, width: w * 0.9, marker: { color: colourOf(g), opacity: 0.35 },
      error_y: { type: "data", array: sd, color: colourOf(g), thickness: 1.4, width: 4 },
      customdata: labels, textposition: "none",
      hovertemplate: `${esc(g)}, %{customdata}: mean %{y:.1f}%<extra></extra>` });
    traces.push({ type: "scatter", mode: "markers", showlegend: false, x: xs, y: ys, text,
      marker: { size: 6, color: colourOf(g), line: { color: dark() ? "#1b1f26" : "#fff", width: 0.8 } },
      hovertemplate: "%{text}: %{y:.1f}%<extra></extra>" });
  });
  const bg = r.motif_background || {};
  traces.push({ type: "scatter", mode: "markers", name: "chance", x: fams.map((_f, i) => i),
    y: fams.map(f => bg[f.key]), text: fams.map(f => `${f.label}: ${f.what}`),
    marker: { symbol: "line-ew", size: 40, color: "#9aa3b2", line: { color: "#9aa3b2", width: 2.5 } },
    hovertemplate: "%{text}<br>chance %{y:.1f}%<extra></extra>" });
  figNote(el, left.size ? `Left out of a motif when they have fewer than ${MIN_SAMPLE_CALLS} calls at its bases: ${[...left].join(", ")}.` : "");
  const L = { ...withLegend(layoutBase()), barmode: "overlay", margin: { l: 64, r: 16, t: 34, b: 56 } };
  L.xaxis = { ...L.xaxis, tickmode: "array", tickvals: fams.map((_f, i) => i), ticktext: labels, range: [-0.6, fams.length - 0.4] };
  L.yaxis = { ...L.yaxis, title: "% of mutations at its bases", range: [0, 105] };
  Plotly.newPlot(el, traces, L, CONFIG);
}

function plotMap() {
  const r = S.result;
  const bySample = {};
  for (const x of r.positional_sample || []) (bySample[x.sample] = bySample[x.sample] || []).push(x);
  const traces = [];
  r.group_order.forEach(g => {
    const rows = r.positional.filter(x => x.group === g).sort((a, b) => a.pos - b.pos);
    const sparse = rows.length < 20;
    // one faint trace per sample, under the group line
    r.samples.filter(s => s.group === g).forEach(s => {
      const pts = (bySample[s.sample] || []).sort((a, b) => a.pos - b.pos);
      if (!pts.length) return;
      traces.push({
        type: "scatter", mode: sparse ? "markers" : "lines", showlegend: false,
        x: pts.map(x => x.pos), y: pts.map(x => x.frequency_pct),
        line: { color: colourOf(g), width: 0.8 }, opacity: 0.45,
        marker: { color: colourOf(g), size: 4 },
        hovertemplate: `${esc(s.sample)} %{x}: %{y:.2f}%<extra></extra>`,
      });
    });
    if (sparse) {                      // a line between four points invents a trend
      const xs = [], ys = [];
      for (const x of rows) { xs.push(x.pos, x.pos, null); ys.push(0, x.frequency_pct, null); }
      traces.push({ type: "scatter", mode: "lines", name: g, x: xs, y: ys,
        line: { color: colourOf(g), width: 2 }, hoverinfo: "skip" });
      traces.push({ type: "scatter", mode: "markers", showlegend: false,
        x: rows.map(x => x.pos), y: rows.map(x => x.frequency_pct),
        marker: { color: colourOf(g), size: 7 },
        hovertemplate: `${g} %{x}: %{y:.2f}%<extra></extra>` });
    } else {
      traces.push({
        type: "scatter", mode: "lines", name: g,
        x: rows.map(x => x.pos), y: rows.map(x => x.frequency_pct),
        line: { color: colourOf(g), width: 1.8 },
        hovertemplate: `${g} %{x}: %{y:.2f}%<extra></extra>`,
      });
    }
  });
  // each motif family's sites as ticks along the axis; WRC/GYW shown, the rest a legend click away
  const sites = r.motif_sites || { WRC_GYW: r.hotspots };
  (r.motif_families || [{ key: "WRC_GYW", label: "WRC/GYW" }]).forEach(f => {
    const xs = sites[f.key] || [];
    if (!xs.length) return;
    traces.push({
      type: "scatter", mode: "markers", name: f.label, x: xs, y: xs.map(() => 0),
      visible: f.key === "WRC_GYW" ? true : "legendonly",
      marker: { color: f.bases === "AT" ? "#c9a227" : "#9aa3b2", size: 4, symbol: "line-ns-open" },
      hovertemplate: `${f.label} %{x}<extra></extra>`,
    });
  });
  const L = { ...withLegend(layoutBase()), hovermode: "closest" };
  L.xaxis.title = "Reference position"; L.yaxis.title = "% of reads mutated";
  L.xaxis.range = [r.params.analysis_window[0] - 2, r.params.analysis_window[1] + 2];
  Plotly.newPlot("f-map", traces, L, CONFIG);
}

/* ── tables ───────────────────────────────────────────────── */
const TABLES = {
  samples: ["Samples", r => r.samples],
  groups: ["Groups and statistics", r => r.groups],
  variants: ["Every mutation", r => r.variants],
  indels: ["Insertions and deletions", r => r.indels || []],
  positional: ["Mutations per position", r => r.positional],
  mutation_load: ["Mutations per read", r => Object.entries(r.mutation_load).flatMap(([s, h]) =>
    Object.entries(h).map(([k, n]) => ({ sample: s, mutations_per_read: +k, reads: n }))
      .sort((a, b) => a.mutations_per_read - b.mutations_per_read))],
  spectrum: ["Substitution spectrum", r => Object.entries(r.spectrum_group).flatMap(([g, c]) =>
    Object.entries(c).map(([k, v]) => ({ group: g, substitution: k, count: v })))],
  animal_tests: ["Group tests (animal as the unit)", r => r.animal_tests || []],
  rs_by_region: ["Replacement and silent, by region", r => r.rs_by_region || []],
  qc: ["QC flags (Sanger)", r => (r.qc || {}).flags || []],
  excluded: ["Clones not analysed (Sanger)", r => (r.qc || {}).excluded || []],
};
function buildTablePicker() {
  $("table-pick").innerHTML = Object.entries(TABLES)
    .map(([k, [label]]) => `<option value="${k}">${esc(label)}</option>`).join("");
  S.table = "samples"; S.sort = null;
  renderTable();
}
function currentRows() {
  let rows = TABLES[$("table-pick").value][1](S.result) || [];
  const q = $("table-filter").value.trim().toLowerCase();
  if (q) rows = rows.filter(r => Object.values(r).some(v => String(v).toLowerCase().includes(q)));
  if (S.sort) {
    rows = [...rows].sort((a, b) => {
      const x = a[S.sort], y = b[S.sort];
      if (x == null) return 1; if (y == null) return -1;
      return (typeof x === "number" && typeof y === "number" ? x - y : String(x).localeCompare(String(y))) * S.sortDir;
    });
  }
  return rows;
}
/* Column names for people. The CSV keeps the raw keys; each header's tooltip gives its key.
   [label, what it means] */
const COLS = {
  sample: ["Sample"], group: ["Group"], control: ["Control"], n: ["Animals"],
  compared_with: ["Compared with", "the control group these P values test against"],
  reads: ["Reads", "reads analysed, after filtering and downsampling"],
  unique_mutations: ["Unique mutations", "each distinct position and base counted once"],
  unique_density_x1e5: ["Unique mutations per 10⁵ bp", "unique mutations per 100,000 analysed bases"],
  one_event_density_x1e5: ["Detection floor per 10⁵ bp", "the density one mutation would give; where a sample with none is drawn"],
  mutation_frequency_x1e3: ["Mutations per 10³ bp", "every mutated base in every read, per 1,000 bases sequenced"],
  mutation_frequency_observed_x1e3: ["Mutations per 10³ observed bp", "as above, over the bases actually read at each position"],
  observed_bases: ["Observed bases"],
  mutated_reads_pct: ["Reads mutated (%)", "reads carrying at least one mutation"],
  mean_mutations_per_read: ["Mutations per read"],
  pct_GC: ["At G:C (%)", "unique mutations at a G or C"],
  pct_transitions: ["Transitions (%)", "A↔G and C↔T, of unique mutations"],
  pct_CT_GA: ["C→T and G→A (%)"],
  pct_GC_in_WRC_GYW: ["G:C in WRC/GYW (%)", "G:C mutations at AID's WRC/GYW hotspots"],
  gc_mutations: ["G:C mutations"],
  pct_AT_in_WA_TW: ["A:T in WA/TW (%)", "A:T mutations at Pol η's WA/TW hotspots"],
  pct_GC_in_WRCH_DGYW: ["G:C in WRCH/DGYW (%)", "G:C mutations at AID's extended WRCH/DGYW hotspots"],
  pct_GC_in_WRCY_RGYW: ["G:C in WRCY/RGYW (%)", "G:C mutations at AID's classic WRCY/RGYW hotspots"],
  pct_GC_in_WGCW: ["G:C in WGCW (%)", "G:C mutations at WGCW, an AID hotspot on both strands"],
  pct_GC_in_SYC_GRS: ["G:C in SYC/GRS (%)", "G:C mutations at AID's SYC/GRS coldspots"],
  pct_AT_in_WATW: ["A:T in WATW (%)", "A:T mutations at WATW, a Pol η hotspot on both strands"],
  at_mutations: ["A:T mutations"],
  analysed_bp: ["Window (bp)", "length of the analysed window"],
  mutations_total: ["Mutations, all reads", "every mutated base in every read"],
  frequency_bases: ["Bases (frequency)", "denominator of the mutation frequency"],
  density_bases: ["Bases (density)", "denominator of the unique-mutation density"],
  model_used: ["Error model"],
  clustered_calls: ["Clustered calls", "calls removed for sitting in a cluster of mismatches"],
  high_fraction_sites_removed: ["Germline-like sites removed", "positions mutated in most reads, set aside as germline"],
  high_fraction_sites: ["Germline-like sites"],
  warnings: ["Warnings"],
  reads_with_junction: ["Reads with a junction"], clones: ["Clones"],
  largest_clone_pct: ["Largest clone (%)", "share of the reads with a junction"],
  clone_weighted_frequency_x1e3: ["Clone-weighted mutations per 10³ bp", "each clone counted once, so one large clone cannot dominate"],
  multi_read_mutations: ["Mutations in 2+ reads"],
  single_read_observed: ["One-read mutations seen"],
  single_read_background: ["One-read background", "one-read events the controls predict for this sample"],
  single_read_estimate: ["One-read mutations, net", "seen minus background"],
  single_read_density_x1e5: ["One-read mutations per 10⁵ bp"],
  single_read_density_ci_low: ["One-read per 10⁵ bp, 95% CI low"],
  single_read_density_ci_high: ["One-read per 10⁵ bp, 95% CI high"],
  combined_density_x1e5: ["All mutations per 10⁵ bp", "called mutations plus the one-read estimate"],
  combined_density_ci_low: ["All per 10⁵ bp, 95% CI low"], combined_density_ci_high: ["All per 10⁵ bp, 95% CI high"],
  combined_density_ci_scope: ["CI covers"], null_from: ["Background from"], null_matching: ["Controls matched by"],
  pos: ["Position"], ref: ["Ref"], alt: ["Alt"], depth: ["Depth"], frequency: ["Frequency", "share of reads carrying it"],
  p_value: ["P"], plus_reads: ["Forward reads"], minus_reads: ["Reverse reads"], context: ["Context"],
  homopolymer: ["Homopolymer run (bp)", "longest run of one base touching the position"], WRC_GYW: ["At WRC/GYW"],
  motifs: ["Motifs", "the hot- and coldspot families the position sits in"],
  mutated_reads: ["Mutated reads"], frequency_pct: ["Reads mutated (%)"],
  mutations_per_read: ["Mutations per read"], substitution: ["Substitution"], count: ["Count"],
  region: ["Region"], replacement: ["Replacement"], silent: ["Silent"], stop: ["Stop"],
  non_coding: ["Non-coding"], rs_ratio: ["R/S", "replacement calls per silent call"],
  metric: ["Measure"], animals_case: ["Animals"], animals_control: ["Control animals"],
  rate_case: ["Rate"], rate_control: ["Control rate"], dispersion: ["Dispersion"], fold: ["Fold"],
  fold_control_zero_one_event: ["Fold, zero control at one event"], p_welch_log: ["P, Welch on log"],
  zero_count_convention: ["Zero-count rule"], p_welch_log_half_event: ["P, log, zeros at half an event"],
  p_welch_log_third_of_smallest: ["P, log, zeros at a third of the smallest"],
  fold_geometric: ["Fold (geometric)"], fold_ci: ["Fold, 95% CI"], p_permutation: ["P, permutation"],
  lrt_negbinomial: ["LRT, negative binomial"], p_negbinomial: ["P, negative binomial"],
  p_pooled_reads: ["P, reads pooled", "treats reads, not animals, as replicates: shown for comparison only"],
};
// the groups table's statistics are "<measure><suffix>"; the control's name is sanitised in the key
const STAT_SUFFIXES = [
  [/^_mean$/, "mean"], [/^_sd$/, "SD"],
  [/^_p_welch_vs_/, "P, Welch"], [/^_p_welch_log_vs_/, "P, Welch on log"],
  [/^_sensitivity_p_welch_log_half_event_vs_/, "P, log, zeros at half an event"],
  [/^_sensitivity_p_welch_log_third_of_smallest_vs_/, "P, log, zeros at a third of the smallest"],
  [/^_p_mwu_vs_/, "P, Mann–Whitney"],
];
function colInfo(c) {
  if (COLS[c]) return { label: COLS[c][0], note: COLS[c][1] || "" };
  const base = Object.keys(COLS).filter(k => c.startsWith(k + "_")).sort((a, b) => b.length - a.length)[0];
  const stat = base && STAT_SUFFIXES.find(([re]) => re.test(c.slice(base.length)));
  if (stat) return { label: `${COLS[base][0]}: ${stat[1]}`, note: COLS[base][1] || "", stat: stat[1], base };
  return { label: c.replace(/_/g, " "), note: "" };
}
// the columns most people want; the rest are one tick away
const KEY_COLS = {
  samples: ["sample", "group", "control", "reads", "unique_mutations", "unique_density_x1e5",
            "mutation_frequency_x1e3", "mutated_reads_pct", "mean_mutations_per_read", "pct_GC_in_WRC_GYW",
            "pct_transitions", "clones", "largest_clone_pct", "warnings"],
  groups: (c) => ["group", "n", "control", "compared_with"].includes(c) || (() => {
    const i = colInfo(c);
    return ["unique_density_x1e5", "mutation_frequency_x1e3", "mutated_reads_pct"].includes(i.base)
      && ["mean", "SD", "P, Welch on log", "P, Mann–Whitney"].includes(i.stat);
  })(),
  animal_tests: ["metric", "group", "control", "animals_case", "animals_control", "rate_case", "rate_control",
                 "fold", "fold_ci", "p_welch_log", "p_permutation", "p_negbinomial"],
};
const isP = (c) => /(^|_)p_|^p$|^p_value$/.test(c);
// one number format per column, chosen from the column's own values, so its decimals line up
function colFormat(c, rows) {
  const vals = rows.map(r => r[c]).filter(v => typeof v === "number" && Number.isFinite(v));
  if (!vals.length) return null;
  if (isP(c)) return (v) => v === 0 ? "&lt;1e-308"    // underflowed: smaller than a double, not zero
    : v < 1e-4 ? v.toExponential(1) : v.toFixed(4);
  if (vals.every(Number.isInteger)) return fmtInt;
  let d;
  if (/(^|_)pct(_|$)|\(%\)/.test(c) || /\(%\)/.test(colInfo(c).label)) d = 1;
  else {
    const nz = vals.filter(v => v !== 0).map(Math.abs).sort((a, b) => a - b);
    const med = nz.length ? nz[Math.floor(nz.length / 2)] : 1;
    d = Math.max(0, Math.min(4, 2 - Math.floor(Math.log10(med))));
  }
  // a value the column's decimals would round to zero keeps its size
  return (v) => v !== 0 && Math.abs(v) < 0.5 * 10 ** -d ? v.toExponential(1)
    : v.toLocaleString("en-GB", { minimumFractionDigits: d, maximumFractionDigits: d });
}
function fmtCell(v, f) {
  if (v === true) return "<span class='pill'>yes</span>";
  if (v === false || v == null) return "";
  if (typeof v === "number") return f ? f(v) : String(v);
  // an interval: two numbers
  if (Array.isArray(v) && v.length === 2 && v.every(x => typeof x === "number"))
    return v.map(x => Math.abs(x) >= 1000 ? fmtInt(Math.round(x)) : String(+x.toPrecision(3))).join("–");
  if (Array.isArray(v)) return esc(v.map(x => x && typeof x === "object"
    ? (x.pos != null ? `${x.pos}${x.alt ? " " + x.alt : ""}${x.frequency != null ? ` (${Math.round(100 * x.frequency)}%)` : ""}` : JSON.stringify(x))
    : x).join(", "));
  if (typeof v === "object") return esc(JSON.stringify(v));
  return esc(v);
}

function renderTable() {
  const rows = currentRows();
  const thead = $("table").tHead, tbody = $("table").tBodies[0];
  if (!rows.length) {
    thead.innerHTML = ""; tbody.innerHTML = "<tr><td class=muted>No rows.</td></tr>";
    $("table-all-box").hidden = true; $("table-count").textContent = ""; return;
  }
  // every column any row has: a control without peers lacks the estimates
  const every = [...new Set(rows.flatMap(r => Object.keys(r)))];
  const key = KEY_COLS[$("table-pick").value];
  const keep = !key ? () => true : typeof key === "function" ? key : (c) => key.includes(c);
  const fewer = every.filter(keep);
  const trimmed = key && fewer.length < every.length;
  $("table-all-box").hidden = !trimmed;
  $("table-all-n").textContent = trimmed ? `(${every.length})` : "";
  const cols = trimmed && !$("table-all").checked ? fewer : every;
  // the control's name in a groups key is sanitised; show the real one
  const ctl = rows.find(r => r.compared_with)?.compared_with;
  const fmts = Object.fromEntries(cols.map(c => [c, colFormat(c, rows)]));
  const num = (c) => rows.some(r => typeof r[c] === "number");
  thead.innerHTML = "<tr>" + cols.map(c => {
    const i = colInfo(c);
    const label = i.stat && i.stat.startsWith("P") && ctl ? `${i.label} vs ${ctl}` : i.label;
    const tip = `${i.note ? i.note + "\n" : ""}CSV column: ${c}`;
    return `<th onclick="sortBy('${c}')" title="${esc(tip)}" class="${num(c) ? "num" : ""}">${esc(label)}${S.sort === c ? (S.sortDir > 0 ? " ▲" : " ▼") : ""}</th>`;
  }).join("") + "</tr>";
  tbody.innerHTML = rows.slice(0, 3000).map(r => "<tr>" + cols.map(c =>
    `<td class="${typeof r[c] === "number" ? "num" : ""}">${fmtCell(r[c], fmts[c])}</td>`).join("") + "</tr>").join("");
  $("table-count").textContent = `${rows.length} row${rows.length === 1 ? "" : "s"}${rows.length > 3000 ? " (showing 3000)" : ""}`;
}
function sortBy(c) { S.sortDir = S.sort === c ? -S.sortDir : 1; S.sort = c; renderTable(); }
function downloadCurrentTable() {
  const map = { samples: "samples.csv", groups: "groups.csv", variants: "variants.csv",
    indels: "indels.csv",
    positional: "positional_frequency.csv", spectrum: "spectrum.csv",
    mutation_load: "mutation_load.csv", animal_tests: "animal_tests.csv", rs_by_region: "rs_by_region.csv",
    qc: "qc_flags.csv", excluded: "excluded.csv" };
  dl(map[$("table-pick").value]);
}

/* ── export ───────────────────────────────────────────────── */
const FIGURES = [
  ["density_by_group", "Mutation density"],
  ["mutation_frequency_by_group", "Mutation frequency"],
  ["mutated_reads_by_group", "% of reads mutated"],
  ["substitution_matrix", "From → to substitutions"],
  ["spectrum_by_group", "Substitution spectrum"],
  ["positional_map", "Position along the amplicon"],
  ["mutation_load", "Mutations per read"],
  ["hotspot_targeting", "Hotspot targeting"],
];
// Detail and diagnostics. Each is drawn only when the run carries what it
// needs, so the Export tab asks the server which ones exist.
const DIAGNOSTIC_FIGURES = [
  ["coverage", "Coverage along the amplicon"],
  ["clone_map", "Clones from the V(D)J junction"],
  ["targeting_profile", "Do the mutations follow the AID targeting model?"],
  ["targeting_spectrum", "Spectrum against the targeting model"],
  ["mutability_hedgehog", "The targeting model (hedgehog plot)"],
  ["spectrum_corrected", "Spectrum, corrected for base composition"],
  ["saturation", "Saturation (is the depth enough?)"],
  ["combined_density", "Combined density with 95% intervals"],
  ["single_read_excess", "Single-read excess over controls"],
  ["indels", "Insertions and deletions"],
  ["positional_by_sample", "Position along the amplicon, per sample"],
  ["quality_histogram", "Base quality: mismatches vs matches"],
  ["phasing_blocks", "Phasing: how reads carry a block"],
  ["qc_waterfall", "Where reads are lost"],
  ["control_background", "The background your controls measured"],
  ["clonality", "Clonal structure"],
  ["volcano", "What the caller kept"],
  ["replicate_concordance", "Do replicates agree?"],
];
function buildExport() {
  $("export-tables").innerHTML = ["samples.csv", "groups.csv", "variants.csv",
    ...((S.result.indels || []).length ? ["indels.csv"] : []),
    ...(S.result.coverage && Object.keys(S.result.coverage).length ? ["coverage.csv"] : []),
    ...(S.result.clones && Object.keys(S.result.clones).length ? ["clones.csv"] : []),
    "positional_frequency.csv", "spectrum.csv", "mutation_load.csv"]
    .map(f => `<button class="btn small" onclick="dl('${f}')">${f}</button>`).join("");
  const avail = new Set(S.result.available_figures || FIGURES.map(([f]) => f));
  const shown = [...FIGURES, ...DIAGNOSTIC_FIGURES].filter(([f]) => avail.has(f));
  $("export-figs").innerHTML = shown.map(([f, label]) => `
    <div class="row" style="margin:4px 0">
      <span style="flex:2 1 220px">${esc(label)}</span>
      ${["png", "svg", "pdf"].map(x => `<button class="btn small" onclick="dl('${f}.${x}')">${x.toUpperCase()}</button>`).join("")}
    </div>`).join("");
}
function dl(what) {
  const a = document.createElement("a");
  a.href = `/api/export/${S.job}/${what}`; a.download = what;
  document.body.appendChild(a); a.click(); a.remove();
}
async function saveLocal() {
  const r = await fetch("/api/save-to/" + S.job, {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ path: $("save-path").value }),
  });
  const j = await r.json();
  $("save-msg").textContent = r.ok ? `Saved ${j.files.length} files to ${j.folder}` : (j.detail || "could not save");
}

matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
  if (S.result) { $("figs").dataset.job = ""; drawFigures(); }
});
const READY = init();


/* ── reads: the alignment, before and after processing ────── */
// Colours run from "analysed" through the stages that drop a read, in the
// order the pipeline applies them, so the legend reads like the pipeline.
const FATE_COLOUR = {
  "kept": "#3F6FA8",
  "kept (short)": "#6FA0D8",
  "not sampled": "#9AA6B4",
  "minor haplotype": "#8E6BB5",
  "alignment too short": "#C98A2E",
  "chimeric": "#B5651D",
  "contaminant": "#C0392B",
  "off target": "#A33E6B",
  "low quality": "#7A8794",
  "too long": "#6B7A5A",
  "no quality segment": "#7A8794",
  "unmapped": "#5D6D7E",
};
// an analysed read: full length, or kept below the minimum span (rescue_short_reads)
const isKept = (fate) => fate === "kept" || fate === "kept (short)";
const RD = { sample: null, from: 1, to: null, offset: 0, limit: 300, data: null,
             index: null, rowH: 7, left: 108 };

async function rdInit() {
  if (RD.index) return;
  const r = await fetch(`/api/alignment/${S.job}`);
  if (!r.ok) { $("rd-status").textContent = "This run did not record read alignments."; return; }
  RD.index = await r.json();
  $("rd-sample").innerHTML = RD.index.samples
    .map(s => `<option value="${esc(s.sample)}">${esc(s.sample)} — ${s.reads.toLocaleString()} reads</option>`).join("");
  const first = RD.index.samples[0];
  if (!first) { $("rd-status").textContent = "No alignments were recorded."; return; }
  RD.sample = first.sample; RD.from = 1; RD.to = first.reference_length;
  $("rd-from").value = RD.from; $("rd-to").value = RD.to;
  rdLoad(true);
}

function rdCurrent() {
  return (RD.index.samples || []).find(s => s.sample === $("rd-sample").value);
}

async function rdLoad(reset) {
  if (!RD.index) return;
  if (reset) RD.offset = 0;
  RD.sample = $("rd-sample").value;
  RD.from = Math.max(1, +$("rd-from").value || 1);
  RD.to = Math.max(RD.from + 1, +$("rd-to").value || (rdCurrent() || {}).reference_length || 1000);
  const show = $("rd-show").value;
  let fates = "";
  if (show === "kept") fates = "kept";
  else if (show === "dropped") fates = Object.keys(FATE_COLOUR).filter(f => !isKept(f)).join(",");
  const q = new URLSearchParams({ start: RD.from, end: RD.to, offset: RD.offset,
                                  limit: RD.limit, ...(fates ? { fates } : {}) });
  const r = await fetch(`/api/alignment/${S.job}/${encodeURIComponent(RD.sample)}?${q}`);
  if (!r.ok) { $("rd-status").textContent = (await r.json()).detail || "could not load"; return; }
  RD.data = await r.json();
  rdDraw();
}

function rdZoom(by) {
  const mid = (RD.from + RD.to) / 2, span = (RD.to - RD.from) * by;
  const len = (rdCurrent() || {}).reference_length || RD.to;
  $("rd-from").value = Math.max(1, Math.round(mid - span / 2));
  $("rd-to").value = Math.min(len, Math.round(mid + span / 2));
  rdLoad(true);
}
function rdAll() {
  $("rd-from").value = 1;
  $("rd-to").value = (rdCurrent() || {}).reference_length || 1000;
  rdLoad(true);
}
function rdPage(dir) {
  const next = RD.offset + dir * RD.limit;
  if (next < 0 || next >= (RD.data ? RD.data.total : 0)) return;
  RD.offset = next; rdLoad(false);
}

function rdDraw() {
  const d = RD.data, cv = $("rd-canvas"), ctx = cv.getContext("2d");
  const reads = d.reads, span = d.end - d.start + 1;
  const width = Math.max(760, $("rd-scroll").clientWidth - 2);
  const plot = width - RD.left - 14;
  const perBase = plot / span;
  const head = 44;
  const dpr = window.devicePixelRatio || 1;
  const height = head + reads.length * RD.rowH + 10;
  cv.width = width * dpr; cv.height = height * dpr;
  cv.style.width = width + "px"; cv.style.height = height + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const css = getComputedStyle(document.body);
  const fg = css.getPropertyValue("--fg") || "#111";
  const bg = css.getPropertyValue("--bg") || "#fff";
  ctx.fillStyle = bg; ctx.fillRect(0, 0, width, height);
  const x = p => RD.left + (p - d.start) * perBase;

  // called mutations, as vertical guides behind everything
  if (perBase >= 1.2) {
    for (const c of d.called) {
      ctx.fillStyle = "rgba(192,57,43,.13)";
      ctx.fillRect(x(c.pos), head - 6, Math.max(1, perBase), height - head);
    }
  }
  // the reference: its sequence when there is room, otherwise a rule
  ctx.font = "11px ui-monospace, Menlo, monospace";
  ctx.fillStyle = fg;
  ctx.globalAlpha = 0.75;
  if (perBase >= 7) {
    for (let p = d.start; p <= d.end; p++) {
      ctx.fillText(d.sequence[p - d.start] || "", x(p) + perBase / 2 - 3, 26);
    }
  } else {
    ctx.fillRect(RD.left, 22, plot, 2);
  }
  ctx.globalAlpha = 1;
  ctx.font = "10px system-ui, sans-serif";
  for (let i = 0; i <= 6; i++) {                       // position ruler
    const p = Math.round(d.start + (span - 1) * i / 6);
    ctx.fillStyle = fg; ctx.globalAlpha = 0.55;
    ctx.fillText(String(p), Math.min(x(p), width - 30), 12);
    ctx.globalAlpha = 1;
  }

  let band = null;
  reads.forEach((rd, i) => {
    const y = head + i * RD.rowH;
    const colour = FATE_COLOUR[rd.fate] || "#888";
    if (rd.fate !== band) {                              // name each band once
      band = rd.fate;
      ctx.fillStyle = fg; ctx.globalAlpha = 0.1;
      ctx.fillRect(0, y, width, 1);
      ctx.globalAlpha = 0.85; ctx.font = "600 9px system-ui, sans-serif";
      ctx.fillStyle = colour;
      ctx.fillText(rd.fate, 4, y + 9);
      ctx.globalAlpha = 1;
    }
    const a = Math.max(RD.left, x(Math.max(rd.start, d.start)));
    const b = Math.min(RD.left + plot, x(Math.min(rd.end || d.end, d.end)));
    ctx.fillStyle = colour;
    ctx.globalAlpha = isKept(rd.fate) ? 0.3 : 0.22;
    if (b > a) ctx.fillRect(a, y + 1, b - a, RD.rowH - 2);   // the read's span
    ctx.globalAlpha = 1;
    for (const [pos, base] of rd.diffs) {                     // its mismatches
      if (pos < d.start || pos > d.end) continue;
      ctx.fillStyle = colour;
      ctx.fillRect(x(pos), y + 1, Math.max(1.4, perBase), RD.rowH - 2);
      if (perBase >= 7) {
        ctx.fillStyle = bg; ctx.font = "9px ui-monospace, monospace";
        ctx.fillText(base, x(pos) + perBase / 2 - 2.5, y + RD.rowH - 1.5);
      }
    }
    for (const [pos, op] of rd.indels) {                      // and its indels
      if (pos < d.start || pos > d.end) continue;
      ctx.fillStyle = op === "D" ? "#111" : "#1B8A5A";
      ctx.fillRect(x(pos) - 1, y, 2.2, RD.rowH);
    }
  });
  // row labels only when the rows are tall enough to read
  if (RD.rowH >= 14) {
    ctx.font = "9px ui-monospace, monospace"; ctx.globalAlpha = 0.7;
    reads.forEach((rd, i) => {
      ctx.fillStyle = fg;
      ctx.fillText(rd.name.slice(0, 14), 4, head + i * RD.rowH + RD.rowH - 1.5);
    });
    ctx.globalAlpha = 1;
  }
  cv.onclick = ev => {
    const r = cv.getBoundingClientRect();
    const i = Math.floor((ev.clientY - r.top - head) / RD.rowH);
    if (i >= 0 && i < reads.length) rdDetail(reads[i]);
  };

  const c = d.counts, kept = (c.kept || 0) + (c["kept (short)"] || 0);
  const total = Object.values(c).reduce((a, b) => a + b, 0);
  $("rd-status").textContent =
    `${RD.sample}: ${total.toLocaleString()} reads recorded, ${kept.toLocaleString()} analysed`
    + (d.truncated ? ` (${d.truncated.toLocaleString()} more not recorded)` : "")
    + ` · showing ${reads.length} of ${d.total.toLocaleString()} over ${d.start}–${d.end}`
    + ` · ${d.called.length} called mutation${d.called.length === 1 ? "" : "s"} in view`;
  $("rd-page").textContent = `reads ${d.offset + 1}–${d.offset + reads.length} of ${d.total.toLocaleString()}`;
  $("rd-legend").innerHTML = Object.entries(c)
    .sort((a, b) => b[1] - a[1])
    .map(([f, n]) => `<span class="hint"><span style="display:inline-block;width:10px;height:10px;
      border-radius:2px;background:${FATE_COLOUR[f] || "#888"};margin-right:5px"></span>${esc(f)} (${n.toLocaleString()})</span>`)
    .join("") + `<span class="hint"><span style="display:inline-block;width:3px;height:10px;
      background:#111;margin-right:5px"></span>deletion</span>
      <span class="hint"><span style="display:inline-block;width:3px;height:10px;
      background:#1B8A5A;margin-right:5px"></span>insertion</span>`;
}

function rdDetail(rd) {
  const box = $("rd-detail");
  box.style.display = "";
  const diffs = rd.diffs.map(([p, b, q]) => `${p}${b} (Q${q})`).join(", ") || "none in this window";
  const ind = rd.indels.map(([p, op, x]) => `${p} ${op === "D" ? "del " + x : "ins " + x}`).join(", ") || "none";
  box.innerHTML = `<h3 style="margin:0 0 6px">${esc(rd.name)}</h3>
    <p class="hint" style="margin:0 0 8px">${esc(rd.fate)} ·
      ${rd.start}–${rd.end} · ${rd.strand === 1 ? "forward" : "reverse"} strand</p>
    ${rd.reason ? `<p class="hint" style="margin:0 0 8px">${esc(rd.reason)}</p>` : ""}
    <p style="margin:0 0 4px"><strong>Mismatches:</strong> ${esc(diffs)}</p>
    <p style="margin:0"><strong>Indels:</strong> ${esc(ind)}</p>`;
}

/* ── class-switch junctions (switch.py) ───────────────────── */
const SW = { reference: null, samples: [], job: null, result: null };
const MH_CLASSES = [["blunt", "#9DB4C0"], ["1-3 nt", "#5B9BD5"], [">=4 nt", "#FF6B57"], ["insertion", "#B07A1E"]];

function swInit() {
  dropZone($("sw-drop-ref"), false, (f) => upload(f).then(async ([x]) => {
    if (!x) return;
    SW.reference = x;
    $("sw-ref-name").textContent = "Reference: " + x.name;
    $("sw-error").innerHTML = "";
    const r = await fetch("/api/fasta-names", { method: "POST",
      headers: { "Content-Type": "application/json" }, body: JSON.stringify({ path: x.path }) });
    const d = await r.json();
    if (!r.ok) return swError(d.detail || "could not read the FASTA");
    $("sw-donor").innerHTML = d.records.map((rec, i) =>
      `<option value="${esc(rec.name)}"${i === 0 ? " selected" : ""}>${esc(rec.name)} (${rec.length} bp)</option>`).join("");
    if (d.records.length < 2) swError("The reference has one record; it needs the donor (Sμ) and at least one acceptor.");
  }));
  dropZone($("sw-drop-samples"), true, (f) => upload(f).then(swAdd), "sw-pick-folder");
}

function swAdd(files) {
  const folders = new Map();
  for (const f of files) {
    const rel = f.rel || "";
    const top = rel.includes("/") ? rel.split("/")[0] : null;
    if (top) {
      if (!folders.has(top)) folders.set(top, (f.folder || f.path).replace(/\/[^/]*$/, ""));
    } else SW.samples.push({ name: strip(f.name), path: f.path, group: "" });
  }
  for (const [name, dir] of folders) SW.samples.push({ name, path: dir, group: "" });
  swRender();
}

function swFromSetup() {
  SW.samples = S.samples.map(s => ({ name: s.name, path: s.path, group: s.group || "" }));
  swRender();
}

function swRender() {
  $("sw-samples").innerHTML = SW.samples.map((s, i) => `
    <div class="sample">
      <span class="nm" title="${esc(s.path)}">${esc(s.name)}</span>
      <input placeholder="group" value="${esc(s.group)}" oninput="SW.samples[${i}].group=this.value;swGroups()">
      <button class="x" title="Remove" onclick="SW.samples.splice(${i},1);swRender()">✕</button>
    </div>`).join("");
  swGroups();
}

function swGroups() {
  const gs = [...new Set(SW.samples.map(s => s.group || s.name))];
  const opts = ["<option value=''>(no comparison)</option>"].concat(gs.map(g => `<option>${esc(g)}</option>`)).join("");
  const keep = [$("sw-case").value, $("sw-control").value];
  $("sw-case").innerHTML = opts; $("sw-control").innerHTML = opts;
  if (gs.includes(keep[0])) $("sw-case").value = keep[0];
  if (gs.includes(keep[1])) $("sw-control").value = keep[1];
}

function swError(msg) { $("sw-error").innerHTML = `<div class="msg err" style="margin-top:8px">${esc(msg)}</div>`; }

async function swRun() {
  $("sw-error").innerHTML = "";
  if (!SW.reference) return swError("Choose the switch-region reference first.");
  if (!SW.samples.length) return swError("Add at least one sample.");
  const c = $("sw-case").value, k = $("sw-control").value;
  if ((c || k) && (!c || !k || c === k)) return swError("Pick two different groups to compare, or neither.");
  $("sw-run").disabled = true; $("sw-status").style.display = "";
  const body = { reference: SW.reference.path, donor: $("sw-donor").value,
    samples: SW.samples.map(s => ({ name: s.name, path: s.path, group: s.group || s.name })),
    compare: { case: c, control: k } };
  const r = await fetch("/api/switch/run", { method: "POST",
    headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  const d = await r.json();
  if (!r.ok) { $("sw-run").disabled = false; return swError(d.detail || "could not start"); }
  SW.job = d.job;
  swPoll();
}

async function swPoll() {
  const j = await (await fetch("/api/job/" + SW.job)).json();
  $("sw-bar").style.width = Math.round(100 * j.progress) + "%";
  $("sw-msg").textContent = j.message;
  if (j.state === "error") { $("sw-run").disabled = false; return swError(j.error); }
  if (j.state !== "done") return setTimeout(swPoll, 700);
  SW.result = await (await fetch("/api/switch/result/" + SW.job)).json();
  $("sw-run").disabled = false;
  swShow();
}

function swTable(el, rows) {
  if (!rows.length) { el.innerHTML = "<tr><td class=muted>No rows.</td></tr>"; return; }
  const cols = [...new Set(rows.flatMap(r => Object.keys(r)))];
  el.innerHTML = "<thead><tr>" + cols.map(c => `<th>${esc(c.replace(/_/g, " "))}</th>`).join("") +
    "</tr></thead><tbody>" + rows.map(r => "<tr>" + cols.map(c => {
      const v = r[c];
      return `<td class="${typeof v === "number" ? "num" : ""}">${v == null ? "" : typeof v === "number" ? fmtNum(c, v) : esc(v)}</td>`;
    }).join("") + "</tr>").join("") + "</tbody>";
}

function swShow() {
  const R = SW.result;
  $("sw-results").hidden = false;
  const names = Object.keys(R.samples);
  swTable($("sw-summary"), names.map(n => {
    const s = R.samples[n], c = s.microhomology_classes_pct || {};
    return { sample: n, group: R.groups[n], reads: s.reads, junction_reads: s.junction_reads,
      unique_junctions: s.unique_junctions,
      acceptors: Object.entries(s.acceptors || {}).map(([k, v]) => `${k} ${v}`).join(", "),
      "blunt %": c["blunt"] || 0, "1-3 nt %": c["1-3 nt"] || 0, ">=4 nt %": c[">=4 nt"] || 0,
      "insertion %": c["insertion"] || 0, mean_microhomology: s.mean_microhomology,
      "ambiguous %": s.ambiguous_pct, mutations_near_join_per_1e3: s.mutations_near_join_per_1e3,
      internal_deletions: s.internal_deletions };
  }));
  $("sw-comp-card").hidden = !R.comparison.length;
  swTable($("sw-comp"), R.comparison);
  const L = { ...withLegend(layoutBase()), barmode: "stack" };
  L.yaxis.title = "% of unique junctions";
  Plotly.newPlot("sw-fig-mh", MH_CLASSES.map(([k, col]) => ({
    type: "bar", name: k, x: names, y: names.map(n => (R.samples[n].microhomology_classes_pct || {})[k] || 0),
    marker: { color: col } })), L, CONFIG);
  const sw = R.junctions.filter(j => j.kind === "switch");
  const L2 = { ...withLegend(layoutBase()), barmode: "overlay" };
  L2.xaxis.title = "Breakpoint in the donor (bp)"; L2.yaxis.title = "Junction reads";
  Plotly.newPlot("sw-fig-bp", names.map((n, i) => ({
    type: "histogram", name: n, x: sw.filter(j => j.sample === n).map(j => j.donor_break),
    opacity: 0.6, nbinsx: 40, marker: { color: PALETTE[i % PALETTE.length] } })), L2, CONFIG);
}

function swDl(what) { location.href = `/api/switch/export/${SW.job}/${what}`; }


/* ── repair routes (repair.py) ────────────────────────────── */
const RP_INDICES = [["pct_AT", "A:T of all"], ["pct_GC_Ts_of_GC", "G:C Ts of G:C"],
  ["pct_GC_Tv_of_GC", "G:C Tv of G:C"], ["pct_CG_of_GC_Tv", "C>G + G>C of Tv"]];
const RP_MIN_MUTATIONS = 20;     // a group with fewer gives proportions too noisy to chart

function resultGroups() { return S.result ? S.result.group_order || [...new Set(S.result.samples.map(s => s.group))] : []; }

function rpInit() {
  const gs = resultGroups();
  const opts = gs.map(g => `<option>${esc(g)}</option>`).join("");
  const keep = [$("rp-case").value, $("rp-control").value];
  $("rp-case").innerHTML = opts;
  $("rp-control").innerHTML = "<option value=''>(none)</option>" + opts;
  const ctrl = S.result.samples.find(s => s.control);
  $("rp-case").value = gs.includes(keep[0]) ? keep[0] : (gs.find(g => !ctrl || g !== ctrl.group) || gs[0]);
  $("rp-control").value = gs.includes(keep[1]) ? keep[1] : "";
  $("rp-refs").innerHTML = gs.map(g => `<label style="margin:0;display:flex;gap:6px;align-items:center;color:var(--fg)">
      <input type="checkbox" style="width:auto" value="${esc(g)}"> ${esc(g)}</label>`).join("");
}

async function rpRun() {
  $("rp-error").innerHTML = "";
  const refs = {};
  for (const c of $("rp-refs").querySelectorAll("input:checked")) refs[c.value] = c.value;
  if (refs[$("rp-case").value]) return $("rp-error").innerHTML = `<div class="msg err" style="margin-top:8px">The group cannot be its own reference.</div>`;
  if (Object.keys(refs).length === 1) return $("rp-error").innerHTML = `<div class="msg err" style="margin-top:8px">A decomposition needs two or more reference groups.</div>`;
  const r = await fetch("/api/repair/" + S.job, { method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ case: $("rp-case").value, control: $("rp-control").value, references: refs }) });
  const d = await r.json();
  if (!r.ok) return $("rp-error").innerHTML = `<div class="msg err" style="margin-top:8px">${esc(d.detail || "failed")}</div>`;
  $("rp-results").hidden = false;
  swTable($("rp-animals"), d.per_animal);
  $("rp-comp-card").hidden = !d.comparison.length;
  swTable($("rp-comp"), d.comparison);
  const total = (g) => d.per_animal.filter(a => a.group === g).reduce((t, a) => t + a.mutations, 0);
  const all = resultGroups(), gs = all.filter(g => total(g) >= RP_MIN_MUTATIONS);
  const left = all.filter(g => !gs.includes(g));
  const L = withLegend(layoutBase()); L.barmode = "group"; L.yaxis.title = "% (mean of animals)";
  if (left.length) L.annotations = [{ text: `Not shown, too few mutations: ${left.map(g => `${g} (${total(g)})`).join(", ")}`,
    xref: "paper", yref: "paper", x: 0, y: -0.22, showarrow: false, xanchor: "left", font: { size: 11 } }];
  Plotly.newPlot("rp-fig", gs.map((g, i) => {
    const rows = d.per_animal.filter(a => a.group === g);
    return { type: "bar", name: g, x: RP_INDICES.map(x => x[1]),
      y: RP_INDICES.map(([k]) => { const v = rows.map(a => a[k]).filter(v => v != null); return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null; }),
      marker: { color: colourOf(g) } };
  }), L, CONFIG);
  $("rp-dec-card").hidden = !d.decomposition;
  if (d.decomposition) {
    const D = d.decomposition, names = Object.keys(D.weights);
    $("rp-dec").innerHTML = `<p class="hint">${esc(D.group)} as a mixture of ${names.map(esc).join(", ")}; fit R² ${D.r2.toFixed(2)}.</p>`
      + (D.warning ? `<div class="msg err" style="margin:6px 0">${esc(D.warning)}</div>` : "");
    const L2 = layoutBase(); L2.yaxis.title = "share of the spectrum"; L2.yaxis.range = [0, 1];
    Plotly.newPlot("rp-dec-fig", [{ type: "bar", x: names, y: names.map(n => D.weights[n]),
      error_y: { type: "data", symmetric: false, array: names.map(n => D.ci[n][1] - D.weights[n]),
                 arrayminus: names.map(n => D.weights[n] - D.ci[n][0]) },
      marker: { color: names.map(n => colourOf(n)) } }], L2, CONFIG);
  }
}

/* ── power (power.py) ─────────────────────────────────────── */
function pwSource() {
  const manual = $("pw-source").value === "manual";
  $("pw-manual").hidden = !manual; $("pw-from-result").hidden = manual;
}

function pwInit() {
  const gs = resultGroups().filter(g => !(S.result.samples.find(s => s.group === g) || {}).control);
  $("pw-group").innerHTML = gs.map(g => `<option>${esc(g)}</option>`).join("");
  $("pw-noresult").hidden = !!S.result;
  if (!S.result) { $("pw-source").value = "manual"; }
  pwSource();
}

async function pwRun() {
  $("pw-error").innerHTML = "";
  const nums = (id) => $(id).value.split(/[,\s]+/).filter(Boolean).map(Number).filter(n => n > 0);
  const body = { fold: +$("pw-fold").value, animals: nums("pw-animals"), reads: nums("pw-reads"),
    test: $("pw-test").value, sims: +$("pw-sims").value };
  if ($("pw-source").value === "result") {
    if (!S.result) return $("pw-error").innerHTML = `<div class="msg err">Run an analysis first, or choose numbers.</div>`;
    Object.assign(body, { job: S.job, group: $("pw-group").value, metric: $("pw-metric").value });
  } else Object.assign(body, { rate: +$("pw-rate").value, cv: +$("pw-cv").value,
    dispersion: +$("pw-disp").value, bp: +$("pw-bp").value });
  $("pw-run").disabled = true; $("pw-msg").textContent = "simulating…";
  const r = await fetch("/api/power", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  const d = await r.json();
  $("pw-run").disabled = false; $("pw-msg").textContent = "";
  if (!r.ok) return $("pw-error").innerHTML = `<div class="msg err">${esc(d.detail || "failed")}</div>`;
  const U = d.used;
  $("pw-used").textContent = `${U.fold}-fold difference; rate ${U.rate.toExponential(2)} per base, animal-to-animal CV ${U.cv_between.toFixed(2)}, clumping ${U.dispersion.toFixed(1)}`
    + (d.pilot ? ` (from ${d.pilot.group}, ${d.pilot.animals} animals${d.pilot.animals < 5 ? ": an estimate from so few animals is itself uncertain" : ""})` : "") + `; ${U.sims} experiments per cell.`;
  const animals = [...new Set(d.rows.map(x => x.animals_per_group))], reads = [...new Set(d.rows.map(x => x.reads_per_animal))];
  const cell = (p) => `<td class="num" style="background:${p >= 0.8 ? "#2d7a4b33" : p >= 0.5 ? "#B07A1E26" : "transparent"}">${p.toFixed(2)}</td>`;
  $("pw-table").innerHTML = "<thead><tr><th>animals per group</th>" + reads.map(r => `<th class="num">${r} reads</th>`).join("") + "</tr></thead><tbody>"
    + animals.map(a => `<tr><td>${a}</td>` + reads.map(r => cell(d.rows.find(x => x.animals_per_group === a && x.reads_per_animal === r).power)).join("") + "</tr>").join("") + "</tbody>";
  $("pw-results").hidden = false;
}


/* ── screen (demux.screen_table) ──────────────────────────── */
let SC = null;
async function scRun() {
  $("sc-error").innerHTML = "";
  const r = await fetch("/api/screen/" + S.job, { method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ metric: $("sc-metric").value }) });
  const d = await r.json();
  if (!r.ok) return $("sc-error").innerHTML = `<div class="msg err" style="margin-top:8px">${esc(d.detail || "failed")}</div>`;
  SC = d;
  const qcut = +$("sc-q").value || 0.1;
  $("sc-results").hidden = false;
  swTable($("sc-table"), d.rows.map(x => ({ condition: x.condition, control: x.control, animals: x.animals,
    fold: x.fold == null ? "" : fmtNum("", x.fold) + (x.fold_control_zero_one_event ? " †" : ""),
    "fold 95% CI": x.fold_ci ? `${fmtNum("", x.fold_ci[0])}–${fmtNum("", x.fold_ci[1])}` : "",
    p: x.p, q: x.q, [`q < ${qcut}`]: x.q != null && x.q < qcut ? "yes" : "" })));
  const pts = d.rows.filter(x => x.fold);
  const lb = d.rows.some(x => x.fold_control_zero_one_event);
  $("sc-error").innerHTML = lb ? `<div class="msg" style="margin-top:8px">The control group has no mutations,
    so folds marked † are against one event per control animal (as the plots show a zero). This is a zero-count convention, not a confidence bound.</div>` : "";
  const hit = (x) => x.q != null && x.q < qcut;
  const L = layoutBase(); L.xaxis.title = "log2 fold vs control"; L.yaxis.title = "−log10 q";
  L.shapes = [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: -Math.log10(qcut), y1: -Math.log10(qcut),
    line: { dash: "dash", color: "#999", width: 1 } }];
  Plotly.newPlot("sc-fig", [{ type: "scatter", mode: "markers+text",
    x: pts.map(x => Math.log2(x.fold)), y: pts.map(x => -Math.log10(Math.max(x.q ?? 1, 1e-6))),
    text: pts.map(x => hit(x) ? x.condition : ""), textposition: "top center",
    hovertext: pts.map(x => `${x.condition}: ${x.fold.toFixed(2)}-fold, q ${x.q == null ? "–" : x.q.toPrecision(2)}`),
    hoverinfo: "text", marker: { size: 10, color: pts.map(x => hit(x) ? "#FF6B57" : "#9DB4C0") } }], L, CONFIG);
}
function scDl() { if (SC) location.href = `/api/screen/${S.job}/screen.csv?metric=${encodeURIComponent(SC.metric)}`; }
