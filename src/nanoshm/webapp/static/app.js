/* nanoshm desktop app: set up a project, run it, read the results, export. */
"use strict";

const S = {                       // everything the page knows
  options: null, reference: null, contaminant: null,
  samples: [],                    // {name, path, group, control}
  job: null, result: null, table: "samples", sort: null, sortDir: 1,
};
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

/* ── tabs ─────────────────────────────────────────────────── */
function showTab(name) {
  for (const t of ["setup", "results", "reads", "tables", "export"]) {
    $("page-" + t).classList.toggle("active", t === name);
    $("tab-" + t).setAttribute("aria-selected", String(t === name));
  }
  if (name === "results") drawFigures();
  if (name === "reads") rdInit();
}
function enableResultTabs(on) {
  for (const t of ["results", "tables", "export"]) $("tab-" + t).disabled = !on;
  $("tab-reads").disabled = !(on && S.result && S.result.has_alignments);
}

/* ── start up ─────────────────────────────────────────────── */
async function init() {
  const r = await fetch("/api/options");
  S.options = await r.json();
  $("version").textContent = "v" + S.options.version;
  fill($("region"), S.options.presets, "jh4_intron");
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
  dropZone($("drop-contam"), false, (f) => upload(f).then(([x]) => {
    S.contaminant = x; $("contam-name").textContent = "Contaminant: " + x.name;
  }));
  dropZone($("drop-samples"), true, (f) => upload(f).then(addSamples), "pick-folder");
  if (new URLSearchParams(location.search).has("desktop")) $("save-local").style.display = "";
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
    const c = (S.options.contaminants || []).find(x => r && r.id.startsWith(x.for));
    if (c && !S.contaminant) {
      S.contaminant = { name: c.name, path: c.path };
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

function clearSamples() { S.samples = []; renderSamples(); }

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
      <label style="margin:0;display:flex;gap:5px;align-items:center;color:var(--fg)">
        <input type="checkbox" style="width:auto" ${s.control ? "checked" : ""}
               onchange="S.samples[${i}].control=this.checked"> control
      </label>
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
  if ($("call-indels").checked) over.exclude_indels = false;
  const saturation = $("saturation").checked;

  const project = {
    name: $("proj-name").value || "nanoshm project",
    region: $("region").value, platform: $("platform").value, mode: $("mode").value,
    reference: S.reference.path,
    contaminants: S.contaminant ? [S.contaminant.path] : [],
    depth: $("depth").value ? +$("depth").value : null,
    saturation,
    browse: $("keep-alignment").checked,
    overrides: over,
    samples: S.samples.map(s => ({ name: s.name, path: s.path, group: s.group || s.name, control: s.control })),
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
  toImageButtonOptions: { scale: 3, filename: "nanoshm" },
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
    const word = strong ? "clearly above" : (l.fold && l.fold >= 2 ? "above" : "not separated from");
    return `<div class="verdict ${strong ? "yes" : "maybe"}">
      <b>${esc(l.g)}</b> is ${word} <b>${esc(l.ctrl)}</b>:
      ${l.a == null ? "—" : l.a.toFixed(1)} vs ${l.b == null ? "—" : l.b.toFixed(2)} per 10⁵ bp
      ${l.fold ? `· <b>${l.fold.toFixed(0)}-fold</b>` : ""}
      ${l.p != null ? `· P = ${l.p < 1e-4 ? l.p.toExponential(1) : l.p.toFixed(4)}
        <span class="hint">(Welch on log₁₀ densities, n = ${l.n} vs ${ctrlGroup ? ctrlGroup.n : "?"})</span>` : ""}
    </div>`;
  }).join("");

  $("headline").innerHTML = `
    <h2>${esc(r.name)}</h2>
    ${verdict || '<p class="hint">No control group, so there is nothing to compare against.</p>'}
    <dl class="kv">
      <dt>Region</dt><dd>${esc(p.region)} · ${esc(p.reference_name)} ·
        window ${p.analysis_window[0]}–${p.analysis_window[1]} of ${p.reference_length} bp</dd>
      <dt>Data</dt><dd>${esc(p.platform)} · mode ${esc(p.mode)} (${esc(p.error_model)})${p.depth ? " · " + p.depth + " reads per sample" : ""}</dd>
      <dt>Samples</dt><dd>${r.samples.length} in ${r.group_order.length} group(s), ${ctrls.length} control(s)</dd>
    </dl>
    ${(r.notes || []).map(n => `<div class="msg err" style="margin-top:10px">${esc(n)}</div>`).join("")}`;
}

function figShell(id, title, note) {
  return `<div class="fig"><div class="plot" id="${id}"></div>
    <div class="cap"><h4>${esc(title)}</h4>
      <span class="hint">${esc(note || "")}</span></div></div>`;
}

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
    + figShell("f-hot", "AID hotspot targeting", "% of G:C mutations in WRC/GYW");
  plotDots("f-density", "unique_density_x1e5", "Unique mutations / 10⁵ bp");
  plotDots("f-freq", "mutation_frequency_x1e3", "Mutations / 10³ bp");
  plotMatrix();
  plotSpectrum();
  plotLoad();
  plotMap();
  plotDots("f-mutated", "mutated_reads_pct", "% of reads mutated");
  plotDots("f-hot", "pct_GC_in_WRC_GYW", "% of G:C mutations at WRC/GYW");
}

function plotDots(el, key, title) {
  const r = S.result, traces = [], shapes = [], annotations = [];
  const per = {};
  r.group_order.forEach(g => {
    per[g] = r.samples.filter(s => s.group === g && s[key] != null)
      .map(s => ({ name: s.sample, ...valueOf(s, key) }));
  });
  const all = Object.values(per).flat().map(x => x.v);
  if (!all.length) { $(el).parentElement.style.display = "none"; return; }
  const pos = all.filter(v => v > 0);
  const log = pos.length > 0 && (Math.max(...pos) / Math.min(...pos) > 8 || all.some(v => v <= 0));
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
    centre[g] = m;
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
        size: 11, color: rows.map(x => x.censored ? "#ffffff" : colourOf(g)),
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
    const st = pFor(r, g, key, log);
    if (!st || idx(st.ctrl) < 0) return;
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
  if (log) L.yaxis.range = [Math.log10(floor / 2), Math.log10(top * Math.pow(2.0, step + 1.6))];
  Plotly.newPlot(el, traces, L, CONFIG);
}

function plotSpectrum() {
  const r = S.result;
  const subs = [];
  for (const a of "ATCG") for (const b of "ATCG") if (a !== b) subs.push(a + ">" + b);
  const counts = {};
  r.group_order.forEach(g => counts[g] = Object.values(r.spectrum_group[g] || {}).reduce((a, b) => a + b, 0));
  // 1 of 4 calls is "25%": a group with a handful of calls has no spectrum
  let shown = r.group_order.filter(g => counts[g] >= 20);
  const dropped = r.group_order.filter(g => counts[g] < 20);
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
  if (dropped.length && dropped.length < r.group_order.length) {
    L.annotations = [{
      xref: "paper", yref: "paper", x: 0, y: -0.26, showarrow: false, font: { size: 11 },
      text: "Not shown, too few calls for a spectrum: " +
        dropped.map(g => `${g} (n = ${counts[g]})`).join(", "),
    }];
    L.margin = { ...L.margin, b: 76 };
  }
  Plotly.newPlot("f-spectrum", traces, L, CONFIG);
}

function plotMatrix() {
  const r = S.result, bb = ["A", "T", "C", "G"];
  const rows = bb.slice().reverse();   // drawn bottom-up, so A ends up on top
  const n = r.group_order.length;
  const traces = [], annotations = [];
  r.group_order.forEach((g, gi) => {
    const c = r.spectrum_group[g] || {};
    const tot = Object.values(c).reduce((a, b) => a + b, 0);
    const pct = (a, b) => a === b ? null : 100 * (c[`${a}>${b}`] || 0) / (tot || 1);
    const z = rows.map(a => bb.map(b => pct(a, b)));
    const zmax = Math.max(20, ...z.flat().map(v => v || 0));
    traces.push({
      type: "heatmap", z, x: bb, y: rows, xaxis: `x${gi + 1}`, yaxis: `y${gi + 1}`,
      colorscale: "Blues", reversescale: true, showscale: false, zmin: 0, zmax,
      xgap: 2, ygap: 2,
      hovertemplate: `${g}  %{y}→%{x}: %{z:.1f}%<extra></extra>`,
    });
    rows.forEach(a => bb.forEach(b => {
      const v = a === b ? null : (c[`${a}>${b}`] || 0);
      annotations.push({
        xref: `x${gi + 1}`, yref: `y${gi + 1}`, x: b, y: a, showarrow: false,
        text: v === null ? "—"
          : `<b>${v}</b><br><span style="font-size:9px">${(pct(a, b)).toFixed(1)}%</span>`,
        font: { size: 12, color: (pct(a, b) || 0) > 0.55 * zmax ? "#fff" : (dark() ? "#dfe4ec" : "#222") },
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
      marker: { colors: LOAD_SHADES, line: { color: dark() ? "#555" : "#333", width: 1 } },
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
  traces.push({
    type: "scatter", mode: "markers", name: "WRC/GYW", x: r.hotspots,
    y: r.hotspots.map(() => 0), marker: { color: "#9aa3b2", size: 4, symbol: "line-ns-open" },
    hovertemplate: "hotspot %{x}<extra></extra>",
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
function fmtNum(col, v) {
  if (Number.isInteger(v)) return v;
  // a p-value that underflowed to zero is not zero, it is smaller than a double
  if (col === "p_value" && v === 0) return "&lt;1e-308";
  if (v !== 0 && (Math.abs(v) < 1e-4 || Math.abs(v) >= 1e6)) return v.toExponential(3);
  return String(parseFloat(v.toPrecision(4)));
}

function renderTable() {
  const rows = currentRows();
  const thead = $("table").tHead, tbody = $("table").tBodies[0];
  if (!rows.length) { thead.innerHTML = ""; tbody.innerHTML = "<tr><td class=muted>No rows.</td></tr>"; return; }
  // every column any row has: a control without peers lacks the estimates
  const cols = [...new Set(rows.flatMap(r => Object.keys(r)))];
  thead.innerHTML = "<tr>" + cols.map(c =>
    `<th onclick="sortBy('${c}')" class="${typeof rows[0][c] === "number" ? "num" : ""}">${esc(c.replace(/_/g, " "))}${S.sort === c ? (S.sortDir > 0 ? " ▲" : " ▼") : ""}</th>`).join("") + "</tr>";
  tbody.innerHTML = rows.slice(0, 3000).map(r => "<tr>" + cols.map(c => {
    let v = r[c];
    if (v === true) v = "<span class='pill'>yes</span>";
    else if (v === false) v = "";
    else if (v == null) v = "";
    else if (typeof v === "number") v = fmtNum(c, v);
    else v = esc(v);
    return `<td class="${typeof r[c] === "number" ? "num" : ""}">${v}</td>`;
  }).join("") + "</tr>").join("");
  $("table-count").textContent = `${rows.length} row${rows.length === 1 ? "" : "s"}${rows.length > 3000 ? " (showing 3000)" : ""}`;
}
function sortBy(c) { S.sortDir = S.sort === c ? -S.sortDir : 1; S.sort = c; renderTable(); }
function downloadCurrentTable() {
  const map = { samples: "samples.csv", groups: "groups.csv", variants: "variants.csv",
    indels: "indels.csv",
    positional: "positional_frequency.csv", spectrum: "spectrum.csv",
    mutation_load: "mutation_load.csv" };
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
init();


/* ── reads: the alignment, before and after processing ────── */
// Colours run from "analysed" through the stages that drop a read, in the
// order the pipeline applies them, so the legend reads like the pipeline.
const FATE_COLOUR = {
  "kept": "#3F6FA8",
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
  else if (show === "dropped") fates = Object.keys(FATE_COLOUR).filter(f => f !== "kept").join(",");
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
    ctx.globalAlpha = rd.fate === "kept" ? 0.3 : 0.22;
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

  const c = d.counts, kept = c.kept || 0;
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
    <p style="margin:0 0 4px"><strong>Mismatches:</strong> ${esc(diffs)}</p>
    <p style="margin:0"><strong>Indels:</strong> ${esc(ind)}</p>`;
}
