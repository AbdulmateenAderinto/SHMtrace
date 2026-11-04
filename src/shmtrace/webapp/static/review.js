/* Check calls: click through the called mutations, see the reads and the evidence behind
   each, and mark it accepted, rejected or unsure. The marks are stored with the run on the
   server (POST /api/review) and go into variants.csv and the Geneious-style tables; nothing
   is re-called. Each call also carries a suggested verdict with its reasons (review.py), and
   the strip at the top draws every analysed read of the sample against the reference, at base
   resolution, scrolling sideways, with the called mutations in colour. */
"use strict";

const RV = { calls: [], shown: [], at: -1, reads: {}, sug: {} };
const SUG = { accept: ["Accept", "✓"], reject: ["Reject", "✗"], check: ["Check by eye", "?"] };
const rvSug = (v) => RV.sug[rvKey(v)] || null;
const MARK_OF = { accept: "accepted", reject: "rejected" };
const RV_MARKS = { accepted: "Accepted", rejected: "Rejected", unsure: "Unsure" };
const rvKey = (v) => `${v.sample}|${v.pos}|${v.ref}>${v.alt}`;
const rvMark = (v) => ((S.result.review || {})[rvKey(v)] || {}).mark || "";
const pct = (x, d = 2) => `${(100 * x).toFixed(d)}%`;
// v.homopolymer is the longest run of one base touching the call; the caller flags runs this long
const rvHomopolymer = (v) => v.homopolymer >= ((S.result.params || {}).homopolymer_flag || 4);

async function rvInit() {
  if (RV.job === S.job) return;
  RV.job = S.job; RV.reads = {}; RV.sug = {}; AL.sample = null;
  const r = await fetch(`/api/suggestions/${S.job}`);
  if (r.ok) RV.sug = await r.json();
  anFill();
  edLoad();
  S.result.review = S.result.review || {};
  RV.calls = S.result.variants || [];
  const counts = {};
  for (const v of RV.calls) counts[v.sample] = (counts[v.sample] || 0) + 1;
  $("rv-sample").innerHTML = `<option value="">Every sample (${fmtInt(RV.calls.length)})</option>`
    + S.result.samples.map(s => `<option value="${esc(s.sample)}">${esc(s.sample)}${s.control ? " (control)" : ""}`
      + ` — ${fmtInt(counts[s.sample] || 0)}</option>`).join("");
  rvList();
}

function rvList() {
  const sample = $("rv-sample").value, show = $("rv-show").value, sort = $("rv-sort").value;
  let xs = RV.calls.filter(v => !sample || v.sample === sample);
  if (show === "todo") xs = xs.filter(v => !rvMark(v));
  else if (show.startsWith("s-")) xs = xs.filter(v => (rvSug(v) || {}).suggestion === show.slice(2));
  else if (show === "differs") xs = xs.filter(v => {
    const m = rvMark(v), s = (rvSug(v) || {}).suggestion;
    return m && MARK_OF[s] && MARK_OF[s] !== m;
  });
  else if (show) xs = xs.filter(v => rvMark(v) === show);
  const order = S.result.samples.map(s => s.sample);
  xs = [...xs].sort(sort === "reads" ? (a, b) => a.reads - b.reads || a.pos - b.pos
    : sort === "p" ? (a, b) => (b.p_value ?? 0) - (a.p_value ?? 0)
    : (a, b) => order.indexOf(a.sample) - order.indexOf(b.sample) || a.pos - b.pos || a.alt.localeCompare(b.alt));
  RV.shown = xs;
  $("rv-list").innerHTML = xs.length ? xs.map((v, i) => `
    <div class="rv-item" data-i="${i}" onclick="rvOpen(${i})">
      <span class="rv-name">${sample ? "" : `<span class="muted">${esc(v.sample)}</span> `}${v.pos} ${esc(v.ref)}→${esc(v.alt)}</span>
      <span class="rv-meta">${fmtInt(v.reads)} read${v.reads === 1 ? "" : "s"} · ${pct(v.frequency)}</span>
      ${rvSugChip(rvSug(v))}${rvChip(rvMark(v))}
    </div>`).join("") : '<p class="hint" style="padding:8px">No calls to show.</p>';
  rvProgress();
  const keep = RV.current ? xs.findIndex(v => rvKey(v) === rvKey(RV.current)) : -1;
  rvOpen(keep >= 0 ? keep : (xs.length ? 0 : -1));
}
const rvChip = (m) => m ? `<span class="rv-chip ${m}">${RV_MARKS[m]}</span>` : '<span class="rv-chip"></span>';
const rvSugChip = (s) => s ? `<span class="rv-sug ${s.suggestion}" title="Suggested: ${SUG[s.suggestion][0]}">${SUG[s.suggestion][1]}</span>` : "";

function rvProgress() {
  const r = S.result.review || {};
  const n = RV.calls.length, done = RV.calls.filter(v => r[rvKey(v)] && r[rvKey(v)].mark).length;
  const by = (m) => RV.calls.filter(v => rvMark(v) === m).length;
  $("rv-progress").textContent = `${fmtInt(done)} of ${fmtInt(n)} reviewed`
    + (done ? ` · ${by("accepted")} accepted, ${by("rejected")} rejected, ${by("unsure")} unsure` : "");
}

async function rvOpen(i) {
  RV.at = i;
  document.querySelectorAll("#rv-list .rv-item").forEach(el => el.classList.toggle("on", +el.dataset.i === i));
  const v = RV.shown[i];
  RV.current = v;
  if (!v) { $("rv-detail").innerHTML = '<p class="hint">Choose a call.</p>'; return; }
  document.querySelector(`#rv-list .rv-item[data-i="${i}"]`)?.scrollIntoView({ block: "nearest" });
  rvDetail(v, null);
  alShow(v);
  const key = rvKey(v);
  if (!RV.reads[key]) {
    const q = new URLSearchParams({ pos: v.pos, alt: v.alt, flank: 20, others: 30 });
    const r = await fetch(`/api/call-reads/${S.job}/${encodeURIComponent(v.sample)}?${q}`);
    RV.reads[key] = r.ok ? await r.json() : { error: (await r.json().catch(() => ({}))).detail || "the reads could not be loaded" };
  }
  if (RV.current === v) rvDetail(v, RV.reads[key]);
}

/* ── reading frame and regions (coding.py): codons, amino acids, replacement or silent ── */
const EFFECT = { silent: "silent", replacement: "replacement", stop: "stop (nonsense)", "non-coding": "outside the coding stretch" };
function rvCodon(v) {
  const c = (S.result.coding || {})[rvKey(v)];
  if (!c) return "";
  const codon = c.codon_ref ? `${esc(c.codon_ref)} → ${esc(c.codon_alt)}, ${esc(c.aa_ref)}${c.codon_number}${esc(c.aa_alt)} · ` : "";
  return `<dt>Codon</dt><dd>${codon}<b>${EFFECT[c.effect]}</b>${c.region ? ` · in ${esc(c.region)}` : ""}</dd>`;
}
function anFill() {
  const a = S.result.annotation;
  $("an-cds").value = a && a.cds ? `${a.cds[0]}-${a.cds[1]}` : "";
  $("an-regions").value = a ? (a.regions || []).map(r => `${r.name} ${r.start}-${r.end}`).join("\n") : "";
  $("an-summary").textContent = a ? [a.cds && `coding ${a.cds[0]}-${a.cds[1]}`, (a.regions || []).map(r => r.name).join(", ")]
    .filter(Boolean).join(" · ") : "none: give them for a coding amplicon (a V region) to see codons, amino acids and replacement or silent calls";
}
async function anApply(clear) {
  const body = clear ? {} : { cds: $("an-cds").value.trim() || null, regions: $("an-regions").value };
  const r = await fetch(`/api/annotation/${S.job}`, { method: "POST", headers: { "Content-Type": "application/json" },
                                                     body: JSON.stringify(body) });
  const d = await r.json();
  if (!r.ok) { $("an-msg").textContent = typeof d.detail === "string" ? d.detail : "could not use that"; return; }
  S.result.annotation = d.annotation; S.result.coding = d.coding; S.result.rs_by_region = d.rs_by_region;
  const n = Object.values(d.coding).filter(c => c.effect !== "non-coding");
  const rr = n.filter(c => c.effect === "replacement").length, ss = n.filter(c => c.effect === "silent").length;
  $("an-msg").textContent = d.annotation ? `${fmtInt(n.length)} calls in the coding stretch: ${fmtInt(rr)} replacement, `
    + `${fmtInt(ss)} silent, ${fmtInt(n.length - rr - ss)} stop. Saved with the project; see also Tables → Replacement and silent.` : "Cleared.";
  anFill(); alSize(); alDraw();
  if (RV.current) rvDetail(RV.current, RV.reads[rvKey(RV.current)]);
}

/* the facts the detail lists; the verdict and its reasons come from the server (review.py) */
function rvContext(v) {
  // the controls' rate for this change; older runs recorded only the highest at the position
  const alt = S.result.control_background_alt;
  const bg = +((alt ? alt[`${v.pos}${v.alt}`] : (S.result.control_background || {})[v.pos]) || 0);
  const elsewhere = RV.calls.filter(x => x.pos === v.pos && x.alt === v.alt && x.sample !== v.sample);
  return { bg, elsewhere };
}

function rvSuggestion(v) {
  const s = rvSug(v);
  if (!s) return "";
  const icon = (r) => r.decisive ? "✗✗" : r.weight > 0 ? "✓".repeat(r.weight) : r.weight < 0 ? "✗".repeat(-r.weight) : "·";
  return `<div class="rv-suggest ${s.suggestion}">
      <div class="rv-suggest-head"><b>Suggested: ${SUG[s.suggestion][0]}</b>
        <span class="hint">score ${s.score > 0 ? "+" : ""}${s.score}; from how SHMtrace processed this call. Your mark decides.</span></div>
      <ul>${s.reasons.map(r => `<li class="${r.decisive || r.weight < 0 ? "neg" : r.weight > 0 ? "pos" : "neu"}">
        <span class="rv-w">${icon(r)}</span>${esc(r.text)}${r.decisive ? " <b>Decisive.</b>" : ""}</li>`).join("")
        || '<li class="neu"><span class="rv-w">·</span>Nothing for or against beyond passing the test.</li>'}</ul>
    </div>`;
}

function rvDetail(v, d) {
  const mark = rvMark(v), note = ((S.result.review || {})[rvKey(v)] || {}).note || "";
  const { bg, elsewhere } = rvContext(v);
  const ctx = String(v.context || ""), mid = Math.floor(ctx.length / 2);
  const ctxHtml = ctx ? `${esc(ctx.slice(0, mid))}<b class="rv-at">${esc(ctx[mid])}</b>${esc(ctx.slice(mid + 1))}` : "—";
  const q = d && d.carrying ? d.carrying.map(r => (r.diffs.find(x => x[0] === v.pos) || [])[2]).filter(x => x != null) : [];
  $("rv-detail").innerHTML = `
    <div class="rv-head">
      <h3>${esc(v.sample)} · position ${v.pos} · ${esc(v.ref)} → ${esc(v.alt)}</h3>
      <div class="rv-buttons">
        ${Object.entries(RV_MARKS).map(([m, label]) => `<button class="btn small ${mark === m ? "" : "ghost"} rv-${m}"
          onclick="rvSet('${m}')">${label} <kbd>${m[0].toUpperCase()}</kbd></button>`).join("")}
        ${mark ? `<button class="btn small ghost" onclick="rvSet('')">Clear</button>` : ""}
      </div>
    </div>
    ${rvSuggestion(v)}
    <input type="text" id="rv-note" placeholder="Note (optional), saved with the mark" value="${esc(note)}"
           onchange="rvSet(rvMark(RV.current), true)">
    <dl class="kv rv-kv">
      <dt>Supporting reads</dt><dd>${fmtInt(v.reads)} of ${fmtInt(v.depth)} (${pct(v.frequency)}) ·
        ${v.plus_reads} forward, ${v.minus_reads} reverse</dd>
      ${rvFullSet(v)}
      <dt>P</dt><dd>${v.p_value == null ? "—" : v.p_value < 1e-4 ? v.p_value.toExponential(1) : v.p_value.toFixed(4)}</dd>
      <dt>Context</dt><dd class="mono">${ctxHtml}${rvHomopolymer(v) ? ` · ${v.homopolymer}-base homopolymer` : ""}${v.motifs ? ` · ${esc(v.motifs)}` : v.WRC_GYW ? " · WRC/GYW" : ""}</dd>
      <dt>Controls here</dt><dd>${bg ? `error rate ${pct(bg, 3)}, about ${(bg * v.depth).toFixed(1)} reads expected` : "no rate recorded"}</dd>
      <dt>Base quality</dt><dd>${q.length ? q.join(", ") : d ? "—" : "loading…"}</dd>
      ${rvCodon(v)}
      <dt>Same change</dt><dd>${elsewhere.length ? elsewhere.map(x => `${esc(x.sample)} (${x.reads})`).join(", ") : "in no other sample"}</dd>
    </dl>
    ${S.result.params.platform === "sanger" ? '<div id="rv-traces" class="rv-traces"></div>' : ""}
    <div id="rv-pile">${d ? rvPile(v, d) : '<p class="hint">Loading the reads…</p>'}</div>`;
  if (d && d.carrying && S.result.params.platform === "sanger") rvTraces(v, d);
}

/* ── Sanger chromatograms: the peaks under a call, clone by clone ── */
const TRACE_COLOUR = { A: "#1b9e3e", C: "#1f6fb4", G: "#d68a00", T: "#c0392b" };
const TRACES = {};
async function traceWindow(sample, clone, pos) {
  const key = `${sample}|${clone}|${pos}`;
  if (!(key in TRACES)) {
    const q = new URLSearchParams({ clone, pos, flank: 10 });
    const r = await fetch(`/api/chromatogram/${S.job}/${encodeURIComponent(sample)}?${q}`);
    TRACES[key] = r.ok ? await r.json() : { error: (await r.json().catch(() => ({}))).detail || "no trace" };
  }
  return TRACES[key];
}
function drawTrace(cv, w, alt, edits = []) {
  // as wide as the row allows, so the call column sits mid-view without scrolling
  const room = (cv.parentElement?.clientWidth || 600) - (cv.previousElementSibling?.offsetWidth || 0) - 16;
  const per = Math.max(12, Math.min(24, Math.floor(room / w.columns.length))), W = w.columns.length * per, H = 96, top = 16, bottom = 18, dpr = window.devicePixelRatio || 1;
  cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px";
  const ctx = cv.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const css = getComputedStyle(document.body), fg = css.getPropertyValue("--fg").trim(), accent = css.getPropertyValue("--accent").trim();
  ctx.clearRect(0, 0, W, H);
  // spread the peaks evenly, one column per reference position, and draw each channel between them
  const cols = w.columns, known = cols.map((c, i) => [i, c.peak]).filter(([, p]) => p != null);
  const xOfSample = (s) => {                     // piecewise-linear: sample index -> pixel
    for (let k = 1; k < known.length; k++) {
      const [i0, p0] = known[k - 1], [i1, p1] = known[k];
      if (s <= p1) return ((i0 + (s - p0) / Math.max(1, p1 - p0) * (i1 - i0)) + 0.5) * per;
    }
    const [iL, pL] = known[known.length - 1], [iF, pF] = known[0];
    return s < pF ? (iF + 0.5 - (pF - s) / per) * per : (iL + 0.5 + (s - pL) / 12) * per;
  };
  const max = Math.max(1, ...Object.values(w.channels).flat());
  const ci = cols.findIndex(c => c.pos === w.pos);
  if (ci >= 0) { ctx.fillStyle = accent; ctx.globalAlpha = 0.16; ctx.fillRect(ci * per, 0, per, H); ctx.globalAlpha = 1; }
  for (const [base, ys] of Object.entries(w.channels)) {
    ctx.strokeStyle = TRACE_COLOUR[base] || fg; ctx.lineWidth = 1.3; ctx.beginPath();
    ys.forEach((y, s) => { const x = xOfSample(s), yy = H - bottom - (y / max) * (H - top - bottom);
                           s ? ctx.lineTo(x, yy) : ctx.moveTo(x, yy); });
    ctx.stroke();
  }
  ctx.font = "11px ui-monospace, Menlo, monospace"; ctx.textAlign = "center";
  cols.forEach((c, i) => {
    if (!c.base) return;
    ctx.fillStyle = TRACE_COLOUR[c.base] || fg; ctx.globalAlpha = c.pos === w.pos ? 1 : 0.8;
    ctx.font = `${c.pos === w.pos ? "700 " : ""}11px ui-monospace, Menlo, monospace`;
    ctx.fillText(c.base, (i + 0.5) * per, H - 5);
    ctx.globalAlpha = 0.5; ctx.fillStyle = fg; ctx.font = "8px system-ui, sans-serif";
    ctx.fillText(String(c.quality ?? ""), (i + 0.5) * per, 9);
    ctx.globalAlpha = 1;
    const e = edits.find(x => x.pos === c.pos);
    if (e) {                                       // a hand edit: the corrected base over the called one
      ctx.fillStyle = "#c0392b"; ctx.font = "700 11px ui-monospace, Menlo, monospace";
      ctx.fillText("✎" + e.to, (i + 0.5) * per, 21);
    }
  });
}
async function rvTraces(v, d) {
  const el = $("rv-traces"); if (!el) return;
  const clones = d.carrying.slice(0, 6).map(r => r.name);
  el.innerHTML = `<p class="hint">Traces of the clone${clones.length === 1 ? "" : "s"} carrying ${esc(v.ref)}${v.pos}${esc(v.alt)}`
    + `${d.carrying.length > 6 ? ` (6 of ${d.carrying.length})` : ""}, along the reference; the shaded column is the call,
      the small numbers the base qualities. A clean call has one tall peak there; a shoulder or a second peak argues against it.</p>`
    + clones.map((c, k) => `<div class="rv-trace"><span class="hint">${esc(c)}</span><canvas id="rv-tr-${k}"></canvas></div>`).join("");
  for (const [k, c] of clones.entries()) {
    const w = await traceWindow(v.sample, c, v.pos);
    if (RV.current !== v) return;
    const cv = $(`rv-tr-${k}`); if (!cv) return;
    if (w.error) { cv.replaceWith(Object.assign(document.createElement("span"), { className: "hint", textContent: ` — ${w.error}` })); continue; }
    drawTrace(cv, w, v.alt, edFor(v.sample, c));
    cv.classList.add("editable");
    cv.title = (w.reverse ? "read against the reference; shown reversed and complemented. " : "")
      + "Click a base to correct its call.";
    cv.onclick = (ev) => edPick(v.sample, c, w, ev, cv);
  }
}

/* ── hand edits to Sanger base calls (edits.py): recorded here, applied by running again ── */
const ED = { list: [], applied: [], form: null };
const edFor = (sample, clone) => ED.list.filter(e => e.sample === sample && e.clone === clone);
async function edLoad() {
  if (S.result.params.platform !== "sanger") { $("ed-bar").hidden = true; return; }
  const r = await fetch(`/api/edits/${S.job}`);
  if (r.ok) edTake(await r.json());
  edBar();
}
const edTake = (d) => { ED.list = d.edits || []; ED.applied = d.applied || []; };
const edKey = (e) => `${e.sample}|${e.clone}|${e.pos}|${e.to}`;
function edBar() {
  const el = $("ed-bar"), applied = new Set(ED.applied.map(edKey));
  const pending = ED.list.filter(e => !applied.has(edKey(e))).length + ED.applied.filter(e => !ED.list.some(x => edKey(x) === edKey(e))).length;
  el.hidden = !ED.list.length && !ED.applied.length;
  if (el.hidden) return;
  el.innerHTML = `<b>Hand edits to base calls: ${fmtInt(ED.list.length)}</b>
    ${pending ? `<span class="hint"> · ${fmtInt(pending)} not yet in the analysis</span>
      <button class="btn small" style="margin-left:8px" onclick="edRerun()">Run again with the edits</button>`
      : '<span class="hint"> · all in this analysis</span>'}
    <ul>${ED.list.map(e => `<li>${esc(e.sample)} ${esc(e.clone)} position ${e.pos}: trace called ${esc(e.called)}, corrected to
      <b>${esc(e.to)}</b>${e.note ? ` — ${esc(e.note)}` : ""} <span class="muted">${esc(e.when)}</span>
      ${applied.has(edKey(e)) ? "" : " <i class='muted'>(not yet applied)</i>"}
      <a href="#" onclick="edSave('${esc(e.sample)}','${esc(e.clone)}',${e.pos},'');return false">remove</a></li>`).join("")}</ul>
    <p class="hint" id="ed-msg" style="margin:0"></p>`;
}
function edPick(sample, clone, w, ev, cv) {
  const box = cv.getBoundingClientRect(), per = box.width / w.columns.length;
  const c = w.columns[Math.floor((ev.clientX - box.left) / per)];
  if (!c || !c.base) return;
  const cur = edFor(sample, clone).find(e => e.pos === c.pos);
  let form = $("ed-form");
  if (!form) { form = document.createElement("div"); form.id = "ed-form"; form.className = "ed-form"; }
  cv.closest(".rv-trace").after(form);
  form.innerHTML = `<span><b>${esc(clone)}</b>, position ${c.pos}: the trace calls <b>${esc(c.base)}</b> (Q${c.quality ?? "?"})${cur ? `, corrected to <b>${esc(cur.to)}</b>` : ""}.</span>
    <span>Correct to ${["A", "C", "G", "T", "N"].filter(b => b !== c.base).map(b =>
      `<button class="btn small ${cur && cur.to === b ? "" : "ghost"}" onclick="edSave('${esc(sample)}','${esc(clone)}',${c.pos},'${b}')">${b}</button>`).join(" ")}</span>
    <input type="text" id="ed-note" placeholder="why (kept with the edit)" value="${esc(cur ? cur.note : "")}">
    ${cur ? `<button class="btn small ghost" onclick="edSave('${esc(sample)}','${esc(clone)}',${c.pos},'')">Remove edit</button>` : ""}
    <button class="btn small ghost" onclick="this.parentElement.remove()">Close</button>`;
}
async function edSave(sample, clone, pos, to) {
  const note = ($("ed-note") || {}).value || "";
  const r = await fetch(`/api/edits/${S.job}`, { method: "POST", headers: { "Content-Type": "application/json" },
                                               body: JSON.stringify({ sample, clone, pos, to, note }) });
  const d = await r.json();
  if (!r.ok) { const f = $("ed-form"); if (f) f.insertAdjacentHTML("beforeend", `<span class="msg err">${esc(d.detail)}</span>`); return; }
  edTake(d); $("ed-form")?.remove(); edBar();
  if (RV.current) rvTraces(RV.current, RV.reads[rvKey(RV.current)]);
}
// a status line under the edits, made if the bar was redrawn without one
function edStatus(text) {
  let el = $("ed-msg");
  if (!el) { $("ed-bar").hidden = false; $("ed-bar").insertAdjacentHTML("beforeend", '<p class="hint" id="ed-msg" style="margin:0"></p>'); el = $("ed-msg"); }
  el.textContent = text;
}
async function edRerun() {
  if (!confirm("Run the project again with the edited base calls? The new run keeps your review marks for calls that "
    + "still exist, the reading frame and regions, and the list of edits; this run stays as it is under Recent projects.")) return;
  const r = await fetch(`/api/edits/${S.job}/rerun`, { method: "POST" });
  const d = await r.json();
  if (!r.ok) { edStatus(d.detail || "could not start"); return; }
  S.job = d.job;
  const tick = async () => {
    const j = await (await fetch("/api/job/" + S.job)).json();
    edStatus(`Running again with the edits: ${j.message} (${Math.round(100 * j.progress)}%)`);
    if (j.state === "done") { RV.job = null; await loadResult(); showTab("calls"); return; }
    if (j.state === "error") { edStatus(j.error); return; }
    setTimeout(tick, 800);
  };
  tick();
}

/* the reads at the call as text: · where a read matches, the base where it does not */
function rvPile(v, d) {
  if (d.error) return `<p class="hint">${esc(d.error)}</p>`;
  const span = d.end - d.start + 1, col = v.pos - d.start;
  const called = new Set((d.called || []).map(c => c.pos));
  const row = (r) => {
    const cells = Array(span).fill(" ");
    for (let i = 0; i < span; i++) { const p = d.start + i; if (p >= r.start && p <= r.end) cells[i] = "·"; }
    const out = cells.map(c => c === " " ? '<i class="off"> </i>' : "·");
    for (const [p, b, q] of r.diffs) {
      const i = p - d.start; if (i < 0 || i >= span) continue;
      const cls = p === v.pos && b === v.alt ? "alt" : "mm";
      out[i] = `<i class="${cls}" style="opacity:${Math.max(0.35, Math.min(1, q / 40)).toFixed(2)}" title="${p}: ${b}, quality ${q}">${b}</i>`;
    }
    for (const [p, op, x] of r.indels || []) {
      if (op === "D") {                 // a deletion of x bases from p
        for (let i = p - d.start; i < p - d.start + (+x || 1); i++) if (i >= 0 && i < span) out[i] = '<i class="ind" title="deletion">-</i>';
      } else {
        const i = p - d.start;
        if (i >= 0 && i < span) out[i] = `<i class="ind" title="insertion of ${esc(x)} before ${p}">^</i>`;
      }
    }
    out[col] = `<span class="rv-col">${out[col]}</span>`;
    return `<div class="rv-row"><span class="rv-strand">${r.strand < 0 ? "←" : "→"}</span>${out.join("")}</div>`;
  };
  const ruler = Array.from({ length: span }, (_, i) => {
    const p = d.start + i;
    return p === v.pos ? `<span class="rv-col">▼</span>` : called.has(p) ? '<i class="cl" title="another call">•</i>' : " ";
  }).join("");
  const refRow = [...d.sequence].map((b, i) => i === col ? `<span class="rv-col"><b>${b}</b></span>` : b).join("");
  return `
    <p class="hint">${fmtInt(d.carrying.length)} of ${fmtInt(d.covering)} reads covering ${v.pos} carry ${esc(v.alt)};
      ${fmtInt(d.others.length)} of the other ${fmtInt(d.others_total)} are shown below them, spread over the pileup.
      Positions ${d.start}–${d.end}. ${d.other_alt_reads ? `${d.other_alt_reads} other reads read a different base here.` : ""}</p>
    <div class="rv-pileup mono">
      <div class="rv-row"><span class="rv-strand"></span>${ruler}</div>
      <div class="rv-row rv-ref"><span class="rv-strand"></span>${refRow}</div>
      ${d.carrying.map(row).join("")}
      ${d.others.length ? `<div class="rv-sep">reads without ${esc(v.alt)} at ${v.pos}</div>` + d.others.map(row).join("") : ""}
    </div>`;
}

async function rvSet(mark, keepPlace) {
  const v = RV.current; if (!v) return;
  const note = ($("rv-note") || {}).value || "";
  const key = rvKey(v);
  const r = await fetch(`/api/review/${S.job}`, { method: "POST", headers: { "Content-Type": "application/json" },
                                                  body: JSON.stringify({ marks: { [key]: { mark, note } } }) });
  if (!r.ok) { $("rv-progress").textContent = (await r.json().catch(() => ({}))).detail || "could not save the mark"; return; }
  if (!mark && !note) delete S.result.review[key]; else S.result.review[key] = { mark, note };
  const el = document.querySelector(`#rv-list .rv-item[data-i="${RV.at}"] .rv-chip`);
  if (el) el.outerHTML = rvChip(mark);
  rvProgress();
  if (keepPlace || !mark) return rvDetail(v, RV.reads[key]);
  // marking a call moves on to the next one, as a reviewer would
  if ($("rv-show").value === "todo") { rvList(); return; }
  rvOpen(Math.min(RV.at + 1, RV.shown.length - 1));
}

async function rvApplySuggestions() {
  const todo = RV.calls.filter(v => !rvMark(v) && MARK_OF[(rvSug(v) || {}).suggestion]);
  if (!todo.length) { $("rv-progress").textContent = "No unreviewed call has an accept or reject suggestion."; return; }
  const acc = todo.filter(v => rvSug(v).suggestion === "accept").length;
  if (!confirm(`Mark ${todo.length} unreviewed calls as suggested: ${acc} accepted and ${todo.length - acc} rejected? `
    + "Calls suggested for checking by eye, and calls you have already marked, are left alone. "
    + "Each gets the note \"as suggested\", so they can be told apart from calls you looked at.")) return;
  const marks = Object.fromEntries(todo.map(v => [rvKey(v), { mark: MARK_OF[rvSug(v).suggestion], note: "as suggested" }]));
  const r = await fetch(`/api/review/${S.job}`, { method: "POST", headers: { "Content-Type": "application/json" },
                                                  body: JSON.stringify({ marks }) });
  if (!r.ok) { $("rv-progress").textContent = (await r.json().catch(() => ({}))).detail || "could not save the marks"; return; }
  Object.assign(S.result.review, marks);
  rvList();
}

// the call over every usable read, where downsampling drew only some of them (project.full_counts)
function rvFullSet(v) {
  const f = (S.result.full_set || {})[v.sample], c = f && f.calls[`${v.pos}${v.alt}`];
  if (!c) return "";
  return `<dt>All usable reads</dt><dd>${fmtInt(c[0])} of ${fmtInt(c[1])}${c[1] ? ` (${pct(c[0] / c[1])})` : ""}
    <span class="muted">· the analysis drew ${fmtInt(f.analysed)} of the ${fmtInt(f.usable)} usable reads;
    the call and its test rest on those; this counts every usable read, for comparison</span></dd>`;
}

/* ── the alignment strip: every analysed read of a sample, base by base, scrolling sideways ── */
// head: ruler (0-16), reference (20-34), coverage (36-52), identity (54-64)
const AL = { sample: null, data: null, rows: [], cw: 11, rh: 13, head: 70, gutter: 232, cov: [], mis: [],
             called: new Map(), focus: null, picked: null };
// the standard genetic code, as coding.py has it
const CODON_JS = (() => { const b = "TCAG", aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", o = {};
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++) o[b[i] + b[j] + b[k]] = aa[16 * i + 4 * j + k];
  return o; })();
const BASE_COLOUR = { A: "#1b9e3e", C: "#1f6fb4", G: "#d68a00", T: "#c0392b" };
// a read's orientation against the reference, as Geneious colours it: forward and reverse
const STRAND = { 1: { arrow: "▶", word: "forward", colour: "#4a90d9" }, [-1]: { arrow: "◀", word: "reverse", colour: "#d9694a" } };
const strandOf = (r) => STRAND[r.strand < 0 ? -1 : 1];

async function alShow(v) {
  AL.focus = v;
  if (AL.sample !== v.sample) {
    AL.sample = v.sample; AL.data = null; AL.picked = null; alReadPanel();
    $("al-title").textContent = `${v.sample}: reads against the reference`;
    $("al-status").textContent = "Loading the reads…";
    const len = S.result.params.reference_length;
    // with the reads downsampling left out, where the run kept them (full_set), shown on request
    const full = (S.result.full_set || {})[v.sample];
    $("al-undrawn-wrap").hidden = !full;
    const q = new URLSearchParams({ start: 1, end: len, fates: full ? "kept,not sampled" : "kept", limit: 10000 });
    const r = await fetch(`/api/alignment/${S.job}/${encodeURIComponent(v.sample)}?${q}`);
    if (AL.sample !== v.sample) return;
    AL.data = r.ok ? await r.json() : null;
    if (AL.data) alTracks();
    AL.called = new Map();
    for (const c of RV.calls.filter(x => x.sample === v.sample)) {
      if (!AL.called.has(c.pos)) AL.called.set(c.pos, new Set());
      AL.called.get(c.pos).add(c.alt);
    }
  }
  if (!AL.data) { $("al-status").textContent = "No reads were recorded for this sample."; return; }
  alSort(); alSize();               // the scroll box must have its full size before it can scroll
  const sc = $("al-scroll");
  sc.scrollLeft = Math.max(0, (v.pos - AL.data.start + 0.5) * AL.cw - (sc.clientWidth - AL.gutter) / 2);
  sc.scrollTop = 0;
  alDraw(); alReadPanel();
}

// coverage and identity per position, over the analysed reads, as Geneious's graphs above an alignment
function alTracks() {
  const d = AL.data, n = d.end - d.start + 1;
  AL.cov = new Array(n).fill(0); AL.mis = new Array(n).fill(0);
  for (const r of d.reads) {
    if (r.fate !== "kept") continue;          // the graphs describe the analysed reads
    for (let p = Math.max(r.start, d.start); p <= Math.min(r.end, d.end); p++) AL.cov[p - d.start]++;
    for (const [p] of r.diffs) if (p >= d.start && p <= d.end) AL.mis[p - d.start]++;
    for (const [p, op, len] of r.indels || []) if (op === "D")
      for (let q = p; q < p + (+len || 1); q++) if (q >= d.start && q <= d.end) AL.mis[q - d.start]++;
  }
  AL.covMax = Math.max(1, ...AL.cov);
}
const alCarries = (r, v) => v && r.diffs.some(d => d[0] === v.pos && d[1] === v.alt);
function alSort() {
  if (!AL.data) return;
  const v = AL.focus, by = $("al-sort").value;
  const calledIn = (r) => r.diffs.filter(d => AL.called.get(d[0])?.has(d[1])).length;
  const undrawn = $("al-undrawn").checked && !$("al-undrawn-wrap").hidden;
  const order = by === "start" ? (a, b) => a.start - b.start
    : by === "most" ? (a, b) => calledIn(b) - calledIn(a) || a.start - b.start
    : by === "strand" ? (a, b) => (b.strand || 1) - (a.strand || 1) || a.start - b.start
    : (a, b) => alCarries(b, v) - alCarries(a, v) || calledIn(b) - calledIn(a) || a.start - b.start;
  // the analysed reads first, then (when asked for) the ones downsampling did not draw
  AL.rows = AL.data.reads.filter(r => r.fate === "kept" || undrawn)
    .sort((a, b) => (a.fate !== "kept") - (b.fate !== "kept") || order(a, b));
  const n = v ? AL.rows.filter(r => alCarries(r, v)).length : 0;
  // what to show: a name search and a filter, both narrowing the reads already ordered
  const q = ($("al-name").value || "").trim().toLowerCase(), f = $("al-filter").value;
  const hasCall = (r) => r.diffs.some(d => AL.called.get(d[0])?.has(d[1]));
  const all = AL.rows.length;
  AL.rows = AL.rows.filter(r => (!q || r.name.toLowerCase().includes(q))
    && (f === "called" ? hasCall(r) : f === "clean" ? !hasCall(r) : f === "focus" ? alCarries(r, v)
      : f === "fwd" ? (r.strand || 1) > 0 : f === "rev" ? r.strand < 0 : true));
  const fwd = AL.rows.filter(r => (r.strand || 1) > 0).length;
  const nk = AL.data.reads.filter(r => r.fate === "kept").length, nu = AL.data.reads.length - nk;
  const shownU = AL.rows.filter(r => r.fate !== "kept").length;
  $("al-status").innerHTML = (AL.rows.length < all ? `Showing ${fmtInt(AL.rows.length)} of ` : "")
    + `${fmtInt(nk)} analysed reads` + (nu ? (undrawn ? ` + ${fmtInt(nu)} left out by downsampling`
      : ` (${fmtInt(nu)} more left out by downsampling can be shown)`) : "")
    + (AL.data.total > AL.data.reads.length ? ` (the first ${fmtInt(AL.data.reads.length)} of ${fmtInt(AL.data.total)})` : "")
    + (shownU ? `, those faded and after the analysed ones` : "")
    + ` · <span class="al-key" style="color:${STRAND[1].colour}">▶ ${fmtInt(fwd)} forward</span>`
    + ` · <span class="al-key" style="color:${STRAND[-1].colour}">◀ ${fmtInt(AL.rows.length - fwd)} reverse</span>`
    + (v ? ` · ${fmtInt(n)} carry ${esc(v.ref)}${v.pos}${esc(v.alt)}` : "")
    + ". Coloured base: a called mutation; grey: a mismatch SHMtrace did not call. Click a read to see it and extract it; "
    + "click a shaded column in the ruler to open its call.";
  $("al-extract-call").disabled = !n;
  $("al-extract-call").textContent = n ? `Extract the ${fmtInt(n)} read${n === 1 ? "" : "s"} carrying this call` : "Extract the reads carrying this call";
}
function alGoTo() {
  const p = +$("al-goto").value, d = AL.data;
  if (!d || !p) return;
  const sc = $("al-scroll"), at = Math.max(d.start, Math.min(d.end, p));
  sc.scrollLeft = Math.max(0, (at - d.start + 0.5) * AL.cw - (sc.clientWidth - AL.gutter) / 2);
  AL.mark = at; alDraw();
}
function alZoom() {
  const sc = $("al-scroll"), mid = (sc.scrollLeft + (sc.clientWidth - AL.gutter) / 2) / AL.cw;
  AL.cw = +$("al-zoom").value;
  alSize(); sc.scrollLeft = Math.max(0, mid * AL.cw - (sc.clientWidth - AL.gutter) / 2); alDraw();
}
function alSize() {
  if (!AL.data) return;
  // two more rows under the identity graph once a reading frame or regions are given
  AL.head = S.result.annotation ? 98 : 70;
  const span = AL.data.end - AL.data.start + 1;
  const inner = $("al-inner");
  inner.style.width = `${AL.gutter + span * AL.cw + 20}px`;
  inner.style.height = `${AL.head + AL.rows.length * AL.rh + 8}px`;
}
let alPending = false;
function alSchedule() { if (!alPending) { alPending = true; requestAnimationFrame(() => { alPending = false; alDraw(); }); } }

function alDraw() {
  const d = AL.data, sc = $("al-scroll"), cv = $("al-canvas");
  if (!d || !$("page-calls").classList.contains("active")) return;
  const W = sc.clientWidth, H = sc.clientHeight, dpr = window.devicePixelRatio || 1;
  if (cv.width !== W * dpr || cv.height !== H * dpr) {
    cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px";
  }
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const css = getComputedStyle(document.body);
  const fg = css.getPropertyValue("--fg").trim() || "#111", panel = css.getPropertyValue("--panel").trim() || "#fff";
  const line = css.getPropertyValue("--line").trim() || "#ddd", accent = css.getPropertyValue("--accent").trim() || "#1F6F8B";
  const muted = css.getPropertyValue("--muted").trim() || "#888";
  ctx.fillStyle = panel; ctx.fillRect(0, 0, W, H);
  const cw = AL.cw, rh = AL.rh, head = AL.head, G = AL.gutter;
  const p0 = d.start + Math.floor(sc.scrollLeft / cw), p1 = Math.min(d.end, p0 + Math.ceil((W - G) / cw) + 1);
  const x = (p) => G + (p - d.start) * cw - sc.scrollLeft;
  const r0 = Math.floor(sc.scrollTop / rh), r1 = Math.min(AL.rows.length, r0 + Math.ceil((H - head) / rh) + 1);
  const y = (i) => head + i * rh - sc.scrollTop;
  const letters = cw >= 7, v = AL.focus, noise = $("al-noise").checked;

  ctx.save(); ctx.beginPath(); ctx.rect(G, 0, W - G, H); ctx.clip();
  // columns holding calls, then the one being looked at
  for (let p = p0; p <= p1; p++) if (AL.called.has(p)) {
    ctx.fillStyle = accent; ctx.globalAlpha = 0.10; ctx.fillRect(x(p), head - 4, cw, H); ctx.globalAlpha = 1;
  }
  if (v && v.pos >= p0 && v.pos <= p1) {
    ctx.strokeStyle = accent; ctx.lineWidth = 1.5; ctx.strokeRect(x(v.pos) + 0.5, 1, cw - 1, H - 2);
  }
  if (AL.mark && AL.mark >= p0 && AL.mark <= p1) {        // where Go to went
    ctx.strokeStyle = muted; ctx.setLineDash([3, 3]); ctx.lineWidth = 1;
    ctx.strokeRect(x(AL.mark) + 0.5, 1, cw - 1, H - 2); ctx.setLineDash([]);
  }
  // the reads: each bar in its orientation's colour, its mismatches on top
  ctx.save(); ctx.beginPath(); ctx.rect(G, head - 2, W - G, H); ctx.clip();
  ctx.font = `${Math.min(12, cw)}px ui-monospace, Menlo, monospace`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
  for (let i = r0; i < r1; i++) {
    const r = AL.rows[i], yy = y(i), st = strandOf(r), fade = r.fate === "kept" ? 1 : 0.4;
    if (AL.picked === r.name) { ctx.fillStyle = accent; ctx.globalAlpha = 0.22; ctx.fillRect(G, yy, W - G, rh); ctx.globalAlpha = 1; }
    const a = Math.max(r.start, p0), b = Math.min(r.end, p1);
    if (b >= a) {
      ctx.fillStyle = st.colour; ctx.globalAlpha = 0.28 * fade; ctx.fillRect(x(a), yy + 1, (b - a + 1) * cw, rh - 2); ctx.globalAlpha = 1;
      // the read's leading end points the way it was sequenced
      ctx.fillStyle = st.colour;
      if (r.strand < 0 && r.start >= p0) { ctx.beginPath(); ctx.moveTo(x(r.start) - 5, yy + rh / 2); ctx.lineTo(x(r.start), yy + 1); ctx.lineTo(x(r.start), yy + rh - 1); ctx.fill(); }
      if (r.strand >= 0 && r.end <= p1) { ctx.beginPath(); ctx.moveTo(x(r.end + 1) + 5, yy + rh / 2); ctx.lineTo(x(r.end + 1), yy + 1); ctx.lineTo(x(r.end + 1), yy + rh - 1); ctx.fill(); }
    }
    for (const [p, base, q] of r.diffs) {
      if (p < p0 || p > p1) continue;
      const called = AL.called.get(p)?.has(base);
      if (!called && !noise) continue;
      ctx.globalAlpha = fade * (called ? 1 : 0.35 + 0.4 * Math.min(1, q / 40));
      ctx.fillStyle = called ? BASE_COLOUR[base] || fg : "#8a919c";
      ctx.fillRect(x(p), yy + 1, cw, rh - 2);
      if (letters) { ctx.fillStyle = "#fff"; ctx.fillText(base, x(p) + cw / 2, yy + rh / 2 + 0.5); }
      ctx.globalAlpha = 1;
    }
    for (const [p, op, len] of r.indels || []) {
      if (p < p0 || p > p1) continue;
      ctx.fillStyle = fg;
      if (op === "D") ctx.fillRect(x(p), yy + rh / 2 - 1, (+len || 1) * cw, 2);
      else { ctx.fillStyle = "#7d3c98"; ctx.fillRect(x(p) - 1, yy, 2, rh); }
    }
  }
  ctx.restore();
  // the reference and a ruler, pinned at the top
  ctx.fillStyle = panel; ctx.fillRect(G, 0, W - G, head - 4);
  ctx.fillStyle = fg; ctx.font = "10px system-ui, sans-serif"; ctx.textAlign = "left"; ctx.textBaseline = "alphabetic";
  const step = cw >= 7 ? 10 : 50;
  for (let p = Math.ceil(p0 / step) * step; p <= p1; p += step) {
    ctx.globalAlpha = 0.6; ctx.fillText(String(p), x(p), 11); ctx.fillRect(x(p), 13, 1, 4); ctx.globalAlpha = 1;
  }
  ctx.font = `${Math.min(12, cw)}px ui-monospace, Menlo, monospace`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
  for (let p = p0; p <= p1; p++) {
    const b = d.sequence[p - d.start] || "";
    if (letters) { ctx.fillStyle = fg; ctx.globalAlpha = AL.called.has(p) ? 1 : 0.7; ctx.fillText(b, x(p) + cw / 2, 28); }
    else { ctx.fillStyle = BASE_COLOUR[b] || fg; ctx.globalAlpha = 0.6; ctx.fillRect(x(p), 24, cw, 6); }
    ctx.globalAlpha = 1;
  }
  for (let p = p0; p <= p1; p++) {                // coverage, then identity: green full, amber, red
    const k = p - d.start, c = AL.cov[k] || 0, id = c ? 1 - (AL.mis[k] || 0) / c : null;
    const hgt = 15 * c / AL.covMax;
    ctx.fillStyle = muted; ctx.globalAlpha = 0.55; ctx.fillRect(x(p), 52 - hgt, cw, hgt); ctx.globalAlpha = 1;
    if (id != null) {
      ctx.fillStyle = id >= 0.99 ? "#2d9a5a" : id >= 0.9 ? "#d6a400" : "#c0392b";
      ctx.globalAlpha = id >= 0.99 ? 0.55 : 0.9; ctx.fillRect(x(p), 55, cw, 9); ctx.globalAlpha = 1;
    }
  }
  const ann = S.result.annotation;
  if (ann) {                                      // amino acids on their codons, then the named regions
    const ref = S.result.reference || "";
    if (ann.cds) {
      ctx.font = "600 10px ui-monospace, Menlo, monospace"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      for (let s = ann.cds[0] + 3 * Math.max(0, Math.floor((p0 - ann.cds[0]) / 3)); s + 2 <= Math.min(ann.cds[1], p1 + 2); s += 3) {
        const aa = CODON_JS[(ref.slice(s - 1, s + 2) || "").toUpperCase()] || "?";
        const k = (s - ann.cds[0]) / 3;
        ctx.fillStyle = aa === "*" ? "#c0392b" : accent; ctx.globalAlpha = k % 2 ? 0.16 : 0.28;
        ctx.fillRect(x(s), 67, 3 * cw, 12); ctx.globalAlpha = 1;
        ctx.fillStyle = fg; if (cw * 3 >= 9) ctx.fillText(aa, x(s) + 1.5 * cw, 73);
      }
    }
    const pal = ["#6a4c93", "#2a7f8f", "#b07a1e", "#3c7d4f", "#c2185b", "#5d6d7e"];
    (ann.regions || []).forEach((rg, i) => {
      const a = Math.max(rg.start, p0), b = Math.min(rg.end, p1);
      if (b < a) return;
      ctx.fillStyle = pal[i % pal.length]; ctx.globalAlpha = 0.75; ctx.fillRect(x(a), 81, (b - a + 1) * cw, 11); ctx.globalAlpha = 1;
      ctx.fillStyle = "#fff"; ctx.font = "600 9px system-ui, sans-serif"; ctx.textAlign = "left"; ctx.textBaseline = "middle";
      ctx.fillText(rg.name, x(a) + 3, 86.5);
    });
  }
  ctx.restore();
  ctx.fillStyle = line; ctx.fillRect(0, head - 5, W, 1);

  // the side column: which read, in which orientation, under which name (stays put while scrolling)
  ctx.fillStyle = panel; ctx.fillRect(0, 0, G, H);
  ctx.fillStyle = muted; ctx.font = "600 10px system-ui, sans-serif"; ctx.textAlign = "left"; ctx.textBaseline = "alphabetic";
  ctx.fillText("ruler", 6, 11); ctx.fillText("reference", 6, 31);
  ctx.fillText(`coverage (max ${fmtInt(AL.covMax || 0)})`, 6, 48); ctx.fillText("identity", 6, 63);
  if (S.result.annotation) { ctx.fillText("amino acids", 6, 76); ctx.fillText("regions", 6, 90); }

  ctx.save(); ctx.beginPath(); ctx.rect(0, head - 2, G, H); ctx.clip();
  ctx.textBaseline = "middle";
  for (let i = r0; i < r1; i++) {
    const r = AL.rows[i], yy = y(i), st = strandOf(r), picked = AL.picked === r.name;
    if (picked) { ctx.fillStyle = accent; ctx.globalAlpha = 0.22; ctx.fillRect(0, yy, G, rh); ctx.globalAlpha = 1; }
    else if (v && alCarries(r, v)) { ctx.fillStyle = accent; ctx.globalAlpha = 0.08; ctx.fillRect(0, yy, G, rh); ctx.globalAlpha = 1; }
    ctx.font = "10px ui-monospace, Menlo, monospace";
    ctx.fillStyle = muted;
    if (r.fate === "kept") ctx.fillText(String(i + 1), 6, yy + rh / 2);
    else { ctx.font = "italic 8px system-ui, sans-serif"; ctx.fillText("not drawn", 4, yy + rh / 2); }
    ctx.fillStyle = st.colour; ctx.font = "11px system-ui, sans-serif"; ctx.fillText(st.arrow, 46, yy + rh / 2);
    ctx.fillStyle = picked ? fg : muted; ctx.font = `${picked ? "600 " : ""}10px ui-monospace, Menlo, monospace`;
    let name = r.name;
    while (name.length > 4 && ctx.measureText(name).width > G - 70) name = name.slice(0, -2);
    ctx.fillText(name === r.name ? name : name + "…", 62, yy + rh / 2);
  }
  ctx.restore();
  ctx.fillStyle = line; ctx.fillRect(G - 1, 0, 1, H);
}

function alAt(ev) {
  const sc = $("al-scroll"), box = $("al-canvas").getBoundingClientRect();
  const cx = ev.clientX - box.left, cy = ev.clientY - box.top;
  const p = AL.data.start + Math.floor((cx - AL.gutter + sc.scrollLeft) / AL.cw);
  const i = Math.floor((cy - AL.head + sc.scrollTop) / AL.rh);
  return { p, inGutter: cx < AL.gutter, inHead: cy < AL.head, i, r: cy >= AL.head ? AL.rows[i] : null };
}

/* the read picked in the strip: its name, orientation and evidence, and ways to take it away */
function alReadPanel() {
  const el = $("al-read");
  const r = AL.picked && AL.data && AL.rows.find(x => x.name === AL.picked);
  if (!r) { el.innerHTML = '<p class="hint">Click a read in the strip to see its name, orientation and mismatches, and to extract it.</p>'; return; }
  const st = strandOf(r), i = AL.rows.indexOf(r) + 1, drawn = r.fate === "kept";
  const calls = r.diffs.filter(d => AL.called.get(d[0])?.has(d[1]));
  const other = r.diffs.filter(d => !AL.called.get(d[0])?.has(d[1]));
  const sanger = S.result.params.platform === "sanger";
  const callLinks = calls.map(([p, b, q]) => {
    const ref = AL.data.sequence[p - AL.data.start] || "?";
    return `<a href="#" onclick="alOpenCall(${p}, '${b}');return false">${p} ${esc(ref)}→${esc(b)}</a> <span class="muted">Q${q}</span>`;
  }).join(", ");
  el.innerHTML = `
    <div class="al-read-head">
      <b>${drawn ? `Read ${fmtInt(i)} of ${fmtInt(AL.rows.length)}` : "A read downsampling did not draw: not part of the analysis"}</b>
      <code class="al-name" title="the read's name in the input file">${esc(r.name)}</code>
      <button class="btn small ghost" onclick="alCopyName()">Copy name</button>
      <button class="btn small ghost" onclick="alExtract([AL.picked], 'fastq')">Extract read (.fastq)</button>
      <button class="btn small ghost" onclick="alExtract([AL.picked], 'fasta')">.fasta</button>
      ${sanger ? `<a class="btn small ghost" href="/api/trace/${S.job}/${encodeURIComponent(AL.sample)}?clone=${encodeURIComponent(r.name)}">Original trace (.ab1)</a>` : ""}
    </div>
    <dl class="kv rv-kv">
      <dt>Orientation</dt><dd><span style="color:${st.colour};font-weight:700">${st.arrow} ${st.word}</span>
        ${r.strand < 0 ? "· sequenced against the reference; reverse-complemented to align" : "· sequenced along the reference"}</dd>
      <dt>Aligned span</dt><dd>${fmtInt(r.start)}–${fmtInt(r.end)} (${fmtInt(r.end - r.start + 1)} bp of the reference)</dd>
      <dt>Called mutations</dt><dd>${calls.length ? callLinks : "none"}</dd>
      <dt>Not called</dt><dd>${other.length ? `${fmtInt(other.length)} mismatch${other.length === 1 ? "" : "es"}: `
        + other.slice(0, 30).map(([p, b, q]) => `${p}${esc(b)} <span class="muted">Q${q}</span>`).join(", ")
        + (other.length > 30 ? " …" : "") : "none"}</dd>
      ${(r.indels || []).length ? `<dt>Indels</dt><dd>${r.indels.map(([p, op, x]) => op === "D" ? `${p} deletion of ${x}` : `${p} insertion of ${esc(x)}`).join(", ")}</dd>` : ""}
    </dl>
    ${sanger && AL.focus && r.start <= AL.focus.pos && AL.focus.pos <= r.end
      ? `<div class="rv-trace"><span class="hint">trace at ${AL.focus.pos}</span><canvas id="al-read-trace"></canvas></div>` : ""}`;
  if (sanger && AL.focus && $("al-read-trace")) {
    const v = AL.focus, name = r.name;
    traceWindow(AL.sample, name, v.pos).then(w => {
      const cv = $("al-read-trace");
      if (!cv || AL.picked !== name) return;
      if (w.error) cv.replaceWith(Object.assign(document.createElement("span"), { className: "hint", textContent: w.error }));
      else drawTrace(cv, w, v.alt);
    });
  }
}
function alOpenCall(p, alt) {
  let i = RV.shown.findIndex(x => x.sample === AL.sample && x.pos === p && x.alt === alt);
  if (i < 0) { $("rv-show").value = ""; $("rv-sample").value = AL.sample; rvList();
               i = RV.shown.findIndex(x => x.sample === AL.sample && x.pos === p && x.alt === alt); }
  if (i >= 0) rvOpen(i);
}
async function alCopyName() {
  try { await navigator.clipboard.writeText(AL.picked); $("al-status").textContent = `Copied ${AL.picked}.`; }
  catch { $("al-status").textContent = "The browser would not copy; select the name and copy it."; }
}
// the reads as they are in the input file, by name (webapp /api/extract-reads)
async function alExtract(names, format = "fastq") {
  const r = await fetch(`/api/extract-reads/${S.job}/${encodeURIComponent(AL.sample)}`, {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ names, format }) });
  if (!r.ok) { $("al-status").textContent = (await r.json().catch(() => ({}))).detail || "could not extract the reads"; return; }
  const blob = await r.blob(), a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = (r.headers.get("Content-Disposition") || "").match(/filename="([^"]+)"/)?.[1] || `reads.${format}`;
  document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 5000);
  const miss = +r.headers.get("X-Reads-Missing") || 0;
  $("al-status").textContent = `Extracted ${r.headers.get("X-Reads-Found")} read(s) to ${a.download}`
    + (miss ? `; ${miss} could not be found in the input file.` : ".");
}
function alExtractCall() {
  const v = AL.focus;
  if (v) alExtract(AL.rows.filter(r => alCarries(r, v)).map(r => r.name));
}

function rvDownload(dropRejected) {
  location.href = `/api/export/${S.job}/geneious.zip${dropRejected ? "?drop_rejected=true" : ""}`;
}

document.addEventListener("keydown", (e) => {
  if (!$("page-calls").classList.contains("active") || e.metaKey || e.ctrlKey || e.altKey) return;
  if (["INPUT", "SELECT", "TEXTAREA"].includes(document.activeElement.tagName)) return;
  const k = e.key.toLowerCase();
  if (k === "arrowdown" || k === "j") { e.preventDefault(); rvOpen(Math.min(RV.at + 1, RV.shown.length - 1)); }
  else if (k === "arrowup" || k === "k") { e.preventDefault(); rvOpen(Math.max(RV.at - 1, 0)); }
  else if (k === "a") rvSet("accepted");
  else if (k === "r") rvSet("rejected");
  else if (k === "u") rvSet("unsure");
  else if (k === "0") rvSet("");
});

$("al-scroll").addEventListener("scroll", alSchedule);
window.addEventListener("resize", alSchedule);
$("al-canvas").addEventListener("mousemove", (ev) => {
  if (!AL.data) return;
  const { p, r, inGutter } = alAt(ev), d = r && r.diffs.find(x => x[0] === p);
  const calls = AL.called.get(p);
  $("al-canvas").title = r && inGutter ? `read ${AL.rows.indexOf(r) + 1}: ${r.name}\n${strandOf(r).word}, ${r.start}–${r.end}`
    : `position ${p} · reference ${AL.data.sequence[p - AL.data.start] || "?"}`
      + (calls ? ` · called: ${[...calls].join(", ")}` : "")
      + (r ? `\nread ${AL.rows.indexOf(r) + 1}: ${r.name} (${strandOf(r).word}, ${r.start}–${r.end})`
        + (d ? `\nreads ${d[1]}, quality ${d[2]}${AL.called.get(p)?.has(d[1]) ? "" : " (not called)"}` : "") : "");
});
$("al-canvas").addEventListener("click", (ev) => {
  if (!AL.data) return;
  const { p, r, inHead, inGutter } = alAt(ev);
  if (inHead && !inGutter && AL.called.has(p)) {          // the ruler: open that column's call
    const alts = [...AL.called.get(p)];
    return alOpenCall(p, alts.includes(AL.focus?.alt) && AL.focus.pos === p ? AL.focus.alt : alts[0]);
  }
  if (r) { AL.picked = AL.picked === r.name ? null : r.name; alDraw(); alReadPanel(); }
});
