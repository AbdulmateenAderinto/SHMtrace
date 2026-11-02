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
      <dt>P</dt><dd>${v.p_value == null ? "—" : v.p_value < 1e-4 ? v.p_value.toExponential(1) : v.p_value.toFixed(4)}</dd>
      <dt>Context</dt><dd class="mono">${ctxHtml}${rvHomopolymer(v) ? ` · ${v.homopolymer}-base homopolymer` : ""}${v.WRC_GYW ? " · WRC/GYW" : ""}</dd>
      <dt>Controls here</dt><dd>${bg ? `error rate ${pct(bg, 3)}, about ${(bg * v.depth).toFixed(1)} reads expected` : "no rate recorded"}</dd>
      <dt>Base quality</dt><dd>${q.length ? q.join(", ") : d ? "—" : "loading…"}</dd>
      <dt>Same change</dt><dd>${elsewhere.length ? elsewhere.map(x => `${esc(x.sample)} (${x.reads})`).join(", ") : "in no other sample"}</dd>
    </dl>
    <div id="rv-pile">${d ? rvPile(v, d) : '<p class="hint">Loading the reads…</p>'}</div>`;
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

/* ── the alignment strip: every analysed read of a sample, base by base, scrolling sideways ── */
const AL = { sample: null, data: null, rows: [], cw: 11, rh: 12, head: 40, called: new Map(), focus: null };
const BASE_COLOUR = { A: "#1b9e3e", C: "#1f6fb4", G: "#d68a00", T: "#c0392b" };

async function alShow(v) {
  AL.focus = v;
  if (AL.sample !== v.sample) {
    AL.sample = v.sample; AL.data = null;
    $("al-title").textContent = `${v.sample}: reads against the reference`;
    $("al-status").textContent = "Loading the reads…";
    const len = S.result.params.reference_length;
    const q = new URLSearchParams({ start: 1, end: len, fates: "kept", limit: 2000 });
    const r = await fetch(`/api/alignment/${S.job}/${encodeURIComponent(v.sample)}?${q}`);
    if (AL.sample !== v.sample) return;
    AL.data = r.ok ? await r.json() : null;
    AL.called = new Map();
    for (const c of RV.calls.filter(x => x.sample === v.sample)) {
      if (!AL.called.has(c.pos)) AL.called.set(c.pos, new Set());
      AL.called.get(c.pos).add(c.alt);
    }
  }
  if (!AL.data) { $("al-status").textContent = "No reads were recorded for this sample."; return; }
  alSort(); alSize();               // the scroll box must have its full size before it can scroll
  const sc = $("al-scroll");
  sc.scrollLeft = Math.max(0, (v.pos - AL.data.start + 0.5) * AL.cw - sc.clientWidth / 2);
  sc.scrollTop = 0;
  alDraw();
}

function alSort() {
  if (!AL.data) return;
  const v = AL.focus, by = $("al-sort").value;
  const carries = (r) => v && r.diffs.some(d => d[0] === v.pos && d[1] === v.alt);
  const calledIn = (r) => r.diffs.filter(d => AL.called.get(d[0])?.has(d[1])).length;
  AL.rows = [...AL.data.reads].sort(by === "start" ? (a, b) => a.start - b.start
    : by === "most" ? (a, b) => calledIn(b) - calledIn(a) || a.start - b.start
    : (a, b) => carries(b) - carries(a) || calledIn(b) - calledIn(a) || a.start - b.start);
  const n = v ? AL.rows.filter(carries).length : 0;
  $("al-status").textContent = `${fmtInt(AL.rows.length)} analysed reads of ${fmtInt(AL.data.total)}`
    + (AL.data.total > AL.rows.length ? " (the first 2,000 by position)" : "")
    + ` over ${AL.data.start}–${AL.data.end}`
    + (v ? ` · ${fmtInt(n)} carry ${v.ref}${v.pos}${v.alt}` : "")
    + ". Coloured: a called mutation; grey: a mismatch SHMtrace did not call (sequencing error or background). "
    + "Shaded columns hold calls; click one to open it.";
}
function alZoom() {
  const sc = $("al-scroll"), mid = (sc.scrollLeft + sc.clientWidth / 2) / AL.cw;
  AL.cw = +$("al-zoom").value;
  alSize(); sc.scrollLeft = Math.max(0, mid * AL.cw - sc.clientWidth / 2); alDraw();
}
function alSize() {
  if (!AL.data) return;
  const span = AL.data.end - AL.data.start + 1;
  const inner = $("al-inner");
  inner.style.width = `${span * AL.cw + 20}px`;
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
  ctx.fillStyle = panel; ctx.fillRect(0, 0, W, H);
  const cw = AL.cw, rh = AL.rh, head = AL.head;
  const p0 = d.start + Math.floor(sc.scrollLeft / cw), p1 = Math.min(d.end, p0 + Math.ceil(W / cw) + 1);
  const x = (p) => (p - d.start) * cw - sc.scrollLeft;
  const r0 = Math.floor(sc.scrollTop / rh), r1 = Math.min(AL.rows.length, r0 + Math.ceil((H - head) / rh) + 1);
  const y = (i) => head + i * rh - sc.scrollTop;
  const letters = cw >= 7, v = AL.focus, noise = $("al-noise").checked;

  // columns holding calls, then the one being looked at
  for (let p = p0; p <= p1; p++) if (AL.called.has(p)) {
    ctx.fillStyle = accent; ctx.globalAlpha = 0.10; ctx.fillRect(x(p), head - 4, cw, H); ctx.globalAlpha = 1;
  }
  if (v && v.pos >= p0 && v.pos <= p1) {
    ctx.strokeStyle = accent; ctx.lineWidth = 1.5; ctx.strokeRect(x(v.pos) + 0.5, 1, cw - 1, H - 2);
  }
  // the reads
  ctx.save(); ctx.beginPath(); ctx.rect(0, head - 2, W, H); ctx.clip();
  ctx.font = `${Math.min(12, cw)}px ui-monospace, Menlo, monospace`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
  for (let i = r0; i < r1; i++) {
    const r = AL.rows[i], yy = y(i);
    const a = Math.max(r.start, p0), b = Math.min(r.end, p1);
    if (b >= a) { ctx.fillStyle = line; ctx.fillRect(x(a), yy + 1, (b - a + 1) * cw, rh - 2); }
    for (const [p, base, q] of r.diffs) {
      if (p < p0 || p > p1) continue;
      const called = AL.called.get(p)?.has(base);
      if (!called && !noise) continue;
      ctx.globalAlpha = called ? 1 : 0.35 + 0.4 * Math.min(1, q / 40);
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
  ctx.fillStyle = panel; ctx.fillRect(0, 0, W, head - 4);
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
  ctx.fillStyle = line; ctx.fillRect(0, head - 5, W, 1);
}

function alAt(ev) {
  const sc = $("al-scroll"), box = $("al-canvas").getBoundingClientRect();
  const p = AL.data.start + Math.floor((ev.clientX - box.left + sc.scrollLeft) / AL.cw);
  const i = Math.floor((ev.clientY - box.top - AL.head + sc.scrollTop) / AL.rh);
  return { p, r: ev.clientY - box.top >= AL.head ? AL.rows[i] : null };
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
  const { p, r } = alAt(ev), d = r && r.diffs.find(x => x[0] === p);
  const calls = AL.called.get(p);
  $("al-canvas").title = `position ${p} · reference ${AL.data.sequence[p - AL.data.start] || "?"}`
    + (calls ? ` · called: ${[...calls].join(", ")}` : "")
    + (r ? `\nread ${r.name} (${r.start}–${r.end})` + (d ? ` reads ${d[1]}, quality ${d[2]}${AL.called.get(p)?.has(d[1]) ? "" : " (not called)"}` : "") : "");
});
$("al-canvas").addEventListener("click", (ev) => {
  if (!AL.data) return;
  const { p, r } = alAt(ev);
  if (!AL.called.has(p)) return;
  const base = r && r.diffs.find(x => x[0] === p)?.[1];
  const alt = base && AL.called.get(p).has(base) ? base : [...AL.called.get(p)][0];
  let i = RV.shown.findIndex(x => x.sample === AL.sample && x.pos === p && x.alt === alt);
  if (i < 0) { $("rv-show").value = ""; $("rv-sample").value = AL.sample; rvList();
               i = RV.shown.findIndex(x => x.sample === AL.sample && x.pos === p && x.alt === alt); }
  if (i >= 0) rvOpen(i);
});
