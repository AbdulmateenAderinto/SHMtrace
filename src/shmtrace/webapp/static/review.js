/* Check calls: click through the called mutations, see the reads and the evidence behind
   each, and mark it accepted, rejected or unsure. The marks are stored with the run on the
   server (POST /api/review) and go into variants.csv and the Geneious-style tables; nothing
   is re-called. */
"use strict";

const RV = { calls: [], shown: [], at: -1, reads: {} };
const RV_MARKS = { accepted: "Accepted", rejected: "Rejected", unsure: "Unsure" };
const rvKey = (v) => `${v.sample}|${v.pos}|${v.ref}>${v.alt}`;
const rvMark = (v) => ((S.result.review || {})[rvKey(v)] || {}).mark || "";
const pct = (x, d = 2) => `${(100 * x).toFixed(d)}%`;
// v.homopolymer is the longest run of one base touching the call; the caller flags runs this long
const rvHomopolymer = (v) => v.homopolymer >= ((S.result.params || {}).homopolymer_flag || 4);

function rvInit() {
  if (RV.job === S.job) return;
  RV.job = S.job; RV.reads = {};
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
      ${rvChip(rvMark(v))}
    </div>`).join("") : '<p class="hint" style="padding:8px">No calls to show.</p>';
  rvProgress();
  const keep = RV.current ? xs.findIndex(v => rvKey(v) === rvKey(RV.current)) : -1;
  rvOpen(keep >= 0 ? keep : (xs.length ? 0 : -1));
}
const rvChip = (m) => m ? `<span class="rv-chip ${m}">${RV_MARKS[m]}</span>` : '<span class="rv-chip"></span>';

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
  const key = rvKey(v);
  if (!RV.reads[key]) {
    const q = new URLSearchParams({ pos: v.pos, alt: v.alt, flank: 20, others: 30 });
    const r = await fetch(`/api/call-reads/${S.job}/${encodeURIComponent(v.sample)}?${q}`);
    RV.reads[key] = r.ok ? await r.json() : { error: (await r.json().catch(() => ({}))).detail || "the reads could not be loaded" };
  }
  if (RV.current === v) rvDetail(v, RV.reads[key]);
}

/* what a careful person would check, worked out from the call and its reads */
function rvFlags(v, d) {
  const flags = [], bg = +((S.result.control_background || {})[v.pos] || 0);
  const sample = S.result.samples.find(s => s.sample === v.sample) || {};
  if (sample.control) flags.push(["warn", "This sample is a control: a call here is background unless it is a real germline difference."]);
  // only telling when the sample has reads on both strands: some sets arrive already oriented
  const bothStrands = RV.calls.some(x => x.sample === v.sample && x.minus_reads > 0)
    && RV.calls.some(x => x.sample === v.sample && x.plus_reads > 0);
  if (bothStrands && v.reads >= 3 && (!v.plus_reads || !v.minus_reads))
    flags.push(["warn", `Every supporting read is on one strand (${v.plus_reads} forward, ${v.minus_reads} reverse).`]);
  if (rvHomopolymer(v)) flags.push(["warn", `Next to a run of ${v.homopolymer} identical bases, where nanopore reads miscall most often.`]);
  if (bg && bg * v.depth >= v.reads)
    flags.push(["warn", `The controls' error rate here (${pct(bg, 3)}) predicts ${(bg * v.depth).toFixed(1)} reads at this depth; the call has ${v.reads}.`]);
  const elsewhere = RV.calls.filter(x => x.pos === v.pos && x.alt === v.alt && x.sample !== v.sample);
  const ctrlNames = new Set(S.result.samples.filter(s => s.control).map(s => s.sample));
  const inCtrl = elsewhere.filter(x => ctrlNames.has(x.sample));
  if (inCtrl.length && !sample.control)
    flags.push(["warn", `Also called in ${inCtrl.length === 1 ? "a control" : inCtrl.length + " controls"} (${inCtrl.map(x => x.sample).join(", ")}).`]);
  if (d && d.carrying) {
    const q = d.carrying.map(r => (r.diffs.find(x => x[0] === v.pos) || [])[2]).filter(x => x != null);
    const low = q.filter(x => x < 20).length;
    if (low) flags.push(["warn", `${low} of ${q.length} supporting bases have quality below 20.`]);
    const nearEnd = d.carrying.filter(r => v.pos - r.start < 10 || r.end - v.pos < 10).length;
    if (nearEnd) flags.push(["warn", `${nearEnd} supporting read${nearEnd === 1 ? " ends" : "s end"} within 10 bp of the call.`]);
  }
  if (v.WRC_GYW) flags.push(["ok", "At a WRC/GYW hotspot, where AID prefers to act."]);
  if (!flags.some(f => f[0] === "warn")) flags.push(["ok", "Nothing in the reads argues against this call."]);
  return { flags, bg, elsewhere };
}

function rvDetail(v, d) {
  const mark = rvMark(v), note = ((S.result.review || {})[rvKey(v)] || {}).note || "";
  const { flags, bg, elsewhere } = rvFlags(v, d);
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
    <div class="rv-flags">${flags.map(([k, t]) => `<div class="msg ${k === "ok" ? "ok" : "warn"}">${esc(t)}</div>`).join("")}</div>
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
