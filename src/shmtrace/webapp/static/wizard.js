/* The set-up guide: four steps over the same state app.js runs a project from.
   Every control app.js reads (#region, #platform, #ref-pick, S.reference, S.samples …)
   is kept, so runProject() is unchanged; this file only adds the steps, the choice
   cards, the checks and the summary. */
"use strict";

const W = { step: 1, ref: null, checks: {}, queue: [], busy: false };

const REGION_TITLES = {
  generic_amplicon: "My own amplicon", jh4_intron: "JH4 intron (mouse)",
  igh_smu: "Sμ switch region (mouse)", jh4_paralogue: "JH4 paralogue",
  n114_vregion: "N114 V region (one hybridoma)", sanger_clones: "Sanger clones",
};
const PLATFORM_TITLES = { ont: "Oxford Nanopore", pacbio: "PacBio HiFi", illumina: "Illumina", sanger: "Sanger" };
const fmt = (n) => Number(n).toLocaleString("en-GB");

/* ── choice cards ─────────────────────────────────────────── */
function wizCards(el, items, current, onPick) {
  el.innerHTML = items.map(([v, title, text, meta]) => `
    <label class="choice${v === current ? " on" : ""}" data-v="${esc(v)}">
      <input type="radio" name="${el.id}" value="${esc(v)}" ${v === current ? "checked" : ""}>
      <b>${esc(title)}</b><span>${esc(text)}</span>${meta ? `<em>${esc(meta)}</em>` : ""}
    </label>`).join("");
  el.querySelectorAll("input").forEach(inp => inp.onchange = () => onPick(inp.value));
}
function presetMeta(name) {
  const d = (S.options.preset_details || {})[name];
  if (!d) return "";
  const parts = [];
  if (d.window) parts.push(`window ${d.window[0]}–${d.window[1]} (${fmt(d.window[1] - d.window[0] + 1)} bp)`);
  else parts.push("whole amplicon");
  parts.push(d.depth ? `${fmt(d.depth)} reads per sample` : "every read");
  return parts.join(" · ");
}
// The general case comes first and is the default; the built-in assays are options, never assumed.
const GENERIC = "generic_amplicon";
const OWN_REFERENCE = [GENERIC, "sanger_clones"];   // regions with no bundled reference
function wizPickRegion(v) {
  $("region").value = v; $("region").onchange && $("region").onchange();
  wizCards($("region-cards"), regionItems(), v, wizPickRegion);
  // until a reference is chosen, step 2 opens on the likely source: your file, or the built-in one
  if (!S.reference) {
    const mode = OWN_REFERENCE.includes(v) ? "own" : "builtin";
    document.querySelector(`input[name="ref-mode"][value="${mode}"]`).checked = true;
    wizRefMode();
  }
  wizRefresh();
}
function wizPickPlatform(v) {
  const was = $("platform").value;
  $("platform").value = v; $("platform").onchange && $("platform").onchange();
  wizCards($("platform-cards"), platformItems(), v, wizPickPlatform);
  if (v === "sanger" && $("region").value !== "sanger_clones" && S.options.presets.includes("sanger_clones")) {
    wizPickRegion("sanger_clones");
  } else if (was === "sanger" && v !== "sanger" && $("region").value === "sanger_clones") {
    wizPickRegion(GENERIC);
  }
  wizRefresh();
}
const REGION_ORDER = [GENERIC, "jh4_intron", "igh_smu", "jh4_paralogue", "n114_vregion", "sanger_clones"];
const firstSentence = (t) => (t.match(/^.*?\.(\s|$)/) || [t])[0].trim();
const regionItems = () => [...S.options.presets]
  .sort((a, b) => (REGION_ORDER.indexOf(a) + 99) % 99 - (REGION_ORDER.indexOf(b) + 99) % 99)
  .map(p => [p, REGION_TITLES[p] || p,
    S.options.preset_help[p] || firstSentence(((S.options.preset_details || {})[p] || {}).description || ""),
    presetMeta(p)]);
const platformItems = () => S.options.platforms.map(p => [p, PLATFORM_TITLES[p] || p,
  (S.options.platform_help || {})[p] || "", ""]);

/* ── reference ────────────────────────────────────────────── */
function wizRefMode() {
  const own = document.querySelector('input[name="ref-mode"]:checked').value === "own";
  $("ref-own-box").hidden = !own; $("ref-builtin-box").hidden = own;
  $("ref-mode-own").classList.toggle("on", own); $("ref-mode-builtin").classList.toggle("on", !own);
}
async function wizCheckRef() {
  const ref = S.reference;
  if (!ref) { W.ref = null; return wizRefresh(); }
  if (W.ref && W.ref.path === ref.path) return wizRefresh();
  W.ref = { path: ref.path, pending: true };
  wizRefresh();
  const r = await fetch("/api/fasta-names", { method: "POST", headers: { "Content-Type": "application/json" },
                                             body: JSON.stringify({ path: ref.path }) });
  const d = await r.json().catch(() => ({}));
  W.ref = r.ok && d.records && d.records.length
    ? { path: ref.path, name: d.records[0].name, length: d.records[0].length, records: d.records.length }
    : { path: ref.path, error: (d.detail || "this file could not be read as a FASTA") };
  wizRefresh();
}

/* ── per-sample read check, one file at a time ────────────── */
function wizQueueChecks() {
  for (const s of S.samples) {
    if (!(s.path in W.checks) && !W.queue.includes(s.path) && typeof s.path === "string") W.queue.push(s.path);
  }
  wizRunQueue();
}
async function wizRunQueue() {
  if (W.busy || !W.queue.length) return;
  W.busy = true;
  const path = W.queue.shift();
  try {
    const r = await fetch("/api/check-reads", { method: "POST", headers: { "Content-Type": "application/json" },
                                               body: JSON.stringify({ path }) });
    W.checks[path] = await r.json();
  } catch (e) { W.checks[path] = { ok: false, error: String(e) }; }
  W.busy = false;
  decorateSamples(); wizRefresh();
  wizRunQueue();
}
function decorateSamples() {
  const rows = $("samples").querySelectorAll(".sample");
  $("sample-head").hidden = !rows.length;
  rows.forEach((row, i) => {
    const s = S.samples[i];
    if (!s) return;
    row.querySelector(".chk")?.remove();
    const c = W.checks[s.path];
    const span = document.createElement("span");
    span.className = "chk";
    if (!c) span.textContent = "checking…";
    else if (!c.ok) { span.textContent = "unreadable"; span.classList.add("bad"); span.title = c.error; }
    else span.textContent = `${fmt(c.reads)}${c.more ? "+" : ""} reads · ${fmt(c.median_length)} bp`;
    row.querySelector(".nm").appendChild(span);
  });
}

/* ── what is wrong, step by step ──────────────────────────── */
function requiredLength() {
  const p = $("platform").value, d = (S.options.preset_details || {})[$("region").value] || {};
  if (p === "illumina") return 100;
  if (p === "sanger") return 200;
  return d.min_align_len || 0;
}
function wizIssues(step) {
  const err = [], warn = [];
  const d = (S.options.preset_details || {})[$("region").value] || {};
  if (step === 1) {
    if (!$("region").value) err.push("Choose the region you amplified.");
    if (!$("platform").value) err.push("Choose how it was sequenced.");
  }
  if (step === 2) {
    if (!S.reference) err.push("Choose a built-in reference or drop your own FASTA.");
    else if (W.ref && W.ref.error) err.push(`The reference could not be read: ${W.ref.error}`);
    else if (W.ref && W.ref.length != null) {
      if (d.window && W.ref.length < d.window[1] && !$("window").value.trim()) {
        warn.push(`The reference is ${fmt(W.ref.length)} bp, but this region analyses positions ${d.window[0]}–${d.window[1]}. `
          + "Pick “My own amplicon” in step 1, or set your own window under More options.");
      }
      if (W.ref.records > 1) warn.push(`The FASTA has ${W.ref.records} records; only the first, ${W.ref.name}, is used.`);
    }
  }
  if (step === 3) {
    const n = S.samples.length;
    if (!n) err.push("Add at least one sample.");
    const names = S.samples.map(s => s.name);
    const dup = names.filter((x, i) => names.indexOf(x) !== i);
    if (dup.length) err.push(`Two samples are called ${[...new Set(dup)].join(", ")}. Rename one file and add it again.`);
    const bad = S.samples.filter(s => W.checks[s.path] && !W.checks[s.path].ok);
    if (bad.length) err.push(`Could not read ${bad.map(s => s.name).join(", ")}. Is it a read file?`);
    const ctrl = S.samples.filter(s => s.control);
    const sanger = $("platform").value === "sanger";
    if (n && !ctrl.length && !sanger) {
      warn.push("No sample is marked as a control. Without one, the background can only be estimated from base "
        + "qualities, which calls more false mutations. Mark your AID-null or empty-vector samples.");
    } else if (ctrl.length === 1 && !sanger) {
      warn.push("Only one control. The background cannot be calibrated from a single control; two or three are much better.");
    }
    const groups = {};
    for (const s of S.samples) (groups[s.group || s.name] ||= []).push(s);
    const mixed = Object.entries(groups).filter(([, xs]) => xs.some(s => s.control) && xs.some(s => !s.control));
    for (const [g] of mixed) warn.push(`Group “${g}” has both controls and samples. Statistics need a group made only of controls.`);
    if (ctrl.length && !Object.values(groups).some(xs => xs.every(s => s.control))) {
      warn.push("No group is made only of controls, so the groups table will have no p-values.");
    }
    const need = requiredLength();
    const short = S.samples.filter(s => W.checks[s.path]?.ok && need && W.checks[s.path].median_length < need);
    if (short.length) warn.push(`Most reads in ${short.map(s => s.name).join(", ")} are shorter than the ${fmt(need)} bp `
      + "this region needs. Check the region and the platform in step 1.");
    const depth = +$("depth").value || d.depth;
    const few = S.samples.filter(s => W.checks[s.path]?.ok && !W.checks[s.path].more && depth && W.checks[s.path].reads < depth);
    if (few.length && !sanger) warn.push(`${few.map(s => `${s.name} (${fmt(W.checks[s.path].reads)})`).join(", ")} `
      + `${few.length === 1 ? "has" : "have"} fewer than ${fmt(depth)} reads in the file. Every sample is analysed at the same depth, `
      + "so the others will be cut to match after filtering.");
    if (W.queue.length || W.busy) warn.push("Still checking the read files…");
  }
  return { err, warn };
}
function allIssues() {
  return [1, 2, 3].map(s => ({ step: s, ...wizIssues(s) }));
}
const stepOk = (s) => !wizIssues(s).err.length;

/* ── navigation ───────────────────────────────────────────── */
function wizGo(step) {
  step = Math.max(1, Math.min(4, step));
  // only forward past steps that are complete
  for (let s = 1; s < step; s++) if (!stepOk(s)) { step = s; break; }
  W.step = step;
  document.querySelectorAll(".wstep").forEach(el => el.hidden = +el.dataset.step !== step);
  window.scrollTo({ top: 0, behavior: "smooth" });
  wizRefresh();
}
function wizRefresh() {
  if (!S.options) return;
  const ok = [1, 2, 3].map(stepOk);
  document.querySelectorAll("#stepper li").forEach(li => {
    const s = +li.dataset.step;
    li.classList.toggle("on", s === W.step);
    li.classList.toggle("done", s < 4 && ok[s - 1] && s !== W.step);
    li.onclick = () => wizGo(s);
  });
  $("wiz-back").style.visibility = W.step > 1 ? "visible" : "hidden";
  $("wiz-next").hidden = W.step === 4;
  if (W.step < 4) {
    const { err, warn } = wizIssues(W.step);
    $("wiz-next").disabled = !!err.length;
    $("wiz-next-hint").textContent = err.length ? err[0] : "";
    $("wiz-issues").innerHTML = warn.map(w => `<div class="msg warn">${esc(w)}</div>`).join("");
  } else {
    $("wiz-next-hint").textContent = "";
    $("wiz-issues").innerHTML = "";
    renderReview();
  }
  if (W.step === 2) renderRefCheck();
}
function renderRefCheck() {
  const r = W.ref;
  $("ref-check").innerHTML = !r ? "" : r.pending ? '<p class="hint">Checking the reference…</p>'
    : r.error ? "" : `<div class="msg ok">${esc(r.name)}: ${fmt(r.length)} bp${r.records > 1 ? `, first of ${r.records} records` : ""}.</div>`;
}

/* ── the summary on step 4 ────────────────────────────────── */
function renderReview() {
  const region = $("region").value, platform = $("platform").value;
  const d = (S.options.preset_details || {})[region] || {};
  const groups = {};
  for (const s of S.samples) (groups[s.group || s.name] ||= []).push(s);
  const ctrlGroups = Object.entries(groups).filter(([, xs]) => xs.length && xs.every(s => s.control)).map(([g]) => g);
  const others = Object.keys(groups).filter(g => !ctrlGroups.includes(g));
  const known = S.samples.map(s => W.checks[s.path]).filter(c => c && c.ok && !c.more).map(c => c.reads);
  const depth = +$("depth").value || d.depth;
  const win = $("window").value.trim() || (d.window ? `${d.window[0]}–${d.window[1]}` : "");
  const span = win ? `positions ${win}`
    : `the whole amplicon${d.end_trim_bp ? ` less ${d.end_trim_bp} bp at each end` : ""}`;
  const rows = [
    ["Region", `${REGION_TITLES[region] || region}: ${span}${depth ? `, up to ${fmt(depth)} reads per sample` : ""}`],
    ["Platform", `${PLATFORM_TITLES[platform] || platform}, calling mode ${$("mode").value}`],
    ["Reference", S.reference ? `${S.reference.name}${W.ref && W.ref.length ? ` (${fmt(W.ref.length)} bp)` : ""}` : "none"],
    ["Set aside", [S.contaminant && `reads that match ${S.contaminant.name}`, S.germlineVcf && `germline variants in ${S.germlineVcf.name}`]
      .filter(Boolean).join("; ") || "nothing"],
    ["Samples", Object.entries(groups).map(([g, xs]) => `${g} (${xs.length}${xs.every(s => s.control) ? ", controls" : ""})`).join(", ") || "none"],
    ["Comparison", ctrlGroups.length && others.length
      ? `${others.join(", ")} against ${ctrlGroups.join(", ")}, each animal or clone counted once (Welch and Mann–Whitney)`
      : "no all-control group, so no p-values"],
  ];
  if (known.length && depth && Math.min(...known) < depth) {
    rows.splice(1, 0, ["Depth", `about ${fmt(Math.min(...known))} reads per sample or fewer (the smallest file), the same for every sample`]);
  }
  const issues = allIssues();
  const errs = issues.flatMap(i => i.err.map(e => [i.step, e]));
  const warns = issues.flatMap(i => i.warn.map(e => [i.step, e]));
  $("review").innerHTML = `<dl class="kv review">${rows.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join("")}</dl>`
    + errs.map(([s, e]) => `<div class="msg err">${esc(e)} <a href="#" onclick="wizGo(${s});return false">Fix in step ${s}</a></div>`).join("")
    + warns.map(([s, e]) => `<div class="msg warn">${esc(e)} <a href="#" onclick="wizGo(${s});return false">Step ${s}</a></div>`).join("")
    + (!errs.length && !warns.length ? '<div class="msg ok">Everything checks out. Press Analyse.</div>' : "");
  $("run-btn").disabled = !!errs.length;
}

/* ── example project ──────────────────────────────────────── */
async function wizExample() {
  const r = await fetch("/api/example");
  if (!r.ok) { $("wiz-issues").innerHTML = `<div class="msg err">${esc((await r.json()).detail)}</div>`; return; }
  const ex = await r.json();
  wizPickRegion(ex.region); wizPickPlatform(ex.platform);
  document.querySelector('input[name="ref-mode"][value="own"]').checked = true; wizRefMode();
  S.reference = ex.reference; $("ref-name").textContent = "Reference: " + ex.reference.name;
  S.contaminant = ex.contaminant; $("contam-name").textContent = "Contaminant: " + ex.contaminant.name;
  S.samples = ex.samples.map(s => ({ ...s, control_kind: s.control ? "aid_null" : "" }));
  $("proj-name").value = ex.name;
  renderSamples(); wizCheckRef();
  wizGo(4);
}

/* ── start ────────────────────────────────────────────────── */
READY.then(() => {
  wizCards($("region-cards"), regionItems(), $("region").value, wizPickRegion);
  wizCards($("platform-cards"), platformItems(), $("platform").value, wizPickPlatform);
  // built-in references: the first option asks for a choice instead of meaning "my own file"
  const first = $("ref-pick").querySelector('option[value=""]');
  if (first) first.textContent = "Choose one…";
  $("ref-pick").addEventListener("change", wizCheckRef);
  wizRefMode();
  // app.js updates these lines when a file is chosen; that is the cue to re-check
  const watch = new MutationObserver(() => { wizCheckRef(); wizRefresh(); });
  for (const id of ["ref-name", "contam-name", "vcf-name", "sheet-msg"]) watch.observe($(id), { childList: true, characterData: true, subtree: true });
  for (const id of ["window", "depth"]) $(id).addEventListener("input", wizRefresh);
  $("mode").addEventListener("change", wizRefresh);
  // Samples added since the last render get their controls guessed straight away, as the guide
  // promises. Only new ones: a box the user unticked stays unticked when more files arrive.
  // (The drop zones hold app.js's addSamples itself, so the hook is here, not there.)
  const CONTROL_NAME = /(^|[^a-z])(ko|knockout|empty|empty_?vector|ev|aicda|aid_?ko|null|control|ctrl|germline|untransfected|mock|uninfected)([^a-z]|$)/i;
  const seen = new Set();
  const rs = renderSamples;
  renderSamples = function () {
    for (const s of S.samples) {
      const key = s.name + "\u0000" + s.path;
      if (!seen.has(key)) { seen.add(key); if (!s.control && CONTROL_NAME.test(s.name)) s.control = true; }
    }
    rs(); decorateSamples(); wizQueueChecks(); wizRefresh();
  };
  wizGo(1);
});
