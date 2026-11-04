/* Saved projects: every finished run is kept on this computer with its review marks
   (session.py), listed on Set up to reopen, and can be saved as one .shmtrace file. */
"use strict";

async function pjRecent() {
  const r = await fetch("/api/projects");
  if (!r.ok) return;
  const d = await r.json(), el = $("pj-recent");
  el.hidden = !d.projects.length;
  if (!d.projects.length) return;
  const when = (s) => new Date(s).toLocaleString("en-GB", { day: "numeric", month: "short", hour: "2-digit", minute: "2-digit" });
  el.innerHTML = `<div class="row" style="align-items:baseline"><b>Recent projects</b>
      <span class="hint" title="${esc(d.folder)}">saved on this computer, with their review marks</span></div>
    <div class="pj-list">${d.projects.slice(0, 8).map(p => `
      <button class="pj-item" onclick="pjOpen('${esc(p.folder)}')" title="${esc(p.folder)}">
        <b>${esc(p.name)}</b>
        <span class="hint">${when(p.saved)} · ${fmtInt(p.samples)} samples · ${fmtInt(p.calls)} calls`
        + `${p.calls ? ` · ${fmtInt(p.reviewed)} reviewed` : ""}${p.platform ? ` · ${esc(p.platform)}` : ""}</span>
      </button>`).join("")}</div>`;
}

async function pjShow(r) {
  if (!r.ok) {
    const d = await r.json().catch(() => ({}));
    showTab("setup");
    $("wiz-issues").innerHTML = `<div class="msg err">${esc(typeof d.detail === "string" ? d.detail : "could not open the project")}</div>`;
    return;
  }
  S.job = (await r.json()).job;
  await loadResult();
}
async function pjOpen(folder) {
  await pjShow(await fetch("/api/projects/open", { method: "POST", headers: { "Content-Type": "application/json" },
                                                   body: JSON.stringify({ folder }) }));
}
async function pjUpload(input) {
  const f = input.files[0]; input.value = "";
  if (!f) return;
  const body = new FormData(); body.append("file", f);
  await pjShow(await fetch("/api/projects/upload", { method: "POST", body }));
}
function pjSave() {
  const a = document.createElement("a");
  a.href = `/api/project-file/${S.job}`;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(pjAfterLoad, 1500);       // the download saved it too (a sample-data run was not kept before)
}
async function pjAfterLoad() {
  const r = await fetch(`/api/job/${S.job}`);
  const j = r.ok ? await r.json() : {};
  $("pj-saved").textContent = j.saved ? "Saved on this computer with its review marks; reopen it from Set up."
    : j.example ? "Sample data: not kept, so it runs afresh each time. Save project file keeps it." : "";
  pjRecent();
}

READY.then(pjRecent);
