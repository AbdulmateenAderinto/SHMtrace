"""HTTP layer for the shmtrace desktop app (and the same thing served locally).

The GUI is a single page talking to this API. A project (many samples, groups,
controls, one reference) is submitted once, runs in a background thread, and
the page polls for progress; the finished result is kept in memory and every
export is a download from it.

Endpoints
  GET  /                      the app
  GET  /api/options           presets, platforms, modes, version
  POST /api/files             upload files, returns server-side paths
  POST /api/local-scan        list a folder on this machine (desktop app)
  POST /api/run               start a project; returns {"job": id}
  GET  /api/job/{id}          {"state","progress","message","error"}
  GET  /api/result/{id}       the full project result (JSON)
  GET  /api/export/{id}/{what} report.html | tables.xlsx | bundle.zip |
                              <table>.csv | <figure>.png|svg|pdf
  POST /api/fasta-names       record names of a FASTA (the switch tab's donor picker)
  POST /api/switch/run        class-switch junctions for several samples; {"job": id}
  GET  /api/switch/result/{id}  per-sample summaries, comparison and every junction
  GET  /api/switch/export/{id}/{what}  summary.csv | comparison.csv | junctions.csv
  POST /api/repair/{id}       repair-route indices per animal, a group comparison and a
                              spectrum decomposition against reference groups
  POST /api/power             power by simulation, from a finished project or given numbers
  POST /api/screen/{id}       every condition against the control, ranked, with BH q-values
                              (pooled runs in a sample sheet are split by barcode at /api/run)
  POST /api/compare           the older two-group comparison (kept)
  POST /api/analyze           the older single-sample endpoint (kept)
"""
from __future__ import annotations

import base64
import json
import shutil
import tempfile
import threading
import uuid
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from .. import __version__, browse, pipeline, project, report, viz
from .. import compare as compare_mod
from ..config import available_presets, load_preset

APP_DIR = Path(__file__).resolve().parent
STATIC = APP_DIR / "static"
MAX_UPLOAD_MB = 2000

app = FastAPI(title=f"shmtrace {__version__}")

UPLOADS = Path(tempfile.mkdtemp(prefix="shmtrace_uploads_"))
JOBS: dict[str, dict] = {}
_LOCK = threading.Lock()


def _finite(x):
    """JSON has no NaN or infinity: an undefined statistic is sent as null."""
    import math
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {k: _finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_finite(v) for v in x]
    if hasattr(x, "item"):                       # numpy scalars
        return _finite(x.item())
    return x


# ---- pages ------------------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
def index():
    return (STATIC / "index.html").read_text()


@app.get("/classic", response_class=HTMLResponse)
def classic():
    return (STATIC / "classic.html").read_text()


def _contaminant_references() -> list[dict]:
    d = APP_DIR.parent / "references"
    return [{"id": f.stem, "name": f.read_text().split("\n", 1)[0].lstrip(">").strip(),
             "path": str(f), "for": f.stem.replace("_contaminant", "")}
            for f in sorted(d.glob("*contaminant*.fasta"))]


def _bundled_references() -> list[dict]:
    d = APP_DIR.parent / "references"
    out = []
    for f in sorted(d.glob("*.fasta")):
        if "contaminant" in f.stem:        # not something to analyse against
            continue
        head = f.read_text().split("\n", 1)[0].lstrip(">").strip()
        length = sum(len(x.strip()) for x in f.read_text().split("\n")[1:]
                     if x and not x.startswith(">"))
        out.append({"id": f.stem, "name": head, "length": length, "path": str(f)})
    return out


@app.get("/api/options")
def options():
    return {
        "version": __version__,
        "presets": available_presets(),
        "references": _bundled_references(),
        "contaminants": _contaminant_references(),
        "platforms": list(getattr(pipeline, "PLATFORMS", ("ont",))),
        "modes": ["recommended", "sensitive", "geneious", "v0.1", "strand", "pbinom", "qcal+filter",
                  "sanger"],
        "preset_help": {
            "jh4_intron": "Mouse JH4 downstream intron, the neutral polyclonal assay.",
            "n114_vregion": "N114 hybridoma V region; the parental sequence is known.",
            "igh_smu": "Mouse IgH Sμ switch region.",
            "generic_amplicon": "Any other region: JH6, Sγ, a V gene, an AID off-target.",
            "sanger_clones": "Sanger traces, one clone per read.",
        },
        "mode_help": {
            "recommended": "Quality-calibrated against your controls. The default.",
            "sensitive": "Recommended at a false-discovery rate of 0.20: finds more mutations seen "
                         "in a single read, but about one call in five is background.",
            "geneious": "Reproduces Geneious's rule (uncorrected p ≤ 0.05).",
            "v0.1": "The previous default: per-position control rates, BH-corrected.",
            "strand": "Strand-aware; each strand tested against its own control rate.",
            "pbinom": "Exact test on base qualities; no controls needed.",
            "qcal+filter": "Recommended, plus a strand-bias filter.",
            "sanger": "Every high-quality mismatch in a clone is a mutation.",
        },
    }


# ---- files ------------------------------------------------------------------
async def _save(upload: UploadFile, dest: Path) -> Path:
    size = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as fh:
        while chunk := await upload.read(1 << 20):
            size += len(chunk)
            if size > MAX_UPLOAD_MB * (1 << 20):
                raise ValueError(f"{upload.filename} is over {MAX_UPLOAD_MB} MB")
            fh.write(chunk)
    return dest


def _safe_rel(rel: str) -> Path:
    """A browser-supplied relative path, stripped of anything that could escape."""
    parts = [p for p in Path(rel.replace("\\", "/")).parts
             if p not in ("", ".", "..", "/") and not p.startswith("~")]
    return Path(*parts[-3:]) if parts else Path("file")


@app.post("/api/files")
async def upload_files(files: list[UploadFile] = File(...),
                       paths: list[str] = Form(default=[])):
    """Store uploads and hand back the paths the project will refer to.

    `paths` carries each file's path inside a dropped folder, so a folder of
    Sanger traces stays a folder on this side and can be read as one sample.
    """
    out = []
    folder = UPLOADS / uuid.uuid4().hex[:8]
    for i, up in enumerate(files):
        rel = _safe_rel(paths[i] if i < len(paths) and paths[i] else (up.filename or "file"))
        dest = await _save(up, folder / rel)
        out.append({"name": dest.name, "rel": str(rel), "path": str(dest),
                    "folder": str(dest.parent) if rel.parent != Path(".") else "",
                    "bytes": dest.stat().st_size})
    return {"files": out}


@app.post("/api/local-scan")
def local_scan(payload: dict):
    """List sequence files in a folder on this machine (desktop app only)."""
    root = Path(payload.get("path", "")).expanduser()
    if not root.is_dir():
        raise HTTPException(400, f"not a folder: {root}")
    exts = set(getattr(__import__("shmtrace.io_utils", fromlist=["x"]), "READ_EXTENSIONS",
                       (".fastq", ".fq", ".fasta", ".fa", ".bam", ".ab1")))
    files = []
    for p in sorted(root.iterdir()):
        if p.is_file() and (p.suffix.lower() in exts
                            or (p.suffix.lower() == ".gz" and p.with_suffix("").suffix.lower() in exts)):
            files.append({"name": p.name, "path": str(p), "bytes": p.stat().st_size})
    return {"folder": str(root), "files": files}


@app.post("/api/samplesheet")
def load_samplesheet(payload: dict):
    """Read a sample sheet on this machine: its samples, warnings and exclusions.

    The sheet's file paths are local (absolute, ~, or relative to the sheet), so
    this is for the desktop app or a server reading its own disk.
    """
    from .. import samplesheet
    path = Path(payload.get("path", "")).expanduser()
    if not path.is_file():
        raise HTTPException(400, f"no sample sheet at {path}")
    try:
        sh = samplesheet.load(str(path))
    except samplesheet.SheetError as e:
        raise HTTPException(400, {"problems": e.problems}) from None
    pooled: dict = {}
    for r in sh.used:
        if r.barcode:
            pooled.setdefault(r.file, {})[r.sample] = [r.barcode, r.barcode2 or None]
    return {"samples": sh.samples(), "warnings": sh.warnings, "excluded": sh.excluded(),
            "pooled": pooled,
            "sheet": {"path": sh.path, "provenance": sh.provenance(), "aliases": sh.aliases(),
                      "excluded": sh.excluded(), "warnings": sh.warnings, "pooled": pooled}}


# ---- running a project -------------------------------------------------------
def _run_job(job_id: str, proj: dict):
    job = JOBS[job_id]

    def progress(frac, msg):
        job["progress"], job["message"] = round(float(frac), 3), str(msg)

    try:
        job["state"] = "running"
        pooled = (proj.get("sample_sheet") or {}).get("pooled") or {}
        if pooled:
            progress(0.02, "splitting the pooled run by barcode")
            _demultiplex(job_id, proj, pooled)
        res = project.run_project(proj, progress)
        if proj.get("_demux_notes"):
            res["notes"] = proj["_demux_notes"] + res.get("notes", [])
            res["demux"] = proj.get("_demux_stats")
        job["result"] = res
        job["state"] = "done"
        job["progress"] = 1.0
        job["message"] = "finished"
    except Exception as exc:
        job["state"] = "error"
        job["error"] = f"{type(exc).__name__}: {exc}"
        job["message"] = "failed"


def _demultiplex(job_id: str, proj: dict, pooled: dict) -> None:
    """Split each pooled FASTQ by inline barcode (demux.py) and point its samples at the parts."""
    from ..demux import demultiplex
    from ..io_utils import read_fasta
    ref = read_fasta(proj["reference"])[0][1].upper()
    opts = proj.get("demux") or {}
    left = (opts.get("left_anchor") or ref[:20]).upper()
    right = (opts.get("right_anchor") or ref[-20:]).upper()
    paths, stats, notes = {}, {}, []
    for k, (f, samples) in enumerate(pooled.items()):
        d = demultiplex(f, {n: (b[0], b[1]) for n, b in samples.items()}, left, right,
                        str(UPLOADS / f"demux_{job_id}_{k}"))
        paths.update(d["paths"])
        st = d["stats"]
        stats[Path(f).name] = st
        notes.append(f"{Path(f).name} split by barcode: {st['reads'] - st['unassigned']} of "
                     f"{st['reads']} reads assigned ("
                     + ", ".join(f"{n} {c}" for n, c in st["per_sample"].items())
                     + f"); {st['unassigned']} without a clear barcode were left out.")
    for smp in proj["samples"]:
        if smp["name"] in paths:
            smp["path"] = paths[smp["name"]]
    proj["_demux_notes"], proj["_demux_stats"] = notes, stats


@app.post("/api/run")
def run(payload: dict):
    proj = payload.get("project") or payload
    for key in ("reference", "samples", "region"):
        if not proj.get(key):
            raise HTTPException(400, f"the project needs '{key}'")
    if len(proj["samples"]) < 1:
        raise HTTPException(400, "add at least one sample")
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        JOBS[job_id] = {"state": "queued", "progress": 0.0, "message": "queued",
                        "error": None, "result": None, "project": proj}
    threading.Thread(target=_run_job, args=(job_id, proj), daemon=True).start()
    return {"job": job_id}


# ---- class-switch junctions (switch.py) ---------------------------------------
@app.post("/api/fasta-names")
def fasta_names(payload: dict):
    from ..io_utils import read_fasta
    path = Path(payload.get("path", "")).expanduser()
    if not path.is_file():
        raise HTTPException(400, f"no FASTA at {path}")
    try:
        recs = read_fasta(str(path))
    except Exception as e:                      # noqa: BLE001 - shown to the user
        raise HTTPException(400, f"could not read {path.name}: {e}") from None
    return {"records": [{"name": n, "length": len(s_)} for n, s_ in recs]}


def _run_switch(job_id: str, req: dict):
    from .. import switch
    job = JOBS[job_id]
    try:
        job["state"] = "running"
        samples = req["samples"]
        per, groups, rows = {}, {}, []
        for i, smp in enumerate(samples):
            job["progress"] = round(0.05 + 0.9 * i / len(samples), 3)
            job["message"] = f"finding junctions in {smp['name']}"
            r = switch.analyse(smp["path"], req["reference"], req.get("donor") or None)
            per[smp["name"]] = r["summary"] | {"reads": r["reads"]}
            groups[smp["name"]] = smp.get("group") or smp["name"]
            for j in r["junctions"]:
                d = switch.asdict(j)
                rows.append({"sample": smp["name"], "group": groups[smp["name"]], **d})
        comp = []
        cmp_ = req.get("compare") or {}
        if cmp_.get("case") and cmp_.get("control"):
            comp = switch.compare_groups(per, groups, cmp_["case"], cmp_["control"])
        job["result"] = {"samples": per, "groups": groups, "junctions": rows,
                         "comparison": comp, "reference": req["reference"],
                         "donor": req.get("donor")}
        job["state"], job["progress"], job["message"] = "done", 1.0, "finished"
    except Exception as exc:                    # noqa: BLE001 - reported to the page
        job["state"], job["error"], job["message"] = "error", f"{type(exc).__name__}: {exc}", "failed"


@app.post("/api/switch/run")
def switch_run(payload: dict):
    if not payload.get("reference"):
        raise HTTPException(400, "choose the switch-region reference FASTA")
    if not payload.get("samples"):
        raise HTTPException(400, "add at least one sample")
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        JOBS[job_id] = {"state": "queued", "progress": 0.0, "message": "queued",
                        "error": None, "result": None, "project": payload}
    threading.Thread(target=_run_switch, args=(job_id, payload), daemon=True).start()
    return {"job": job_id}


@app.get("/api/switch/result/{job_id}")
def switch_result(job_id: str):
    return JSONResponse(_finite(json.loads(json.dumps(_result(job_id), default=str))))


@app.get("/api/switch/export/{job_id}/{what}")
def switch_export(job_id: str, what: str):
    res = _result(job_id)
    if what == "junctions.csv":
        rows = res["junctions"]
    elif what == "summary.csv":
        rows = []
        for n, sm in res["samples"].items():
            flat = {"sample": n, "group": res["groups"].get(n, n)}
            for k, v in sm.items():
                if isinstance(v, dict):
                    flat.update({f"{k}:{kk}": vv for kk, vv in v.items()})
                else:
                    flat[k] = v
            rows.append(flat)
    elif what == "comparison.csv":
        rows = res["comparison"]
    else:
        raise HTTPException(404, f"cannot export {what}")
    cols = list(dict.fromkeys(k for r in rows for k in r)) or ["sample"]
    return Response(report._rows_csv(rows, cols) or ",".join(cols) + "\n", media_type="text/csv",
                    headers={"Content-Disposition": f"attachment; filename=switch_{what}"})


# ---- repair routes (repair.py) and power (power.py) ---------------------------
@app.post("/api/repair/{job_id}")
def repair_routes(job_id: str, payload: dict):
    from .. import repair
    res = _result(job_id)
    groups = {s["sample"]: s["group"] for s in res["samples"]}
    by_animal: dict = {}
    for v in res["variants"]:
        by_animal.setdefault(v["sample"], []).append((v["pos"], v["ref"], v["alt"]))
    per = []
    for name, g in groups.items():
        i = repair.indices(by_animal.get(name, []), res["reference"])
        i.pop("_counts", None)
        per.append({"sample": name, "group": g, **i})
    out: dict = {"per_animal": per, "comparison": [], "decomposition": None}
    case, control = payload.get("case"), payload.get("control")
    if case and control and case != control:
        out["comparison"] = repair.compare_indices(res, case, control)
    refs = {k: v for k, v in (payload.get("references") or {}).items() if v}
    if case and refs:
        spectra = {n: repair.group_spectra(res, g) for n, g in refs.items()}
        if any(not sum(sum(x.values()) for x in sp) for sp in spectra.values()):
            raise HTTPException(400, "a reference group has no mutations to build a spectrum from")
        d = repair.decompose_animals(repair.group_spectra(res, case), spectra,
                                     boot=int(payload.get("boot", 499)))
        sparse = {n: sum(sum(x.values()) for x in sp) for n, sp in spectra.items()}
        sparse = {n: k for n, k in sparse.items() if k < 20}
        out["decomposition"] = {"group": case, "weights": d["weights"], "ci": d["ci"],
                                "r2": d["r2"], "residual": d["residual"],
                                "warning": ("Reference spectra built from fewer than 20 mutations "
                                            "make the fit unreliable: " + ", ".join(
                                                f"{n} ({k})" for n, k in sparse.items()))
                                if sparse else ""}
    return JSONResponse(_finite(out))


@app.post("/api/power")
def power_table(payload: dict):
    from .. import power
    metric = payload.get("metric") or "frequency"
    pilot = None
    if payload.get("job"):
        try:
            pilot = power.pilot(_result(payload["job"]), payload.get("group"), metric)
        except ValueError as e:
            raise HTTPException(400, str(e)) from None
        rate, cv, disp, bp = pilot["rate"], pilot["cv_between"], pilot["dispersion"], pilot["bp"]
    else:
        try:
            rate, cv = float(payload["rate"]), float(payload["cv"])
            disp, bp = float(payload.get("dispersion") or 1), int(payload["bp"])
        except (KeyError, TypeError, ValueError):
            raise HTTPException(400, "give a finished project and group, or rate, cv and bp") from None
    if rate <= 0:
        raise HTTPException(400, "the pilot rate is zero; pick a group with mutations")
    fold = float(payload.get("fold") or 2)
    animals = [int(x) for x in payload.get("animals") or [3, 4, 5, 6, 8, 10]][:12]
    reads = [int(x) for x in payload.get("reads") or [125, 250, 500, 1000]][:8]
    sims = max(100, min(int(payload.get("sims") or 1000), 5000))
    rows = power.table(rate, fold, animals, reads, bp, cv, disp, sims=sims,
                       alpha=float(payload.get("alpha") or 0.05),
                       test=payload.get("test") or "welch")
    return JSONResponse(_finite({"pilot": pilot, "used": {"rate": rate, "cv_between": cv,
                                                          "dispersion": disp, "bp": bp, "fold": fold,
                                                          "metric": metric, "sims": sims},
                                 "rows": rows}))


@app.post("/api/screen/{job_id}")
def screen(job_id: str, payload: dict):
    from ..demux import screen_table
    res = _result(job_id)
    metric = payload.get("metric") or "unique_density"
    rows = screen_table(res, metric)
    if not rows:
        raise HTTPException(400, "no group comparisons: the analysis needs a control group "
                                 "and at least one other group")
    return JSONResponse(_finite({"metric": metric, "rows": rows}))


@app.get("/api/screen/{job_id}/screen.csv")
def screen_csv(job_id: str, metric: str = "unique_density"):
    from ..demux import screen_table
    rows = [{**{k: v for k, v in r.items() if k != "fold_ci"},
             "fold_ci_low": (r.get("fold_ci") or [None, None])[0],
             "fold_ci_high": (r.get("fold_ci") or [None, None])[1]}
            for r in screen_table(_result(job_id), metric)]
    cols = ["condition", "control", "animals", "fold", "fold_ci_low", "fold_ci_high", "p", "q"]
    return Response(report._rows_csv(rows, cols) or ",".join(cols) + "\n", media_type="text/csv",
                    headers={"Content-Disposition": "attachment; filename=screen.csv"})


@app.get("/api/job/{job_id}")
def job_status(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        raise HTTPException(404, "no such job")
    return {k: job[k] for k in ("state", "progress", "message", "error")}


def _result(job_id: str) -> dict:
    job = JOBS.get(job_id)
    if not job:
        raise HTTPException(404, "no such job")
    if job["state"] != "done":
        raise HTTPException(409, f"job is {job['state']}")
    return job["result"]


@app.get("/api/result/{job_id}")
def result(job_id: str):
    return JSONResponse(_finite(json.loads(json.dumps(project.public(_result(job_id)), default=str))))


@app.get("/api/export/{job_id}/{what}")
def export(job_id: str, what: str):
    res = _result(job_id)
    if what == "report.html":
        return HTMLResponse(report.html_report(res),
                            headers={"Content-Disposition": "attachment; filename=report.html"})
    if what == "bundle.zip":
        return Response(report.zip_bytes(res), media_type="application/zip",
                        headers={"Content-Disposition": "attachment; filename=shmtrace_results.zip"})
    if what == "tables.xlsx":
        try:
            blob = report.workbook(res)
        except ImportError:
            raise HTTPException(501, "Excel export needs openpyxl") from None
        return Response(blob, media_type=(
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
            headers={"Content-Disposition": "attachment; filename=shmtrace_tables.xlsx"})
    if what == "airr_rearrangements.tsv":
        from .. import airr
        rows = airr.project_rows(res)
        if not rows:
            raise HTTPException(404, "no analysed reads to export")
        return Response(airr.write_tsv(rows), media_type="text/tab-separated-values",
                        headers={"Content-Disposition": f"attachment; filename={what}"})
    if what.endswith(".csv"):
        tables = report.tables(res)
        if what not in tables:
            raise HTTPException(404, f"no table {what}")
        return Response(tables[what], media_type="text/csv",
                        headers={"Content-Disposition": f"attachment; filename={what}"})
    if what.rsplit(".", 1)[-1] in ("png", "svg", "pdf"):
        fmt = what.rsplit(".", 1)[-1]
        figs = report.figures(res, fmt=fmt)
        if what not in figs:
            raise HTTPException(404, f"no figure {what}")
        mime = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}[fmt]
        return Response(figs[what], media_type=mime,
                        headers={"Content-Disposition": f"attachment; filename={what}"})
    raise HTTPException(404, f"cannot export {what}")


@app.get("/api/alignment/{job_id}/{sample}")
def alignment(job_id: str, sample: str, start: int = 1, end: int | None = None,
              offset: int = 0, limit: int = 400, fates: str | None = None):
    """A window of one sample's read alignment, for the Reads tab.

    `fates` is a comma-separated filter ("kept,contaminant"); omit it for all.
    """
    res = _result(job_id)
    data = (res.get("_alignments") or {}).get(sample)
    if data is None:
        raise HTTPException(404, "this run did not record read alignments; "
                                 "tick 'Keep the read alignment' and run again")
    want = {f for f in fates.split(",") if f} if fates else None
    return browse.window(data, start=max(1, start), end=end, fates=want,
                         offset=max(0, offset), limit=min(2000, max(1, limit)))


@app.get("/api/alignment/{job_id}")
def alignment_index(job_id: str):
    """Which samples have a recorded alignment, and how big each is."""
    res = _result(job_id)
    al = res.get("_alignments") or {}
    return {"samples": [{"sample": k, "reads": len(v["reads"]), "counts": v["counts"],
                         "truncated": v["truncated"],
                         "reference_length": len(v["reference"])}
                        for k, v in al.items()],
            "fates": [{"name": n, "description": d} for n, d in browse.FATES]}


@app.get("/api/figure/{job_id}/{name}")
def figure(job_id: str, name: str):
    """A figure for display in the page (not a download)."""
    figs = report.figures(_result(job_id), fmt="png")
    if name not in figs:
        raise HTTPException(404, f"no figure {name}")
    return Response(figs[name], media_type="image/png",
                    headers={"Cache-Control": "no-store"})


@app.post("/api/save-to/{job_id}")
def save_to(job_id: str, payload: dict):
    """Desktop app: write every export into a folder the user picked."""
    out = Path(payload.get("path", "")).expanduser()
    written = report.write_all(_result(job_id), str(out))
    return {"folder": str(out), "files": [Path(p).name for p in written]}


# ---- the older endpoints, kept so existing scripts and the classic page work --
@app.post("/analyze")
async def analyze(
    fastq: UploadFile = File(...),
    reference: UploadFile = File(...),
    contaminant: UploadFile | None = File(None),
    null: list[UploadFile] | None = File(None),
    preset: str = Form("jh4_intron"),
    mode: str = Form("recommended"),
    platform: str = Form("ont"),
    downsample_to: str = Form(""),
    max_ee_rate: str = Form(""),
    max_pvalue: str = Form(""),
    germline_freq_threshold: str = Form(""),
    end_trim_bp: str = Form(""),
    analysis_window: str = Form(""),
):
    work = Path(tempfile.mkdtemp(prefix="shm_"))
    try:
        fq = await _save(fastq, work / "reads.fastq")
        ref = await _save(reference, work / "reference.fasta")
        contams = []
        if contaminant is not None:
            contams = [str(await _save(contaminant, work / "contaminant.fasta"))]

        def num(s, cast):
            return cast(s) if str(s).strip() else None

        null_paths = []
        for i, up in enumerate(null or []):
            if up.filename:
                null_paths.append(str(await _save(up, work / f"null_{i}.fastq")))

        overrides = dict(
            downsample_to=num(downsample_to, int),
            max_ee_rate=num(max_ee_rate, float),
            max_pvalue=num(max_pvalue, float),
            germline_freq_threshold=num(germline_freq_threshold, float),
            end_trim_bp=num(end_trim_bp, int),
        )
        if analysis_window.strip():
            a, b = analysis_window.replace(" ", "").split("-")
            overrides["analysis_window"] = (int(a), int(b))
        if hasattr(pipeline, "build_params"):
            p = pipeline.build_params(preset, platform=platform, mode=mode,
                                      **{k: v for k, v in overrides.items() if v is not None})
        else:
            overrides.update(pipeline.mode_overrides(mode))
            p = load_preset(preset, **overrides)
        null_rates = s_null = known = None
        if null_paths:
            null_rates, s_null, known = pipeline.fit_null(null_paths, str(fq), str(ref),
                                                          contams, p)
            if p.error_model == "geneious":
                null_rates = None
        result = pipeline.run(str(fq), str(ref), contams, p, null_rates, known)
        if s_null:
            result["summary"]["null_model"] = s_null
        variants = [v.__dict__ for v in result["variants"][:5000]]
        plots = {k: "data:image/png;base64," + base64.b64encode(v).decode()
                 for k, v in viz.render_all(result).items()}
        return JSONResponse({"summary": result["summary"],
                             "variants": variants, "plots": plots})
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)
    finally:
        shutil.rmtree(work, ignore_errors=True)


@app.post("/compare")
async def compare_groups(
    control: list[UploadFile] = File(...),
    case: list[UploadFile] = File(...),
    reference: UploadFile = File(...),
    contaminant: UploadFile | None = File(None),
    preset: str = Form("jh4_intron"),
    subregion: str = Form(""),
    mode: str = Form("recommended"),
):
    work = Path(tempfile.mkdtemp(prefix="shm_cmp_"))
    try:
        ctrl_paths, case_paths = [], []
        for i, up in enumerate(control):
            ctrl_paths.append(str(await _save(up, work / f"ctrl_{i}.fastq")))
        for i, up in enumerate(case):
            case_paths.append(str(await _save(up, work / f"case_{i}.fastq")))
        ref = str(await _save(reference, work / "reference.fasta"))
        contams = []
        if contaminant is not None:
            contams = [str(await _save(contaminant, work / "contaminant.fasta"))]

        sub = None
        if subregion.strip():
            a, b = subregion.replace(" ", "").split("-")
            sub = list(range(int(a), int(b) + 1))

        p = load_preset(preset, **pipeline.mode_overrides(mode))
        res = compare_mod.compare(ctrl_paths, case_paths, ref, contams, p, sub)
        plots = {k: "data:image/png;base64," + base64.b64encode(v).decode()
                 for k, v in viz.comparison_plots(res).items()}
        payload = {k: v for k, v in res.items() if not k.startswith("_")}
        payload["markdown"] = compare_mod.report_markdown(res)
        return JSONResponse({"report": payload, "plots": plots})
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)
    finally:
        shutil.rmtree(work, ignore_errors=True)


@app.get("/api/reference/{name}")
def builtin_reference(name: str):
    """Bundled reference sequences, so common regions need no file."""
    path = APP_DIR.parent / "references" / f"{name}.fasta"
    if not path.exists():
        raise HTTPException(404, f"no bundled reference {name}")
    return FileResponse(path, media_type="text/plain")


app.mount("/static", StaticFiles(directory=STATIC), name="static")
