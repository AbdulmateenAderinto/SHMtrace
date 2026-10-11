"""HTTP layer for the nanoshm desktop app (and the same thing served locally).

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

from .. import __version__, pipeline, project, report, viz
from .. import compare as compare_mod
from ..config import available_presets, load_preset

APP_DIR = Path(__file__).resolve().parent
STATIC = APP_DIR / "static"
MAX_UPLOAD_MB = 2000

app = FastAPI(title=f"nanoshm {__version__}")

UPLOADS = Path(tempfile.mkdtemp(prefix="nanoshm_uploads_"))
JOBS: dict[str, dict] = {}
_LOCK = threading.Lock()


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
        "modes": ["recommended", "geneious", "v0.1", "strand", "pbinom", "qcal+filter",
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
    exts = set(getattr(__import__("nanoshm.io_utils", fromlist=["x"]), "READ_EXTENSIONS",
                       (".fastq", ".fq", ".fasta", ".fa", ".bam", ".ab1")))
    files = []
    for p in sorted(root.iterdir()):
        if p.is_file() and (p.suffix.lower() in exts
                            or (p.suffix.lower() == ".gz" and p.with_suffix("").suffix.lower() in exts)):
            files.append({"name": p.name, "path": str(p), "bytes": p.stat().st_size})
    return {"folder": str(root), "files": files}


# ---- running a project -------------------------------------------------------
def _run_job(job_id: str, proj: dict):
    job = JOBS[job_id]

    def progress(frac, msg):
        job["progress"], job["message"] = round(float(frac), 3), str(msg)

    try:
        job["state"] = "running"
        res = project.run_project(proj, progress)
        job["result"] = res
        job["state"] = "done"
        job["progress"] = 1.0
        job["message"] = "finished"
    except Exception as exc:
        job["state"] = "error"
        job["error"] = f"{type(exc).__name__}: {exc}"
        job["message"] = "failed"


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
    return JSONResponse(json.loads(json.dumps(project.public(_result(job_id)), default=str)))


@app.get("/api/export/{job_id}/{what}")
def export(job_id: str, what: str):
    res = _result(job_id)
    if what == "report.html":
        return HTMLResponse(report.html_report(res),
                            headers={"Content-Disposition": "attachment; filename=report.html"})
    if what == "bundle.zip":
        return Response(report.zip_bytes(res), media_type="application/zip",
                        headers={"Content-Disposition": "attachment; filename=nanoshm_results.zip"})
    if what == "tables.xlsx":
        try:
            blob = report.workbook(res)
        except ImportError:
            raise HTTPException(501, "Excel export needs openpyxl") from None
        return Response(blob, media_type=(
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
            headers={"Content-Disposition": "attachment; filename=nanoshm_tables.xlsx"})
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
