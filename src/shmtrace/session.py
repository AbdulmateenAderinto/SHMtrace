"""Saved projects: a finished run, its inputs and its review marks, kept so it can be reopened.

A run otherwise lives only while the app is open. A saved project holds what the app
shows and exports (the result without the bulk per-read detail), the analysed reads the
Check calls strip draws, the project as it was run (so reads can still be extracted from
the input files, while those stay where they were) and the review marks.

On disk a project is a folder:

    manifest.json   format, versions, name, when it was saved, sample and call counts
    project.json    the project as run: reference, samples and their input paths, settings
    result.json.gz  the result, with the analysed reads under "_call_reads"
    review.json     the review marks, rewritten on its own each time a mark changes

and it travels as a `.shmtrace` file, a zip of that folder. Reopened, everything the app
shows and exports comes back, except what needs the per-read detail a run keeps only in
memory: the AIRR table and the diagnostic figures that draw individual reads.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import re
import sys
import time
import zipfile
from pathlib import Path

FORMAT = "shmtrace-project"
VERSION = 1
# kept beyond the public result: what Check calls draws, and the dropped reads when kept
KEEP_PRIVATE = ("_call_reads", "_alignments")


def home() -> Path:
    """Where saved projects live: SHMTRACE_HOME, else the platform's application-data folder."""
    if os.environ.get("SHMTRACE_HOME"):
        return Path(os.environ["SHMTRACE_HOME"]).expanduser() / "projects"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "SHMtrace" / "projects"
    if os.name == "nt" and os.environ.get("APPDATA"):
        return Path(os.environ["APPDATA"]) / "SHMtrace" / "projects"
    return Path.home() / ".shmtrace" / "projects"


def _jsonable(res: dict) -> dict:
    out = {k: v for k, v in res.items() if not k.startswith("_") or k in KEEP_PRIVATE}
    return json.loads(json.dumps(out, default=str))


def _manifest(res: dict, name: str) -> dict:
    from . import __version__
    review = res.get("review") or {}
    return {"format": FORMAT, "version": VERSION, "shmtrace": __version__, "name": name,
            "saved": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "samples": len(res.get("samples") or []), "calls": len(res.get("variants") or []),
            "reviewed": sum(1 for v in review.values() if v.get("mark")),
            "region": (res.get("params") or {}).get("region"),
            "platform": (res.get("params") or {}).get("platform")}


def folder_name(name: str, key: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name or "project").strip("_")[:60] or "project"
    return f"{time.strftime('%Y-%m-%d')}_{stem}_{key}"


def write(folder: Path, res: dict, proj: dict) -> Path:
    """Save a run as a project folder (created if need be); returns the folder."""
    folder.mkdir(parents=True, exist_ok=True)
    name = res.get("name") or proj.get("name") or "project"
    (folder / "project.json").write_text(json.dumps(proj, indent=1, default=str))
    with gzip.open(folder / "result.json.gz", "wt") as fh:
        json.dump(_jsonable({k: v for k, v in res.items() if k != "review"}), fh)
    write_review(folder, res)
    return folder


def write_review(folder: Path, res: dict) -> None:
    """Rewrite just the marks and the manifest: cheap enough to do on every mark."""
    name = res.get("name") or "project"
    tmp = folder / "review.json.tmp"
    tmp.write_text(json.dumps(res.get("review") or {}, indent=1))
    tmp.replace(folder / "review.json")
    (folder / "manifest.json").write_text(json.dumps(_manifest(res, name), indent=1))


def read(folder: Path) -> tuple[dict, dict]:
    """(result, project) from a project folder."""
    man = json.loads((folder / "manifest.json").read_text())
    _check(man)
    with gzip.open(folder / "result.json.gz", "rt") as fh:
        res = json.load(fh)
    review = folder / "review.json"
    res["review"] = json.loads(review.read_text()) if review.is_file() else {}
    proj = json.loads((folder / "project.json").read_text())
    return res, proj


def _check(man: dict) -> None:
    if man.get("format") != FORMAT:
        raise ValueError("this is not a saved SHMtrace project")
    if int(man.get("version", 0)) > VERSION:
        raise ValueError("this project was saved by a newer SHMtrace; update the app to open it")


def pack(folder: Path) -> bytes:
    """The project folder as one .shmtrace file (a zip)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name in ("manifest.json", "project.json", "result.json.gz", "review.json"):
            if (folder / name).is_file():
                z.write(folder / name, name)
    return buf.getvalue()


def unpack(blob: bytes, into: Path) -> Path:
    """A .shmtrace file back into a project folder under `into`; returns the folder."""
    try:
        z = zipfile.ZipFile(io.BytesIO(blob))
    except zipfile.BadZipFile:
        raise ValueError("this is not a saved SHMtrace project") from None
    names = set(z.namelist())
    if not {"manifest.json", "project.json", "result.json.gz"} <= names:
        raise ValueError("this is not a saved SHMtrace project")
    man = json.loads(z.read("manifest.json"))
    _check(man)
    folder = into / folder_name(man.get("name", "project"), f"opened{int(time.time())}")
    folder.mkdir(parents=True, exist_ok=True)
    for name in ("manifest.json", "project.json", "result.json.gz", "review.json"):
        if name in names:
            (folder / name).write_bytes(z.read(name))
    return folder


def recent(root: Path | None = None, limit: int = 20) -> list[dict]:
    """Saved projects, newest first, from their manifests."""
    root = root or home()
    out = []
    for f in root.glob("*/manifest.json") if root.is_dir() else []:
        try:
            man = json.loads(f.read_text())
            _check(man)
        except (ValueError, OSError):
            continue
        out.append({**man, "folder": f.parent.name})
    out.sort(key=lambda m: m.get("saved", ""), reverse=True)
    return out[:limit]
