"""Saved projects (session.py): a run, its inputs and its marks survive the app closing."""
import gzip
import io
import json
import zipfile

import pytest

from shmtrace import session


def _res():
    return {"name": "JH4 Nov", "params": {"region": "jh4_intron", "platform": "ont"},
            "samples": [{"sample": "WT1", "group": "Wt"}, {"sample": "KO1", "group": "KO", "control": True}],
            "variants": [{"sample": "WT1", "pos": 12, "ref": "C", "alt": "T", "reads": 3}],
            "review": {"WT1|12|C>T": {"mark": "accepted", "note": "clean"}},
            "_call_reads": {"WT1": {"reads": [{"name": "r1", "fate": "kept", "strand": -1, "start": 1, "end": 30,
                                               "diffs": [[12, "T", 40]]}]}},
            "_details": {"WT1": {"_aligned": [object()]}}, "_suggestions": {"x": 1}}


def test_a_project_round_trips_with_its_marks_and_reads_but_not_the_in_memory_detail(tmp_path):
    proj = {"name": "JH4 Nov", "samples": [{"name": "WT1", "path": "/data/WT1.fastq"}]}
    folder = session.write(tmp_path / "p", _res(), proj)
    res, p = session.read(folder)
    assert p == proj
    assert res["review"] == {"WT1|12|C>T": {"mark": "accepted", "note": "clean"}}
    assert res["_call_reads"]["WT1"]["reads"][0]["strand"] == -1
    assert "_details" not in res and "_suggestions" not in res
    man = json.loads((folder / "manifest.json").read_text())
    assert man["format"] == session.FORMAT and man["reviewed"] == 1 and man["calls"] == 1
    with gzip.open(folder / "result.json.gz", "rt") as fh:
        assert "review" not in json.load(fh)        # the marks live in review.json alone


def test_marks_are_rewritten_on_their_own(tmp_path):
    res = _res()
    folder = session.write(tmp_path / "p", res, {})
    before = (folder / "result.json.gz").stat().st_mtime_ns
    res["review"]["WT1|12|C>T"] = {"mark": "rejected", "note": ""}
    session.write_review(folder, res)
    assert session.read(folder)[0]["review"]["WT1|12|C>T"]["mark"] == "rejected"
    assert (folder / "result.json.gz").stat().st_mtime_ns == before


def test_a_shmtrace_file_packs_and_unpacks_and_refuses_anything_else(tmp_path):
    folder = session.write(tmp_path / "p", _res(), {"name": "x"})
    blob = session.pack(folder)
    assert sorted(zipfile.ZipFile(io.BytesIO(blob)).namelist()) == [
        "manifest.json", "project.json", "result.json.gz", "review.json"]
    again = session.unpack(blob, tmp_path / "home")
    assert session.read(again)[0]["review"] == _res()["review"]
    with pytest.raises(ValueError, match="not a saved SHMtrace project"):
        session.unpack(b"not a zip", tmp_path / "home")
    newer = io.BytesIO()
    with zipfile.ZipFile(newer, "w") as z:
        for n in ("manifest.json", "project.json", "result.json.gz"):
            z.writestr(n, json.dumps({"format": session.FORMAT, "version": 99}) if n == "manifest.json" else "{}")
    with pytest.raises(ValueError, match="newer SHMtrace"):
        session.unpack(newer.getvalue(), tmp_path / "home")


def test_recent_lists_newest_first(tmp_path):
    a = session.write(tmp_path / "a", {**_res(), "name": "first"}, {})
    b = session.write(tmp_path / "b", {**_res(), "name": "second"}, {})
    for f, when in ((a, "2026-11-01T10:00:00"), (b, "2026-11-02T10:00:00")):
        m = json.loads((f / "manifest.json").read_text()); m["saved"] = when
        (f / "manifest.json").write_text(json.dumps(m))
    assert [p["name"] for p in session.recent(tmp_path)] == ["second", "first"]


fastapi = pytest.importorskip("fastapi")


def test_the_app_saves_a_run_reopens_it_and_keeps_marks_made_after(tmp_path):
    from fastapi.testclient import TestClient

    from shmtrace.webapp import main
    c = TestClient(main.app)
    main.JOBS["t-sess"] = {"state": "done", "result": _res(), "progress": 1, "message": "", "error": None,
                           "project": {"name": "JH4 Nov", "samples": []}}
    blob = c.get("/api/project-file/t-sess")
    assert blob.status_code == 200 and blob.headers["content-disposition"].endswith('JH4_Nov.shmtrace"')
    assert c.get("/api/job/t-sess").json()["saved"] is True
    c.post("/api/review/t-sess", json={"marks": {"WT1|12|C>T": {"mark": "unsure"}}})
    listed = c.get("/api/projects").json()["projects"]
    assert listed[0]["name"] == "JH4 Nov"
    job = c.post("/api/projects/open", json={"folder": listed[0]["folder"]}).json()["job"]
    assert c.get(f"/api/review/{job}").json()["WT1|12|C>T"]["mark"] == "unsure"
    assert c.get(f"/api/result/{job}").json()["name"] == "JH4 Nov"
    up = c.post("/api/projects/upload", files={"file": ("x.shmtrace", blob.content)}).json()["job"]
    assert c.get(f"/api/review/{up}").json()["WT1|12|C>T"]["mark"] == "accepted"   # as the file was saved
    assert c.post("/api/projects/open", json={"folder": "../etc"}).status_code == 404
    assert c.post("/api/projects/upload", files={"file": ("x.shmtrace", b"nope")}).status_code == 400
    for k in ("t-sess", job, up):
        main.JOBS.pop(k, None)


def test_a_sample_data_run_is_not_kept_until_asked(tmp_path, monkeypatch):
    """The bundled sample data is for running in front of people, so a run of it is not saved on
    its own; Save project file still keeps it, and any other run is saved as before."""
    from fastapi.testclient import TestClient

    from shmtrace import project
    from shmtrace.webapp import main
    monkeypatch.setattr(project, "run_project", lambda proj, progress=None: _res())
    c = TestClient(main.app)
    for job, example in (("t-demo", True), ("t-own", False)):
        main.JOBS[job] = {"state": "queued", "progress": 0, "message": "", "error": None, "result": None,
                          "project": {"name": "JH4 Nov", "samples": [], "example": example}}
        main._run_job(job, main.JOBS[job]["project"])
    demo, own = c.get("/api/job/t-demo").json(), c.get("/api/job/t-own").json()
    assert demo["state"] == "done" and demo["saved"] is False and demo["example"] is True
    assert own["saved"] is True and own["example"] is False
    assert c.get("/api/project-file/t-demo").status_code == 200
    assert c.get("/api/job/t-demo").json()["saved"] is True
    for k in ("t-demo", "t-own"):
        main.JOBS.pop(k, None)
