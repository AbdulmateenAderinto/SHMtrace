"""The set-up guide's server side: preset details, the example project and the read check."""
from pathlib import Path

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from shmtrace.webapp.main import app  # noqa: E402

client = TestClient(app)
ROOT = Path(__file__).resolve().parents[1]


def test_options_describe_every_preset_and_platform():
    o = client.get("/api/options").json()
    assert set(o["preset_details"]) == set(o["presets"])
    assert o["preset_details"]["jh4_intron"]["window"] == [101, 629]
    assert set(o["platform_help"]) == set(o["platforms"])


@pytest.mark.skipif(not (ROOT / "examples_compare" / "samples.csv").is_file(), reason="no example data")
def test_example_project_is_three_against_three_with_readable_files():
    ex = client.get("/api/example").json()
    assert ex["region"] == "jh4_intron"
    assert sum(s["control"] for s in ex["samples"]) == 3
    assert all(Path(s["path"]).is_file() for s in ex["samples"])
    check = client.post("/api/check-reads", json={"path": ex["samples"][0]["path"]}).json()
    assert check["ok"] and check["reads"] == 500 and check["median_length"] > 1000


def test_read_check_reports_a_missing_file_instead_of_failing(tmp_path):
    r = client.post("/api/check-reads", json={"path": str(tmp_path / "missing.fastq")}).json()
    assert r["ok"] is False and r["error"]
