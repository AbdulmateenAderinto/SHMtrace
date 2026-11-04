"""Sanger traces around a position (chromatogram.py), on a clone from the bundled Sanger example."""
from pathlib import Path

import pytest

from shmtrace import chromatogram

ROOT = Path(__file__).resolve().parents[1]
CLONE = ROOT / "examples_sanger/traces/berrumar_391861_chromats/7-2.1-7-user_added.ab1"
REF = str(ROOT / "src/shmtrace/references/mouse_jh4_intron.fasta")
pytestmark = pytest.mark.skipif(not CLONE.is_file(), reason="no Sanger example")


def test_a_reverse_clone_is_shown_along_the_reference_with_its_called_base():
    w = chromatogram.window(str(CLONE), REF, 104, flank=5)
    assert w["reverse"] is True                       # this clone reads against the reference
    cols = {c["pos"]: c for c in w["columns"]}
    assert cols[104]["base"] == "A"                    # the G104A the analysis called
    peaks = [c["peak"] for c in w["columns"]]
    assert peaks == sorted(peaks)                      # reversed, so peaks run left to right
    assert set(w["channels"]) == set("ACGT") and all(len(v) == w["length"] for v in w["channels"].values())
    # the call's own channel has the tallest signal at its peak
    at = cols[104]["peak"]
    assert max(w["channels"], key=lambda b: w["channels"][b][at]) == "A"


def test_a_position_the_clone_does_not_cover_is_refused():
    with pytest.raises(ValueError, match="does not cover"):
        chromatogram.window(str(CLONE), REF, 1225, flank=3)
