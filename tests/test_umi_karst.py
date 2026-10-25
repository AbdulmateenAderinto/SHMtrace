"""UMI consensus on real nanopore UMI reads (Karst et al. 2021), where the data are present.

The reads and references are downloaded separately (docs/calibration/UMI.md); elsewhere
the test is skipped.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path.home() / "Nanoshm/public_data/karst_umi"
ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(not (DATA / "ont_r10_zymo_rrna.fq.gz").exists(),
                                reason="Karst UMI test data not on this machine")


def test_consensus_removes_substitution_errors_on_real_umi_reads():
    subprocess.run([sys.executable, str(ROOT / "docs/calibration/umi_karst.py"), str(DATA)],
                   check=True, capture_output=True)
    s = json.load(open(ROOT / "docs/calibration/umi_karst/results.json"))["summary"]
    assert s["molecules_3plus_reads"] >= 90
    assert s["raw_read_error_pct"] > 5 and s["error_reduction_fold"] > 30
    deep = s["by_reads_per_molecule"]["20-+"]
    assert deep["substitution_error_pct"] < 0.02
    assert s["consensus_matching_group_operon_pct"] >= 99
    # dual 18-nt UMIs: the directional step (+2, the default) merges no different
    # species and leaves fewer error-prone split-off groups than grouping without it
    d = s["directional"]
    assert d["+2"]["groups_mixing_species_pct"] == 0
    assert d["+2"]["consensus_substitution_error_pct"] < d["off"]["consensus_substitution_error_pct"]
