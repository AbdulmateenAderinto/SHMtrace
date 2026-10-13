"""J_H4 and N114 manuscript projects re-run with the 5-mer targeting model (S5F index per sample),
clones and the single-read estimate, for the methods paper. Same reads and settings as
jh4_shmtrace.py / n114_shmtrace.py (delivered reads, 500 per library).

    python docs/calibration/crm/jh4_extras.py
Writes out/extras_<jh4|jh4_full|n114>/project.json.
"""
import json
import sys
from pathlib import Path

from shmtrace import project

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402

S5F = str(Path.home() / "Nanoshm/models/HH_S5F.csv")
OUT = Path(__file__).parent / "out"
RUNS = {"jh4": ("jh4", "jh4_intron", {}), "jh4_full": ("jh4", "jh4_intron", {"analysis_window": [90, 1230]}),
        "n114": ("n114", "n114_vregion", {})}
for name, (ds_key, region, ov) in RUNS.items():
    ds = DATASETS[ds_key]
    proj = {"name": name, "region": region, "reference": ds["reference"], "contaminants": ds.get("contaminants", []),
            "mode": "recommended", "depth": 500, "overrides": ov, "targeting_model": S5F,
            "samples": [{"name": n, "path": s["fastq"], "group": "control" if s["control"] else "case",
                         "control": s["control"]} for n, s in ds["samples"].items()]}
    res = project.run_project(proj)
    (OUT / f"extras_{name}").mkdir(parents=True, exist_ok=True)
    json.dump(project.public(res), open(OUT / f"extras_{name}/project.json", "w"), default=str)
    for s in res["samples"]:
        print(name, s["sample"], s["unique_mutations"], s.get("s5f_index"), s.get("s5f_p"), s.get("clones"),
              s.get("largest_clone_pct"), s.get("clone_weighted_frequency_x1e3"), s.get("mutation_frequency_x1e3"),
              s.get("single_read_estimate"), s.get("combined_density_x1e5"), flush=True)
