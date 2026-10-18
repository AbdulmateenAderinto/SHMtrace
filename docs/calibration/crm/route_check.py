"""Major issue 1: does the production route (project.run_project) reproduce the benchmark route?

    python crm/route_check.py        -> crm/out/route_check.json

Figure 4 and Table S2 are built by benchmark_real.py, which drives pipeline.analyse directly so that every
caller sees exactly the same reads. The product itself runs project.run_project. The referee asked which
route produced each Figure 4 library and how the two differ, so this runs the production route over the same
datasets and compares them library by library.
"""
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(CRM))
from shmtrace import project  # noqa: E402

from datasets import DATASETS  # noqa: E402

REGION = {"jh4": "jh4_intron", "n114": "n114_vregion", "su_lps": "generic_amplicon",
          "su_ch12": "generic_amplicon", "pim1": "generic_amplicon",
          "ot1": "generic_amplicon", "ot2": "generic_amplicon"}


def run(key):
    ds = DATASETS[key]
    proj = {"name": key, "region": REGION[key], "reference": ds["reference"],
            "contaminants": ds.get("contaminants", []), "mode": "recommended", "depth": ds["depth"],
            "overrides": dict(ds.get("overrides") or {}),
            "samples": [{"name": n, "path": s["fastq"], "group": "control" if s["control"] else "case",
                         "control": s["control"]} for n, s in ds["samples"].items()]}
    return {s["sample"]: s for s in project.public(project.run_project(proj))["samples"]}


def main(keys):
    out = json.loads((CRM / "out/route_check.json").read_text()) if (CRM / "out/route_check.json").exists() else {}
    for key in keys:
        f = CRM / f"out/benchmark_real_{key}.json"
        if not f.exists():
            print(f"{key}: no benchmark file"); continue
        bench = {r["sample"]: r for r in json.load(open(f))["rows"] if r["mode"] == "recommended"}
        prod = run(key)
        rows, same = {}, 0
        for s, b in bench.items():
            p = prod.get(s)
            if p is None:
                rows[s] = {"benchmark_calls": b["calls"], "production_calls": None}
                continue
            agree = b["calls"] == p["unique_mutations"]
            same += agree
            rows[s] = {"benchmark_calls": b["calls"], "production_calls": p["unique_mutations"],
                       "agree": agree, "benchmark_bases": int(b["bases"]),
                       "production_bases": int(p["frequency_bases"] or 0)}
        out[key] = {"libraries": len(rows), "calls_identical": same, "rows": rows}
        (CRM / "out/route_check.json").write_text(json.dumps(out, indent=1))
        print(f"{key}: {same}/{len(rows)} libraries identical")
        for s, r in rows.items():
            if not r.get("agree"):
                print(f"    DIFFERS {s}: benchmark {r['benchmark_calls']} vs production {r['production_calls']}")


if __name__ == "__main__":
    main(sys.argv[1:] or list(REGION))
