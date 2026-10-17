"""Depth series for Figure S5 with repeated draws, so the spread of one draw is visible.

    python crm/depth_draws.py {jh4|n114} [draws]      -> crm/out/depth_draws_<assay>.json

docs/calibration/depth_shmtrace.py runs each depth once, on SHMtrace's fixed seed. A single draw cannot show
how much of the curve is the depth and how much is which reads were drawn, so this repeats the whole analysis
at each depth with a different draw (downsample_seed), controls redrawn and recalibrated at the same depth
each time, exactly as a real run at that depth would be.

Zeros are kept as zeros here: a library with no calls has a density of 0, which has no logarithm and cannot
enter a geometric mean unmodified. The figure states how it is handled rather than letting it vanish.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from depth_shmtrace import ASSAYS, DEPTHS  # noqa: E402  (its assay definitions, unchanged)

from shmtrace import project  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"


def main():
    name = sys.argv[1]
    draws = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    a = ASSAYS[name]
    out = {"assay": name, "draws": draws, "depths": DEPTHS, "runs": []}
    for d in DEPTHS:
        for seed in range(draws):
            proj = {"name": f"{name} depth {d} draw {seed}", "region": a["region"],
                    "reference": str(a["reference"]), "contaminants": a["contaminants"],
                    "mode": "recommended", "depth": d, "overrides": {"downsample_seed": seed},
                    "samples": [{"name": n, "path": str(a["dir"] / f"{f}.fastq"), "group": g,
                                 "control": g == a["control"]} for n, (f, g) in a["libs"].items()]}
            res = project.run_project(proj)
            for s in res["samples"]:
                out["runs"].append({"depth": d, "draw": seed, "sample": s["sample"], "group": s["group"],
                                    "reads": s["reads"], "calls": s["unique_mutations"],
                                    "analysed_bp": s["analysed_bp"], "density": s["unique_density_x1e5"]})
            print(d, seed, {s["sample"]: s["unique_mutations"] for s in res["samples"]}, flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump(out, open(OUT / f"depth_draws_{name}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
