"""Sensitivity of the control-based threshold to how the decoy tail is summarised.

    python crm/fdr_estimator.py [reps]      -> crm/out/fdr_estimator.json

The review asks for sensitivity to the handling of the zero-decoy tail rather than taking a target-decoy
convention on trust. SHMtrace averages the held-out controls' candidate counts at each threshold
(`control_decoy_estimate = "mean"`, an averaging strategy of the kind Keich and Noble recommend to reduce
the variance of a small-set target-decoy estimate). The conservative alternative already in the code is
"max": the threshold must satisfy the *noisiest* control rather than the average one.

This runs the same fresh libraries under both, so the cost in recall and the gain in false calls can be read
off directly rather than argued. The "mean" rows come from fdr_validation.py; only "max" is run here, on the
cells where the two can differ.
"""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fdr_validation as fv  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cells = [(signal, noise, 7000 + rep)
             for signal in ("null", "polyclonal")
             for noise in ("none", "one_noisy", "case_only")
             for rep in range(reps)]
    print(f"{len(cells)} runs (estimate=max)", flush=True)
    rows, workdir = [], tempfile.mkdtemp(prefix="shmtrace_fdrmax_")
    for i, (signal, noise, seed) in enumerate(cells, 1):
        r = fv.run_cell(signal, 3, noise, True, 500, seed, workdir, estimate="max")
        rows.append(r)
        print(f"[{i}/{len(cells)}] {signal} {noise} max -> {r['discoveries']} calls, "
              f"{r['false_calls']} false, FDP {r['fdp']}", flush=True)
        OUT.mkdir(exist_ok=True)
        json.dump({"runs": rows}, open(OUT / "fdr_estimator.json", "w"), indent=1)
    print("written", OUT / "fdr_estimator.json")


if __name__ == "__main__":
    main()
