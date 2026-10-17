"""More replicates of the all-null cells, where the measure the review asks for needs them.

    python crm/fdr_null_reps.py [reps]     -> crm/out/fdr_null_reps.json

"The probability of any discovery in an all-null library" is the measure that decides whether the threshold
controls anything, and it is the one five runs cannot pin down: none out of five gives a 95% interval of
0 to 52%. These are the same cells as fdr_validation.py's all-null arm, run many more times so the interval
is worth quoting.
"""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fdr_validation as fv  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    cells = [(noise, excl, 9000 + rep)
             for noise in ("none", "one_noisy")
             for excl in (True, False)
             for rep in range(reps)]
    print(f"{len(cells)} runs", flush=True)
    rows, workdir = [], tempfile.mkdtemp(prefix="shmtrace_fdrnull_")
    for i, (noise, excl, seed) in enumerate(cells, 1):
        r = fv.run_cell("null", 3, noise, excl, 500, seed, workdir)
        rows.append(r)
        print(f"[{i}/{len(cells)}] {noise} excl={excl} -> {r['discoveries']} calls, "
              f"set aside {r['decoys_set_aside']}", flush=True)
        OUT.mkdir(exist_ok=True)
        json.dump({"runs": rows}, open(OUT / "fdr_null_reps.json", "w"), indent=1)
    print("written", OUT / "fdr_null_reps.json")


if __name__ == "__main__":
    main()
