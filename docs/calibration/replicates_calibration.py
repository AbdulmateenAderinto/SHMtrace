"""Are the animal-level tests' p-values honest? Type I error and power by simulation.

    python docs/calibration/replicates_calibration.py --out replicates_calibration

Animals' true rates are log-normal around their group's rate (between-animal CV as in
the real J꜀4 wild types, about 1); each animal's count is negative binomial on top
(reads within an animal share clones). Groups of 3 and 5 animals, 262,913 bases each
(500 reads x 526 bp). With no group difference, a test at alpha = 0.05 should reject
5% of the time; with a 2-fold or 4-fold difference, the rejection rate is the power.
"""
import argparse
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from shmtrace import replicates

RATE = 3e-3            # mutations per base, a wild-type-like level
N_BASES = 262913


def one(args):
    n_animals, fold, cv, seed = args
    rng = np.random.default_rng(seed)
    sd = math.sqrt(math.log(1 + cv ** 2))

    def group(rate):
        r = rate * np.exp(rng.normal(-sd * sd / 2, sd, n_animals))
        lam = r * N_BASES
        k = rng.negative_binomial(5.0, 5.0 / (5.0 + lam))      # within-animal clumping
        return [(int(x), N_BASES) for x in k]
    c = replicates.compare(group(RATE), group(RATE * fold), boot=199, seed=seed)
    return {"n": n_animals, "fold": fold, "cv": cv,
            **{k: c.get(k) for k in ("p_welch_log", "p_negbinomial", "p_permutation",
                                     "p_pooled_reads")}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="replicates_calibration")
    ap.add_argument("--sims", type=int, default=300)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tasks = [(n, f, 1.0, 1000 * i + 7 * n + int(10 * f))
             for n in (3, 5) for f in (1.0, 2.0, 4.0) for i in range(a.sims)]
    with ProcessPoolExecutor(a.jobs) as ex:
        rows = list(ex.map(one, tasks, chunksize=8))
    json.dump(rows, open(os.path.join(a.out, "runs.json"), "w"))
    lines = ["| Animals per group | True fold | Welch on log10 | Negative binomial | Permutation | "
             "Pooled reads (pseudoreplication) |", "|---:|---:|---:|---:|---:|---:|"]
    for n in (3, 5):
        for f in (1.0, 2.0, 4.0):
            rs = [r for r in rows if r["n"] == n and r["fold"] == f]

            def rate(k, rs=rs):
                v = [r[k] for r in rs if r[k] is not None]
                return f"{sum(x < 0.05 for x in v) / len(v):.3f}" if v else "–"
            lines.append(f"| {n} | {f:g} | {rate('p_welch_log')} | {rate('p_negbinomial')} | "
                         f"{rate('p_permutation')} | {rate('p_pooled_reads')} |")
    md = ("Share of simulated experiments with p < 0.05 (fold 1: false positives, should be "
          "about 0.05; otherwise power).\n\n" + "\n".join(lines) + "\n")
    open(os.path.join(a.out, "report.md"), "w").write(md)
    print(md)


if __name__ == "__main__":
    main()
