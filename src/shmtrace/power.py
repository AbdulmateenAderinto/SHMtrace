"""How many animals, and how many reads per animal, does an SHM comparison need?

Power by simulation, with the animal as the unit (replicates.py). Three numbers
describe a design, and each can be taken from a pilot project:

  rate        the control group's mutations per base (or its unique-mutation density)
  cv_between  how much animals of one group differ: the coefficient of variation of
              their rates (log-normal). In the J꜀4 wild types it is about 1, so animal
              number usually matters more than read depth.
  dispersion  how clumped mutations are within an animal (reads share clones): the
              variance of mutations per read over their mean (index of dispersion).
              The count for an animal is negative binomial with that variance.

Each simulated experiment draws `animals` per group, a case group at `fold` times the
control rate, and applies the test (Welch on log10 rates by default, the manuscript's;
"permutation" is exact but cannot go below p = 0.1 with 3 vs 3).

    shmtrace-power --pilot project.json --group WT --fold 0.5 --animals 3,4,5,6,8,10 \\
        --reads 125,250,500,1000
    shmtrace-power --rate 0.003 --cv 1.0 --dispersion 6 --bp 529 --fold 2
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats

from .replicates import ZERO_EVENTS, permutation_p

METRICS = {
    # metric: (per-animal rate from a sample row, per-read clumping applies)
    "frequency": (lambda s: s["mutations_total"] / s["frequency_bases"]
                  if s.get("frequency_bases") else None, True),
    "density": (lambda s: s["unique_mutations"] / s["density_bases"]
                if s.get("density_bases") else None, False),
    "clone_weighted": (lambda s: (s.get("clone_weighted_frequency_x1e3") or 0) / 1e3
                       if s.get("clone_weighted_frequency_x1e3") is not None else None, False),
}


def pilot(res: dict, group: str, metric: str = "frequency") -> dict:
    """rate, cv_between, dispersion and analysed bp from a project result's group.

    metric: "frequency" (mutations per base, read-weighted), "density" (unique
    mutations per base sequenced) or "clone_weighted" (per V(D)J clone). Animals
    differ least in the last two when a clone dominates some of them.
    """
    rows = [s for s in res["samples"] if s["group"] == group]
    if not rows:
        raise ValueError(f"no group {group!r} in the result")
    get, clumped = METRICS[metric]
    rates = [r for r in (get(s) for s in rows) if r is not None]
    pos = [r for r in rates if r > 0]
    if len(pos) >= 2:
        sd = float(np.std(np.log(pos), ddof=1))
        cv = math.sqrt(math.exp(sd * sd) - 1)
    else:
        cv = 1.0
    disp = []
    for s in (rows if clumped else []):
        hist = (res.get("mutation_load") or {}).get(s["sample"]) or {}
        n = sum(hist.values())
        if n:
            m = sum(int(k) * v for k, v in hist.items()) / n
            v = sum(v * (int(k) - m) ** 2 for k, v in hist.items()) / max(n - 1, 1)
            if m > 0:
                disp.append(v / m)
    return {"group": group, "metric": metric, "animals": len(rows), "rate": float(np.mean(rates)) if rates else 0.0,
            "cv_between": cv, "dispersion": float(np.mean(disp)) if disp else 1.0,
            "bp": rows[0].get("analysed_bp"), "reads": int(np.mean([s["reads"] for s in rows]))}


def simulate(rate: float, fold: float, animals: int, reads: int, bp: int, cv_between: float,
             dispersion: float, sims: int = 2000, alpha: float = 0.05, test: str = "welch",
             seed: int = 1) -> float:
    """Share of simulated experiments in which `test` finds the difference at `alpha`."""
    rng = np.random.default_rng(seed)
    sd = math.sqrt(math.log(1 + cv_between ** 2))
    n_bases = reads * bp
    extra = max(dispersion - 1.0, 1e-9)     # NB: variance = mean x dispersion

    def counts(r):
        true = r * np.exp(rng.normal(-sd * sd / 2, sd, animals))
        lam = true * n_bases
        size = lam / extra
        return rng.negative_binomial(np.maximum(size, 1e-9), size / (size + lam))
    hits = 0
    for _ in range(sims):
        a, b = counts(rate), counts(rate * fold)
        # the zero-count rule of the test being powered: one event on the animal's
        # own exposure (replicates.ZERO_EVENTS)
        la = np.log10(np.maximum(a, ZERO_EVENTS) / n_bases)
        lb = np.log10(np.maximum(b, ZERO_EVENTS) / n_bases)
        if test == "permutation":
            p = permutation_p(list(lb), list(la))
        else:
            p = float(stats.ttest_ind(lb, la, equal_var=False).pvalue) if animals >= 2 else 1.0
        hits += p < alpha
    return hits / sims


def table(rate, fold, animals_list, reads_list, bp, cv_between, dispersion, **kw) -> list[dict]:
    return [{"animals_per_group": n, "reads_per_animal": r,
             "power": simulate(rate, fold, n, r, bp, cv_between, dispersion, **kw)}
            for n in animals_list for r in reads_list]


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="shmtrace-power", description=__doc__.split("\n\n")[0])
    ap.add_argument("--pilot", help="project.json from a pilot run")
    ap.add_argument("--group", help="the pilot group to take rate, variation and clumping from")
    ap.add_argument("--metric", default="frequency", choices=tuple(METRICS),
                    help="the measure to be compared")
    ap.add_argument("--rate", type=float)
    ap.add_argument("--cv", type=float, help="between-animal coefficient of variation")
    ap.add_argument("--dispersion", type=float, help="variance / mean of mutations per read")
    ap.add_argument("--bp", type=int, help="analysed bases per read")
    ap.add_argument("--fold", type=float, required=True, help="the case group's rate / control's")
    ap.add_argument("--animals", default="3,4,5,6,8,10")
    ap.add_argument("--reads", default="125,250,500,1000")
    ap.add_argument("--test", default="welch", choices=("welch", "permutation"))
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--sims", type=int, default=2000)
    a = ap.parse_args(argv)
    if a.pilot:
        pl = pilot(json.load(open(a.pilot)), a.group, a.metric)
        print("pilot:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in pl.items()})
        rate, cv, disp, bp = pl["rate"], pl["cv_between"], pl["dispersion"], pl["bp"]
    else:
        rate, cv, disp, bp = a.rate, a.cv, a.dispersion, a.bp
    if None in (rate, cv, disp, bp):
        raise SystemExit("give --pilot and --group, or all of --rate --cv --dispersion --bp")
    animals = [int(x) for x in a.animals.split(",")]
    reads = [int(x) for x in a.reads.split(",")]
    rows = table(rate, a.fold, animals, reads, bp, cv, disp, sims=a.sims, alpha=a.alpha,
                 test=a.test)
    print(f"\nPower to detect a {a.fold:g}-fold difference ({a.test}, alpha {a.alpha:g})\n")
    print("animals/group | " + " | ".join(f"{r} reads" for r in reads))
    for n in animals:
        print(f"{n:13d} | " + " | ".join(f"{x['power']:.2f}".rjust(len(f'{r} reads'))
                                        for x, r in zip([x for x in rows if x['animals_per_group'] == n],
                                                        reads, strict=True)))


if __name__ == "__main__":
    main()
