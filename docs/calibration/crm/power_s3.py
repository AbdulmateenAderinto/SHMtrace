"""Power for Figure S3B: animals per group needed for a 2-fold change, with the pilot CV's own uncertainty.

    python crm/power_s3.py      -> crm/out/power_s3.json

The pilot is the three J_H4 wild-type mice at 500 reads (the J_H4 project). For each measure the between-animal
CV is estimated from those three mice, and so is very uncertain: with n = 3 the log-scale variance has 2 degrees
of freedom, and its 95% chi-square interval spans a wide range of CVs. Power is simulated with shmtrace.power
(Welch on log10 rates, alpha 0.05, two-sided, 1,000 experiments per cell) at the CV estimate and at both ends of
that interval, out to the group size at which the estimate reaches 80% power. With three mice the interval's
upper end is unbounded in practice (a CV above 10^5), so the figure also simulates twice the estimate as a
stated, arbitrary pessimistic case.
"""
import json
import math
from pathlib import Path

from scipy import stats

from shmtrace import power

CRM = Path(__file__).resolve().parent
CAL = CRM.parent
PILOT = CAL / "figures/jh4_shmtrace/project.json"
ANIMALS = [3, 5, 8, 10, 15, 20, 30, 40, 60]
MEASURES = {"density": "unique density", "clone_weighted": "clone-weighted frequency",
            "frequency": "read-weighted frequency"}
FOLD, READS, SIMS, ALPHA = 2.0, 500, 1000, 0.05


def cv_interval(cv: float, n: int) -> tuple[float, float]:
    """95% interval for a log-normal CV estimated from n animals (chi-square on the log-scale variance)."""
    s2 = math.log(1 + cv * cv)
    df = n - 1
    lo, hi = df * s2 / stats.chi2.ppf(0.975, df), df * s2 / stats.chi2.ppf(0.025, df)
    return math.sqrt(math.exp(lo) - 1), math.sqrt(math.exp(hi) - 1)


def main():
    res = json.load(open(PILOT))
    out = {"pilot": str(PILOT.relative_to(CAL)), "fold": FOLD, "reads": READS, "sims": SIMS, "alpha": ALPHA,
           "test": "Welch on log10 rates, two-sided", "animals": ANIMALS, "measures": {}}
    for metric, label in MEASURES.items():
        pl = power.pilot(res, "Wt", metric)
        lo, hi = cv_interval(pl["cv_between"], pl["animals"])
        row = {"label": label, "cv": pl["cv_between"], "cv_ci": [lo, hi], "rate": pl["rate"],
               "dispersion": pl["dispersion"], "pilot_animals": pl["animals"]}
        for name, cv in (("power", pl["cv_between"]), ("power_cv_low", lo), ("power_cv_double", 2 * pl["cv_between"])):
            row[name] = [power.simulate(pl["rate"], FOLD, n, READS, pl["bp"], cv, pl["dispersion"], sims=SIMS,
                                        alpha=ALPHA, seed=11) for n in ANIMALS]
        out["measures"][metric] = row
        print(label, f"CV {pl['cv_between']:.2f} [{lo:.2f}-{hi:.2f}]", [round(x, 2) for x in row["power"]])
    json.dump(out, open(CRM / "out/power_s3.json", "w"), indent=1)


if __name__ == "__main__":
    main()
