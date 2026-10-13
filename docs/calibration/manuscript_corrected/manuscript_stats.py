"""Statistics quoted in the manuscript (Ver17), from sanger_calls/corrected_sanger.json and the
SHMtrace ONT re-analysis in ../figures/jh4_shmtrace/project.json. Run from this folder: python manuscript_stats.py

The spectrum statistics come from shmtrace.concordance (pooled Pearson r with Fisher CI, slope,
repeated-measures platform x type interaction, paired t per type with Sidak, G:C Ts/Tv).
Densities: Welch t on log10, group geometric means with 20,000-resample bootstrap CIs.
"""
import json
import random

import numpy as np
from scipy import stats

from shmtrace import concordance

cs = {k: v for k, v in json.load(open("sanger_calls/corrected_sanger.json")).items() if isinstance(v, dict)}
# ONT: the SHMtrace re-analysis (../figures/jh4_shmtrace/project.json, jh4_shmtrace.py)
JH4 = json.load(open("../figures/jh4_shmtrace/project.json"))
S = ["WT2", "WT3", "WT4"]


def counts_from_calls(vs):
    out = {}
    for _p, f, t, _n in vs:
        out[f"{f}>{t}"] = out.get(f"{f}>{t}", 0) + 1
    return out


sp = concordance.spectrum({k: counts_from_calls(cs[f"Wt {k[2]}"]["variants"]) for k in S},
                          {k: JH4["spectrum_sample"][k] for k in S})
pr = sp["paired"]
print(f"Pearson r {sp['pearson_r']:.5f} [{sp['pearson_ci'][0]:.2f}-{sp['pearson_ci'][1]:.2f}] "
      f"slope {sp['slope']:.2f}; mean|diff| {sp['mean_abs_diff_pp']:.2f} max {sp['max_abs_diff_pp']:.2f}")
print(f"interaction F({pr['interaction_df'][0]},{pr['interaction_df'][1]})={pr['interaction_F']:.2f} "
      f"p={pr['interaction_p']:.3f}")
best = min(pr["per_type"], key=lambda x: x["p_sidak"])
print("min Sidak", round(best["p_sidak"], 3), best["type"])
print("Ts/Tv adj p", [round(pr["ts_tv"][k]["p_sidak"], 3) for k in ("gc_transitions", "gc_transversions")])

dens = {k: 1e5 * max(v["unique"], 1) / v["bases"] for k, v in cs.items()}
wt = [dens[f"Wt {i}"] for i in (2, 3, 4)]
ko = [dens[f"Aicda−/− {i}"] for i in (1, 2, 3)]
rng = random.Random(1)


def geo(x):
    return 10 ** np.mean(np.log10(x))


def ci(x, n=20000):
    v = sorted(geo([rng.choice(x) for _ in x]) for _ in range(n))
    return v[int(.025 * n)], v[int(.975 * n) - 1]


print("densities", {k: round(v, 1) for k, v in dens.items()})
print(f"Wt geo {geo(wt):.0f} {[round(x) for x in ci(wt)]}  KO geo {geo(ko):.1f} "
      f"{[round(x, 1) for x in ci(ko)]}  fold {geo(wt) / geo(ko):.1f}  "
      f"P {stats.ttest_ind(np.log10(wt), np.log10(ko), equal_var=False).pvalue:.4f}")
loads = [x for k, v in cs.items() if k.startswith("Wt") for x in v["loads"]]
print("pie n", len(loads), [sum(1 for x in loads if a <= x <= b)
                            for a, b in ((0, 0), (1, 3), (4, 10), (11, 20), (21, 10 ** 6))])
print("Sanger unique pooled", sp["calls_a"], "ONT", sp["calls_b"], "| clones",
      {k: v["clones"] for k, v in cs.items()}, "| effective", {k: round(v["bases"] / 529, 1) for k, v in cs.items()})


# ONT densities (SHMtrace; a mouse with no calls at its one-event density, the plots' convention)
JS = {x["sample"]: x for x in JH4["samples"]}
od = {k: (JS[k]["unique_density_x1e5"] or JS[k]["one_event_density_x1e5"]) for k in JS}
owt = [od[k] for k in ("WT2", "WT3", "WT4")]
oko = [od[k] for k in ("KO1", "KO2", "KO3")]
print("ONT densities", {k: round(v, 2) for k, v in od.items()}, "calls",
      {k: JS[k]["unique_mutations"] for k in JS})
print(f"ONT Wt geo {geo(owt):.0f} {[round(x) for x in ci(owt)]}  KO geo {geo(oko):.2f} "
      f"{[round(x, 2) for x in ci(oko)]}  fold {geo(owt) / geo(oko):.0f}  "
      f"P {stats.ttest_ind(np.log10(owt), np.log10(oko), equal_var=False).pvalue:.4f}")
ab = {t["metric"]: t for t in JH4["animal_tests"]}["unique_density"]
print("SHMtrace NB test (unique density):", {k: ab[k] for k in ("fold", "p_negbinomial", "p_permutation")})
