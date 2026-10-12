"""Statistics quoted in the manuscript (Ver17), from sanger_calls/corrected_sanger.json and the
ONT matrices in scripts/figure3_concordance.py. Run from this folder: python manuscript_stats.py

The spectrum statistics come from nanoshm.concordance (pooled Pearson r with Fisher CI, slope,
repeated-measures platform x type interaction, paired t per type with Sidak, G:C Ts/Tv).
Densities: Welch t on log10, group geometric means with 20,000-resample bootstrap CIs.
"""
import json
import random

import numpy as np
from scipy import stats

from nanoshm import concordance

cs = {k: v for k, v in json.load(open("sanger_calls/corrected_sanger.json")).items() if isinstance(v, dict)}
src = open("scripts/figure3_concordance.py").read()
ns = {"np": np}
exec('BASES=["A","T","C","G"]; IDX={b:i for i,b in enumerate(BASES)}\ndef m(r): return np.array(r,float)', ns)
exec(src[src.index('ont={"WT2"'):src.index('S=["WT2"')], ns)
ont, IDX = ns["ont"], ns["IDX"]
S = ["WT2", "WT3", "WT4"]


def counts_from_matrix(M):          # rows = to-base, columns = from-base
    return {f"{f}>{t}": M[IDX[t]][IDX[f]] for f in "ACGT" for t in "ACGT" if f != t}


def counts_from_calls(vs):
    out = {}
    for _p, f, t, _n in vs:
        out[f"{f}>{t}"] = out.get(f"{f}>{t}", 0) + 1
    return out


sp = concordance.spectrum({k: counts_from_calls(cs[f"Wt {k[2]}"]["variants"]) for k in S},
                          {k: counts_from_matrix(ont[k]) for k in S})
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
