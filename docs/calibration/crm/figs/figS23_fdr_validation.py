"""Figure S23: the whole calling procedure tested on fresh simulations, with the noisy-control rule on and off.

A  all-null libraries: the chance of any discovery at all, by noise condition, rule on and off
B  all-null libraries: mean false calls per library, the measure that hides A
C  libraries with planted mutations: the false-discovery proportion of each library, against the 5% target
D  what the rule costs and buys: control number and depth, under a heterogeneous null

Data: out/fdr_validation.json (fdr_validation.py), which runs project.run_project itself, so every number
here comes from the complete procedure: margin, threshold, and the noisy-decoy rule.
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from scipy import stats
from style import CTRL, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

LABEL = {"none": "unmodified\nbackground", "one_noisy": "one noisy\ncontrol",
         "case_only": "case-only\nartifact", "shared": "shared\nartifact"}
ORDER = ("none", "one_noisy", "case_only", "shared")
ON, OFF = "#1f6f8b", "#c08a3e"          # rule on / rule off


def clopper_pearson(k, n):
    if not n:
        return 0.0, 0.0, 0.0
    lo, hi = stats.beta.ppf(0.025, k, n - k + 1) if k else 0.0, stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return k / n, lo, hi


def read_runs(path):
    data = json.load(open(path))
    if data.get("revision_provenance", {}).get("complete") is False:
        raise RuntimeError(f"{path.name} is still being regenerated")
    return data["runs"]


def pick(runs, **kw):
    return [r for r in runs if all(r[k] == v for k, v in kw.items())]


def paired(ax, runs, signal, value, ylabel, title, as_rate=False, depth=500, controls=3):
    xs = np.arange(len(ORDER))
    for j, (excl, col, lab) in enumerate(((True, ON, "rule on"), (False, OFF, "rule off"))):
        ys, los, his = [], [], []
        for noise in ORDER:
            rs = pick(runs, signal=signal, noise=noise, exclude_noisy=excl, depth=depth, controls=controls)
            if as_rate:
                k = sum(1 for r in rs if r["any_discovery"])
                p, lo, hi = clopper_pearson(k, len(rs))
                ys.append(100 * p)
                los.append(100 * (p - lo))
                his.append(100 * (hi - p))
            else:
                v = [r[value] for r in rs]
                ys.append(np.mean(v) if v else 0)
                los.append(0)
                his.append(np.std(v) / max(1, len(v)) ** 0.5 if v else 0)
        ax.bar(xs + (j - 0.5) * 0.38, ys, 0.36, color=col, label=lab)
        ax.errorbar(xs + (j - 0.5) * 0.38, ys, yerr=[los, his], fmt="none", ecolor=INK, lw=0.6, capsize=1.2)
    ax.set_xticks(xs, [LABEL[n] for n in ORDER], fontsize=6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=6.6)


def main():
    setup()
    runs = read_runs(DATA / "fdr_validation.json")
    # the all-null cells are repeated many more times (fdr_null_reps.py): with five runs and no discovery
    # the interval on "any discovery" runs to 52%, which is not worth quoting
    extra = DATA / "fdr_null_reps.json"
    if extra.exists():
        runs = runs + read_runs(extra)
    est = DATA / "fdr_estimator.json"
    max_runs = read_runs(est) if est.exists() else []
    fig = plt.figure(figsize=(FULL, 134 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.105, right=0.985, top=0.90, bottom=0.13, wspace=0.43,
                           hspace=0.78)
    # A: the chance of any discovery in a library with nothing to find
    ax = fig.add_subplot(gs[0, 0])
    paired(ax, runs, "null", None, "Libraries with any call (%)",
           "All-null resamples: any discovery", as_rate=True)
    ax.axhline(5, color=SOFT, lw=0.6, ls=":")
    ax.text(3.45, 6, "5%", fontsize=6, color=MID, ha="right")
    ax.set_ylim(0, 112)
    ax.legend(fontsize=6, loc="upper left", handlelength=1.2)
    letter(ax, "A", dx=-0.22)
    # B: the same libraries, summarised the way the manuscript used to
    ax = fig.add_subplot(gs[0, 1])
    paired(ax, runs, "null", "false_calls", "False calls per library",
           "Mean false calls in the same resamples")
    ax.legend(fontsize=6, loc="upper left", handlelength=1.2)
    letter(ax, "B", dx=-0.22)
    # C: libraries with real mutations: each library's own false-discovery proportion
    ax = fig.add_subplot(gs[1, 0])
    xs = np.arange(len(ORDER))
    for j, (excl, col, _lab) in enumerate(((True, ON, "rule on"), (False, OFF, "rule off"))):
        for i, noise in enumerate(ORDER):
            rs = pick(runs, signal="polyclonal", noise=noise, exclude_noisy=excl, depth=500, controls=3)
            v = [100 * r["fdp"] for r in rs]
            x = i + (j - 0.5) * 0.38
            if v:
                ax.scatter(x + np.linspace(-0.08, 0.08, len(v)), v, s=5, color=col, linewidths=0, zorder=3)
                ax.plot([x - 0.13, x + 0.13], [np.mean(v)] * 2, color=INK, lw=1.1, zorder=4)
    ax.axhline(5, color=CTRL, lw=0.7, ls="--")
    ax.text(3.45, 6.5, "5% target", fontsize=6, color=CTRL, ha="right")
    ax.set_xticks(xs, [LABEL[n] for n in ORDER], fontsize=6)
    ax.set_ylim(0, None)
    ax.set_ylabel("False-discovery proportion (%)")
    ax.set_title("Mutation-bearing resamples:\nper-library false-discovery proportion", fontsize=6.6)
    ax.scatter([], [], s=5, color=ON, label="rule on")
    ax.scatter([], [], s=5, color=OFF, label="rule off")
    ax.legend(fontsize=6, loc="upper left", handlelength=1, handletextpad=0.3)
    letter(ax, "C", dx=-0.22)
    # D: how many controls, and at what depth, the rule can act
    ax = fig.add_subplot(gs[1, 1])
    cells = [("2 controls", {"controls": 2, "depth": 500}), ("3 controls", {"controls": 3, "depth": 500}),
             ("250 reads", {"controls": 3, "depth": 250}), ("1,000 reads", {"controls": 3, "depth": 1000})]
    xs = np.arange(len(cells))
    for j, (excl, col, lab) in enumerate(((True, ON, "rule on"), (False, OFF, "rule off"))):
        fired, calls = [], []
        for _label, kw in cells:
            rs = pick(runs, signal="null", noise="one_noisy", exclude_noisy=excl, **kw)
            fired.append(100 * np.mean([bool(r["decoys_set_aside"]) for r in rs]) if rs else 0)
            calls.append(np.mean([r["false_calls"] for r in rs]) if rs else 0)
        ax.bar(xs + (j - 0.5) * 0.38, fired, 0.36, color=col, label=lab)
    ax.set_xticks(xs, [c[0] for c in cells], fontsize=6)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Runs with a control\nexcluded (%)")
    ax.set_title("When the exclusion rule activates", fontsize=6.6)
    ax.legend(fontsize=6, loc="upper left", handlelength=1.2)
    letter(ax, "D", dx=-0.22)
    if max_runs:
        rows = []
        for noise in ("none", "one_noisy", "case_only"):
            mean_rs = pick(runs, signal="polyclonal", noise=noise, exclude_noisy=True, depth=500, controls=3)
            max_rs = [r for r in max_runs if r["signal"] == "polyclonal" and r["noise"] == noise]
            if mean_rs and max_rs:
                rows.append((noise,
                             float(np.mean([r["true_calls"] for r in mean_rs])),
                             float(np.mean([r["true_calls"] for r in max_rs])),
                             100 * float(np.mean([r["fdp"] for r in mean_rs])),
                             100 * float(np.mean([r["fdp"] for r in max_rs]))))
        json.dump([{"noise": n, "true_mean": a, "true_max": b, "fdp_mean_pct": c, "fdp_max_pct": d}
                   for n, a, b, c, d in rows], open(DATA / "figS15_estimator.json", "w"), indent=1)
        print("estimator (mean vs max):")
        for n, a, b, c, d in rows:
            print(f"  {n:10s} true {a:6.1f} -> {b:6.1f}   FDP {c:4.1f}% -> {d:4.1f}%")
    save(fig, "FigureS23")

    # the numbers the legend and the response quote
    out = {}
    for signal in ("null", "polyclonal"):
        for noise in ORDER:
            for excl in (True, False):
                rs = pick(runs, signal=signal, noise=noise, exclude_noisy=excl, depth=500, controls=3)
                if not rs:
                    continue
                fdp = [r["fdp"] for r in rs]
                tp = sum(r["true_calls"] for r in rs)
                fp = sum(r["false_calls"] for r in rs)
                out[f"{signal}/{noise}/{'on' if excl else 'off'}"] = {
                    "runs": len(rs),
                    "any_discovery_pct": round(100 * np.mean([r["any_discovery"] for r in rs]), 1),
                    "mean_false_calls": round(float(np.mean([r["false_calls"] for r in rs])), 2),
                    "mean_fdp_pct": round(100 * float(np.mean(fdp)), 1),
                    "max_fdp_pct": round(100 * float(np.max(fdp)), 1),
                    "pooled_fdr_pct": round(100 * fp / (tp + fp), 1) if tp + fp else 0.0,
                    "set_aside_pct": round(100 * np.mean([bool(r["decoys_set_aside"]) for r in rs]), 1),
                    "mean_true_calls": round(float(np.mean([r["true_calls"] for r in rs])), 1),
                }
    json.dump(out, open(DATA / "figS15_numbers.json", "w"), indent=1)
    for k, v in out.items():
        print(f"{k:34s} any {v['any_discovery_pct']:5.1f}%  meanFDP {v['mean_fdp_pct']:5.1f}%  "
              f"pooled {v['pooled_fdr_pct']:5.1f}%  true {v['mean_true_calls']:6.1f}  "
              f"set aside {v['set_aside_pct']:5.1f}%")


if __name__ == "__main__":
    main()
