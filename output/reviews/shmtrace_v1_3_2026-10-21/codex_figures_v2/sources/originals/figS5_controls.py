"""Figure S5: how many AID-null controls, and which (controls_design.py).

A  recall of planted lineage mutations by number of controls (J_H4, the whole run_project procedure)
B  false calls per library by number of controls, in lineage-bearing and all-null libraries
C  one control carrying recurrent artifacts: false calls with the noisy-control rule off and on, beside the clean set
Writes out/controls_design.json, the numbers the text quotes.
"""
import json
import random
from collections import defaultdict

from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import CASE, CASEF, CTRL, DATA, FULL, INK, MID, MM, letter, save, setup

NS = (1, 2, 3, 4)
ARM_STYLE = {"matched (2)": ("#274b7e", "o"), "other experiment (2)": ("#a3621b", "s"),
             "other experiment (3)": ("#a3621b", "D"), "pooled (5)": ("#1f6f8b", "^")}


def summary():
    runs = json.load(open(DATA / "controls_design_jh4.json"))["runs"]
    cells = defaultdict(list)
    for r in runs:
        cells[(r["signal"], r["controls"], r["noise"], r["exclude_noisy"])].append(r)

    def cell(rs):
        calls = sum(r["discoveries"] for r in rs)
        planted = [r for r in rs if r["planted"]]
        return {"runs": len(rs), "false": [r["false_calls"] for r in rs],
                "recall": [r["true_calls"] / r["planted"] for r in planted],
                "fdp_pooled": sum(r["false_calls"] for r in rs) / calls if calls else 0.0,
                "any_call": sum(r["any_discovery"] for r in rs),
                "set_aside": sum(bool(r["decoys_set_aside"]) for r in rs)}

    jh4 = {f"{s}|{n}|{noise}|{'on' if ex else 'off'}": cell(rs) for (s, n, noise, ex), rs in cells.items()}
    su = json.load(open(DATA / "controls_design_su.json"))
    by = defaultdict(lambda: defaultdict(list))
    for r in su["rows"]:
        by[r["arm"]][r["case"]].append(r)
    arms = {}
    for arm, cases in by.items():
        rs = [r for v in cases.values() for r in v]
        tp, fp = sum(r["true_calls"] for r in rs), sum(r["false_calls"] for r in rs)
        per_case = {c: {"recall": sum(r["true_calls"] for r in v) / sum(r["planted"] for r in v),
                        "false_per_library": sum(r["false_calls"] for r in v) / len(v)} for c, v in cases.items()}
        arms[arm] = {"libraries": len(rs), "recall": tp / sum(r["planted"] for r in rs),
                     "false_per_library": fp / len(rs), "fdp_pooled": fp / (tp + fp) if tp + fp else 0.0,
                     "per_case": per_case}
    return {"jh4": jh4, "su": arms}


def strip(ax, x, vals, face, edge, rng, spread=0.13, d=3.2, marker="o"):
    xs = [x + rng.uniform(-spread, spread) for _ in vals]
    ax.scatter(xs, vals, s=d * d, facecolor=face, edgecolor=edge, linewidth=0.5, marker=marker, zorder=3)


def mean_bar(ax, x, vals, color, half=0.24):
    if vals:
        m = sum(vals) / len(vals)
        ax.plot([x - half, x + half], [m, m], color=color, linewidth=1.3, solid_capstyle="butt", zorder=4)


def main():
    setup()
    s = summary()
    json.dump(s, open(DATA / "controls_design.json", "w"), indent=1)
    J = s["jh4"]
    rng = random.Random(1)
    fig = plt.figure(figsize=(FULL, 64 * MM))
    gs = gridspec.GridSpec(1, 3, figure=fig, left=0.065, right=0.99, top=0.88, bottom=0.2, wspace=0.42,
                           width_ratios=[1, 1, 1.1])

    # A: recall by number of controls
    ax = fig.add_subplot(gs[0])
    for i, n in enumerate(NS):
        v = [100 * x for x in J[f"polyclonal|{n}|none|on"]["recall"]]
        strip(ax, i, v, CASEF, CASE, rng)
        mean_bar(ax, i, v, CASE)
    ax.set_xticks(range(len(NS)), [str(n) for n in NS]); ax.set_xlim(-0.6, len(NS) - 0.4)
    ax.set_ylim(50, 100); ax.set_yticks([50, 75, 100])
    ax.set_xlabel("AID-null controls"); ax.set_ylabel("Planted mutations recalled (%)")
    ax.text(0, 59, "BH, no\ndecoys", ha="center", va="top", color=MID, fontsize=6)
    letter(ax, "A", dx=-0.42)

    # B: false calls per library by number of controls; all-null libraries drawn open on the same axis
    ax = fig.add_subplot(gs[1])
    for i, n in enumerate(NS):
        v = J[f"polyclonal|{n}|none|on"]["false"]
        strip(ax, i - 0.17, v, CASEF, CASE, rng, spread=0.08)
        mean_bar(ax, i - 0.17, v, CASE, half=0.14)
        z = J[f"null|{n}|none|on"]["false"]
        strip(ax, i + 0.2, z, "white", CTRL, rng, spread=0.08)
        ax.text(i, 31, f"{100 * J[f'polyclonal|{n}|none|on']['fdp_pooled']:.1f}%", ha="center", va="bottom",
                color=INK, fontsize=6)
    ax.text(-0.55, 34, "false-discovery proportion:", color=MID, fontsize=6, va="bottom")
    ax.set_xticks(range(len(NS)), [str(n) for n in NS]); ax.set_xlim(-0.6, len(NS) - 0.4)
    ax.set_ylim(-1, 37); ax.set_yticks([0, 10, 20, 30])
    ax.set_xlabel("AID-null controls"); ax.set_ylabel("False calls per library")
    ax.scatter([], [], s=10, facecolor=CASEF, edgecolor=CASE, linewidth=0.5, label="lineage library")
    ax.scatter([], [], s=10, facecolor="white", edgecolor=CTRL, linewidth=0.5, label="all-null library")
    ax.legend(loc="upper left", handletextpad=0.1, borderaxespad=0.1, labelspacing=0.3, bbox_to_anchor=(-0.02, 0.83))
    letter(ax, "B", dx=-0.42)

    # C: one noisy control, rule off and on, beside the clean set (lineage libraries)
    ax = fig.add_subplot(gs[2])
    sets = [("clean", "none|on"), ("off", "one_noisy|off"), ("on", "one_noisy|on")]
    XC = (0, 1.4, 2.4)                    # clean, then the noisy pair
    for j, n in enumerate((3, 4)):
        for k, (_lab, key) in enumerate(sets):
            x = j * 4 + XC[k]
            v = J[f"polyclonal|{n}|{key}"]["false"]
            face, edge = (CASEF, CASE) if k == 0 else ("#e9c49d", "#a3621b")
            strip(ax, x, v, face, edge, rng, spread=0.12)
            mean_bar(ax, x, v, edge, half=0.22)
        ax.text(j * 4 + 1.2, 23.5, f"{n} controls", ha="center", va="bottom", color=INK, fontsize=6)
    ax.set_xticks([j * 4 + XC[k] for j in range(2) for k in range(3)], [lab for _ in range(2) for lab, _k in sets])
    ax.set_xlabel("One noisy control: rule off, on")
    ax.set_xlim(-0.6, 7.0); ax.set_ylim(-1, 27); ax.set_yticks([0, 10, 20])
    ax.set_ylabel("False calls per library")
    letter(ax, "C", dx=-0.42)

    save(fig, "FigureS5_controls")


if __name__ == "__main__":
    main()
