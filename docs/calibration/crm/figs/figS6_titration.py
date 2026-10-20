"""Figure S6: real hypermutation diluted into AID-null reads (titration.py).

A  J_H4: mutation frequency against the frequency expected from the donor and diluent alone, every fraction
B  N114, the same (parental reference: ancestral_reference on)
C  unique mutations called against the share of donor reads, every donor: detection is all-or-nothing, and where it
   starts depends on how many mutations the donor carries
Writes out/titration_summary.json, the numbers the text quotes.
"""
import json
from collections import defaultdict

from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

DONOR_STYLE = {"WT2": ("#274b7e", "o"), "WT3": ("#6d9ad0", "s"), "WT4": ("#1f6f8b", "^"),
               "C1": ("#6a4c93", "o"), "D3": ("#b9a8d3", "s"), "D5": ("#a3621b", "^")}
ZERO_F = 0.003                        # where a frequency of 0 is drawn on the log axes (per 10^3 bp)
ZERO_U = 0.5                          # where 0 unique mutations are drawn


def summary():
    d = json.load(open(DATA / "titration.json"))
    by = defaultdict(lambda: defaultdict(list))
    for r in d["rows"]:
        by[(r["locus"], r["donor"])][r["fraction"]].append(r)
    out = {}
    for (locus, donor), fr in by.items():
        mean = lambda f, m: sum(r[m] for r in fr[f]) / len(fr[f])  # noqa: E731
        f1, f0 = mean(1.0, "mutation_frequency_x1e3"), mean(0.0, "mutation_frequency_x1e3")
        pts = []
        for f in sorted(fr):
            vals = [r["mutation_frequency_x1e3"] for r in fr[f]]
            pts.append({"fraction": f, "expected": f * f1 + (1 - f) * f0, "observed": vals,
                        "unique": [r["unique_mutations"] for r in fr[f]]})
        # the smallest fraction from which every replicate at every larger fraction made calls
        detect = next((p["fraction"] for i, p in enumerate(pts) if p["fraction"] > 0
                       and all(min(q["unique"]) > 0 for q in pts[i:])), None)
        ratios = [sum(p["observed"]) / len(p["observed"]) / p["expected"] for p in pts if p["fraction"] >= 0.25
                  and p["fraction"] < 1]
        out[f"{locus}|{donor}"] = {"locus": locus, "donor": donor, "points": pts, "detected_from": detect,
                                   "unique_at_full": mean(1.0, "unique_mutations"), "ratio_25_50": ratios}
    default = DATA / "titration_default_n114.json"
    if default.exists():
        dd = defaultdict(list)
        for r in json.load(open(default))["rows"]:
            if r["locus"] == "n114":
                dd[(r["donor"], r["fraction"])].append(r["mutation_frequency_x1e3"])
        out["n114_default"] = {f"{k[0]}|{k[1]}": sum(v) / len(v) for k, v in dd.items()}
    return out


def freq_panel(ax, s, locus, title):
    ax.plot([ZERO_F, 50], [ZERO_F, 50], color=SOFT, linestyle=(0, (3, 2)), linewidth=0.6, zorder=1)
    for v in (x for x in s.values() if isinstance(x, dict) and x.get("locus") == locus):
        col, mk = DONOR_STYLE[v["donor"]]
        xs, ys, lo, hi = [], [], [], []
        for p in v["points"]:
            if p["fraction"] == 0:
                continue
            m = sum(p["observed"]) / len(p["observed"])
            xs.append(p["expected"]); ys.append(max(m, ZERO_F))
            lo.append(max(min(p["observed"]), ZERO_F)); hi.append(max(max(p["observed"]), ZERO_F))
        ax.vlines(xs, lo, hi, color=col, linewidth=0.6, zorder=2)
        ax.plot(xs, ys, color=col, linewidth=0.8, marker=mk, markersize=3.2, markeredgewidth=0, zorder=3,
                label=v["donor"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.01, 40); ax.set_ylim(0.002, 40)
    ticks = [0.01, 0.1, 1, 10]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks]); ax.set_xticks([], minor=True)
    ax.set_yticks([ZERO_F] + ticks, ["0"] + [f"{t:g}" for t in ticks]); ax.set_yticks([], minor=True)
    ax.set_xlabel("Expected (per 10$^{3}$ bp)"); ax.set_ylabel("Observed mutation frequency\n(per 10$^{3}$ bp)")
    ax.set_title(title, fontsize=7, loc="left", pad=8)
    ax.legend(loc="upper left", handlelength=1.4, borderaxespad=0.2, labelspacing=0.3)


def main():
    setup()
    s = summary()
    json.dump(s, open(DATA / "titration_summary.json", "w"), indent=1)
    fig = plt.figure(figsize=(FULL, 62 * MM))
    gs = gridspec.GridSpec(1, 3, figure=fig, left=0.095, right=0.985, top=0.88, bottom=0.2, wspace=0.45)
    ax = fig.add_subplot(gs[0]); freq_panel(ax, s, "jh4", "J$_H$4"); letter(ax, "A", dx=-0.42)
    ax = fig.add_subplot(gs[1]); freq_panel(ax, s, "n114", "N114"); letter(ax, "B", dx=-0.42)

    # C: unique mutations called by share of donor reads
    ax = fig.add_subplot(gs[2])
    for v in (x for x in s.values() if isinstance(x, dict) and "points" in x):
        col, mk = DONOR_STYLE[v["donor"]]
        pts = [p for p in v["points"] if p["fraction"] > 0]
        xs = [100 * p["fraction"] for p in pts]
        ax.plot(xs, [max(sum(p["unique"]) / len(p["unique"]), ZERO_U) for p in pts], color=col, linewidth=0.8,
                marker=mk, markersize=3.2, markeredgewidth=0, linestyle="-" if v["locus"] == "jh4" else (0, (3, 1.5)),
                label=v["donor"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([1, 2, 5, 10, 25, 50, 100], ["1", "2", "5", "10", "25", "50", "100"]); ax.set_xticks([], minor=True)
    ax.set_yticks([ZERO_U, 1, 10, 100], ["0", "1", "10", "100"]); ax.set_yticks([], minor=True)
    ax.set_ylim(0.35, 3000)
    ax.set_xlabel("Donor reads in the library (%)"); ax.set_ylabel("Unique mutations called")
    ax.text(0.03, 0.97, "colours as in A and B\nsolid, J$_H$4; dashed, N114", transform=ax.transAxes, ha="left",
            va="top", color=MID, fontsize=6)
    letter(ax, "C", dx=-0.32)
    for a in fig.axes:
        a.tick_params(labelcolor=INK)
    save(fig, "FigureS6_titration")


if __name__ == "__main__":
    main()
