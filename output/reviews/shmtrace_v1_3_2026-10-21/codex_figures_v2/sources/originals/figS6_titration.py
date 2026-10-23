"""Figure S6: real hypermutation diluted into AID-null reads (titration.py).

A, B  mutation frequency observed against the frequency expected from the donor and diluent alone (J_H4; N114,
      which uses the parental reference and so runs with ancestral_reference)
C, D  unique mutations called against the share of donor reads, per donor

Zero is not a point on a logarithmic axis, so draws that produced no call are not placed on one. They sit as open
markers in a separate "no call" band below a break, the connecting line stops at the cliff rather than running
through it, and each of the three draws is drawn in its own right instead of as a range bar: with n = 3 a bar
reaching from zero to a real value is a scale artifact, not a measurement. The band is positioned from the data,
so it never overlaps a real value.

Writes out/titration_summary.json, the numbers the text quotes.
"""
import json
from collections import defaultdict

from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

DONOR_STYLE = {"WT2": ("#274b7e", "o"), "WT3": ("#6d9ad0", "s"), "WT4": ("#1f6f8b", "^"),
               "C1": ("#6a4c93", "o"), "D3": ("#b9a8d3", "s"), "D5": ("#a3621b", "^")}
LOCUS = {"jh4": "J$_H$4", "n114": "N114 V$_H$"}


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


def band_geometry(values):
    """Where the no-call strip sits, given every positive value a panel will draw.

    Placed below the smallest real value so the strip can never overlap one: the break sits a
    factor below the data, and the strip below that again. Returns (floor, strip_lo, strip_hi,
    break, marker y)."""
    v0 = min(values)
    brk = v0 / 2.2
    hi, lo = brk / 1.25, brk / 3.0
    return lo / 1.5, lo, hi, brk, (lo * hi) ** 0.5


def no_call_band(ax, lo, hi, brk, zero_at):
    ax.axhspan(lo, hi, color="#eef1f4", zorder=0, linewidth=0)
    ax.axhline(brk, color=SOFT, linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
    # right-aligned inside the strip: the no-call markers always sit at the low end of x
    ax.text(0.985, zero_at, "no call", fontsize=6, color=MID, va="center", ha="right",
            transform=ax.get_yaxis_transform(), zorder=4)


def series(ax, xs, means, draws, zero_at, col, mk, label):
    """One donor: the three draws behind, the mean in front, the line joining only consecutive
    points that were called, and open markers in the strip wherever a draw found nothing."""
    for x, dr in zip(xs, draws):
        for v in dr:
            called = v > 0
            ax.plot([x], [v if called else zero_at], marker=mk, markersize=2.4 if called else 3.0,
                    markeredgewidth=0 if called else 0.7, color=col, alpha=0.45 if called else 0.85,
                    markerfacecolor=col if called else "white", linestyle="none", zorder=2)
    run_x, run_y = [], []
    for x, m in zip(xs, means):
        if m > 0:
            run_x.append(x); run_y.append(m)
            continue
        if run_x:
            ax.plot(run_x, run_y, color=col, linewidth=0.9, zorder=3)
        run_x, run_y = [], []
    if run_x:
        ax.plot(run_x, run_y, color=col, linewidth=0.9, zorder=3)
    mx = [x for x, m in zip(xs, means) if m > 0]
    my = [m for m in means if m > 0]
    ax.plot(mx, my, color=col, linewidth=0, marker=mk, markersize=3.4, markeredgewidth=0,
            linestyle="none", zorder=4, label=label)


def panel(ax, s, locus, kind):
    """kind "freq": observed against expected frequency. kind "unique": calls against donor share."""
    pts = [v for v in s.values() if isinstance(v, dict) and v.get("locus") == locus]
    key = "observed" if kind == "freq" else "unique"
    rows = {v["donor"]: [p for p in v["points"] if p["fraction"] > 0] for v in pts}
    xs = {d: [(p["expected"] if kind == "freq" else 100 * p["fraction"]) for p in q] for d, q in rows.items()}
    means = {d: [sum(p[key]) / len(p[key]) for p in q] for d, q in rows.items()}
    draws = {d: [p[key] for p in q] for d, q in rows.items()}
    positive = [v for d in rows for v in means[d] + [x for dr in draws[d] for x in dr] if v > 0]
    floor, lo, hi, brk, zero_at = band_geometry(positive)
    no_call_band(ax, lo, hi, brk, zero_at)
    if kind == "freq":
        # the equality line is a reference for the points that were called, and spans only those
        ex = [x for d in rows for x, m in zip(xs[d], means[d]) if m > 0]
        ax.plot([min(ex), max(ex)], [min(ex), max(ex)], color=SOFT, linestyle=(0, (3, 2)),
                linewidth=0.6, zorder=1)
    for v in pts:
        col, mk = DONOR_STYLE[v["donor"]]
        series(ax, xs[v["donor"]], means[v["donor"]], draws[v["donor"]], zero_at, col, mk, v["donor"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(floor, max(positive) * 2.6)
    if kind == "freq":
        allx = [x for d in rows for x in xs[d]]
        ax.set_xlim(min(allx) / 1.6, max(allx) * 1.8)
        ticks = [0.01, 0.1, 1, 10]
        ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        yt = [t for t in ticks if t > brk]      # no tick may sit below the no-call strip
        ax.set_yticks(yt, [f"{t:g}" for t in yt])
        ax.set_xlabel("Expected (per 10$^{3}$ bp)")
        ax.set_ylabel("Observed mutation frequency\n(per 10$^{3}$ bp)")
    else:
        ax.set_xlim(0.8, 160)
        ax.set_xticks([1, 2, 5, 10, 25, 50, 100], ["1", "2", "5", "10", "25", "50", "100"])
        yt = [t for t in (1, 10, 100) if t > brk]
        ax.set_yticks(yt, [f"{t:g}" for t in yt])
        ax.set_xlabel("Donor reads in the library (%)")
        ax.set_ylabel("Unique mutations called")
    ax.set_xticks([], minor=True); ax.set_yticks([], minor=True)
    ax.set_title(LOCUS[locus], fontsize=7, loc="left", pad=8)
    ax.legend(loc="upper left", handlelength=1.2, borderaxespad=0.3, labelspacing=0.25,
              handletextpad=0.5, borderpad=0.4, framealpha=0.9, edgecolor="none")


def main():
    setup()
    s = summary()
    json.dump(s, open(DATA / "titration_summary.json", "w"), indent=1)
    fig = plt.figure(figsize=(FULL, 108 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.085, right=0.985, top=0.93, bottom=0.085,
                           wspace=0.28, hspace=0.45)
    for row, kind in enumerate(("freq", "unique")):
        for col, locus in enumerate(("jh4", "n114")):
            ax = fig.add_subplot(gs[row, col])
            panel(ax, s, locus, kind)
            letter(ax, "ABCD"[row * 2 + col], dx=-0.22)
    for a in fig.axes:
        a.tick_params(labelcolor=INK)
    save(fig, "FigureS6_titration")


if __name__ == "__main__":
    main()
