"""Figure S11: SHMtrace against the Immcantation route on the same N114 reads (immcantation_compare.py).

A  mutation frequency over the V segment per library: IgBLAST + SHazaM with the IMGT germline, with the parental
   V as germline, and with bases below Q20 masked first; SHMtrace's calls over the same segment
B  separation: geometric mean of the pMX-Aicda libraries over that of the empty-vector libraries
A library with no SHMtrace call is drawn at one event on its own exposure (a triangle), as in Figure 4.
Writes out/immcantation_summary.json, the numbers the text quotes.
"""
import json
import math

from matplotlib import gridspec
from matplotlib import pyplot as plt
from style import SCRATCH  # noqa: E402
from style import CASE, CASEF, CTRL, DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

ROUTES = [("imgt", "Immcantation,\nIMGT germline"), ("parental", "Immcantation,\nparental germline"),
          ("parental_q20", "Immcantation,\nparental, Q<20\nmasked"), ("shmtrace", "SHMtrace")]
# what each route uses as reference, which bases it keeps and what it counts: the differences behind panel A
RULES = [("IMGT IGHV1-47", "all bases", "every mismatch"), ("parental V", "all bases", "every mismatch"),
         ("parental V", "Q ≥ 20 only", "every mismatch"), ("parental V", "control-calibrated", "called variants")]


def summary():
    d = json.load(open(DATA / "immcantation_compare.json"))
    lo, hi = d["v_region_reference"]
    vlen = hi - lo + 1
    out = {"v_call": d["v_call"], "v_identity": d["v_identity_parental_to_imgt"], "v_length": vlen, "routes": {}}
    for key, _lab in ROUTES:
        per = {}
        for r in d["rows"]:
            v = r[key]
            floor = 1 / (r["reads"] * vlen)               # one event on the library's own exposure
            per[r["library"]] = {"group": r["group"], "value": v, "plotted": max(v, floor), "zero": v == 0}
        geo = lambda g: math.exp(sum(math.log(x["plotted"]) for x in per.values() if x["group"] == g)  # noqa: E731
                                 / sum(1 for x in per.values() if x["group"] == g))
        out["routes"][key] = {"libraries": per, "aid_geo": geo("pMX-Aicda"), "empty_geo": geo("empty vector"),
                              "fold": geo("pMX-Aicda") / geo("empty vector")}
    return out


def main():
    setup()
    s = summary()
    json.dump(s, open(SCRATCH / "immcantation_summary.json", "w"), indent=1)
    fig = plt.figure(figsize=(FULL, 62 * MM))
    fig.set_size_inches(FULL, 70 * MM)
    gs = gridspec.GridSpec(1, 2, figure=fig, left=0.08, right=0.985, top=0.9, bottom=0.31, wspace=0.42,
                           width_ratios=[2.5, 1])
    zero_y = 3e-3                                       # the '0' strip, below an axis break

    ax = fig.add_subplot(gs[0])
    for i, (key, _lab) in enumerate(ROUTES):
        libs = s["routes"][key]["libraries"]
        for j, g in enumerate(("pMX-Aicda", "empty vector")):
            vals = [v for v in libs.values() if v["group"] == g]
            for k, v in enumerate(vals):
                x = i + (j - 0.5) * 0.36 + (k - 1) * 0.07
                face, edge = (CASEF, CASE) if j == 0 else ("white", CTRL)
                ax.scatter([x], [zero_y if v["zero"] else 1e3 * v["plotted"]], s=14, facecolor=face,
                           edgecolor=edge, linewidth=0.6, zorder=3)
            geo = s["routes"][key]["aid_geo" if j == 0 else "empty_geo"]
            ax.plot([i + (j - 0.5) * 0.36 - 0.13, i + (j - 0.5) * 0.36 + 0.13], [1e3 * geo] * 2,
                    color=CASE if j == 0 else CTRL, linewidth=1.2, solid_capstyle="butt", zorder=4)
    ax.axvline(2.5, color=SOFT, linewidth=0.5)
    ax.set_yscale("log"); ax.set_ylim(0.0018, 30)
    ax.axhspan(0.0018, 0.0048, color=SOFT, alpha=0.15, lw=0)
    for dy in (0.0052, 0.006):
        ax.plot([-0.012, 0.012], [dy] * 2, transform=ax.get_yaxis_transform(), color=INK, lw=0.6, clip_on=False)
    ax.set_yticks([zero_y, 0.01, 0.1, 1, 10], ["0", "0.01", "0.1", "1", "10"])
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xticks(range(len(ROUTES)), [lab for _k, lab in ROUTES])
    for i, rule in enumerate(RULES):
        ax.text(i, -0.3, "\n".join(rule), transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6, color=MID, linespacing=1.15)
    ax.text(-0.6, -0.3, "reference\nbases kept\ncounted", transform=ax.get_xaxis_transform(), ha="right",
            va="top", fontsize=6, color=MID, style="italic", linespacing=1.15)
    ax.set_xlim(-0.6, len(ROUTES) - 0.4)
    ax.set_ylabel("Mutation frequency over V\n(per 10$^{3}$ bp)")
    ax.scatter([], [], s=14, facecolor=CASEF, edgecolor=CASE, linewidth=0.6, label="pMX-$\\it{Aicda}$")
    ax.scatter([], [], s=14, facecolor="white", edgecolor=CTRL, linewidth=0.6, label="empty vector")
    ax.plot([], [], color=INK, lw=1.2, label="geometric mean (zeros at one event)")
    ax.legend(loc="lower left", handletextpad=0.3, borderaxespad=0.2, labelspacing=0.3, bbox_to_anchor=(0, 0.12))
    ax.set_title("Absolute frequencies on the same reads", fontsize=7)
    letter(ax, "A", dx=-0.13)

    ax = fig.add_subplot(gs[1])
    folds = [s["routes"][k]["fold"] for k, _l in ROUTES]
    cols = ["#b8c0c8", "#b8c0c8", "#b8c0c8", "#1f6f8b"]
    ax.bar(range(len(ROUTES)), folds, color=cols, width=0.7)
    for i, f in enumerate(folds):
        ax.text(i, f * 1.25, f"{f:.0f}×", ha="center", va="bottom", color=INK, fontsize=6)
    ax.set_yscale("log"); ax.set_ylim(1, 2000)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xticks(range(len(ROUTES)), ["IMGT", "parental", "Q<20\nmasked", "SHM-\ntrace"], rotation=0)
    ax.set_ylabel("pMX-$\\it{Aicda}$ / empty vector\n(geometric means)")
    ax.set_title("Ratio (secondary)", fontsize=7)
    ax.text(0.02, 0.98, "grey: Immcantation routes;\na near-zero denominator\ninflates a ratio (see A)",
            transform=ax.transAxes, ha="left", va="top", color=MID, fontsize=6)
    letter(ax, "B", dx=-0.42)
    for a in fig.axes:
        a.tick_params(labelcolor=INK)
    save(fig, "FigureS12_immcantation")


if __name__ == "__main__":
    main()
