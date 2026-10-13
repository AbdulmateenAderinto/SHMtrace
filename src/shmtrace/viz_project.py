"""Project figures: many samples in groups, drawn the way the paper draws them.

The styles follow the manuscript's Figures 2 to 4: a log-scale dot plot with a
geometric mean and a significance bracket, a from->to substitution matrix, the
twelve substitution classes with every sample shown, a positional map with one
thin line per sample under the group mean, and a mutations-per-read donut.
"""
from __future__ import annotations

import io
import math
import random
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# the manuscript's two-genotype palette, then distinct colours for further groups
MS_BLUE = "#1F6FB4"
MS_RED = "#C0392B"
GROUP_COLOURS = [MS_BLUE, MS_RED, "#3C7D4F", "#B07A1E", "#6A4C93", "#2A7F8F",
                 "#8C564B", "#C2185B", "#5D6D7E", "#7F8C1E"]
SUBS12 = [f"{a}>{b}" for a in "ATCG" for b in "ATCG" if a != b]


def _save(fig, fmt: str) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def _group_colour(res: dict) -> dict[str, str]:
    """Blue for the test groups, red for the control, in the manuscript's order."""
    order = res["group_order"]
    ctrl = {g["group"] for g in res.get("groups") or [] if g.get("control")}
    out, spare = {}, [c for c in GROUP_COLOURS if c not in (MS_BLUE, MS_RED)]
    for g in order:
        if g in ctrl:
            out[g] = MS_RED
        elif MS_BLUE not in out.values():
            out[g] = MS_BLUE
        else:
            out[g] = spare.pop(0) if spare else MS_BLUE
    return out


def _geo(vals: list[float], floor: float) -> tuple[float, float, float]:
    """Geometric mean and the geometric SD interval, zeros raised to the floor."""
    xs = [math.log10(max(v, floor)) for v in vals]
    m = sum(xs) / len(xs)
    sd = (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5 if len(xs) > 1 else 0.0
    return 10 ** m, 10 ** (m - sd), 10 ** (m + sd)


def _stat_for(res: dict, group: str, metric: str,
              log: bool = False) -> tuple[float | None, str | None]:
    """The Welch p against the control group, and the control's name.

    On a log axis the p that belongs with the plot is the one computed on log10
    densities, the same transform the eye is reading.
    """
    want = f"{metric}_p_welch_log_vs_" if log else f"{metric}_p_welch_vs_"
    for row in res.get("groups") or []:
        if row["group"] != group:
            continue
        ctrl = row.get("compared_with") or None
        if not ctrl:
            return None, None
        for k, v in row.items():
            if k.startswith(want):
                return v, ctrl
        for k, v in row.items():   # log p missing (percentages), fall back
            if k.startswith(f"{metric}_p_welch_vs_"):
                return v, ctrl
    return None, None


def group_dotplot(res: dict, metric: str, label: str, fmt: str = "png") -> bytes:
    """One dot per sample on a log axis, with the geometric mean and a p-value bracket.

    This is the manuscript's Figure 2A: the reader should be able to count the
    animals, see the spread, and read the group difference off the bracket.
    """
    order = res["group_order"]
    col = _group_colour(res)
    # a zero density is censored to the one-event scale, as the paper does, so
    # an empty sample still appears and still counts towards the group
    vals, censored = {}, set()
    for g in order:
        vals[g] = []
        for s in res["samples"]:
            if s["group"] != g or s[metric] is None:
                continue
            v, cen = s[metric], False
            if v == 0 and metric == "unique_density_x1e5" and s.get("one_event_density_x1e5"):
                v, cen = s["one_event_density_x1e5"], True
            vals[g].append(v)
            if cen:
                censored.add((g, len(vals[g]) - 1))
    allv = [v for vs in vals.values() for v in vs]
    if not allv:
        fig, ax = plt.subplots(figsize=(3.4, 3.4))
        ax.text(0.5, 0.5, "nothing to plot", ha="center", va="center", color="#777")
        ax.axis("off")
        return _save(fig, fmt)

    pos = [v for v in allv if v > 0]
    # a sample with no mutations cannot go on a log axis, so it sits one step
    # below the smallest real value and is drawn open, never silently dropped
    floor = (min(pos) / 3) if pos else 1e-3
    log = bool(pos) and (max(pos) / min(pos) > 8 or 0 in allv)

    fig, ax = plt.subplots(figsize=(max(2.8, 1.15 * len(order) + 1.1), 3.8))
    rng = random.Random(0)
    for i, g in enumerate(order):
        vs = vals[g]
        if not vs:
            continue
        if log:
            m, lo, hi = _geo(vs, floor)
        else:
            m = sum(vs) / len(vs)
            sd = ((sum((v - m) ** 2 for v in vs) / (len(vs) - 1)) ** 0.5) if len(vs) > 1 else 0.0
            lo, hi = m - sd, m + sd
        ax.hlines(m, i - 0.26, i + 0.26, color=col[g], lw=2.4, zorder=4)
        ax.vlines(i, lo, hi, color=col[g], lw=1.4, zorder=3)
        for cap in (lo, hi):
            ax.hlines(cap, i - 0.12, i + 0.12, color=col[g], lw=1.4, zorder=3)
        for n, v in enumerate(vs):
            y = max(v, floor) if log else v
            open_dot = log and ((g, n) in censored or v <= 0)
            ax.plot(i + rng.uniform(-0.13, 0.13), y, "o", ms=7, zorder=5,
                    color="white" if open_dot else col[g],
                    markeredgecolor=col[g], markeredgewidth=1.4)

    if log:
        ax.set_yscale("log")
        ax.grid(axis="y", which="major", color="#DDD", ls="--", lw=0.8)
        ax.set_axisbelow(True)

    # significance brackets against the control group
    top = max(max(vs) for vs in vals.values() if vs)
    step = 0
    for i, g in enumerate(order):
        p, ctrl = _stat_for(res, g, metric, log=log)
        if p is None or ctrl not in order:
            continue
        j = order.index(ctrl)
        a, b = sorted((i, j))
        y = (top * (1.9 ** (1 + step))) if log else top * (1.12 + 0.12 * step)
        ax.plot([a, a, b, b], [y / 1.25 if log else y - top * 0.03, y, y,
                               y / 1.25 if log else y - top * 0.03],
                color="#333", lw=1.0, clip_on=False)
        gm_i = _geo(vals[g], floor)[0] if log else sum(vals[g]) / len(vals[g])
        gm_j = _geo(vals[ctrl], floor)[0] if log else sum(vals[ctrl]) / len(vals[ctrl])
        fold = (gm_i / gm_j) if gm_j else None
        ptxt = f"{p:.4f}" if p >= 1e-4 else f"{p:.1e}"
        txt = f"P = {ptxt}" + (f" (~{fold:.0f}×)" if fold and fold >= 1.5 else "")
        ax.text((a + b) / 2, y * 1.08 if log else y + top * 0.02, txt, ha="center",
                va="bottom", fontsize=9.5, color="#333")
        step += 1

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=30 if len(order) > 4 else 0,
                       ha="right" if len(order) > 4 else "center")
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylabel(label)
    ax.spines[["top", "right"]].set_visible(False)
    return _save(fig, fmt)


MATRIX_BASES = "ATCG"   # the order the manuscript's Figure 3A uses


def group_matrix(res: dict, fmt: str = "png") -> bytes:
    """From (rows) to (columns) substitution counts per group, the Figure 3A matrix.

    The right-hand strip is each starting base's share of the group's mutations,
    which is where AID's G:C preference shows up at a glance.
    """
    order = res["group_order"]
    bb = MATRIX_BASES
    fig, axes = plt.subplots(1, len(order), figsize=(3.6 * len(order), 3.3), squeeze=False)
    for ax, g in zip(axes[0], order, strict=True):
        c = res["spectrum_group"].get(g, {})
        tot = sum(c.values())
        grid = [[(c.get(f"{a}>{b}", 0) if a != b else None) for b in bb] for a in bb]
        pct = [[0.0 if v is None else 100 * v / (tot or 1) for v in row] for row in grid]
        rowpct = [sum(r) for r in pct]
        vmax = max(20.0, max(max(r) for r in pct))
        # the matrix, then a one-column strip of row totals, drawn beside it
        ax.imshow([r + [0.0] for r in pct], cmap="Blues", vmin=0, vmax=vmax)
        ax.imshow([[v] for v in rowpct], cmap="Blues", vmin=0, vmax=max(vmax, max(rowpct)),
                  extent=(4.6, 5.4, 3.5, -0.5))
        for i, a in enumerate(bb):
            for j, b in enumerate(bb):
                if a == b:
                    ax.text(j, i, "\u2014", ha="center", va="center", color="#999")
                    continue
                dark = pct[i][j] > 0.55 * vmax
                ax.text(j, i - 0.12, str(grid[i][j]), ha="center", va="center", fontsize=10,
                        fontweight="bold", color="white" if dark else "#222")
                ax.text(j, i + 0.26, f"{pct[i][j]:.1f}%", ha="center", va="center",
                        fontsize=7.5, color="white" if dark else "#555")
            ax.text(5.0, i, f"{rowpct[i]:.1f}", ha="center", va="center", fontsize=10,
                    fontweight="bold",
                    color="white" if rowpct[i] > 0.55 * max(vmax, max(rowpct)) else "#222")
        ax.set_xticks(list(range(4)) + [5.0], list(bb) + ["%"])
        ax.set_yticks(range(4), list(bb))
        ax.set_xlim(-0.5, 5.6)
        ax.set_ylim(3.5, -0.5)
        ax.set_xlabel("To")
        ax.set_ylabel("From")
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_title(f"{g}  (n = {tot})", pad=26, fontsize=11, fontweight="bold")
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.text(0.5, -0.02, "cell: unique mutations, and % of the group's total; "
                         "right-hand strip: % starting from that base",
             ha="center", fontsize=9, color="#555")
    return _save(fig, fmt)


def group_spectrum(res: dict, fmt: str = "png") -> bytes:
    """The twelve substitution classes: group mean, SD, and every sample as a point."""
    # a group with a handful of calls has no spectrum: 1 of 4 is "25%", which
    # would tower over a real 10%. Those groups are named, not drawn.
    MIN_CALLS = 20
    counts = {g: sum((res["spectrum_group"].get(g) or {}).values())
              for g in res["group_order"]}
    order = [g for g in res["group_order"] if counts[g] >= MIN_CALLS]
    dropped = [g for g in res["group_order"] if counts[g] < MIN_CALLS]
    if not order:                       # nothing passes: show everything instead
        order, dropped = res["group_order"], []
    col = _group_colour(res)
    per = {}   # group -> class -> list of per-sample percentages
    for s in res["samples"]:
        if s["group"] not in order:
            continue
        c = res["spectrum_sample"].get(s["sample"], {})
        tot = sum(c.values())
        if not tot:
            continue
        per.setdefault(s["group"], {k: [] for k in SUBS12})
        for k in SUBS12:
            per[s["group"]][k].append(100 * c.get(k, 0) / tot)

    fig, ax = plt.subplots(figsize=(9.2, 3.6))
    w = 0.8 / max(1, len(order))
    rng = random.Random(1)
    for j, g in enumerate(order):
        xs = [i + j * w - 0.4 + w / 2 for i in range(12)]
        vals = per.get(g)
        if vals:
            means = [sum(vals[k]) / len(vals[k]) for k in SUBS12]
            sds = [((sum((v - m) ** 2 for v in vals[k]) / (len(vals[k]) - 1)) ** 0.5)
                   if len(vals[k]) > 1 else 0.0 for k, m in zip(SUBS12, means, strict=True)]
        else:   # a group with no called mutations still gets its legend entry
            c = res["spectrum_group"].get(g, {})
            tot = sum(c.values()) or 1
            means = [100 * c.get(k, 0) / tot for k in SUBS12]
            sds = [0.0] * 12
        ax.bar(xs, means, width=w * 0.92, color=col[g], label=f"{g}", zorder=2)
        ax.errorbar(xs, means, yerr=[[0] * 12, sds], fmt="none", ecolor=col[g],
                    elinewidth=1.1, capsize=2.5, zorder=3)
        if vals:
            for x, k in zip(xs, SUBS12, strict=True):
                for v in vals[k]:
                    ax.plot(x + rng.uniform(-w * 0.18, w * 0.18), v, "o", ms=3.4,
                            mfc="white", mec=col[g], mew=0.9, zorder=4)
    ax.set_xticks(range(12), [s.replace(">", "→") for s in SUBS12])
    ax.set_ylabel("% of unique mutations")
    ax.grid(axis="y", color="#DDD", ls="--", lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, ncol=min(4, len(order)),
              loc="upper center", bbox_to_anchor=(0.5, 1.16))
    ax.spines[["top", "right"]].set_visible(False)
    if dropped:
        names = ", ".join(f"{g} (n = {counts[g]})" for g in dropped)
        fig.text(0.5, -0.17, f"Not shown, too few calls for a spectrum: {names}",
                 ha="center", fontsize=9, color="#777")
    return _save(fig, fmt)


def group_positional(res: dict, fmt: str = "png") -> bytes:
    """Where the mutations sit: the group mean in bold over one thin line per sample."""
    order = res["group_order"]
    col = _group_colour(res)
    lo, hi = res["params"]["analysis_window"]
    by_sample = {}
    for r in res.get("positional_sample") or []:
        by_sample.setdefault(r["sample"], []).append((r["pos"], r["frequency_pct"]))
    gmap = {s["sample"]: s["group"] for s in res["samples"]}

    fig, axes = plt.subplots(len(order), 1, figsize=(9.4, 1.9 * len(order) + 0.6),
                             sharex=True, squeeze=False)
    for ax, g in zip(axes[:, 0], order, strict=True):
        rows = sorted((r["pos"], r["frequency_pct"]) for r in res["positional"]
                      if r["group"] == g)
        mine = {n: sorted(pts) for n, pts in by_sample.items() if gmap.get(n) == g}
        ymax = max([f for _, f in rows] + [f for pts in mine.values() for _, f in pts]
                   or [0]) or 1.0
        # with only a handful of mutated positions a line would invent a trend
        # between them, so a sparse group is drawn as stems
        sparse = len(rows) < 20
        for pts in mine.values():
            if sparse:
                ax.plot([p for p, _ in pts], [f for _, f in pts], "o", ms=2.6,
                        color=col[g], alpha=0.5, zorder=2)
            else:
                ax.plot([p for p, _ in pts], [f for _, f in pts], color=col[g], lw=0.7,
                        alpha=0.4, zorder=2)
        if sparse:
            ax.vlines([p for p, _ in rows], 0, [f for _, f in rows], color=col[g],
                      lw=1.8, zorder=3)
            ax.plot([p for p, _ in rows], [f for _, f in rows], "o", ms=4,
                    color=col[g], zorder=4)
        else:
            ax.plot([p for p, _ in rows], [f for _, f in rows], color=col[g], lw=1.8,
                    zorder=3)
        ax.plot(res["hotspots"], [-0.07 * ymax] * len(res["hotspots"]), "|",
                color="#AAA", ms=4, zorder=1)
        ax.set_ylim(-0.13 * ymax, 1.1 * ymax)
        # the group's name above its panel, so it never meets the shared y-axis label
        ax.set_title(g, loc="left", fontsize=10, fontweight="bold", pad=3)
        ax.set_xlim(lo - 2, hi + 2)
        ax.grid(axis="y", color="#EEE", lw=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1, 0].set_xlabel("Position in the amplicon (bp).  "
                           "Grey ticks: WRC/GYW hotspots.  "
                           "Faint traces: one sample each.  "
                           "Each panel has its own y scale")
    fig.supylabel("% of reads mutated at that position", fontsize=10)
    return _save(fig, fmt)


LOAD_CLASSES = [("0", 0, 0), ("1–3", 1, 3), ("4–10", 4, 10),
                ("11–20", 11, 20), ("20+", 21, 10 ** 9)]
LOAD_SHADES = ["#FFFFFF", "#D6E6F4", "#A8CBE8", "#5B9BD5", "#1F4E79"]


def group_load(res: dict, fmt: str = "png") -> bytes:
    """Mutations per read as a donut per group, the manuscript's Figure 3E."""
    order = res["group_order"]
    gmap = {s["sample"]: s["group"] for s in res["samples"]}
    pooled = {g: Counter() for g in order}
    for s, hist in res["mutation_load"].items():
        for k, n in hist.items():
            pooled[gmap[s]][int(k)] += n

    fig, axes = plt.subplots(1, len(order), figsize=(2.9 * len(order) + 1.4, 3.4),
                             squeeze=False)
    for ax, g in zip(axes[0], order, strict=True):
        counts = [sum(n for k, n in pooled[g].items() if a <= k <= b)
                  for _, a, b in LOAD_CLASSES]
        tot = sum(counts)
        if not tot:
            ax.text(0.5, 0.5, "no reads", ha="center", va="center", color="#777")
            ax.axis("off")
            continue
        wedges, _ = ax.pie(counts, colors=LOAD_SHADES, startangle=90, counterclock=False,
                           wedgeprops={"width": 0.42, "edgecolor": "#333", "linewidth": 0.8})
        for n, (wedge, c) in enumerate(zip(wedges, counts, strict=True)):
            pct = 100 * c / tot
            if pct < 3:
                continue
            ang = math.radians((wedge.theta1 + wedge.theta2) / 2)
            ax.text(0.79 * math.cos(ang), 0.79 * math.sin(ang), f"{pct:.0f}%",
                    ha="center", va="center", fontsize=9.5, fontweight="bold",
                    color="white" if n >= 3 else "#222")
        ax.text(0, 0, f"n = {tot:,}", ha="center", va="center", fontsize=10)
        ax.set_title(g, fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=s, ec="#333", lw=0.6)
               for s in LOAD_SHADES]
    fig.legend(handles, [c[0] for c in LOAD_CLASSES], title="Mutations\nper read",
               frameon=False, loc="center right", fontsize=9)
    return _save(fig, fmt)


def project_figures(res: dict, fmt: str = "png") -> dict[str, bytes]:
    out = {
        f"density_by_group.{fmt}": group_dotplot(
            res, "unique_density_x1e5", "Unique mutations per 10⁵ bp", fmt),
        f"mutation_frequency_by_group.{fmt}": group_dotplot(
            res, "mutation_frequency_x1e3", "Mutations per 10³ bp", fmt),
        f"mutated_reads_by_group.{fmt}": group_dotplot(
            res, "mutated_reads_pct", "% of reads mutated", fmt),
        f"substitution_matrix.{fmt}": group_matrix(res, fmt),
        f"spectrum_by_group.{fmt}": group_spectrum(res, fmt),
        f"positional_map.{fmt}": group_positional(res, fmt),
        f"mutation_load.{fmt}": group_load(res, fmt),
    }
    if any(s["pct_GC_in_WRC_GYW"] is not None for s in res["samples"]):
        out[f"hotspot_targeting.{fmt}"] = group_dotplot(
            res, "pct_GC_in_WRC_GYW", "% of G:C mutations at WRC/GYW", fmt)
    return out


# ---- comparing two runs of the same samples (platforms, callers, depths) -----
MS_PINK = "#C2548A"


def spectrum_concordance(res_a: dict, res_b: dict, group: str,
                         label_a: str = "A", label_b: str = "B",
                         fmt: str = "png") -> bytes:
    """The manuscript's Figure 3C: one run's spectrum against another's.

    Each point is one of the twelve substitution classes, as a percentage of
    that run's calls in `group`. The dashed line is identity, so a point above
    it is a class the second run finds more often. Pearson r and the
    least-squares slope are printed on the plot; with twelve dependent
    categories they describe agreement, they do not test it.
    """
    ca, cb = res_a["spectrum_group"].get(group, {}), res_b["spectrum_group"].get(group, {})
    ta, tb = sum(ca.values()) or 1, sum(cb.values()) or 1
    xs = [100 * ca.get(k, 0) / ta for k in SUBS12]
    ys = [100 * cb.get(k, 0) / tb for k in SUBS12]
    n = len(SUBS12)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    r = sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")
    slope = sxy / sxx if sxx else float("nan")

    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    top = max(xs + ys) * 1.18 or 1.0
    ax.plot([0, top], [0, top], "--", color="#999", lw=1.1, zorder=1)
    ax.scatter(xs, ys, s=54, color=MS_PINK, zorder=3)
    for k, x, y in zip(SUBS12, xs, ys, strict=True):
        ax.annotate(k.replace(">", "→"), (x, y), textcoords="offset points",
                    xytext=(7, -3), fontsize=9, color="#333")
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel(f"{label_a} (% of calls)")
    ax.set_ylabel(f"{label_b} (% of calls)")
    ax.text(0.97, 0.06, f"Pearson r = {r:.2f}\nslope = {slope:.2f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=11, color=MS_BLUE, fontweight="bold")
    ax.set_title(f"{group}: n = {ta} and {tb} unique calls", fontsize=10, color="#555")
    ax.spines[["top", "right"]].set_visible(False)
    return _save(fig, fmt)


def compare_dotplots(runs: list[tuple[str, dict]], metric: str, label: str,
                     fmt: str = "png") -> bytes:
    """One `group_dotplot` panel per run, side by side on their own scales.

    This is the manuscript's Figure 2A and 2B: the same genotypes measured two
    ways. The scales are separate because the two assays count different things
    (reads against clones), so only the separation within a panel is comparable.
    """
    import matplotlib.image as mpimg
    panels = [(name, group_dotplot(res, metric, label if i == 0 else "", "png"))
              for i, (name, res) in enumerate(runs)]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 4.2), squeeze=False)
    for ax, (name, png) in zip(axes[0], panels, strict=True):
        ax.imshow(mpimg.imread(io.BytesIO(png)))
        ax.set_title(name, fontweight="bold")
        ax.axis("off")
    return _save(fig, fmt)
