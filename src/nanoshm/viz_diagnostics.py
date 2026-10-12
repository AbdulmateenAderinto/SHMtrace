"""Diagnostic and detail figures, beyond the manuscript's seven.

Three kinds:

* detail the project already computes but never drew: the saturation curve,
  the single-read excess, indels, per-sample positional maps;
* diagnostics that make the method auditable rather than taken on trust: the
  base-quality histogram the quality-calibrated model keys on, the read-by-site
  carriage behind a phasing decision, the read-attrition waterfall, the fitted
  control background, and the caller's own decision boundary;
* biology the tables hold but no figure shows: clonal structure and the
  agreement between replicates of a group.

Every function takes the dict `project.run_project` returns and gives back
image bytes, as `viz_project` does. Functions return None when the result does
not carry what they need, and `diagnostic_figures` simply skips those.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .viz_project import MS_BLUE, MS_RED, SUBS12, _group_colour, _save

GREY = "#8A8F98"
FAINT = "#D7DBE0"


def _details(res: dict) -> dict:
    return res.get("_details") or {}


def _groups_of(res: dict) -> dict:
    return {s["sample"]: s["group"] for s in res["samples"]}


def _ax(figsize=(7.2, 4.2)):
    fig, ax = plt.subplots(figsize=figsize)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return fig, ax


# ---- 1. saturation -----------------------------------------------------------

def saturation(res: dict, fmt: str = "png") -> bytes | None:
    """Unique mutations against sequencing depth, one line per sample.

    Density saturates with depth, so a density is only meaningful at a stated
    depth. The curve shows whether the analysed depth sits on the plateau or
    on the steep part, where another thousand reads would change the answer.
    """
    curves = {n: d.get("rarefaction") for n, d in _details(res).items()}
    curves = {n: c for n, c in curves.items() if c}
    if not curves:
        return None
    colour = _group_colour(res)
    groups = _groups_of(res)
    fig, ax = _ax()
    ends: list = []
    for name, rows in sorted(curves.items()):
        xs = [r["depth"] for r in rows]
        ys = [r["unique_variants_mean"] for r in rows]
        sd = [r["unique_variants_sd"] for r in rows]
        c = colour.get(groups.get(name), MS_BLUE)
        ax.errorbar(xs, ys, yerr=sd, color=c, lw=1.4, marker="o", ms=3,
                    capsize=2, alpha=0.85)
        if xs:
            ends.append((ys[-1], xs[-1], name, c))
    # labels at the right-hand end, nudged apart when curves finish together
    span = max((y for y, *_ in ends), default=1) or 1
    last = None
    for y, x, name, c in sorted(ends):
        y_lab = y if last is None else max(y, last + 0.045 * span)
        ax.annotate(name, (x, y), fontsize=7, color=c, va="center",
                    xytext=(7, 0), textcoords="offset points",
                    xycoords="data", annotation_clip=False)
        last = y_lab
    depth = (res.get("params") or {}).get("depth")
    if depth:
        ax.axvline(depth, color=GREY, ls="--", lw=1)
        ax.annotate(f"analysed at {depth:,}", (depth, ax.get_ylim()[1]), fontsize=7,
                    color=GREY, rotation=90, va="top", xytext=(-11, -4),
                    textcoords="offset points")
    ax.set_xlabel("Reads sampled")
    ax.set_ylabel("Unique mutations")
    ax.set_title("Saturation: is the depth enough?", loc="left", fontsize=11)
    ax.margins(x=0.12)
    return _save(fig, fmt)


# ---- 2. single-read excess ---------------------------------------------------

def single_read_excess(res: dict, fmt: str = "png") -> bytes | None:
    """Single-read events observed against the control background, by class.

    A mutation in one read cannot be confirmed site by site, so it is estimated
    in aggregate. This is that estimate made visible: for each substitution
    class, what the sample shows and what the controls predict. The gap is the
    signal; the controls' own bar is the price of admission.
    """
    sr = {k: v for k, v in (res.get("single_read") or {}).items() if v}
    if not sr:
        return None
    groups = _groups_of(res)
    ctrl = {g["group"] for g in res.get("groups") or [] if g.get("control")}
    obs, exp = Counter(), Counter()
    for name, est in sr.items():
        if groups.get(name) in ctrl:
            continue
        obs.update(est.get("observed_by_class") or {})
        for k, v in (est.get("expected_by_class") or {}).items():
            exp[k] += v
    if not obs and not exp:
        return None
    fig, ax = _ax((7.6, 4.0))
    x = range(len(SUBS12))
    ax.bar([i - 0.2 for i in x], [obs.get(k, 0) for k in SUBS12], width=0.4,
           color=MS_BLUE, label="observed in samples")
    ax.bar([i + 0.2 for i in x], [exp.get(k, 0) for k in SUBS12], width=0.4,
           color=FAINT, edgecolor=GREY, label="expected from controls")
    ax.set_xticks(list(x))
    ax.set_xticklabels([s.replace(">", "→") for s in SUBS12], fontsize=8)
    ax.set_ylabel("Single-read events")
    ax.set_title("Single-read mutations against the control background",
                 loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, fmt)


# ---- 3. indels ---------------------------------------------------------------

def indel_figure(res: dict, fmt: str = "png") -> bytes | None:
    """Indel size, direction and position. Only drawn when indels were called."""
    rows = res.get("indels") or []
    if not rows:
        return None
    groups = _groups_of(res)
    colour = _group_colour(res)
    order = [g for g in res["group_order"] if any(groups.get(r["sample"]) == g for r in rows)]
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.8))
    for ax in axes:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    ax = axes[0]                                   # size distribution
    bins = [1, 2, 3, 4, 6, 11, 21, 51, 10 ** 6]
    labels = ["1", "2", "3", "4-5", "6-10", "11-20", "21-50", ">50"]
    for g in order:
        lens = [r["length"] for r in rows if groups.get(r["sample"]) == g]
        counts = [sum(1 for v in lens if bins[i] <= v < bins[i + 1]) for i in range(len(labels))]
        n = sum(counts) or 1
        ax.plot(labels, [100 * c / n for c in counts], marker="o", ms=4,
                color=colour.get(g, MS_BLUE), label=f"{g} (n={n})")
    ax.set_xlabel("Indel length (bp)")
    ax.set_ylabel("% of the group's indels")
    ax.set_title("Size", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]                                   # insertions against deletions
    # every group, including one with no indels at all: a flat zero next to a
    # tall bar is the comparison
    all_groups = list(res["group_order"])
    n_per = Counter(s["sample"] for s in res["samples"])
    for i, g in enumerate(all_groups):
        members = [s["sample"] for s in res["samples"] if s["group"] == g] or [None]
        per = {m: [0, 0] for m in members}
        for r in rows:
            if groups.get(r["sample"]) == g:
                per.setdefault(r["sample"], [0, 0])[0 if r["kind"] == "del" else 1] += 1
        n = max(1, len(members))
        ax.bar(i - 0.18, sum(v[0] for v in per.values()) / n, width=0.34,
               color=colour.get(g, MS_BLUE))
        ax.bar(i + 0.18, sum(v[1] for v in per.values()) / n, width=0.34,
               color=colour.get(g, MS_BLUE), alpha=0.45)
    ax.set_xticks(range(len(all_groups)))
    ax.set_xticklabels(all_groups, fontsize=9)
    ax.set_ylabel("Per sample (mean)")
    ax.set_title("Deletions (solid) and insertions (pale)", loc="left", fontsize=10)

    ax = axes[2]                                   # position along the amplicon
    for g in order:
        pos = [r["pos"] for r in rows if groups.get(r["sample"]) == g]
        if pos:
            ax.hist(pos, bins=28, histtype="step", lw=1.5,
                    color=colour.get(g, MS_BLUE), label=g)
    ax.set_xlabel("Reference position")
    ax.set_ylabel("Indels")
    ax.set_title("Position", loc="left", fontsize=10)
    fig.suptitle("Insertions and deletions", x=0.007, ha="left", fontsize=11)
    fig.tight_layout()
    return _save(fig, fmt)


# ---- 4. per-sample positional map --------------------------------------------

def positional_by_sample(res: dict, fmt: str = "png") -> bytes | None:
    """One panel per sample, so a group mean cannot hide a single animal."""
    rows = res.get("positional_sample") or []
    if not rows:
        return None
    groups = _groups_of(res)
    colour = _group_colour(res)
    names = [s["sample"] for s in res["samples"]]
    per = defaultdict(list)
    for r in rows:
        per[r["sample"]].append((r["pos"], r["frequency_pct"]))
    n = len(names)
    ncol = 3
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 1.9 * nrow),
                             sharex=True, sharey=True, squeeze=False)
    hi = max((v for vals in per.values() for _p, v in vals), default=1.0)
    win = (res.get("params") or {}).get("analysis_window") or [None, None]
    for i, name in enumerate(names):
        ax = axes[i // ncol][i % ncol]
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        pts = sorted(per.get(name, []))
        if pts:
            ax.vlines([p for p, _v in pts], 0, [v for _v2, v in pts],
                      color=colour.get(groups.get(name), MS_BLUE), lw=0.7)
        ax.set_title(f"{name}  ({groups.get(name)})", loc="left", fontsize=9)
        ax.set_ylim(0, hi * 1.08)
        if win[0]:
            ax.set_xlim(win[0], win[1])
    for i in range(n, nrow * ncol):
        axes[i // ncol][i % ncol].axis("off")
    fig.supxlabel("Reference position", fontsize=9)
    fig.supylabel("% of reads mutated", fontsize=9)
    fig.tight_layout()
    return _save(fig, fmt)


# ---- 5. combined density with intervals --------------------------------------

def combined_density(res: dict, fmt: str = "png") -> bytes | None:
    """Per-site calls plus the estimated single-read mutations, with 95% intervals.

    The plain density plot shows point estimates of what the test called. This
    adds the single-read estimate and, more importantly, its uncertainty.
    """
    samples = [s for s in res["samples"] if s.get("combined_density_x1e5") is not None]
    if not samples:
        return None
    colour = _group_colour(res)
    order = [g for g in res["group_order"] if any(s["group"] == g for s in samples)]
    fig, ax = _ax((max(6.0, 1.1 * len(samples)), 4.2))
    xs, labels = [], []
    for i, s in enumerate(sorted(samples, key=lambda s: (order.index(s["group"]), s["sample"]))):
        y = s["combined_density_x1e5"]
        lo = max(0.0, s.get("combined_density_ci_low") or y)
        hi = s.get("combined_density_ci_high") or y
        c = colour.get(s["group"], MS_BLUE)
        ax.errorbar(i, y, yerr=[[y - lo], [hi - y]], fmt="o", ms=6, color=c,
                    capsize=3, lw=1.3)
        ax.plot(i, s["unique_density_x1e5"], marker="_", ms=11, color=GREY, zorder=3)
        xs.append(i)
        labels.append(s["sample"])
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mutations per 10⁵ bp")
    ax.set_title("Combined density with 95% intervals "
                 "(grey dash: per-site calls only)", loc="left", fontsize=11)
    return _save(fig, fmt)


# ---- 6. the base-quality histogram -------------------------------------------

def quality_histogram(res: dict, fmt: str = "png") -> bytes | None:
    """Base quality of mismatching against matching bases.

    This is the assumption the quality-calibrated model rests on, drawn rather
    than asserted. Sequencing error piles up at low quality; a real mutation is
    read at ordinary quality, like any other base. The further apart the two
    distributions, the more a quality-aware test can do.
    """
    mm, ref = Counter(), Counter()
    for d in _details(res).values():
        c = (d.get("summary") or {}).get("calling") or {}
        for q, n in (c.get("variant_base_quality_hist") or {}).items():
            mm[int(q)] += n
        for q, n in (c.get("reference_base_quality_hist") or {}).items():
            ref[int(q)] += n
    if not mm:
        return None
    fig, ax = _ax((7.4, 4.0))
    qs = range(0, 61)
    tm, tr = sum(mm.values()) or 1, sum(ref.values()) or 1
    ax.fill_between(list(qs), [100 * ref.get(q, 0) / tr for q in qs], step="mid",
                    color=GREY, alpha=0.25, label=f"matches the reference (n={tr:,})")
    ax.step(list(qs), [100 * mm.get(q, 0) / tm for q in qs], where="mid",
            color=MS_RED, lw=1.6, label=f"mismatches (n={tm:,})")
    top = max(mm.values()) / tm * 100 if mm else 1
    if mm.get(60) or mm.get(50):
        cap = 60 if mm.get(60) else 50
        ax.annotate("basecaller quality cap", (cap, 100 * mm[cap] / tm), fontsize=7,
                    color=GREY, ha="right", va="bottom",
                    xytext=(-6, 2), textcoords="offset points")
    ax.set_ylim(0, top * 1.18)
    minq = next((((d.get("summary") or {}).get("calling") or {}).get("min_base_quality")
                 for d in _details(res).values()), None)
    if minq:
        ax.axvline(minq - 0.5, color=GREY, ls="--", lw=1)
        ax.annotate(f"min_base_quality {minq}", (minq, ax.get_ylim()[1]), fontsize=7,
                    color=GREY, va="top", xytext=(3, -4), textcoords="offset points")
    hiq = max(max(mm, default=0), max(ref, default=0))
    ax.set_xlim(-1, hiq + 1)
    ax.set_xlabel("Base quality (Phred)")
    ax.set_ylabel("% of that population's bases")
    ax.set_title("What the quality-calibrated model keys on", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, fmt)


# ---- 7. phasing ---------------------------------------------------------------

def phasing_blocks(res: dict, fmt: str = "png") -> bytes | None:
    """How many of a block's sites each read carries.

    Somatic hypermutation does not co-occur: every read is a different B cell,
    so reads carry their own scattered mutations and the distribution is
    graded. A second haplotype, paralogue or contaminant is all-or-nothing, so
    the distribution is U-shaped. That shape, not the frequency, is what
    decides whether a block is dropped.
    """
    blocks = []
    for name, d in _details(res).items():
        for b in ((d.get("summary") or {}).get("reads") or {}).get("blocks") or []:
            if b.get("carried_hist"):
                blocks.append((name, b))
    if not blocks:
        return None
    blocks = blocks[:6]
    fig, axes = plt.subplots(1, len(blocks), figsize=(3.5 * len(blocks), 3.4),
                             squeeze=False)
    for ax, (name, b) in zip(axes[0], blocks, strict=True):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        hist = {int(k): v for k, v in b["carried_hist"].items()}
        ks = list(range(0, b["n_sites"] + 1))
        drop = b.get("action") == "drop"
        ax.bar(ks, [hist.get(k, 0) for k in ks],
               color=MS_RED if drop else MS_BLUE, alpha=0.85)
        ax.set_title(f"{name}: {b['n_sites']} sites\n{b.get('action')} — "
                     f"{b.get('why', '')[:42]}", loc="left", fontsize=8)
        ax.set_xlabel("Block sites carried by a read")
        ax.set_ylabel("Reads")
    fig.suptitle("Phasing: all-or-nothing means a haplotype, graded means hypermutation",
                 x=0.007, ha="left", fontsize=11)
    fig.tight_layout()
    return _save(fig, fmt)


# ---- 8. read attrition --------------------------------------------------------

_STEPS = [("reads_in", "input"), ("passed_prefilter", "quality + length"),
          ("aligned_before_phasing", "on target"), ("aligned_to_target", "after phasing"),
          ("downsampled_to", "analysed")]


def qc_waterfall(res: dict, fmt: str = "png") -> bytes | None:
    """Reads surviving each stage, per sample.

    A library that is largely contaminant or largely low quality should be
    visible before anyone reads a p-value.
    """
    det = _details(res)
    if not det:
        return None
    names = [s["sample"] for s in res["samples"] if s["sample"] in det]
    rows = []
    for n in names:
        r = (det[n].get("summary") or {}).get("reads") or {}
        rows.append([r.get(k) for k, _lab in _STEPS])
    if not any(any(v for v in r) for r in rows):
        return None
    groups = _groups_of(res)
    colour = _group_colour(res)
    labels = [lab for _k, lab in _STEPS]
    fig, ax = _ax((7.4, 4.2))
    worst = 100.0
    for n, vals in zip(names, rows, strict=True):
        first = next((v for v in vals if v), 0) or 1
        xs = [i for i, v in enumerate(vals) if v is not None]
        ys = [100 * vals[i] / first for i in xs]
        if not xs:
            continue
        worst = min(worst, min(ys))
        c = colour.get(groups.get(n), MS_BLUE)
        ax.plot(xs, ys, marker="o", ms=4, lw=1.3, color=c, alpha=0.85)
        ax.annotate(n, (xs[-1], ys[-1]), fontsize=7, color=c, va="center",
                    xytext=(6, 0), textcoords="offset points", annotation_clip=False)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("% of input reads surviving")
    ax.set_ylim(max(0, worst - 8), 103)
    ax.margins(x=0.16)
    ax.set_title("Where reads are lost", loc="left", fontsize=11)
    return _save(fig, fmt)


# ---- 9. the fitted control background ----------------------------------------

def control_background(res: dict, fmt: str = "png") -> bytes | None:
    """The error rate the controls measured, position by position.

    Every call is judged against this. Spikes are positions where the controls
    themselves err often, and a call there needs much more support. A flat,
    low line means the controls are doing their job.
    """
    bg = res.get("control_background") or {}
    if not bg:
        return None
    xs = sorted(int(k) for k in bg)
    ys = [bg[str(x)] for x in xs]
    fig, ax = _ax((7.6, 3.6))
    ax.plot(xs, ys, lw=0.8, color=MS_RED)
    ax.fill_between(xs, 0, ys, color=MS_RED, alpha=0.18)
    med = sorted(ys)[len(ys) // 2] if ys else 0
    ax.axhline(med, color=GREY, ls="--", lw=1)
    ax.annotate(f"median {med:.2g}", (xs[0], med), fontsize=7, color=GREY,
                va="bottom", xytext=(2, 2), textcoords="offset points")
    ax.set_yscale("log")
    ax.set_xlabel("Reference position")
    ax.set_ylabel("Fitted error rate")
    ax.set_title("The background your controls measured", loc="left", fontsize=11)
    return _save(fig, fmt)


# ---- 10. clonality ------------------------------------------------------------

def clonality(res: dict, fmt: str = "png") -> bytes | None:
    """How concentrated a sample's mutated reads are in shared profiles.

    Reads are grouped by the exact set of called positions they carry. Under
    independent hypermutation almost every mutated read is its own profile; an
    expanded clone shows up as one profile carried by many reads. The bar is
    the share of mutated reads in the largest profile.
    """
    det = _details(res)
    if not det:
        return None
    groups = _groups_of(res)
    colour = _group_colour(res)
    names, top, uniq = [], [], []
    for s in res["samples"]:
        rows = (det.get(s["sample"]) or {}).get("per_read") or []
        prof = Counter(tuple(r["positions"]) for r in rows if r["positions"])
        if not prof:
            continue
        n = sum(prof.values())
        names.append(s["sample"])
        top.append(100 * prof.most_common(1)[0][1] / n)
        uniq.append(100 * sum(1 for v in prof.values() if v == 1) / n)
    if not names:
        return None
    fig, ax = _ax((max(6.2, 1.05 * len(names)), 4.0))
    for i, n in enumerate(names):
        c = colour.get(groups.get(n), MS_BLUE)
        ax.bar(i - 0.19, top[i], width=0.36, color=c)
        ax.bar(i + 0.19, uniq[i], width=0.36, color=c, alpha=0.4)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("% of mutated reads")
    ax.set_title("Clonal structure: largest shared profile (solid), "
                 "one-off profiles (pale)", loc="left", fontsize=11)
    return _save(fig, fmt)


# ---- 11. the decision boundary ------------------------------------------------

def volcano(res: dict, fmt: str = "png") -> bytes | None:
    """Every candidate the caller tested: frequency against p-value.

    Filled points were called, open points were not. It shows what the test
    threw away, which is otherwise invisible.
    """
    kept_x, kept_y, drop_x, drop_y = [], [], [], []
    for d in _details(res).values():
        for freq, pv, called in ((d.get("summary") or {}).get("calling") or {}).get("_candidates") or []:
            y = -math.log10(max(pv, 1e-300)) if pv is not None else 0.0
            (kept_x if called else drop_x).append(max(freq, 1e-6))
            (kept_y if called else drop_y).append(y)
    if not kept_x and not drop_x:
        return None
    # a handful of p = 0 calls would squash everything else, so clip the axis
    allv = sorted(kept_y + drop_y)
    cap = allv[int(0.99 * (len(allv) - 1))] if allv else 1.0
    cap = max(cap, 10.0)
    n_above = sum(1 for v in kept_y + drop_y if v > cap)
    ky = [min(v, cap) for v in kept_y]
    dy = [min(v, cap) for v in drop_y]
    fig, ax = _ax((7.0, 4.4))
    ax.scatter(drop_x, dy, s=7, facecolors="none", edgecolors=GREY,
               lw=0.5, alpha=0.5, label=f"not called ({len(drop_x):,})")
    ax.scatter(kept_x, ky, s=9, color=MS_BLUE, alpha=0.7,
               label=f"called ({len(kept_x):,})")
    if kept_y:
        # the largest p that still got called: the effective cutoff
        edge = min(kept_y)
        ax.axhline(edge, color=MS_RED, ls="--", lw=1)
        ax.annotate(f"weakest call kept, p = {10 ** -edge:.1g}", (ax.get_xlim()[0], edge),
                    fontsize=7, color=MS_RED, va="bottom",
                    xytext=(3, 2), textcoords="offset points")
    if n_above:
        ax.annotate(f"{n_above} point(s) above the axis", (0.99, 0.99), fontsize=7,
                    color=GREY, ha="right", va="top", xycoords="axes fraction")
    ax.set_xscale("log")
    ax.set_xlabel("Variant frequency in the column")
    ax.set_ylabel("−log₁₀ p")
    ax.set_title("What the caller kept, and what it did not", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    return _save(fig, fmt)


# ---- 12. replicate agreement ---------------------------------------------------

def replicate_concordance(res: dict, fmt: str = "png") -> bytes | None:
    """Per-position mutated-read counts, each sample against its group's others.

    Replicates of one group should agree position by position. Scatter far from
    the diagonal, or one animal carrying all the signal, is worth knowing.
    """
    rows = res.get("positional_sample") or []
    if not rows:
        return None
    groups = _groups_of(res)
    by = defaultdict(lambda: defaultdict(float))
    for r in rows:
        by[r["sample"]][r["pos"]] = r["frequency_pct"]
    colour = _group_colour(res)
    pts = defaultdict(lambda: ([], []))
    for g in res["group_order"]:
        members = [s for s in by if groups.get(s) == g]
        if len(members) < 2:
            continue
        positions = {p for m in members for p in by[m]}
        for m in members:
            others = [o for o in members if o != m]
            for pos in positions:
                mean_other = sum(by[o].get(pos, 0.0) for o in others) / len(others)
                pts[g][0].append(max(by[m].get(pos, 0.0), 1e-3))
                pts[g][1].append(max(mean_other, 1e-3))
    if not pts:
        return None
    fig, ax = _ax((5.4, 5.0))
    for g, (x, y) in pts.items():
        ax.scatter(x, y, s=9, alpha=0.45, color=colour.get(g, MS_BLUE), label=g)
    lim = max(max(x) for x, _y in pts.values()) * 1.3
    ax.plot([1e-3, lim], [1e-3, lim], color=GREY, ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("This sample, % of reads mutated")
    ax.set_ylabel("Mean of the group's other samples")
    ax.set_title("Do replicates agree position by position?", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, fmt)


def _any_detail(res: dict, *path) -> bool:
    for d in _details(res).values():
        cur = d
        for key in path:
            cur = (cur or {}).get(key) if isinstance(cur, dict) else None
        if cur:
            return True
    return False


def _has_blocks(res: dict) -> bool:
    return any(b.get("carried_hist")
               for d in _details(res).values()
               for b in ((d.get("summary") or {}).get("reads") or {}).get("blocks") or [])


def _has_clones(res: dict) -> bool:
    return any(r.get("positions")
               for d in _details(res).values() for r in (d.get("per_read") or []))


def _has_replicates(res: dict) -> bool:
    if not res.get("positional_sample"):
        return False
    counts = Counter(s["group"] for s in res.get("samples", []))
    return any(v >= 2 for v in counts.values())


# (file name, function, label, a cheap test of whether the data is there).
# The predicate exists so `project.public` can list the drawable figures for
# the browser without rendering them all first.
FIGURES = [
    ("saturation", saturation, "Saturation (depth)",
     lambda r: _any_detail(r, "rarefaction")),
    ("single_read_excess", single_read_excess, "Single-read excess",
     lambda r: bool(r.get("single_read"))),
    ("indels", indel_figure, "Insertions and deletions",
     lambda r: bool(r.get("indels"))),
    ("positional_by_sample", positional_by_sample, "Position, per sample",
     lambda r: bool(r.get("positional_sample"))),
    ("combined_density", combined_density, "Combined density with intervals",
     lambda r: any(s.get("combined_density_x1e5") is not None for s in r.get("samples", []))),
    ("quality_histogram", quality_histogram, "Base-quality histogram",
     lambda r: _any_detail(r, "summary", "calling", "variant_base_quality_hist")),
    ("phasing_blocks", phasing_blocks, "Phasing blocks", _has_blocks),
    ("qc_waterfall", qc_waterfall, "Where reads are lost",
     lambda r: _any_detail(r, "summary", "reads")),
    ("control_background", control_background, "Control background",
     lambda r: bool(r.get("control_background"))),
    ("clonality", clonality, "Clonal structure", _has_clones),
    ("volcano", volcano, "Caller decision boundary",
     lambda r: _any_detail(r, "summary", "calling", "_candidates")),
    ("replicate_concordance", replicate_concordance, "Replicate agreement",
     _has_replicates),
]


def drawable(res: dict) -> list[str]:
    """Names of the diagnostics this result can support, without drawing them."""
    out = []
    for name, _fn, _label, ok in FIGURES:
        try:
            if ok(res):
                out.append(name)
        except Exception:
            pass
    return out


def diagnostic_figures(res: dict, fmt: str = "png") -> dict[str, bytes]:
    """Every diagnostic the result can support, skipping those it cannot."""
    out = {}
    for name, fn, _label, _ok in FIGURES:
        try:
            blob = fn(res, fmt)
        except Exception:                      # a figure must never fail a run
            blob = None
        if blob:
            out[f"{name}.{fmt}"] = blob
    return out
