"""How SHMtrace turns reads into calls: the decisions, the equations, and real data at each step.

Six stages, each with what it does, the rule it applies, and a panel of data from the J_H4
libraries showing that stage working. Everything plotted is read from out/intuition_jh4.json
and the project outputs; nothing is drawn by hand.
"""
import json

import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch
from style import CASE, CASEF, CTRL, CTRLF, DATA, FULL, INK, MID, MM, SOFT, save, setup

D = json.load(open(DATA / "intuition_jh4.json"))
ACCT = json.load(open(DATA.parents[1] / "figures/jh4_shmtrace/table_s1.json"))
ACC = "#1f6f8b"
BOX = "#f3f6f8"

STAGES = [
    ("1", "Read preparation",
     "Every read is gated before it is used: too long for the amplicon, or noisier on average than\n"
     "one error in a hundred bases, and it is dropped.",
     r"$\bar{E}=\frac{1}{L}\sum_{i=1}^{L}10^{-Q_i/10}\;\leq\;0.01$"),
    ("2", "Competitive alignment",
     "Each read is aligned against the target, its paralogue and any contaminant at once, and goes to\n"
     "whichever it matches best. A paralogue can be a target of its own (run_targets).",
     r"$t^{*}(r)=\arg\max_{t}\;\mathrm{score}(r,\,t)$"     "   over the target, its paralogue and every contaminant"),
    ("3", "Background measured in the controls",
     "The error rate is not assumed from the quality scores: it is measured in the AID-null libraries,\n"
     "per position, substitution, strand and quality, then inflated by one factor fitted by leaving\n"
     "each control out in turn.",
     r"$\hat{q}=\kappa\,\frac{\#\{\mathrm{miscalls}\;x\to a,\;\mathrm{strand}\;s,\;Q\geq q\}}"
     r"{\#\{\mathrm{control\;bases\;at}\;x,\;s,\;Q\geq q\}}$"),
    ("4", "One test per candidate",
     "Supporting reads are ranked by quality and every depth of evidence is tested against that\n"
     "measured rate; the best test is kept and corrected for having tried them all. A change carried\n"
     "by half the reads is germline, not a mutation.",
     r"$\lambda_j=\sum_{s}n_s\,\hat{q}_s(\geq q_{(j)}),\quad"
     r"p=k\,\min_{j\leq k}P\!\left(X\geq j\right),\;X\sim\mathrm{Poisson}(\lambda_j)$"),
    ("5", "Threshold set by the controls themselves",
     "Each control is called against the other controls, giving decoy p-values. The p-value threshold\n"
     "is the largest at which the false calls the decoys predict are at most 5% of the calls kept.",
     r"$\widehat{FP}(t)=\frac{1}{|D|}\sum_{d\in D}\left|\{p^{d}\leq t\}\right|\max\!\left(1,\frac{B}{B_d}\right),"
     r"\quad t^{*}=\max\{t:\widehat{FP}(t)\leq\alpha\left|\{p\leq t\}\right|\}$"),
    ("6", "One number per animal",
     "Calls are summarised per library and compared between groups on the log scale, with one value\n"
     "per animal, so clones of one lineage cannot stand in for independent mice.",
     r"$\rho=\frac{\#\,\mathrm{unique\;calls}}{\#\,\mathrm{analysed\;bases}}\times10^{5},"
     r"\quad \mathrm{Welch\;}t\mathrm{\;on\;}\log_{10}\rho,\;n=\mathrm{animals}$"),
]


def _stage_box(ax, num, title, body, eq):
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.003, 0.03), 0.994, 0.95, boxstyle="round,pad=0.004,rounding_size=0.02",
                                transform=ax.transAxes, facecolor=BOX, edgecolor="#d9e2e8", lw=0.7, zorder=0))
    lines = body.count("\n") + 1
    ax.text(0.022, 0.93, num, transform=ax.transAxes, fontsize=10.5, fontweight="bold", color=ACC, va="top")
    ax.text(0.062, 0.935, title, transform=ax.transAxes, fontsize=7.4, fontweight="bold", color=INK, va="top")
    ax.text(0.062, 0.745, body, transform=ax.transAxes, fontsize=5.9, color=MID, va="top", linespacing=1.4)
    # mathtext with fractions and sums is tall, so the equation is centred well clear of the text
    ax.text(0.062, max(0.18, 0.32 - 0.055 * lines), eq, transform=ax.transAxes, fontsize=6.4, color=INK, va="center")


def panel1(ax):
    """The fate of every delivered read, measured in the six J_H4 libraries."""
    parts = [("too_long", "over-length", "#5f6b75"), ("low_quality", "expected error > 0.01", "#9aa5b1"),
             ("not_target", "unmapped or paralogue", "#d9952b"), ("chimeric", "chimeric", "#c0392b"),
             ("short_alignment", "alignment too short", "#e6b98a"),
             ("paralogue_alleles", "paralogue alleles", "#6a4c93"), ("usable", "usable", ACC)]
    libs = list(ACCT)
    for i, lib in enumerate(libs):
        a, left = ACCT[lib], 0
        for key, lab, col in parts:
            v = 100 * a[key] / a["delivered"]
            ax.bar(i, v, bottom=left, color=col, width=0.74, lw=0, label=lab if i == 0 else None)
            left += v
    ax.set_xticks(range(len(libs)), libs, fontsize=4.6, rotation=35, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of delivered reads", fontsize=5.6)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=4.1, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=0.8,
              handletextpad=0.3, labelspacing=0.2, columnspacing=0.8)
    ax.set_title("The fate of every read", fontsize=6.2, pad=4)


def alignment_points():
    """Every read of one library scored against the target and against the paralogue (mappy).

    Cached in out/, and built here when it is not there, so the figure needs nothing but the
    bundled example library and the two references.
    """
    cache = DATA / "engine_alignment.json"
    if cache.exists():
        return json.load(open(cache))
    import mappy

    from shmtrace.io_utils import read_fasta, read_fastq
    root = DATA.parents[3]
    refs = {"target": read_fasta(str(root / "src/shmtrace/references/mouse_jh4_intron.fasta"))[0][1],
            "paralogue": dict(read_fasta(str(root / "src/shmtrace/references/mouse_jh4_contaminant.fasta")))
            ["JH4_paralog_sanger_consensus"]}
    idx = {k: mappy.Aligner(seq=v, preset="map-ont") for k, v in refs.items()}

    def ident(a, seq):
        h = max(a.map(seq), key=lambda x: x.mlen, default=None)
        return h.mlen / h.blen if h else 0.0

    out: dict = {"target": [], "paralogue": []}
    for _n, seq, _q in list(read_fastq(str(root / "examples/WT3_JH4_sub500.fastq")))[:400]:
        a, b = ident(idx["target"], seq), ident(idx["paralogue"], seq)
        out["target" if a >= b else "paralogue"].append([round(a, 4), round(b, 4)])
    DATA.mkdir(exist_ok=True)
    json.dump(out, open(cache, "w"))
    return out


def panel2(ax):
    """Identity to the target against identity to the paralogue, per read."""
    pts = alignment_points()
    for lab, key, col, f in (("to target", "target", CASE, CASEF), ("to paralogue", "paralogue", "#7a4f9a", "#bda3d4")):
        xy = pts[key]
        ax.scatter([q[0] for q in xy], [q[1] for q in xy], s=3.5, facecolor=f, edgecolor=col, lw=0.25, label=lab)
    lim = [min(q[0] for v in pts.values() for q in v) - 0.01, 1.005]
    ax.plot(lim, lim, color=SOFT, lw=0.6, ls="--")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("identity to target", fontsize=5.6)
    ax.set_ylabel("identity to paralogue", fontsize=5.6)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=4.8, handletextpad=0.2, borderaxespad=0.1, loc="upper left")
    ax.set_title("Each read goes to its best match", fontsize=6.2, pad=4)


def panel3(ax):
    """Reported base quality against the error actually measured in the controls."""
    q = D["quality"]
    xs, ys = [], []
    for k, (mm, bases) in sorted(q.items(), key=lambda kv: int(kv[0])):
        if bases >= 2000:
            xs.append(int(k))
            ys.append(-10 * np.log10(max(mm, 0.5) / bases))
    ax.plot(xs, ys, "o-", ms=2.2, lw=0.9, color=CTRL, markerfacecolor=CTRLF, markeredgecolor=CTRL, mew=0.5,
            label="measured in controls")
    lim = [min(xs), max(xs)]
    ax.plot(lim, lim, color=INK, lw=0.7, ls="--", label="what the score promises")
    ax.set_xlabel("reported quality", fontsize=5.6)
    ax.set_ylabel("measured quality", fontsize=5.6)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=4.8, handletextpad=0.3, borderaxespad=0.1, loc="upper left")
    ax.set_title("Quality scores are optimistic", fontsize=6.2, pad=4)


def panel4(ax):
    """Two look-alike candidates: the one the controls clear, and the one they do not."""
    sites = D["anatomy"]["sites"]
    labels, case_v, ctl_v, called = [], [], [], []
    for key in ("true", "error"):
        s = sites.get(key)
        if not s:
            continue
        libs = s["libraries"]
        case = next(v for k, v in libs.items() if k == D["anatomy"]["case"])
        ctls = [v for k, v in libs.items() if k != D["anatomy"]["case"]]
        labels.append(f"{s['ref']}{s['pos']}{s['alt']}\n{s['context']}")
        case_v.append(case["alt_plus"] + case["alt_minus"])
        ctl_v.append(np.mean([c["alt_plus"] + c["alt_minus"] for c in ctls]))
        called.append(s["called"])
    x = np.arange(len(labels))
    ax.bar(x - 0.18, case_v, width=0.34, color=CASEF, edgecolor=CASE, lw=0.6, label="wild type")
    ax.bar(x + 0.18, ctl_v, width=0.34, color=CTRLF, edgecolor=CTRL, lw=0.6, label="$Aicda^{-/-}$")
    for i, c in enumerate(called):
        ax.text(i, max(case_v[i], ctl_v[i]) + 1.4, "called" if c else "not called", ha="center",
                fontsize=5.2, color=CASE if c else MID, fontweight="bold" if c else "normal")
    ax.set_xticks(x, labels, fontsize=5)
    ax.set_ylim(0, max(case_v + ctl_v) * 1.45)
    ax.set_ylabel("reads with the change", fontsize=5.6)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=4.8, handletextpad=0.3, borderaxespad=0.1, loc="upper right")
    ax.set_title("Same depth, opposite verdicts", fontsize=6.2, pad=4)


def panel5(ax):
    """Candidates kept against p-value, for a wild-type library and for the decoys."""
    case = D["anatomy"]["case"]
    t = D["anatomy"]["threshold"]
    ps = np.sort(np.array(D["candidates"][case]["p"]))
    grid = np.logspace(-12, 0, 160)
    ax.plot(grid, [np.searchsorted(ps, g, "right") for g in grid], color=CASE, lw=1.1, label="wild type")
    for i, (_name, d) in enumerate(sorted(D["decoys"].items())):
        dp = np.sort(np.array(d[0] if isinstance(d, list) and d and isinstance(d[0], list) else d))
        ax.plot(grid, [np.searchsorted(dp, g, "right") for g in grid], color=CTRL, lw=0.7, alpha=0.75,
                label="held-out controls (decoys)" if i == 0 else None)
    ax.axvline(t, color=INK, lw=0.8, ls="--")
    ax.text(t, ax.get_ylim()[1] * 0.97, f"  $t^*$ = {t:.3f}", fontsize=5, color=INK, va="top")
    ax.set_xscale("log")
    ax.set_xlabel("p-value threshold", fontsize=5.6)
    ax.set_ylabel("candidates kept", fontsize=5.6)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=4.8, handletextpad=0.3, borderaxespad=0.1, loc="upper left")
    ax.set_title("The controls set the threshold", fontsize=6.2, pad=4)


def panel6(ax):
    """Unique-variant density per animal, wild type against Aicda-/-."""
    calls = D["calls"]
    win = D["window"]
    bases = (win[1] - win[0] + 1) * D["depth"]
    groups = [("wild type", D["cases"], CASE, CASEF), ("$Aicda^{-/-}$", D["controls"], CTRL, CTRLF)]
    for i, (_lab, names, col, face) in enumerate(groups):
        vals = [max(len(calls[n]), 0.5) / bases * 1e5 for n in names if n in calls]
        ax.scatter(np.full(len(vals), i) + np.linspace(-0.09, 0.09, len(vals)), vals, s=14,
                   facecolor=face, edgecolor=col, lw=0.6, zorder=3)
        g = 10 ** np.mean(np.log10(vals))
        ax.plot([i - 0.22, i + 0.22], [g, g], color=col, lw=1.3)
    ax.set_yscale("log")
    ax.set_xticks(range(2), [g[0] for g in groups], fontsize=5.6)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("unique variants / 10$^5$ bases", fontsize=5.6)
    ax.tick_params(labelsize=5)
    ax.set_title("One point per animal", fontsize=6.2, pad=4)


def main():
    setup()
    fig = plt.figure(figsize=(FULL, 210 * MM))
    gs = gridspec.GridSpec(6, 2, figure=fig, left=0.012, right=0.955, top=0.93, bottom=0.028,
                           width_ratios=[2.2, 1], wspace=0.28, hspace=0.6)
    fig.text(0.012, 0.978, "How SHMtrace turns reads into calls", fontsize=10, fontweight="bold", color=INK)
    fig.text(0.012, 0.9655, "Every rule the engine applies, and the same rule working on the J$_H$4 libraries "
                            "(3 wild-type and 3 $Aicda^{-/-}$ mice, 500 reads each).",
             fontsize=6.2, color=MID)
    panels = [panel1, panel2, panel3, panel4, panel5, panel6]
    for i, ((num, title, body, eq), panel) in enumerate(zip(STAGES, panels, strict=True)):
        _stage_box(fig.add_subplot(gs[i, 0]), num, title, body, eq)
        ax = fig.add_subplot(gs[i, 1])
        panel(ax)
    save(fig, "Figure_Engine")


if __name__ == "__main__":
    main()
