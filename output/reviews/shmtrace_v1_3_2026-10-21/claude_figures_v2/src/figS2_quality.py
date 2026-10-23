"""Figure S2: reported base quality against the error measured in the AID-null controls, per run.

Each point is the mismatch rate among control bases of one reported quality (germline sites excluded; bins of
two quality values with at least 2,000 bases); the diagonal is the rate the quality score promises,
10^(-Q/10). The four loci were sequenced in three runs (J_H4; N114; Smu and Pim1 from PRJNA1481264).
Data: out/intuition_<dataset>.json (intuition_data.py).
"""
import json

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from style import SCRATCH  # noqa: E402
from style import DATA, FULL, INK, MID, MM, SOFT, letter, save, setup

SETS = {"jh4": ("J$_H$4 (this study)", "#274b7e", "o"), "n114": ("N114 V$_H$ (this study)", "#6d9ad0", "s"),
        "su_lps": ("Sμ (PRJNA1481264)", "#1f6f8b", "^"), "pim1": ("$Pim1$ (PRJNA1481264)", "#9a3439", "D")}
BIN, MIN_BASES = 2, 2000


def cp(m, n):
    """95% Clopper-Pearson interval for m mismatches in n bases."""
    lo = stats.beta.ppf(0.025, m, n - m + 1) if m else 0.0
    hi = stats.beta.ppf(0.975, m + 1, n - m) if m < n else 1.0
    return lo, hi


def main():
    setup()
    # full width, the two panels side by side, in the style of the other figures (it was drawn 85 mm wide and
    # placed at page width, nearly twice its size)
    fig = plt.figure(figsize=(FULL, 72 * MM))
    ax = fig.add_axes([0.09, 0.17, 0.5, 0.64])
    axn = fig.add_axes([0.71, 0.17, 0.27, 0.64])
    out = {}
    for key, (label, col, mk) in SETS.items():
        q = json.load(open(DATA / f"intuition_{key}.json"))["quality"]
        mm, n = {}, {}
        for qq, (m, b) in q.items():
            k = int(qq) // BIN * BIN
            mm[k] = mm.get(k, 0) + m
            n[k] = n.get(k, 0) + b
        ks = [k for k in sorted(n) if n[k] >= MIN_BASES]
        seen = [k for k in ks if mm[k]]
        xs, ys = [k + BIN / 2 for k in seen], [mm[k] / n[k] for k in seen]
        ax.plot(xs, ys, marker=mk, ms=2.6, lw=0.8, color=col, label=label)
        ci = [cp(mm[k], n[k]) for k in ks]
        ax.fill_between([k + BIN / 2 for k in ks], [max(c[0], 1e-6) for c in ci], [c[1] for c in ci], color=col,
                        alpha=0.12, lw=0, step=None)
        # bins with no mismatch: only an upper bound
        zero = [k for k in ks if not mm[k]]
        ax.scatter([k + BIN / 2 for k in zero], [cp(0, n[k])[1] for k in zero], marker="v", s=7, color=col,
                   linewidths=0, zorder=3)
        axn.plot([k + BIN / 2 for k in ks], [n[k] for k in ks], marker=mk, ms=1.8, lw=0.6, color=col)
        tot_m, tot_b = sum(m for qq, (m, b) in q.items() if int(qq) >= 20), sum(b for qq, (m, b) in q.items()
                                                                                 if int(qq) >= 20)
        out[key] = tot_m / tot_b
    qs = np.arange(0, 51)
    ax.plot(qs, 10 ** (-qs / 10), color=INK, lw=0.7, ls="--", label="what the quality score promises")
    ax.set_yscale("log")
    ax.set_xlim(0, 52)
    ax.set_ylim(1e-6, 1)
    ax.set_ylabel("Mismatch rate in AID-null controls")
    ax.set_xlabel("Reported base quality")
    ax.legend(fontsize=6, loc="lower left", handlelength=1.6, borderaxespad=0.3)
    ax.set_title("Reported quality against the error the controls show", fontsize=7, pad=10)
    ax.text(0.0, 1.012, "shading, 95% Clopper–Pearson interval; ▼, a bin with no mismatch (its upper bound)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=MID)
    letter(ax, "A", dx=-0.13, dy=1.08)
    axn.set_yscale("log")
    axn.set_xlim(0, 52)
    axn.set_ylabel("Control bases per bin")
    axn.set_xlabel("Reported base quality")
    axn.set_title("How many bases each bin holds", fontsize=7, pad=10)
    axn.text(0.0, 1.012, f"bins of {BIN} quality values; shown with ≥ {MIN_BASES:,} bases", transform=axn.transAxes,
             ha="left", va="bottom", fontsize=6, color=MID)
    letter(axn, "B", dx=-0.3, dy=1.08)
    save(fig, "FigureS2")
    json.dump({"error_q20plus": out}, open(SCRATCH / "figS2_numbers.json", "w"), indent=1)


if __name__ == "__main__":
    main()
