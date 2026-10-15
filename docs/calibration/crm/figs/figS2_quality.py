"""Figure S2: reported base quality against the error measured in the AID-null controls, per run.

Each point is the mismatch rate among control bases of one reported quality (germline sites excluded; bins of
two quality values with at least 2,000 bases); the diagonal is the rate the quality score promises,
10^(-Q/10). The four loci were sequenced in three runs (J_H4; N114; Smu and Pim1 from PRJNA1481264).
Data: out/intuition_<dataset>.json (intuition_data.py).
"""
import json

import numpy as np
from matplotlib import pyplot as plt
from style import DATA, HALF, INK, MM, SOFT, save, setup

SETS = {"jh4": ("J$_H$4 (this study)", "#274b7e", "o"), "n114": ("N114 (this study)", "#6d9ad0", "s"),
        "su_lps": ("Sμ (PRJNA1481264)", "#1f6f8b", "^"), "pim1": ("$Pim1$ (PRJNA1481264)", "#9a3439", "D")}
BIN, MIN_BASES = 2, 2000


def main():
    setup()
    fig, ax = plt.subplots(figsize=(HALF, 70 * MM))
    fig.subplots_adjust(left=0.17, right=0.97, top=0.9, bottom=0.16)
    out = {}
    for key, (label, col, mk) in SETS.items():
        q = json.load(open(DATA / f"intuition_{key}.json"))["quality"]
        mm, n = {}, {}
        for qq, (m, b) in q.items():
            k = int(qq) // BIN * BIN
            mm[k] = mm.get(k, 0) + m
            n[k] = n.get(k, 0) + b
        ks = [k for k in sorted(n) if n[k] >= MIN_BASES and mm[k]]
        xs, ys = [k + BIN / 2 for k in ks], [mm[k] / n[k] for k in ks]
        ax.plot(xs, ys, marker=mk, ms=2.6, lw=0.8, color=col, label=label)
        tot_m, tot_b = sum(m for qq, (m, b) in q.items() if int(qq) >= 20), sum(b for qq, (m, b) in q.items()
                                                                                 if int(qq) >= 20)
        out[key] = tot_m / tot_b
    qs = np.arange(0, 51)
    ax.plot(qs, 10 ** (-qs / 10), color=INK, lw=0.7, ls="--", label="what the quality score promises")
    ax.set_yscale("log")
    ax.set_xlim(0, 50)
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel("Reported base quality")
    ax.set_ylabel("Mismatch rate in AID-null controls")
    ax.legend(fontsize=5.2, loc="lower left", handlelength=1.6, borderaxespad=0.3)
    ax.grid(True, which="major", color=SOFT, lw=0.2, alpha=0.5)
    save(fig, "FigureS2")
    json.dump({"error_q20plus": out}, open(DATA / "figS8_numbers.json", "w"), indent=1)


if __name__ == "__main__":
    main()
