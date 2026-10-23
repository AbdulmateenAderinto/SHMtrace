"""Explicit, reproducible edits to copied figure builders only."""
from pathlib import Path
P=Path(__file__).parent

def change(name,fn):
 p=P/name;s=p.read_text();p.write_text(fn(s))

def fig2(s):
 start=s.index('    # the range in which the somatic callers')
 end=s.index('    ax.axvline(50',start)
 s=s[:start]+s[end:]
 s=s.replace('    from fig3_calibration import by_fraction\n','')
 s=s.replace('Control FDR,','Empirical threshold,')
 s=s.replace('ax.legend(fontsize=6, loc="upper left", handlelength=1.4, borderaxespad=0.2)', 'ax.legend(fontsize=6.5, loc="upper left", handlelength=1.3, borderaxespad=0.2)')
 return s
change('fig2_call.py',fig2)

def fig5(s):
 start=s.index('    # C: the same comparison')
 end=s.index('    save(fig, "Figure5")',start)
 s=s[:start]+'''    # C: paired biological observations; one marker shape per mouse.
    ax = fig.add_subplot(gs[0, 2:])
    pm = json.load(open(DATA / "fig5_per_mouse.json"))
    marks = {"WT2": "o", "WT3": "s", "WT4": "^"}
    xs = np.arange(len(pm["classes"]))
    for i, m in enumerate(pm["mice"]):
        off = (i - 1) * 0.08
        for k in xs:
            va = pm["share"]["sanger"][m][k]
            vb = pm["share"]["ont"][m][k]
            ax.plot([k - 0.18 + off, k + 0.18 + off], [va, vb], color=SOFT, lw=0.6, alpha=0.65)
        for j, (plat, col) in enumerate((("sanger", SANG_C), ("ont", ONT_C))):
            ax.scatter(xs + (j - 0.5) * 0.36 + off,
                pm["share"][plat][m], s=19, marker=marks[m], facecolor=col,
                edgecolor="white", lw=0.5, zorder=3)
    ax.set_xticks(xs, [c.replace(" ", "\\n", 1) for c in pm["classes"]])
    ax.set_xlim(-0.6, len(xs)-0.4); ax.set_ylim(0, 62)
    ax.set_ylabel("Share of each mouse's calls (%)")
    ax.set_title("Paired spectra in the same three mice", fontsize=7)
    from matplotlib.lines import Line2D
    hs=[Line2D([],[],marker="o",ls="",color=SANG_C,label="Sanger"),
        Line2D([],[],marker="o",ls="",color=ONT_C,label="Nanopore")]
    hs += [Line2D([],[],marker=mk,ls="",color=INK,label=f"Wt {m[-1]}") for m,mk in marks.items()]
    ax.legend(handles=hs,fontsize=6.5,ncol=3,loc="upper left",columnspacing=0.8,handletextpad=0.3)
    letter(ax,"C",dx=-0.1)
''' +s[end:]
 s=s.replace('color=SANG_C, lw=1.1, marker="o", ms=2.2','color=SANG_C, lw=0.8, marker="o", ms=1.8')
 s=s.replace('color=ONT_C, lw=1.1, marker="o", ms=2.2','color=ONT_C, lw=0.8, marker="o", ms=1.8')
 s=s.replace('label="sites"','label="motif availability"')
 s=s.replace('ax.bar(x, np.mean(vals), 0.8, color=col, alpha=0.85)','ax.bar(x, np.mean(vals), 0.8, color=col, alpha=0.18)')
 s=s.replace('s=8, color=INK','s=14, color=col')
 return s
change('fig5_sanger.py',fig5)

def fig6(s):
 s=s.replace('color="#6d9ad0", lw=0.6','color="#6d9ad0", lw=0.6, alpha=0.4')
 s=s.replace('"Per position, relative to mean"','"Events/position ÷ mean across quintiles"')
 return s
change('fig6_longreads.py',fig6)
# Keep an obvious genotype colour distinction without red-green pairing.
change('figS14_sanger_cohort2.py',lambda s:s.replace('"#4f7a3a", "#9cc28a"','"#6a4c93", "#b9a3d6"'))
# No arbitrary positive floor for a genuinely zero observed FDP.
change('figS4_fdr_calibration.py',lambda s:s.replace('max(v["fdp_pooled"], ZERO)','v["fdp_pooled"]')
 .replace('ax.set_yscale("log"); ax.set_ylim(0.0022, 0.6)',
 'ax.set_yscale("symlog", linthresh=0.01, linscale=0.5); ax.set_ylim(-0.001, 0.6)')
 .replace('ax.set_yticks([ZERO] + yt','ax.set_yticks([0] + yt'))
