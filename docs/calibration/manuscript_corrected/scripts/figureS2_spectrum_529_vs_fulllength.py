#!/usr/bin/env python3
# =============================================================================
# figureS2_spectrum_529_vs_fulllength.py  ->  Supplementary Figure S2
# -----------------------------------------------------------------------------
# WT ONT JH4 spectrum, 529 bp window (n=571) vs reliable 1,141-bp interval (n=739).
# Panels A/B: From->To matrices with shared colour scale + shared colorbar
# (styled to match Figure 3A). Panel C: agreement scatter (Pearson r).
# Matrices are the July-13 phased sub500 values, inlined for portability.
# Run:  python figureS2_spectrum_529_vs_fulllength.py
# =============================================================================
import os
from config import OUT_DIR
# -*- coding: utf-8 -*-
import numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":12,"savefig.dpi":600,"figure.facecolor":"white"})
SEQ=LinearSegmentedColormap.from_list("seqblue",["#F5F8FC","#C4D8EE","#6FA0CE","#2F6BA6","#173D63"])
TEAL="#CC79A7"; BLUE="#0072B2"  # ONT accent = reddish-purple; stat text = blue
BASES=["A","T","C","G"]; IDX={b:i for i,b in enumerate(BASES)}
MUT=["A>T","A>C","A>G","T>A","T>C","T>G","C>A","C>T","C>G","G>A","G>T","G>C"]
def c(M,f,t): return M[IDX[t]][IDX[f]]
# SHMtrace calls, Wt mice pooled: the 529-bp window (jh4_shmtrace.py) and positions
# 90-1,230 (jh4_shmtrace.py --full-length); previously Geneious July-13 values
import json
from pathlib import Path
_CAL=Path(__file__).resolve().parents[2]
def _pooled(run):
    M=np.zeros((4,4))
    for k,v in json.load(open(_CAL/f"figures/{run}/project.json"))["spectrum_group"]["Wt"].items():
        M[IDX[k[2]]][IDX[k[0]]]=v
    return M
win=_pooled("jh4_shmtrace"); full=_pooled("jh4_shmtrace_full")
NW,NF=int(win.sum()),int(full.sum())
# shared vmax across both matrices (as in Figure 3A)
MTw=win.T; MTf=full.T
gvmax=max(100*MTw.max()/MTw.sum(),100*MTf.max()/MTf.sum())

# Reviewer request: place panels A, B, C beside each other in a single row
# (each panel smaller) rather than stacking the scatter beneath the matrices.
fig=plt.figure(figsize=(14.8,5.4))
outer=fig.add_gridspec(1,2,width_ratios=[2.0,1.03],wspace=0.16,
                       left=0.035,right=0.985,top=0.92,bottom=0.13)
# A and B share a block so the colour bar can sit as a tight row directly
# beneath BOTH matrices; the narrower B cell keeps the two at identical scale.
gsAB=outer[0,0].subgridspec(2,2,height_ratios=[18,1.0],hspace=0.02,wspace=0.02,
                            width_ratios=[8.1,7.4])
def heat(ax,M,title,show_from=True):
    MT=M.T; P=100*MT/MT.sum(); rowpct=100*MT.sum(1)/MT.sum()
    for i in range(4):
        for j in range(4):
            x,y=j,3-i
            if i==j: ax.add_patch(Rectangle((x,y),1,1,fc="white",ec="white")); ax.text(x+.5,y+.5,"—",ha="center",va="center",color="#BBB",fontsize=13); continue
            v=P[i,j]; ax.add_patch(Rectangle((x,y),1,1,fc=SEQ(v/gvmax),ec="white",lw=2))
            tc="white" if v/gvmax>0.55 else "#1a1a1a"
            ax.text(x+.5,y+.60,f"{int(MT[i,j])}",ha="center",va="center",color=tc,fontsize=13,fontweight="bold")
            ax.text(x+.5,y+.30,f"{v:.1f}%",ha="center",va="center",color=tc,fontsize=8)
    for i in range(4): ax.add_patch(Rectangle((5,3-i),1,1,fc=SEQ(rowpct[i]/max(rowpct)*.6),ec="white",lw=2)); ax.text(5.5,3-i+.5,f"{rowpct[i]:.1f}",ha="center",va="center",fontsize=11,fontweight="bold")
    for j,b in enumerate(BASES): ax.text(j+.5,4.16,b,ha="center",fontsize=12,fontweight="bold")
    ax.text(5.5,4.16,"%",ha="center",fontsize=12,fontweight="bold")
    # "To" and the panel title are both va="bottom" so the gap between them is
    # measured from the glyph box, not the baseline. With the default (baseline)
    # anchoring the title's parenthesis descenders ran into the "To" cap height.
    ax.text(1.5,4.52,"To",ha="center",va="bottom",fontsize=11,fontweight="bold")
    for i,b in enumerate(BASES): ax.text(-0.55,3-i+.5,b,ha="center",va="center",fontsize=12,fontweight="bold")
    if show_from:
        ax.text(-1.35,2,"From",ha="center",va="center",fontsize=12,fontweight="bold",rotation=90)
    lo=-1.7 if show_from else -1.0
    ax.text((lo+6.4)/2,5.06,title,ha="center",va="bottom",fontsize=13,fontweight="bold")
    ax.set_xlim(lo,6.4); ax.set_ylim(-0.5,5.55); ax.set_aspect("equal"); ax.axis("off")
axA=fig.add_subplot(gsAB[0,0]); heat(axA,win,f"529 bp window (n = {NW:,})")
axB=fig.add_subplot(gsAB[0,1]); heat(axB,full,f"Reliable 1,141-bp interval (n = {NF:,})",show_from=False)
# colour bar as a tight row spanning beneath both matrices
cax=fig.add_subplot(gsAB[1,:])
sm=plt.cm.ScalarMappable(cmap=SEQ,norm=plt.Normalize(0,gvmax))
cb=fig.colorbar(sm,cax=cax,orientation="horizontal"); cb.set_label("% of total mutations",fontsize=9); cb.ax.tick_params(labelsize=8)

# ---- third column: agreement scatter C ----
axC=fig.add_subplot(outer[0,1])
def props(M): tot=M.sum(); return np.array([100*c(M,*mt.split(">"))/tot for mt in MUT])
wp=props(win); fp=props(full); r=np.corrcoef(wp,fp)[0,1]
lim=max(wp.max(),fp.max())*1.12
axC.plot([0,lim],[0,lim],ls="--",color="#999",lw=1.3)
axC.scatter(wp,fp,s=85,facecolor=TEAL,edgecolor="white",lw=0.9,zorder=3)
# Label all 12 substitution types; manual offsets and leader lines keep the
# tightly clustered low-frequency points legible at manuscript scale.
label_offsets={
    "A>T":(8, 8), "A>C":(-32, 8), "A>G":(8, 7),
    "T>A":(-31,-11), "T>C":(8,-8), "T>G":(8, 8),
    "C>A":(-30,-12), "C>T":(-32, 8), "C>G":(-30, 8),
    "G>A":(8,-12), "G>T":(8,-11), "G>C":(8,-10),
}
for i,t in enumerate(MUT):
    axC.annotate(t.replace(">","→"),(wp[i],fp[i]),fontsize=8.5,
                 xytext=label_offsets[t],textcoords="offset points",
                 fontweight="bold",color="#333",
                 arrowprops=dict(arrowstyle="-",color="#777",lw=.55,
                                 shrinkA=1,shrinkB=4))
axC.set_xlim(0,lim); axC.set_ylim(0,lim); axC.set_aspect("equal")
axC.set_xlabel("529 bp window (% of total mutations)",fontsize=11); axC.set_ylabel("Reliable 1,141-bp interval (% of total mutations)",fontsize=11)
axC.text(lim*0.97,lim*0.05,f"Pearson r = {r:.3f}",fontsize=12.5,color=BLUE,fontweight="bold",ha="right",va="bottom")
for s in ["top","right"]: axC.spines[s].set_visible(False)
# Panel letters on ONE shared baseline (figure coords) so A, B and C line up
_tops=[axA.get_position().y1,axB.get_position().y1,axC.get_position().y1]
_ylet=max(_tops)+0.015
for _ax,_lab,_dx in [(axA,"A",-0.030),(axB,"B",-0.020),(axC,"C",-0.055)]:
    fig.text(_ax.get_position().x0+_dx,_ylet,_lab,fontsize=22,fontweight="bold",va="top",ha="left")
# ---- verify the title clears the "To" label on both panels -------------------
fig.canvas.draw()
_r=fig.canvas.get_renderer()
for _ax,_name in [(axA,"A"),(axB,"B")]:
    _t=[t for t in _ax.texts if t.get_text()=="To"][0]
    _ti=[t for t in _ax.texts if t.get_text().startswith(("529","Reliable"))][0]
    _gap=_ti.get_window_extent(_r).y0-_t.get_window_extent(_r).y1
    _pt=_gap*72/fig.dpi
    print(f"panel {_name}: title/To gap = {_pt:.1f} pt", "OK" if _pt>=4 else "TOO TIGHT")

OUT=OUT_DIR+os.sep
for ext in ("png","pdf"):
    fig.savefig(OUT+f"FigureS_spectrum_529_vs_fulllength.{ext}",dpi=600,
                bbox_inches="tight",pad_inches=.04,facecolor="white")
from PIL import Image as _Image
_p=OUT+"FigureS_spectrum_529_vs_fulllength.png"
_im=_Image.open(_p).convert("RGB"); _w,_h=_im.size; _strip=int(_h*0.12)
_c=_Image.new("RGB",(_w,_h+_strip),"white"); _c.paste(_im,(0,_strip)); _c.save(_p)
print("top guard added: %dpx"%_strip)
print(f"S2 saved r={r:.3f} n {NW} {NF} maxdiff {abs(wp-fp).max():.2f}")
