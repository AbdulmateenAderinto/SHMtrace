#!/usr/bin/env python3
# =============================================================================
# figureS7_recurrence_weighted.py  ->  Supplementary Figure S7
# -----------------------------------------------------------------------------
# Recurrence-weighted N114 V-region substitution spectra (pMX-Aicda vs empty
# vector): each substitution is tallied by the number of supporting reads,
# rather than counted once per unique position. This is the companion to the
# unique-mutation matrices in Figure 4B/C and shows that read-level weighting
# preserves the same G:C-biased, AID-catalysed signature (dominant C->T / G->A).
# Matrices are the audited October C1, D3, D5 sub500 exports (M[to][from]).
# Run:  python figureS7_recurrence_weighted.py
# =============================================================================
import os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from config import OUT_DIR
mpl.rcParams.update({"font.family":"DejaVu Sans","savefig.dpi":600,"font.size":12,
    "figure.facecolor":"white","text.color":"#1a1a1a",
    "mathtext.fontset":"custom","mathtext.rm":"DejaVu Sans","mathtext.bf":"DejaVu Sans:bold",
    "mathtext.it":"DejaVu Sans:italic"})   # regular-weight italic gene name (matches Figure 4)
SEQ=LinearSegmentedColormap.from_list("seqblue",["#F5F8FC","#C4D8EE","#6FA0CE","#2F6BA6","#173D63"])
BASES=["A","T","C","G"]
def m(r): return np.array(r,float)
# recurrence-weighted matrices (each unique mutation x supporting reads)
aid_w=m([[0,2,120,577],[2,0,470,72],[1,0,0,101],[1,1,128,0]])  # n = 1,475
emp_w=m([[0,0,8,0],[0,0,0,9],[0,0,0,0],[2,0,0,0]])             # n = 19

def heat(ax,M,label,n,gvmax,show_from=True):
    MT=M.T; P=100*MT/MT.sum(); rowpct=100*MT.sum(1)/MT.sum()
    for i in range(4):
        for j in range(4):
            x,y=j,3-i
            if i==j:
                ax.add_patch(Rectangle((x,y),1,1,fc="white",ec="white"))
                ax.text(x+.5,y+.5,"—",ha="center",va="center",color="#AFAFAF",fontsize=14); continue
            v=P[i,j]; ax.add_patch(Rectangle((x,y),1,1,fc=SEQ(v/gvmax),ec="white",lw=2))
            tc="white" if v/gvmax>0.55 else "#1a1a1a"
            ax.text(x+.5,y+.60,f"{int(MT[i,j])}",ha="center",va="center",color=tc,fontsize=12,fontweight="bold")
            ax.text(x+.5,y+.30,f"{v:.1f}%",ha="center",va="center",color=tc,fontsize=7.5)
    for i in range(4):
        ax.add_patch(Rectangle((5,3-i),1,1,fc=SEQ(rowpct[i]/max(rowpct)*.6),ec="white",lw=2))
        ax.text(5.5,3-i+.5,f"{rowpct[i]:.1f}",ha="center",va="center",fontsize=10.5,fontweight="bold")
    for j,b in enumerate(BASES): ax.text(j+.5,4.16,b,ha="center",va="bottom",fontsize=11,fontweight="bold")
    ax.text(5.5,4.16,"%",ha="center",va="bottom",fontsize=11,fontweight="bold")
    ax.text(1.5,4.72,"To",ha="center",va="bottom",fontsize=11,fontweight="bold")
    for i,b in enumerate(BASES): ax.text(-0.55,3-i+.5,b,ha="center",va="center",fontsize=11,fontweight="bold")
    # "From" only on the left matrix; dropping it on the right closes the gap
    if show_from:
        ax.text(-1.3,2,"From",ha="center",va="center",fontsize=11,fontweight="bold",rotation=90)
    lo=-1.6 if show_from else -0.95
    ax.text((lo+6.4)/2,-0.7,f"{label} (n = {n:,})",ha="center",va="top",fontsize=11.5)
    ax.set_xlim(lo,6.4); ax.set_ylim(-1.15,4.75); ax.set_aspect("equal"); ax.axis("off")

# ---- canvas: matrices tightened together with the colour bar as a single row
#      directly beneath BOTH of them (matches Figures 3A and 4B/C) -------------
gv=max(100*aid_w.T.max()/aid_w.T.sum(),100*emp_w.T.max()/emp_w.T.sum())
fig=plt.figure(figsize=(10.4,5.0))
# B's x-range is narrower (no 'From' label), so give it a proportionally
# narrower cell to keep both matrices at exactly the same scale.
gs=fig.add_gridspec(2,2,height_ratios=[18,1.0],hspace=0.02,wspace=0.02,
                    width_ratios=[8.0,7.35],
                    left=0.045,right=0.985,top=0.92,bottom=0.135)
axA=fig.add_subplot(gs[0,0]); heat(axA,aid_w,"pMX-$\\it{Aicda}$",1475,gv)
axB=fig.add_subplot(gs[0,1]); heat(axB,emp_w,"pMX-Empty vector",19,gv,show_from=False)
sm=plt.cm.ScalarMappable(cmap=SEQ,norm=plt.Normalize(0,gv))
cax=fig.add_subplot(gs[1,:])
cb=fig.colorbar(sm,cax=cax,orientation="horizontal")
cb.ax.tick_params(labelsize=7.5)   # unit is stated in the note below
# Descriptive note directly under the colour bar.
_pA=axA.get_position(); _pB=axB.get_position()
fig.text(_pA.x0+(_pB.x1-_pA.x0)*0.5,0.065,
         "Read-weighted substitutions: each unique mutation tallied by its number of "
         "supporting reads; cell = count and % of total.",
         fontsize=8.2,va="top",ha="center")
# Panel letters on ONE shared baseline (figure coords) so A and B line up
_ylet=max(axA.get_position().y1,axB.get_position().y1)+0.012
for _ax,_lab,_dx in [(axA,"A",-0.012),(axB,"B",-0.012)]:
    fig.text(_ax.get_position().x0+_dx,_ylet,_lab,fontsize=20,fontweight="bold",va="top",ha="left")
OUT=OUT_DIR+os.sep
for ext in ("png","pdf"):
    fig.savefig(OUT+f"FigureS7_recurrence_weighted.{ext}",dpi=600,bbox_inches="tight",pad_inches=.04,facecolor="white")
# Renderers (LibreOffice/Word) shave a few pixels off the top of a wide inline
# image, which was clipping the bold panel letters. Bake a small white margin
# into the raster so the letters can never be cropped.
from PIL import Image as _Image
_p=OUT+"FigureS7_recurrence_weighted.png"
_im=_Image.open(_p).convert("RGB"); _w,_h=_im.size; _strip=int(_h*0.10)
_c=_Image.new("RGB",(_w,_h+_strip),"white"); _c.paste(_im,(0,_strip)); _c.save(_p)
print("top guard added: %dpx"%_strip)
print("S7 recurrence-weighted written (Aicda n=1,475 | Empty n=19)")
