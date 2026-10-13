#!/usr/bin/env python3
# =============================================================================
# figure4_hybridoma_combined.py  ->  Figure 4 (combines former Fig 4 + Fig 5)
# -----------------------------------------------------------------------------
# N114 hybridoma V-region SHM, pMX-Aicda vs pMX-Empty vector, by ONT.
#   A  unique-variant density (dot plot, geometric mean, n=3 subclones, Welch t)
#   B  pMX-Aicda pooled UNIQUE-mutation spectrum matrix
#   C  pMX-Empty pooled UNIQUE-mutation spectrum matrix
# Square 2x2 layout, panel letters only (no descriptive titles on top of panels),
# minimal white space. The recurrence-weighted matrices (old Fig 5B) move to the
# supplement. Values: SHMtrace re-analysis of the deposited N114 reads (500 reads per
# subclone); previously the Geneious sub500 exports.
# Run:  python figure4_hybridoma_combined.py
# =============================================================================
import os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats as _st
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from config import OUT_DIR
mpl.rcParams.update({"font.family":"DejaVu Sans","savefig.dpi":600,"font.size":12,
    "axes.edgecolor":"black","xtick.color":"black","ytick.color":"black",
    "text.color":"#1a1a1a","axes.labelcolor":"#1a1a1a",
    "mathtext.fontset":"custom","mathtext.rm":"DejaVu Sans","mathtext.bf":"DejaVu Sans:bold",
    "mathtext.it":"DejaVu Sans:italic"})   # regular-weight italic gene name (consistent with surrounding text)
BLUE=("#6D9AD0","#274B7E"); RED=("#D07A7E","#9A3439"); GRID="#C9C9C9"
SEQ=LinearSegmentedColormap.from_list("seqblue",["#F5F8FC","#C4D8EE","#6FA0CE","#2F6BA6","#173D63"])
BASES=["A","T","C","G"]

# ---- data: the SHMtrace re-analysis of the deposited N114 reads (500 reads per
# subclone, recommended quality-calibrated mode; docs/calibration/n114_shmtrace.py) ----
import json as _json
from pathlib import Path as _Path
_RES = _json.load(open(_Path(__file__).resolve().parents[2] / "figures/n114_shmtrace/project.json"))
_S = {s["sample"]: s for s in _RES["samples"]}
AID_IDS, EMP_IDS = ["C1", "D3", "D5"], ["C3", "C4", "D1"]
AID_N = [_S[k]["unique_mutations"] for k in AID_IDS]
EMP_N = [_S[k]["unique_mutations"] for k in EMP_IDS]
# a subclone with no calls is drawn (and tested) at its one-event density: the left-
# censoring convention, so the fold is a lower bound
def _dens(ids):
    out, cens = [], []
    for k in ids:
        d, one = _S[k]["unique_density_x1e5"], _S[k]["one_event_density_x1e5"]
        out.append(d if d > 0 else one)
        cens.append(d <= 0)
    return out, cens
aicda, aicda_cens = _dens(AID_IDS)
empty, empty_cens = _dens(EMP_IDS)
def _matrix(group):
    c = _RES["spectrum_group"].get(group, {})
    M = np.zeros((4, 4))
    for f in BASES:
        for t in BASES:
            if f != t:
                M[BASES.index(t)][BASES.index(f)] = c.get(f"{f}>{t}", 0)
    return M
aid_u = _matrix("pMX-Aicda")     # M[to][from]
emp_u = _matrix("pMX-Empty")
N_AID, N_EMP = int(aid_u.sum()), int(emp_u.sum())

def geo(x):
    x=np.asarray(x,float); return 10**np.mean(np.log10(x)),10**np.std(np.log10(x),ddof=1)
# statistics computed from the data, never typed in by hand
PVAL=f"{_st.ttest_ind(np.log10(aicda),np.log10(empty),equal_var=False).pvalue:.4f}"
FOLD=round(geo(aicda)[0]/geo(empty)[0])
FOLD_TXT=("≥ " if any(empty_cens) else "~")+f"{FOLD}×"

def dotplot(ax):
    rng=np.random.default_rng(4)
    for i,(d,cens,(fill,edge)) in enumerate([(aicda,aicda_cens,BLUE),(empty,empty_cens,RED)]):
        gm,gs=geo(d)
        ax.plot([i,i],[gm/gs,gm*gs],color=edge,lw=1.6,zorder=3)
        for yy in (gm/gs,gm*gs): ax.plot([i-.07,i+.07],[yy,yy],color=edge,lw=1.6,zorder=3)
        ax.plot([i-.26,i+.26],[gm,gm],color=edge,lw=3.0,zorder=5,solid_capstyle="round")
        xs=i+(rng.random(len(d))-.5)*0.16
        ax.scatter(xs,d,facecolors=["white" if c else fill for c in cens],edgecolors=edge,s=62,linewidths=1.5,zorder=6)
    ax.set_yscale("log"); ax.set_ylim(0.1,150); ax.set_xlim(-0.7,1.7)
    if all(empty_cens):
        ax.text(1,min(empty)*0.62,"no calls\n(one-event scale)",ha="center",va="top",fontsize=7.5,color="#7A2E33")
    ax.set_xticks([0,1]); ax.set_xticklabels(["pMX-$\\it{Aicda}$","pMX-Empty vector"],rotation=25,ha="right",fontsize=9)
    ax.set_ylabel("Mutation density, N114 V region\n(×10$^{-5}$ mut/bp; 615 bp)",fontsize=9)
    ax.set_axisbelow(True); ax.yaxis.grid(True,linestyle=(0,(5,4)),color=GRID,linewidth=1.1,zorder=0)
    for sp in ("left","bottom"): ax.spines[sp].set_linewidth(1.9)
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(width=1.5,length=5,labelsize=9)
    u1=geo(aicda); u2=geo(empty); upper=max(u1[0]*u1[1],u2[0]*u2[1],max(aicda),max(empty)); yb=upper*1.45
    ax.plot([0,0,1,1],[yb,yb*1.4,yb*1.4,yb],color="black",lw=1.5,zorder=8)
    ax.text(.5,yb*1.75,f"$P$ = {PVAL} ({FOLD_TXT})",ha="center",fontsize=9.5)

def heat(ax,M,label,n,gvmax,show_from=True):
    MT=M.T; tot=MT.sum() or 1.0      # a group with no calls: every cell 0
    P=100*MT/tot; rowpct=100*MT.sum(1)/tot
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
        ax.add_patch(Rectangle((5,3-i),1,1,fc=SEQ(rowpct[i]/(max(rowpct) or 1)*.6),ec="white",lw=2))
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

# ---- canvas (reviewer request): all three panels on ONE level, A beside B and C
#      and made smaller; colour bar + note beneath the two matrices -------------
gv=max(100*aid_u.T.max()/(aid_u.T.sum() or 1),100*emp_u.T.max()/(emp_u.T.sum() or 1))
fig=plt.figure(figsize=(13.6,5.2))
outer=fig.add_gridspec(1,2,width_ratios=[0.52,2.0],wspace=0.13,
                       left=0.05,right=0.985,top=0.93,bottom=0.115)
# A: dot plot, narrower than the matrices so it does not dominate the row
axA=fig.add_subplot(outer[0,0]); dotplot(axA)
# B, C: matrices sharing one block, with the colour bar as a tight row directly
# beneath BOTH of them (mirrors the Figure 3A treatment)
# C's x-range is narrower (no 'From' label), so give it a proportionally
# narrower cell to keep both matrices at exactly the same scale.
gsBC=outer[0,1].subgridspec(2,2,height_ratios=[18,1.0],hspace=0.02,wspace=0.02,
                            width_ratios=[8.0,7.35])
axB=fig.add_subplot(gsBC[0,0]); heat(axB,aid_u,"pMX-$\\it{Aicda}$",N_AID,gv)
axC=fig.add_subplot(gsBC[0,1]); heat(axC,emp_u,"pMX-Empty vector",N_EMP,gv,show_from=False)
sm=plt.cm.ScalarMappable(cmap=SEQ,norm=plt.Normalize(0,gv))
cax=fig.add_subplot(gsBC[1,:])
cb=fig.colorbar(sm,cax=cax,orientation="horizontal")
cb.ax.tick_params(labelsize=7.5)   # unit is stated in the note below
# Descriptive note directly under the colour bar.
_pB=axB.get_position(); _pC=axC.get_position()
fig.text(_pB.x0+(_pC.x1-_pB.x0)*0.5,0.055,
         "Pooled unique substitutions across three pMX-$\\it{Aicda}$ subclones (C1, D3, D5) "
         "vs three empty-vector controls; cell = count and % of total.",
         fontsize=8.2,va="top",ha="center")
# panel letters
# One shared baseline (figure coords) keeps A, B and C exactly in line
_ylet=max(axA.get_position().y1,axB.get_position().y1,axC.get_position().y1)+0.012
for _ax,_lab,_dx in [(axA,"A",-0.038),(axB,"B",-0.012),(axC,"C",-0.012)]:
    fig.text(_ax.get_position().x0+_dx,_ylet,_lab,fontsize=20,fontweight="bold",va="top",ha="left")
OUT=OUT_DIR+os.sep
for ext in ("png","pdf"):
    fig.savefig(OUT+f"Figure4_hybridoma_combined.{ext}",dpi=600,bbox_inches="tight",pad_inches=.04,facecolor="white")
print("pMX-Aicda GM %.1f | Empty GM %.2f%s | fold %s | P %s | calls %s / %s"%(geo(aicda)[0],geo(empty)[0],
      " (all at the one-event density)" if all(empty_cens) else "",FOLD_TXT,PVAL,AID_N,EMP_N))
# Renderers (LibreOffice/Word) shave a few pixels off the top of a wide inline
# image, which was clipping the bold panel letters. Bake a small white margin
# into the raster so the letters can never be cropped.
from PIL import Image as _Image
_p=OUT+"Figure4_hybridoma_combined.png"
_im=_Image.open(_p).convert("RGB"); _w,_h=_im.size; _strip=int(_h*0.10)
_c=_Image.new("RGB",(_w,_h+_strip),"white"); _c.paste(_im,(0,_strip)); _c.save(_p)
print("top guard added: %dpx"%_strip)
print("Figure 4 (combined) written to",OUT)
