#!/usr/bin/env python3
# =============================================================================
# figure2_jh4_combined.py  ->  Figure 2 (combines former Fig 2 + Suppl. Fig S1A)
# -----------------------------------------------------------------------------
# JH4-intron SHM, wild-type vs Aicda-/-, promoted to a single main figure:
#   A  ONT   unique-variant density (529 bp window, WT vs Aicda-/-)
#   B  Sanger unique-variant density (same mice)
#   C  full-length per-position SHM density + read coverage (WT; 90-1,230 bp)
# Panel C is the former Supplementary Figure S1A, moved up per reviewer request.
# Panel letters only (no descriptive titles on top of panels); minimal white
# space; A and B share one y-axis so the platforms are directly comparable.
# Run:  python figure2_jh4_combined.py
# =============================================================================
import os, csv, glob, statistics as st
from scipy import stats as _st
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
import pysam
from config import OUT_DIR, DATA_JH4
mpl.rcParams.update({"font.family":"DejaVu Sans","savefig.dpi":600,"font.size":12,
    "axes.edgecolor":"black","xtick.color":"black","ytick.color":"black",
    "text.color":"#1a1a1a","axes.labelcolor":"#1a1a1a"})
BLUE=("#6D9AD0","#274B7E"); RED=("#D07A7E","#9A3439"); GRID="#C9C9C9"
BLUEL="#2C6FB0"; GREY="#9AA3AD"; SHADE="#DCE7F2"

# ---- panels A/B: per-mouse unique-variant density (x10^-5 / bp), n=3 --------
# Densities are derived from the audited unique-variant counts rather than being
# typed in pre-rounded: rounding the KO values to one decimal shifted the plotted
# points and changed the recomputed fold-change from 47.0 to 46.4.
#   ONT sub500, 529 bp window: Wt 161/342/68, Aicda-/- 2/6/3 unique variants
#   Sanger, 529 bp window:     Wt 189/227/99 over 37/36/38 independent clones;
#                              Aicda-/- 0/1/1 over 10/9/9 clones (a zero-mutation
#                              sample is placed at its own 1-mutation floor).
def _dens(counts, denoms, bp=529): return [c/(d*bp)*1e5 for c, d in zip(counts, denoms)]
ont ={"WT":_dens([161,342,68],[500]*3), "KO":_dens([2,6,3],[500]*3)}
# CORRECTED SANGER (nanoshm, organised .ab1 traces; see SANGER_VS_NANOPORE.md):
# unique mutations over the bases each clone read, expressed as full 529-bp clones.
import json as _json
from config import DATA_SANGER as _DS
_cs=_json.load(open(os.path.join(_DS,"corrected_sanger.json")))
def _c(k): return _cs[k]["unique"], _cs[k]["bases"]/529
_wt=[_c(k) for k in ("Wt 2","Wt 3","Wt 4")]; _ko=[_c(k) for k in ("Aicda−/− 1","Aicda−/− 2","Aicda−/− 3")]
# a zero-mutation sample sits at its own one-mutation floor, as before
sang={"WT":_dens([max(u,1) for u,_ in _wt],[n for _,n in _wt]), "KO":_dens([max(u,1) for u,_ in _ko],[n for _,n in _ko])}
YLIM=(0.4,8000)
YLAB=("Mutation density, J$_H$4 region\nof GC PP B cells (×10$^{-5}$ mut/bp)")

def geo(x):
    x=np.asarray(x,float); return 10**np.mean(np.log10(x)),10**np.std(np.log10(x),ddof=1)
def welch(a,b):
    """Two-tailed Welch t-test on log10 densities; returns (p-string, fold)."""
    p=_st.ttest_ind(np.log10(a),np.log10(b),equal_var=False).pvalue
    return f"{p:.4f}", round(geo(a)[0]/geo(b)[0])
def fstyle(ax):
    ax.set_axisbelow(True); ax.yaxis.grid(True,linestyle=(0,(5,4)),color=GRID,linewidth=1.1,zorder=0)
    for sp in ("left","bottom"): ax.spines[sp].set_linewidth(1.9)
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(width=1.5,length=5,labelsize=9)
def freqpanel(ax,g1,g2,pval,fold,ylab=None):
    rng=np.random.default_rng(4)
    for i,(d,(fill,edge)) in enumerate([(g1,BLUE),(g2,RED)]):
        gm,gs=geo(d)
        ax.plot([i,i],[gm/gs,gm*gs],color=edge,lw=1.6,zorder=3)
        for yy in (gm/gs,gm*gs): ax.plot([i-.07,i+.07],[yy,yy],color=edge,lw=1.6,zorder=3)
        ax.plot([i-.26,i+.26],[gm,gm],color=edge,lw=3.0,zorder=5,solid_capstyle="round")
        xs=i+(rng.random(len(d))-.5)*0.16
        ax.scatter(xs,d,facecolors=fill,edgecolors=edge,s=58,linewidths=1.5,zorder=6)
    ax.set_yscale("log"); ax.set_ylim(*YLIM); ax.set_xlim(-0.7,1.7)
    ax.set_xticks([0,1]); ax.set_xticklabels(["$\\it{Wt}$","$\\it{Aicda}^{-/-}$"],rotation=25,ha="right",fontsize=9)
    fstyle(ax)
    if ylab: ax.set_ylabel(ylab,fontsize=9)
    u1=geo(g1); u2=geo(g2); upper=max(u1[0]*u1[1],u2[0]*u2[1],max(g1),max(g2)); yb=upper*1.75
    ax.plot([0,0,1,1],[yb,yb*1.28,yb*1.28,yb],color="black",lw=1.5,zorder=8)
    ax.text(.5,yb*1.45,f"$P$ = {pval} (~{fold}×)",ha="center",fontsize=9)

# ---- panel C: full-length per-position density + coverage (WT) -------------
BASES=set("ACGT"); REFLEN=1230; LO,HI=90,1230; DIR=DATA_JH4+os.sep
WT=[("WT2","DHJGKB_2_Wildtype_2_C57BL6_sub500"),
    ("WT3","DHJGKB_3_Wildtype_3_C57BL6_sub500"),
    ("WT4","DHJGKB_4_Wildtype_4_C57BL6_sub500")]
def muts(stem):
    rows=list(csv.DictReader(open(glob.glob(DIR+stem+"*Annotations.csv")[0],encoding='latin-1'))); out=[]
    for r in rows:
        if 'SNP' not in r.get('Polymorphism Type',''): continue
        ref=(r.get('Reference Nucleotide(s)') or '').strip().upper(); var=(r.get('Variant Nucleotide(s)') or '').strip().upper()
        if not(len(ref)==1 and len(var)==1 and ref in BASES and var in BASES and ref!=var): continue
        try: vf=float((r.get('Variant Frequency','') or '').replace('%','').strip())
        except: vf=0
        if vf>=50: continue
        try: pv=float((r.get('Variant P-Value (approximate)','') or '').strip())
        except: pv=1
        if pv>0.05: continue
        p=int(float(r['Minimum']))
        if LO<=p<=HI: out.append(p)
    return out
def depth(stem):
    bam=pysam.AlignmentFile(glob.glob(DIR+stem+"*.bam")[0],"rb"); d=[0]*(REFLEN+2)
    for col in bam.pileup(truncate=True,min_base_quality=0):
        pp=col.reference_pos+1
        if 1<=pp<=REFLEN: d[pp]=col.nsegments
    bam.close(); return d
BIN=30; nb=(HI-LO)//BIN+1; bcent=[LO+i*BIN+BIN/2 for i in range(nb)]
dens=[]; cov=[]
for _,stem in WT:
    pos=muts(stem); dp=depth(stem)
    mb=[0]*nb; cb=[[] for _ in range(nb)]
    for pp in pos: mb[min((pp-LO)//BIN,nb-1)]+=1
    for pp in range(LO,HI+1): cb[min((pp-LO)//BIN,nb-1)].append(dp[pp])
    cm=[sum(c)/len(c) if c else 0 for c in cb]
    dens.append([mb[i]/(cm[i]*BIN)*1e5 if cm[i]>0 else 0 for i in range(nb)]); cov.append(cm)
meandens=[st.mean(d[i] for d in dens) for i in range(nb)]
meancov=[st.mean(c[i] for c in cov) for i in range(nb)]

# ---- layout (reviewer request): A, B and C all on ONE level, with the two dot
#      plots kept narrow so the full-length map C has room ---------------------
fig=plt.figure(figsize=(14.4,4.6))
# a narrow spacer column keeps C's y-axis label clear of panel B
gs=fig.add_gridspec(1,4,width_ratios=[0.46,0.46,0.22,1.75],wspace=0.10,
                    left=0.055,right=0.945,top=0.905,bottom=0.20)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1],sharey=axA)
# statistics are computed from the data, never typed in by hand
freqpanel(axA,ont["WT"],ont["KO"],*welch(ont["WT"],ont["KO"]),ylab=YLAB)
freqpanel(axB,sang["WT"],sang["KO"],*welch(sang["WT"],sang["KO"]))
axB.tick_params(axis="y",labelleft=False)
axA.set_xlabel("Oxford Nanopore (ONT)",fontsize=9.5,labelpad=6)
axB.set_xlabel("Sanger",fontsize=9.5,labelpad=6)

axC=fig.add_subplot(gs[0,3]); axC.axvspan(LO,540,color=SHADE,zorder=0)
for d in dens: axC.plot(bcent,d,color=BLUEL,lw=0.9,alpha=0.35)
axC.plot(bcent,meandens,color=BLUEL,lw=2.6,label="Wt mean SHM density")
axc=axC.twinx(); axc.plot(bcent,meancov,color=GREY,lw=1.6,ls=(0,(4,2)))
axc.set_ylabel("Read coverage",color=GREY,fontsize=9.5); axc.tick_params(axis='y',colors=GREY); axc.set_ylim(0,600)
axC.set_xlim(LO,HI); axC.set_ylim(0,max(meandens)*1.30)
axC.set_xlabel("Position in J$_H$4 downstream intron amplicon (bp)",fontsize=9.5)
axC.set_ylabel("Mutation density\n(×10$^{-5}$ mut/bp)",color=BLUEL,fontsize=9.5); axC.tick_params(axis='y',colors=BLUEL)
axC.annotate("overlap with ~529 bp Sanger window",xy=(315,max(meandens)*1.17),
             ha="center",va="center",fontsize=8.5,color="#33608F",fontweight="bold",
             bbox=dict(facecolor="white",edgecolor="none",alpha=.82,pad=2),zorder=6)
axC.legend(loc="upper right",frameon=False,fontsize=8.5); axC.spines["top"].set_visible(False)

# Panel letters on ONE shared baseline in figure coordinates, sitting in the
# clear band above the axes. Placing them inside each axes previously put B
# on top of its own P-value label and C on top of panel C's frame.
_ylet=max(axA.get_position().y1,axB.get_position().y1,axC.get_position().y1)+0.052
for _ax,_lab,_dx in [(axA,"A",-0.048),(axB,"B",-0.022),(axC,"C",-0.052)]:
    fig.text(_ax.get_position().x0+_dx,_ylet,_lab,fontsize=20,fontweight="bold",
             va="top",ha="left")
OUT=OUT_DIR+os.sep
for ext in ("png","pdf"):
    fig.savefig(OUT+f"Figure2_JH4_combined.{ext}",dpi=600,bbox_inches="tight",pad_inches=.10,facecolor="white")
# Renderers (LibreOffice/Word) shave a few pixels off the top of a wide inline
# image, which was clipping the bold panel letters. Bake a small white margin
# into the raster so the letters can never be cropped.
from PIL import Image as _Image
_p=OUT+"Figure2_JH4_combined.png"
_im=_Image.open(_p).convert("RGB"); _w,_h=_im.size; _strip=int(_h*0.14)
_c=_Image.new("RGB",(_w,_h+_strip),"white"); _c.paste(_im,(0,_strip)); _c.save(_p)
print("top guard added: %dpx"%_strip)
print("Figure 2: ONT %s | Sanger %s"%(welch(ont["WT"],ont["KO"]),welch(sang["WT"],sang["KO"])))
print("Figure 2 (combined) written; WT full-length bins=%d"%nb)
