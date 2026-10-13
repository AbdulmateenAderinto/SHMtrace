#!/usr/bin/env python3
# =============================================================================
# figure3_concordance.py  ->  Figure 3 (composite A-E)
# -----------------------------------------------------------------------------
# ONT-vs-Sanger JH4 spectrum concordance (Wt, n=3 paired mice). Panels:
#   A  pooled From->To substitution matrices (Sanger | ONT), shared colour scale
#   B  12 substitution-type proportions (paired tests + Sidak, adj p shown)
#   C  Sanger-vs-ONT agreement of the 12 proportions (Pearson r, slope)
#   D  G:C transition/transversion proportions (paired tests + Sidak)
#   E  mutations-per-sequence distribution (Sanger vs ONT)
# ONT per-mouse matrices are the July-13 phased sub500 values (529 bp window),
# inlined for portability. Statistics were recalculated as mouse-matched tests.
# from the original analysis. Recompute the ONT matrices with recompute_jh4.py.
# Run:  python figure3_concordance.py
# =============================================================================
import os
from config import OUT_DIR
# -*- coding: utf-8 -*-
import numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Circle, Patch, Rectangle
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":12,"axes.edgecolor":"black",
 "xtick.color":"black","ytick.color":"black","text.color":"#1a1a1a","axes.labelcolor":"#1a1a1a",
 "savefig.dpi":600,"figure.facecolor":"white"})
BLUE=("#6D9AD0","#274B7E"); RED=("#D07A7E","#9A3439")
# Sanger vs ONT platform palette (colour-blind-safe: blue + reddish-purple, Okabe–Ito)
BLK=("#0072B2","#005A8D")     # "BLK" name retained; now Sanger = blue
TEALP=("#CC79A7","#A85C86")   # "TEALP" name retained; now ONT = reddish-purple
SEQ=LinearSegmentedColormap.from_list("seqblue",["#F5F8FC","#C4D8EE","#6FA0CE","#2F6BA6","#173D63"])
# Reference-style sequential blue/purple mutation bins for the two donuts.
PIE=['#FFFFFF','#DCEAF7','#A8CBE8','#5A9BD5','#1F4E79']; GRID="#C9C9C9"
BASES=["A","T","C","G"]; IDX={b:i for i,b in enumerate(BASES)}
MUT=["A>T","A>C","A>G","T>A","T>C","T>G","C>A","C>T","C>G","G>A","G>T","G>C"]
P12=None   # computed below from the matrices (paired t, Šídák over 12 types)
def m(r): return np.array(r,float)
def c(M,f,t): return M[IDX[t]][IDX[f]]
# CORRECTED SANGER per-mouse matrices M[to][from], built from the organised .ab1
# traces by shmtrace (docs/calibration/SANGER_VS_NANOPORE.md)
import json as _json
from config import DATA_SANGER as _DS
_cs=_json.load(open(os.path.join(_DS,"corrected_sanger.json")))
def _mat(vs):
    M=np.zeros((4,4))
    for _p,f,t,_n in vs: M[IDX[t]][IDX[f]]+=1
    return M
sanger={k.replace("Wt ","WT"):_mat(v["variants"]) for k,v in _cs.items() if k.startswith("Wt")}
_loads=[x for k,v in _cs.items() if k.startswith("Wt") for x in v["loads"]]
SANGER_PIE=[sum(1 for x in _loads if a<=x<=b) for a,b in ((0,0),(1,3),(4,10),(11,20),(21,10**6))]
SANGER_N=len(_loads)
# ONT per-mouse matrices M[to][from]: SHMtrace calls in the 529-bp window, 500 reads per
# mouse (docs/calibration/jh4_shmtrace.py; previously the Geneious July 13 exports)
from pathlib import Path as _Path
_JH4=_json.load(open(_Path(__file__).resolve().parents[2]/"figures/jh4_shmtrace/project.json"))
def _ont_mat(k):
    M=np.zeros((4,4))
    for ft,n in _JH4["spectrum_sample"][k].items(): M[IDX[ft[2]]][IDX[ft[0]]]=n
    return M
ont={k:_ont_mat(k) for k in ("WT2","WT3","WT4")}
# per-read called mutations (ONT donut, panel E): the same reads and calls
_oload={}
for k in ("WT2","WT3","WT4"):
    for x,n in _JH4["mutation_load"][k].items(): _oload[int(x)]=_oload.get(int(x),0)+n
ONT_PIE=[sum(n for x,n in _oload.items() if a<=x<=b) for a,b in ((0,0),(1,3),(4,10),(11,20),(21,10**6))]
ONT_N=sum(_oload.values())
S=["WT2","WT3","WT4"]; SM=[sanger[s] for s in S]; OM=[ont[s] for s in S]
from scipy import stats as _st
def _fmtp(p): return ">0.9999" if p>0.9999 else f"{p:.4f}"
def _sidak(p,k): return 1-(1-p)**k
def _wilson(k,n,z=1.96):
    if n==0: return (0.0,0.0)
    ph=k/n; d=1+z*z/n; cc=(ph+z*z/(2*n))/d; h=z*np.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d
    return 100*(cc-h),100*(cc+h)
def style(ax):
    ax.set_axisbelow(True); ax.yaxis.grid(True,linestyle=(0,(5,4)),color=GRID,linewidth=1.1,zorder=0)
    for sp in ('left','bottom'): ax.spines[sp].set_linewidth(1.9); ax.spines[sp].set_color('black')
    ax.spines[['top','right']].set_visible(False); ax.tick_params(width=1.7,length=6,labelsize=12)
def pts(ax,xc,vals,color,r=0.052,jit=0.02,seed=0,s=62,mec=None):
    ec=mec if mec else color; rng=np.random.default_rng(seed)
    for v in vals:
        x=xc+(rng.random()-.5)*jit
        ax.scatter(x,v,marker='o',facecolors='white',edgecolors=ec,s=s,linewidths=1.6,zorder=6)
        ax.plot([x-r,x+r],[v,v],color=ec,lw=1.4,zorder=7)
def bracket(ax,x1,x2,y,p,dh=0.5,fs=9.5):
    ax.plot([x1,x1,x2,x2],[y,y+dh,y+dh,y],color="black",lw=1.2,zorder=8)
    ax.text((x1+x2)/2,y+dh*1.25,p,ha="center",va="bottom",fontsize=fs)
def props12(mats):
    A=np.zeros((12,len(mats)))
    for j,M in enumerate(mats):
        tot=M.sum()
        for i,mt in enumerate(MUT): f,t=mt.split(">"); A[i,j]=100*c(M,f,t)/tot
    return A
def gctstv(mats):
    ts=[];tv=[]
    for M in mats:
        tot=M.sum(); ts.append(100*(c(M,"G","A")+c(M,"C","T"))/tot)
        tv.append(100*(c(M,"G","T")+c(M,"G","C")+c(M,"C","A")+c(M,"C","G"))/tot)
    return np.array(ts),np.array(tv)

# Compact three-row narrative. The full-width matrix and 12-type panels retain
# their readable scales, while C, D, and the two-donut E panel share one dense
# bottom row. This removes the empty half-row areas created by the former
# four-subfigure stack and lets Word display the composite at full text width.
fig=plt.figure(figsize=(15.2,10.6))
outer=fig.add_gridspec(
    2,1,height_ratios=[1.02,1.0],hspace=.20,
    left=.045,right=.985,top=.965,bottom=.055
)
# Reviewer request: panels A and B sit on the same level, with the two matrices
# pulled together (no dead space between them) and B narrowed to fit alongside.
top=outer[0].subgridspec(1,2,width_ratios=[1.30,1.0],wspace=.10)
# ---------- Panel A heatmaps (styled to match Supplementary Figure S2) ----------
# The right matrix has a narrower x-range (no rotated "From" label), so it needs
# a proportionally narrower cell. With aspect="equal" and equal cells it would
# otherwise render at a larger scale than the left matrix, which also pushed its
# title to a different height. Ratios = each panel's x-span: 5.48-(-1.35)=6.83
# and 5.48-(-0.72)=6.20.
gsA=top[0,0].subgridspec(2,2,height_ratios=[18,1.1],hspace=.02,wspace=.01,
                         width_ratios=[6.83,6.20])
axes=[fig.add_subplot(gsA[0,0]),fig.add_subplot(gsA[0,1])]
def heat(ax,M,title,show_from=True):
    # From = rows, To = columns (identical layout to Suppl. Fig. S2). The blue
    # summary column gives the per-From-base total; a shared colour bar sits
    # below both matrices.
    MT=M.T; P=100*MT/MT.sum(); rowpct=100*MT.sum(1)/MT.sum()
    for i in range(4):
        for j in range(4):
            x,y=j,3-i
            if i==j:
                ax.add_patch(Rectangle((x,y),1,1,fc="white",ec="white"))
                ax.text(x+.5,y+.5,"—",ha="center",va="center",color="#AFAFAF",fontsize=12)
                continue
            v=P[i,j]
            ax.add_patch(Rectangle((x,y),1,1,fc=SEQ(v/gvmax),ec="white",lw=2))
            tc="white" if v/gvmax>0.55 else "#1a1a1a"
            ax.text(x+.5,y+.60,f"{int(MT[i,j])}",ha="center",va="center",
                    color=tc,fontsize=12,fontweight="bold")
            ax.text(x+.5,y+.30,f"{v:.1f}%",ha="center",va="center",
                    color=tc,fontsize=7.4)
    sx=4.28
    for i in range(4):
        ax.add_patch(Rectangle((sx,3-i),1,1,fc=SEQ(rowpct[i]/max(rowpct)*.6),ec="white",lw=2))
        ax.text(sx+.5,3-i+.5,f"{rowpct[i]:.1f}",ha="center",va="center",
                fontsize=10,fontweight="bold")
    for j,b in enumerate(BASES):
        ax.text(j+.5,4.08,b,ha="center",va="bottom",fontsize=11.5,fontweight="bold")
    ax.text(sx+.5,4.08,"%",ha="center",va="bottom",fontsize=11,fontweight="bold")
    ax.text(1.5,4.58,"To",ha="center",va="bottom",fontsize=11,fontweight="bold")
    for i,b in enumerate(BASES):
        ax.text(-.42,3-i+.5,b,ha="center",va="center",fontsize=11.5,fontweight="bold")
    # The rotated "From" label is only drawn on the left matrix; dropping it on the
    # right one removes the dead margin that separated the two panels.
    if show_from:
        ax.text(-1.02,2,"From",ha="center",va="center",fontsize=11,fontweight="bold",rotation=90)
    lo=-1.35 if show_from else -0.72
    ax.text((lo+5.48)/2,5.18,title,ha="center",va="bottom",fontsize=13,fontweight="bold")
    ax.set_xlim(lo,5.48); ax.set_ylim(-.05,5.55); ax.set_aspect("equal"); ax.axis("off")
MTs=sum(SM); MTo=sum(OM)
gvmax=max(100*MTs.max()/MTs.sum(),100*MTo.max()/MTo.sum())
heat(axes[0],MTs,"$\\mathbfit{Wt}$ (Sanger)"); heat(axes[1],MTo,"$\\mathbfit{Wt}$ (ONT)",show_from=False)
sm=plt.cm.ScalarMappable(cmap=SEQ,norm=plt.Normalize(0,gvmax))
cax=fig.add_subplot(gsA[1,:]); cb=fig.colorbar(sm,cax=cax,orientation="horizontal")
cb.set_label("% of total mutations",fontsize=9); cb.ax.tick_params(labelsize=7.5)

# ---------- Panel B: 12-type spectrum ----------
# Constrain B to ~62% width (centred) so the 12-category bar chart is not
# stretched horizontally across the full figure width.
axB=fig.add_subplot(top[0,1])
style(axB); g1=props12(SM); g2=props12(OM); x=np.arange(12); off=0.205; bw=0.36
# CORRECTED: bars are each platform's pooled proportion (all unique mutations of
# the three mice) with a 95% Wilson interval; a mean of three mice let Sanger Wt2
# (6 mutations) swing a type by tens of points. Points: each mouse, area by its
# number of mutations. Labels: paired t on the per-mouse proportions, Šídák x12.
P12=[_fmtp(_sidak(float(_st.ttest_rel(g1[i],g2[i]).pvalue),12)) for i in range(12)]
_hi=[]
for mats,arr,o,(fill,edge),lab in [(SM,g1,-off,BLK,"Sanger"),(OM,g2,off,TEALP,"ONT")]:
    T=sum(mats); N=T.sum(); k=np.array([c(T,*mt.split(">")) for mt in MUT])
    pooled=100*k/N; ci=np.array([_wilson(kk,N) for kk in k]); _hi.append(ci[:,1])
    axB.bar(x+o,pooled,bw,yerr=[pooled-ci[:,0],ci[:,1]-pooled],capsize=3,color=fill,edgecolor=edge,linewidth=1.7,error_kw=dict(ecolor=edge,lw=1.4),label=f"{lab} (n = {int(N)})",zorder=2)
    nm=[m.sum() for m in mats]; sz=[max(8,46*np.sqrt(n_/max(nm))) for n_ in nm]
    for i in range(12):
        for kk in range(3):
            axB.scatter(x[i]+o+(kk-1)*.06,min(arr[i][kk],25.5),s=sz[kk],facecolor="white",edgecolor=edge,lw=1.2,zorder=6)
for i in range(12):
    ytop=max(_hi[0][i],_hi[1][i])+0.7
    bracket(axB,x[i]-off,x[i]+off,ytop,P12[i],dh=0.5,fs=6.2)
axB.set_xticks(x); axB.set_xticklabels([q.replace(">","→") for q in MUT],fontsize=8.6,rotation=45,ha="right")
axB.set_ylabel("Proportion of total mutations (%)",fontsize=9.5); axB.set_ylim(0,26)
axB.legend(loc="upper center",bbox_to_anchor=(.5,1.06),ncol=2,frameon=False,fontsize=9)
axB.margins(x=.02)
# A and B share one baseline in figure coords so they line up exactly
_topAB=max(axes[0].get_position().y1,axes[1].get_position().y1,axB.get_position().y1)+0.012
fig.text(axes[0].get_position().x0-0.030,_topAB,"A",fontsize=22,fontweight="bold",va="top",ha="left")
fig.text(axB.get_position().x0-0.055,_topAB,"B",fontsize=22,fontweight="bold",va="top",ha="left")

# ---------- Panel C agreement + Panel D G:C Ts/Tv ----------
gsBottom=outer[1].subgridspec(1,3,width_ratios=[1.08,.82,2.10],wspace=.25)
axC=fig.add_subplot(gsBottom[0,0]); axD=fig.add_subplot(gsBottom[0,1])
def props_flat(M):
    tot=M.sum(); return np.array([100*c(M,*mt.split(">"))/tot for mt in MUT])
sp=props_flat(sum(SM)); op=props_flat(sum(OM))
r=np.corrcoef(sp,op)[0,1]; slope,intc=np.polyfit(sp,op,1)
lim=max(sp.max(),op.max())*1.12
axC.plot([0,lim],[0,lim],ls="--",color="#999",lw=1.3,zorder=1)
axC.scatter(sp,op,s=85,facecolor=TEALP[1],edgecolor="white",lw=0.9,zorder=3)
# Structured labels: every leader is horizontal and 26 points long. Types are
# assigned to the left or right column to preserve vertical separation.
label_side={
    "A>T":1, "A>C":-1, "A>G":1,
    "T>A":1, "T>C":1, "T>G":-1,
    "C>A":1, "C>T":1, "C>G":1,
    "G>A":-1, "G>T":-1, "G>C":1,
}
for i,t in enumerate(MUT):
    dx=20*label_side[t]
    axC.annotate(t.replace(">","→"),(sp[i],op[i]),fontsize=8.0,
                 xytext=(dx,0),textcoords="offset points",
                 ha="right" if dx<0 else "left",
                 va="center",
                 color="#333",fontweight="bold",
                 arrowprops=dict(arrowstyle="-",color="#777",lw=.55,
                                 shrinkA=2,shrinkB=4),annotation_clip=False)
axC.set_xlim(0,lim); axC.set_ylim(0,lim); axC.set_aspect("equal")
axC.set_xlabel("Sanger (% of total mutations)",fontsize=9.5); axC.set_ylabel("ONT (% of total mutations)",fontsize=9.5)
axC.text(.97,.055,f"Pearson r = {r:.2f}\nslope = {slope:.2f}",
         transform=axC.transAxes,fontsize=10.5,color="#0072B2",fontweight="bold",ha="right",va="bottom")
for s2 in ["top","right"]: axC.spines[s2].set_visible(False)
axC.tick_params(labelsize=9,width=1.3,length=4)
axC.spines['left'].set_linewidth(1.4); axC.spines['bottom'].set_linewidth(1.4)

st,sv=gctstv(SM); ot,ov=gctstv(OM); style(axD)
x=np.arange(2); off=0.20; bw=0.30
for mats,(tt,vv),o,(fill,edge),lab in [(SM,(st,sv),-off,BLK,"Sanger"),(OM,(ot,ov),off,TEALP,"ONT")]:
    V=np.array([tt,vv]); T=sum(mats); N=T.sum()
    k=np.array([c(T,"G","A")+c(T,"C","T"), c(T,"G","T")+c(T,"G","C")+c(T,"C","A")+c(T,"C","G")])
    mn=100*k/N; ci=np.array([_wilson(kk,N) for kk in k])
    axD.bar(x+o,mn,bw,yerr=[mn-ci[:,0],ci[:,1]-mn],capsize=4,color=fill,edgecolor=edge,linewidth=1.8,error_kw=dict(ecolor=edge,lw=1.4),label=lab,zorder=2)
    mec=edge
    for i in range(2): pts(axD,x[i]+o,V[i],edge,r=0.055,jit=0.035,seed=2,s=44,mec=mec)
_ptt=[_sidak(float(_st.ttest_rel(a_,b_).pvalue),2) for a_,b_ in ((st,ot),(sv,ov))]
for i,pt in enumerate([f"P = {p:.4f}" for p in _ptt]):
    pair_top=max((st if i==0 else sv).max(),(ot if i==0 else ov).max())+1.0
    bracket(axD,x[i]-off,x[i]+off,pair_top,pt,dh=0.7,fs=9.5)
axD.set_xticks(x); axD.set_xticklabels(["G:C Ts","G:C Tv"],fontsize=12)
axD.set_ylabel("Proportion of total mutations (%)",fontsize=9.5)
axD.set_ylim(0,42); axD.legend(loc="upper center",bbox_to_anchor=(.5,1.16),ncol=2,frameon=False,fontsize=10)
axD.set_box_aspect(1.08); axD.set_anchor("C")

# ---------- Panel E: compact, reference-style mutations-per-sequence donuts ----------
gsE=gsBottom[0,2].subgridspec(1,3,width_ratios=[1,1,.58],wspace=.02)
axd0=fig.add_subplot(gsE[0,0]); axd1=fig.add_subplot(gsE[0,1]); axleg=fig.add_subplot(gsE[0,2])
axd=[axd0,axd1]
LAB=["0","1–3","4–10","11–20","20+"]
def label_colour(fill):
    """Choose readable percentage text for light and dark donut wedges."""
    r,g,b=to_rgb(fill)
    luminance=.2126*r+.7152*g+.0722*b
    return "white" if luminance<.50 else "#1a1a1a"
# The Sanger distribution contains 92 evaluable C57BL/6-background sequences
# plus 19 largely unmutated Agouti-background sequences. The latter contribute
# to the sequence denominator and the zero-mutation class, but not to the
# pooled mutation spectrum. Counts are therefore 50, 12, 15, 19, and 15 of 111.
# Panel E: Sanger clones (corrected_sanger.json) and ONT reads (SHMtrace per-read calls).
for a,vals,N,lab in [(axd[0],SANGER_PIE,SANGER_N,"$\\mathbfit{Wt}$ (Sanger)"),(axd[1],ONT_PIE,ONT_N,"$\\mathbfit{Wt}$ (ONT)")]:
    v=np.array(vals,float)
    radius=1.0; width=.58; inner=radius-width
    wed,_=a.pie(v,colors=PIE,startangle=90,counterclock=False,radius=radius,
                wedgeprops=dict(width=width,edgecolor="none"))
    # The reference figure emphasizes the donut silhouette rather than every
    # internal wedge boundary: one strong outer ring and one strong inner ring.
    a.add_patch(Circle((0,0),radius,fill=False,ec="black",lw=2.5,zorder=5))
    a.add_patch(Circle((0,0),inner,fill=False,ec="black",lw=2.0,zorder=5))
    for k,(wd,pp) in enumerate(zip(wed,100*v/v.sum())):
        if pp<1: continue
        ang=np.deg2rad((wd.theta1+wd.theta2)/2); cx,cy=np.cos(ang),np.sin(ang)
        if pp>=5:
            # All internal percentages use one radial baseline. This keeps the
            # labels visually aligned instead of pushing the 34%/51% values
            # toward the outside ring.
            rr=.75
            a.text(rr*cx,rr*cy,f"{pp:.0f}%",ha="center",va="center",
                   fontsize=12.6,fontweight="bold",color=label_colour(PIE[k]))
        else:
            tc=PIE[k] if k>0 else "#1a1a1a"
            # The lone 3% label is centred above its slice with a short vertical
            # leader, avoiding the slightly diagonal placement in the old panel.
            a.annotate(f"{pp:.0f}%",xy=(0,.98),xytext=(0,1.18),
                       ha="center",va="center",fontsize=12.0,fontweight="bold",fontstyle="italic",color=tc,
                       arrowprops=dict(arrowstyle="-",color=tc,lw=1.0),annotation_clip=False)
    a.text(0,0,f"$n$ = {N:,}",ha="center",va="center",fontsize=10.5)
    # Both method labels share the same baseline and identical axes extents.
    a.text(0,-1.16,lab,ha="center",va="top",fontsize=12.8,fontweight="bold")
    a.set_aspect("equal"); a.set_xlim(-1.10,1.10); a.set_ylim(-1.27,1.24); a.axis("off")
leg=[Patch(facecolor=cc,edgecolor="black",label=l,linewidth=1.2) for cc,l in zip(PIE,LAB)]
axleg.axis("off")
axleg.legend(handles=leg,title="Mutations\nper sequence",loc="center left",
             bbox_to_anchor=(-.05,.50),frameon=False,fontsize=10.5,title_fontsize=11.0,handlelength=1.65)
axd0.text(1.01,1.03,"$J_{H4}$",transform=axd0.transAxes,fontsize=12.5,
          fontstyle="italic",ha="center",va="bottom")

# Reserve a clean, shared header lane for the bottom-row panel letters.  Using
# figure coordinates keeps C/D/E on one baseline and prevents them from
# colliding with C's top y tick, D's y-axis/legend, or E's donut annotations.
bottom_header_y=max(axC.get_position().y1,axD.get_position().y1,
                    axd0.get_position().y1)+.018
for letter,axis,xpad in [("C",axC,.032),("D",axD,.030),("E",axd0,.030)]:
    fig.text(axis.get_position().x0-xpad,bottom_header_y,letter,
             fontsize=22,fontweight="bold",ha="left",va="bottom")
OUT=OUT_DIR+os.sep
for ext in ("png","pdf"):
    fig.savefig(OUT+f"Figure3_concordance.{ext}",dpi=600,
                bbox_inches="tight",pad_inches=.04,facecolor="white")
print(f"Fig3 saved. ONT n={int(sum(OM).sum())} r={r:.3f} slope={slope:.2f} ONT pie {ONT_PIE} of {ONT_N}")
# Renderers (LibreOffice/Word) shave a few pixels off the top of a wide inline
# image, which was clipping the bold panel letters. Bake a small white margin
# into the raster so the letters can never be cropped.
from PIL import Image as _Image
_p=OUT+"Figure3_concordance.png"
_im=_Image.open(_p).convert("RGB"); _w,_h=_im.size; _strip=int(_h*0.10)
_c=_Image.new("RGB",(_w,_h+_strip),"white"); _c.paste(_im,(0,_strip)); _c.save(_p)
print("top guard added: %dpx"%_strip)
