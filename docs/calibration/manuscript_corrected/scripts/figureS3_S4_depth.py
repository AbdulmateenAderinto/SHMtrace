#!/usr/bin/env python3
# =============================================================================
# figureS3_S4_depth.py  ->  Supplementary Figures S3 and S4
# -----------------------------------------------------------------------------
# Down-sampling analysis rebuilt from the depth-subsampled Geneious variant
# exports (sub30..sub1000) for the JH4 intron (WT/Aicda-/-) and the N114 V
# region (pMX-Aicda / pMX-Empty vector).
#   S3  A/C  unique variants detected vs depth   (variant discovery)
#       B/D  mutation density vs depth              (per-read, depth-dependent)
#   S4  A/B  per-sample discovery normalised to its own maximum (saturation)
# Filters match the pipeline: single-nt substitutions, Variant P-Value <= 0.05,
# germline/parental removed (Variant Frequency < 50%), restricted to the
# reliably aligned window (JH4 101-629 = 529 bp; V region 21-635 = 615 bp).
# Labels use "mutation density" (x10^-5 mut/bp) and the "pMX-Empty vector"
# naming, consistent with Figures 2 and 4.
# Run:  python figureS3_S4_depth.py
# =============================================================================
import os, csv, glob, math, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from config import OUT_DIR
mpl.rcParams.update({"font.family":"DejaVu Sans","savefig.dpi":600,"font.size":12,
    "axes.edgecolor":"black","xtick.color":"black","ytick.color":"black",
    "text.color":"#1a1a1a","axes.labelcolor":"#1a1a1a",
    # Regular-weight italic gene names, matching Figures 2/4 and Supplementary S5.
    # (These previously rendered bold-italic, which was the only place in the
    # figure set where "Wt", "Aicda-/-" and "pMX-Aicda" were heavier than the
    # surrounding text.)
    "mathtext.fontset":"custom","mathtext.rm":"DejaVu Sans","mathtext.bf":"DejaVu Sans:bold",
    "mathtext.it":"DejaVu Sans:italic","mathtext.cal":"DejaVu Sans:italic"})

HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
DIR_JH4=os.path.join(ROOT,"July 13th JH4")
DIR_HYB=os.path.join(ROOT,"July 10th ","Hybridoma")
DEPTHS=[30,60,125,250,500,1000]
BLUE="#0072B2"; BLUE_E="#005A8D"; RED="#C0433F"; RED_E="#8E2F2C"; GRID="#C9C9C9"; VLINE="#9AA3AC"

def pf(x):
    if x is None: return None
    x=str(x).replace('%','').replace(',','').strip()
    try: return float(x)
    except: return None

def uniq_count(path, lo, hi):
    if not path or not os.path.exists(path): return None
    n=0
    for r in csv.DictReader(open(path,encoding='utf-8-sig')):
        if len(r.get('Reference Nucleotide(s)','') or '')!=1 or len(r.get('Variant Nucleotide(s)','') or '')!=1: continue
        p=pf(r.get('Variant P-Value (approximate)'));  vf=pf(r.get('Variant Frequency'))
        if p is None or p>0.05 or vf is None or vf>=50: continue
        pos=pf(r.get('Minimum')) or pf(r.get('Maximum'))
        if pos is None or pos<lo or pos>hi: continue
        n+=1
    return n

def find(dirpath, token, depth, suffix):
    g=glob.glob(os.path.join(dirpath, f"*{token}*sub{depth} {suffix}"))
    if len(g)!=1:
        raise RuntimeError(f"token '{token}' matched {len(g)} files at depth {depth} (must be exactly 1): "
                           +", ".join(os.path.basename(x) for x in g))
    return g[0]

# subclone / mouse tokens - fully specified so each matches EXACTLY ONE file.
# The three independent AID subclones (C1, D3, D5) are pinned to the SAME
# sequencing runs used in Figure 4 (C1 = the October 10.06.24 run).
JH4={"WT":["Wildtype_2","Wildtype_3","Wildtype_4"], "KO":["AID_KO_1","AID_KO_2","AID_KO_3"]}
HYB={"AID":["C1_AID_10.06.24","D3_AID_10.01","D5_AID_6.28"],
     "EMP":["C3_Empty_March_2024_1","C4_Empty_March_2024","D1_Empty_March_2024"]}
# Note: the three AID files above are exactly those used in Figure 4
# (C1 October 10.06.24 = 80, D3 10.01 = 71, D5 6.28 = 147 unique variants at 500 reads).
# C1 uses the October run (all subclones sequenced in the same batch; March C1 archived).
JH4_SUF="assembled to Jh4 Annotations.csv"; HYB_SUF="assembled to seq_23353 Annotations.csv"

def series(system):
    """returns dict grp -> list-over-depth of [per-sample counts]"""
    out={}
    if system=="JH4":
        lo,hi,bp,d_dir,suf,groups=101,629,529,DIR_JH4,JH4_SUF,JH4
    else:
        lo,hi,bp,d_dir,suf,groups=21,635,615,DIR_HYB,HYB_SUF,HYB
    for grp,toks in groups.items():
        perdepth=[]
        for d in DEPTHS:
            vals=[uniq_count(find(d_dir,t,d,suf),lo,hi) for t in toks]
            perdepth.append([v for v in vals if v is not None])
        out[grp]=perdepth
    return out, bp

def geo(xs, floor=1.0):
    xs=[max(x,floor) for x in xs]; m=np.mean(np.log10(xs)); s=np.std(np.log10(xs),ddof=1) if len(xs)>1 else 0
    return 10**m, 10**s

def plot_depth(ax, counts_by_depth, color, edge, label, bp=None):
    x=np.arange(len(DEPTHS)); gm=[]; lo=[]; hi=[]
    for i,cs in enumerate(counts_by_depth):
        if bp is None:
            # A zero-variant sample receives a one-variant pseudocount only for
            # the geometric summary; the plotted per-sample counts are otherwise
            # unchanged.
            vals=cs
            g,s=geo(vals)
        else:
            # Apply the same one-variant pseudocount BEFORE converting counts to
            # density. Flooring the already-converted density at 1 x 10^-5 would
            # be depth-dependent and would not match the stated method.
            vals=[max(c,1)/(DEPTHS[i]*bp)*1e5 for c in cs]
            g,s=geo(vals,floor=1e-300)
        gm.append(g); lo.append(g/s); hi.append(g*s)
    gm=np.array(gm); lo=np.array(lo); hi=np.array(hi)
    ax.errorbar(x,gm,yerr=[gm-lo,hi-gm],fmt='-o',color=color,mfc='white',mec=edge,mew=1.8,
                ms=9,lw=2.0,capsize=4,ecolor=edge,elinewidth=1.6,label=label,zorder=5)
    return gm

def style(ax, ylab, title):
    ax.set_yscale("log"); ax.set_xticks(range(len(DEPTHS))); ax.set_xticklabels(DEPTHS)
    ax.set_xlabel("Reads subsampled"); ax.set_ylabel(ylab)
    ax.axvline(DEPTHS.index(500),color=VLINE,ls="--",lw=1.6,zorder=1)
    ax.set_title(title,loc="left",fontweight="bold",fontsize=13,pad=8)
    ax.yaxis.grid(True,ls=(0,(5,4)),color=GRID,lw=0.9,zorder=0); ax.set_axisbelow(True)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

# ---------------------------------------------------------------- data --------
jh4,bpJ=series("JH4"); hyb,bpV=series("HYB")
# Units follow the terminology applied throughout the manuscript (Figures 2/4):
# "mutation density" in units of x10^-5 mutations per base pair.
YDENS="Mutation density\n(×10$^{-5}$ mut/bp)"
YCOUNT="Unique variants detected"
# One legend per row only: A/B share the Wt vs Aicda-/- key and C/D share the
# pMX key, so repeating it in all four panels was redundant clutter that also
# forced a different placement in each panel.
# Anchored just below mid-height on the left: on both log axes this band sits in
# the clear gap between the upper (blue) and lower (red) series in A and C.
LEGKW=dict(loc="center left",bbox_to_anchor=(0.12,0.36))

# ===================== Supplementary Figure S3 ================================
fig,axes=plt.subplots(2,2,figsize=(13,9.2))
plot_depth(axes[0,0],jh4["WT"],BLUE,BLUE_E,"$\\it{Wt}$")
plot_depth(axes[0,0],jh4["KO"],RED,RED_E,"$\\it{Aicda}^{-/-}$")
style(axes[0,0],YCOUNT,"$\\mathbfit{J_{H}4}$ intron"); axes[0,0].legend(frameon=False,fontsize=11,**LEGKW)
plot_depth(axes[0,1],jh4["WT"],BLUE,BLUE_E,"$\\it{Wt}$",bp=bpJ)
plot_depth(axes[0,1],jh4["KO"],RED,RED_E,"$\\it{Aicda}^{-/-}$",bp=bpJ)
style(axes[0,1],YDENS,"$\\mathbfit{J_{H}4}$ intron")
plot_depth(axes[1,0],hyb["AID"],BLUE,BLUE_E,"pMX-$\\it{Aicda}$")
plot_depth(axes[1,0],hyb["EMP"],RED,RED_E,"pMX-Empty vector")
style(axes[1,0],YCOUNT,"N114 V region"); axes[1,0].legend(frameon=False,fontsize=11,**LEGKW)
plot_depth(axes[1,1],hyb["AID"],BLUE,BLUE_E,"pMX-$\\it{Aicda}$",bp=bpV)
plot_depth(axes[1,1],hyb["EMP"],RED,RED_E,"pMX-Empty vector",bp=bpV)
style(axes[1,1],YDENS,"N114 V region")
fig.tight_layout(rect=(0,0,1,0.985))
# Panel letters on ONE shared baseline per row (figure coords), matching the
# treatment used in Figures 3/4 and Supplementary S1/S5.
for _row,_labs in ((0,"AB"),(1,"CD")):
    _ytop=max(axes[_row,c].get_position().y1 for c in (0,1))+0.022
    for _c,_lab in enumerate(_labs):
        fig.text(axes[_row,_c].get_position().x0-0.052,_ytop,_lab,
                 fontsize=20,fontweight="bold",va="top",ha="left")
for ext in ("png","pdf"):
    fig.savefig(os.path.join(OUT_DIR,f"FigureS3_depth_discovery_density.{ext}"),dpi=600,bbox_inches="tight",facecolor="white")
plt.close(fig)
# Renderers (LibreOffice/Word) shave a few pixels off the top of a wide inline
# image, which was clipping the bold panel letters. Bake a small white margin
# into the raster so the letters can never be cropped.
from PIL import Image as _Image
_p=os.path.join(OUT_DIR,"FigureS3_depth_discovery_density.png")
_im=_Image.open(_p).convert("RGB"); _w,_h=_im.size; _strip=int(_h*0.045)
_c=_Image.new("RGB",(_w,_h+_strip),"white"); _c.paste(_im,(0,_strip)); _c.save(_p)
print("S6 top guard added: %dpx"%_strip)

# ===================== Supplementary Figure S4 ===============================
def saturation(counts_by_depth):
    # per sample: value at each depth / that sample's own max, in %
    nS=min(len(c) for c in counts_by_depth)
    perc=[]
    for s in range(nS):
        vals=[counts_by_depth[i][s] for i in range(len(DEPTHS))]
        mx=max(vals) or 1
        perc.append([100*v/mx for v in vals])
    arr=np.array(perc); return arr.mean(0), arr.std(0,ddof=1), arr

fig2,ax2=plt.subplots(1,2,figsize=(13,5.2))
for ax,(cbd,ttl,leg) in zip(ax2,[(jh4["WT"],"A   $\\mathbfit{J_{H}4}$ intron (wild-type)","$\\it{Wt}$ mean ± SD (n = 3)"),
                                  (hyb["AID"],"B   N114 V region (pMX-$\\mathbfit{Aicda}$)","pMX-$\\it{Aicda}$ mean ± SD (n = 3)")]):
    mean,sd,_=saturation(cbd); x=np.arange(len(DEPTHS))
    ax.errorbar(x,mean,yerr=sd,fmt='-o',color=BLUE,mfc='white',mec=BLUE_E,mew=1.8,ms=9,lw=2.2,capsize=4,ecolor=BLUE_E,elinewidth=1.6,label=leg,zorder=5)
    ax.axhline(100,color="#9AA3AC",ls=":",lw=1.3,zorder=1)
    ax.axvline(DEPTHS.index(500),color=VLINE,ls="--",lw=1.6,zorder=1)
    p500=mean[DEPTHS.index(500)]
    ax.text(0.12,90,f"{p500:.0f}% recovered\nat 500 reads",color=BLUE_E,fontsize=13,fontweight="bold",ha="left",va="top")
    ax.set_xticks(x); ax.set_xticklabels(DEPTHS); ax.set_xlabel("Reads subsampled")
    ax.set_ylabel("Variants recovered (% of maximum)"); ax.set_ylim(0,112)
    _letter,_,_rest=ttl.partition("   ")
    ax.text(-0.12,1.06,_letter.strip(),transform=ax.transAxes,fontsize=20,
            fontweight="bold",va="top",ha="left")
    ax.set_title(_rest.strip(),loc="left",fontweight="bold",fontsize=13)
    ax.legend(frameon=False,fontsize=11,loc="lower right")
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig2.tight_layout()
for ext in ("png","pdf"):
    fig2.savefig(os.path.join(OUT_DIR,f"FigureS4_saturation.{ext}"),dpi=600,bbox_inches="tight",facecolor="white")
plt.close(fig2)
print("S3 + S4 written to",OUT_DIR)  # S5 (Sanger-vs-Nanopore) removed per author request
