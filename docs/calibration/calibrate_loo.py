"""Leave-one-out calibration of the variant caller against AID-null controls."""
import sys, math
sys.path.insert(0,"/sessions/festive-dreamy-davinci/mnt/nanoshm/src")
import numpy as np
from collections import defaultdict
from nanoshm.config import load_preset
from nanoshm.io_utils import read_fastq, read_fasta
from nanoshm.reads_filter import prefilter, downsample
from nanoshm.align import build_index, align_reads

B="/sessions/festive-dreamy-davinci/mnt/nanoshm/examples_compare"
R=f"{B}/JH4_reference.fasta"; C=f"{B}/JH4_contaminant.fasta"
BASES=("A","C","G","T"); LO,HI=100,628
p=load_preset("jh4_intron")
ref=dict(read_fasta(R))["JH4_C57BL6"]

def pileup(tag):
    reads=list(read_fastq(f"{B}/{tag}.fastq")); kept,_=prefilter(reads,p); kept,_=downsample(kept,p)
    al,tn,cn=build_index([R],[C])
    counts=defaultdict(lambda: defaultdict(int)); depth=defaultdict(int); lam=defaultdict(lambda: defaultdict(float))
    for r in align_reads(kept,al,tn,cn,p):
        if r.ctg not in tn: continue
        for rp,b,q in r.pairs:
            if b not in BASES: continue
            counts[rp][b]+=1; depth[rp]+=1
            e3=(10**(-q/10.0))/3
            for a in BASES:
                if a!=b: lam[rp][a]+=e3
    return counts,depth,lam

PU={t:pileup(t) for t in ("KO1","KO2","KO3","WT2","WT3","WT4")}
print("pileups built", flush=True)

def slots(counts,depth):
    for pos in range(LO,HI+1):
        d=depth.get(pos,0)
        if d==0: continue
        rb=ref[pos].upper()
        if rb not in BASES: continue
        for alt in BASES:
            if alt==rb: continue
            k=counts.get(pos,{}).get(alt,0)
            if k==0: continue
            if k/d>=0.5: continue                 # germline rule
            yield pos,rb,alt,k,d

def pois_sf(k,l):
    if k<=0: return 1.0
    if l<=0: return 0.0
    t=math.exp(-l); c=t
    for i in range(1,k): t*=l/i; c+=t
    return max(0.0,1.0-c)

def binom_sf(k,n,q):
    """P(X>=k), X~Bin(n,q)  via survival using regularized incomplete beta."""
    from scipy.stats import binom
    return float(binom.sf(k-1,n,q))

def bh(pvals, alpha=0.05):
    """Benjamini-Hochberg: returns boolean keep-mask."""
    m=len(pvals)
    if m==0: return []
    order=np.argsort(pvals); keep=np.zeros(m,bool); thresh=0
    for rank,idx in enumerate(order,1):
        if pvals[idx] <= alpha*rank/m: thresh=rank
    keep[order[:thresh]]=True
    return keep

def empirical_rate(train_tags, alpha=0.5, beta=500.0):
    """Shrunk per-(pos,alt) error rate from AID-null training samples."""
    num=defaultdict(float); den=defaultdict(float)
    for t in train_tags:
        c,d,_=PU[t]
        for pos in range(LO,HI+1):
            dep=d.get(pos,0)
            if dep==0: continue
            rb=ref[pos].upper()
            if rb not in BASES: continue
            for a in BASES:
                if a==rb: continue
                k=c.get(pos,{}).get(a,0)
                if k/dep>=0.5: continue           # germline, not error
                num[(pos,a)]+=k; den[(pos,a)]+=dep
    return {key:(num[key]+alpha)/(den[key]+alpha+beta) for key in den}

def call(tag, strategy, eps=None, min_reads=1):
    c,d,l=PU[tag]; cand=list(slots(c,d))
    if not cand: return 0,0.0
    ps=[]
    for pos,rb,alt,k,dep in cand:
        if strategy.startswith("empirical"):
            q=eps.get((pos,alt), 0.001)
            ps.append(binom_sf(k,dep,q))
        else:
            lam=l[pos].get(alt,0.0) or dep*p.per_base_error/3
            ps.append(pois_sf(k,lam))
    ps=np.array(ps)
    if "bh" in strategy: keep=bh(ps,0.05)
    else: keep=ps<=0.05
    keep=keep & np.array([k>=min_reads for _,_,_,k,_ in cand])
    n=int(keep.sum())
    covered=[pos for pos in range(LO,HI+1) if d.get(pos,0)>0]
    mean_depth=sum(d.get(x,0) for x in covered)/max(1,len(covered))
    dens=n/(len(covered)*mean_depth)*1e5
    return n, dens

KO=("KO1","KO2","KO3"); WT=("WT2","WT3","WT4")
STRATS=[("A  baseline (current)","poisson",1),
        ("B  + min 3 reads","poisson",3),
        ("C  + BH FDR 5%","poisson_bh",1),
        ("D  empirical null + BH","empirical_bh",1),
        ("E  empirical + BH + min3","empirical_bh",3)]
print(f"\n{'strategy':26} {'KO bg (held out)':>18} {'WT signal':>12} {'ratio':>7}")
for label,strat,mr in STRATS:
    ko_d=[]; wt_d=[]
    for held in KO:                                  # leave-one-out over the AID-null mice
        train=[t for t in KO if t!=held]
        eps=empirical_rate(train) if strat.startswith("empirical") else None
        ko_d.append(call(held,strat,eps,mr)[1])
        if held==KO[0]:
            for w in WT: wt_d.append(call(w,strat,eps,mr)[1])
    k=np.mean(ko_d); w=np.mean(wt_d)
    print(f"{label:26} {k:12.2f} x1e-5 {w:9.1f} {w/k:7.1f}x")
print(f"\npaper target: KO background ~1.2 x1e-5, WT ~59 x1e-5  (ratio ~49x)")
