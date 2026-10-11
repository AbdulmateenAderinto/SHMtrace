import json, statistics as st, random
from collections import defaultdict
from dataclasses import replace
from nanoshm import pipeline, motifs, null_model
from nanoshm.config import load_preset
from nanoshm.io_utils import read_fasta, read_fastq, phred
from nanoshm.align import build_index
from nanoshm.variants import _geneious_pvalue
E="examples_compare/"; REF=E+"JH4_reference.fasta"; CON=E+"JH4_contaminant.fasta"
ref=read_fasta(REF)[0][1]; lo,hi=100,628
S=["WT2","WT3","WT4","KO1","KO2","KO3"]
p=load_preset("jh4_intron", phase_mode="off")
al,tn,cn=build_index([REF],[CON])
data={}
for s in S:
    reads,_=pipeline.prepare_reads(E+s+".fastq",REF,[CON],p)
    reads=reads[:500]
    strand={}
    seqs={n:sq for n,sq,q in read_fastq(E+s+".fastq")}
    for r in reads:
        h=[x for x in al.map(seqs[r.qname]) if x.is_primary]
        strand[r.qname]=h[0].strand if h else 0
    data[s]=(reads,strand)
# leave-one-out nulls for recommended
res={}
for s in S:
    reads,strand=data[s]
    ctrl=[c for c in ["KO1","KO2","KO3"] if c!=s]
    rates,_=null_model.fit_aligned([data[c][0] for c in ctrl],ref,p,lo,hi)
    rec=pipeline.analyse(reads,REF,p,{},null_rates=rates,depth=500)
    res[s]={(v.pos-1,v.alt) for v in rec["variants"]}
out={}
gc_sites=[i for i in range(lo,hi+1) if ref[i] in "GC"]; at_sites=[i for i in range(lo,hi+1) if ref[i] in "AT"]
def hot(pos):
    b=ref[pos]
    if b in "GC": return motifs.is_site("WRC",ref,pos) or motifs.is_site("GYW",ref,pos)
    return motifs.is_site("WA",ref,pos) or motifs.is_site("TW",ref,pos)
exp_gc=sum(hot(i) for i in gc_sites)/len(gc_sites); exp_at=sum(hot(i) for i in at_sites)/len(at_sites)
print("expected hotspot fraction if random: G:C %.2f  A:T %.2f"%(exp_gc,exp_at))
cand_by={}
for s in S:
    reads,strand=data[s]
    cnt=defaultdict(int); qs=defaultdict(float); dep=defaultdict(int); plus=defaultdict(int); depplus=defaultdict(int)
    carriers=defaultdict(list)
    readQ={r.qname: st.mean(q for _,_,q in r.pairs) for r in reads}
    for r in reads:
        for pos,b,q in r.pairs:
            if lo<=pos<=hi and b in "ACGT":
                dep[pos]+=1; depplus[pos]+= strand[r.qname]==1
                if b!=ref[pos]:
                    cnt[(pos,b)]+=1; qs[(pos,b)]+=q; plus[(pos,b)]+= strand[r.qname]==1; carriers[(pos,b)].append(r.qname)
    cands={}
    for k,c in cnt.items():
        if c/dep[k[0]]>=0.5: continue
        mq=qs[k]/c
        cands[k]=dict(k=c,q=mq,g=_geneious_pvalue(c,dep[k[0]],mq)<=0.05,r=k in res[s],plus=plus[k]/c,bgplus=depplus[k[0]]/dep[k[0]],rq=st.mean(readQ[n] for n in carriers[k]))
    cand_by[s]=cands
    allQ=[q for r in reads for _,_,q in r.pairs]
    out[s]=dict(read_meanQ=round(st.median(readQ.values()),1), base_meanQ=round(st.mean(allQ),1), plus_frac=round(sum(v==1 for v in strand.values())/len(strand),2))
for s in S: print(s,out[s])
def summ(lst,label):
    if not lst: print(f"  {label}: n=0"); return
    gc=[k for k in lst if ref[k[0]] in "GC"]; at=[k for k in lst if ref[k[0]] in "AT"]
    fg=sum(hot(k[0]) for k in gc)/len(gc) if gc else float('nan'); fa=sum(hot(k[0]) for k in at)/len(at) if at else float('nan')
    print(f"  {label}: n={len(lst)}  G:C at WRC/GYW {fg:.2f} (n={len(gc)})  A:T at WA/TW {fa:.2f} (n={len(at)})  G:C share {len(gc)/len(lst):.2f}")
print("\nHOTSPOT TARGETING")
pool=defaultdict(list)
for s in S:
    grp="WT" if s.startswith("WT") else "KO"
    c=cand_by[s]
    for k,v in c.items():
        pool[grp+" all candidates"].append(k)
        pool[grp+(" geneious kept" if v["g"] else " geneious rejected")].append(k)
        pool[grp+(" recommended kept" if v["r"] else " recommended rejected")].append(k)
        if v["k"]==1: pool[grp+(" single-read kept(geneious)" if v["g"] else " single-read rejected(geneious)")].append(k)
for lab in sorted(pool): summ(pool[lab],lab)
print("\nRECURRENCE at KO error-prone slots (candidate in >=2 of 3 KO)")
kc=defaultdict(int)
for s in ["KO1","KO2","KO3"]:
    for k in cand_by[s]: kc[k]+=1
hotslots={k for k,n in kc.items() if n>=2}
nslots=sum(3 for i in range(lo,hi+1))
print("  KO error-prone slots:",len(hotslots),"of",nslots, f"({len(hotslots)/nslots:.3f})")
for s in ["WT2","WT3","WT4"]:
    c=cand_by[s]
    for lab,sel in [("geneious kept",[k for k,v in c.items() if v["g"]]),("recommended kept",[k for k,v in c.items() if v["r"]]),("rejected",[k for k,v in c.items() if not v["g"] and not v["r"]])]:
        print(f"  {s} {lab}: {sum(k in hotslots for k in sel)}/{len(sel)} = {sum(k in hotslots for k in sel)/max(1,len(sel)):.2f}")
print("\nSINGLE-READ candidates: alt-base quality and carrying-read quality")
for s in S:
    sr=[v for v in cand_by[s].values() if v["k"]==1]
    qs=[v["q"] for v in sr]
    print(f"  {s}: n={len(sr)} median altQ={st.median(qs):.0f}  frac>=Q40={sum(q>=40 for q in qs)/len(qs):.2f}  median carrying-read meanQ={st.median(v['rq'] for v in sr):.1f}")
print("\nSTRAND: multi-read (>=3) calls, |plus frac - column plus frac|")
for s in S:
    mr=[v for v in cand_by[s].values() if v["k"]>=3]
    for lab,sel in [("kept",[v for v in mr if v["g"] or v["r"]]),("rejected",[v for v in mr if not(v["g"] or v["r"])])]:
        if sel:
            d=[abs(v["plus"]-v["bgplus"]) for v in sel]; ext=sum(v["plus"] in (0.0,1.0) for v in sel)
            print(f"  {s} {lab}: n={len(sel)} mean|dev|={st.mean(d):.2f} all-one-strand={ext}")
print("\nDISPUTED: geneious kept but recommended rejected, and vice versa")
for grp in ["WT","KO"]:
    a=[];b=[];both=[]
    for s in S:
        if not s.startswith(grp): continue
        for k,v in cand_by[s].items():
            if v["g"] and not v["r"]: a.append(k)
            if v["r"] and not v["g"]: b.append(k)
            if v["r"] and v["g"]: both.append(k)
    summ(both,grp+" kept by both"); summ(a,grp+" geneious-only"); summ(b,grp+" recommended-only")
    if grp=="WT":
        ks=[cand_by[s][k]["k"] for s in S if s.startswith("WT") for k,v in cand_by[s].items() if v["g"] and not v["r"]]
        from collections import Counter; print("   geneious-only support:",sorted(Counter(min(x,5) for x in ks).items()))
