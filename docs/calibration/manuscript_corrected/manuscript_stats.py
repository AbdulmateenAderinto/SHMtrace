"""Statistics quoted in the manuscript (Ver17), from sanger_calls/corrected_sanger.json and the
ONT matrices in scripts/figure3_concordance.py. Run from this folder: python manuscript_stats.py

Spectrum: pooled Pearson r with Fisher CI, slope; two-way repeated-measures platform x type
interaction; paired t per type with Sidak adjustment; G:C Ts/Tv paired t (Sidak, 2 tests).
Density: Welch t on log10, group geometric means with 20,000-resample bootstrap CIs.
"""
import json, numpy as np, math, random
from scipy import stats
cs={k:v for k,v in json.load(open("sanger_calls/corrected_sanger.json")).items() if isinstance(v,dict)}
src=open("scripts/figure3_concordance.py").read()
ns={"np":np}
exec('BASES=["A","T","C","G"]; IDX={b:i for i,b in enumerate(BASES)}\ndef m(r): return np.array(r,float)', ns)
i=src.index('ont={"WT2"'); j=src.index('S=["WT2"')
exec(src[i:j], ns); ont=ns["ont"]; IDX=ns["IDX"]
MUT=["A>T","A>C","A>G","T>A","T>C","T>G","C>A","C>T","C>G","G>A","G>T","G>C"]
def mat(vs):
    M=np.zeros((4,4))
    for _p,f,t,_n in vs: M[IDX[t]][IDX[f]]+=1
    return M
S=["WT2","WT3","WT4"]
SM=[mat(cs[f"Wt {k[2]}"]["variants"]) for k in S]; OM=[ont[k] for k in S]
c=lambda M,f,t: M[IDX[t]][IDX[f]]
def props(M): tot=M.sum(); return np.array([100*c(M,*t.split(">"))/tot for t in MUT])
ps,po=props(sum(SM)),props(sum(OM))
r,p=stats.pearsonr(ps,po); z=math.atanh(r); se=1/math.sqrt(9)
print(f"Pearson r {r:.5f} [{math.tanh(z-1.96*se):.2f}-{math.tanh(z+1.96*se):.2f}] slope {np.polyfit(ps,po,1)[0]:.2f}; mean|diff| {np.mean(abs(ps-po)):.2f} max {np.max(abs(ps-po)):.2f}")
Ps=np.array([props(M) for M in SM]); Po=np.array([props(M) for M in OM])
y=np.stack([Ps,Po],axis=1); Sn=3; g=y.mean()
m_pt,m_ps,m_ts=y.mean(0),y.mean(2),y.mean(1); m_p,m_t,m_s=y.mean((0,2)),y.mean((0,1)),y.mean((1,2))
ss_pt=Sn*((m_pt-m_p[:,None]-m_t[None,:]+g)**2).sum()
res=y-m_pt[None]-m_ps[:,:,None]-m_ts[:,None,:]+m_p[None,:,None]+m_t[None,None,:]+m_s[:,None,None]-g
F=(ss_pt/11)/((res**2).sum()/22); print(f"interaction F(11,22)={F:.2f} p={stats.f.sf(F,11,22):.3f}")
adj=[1-(1-stats.ttest_rel(Ps[:,i],Po[:,i]).pvalue)**12 for i in range(12)]
print("min Sidak",round(min(adj),3),MUT[int(np.argmin(adj))])
def gctstv(M):
    tot=M.sum(); return 100*(c(M,"G","A")+c(M,"C","T"))/tot, 100*(c(M,"G","T")+c(M,"G","C")+c(M,"C","A")+c(M,"C","G"))/tot
st=np.array([gctstv(M) for M in SM]); ot=np.array([gctstv(M) for M in OM])
print("Ts/Tv adj p",[round(1-(1-stats.ttest_rel(st[:,k],ot[:,k]).pvalue)**2,3) for k in (0,1)])
dens={k:1e5*max(v["unique"],1)/v["bases"] for k,v in cs.items()}
wt=[dens[f"Wt {i}"] for i in (2,3,4)]; ko=[dens[f"Aicda−/− {i}"] for i in (1,2,3)]
geo=lambda x: 10**np.mean(np.log10(x)); rng=random.Random(1)
def ci(x,n=20000):
    v=sorted(geo([rng.choice(x) for _ in x]) for _ in range(n)); return v[int(.025*n)],v[int(.975*n)-1]
print("densities",{k:round(v,1) for k,v in dens.items()})
print(f"Wt geo {geo(wt):.0f} {[round(x) for x in ci(wt)]}  KO geo {geo(ko):.1f} {[round(x,1) for x in ci(ko)]}  fold {geo(wt)/geo(ko):.1f}  P {stats.ttest_ind(np.log10(wt),np.log10(ko),equal_var=False).pvalue:.4f}")
loads=[x for k,v in cs.items() if k.startswith("Wt") for x in v["loads"]]
print("pie n",len(loads),[sum(1 for x in loads if a<=x<=b) for a,b in ((0,0),(1,3),(4,10),(11,20),(21,10**6))])
print("Sanger unique pooled",int(sum(SM).sum()),"ONT",int(sum(OM).sum()),"| clones",{k:v["clones"] for k,v in cs.items()},"| effective",{k:round(v["bases"]/529,1) for k,v in cs.items()})
