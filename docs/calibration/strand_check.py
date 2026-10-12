"""Strand support of shmtrace calls on the deposited JH4 sub500 FASTQs (mixed orientation)."""
import csv, mappy
from collections import defaultdict
from dataclasses import replace
from scipy.stats import binomtest
from shmtrace import pipeline, config
rows=list(csv.DictReader(open("jh4_sheet.csv")))
REF="deposit_Ver9/02_JH4_ONT/Reference_Sequences/jh4-base.fasta"; CON=["deposit_Ver9/02_JH4_ONT/Reference_Sequences/Contaminant_JH4_ref.fasta"]
al=mappy.Aligner(REF,preset="map-ont")
p=replace(config.load_preset("jh4_intron"),**pipeline.mode_overrides("geneious"))
for r in rows:
    strand={n:(next(al.map(s),None) or type("x",(),{"strand":0})).strand for n,s,q in mappy.fastx_read(r["fastq"])}
    res=pipeline.run(r["fastq"],REF,CON,p)
    reads=res["_aligned_all"]
    sup=defaultdict(lambda:[0,0]); dep=defaultdict(lambda:[0,0])
    for a in reads:
        sd=0 if strand[a.qname]==1 else 1
        for rp,b,q in a.pairs:
            dep[rp+1][sd]+=1; sup[(rp+1,b)][sd]+=1
    multi=biased=0; single=[0,0]
    for v in res["variants"]:
        pl,mi=sup[(v.pos,v.alt)]; dp,dm=dep[v.pos]
        if pl+mi==1: single[0 if pl else 1]+=1
        else:
            multi+=1
            if binomtest(pl,pl+mi,dp/(dp+dm)).pvalue<0.01: biased+=1
    fr=sum(strand.values().__iter__().__class__.__mro__ and [1 for x in strand.values() if x==1])/len(strand)
    print(f"{r['sample']}: calls {len(res['variants'])}, plus-strand reads {fr:.2f}, multi-read calls {multi} (strand-biased p<0.01: {biased}), single-read calls +/-: {single[0]}/{single[1]}")
