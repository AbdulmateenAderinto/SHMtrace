import sys
sys.path.insert(0,"/sessions/festive-dreamy-davinci/mnt/nanoshm/src")
import numpy as np
from collections import defaultdict
from nanoshm.config import load_preset
from nanoshm.io_utils import read_fastq, read_fasta
from nanoshm.reads_filter import prefilter
from nanoshm.align import build_index, align_reads
D="/sessions/festive-dreamy-davinci/mnt/June 2nd analyses/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ"
REF="/sessions/festive-dreamy-davinci/mnt/nanoshm/examples_compare/JH4_reference.fasta"
CON="/sessions/festive-dreamy-davinci/mnt/nanoshm/examples_compare/JH4_contaminant.fasta"
p=load_preset("jh4_intron"); p.downsample_to=None
ref=dict(read_fasta(REF))["JH4_C57BL6"]
TOP=[(216,'A'),(150,'A'),(110,'G'),(220,'C'),(160,'T'),(165,'G'),(164,'A'),(158,'G'),
     (145,'G'),(142,'A'),(276,'G'),(260,'G'),(230,'G'),(299,'A'),(133,'C'),(212,'G')]  # 0-based
for tag,stem in [("KO1","DHJGKB_5_AID_KO_1"),("WT3","DHJGKB_3_Wildtype_3")]:
    reads=list(read_fastq(f"{D}/{stem}.fastq")); kept,_=prefilter(reads,p)
    al,tn,cn=build_index([REF],[CON])
    aligned=[a for a in align_reads(kept,al,tn,cn,p) if a.ctg in tn]
    carries=defaultdict(set)     # read index -> set of top slots it carries
    covers=defaultdict(set)
    for i,r in enumerate(aligned):
        seen={rp:b for rp,b,_q in r.pairs}
        for pos,alt in TOP:
            if pos in seen:
                covers[i].add((pos,alt))
                if seen[pos]==alt: carries[i].add((pos,alt))
    n=[len(carries[i]) for i in range(len(aligned)) if len(covers[i])>=14]
    n=np.array(n)
    print(f"\n{tag}: {len(n)} reads covering >=14 of the 16 hot slots")
    print(f"   reads carrying 0 of them : {(n==0).sum():5d} ({100*(n==0).mean():.1f}%)")
    print(f"   reads carrying 1-3       : {((n>=1)&(n<=3)).sum():5d} ({100*((n>=1)&(n<=3)).mean():.1f}%)")
    print(f"   reads carrying >=12      : {(n>=12).sum():5d} ({100*(n>=12).mean():.1f}%)")
    print(f"   reads carrying all 16    : {(n==16).sum():5d} ({100*(n==16).mean():.1f}%)")
    print(f"   mean carried = {n.mean():.2f} of 16")
