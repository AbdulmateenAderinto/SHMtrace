"""Build per-position, per-alt pileups for every sample at full depth."""
import json, os, pickle, sys
sys.path.insert(0, "/sessions/festive-dreamy-davinci/mnt/shmtrace/src")
from collections import defaultdict
from shmtrace.config import load_preset
from shmtrace.io_utils import read_fastq, read_fasta
from shmtrace.reads_filter import prefilter
from shmtrace.align import build_index, align_reads

D = "/sessions/festive-dreamy-davinci/mnt/June 2nd analyses/submission_v4/ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ"
REF = "/sessions/festive-dreamy-davinci/mnt/shmtrace/examples_compare/JH4_reference.fasta"
CON = "/sessions/festive-dreamy-davinci/mnt/shmtrace/examples_compare/JH4_contaminant.fasta"
SAMPLES = {
    "WT2": "DHJGKB_2_Wildtype_2", "WT3": "DHJGKB_3_Wildtype_3", "WT4": "DHJGKB_4_Wildtype_4",
    "KO1": "DHJGKB_5_AID_KO_1",   "KO2": "DHJGKB_6_AID_KO_2",   "KO3": "DHJGKB_7_AID_KO_3",
}
BASES = ("A","C","G","T")
p = load_preset("jh4_intron")
p.downsample_to = None                     # full depth: more power to see the error structure
ref = dict(read_fasta(REF))["JH4_C57BL6"]

out = {}
for tag, stem in SAMPLES.items():
    reads = list(read_fastq(f"{D}/{stem}.fastq"))
    kept, s = prefilter(reads, p)
    aligner, tnames, cnames = build_index([REF], [CON])
    aligned = [a for a in align_reads(kept, aligner, tnames, cnames, p) if a.ctg in tnames]
    counts = defaultdict(lambda: defaultdict(int)); depth = defaultdict(int)
    qsum   = defaultdict(float)
    for r in aligned:
        for rpos, base, q in r.pairs:
            if base not in BASES: continue
            counts[rpos][base] += 1; depth[rpos] += 1
            qsum[rpos] += 10 ** (-q/10.0)
    out[tag] = {
        "counts": {k: dict(v) for k, v in counts.items()},
        "depth": dict(depth), "qsum": dict(qsum),
        "reads_in": s["reads_in"], "passed": s["passed_prefilter"], "aligned": len(aligned),
    }
    print(f"{tag}: in={s['reads_in']} passed={s['passed_prefilter']} aligned={len(aligned)}", flush=True)

pickle.dump({"ref": ref, "samples": out}, open("/sessions/festive-dreamy-davinci/mnt/outputs/calib/pileups.pkl","wb"))
print("saved pileups.pkl")
