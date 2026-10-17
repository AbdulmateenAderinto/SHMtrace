"""How SHMtrace merges Illumina read pairs, measured on two of the public libraries, for Figure S12.

    python docs/calibration/crm/scicanvas/s12_illumina.py   -> crm/out/s12_illumina.json

The runs record how many reads were analyzed but not how the pairs merged, so this re-merges the first 20,000
pairs of one wild-type J_H4 mouse (PRJNA973112) and one AID-expressing CH12F3 replicate (PRJNA918596) with
io_utils.merge_pair, the function the pipeline uses, and reports the share merged, the overlap length, and how
often the two reads disagree inside the overlap. It also reports the
distinct base qualities in the reads: these two SRA files carry a flat Q30 on every base.
"""
import gzip
import json
import statistics
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CRM.parents[2] / "src"))
from shmtrace import io_utils  # noqa: E402

PUB = Path.home() / "Nanoshm/public_data"
N = 20_000
LIBS = {"jh4_wt": (PUB / "PRJNA973112_sub", None), "ch12_aid": (PUB / "PRJNA918596_ch12/pooled", "AID_r1")}


def fastq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip()
            fh.readline()
            yield h.strip(), s, fh.readline().strip()


def paths(key):
    d, stem = LIBS[key]
    if stem is None:   # the first wild-type mouse in the sample list
        stem = next(r for m, r in (ln.split() for ln in open(d / "mouse_runs.txt")) if "_WT" in m)
    return d / f"{stem}_1.fastq.gz", d / f"{stem}_2.fastq.gz", stem


def measure(key):
    r1p, r2p, stem = paths(key)
    n = merged = 0
    ov, dis, l1, l2, lm, quals = [], [], [], [], [], set()
    for (_h, s1, q1), (_h2, s2, q2) in zip(fastq(r1p), fastq(r2p), strict=False):
        n += 1
        o = io_utils._best_offset(s1.upper(), io_utils._revcomp(s2.upper()), 20, 0.1, 12)
        m = io_utils.merge_pair(s1, q1, s2, q2)
        if m is not None:
            merged += 1
            a, b = max(0, o), min(len(s1), o + len(s2))
            rc2 = io_utils._revcomp(s2.upper())
            d = sum(x != y for x, y in zip(s1.upper()[a:b], rc2[a - o:b - o], strict=True))
            ov.append(b - a)
            dis.append(d / (b - a))
            l1.append(len(s1))
            l2.append(len(s2))
            lm.append(len(m[0]))
            quals.update(q1)
        if n == N:
            break
    return {"library": stem, "pairs": n, "merged_pct": 100 * merged / n,
            "r1_len_median": statistics.median(l1), "r2_len_median": statistics.median(l2),
            "merged_len_median": statistics.median(lm), "overlap_median": statistics.median(ov),
            "overlap_iqr": statistics.quantiles(ov, n=4)[::2],
            "disagree_per_overlap_base_pct": 100 * statistics.fmean(dis),
            "base_qualities_r1": sorted(ord(c) - 33 for c in quals)}


def main():
    out = {k: measure(k) for k in LIBS}
    (CRM / "out").mkdir(exist_ok=True)
    json.dump(out, open(CRM / "out/s12_illumina.json", "w"), indent=1)
    for k, v in out.items():
        print(k, {x: (round(y, 2) if isinstance(y, float) else y) for x, y in v.items()})


if __name__ == "__main__":
    main()
