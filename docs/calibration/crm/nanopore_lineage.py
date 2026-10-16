"""Which mouse a Sanger clone's lineage appears in, checked against the nanopore libraries.

    python docs/calibration/crm/nanopore_lineage.py <fastq-dir> [sheet]

The J_H4 nanopore run covers the same mice as the Sanger studies, so a clone that carries several
mutations should be found again among the reads of its own mouse: somatic mutations are not
inherited, so no other mouse carries the same combination except by coincidence. Each clone's
mutations (germline removed, as the Sanger run removes them) are looked for in every read of every
library; a library is reported when its reads carry THRESHOLD or more of them, which pairs of
mutations at the intron's hotspots do not reach.

Writes out/nanopore_lineage.json: per clone, the reads of each library sharing 3 or more of its
mutations. Clones with fewer than THRESHOLD mutations cannot be placed and are reported as such.
"""
import collections
import json
import re
import sys
from pathlib import Path

import mappy

CRM = Path(__file__).resolve().parent
REF = CRM.parents[2] / "src/shmtrace/references/mouse_jh4_intron.fasta"
OUT = CRM / "out/nanopore_lineage.json"
LO, HI = 100, 628
MIN_Q = 15            # a mismatch in a read counts only at this base quality
THRESHOLD = 3         # shared mutations before a read is called the same lineage
COMP = str.maketrans("ACGTN", "TGCAN")


def read_subs(fastq, ref, aligner):
    """Every read's substitutions against the reference, over the Sanger window."""
    out = []
    for _name, seq, qual in mappy.fastx_read(str(fastq)):
        h = next(iter(aligner.map(seq, cs=True)), None)
        if h is None or h.mlen / h.blen < 0.85 or h.r_st > LO or h.r_en < HI:
            continue
        if h.strand == -1:
            seq, qual = seq.translate(COMP)[::-1], qual[::-1]
            h = next(iter(aligner.map(seq, cs=True)), None)
            if h is None or h.strand != 1:
                continue
        r, q, subs = h.r_st, h.q_st, {}
        for op, val in re.findall(r"([:*+\-~])([0-9]+|[a-z]+)", h.cs):
            if op == ":":
                r += int(val)
                q += int(val)
            elif op == "*":
                if ord(qual[q]) - 33 >= MIN_Q:
                    subs[r] = val[1].upper()
                r += 1
                q += 1
            elif op == "+":
                q += len(val)
            elif op == "-":
                r += len(val)
        out.append(subs)
    return out


def main(fastq_dir, clones):
    """clones: {clone id: (sample, [(pos, base), ...])}, germline already removed."""
    ref = "".join(ln.strip() for ln in open(REF) if not ln.startswith(">")).upper()
    aligner = mappy.Aligner(seq=ref, preset="map-ont")
    libs = {}
    for f in sorted(Path(fastq_dir).glob("*.fastq")):
        name = re.sub(r"^[A-Z0-9]+_\d+_", "", f.stem)
        libs[name] = read_subs(f, ref, aligner)
        print(f"{name}: {len(libs[name])} reads spanning {LO}-{HI}")
    rows = []
    for cid, (sample, muts) in sorted(clones.items()):
        row = {"clone": cid, "sample": sample, "mutations": len(muts)}
        if len(muts) < THRESHOLD:
            row["verdict"] = "too few mutations to place"
        else:
            hits = {}
            for lib, reads in libs.items():
                shared = collections.Counter(sum(1 for p, b in muts if r.get(p) == b) for r in reads)
                n = sum(v for k, v in shared.items() if k >= THRESHOLD)
                if n:
                    hits[lib] = {"reads": n, "best": max(shared), "of": len(reads)}
            row["libraries"] = hits
            row["verdict"] = max(hits, key=lambda k: hits[k]["best"]) if hits else "not found in any library"
        rows.append(row)
    OUT.parent.mkdir(exist_ok=True)
    json.dump(rows, OUT.open("w"), indent=1)
    placed = sum(1 for r in rows if r["verdict"] in libs)
    print(f"{placed} of {len(rows)} clones placed; out/{OUT.name}")
    return rows


if __name__ == "__main__":
    main(sys.argv[1], json.load(open(sys.argv[2])) if len(sys.argv) > 2 else {})
