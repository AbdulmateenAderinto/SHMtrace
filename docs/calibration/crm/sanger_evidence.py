"""The Mkln1 J_H4 clones pruned and moved on sequence evidence alone: a sensitivity analysis.

    python docs/calibration/crm/sanger_evidence.py <nanopore-fastq-dir>

The primary analysis (sanger_cohort2.py mkln) takes every clone as the original project filed it. This
one asks whether the result depends on clones whose filing the sequences contradict. Its rules were
fixed before it was run, and none of them refers to the published numbers:

1. Paralogue amplicons are set aside (SHMtrace's competitive alignment, as in the primary analysis).
2. Byte-identical copies of a trace count once (the sample sheet, as in the primary analysis).
3. A lineage filed under two mice (SHMtrace's duplicate_clone and shared_mutations checks) is tested
   against the J_H4 nanopore libraries of the same mice. Its shared mutations are looked for in every
   read, leaving out positions where the Aicda-/- reads themselves mismatch (more than 0.5% of them),
   which are nanopore errors and not markers. A read belongs to the lineage when it carries more than
   half of those mutations and at least three. If exactly one library has three or more such reads,
   and it is a library of a mouse in this study, every clone of the lineage moves to that mouse, even
   one it was not filed under. If no library does, or more than one, or the one is an Aicda-/- or
   another study's library, the lineage's clones leave the per-mouse counts.
4. A clone filed under an Aicda-/- mouse that SHMtrace flags as carrying mutations beyond its base
   qualities (control_mutations) is set aside: an AID-null mouse cannot have made it.
5. Clones the original project never filed into a mouse folder stay out, as in the primary analysis.

Writes samplesheets/jh4_mkln_barbulescu2025_evidence.csv (each change recorded in the clone's note,
each removal in its exclude column) and out/sanger_mkln_evidence.json, and prints each decision.
"""
import collections
import csv
import json
import sys
from pathlib import Path

import mappy

CRM = Path(__file__).resolve().parent
sys.path.insert(0, str(CRM))
import nanopore_lineage as nl  # noqa: E402
import sanger_cohort2 as sc  # noqa: E402

CAL = CRM.parent
SHEET = CAL / sc.STUDIES["mkln"]
OUT_SHEET = CAL / "samplesheets/jh4_mkln_barbulescu2025_evidence.csv"
LIBRARY_OF = {"MKLN_KO_2": "MKLN_K2", "MKLN_KO_3": "MKLN_K3", "MKLN_KO_4": "MKLN_K4",
              "Wildtype_2": "WT2", "Wildtype_3": "WT3", "Wildtype_4": "WT4"}
AID_NULL = ("AID_KO_1", "AID_KO_2", "AID_KO_3")
OTHER = ("MAEA_KO_1",)          # another study's mouse: a claim from it makes a lineage ambiguous
NOISE = 0.005          # a position mismatched in more of the Aicda-/- reads than this is not a marker
MIN_SHARED = 3         # mutations a read must carry, at the least
MIN_READS = 3          # reads a library must have before it claims the lineage


def lineages(qc, sample_of):
    """Clones filed under two or more mice that SHMtrace ties together, with their shared mutations."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    shared = collections.defaultdict(set)
    for r in qc["pairs"]:
        a, b = r["a"], r["b"]
        if sample_of.get(a) == sample_of.get(b):
            continue
        parent[find(a)] = find(b)
        for p, base in r.get("shared_mutations", []):
            shared[a].add((p, base))
            shared[b].add((p, base))
    comps = collections.defaultdict(list)
    for c in parent:
        comps[find(c)].append(c)
    out = []
    for cs in comps.values():
        muts = set().union(*(shared[c] for c in cs))
        out.append({"clones": sorted(cs), "samples": sorted({sample_of[c] for c in cs}), "mutations": sorted(muts)})
    return out


def judge(lineage, libs, noisy):
    """Which library, if any, the lineage belongs to, with the evidence."""
    muts = [(p - 1, b) for p, b in lineage["mutations"] if (p - 1, b) not in noisy]   # 1-based -> 0-based
    need = max(MIN_SHARED, len(muts) // 2 + 1)
    reads = {lib: sum(1 for r in rs if sum(1 for p, b in muts if r.get(p) == b) >= need) for lib, rs in libs.items()}
    claims = sorted(lib for lib, n in reads.items() if n >= MIN_READS)
    if len(muts) < MIN_SHARED:
        verdict, reason = None, f"only {len(muts)} usable mutations"
    elif len(claims) == 1 and claims[0] in LIBRARY_OF:
        verdict, reason = LIBRARY_OF[claims[0]], f"{reads[claims[0]]} reads of {claims[0]} carry {need}+ of {len(muts)}"
    elif not claims:
        verdict, reason = None, f"no library has {MIN_READS} reads carrying {need}+ of {len(muts)}"
    else:
        verdict, reason = None, "claimed by " + ", ".join(f"{c} ({reads[c]})" for c in claims)
    return verdict, reason, {k: v for k, v in reads.items() if v}


def main(fastq_dir):
    asfiled = json.load(open(CRM / "out/sanger_mkln.json"))["project"]
    rows = list(csv.DictReader(open(SHEET)))
    header = list(rows[0])
    sample_of = {r["clone_id"]: r["sample"] for r in rows}

    ref = "".join(ln.strip() for ln in open(nl.REF) if not ln.startswith(">")).upper()
    aligner = mappy.Aligner(seq=ref, preset="map-ont")
    libs = {}
    for f in sorted(Path(fastq_dir).glob("*.fastq")):
        name = f.stem.split("_", 2)[2]
        if name in LIBRARY_OF or name in AID_NULL or name in OTHER:
            libs[name] = nl.read_subs(f, ref, aligner)
    n_null = sum(len(libs[k]) for k in AID_NULL)
    null_count = collections.Counter(m for k in AID_NULL for r in libs[k] for m in r.items())
    noisy = {m for m, n in null_count.items() if n / n_null > NOISE}

    moves, drops = {}, {}
    for lin in lineages(asfiled["qc"], sample_of):
        verdict, reason, reads = judge(lin, libs, noisy)
        print(f"lineage {lin['clones']} filed under {lin['samples']}: "
              f"{'-> ' + verdict if verdict else 'left out'} ({reason}; reads {reads})")
        for c in lin["clones"]:
            if verdict:
                if sample_of[c] != verdict:
                    moves[c] = (verdict, f"moved from {sample_of[c]}: lineage of {len(lin['clones'])} clones; {reason}")
            else:
                drops[c] = f"lineage filed under {' and '.join(lin['samples'])}; {reason}"
    for f in asfiled["qc"]["flags"]:
        if f["check"] == "control_mutations":
            drops[f["clone"]] = "mutated clone filed under an Aicda-/- mouse: " + f["detail"]
            print(f"control clone {f['clone']} ({f['sample']}) left out")

    group_of = {r["sample"]: (r["group"], r["control"]) for r in rows}
    for r in rows:
        c = r["clone_id"]
        if c in moves:
            r["sample"], why = moves[c]
            r["group"], r["control"] = group_of[r["sample"]]       # a clone takes its new mouse's genotype
            r["note"] = (r["note"] + "; " if r["note"] else "") + "evidence set: " + why
        elif c in drops and not r["exclude"]:
            r["exclude"] = "evidence set: " + drops[c]
    with open(OUT_SHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, header)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(moves)} clones moved, {len(drops)} left out; {OUT_SHEET.name}")
    sc.STUDIES["mkln_evidence"] = str(OUT_SHEET.relative_to(CAL))
    sc.run("mkln_evidence")


if __name__ == "__main__":
    main(sys.argv[1])
