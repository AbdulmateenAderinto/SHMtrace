"""Class-switch junctions in a public PacBio dataset (PRJNA831666), analysed with shmtrace.switch.

    python crm/switch_pacbio.py            -> crm/out/switch_pacbio/

Vincendeau et al. 2022 (Nat Commun 13:3707) sequenced Sμ-Sγ1 switch products from primary B cells of
wild-type, Shld1-knockout, Xlf-knockout and double-knockout mice, three biological replicates each. The study
compares repair genotypes against wild type; it has no AID-null arm, so it tests the junction module and not
the substitution caller's control calibration.

Read accounting matters here and is reported per sample: the libraries carry a large population of short
off-target molecules (median about 350 bp) alongside the switch products (median about 1.3 kb), so the
fraction of reads that yield a junction is not a switching efficiency. Neither is the junction count: long-range
PCR and sequencing both condition on recoverable products. The unit of analysis is the biological sample, and
the quantities compared are shapes of the junctions (microhomology classes, mean microhomology, mutations near
the join), not their number.

Reads are high-accuracy despite the archive's `_subreads` file names: median identity against the bundled Sμ
fragment is 99.8% over alignments of 600 bp or more, so `map-hifi` is used rather than the module's `map-ont`
default.

`switch.py` reports an `ambiguous` flag for junctions whose acceptor part also aligns elsewhere. On these data
minimap2 returns no secondary alignments for junction reads, so that flag is always False and its 0% must not
be read as "breakpoints are unambiguous". This script measures the same thing directly instead, by re-aligning
each junction's acceptor-side segment on its own and counting near-equal placements.
"""
import csv
import json
import statistics
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
import mappy as mp  # noqa: E402

from shmtrace import switch  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

DATA = Path.home() / "Nanoshm/public_data/PRJNA831666"
REFDIR = CRM / "out/switch_ref"
OUT = CRM / "out/switch_pacbio"
PRESET = "map-hifi"
# PRJNA831666: run -> (genotype, replicate). Genotypes are the submitters' own sample descriptions.
SAMPLES = [("SRR18913496", "WT", 1), ("SRR18913495", "WT", 2), ("SRR18913492", "WT", 3),
           ("SRR18913491", "Shld1", 1), ("SRR18913490", "Shld1", 2), ("SRR18913489", "Shld1", 3),
           ("SRR18913488", "Xlf", 1), ("SRR18913487", "Xlf", 2), ("SRR18913486", "Xlf", 3),
           ("SRR18913485", "Shld1_Xlf", 1), ("SRR18913494", "Shld1_Xlf", 2), ("SRR18913493", "Shld1_Xlf", 3)]


def acceptor_ambiguity(junctions, seqs, acc_ref, seg=60, slack=5):
    """Share of junctions whose acceptor-side segment has more than one near-equal placement.

    switch.py's own flag depends on minimap2 reporting secondary alignments, which it does not do here. This
    aligns the segment on its own, with secondary alignments kept, and counts placements scoring within
    `slack` matched bases of the best one.
    """
    # a short query needs the short-read preset: map-hifi returns nothing for a 60-bp segment
    al = mp.Aligner(seq=acc_ref, preset="sr", best_n=50)
    amb = n = 0
    for j in junctions:
        s = seqs.get(j.read)
        if not s:
            continue
        # the read was oriented to the donor inside call_read; use the stored orientation
        start = s[1]
        seg_seq = s[0][start:start + seg]
        if len(seg_seq) < seg:
            continue
        hits = sorted((h.mlen for h in al.map(seg_seq)), reverse=True)
        if not hits:
            continue
        n += 1
        amb += sum(1 for m in hits if m >= hits[0] - slack) > 1
    return {"checked": n, "ambiguous": amb, "ambiguous_pct": round(100 * amb / n, 1) if n else None}


def analyse_sample(run, ref_path, refs, donor):
    """switch.analyse, plus read accounting and the direct ambiguity measure."""
    path = DATA / f"{run}_subreads.fastq.gz"
    al = mp.Aligner(str(ref_path), preset=PRESET, best_n=10)
    lens_all, lens_junc, junctions, seqs = [], [], [], {}
    mapped = 0
    for name, seq, qual in read_fastq(str(path)):
        lens_all.append(len(seq))
        hits = [h for h in al.map(seq) if h.q_en - h.q_st >= 60]
        if hits:
            mapped += 1
        j = switch.call_read(name, seq, qual, al, refs, donor)
        if j:
            junctions.append(j)
            lens_junc.append(len(seq))
            d0 = max((h for h in al.map(seq) if h.ctg == donor), key=lambda h: h.mlen, default=None)
            oriented = mp.revcomp(seq) if (d0 and d0.strand == -1) else seq
            # where the acceptor part starts in the oriented read, for the ambiguity probe
            acc = [h for h in al.map(oriented) if h.ctg != donor and h.strand == 1]
            if acc:
                seqs[name] = (oriented, max(acc, key=lambda h: h.mlen).q_st)
    summary = switch.summarise(junctions)
    uniq = {j.key: j for j in switch.unique_junctions([j for j in junctions if j.kind == "switch"])}
    return {"run": run, "reads": len(lens_all), "reads_mapped_to_switch_refs": mapped,
            "read_len_median_all": statistics.median(lens_all) if lens_all else None,
            "read_len_median_junction": statistics.median(lens_junc) if lens_junc else None,
            "reads_per_unique_junction": round(summary["junction_reads"] / summary["unique_junctions"], 2)
            if summary["unique_junctions"] else None,
            **summary,
            "acceptor_ambiguity_direct": acceptor_ambiguity(list(uniq.values()), seqs, refs["Sg1"])}, junctions


def main():
    ref_path = REFDIR / "mouse_Smu_Sg1.fasta"
    if not ref_path.exists():
        sys.exit("run crm/build_switch_reference.py first")
    recs = read_fasta(str(ref_path))
    refs = {n.split()[0]: s.upper() for n, s in recs}
    donor = "Smu"
    OUT.mkdir(parents=True, exist_ok=True)
    per_sample, groups = {}, {}
    for run, geno, rep in SAMPLES:
        name = f"{geno}_{rep}"
        res, junctions = analyse_sample(run, ref_path, refs, donor)
        per_sample[name] = res
        groups[name] = geno
        if junctions:
            switch.write_csv(junctions, str(OUT / f"junctions_{name}.csv"))
        a = res["acceptor_ambiguity_direct"]
        print(f"{name:12s} {run} reads {res['reads']:>6} mapped {res['reads_mapped_to_switch_refs']:>5} "
              f"junction reads {res['junction_reads']:>5} unique {res['unique_junctions']:>4} "
              f"reads/junction {res['reads_per_unique_junction']} mean MH {res['mean_microhomology']} "
              f"ambiguous(direct) {a['ambiguous_pct']}%")
    comparisons = []
    for case in ("Shld1", "Xlf", "Shld1_Xlf"):
        comparisons += switch.compare_groups(per_sample, groups, case, "WT")
    json.dump({"dataset": "PRJNA831666 (Vincendeau et al. 2022, Nat Commun 13:3707)",
               "reference": json.load(open(REFDIR / "reference_provenance.json")),
               "preset": PRESET, "unit": "biological sample; unique junctions within a sample",
               "samples": per_sample, "groups": groups, "comparisons": comparisons},
              open(OUT / "switch_pacbio.json", "w"), indent=1)
    with open(OUT / "comparisons.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comparisons[0]))
        w.writeheader()
        w.writerows(comparisons)
    print()
    for r in comparisons:
        if r.get("p_welch") is not None:
            print(f"{r['group']:10s} vs WT  {r['feature']:28s} {r['mean_group']:>7} vs {r['mean_control']:>7} "
                  f"Welch p={r['p_welch']:.3f} perm p={r['p_permutation']:.3f}")


if __name__ == "__main__":
    main()
