"""Targeted AID mutagenesis on an EGFP reporter (Ma et al. 2016, SRP075898).

    python crm/tam_reporter.py        -> crm/out/tam/

TAM tethers a cytidine deaminase to dCas9 and aims it with sgRNAs, so the question here is not whether a
deaminase mutates - the paper already shows that - but whether the calls SHMtrace makes are *where the guide
points*. That is a different test from the untargeted reporter already in the manuscript.

The reporter arm has four libraries: a mock with no editor, AID alone, dCas9-AIDX with one sgRNA against
EGFP, and dCas9-AIDX with a pool of sgRNAs. One library per condition, so nothing here is replicated and no
group statistic is computed; the endpoint is where the calls fall along the amplicon.

Read handling needed care. The archive labels SRR3608504 (mock) SINGLE, but it exports three FASTQ files:
spots 1-809,813 are a paired block (`_1`/`_2`) and spots 809,814-1,292,558 are a single-end block (the
unsuffixed file), together matching the run's stated read count. The other three libraries are single-end
only. This uses the single-end block of the mock, so every library contributes one read of the same kind;
the paired block is a different preparation and is left out rather than mixed in.

The off-target arm of the same project is not analysed here: it is genome-wide, fragmented material that
would have to be forced through amplicon-shaped gates, and `DATASET_LEDGER.md` records that as a separate
capability rather than pretending the amplicon path covers it.
"""
import collections
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
from shmtrace import io_utils, project  # noqa: E402

DATA = Path.home() / "Nanoshm/public_data/SRP075898"
OUT = CRM / "out/tam"
DEPTH = 50_000
# run -> (name, what it is, is the no-editor control)
REPORTER = [("SRR3608504", "mock", "no editor", True),
            ("SRR3615407", "AIDx", "AID alone, untethered", False),
            ("SRR3615409", "dCas9_AIDX_sgRNA", "dCas9-AIDX with one sgRNA against EGFP", False),
            ("SRR3615429", "dCas9_AIDX_pool", "dCas9-AIDX with a pool of sgRNAs", False)]


AMPLICON_STARTS = ["ACGTCTATATCATGGCCGAC", "CCCCTTGAACCTCCTCGTTC", "CCTGACCTACGGCGTGCAGT",
                   "GCCCTATGTCCACTTCAGGA", "TTTTCTCCCATTCCTGCGCT"]


def build_references(run, out_path, per_amplicon=6000):
    """One consensus per amplicon, from the no-editor reads.

    The reporter libraries are multiplexed: five amplicons share each file, each with its own start. A single
    position-wise consensus over all reads would average different sequences together, so reads are first
    grouped by their exact 20-base start - safe in the no-editor library, where almost nothing is mutated -
    and each group is collapsed on its own.
    """
    cols = {k: collections.defaultdict(collections.Counter) for k in AMPLICON_STARTS}
    n = collections.Counter()
    for _nm, seq, _q in io_utils.read_fastq(str(DATA / f"{run}.fastq.gz")):
        k = seq[:20]
        if k not in cols or len(seq) < 250:
            continue
        if n[k] >= per_amplicon:
            continue
        for i, b in enumerate(seq[:300]):
            cols[k][i][b] += 1
        n[k] += 1
        if all(n[a] >= per_amplicon for a in AMPLICON_STARTS):
            break
    refs, weakest = {}, {}
    with open(out_path, "w") as fh:
        for a in AMPLICON_STARTS:
            c = cols[a]
            if not c:
                continue
            seq = "".join(c[i].most_common(1)[0][0] for i in sorted(c))
            refs[a] = seq
            weakest[a] = min(c[i].most_common(1)[0][1] / sum(c[i].values()) for i in c)
            fh.write(f">amp_{a}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
    return refs, dict(n), weakest


def main():
    missing = [r for r, *_ in REPORTER if not (DATA / f"{r}.fastq.gz").exists()]
    if missing:
        sys.exit(f"missing reads for {missing}")
    OUT.mkdir(parents=True, exist_ok=True)
    refs, counts, weakest = build_references("SRR3608504", OUT / "amplicons.fasta")
    print(f"{len(refs)} amplicons from the no-editor library; reads used {counts}", flush=True)
    for a, w in weakest.items():
        print(f"   amp_{a}: {len(refs[a])} bp, weakest consensus position {100 * w:.1f}%", flush=True)
    what = {name: w for _r, name, w, _c in REPORTER}
    out = {"dataset": "SRP075898 (Ma et al. 2016, Nat Methods)",
           "read_handling": "single-end block only; the mock's paired block is excluded",
           "replication": "one library per condition; no group statistic is computed",
           "amplicons": {f"amp_{a}": {"bp": len(s_), "weakest_consensus_pct": round(100 * weakest[a], 1)}
                         for a, s_ in refs.items()},
           "depth": DEPTH, "results": {}}
    # one project per amplicon; the other amplicons are competing references, which is how the pipeline
    # assigns a read to the sequence it actually came from
    for a in refs:
        single = OUT / f"ref_{a}.fasta"
        single.write_text(f">amp_{a}\n" + "\n".join(refs[a][i:i + 70] for i in range(0, len(refs[a]), 70))
                          + "\n")
        others = OUT / f"others_{a}.fasta"
        others.write_text("".join(f">amp_{b}\n{refs[b]}\n" for b in refs if b != a))
        samples = [{"name": name, "path": str(DATA / f"{run}.fastq.gz"), "group": name, "control": ctrl}
                   for run, name, _w, ctrl in REPORTER]
        cfg = {"name": f"TAM amp_{a}", "region": "generic_amplicon", "platform": "illumina",
               "reference": str(single), "contaminants": [str(others)], "depth": DEPTH, "samples": samples,
               "overrides": {"analysis_window": [10, len(refs[a]) - 10], "end_trim_bp": 0,
                             "min_align_len": 200, "max_total_len": 400, "exclude_indels": True, "jobs": 1}}
        pub = project.public(project.run_project(cfg))
        rows = {}
        for s_ in pub["samples"]:
            calls = sorted(v["pos"] for v in pub["variants"] if v["sample"] == s_["sample"])
            rows[s_["sample"]] = {"condition": what[s_["sample"]], "reads": s_["reads"],
                                  "unique_mutations": s_["unique_mutations"],
                                  "mutation_frequency_x1e3": s_["mutation_frequency_x1e3"],
                                  "pct_CT_GA": s_.get("pct_CT_GA"), "called_positions": calls}
            if len(calls) >= 4:
                q1, q3 = calls[len(calls) // 4], calls[(3 * len(calls)) // 4]
                rows[s_["sample"]]["interquartile_span_bp"] = q3 - q1
                rows[s_["sample"]]["median_position"] = calls[len(calls) // 2]
        out["results"][f"amp_{a}"] = rows
        print(f"\namp_{a} ({len(refs[a])} bp)", flush=True)
        for name, v in rows.items():
            span = v.get("interquartile_span_bp")
            print(f"   {name:20s} reads {v['reads']:>6} unique {v['unique_mutations']:>4} "
                  f"freq/1e3 {v['mutation_frequency_x1e3']:>8} C>T/G>A {str(v['pct_CT_GA']):>6}"
                  + (f"  middle half within {span} bp" if span else ""), flush=True)
    json.dump(out, open(OUT / "tam_reporter.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
