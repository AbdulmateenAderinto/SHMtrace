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


def raw_rates(refs, min_q=30, limit=150_000):
    """Mismatches per 1,000 bases against each amplicon's own reference, with no model in the way.

    This exists because the calibrated path cannot be used on this arm: it has one no-editor library, and the
    error model is fit from the controls, so with a single control it has nothing to calibrate against. A
    plain rate is interpretable without that machinery, and is reported as what it is - a descriptive rate
    that still contains sequencing error, not a mutation call.
    """
    out = {}
    for run, name, _w, _c in REPORTER:
        mism, bases = collections.Counter(), collections.Counter()
        n = 0
        for _nm, seq, qual in io_utils.read_fastq(str(DATA / f"{run}.fastq.gz")):
            n += 1
            if n > limit:
                break
            k = seq[:20]
            ref = refs.get(k)
            if ref is None:
                continue
            for i, (b, q) in enumerate(zip(seq, qual, strict=False)):
                if i >= len(ref) or ord(q) - 33 < min_q:
                    continue
                bases[k] += 1
                if b != ref[i]:
                    mism[k] += 1
        out[name] = {f"amp_{k}": round(1000 * mism[k] / bases[k], 3)
                     for k in refs if bases[k] > 10_000}
    return out


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
           "amplicons": {f"amp_{a_}": {"bp": len(s_), "weakest_consensus_pct": round(100 * weakest[a_], 1)}
                         for a_, s_ in refs.items()},
           "depth": DEPTH}
    raw = raw_rates(refs)
    print("\nmismatches per 1,000 bases at Q>=30, against each amplicon's own reference:", flush=True)
    amps = sorted({a for v in raw.values() for a in v})
    print("   " + "library".ljust(20) + "".join(a[:14].ljust(16) for a in amps), flush=True)
    for name in (n for _r, n, _w, _c in REPORTER):
        print("   " + name.ljust(20) + "".join(str(raw[name].get(a, "-")).ljust(16) for a in amps), flush=True)
    # the calibrated path, for contrast, on the amplicon with the clearest raw signal
    target = max(amps, key=lambda a: raw["dCas9_AIDX_sgRNA"].get(a, 0) - raw["mock"].get(a, 0))
    key = target.replace("amp_", "")
    # the single reference and its competitors, written for this run
    (OUT / f"ref_{key}.fasta").write_text(
        f">amp_{key}\n" + "\n".join(refs[key][i:i + 70] for i in range(0, len(refs[key]), 70)) + "\n")
    (OUT / f"others_{key}.fasta").write_text(
        "".join(f">amp_{b_}\n{refs[b_]}\n" for b_ in refs if b_ != key))
    samples = [{"name": name, "path": str(DATA / f"{run}.fastq.gz"), "group": name, "control": ctrl}
               for run, name, _w, ctrl in REPORTER]
    cfg = {"name": f"TAM {target}", "region": "generic_amplicon", "platform": "illumina",
           "reference": str(OUT / f"ref_{key}.fasta"), "contaminants": [str(OUT / f"others_{key}.fasta")],
           "depth": DEPTH, "samples": samples,
           "overrides": {"analysis_window": [10, len(refs[key]) - 10], "end_trim_bp": 0,
                         "min_align_len": 200, "max_total_len": 400, "exclude_indels": True, "jobs": 1}}
    pub = project.public(project.run_project(cfg))
    calibrated = {s_["sample"]: {"condition": what[s_["sample"]], "reads": s_["reads"],
                                 "unique_mutations": s_["unique_mutations"],
                                 "mutation_frequency_x1e3": s_["mutation_frequency_x1e3"]}
                  for s_ in pub["samples"]}
    out["raw_mismatch_per_1e3"] = raw
    out["calibrated_run"] = {"amplicon": target, "controls": 1,
                             "note": "one no-editor library, so the error model has a single control; "
                                     "reported for contrast with the raw rates, not as a calibrated result",
                             "samples": calibrated}
    print(f"\ncalibrated path on {target} (one control):", flush=True)
    for name, v in calibrated.items():
        print(f"   {name:20s} unique {v['unique_mutations']:>4}  freq/1e3 {v['mutation_frequency_x1e3']}",
              flush=True)
    json.dump(out, open(OUT / "tam_reporter.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
