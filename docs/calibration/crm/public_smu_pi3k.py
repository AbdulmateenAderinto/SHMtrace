"""Sμ mutation libraries of GSE77788 (Compagno et al. 2017), re-analysed with the project workflow.

    python crm/public_smu_pi3k.py           -> crm/out/public_smu/

Primary B cells from wild-type and *Aicda*-null mice were treated with DMSO, idelalisib or duvelisib, three
replicates per group and treatment, and the Sμ region was amplified and sequenced. The published claim is that
PI3Kδ blockade raises AID-dependent mutation; this asks whether the control-calibrated workflow recovers a
treatment effect while the *Aicda*-null libraries stay at background.

Read handling, which decides whether the comparison is fair at all:

* The deposit mixes layouts, and the mix is not random: all nine *Aicda*-null libraries are paired, while six
  of the nine wild-type libraries are 300-bp single-end. Analysed as deposited, genotype would be confounded
  with library preparation.
* The paired mates do not overlap. R1 and the single reads cover one end of the 651-bp amplicon (reference
  positions 401-650 and 351-650), R2 covers the other (1-250), with a gap between. Merging them is therefore
  impossible: `io_utils.merge_pairs` keeps only 9-21% of pairs, and those few are the shortest molecules, a
  length-biased subset that would differ systematically between the layouts.
* So every library is analysed the same way: R1 alone for the paired libraries, the single file for the
  single-end ones, over the stretch all of them read. Nothing is merged and no library contributes its mate.

The window is reference positions 410-625: the common covered stretch, less the primer at the 3' end of the
reference and the last bases each read reaches. The reference is the bundled 651-bp Sμ fragment, which is
exactly this amplicon (its ends match where the reads start and stop).

Controls are condition-matched: each treatment is its own project, calibrated on the *Aicda*-null libraries of
that same treatment, so a drug effect on the background cannot be mistaken for a drug effect on mutation.

These libraries are sequenced far beyond one read per molecule, so read-weighted mutation frequency is the
estimand; unique-variant density is depth-dependent and is reported only alongside it.
"""
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
from shmtrace import project, report  # noqa: E402

DATA = Path.home() / "Nanoshm/public_data/GSE77788"
REF = REPO / "src/shmtrace/references/mouse_igh_smu_ch12.fasta"
S5F = str(Path.home() / "Nanoshm/models/HH_S5F.csv")
OUT = CRM / "out/public_smu"
WINDOW = [410, 625]
DEPTH = 20_000
# GSE77788: run -> (genotype, treatment, replicate, layout). Layout is from the ENA records and is confirmed
# by the files themselves: the single-end runs carry one 300-bp read, the paired runs two 250-bp mates.
RUNS = [("SRR3160391", "WT", "DMSO", 1, "single"), ("SRR3160392", "WT", "DMSO", 2, "single"),
        ("SRR3160393", "WT", "DMSO", 3, "paired"), ("SRR3160394", "WT", "Idelalisib", 1, "single"),
        ("SRR3160395", "WT", "Idelalisib", 2, "single"), ("SRR3160396", "WT", "Idelalisib", 3, "paired"),
        ("SRR3160397", "WT", "Duvelisib", 1, "single"), ("SRR3160398", "WT", "Duvelisib", 2, "single"),
        ("SRR3160399", "WT", "Duvelisib", 3, "paired"), ("SRR3160400", "AID_KO", "DMSO", 1, "paired"),
        ("SRR3160401", "AID_KO", "DMSO", 2, "paired"), ("SRR3160402", "AID_KO", "DMSO", 3, "paired"),
        ("SRR3160403", "AID_KO", "Idelalisib", 1, "paired"), ("SRR3160404", "AID_KO", "Idelalisib", 2, "paired"),
        ("SRR3160405", "AID_KO", "Idelalisib", 3, "paired"), ("SRR3160406", "AID_KO", "Duvelisib", 1, "paired"),
        ("SRR3160407", "AID_KO", "Duvelisib", 2, "paired"), ("SRR3160408", "AID_KO", "Duvelisib", 3, "paired")]


def read_path(run, layout):
    """R1 alone for a paired run; the single file otherwise. Never the pair: the mates do not overlap."""
    return str(DATA / (f"{run}_1.fastq.gz" if layout == "paired" else f"{run}.fastq.gz"))


def build(treatment):
    samples = [{"name": f"{g}_{t}_{r}", "path": read_path(run, lay), "group": g, "control": g == "AID_KO"}
               for run, g, t, r, lay in RUNS if t == treatment]
    return {"name": f"GSE77788 Smu, {treatment}", "region": "generic_amplicon", "platform": "illumina",
            "reference": str(REF), "depth": DEPTH, "targeting_model": S5F, "samples": samples,
            "overrides": {"analysis_window": WINDOW, "end_trim_bp": 0, "min_align_len": 150,
                          "max_total_len": 400, "exclude_indels": True, "jobs": 0}}


# The informative half of the amplicon. AID acts on the Sμ repeats, which lie in reference positions 1-327;
# only R2 reads reach them. Every R2 file belongs to a paired library, so a contrast built from R2 alone has no
# layout confound at all - but only one wild-type replicate per treatment is paired, so it cannot test the drug
# effect. It can test the genotype effect, pooling the three paired wild-type libraries against the nine
# *Aicda*-null ones, with treatment left as a nuisance factor.
R2_WINDOW = [20, 240]


def build_repeats():
    samples = []
    for run, g, t, r, lay in RUNS:
        if lay != "paired":
            continue
        samples.append({"name": f"{g}_{t}_{r}", "path": str(DATA / f"{run}_2.fastq.gz"), "group": g,
                        "control": g == "AID_KO"})
    return {"name": "GSE77788 Smu repeats, paired libraries only (R2)", "region": "generic_amplicon",
            "platform": "illumina", "reference": str(REF), "depth": DEPTH, "targeting_model": S5F,
            "samples": samples,
            "overrides": {"analysis_window": R2_WINDOW, "end_trim_bp": 0, "min_align_len": 150,
                          "max_total_len": 400, "exclude_indels": True, "jobs": 0}}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    layouts = {f"{g}_{t}_{r}": lay for run, g, t, r, lay in RUNS}
    runs = {f"{g}_{t}_{r}": run for run, g, t, r, lay in RUNS}
    allout = {"dataset": "GSE77788 / PRJNA311552 (Compagno et al. 2017, Nature 542:489)",
              "window_reference_positions": WINDOW, "depth_per_library": DEPTH,
              "read_handling": "R1 only for paired runs; single file otherwise; no merging (mates do not overlap)",
              "treatments": {}}
    for treatment in ("DMSO", "Idelalisib", "Duvelisib"):
        res = project.run_project(build(treatment))
        out = OUT / treatment
        out.mkdir(parents=True, exist_ok=True)
        report.write_all(res, str(out))
        json.dump(project.public(res), open(out / "project.json", "w"), default=str)
        rows = []
        for r in sorted(res["samples"], key=lambda r: (r["group"], r["sample"])):
            rows.append({"sample": r["sample"], "run": runs[r["sample"]], "layout": layouts[r["sample"]],
                         "group": r["group"], "reads": r["reads"], "unique_mutations": r["unique_mutations"],
                         "mutation_frequency_x1e3": r["mutation_frequency_x1e3"],
                         "unique_density_x1e5": r.get("unique_density_x1e5"),
                         "pct_CT_GA": r.get("pct_CT_GA"), "s5f_index": r.get("s5f_index"),
                         "analysed_bp": r.get("analysed_bp")})
            print(f"{treatment:11s} {r['sample']:18s} {layouts[r['sample']]:6s} reads {r['reads']:>6} "
                  f"unique {r['unique_mutations']:>4} freq/1e3 {r['mutation_frequency_x1e3']}", flush=True)
        allout["treatments"][treatment] = {"samples": rows,
                                           "group_tests": project.public(res).get("group_tests")}
    # the repeat-region genotype contrast, paired libraries only
    res = project.run_project(build_repeats())
    out = OUT / "repeats_paired"
    out.mkdir(parents=True, exist_ok=True)
    report.write_all(res, str(out))
    json.dump(project.public(res), open(out / "project.json", "w"), default=str)
    rows = []
    for r in sorted(res["samples"], key=lambda r: (r["group"], r["sample"])):
        rows.append({"sample": r["sample"], "group": r["group"], "reads": r["reads"],
                     "unique_mutations": r["unique_mutations"],
                     "mutation_frequency_x1e3": r["mutation_frequency_x1e3"],
                     "unique_density_x1e5": r.get("unique_density_x1e5"),
                     "pct_CT_GA": r.get("pct_CT_GA"), "s5f_index": r.get("s5f_index")})
        print(f"repeats(R2)  {r['sample']:22s} reads {r['reads']:>6} unique {r['unique_mutations']:>4} "
              f"freq/1e3 {r['mutation_frequency_x1e3']}", flush=True)
    allout["repeats_paired"] = {"window_reference_positions": R2_WINDOW, "samples": rows,
                                "group_tests": project.public(res).get("group_tests"),
                                "note": "genotype contrast only; treatment is a nuisance factor here"}
    json.dump(allout, open(OUT / "public_smu.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
