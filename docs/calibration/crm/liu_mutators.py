"""Different cytidine deaminases on one GFP substrate (Liu et al. 2018, SRP159221).

    python crm/liu_mutators.py        -> crm/out/liu/

The deposit targets the same 238-bp EGFP amplicon with a panel of editors delivered the same way, three
replicates each, alongside a no-editor control. That makes it the contrasting-mutator test: whether the
substitution spectrum SHMtrace reports separates AID from an APOBEC on identical sequence and identical
callable positions, rather than only separating "editor" from "no editor".

The contrast is fixed before looking at any result: dCas9 with DBE-AIDmono against dCas9 with DBE-humanA3a,
both in wild-type 293T and both in the UNG/MSH2-null background, with the NC arm as the control. The endpoint
is the share of calls that are C->T or G->A and the share at WRC/GYW, computed over the positions all arms
cover. APOBEC arms are labelled as such and never pooled with AID.

Two cautions the data impose. The archive holds 193 run records, not 193 experiments: a group of three is
three replicates of one condition, and single-guide arms (one run each) are not replicated and are left out.
And the motif share is descriptive only - it is reported for whatever was called, and is never used to decide
what to call, so it cannot both filter the calls and then testify that they look like AID.
"""
import collections
import csv
import gzip
import json
import re
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
from shmtrace import io_utils, project  # noqa: E402

DATA = Path.home() / "Nanoshm/public_data/SRP159221"
TSV = Path("/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab"
           "/nanoshm/output/reviews/SHMtrace_dataset_search_2026-10-17/SRP159221.tsv")
OUT = CRM / "out/liu"
# the prespecified arms: label -> (group, mutator family). Each has three replicate runs.
ARMS = {
    "293T_puroGFP_NC": ("NC", "none"),
    # the archive's labels are inconsistent about separators, so they are written out exactly as deposited
    "293TpuroGFP_dCas9_DBE_AID_mono": ("AID_wt293T", "AID"),
    "293TpuroGFP_dCas9_DBE-humanA3a": ("A3A_wt293T", "APOBEC"),
    "293T_GFP_UNG-/-MSH2-/-_dCas9_DBE-AIDmono": ("AID_ungmsh2", "AID"),
    "293T_GFP_UNG-/-MSH2-/-_dCas9_DBE-humanA3a": ("A3A_ungmsh2", "APOBEC"),
}
DEPTH = 20_000


def label(title):
    t = re.sub(r"^.*?amplicon seq of ", "", title)
    t = re.sub(r"\s*(Homo sapiens|Mus musculus).*$", "", t)
    return re.sub(r"[-_ ]?(R|replicate)\s*\d*$", "", t).strip()


def runs_by_arm():
    out = collections.defaultdict(list)
    for x in csv.DictReader(open(TSV), delimiter="\t"):
        out[label(x["experiment_title"])].append(x["run_accession"])
    return out


def build_reference(nc_runs, out_path, min_frac=0.9):
    """Consensus of the merged no-editor reads: the construct as the experiment actually carries it."""
    cols = collections.defaultdict(collections.Counter)
    n = 0
    for run in nc_runs:
        for _nm, seq, _q in io_utils.merge_pairs(str(DATA / f"{run}_1.fastq.gz"),
                                                 str(DATA / f"{run}_2.fastq.gz")):
            if len(seq) != 238:
                continue
            for i, b in enumerate(seq):
                cols[i][b] += 1
            n += 1
            if n >= 30_000:
                break
        if n >= 30_000:
            break
    ref = "".join(cols[i].most_common(1)[0][0] for i in range(len(cols)))
    weakest = min((cols[i].most_common(1)[0][1] / sum(cols[i].values())) for i in range(len(cols)))
    out_path.write_text(">EGFP_amplicon consensus of the no-editor arm\n"
                        + "\n".join(ref[i:i + 70] for i in range(0, len(ref), 70)) + "\n")
    return ref, n, weakest


def main():
    if not TSV.exists():
        sys.exit(f"missing run table {TSV}")
    OUT.mkdir(parents=True, exist_ok=True)
    by_arm = runs_by_arm()
    missing = [a for a in ARMS if a not in by_arm]
    if missing:
        sys.exit(f"arms not in the run table: {missing}")
    nc_runs = by_arm["293T_puroGFP_NC"]
    ref, n_used, weakest = build_reference(nc_runs, OUT / "egfp_amplicon.fasta")
    print(f"reference: {len(ref)} bp from {n_used} merged no-editor reads; "
          f"weakest consensus position {100 * weakest:.1f}%")
    samples, groups = [], {}
    for arm, (grp, family) in ARMS.items():
        for i, run in enumerate(sorted(by_arm[arm]), 1):
            name = f"{grp}_{i}"
            r1, r2 = DATA / f"{run}_1.fastq.gz", DATA / f"{run}_2.fastq.gz"
            if not r1.exists() or not r2.exists():
                print(f"   skipping {name}: reads missing")
                continue
            samples.append({"name": name, "path": f"{r1}+{r2}", "group": grp, "control": grp == "NC"})
            groups[name] = {"arm": arm, "run": run, "family": family}
    cfg = {"name": "Liu 2018 GFP editors", "region": "generic_amplicon", "platform": "illumina",
           "reference": str(OUT / "egfp_amplicon.fasta"), "depth": DEPTH, "samples": samples,
           "overrides": {"analysis_window": [10, len(ref) - 10], "end_trim_bp": 0, "min_align_len": 150,
                         "max_total_len": 400, "exclude_indels": True, "jobs": 1}}
    pub = project.public(project.run_project(cfg))
    rows = {}
    for s in pub["samples"]:
        rows[s["sample"]] = {**groups[s["sample"]], "reads": s["reads"],
                             "unique_mutations": s["unique_mutations"],
                             "mutation_frequency_x1e3": s["mutation_frequency_x1e3"],
                             "pct_CT_GA": s.get("pct_CT_GA"),
                             "pct_GC_in_WRC_GYW": s.get("pct_GC_in_WRC_GYW"),
                             "pct_GC": s.get("pct_GC"),
                             "removed_by_frequency_rule": s.get("high_fraction_sites_removed")}
    json.dump({"dataset": "SRP159221 (Liu et al. 2018, Cell Rep)",
               "reference": "238-bp EGFP amplicon, consensus of the no-editor arm",
               "prespecified_contrast": "DBE-AIDmono against DBE-humanA3a, in wild-type and UNG/MSH2-null 293T",
               "depth": DEPTH, "samples": rows}, open(OUT / "liu_mutators.json", "w"), indent=1, default=str)
    print(f"{'sample':16s} {'family':8s} {'reads':>6} {'unique':>7} {'freq/1e3':>9} {'C>T/G>A%':>9} {'WRC/GYW%':>9}")
    for name, v in sorted(rows.items()):
        print(f"{name:16s} {v['family']:8s} {v['reads']:>6} {v['unique_mutations']:>7} "
              f"{v['mutation_frequency_x1e3']:>9} {str(v['pct_CT_GA']):>9} {str(v['pct_GC_in_WRC_GYW']):>9}")


if __name__ == "__main__":
    main()
