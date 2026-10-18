"""AID-diversified MAGE-A1 TCR (Bassan et al. 2021), called against its own ancestral construct.

    python crm/tcr_bassan.py        -> crm/out/tcr/tcr_bassan.json

Run `crm/tcr_prepare.py` first: it reconstructs the 1,785-bp `hT27` reference exactly from the deposit's
CIGAR and MD tags and writes the reads as FASTQ beside the deposit.

What this is, and is not. The construct that went into the cells is the reference, so a variant that
selection drove high is the result rather than germline: the run uses `ancestral_reference`. Only the h5
series has a sequenced parental sample (`1_h5_wt`); h7, h8 and h12 have none, so `1_h5_wt` is the *available
baseline* for the experiment and not a matched AID-off control for each founder. There is one baseline
library, so this cannot establish a false-discovery rate; it is a recovery and linkage demonstration, and the
baseline's own call level is reported beside every sample as the honest comparison.

`11_h12_HA` holds 19 reads and is excluded, with its count reported rather than hidden.

Phasing is off. SHMtrace's phasing exists to drop minor contaminating haplotypes, which is the wrong
operation for a population that is diversified on purpose; linkage is instead counted here from the per-read
calls, over reads that cover every site in the combination, so a combination is never inferred from reads
that could not have shown it.

Selected endpoint frequencies are not intrinsic mutation rates, and founder descendants are correlated, so
shared combinations are not evidence of an independent evolutionary history.
"""
import collections
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
from shmtrace import project  # noqa: E402

DATA = Path.home() / "Nanoshm/public_data/TCR_figshare_12803501"
OUT = CRM / "out/tcr"
WINDOW = [30, 1755]
DEPTH = 6000
MIN_READS = 1000          # a library below this cannot support a population estimate
BASELINE = "h5_wt"
SAMPLES = [("1_h5_wt", "h5", "wt"), ("2_h5_MHA", "h5", "MHA"), ("3_h5_HA", "h5", "HA"),
           ("4_h5_HAHA", "h5", "HAHA"), ("5_h7_HA", "h7", "HA"), ("6_h7_HAHA", "h7", "HAHA"),
           ("7_h8_MHA", "h8", "MHA"), ("8_h8_HA", "h8", "HA"), ("9_h8_HAHE", "h8", "HAHE"),
           ("10_h12_MHA", "h12", "MHA"), ("11_h12_HA", "h12", "HA"), ("12_h12_HAHE", "h12", "HAHE")]


def usable():
    keep, dropped = [], {}
    for stem, founder, cond in SAMPLES:
        fq = DATA / "fastq" / f"{stem}.fastq"
        n = sum(1 for _ in open(fq)) // 4 if fq.exists() else 0
        if n < MIN_READS:
            dropped[stem] = n
        else:
            keep.append((stem, founder, cond, n))
    return keep, dropped


def build(keep):
    samples = []
    for stem, founder, cond, _n in keep:
        name = f"{founder}_{cond}"
        samples.append({"name": name, "path": str(DATA / "fastq" / f"{stem}.fastq"),
                        "group": founder, "control": name == BASELINE})
    return {"name": "Bassan 2021 MAGE-A1 TCR", "region": "generic_amplicon", "platform": "pacbio",
            "reference": str(OUT / "hT27.fasta"), "depth": DEPTH, "samples": samples,
            "overrides": {"analysis_window": WINDOW, "end_trim_bp": 0, "min_align_len": 1200,
                          "max_total_len": 2500, "exclude_indels": True, "ancestral_reference": True,
                          "phase_mode": "off", "jobs": 1}}


def linkage(details, called, top=6):
    """Combinations of called variants carried together, counted only over reads covering every site."""
    per_read = details.get("per_read_mutations") or []
    sites = sorted(called)
    if not sites or not per_read:
        return []
    combos = collections.Counter()
    covered = 0
    for read in per_read:
        muts = {(m["pos"], m["alt"]) for m in read.get("mutations", [])} if isinstance(read, dict) else set()
        span = read.get("covered") if isinstance(read, dict) else None
        if span is None:
            continue
        lo, hi = span
        here = [s for s in sites if lo <= s[0] <= hi]
        if len(here) < len(sites):
            continue
        covered += 1
        combos[tuple(sorted(muts & set(sites)))] += 1
    return [{"combination": [f"{p}{a}" for p, a in c], "reads": n,
             "of_reads_covering_all_sites": covered}
            for c, n in combos.most_common(top) if c]


def policy_sensitivity(keep, depth=2000):
    """What the ancestral_reference policy changes on this deposit, rather than in principle.

    The prompt asks for a bounded sensitivity analysis, and the answer matters: if every selected variant
    sits below the 50% rule, the policy is a precaution that changes nothing here, and saying otherwise
    would overstate it.
    """
    out = {}
    for policy in (True, False):
        cfg = build(keep)
        cfg["depth"] = depth
        cfg["overrides"]["ancestral_reference"] = policy
        pub = project.public(project.run_project(cfg))
        out["on" if policy else "off"] = {
            s["sample"]: {"unique_mutations": s["unique_mutations"],
                          "removed_by_frequency_rule": s.get("high_fraction_sites_removed")}
            for s in pub["samples"]}
    return {"depth": depth, "runs": out}


def main():
    if not (OUT / "hT27.fasta").exists():
        sys.exit("run crm/tcr_prepare.py first")
    keep, dropped = usable()
    res = project.run_project(build(keep))
    pub = project.public(res)
    rows = {}
    for s in pub["samples"]:
        called = {(v["pos"], v["alt"]) for v in pub["variants"] if v["sample"] == s["sample"]}
        rows[s["sample"]] = {"reads": s["reads"], "unique_mutations": s["unique_mutations"],
                             "mutation_frequency_x1e3": s["mutation_frequency_x1e3"],
                             "pct_CT_GA": s.get("pct_CT_GA"), "analysed_bp": s.get("analysed_bp"),
                             "called_positions": sorted(p for p, _a in called)}
    base = rows.get(BASELINE, {})
    out = {"dataset": "figshare 10.6084/m9.figshare.12803501 (Bassan et al. 2021, Eur J Immunol)",
           "reference": "hT27, 1785 bp, reconstructed from the deposit (see reference_provenance.json)",
           "policy": "ancestral_reference; phasing off; one baseline library, so no FDR is claimed",
           "window": WINDOW, "depth": DEPTH, "baseline_sample": BASELINE,
           "excluded_low_yield": dropped, "samples": rows,
           "baseline_unique_mutations": base.get("unique_mutations")}
    sens = policy_sensitivity([k for k in keep if k[0] in ("1_h5_wt", "6_h7_HAHA", "7_h8_MHA")])
    out["policy_sensitivity"] = sens
    json.dump(out, open(OUT / "tcr_bassan.json", "w"), indent=1, default=str)
    on, off = sens["runs"]["on"], sens["runs"]["off"]
    diff = {k: (on[k]["unique_mutations"], off[k]["unique_mutations"]) for k in on
            if on[k]["unique_mutations"] != off[k]["unique_mutations"]}
    print(f"\npolicy sensitivity at depth {sens['depth']}: "
          + ("no sample changed" if not diff else
             ", ".join(f"{k} {a} with the policy vs {b} without" for k, (a, b) in diff.items())))
    print(f"excluded for low yield: {dropped}")
    print(f"{'sample':12s} {'reads':>6} {'unique':>7} {'freq/1e3':>9}  C>T/G>A%")
    for name, v in sorted(rows.items()):
        mark = "  <- baseline" if name == BASELINE else ""
        print(f"{name:12s} {v['reads']:>6} {v['unique_mutations']:>7} "
              f"{v['mutation_frequency_x1e3']:>9}  {v['pct_CT_GA']}{mark}")


if __name__ == "__main__":
    main()
