"""What the 50% germline rule does to an acquired variant that selection has driven high.

    python crm/germline_policy_demo.py        -> crm/out/germline_policy/demo.json

The engineered-diversification datasets (an AID-diversified TCR, base-editor antibody maturation) differ from
a physiological SHM experiment in one way that matters here: the reference is the *known ancestor*, the
construct that was put into the cells, and a mutation acquired during diversification can be driven to high
frequency, or to fixation, by selection. SHMtrace treats a change carried by at least
`germline_freq_threshold` (0.50) of a library's reads as germline and removes it.

This builds a controlled library in which four variants sit at 20%, 50%, 90% and 100% of reads, with a
matched control library carrying none, and reports which survive. It then repeats the run with the threshold
raised above 1, the obvious way to "turn the rule off", and shows what that does to the error model - which is
why raising the threshold is not a safe interface. Nothing here is downloaded; the reads are synthetic.
"""
import json
import random
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
from shmtrace import project  # noqa: E402

OUT = CRM / "out/germline_policy"
SEED = 11
LEN, DEPTH = 600, 400
# acquired variants: (1-based position, share of reads carrying it). The alt is chosen per position as the
# next base after the reference base, so it can never coincide with the reference - a mistake that silently
# plants nothing and looks like a caller failure.
ACQUIRED = [(120, 0.20), (240, 0.50), (360, 0.90), (480, 1.00)]
# a second series around fixation, to separate the frequency rule from what happens at 100%
NEAR_FIXED = [(120, 0.95), (240, 0.98), (360, 0.995), (480, 1.00)]


def alt_for(ref, pos):
    """A base that differs from the reference at this 1-based position."""
    return "ACGT"[("ACGT".index(ref[pos - 1]) + 1) % 4]


def make_reference(rng):
    return "".join(rng.choice("ACGT") for _ in range(LEN))


def write_fastq(path, ref, variants, n, rng, err=0.001):
    """n reads of `ref`, with each (pos, frac) carried by that fraction, plus a low error rate."""
    carriers = {}
    for pos, frac in variants:
        k = round(frac * n)
        carriers[pos] = (alt_for(ref, pos), set(rng.sample(range(n), k)))
    with open(path, "w") as fh:
        for i in range(n):
            seq = list(ref)
            for pos, (alt, who) in carriers.items():
                if i in who:
                    seq[pos - 1] = alt
            for j in range(LEN):                       # sequencing error, uniform and rare
                if rng.random() < err:
                    seq[j] = rng.choice([b for b in "ACGT" if b != seq[j]])
            fh.write(f"@r{i}\n{''.join(seq)}\n+\n{'I' * LEN}\n")


def build(tmp, threshold):
    """threshold None = defaults; a float sets germline_freq_threshold; "policy" = ancestral_reference."""
    samples = [{"name": "case", "path": str(tmp / "case.fastq"), "group": "case", "control": False}]
    for c in ("ctrl1", "ctrl2", "ctrl3"):
        samples.append({"name": c, "path": str(tmp / f"{c}.fastq"), "group": "control", "control": True})
    ov = {"analysis_window": [30, LEN - 30], "end_trim_bp": 0, "min_align_len": 300,
          "max_total_len": 900, "exclude_indels": True, "jobs": 1}
    if threshold == "policy":
        ov["ancestral_reference"] = True
    elif threshold is not None:
        ov["germline_freq_threshold"] = threshold
    return {"name": f"germline demo (threshold={threshold})", "region": "generic_amplicon",
            "platform": "ont", "reference": str(tmp / "ref.fasta"), "depth": DEPTH,
            "samples": samples, "overrides": ov}


def run(tmp, threshold, acquired):
    res = project.run_project(build(tmp, threshold))
    pub = project.public(res)
    case = next(s for s in pub["samples"] if s["sample"] == "case")
    called = {v["pos"]: v["alt"] for v in pub["variants"] if v["sample"] == "case"}
    ctrl_calls = {s["sample"]: s["unique_mutations"] for s in pub["samples"] if s["sample"] != "case"}
    return {"threshold": threshold, "case_unique_mutations": case["unique_mutations"],
            "acquired_recovered": {f"pos {pos} at {100 * frac:g}% of reads": (pos in called)
                                   for pos, frac in acquired},
            "germline_removed": case.get("high_fraction_sites"),
            "control_unique_mutations": ctrl_calls,
            "case_total_called_positions": len(called)}


def series(name, acquired, thresholds, ref, rng):
    tmp = OUT / f"reads_{name}"
    tmp.mkdir(exist_ok=True)
    (tmp / "ref.fasta").write_text(">demo_ancestor\n"
                                   + "\n".join(ref[i:i + 70] for i in range(0, LEN, 70)) + "\n")
    write_fastq(tmp / "case.fastq", ref, acquired, DEPTH, rng)
    for c in ("ctrl1", "ctrl2", "ctrl3"):
        write_fastq(tmp / f"{c}.fastq", ref, [], DEPTH, rng)
    runs = []
    for th in thresholds:
        runs.append(run(tmp, th, acquired))
    return runs


def main():
    rng = random.Random(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    ref = make_reference(random.Random(SEED))
    out = {"spread": {"acquired": [{"pos": p, "alt": alt_for(ref, p), "share": f} for p, f in ACQUIRED],
                      "runs": series("spread", ACQUIRED, [None, 1.01, "policy"], ref, rng)},
           "near_fixed": {"acquired": [{"pos": p, "alt": alt_for(ref, p), "share": f} for p, f in NEAR_FIXED],
                          "runs": series("fixed", NEAR_FIXED, [None, "policy"], ref, rng)}}
    json.dump(out, open(OUT / "demo.json", "w"), indent=1)
    for key in ("spread", "near_fixed"):
        for r in out[key]["runs"]:
            lab = ("default (0.50)" if r["threshold"] is None else
                   "ancestral_reference policy" if r["threshold"] == "policy" else
                   f"threshold {r['threshold']}")
            print(f"\n[{key}] {lab}:")
            for k, v in r["acquired_recovered"].items():
                print(f"   acquired at {k:30s} {'recovered' if v else 'NOT recovered'}")
            print(f"   removed by the germline rule: {r['germline_removed']}")


if __name__ == "__main__":
    main()
