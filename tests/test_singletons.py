"""Single-read mutations estimated against a control background (singletons.py)."""
import random

from nanoshm import project, singletons
from nanoshm.align import AlignedRead

REF = "ACGTTGCAAGCTAGCTAGGATCCATGCATGCAAATTTGGGCCCATATCGCGATCGATCGTAGCTAGCATGCA" * 3


def _reads(n, ref, rng, err=0.0, err_q=10, mutations=(), mut_q=40, base_q=40):
    """n reads over the whole reference; random errors at err_q, and each
    (pos, alt) in `mutations` written into one read at mut_q."""
    out = []
    plant = {i: m for i, m in enumerate(mutations)}
    for i in range(n):
        pairs = []
        for pos, b in enumerate(ref):
            base, q = b, base_q
            if rng.random() < err:
                base, q = rng.choice([x for x in "ACGT" if x != b]), err_q
            pairs.append((pos, base, q))
        if i in plant:
            pos, alt = plant[i]
            pairs[pos] = (pos, alt, mut_q)
        out.append(AlignedRead(f"r{i}", "ref", 0, len(ref), pairs, "", 1 if i % 2 else -1))
    return out


def _mutations(ref, k, rng):
    sites = rng.sample(range(10, len(ref) - 10), k)
    return [(s, rng.choice([x for x in "ACGT" if x != ref[s]])) for s in sites]


def test_controls_against_each_other_estimate_about_zero():
    rng = random.Random(1)
    lo, hi = 0, len(REF) - 1
    libs = [singletons.tally(_reads(300, REF, rng, err=0.01), REF, lo, hi) for _ in range(3)]
    est = singletons.estimate(libs[0], libs[1:], 35, 0.5, REF)
    # low-quality errors never reach the Q35 cut, so the background is empty
    assert est["observed"] == 0 and est["estimate"] == 0


def test_planted_single_read_mutations_are_recovered():
    rng = random.Random(2)
    lo, hi = 0, len(REF) - 1
    ctrls = [singletons.tally(_reads(300, REF, rng, err=0.01), REF, lo, hi) for _ in range(2)]
    muts = _mutations(REF, 40, rng)
    sample = singletons.tally(_reads(300, REF, rng, err=0.01, mutations=muts), REF, lo, hi)
    est = singletons.estimate(sample, ctrls, 35, 0.5, REF)
    assert est["observed"] == 40
    assert est["ci_low"] <= 40 <= est["ci_high"]


def test_low_quality_reads_of_the_same_alt_do_not_hide_a_mutation():
    # one high-quality read carries the mutation; two low-quality reads show
    # the same base by error. It is still one high-quality event.
    rng = random.Random(3)
    reads = _reads(50, REF, rng)
    pos, alt = 30, "A" if REF[30] != "A" else "C"
    for i, q in ((0, 45), (1, 8), (2, 9)):
        pairs = list(reads[i].pairs)
        pairs[pos] = (pos, alt, q)
        reads[i] = AlignedRead(reads[i].qname, "ref", 0, len(REF), pairs, "", 1)
    t = singletons.tally(reads, REF, 0, len(REF) - 1)
    assert [(e[0], e[1]) for e in singletons.count(t, 35, set())["events"]] == [(pos, alt)]


def test_called_events_are_not_counted_twice():
    rng = random.Random(4)
    muts = _mutations(REF, 5, rng)
    t = singletons.tally(_reads(100, REF, rng, mutations=muts), REF, 0, len(REF) - 1)
    assert len(singletons.count(t, 35, set(), called={muts[0]})["events"]) == 4


def test_efficiency_corrects_for_bases_under_the_cut():
    # every read base, and so every real mutation, is Q30: a Q35 cut would see
    # none of them, so the cut steps down until it keeps the bases
    rng = random.Random(5)
    lo, hi = 0, len(REF) - 1
    ctrls = [singletons.tally(_reads(200, REF, rng, base_q=30), REF, lo, hi) for _ in range(2)]
    muts = _mutations(REF, 20, rng)
    sample = singletons.tally(_reads(200, REF, rng, mutations=muts, mut_q=30, base_q=30),
                              REF, lo, hi)
    est = singletons.estimate(sample, ctrls, 35, 0.5, REF)
    assert est["q_min"] == 30 and est["q_min_requested"] == 35
    assert est["estimate"] == 20


def test_project_reports_combined_density(wt_fastq, ko_fastq, reference, contaminant):
    proj = {"name": "t", "region": "jh4_intron", "reference": str(reference),
            "contaminants": [str(contaminant)],
            "samples": [{"name": "WT", "path": str(wt_fastq), "group": "WT"},
                        {"name": "KO", "path": str(ko_fastq), "group": "KO", "control": True}]}
    res = project.run_project(proj)
    wt = next(r for r in res["samples"] if r["sample"] == "WT")
    assert wt["combined_density_x1e5"] >= wt["unique_density_x1e5"]
    assert wt["combined_density_ci_low"] <= wt["combined_density_x1e5"] \
        <= wt["combined_density_ci_high"]
    assert "WT" in res["single_read"]
    # one control has no peers to be measured against
    assert "combined_density_x1e5" not in next(r for r in res["samples"] if r["sample"] == "KO")
