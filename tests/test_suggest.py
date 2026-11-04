"""Suggested verdicts for called mutations (review.py): the decisive rules, the score, and what is left out."""
from shmtrace import review


def _call(sample, pos, alt, reads, depth=500, p=1e-8, plus=None, minus=0, hp=1, hot=False):
    return {"sample": sample, "group": "", "pos": pos, "ref": "A", "alt": alt, "reads": reads, "depth": depth,
            "frequency": reads / depth, "p_value": p, "plus_reads": reads if plus is None else plus,
            "minus_reads": minus, "context": "AAAAA", "homopolymer": hp, "WRC_GYW": hot}


def _reads(sample, pos, alt, quals, start=1, end=200):
    return [{"name": f"{sample}{pos}{i}", "fate": "kept", "start": start, "end": end, "strand": 1,
             "diffs": [[pos, alt, q]], "indels": []} for i, q in enumerate(quals)]


def _res(variants, reads, background=None, samples=None):
    return {"samples": samples or [{"sample": "WT1"}, {"sample": "KO1", "control": True}],
            "variants": variants, "control_background_alt": background or {},
            "_call_reads": {s: {"reads": rs} for s, rs in reads.items()}, "params": {}}


def test_a_well_supported_call_is_suggested_for_acceptance():
    v = _call("WT1", 50, "G", 5)
    s = review.suggest(_res([v], {"WT1": _reads("WT1", 50, "G", [35, 38, 40, 41, 36])}, {"50G": 0.0002}))
    s = s["WT1|50|A>G"]
    assert s["suggestion"] == "accept" and s["score"] >= 2
    assert any("times the reads" in r["text"] for r in s["reasons"])


def test_decisive_reasons_reject_whatever_the_score():
    ctl = _call("KO1", 60, "T", 6)
    bg = _call("WT1", 70, "T", 2)
    shared = _call("WT1", 60, "T", 8)
    res = _res([ctl, bg, shared], {"WT1": _reads("WT1", 60, "T", [40] * 8) + _reads("WT1", 70, "T", [40, 40])},
               {"70T": 0.01})
    s = review.suggest(res)
    assert s["KO1|60|A>T"]["suggestion"] == "reject"                 # a call in a control
    assert s["WT1|70|A>T"]["suggestion"] == "reject"                 # 0.01 x 500 = 5 expected, 2 seen
    assert s["WT1|60|A>T"]["suggestion"] == "reject"                 # the control has it too
    assert all(any(r["decisive"] for r in s[k]["reasons"]) for k in s)


def test_low_quality_and_read_ends_count_against_and_one_read_calls_wait_for_an_eye():
    low = _call("WT1", 80, "C", 3, p=0.01)
    ends = _call("WT1", 90, "C", 1, p=0.01)
    reads = _reads("WT1", 80, "C", [8, 12, 15]) + _reads("WT1", 90, "C", [38], start=85, end=200)
    s = review.suggest(_res([low, ends], {"WT1": reads}))
    assert s["WT1|80|A>C"]["suggestion"] == "reject"
    assert any("end within" in r["text"] for r in s["WT1|90|A>C"]["reasons"])
    assert s["WT1|90|A>C"]["suggestion"] == "check"


def test_the_hotspot_and_other_samples_do_not_steer_the_suggestion():
    a = _call("WT1", 50, "G", 4, hot=False)
    b = _call("WT1", 50, "G", 4, hot=True)
    reads = {"WT1": _reads("WT1", 50, "G", [40] * 4)}
    other = _call("WT2", 50, "G", 4)
    samples = [{"sample": "WT1"}, {"sample": "WT2"}, {"sample": "KO1", "control": True}]
    plain = review.suggest(_res([a], reads, samples=samples))["WT1|50|A>G"]
    hot = review.suggest(_res([b], reads, samples=samples))["WT1|50|A>G"]
    shared = review.suggest(_res([a, other], reads, samples=samples))["WT1|50|A>G"]
    assert plain == hot == shared


def test_one_strand_counts_only_when_the_library_has_both():
    one = _call("WT1", 50, "G", 4, plus=4, minus=0)
    both = _call("WT1", 60, "G", 4, plus=2, minus=2)
    reads = {"WT1": _reads("WT1", 50, "G", [40] * 4) + _reads("WT1", 60, "G", [40] * 4)}
    s = review.suggest(_res([one, both], reads))
    assert any("one strand" in r["text"] for r in s["WT1|50|A>G"]["reasons"])
    alone = review.suggest(_res([one], {"WT1": reads["WT1"][:4]}))
    assert not any("one strand" in r["text"] for r in alone["WT1|50|A>G"]["reasons"])
