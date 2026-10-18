"""The opt-in ancestral-reference policy (config.ancestral_reference).

When the reference is the construct that was put into the cells, a mutation acquired during diversification
can be driven high, or to fixation, by selection. The default rule treats a change carried by at least
`germline_freq_threshold` of a library's reads as germline and removes it, which is right for a physiological
SHM assay and wrong here. These tests pin both behaviours: the default must keep removing such changes, and
the policy must keep them while leaving everything else alone.
"""
import random

import pytest

pytest.importorskip("mappy")

from shmtrace import project  # noqa: E402

LEN, DEPTH = 500, 300
SEED = 7


def _ref(rng):
    return "".join(rng.choice("ACGT") for _ in range(LEN))


def _alt(ref, pos):
    """A base that differs from the reference at this 1-based position."""
    return "ACGT"["ACGT".index(ref[pos - 1]) + 1 & 3]


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _fastq(path, ref, variants, n, rng, err=0.0005, revcomp=False, cover=None):
    """n reads of `ref`; each (pos, frac) is carried by that fraction of them."""
    carriers = {pos: (_alt(ref, pos), set(rng.sample(range(n), round(frac * n))))
                for pos, frac in variants}
    lo, hi = cover or (0, LEN)
    with open(path, "w") as fh:
        for i in range(n):
            seq = list(ref)
            for pos, (alt, who) in carriers.items():
                if i in who:
                    seq[pos - 1] = alt
            for j in range(LEN):
                if rng.random() < err:
                    seq[j] = rng.choice([b for b in "ACGT" if b != seq[j]])
            s = "".join(seq)[lo:hi]
            if revcomp:
                s = _revcomp(s)
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


def _project(tmp, case_variants, ctrl_variants=(), ancestral=False, revcomp=False, cover=None):
    rng = random.Random(SEED)
    ref = _ref(random.Random(SEED))
    (tmp / "ref.fasta").write_text(">ancestor\n" + ref + "\n")
    _fastq(tmp / "case.fastq", ref, case_variants, DEPTH, rng, revcomp=revcomp, cover=cover)
    for c in ("c1", "c2", "c3"):
        _fastq(tmp / f"{c}.fastq", ref, ctrl_variants, DEPTH, rng)
    samples = [{"name": "case", "path": str(tmp / "case.fastq"), "group": "case", "control": False}]
    samples += [{"name": c, "path": str(tmp / f"{c}.fastq"), "group": "ctrl", "control": True}
                for c in ("c1", "c2", "c3")]
    ov = {"analysis_window": [40, LEN - 40], "end_trim_bp": 0, "min_align_len": 200,
          "max_total_len": 800, "exclude_indels": True, "jobs": 1}
    if ancestral:
        ov["ancestral_reference"] = True
    cfg = {"name": "t", "region": "generic_amplicon", "platform": "ont", "reference": str(tmp / "ref.fasta"),
           "depth": DEPTH, "samples": samples, "overrides": ov}
    pub = project.public(project.run_project(cfg))
    called = {v["pos"] for v in pub["variants"] if v["sample"] == "case"}
    case = next(s for s in pub["samples"] if s["sample"] == "case")
    return called, case, pub


def test_default_removes_acquired_variants_at_or_above_the_threshold(tmp_path):
    """The physiological default is unchanged: >= 50% is germline, and the removal is reported."""
    called, case, _ = _project(tmp_path, [(120, 0.20), (240, 0.50), (360, 0.90), (420, 1.00)])
    assert 120 in called                       # below the threshold, kept
    for pos in (240, 360, 420):
        assert pos not in called
    removed = {s["pos"] for s in case["high_fraction_sites"]}
    assert {240, 360, 420} <= removed          # and every removal is visible, not silent


def test_ancestral_reference_keeps_them_including_fixation(tmp_path):
    called, case, _ = _project(tmp_path, [(120, 0.20), (240, 0.50), (360, 0.90), (420, 1.00)],
                               ancestral=True)
    assert {120, 240, 360, 420} <= called
    assert not case["high_fraction_sites"]


def test_baseline_mismatch_is_not_called_either_way(tmp_path):
    """A fixed difference from the deposited reference that the controls share is not a mutation.

    The control-calibrated test handles this, not the germline rule, so the policy must not change it.
    """
    for ancestral in (False, True):
        called, _case, _ = _project(tmp_path, [(300, 1.00), (120, 0.30)], ctrl_variants=[(300, 1.00)],
                                    ancestral=ancestral)
        assert 300 not in called, f"baseline mismatch called with ancestral={ancestral}"
        assert 120 in called


def test_policy_holds_on_reverse_complement_reads(tmp_path):
    called, _case, _ = _project(tmp_path, [(240, 0.90), (420, 1.00)], ancestral=True, revcomp=True)
    assert {240, 420} <= called


def test_policy_holds_with_partial_coverage(tmp_path):
    """Reads covering only part of the amplicon still yield the high-frequency acquired variants."""
    called, _case, _ = _project(tmp_path, [(240, 0.90), (420, 1.00)], ancestral=True, cover=(100, LEN))
    assert 240 in called


def test_default_is_off(tmp_path):
    """The option defaults to False, so presets and existing projects are untouched."""
    from shmtrace.config import Params
    assert Params().ancestral_reference is False
