"""Strand-aware and quality-calibrated nulls, exact tests, strand filter, homopolymers."""
import random
from dataclasses import replace

import pytest

from shmtrace.align import AlignedRead
from shmtrace.config import Params
from shmtrace.null_model import NullModel, estimate_rho, fit_aligned, load_panel, save_panel
from shmtrace.pipeline import mode_overrides
from shmtrace.variants import (
    _betabinom_sf,
    _binom_sf,
    _pbinom_sf,
    _poisson_sf,
    call_variants,
    homopolymer_run,
)

REF = "ACGTTGCAAGCTAGCTTACGGATCCATGCAAATTTCGCGATCGTAGCTAGGCTAACGTTAGC"
A10 = "A" if REF[10] != "A" else "G"          # alternate bases at the two test sites
A30 = "C" if REF[30] != "C" else "G"


def _ref(tmp_path):
    f = tmp_path / "ref.fasta"
    f.write_text(f">ref\n{REF}\n")
    return str(f)


def _reads(n, strand_of=lambda i: 1 if i % 2 else -1, edits=None, q=30, qmap=None):
    """n full-length reads; edits = {read index: {pos0: base}}."""
    edits = edits or {}
    out = []
    for i in range(n):
        seq = list(REF)
        for pos, b in edits.get(i, {}).items():
            seq[pos] = b
        pairs = [(j, seq[j], (qmap or {}).get((i, j), q)) for j in range(len(REF))]
        out.append(AlignedRead(f"r{i}", "ref", 0, len(REF), pairs, "".join(seq), strand_of(i)))
    return out


def _edits(*groups):
    """Merge (read indices, pos0, base) groups into {read: {pos0: base}}."""
    out = {}
    for idx, pos, base in groups:
        for i in idx:
            out.setdefault(i, {})[pos] = base
    return out


def _params(**kw):
    return replace(Params(), downsample_to=None, phase_mode="off", **kw)


# ---- exact tests -------------------------------------------------------------
def test_poisson_binomial_matches_binomial_when_probabilities_are_equal():
    assert _pbinom_sf(3, [0.01] * 200) == pytest.approx(_binom_sf(3, 200, 0.01), rel=1e-6)


def test_poisson_is_a_close_approximation_for_small_error():
    probs = [0.001] * 500
    assert _pbinom_sf(2, probs) == pytest.approx(_poisson_sf(2, 0.5), rel=0.05)


def test_beta_binomial_is_binomial_without_overdispersion_and_wider_with_it():
    assert _betabinom_sf(4, 500, 0.002, 0.0) == pytest.approx(_binom_sf(4, 500, 0.002))
    assert _betabinom_sf(4, 500, 0.002, 0.01) > _betabinom_sf(4, 500, 0.002, 0.0)


# ---- homopolymers ------------------------------------------------------------
def test_homopolymer_run_counts_runs_touching_the_position():
    seq = "ACAAAAGT"
    assert homopolymer_run(seq, 3) == 4          # inside AAAA
    assert homopolymer_run(seq, 6) == 4          # G right after the run
    assert homopolymer_run("ACGTACGT", 2) == 1


# ---- strand-aware null ---------------------------------------------------------
def test_fit_records_strand_specific_rates(tmp_path):
    # position 10 miscalled as T on the minus strand only, in both control files
    files = [_reads(200, edits={i: {10: A10} for i in range(0, 40, 2)}) for _ in range(2)]
    null, stats = fit_aligned(files, REF, _params(), 0, len(REF) - 1, calibrate=False)
    plus, minus = null.strand_rates[(10, A10, 1)], null.strand_rates[(10, A10, -1)]
    assert minus > 10 * plus
    assert stats["strand_specific_error_sites"] >= 1


def test_strand_model_discounts_an_error_seen_on_one_strand_in_controls(tmp_path):
    ref = _ref(tmp_path)
    ctrl = [_reads(200, edits={i: {10: A10} for i in range(0, 40, 2)}) for _ in range(2)]
    null, _ = fit_aligned(ctrl, REF, _params(), 0, len(REF) - 1, calibrate=False)
    # sample: same minus-strand error, plus a real mutation at 30 on both strands
    sample = _reads(200, edits=_edits((range(0, 40, 2), 10, A10), (range(100, 120), 30, A30)))
    p = _params(error_model="strand", multiple_testing="none", max_pvalue=1e-3)
    calls = {(v.pos, v.alt) for v in call_variants(sample, ref, "ref", p, null)[0]}
    assert (31, A30) in calls
    assert (11, A10) not in calls


def test_strand_model_without_controls_falls_back(tmp_path):
    ref = _ref(tmp_path)
    p = _params(error_model="strand")
    _v, stats = call_variants(_reads(50), ref, "ref", p, None)
    assert stats["error_model"] == "read_quality"
    assert "notes" in stats


def test_overdispersion_is_zero_for_identical_libraries_and_positive_for_variable_ones():
    same = [{("k",): (5, 1000)}, {("k",): (5, 1000)}]
    varied = [{("k",): (0, 1000)}, {("k",): (40, 1000)}]
    assert estimate_rho(same, [("k",)]) == 0.0
    assert estimate_rho(varied, [("k",)]) > 0.0


# ---- quality-calibrated null ---------------------------------------------------
def test_qcal_trusts_high_quality_bases_only_as_far_as_controls_justify(tmp_path):
    ref = _ref(tmp_path)
    # controls: low-quality A>G miscalls everywhere (Q8), never at high quality
    rng = random.Random(0)
    ctrl = []
    for _f in range(2):
        edits, qmap = {}, {}
        for i in range(300):
            if rng.random() < 0.3:
                pos = rng.randrange(5, len(REF) - 5)
                if REF[pos] == "A":
                    edits[i] = {pos: "G"}
                    qmap[(i, pos)] = 8
        ctrl.append(_reads(300, edits=edits, qmap=qmap))
    null, _ = fit_aligned(ctrl, REF, _params(), 0, len(REF) - 1, calibrate=False)
    p = _params(error_model="qcal", multiple_testing="none", max_pvalue=1e-3)
    pos = REF.index("AAA") + 1                    # an A
    hi_q = _reads(300, edits={i: {pos: "G"} for i in (1, 2, 4)},
                  qmap={(i, pos): 45 for i in (1, 2, 4)})
    lo_q = _reads(300, edits={i: {pos: "G"} for i in (1, 2, 4)},
                  qmap={(i, pos): 8 for i in (1, 2, 4)})
    assert (pos + 1, "G") in {(v.pos, v.alt) for v in call_variants(hi_q, ref, "ref", p, null)[0]}
    assert (pos + 1, "G") not in {(v.pos, v.alt)
                                  for v in call_variants(lo_q, ref, "ref", p, null)[0]}


def test_qscale_inflates_qcal_rates():
    null = NullModel(qtail={("A", "G", 1): [1e-5] * 61})
    base = null.qrate(5, "A", "G", 1, 40, 500)
    null.qscale = 10
    assert null.qrate(5, "A", "G", 1, 40, 500) == pytest.approx(10 * base)


def test_panel_round_trip(tmp_path):
    null = NullModel({(3, "T"): 0.01}, {(3, "T", 1): 0.02}, 0.003,
                     qtail={("C", "T", 1): [0.1] * 61}, site_quals={(3, "T", 1): [30, 12]},
                     site_depth={(3, 1): 400}, qscale=5.0)
    path = tmp_path / "panel.json"
    save_panel(null, str(path), basecaller="sup@v5")
    back = load_panel(str(path))
    assert back == null
    assert back.strand_rates == null.strand_rates and back.rho == null.rho
    # a site's control qualities are held in ascending order, so qrate can bisect
    assert back.qscale == 5.0 and back.site_quals == {(3, "T", 1): [12, 30]}
    assert back.meta["basecaller"] == "sup@v5"


# ---- strand filter -------------------------------------------------------------
def test_strand_bias_filter_drops_one_strand_support(tmp_path):
    ref = _ref(tmp_path)
    edits = _edits((range(1, 40, 2), 10, A10),        # 20 plus-strand reads only
                   (range(0, 40), 30, A30))           # 20 + 20
    p = _params(multiple_testing="none", strand_filter="bias")
    v, stats = call_variants(_reads(200, edits=edits), ref, "ref", p, None)
    calls = {(x.pos, x.alt) for x in v}
    assert (31, A30) in calls and (11, A10) not in calls
    assert stats["strand_filter_removed"] >= 1


def test_both_strands_filter_needs_a_read_on_each_strand(tmp_path):
    ref = _ref(tmp_path)
    edits = _edits(((1, 3, 5), 10, A10), ((1, 2, 4), 30, A30))
    p = _params(multiple_testing="none", max_pvalue=1.0, strand_filter="both")
    calls = {(x.pos, x.alt) for x in call_variants(_reads(100, edits=edits), ref, "ref", p)[0]}
    assert (31, A30) in calls and (11, A10) not in calls


def test_mode_overrides_accept_filter_suffix():
    assert mode_overrides("qcal+filter") == {"error_model": "qcal", "strand_filter": "bias"}
    with pytest.raises(ValueError):
        mode_overrides("nonsense")
