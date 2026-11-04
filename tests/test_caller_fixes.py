"""The three opt-in caller fixes from docs/calibration/ACCURACY_BREADTH.md.

Each is off by default, so the rest of the suite covers the unchanged behaviour;
these tests show that each one does what it claims when it is switched on.

  read_end_margin_bp  the columns where reads begin are covered by fewer mutant
                      reads than mutation-free ones, because a mismatch near a
                      read's end is cheaper for the aligner to clip than to
                      align. Dropping the same number of columns from every read
                      takes that deficit out of the denominator too.
  mask_control_het    a second allele in a control sits just under the 50%
                      germline cut, so it stays in the test; the mask takes those
                      (position, allele) pairs out of the samples' candidates and
                      the controls' decoys alike.
  flat qualities      a library whose bases all carry one quality leaves the
                      single-read estimate nothing to separate; it is flagged.
"""
import random

from shmtrace import null_model, singletons
from shmtrace.align import AlignedRead, pair_span
from shmtrace.config import Params
from shmtrace.singletons import FLAT_QUALITY_SHARE
from shmtrace.variants import call_variants

REF = ("ACGTTGCAAGCTAGCTAGGATCCATGCATGCAAATTTGGGCCCATATCGCGATCGATCGTAGCTAGCATGCA"
       "GGATCCTTAAGGCCATCGATCGAT")          # 96 bp
MARGIN = 8


def _ref_fasta(tmp_path):
    fa = tmp_path / "ref.fasta"
    fa.write_text(f">ref\n{REF}\n")
    return str(fa)


def _read(name, start, end, muts=(), q=40, strand=1):
    """A read aligned over [start, end), reference bases except at `muts`."""
    alt = dict(muts)
    pairs = [(i, alt.get(i, REF[i]), q) for i in range(start, end)]
    return AlignedRead(name, "ref", start, end, pairs, "", strand)


def _params(**kw):
    base = dict(analysis_window=(1, len(REF)), error_model="auto", multiple_testing="none",
                max_pvalue=1.0, min_variant_reads=1, min_base_quality=20)
    base.update(kw)
    return Params(**base)


def test_pair_span_drops_the_same_number_of_columns_at_each_read_end():
    r = _read("r", 0, len(REF))
    assert pair_span(r, 0) == (0, len(REF))
    assert pair_span(r, MARGIN) == (MARGIN, len(REF) - MARGIN)
    # a read with nothing left over contributes no bases at all, rather than
    # a negative span that would silently wrap round
    short = _read("s", 0, 2 * MARGIN)
    assert pair_span(short, MARGIN) == (0, 0)
    assert pair_span(_read("t", 0, 2 * MARGIN + 1), MARGIN) == (MARGIN, MARGIN + 1)


def test_the_margin_takes_the_columns_where_reads_begin_out_of_the_denominator(tmp_path):
    """Reads that begin inside the window leave those columns biased.

    Every read here begins at column 0 and the library carries one single-read
    mutation at each of the 96 positions, so the arithmetic is exact. A molecule
    whose mutation falls in its first 8 aligned columns is clipped back past it,
    as minimap2 clips a mismatch near a read's end: the mutation is gone and so
    are that read's bases at those columns. The frequency the caller reports is
    then low, because the clipped reads are missing from the numerator while
    their mutation-free neighbours are still in the denominator.
    """
    fa = _ref_fasta(tmp_path)
    n = len(REF)
    reads = []
    for pos in range(n):
        alt = "A" if REF[pos] != "A" else "C"
        strand = 1 if pos % 2 else -1
        if pos < MARGIN:
            # clipped: the alignment starts past the mutation, which is lost
            reads.append(_read(f"c{pos}", MARGIN, n, strand=strand))
        else:
            reads.append(_read(f"m{pos}", 0, n, [(pos, alt)], strand=strand))
    reads += [_read(f"r{i}", 0, n, strand=1 if i % 2 else -1) for i in range(1000 - n)]
    truth = n / (1000 * n)

    def frequency(margin):
        variants, stats = call_variants(reads, fa, "ref", _params(read_end_margin_bp=margin))
        copies = sum(v.supporting_reads for v in variants)
        return copies / stats["analysed_bases"]

    plain, trimmed = frequency(0), frequency(MARGIN)
    # 88 of 96 mutations survive the clip, in 95,936 bases: 0.917 of the truth
    assert plain < 0.93 * truth
    # with the margin, 80 mutations in 79,936 bases: the biased columns have left
    # the analysis on both sides of the ratio
    assert abs(trimmed - truth) < 0.01 * truth


def test_the_margin_warns_that_the_fixed_window_frequency_reads_low(tmp_path):
    """The fixed-window frequency credits every read with the whole window, so trimming
    columns takes mutations out of its numerator and nothing out of its denominator."""
    fa = _ref_fasta(tmp_path)
    reads = [_read(f"r{i}", 0, len(REF), strand=1 if i % 2 else -1) for i in range(20)]
    _v, off = call_variants(reads, fa, "ref", _params())
    assert not any("read_end_margin_bp" in n for n in off.get("notes", []))
    _v, on = call_variants(reads, fa, "ref", _params(read_end_margin_bp=MARGIN))
    note = next(n for n in on["notes"] if "read_end_margin_bp" in n)
    assert "mutation_frequency_observed_x1e3" in note
    assert on["read_end_margin_bp"] == MARGIN


def test_the_margin_does_not_move_a_call_away_from_the_read_ends(tmp_path):
    fa = _ref_fasta(tmp_path)
    mid = len(REF) // 2
    alt = "A" if REF[mid] != "A" else "C"
    reads = [_read(f"m{i}", 0, len(REF), [(mid, alt)], strand=1 if i % 2 else -1)
             for i in range(5)]
    reads += [_read(f"r{i}", 0, len(REF), strand=1 if i % 2 else -1) for i in range(95)]
    for margin in (0, MARGIN):
        variants, _stats = call_variants(reads, fa, "ref", _params(read_end_margin_bp=margin))
        called = [(v.pos, v.alt, v.supporting_reads, v.depth) for v in variants]
        assert called == [(mid + 1, alt, 5, 100)]


def _het_controls(het_pos, alt, n=100, fraction=0.46):
    """Three control libraries each carrying a second allele at `het_pos`."""
    libs = []
    for c in range(3):
        reads = []
        for i in range(n):
            carries = i < round(fraction * n)
            reads.append(_read(f"c{c}_{i}", 0, len(REF),
                               [(het_pos, alt)] if carries else (),
                               strand=1 if i % 2 else -1))
        libs.append(reads)
    return libs


def test_the_controls_second_allele_is_recorded_only_when_the_mask_is_on():
    het_pos = 30
    alt = "A" if REF[het_pos] != "A" else "C"
    libs = _het_controls(het_pos, alt)
    off, _s = null_model.fit_aligned(libs, REF, _params(), 0, len(REF) - 1, calibrate=False)
    assert off.het_sites == set()
    on, stats = null_model.fit_aligned(libs, REF, _params(mask_control_het=True),
                                       0, len(REF) - 1, calibrate=False)
    assert on.het_sites == {(het_pos, alt)}
    assert stats["control_het_sites"] == 1
    assert stats["control_het_alleles"] == [(het_pos + 1, alt)]
    # the frequency rule still owns anything at 50% or more: that is germline,
    # counted as such, and not reported as a second allele
    fixed = _het_controls(het_pos, alt, fraction=0.9)
    high, hstats = null_model.fit_aligned(fixed, REF, _params(mask_control_het=True),
                                          0, len(REF) - 1, calibrate=False)
    assert high.het_sites == set() and hstats["germline_sites_excluded"] > 0


def test_a_site_only_error_reaches_is_not_taken_for_a_second_allele():
    # 2% of reads mismatching is error, well below control_het_min_frequency
    libs = _het_controls(40, "A" if REF[40] != "A" else "C", fraction=0.02)
    null, _s = null_model.fit_aligned(libs, REF, _params(mask_control_het=True),
                                      0, len(REF) - 1, calibrate=False)
    assert null.het_sites == set()


def test_the_masked_allele_is_not_tested_in_the_sample_but_a_mutation_still_is(tmp_path):
    fa = _ref_fasta(tmp_path)
    het_pos, mut_pos = 30, 60
    het_alt = "A" if REF[het_pos] != "A" else "C"
    mut_alt = "A" if REF[mut_pos] != "A" else "C"
    libs = _het_controls(het_pos, het_alt)
    p = _params(mask_control_het=True)
    null, _s = null_model.fit_aligned(libs, REF, p, 0, len(REF) - 1, calibrate=False)

    # the sample carries the same second allele, at a slightly different ratio,
    # plus five reads of a real mutation elsewhere
    sample = []
    for i in range(100):
        muts = [(het_pos, het_alt)] if i < 48 else []
        if 50 <= i < 55:
            muts = muts + [(mut_pos, mut_alt)]
        sample.append(_read(f"s{i}", 0, len(REF), muts, strand=1 if i % 2 else -1))

    variants, stats = call_variants(sample, fa, "ref", p, null)
    assert (het_pos + 1, het_alt) not in [(v.pos, v.alt) for v in variants]
    assert (mut_pos + 1, mut_alt) in [(v.pos, v.alt) for v in variants]
    assert stats["control_het_removed"] == 1
    assert any("second allele" in n for n in stats["notes"])

    # without the mask the same site is a candidate and is tested
    off, _s2 = null_model.fit_aligned(libs, REF, _params(), 0, len(REF) - 1, calibrate=False)
    v_off, s_off = call_variants(sample, fa, "ref", _params(), off)
    assert (het_pos + 1, het_alt) in [(v.pos, v.alt) for v in v_off]
    assert s_off["control_het_removed"] == 0


def test_an_ancestral_reference_keeps_its_high_frequency_changes():
    """With the reference declared the construct that went into the cells, a change
    driven high is the experiment's subject, so neither rule may remove it."""
    het_pos = 30
    alt = "A" if REF[het_pos] != "A" else "C"
    libs = _het_controls(het_pos, alt)
    p = _params(mask_control_het=True, ancestral_reference=True)
    null, stats = null_model.fit_aligned(libs, REF, p, 0, len(REF) - 1, calibrate=False)
    assert null.het_sites == set() and stats["germline_sites_excluded"] == 0


def test_a_saved_panel_carries_the_masked_alleles(tmp_path):
    het_pos = 30
    alt = "A" if REF[het_pos] != "A" else "C"
    libs = _het_controls(het_pos, alt)
    null, _s = null_model.fit_aligned(libs, REF, _params(mask_control_het=True),
                                      0, len(REF) - 1, calibrate=False)
    path = tmp_path / "panel.json"
    null_model.save_panel(null, str(path))
    assert null_model.load_panel(str(path)).het_sites == {(het_pos, alt)}


def test_the_qscale_search_stops_at_its_published_ceiling_by_default():
    """A project needing more than 50-fold inflation used to be given 50 and reported no
    differently from one that fitted. The grid now goes further, but only when asked: a
    higher factor is a stricter test, and raising it costs real calls."""
    grid = null_model.QSCALE_GRID
    assert grid[:7] == (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0)   # unchanged prefix
    assert grid[6] == null_model.QSCALE_SUSPECT
    assert list(grid) == sorted(grid) and grid[-1] > 50.0
    # what the default cap admits is exactly the old grid, so no published run moves
    assert Params().qscale_max == 50.0
    allowed = [c for c in grid if c <= Params().qscale_max]
    assert allowed == [1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0]
    assert [c for c in grid if c <= 1000.0] == list(grid)


def test_an_unsatisfiable_calibration_says_so_instead_of_reporting_a_fit():
    """Two controls that disagree, with the cap left at its smallest factor: the search
    runs out and must warn, rather than returning the ceiling as if it had worked."""
    rng = random.Random(5)
    noisy, clean = [], []
    for i in range(120):
        # the noisy library miscalls at high quality all over the reference; the clean one
        # does not, so no inflation this search is allowed to try reconciles them
        pairs = [(pos, rng.choice([b for b in "ACGT" if b != REF[pos]]) if rng.random() < 0.04
                  else REF[pos], 40) for pos in range(len(REF))]
        noisy.append(AlignedRead(f"n{i}", "ref", 0, len(REF), pairs, "", 1 if i % 2 else -1))
        clean.append(_read(f"c{i}", 0, len(REF), strand=1 if i % 2 else -1))
    p = _params(error_model="qcal", qscale_max=1.0)
    factor, stats = null_model.calibrate_qscale([noisy, clean], REF, "ref", p, 0, len(REF) - 1)
    assert factor == 1.0                       # the only factor the cap admits
    assert "grid maximum" in stats["qscale_source"]
    assert "does not fit this data" in stats["qscale_warning"]
    assert "raise qscale_max" in stats["qscale_warning"]


def _flat_tally(q, n=200):
    reads = [_read(f"f{i}", 0, len(REF), q=q, strand=1 if i % 2 else -1) for i in range(n)]
    return singletons.tally(reads, REF, 0, len(REF) - 1)


def test_one_quality_on_every_base_is_reported_as_flat():
    q = singletons.quality_resolution([_flat_tally(30) for _ in range(3)])
    assert q["flat"] and q["top_share"] == 1.0 and q["top_quality"] == 30 and q["levels"] == 1


def test_two_levels_are_flat_even_when_neither_dominates():
    """The shape the two public Illumina sets have: Q30 where one mate covered the
    base and Q60 where both did. Neither value says how right the base is."""
    merged = [_read(f"a{i}", 0, len(REF), q=60, strand=1 if i % 2 else -1) for i in range(176)]
    single = [_read(f"b{i}", 0, len(REF), q=30, strand=1 if i % 2 else -1) for i in range(24)]
    q = singletons.quality_resolution([singletons.tally(merged + single, REF, 0, len(REF) - 1)])
    assert q["levels"] == 2 and q["top_share"] < FLAT_QUALITY_SHARE and q["flat"]


def test_a_spread_of_qualities_is_not_flat():
    rng = random.Random(3)
    reads = [_read(f"v{i}", 0, len(REF), strand=1 if i % 2 else -1) for i in range(200)]
    for r in reads:
        r.pairs = [(pos, b, rng.choice([20, 25, 30, 35, 40])) for pos, b, _q in r.pairs]
    q = singletons.quality_resolution([singletons.tally(reads, REF, 0, len(REF) - 1)])
    assert not q["flat"] and q["levels"] == 5 and q["top_share"] < 0.3


def test_the_single_read_estimate_carries_the_flat_quality_warning():
    # flat Q30: the cut steps down to it, keeps every base, and separates nothing
    sample, controls = _flat_tally(30), [_flat_tally(30) for _ in range(3)]
    est = singletons.estimate(sample, controls, 35, 0.5, REF)
    assert est["quality_flat"] is True
    assert est["quality_levels"] == 1
    assert est["efficiency"] == 1.0
    assert "cannot separate a real base from an error" in est["warning"]


def test_a_library_with_real_qualities_is_not_warned_about():
    """A basecaller's own spread of qualities, as the nanopore libraries have."""
    rng = random.Random(11)
    levels = [20, 25, 30, 33, 36, 38, 40, 42]

    def lib(tag):
        reads = []
        for i in range(200):
            pairs = [(pos, b, rng.choice(levels)) for pos, b in enumerate(REF)]
            reads.append(AlignedRead(f"{tag}{i}", "ref", 0, len(REF), pairs, "",
                                     1 if i % 2 else -1))
        return singletons.tally(reads, REF, 0, len(REF) - 1)

    est = singletons.estimate(lib("s"), [lib(f"c{c}") for c in range(3)], 35, 0.5, REF)
    assert est["quality_levels"] == len(levels)
    assert "quality_flat" not in est and "warning" not in est
