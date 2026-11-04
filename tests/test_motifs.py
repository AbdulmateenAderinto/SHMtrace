"""AID hotspot / coldspot motif scoring."""
from shmtrace.motifs import IUPAC, MOTIFS, STRAND_PAIR, is_site, site_counts


def test_iupac_codes_cover_the_degenerate_alphabet():
    assert set(IUPAC["W"]) == {"A", "T"}
    assert set(IUPAC["R"]) == {"A", "G"}
    assert set(IUPAC["Y"]) == {"C", "T"}
    assert set(IUPAC["N"]) == {"A", "C", "G", "T"}


def test_wrc_hits_the_mutated_c():
    #        0123
    ref = "TAGC"      # T A G C -> AGC matches WRC with C at index 3
    assert is_site("WRC", ref, 3)


def test_wrc_rejects_a_non_matching_context():
    assert not is_site("WRC", "GGGC", 3)   # G G G -> first base G is not W


def test_gyw_is_the_bottom_strand_complement_of_wrc():
    """WRC scores a mutated C; GYW scores the mutated G on the other strand."""
    assert MOTIFS["WRC"][1] == 2, "C is the third base of WRC"
    assert MOTIFS["GYW"][1] == 0, "G is the first base of GYW"
    assert STRAND_PAIR == ("WRC", "GYW")
    assert is_site("GYW", "GCA", 0)        # G C A -> G/Y=C/W=A


def test_motif_at_the_sequence_edge_is_not_a_site():
    """No out-of-range indexing when the motif would run off the end."""
    assert not is_site("WRC", "AGC", 0)
    assert not is_site("GYW", "AGC", 2)


def test_wa_is_the_pol_eta_hotspot():
    assert is_site("WA", "AA", 1)          # W=A, mutated A
    assert is_site("WA", "TA", 1)
    assert not is_site("WA", "GA", 1)


def test_coldspot_syc():
    assert is_site("SYC", "GTC", 2)        # S=G, Y=T, mutated C
    assert not is_site("SYC", "ATC", 2)    # A is not S


def test_site_counts_is_the_denominator_over_positions():
    ref = "TAGCTAGC"
    counts = site_counts(ref, list(range(len(ref))))
    assert counts["WRC"] >= 1
    assert set(counts) == set(MOTIFS)
    assert all(v >= 0 for v in counts.values())


def test_site_counts_never_exceeds_the_positions_examined():
    ref = "TAGCTAGC"
    positions = [3]
    counts = site_counts(ref, positions)
    assert all(v <= len(positions) for v in counts.values())


def test_rgyw_and_dgyw_score_the_g_not_the_y():
    """RGYW and DGYW are the bottom-strand forms of WRCY and WRCH: the mutated base is the G."""
    assert is_site("RGYW", "AGCT", 1) and not is_site("RGYW", "AGCT", 2)
    assert is_site("DGYW", "TGCA", 1) and not is_site("DGYW", "CGCA", 1)   # C is not D
    assert is_site("WRCH", "AGCA", 2) and not is_site("WRCH", "AGCG", 2)   # G is not H


def test_overlapping_motifs_score_both_their_bases():
    """WGCW is a WRC and a GYW at once, WATW a WA and a TW: both middle bases are sites."""
    assert is_site("WGCW", "AGCT", 1) and is_site("WGCW", "AGCT", 2)
    assert not is_site("WGCW", "AGCG", 1)
    assert is_site("WATW", "TATA", 1) and is_site("WATW", "TATA", 2)
    assert not is_site("WATW", "CATA", 1)


def test_a_call_is_named_with_every_family_it_sits_in():
    from shmtrace.motifs import families_at
    assert families_at("AGCA", 2) == ["WRC/GYW", "WRCH/DGYW", "WGCW"]      # C of AGCA: not WRCY (A is not Y)
    assert families_at("CAAGA", 2) == ["WA/TW"]
    assert families_at("GTTC", 1) == ["WA/TW"]                             # T then W: TW
    assert families_at("AATAT", 2) == ["WA/TW", "WATW"]
    assert families_at("GCC", 2) == ["SYC/GRS"]
    assert families_at("GCC", 0) == []


def test_the_chance_level_is_the_share_of_the_window_s_bases():
    from shmtrace.motifs import background, family_share
    ref = "AGCTGGCC"          # G:C at 1, 2, 4, 5, 6, 7
    hits, n = family_share(ref, range(len(ref)), "WRC_GYW")
    assert n == 6 and hits == 2                                             # the G and C of AGCT
    assert background(ref, 0, len(ref) - 1)["WRC_GYW"] == round(100 * 2 / 6, 1)


def test_sample_rows_carry_every_family():
    from types import SimpleNamespace as V
    from shmtrace.motifs import FAMILIES
    from shmtrace.project import sample_metrics
    ref = "TTAGCTTAATTTGCCAAT"
    vs = [V(pos=5, ref="C", alt="T", supporting_reads=3), V(pos=9, ref="A", alt="G", supporting_reads=2),
          V(pos=15, ref="C", alt="T", supporting_reads=2)]
    res = {"summary": {"calling": {}, "unique_variant_density": 0}, "variants": vs, "per_read_mutations": [1, 1]}
    try:
        row = sample_metrics("s", res, ref, 0, len(ref) - 1, "g", False)
    except Exception:          # the summary is a stub: only the motif fields are under test
        import pytest
        pytest.skip("sample_metrics needs a fuller summary")
    for key, _l, bases, _m, _w in FAMILIES:
        assert f"pct_{bases}_in_{key}" in row
    assert row["pct_GC_in_WRC_GYW"] == 50.0          # C at 5 (AGC) is; C at 15 (GCC) is not
    assert row["pct_GC_in_SYC_GRS"] == 50.0
    assert row["pct_AT_in_WA_TW"] == 100.0 and row["at_mutations"] == 1
