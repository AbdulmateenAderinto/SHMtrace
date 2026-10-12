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
