import random

from shmtrace import clones
from shmtrace import lineage_sim as L


def test_clone_sizes_add_up_and_the_dominant_clone_comes_first():
    rng = random.Random(0)
    s = L.clone_sizes(500, 40, 1.0, 0.6, rng)
    assert sum(s) == 500
    assert len(s) == 40
    assert s[0] == 300
    assert all(x >= 1 for x in s)
    s = L.clone_sizes(500, 40, 1.0, 0.0, rng)
    assert sum(s) == 500 and s == sorted(s, reverse=True)


def test_junctions_are_kept_apart():
    j = L.random_junctions(30, 20, random.Random(1), min_dist=8)
    assert len(set(j)) == 30
    assert all(clones._edit(a, b, 7) >= 8 for i, a in enumerate(j) for b in j[i + 1:])


def test_lineage_mutations_are_shared_down_the_tree():
    ref = "ACGTGACTGGCTACGATCGGCTAGCTAGGCTAACGTAGCTAGCATCGACTAGCTAGCATCGATCGATCGA" * 4
    rng = random.Random(3)
    sc = {"clones": 5, "zipf": 1.0, "dominant": 0.0, "trunk": 6.0, "branch": 2.0,
          "private": 0.0}
    t = L.design(60, ref, 10, len(ref) - 10, sc, rng)
    by_clone = {}
    for r, c in t["read_clone"].items():
        by_clone.setdefault(c, []).append(t["read_muts"].get(r, set()))
    for sets in by_clone.values():
        trunk = set.intersection(*sets)
        # with no private mutations, every read of a clone carries its trunk
        # and the union is trunk + branches only
        assert all(trunk <= s for s in sets)
    # different clones draw their trunks independently
    assert len({frozenset(set.intersection(*v)) for v in by_clone.values()}) > 1


def test_the_dominant_clone_can_be_left_unmutated():
    ref = "ACGTGACTGGCTACGATCGGCTAGCTAGGCTAACGTAGCTAGCATCGACT" * 6
    sc = {"clones": 10, "zipf": 1.0, "dominant": 0.5, "dominant_load": 0.0,
          "trunk": 5.0, "branch": 1.0, "private": 0.0}
    t = L.design(100, ref, 5, len(ref) - 5, sc, random.Random(4))
    dom = [r for r, c in t["read_clone"].items() if c == 0]
    assert len(dom) == 50
    assert not any(t["read_muts"].get(r) for r in dom)


def test_adjusted_rand():
    assert L.adjusted_rand([0, 0, 1, 1], [5, 5, 7, 7]) == 1.0      # labels are arbitrary
    assert L.adjusted_rand([0, 0, 1, 1], [0, 1, 0, 1]) < 0
    assert L.adjusted_rand([0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 2]) > 0.5


def test_qvalues_match_benjamini_hochberg():
    q = L.qvalues([0.01, 0.04, 0.03, 0.5])
    assert [round(x, 4) for x in q] == [0.04, 0.0533, 0.0533, 0.5]
    # untested hypotheses enlarge m
    assert L.qvalues([0.01], m=10) == [0.1]


def test_roc_counts_true_and_false_calls_at_each_threshold():
    calling = {"multiple_testing": "none",
               "_candidate_keys": [(11, "T", 1e-7), (21, "A", 0.03), (31, "G", 0.3)]}
    support = {(10, "T"): 5, (20, "A"): 2, (40, "C"): 3}
    pts = {r["alpha"]: r for r in L.roc(calling, support)}
    assert pts[1e-6]["true_calls"] == 1 and pts[1e-6]["false_calls"] == 0
    assert pts[0.05]["true_calls"] == 2 and pts[0.05]["recall"] == round(2 / 3, 4)
    assert pts[0.5]["false_calls"] == 1


def test_bushy_trees_nest_subclones_inside_their_parents():
    ref = "ACGTGACTGGCTACGATCGGCTAGCTAGGCTAACGTAGCTAGCATCGACTAGCTAGCATCGATCGATCGA" * 4
    sc = {"clones": 3, "zipf": 1.0, "dominant": 0.0, "trunk": 2.0, "branch": 1.5,
          "private": 0.0, "levels": 3}
    t = L.design(90, ref, 10, len(ref) - 10, sc, random.Random(4))
    assert all(x["subclones"] >= 2 for x in t["lineages"] if x["reads"] >= 4)
    # reads of one clone carry nested sets: any two either share a subclone's mutations or not
    by_clone = {}
    for r, c in t["read_clone"].items():
        by_clone.setdefault(c, []).append(t["read_muts"].get(r, set()))
    assert any(len({frozenset(s) for s in v}) > 2 for v in by_clone.values())


def test_amplify_copies_a_molecules_mutations_and_keeps_pcr_errors_out_of_the_truth():
    ref = "ACGTGACTGGCTACGATCGGCTAGCTAGGCTAACGTAGCTAGCATCGACTAGCTAGCATCGATCGATCGA" * 4
    rng = random.Random(5)
    t = L.design(40, ref, 10, len(ref) - 10, L.SCENARIOS["polyclonal"], rng)
    a = L.amplify(t, 200, rng, ref, 10, len(ref) - 10, skew=1.0, early_error=0.01)
    assert len(a["read_clone"]) == 200
    for j, m in a["read_molecule"].items():
        assert a["read_muts"][j] == t["read_muts"].get(m, set())
        assert a["read_clone"][j] == t["read_clone"][m]
        assert not {p for p, _ in a["read_noise"][j]} & {p for p, _ in a["read_muts"][j]}
    # copies of one molecule share its PCR errors
    mols = {}
    for j, m in a["read_molecule"].items():
        mols.setdefault(m, set()).add(frozenset(a["read_noise"][j]))
    assert all(len(v) == 1 for v in mols.values())
    assert sum(map(len, a["read_noise"].values())) > 0


def test_quality_map_moves_qualities_by_quantile():
    ref_q = [40] * 90 + [20] * 10
    alt_q = [38] * 80 + [15] * 20
    m = L.quality_map(ref_q, alt_q)
    assert m[40] == 38 and m[20] == 15
    assert L.quality_map([], alt_q) == {}
