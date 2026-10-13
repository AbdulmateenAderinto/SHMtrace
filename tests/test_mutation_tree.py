"""shmtrace.mutation_tree: nesting of mutations read from the reads that carry them."""
from shmtrace import mutation_tree


def test_nested_lineage():
    # founder F in every read; branch A in reads 1-6, with a later B in reads 1-3;
    # a second branch C in reads 7-10; read 11 unmutated
    reads = {f"r{i}": {(10, "T")} for i in range(1, 11)}
    reads["r11"] = set()
    for i in range(1, 7):
        reads[f"r{i}"].add((20, "G"))
    for i in range(1, 4):
        reads[f"r{i}"].add((30, "A"))
    for i in range(7, 11):
        reads[f"r{i}"].add((40, "C"))
    reads["r5"].add((50, "A"))                      # one read: left out of the tree
    t = mutation_tree.build(reads)
    by_mut = {n.mutations[0]: n for n in t["nodes"]}
    founder = by_mut[(10, "T")]
    assert founder.parent is None
    assert by_mut[(20, "G")].parent == founder.id
    assert by_mut[(40, "C")].parent == founder.id
    assert by_mut[(30, "A")].parent == by_mut[(20, "G")].id
    assert mutation_tree.depth(t, by_mut[(30, "A")].id) == 3
    assert t["reads_at"][0] == 1 and t["left_out"] == 1 and t["fit"] == 1.0


def test_same_reads_merge_into_one_node():
    reads = {f"r{i}": {(5, "A"), (9, "G")} for i in range(5)}
    t = mutation_tree.build(reads)
    assert len(t["nodes"]) == 1 and t["nodes"][0].mutations == [(5, "A"), (9, "G")]


def test_four_gamete_separates_nested_from_recurrent():
    # nested: B inside A; and C apart from both: no pair shows all four combinations
    nested = {f"r{i}": set() for i in range(12)}
    for i in range(8):
        nested[f"r{i}"].add((1, "A"))
    for i in range(4):
        nested[f"r{i}"].add((2, "G"))
    for i in range(8, 11):
        nested[f"r{i}"].add((3, "T"))
    assert mutation_tree.four_gamete(nested)["incompatible"] == 0
    # recurrent: two mutations combined freely across reads
    rec = {f"r{i}": set() for i in range(12)}
    for i in (0, 1, 2, 3, 4, 5):
        rec[f"r{i}"].add((1, "A"))
    for i in (3, 4, 5, 6, 7, 8):
        rec[f"r{i}"].add((2, "G"))
    fg = mutation_tree.four_gamete(rec)
    assert fg["incompatible"] == 1 and fg["incompatible_pct"] == 100.0
