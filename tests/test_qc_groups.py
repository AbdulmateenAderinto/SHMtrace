"""qc.duplicate_groups: pairwise duplicate flags become one row per sequence, with what ties it together."""
from shmtrace import qc


def _clone(sample, muts):
    return {"sample": sample, "muts": set(muts), "load": len(muts)}


def _flag(a, b):
    return {"check": "duplicate_clone", "clone": a, "other_clone": b}


def test_pairs_become_one_group_across_animals():
    clones = {"a": _clone("K2", [(143, "A"), (166, "A"), (259, "G")]),
              "b": _clone("K2", [(143, "A"), (166, "A"), (259, "G")]),
              "c": _clone("WT4", [(143, "A"), (166, "A"), (259, "G"), (614, "T")])}
    groups = qc.duplicate_groups([_flag("a", "c"), _flag("b", "c")], clones)
    assert len(groups) == 1
    g = groups[0]
    assert g["clones"] == ["a", "b", "c"] and g["samples"] == ["K2", "WT4"]
    assert g["severity"] == "problem"                      # two animals cannot make one clone
    # what every member carries, so the lineage can be looked for elsewhere (a nanopore library)
    assert g["shared_mutations"] == [[143, "A"], [166, "A"], [259, "G"]]


def test_one_animal_is_only_a_check():
    clones = {"a": _clone("W3", [(10, "T"), (20, "C")]), "b": _clone("W3", [(10, "T"), (20, "C")])}
    g = qc.duplicate_groups([_flag("a", "b")], clones)[0]
    assert g["severity"] == "check" and g["shared_mutations"] == [[10, "T"], [20, "C"]]
