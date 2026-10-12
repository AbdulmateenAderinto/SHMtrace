"""Switch junctions (switch.py) on simulated recombined reads with known joins."""
import random

import pytest

pytest.importorskip("mappy")

from nanoshm import switch  # noqa: E402


def _region(rng, n, motif):
    """A switch-like region: random sequence studded with its repeat motif."""
    out = []
    while sum(map(len, out)) < n:
        out.append(motif if rng.random() < 0.25 else "".join(rng.choice("ACGT") for _ in range(12)))
    return "".join(out)[:n]


def _noisy(s, rng, sub=0.01, indel=0.01):
    out = []
    for b in s:
        u = rng.random()
        if u < indel / 2:
            continue
        if u < indel:
            out.append(rng.choice("ACGT"))
        out.append(rng.choice([x for x in "ACGT" if x != b]) if rng.random() < sub else b)
    return "".join(out)


def _true_mh(donor, acc, bd, ba):
    k = 0
    while k < 30 and donor[bd - 1 - k] == acc[ba - 1 - k]:
        k += 1
    j = 0
    while j < 30 and donor[bd + j] == acc[ba + j]:
        j += 1
    return k + j


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    rng = random.Random(7)
    d = tmp_path_factory.mktemp("sw")
    donor = _region(rng, 3000, "GAGCTGAGCTGGGGTGAGCT")
    acc = _region(rng, 3000, "GGGACCAGTCCTAGCAGCTGT")
    (d / "ref.fa").write_text(f">Smu\n{donor}\n>Sg1\n{acc}\n")
    truth, reads = {}, []
    for i in range(120):
        bd = rng.randint(1200, 2200)
        kind = i % 3
        if kind == 0:                                   # long microhomology: find a match
            k = rng.randint(4, 8)
            mer = donor[bd - k:bd]
            hits = [x + k for x in range(800, 2200) if acc[x:x + k] == mer]
            if not hits:
                continue
            ba, ins = rng.choice(hits), ""
        else:
            ba = rng.randint(800, 2200)
            ins = "".join(rng.choice("ACGT") for _ in range(rng.randint(2, 5))) if kind == 1 else ""
        clean = donor[bd - 700:bd] + ins + acc[ba:ba + 700]
        name = f"r{i}"
        truth[name] = {"bd": bd, "ba": ba, "ins": ins,
                       "mh": 0 if ins else _true_mh(donor, acc, bd, ba)}
        s = _noisy(clean, rng)
        if rng.random() < 0.5:
            s = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        reads.append((name, s))
    # a few reads that are donor only (no switch)
    for i in range(10):
        a = rng.randint(0, 1500)
        reads.append((f"n{i}", _noisy(donor[a:a + 1400], rng)))
    with open(d / "reads.fq", "w") as fh:
        for n, s in reads:
            fh.write(f"@{n}\n{s}\n+\n{'5' * len(s)}\n")
    return switch.analyse(str(d / "reads.fq"), str(d / "ref.fa"), "Smu"), truth


def test_switched_reads_are_found_and_donor_only_reads_are_not(sim):
    res, truth = sim
    found = {j.read for j in res["junctions"] if j.kind == "switch"}
    assert len(found & set(truth)) >= 0.95 * len(truth)
    assert not {f"n{i}" for i in range(10)} & found


def test_breakpoints_are_placed_within_the_microhomology(sim):
    res, truth = sim
    ok = n = 0
    for j in res["junctions"]:
        t = truth.get(j.read)
        if not t:
            continue
        n += 1
        slack = t["mh"] + 3
        ok += abs(j.donor_break - t["bd"]) <= slack and abs(j.acceptor_break - 1 - t["ba"]) <= slack
    assert ok >= 0.9 * n


def test_microhomology_classes_are_recovered(sim):
    res, truth = sim
    agree = n = 0
    for j in res["junctions"]:
        t = truth.get(j.read)
        if not t:
            continue
        n += 1
        agree += switch.mh_class(j.microhomology, j.insertion) == switch.mh_class(t["mh"], t["ins"])
    assert agree >= 0.75 * n


def test_copies_of_one_junction_count_once():
    j = [switch.Junction("a", "Smu", 1500, "Sg1", 900, "switch", 2, "", False, 0, 100),
         switch.Junction("b", "Smu", 1501, "Sg1", 901, "switch", 1, "", False, 0, 100),
         switch.Junction("c", "Smu", 1500, "Sg1", 900, "switch", 2, "", False, 0, 100),
         switch.Junction("d", "Smu", 1700, "Sg1", 900, "switch", 0, "", False, 0, 100)]
    u = switch.unique_junctions(j)
    assert len(u) == 2 and any(x.donor_break == 1500 and x.microhomology == 2 for x in u)
