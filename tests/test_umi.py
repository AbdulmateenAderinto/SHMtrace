"""UMI extraction, grouping and consensus (umi.py), and an end-to-end simulation."""
import random
from dataclasses import replace

import pytest

from shmtrace import umi
from shmtrace.align import AlignedRead
from shmtrace.config import Params

REF = "ACGTTGCAAGGCTTACCGATGCATGCAAGTCCGATTACGGATCCATGGTACCAGTTCAGGCATGCTAGCTAGGATCCAT" * 3


def test_semiglobal_finds_an_anchor_with_errors():
    text = "NNNNNNNNACGTTGCAAGGCTTACCGATNNNN"
    ed, s, e = umi._semiglobal(text, "ACGTTGCAAGCTTACCGAT")     # one base deleted
    assert ed == 1 and s == 8 and e == 28


def test_left_and_right_umis_are_read_next_to_the_reference_ends():
    left, right = "AAAACCCCGGGGTTTTAC", "TTTTGGGGCCCCAAAAGT"
    seq = left + REF + right
    r = AlignedRead("r", "ref", 0, len(REF), [(i, b, 30) for i, b in enumerate(REF)], seq, 1,
                    [], len(left), len(left) + len(REF))
    p = replace(Params(), umi="flank", umi_ends="both")
    assert umi.extract(r, REF, p) == left + right
    assert umi.extract(r, REF, replace(p, umi_ends="right")) == right


def test_header_umis():
    r = AlignedRead("read7_ACGTACGTAC", "ref", 0, 10, [], "", 1)
    assert umi.extract(r, REF, replace(Params(), umi="header")) == "ACGTACGTAC"


def test_grouping_tolerates_errors_but_keeps_molecules_apart():
    rng = random.Random(1)
    mols = ["".join(rng.choice("ACGT") for _ in range(18)) for _ in range(30)]
    keys = {}
    for m, u in enumerate(mols):
        for c in range(5):
            s = list(u)
            for _ in range(c % 3):                      # 0-2 substitutions per copy
                i = rng.randrange(18)
                s[i] = rng.choice([x for x in "ACGT" if x != s[i]])
            keys[f"m{m}c{c}"] = "".join(s)
    g = umi.group(keys, 3)
    for m in range(30):
        assert len({g[f"m{m}c{c}"] for c in range(5)}) == 1
    assert len(set(g.values())) == 30


def test_a_truncated_or_indel_umi_joins_its_centroid():
    # one base short (cut at the read end), one deleted inside, one inserted inside:
    # each is within the tolerance of the centroid and must not become a molecule of
    # its own, though the pieces it shares are shifted or of another length. Each is
    # seen 3 times, so no directional step could rescue it.
    keys = {f"c{i}": "ACGTACGTACGT" for i in range(4)}
    for v in ("ACGTACGTACG", "ACGTACTACGT", "ACGTAACGTACGT"):
        keys.update({f"{v}{i}": v for i in range(3)})
    assert len(set(umi.group(keys, 2).values())) == 1
    rng = random.Random(4)
    u = "".join(rng.choice("ACGT") for _ in range(18))
    keys = {f"c{i}": u for i in range(4)}
    for v in (u[2:], u[:9] + u[10:], u[:-1]):
        keys.update({f"{v}{i}": v for i in range(3)})
    assert len(set(umi.group(keys, 3).values())) == 1


def test_a_centroid_too_short_to_cut_is_still_found():
    keys = {"a": "ACGTACG", "b": "ACGTACG", "c": "ACGTTCG", "d": "ACGTACGTAC"}
    g = umi.group(keys, 3)            # pieces from the commonest length 7 are 1 nt
    assert g["c"] == g["a"]
    keys = {f"s{i}": "ACGTA" for i in range(3)} | {f"l{i}": "TTTTTTTTTTTT" for i in range(5)}
    keys["s_err"] = "ACGTT"
    g = umi.group(keys, 2)            # pieces of 4 from length 12; "ACGTA" can't be cut
    assert g["s_err"] == g["s0"] and g["s0"] != g["l0"]


def test_directional_step_is_off_unless_asked():
    big = "AAAACCCCGGGGTTTTAC"
    far = "AAAACCCCGGGGTTTTAC"[:13] + "ACGAG"      # 4 edits: beyond max_edits 3
    keys = {f"b{i}": big for i in range(6)} | {"rare": far}
    assert umi.group(keys, 3)["rare"] != umi.group(keys, 3)["b0"]
    g = umi.group(keys, 3, directional=1)
    assert g["rare"] == g["b0"]


def test_consensus_outvotes_errors_and_keeps_a_shared_mutation():
    ref = REF[:60]
    reads = []
    for k in range(5):
        pairs = []
        for i, b in enumerate(ref):
            if i == 30:
                b = "A" if ref[i] != "A" else "C"           # the molecule's mutation
            elif i == 10 + k:
                b = "T" if ref[i] != "T" else "G"           # a different error in each copy
            pairs.append((i, b, 20))
        reads.append(AlignedRead(f"c{k}", "ref", 0, 60, pairs, "", 1))
    c = umi._consensus(reads, ref, "umi0")
    got = {rp: (b, q) for rp, b, q in c.pairs}
    assert got[30][0] != ref[30] and got[30][1] >= 40       # kept, confidently
    assert all(got[10 + k][0] == ref[10 + k] for k in range(5))   # errors outvoted


def test_umi_consensus_end_to_end(examples):
    """A simulated UMI library from real control reads: every planted mutation is
    found, single-molecule ones included, and nothing else is called."""
    pytest.importorskip("mappy")
    import os
    import tempfile

    from shmtrace import lineage_sim as L
    from shmtrace import pipeline, project
    from shmtrace.io_utils import read_fasta, read_fastq
    d = examples.parent / "examples_compare"
    ref_f, con = str(d / "JH4_reference.fasta"), [str(d / "JH4_contaminant.fasta")]
    p = project.build_params({"region": "jh4_intron"})
    ref = read_fasta(ref_f)[0][1]
    lo, hi = pipeline._window(p, ref)
    rng = random.Random(3)
    tmp = tempfile.mkdtemp()

    def lib(src, scen, name):
        aligned, _ = pipeline.prepare_reads(src, ref_f, con, p)
        ok = {r.qname for r in aligned if r.r_en >= len(ref) - 3}
        reads = [x for x in read_fastq(src) if x[0] in ok]
        t = L.design(len(reads) // 6, ref, lo, hi, L.SCENARIOS[scen], rng, None,
                     p.junction_length, p.junction_ref_start - 1)
        out, carried = L.build_umi_library(reads, ref_f, t, p, rng, ends="right")
        f = os.path.join(tmp, name)
        with open(f, "w") as fh:
            for n, s, q in out:
                fh.write(f"@{n}\n{s}\n+\n{q}\n")
        return f, carried
    sim, carried = lib(str(d / "KO3.fastq"), "polyclonal", "sim.fq")
    nulls = [lib(str(d / f"KO{i}.fastq"), "null", f"n{i}.fq")[0] for i in (1, 2)]
    r = project.run_project({
        "name": "u", "region": "jh4_intron", "reference": ref_f, "contaminants": con,
        "overrides": {"umi": "flank", "umi_ends": "right"},
        "samples": [{"name": "SIM", "path": sim, "group": "S"}]
        + [{"name": f"C{i}", "path": f, "group": "C", "control": True} for i, f in enumerate(nulls)]})
    row = next(s for s in r["samples"] if s["sample"] == "SIM")
    sc = L.score_umi(r["_details"]["SIM"], row, carried, lo, hi)
    assert sc["molecules_analysed"] >= 15
    assert sc["group_purity"] >= 0.98
    assert sc["recall"] >= 0.95 and sc["false_calls"] <= 1


def test_anchor_search_accepts_degenerate_primer_bases():
    ed, s, e = umi._semiglobal("NNNNAGAGTTTGATCATGGCTCAGNNN", "AGRGTTYGATYMTGGCTCAG")
    assert ed == 0 and s == 4 and e == 24


def test_directional_step_is_on_for_dual_umis_only():
    p = replace(Params(), umi="flank")
    assert umi.directional(replace(p, umi_ends="both")) == 2
    assert umi.directional(replace(p, umi_ends="right")) == 0
    assert umi.directional(replace(p, umi="header")) == 0
    assert umi.directional(replace(p, umi_ends="right", umi_directional=1)) == 1
