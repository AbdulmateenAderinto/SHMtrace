"""Two targets in one project: each read goes to the locus it matches best.

A J_H4 amplicon also carries a J_H4-like paralogue, and a paralogue is worth measuring
rather than discarding. Here a diverged copy of the reference stands in for it, so the
number of reads written from each locus is known.
"""
import random

import pytest

pytest.importorskip("mappy", reason="alignment backend not installed")

from shmtrace import project  # noqa: E402
from shmtrace.io_utils import read_fasta  # noqa: E402

N_TARGET, N_PARALOG = 24, 9
SWAP = {"A": "G", "G": "A", "C": "T", "T": "C"}


def _fastq(path, ref_a, ref_b, n_a, n_b, muts, rng):
    """`muts` is [(position, how many of the target reads carry it), ...]."""
    with open(path, "w") as fh:
        for i in range(n_a + n_b):
            seq = list(ref_a if i < n_a else ref_b)
            if i < n_a:
                for pos, k in muts:                    # a minority of reads, or it is germline
                    if i < k:
                        seq[pos] = SWAP[seq[pos]]
            s = "".join(seq)
            fh.write(f"@{'t' if i < n_a else 'p'}{i}\n{s}\n+\n{'I' * len(s)}\n")


@pytest.fixture(scope="module")
def proj(tmp_path_factory, reference):
    d = tmp_path_factory.mktemp("targets")
    ref = read_fasta(str(reference))[0][1].upper()
    para = list(ref)
    for i in range(0, len(para), 13):                   # ~8% divergent: a clear second locus
        para[i] = SWAP[para[i]]
    a, b = d / "target.fasta", d / "paralog.fasta"
    a.write_text(f">target\n{ref}\n")
    b.write_text(f">paralog\n{''.join(para)}\n")
    rng = random.Random(7)
    samples = []
    for name, muts in (("S1", [(300, 4), (420, 3)]), ("S2", [(300, 5), (505, 4)]), ("C1", []), ("C2", [])):
        fq = d / f"{name}.fastq"
        _fastq(fq, ref, "".join(para), N_TARGET, N_PARALOG, muts, rng)
        samples.append({"name": name, "path": str(fq), "group": "KO" if name.startswith("C") else "WT",
                        "control": name.startswith("C")})
    return {"name": "two targets", "region": "jh4_intron", "platform": "ont", "samples": samples,
            "targets": [{"name": "target", "reference": str(a)}, {"name": "paralog", "reference": str(b)}]}


@pytest.fixture(scope="module")
def res(proj):
    return project.run_targets(proj)


def test_every_read_goes_to_the_locus_it_came_from(res):
    assert res["targets"] == ["target", "paralog"]
    for row in res["samples"]:
        assert row["per_target"]["target"]["clones"] == N_TARGET, row["sample"]
        assert row["per_target"]["paralog"]["clones"] == N_PARALOG, row["sample"]
    assert res["shared_clones"] == {}
    assert not any("aligned equally well" in n for n in res["notes"])


def test_each_target_is_an_ordinary_project_result(res):
    for name in res["targets"]:
        r = res["by_target"][name]
        assert [s["sample"] for s in r["samples"]] == ["S1", "S2", "C1", "C2"]
        assert r["params"]["reference_name"] == name
    calls = {s["sample"]: s["unique_mutations"] for s in res["by_target"]["target"]["samples"]}
    assert calls["S1"] and calls["S2"], "planted mutations should be called at the target"
    assert not any(s["unique_mutations"] for s in res["by_target"]["paralog"]["samples"]), \
        "the paralogue reads carry no mutations of their own"


def test_pooling_the_paralogue_dilutes_the_frequency(res):
    row = next(r for r in res["samples"] if r["sample"] == "S1")
    target_only = row["per_target"]["target"]["frequency_x1e4"]
    assert row["pooled_frequency_x1e4"] < target_only
    assert row["clones_total"] == N_TARGET + N_PARALOG
    assert row["unique_mutations_total"] == row["per_target"]["target"]["unique_mutations"]


def test_public_targets_drops_the_per_sample_detail(res):
    pub = project.public_targets(res)
    assert all(not any(k.startswith("_") for k in r) for r in pub["by_target"].values())


def test_one_target_is_run_project(proj):
    with pytest.raises(ValueError, match="at least two targets"):
        project.run_targets({**proj, "targets": proj["targets"][:1]})
