"""Quantification against known truth (docs/calibration/QUANT_ACCURACY.md)."""
import pytest

pytest.importorskip("mappy", reason="alignment backend not installed")

from nanoshm import lineage_sim  # noqa: E402
from nanoshm.align import AlignedRead, drop_contaminant_alleles  # noqa: E402


def _read(name, calls):
    return AlignedRead(name, "ref", 0, 100, [(i, calls.get(i, "A"), 40) for i in range(100)], "")


def test_reads_carrying_contaminant_alleles_are_dropped_but_germline_is_not():
    sites = {i: frozenset("C") for i in range(0, 60, 3)}           # 20 diagnostic sites
    germ = {0: "C", 3: "C"}                                        # every read: germline
    reads = [_read(f"r{i}", germ) for i in range(40)]
    reads += [_read("mutated", {**germ, 6: "C", 9: "C", 50: "T"})]  # 2 contaminant bases
    reads += [_read("paralogue", {i: "C" for i in range(0, 60, 3)})]
    reads += [_read("chimera", {**germ, **{i: "C" for i in range(30, 60, 3)}})]
    fates = {}
    kept, dropped = drop_contaminant_alleles(reads, sites, 10, 0.2, fates)
    assert dropped == 2 and set(fates) == {"paralogue", "chimera"}
    assert "mutated" in {r.qname for r in kept}


def test_densities_match_the_truth_when_the_library_is_deeper_than_the_analysis(examples):
    """Built from 497 control reads and analysed at 100: the single-read estimate
    must count events in the analysed reads only. When it counted the whole library
    this ratio was 1.7-2.2; it is now 0.95-1.04."""
    d = examples.parent / "examples_compare"
    r = lineage_sim.run_one(str(d / "KO3.fastq"), [str(d / "KO1.fastq"), str(d / "KO2.fastq")],
                            str(d / "JH4_reference.fasta"), [str(d / "JH4_contaminant.fasta")],
                            "jh4_intron", lineage_sim.SCENARIOS["low"], seed=12, depth=100)
    assert r["est_read_weighted_x1e3"] / r["true_read_weighted_x1e3"] == pytest.approx(1, abs=0.1)
    assert r["est_combined_density_x1e5"] / r["true_unique_density_x1e5"] == pytest.approx(1, abs=0.15)
