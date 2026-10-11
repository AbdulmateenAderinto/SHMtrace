"""Geneious compatibility: the p-value, the export parser and the benchmark."""
from pathlib import Path

import pytest

from nanoshm import pipeline
from nanoshm.benchmark import read_geneious_export, read_sheet, site_concordance
from nanoshm.config import load_preset
from nanoshm.variants import _geneious_pvalue

ROOT = Path(__file__).resolve().parents[1]
CMP = ROOT / "examples_compare"


def test_pvalue_reproduces_the_geneious_manual_example():
    """Three Q20 reads, two carry the variant: 3C2 x 0.01^2 x 0.99 + 0.01^3."""
    assert _geneious_pvalue(2, 3, 20.0) == pytest.approx(0.000298, rel=1e-3)


def test_single_read_needs_a_high_quality_base_to_pass():
    """At 500 reads one supporting base passes p <= 0.05 only above ~Q40."""
    assert _geneious_pvalue(1, 500, 45.0) < 0.05
    assert _geneious_pvalue(1, 500, 30.0) > 0.05


def test_export_parser_reads_geneious_annotation_tables(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text(
        "Name,Minimum,Maximum,Change,Polymorphism Type,Variant Frequency\n"
        "G -> A,120,120,G -> A,SNP (transition),0.4%\n"
        "T -> C,1,000,1,000,T -> C,SNP (transition),0.2%\n".replace("1,000,1,000", "1000,1000")
        + "Deletion,300,301,-AT,Deletion,0.2%\n"
    )
    assert read_geneious_export(str(f)) == {(120, "G", "A"), (1000, "T", "C")}
    assert read_geneious_export(str(f), offset=-100) == {(20, "G", "A"), (900, "T", "C")}


def test_export_parser_reads_nanoshm_variant_tables(tmp_path):
    f = tmp_path / "v.csv"
    f.write_text("position,ref,alt,supporting_reads\n5,A,G,1\n")
    assert read_geneious_export(str(f)) == {(5, "A", "G")}


def test_site_concordance():
    c = site_concordance({(1, "A", "G"), (2, "C", "T")}, {(2, "C", "T"), (3, "G", "A")})
    assert c["shared"] == 1 and c["jaccard"] == pytest.approx(1 / 3, abs=1e-3)
    assert c["precision_vs_geneious"] == 0.5 and c["recall_vs_geneious"] == 0.5


def test_sample_sheet_resolves_paths_and_flags():
    rows = read_sheet(str(CMP / "samples.csv"))
    assert [r["sample"] for r in rows] == ["WT2", "WT3", "WT4", "KO1", "KO2", "KO3"]
    assert [r["null"] for r in rows] == [False] * 3 + [True] * 3
    assert Path(rows[0]["fastq"]).exists()
    assert rows[1]["geneious_calls"] == 342


@pytest.mark.parametrize("sample,geneious,tol", [("WT3", 342, 15), ("KO3", 3, 3)])
def test_geneious_mode_reproduces_the_manuscript_counts(sample, geneious, tol):
    """Manuscript Table S1, Geneious Prime 2026.1.2, 500 reads, 529-bp window."""
    p = load_preset("jh4_intron", **pipeline.mode_overrides("geneious"))
    res = pipeline.run(str(CMP / f"{sample}.fastq"), str(CMP / "JH4_reference.fasta"),
                       [str(CMP / "JH4_contaminant.fasta")], p)
    assert res["summary"]["calling"]["unique_variants"] == pytest.approx(geneious, abs=tol)


def test_export_parser_reapplies_manuscript_filters(tmp_path):
    """Exports list every candidate; keep p <= 0.05, frequency < 50%, in window."""
    f = tmp_path / "g.csv"
    f.write_text(
        "Minimum,Change,Polymorphism Type,Variant P-Value (approximate),Variant Frequency\n"
        "120,G -> A,SNP (transition),1.0E-3,0.4%\n"
        "130,C -> T,SNP (transition),0.2,0.2%\n"
        "140,A -> G,SNP (transition),1.0E-6,97%\n"
        "700,T -> C,SNP (transition),1.0E-4,1%\n"
    )
    assert read_geneious_export(str(f), window=(101, 629)) == {(120, "G", "A")}
    assert len(read_geneious_export(str(f), max_p=None, max_freq=None)) == 4
