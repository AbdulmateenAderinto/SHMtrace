"""Kinds of control, the baseline fallback, quality-less controls and germline VCFs, through run_project."""
from pathlib import Path

import pytest

from shmtrace import project

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples_compare"
pytestmark = pytest.mark.skipif(not (EX / "KO3.fastq").exists(), reason="example libraries not present")


def _run(samples, **extra):
    return project.run_project({
        "name": "t", "region": "jh4_intron", "reference": str(EX / "JH4_reference.fasta"),
        "contaminants": [str(EX / "JH4_contaminant.fasta")], "samples": samples, **extra})


def _wt3():
    return {"name": "WT3", "path": str(EX / "WT3.fastq"), "group": "WT"}


def _ko(k, **kw):
    return {"name": k, "path": str(EX / f"{k}.fastq"), "group": "KO", "control": True, **kw}


def test_a_baseline_is_a_sample_when_a_mutation_free_control_exists():
    res = _run([_wt3(), _ko("KO1", control_kind="aid_null"), _ko("KO2", control_kind="aid_null"),
                _ko("KO3", control_kind="baseline")])
    s = {x["sample"]: x for x in res["samples"]}
    assert s["KO1"]["control"] and s["KO2"]["control"] and not s["KO3"]["control"]
    # sharing the controls' group, it gets its own, so the KO group stays all controls
    assert s["KO3"]["group"] == "KO (baseline)"
    assert any(g["group"] == "KO" and g["control"] for g in res["groups"])
    assert any("KO3 is a baseline" in n and "group of its own" in n for n in res["notes"])


def test_baselines_are_the_controls_when_nothing_else_is_given():
    res = _run([_wt3()] + [_ko(k, control_kind="baseline") for k in ("KO1", "KO2", "KO3")])
    assert all(x["control"] for x in res["samples"] if x["sample"].startswith("KO"))
    assert any(n.startswith("The controls are baselines") for n in res["notes"])


def test_controls_without_qualities_make_every_sample_quality_free(tmp_path):
    fa = []
    for k in ("KO1", "KO2", "KO3"):
        out = tmp_path / f"{k}.fasta"
        with open(EX / f"{k}.fastq") as src, open(out, "w") as dst:
            for i, line in enumerate(src):
                if i % 4 == 0:
                    dst.write(">" + line[1:])
                elif i % 4 == 1:
                    dst.write(line)
        fa.append({"name": k, "path": str(out), "group": "KO", "control": True})
    with_q = _run([_wt3()] + [_ko(k) for k in ("KO1", "KO2", "KO3")])
    no_q = _run([_wt3()] + fa)
    assert any("no base qualities" in n for n in no_q["notes"])
    assert not any("no base qualities" in n for n in with_q["notes"])
    calls = lambda r: next(x["unique_mutations"] for x in r["samples"] if x["sample"] == "WT3")   # noqa: E731
    # conservative, never more calls than with the controls' real qualities
    assert calls(no_q) <= calls(with_q)


def test_a_germline_vcf_removes_the_sites_it_lists(tmp_path):
    base = _run([_wt3()] + [_ko(k) for k in ("KO1", "KO2", "KO3")])
    wt3 = [v for v in base["variants"] if v["sample"] == "WT3"][:2]
    contig = (EX / "JH4_reference.fasta").read_text().split("\n", 1)[0][1:].split()[0]
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWT3\n"
    rec = lambda chrom: "".join(f"{chrom}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t30\tPASS\t.\tGT\t0/1\n"  # noqa: E731
                                for v in wt3)
    # a VCF on another contig name removes nothing and says so
    wrong = tmp_path / "wrong.vcf"
    wrong.write_text(head + rec("chr12"))
    res = _run([{**_wt3(), "germline_vcf": str(wrong)}] + [_ko(k) for k in ("KO1", "KO2", "KO3")])
    assert any("wrong.vcf" in n and "same FASTA" in n for n in res["notes"])
    vcf = tmp_path / "wt3.vcf"
    vcf.write_text(head + rec(contig))
    res = _run([{**_wt3(), "germline_vcf": str(vcf)}] + [_ko(k) for k in ("KO1", "KO2", "KO3")])
    before = next(x["unique_mutations"] for x in base["samples"] if x["sample"] == "WT3")
    after = next(x for x in res["samples"] if x["sample"] == "WT3")
    assert after["unique_mutations"] == before - 2
    assert "germline VCF wt3.vcf" in after["warnings"]
    gone = {(v["pos"], v["alt"]) for v in wt3}
    assert not gone & {(v["pos"], v["alt"]) for v in res["variants"] if v["sample"] == "WT3"}
