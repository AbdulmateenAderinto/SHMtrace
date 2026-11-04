"""Length gates scaled to the reference (Params.min_align_frac / max_total_frac)."""
from shmtrace import pipeline
from shmtrace.config import (
    ADAPTER_ALLOWANCE_BP,
    Params,
    available_presets,
    load_preset,
    scale_lengths,
)
from shmtrace.io_utils import read_fasta
from shmtrace.project import build_params


def test_only_generic_amplicon_scales():
    for name in available_presets():
        p = load_preset(name)
        if name == "generic_amplicon":
            assert (p.min_align_frac, p.max_total_frac) == (0.65, 1.5)
        else:
            assert p.min_align_frac is None and p.max_total_frac is None
            assert scale_lengths(p, 900) is p          # nothing to scale, nothing changes


def test_scaled_to_the_reference():
    q = scale_lengths(Params(min_align_frac=0.65, max_total_frac=1.5), 1230)
    assert q.min_align_len == 800                     # J꜀4's own preset uses 800
    assert q.max_total_len == 1845 + ADAPTER_ALLOWANCE_BP
    assert q.min_align_frac is None and q.max_total_frac is None
    assert scale_lengths(q, 100) is q                 # scaling happens once


def test_illumina_and_sanger_use_their_lengths_as_given():
    for plat, lo, hi in (("illumina", 100, 20000), ("sanger", 200, 20000)):
        p = scale_lengths(pipeline.build_params("generic_amplicon", platform=plat), 900)
        assert (p.min_align_len, p.max_total_len) == (lo, hi)
    p = scale_lengths(pipeline.build_params("generic_amplicon", platform="pacbio"), 900)
    assert (p.min_align_len, p.max_total_len) == (585, 1450)


def test_a_length_set_outright_wins():
    p = scale_lengths(build_params({"region": "generic_amplicon",
                                    "overrides": {"min_align_len": 150}}), 900)
    assert (p.min_align_len, p.max_total_len) == (150, 1450)   # the other one still scales
    p = scale_lengths(build_params({"region": "generic_amplicon",
                                    "overrides": {"min_align_len": 150, "max_total_len": 400}}), 900)
    assert (p.min_align_len, p.max_total_len) == (150, 400)
    p = scale_lengths(load_preset("generic_amplicon", max_total_len=700), 900)
    assert (p.min_align_len, p.max_total_len) == (585, 700)


def test_a_run_filters_by_the_scaled_lengths(tmp_path, examples):
    ref = str(examples / "JH4_reference.fasta")           # 1,230 bp: 800 and 1,945 bp
    seq = read_fasta(ref)[0][1]
    long = seq + "ACGT" * 200                              # 2,030 bp raw
    q = chr(33 + 40)
    fq = tmp_path / "r.fastq"
    fq.write_text("".join(f"@{n}\n{s}\n+\n{q * len(s)}\n" for n, s in
                          [("full", seq[40:900]), ("part", seq[100:700]), ("long", long)]))
    p = build_params({"region": "generic_amplicon"})
    fates: dict = {}
    kept, _st = pipeline.prepare_reads(str(fq), ref, [], p, fates=fates)
    assert [a.qname for a in kept] == ["full"]
    assert fates["part"]["fate"] == "alignment too short"   # 600 bp: kept under the old 200
    assert fates["long"]["fate"] == "too long"
    res = pipeline.run(str(fq), ref, [], p)
    assert res["summary"]["params"]["min_align_len"] == 800
    assert res["summary"]["params"]["max_total_len"] == 1945
