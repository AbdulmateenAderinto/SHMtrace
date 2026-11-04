"""Short reads kept when they align well (Params.rescue_short_reads), off by default."""
import pickle
import random
from dataclasses import fields
from types import SimpleNamespace

import pytest

from shmtrace import browse, pipeline
from shmtrace.align import (
    AlignedRead,
    diagnostic_sites,
    drop_unspecific_short_reads,
    record_sites,
    short_read_flaw,
)
from shmtrace.config import Params, available_presets, load_preset
from shmtrace.io_utils import read_fasta
from shmtrace.project import build_params

REFS = pipeline.__file__.rsplit("/", 1)[0] + "/references/"


def _fastq(path, reads):
    path.write_text("".join(f"@{n}\n{s}\n+\n{chr(33 + 40) * len(s)}\n" for n, s in reads))
    return str(path)


def _mutate(seq, rate, seed=2):
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    rng = random.Random(seed)
    return "".join(swap[b] if rng.random() < rate else b for b in seq)


def test_off_in_every_preset_and_platform():
    assert Params().rescue_short_reads is False
    for name in available_presets():
        assert load_preset(name).rescue_short_reads is False
    for plat in ("ont", "pacbio", "illumina", "sanger"):
        assert pipeline.build_params("jh4_intron", platform=plat).rescue_short_reads is False


@pytest.fixture(scope="module")
def jh4(examples):
    ref = str(examples / "JH4_reference.fasta")
    return ref, read_fasta(ref)[0][1]


def test_a_clean_short_read_is_kept_only_when_asked(tmp_path, jh4):
    ref, seq = jh4
    fq = _fastq(tmp_path / "r.fastq", [("full", seq[40:900]), ("short", seq[100:600])])
    off: dict = {}
    kept, st = pipeline.prepare_reads(fq, ref, [], build_params({"region": "jh4_intron"}), fates=off)
    assert [a.qname for a in kept] == ["full"]
    assert off["short"]["fate"] == "alignment too short" and "reason" not in off["short"]
    assert "short_reads_kept" not in st                  # the off path reports as it always did

    on: dict = {}
    p = build_params({"region": "jh4_intron", "overrides": {"rescue_short_reads": True}})
    kept, st = pipeline.prepare_reads(fq, ref, [], p, fates=on)
    assert sorted(a.qname for a in kept) == ["full", "short"]
    assert [a.qname for a in kept if getattr(a, "short", False)] == ["short"]
    assert st["short_reads_kept"] == 1 and "short" not in on


def test_each_test_names_why_a_short_read_was_not_kept(tmp_path, jh4):
    ref, seq = jh4
    junk = "".join(random.Random(1).choice("ACGT") for _ in range(400))
    fq = _fastq(tmp_path / "r.fastq", [
        ("stub", seq[300:460]),                  # spans less than half of min_align_len
        ("outside", seq[700:1200]),              # long enough, but not in the window 101-629
        ("junk", seq[100:600] + junk),           # a long unaligned stretch
        ("divergent", _mutate(seq[100:600], 0.08)),  # 8% substitutions: a paralogue, not SHM
    ])
    p = build_params({"region": "jh4_intron", "overrides": {"rescue_short_reads": True}})
    fates: dict = {}
    kept, _st = pipeline.prepare_reads(fq, ref, [], p, fates=fates)
    assert kept == []
    assert all(fates[n]["fate"] == "alignment too short" for n in fates)
    assert fates["stub"]["reason"].startswith("aligned span 160 bp")
    assert fates["outside"]["reason"].startswith("covers 0 bp of the analysis window")
    assert "of the read aligned" in fates["junk"]["reason"]
    assert "divergent" in fates["divergent"]["reason"]


def test_flaw_checks_in_order():
    p = Params(rescue_short_reads=True, min_align_len=800, analysis_window=(101, 629))
    hit = dict(r_st=100, r_en=600, q_st=0, q_en=520, NM=5, blen=500, ctg_len=1230)
    assert short_read_flaw(SimpleNamespace(**hit), 520, 1, p) is None
    assert "primary" in short_read_flaw(SimpleNamespace(**hit), 520, 2, p)
    assert "divergent" in short_read_flaw(SimpleNamespace(**{**hit, "NM": 30}), 520, 1, p)
    assert "aligned" in short_read_flaw(SimpleNamespace(**hit), 800, 1, p)
    # window 101-629 is 529 bp: half of it (264.5 bp) must be covered
    assert "window" in short_read_flaw(SimpleNamespace(**{**hit, "r_st": 370, "r_en": 790}), 520, 1, p)
    assert "span" in short_read_flaw(SimpleNamespace(**{**hit, "r_en": 450}), 520, 1, p)


def _read(name, lo, hi, ref, alt=None, short=True):
    pairs = [(i, (alt or {}).get(i, ref[i]), 40) for i in range(lo, hi)]
    r = AlignedRead(name, "t", lo, hi, pairs, "")
    if short:
        r.short = True
    return r


def test_a_short_read_must_show_it_is_target_over_each_contaminant():
    ref = "A" * 300
    # record C covers 0-200 and differs from the reference at 10 sites; D covers 200-300 at one
    records = [("C", [(0, 200)], {k: {"G"} for k in range(5, 200, 20)}),
               ("D", [(200, 300)], {250: {"T"}})]
    reads = [
        _read("target", 0, 200, ref),                                   # 10 reference bases
        _read("mixed", 0, 200, ref, {k: "G" for k in range(5, 100, 20)}),  # 5 v 5
        _read("blind", 100, 240, ref),                                  # 5 of C, none of D under it
        _read("full", 0, 300, ref, short=False),                        # not judged here
    ]
    fates: dict = {}
    kept, n = drop_unspecific_short_reads(reads, records, ref, 5, 0.5, fates)
    assert [r.qname for r in kept] == ["target", "full"] and n == 2
    assert "not told from C" in fates["mixed"]["reason"]
    assert fates["blind"]["fate"] == "alignment too short"


def test_sites_most_of_the_sample_shares_are_not_evidence():
    ref = "A" * 200
    records = [("C", [(0, 200)], {k: {"G"} for k in range(5, 200, 20)})]
    shared = {5: "G", 25: "G", 45: "G"}                 # germline the sample shares with C
    reads = [_read(f"r{i}", 0, 200, ref, shared, short=False) for i in range(5)]
    reads.append(_read("short", 0, 200, ref, shared))
    kept, n = drop_unspecific_short_reads(reads, records, ref, 7, 0.5)
    assert n == 0 and len(kept) == 6                    # 7 reference bases, 0 counted against
    kept, n = drop_unspecific_short_reads(reads, records, ref, 8, 0.5)
    assert n == 1


def test_record_sites_add_up_to_the_diagnostic_sites():
    ref, con = REFS + "mouse_jh4_intron.fasta", [REFS + "mouse_jh4_contaminant.fasta"]
    recs = record_sites(ref, con)
    assert [n for n, _s, _x in recs] == ["Contaminant_JH4", "JH4_paralog_sanger_consensus"]
    union: dict = {}
    for _n, _spans, sites in recs:
        for k, b in sites.items():
            union.setdefault(k, set()).update(b)
    assert union == {k: set(v) for k, v in diagnostic_sites(ref, con).items()}


def test_against_the_bundled_contaminants_a_clean_short_read_stays(tmp_path):
    ref = REFS + "mouse_jh4_intron.fasta"
    seq = read_fasta(ref)[0][1]
    fq = _fastq(tmp_path / "r.fastq", [("short", seq[100:600])])
    p = build_params({"region": "jh4_intron", "overrides": {"rescue_short_reads": True}})
    kept, st = pipeline.prepare_reads(fq, ref, [REFS + "mouse_jh4_contaminant.fasta"], p)
    assert [a.qname for a in kept] == ["short"] and st["short_reads_unspecific"] == 0


def test_the_flag_survives_a_worker_process():
    r = _read("s", 0, 50, "A" * 50)
    back = pickle.loads(pickle.dumps(r))
    assert back.short is True and back == r
    assert "short" not in {f.name for f in fields(AlignedRead)}


def test_the_reads_tab_names_the_fate_and_counts_it_as_analysed(tmp_path, jh4):
    ref, seq = jh4
    full, short = _read("full", 40, 900, seq.upper(), short=False), _read("short", 100, 600, seq.upper())
    data = browse.alignment([full, short], {"stub": {"fate": "alignment too short", "r_st": 300,
                                                     "r_en": 460, "reason": "aligned span 160 bp"}},
                            ref)
    fates = {r["name"]: r["fate"] for r in data["reads"]}
    assert fates == {"full": "kept", "short": "kept (short)", "stub": "alignment too short"}
    assert "kept (short)" in browse.FATE_ORDER and "kept (short)" in browse.KEPT
    assert next(r for r in data["reads"] if r["name"] == "stub")["reason"] == "aligned span 160 bp"
    only = browse.window(data, 1, len(seq), fates={"kept"})
    assert {r["name"] for r in only["reads"]} == {"full", "short"}
    assert browse.call_reads(data, 300, "G")["covering"] == 2


def test_the_frequency_note_and_the_sample_count(tmp_path, jh4):
    ref, seq = jh4
    reads = [(f"f{i}", seq[40:900]) for i in range(6)] + [("short", seq[100:600])]
    fq = _fastq(tmp_path / "r.fastq", reads)
    p = build_params({"region": "jh4_intron", "overrides": {"rescue_short_reads": True}})
    aligned, st = pipeline.prepare_reads(fq, ref, [], p)
    res = pipeline.analyse(aligned, ref, p, st)
    assert res["summary"]["reads"]["short_reads_analysed"] == 1
    assert any("rescue_short_reads" in n for n in res["summary"]["calling"]["notes"])


@pytest.mark.parametrize("platform,span", [("sanger", (150, 330)), ("illumina", (150, 230))])
def test_sanger_and_illumina_keep_their_short_reads_too(tmp_path, jh4, platform, span):
    ref, seq = jh4
    fq = _fastq(tmp_path / "r.fastq", [("short", seq[span[0]:span[1]])])
    p = build_params({"region": "jh4_intron", "platform": platform,
                      "overrides": {"rescue_short_reads": True}})
    assert span[1] - span[0] < p.min_align_len
    kept, st = pipeline.prepare_reads(fq, ref, [], p)
    assert [a.qname for a in kept] == ["short"] and getattr(kept[0], "short", False)
    assert st["short_reads_kept"] == 1
