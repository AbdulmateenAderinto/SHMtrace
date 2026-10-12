"""Sanger, Illumina and PacBio input and calling, on synthetic data with known truth."""
import gzip
import random
import struct
from dataclasses import replace

import pytest

from shmtrace import pipeline
from shmtrace.config import Params
from shmtrace.io_utils import (
    FASTA_DEFAULT_Q,
    find_pairs,
    merge_pair,
    merge_pairs,
    mott_trim,
    parse_abif,
    read_fastq,
)

pytest.importorskip("mappy", reason="alignment backend not installed")

_COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(_COMP)[::-1]


def qchars(qs):
    return "".join(chr(q + 33) for q in qs)


def random_ref(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, sites):
    s = list(seq)
    for pos0, alt in sites.items():
        s[pos0] = alt
    return "".join(s)


def alt_of(base, k=0):
    return [b for b in "ACGT" if b != base][k]


def write_abif(path, seq, quals, sample="clone"):
    """Minimal ABIF file: PBAS2, PCON2 and SMPL1 tags behind a 128-byte header."""
    entries = [(b"PBAS", 2, 2, 1, seq.encode()),
               (b"PCON", 2, 2, 1, bytes(quals)),
               (b"SMPL", 1, 18, 1, bytes([len(sample)]) + sample.encode())]
    body, offsets = b"", []
    for *_x, data in entries:
        offsets.append(128 + len(body))
        body += data
    dir_off = 128 + len(body)
    directory = b""
    for (name, num, etype, esize, data), off in zip(entries, offsets, strict=True):
        n = len(data) // esize
        field = data.ljust(4, b"\0") if len(data) <= 4 else struct.pack(">i", off)
        directory += struct.pack(">4sihhii", name, num, etype, esize, n, len(data)) + field \
            + struct.pack(">i", 0)
    header = struct.pack(">4sH4sihhiiii", b"ABIF", 101, b"tdir", 1, 1023, 28,
                         len(entries), len(directory), dir_off, 0)
    path.write_bytes(header.ljust(128, b"\0") + body + directory)


# ------------------------------------------------------------------ input --

def test_fasta_reads_get_the_documented_constant_quality(tmp_path):
    f = tmp_path / "clones.fa"
    f.write_text(">c1 desc\nACGT\nAC\n>c2\nTTTT\n")
    recs = list(read_fastq(str(f)))
    assert [(n, s) for n, s, _q in recs] == [("c1", "ACGTAC"), ("c2", "TTTT")]
    assert recs[0][2] == chr(FASTA_DEFAULT_Q + 33) * 6


def test_gzipped_fasta_and_seq(tmp_path):
    g = tmp_path / "r.fasta.gz"
    with gzip.open(g, "wt") as fh:
        fh.write(">x\nACGTN\n")
    assert list(read_fastq(str(g)))[0][:2] == ("x", "ACGTN")
    s = tmp_path / "clone_07.seq"
    s.write_text("acgt acgt\nTTGA\n")
    assert list(read_fastq(str(s)))[0][:2] == ("clone_07", "ACGTACGTTTGA")


def test_abif_parser_reads_basecalls_and_qualities(tmp_path):
    f = tmp_path / "A01_clone1.ab1"
    write_abif(f, "ACGTTGCA", [10, 20, 30, 40, 40, 30, 20, 10])
    tags = parse_abif(f.read_bytes())
    assert tags[("PBAS", 2)] == b"ACGTTGCA"
    (name, seq, qual), = list(read_fastq(str(f)))
    assert (name, seq) == ("A01_clone1", "ACGTTGCA")
    assert [ord(c) - 33 for c in qual] == [10, 20, 30, 40, 40, 30, 20, 10]


def test_abif_short_tags_are_stored_inline(tmp_path):
    f = tmp_path / "tiny.abi"
    write_abif(f, "ACG", [30, 31, 32])               # <= 4 bytes: inside the entry
    assert list(read_fastq(str(f))) == [("tiny", "ACG", qchars([30, 31, 32]))]


def test_directory_is_one_sample_and_skips_non_read_files(tmp_path):
    d = tmp_path / "clones"
    d.mkdir()
    write_abif(d / "c1.ab1", "AAAACCCC", [40] * 8)
    (d / "c2.seq").write_text("GGGGTTTT\n")
    (d / "c3.fastq").write_text("@c3\nACAC\n+\nIIII\n")
    (d / "notes.txt").write_text("sent to the core on Monday\n")
    (d / ".DS_Store").write_bytes(b"\0\1junk")
    names = sorted(n for n, _s, _q in read_fastq(str(d)))
    assert names == ["c1", "c2", "c3"]


def test_mott_trim_keeps_the_high_quality_core():
    q = qchars([5] * 10 + [40] * 50 + [8] * 15)
    assert mott_trim(q, 0.01) == (10, 60)
    assert mott_trim(qchars([3] * 10), 0.01) == (0, 0)


# ------------------------------------------------------------- read pairs --

def _pair(frag, l1, l2, q=35):
    r1 = frag[:l1]
    r2 = rc(frag)[:l2]
    return r1, qchars([q] * len(r1)), r2, qchars([q] * len(r2))


def test_merge_pair_reconstructs_the_fragment():
    frag = random_ref(400, 3)
    r1, q1, r2, q2 = _pair(frag, 250, 250)
    seq, qual = merge_pair(r1, q1, r2, q2)
    assert seq == frag
    assert len(qual) == len(seq)
    # overlap bases carry the summed (capped) quality
    assert ord(qual[200]) - 33 == 60 and ord(qual[10]) - 33 == 35


def test_merge_pair_resolves_a_disagreement_by_quality():
    frag = random_ref(300, 4)
    r1, q1, r2, q2 = _pair(frag, 200, 200)
    bad = 150                                       # inside the 100-bp overlap
    r1 = r1[:bad] + alt_of(r1[bad]) + r1[bad + 1:]
    q1 = q1[:bad] + qchars([12]) + q1[bad + 1:]
    seq, qual = merge_pair(r1, q1, r2, q2)
    assert seq == frag
    assert ord(qual[bad]) - 33 == 35 - 12


def test_merge_pair_trims_read_through_and_rejects_non_overlapping():
    frag = random_ref(120, 5)
    adapter = "AGATCGGAAGAGC" * 3
    r1 = (frag + adapter)[:150]
    r2 = (rc(frag) + adapter)[:150]
    seq, _q = merge_pair(r1, qchars([35] * 150), r2, qchars([35] * 150))
    assert seq == frag
    a, b = random_ref(150, 6), random_ref(150, 7)
    assert merge_pair(a, qchars([35] * 150), b, qchars([35] * 150)) is None


def test_merge_pairs_files_and_pair_spec(tmp_path):
    frags = [random_ref(380, 10 + i) for i in range(5)]
    r1f, r2f = tmp_path / "s_R1_001.fastq.gz", tmp_path / "s_R2_001.fastq.gz"
    with gzip.open(r1f, "wt") as f1, gzip.open(r2f, "wt") as f2:
        for i, fr in enumerate(frags):
            r1, q1, r2, q2 = _pair(fr, 250, 250)
            f1.write(f"@m{i}/1\n{r1}\n+\n{q1}\n")
            f2.write(f"@m{i}/2\n{r2}\n+\n{q2}\n")
    stats = {}
    merged = list(merge_pairs(str(r1f), str(r2f), stats=stats))
    assert [s for _n, s, _q in merged] == frags
    assert merged[0][0] == "m0"
    assert stats == {"pairs_in": 5, "pairs_merged": 5, "pairs_unmerged": 0}
    via_spec = list(read_fastq(f"{r1f}+{r2f}"))
    assert via_spec == merged
    via_dir = list(read_fastq(str(tmp_path)))          # _R1/_R2 detected in a folder
    assert via_dir == merged


def test_find_pairs_naming_conventions():
    pairs, singles = find_pairs(["a/x_S1_L001_R1_001.fastq.gz", "a/x_S1_L001_R2_001.fastq.gz",
                                 "a/y_1.fq", "a/y_2.fq", "a/z.fastq"])
    assert pairs == [("a/x_S1_L001_R1_001.fastq.gz", "a/x_S1_L001_R2_001.fastq.gz"),
                     ("a/y_1.fq", "a/y_2.fq")]
    assert singles == ["a/z.fastq"]


# --------------------------------------------------------------- platforms --

def test_platform_overrides_contract():
    jh4 = pipeline.build_params("jh4_intron", platform="ont")
    assert jh4.platform == "ont" and jh4.minimap2_preset == "map-ont"
    assert jh4.min_align_len == 800 and jh4.error_model == "qcal"
    pb = pipeline.build_params("jh4_intron", platform="pacbio")
    assert pb.minimap2_preset == "map-hifi" and pb.min_align_len == 800
    il = pipeline.build_params("jh4_intron", platform="illumina")
    assert il.minimap2_preset == "sr" and il.max_ee_rate == 0.02 and il.min_align_len == 100
    sg = pipeline.build_params("jh4_intron", platform="sanger")
    assert (sg.error_model, sg.downsample_to, sg.phase_mode, sg.quality_trim,
            sg.min_base_quality) == ("sanger", None, "off", "mott", 20)
    assert sg.analysis_window == (101, 629)          # the assay window survives
    for plat in pipeline.PLATFORMS:
        assert pipeline.platform_overrides(plat)["platform"] == plat
    assert pipeline.platform_overrides("Nanopore")["platform"] == "ont"
    with pytest.raises(ValueError):
        pipeline.platform_overrides("454")
    assert pipeline.mode_overrides("recommended", "sanger")["error_model"] == "sanger"
    assert pipeline.mode_overrides("sanger")["error_model"] == "sanger"


SANGER_LO_Q, SANGER_HI_Q = 8, 40


def _sanger_clone(ref, sites, low_q_at=()):
    seq = mutate(ref, sites)
    quals = [SANGER_LO_Q] * 25 + [SANGER_HI_Q] * (len(ref) - 55) + [SANGER_LO_Q] * 30
    for pos0 in low_q_at:
        quals[pos0] = 12
    return seq, quals


def _sanger_params(**kw):
    base = Params(name="synthetic", analysis_window=None)
    p = replace(base, **pipeline.platform_overrides("sanger", base))
    return replace(p, **kw)


@pytest.fixture(scope="module")
def sanger_case(tmp_path_factory):
    d = tmp_path_factory.mktemp("sanger")
    ref = random_ref(700, 42)
    (d / "ref.fasta").write_text(f">amplicon\n{ref}\n")
    clones_dir = d / "clones"
    clones_dir.mkdir()
    germ = {300: alt_of(ref[300])}                   # carried by every clone
    planted = {
        "c1": {100: alt_of(ref[100]), 250: alt_of(ref[250], 1)},
        "c2": {100: alt_of(ref[100])},
        "c3": {400: alt_of(ref[400], 2), 401: alt_of(ref[401]), 550: alt_of(ref[550])},
        "c4": {},
        "c5": {200: alt_of(ref[200])},
        "c6": {610: alt_of(ref[610], 1)},
    }
    low_q = {"c4": {150: alt_of(ref[150])}}          # a Q12 mismatch: not a call
    for i, (name, sites) in enumerate(planted.items()):
        extra = low_q.get(name, {})
        seq, quals = _sanger_clone(ref, {**sites, **germ, **extra}, low_q_at=extra)
        if i % 2:
            seq, quals = rc(seq), quals[::-1]
        write_abif(clones_dir / f"{name}.ab1", seq, quals)
    return d, ref, planted, germ


def test_sanger_calls_every_quality_substitution_per_clone(sanger_case):
    d, ref, planted, germ = sanger_case
    p = _sanger_params()
    res = pipeline.run(str(d / "clones"), str(d / "ref.fasta"), [], p)
    called = {(v.pos, v.alt) for v in res["variants"]}
    expected = {(pos0 + 1, alt) for sites in planted.values() for pos0, alt in sites.items()}
    assert called == expected
    s = res["summary"]
    assert s["calling"]["error_model"] == "sanger_threshold"
    assert s["calling"]["germline_positions_removed"] == 1      # the shared one
    assert s["reads"]["aligned_to_target"] == 6
    assert s["reads"]["reads_trimmed"] == 6
    per_clone = {r["read"]: r["mutations"] for r in res["per_read"]}
    assert per_clone == {name: len(sites) for name, sites in planted.items()}
    assert res["per_read_mutations"] == [r["mutations"] for r in res["per_read"]]
    ml = s["mutation_load"]
    assert ml["reads"] == 6 and ml["mutated_reads"] == 5
    assert ml["mutations_total"] == 8
    # trimmed core (positions 26..670) minus the Q12 base in c4
    assert ml["bases_sequenced"] == 6 * (700 - 55) - 1
    assert ml["mutation_frequency_per_bp"] == pytest.approx(8 / (6 * 645 - 1))
    v100 = next(v for v in res["variants"] if v.pos == 101)
    assert v100.supporting_reads == 2


def test_sanger_germline_rule_is_off_with_few_clones(sanger_case, tmp_path):
    d, ref, planted, germ = sanger_case
    few = tmp_path / "few"
    few.mkdir()
    for name in ("c1", "c2"):
        (few / f"{name}.ab1").write_bytes((d / "clones" / f"{name}.ab1").read_bytes())
    res = pipeline.run(str(few), str(d / "ref.fasta"), [], _sanger_params())
    called = {(v.pos, v.alt) for v in res["variants"]}
    assert (301, germ[300]) in called                # cannot be told from germline
    assert any("germline filter disabled" in n for n in res["summary"]["calling"]["notes"])


def test_sanger_fasta_clones_via_mode(sanger_case, tmp_path):
    """FASTA clones (constant Q30) through mode 'sanger' on a default Params."""
    d, ref, planted, germ = sanger_case
    fa = tmp_path / "clones.fasta"
    fa.write_text("".join(f">{n}\n{mutate(ref, {**s, **germ})[25:670]}\n"
                          for n, s in planted.items()))
    p = replace(Params(name="x", analysis_window=None, min_align_len=200,
                       downsample_to=None), **pipeline.mode_overrides("sanger"))
    res = pipeline.run(str(fa), str(d / "ref.fasta"), [], p)
    expected = {(pos0 + 1, alt) for sites in planted.values() for pos0, alt in sites.items()}
    assert {(v.pos, v.alt) for v in res["variants"]} == expected


# ------------------------------------------------- Illumina-like / PacBio --

def _noisy(seq, quals, rng):
    out = []
    for b, q in zip(seq, quals, strict=True):
        if rng.random() < 10 ** (-q / 10):
            b = rng.choice([x for x in "ACGT" if x != b])
        out.append(b)
    return "".join(out)


@pytest.fixture(scope="module")
def amplicon(tmp_path_factory):
    d = tmp_path_factory.mktemp("amplicon")
    ref = random_ref(450, 7)
    (d / "ref.fasta").write_text(f">amp\n{ref}\n")
    planted = {60: 0.10, 120: 0.08, 200: 0.15, 260: 0.06, 330: 0.12, 400: 0.10}
    alts = {pos: alt_of(ref[pos]) for pos in planted}
    return d, ref, planted, alts


def _molecules(ref, planted, alts, n, rng):
    for _ in range(n):
        yield mutate(ref, {pos: alts[pos] for pos, f in planted.items() if rng.random() < f})


def test_illumina_pairs_through_pipeline_run(amplicon):
    d, ref, planted, alts = amplicon
    rng = random.Random(11)
    q = [36] * 220 + [18] * 30                       # Illumina-like decaying tail
    r1f, r2f = d / "lib_R1.fastq.gz", d / "lib_R2.fastq.gz"
    with gzip.open(r1f, "wt") as f1, gzip.open(r2f, "wt") as f2:
        for i, mol in enumerate(_molecules(ref, planted, alts, 400, rng)):
            r1 = _noisy(mol[:250], q, rng)
            r2 = _noisy(rc(mol)[:250], q, rng)
            f1.write(f"@m{i} 1:N:0\n{r1}\n+\n{qchars(q)}\n")
            f2.write(f"@m{i} 2:N:0\n{r2}\n+\n{qchars(q)}\n")
    base = Params(name="amp", analysis_window=None, downsample_to=None)
    expected = {(pos + 1, alts[pos]) for pos in planted}
    for mode in ("recommended", "pbinom"):
        p = replace(base, **pipeline.platform_overrides("illumina", base),
                    **pipeline.mode_overrides(mode, "illumina"))
        res = pipeline.run(f"{r1f}+{r2f}", str(d / "ref.fasta"), [], p)
        called = {(v.pos, v.alt) for v in res["variants"]}
        assert res["summary"]["reads"]["aligned_to_target"] >= 390
        assert expected <= called, f"{mode}: missed {expected - called}"
        assert len(called - expected) <= 1, f"{mode}: false calls {called - expected}"
        assert len(res["per_read_mutations"]) == res["summary"]["reads"]["aligned_to_target"]


def test_pacbio_hifi_like_reads(amplicon):
    d, ref, planted, alts = amplicon
    rng = random.Random(12)
    fq = d / "hifi.fastq"
    with open(fq, "w") as fh:
        for i, mol in enumerate(_molecules(ref, planted, alts, 300, rng)):
            q = [rng.choice([30, 35, 40, 45]) for _ in mol]
            s = _noisy(mol, q, rng)
            if i % 2:
                s, q = rc(s), q[::-1]
            fh.write(f"@ccs{i}\n{s}\n+\n{qchars(q)}\n")
    base = Params(name="amp", analysis_window=None, min_align_len=400, max_total_len=600,
                  downsample_to=None)
    p = replace(base, **pipeline.platform_overrides("pacbio", base),
                **pipeline.mode_overrides("recommended", "pacbio"))
    res = pipeline.run(str(fq), str(d / "ref.fasta"), [], p)
    called = {(v.pos, v.alt) for v in res["variants"]}
    expected = {(pos + 1, alts[pos]) for pos in planted}
    assert expected <= called and len(called - expected) <= 1


def test_a_one_clone_sanger_sample_borrows_the_projects_germline(sanger_case, tmp_path):
    # one clone cannot tell its germline differences from mutations; the other
    # mice in the project can. On the real J꜀4 clones a single Aicda-/- clone
    # reported the strain's five germline differences as five mutations.
    from shmtrace import project
    d, ref, planted, germ = sanger_case
    many = tmp_path / "many.fasta"
    many.write_text("".join(f">{n}\n{mutate(ref, {**s, **germ})[25:670]}\n"
                            for n, s in planted.items()))
    one = tmp_path / "one.fasta"
    one.write_text(f">solo\n{mutate(ref, {**germ, 500: alt_of(ref[500])})[25:670]}\n")
    proj = {"name": "g", "region": "generic_amplicon", "platform": "sanger",
            "reference": str(d / "ref.fasta"), "contaminants": [],
            "overrides": {"min_align_len": 200},
            "samples": [{"name": "many", "path": str(many), "group": "A"},
                        {"name": "one", "path": str(one), "group": "B"}]}
    res = project.run_project(proj)
    called = {(v.pos, v.alt) for v in res["_details"]["one"]["variants"]}
    assert called == {(501, alt_of(ref[500]))}
    assert (301, germ[300]) not in called
    row = next(s for s in res["samples"] if s["sample"] == "one")
    assert "germline difference(s) removed using the project's other clones" in row["warnings"]
    # the well-sampled mouse is untouched
    many_called = {(v.pos, v.alt) for v in res["_details"]["many"]["variants"]}
    assert many_called == {(p0 + 1, a) for s in planted.values() for p0, a in s.items()}
