"""Fitting an empirical error model, and refusing to use it across runs."""
from shmtrace.null_model import check_same_run, run_token


def _fastq(tmp_path, name, read_id):
    f = tmp_path / name
    f.write_text(f"@{read_id}\nACGT\n+\nIIII\n")
    return str(f)


def test_run_token_parses_the_provider_read_id(tmp_path):
    assert run_token(_fastq(tmp_path, "a.fastq", "0_DHJGKB_5")) == "DHJGKB"
    assert run_token(_fastq(tmp_path, "b.fastq", "17_ZT44NG_10")) == "ZT44NG"


def test_run_token_returns_none_on_unfamiliar_ids(tmp_path):
    assert run_token(_fastq(tmp_path, "c.fastq", "read1")) is None


def test_same_run_is_accepted(tmp_path):
    ctrl = _fastq(tmp_path, "ctrl.fastq", "0_RUNA_1")
    tgt = _fastq(tmp_path, "tgt.fastq", "0_RUNA_2")
    ok, why = check_same_run([ctrl], tgt)
    assert ok
    assert "RUNA" in why


def test_cross_run_is_flagged(tmp_path):
    """Nanopore error is run-specific; a null from another flow cell is not valid."""
    ctrl = _fastq(tmp_path, "ctrl.fastq", "0_RUNA_1")
    tgt = _fastq(tmp_path, "tgt.fastq", "0_RUNB_2")
    ok, why = check_same_run([ctrl], tgt)
    assert not ok
    assert "RUNA" in why and "RUNB" in why


def test_mixed_control_runs_are_flagged(tmp_path):
    a = _fastq(tmp_path, "a.fastq", "0_RUNA_1")
    b = _fastq(tmp_path, "b.fastq", "0_RUNB_1")
    tgt = _fastq(tmp_path, "t.fastq", "0_RUNA_2")
    ok, _why = check_same_run([a, b], tgt)
    assert not ok


def test_unrecognised_ids_do_not_block(tmp_path):
    """Other providers name reads differently; warn-only, never a hard failure."""
    ctrl = _fastq(tmp_path, "ctrl.fastq", "read1")
    tgt = _fastq(tmp_path, "tgt.fastq", "read2")
    ok, why = check_same_run([ctrl], tgt)
    assert ok
    assert "cannot verify" in why


def test_run_token_reads_the_minknow_runid_tag(tmp_path):
    f = tmp_path / "m.fastq"
    f.write_text("@3f2a-read runid=abc123def ch=5 start_time=x\nACGT\n+\nIIII\n")
    assert run_token(str(f)) == "abc123def"


def test_explicit_run_ids_override_read_names(tmp_path):
    ctrl = _fastq(tmp_path, "c.fastq", "0_AAAAAA_1")
    tgt = _fastq(tmp_path, "t.fastq", "0_BBBBBB_2")
    ok, _ = check_same_run([ctrl], tgt, sample_run="R1", control_runs=["R1"])
    assert ok
    ok, _ = check_same_run([ctrl], tgt, sample_run="R1", control_runs=["R2"])
    assert not ok
