"""AID targeting against a 5-mer model (targeting.py), on a small synthetic model."""
import itertools
import random

import pytest

from shmtrace import project
from shmtrace import targeting as tg


@pytest.fixture
def model_csv(tmp_path):
    """Every 5-mer, with WRCY / RGYW hotspots 20x as mutable as the rest."""
    rows = ["fivemer,mutability,to_A,to_C,to_G,to_T"]
    for k in ("".join(t) for t in itertools.product("ACGT", repeat=5)):
        hot = tg.matches(k, "WRCY") or tg.matches(k, "RGYW")
        m = 0.02 if hot else 0.001
        probs = {b: (0.0 if b == k[2] else 1 / 3) for b in "ACGT"}
        if k[2] == "C":
            probs = {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}      # C>T favoured
        rows.append(f"{k},{m}," + ",".join(str(probs[b]) for b in "ACGT"))
    rows.append("NNNNN,0.5,NA,NA,NA,NA")                          # skipped
    f = tmp_path / "model.csv"
    f.write_text("\n".join(rows) + "\n")
    return str(f)


def test_load_skips_ambiguous_fivemers(model_csv):
    m = tg.load(model_csv)
    assert len(m.mutability) == 1024
    assert "NNNNN" not in m.mutability


def test_motif_centres():
    # WRCY: the C is the 5-mer's centre, motif in positions 0-3
    assert tg.matches("AGCTA", "WRCY")
    assert not tg.matches("CGCTA", "WRCY")     # C is not W
    assert not tg.matches("AGCGA", "WRCY")     # G is not Y
    # RGYW: the G is the centre, motif in positions 1-4
    assert tg.matches("AAGCT", "RGYW")
    assert not tg.matches("ACGCT", "RGYW")     # C is not R
    # coldspots
    assert tg.matches("GCCAA", "SYC")
    assert not tg.matches("ACCAA", "SYC")      # A is not S
    assert tg.matches("AAGAC", "GRS")


def test_mutations_on_hotspots_score_above_one(model_csv):
    m = tg.load(model_csv)
    rng = random.Random(1)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    prof = tg.window_profile(m, ref, 10, 590)
    hot = [i for i, v in prof["mutability"].items() if v > prof["mean"]]
    on_hot = tg.targeting_index(rng.sample(hot, 30), prof)
    assert on_hot["index"] > 2 and on_hot["p"] < 0.01
    anywhere = tg.targeting_index(rng.sample(list(prof["mutability"]), 30), prof, seed=3)
    assert 0.3 < anywhere["index"] < 2.5 and anywhere["p"] > 0.01


def test_index_is_withheld_below_the_minimum(model_csv):
    m = tg.load(model_csv)
    ref = "AGCTA" * 60
    prof = tg.window_profile(m, ref, 5, 290)
    few = tg.targeting_index([7, 12, 17], prof)
    assert few["index"] is None and few["n"] == 3 and "fewer than" in few["note"]


def test_expected_spectrum_follows_the_substitution_model(model_csv):
    m = tg.load(model_csv)
    ref = "AGCTA" * 60
    prof = tg.window_profile(m, ref, 5, 290)
    exp = tg.expected_spectrum(m, ref, prof)
    assert abs(sum(exp.values()) - 1) < 1e-4      # six-decimal rounding
    assert exp["C>T"] > exp["C>A"] and exp["C>T"] > exp["C>G"]


def test_project_reports_targeting(examples, model_csv):
    res = project.run_project({
        "name": "t", "region": "jh4_intron",
        "reference": str(examples / "JH4_reference.fasta"),
        "contaminants": [str(examples / "JH4_contaminant.fasta")],
        "targeting_model": model_csv,
        "samples": [{"name": "WT", "path": str(examples / "WT3_JH4_sub500.fastq"), "group": "WT"},
                    {"name": "KO", "path": str(examples / "KO3_JH4_sub500.fastq"),
                     "group": "KO", "control": True}]})
    t = res["targeting"]
    assert t["model"] == "model" and t["groups"]["WT"]["n"] > 10
    wt = next(s for s in res["samples"] if s["sample"] == "WT")
    assert wt["s5f_index"] is not None
    assert all("s5f_mutability_rel" in v for v in res["variants"])
