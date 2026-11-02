"""Codons, amino acids and R/S by region (coding.py)."""
import pytest

from shmtrace import coding

#            codon 1  2   3   4
REF = "GG" + "GAT" "TGG" "TAC" "CTG" + "AA"     # coding 3-14: D W Y L


def test_substitutions_are_placed_in_their_codon_and_classed():
    cds = [3, 14]
    assert coding.classify(REF, cds, 3, "A") == {"codon_number": 1, "codon_ref": "GAT", "codon_alt": "AAT",
                                                 "aa_ref": "D", "aa_alt": "N", "effect": "replacement"}
    assert coding.classify(REF, cds, 5, "C")["effect"] == "silent"         # GAT -> GAC, D
    assert coding.classify(REF, cds, 8, "A")["effect"] == "stop"           # TGG -> TGA
    assert coding.classify(REF, cds, 1, "A") == {"effect": "non-coding"}
    assert coding.classify(REF, None, 5, "C") == {"effect": "non-coding"}


def test_annotations_are_checked_and_regions_read_from_lines():
    a = coding.parse({"cds": "3-14", "regions": "CDR1 6-8\nFR1 3-5"}, len(REF))
    assert a == {"cds": [3, 14], "regions": [{"name": "FR1", "start": 3, "end": 5}, {"name": "CDR1", "start": 6, "end": 8}]}
    assert coding.parse({}, 20) is None
    with pytest.raises(ValueError, match="within"):
        coding.parse({"cds": "3-99"}, len(REF))
    with pytest.raises(ValueError, match="name and a span"):
        coding.parse({"regions": "CDR1"}, len(REF))


def test_rs_counts_per_sample_and_region():
    res = {"reference": REF, "samples": [{"sample": "WT1", "group": "Wt"}],
           "variants": [{"sample": "WT1", "pos": 3, "ref": "G", "alt": "A"},      # FR1 replacement
                        {"sample": "WT1", "pos": 5, "ref": "T", "alt": "C"},      # FR1 silent
                        {"sample": "WT1", "pos": 8, "ref": "G", "alt": "A"},      # CDR1 stop
                        {"sample": "WT1", "pos": 1, "ref": "G", "alt": "A"}],     # outside
           "annotation": {"cds": [3, 14], "regions": [{"name": "FR1", "start": 3, "end": 5},
                                                      {"name": "CDR1", "start": 6, "end": 8}]}}
    coding.apply(res)
    rows = {r["region"]: r for r in res["rs_by_region"]}
    assert rows["all"]["replacement"] == 1 and rows["all"]["silent"] == 1 and rows["all"]["stop"] == 1
    assert rows["all"]["non_coding"] == 1 and rows["all"]["rs_ratio"] == 1.0
    assert rows["FR1"]["replacement"] == 1 and rows["FR1"]["silent"] == 1
    assert rows["CDR1"]["stop"] == 1 and rows["CDR1"]["rs_ratio"] is None
    assert res["coding"]["WT1|8|G>A"]["region"] == "CDR1"
    res["annotation"] = None
    coding.apply(res)
    assert "coding" not in res and "rs_by_region" not in res
