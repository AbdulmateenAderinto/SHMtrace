"""Batch behaviour: shared haplotype evidence and a common analysis depth."""
from shmtrace import batch, pipeline
from shmtrace.config import load_preset

BLOCK = [(10, "T"), (20, "T"), (30, "T"), (40, "T")]
CLONE = [(110, "C"), (120, "C"), (130, "C"), (140, "C")]


def test_a_block_in_two_animals_is_known_but_a_private_clone_is_not():
    known = batch.recurring_blocks({"m1": [BLOCK, CLONE], "m2": [BLOCK]}, 4)
    assert known == [BLOCK]


def test_analyse_warns_when_a_sample_is_below_the_target_depth(ko_fastq, reference, contaminant):
    p = load_preset("jh4_intron")
    aligned, st = pipeline.prepare_reads(str(ko_fastq), str(reference), [str(contaminant)], p)
    res = pipeline.analyse(aligned, str(reference), p, st, depth=len(aligned) + 50)
    assert res["summary"]["warnings"], "a depth shortfall must be reported"


def test_batch_analyses_every_sample_at_one_depth(wt_fastq, ko_fastq, reference, contaminant):
    p = load_preset("jh4_intron")
    b = batch.prepare_batch({"wt": str(wt_fastq), "ko": str(ko_fastq)},
                            str(reference), [str(contaminant)], p, controls={"ko"})
    assert b["depth"] <= min(len(r) for r in b["reads"].values())


def test_rarefaction_rises_with_depth(wt_fastq, reference, contaminant):
    p = load_preset("jh4_intron")
    aligned, _ = pipeline.prepare_reads(str(wt_fastq), str(reference), [str(contaminant)], p)
    rows = pipeline.rarefaction(aligned, str(reference), p, [60, 250], draws=2)
    assert rows[0]["unique_variants_mean"] < rows[1]["unique_variants_mean"]
