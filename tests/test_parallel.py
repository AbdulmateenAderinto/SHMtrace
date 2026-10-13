"""Per-library work on several processes gives the same results, in order."""
from dataclasses import replace
from pathlib import Path

import pytest

from shmtrace import batch, project
from shmtrace.parallel import pmap, resolve

EX = Path(__file__).resolve().parents[1] / "examples_compare"


def test_pmap_keeps_order_and_falls_back_to_one_process():
    args = [(i, i + 1) for i in range(6)]
    assert pmap(pow, args, 1) == pmap(pow, args, 3) == [i ** (i + 1) for i in range(6)]
    assert resolve(0) >= 1 and resolve(2) == 2


@pytest.mark.skipif(not (EX / "KO1.fastq").exists(), reason="example libraries not present")
def test_read_preparation_does_not_depend_on_the_number_of_processes():
    p = project.build_params({"region": "jh4_intron", "samples": []})
    samples = {k: str(EX / f"{k}.fastq") for k in ("WT3", "KO1", "KO2")}
    out = []
    for jobs in (1, 3):
        b = batch.prepare_batch(samples, str(EX / "JH4_reference.fasta"), [str(EX / "JH4_contaminant.fasta")],
                                replace(p, jobs=jobs), controls={"KO1", "KO2"})
        out.append({k: [(a.qname, a.pairs) for a in v] for k, v in b["reads"].items()})
    assert out[0] == out[1]
