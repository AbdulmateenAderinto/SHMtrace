"""Read pre-processing: optional quality trim, expected-error + length gate,
then depth downsampling."""
from __future__ import annotations

import random

from .config import Params
from .io_utils import expected_error_rate, mott_trim


def prefilter(reads, p: Params) -> tuple[list[tuple[str, str, str]], dict]:
    """Apply the EE-rate and raw-length filters (Step 2 of the pipeline).

    Alignment-length is enforced later against the reference; here we only need
    the read to be long enough and clean enough to be worth aligning. With
    p.quality_trim == "mott" each read is first trimmed to its best-quality
    segment (Sanger traces), and the gates apply to the trimmed read.
    """
    if p.quality_trim not in ("off", "mott"):
        raise ValueError(f"unknown quality_trim {p.quality_trim!r}; use 'off' or 'mott'")
    kept, n_in, n_ee, n_len, n_trim, n_empty = [], 0, 0, 0, 0, 0
    for qname, seq, qual in reads:
        n_in += 1
        if p.quality_trim == "mott":
            a, b = mott_trim(qual, p.trim_cutoff)
            if b - a < len(seq):
                n_trim += 1
            if b <= a:
                n_empty += 1
                continue
            seq, qual = seq[a:b], qual[a:b]
        if len(seq) > p.max_total_len:
            n_len += 1
            continue
        if expected_error_rate(qual) > p.max_ee_rate:
            n_ee += 1
            continue
        kept.append((qname, seq, qual))
    stats = {"reads_in": n_in, "dropped_length": n_len, "dropped_ee": n_ee,
             "passed_prefilter": len(kept)}
    if p.quality_trim != "off":
        stats.update(reads_trimmed=n_trim, dropped_no_quality_segment=n_empty)
    return kept, stats


def downsample(reads: list, p: Params) -> tuple[list, dict]:
    if p.downsample_to is None or len(reads) <= p.downsample_to:
        return reads, {"downsampled_to": len(reads), "downsample_applied": False}
    rng = random.Random(p.downsample_seed)
    sample = rng.sample(reads, p.downsample_to)
    return sample, {"downsampled_to": p.downsample_to, "downsample_applied": True}
