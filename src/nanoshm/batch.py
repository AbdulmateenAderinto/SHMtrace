"""Run several samples together: shared haplotype evidence and a common depth.

Two things only make sense across a set of samples rather than one at a time.

1. **Which co-occurring blocks are haplotypes.** A block that recurs in more
   than one animal, or appears in a mutation-free control, is a haplotype,
   paralogue or contaminant. A block confined to one animal may be an expanded
   clone carrying real SHM. So blocks are detected in every file first, and only
   recurring or control-supported blocks are dropped (phase_mode="auto").

2. **A common analysis depth.** Unique-variant density saturates with depth,
   so samples are only comparable when analysed at the same depth. The common
   depth is the smaller of `downsample_to` and the shallowest sample's usable
   reads after phasing.
"""
from __future__ import annotations

from dataclasses import replace

from . import pipeline
from .config import Params
from .io_utils import read_fasta
from .phasing import matches_known, partition


def recurring_blocks(blocks_by_sample: dict[str, list], min_shared: int) -> list:
    """Blocks seen in at least two different samples."""
    out = []
    names = list(blocks_by_sample)
    for a in names:
        others = [b for n in names if n != a for b in blocks_by_sample[n]]
        for block in blocks_by_sample[a]:
            if matches_known(block, others, min_shared) and not matches_known(block, out, min_shared):
                out.append(block)
    return out


def prepare_batch(samples: dict[str, str], reference: str, contaminants: list[str],
                  p: Params, control_blocks: list | None = None,
                  depth: int | None = None, controls: set | None = None) -> dict:
    """Align and phase every sample; return reads, stats and the common depth.

    `samples` maps a sample name to its FASTQ. `control_blocks` are blocks
    found in mutation-free controls (from null_model.fit), which count as known
    regardless of recurrence. Samples named in `controls` are mutation-free:
    every block in them is dropped and counts as known for the others.
    """
    from .null_model import control_params
    controls = set(controls or ())
    ref_seq = read_fasta(reference)[0][1]
    lo, hi = pipeline._window(p, ref_seq)
    auto = p.phase_reads and p.phase_mode == "auto"
    detect = replace(p, phase_mode="flag") if auto else p

    raw, stats, blocks = {}, {}, {}
    for name, fq in samples.items():
        pp = control_params(p) if name in controls else detect
        aligned, st = pipeline.prepare_reads(fq, reference, contaminants, pp)
        raw[name], stats[name] = aligned, st
        blocks[name] = [[(pos - 1, alt) for pos, alt in b["sites"]] for b in st.get("blocks", [])]

    known = list(control_blocks or [])
    for name in controls:
        known += blocks.get(name, [])
    if auto:
        known += recurring_blocks(blocks, p.phase_min_sites)
        for name in samples:
            if name in controls:
                continue
            raw[name], st = partition(raw[name], ref_seq, p, lo, hi, known)
            stats[name] = {**{k: v for k, v in stats[name].items()
                              if k not in st}, **st}

    usable = min((len(a) for a in raw.values()), default=0)
    if depth is None:
        depth = usable if p.downsample_to is None else min(p.downsample_to, usable)
    return {"reads": raw, "stats": stats, "known_blocks": known, "depth": depth,
            "usable_min": usable}
