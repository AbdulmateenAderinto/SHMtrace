"""Competitive alignment with minimap2 (via the `mappy` Python binding).

Replaces the Geneious mapper. Target reference(s) and contaminant reference(s)
are indexed together;
each read is assigned to its best-scoring hit. Reads that
map best to a contaminant, that align too short, or that are chimeric
(supplementary alignments) are dropped. Retained reads are returned with the
per-reference-position aligned pairs needed for pileup-based variant calling.
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass

import mappy as mp

from .config import Params
from .io_utils import phred, read_fasta

_COMP = str.maketrans("ACGTN", "TGCAN")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass
class AlignedRead:
    qname: str
    ctg: str                 # reference contig the read was assigned to
    r_st: int                # 0-based start on reference
    r_en: int
    pairs: list              # list of (rpos, base, qval) at aligned columns
    seq: str
    strand: int = 1          # +1 read runs along the reference, -1 reverse-complemented


def build_index(target_fastas: list[str], contaminant_fastas: list[str],
                preset: str = "map-ont"):
    """Write a combined reference and return (aligner, target_names, contam_names)."""
    targets, contams = [], []
    for p in target_fastas:
        targets += read_fasta(p)
    for p in contaminant_fastas or []:
        contams += read_fasta(p)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    for name, seq in targets + contams:
        tmp.write(f">{name}\n{seq}\n")
    tmp.close()
    aligner = mp.Aligner(tmp.name, preset=preset)
    if not aligner:
        raise RuntimeError("Failed to build minimap2 index")
    os.unlink(tmp.name)
    return aligner, {n for n, _ in targets}, {n for n, _ in contams}


def align_reads(reads, aligner, target_names, contam_names, p: Params):
    """reads = iterable of (qname, seq, qual). Yields AlignedRead for kept reads.

    Returns a stats dict via the .stats attribute on the generator's closure —
    here we just yield and let the caller tally
    see pipeline.py.
    """
    for qname, seq, qual in reads:
        best = None
        n_primary = 0
        for hit in aligner.map(seq):
            if hit.is_primary:
                n_primary += 1
                if best is None or hit.mlen > best.mlen:
                    best = hit
        if best is None:
            continue                                  # unmapped
        if p.drop_supplementary and n_primary > 1:
            continue                                  # chimera / structural anomaly
        if best.ctg in contam_names:
            continue                                  # best hit is a contaminant
        if best.ctg not in target_names:
            continue
        if (best.r_en - best.r_st) < p.min_align_len:
            continue                                  # aligned span too short
        # minimap2 CIGAR is reference-forward. For a reverse-strand hit the
        # query bases/qualities must come from the reverse-complemented read,
        # with query coordinates remapped onto it.
        if best.strand == -1:
            wseq = _revcomp(seq)
            wquals = phred(qual)[::-1]
            q_start = len(seq) - best.q_en
        else:
            wseq = seq
            wquals = phred(qual)
            q_start = best.q_st
        pairs = _aligned_pairs(best, wseq, wquals, q_start)
        yield AlignedRead(qname, best.ctg, best.r_st, best.r_en, pairs, wseq, best.strand)


def _aligned_pairs(hit, seq: str, quals: list[int], q_start: int):
    """Walk the CIGAR to produce (rpos, base, qval) at aligned ref columns.

    `seq`/`quals` are already in reference-forward orientation and `q_start` is
    the query offset within them. Insertions/deletions are skipped. Only aligned
    columns (M/=/X) yield a base + its Phred quality for quality-aware pileup.
    """
    pairs = []
    qpos, rpos = q_start, hit.r_st
    for length, op in hit.cigar:            # op: 0=M,1=I,2=D,3=N,4=S,7==,8=X
        if op in (0, 7, 8):                 # alignment match/mismatch
            for k in range(length):
                pairs.append((rpos + k, seq[qpos + k], quals[qpos + k]))
            qpos += length
            rpos += length
        elif op == 1:                       # insertion to reference
            qpos += length
        elif op in (2, 3):                  # deletion / skip in read
            rpos += length
        elif op == 4:                       # soft clip
            qpos += length
    return pairs
