"""Competitive alignment with minimap2 (via the `mappy` Python binding).

Replaces the Geneious mapper. Target reference(s) and contaminant reference(s)
are indexed together;
each read is assigned to its best-scoring hit. Reads that
map best to a contaminant, that align too short, or that are chimeric
(supplementary alignments) are dropped. Retained reads are returned with the
per-reference-position aligned pairs needed for pileup-based variant calling.
"""
from __future__ import annotations

import functools
import os
import tempfile
from dataclasses import dataclass, field

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
    # (rpos, "I" or "D", bases, quality): an insertion before reference
    # position rpos, or a deletion of the reference bases starting there, as
    # aligned (indels.normalise left-aligns them). Quality is the weakest
    # inserted base, or the weaker of the two bases flanking a deletion.
    indels: list = field(default_factory=list)
    # where the alignment starts and ends in `seq`, so the unaligned 5' sequence
    # (for J꜀4, the V(D)J junction) and 3' sequence (a UMI) can be read back
    q_st: int = 0
    q_en: int = -1


def diagnostic_sites(reference: str, contaminant_fastas: list[str],
                     preset: str = "map-ont") -> dict[int, set[str]]:
    """Positions (0-based, on the reference) where a contaminant differs from it,
    with the contaminant's base(s) there. Only aligned columns; indels are left out.

    These identify a read from a contaminant or paralogue that the competitive
    alignment kept: one longer than the contaminant record (it scores higher on
    the full-length reference), or a PCR chimera with a contaminant stretch.
    """
    return dict(_diagnostic_sites(reference, tuple(contaminant_fastas or ()), preset))


@functools.lru_cache(maxsize=16)
def _diagnostic_sites(reference: str, contaminant_fastas: tuple, preset: str):
    ref = read_fasta(reference)[0][1].upper()
    aligner = mp.Aligner(seq=ref, preset=preset)
    sites: dict[int, set[str]] = {}
    for f in contaminant_fastas:
        for _name, con in read_fasta(f):
            con = con.upper()
            for h in aligner.map(con):
                if not h.is_primary:
                    continue
                q = con if h.strand == 1 else _revcomp(con)
                rp, qp = h.r_st, (h.q_st if h.strand == 1 else len(con) - h.q_en)
                for length, op in h.cigar:
                    if op in (0, 7, 8):
                        for k in range(length):
                            b = q[qp + k]
                            if b in "ACGT" and b != ref[rp + k]:
                                sites.setdefault(rp + k, set()).add(b)
                        rp += length
                        qp += length
                    elif op == 1:
                        qp += length
                    elif op in (2, 3):
                        rp += length
    return tuple((k, frozenset(v)) for k, v in sorted(sites.items()))


def drop_contaminant_alleles(aligned: list, sites: dict, min_alleles: int,
                             germline_frac: float = 0.2, fates: dict | None = None):
    """Drop reads carrying `min_alleles` or more contaminant bases at diagnostic sites.

    A site where `germline_frac` or more of this sample's reads carry the
    contaminant's base is a germline difference the sample shares with the
    contaminant, not evidence, and is left out. Returns (kept, dropped count).
    """
    if not sites or min_alleles <= 0:
        return aligned, 0
    cov: dict[int, int] = {}
    car: dict[int, int] = {}
    for r in aligned:
        for rp, b, _q in r.pairs:
            if rp in sites:
                cov[rp] = cov.get(rp, 0) + 1
                car[rp] = car.get(rp, 0) + (b in sites[rp])
    use = {s: v for s, v in sites.items() if cov.get(s) and car.get(s, 0) / cov[s] < germline_frac}
    kept, dropped = [], 0
    for r in aligned:
        hits = sum(1 for rp, b, _q in r.pairs if rp in use and b in use[rp])
        if hits >= min_alleles:
            dropped += 1
            if fates is not None:
                fates[r.qname] = {"fate": "contaminant alleles", "alleles": hits}
            continue
        kept.append(r)
    return kept, dropped


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


def align_reads(reads, aligner, target_names, contam_names, p: Params,
                fates: dict | None = None):
    """reads = iterable of (qname, seq, qual). Yields AlignedRead for kept reads.

    Returns a stats dict via the .stats attribute on the generator's closure —
    here we just yield and let the caller tally
    see pipeline.py.

    Pass `fates` to record why each read was dropped, keyed by read name. A
    read that aligned but was rejected is recorded with its alignment, so the
    read browser can show it in place against the reference.
    """
    def note(name, why, hit=None, seq=None, qual=None):
        if fates is None:
            return
        rec = {"fate": why}
        if hit is not None:
            rec.update(ctg=hit.ctg, r_st=hit.r_st, r_en=hit.r_en, strand=hit.strand)
            if seq is not None:
                rec["pairs"] = _oriented_pairs(hit, seq, qual)
        fates[name] = rec

    for qname, seq, qual in reads:
        best = None
        n_primary = 0
        for hit in aligner.map(seq):
            if hit.is_primary:
                n_primary += 1
                if best is None or hit.mlen > best.mlen:
                    best = hit
        if best is None:
            note(qname, "unmapped")
            continue
        if p.drop_supplementary and n_primary > 1:
            note(qname, "chimeric", best, seq, qual)
            continue                                  # chimera / structural anomaly
        if best.ctg in contam_names:
            note(qname, "contaminant", best)
            continue                                  # best hit is a contaminant
        if best.ctg not in target_names:
            note(qname, "off target", best)
            continue
        if (best.r_en - best.r_st) < p.min_align_len:
            note(qname, "alignment too short", best, seq, qual)
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
        pairs, indels = _aligned_pairs(best, wseq, wquals, q_start, with_indels=True)
        yield AlignedRead(qname, best.ctg, best.r_st, best.r_en, pairs, wseq, best.strand,
                          indels, q_start, q_start + (best.q_en - best.q_st))


def _oriented_pairs(hit, seq: str, qual: str):
    """The aligned (ref position, base, quality) triples for a hit, any strand."""
    if hit.strand == -1:
        wseq, wquals, q_start = _revcomp(seq), phred(qual)[::-1], len(seq) - hit.q_en
    else:
        wseq, wquals, q_start = seq, phred(qual), hit.q_st
    return _aligned_pairs(hit, wseq, wquals, q_start)


def _aligned_pairs(hit, seq: str, quals: list[int], q_start: int, with_indels: bool = False):
    """Walk the CIGAR to produce (rpos, base, qval) at aligned ref columns.

    `seq`/`quals` are already in reference-forward orientation and `q_start` is
    the query offset within them. Only aligned columns (M/=/X) yield a base +
    its Phred quality for quality-aware pileup. With `with_indels`, insertions
    and deletions are also returned, as AlignedRead.indels describes.
    """
    pairs, indels = [], []
    qpos, rpos = q_start, hit.r_st
    for length, op in hit.cigar:            # op: 0=M,1=I,2=D,3=N,4=S,7==,8=X
        if op in (0, 7, 8):                 # alignment match/mismatch
            for k in range(length):
                pairs.append((rpos + k, seq[qpos + k], quals[qpos + k]))
            qpos += length
            rpos += length
        elif op == 1:                       # insertion to reference
            if with_indels:
                indels.append((rpos, "I", seq[qpos:qpos + length],
                               min(quals[qpos:qpos + length])))
            qpos += length
        elif op == 2:                       # deletion in read
            if with_indels:
                flank = [quals[i] for i in (qpos - 1, qpos) if 0 <= i < len(quals)]
                indels.append((rpos, "D", length, min(flank) if flank else 0))
            rpos += length
        elif op == 3:                       # skip
            rpos += length
        elif op == 4:                       # soft clip
            qpos += length
    return (pairs, indels) if with_indels else pairs


def trim_messy_ends(read: AlignedRead, ref: str, window: int,
                    max_mismatch: int) -> AlignedRead | None:
    """Cut an aligned Sanger read back to where it is clean, from each end.

    The start of the clean part is the first aligned column from which `window`
    consecutive columns hold no ambiguity code, no indel and at most
    `max_mismatch` mismatches against the reference; the end is found the same
    way from the other side. Mutations inside the clean part are kept. Returns
    None if no clean stretch exists.
    """
    pairs = read.pairs
    n = len(pairs)
    if window <= 0:
        return read
    if n < window:
        return None
    indel_at = {rpos for rpos, *_ in read.indels}

    def clean(i: int) -> bool:
        seg = pairs[i:i + window]
        lo_r, hi_r = seg[0][0], seg[-1][0]
        if hi_r - lo_r != window - 1:
            return False                      # a deletion inside the run
        if any(lo_r <= r <= hi_r + 1 for r in indel_at):
            return False
        mm = 0
        for rpos, base, _q in seg:
            if base not in "ACGT":
                return False
            mm += base != ref[rpos].upper()
        return mm <= max_mismatch

    start = next((i for i in range(n - window + 1) if clean(i)), None)
    if start is None:
        return None
    end = next(j for j in range(n - window, start - 1, -1) if clean(j)) + window - 1
    kept = pairs[start:end + 1]
    lo_r, hi_r = kept[0][0], kept[-1][0]
    indels = [e for e in read.indels if lo_r < e[0] <= hi_r]
    return AlignedRead(read.qname, read.ctg, lo_r, hi_r + 1, kept, read.seq, read.strand,
                       indels)
