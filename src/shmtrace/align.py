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
from array import array
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

    def compact(self) -> AlignedRead:
        """Hold the aligned columns packed, about 6 bytes a base in place of about 90.

        A list of one (rpos, base, qval) tuple per aligned base is what every pileup reads,
        and a project keeps every usable read of every library at once: 15 Illumina libraries
        of 100,000 reads came to 44 GB that way, and the garbage collector walking those
        tuples is what then stalls a run. Packed, `pairs` builds the list afresh on each use
        (about as long as one pass over it), so a change made to that list is not kept: a
        read is compacted only where nothing edits its columns. Columns that cannot be packed
        (an empty base, a position or quality out of range) are left as they are.
        """
        if "_pairs" not in self.__dict__:
            return self
        cols = self.__dict__["_pairs"]
        try:
            rpos = array("i", [c[0] for c in cols])
            bases = "".join(c[1] for c in cols)
            quals = bytes(c[2] for c in cols)
            if len(bases) != len(cols) or not all(type(c[0]) is int and type(c[2]) is int for c in cols):
                raise ValueError
        except (ValueError, TypeError, OverflowError, IndexError):
            return self
        del self.__dict__["_pairs"]
        self.__dict__["_packed"] = (rpos, bases, quals)
        return self

    def __reduce__(self):
        # Reads cross process boundaries in bulk (parallel.py), and a list of one tuple per
        # aligned base is slow to pickle: pack the columns as two arrays and a string instead.
        # The rebuilt read is equal to this one, and arrives compacted; anything unusual is
        # pickled as it is. Attributes set on a read after it was built (umi.py's `members`,
        # the reads a consensus molecule came from) travel as pickle state, so a worker
        # process returns the same reads as a serial run (audit F15).
        extra = {k: v for k, v in self.__dict__.items() if k not in _OWN} or None
        packed = self.__dict__.get("_packed")
        if packed is None:
            packed = AlignedRead(self.qname, self.ctg, 0, 0, self.pairs, "").compact().__dict__.get("_packed")
        if packed is None:
            return (AlignedRead, (self.qname, self.ctg, self.r_st, self.r_en, self.pairs, self.seq,
                                  self.strand, self.indels, self.q_st, self.q_en), extra)
        rpos, bases, quals = packed
        return (_unpack_read, (self.qname, self.ctg, self.r_st, self.r_en, rpos.tobytes(), bases, quals,
                               self.seq, self.strand, self.indels, self.q_st, self.q_en), extra)


def _get_pairs(self) -> list:
    cols = self.__dict__.get("_pairs")
    if cols is not None:
        return cols
    rpos, bases, quals = self.__dict__["_packed"]
    return list(zip(rpos.tolist(), bases, quals))


def _set_pairs(self, cols) -> None:
    self.__dict__.pop("_packed", None)
    self.__dict__["_pairs"] = cols


# `pairs` is a property set after the dataclass is made, so the generated __init__, __eq__ and
# __repr__ go through it: a read built from a list keeps that list; a compacted one rebuilds it
AlignedRead.pairs = property(_get_pairs, _set_pairs, doc="list of (rpos, base, qval) at aligned columns")
_OWN = set(AlignedRead.__dataclass_fields__) | {"_pairs", "_packed"}


def _unpack_read(qname, ctg, r_st, r_en, rpos, bases, quals, seq, strand, indels, q_st, q_en):
    pos = array("i")
    pos.frombytes(rpos)
    r = AlignedRead(qname, ctg, r_st, r_en, None, seq, strand, indels, q_st, q_en)
    r.__dict__.pop("_pairs")
    r.__dict__["_packed"] = (pos, bases, quals)
    return r


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
    sites: dict[int, set[str]] = {}
    for _name, _spans, rec_sites in _record_sites(reference, contaminant_fastas, preset):
        for k, b in rec_sites:
            sites.setdefault(k, set()).update(b)
    return tuple((k, frozenset(v)) for k, v in sorted(sites.items()))


def record_sites(reference: str, contaminant_fastas: list[str],
                 preset: str = "map-ont") -> list[tuple[str, list, dict[int, set[str]]]]:
    """`diagnostic_sites` one contaminant record at a time: (record name, the
    half-open reference spans it aligns to, {position: its base(s) there}).

    The short-read check needs these apart: a site where only the paralogue
    differs says nothing about whether a read is the deposit's contaminant.
    """
    return [(n, list(spans), {k: set(b) for k, b in s})
            for n, spans, s in _record_sites(reference, tuple(contaminant_fastas or ()), preset)]


@functools.lru_cache(maxsize=16)
def _record_sites(reference: str, contaminant_fastas: tuple, preset: str):
    ref = read_fasta(reference)[0][1].upper()
    aligner = mp.Aligner(seq=ref, preset=preset)
    out = []
    for f in contaminant_fastas:
        for name, con in read_fasta(f):
            con = con.upper()
            sites: dict[int, set[str]] = {}
            spans = []
            for h in aligner.map(con):
                if not h.is_primary:
                    continue
                spans.append((h.r_st, h.r_en))
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
            out.append((name, tuple(spans),
                        tuple((k, frozenset(v)) for k, v in sorted(sites.items()))))
    return tuple(out)


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


def drop_unspecific_short_reads(aligned: list, records: list, ref: str, min_margin: int,
                                germline_frac: float = 0.5, fates: dict | None = None):
    """Drop rescued short reads (Params.rescue_short_reads) that do not show they are target.

    For each contaminant record whose span the read overlaps, count the record's
    diagnostic sites where the read has the reference base (t) and where it has the
    record's (c); the read stays only if t - c >= `min_margin` for every such record.
    A record with no diagnostic site under the read gives t = c = 0: there the read
    cannot be told from that contaminant, so it goes. Sites `germline_frac` or more of
    the sample's reads share with the record are left out, as drop_contaminant_alleles
    does. Full-length reads are not judged here. Returns (kept, dropped count).
    """
    short = [r for r in aligned if getattr(r, "short", False)]
    if not short or not records:
        return aligned, 0
    ref = ref.upper()
    union: dict[int, list[int]] = {}
    for i, (_n, _spans, sites) in enumerate(records):
        for k in sites:
            union.setdefault(k, []).append(i)
    cov: dict[tuple[int, int], int] = {}
    car: dict[tuple[int, int], int] = {}
    for r in aligned:
        for rp, b, _q in r.pairs:
            for i in union.get(rp, ()):
                cov[i, rp] = cov.get((i, rp), 0) + 1
                car[i, rp] = car.get((i, rp), 0) + (b in records[i][2][rp])
    shared = {key for key, n in cov.items() if car.get(key, 0) / n >= germline_frac}
    bad = set()
    for r in short:
        bases = {rp: b for rp, b, _q in r.pairs}
        for i, (name, spans, sites) in enumerate(records):
            if not any(st < r.r_en and r.r_st < en for st, en in spans):
                continue
            t = c = 0
            for k, alts in sites.items():
                b = bases.get(k)
                if b is None or (i, k) in shared:
                    continue
                c += b in alts
                t += b == ref[k]
            if t - c < min_margin:
                bad.add(r.qname)
                if fates is not None:
                    fates[r.qname] = {
                        "fate": "alignment too short", "ctg": r.ctg, "r_st": r.r_st,
                        "r_en": r.r_en, "strand": r.strand, "pairs": r.pairs,
                        "reason": f"not told from {name}: {t} reference and {c} contaminant "
                                  f"base(s) at its diagnostic sites, margin below {min_margin}"}
                break
    return [r for r in aligned if r.qname not in bad], len(bad)


def short_read_flaw(hit, read_len: int, n_primary: int, p: Params) -> str | None:
    """Why an alignment spanning less than p.min_align_len is not kept as a short read,
    or None when it is (Params.rescue_short_reads). The target-over-contaminant test
    needs the whole sample and comes later (drop_unspecific_short_reads)."""
    frac = p.short_read_min_fraction
    span = hit.r_en - hit.r_st
    if span < frac * p.min_align_len:
        return f"aligned span {span} bp, below {frac * p.min_align_len:.0f}"
    win = p.analysis_window
    lo, hi = (win[0] - 1, win[1]) if win else (0, hit.ctg_len)
    need = frac * min(p.min_align_len, hi - lo)
    inside = max(0, min(hit.r_en, hi) - max(hit.r_st, lo))
    if inside < need:
        return f"covers {inside} bp of the analysis window, below {need:.0f}"
    cover = (hit.q_en - hit.q_st) / read_len if read_len else 0.0
    if cover < p.short_read_min_query_cover:
        return f"{100 * cover:.0f}% of the read aligned, below {100 * p.short_read_min_query_cover:.0f}%"
    div = hit.NM / hit.blen if hit.blen else 1.0
    if div > p.short_read_max_divergence:
        return f"{100 * div:.1f}% divergent from the reference, above {100 * p.short_read_max_divergence:.1f}%"
    if n_primary > 1:
        return "more than one primary alignment"
    return None


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
    def note(name, why, hit=None, seq=None, qual=None, reason=None):
        if fates is None:
            return
        rec = {"fate": why}
        if hit is not None:
            rec.update(ctg=hit.ctg, r_st=hit.r_st, r_en=hit.r_en, strand=hit.strand)
            if seq is not None:
                rec["pairs"] = _oriented_pairs(hit, seq, qual)
        if reason:
            rec["reason"] = reason
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
        short = (best.r_en - best.r_st) < p.min_align_len
        if short:
            why = (short_read_flaw(best, len(seq), n_primary, p) if p.rescue_short_reads
                   else None)
            if not p.rescue_short_reads or why:
                note(qname, "alignment too short", best, seq, qual, why)
                continue                              # aligned span too short
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
        read = AlignedRead(qname, best.ctg, best.r_st, best.r_en, pairs, wseq, best.strand,
                           indels, q_start, q_start + (best.q_en - best.q_st)).compact()
        if short:
            read.short = True       # kept below min_align_len (rescue_short_reads)
        yield read


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


def pair_span(read, margin: int) -> tuple[int, int]:
    """The half-open range of `read.pairs` a pileup may count, less `margin`
    aligned columns at each end of the read (Params.read_end_margin_bp).

    minimap2 soft-clips a mismatch near a read's end rather than align it, so the
    columns where reads begin and end are covered by fewer mutant reads than
    mutation-free ones and read a lower frequency than the library carries.
    Dropping the same number of columns from every read takes that deficit out of
    the depth as well as out of the counts, at the cost of coverage there. A read
    with nothing left over contributes no bases at all.

    Every pileup that feeds a call applies this: the sample's (variants), the
    controls' background and their decoys (null_model), the project-wide tallies
    (project) and the single-read estimate (singletons). With margin 0 the span
    is the whole read, as every published number had it.
    """
    n = len(read.pairs)
    if margin <= 0:
        return 0, n
    return (margin, n - margin) if n > 2 * margin else (0, 0)


def pair_window(read, margin: int):
    """(the aligned columns a pileup may count, the index of the first in read.pairs).

    `pair_span` as the pileups use it. With margin 0 this is the read's own list,
    not a copy of it, so the hot loops pay nothing for the option being off.
    """
    if margin <= 0:
        return read.pairs, 0
    i0, i1 = pair_span(read, margin)
    return read.pairs[i0:i1], i0


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
