"""Insertions and deletions, tested against the controls' own indel rate.

Nanopore reads carry about one indel error each, mostly one base long and
nearly half of them in homopolymers, so no quality-based null can separate a
real indel from noise. What makes it possible is that these errors are
systematic: the same deletion in the same homopolymer recurs at a similar rate
in every library. A mutation-free control therefore measures, for each exact
indel, how often it appears by error, and a sample is tested against that.

  rate(indel) = max(class rate, site rate shrunk toward the class rate) x c

The class is (insertion or deletion, length up to 4, homopolymer run up to 5):
an indel the controls never showed still gets the rate of its kind. c inflates
the rates for between-library variation and is fitted as the quality-calibrated
substitution model's is: each control is tested against the others, and c is
the smallest value at which no more than 1% of its indels reach p < 0.01.
Candidates are corrected by Benjamini-Hochberg.

Indels are counted only at quality >= p.indel_min_quality (the weakest
inserted base, or the weaker base flanking a deletion): nanopore indel errors
mostly come with low-quality flanks and a real indel is read cleanly. One- and
two-read indels are not called (p.indel_min_reads): at about 500 reads the
controls show around 225 distinct one-read indels per library, and a few real
two-read indel molecules that are not AID's.
Under the Sanger platform (one read per clone, low indel error) every indel
whose bases reach the quality threshold is called, as substitutions are.

See docs/calibration/INDELS.md.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass

from scipy.stats import binom

from .variants import _benjamini_hochberg, homopolymer_run

QSCALE_GRID = (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0)
DEFAULT_QSCALE = 3.0      # used with one control; measured on the J꜀4 knockouts
SITE_PRIOR_TRIALS = 50.0  # pseudo-reads pulling a site's rate toward its class


@dataclass
class Indel:
    pos: int               # 1-based: first deleted base, or the base an insertion precedes
    kind: str              # "ins" or "del"
    seq: str               # inserted or deleted bases
    supporting_reads: int
    depth: int
    frequency: float
    p_value: float
    plus_reads: int = 0
    minus_reads: int = 0
    homopolymer: int = 1   # identical-base run at the site in the reference

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def frameshift(self) -> bool:
        return len(self.seq) % 3 != 0


def normalise(ref: str, rpos: int, op: str, x) -> tuple[int, str, str]:
    """Left-align an indel so the same event lines up in every read.

    An aligner may place a deletion anywhere along a homopolymer or tandem
    repeat; shifting it as far left as the sequence allows gives one
    representation. Returns (0-based position, "ins"/"del", bases).
    """
    if op == "D":
        n = int(x)
        while rpos > 0 and ref[rpos - 1] == ref[rpos + n - 1]:
            rpos -= 1
        return rpos, "del", ref[rpos:rpos + n]
    s = str(x)
    while rpos > 0 and s and ref[rpos - 1] == s[-1]:
        s = s[-1] + s[:-1]
        rpos -= 1
    return rpos, "ins", s


def context(ref: str, pos0: int, kind: str, seq: str) -> tuple:
    hp = homopolymer_run(ref, pos0 if kind == "del" else max(0, pos0 - 1))
    return kind, min(len(seq), 4), min(hp, 5)


class Pileup:
    """Indel events in one library, and how many reads span each."""

    def __init__(self, aligned, ref: str, lo: int, hi: int, min_quality: int = 0):
        self.ref = ref
        self.events: dict = defaultdict(list)       # (pos0, kind, seq) -> [(q, strand)]
        starts, ends = [], []
        for r in aligned:
            starts.append(r.r_st)
            ends.append(r.r_en)
            s = getattr(r, "strand", 1)
            for rpos, op, x, q in getattr(r, "indels", ()):
                if q < min_quality:
                    continue
                pos0, kind, seq = normalise(ref, rpos, op, x)
                if lo <= pos0 <= hi:
                    self.events[(pos0, kind, seq)].append((q, s))
        self.starts, self.ends = sorted(starts), sorted(ends)
        self.n = len(starts)

    def depth(self, pos0: int, kind: str, seq: str) -> int:
        """Reads that cover the base before the event and the base after it."""
        first = pos0 - 1
        last = pos0 + (len(seq) if kind == "del" else 0)
        started = bisect_right(self.starts, first)          # r_st <= first
        ended = bisect_left(self.ends, last + 1)             # r_en < last + 1
        return max(0, started - ended)


class IndelNull:
    """Per-indel and per-class error rates from control pileups."""

    def __init__(self, piles: list[Pileup], ref: str):
        self.ref = ref
        self.site_k: dict = defaultdict(int)
        self.site_n: dict = defaultdict(int)
        cls_k: dict = defaultdict(int)
        self.qscale = 1.0
        keys = set()
        for pl in piles:
            keys |= set(pl.events)
        for pl in piles:
            for key in keys:
                n = pl.depth(*key)
                k = len(pl.events.get(key, ()))
                if n and k / n >= 0.5:
                    continue                    # a fixed difference, not error
                self.site_k[key] += k
                self.site_n[key] += n
                cls_k[context(ref, *key)] += k
        # the class denominator: read-positions at which such an event could
        # occur, i.e. analysed reads x positions of that context
        n_reads = sum(pl.n for pl in piles)
        positions: dict = defaultdict(int)
        for pos0 in range(1, len(ref)):
            for kind in ("del", "ins"):
                for ln in (1, 2, 3, 4):
                    positions[context(ref, pos0, kind, "N" * ln)] += 1
        self.cls_rate = {c: (cls_k.get(c, 0) + 0.5) / (positions[c] * max(1, n_reads) + 1.0)
                         for c in positions}

    def rate(self, key: tuple) -> float:
        g = self.cls_rate.get(context(self.ref, *key), 1e-6)
        k, n = self.site_k.get(key, 0), self.site_n.get(key, 0)
        site = (k + SITE_PRIOR_TRIALS * g) / (n + SITE_PRIOR_TRIALS) if n else g
        return min(1.0, self.qscale * max(g, site))


def _binom_sf(k: int, n: int, q: float) -> float:
    return float(binom.sf(k - 1, n, min(max(q, 1e-12), 1.0))) if k > 0 else 1.0


def call(pile: Pileup, null: IndelNull | None, p,
         sanger: bool = False) -> tuple[list[Indel], dict]:
    """Indels in one library: tested against `null`, or (Sanger) thresholded."""
    ref = pile.ref
    cands, pvals = [], []
    germline = 0
    for key, reads in pile.events.items():
        pos0, kind, seq = key
        if pos0 < p.end_trim_bp or pos0 >= len(ref) - p.end_trim_bp:
            continue
        d = pile.depth(*key)
        k = len(reads)
        if not d or k > d or (not sanger and k < p.indel_min_reads):
            continue
        if not sanger and k / d >= p.germline_freq_threshold:
            germline += 1
            continue
        if sanger:
            pv = 0.0
        elif null is None:
            continue
        else:
            pv = _binom_sf(k, d, null.rate(key))
        plus = sum(1 for _q, s in reads if s == 1)
        cands.append(Indel(pos0 + 1, kind, seq, k, d, k / d, pv, plus, k - plus,
                           homopolymer_run(ref, pos0)))
        pvals.append(pv)
    if sanger:
        keep = [True] * len(cands)
    elif p.multiple_testing == "bh":
        keep = _benjamini_hochberg(pvals, p.fdr_alpha)
    else:
        keep = [pv <= p.max_pvalue for pv in pvals]
    called = sorted((c for c, kp in zip(cands, keep, strict=True) if kp),
                    key=lambda v: (v.pos, v.kind, v.seq))
    stats = {"indel_candidates": len(cands), "indels_called": len(called),
             "indel_germline_removed": germline,
             "indel_model": "sanger_threshold" if sanger else
             ("control_rate" if null is not None else "none")}
    return called, stats


def calibrate(piles: list[Pileup], ref: str, p, alpha: float = 0.01) -> tuple[float, dict]:
    """Smallest rate inflation at which controls tested against each other look clean."""
    if len(piles) < 2:
        return DEFAULT_QSCALE, {"indel_qscale_source": "default (one control file)"}
    folds = [(t, IndelNull([o for j, o in enumerate(piles) if j != i], ref))
             for i, t in enumerate(piles)]
    trace = {}
    for c in QSCALE_GRID:
        small = m = 0
        for test, null in folds:
            null.qscale = c
            for key, reads in test.events.items():
                d = test.depth(*key)
                if not d or len(reads) / d >= p.germline_freq_threshold:
                    continue
                m += 1
                small += _binom_sf(len(reads), d, null.rate(key)) < alpha
        trace[str(c)] = (small, m)
        if small <= alpha * m:
            return c, {"indel_qscale_source": "leave-one-control-out", "indel_qscale_trace": trace}
    return QSCALE_GRID[-1], {"indel_qscale_source": "leave-one-control-out (grid maximum)",
                             "indel_qscale_trace": trace}
