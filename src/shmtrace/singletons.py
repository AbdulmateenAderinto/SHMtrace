"""Single-read mutations, estimated in aggregate against a control background.

A mutation carried by one read cannot be told from a sequencing error by any
per-site test: at 500 reads the controls show the same kind of single
high-quality miscall at some position or other in every library
(docs/calibration/ERROR_MODELS.md). Yet single-read mutations are most of what
the paper and SHMTool count as unique mutations in a polyclonal assay.

So this module does not call them one by one. It counts them:

  observed   substitutions carried by exactly one read at quality >= q_min,
             and not already called by the per-site test
  expected   the same count in mutation-free controls, per base sequenced,
             times the bases sequenced in this sample
  efficiency the share of this sample's correctly read bases with quality
             >= q_min: a real mutated base is read at ordinary quality, so this
             is the share of real single-read mutations the quality cut keeps
  estimate   (observed - expected) / efficiency

Errors cluster at low quality and real bases do not, so above q_min the
background is a few events per library and the excess is the signal. The
estimate should not depend on q_min; that it does not, between Q30 and Q40 on
JH4, is the check that the efficiency correction is right
(docs/calibration/SINGLE_READ.md).

Samples and controls must be at the same depth: whether a mutation is seen in
exactly one read depends on how many reads there are.
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field

BASES = "ACGT"
# a site where a control shows an alt in at least this share of reads is a
# haplotype or contaminant there, not error, and is left out (as in null_model)
POLYMORPHIC_IN_CONTROLS = 0.02
# the quality cut is lowered (in steps of 5, not below MIN_Q) until it keeps at
# least this share of every library's correctly read bases
MIN_EFFICIENCY = 0.5
MIN_Q = 20


def _list_by_key():
    # module-level factories, not lambdas, so a Tally can be pickled (parallel.py)
    return defaultdict(list)


def _count_by_key():
    return defaultdict(int)


@dataclass
class Tally:
    """Per-position counts from one library's aligned reads, inside the window."""
    depth: dict = field(default_factory=lambda: defaultdict(int))
    # pos -> alt -> [(quality, strand)] of every read showing that alt
    alts: dict = field(default_factory=lambda: defaultdict(_list_by_key))
    # pos -> qualities of the bases matching the reference, as a histogram
    ref_quals: dict = field(default_factory=lambda: defaultdict(_count_by_key))


def tally(aligned, ref: str, lo: int, hi: int) -> Tally:
    t = Tally()
    for r in aligned:
        s = getattr(r, "strand", 1)
        for pos, base, q in r.pairs:
            if not (lo <= pos <= hi) or base not in BASES:
                continue
            rb = ref[pos].upper()
            if rb not in BASES:
                continue
            t.depth[pos] += 1
            if base == rb:
                t.ref_quals[pos][q] += 1
            else:
                t.alts[pos][base].append((q, s))
    return t


def excluded_sites(sample: Tally, controls: list[Tally], germline_freq: float) -> set[int]:
    """Positions where a non-reference base is not error: germline or polymorphic."""
    out = set()
    for pos, alts in sample.alts.items():
        d = sample.depth[pos]
        if any(len(v) / d >= germline_freq for v in alts.values()):
            out.add(pos)
    for c in controls:
        for pos, alts in c.alts.items():
            d = c.depth[pos]
            if any(len(v) / d >= POLYMORPHIC_IN_CONTROLS for v in alts.values()):
                out.add(pos)
    return out


def count(t: Tally, q_min: int, skip: set[int], called: set | None = None) -> dict:
    """Single-read events at quality >= q_min, and the bases they arose from.

    An event is a (position, alt) carried by exactly one read at quality
    >= q_min. Low-quality reads showing the same alt are ignored: at nanopore
    error rates a column of 500 reads often holds one or two of them, and
    requiring the alt in exactly one read overall loses about 15% of real
    single-read mutations in the spike-in. Events the per-site test already
    called are left out, so they are not counted twice.
    """
    called = called or set()
    events = []
    for pos, alts in t.alts.items():
        if pos in skip:
            continue
        for alt, reads in alts.items():
            if (pos, alt) in called:
                continue
            hq = [r for r in reads if r[0] >= q_min]
            if len(hq) == 1:
                events.append((pos, alt, hq[0][0], hq[0][1]))
    bases = sum(d for pos, d in t.depth.items() if pos not in skip)
    ref_all = ref_hi = 0
    for pos, hist in t.ref_quals.items():
        if pos in skip:
            continue
        for q, n in hist.items():
            ref_all += n
            if q >= q_min:
                ref_hi += n
    return {"events": events, "bases": bases,
            "efficiency": ref_hi / ref_all if ref_all else 0.0}


def estimate(sample: Tally, controls: list[Tally], q_min: int, germline_freq: float,
             ref: str, called: set | None = None, z: float = 1.96) -> dict | None:
    """Background-subtracted count of single-read mutations, with a 95% interval.

    Variance: the sample's own count is Poisson in the signal and overdispersed
    in the background (phi, from the spread between control libraries when
    there are three or more), and the background rate carries the uncertainty
    of the control count it was measured from.
    """
    if not controls:
        return None
    skip = excluded_sites(sample, controls, germline_freq)
    # Older basecallers cap quality lower; if the cut keeps under half of any
    # library's correct bases, step it down so the estimate is not built on a
    # sliver of the reads
    asked = q_min
    while True:
        s = count(sample, q_min, skip, called)
        cs = [count(c, q_min, skip) for c in controls]
        if q_min <= MIN_Q or min([s["efficiency"]] + [c["efficiency"] for c in cs]) \
                >= MIN_EFFICIENCY:
            break
        q_min = max(MIN_Q, q_min - 5)
    kc = sum(len(c["events"]) for c in cs)
    nc = sum(c["bases"] for c in cs)
    if not nc or not s["bases"] or s["efficiency"] <= 0:
        return None
    rate = kc / nc
    obs = len(s["events"])
    exp = s["bases"] * rate
    phi = 1.0
    if len(cs) >= 3 and kc:
        x2 = sum((len(c["events"]) - c["bases"] * rate) ** 2 / (c["bases"] * rate)
                 for c in cs if c["bases"])
        phi = max(1.0, x2 / (len(cs) - 1))
    eff = s["efficiency"]
    est = (obs - exp) / eff
    var = obs + (phi - 1.0) * exp + phi * exp * exp / max(kc, 1)
    se = math.sqrt(max(var, 0.0)) / eff
    # the spectrum of the excess: observed single-read events by class, and the
    # background expected for each class from the controls' own events
    obs_cls = defaultdict(int)
    for pos, alt, _q, _st in s["events"]:
        obs_cls[f"{ref[pos].upper()}>{alt}"] += 1
    ctrl_cls = defaultdict(int)
    for c in cs:
        for pos, alt, _q, _st in c["events"]:
            ctrl_cls[f"{ref[pos].upper()}>{alt}"] += 1
    exp_cls = {k: round(v * s["bases"] / nc, 3) for k, v in ctrl_cls.items()}
    return {
        "q_min": q_min,
        "q_min_requested": asked,
        "observed": obs,
        "expected": round(exp, 3),
        "efficiency": round(eff, 4),
        "estimate": round(est, 2),
        "ci_low": round(est - z * se, 2),
        "ci_high": round(est + z * se, 2),
        "overdispersion": round(phi, 3),
        "bases": s["bases"],
        "control_events": kc,
        "control_bases": nc,
        "sites_excluded": len(skip),
        "observed_by_class": dict(obs_cls),
        "expected_by_class": exp_cls,
    }
