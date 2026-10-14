"""Base-quality recalibration from a run's own mutation-free controls (opt-in, Params.recalibrate).

A nanopore base's reported quality says less than its neighbours do: on the J_H4 and N114
AID-null controls, a base reported at Q15-34 is miscalled 4-8 times as often when a base within
two positions is below Q10 as when every neighbour is at Q20 or more. This module measures, for
each reported-quality band, how much more or less often bases are miscalled in each neighbourhood
than in the band as a whole, and shifts each base's quality by that amount (a shift, as in GATK's
BQSR, not a remap, so the basecaller's own ranking within a band is kept). Shrinkage is counted in
expected errors, so a sparse cell is not promoted for showing no errors.

The table is fitted on the same control reads as the null it belongs to (null_model.fit_aligned
stores it as `NullModel.recal`), and both the null's control qualities and the tested sample's
qualities are read through it, so a control tested against the others (leave-one-control-out
calibration, decoys) is always recalibrated with a table that never saw its own reads. Only the
quality-calibrated test (error_model "qcal") reads effective qualities.

It is off by default because it does not help measurably on R10.4.1 SUP reads: on the planted-
mutation benchmark (docs/calibration/RECAL.md) recall and false calls moved by about 0.02 and two
calls per library, in both directions. The fitted shifts are 0 at Q35+, where most calls are
made, and the per-site control rates already absorb the rest.
"""
from __future__ import annotations

from collections import defaultdict
from math import log10

QMAX = 60
# neighbourhood: the lowest reported quality within FLANK aligned bases on either side
FLANK = 2
NB_EDGES = (10, 20)            # bins: < 10, 10-19, >= 20
N_NB = len(NB_EDGES) + 1
BASES = ("A", "C", "G", "T")
# pseudo-ERRORS pulling a (quality band, neighbourhood) cell toward its band's rate. Counted in
# errors, not bases: at Q35+ errors are ~1 in 10^4, so a cell of a few thousand bases with no
# error says nothing, and a prior in bases would let it look error-free and be promoted.
PRIOR_ERRORS = 20.0
# reported-quality bands within which the neighbourhood shift is measured
Q_BANDS = ((0, 14), (15, 24), (25, 34), (35, QMAX))


def _nb_bin(qmin: int) -> int:
    for k, edge in enumerate(NB_EDGES):
        if qmin < edge:
            return k
    return N_NB - 1


def neighbour_bins(pairs) -> list[int]:
    """For each aligned base of a read, the bin of the lowest quality among its neighbours."""
    qs = [q for _r, _b, q in pairs]
    n = len(qs)
    out = []
    for i in range(n):
        a, b = max(0, i - FLANK), min(n, i + FLANK + 1)
        nb = [qs[j] for j in range(a, b) if j != i]
        out.append(_nb_bin(min(nb)) if nb else N_NB - 1)
    return out


class Recal:
    """Effective quality for each (reported quality, neighbourhood bin)."""

    def __init__(self, table: list[list[int]], meta: dict | None = None):
        self.table = table                  # table[q][nb] -> effective quality
        self.meta = dict(meta or {})

    def quals(self, pairs) -> list[int]:
        t = self.table
        return [t[min(QMAX, max(0, int(q)))][nb] for (_r, _b, q), nb in zip(pairs, neighbour_bins(pairs),
                                                                              strict=True)]

    def to_json(self) -> dict:
        return {"table": self.table, "meta": self.meta}

    @classmethod
    def from_json(cls, data: dict | None) -> Recal | None:
        return cls(data["table"], data.get("meta")) if data else None


def fit(aligned_per_file: list[list], ref: str, lo: int, hi: int, germline_frac: float) -> Recal:
    """Fit the table on control reads. Sites where `germline_frac` or more of the pooled reads
    differ from the reference are fixed differences, not errors, and are left out."""
    ref_u = ref.upper()
    cnt: dict = defaultdict(lambda: defaultdict(int))
    dep: dict = defaultdict(int)
    for aligned in aligned_per_file:
        for r in aligned:
            for rpos, base, _q in r.pairs:
                if base in BASES and lo <= rpos <= hi:
                    cnt[rpos][base] += 1
                    dep[rpos] += 1
    germ = {(rp, a) for rp, c in cnt.items() for a, k in c.items()
            if a != ref_u[rp] and k / dep[rp] >= germline_frac}
    err = [[0] * N_NB for _ in range(QMAX + 1)]
    tot = [[0] * N_NB for _ in range(QMAX + 1)]
    for aligned in aligned_per_file:
        for r in aligned:
            for (rpos, base, q), nb in zip(r.pairs, neighbour_bins(r.pairs), strict=True):
                if base not in BASES or not lo <= rpos <= hi or (rpos, base) in germ:
                    continue
                q = min(QMAX, max(0, int(q)))
                tot[q][nb] += 1
                err[q][nb] += base != ref_u[rpos]
    # The reported quality keeps its resolution: each base is shifted by how much more or less
    # often bases of its quality band are miscalled in its neighbourhood than in the band as a
    # whole (a shift, as in BQSR, not a remap - a remap collapses the high qualities, which the
    # controls are too few to tell apart, and loses the basecaller's own ranking there).
    shift = []
    for band_lo, band_hi in Q_BANDS:
        e_band = sum(err[q][nb] for q in range(band_lo, band_hi + 1) for nb in range(N_NB))
        n_band = sum(tot[q][nb] for q in range(band_lo, band_hi + 1) for nb in range(N_NB))
        rate = (e_band + 1.0) / (n_band + 1.0)
        row = []
        for nb in range(N_NB):
            e = sum(err[q][nb] for q in range(band_lo, band_hi + 1))
            n = sum(tot[q][nb] for q in range(band_lo, band_hi + 1))
            prior = PRIOR_ERRORS / rate                   # bases holding PRIOR_ERRORS errors at the band rate
            cell = (e + PRIOR_ERRORS) / (n + prior)
            row.append(round(10 * log10(rate / cell)))
        shift.append(row)
    table = []
    for q in range(QMAX + 1):
        band = next(i for i, (a, b_) in enumerate(Q_BANDS) if a <= q <= b_)
        table.append([min(QMAX, max(0, q + shift[band][nb])) for nb in range(N_NB)])
    meta = {"shift_db": shift, "bases": sum(map(sum, tot)), "errors": sum(map(sum, err)), "germline_sites": len(germ)}
    return Recal(table, meta)
