"""Pileup-based SNV calling that reproduces the manuscript's variant rules.

Rules (mirrors the Geneious workflow used in the paper):
  * SNVs only (indels excluded when p.exclude_indels).
  * A variant = one unique substitution at one reference position, counted once
    regardless of how many reads carry it (recurrence is reported separately as
    the supporting-read count / raw frequency).
  * Positions where an alternate allele reaches >= germline_freq_threshold are
    treated as germline/parental/fixed and removed (the "< 50%" rule).
  * Each candidate is tested against a sequencing-error null (per-column
    quality, or an empirical per-position rate measured from mutation-free
    control reads) and the p-values are corrected across all candidates by
    Benjamini-Hochberg. See docs/calibration/FINDINGS.md.
  * With error_model="geneious" the test instead follows Geneious's approximate
    variant p-value: a binomial tail at the mean quality of the variant bases
    (see _geneious_pvalue). Used uncorrected at max_pvalue, this is the
    manuscript's calling rule.
  * With error_model="strand" each strand is tested against its own
    control-derived error rate with beta-binomial overdispersion, and the two
    strands' p-values are combined (null_model.NullModel).
  * With error_model="sanger" (one read = one clone) every substitution with
    base quality >= min_base_quality is called without a test; the germline
    rule applies only where enough clones cover the position.
  * Bases below min_base_quality are left out of the pileup.
  * An optional strand filter drops calls whose support is lopsided between
    strands, and calls in homopolymers can be flagged or masked.
  * Variants within end_trim_bp of the reference ends are dropped (N114 = 20).
  * Restricted to analysis_window if given.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import exp, lgamma

import numpy
from scipy.stats import betabinom, binom, chi2, poisson

from .config import Params
from .io_utils import read_fasta

BASES = ("A", "C", "G", "T")


@dataclass
class Variant:
    pos: int          # 1-based reference position
    ref: str
    alt: str
    supporting_reads: int
    depth: int
    frequency: float
    p_value: float
    plus_reads: int = 0       # supporting reads sequenced along the reference
    minus_reads: int = 0      # supporting reads sequenced reverse-complemented
    strand_bias_p: float = 1.0
    homopolymer: int = 1      # longest identical-base run touching the position
    context: str = ""         # reference 5-mer centred on the position


def _poisson_sf(k: int, lam: float) -> float:
    """P(X >= k) for X~Poisson(lam). Used as the false-call null per column/allele.

    lam = expected number of erroneous reads showing this alt, summed from the
    per-base error probabilities of the reads at this column (a Poisson-binomial
    well approximated by Poisson for small per-read error). A real mutation is
    supported by more high-quality reads than sequencing error predicts.
    """
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    # P(X <= k-1)
    term = exp(-lam)
    cdf = term
    for i in range(1, k):
        term *= lam / i
        cdf += term
    return max(0.0, 1.0 - cdf)


def _binom_sf(k: int, n: int, q: float) -> float:
    """P(X >= k) for X ~ Binomial(n, q). Used with an empirical error rate.

    Summed in log space so that the small q and large n typical here do not
    underflow.
    """
    if k <= 0:
        return 1.0
    if q <= 0.0:
        return 0.0
    if q >= 1.0 or k > n:
        return 1.0 if k <= 0 else (1.0 if q >= 1.0 else 0.0)
    from math import expm1, log
    log_q, log_1q = log(q), log(1.0 - q)
    total = 0.0
    for i in range(0, k):
        logc = lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1)
        total += exp(logc + i * log_q + (n - i) * log_1q)
    return max(0.0, min(1.0, -expm1(log(total)) if 0 < total < 1 else 1.0 - total))


def _geneious_pvalue(k: int, depth: int, mean_alt_q: float) -> float:
    """Geneious's approximate variant p-value.

    The Geneious manual describes the p-value as the probability of sequencing
    error producing the observed variant bases, computed by averaging the
    quality of the bases that carry the variant and taking a binomial tail:
    P(X >= k), X ~ Binomial(depth, 10^(-Q_mean/10)). Its worked example (three
    Q20 reads, two carrying the variant) gives 3C2 x 0.01^2 x 0.99 + 0.01^3.
    Note the per-base error is not divided by 3.

    On the bundled JH4 libraries this reproduces the manuscript's Geneious call
    counts to within a few calls per sample (docs/calibration/GENEIOUS.md).
    """
    if k <= 0:
        return 1.0
    e = 10 ** (-mean_alt_q / 10.0)
    return float(binom.sf(k - 1, depth, e))


def _benjamini_hochberg(pvals: list[float], alpha: float) -> list[bool]:
    """Return a keep-mask controlling the false discovery rate at `alpha`.

    The caller tests (analysed positions x 3) hypotheses per sample. For the
    529 bp JH4 window that is 1,587 tests, at which an uncorrected p <= 0.05
    expects ~79 false positives, comparable to the entire observed background.
    """
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    cutoff_rank = 0
    for rank, idx in enumerate(order, start=1):
        if pvals[idx] <= alpha * rank / m:
            cutoff_rank = rank
    keep = [False] * m
    for idx in order[:cutoff_rank]:
        keep[idx] = True
    return keep


def _pbinom_sf(k: int, probs: list[float]) -> float:
    """Exact P(X >= k) for X = sum of independent Bernoulli(probs) (LoFreq's test).

    Dynamic programming over reads, keeping only P(X = 0..k-1), so the cost is
    reads x k. The Poisson version (_poisson_sf) is its approximation.
    """
    if k <= 0:
        return 1.0
    dist = [1.0] + [0.0] * (k - 1)            # P(X = j) for j < k
    for e in probs:
        for j in range(k - 1, 0, -1):
            dist[j] = dist[j] * (1.0 - e) + dist[j - 1] * e
        dist[0] *= (1.0 - e)
    return max(0.0, min(1.0, 1.0 - sum(dist)))


def _betabinom_sf(k: int, n: int, mu: float, rho: float) -> float:
    """P(X >= k), X ~ BetaBinomial(n, mean mu, overdispersion rho); binomial when rho=0."""
    if k <= 0:
        return 1.0
    if n <= 0 or k > n:
        return 0.0
    mu = min(max(mu, 1e-12), 1 - 1e-12)
    if rho <= 0:
        return _binom_sf(k, n, mu)
    s = (1.0 - rho) / rho                      # alpha + beta
    return float(betabinom.sf(k - 1, n, mu * s, (1.0 - mu) * s))


def _fisher(pvals: list[float]) -> float:
    pv = [min(1.0, max(p, 1e-300)) for p in pvals]
    from math import log
    stat = -2.0 * sum(log(p) for p in pv)
    return float(chi2.sf(stat, 2 * len(pv)))


def homopolymer_run(ref: str, pos0: int) -> int:
    """Longest run of identical bases that includes pos0 or ends next to it.

    A substitution beside a run (e.g. the base after AAAA) is as exposed to
    homopolymer miscalls as one inside it, so neighbours' runs count too.
    """
    def run_at(i: int) -> int:
        if i < 0 or i >= len(ref):
            return 0
        b = ref[i]
        lo = i
        while lo > 0 and ref[lo - 1] == b:
            lo -= 1
        hi = i
        while hi < len(ref) - 1 and ref[hi + 1] == b:
            hi += 1
        return hi - lo + 1
    return max(run_at(pos0), run_at(pos0 - 1), run_at(pos0 + 1))


def call_variants(aligned_reads, target_ref_fasta: str, target_name: str,
                  p: Params, null_rates: dict | None = None) -> tuple[list[Variant], dict]:
    ref = dict(read_fasta(target_ref_fasta))[target_name]
    reflen = len(ref)
    model = p.error_model
    strand_null = getattr(null_rates, "has_strand", False)
    notes = []
    if model == "strand" and not strand_null:
        notes.append("error_model 'strand' needs control reads (--null or a panel); "
                     "used the read-quality model instead")
        model = "auto"
    if model == "qcal" and not getattr(null_rates, "has_qcal", False):
        notes.append("error_model 'qcal' needs control reads (--null or a panel); "
                     "used the read-quality model instead")
        model = "auto"
    minq = p.min_base_quality
    if model == "sanger" and minq <= 0:
        minq = 20

    win = p.analysis_window
    lo = (win[0] - 1) if win else 0
    hi = (win[1] - 1) if win else reflen - 1

    counts = defaultdict(lambda: defaultdict(int))
    err_lambda = defaultdict(lambda: defaultdict(float))   # expected false reads / alt
    depth = defaultdict(int)
    s_counts = defaultdict(int)                             # (pos, base, strand) -> reads
    s_depth = defaultdict(int)                              # (pos, strand) -> reads
    s_lambda = defaultdict(float)                           # (pos, alt, strand) -> expected errors
    alt_qsum = defaultdict(lambda: defaultdict(int))        # pos -> base -> summed q
    alt_quals = defaultdict(list)                           # (pos, base, strand) -> qualities
    qual_hist: dict = defaultdict(int)                      # quality -> mismatching bases
    ref_qual_hist: dict = defaultdict(int)                  # quality -> matching bases
    col_probs = defaultdict(list) if model == "pbinom" else None
    n_clones = 0
    for r in aligned_reads:
        if r.ctg != target_name:
            continue
        n_clones += 1
        s = getattr(r, "strand", 1)
        for rpos, base, q in r.pairs:
            if base not in BASES or q < minq:
                continue
            counts[rpos][base] += 1
            depth[rpos] += 1
            if base == ref[rpos].upper() and lo <= rpos <= hi:
                ref_qual_hist[min(int(q), 60)] += 1
            s_counts[(rpos, base, s)] += 1
            s_depth[(rpos, s)] += 1
            alt_qsum[rpos][base] += q
            if base != ref[rpos].upper() and lo <= rpos <= hi:
                # every mismatching base's quality, whatever the model: this is
                # what the quality-calibrated test keys on
                qual_hist[min(int(q), 60)] += 1
            if model in ("qcal", "sanger") and base != ref[rpos]:
                alt_quals[(rpos, base, s)].append(q)
            # per-read chance of an error landing on each of the 3 non-observed
            # bases ~ (10^(-Q/10))/3; contributes to that alt's false-call rate.
            e3 = (10 ** (-q / 10.0)) / 3.0
            if col_probs is not None:
                col_probs[rpos].append((base, e3))
            for alt in BASES:
                if alt != base:
                    err_lambda[rpos][alt] += e3
                    s_lambda[(rpos, alt, s)] += e3

    candidates: list[Variant] = []
    pvals: list[float] = []
    germline_removed = 0
    n_hypotheses = 0
    germline_on = True
    germline_min_depth = 0
    if model == "sanger":
        germline_min_depth = p.sanger_germline_min_clones
        if n_clones < p.sanger_germline_min_clones:
            germline_on = False
            notes.append(f"germline filter disabled: only {n_clones} clone(s) "
                         f"(< {p.sanger_germline_min_clones}), so a shared mutation "
                         "cannot be told from germline")
    for rpos in range(lo, hi + 1):
        d = depth.get(rpos, 0)
        if d == 0:
            continue
        if rpos < p.end_trim_bp or rpos >= reflen - p.end_trim_bp:
            continue
        refbase = ref[rpos].upper()
        if refbase not in BASES:
            continue
        n_hypotheses += 3
        for alt in BASES:
            if alt == refbase:
                continue
            c = counts[rpos].get(alt, 0)
            if c < p.min_variant_reads:
                continue
            freq = c / d
            if (germline_on and d >= germline_min_depth
                    and freq >= p.germline_freq_threshold):   # fixed/germline -> not SHM
                germline_removed += 1
                continue
            if model == "sanger":
                # informational: error probability of the best supporting base
                qs = [q for st in (1, -1) for q in alt_quals.get((rpos, alt, st), ())]
                pv = 10 ** (-max(qs) / 10.0) if qs else 1.0
            elif model == "geneious":
                pv = _geneious_pvalue(c, d, alt_qsum[rpos][alt] / c)
            elif model == "strand":
                pv = _strand_pvalue(rpos, alt, d, s_counts, s_depth, s_lambda,
                                    null_rates, p)
            elif model == "qcal":
                pv = _qcal_pvalue(rpos, refbase, alt, alt_quals, s_depth, null_rates, p)
            elif null_rates is not None and (rpos, alt) in null_rates:
                # empirical per-position error measured from control reads
                q = null_rates[(rpos, alt)]
                if p.null_rate_floor:
                    lam = err_lambda[rpos].get(alt, 0.0) or (d * p.per_base_error / 3.0)
                    q = max(q, lam / d)
                pv = _binom_sf(c, d, q)
            elif model == "pbinom":
                # LoFreq: every read at the column could have been miscalled
                # as `alt`, including the reads that show it
                pv = _pbinom_sf(c, [e for _b, e in col_probs[rpos]])
            else:
                lam = err_lambda[rpos].get(alt, 0.0) or (d * p.per_base_error / 3.0)
                pv = _poisson_sf(c, lam)
            plus = s_counts.get((rpos, alt, 1), 0)
            minus = s_counts.get((rpos, alt, -1), 0)
            dp, dm = s_depth.get((rpos, 1), 0), s_depth.get((rpos, -1), 0)
            sb = 1.0
            if plus + minus >= p.strand_bias_min_reads and dp and dm:
                from scipy.stats import binomtest
                sb = float(binomtest(plus, plus + minus, dp / (dp + dm)).pvalue)
            candidates.append(Variant(rpos + 1, refbase, alt, c, d, freq, pv, plus, minus, sb,
                                      homopolymer_run(ref, rpos),
                                      ref[max(0, rpos - 2):rpos + 3]))
            pvals.append(pv)

    if model == "sanger":
        keep = [True] * len(pvals)               # quality threshold already applied
    elif p.multiple_testing == "bh":
        if p.bh_family == "all":
            # untested (unseen) alternatives have p = 1: they only enlarge m
            pad = max(0, n_hypotheses - len(pvals))
            keep = _benjamini_hochberg(pvals + [1.0] * pad, p.fdr_alpha)[:len(pvals)]
        elif p.bh_family == "candidates":
            keep = _benjamini_hochberg(pvals, p.fdr_alpha)
        else:
            raise ValueError(f"unknown bh_family {p.bh_family!r}; use 'candidates' or 'all'")
    else:
        keep = [pv <= p.max_pvalue for pv in pvals]
    variants = [v for v, k in zip(candidates, keep, strict=True) if k]

    called = len(variants)
    variants, strand_removed = _strand_filter(variants, p, s_depth)
    hp_flagged = sum(1 for v in variants if v.homopolymer >= p.homopolymer_flag)
    hp_removed = 0
    if p.mask_homopolymers:
        hp_removed = hp_flagged
        variants = [v for v in variants if v.homopolymer < p.homopolymer_flag]

    win_positions = range(lo, hi + 1)
    covered = [pos for pos in win_positions if depth.get(pos, 0) > 0]
    win_depth_sum = sum(depth.get(pos, 0) for pos in win_positions)
    plus_reads = sum(r.strand == 1 for r in aligned_reads if hasattr(r, "strand"))
    stats = {
        "reference": target_name,
        "reference_length": reflen,
        "analysed_positions": len(covered),
        "mean_depth": round(win_depth_sum / max(1, len(covered)), 1),
        "germline_positions_removed": germline_removed,
        "candidates_tested": len(candidates),
        "hypotheses": n_hypotheses if p.bh_family == "all" else len(candidates),
        "multiple_testing": "none" if model == "sanger" else p.multiple_testing,
        "error_model": ("sanger_threshold" if model == "sanger"
                        else "geneious" if model == "geneious"
                        else "strand_betabinomial" if model == "strand"
                        else "quality_calibrated" if model == "qcal"
                        else "empirical" if null_rates else
                        "poisson_binomial" if model == "pbinom" else "read_quality"),
        "passed_test": called,
        "strand_filter": p.strand_filter,
        "strand_filter_removed": strand_removed,
        "homopolymer_flagged": hp_flagged,
        "homopolymer_removed": hp_removed,
        "plus_strand_read_fraction": round(plus_reads / len(aligned_reads), 3)
        if aligned_reads else None,
        "unique_variants": len(variants),
        "min_base_quality": minq,
        "variant_base_quality_hist": dict(sorted(qual_hist.items())),
        "reference_base_quality_hist": dict(sorted(ref_qual_hist.items())),
    }
    if model == "sanger":
        stats["clones"] = n_clones
        stats["germline_filter"] = ("on" if germline_on else "off")
    if model == "strand":
        stats["overdispersion_rho"] = round(null_rates.rho, 6)
    if notes:
        stats["notes"] = notes
    # every candidate the test saw, kept or not: the decision boundary a volcano
    # plot draws. Underscored so write_outputs keeps it out of summary.json.
    stats["_candidates"] = [(round(v.frequency, 6), v.p_value, bool(k))
                            for v, k in zip(candidates, keep, strict=True)]
    return variants, stats


def _strand_pvalue(rpos, alt, d, s_counts, s_depth, s_lambda, null, p) -> float:
    """Test each strand against its own control rate, then combine.

    The rate for a strand comes from the control reads on that strand
    (null.strand_rates), never below what that strand's read qualities alone
    predict when null_rate_floor is set. A strand with no reads contributes
    nothing.
    """
    per = []
    for s in (1, -1):
        n = s_depth.get((rpos, s), 0)
        if n == 0:
            continue
        k = s_counts.get((rpos, alt, s), 0)
        mu = null.strand_rates.get((rpos, alt, s))
        if mu is None:
            mu = null.get((rpos, alt), p.per_base_error / 3.0)
        if p.null_rate_floor:
            lam = s_lambda.get((rpos, alt, s), 0.0) or (n * p.per_base_error / 3.0)
            mu = max(mu, lam / n)
        per.append(_betabinom_sf(k, n, mu, null.rho))
    if not per:
        return 1.0
    if len(per) == 1:
        return per[0]
    if p.strand_combine == "both":
        return max(per)
    return _fisher(per)


def _qcal_pvalue(rpos, refbase, alt, alt_quals, s_depth, null, p) -> float:
    """Quality-calibrated test against control-measured miscall rates.

    Geneious's insight is that a high-quality variant base is more believable
    than a low-quality one; its weakness is that it trusts the basecaller's
    quality scores as stated. Here the trust is measured: in mutation-free
    control reads, how often is this reference base miscalled as `alt`, on this
    strand, with at least this quality (null.qrate)? For the supporting reads
    ranked by quality q(1) >= q(2) >= ..., the j-th test asks how likely j or
    more such miscalls with quality >= q(j) are among this column's reads,
    X ~ Poisson(sum over strands of reads x rate). The smallest of these
    p-values, times the number of tests, is the call's p-value.
    """
    quals = []
    for s in (1, -1):
        quals += [(q, s) for q in alt_quals.get((rpos, alt, s), ())]
    if not quals:
        return 1.0
    quals.sort(reverse=True)
    # Every j is tested, but the Poisson tails are evaluated in one vectorised
    # call. Summing each tail by hand costs O(j), so testing all k reads cost
    # O(k^2); at amplicon depth k reaches the thousands (8,967 at one Ramos
    # site) and that dominated the run. The regularised incomplete gamma gives
    # the same tail in constant time per j.
    lams = []
    for qj, _s in quals:
        lam = 0.0
        for s in (1, -1):
            n = s_depth.get((rpos, s), 0)
            if n:
                lam += n * null.qrate(rpos, refbase, alt, s, qj, p.null_prior_beta)
        lams.append(lam)
    js = numpy.arange(len(lams))          # sf(j-1, lam) = P(X >= j)
    best = float(poisson.sf(js, numpy.asarray(lams)).min())
    return min(1.0, best * len(quals))


def _strand_filter(variants: list[Variant], p: Params, s_depth) -> tuple[list[Variant], int]:
    if p.strand_filter == "off":
        return variants, 0
    kept = []
    for v in variants:
        if p.strand_filter == "both":
            ok = v.plus_reads > 0 and v.minus_reads > 0
        elif p.strand_filter == "bias":
            ok = v.strand_bias_p >= p.strand_bias_alpha
        else:
            raise ValueError(f"unknown strand_filter {p.strand_filter!r}")
        if ok:
            kept.append(v)
    return kept, len(variants) - len(kept)
