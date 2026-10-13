"""Group comparisons with the animal as the unit.

Reads and clones from one mouse are not independent: a clone's reads share its
mutations, and mice differ from each other far more than reads within a mouse do.
Pooling every read of a group and testing reads against reads treats thousands of
correlated observations as independent and gives p-values that are far too small
(pseudoreplication). Every test here uses one value per animal.

For a count-based rate (mutations per base sequenced), each animal contributes
(k mutations, n bases). Per group:

  negative binomial
                  counts k ~ NB(mean = rate x n, variance mean + phi mean^2): the animals'
                  rates vary around their group's rate with one dispersion phi shared by
                  the groups (a constant coefficient of variation, as on the log scale),
                  fitted under the full model as DESeq2 and edgeR do. The group effect
                  is a likelihood-ratio test with the p-value by parametric bootstrap
                  (small groups make the chi-square approximation unsafe).
  permutation     exact (or 20,000-draw) permutation of animals between the groups on
                  log10 rates (an animal with no mutations counts half an event): no distributional assumption at all; with 3 vs 3 the
                  smallest two-sided p is 0.1, which is the honest limit of 6 animals.
  Welch on log10  the manuscript's test.
  fold            the ratio of the groups' fitted rates, and the ratio of geometric
                  means with a 95% interval on the log scale (Welch-Satterthwaite).
  pooled reads    the two-proportion test on summed counts, reported only to show how
                  much smaller pseudoreplication would make the p-value.
"""
from __future__ import annotations

import itertools
import math
import random

import numpy as np
from scipy import optimize, stats
from scipy.special import gammaln


def _nb_ll(k, mu, phi):
    """Negative binomial log-likelihood: mean mu, variance mu + phi mu^2."""
    k, mu = np.asarray(k, float), np.maximum(np.asarray(mu, float), 1e-300)
    phi = max(phi, 1e-10)
    r = 1 / phi
    return float(np.sum(gammaln(k + r) - gammaln(r) - gammaln(k + 1)
                        + r * np.log(r / (r + mu)) + k * np.log(mu / (r + mu))))


def _group_rate(k, n, phi):
    """MLE of a group's rate given phi (1-D, on the log scale)."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    if k.sum() == 0:
        return 0.0
    m0 = math.log(k.sum() / n.sum())
    res = optimize.minimize_scalar(lambda b: -_nb_ll(k, n * math.exp(b), phi),
                                   bounds=(m0 - 10, m0 + 10), method="bounded")
    return math.exp(res.x)


def fit_nb(groups: list[tuple[list, list]], common: bool = False,
           phi: float | None = None) -> tuple[list[float], float, float]:
    """Rates per group (one shared rate if `common`) and a shared dispersion phi
    (estimated unless given): returns (rates, phi, loglik)."""
    ks = [np.asarray(k, float) for k, _ in groups]
    ns = [np.asarray(n, float) for _, n in groups]
    if common:
        ks, ns = [np.concatenate(ks)], [np.concatenate(ns)]

    def ll_at(ph):
        rates = [_group_rate(k, n, ph) for k, n in zip(ks, ns, strict=True)]
        return sum(_nb_ll(k, n * r, ph) for k, n, r in zip(ks, ns, rates, strict=True)), rates
    if phi is None:
        res = optimize.minimize_scalar(lambda lp: -ll_at(math.exp(lp))[0],
                                       bounds=(-12, 4), method="bounded")
        phi = math.exp(res.x)
    ll, rates = ll_at(phi)
    if common:
        rates = rates * len(groups)
    return rates, phi, ll


def _log_rates(k, n):
    """log10 rates; an animal with no mutations counts half an event (continuity
    correction), so its log is finite and set by its own depth, not by other animals."""
    return [math.log10(max(ki, 0.5) / ni) for ki, ni in zip(k, n, strict=True) if ni]


def permutation_p(a: list[float], b: list[float], max_exact: int = 20000,
                  seed: int = 1) -> float:
    """Two-sided permutation p for the difference in means of a and b."""
    allv = a + b
    obs = abs(np.mean(a) - np.mean(b))
    na = len(a)
    combos = math.comb(len(allv), na)
    if combos <= max_exact:
        hits = total = 0
        for idx in itertools.combinations(range(len(allv)), na):
            s = set(idx)
            x = [allv[i] for i in idx]
            y = [allv[i] for i in range(len(allv)) if i not in s]
            total += 1
            hits += abs(np.mean(x) - np.mean(y)) >= obs - 1e-12
        return float(hits / total)
    rng = random.Random(seed)
    hits = 0
    for _ in range(max_exact):
        rng.shuffle(allv)
        hits += abs(np.mean(allv[:na]) - np.mean(allv[na:])) >= obs - 1e-12
    return float((hits + 1) / (max_exact + 1))


def compare(ctrl: list[tuple[int, int]], case: list[tuple[int, int]], boot: int = 499,
            seed: int = 1) -> dict:
    """Compare two groups of animals, each (k mutations, n bases)."""
    kc, nc = [x[0] for x in ctrl], [x[1] for x in ctrl]
    kx, nx = [x[0] for x in case], [x[1] for x in case]
    out: dict = {"animals_case": len(case), "animals_control": len(ctrl)}
    (mx, mc), phi, ll1 = fit_nb([(kx, nx), (kc, nc)])
    out.update(rate_case=mx, rate_control=mc, dispersion=phi,
               fold=(mx / mc) if mc > 0 else (mx / (0.5 / sum(nc)) if mx > 0 else 1.0),
               # a control group with no mutations at all: the fold is against half an event
               # over all its bases, so it is a lower bound
               fold_is_lower_bound=bool(mc <= 0 < mx))
    lc, lx = _log_rates(kc, nc), _log_rates(kx, nx)
    if len(lc) >= 2 and len(lx) >= 2:
        t = stats.ttest_ind(lx, lc, equal_var=False)
        out["p_welch_log"] = float(t.pvalue)
        # fold interval on the log scale (Welch-Satterthwaite df)
        va, vb = np.var(lx, ddof=1) / len(lx), np.var(lc, ddof=1) / len(lc)
        se = math.sqrt(va + vb)
        df = (va + vb) ** 2 / ((va ** 2 / (len(lx) - 1) if va else 0) + (vb ** 2 / (len(lc) - 1) if vb else 0)
                               or 1e-12)
        d = np.mean(lx) - np.mean(lc)
        q = stats.t.ppf(0.975, max(df, 1))
        out["fold_geometric"] = 10 ** d
        out["fold_ci"] = [10 ** (d - q * se), 10 ** (d + q * se)]
    out["p_permutation"] = permutation_p(lx, lc) if lc and lx else None
    if len(ctrl) >= 2 and len(case) >= 2:
        # negative binomial with one dispersion for both groups, fitted under the full
        # model (as DESeq2 and edgeR do); the p-value by parametric bootstrap under a
        # common rate with that dispersion, because 2-10 animals per group are too few
        # for the chi-square approximation
        (m0, _), _p0, ll0 = fit_nb([(kx, nx), (kc, nc)], common=True, phi=phi)
        stat = max(0.0, 2 * (ll1 - ll0))
        rng = np.random.default_rng(seed)
        hits = 0
        allk_n = np.asarray(list(nx) + list(nc), float)
        for _ in range(boot):
            lam = allk_n * m0
            sim = rng.negative_binomial(1 / phi, (1 / phi) / ((1 / phi) + lam)) if phi > 1e-10 \
                else rng.poisson(lam)
            sx, sc = list(sim[:len(nx)]), list(sim[len(nx):])
            _r, ph_b, l1 = fit_nb([(sx, nx), (sc, nc)])
            _r0, _p, l0 = fit_nb([(sx, nx), (sc, nc)], common=True, phi=ph_b)
            hits += 2 * (l1 - l0) >= stat - 1e-9
        out["lrt_negbinomial"] = stat
        out["p_negbinomial"] = (hits + 1) / (boot + 1)
    # what pooling reads would claim (pseudoreplication), for contrast only
    K, N = sum(kx) + sum(kc), sum(nx) + sum(nc)
    if 0 < K < N:
        p_pool = K / N
        se = math.sqrt(p_pool * (1 - p_pool) * (1 / sum(nx) + 1 / sum(nc)))
        z = (sum(kx) / sum(nx) - sum(kc) / sum(nc)) / se if se else 0.0
        out["p_pooled_reads"] = float(2 * stats.norm.sf(abs(z)))
    return out


METRICS = {"mutation_frequency": ("mutations_total", "frequency_bases"),
           "unique_density": ("unique_mutations", "density_bases")}


def group_tests(samples: list[dict], boot: int = 499) -> list[dict]:
    """Every non-control group against the control group, per count-based metric."""
    groups: dict = {}
    for r in samples:
        groups.setdefault(r["group"], []).append(r)
    ctrl = next((g for g, rs in groups.items() if all(r["control"] for r in rs)), None)
    if ctrl is None:
        return []
    rows = []
    for g, rs in groups.items():
        if g == ctrl:
            continue
        for metric, (kk, nk) in METRICS.items():
            case = [(int(r[kk]), int(r[nk])) for r in rs if r.get(nk)]
            con = [(int(r[kk]), int(r[nk])) for r in groups[ctrl] if r.get(nk)]
            if not case or not con:
                continue
            c = compare(con, case, boot=boot)
            rows.append({"metric": metric, "group": g, "control": ctrl, **{
                k: (round(v, 6) if isinstance(v, float) else v) for k, v in c.items()}})
    return rows
