"""Empirical per-position error model fitted from mutation-free control reads.

The quality-derived null in `variants.py` assumes sequencing errors are
independent across reads: lambda = sum(10^(-Q/10))/3 over the column. Nanopore
error is not independent. It is context-driven, so the same position errs in
many reads at once, and measured against Aicda-/- controls the real rate runs
about 3x higher than the quality model predicts
(docs/calibration/FINDINGS.md).

Given control reads that should carry no somatic mutations (an AID-deficient
animal, a germline sample, an unmutated cell line), we can measure the error
rate at each position and alt base directly and test candidates against that
instead.

Rates are shrunk toward a weak Beta prior so that a position seen at low depth
in the control does not produce an overconfident rate of exactly zero, which
would make every observation there look infinitely significant.
"""
from __future__ import annotations

import gzip
import os
import re
import tempfile
from collections import defaultdict
from dataclasses import replace

from .align import align_reads, build_index
from .config import Params
from .io_utils import read_fasta, read_fastq
from .phasing import partition
from .reads_filter import prefilter

BASES = ("A", "C", "G", "T")
QMAX = 60
# a site where controls show the alt in at least this share of reads is
# polymorphic there (a haplotype or contaminant), so it is kept out of the
# pooled quality model and judged on its own rate
POLYMORPHIC_IN_CONTROLS = 0.02


def _first_header(fastq: str) -> str | None:
    from .io_utils import is_alignment_file, sniff_format
    if os.path.isdir(fastq) or not os.path.exists(fastq) or is_alignment_file(fastq) \
            or sniff_format(fastq) != "fastq":
        # directory, R1+R2 pair spec, FASTA, trace, BAM: the first read's name
        try:
            return next(read_fastq(fastq))[0]
        except Exception:
            return None
    try:
        opener = gzip.open if fastq.endswith(".gz") else open
        with opener(fastq, "rt") as fh:
            line = fh.readline()
        return line[1:].strip() if line.startswith("@") else None
    except Exception:
        return None


def run_token(fastq: str) -> str | None:
    """Best-effort sequencing-run identifier from the first read header.

    Tries, in order: the `runid=` tag that MinKNOW and Dorado write into FASTQ
    headers, then the Plasmidsaurus read-name convention `<n>_<RUN>_<sample>`.
    Used only to warn when a null is fitted on a different run from the sample
    it is applied to; nanopore error is run-specific, so a cross-run null is not
    trustworthy. When neither is present, pass the run explicitly (CLI
    --run-id / --null-run-id, or the `run` column of a sample sheet).
    """
    header = _first_header(fastq)
    if not header:
        return None
    m = re.search(r"(?:^|\s)runid=([A-Za-z0-9_-]+)", header)
    if m:
        return m.group(1)
    parts = header.split()[0].split("_")
    return parts[1] if len(parts) >= 3 else None


def check_same_run(control_fastqs: list[str], target_fastq: str | None,
                   sample_run: str | None = None,
                   control_runs: list[str] | None = None) -> tuple[bool, str]:
    """Warn if the null's controls do not come from the target's run.

    Explicit run identifiers win over anything parsed from read headers.
    """
    tgt = sample_run or (run_token(target_fastq) if target_fastq else None)
    ctrl = set(control_runs) if control_runs else {run_token(f) for f in control_fastqs}
    ctrl.discard(None)
    if tgt is None or not ctrl:
        return True, "run identifiers not recognised; cannot verify"
    if ctrl == {tgt}:
        return True, f"same run ({tgt})"
    return False, (f"control run(s) {sorted(ctrl)} differ from sample run {tgt}. "
                   "Nanopore error is run-specific, so this null may not describe "
                   "this sample. Prefer a control sequenced on the same flow cell.")


def fit(control_fastqs: list[str], reference: str, contaminants: list[str],
        p: Params, target_name: str | None = None) -> tuple[dict, dict]:
    """Measure per-(position, alt) error from control reads.

    Returns (rates, stats). `rates` maps (0-based position, alt base) to a
    shrunk error probability, suitable as the null for a binomial test.
    """
    refs = dict(read_fasta(reference))
    target_name = target_name or read_fasta(reference)[0][0]
    ref = refs[target_name]
    lo = (p.analysis_window[0] - 1) if p.analysis_window else 0
    hi = (p.analysis_window[1] - 1) if p.analysis_window else len(ref) - 1

    # Controls carry no SHM, so any co-occurring block in them is a haplotype,
    # paralogue or contaminant: always drop it here, and hand the blocks back so
    # the same blocks can be dropped from the samples (phase_mode="auto").
    pc = control_params(p)
    blocks: list = []
    per_file = []
    aligner, tnames, cnames = build_index([reference], contaminants, p.minimap2_preset)
    for fq in control_fastqs:
        reads = list(read_fastq(fq))
        kept, _ = prefilter(reads, p)
        aligned = [a for a in align_reads(kept, aligner, tnames, cnames, p) if a.ctg == target_name]
        aligned, s_ph = partition(aligned, ref, pc, lo, hi)
        blocks += [[(pos - 1, alt) for pos, alt in b["sites"]] for b in s_ph.get("blocks", [])]
        per_file.append(aligned)
    rates, stats = fit_aligned(per_file, ref, p, lo, hi, ref_name=target_name)
    stats["control_files"] = len(control_fastqs)
    stats["haplotype_blocks_found"] = len(blocks)
    stats["_blocks"] = blocks
    return rates, stats


def control_params(p: Params) -> Params:
    """Phasing settings for mutation-free controls: drop every block."""
    return replace(p, phase_mode="drop") if p.phase_reads and p.phase_mode != "off" else p


class NullModel(dict):
    """Per-(position, alt) error rates from controls, plus a strand-aware model.

    Behaves as the plain {(pos0, alt): rate} dict the binomial test has always
    used. On top of that it carries:

    * `strand_rates`: {(pos0, alt, strand): rate}, strand +1 or -1. Nanopore
      errors are often specific to the strand a read was sequenced on (the
      basecaller sees a different k-mer context), so an error that is common on
      one strand can be rare on the other.
    * `rho`: beta-binomial overdispersion. Errors at a position vary between
      libraries more than binomial sampling allows; rho is fitted from the
      spread between control files (0 means binomial).
    * `meta`: provenance (control files, run, basecaller) for saved panels.
    """

    def __init__(self, rates=None, strand_rates=None, rho=0.0, meta=None,
                 qtail=None, site_quals=None, site_depth=None, qscale=1.0):
        super().__init__(rates or {})
        self.strand_rates = dict(strand_rates or {})
        self.rho = float(rho)
        self.meta = dict(meta or {})
        # quality-calibrated model (error_model="qcal"):
        # qtail[(ref, alt, strand)][q] = share of control bases at `ref` sites on
        #   that strand miscalled as `alt` with quality >= q, pooled over sites
        # site_quals[(pos0, alt, strand)] = qualities of those miscalls at one site
        # site_depth[(pos0, strand)] = control bases observed there
        self.qtail = {k: list(v) for k, v in (qtail or {}).items()}
        self.site_quals = {k: sorted(v, reverse=True) for k, v in (site_quals or {}).items()}
        self.site_depth = dict(site_depth or {})
        # between-library inflation of the quality-calibrated rates, measured by
        # testing each control against the others (calibrate_qscale)
        self.qscale = float(qscale)

    @property
    def has_qcal(self) -> bool:
        return bool(self.qtail)

    def qrate(self, pos0: int, ref: str, alt: str, strand: int, q: int,
              prior_trials: float) -> float:
        """Rate of `alt` miscalls with quality >= q at this site and strand.

        The pooled rate across all sites is the floor. A site whose controls
        show more high-quality miscalls than that gets its own, shrunk toward the
        pooled rate with `prior_trials` pseudo-observations.
        """
        q = max(0, min(QMAX, int(q)))
        tail = self.qtail.get((ref, alt, strand)) or self.qtail.get((ref, alt, -strand))
        g = tail[q] if tail else 0.0
        n = self.site_depth.get((pos0, strand), 0)
        k = sum(1 for x in self.site_quals.get((pos0, alt, strand), ()) if x >= q)
        site = (k + prior_trials * g) / (n + prior_trials) if n else g
        return min(1.0, self.qscale * max(g, site))

    @property
    def has_strand(self) -> bool:
        return bool(self.strand_rates)

    def to_json(self) -> dict:
        return {"format": "nanoshm-null-panel", "version": 2, "rho": self.rho,
                "qscale": self.qscale,
                "meta": self.meta,
                "qtail": [[k[0], k[1], k[2], v] for k, v in sorted(self.qtail.items())],
                "site_quals": [[k[0], k[1], k[2], v]
                               for k, v in sorted(self.site_quals.items())],
                "site_depth": [[k[0], k[1], v] for k, v in sorted(self.site_depth.items())],
                "rates": [[k[0], k[1], v] for k, v in sorted(self.items())],
                "strand_rates": [[k[0], k[1], k[2], v]
                                 for k, v in sorted(self.strand_rates.items())]}

    @classmethod
    def from_json(cls, data: dict) -> NullModel:
        if data.get("format") != "nanoshm-null-panel":
            raise ValueError("not a nanoshm null panel")
        return cls({(int(a), b): float(v) for a, b, v in data["rates"]},
                   {(int(a), b, int(s)): float(v) for a, b, s, v in data.get("strand_rates", [])},
                   data.get("rho", 0.0), data.get("meta"),
                   {(a, b, int(s)): v for a, b, s, v in data.get("qtail", [])},
                   {(int(a), b, int(s)): v for a, b, s, v in data.get("site_quals", [])},
                   {(int(a), int(s)): v for a, s, v in data.get("site_depth", [])},
                   data.get("qscale", 1.0))


def save_panel(null: NullModel, path: str, **meta) -> None:
    """Write a fitted null to JSON so it can be reused on later runs.

    A panel describes one chemistry, basecaller model and amplicon. Reuse it only
    for the same combination, and prefer a same-run control when there is one.
    """
    import json
    null.meta.update({k: v for k, v in meta.items() if v is not None})
    with open(path, "w") as fh:
        json.dump(null.to_json(), fh)


def load_panel(path: str) -> NullModel:
    import json
    with open(path) as fh:
        return NullModel.from_json(json.load(fh))


def estimate_rho(per_file: list[dict], keys) -> float:
    """Method-of-moments beta-binomial overdispersion, pooled over sites.

    For each site, the spread of the per-file error counts around the pooled
    rate is compared with what binomial sampling alone would give:
    Var(k_i) = n_i mu (1 - mu) (1 + (n_i - 1) rho). Needs two or more files.
    """
    if len(per_file) < 2:
        return 0.0
    num = den = 0.0
    for key in keys:
        obs = [(f[key][0], f[key][1]) for f in per_file if key in f and f[key][1] > 0]
        if len(obs) < 2:
            continue
        n_tot = sum(n for _k, n in obs)
        mu = sum(k for k, _n in obs) / n_tot
        if mu <= 0 or mu >= 1:
            continue
        v = mu * (1 - mu)
        for k, n in obs:
            num += (k - n * mu) ** 2 - n * v
            den += n * (n - 1) * v
    if den <= 0:
        return 0.0
    return float(min(0.1, max(0.0, num / den)))


QSCALE_GRID = (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0)
DEFAULT_QSCALE = 10.0   # used with a single control file; see calibrate_qscale


def calibrate_qscale(aligned_per_file: list[list], ref: str, ref_name: str, p: Params,
                     lo: int, hi: int, alpha: float = 0.01) -> tuple[float, dict]:
    """Inflate the quality-calibrated rates until controls look like controls.

    Pooled rates are fitted from some libraries and applied to another, and
    nanopore error varies between libraries more than within one. Each control
    is tested against a null fitted from the other controls; the smallest
    factor in QSCALE_GRID at which no more than `alpha` of its candidates reach
    p < alpha is used. On the JH4 Aicda-/- libraries this lands at 10.
    """
    from .variants import call_variants
    if len(aligned_per_file) < 2:
        return DEFAULT_QSCALE, {"qscale_source": "default (one control file)"}
    folds = []
    for i, test in enumerate(aligned_per_file):
        others = [a for j, a in enumerate(aligned_per_file) if j != i]
        folds.append((test, fit_aligned(others, ref, p, lo, hi, calibrate=False)[0]))
    pp = replace(p, error_model="qcal", multiple_testing="none", max_pvalue=1.0,
                 strand_filter="off", mask_homopolymers=False)
    trace = {}
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    tmp.write(f">{ref_name}\n{ref}\n")
    tmp.close()
    try:
        for c in QSCALE_GRID:
            m = small = 0
            for test, null in folds:
                null.qscale = c
                vs, _st = call_variants(test, tmp.name, ref_name, pp, null)
                m += len(vs)
                small += sum(v.p_value < alpha for v in vs)
            trace[c] = (small, m)
            if small <= alpha * m:
                return c, {"qscale_source": "leave-one-control-out",
                           "qscale_trace": {str(k): v for k, v in trace.items()}}
    finally:
        os.remove(tmp.name)
    return QSCALE_GRID[-1], {"qscale_source": "leave-one-control-out (grid maximum)",
                             "qscale_trace": {str(k): v for k, v in trace.items()}}


def fit_aligned(aligned_per_file: list[list], ref: str, p: Params,
                lo: int, hi: int, calibrate: bool = True,
                ref_name: str = "ref") -> tuple[NullModel, dict]:
    """Fit per-(position, alt) and per-strand rates from aligned, phased control reads."""
    alt_counts: dict = defaultdict(float)
    trials: dict = defaultdict(float)
    s_alt: dict = defaultdict(float)
    s_trials: dict = defaultdict(float)
    per_file_strand: list[dict] = []
    site_quals: dict = defaultdict(list)
    site_depth: dict = defaultdict(int)
    q_hist: dict = defaultdict(lambda: [0] * (QMAX + 1))
    ref_obs: dict = defaultdict(int)
    class_k: dict = defaultdict(int)          # (ref, alt, strand) pooled, non-polymorphic
    class_n: dict = defaultdict(int)
    site_ref: dict = {}
    stats = {"control_files": len(aligned_per_file), "control_reads_used": 0,
             "germline_sites_excluded": 0}
    for aligned in aligned_per_file:
        stats["control_reads_used"] += len(aligned)
        counts: dict = defaultdict(lambda: defaultdict(int))
        depth: dict = defaultdict(int)
        sc: dict = defaultdict(int)
        sd: dict = defaultdict(int)
        quals: dict = defaultdict(list)          # (pos, base, strand) -> qualities
        for r in aligned:
            s = getattr(r, "strand", 1)
            for rpos, base, q in r.pairs:
                if base in BASES:
                    counts[rpos][base] += 1
                    depth[rpos] += 1
                    sc[(rpos, base, s)] += 1
                    sd[(rpos, s)] += 1
                    if lo <= rpos <= hi and base != ref[rpos].upper():
                        quals[(rpos, base, s)].append(min(QMAX, q))
        this_file: dict = {}
        for pos in range(lo, hi + 1):
            d = depth.get(pos, 0)
            if d == 0:
                continue
            refbase = ref[pos].upper()
            if refbase not in BASES:
                continue
            for alt in BASES:
                if alt == refbase:
                    continue
                k = counts[pos].get(alt, 0)
                if k / d >= p.germline_freq_threshold:
                    # a fixed difference from the reference, not error
                    stats["germline_sites_excluded"] += 1
                    continue
                alt_counts[(pos, alt)] += k
                trials[(pos, alt)] += d
                polymorphic = k / d >= POLYMORPHIC_IN_CONTROLS
                for s in (1, -1):
                    n_s = sd.get((pos, s), 0)
                    if n_s == 0:
                        continue
                    k_s = sc.get((pos, alt, s), 0)
                    s_alt[(pos, alt, s)] += k_s
                    s_trials[(pos, alt, s)] += n_s
                    this_file[(pos, alt, s)] = (k_s, n_s)
                    site_ref[(pos, alt, s)] = refbase
                    if not polymorphic:
                        class_k[(refbase, alt, s)] += k_s
                        class_n[(refbase, alt, s)] += n_s
                    qs = quals.get((pos, alt, s), [])
                    site_quals[(pos, alt, s)].extend(qs)
                    if not polymorphic:
                        for q in qs:
                            q_hist[(refbase, alt, s)][q] += 1
            for s in (1, -1):
                n_s = sd.get((pos, s), 0)
                site_depth[(pos, s)] += n_s
                ref_obs[(refbase, s)] += n_s
        per_file_strand.append(this_file)

    a, b = p.null_prior_alpha, p.null_prior_beta
    rates = {key: (alt_counts[key] + a) / (trials[key] + a + b) for key in trials}
    # Strand rates are shrunk toward the pooled rate for that base change and
    # strand, with p.strand_prior_trials pseudo-reads, and never below it. A
    # light prior lets a site the controls show erring often (a haplotype, a
    # systematic miscall) keep its high rate; the floor stops a site that just
    # happened to show no errors from looking error-free.
    B = p.strand_prior_trials
    cls = {key: (class_k[key] + a) / (class_n[key] + a + b) for key in class_n}
    srates = {}
    for key in s_trials:
        g = cls.get((site_ref[key], key[1], key[2]), a / (a + b))
        srates[key] = max(g, (s_alt[key] + B * g) / (s_trials[key] + B))
    rho = estimate_rho(per_file_strand, list(s_trials))
    qtail = {}
    for (rb, alt, s), hist in q_hist.items():
        n = ref_obs.get((rb, s), 0)
        if n == 0:
            continue
        cum, tail = 0, [0.0] * (QMAX + 1)
        for q in range(QMAX, -1, -1):
            cum += hist[q]
            # half a pseudo-count keeps a never-seen quality tail above zero
            tail[q] = (cum + 0.5) / (n + 1.0)
        qtail[(rb, alt, s)] = tail
    for rb in BASES:                     # classes never miscalled in the controls
        for alt in BASES:
            for s in (1, -1):
                if alt != rb and (rb, alt, s) not in qtail and ref_obs.get((rb, s)):
                    qtail[(rb, alt, s)] = [0.5 / (ref_obs[(rb, s)] + 1.0)] * (QMAX + 1)
    null = NullModel(rates, srates, rho, qtail=qtail,
                     site_quals={k: v for k, v in site_quals.items() if v},
                     site_depth=site_depth)
    if calibrate:
        null.qscale, s_cal = calibrate_qscale(aligned_per_file, ref, ref_name, p, lo, hi)
        stats.update(s_cal)
        stats["qscale"] = null.qscale
    stats["sites_modelled"] = len(rates)
    stats["strand_sites_modelled"] = len(srates)
    stats["overdispersion_rho"] = round(rho, 6)
    if rates:
        vals = sorted(rates.values())
        stats["median_rate"] = round(vals[len(vals) // 2], 6)
        stats["max_rate"] = round(vals[-1], 6)
    if srates:
        # how different the two strands are: sites where one strand's rate is
        # at least 3x the other's (only where either strand saw an error)
        lop = 0
        both = 0
        for (pos, alt, s), v in srates.items():
            if s != 1 or (pos, alt, -1) not in srates:
                continue
            w = srates[(pos, alt, -1)]
            if s_alt[(pos, alt, 1)] + s_alt[(pos, alt, -1)] == 0:
                continue
            both += 1
            if max(v, w) >= 3 * min(v, w):
                lop += 1
        stats["strand_specific_error_sites"] = lop
        stats["error_sites_seen"] = both
    return null, stats
