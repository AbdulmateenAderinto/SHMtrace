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
from bisect import bisect_left
from collections import defaultdict
from dataclasses import replace

from .align import pair_window
from .config import Params
from .io_utils import read_fasta, read_fastq
from .variants import gate_quality

BASES = ("A", "C", "G", "T")
QMAX = 60
# a site where controls show the alt in at least this share of reads is
# polymorphic there (a haplotype or contaminant), so it is kept out of the
# pooled quality model and judged on its own rate
POLYMORPHIC_IN_CONTROLS = 0.02


def _usable(read, margin: int):
    """The read's aligned columns the pileup counts (align.pair_window)."""
    return pair_window(read, margin)[0]


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
    from .pipeline import prepare_reads
    records = read_fasta(reference)
    target_name = target_name or records[0][0]
    if target_name != records[0][0]:
        # prepare_reads, like the sample route, analyses the first reference record
        raise ValueError(f"the null is fitted on the first reference record ({records[0][0]}), "
                         f"not {target_name}")
    ref = records[0][1]
    lo = (p.analysis_window[0] - 1) if p.analysis_window else 0
    hi = (p.analysis_window[1] - 1) if p.analysis_window else len(ref) - 1

    # Controls carry no SHM, so any co-occurring block in them is a haplotype,
    # paralogue or contaminant: always drop it here, and hand the blocks back so
    # the same blocks can be dropped from the samples (phase_mode="auto").
    # Controls are prepared exactly as samples are (pipeline.prepare_reads: the
    # paralog-allele filter, UMI consensus, ignore_qualities, Sanger trimming),
    # so the null describes the same observational unit (audit F07). Only the
    # phasing policy differs, by design.
    pc = control_params(p)
    blocks: list = []
    per_file = []
    for fq in control_fastqs:
        aligned, s_prep = prepare_reads(fq, reference, contaminants, pc)
        aligned = [a for a in aligned if a.ctg == target_name]
        blocks += [[(pos - 1, alt) for pos, alt in b["sites"]] for b in s_prep.get("blocks", [])]
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
        # ascending, so qrate can bisect instead of counting
        self.site_quals = {k: sorted(v) for k, v in (site_quals or {}).items()}
        self.site_depth = dict(site_depth or {})
        # between-library inflation of the quality-calibrated rates, measured by
        # testing each control against the others (calibrate_qscale)
        self.qscale = float(qscale)
        # control-based FDR (multiple_testing="control"): for each control, the
        # p-values of its candidates when it was called against the other
        # controls at the samples' depth, and the bases it had in the window.
        # {control name: (sorted p-values, bases)}
        self.decoys: dict = {}
        # base-quality recalibration fitted on the same control reads (recal.py), or None
        self.recal = None
        # (pos0, alt) a control library carries as a second allele: an allele fraction
        # between control_het_min_frequency and germline_freq_threshold, which the
        # frequency rule is too high to catch (Params.mask_control_het). Empty unless
        # that is on. variants._candidate_exclusion drops these from the samples'
        # candidates and from the decoys alike.
        self.het_sites: set = set()

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
        site = self.site_quals.get((pos0, alt, strand))
        k = (len(site) - bisect_left(site, q)) if site else 0
        site = (k + prior_trials * g) / (n + prior_trials) if n else g
        return min(1.0, self.qscale * max(g, site))

    @property
    def has_strand(self) -> bool:
        return bool(self.strand_rates)

    def to_json(self) -> dict:
        return {"format": "shmtrace-null-panel", "version": 3, "rho": self.rho,
                "qscale": self.qscale,
                "decoys": [[name, list(pvalues), bases]
                           for name, (pvalues, bases) in sorted(self.decoys.items())],
                "meta": self.meta,
                "qtail": [[k[0], k[1], k[2], v] for k, v in sorted(self.qtail.items())],
                "site_quals": [[k[0], k[1], k[2], v]
                               for k, v in sorted(self.site_quals.items())],
                "site_depth": [[k[0], k[1], v] for k, v in sorted(self.site_depth.items())],
                "rates": [[k[0], k[1], v] for k, v in sorted(self.items())],
                "strand_rates": [[k[0], k[1], k[2], v]
                                 for k, v in sorted(self.strand_rates.items())],
                "het_sites": [[pos, alt] for pos, alt in sorted(self.het_sites)],
                "recal": self.recal.to_json() if self.recal is not None else None}

    @classmethod
    def from_json(cls, data: dict) -> NullModel:
        if data.get("format") != "shmtrace-null-panel":
            raise ValueError("not a shmtrace null panel")
        null = cls({(int(a), b): float(v) for a, b, v in data["rates"]},
                   {(int(a), b, int(s)): float(v) for a, b, s, v in data.get("strand_rates", [])},
                   data.get("rho", 0.0), data.get("meta"),
                   {(a, b, int(s)): v for a, b, s, v in data.get("qtail", [])},
                   {(int(a), b, int(s)): v for a, b, s, v in data.get("site_quals", [])},
                   {(int(a), int(s)): v for a, s, v in data.get("site_depth", [])},
                   data.get("qscale", 1.0))._with_recal(data.get("recal"))
        # Version 1/2 panels remain readable as rate-only panels. The caller
        # reports the BH fallback if control thresholding has no saved decoys.
        null.decoys = {name: (sorted(float(x) for x in pvalues), int(bases))
                       for name, pvalues, bases in data.get("decoys", [])}
        null.het_sites = {(int(pos), alt) for pos, alt in data.get("het_sites", [])}
        return null

    def _with_recal(self, data):
        from .recal import Recal
        self.recal = Recal.from_json(data)
        return self


def save_panel(null: NullModel, path: str, **meta) -> None:
    """Write rates, calibration and decoys for reproducible reuse.

    Restoring the state does not establish calibration across sequencing runs.
    Use matched same-run controls; chemistry, basecaller, reference, window and
    preparation settings must also be compatible. Compatibility is not enforced
    by this serializer.
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


# The inflation factors calibrate_qscale may search, smallest first. Params.qscale_max
# caps which of them it actually tries, and defaults to 50.0 -- the grid's whole extent
# until October 2026, and the published behaviour. The steps past 50 exist because 50 was
# both the last step and the ceiling: a project that needed more was given 50 silently and
# called almost nothing. Raising qscale_max lets such a project say how far off it is, at
# the cost of being stricter than the published runs were; an ONT library analysed well
# below its depth loses a few true calls that way, so it is not the default.
QSCALE_GRID = (1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0)
# above this the quality model is being stretched far beyond anything the nanopore
# calibration needed, and the calls should be treated as unreliable rather than reported
# plainly
QSCALE_SUSPECT = 50.0
DEFAULT_QSCALE = 10.0   # used with a single control file; see calibrate_qscale


def window_bases(aligned, lo: int, hi: int, minq: int = 0, margin: int = 0) -> int:
    """Bases read inside the window (the denominator the decoy counts are scaled by).

    `minq` is the caller's base-quality gate, so the decoys' total counts the same
    bases as the sample's (variants.call_variants, win_bases). `margin` likewise
    leaves out the read-end columns the pileup leaves out."""
    return sum(1 for r in aligned for rpos, base, q in _usable(r, margin)
               if base in BASES and q >= minq and lo <= rpos <= hi)


def _fold(i, aligned_per_file, ref, p, lo, hi, depth, null=None):
    """Control i held out: its reads (at the samples' depth, so its candidates, the
    decoys of multiple_testing="control", are what a sample would show) and the null
    fitted from the other controls."""
    test = aligned_per_file[i]
    others = [a for j, a in enumerate(aligned_per_file) if j != i]
    if depth:
        from .reads_filter import downsample
        test = downsample(test, replace(p, downsample_to=depth))[0]
    return test, null if null is not None else fit_aligned(others, ref, p, lo, hi, calibrate=False)[0]


def _held_out_pvalues(test, ref_fasta, ref_name, pp, null, grid=QSCALE_GRID):
    from .variants import qcal_pvalues_by_qscale
    return qcal_pvalues_by_qscale(test, ref_fasta, ref_name, pp, null, grid)


def calibrate_qscale(aligned_per_file: list[list], ref: str, ref_name: str, p: Params,
                     lo: int, hi: int, alpha: float = 0.01, depth: int | None = None,
                     names: list[str] | None = None, fold_nulls: list | None = None) -> tuple[float, dict]:
    """Inflate the quality-calibrated rates until controls look like controls.

    Pooled rates are fitted from some libraries and applied to another, and
    nanopore error varies between libraries more than within one. Each control
    is tested against a null fitted from the other controls; the smallest
    factor in QSCALE_GRID at which no more than `alpha` of its candidates reach
    p < alpha is used. On the JH4 Aicda-/- libraries this lands at 10.

    The search is capped at `p.qscale_max` (default 50.0, the published ceiling).
    A factor above QSCALE_SUSPECT, and the cap in particular, means the controls
    could not be made to look like controls by any inflation the search tried;
    `qscale_warning` says so, and the project carries it as a note. The calls from
    such a run are not trustworthy, whichever way the search ended.
    """
    if len(aligned_per_file) < 2:
        return DEFAULT_QSCALE, {"qscale_source": "default (one control file)"}
    from .parallel import pmap
    # the folds are fitted here, not on worker processes: a fit takes seconds, but sending a
    # deep library's reads to another process costs about as long (fold_nulls, when given,
    # are the same fits made by the caller and are reused)
    folds = [_fold(i, aligned_per_file, ref, p, lo, hi, depth,
                   fold_nulls[i] if fold_nulls is not None else None) for i in range(len(aligned_per_file))]
    names = list(names) if names else [f"control_{i + 1}" for i in range(len(folds))]
    pp = replace(p, error_model="qcal", multiple_testing="none", max_pvalue=1.0,
                 strand_filter="off", mask_homopolymers=False)
    trace = {}
    grid = [c for c in QSCALE_GRID if c <= p.qscale_max] or [QSCALE_GRID[0]]
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    tmp.write(f">{ref_name}\n{ref}\n")
    tmp.close()
    try:
        decoys = {}
        bases = [window_bases(test, lo, hi, gate_quality(p), p.read_end_margin_bp)
                 for test, _null in folds]
        # each held-out control is tallied once and tested at every margin (on its own
        # process, parallel.py); the margins are then read in order, as before
        grids = pmap(_held_out_pvalues, [(test, tmp.name, ref_name, pp, null, grid)
                                         for test, null in folds], p.jobs)
        for c in grid:
            m = small = 0
            decoys = {}
            pvs = [g[c] for g in grids]
            for name, pv, nb in zip(names, pvs, bases, strict=True):
                m += len(pv)
                small += sum(x < alpha for x in pv)
                decoys[name] = (pv, nb)
            trace[c] = (small, m)
            if small <= alpha * m:
                out = {"qscale_source": "leave-one-control-out",
                       "qscale_trace": {str(k): v for k, v in trace.items()},
                       "decoys": decoys}
                if c > QSCALE_SUSPECT:
                    out["qscale_warning"] = (
                        f"the quality-calibrated error rates had to be inflated {c:g}-fold before "
                        "the controls stopped looking like samples (the nanopore libraries this was "
                        "calibrated on need 10-fold). The reported base qualities do not describe "
                        "this data's errors; treat the calls as unreliable and compare libraries "
                        "only with the same controls")
                return c, out
    finally:
        os.remove(tmp.name)
    return grid[-1], {"qscale_source": "leave-one-control-out (grid maximum)",
                      "qscale_trace": {str(k): v for k, v in trace.items()},
                      "qscale_warning": (
                          f"no inflation up to {grid[-1]:g}-fold made the controls look like controls, "
                          "so the rates are inflated to that maximum and little or nothing can be "
                          "called. The quality model does not fit this data"
                          + ("; raise qscale_max to see how far off it is"
                             if grid[-1] < QSCALE_GRID[-1] else "")),
                      "decoys": decoys}


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
    het_sites: set = set()
    stats = {"control_files": len(aligned_per_file), "control_reads_used": 0,
             "germline_sites_excluded": 0}
    recal = None
    if p.recalibrate and p.error_model == "qcal":
        from .recal import fit as fit_recal
        recal = fit_recal(aligned_per_file, ref, lo, hi, p.germline_freq_threshold)
        stats["recalibration"] = recal.meta
    # Bases below the caller's quality gate never reach a sample's pileup, so they
    # are left out of the null's counts and depths too: every rate, and every
    # quality tail, is conditional on a base passing the gate (audit F10).
    minq = gate_quality(p)
    for aligned in aligned_per_file:
        stats["control_reads_used"] += len(aligned)
        counts: dict = defaultdict(lambda: defaultdict(int))
        depth: dict = defaultdict(int)
        sc: dict = defaultdict(int)
        sd: dict = defaultdict(int)
        quals: dict = defaultdict(list)          # (pos, base, strand) -> qualities
        ref_u = ref.upper()
        for r in aligned:
            s = getattr(r, "strand", 1)
            eq = recal.quals(r.pairs) if recal is not None else None
            cols, off = pair_window(r, p.read_end_margin_bp)
            for i, (rpos, base, q) in enumerate(cols, off):
                # only the window is read below, so bases outside it are skipped
                if base in BASES and q >= minq and lo <= rpos <= hi:
                    counts[rpos][base] += 1
                    depth[rpos] += 1
                    sc[(rpos, base, s)] += 1
                    sd[(rpos, s)] += 1
                    if base != ref_u[rpos]:
                        quals[(rpos, base, s)].append(min(QMAX, eq[i] if eq is not None else q))
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
                    # A fixed difference from the reference, not error. With a declared ancestral
                    # reference it is kept, as a position-specific rate: the controls are the
                    # un-diversified baseline, so a difference they carry is baseline material and the
                    # case must not be called for carrying it too.
                    if not p.ancestral_reference:
                        stats["germline_sites_excluded"] += 1
                        continue
                elif (p.mask_control_het and not p.ancestral_reference
                      and d >= p.control_het_min_depth
                      and k / d >= p.control_het_min_frequency):
                    # Below the germline cut but far above error: a second allele this
                    # control carries. One control is enough to establish it, so that a
                    # variant private to one library is caught too. The site still gets
                    # its rate fitted below, so the background model is unchanged; what
                    # the mask does is stop the caller testing there at all.
                    # Not under a declared ancestral reference, where a change driven
                    # high is the point of the experiment.
                    het_sites.add((pos, alt))
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
    null.recal = recal
    null.het_sites = het_sites
    if p.mask_control_het:
        stats["control_het_sites"] = len(het_sites)
        # (1-based position, allele) for the report
        stats["control_het_alleles"] = sorted((pos + 1, alt) for pos, alt in het_sites)
    if calibrate:
        null.qscale, s_cal = calibrate_qscale(aligned_per_file, ref, ref_name, p, lo, hi,
                                              depth=p.downsample_to)
        null.decoys = s_cal.pop("decoys", {})
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
