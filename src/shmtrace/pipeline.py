"""End-to-end SHM analysis: FASTQ + reference(+contaminants) -> variants + metrics.

This is the single entry point the CLI and the web backend both call. It is the
open-source replacement for the Geneious "assemble + find variants" step.

    python -m shmtrace.pipeline \
        --fastq reads.fastq --reference ref.fasta \
        --contaminant contam.fasta --preset presets/jh4_intron.yaml \
        --out results/

Everything is deterministic given a seed, so results are reproducible.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import statistics
from dataclasses import replace

from . import __version__, metrics
from .align import (
    align_reads,
    build_index,
    diagnostic_sites,
    drop_contaminant_alleles,
    trim_messy_ends,
)
from .config import Params, load_preset, to_dict
from .io_utils import read_fasta, read_fastq
from .phasing import partition
from .reads_filter import downsample, prefilter
from .variants import call_variants


def prepare_reads(fastq: str, reference: str, contaminants: list[str], p: Params,
                  known_blocks: list | None = None,
                  fates: dict | None = None, detection: dict | None = None) -> tuple[list, dict]:
    """Filter, align and phase one sample. Returns (aligned reads, read stats).

    Split from `analyse` so that a batch can see every sample's usable depth
    before choosing a common depth for all of them.

    Pass `fates` to record what became of every input read, which is what the
    read browser draws. Reads that survive to the end are not recorded here;
    the caller knows them already. Pass `detection` (an empty dict) to receive the
    phasing blocks detected in the returned reads (phasing.partition's `keep`).
    """
    ref_seq = dict(read_fasta(reference))[read_fasta(reference)[0][0]]
    lo, hi = _window(p, ref_seq)
    reads = list(read_fastq(fastq))
    kept, s_pre = prefilter(reads, p, fates)
    if p.ignore_qualities:
        from .io_utils import FASTA_DEFAULT_Q
        q = chr(33 + FASTA_DEFAULT_Q)
        kept = [(n, s, q * len(s)) for n, s, _q in kept]
    aligner, target_names, contam_names = build_index([reference], contaminants,
                                                      p.minimap2_preset)
    aligned = [a for a in align_reads(kept, aligner, target_names, contam_names, p, fates)
               if a.ctg in target_names]
    if contaminants and p.contaminant_min_alleles > 0:
        sites = diagnostic_sites(reference, contaminants, p.minimap2_preset)
        aligned, n_drop = drop_contaminant_alleles(aligned, sites, p.contaminant_min_alleles,
                                                   p.contaminant_germline_frac, fates)
        s_pre["dropped_contaminant_alleles"] = n_drop
    if p.umi != "off":
        # one consensus read per molecule (umi.py); everything after counts molecules
        from .umi import consensus_reads
        aligned, s_umi = consensus_reads(aligned, ref_seq, p, fates)
        s_pre.update(s_umi)
    if p.platform == "sanger" and p.sanger_clean_window:
        # cut each clone back to where its sequence is clean (Params docs)
        trimmed = [trim_messy_ends(a, ref_seq, p.sanger_clean_window,
                                   p.sanger_clean_max_mismatch) for a in aligned]
        s_pre["clones_without_clean_stretch"] = sum(t is None for t in trimmed)
        if fates is not None:
            for a, t in zip(aligned, trimmed, strict=True):
                if t is None:
                    fates[a.qname] = {"fate": "no clean stretch"}
        aligned = [t for t in trimmed if t is not None]
    s_pre["aligned_before_phasing"] = len(aligned)
    # Remove a co-occurring minor haplotype before depth standardisation, so the
    # retained depth is all target population. Phasing before downsampling also
    # gives the site-detection step the full read set to work with.
    before = {a.qname for a in aligned} if fates is not None else set()
    aligned, s_phase = partition(aligned, ref_seq, p, lo, hi, known_blocks, keep=detection)
    if fates is not None:
        survived = {a.qname for a in aligned}
        for a in _dropped_by_phasing(before, survived):
            fates[a] = {"fate": "minor haplotype"}
    return aligned, {**s_pre, **s_phase}


def _dropped_by_phasing(before: set, after: set) -> set:
    return before - after


def analyse(aligned: list, reference: str, p: Params, read_stats: dict | None = None,
            null_rates: dict | None = None, depth: int | None = None,
            germline: set | None = None, drawn: dict | None = None,
            germline_source: str | None = None) -> dict:
    """Downsample, call variants and summarise already-prepared reads.

    `depth` overrides p.downsample_to, e.g. a batch's common depth. If fewer
    reads are available than the target the sample is analysed at what it has
    and flagged: unique-variant density saturates with depth, so densities from
    different depths are not comparable.

    `germline` is a set of (0-based position, alt) known to be germline from
    elsewhere, removed after calling. A Sanger sample of one or two clones
    cannot tell its own germline differences from mutations, so a project
    supplies the ones its other clones show (project.sanger_germline).

    `drawn` is the record downsample() returned when the caller has already drawn
    `aligned` to the target depth (project.py draws before sending reads to a worker,
    so only the analysed reads are pickled); the reads are then used as given.
    """
    target_name = read_fasta(reference)[0][0]
    target = depth if depth is not None else p.downsample_to
    pd = replace(p, downsample_to=target)
    if drawn is None:
        aligned, s_down = downsample(aligned, pd)
    else:
        s_down = dict(drawn)
    warnings = []
    if target is not None and len(aligned) < target:
        warnings.append(f"only {len(aligned)} reads available for a target depth of {target}; "
                        "densities from different depths are not comparable")
    s_down["depth_target"] = target
    s_align = {"aligned_to_target": len(aligned)}

    variants, s_var = call_variants(aligned, reference, target_name, p, null_rates)
    if germline:
        kept = [v for v in variants if (v.pos - 1, v.alt) not in germline]
        if len(kept) < len(variants):
            s_var["germline_from_project"] = len(variants) - len(kept)
            s_var["germline_positions_removed"] += len(variants) - len(kept)
            s_var["unique_variants"] = len(kept)
            warnings.append(f"{len(variants) - len(kept)} germline difference(s) removed using "
                            + (germline_source or "the project's other clones; this sample has too "
                               "few clones to tell them apart itself"))
        variants = kept

    # per-read SHM load: how many called-variant sites each read carries, and
    # how many usable bases it contributed (for mutations per bp sequenced;
    # with Sanger, one read = one clone)
    per_read = per_read_table(aligned, variants, p, len(read_fasta(reference)[0][1]))
    per_read_mutations = [r["mutations"] for r in per_read]

    density = metrics.unique_variant_density(
        variants, s_var["analysed_positions"], s_var["mean_depth"], p)
    summary = {
        "shmtrace_version": __version__,
        "params": to_dict(p),
        "reads": {**(read_stats or {}), **s_down, **s_align},
        "calling": s_var,
        "unique_variant_density": round(density, 3),
        "spectrum_unique": metrics.spectrum_matrix(variants, weighted=False),
        "spectrum_weighted": metrics.spectrum_matrix(variants, weighted=True),
        "gc_bias_unique": metrics.gc_at_summary(variants, weighted=False),
        "mutation_load": mutation_load_summary(per_read),
        "warnings": warnings,
    }
    return {"summary": summary, "variants": variants,
            "per_read_mutations": per_read_mutations, "per_read": per_read,
            "_aligned": aligned}


def per_read_table(aligned: list, variants: list, p: Params, reflen: int) -> list[dict]:
    """One row per analysed read (Sanger: per clone), in the order of `aligned`.

    read:      read name (Sanger: the trace file's name)
    strand:    +1 along the reference, -1 reverse-complemented
    mutations: called-variant sites the read carries (base quality >=
               min_base_quality, inside the analysis window)
    bases:     reference positions the read covers inside the analysed window
               (end trims applied) with an A/C/G/T call of quality >=
               min_base_quality: the denominator for mutation frequency
    positions: 1-based positions of the read's mutations
    """
    lo = (p.analysis_window[0] - 1) if p.analysis_window else 0
    hi = (p.analysis_window[1] - 1) if p.analysis_window else reflen - 1
    lo, hi = max(lo, p.end_trim_bp), min(hi, reflen - 1 - p.end_trim_bp)
    minq = p.min_base_quality
    if p.error_model == "sanger" and minq <= 0:
        minq = 20
    variant_sites = {(v.pos - 1, v.alt) for v in variants}
    rows = []
    for r in aligned:
        n_mut, n_bases, where = 0, 0, []
        for rpos, base, q in r.pairs:
            if q < minq or base not in "ACGT" or rpos < lo or rpos > hi:
                continue
            n_bases += 1
            if (rpos, base) in variant_sites:
                n_mut += 1
                where.append(rpos + 1)
        rows.append({"read": r.qname, "strand": getattr(r, "strand", 1),
                     "mutations": n_mut, "bases": n_bases, "positions": where})
    return rows


def mutation_load_summary(per_read: list[dict]) -> dict:
    """Clone/read-level load: mutations per read and per base sequenced."""
    n = len(per_read)
    muts = [r["mutations"] for r in per_read]
    bases = sum(r["bases"] for r in per_read)
    total = sum(muts)
    return {
        "reads": n,
        "mutated_reads": sum(1 for m in muts if m > 0),
        "fraction_mutated": round(sum(1 for m in muts if m > 0) / n, 4) if n else None,
        "mutations_total": total,
        "mean_mutations_per_read": round(total / n, 3) if n else None,
        "median_mutations_per_read": statistics.median(muts) if muts else None,
        "bases_sequenced": bases,
        "mutation_frequency_per_bp": (total / bases) if bases else None,
    }


def run(fastq: str, reference: str, contaminants: list[str], p: Params,
        null_rates: dict | None = None, known_blocks: list | None = None,
        depth: int | None = None) -> dict:
    aligned, s_reads = prepare_reads(fastq, reference, contaminants, p, known_blocks)
    result = analyse(aligned, reference, p, s_reads, null_rates, depth)
    result["_aligned_all"] = aligned
    return result


def rarefaction(aligned: list, reference: str, p: Params, depths: list[int],
                draws: int = 5, null_rates: dict | None = None) -> list[dict]:
    """Unique variants and density at several depths, several random draws each.

    The manuscript's density metric falls as its distinct-call numerator
    plateaus, so a single density is only meaningful at a stated depth. This
    gives the whole curve.
    """
    target_name = read_fasta(reference)[0][0]
    rows = []
    for d in depths:
        if d > len(aligned):
            continue
        counts, dens = [], []
        for seed in range(draws):
            sub = random.Random(seed).sample(aligned, d)
            variants, s_var = call_variants(sub, reference, target_name, p, null_rates)
            counts.append(len(variants))
            dens.append(metrics.unique_variant_density(
                variants, s_var["analysed_positions"], s_var["mean_depth"], p))
        rows.append({"depth": d, "draws": draws,
                     "unique_variants_mean": round(statistics.mean(counts), 2),
                     "unique_variants_sd": round(statistics.pstdev(counts), 2),
                     "density_mean": round(statistics.mean(dens), 3),
                     "density_sd": round(statistics.pstdev(dens), 3)})
    return rows


def _window(p: Params, ref_seq: str) -> tuple[int, int]:
    lo = (p.analysis_window[0] - 1) if p.analysis_window else 0
    hi = (p.analysis_window[1] - 1) if p.analysis_window else len(ref_seq) - 1
    return lo, hi


def write_outputs(result: dict, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        # underscored keys are bulk data for figures, not part of the summary
        summary = dict(result["summary"])
        summary["calling"] = {k: v for k, v in summary["calling"].items()
                              if not k.startswith("_")}
        json.dump(summary, fh, indent=2)
    with open(os.path.join(out_dir, "variants.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "ref", "alt", "supporting_reads", "depth",
                    "frequency", "p_value", "plus_strand_reads", "minus_strand_reads",
                    "strand_bias_p", "homopolymer", "context"])
        for v in result["variants"]:
            w.writerow([v.pos, v.ref, v.alt, v.supporting_reads, v.depth,
                        f"{v.frequency:.4f}", f"{v.p_value:.3e}", v.plus_reads,
                        v.minus_reads, f"{v.strand_bias_p:.3g}", v.homopolymer, v.context])


MODES = {
    # quality-calibrated null from same-run controls, BH; without controls it
    # falls back to the read-quality test. Chosen on the 2026-10-11 benchmark:
    # Geneious-level WT signal with a quarter of its control background
    # (docs/calibration/ERROR_MODELS.md)
    "recommended": {"error_model": "qcal", "multiple_testing": "control"},
    # recommended at a false-discovery rate of 0.20 instead of 0.05: more of the mutations seen in
    # one read, at about four times the background. On held-out controls with planted mutations
    # (docs/calibration/crm/spikein.py, profile recommended_fdr20) one-read recall rose from 30% to
    # 72% (JH4) and 58% to 87% (N114); the false-call share was 0.18 and 0.16
    "sensitive": {"error_model": "qcal", "multiple_testing": "control", "fdr_alpha": 0.20},
    # the v0.1 default: per-position empirical null (or read quality), BH
    "v0.1": {},
    # the manuscript's rule: Geneious approximate p-value, uncorrected, no phasing
    "geneious": {"error_model": "geneious", "multiple_testing": "none", "phase_mode": "off"},
    # quality-calibrated null from controls, BH
    "qcal": {"error_model": "qcal"},
    # strand-aware beta-binomial null from controls, BH
    "strand": {"error_model": "strand"},
    # exact Poisson-binomial on read qualities, BH
    "pbinom": {"error_model": "pbinom"},
    # Sanger clones: every substitution with base Q >= 20 in the trimmed read,
    # no statistical test, germline rule only with >= 4 clones, no phasing
    "sanger": {"error_model": "sanger", "multiple_testing": "none", "phase_mode": "off",
               "min_base_quality": 20, "min_variant_reads": 1},
}

PLATFORMS = ("ont", "pacbio", "illumina", "sanger")
_PLATFORM_ALIASES = {"ont": "ont", "nanopore": "ont", "oxford nanopore": "ont",
                     "pacbio": "pacbio", "hifi": "pacbio", "pacbio hifi": "pacbio",
                     "illumina": "illumina", "sanger": "sanger"}


def normalise_platform(platform: str | None) -> str:
    """Canonical platform name ("ont", "pacbio", "illumina", "sanger")."""
    key = (platform or "ont").strip().lower()
    if key not in _PLATFORM_ALIASES:
        raise ValueError(f"unknown platform {platform!r}; choose from {', '.join(PLATFORMS)}")
    return _PLATFORM_ALIASES[key]


def platform_overrides(platform: str, base: Params | None = None) -> dict:
    """Params overrides that adapt an assay preset to a sequencing platform.

    Long-read platforms keep the preset's amplicon-specific length filters.
    Platforms whose reads are shorter than a long-read amplicon (Illumina,
    Sanger) lower min_align_len, to at most 100 / 200 bp (the preset's own
    value when that is already smaller, if `base` is given).

    ont      map-ont; the preset as written
    pacbio   map-hifi; the preset's EE gate (HiFi reads clear it easily);
             Benjamini-Hochberg over every position x alt (bh_family="all")
    illumina sr; EE gate relaxed to 0.02 (poor read tails are common but the
             merged overlap is high quality); bh_family="all"
    sanger   map-ont (tolerant of a heavily mutated clone, works on single
             500-1000 bp reads); Mott trim at Q20, EE gate 0.05 on the trimmed
             read, no downsampling, no phasing, and the "sanger" calling mode
             (every substitution with base Q >= 20, no test)
    """
    plat = normalise_platform(platform)
    out: dict = {"platform": plat}
    if plat == "ont":
        out.update(minimap2_preset="map-ont")
    elif plat == "pacbio":
        out.update(minimap2_preset="map-hifi", bh_family="all")
    elif plat == "illumina":
        # Illumina qualities top out near Q40 (binned at 37 on NovaSeq), so the
        # single-read cut sits lower than on nanopore
        out.update(minimap2_preset="sr", bh_family="all", single_read_min_quality=30,
                   max_ee_rate=max(0.02, base.max_ee_rate if base else 0.0),
                   min_align_len=min(100, base.min_align_len) if base else 100)
    elif plat == "sanger":
        out.update(minimap2_preset="map-ont", quality_trim="mott", trim_cutoff=0.01,
                   max_ee_rate=0.05,
                   min_align_len=min(200, base.min_align_len) if base else 200,
                   max_total_len=max(3000, base.max_total_len) if base else 3000,
                   downsample_to=None, single_read_min_quality=None, **MODES["sanger"])
    return out


def build_params(preset: str, platform: str | None = None, mode: str | None = "recommended",
                 **overrides) -> Params:
    """Preset, then platform_overrides, then mode_overrides, then `overrides`.

    The one-call way to get Params for a run. Overrides that are None are
    ignored. Under platform "sanger" the "recommended" mode means "sanger".
    """
    base = load_preset(preset)
    plat = normalise_platform(platform or base.platform)
    ov = platform_overrides(plat, base)
    ov.update(mode_overrides(mode, plat))
    ov.update({k: v for k, v in overrides.items() if v is not None})
    if isinstance(ov.get("analysis_window"), list):
        ov["analysis_window"] = tuple(ov["analysis_window"])
    return replace(base, **ov)


def mode_overrides(mode: str | None, platform: str | None = None) -> dict:
    """Named calling profiles (MODES). Append "+filter" for the strand-bias filter.

    With platform="sanger", "recommended" (and None) select the "sanger" mode,
    since a statistical test is meaningless on one read per clone.

    "recommended": read-quality or empirical null, BH correction, phasing
                   that only drops blocks seen in controls (the defaults)
    "sensitive":   recommended at a false-discovery rate of 0.20: more one-read
                   mutations, about one call in five from the background
    "geneious":    the manuscript's rule: Geneious's approximate p-value,
                   uncorrected at max_pvalue, no phasing
    "qcal", "strand", "pbinom": the corresponding error_model with BH
    "sanger":      one read per clone; every substitution with base Q >= 20,
                   no test, no phasing
    """
    sanger = platform is not None and normalise_platform(platform) == "sanger"
    if mode is None:
        return dict(MODES["sanger"]) if sanger else {}
    base, _, extra = mode.partition("+")
    if sanger and base in ("recommended", "sensitive"):
        base = "sanger"
    if base not in MODES or extra not in ("", "filter"):
        raise ValueError(f"unknown mode {mode!r}; choose from {', '.join(MODES)} "
                         "(optionally with +filter)")
    out = dict(MODES[base])
    if extra == "filter":
        out["strand_filter"] = "bias"
    return out


def fit_null(null_fastqs: list[str], target_fastq: str | None, reference: str,
             contaminants: list[str], p: Params, force: bool = False,
             sample_run: str | None = None, null_runs: list[str] | None = None):
    """Fit the empirical null if the controls are usable. Returns (rates, stats, blocks)."""
    from . import null_model
    same_run, why = null_model.check_same_run(null_fastqs, target_fastq,
                                              sample_run=sample_run, control_runs=null_runs)
    if not same_run:
        print(f"[null] WARNING: {why}")
        if not force:
            print("[null] not applied. Re-run with --force-null to override, or use a "
                  "same-run control. Falling back to the read-quality model.")
            return None, {"applied": False, "same_run": False, "why": why}, None
    rates, s_null = null_model.fit(null_fastqs, reference, contaminants, p)
    s_null.update(applied=True, same_run=same_run, why=why,
                  provenance=[os.path.basename(f) for f in null_fastqs])
    print(f"[null] fitted from {s_null['control_files']} control file(s) ({why}), "
          f"{s_null['control_reads_used']} reads, {s_null['sites_modelled']} sites "
          f"(median rate {s_null.get('median_rate')})")
    blocks = s_null.pop("_blocks", None)
    return rates, s_null, blocks


def write_per_read(result: dict, out_dir: str):
    """per_read.csv: one row per analysed read (Sanger: per clone)."""
    with open(os.path.join(out_dir, "per_read.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["read", "strand", "mutations", "bases", "positions"])
        for r in result.get("per_read", []):
            w.writerow([r["read"], r["strand"], r["mutations"], r["bases"],
                        ";".join(str(x) for x in r["positions"])])


def write_rarefaction(rows: list[dict], out_dir: str):
    with open(os.path.join(out_dir, "rarefaction.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["depth"])
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description="Open-source SHM variant pipeline")
    ap.add_argument("--fastq", required=True,
                    help="reads of one sample: FASTQ, FASTA, .ab1/.abi, .seq, BAM (any text "
                         "format may be gzipped), a directory (all files = one sample, e.g. "
                         "Sanger clones), or Illumina pairs as R1.fastq.gz+R2.fastq.gz")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--contaminant", action="append", default=[])
    ap.add_argument("--preset", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--platform", default="ont",
                    help="sequencing platform: ont (default), pacbio, illumina or sanger; "
                         "sets the aligner preset and read filters (see README 'Platforms')")
    ap.add_argument("--mode", default="recommended",
                    help="calling profile: recommended (default; = qcal, or sanger under "
                         "--platform sanger), sensitive (recommended at FDR 0.20), v0.1, geneious (the manuscript's rule), qcal, "
                         "strand, pbinom or sanger; add +filter for "
                         "the strand-bias filter, e.g. qcal+filter")
    ap.add_argument("--min-base-quality", type=int,
                    help="ignore bases below this Phred in the pileup (sanger default 20)")
    # a few common overrides exposed on the CLI (the web UI exposes all of them)
    ap.add_argument("--downsample-to", type=int)
    ap.add_argument("--max-ee-rate", type=float)
    ap.add_argument("--end-trim-bp", type=int)
    ap.add_argument("--plots", action="store_true", help="also write PNG figures")
    ap.add_argument("--null", action="append", default=[], metavar="FASTQ",
                    help="mutation-free control reads (e.g. AID-deficient) used to measure "
                         "the per-position error rate; repeatable. Strongly recommended.")
    ap.add_argument("--phase-mode", choices=["auto", "drop", "flag", "off"],
                    help="minor-haplotype handling (default auto: drop only blocks also "
                         "seen in the --null controls)")
    ap.add_argument("--no-phasing", action="store_true", help="same as --phase-mode off")
    ap.add_argument("--fdr", type=float, metavar="ALPHA",
                    help="false discovery rate for Benjamini-Hochberg (default 0.05)")
    ap.add_argument("--run-id", help="sequencing run of --fastq, if read names do not say")
    ap.add_argument("--null-run-id", help="sequencing run of the --null files")
    ap.add_argument("--force-null", action="store_true",
                    help="use --null control reads even if they come from a different "
                         "sequencing run (not recommended)")
    ap.add_argument("--no-correction", action="store_true",
                    help="disable multiple-testing correction (not recommended)")
    ap.add_argument("--error-model", choices=["auto", "qcal", "strand", "pbinom", "geneious"],
                    help="variant test (default from --mode). 'strand': strand-aware "
                         "beta-binomial null from --null controls or --null-panel")
    ap.add_argument("--strand-filter", choices=["off", "bias", "both"],
                    help="drop calls with lopsided strand support ('bias') or without a "
                         "read on each strand ('both')")
    ap.add_argument("--mask-homopolymers", action="store_true",
                    help="drop calls touching a homopolymer run (flagged either way)")
    ap.add_argument("--null-panel", metavar="JSON",
                    help="saved background panel (from --save-null-panel) used when no "
                         "same-run control is available")
    ap.add_argument("--save-null-panel", metavar="JSON",
                    help="write the null fitted from --null controls to this file")
    ap.add_argument("--rarefaction", metavar="DEPTHS",
                    help="comma-separated depths, e.g. 30,60,125,250,500; writes "
                         "rarefaction.csv (5 random draws per depth)")
    args = ap.parse_args()

    overrides = dict(downsample_to=args.downsample_to, max_ee_rate=args.max_ee_rate,
                     end_trim_bp=args.end_trim_bp, fdr_alpha=args.fdr,
                     phase_mode=args.phase_mode, min_base_quality=args.min_base_quality)
    if args.no_phasing:
        overrides["phase_mode"] = "off"
    if args.no_correction:
        overrides["multiple_testing"] = "none"
    if args.error_model:
        overrides["error_model"] = args.error_model
    if args.strand_filter:
        overrides["strand_filter"] = args.strand_filter
    if args.mask_homopolymers:
        overrides["mask_homopolymers"] = True
    p = build_params(args.preset, platform=args.platform, mode=args.mode, **overrides)

    null_rates = s_null = known = None
    if args.null:
        null_rates, s_null, known = fit_null(
            args.null, args.fastq, args.reference, args.contaminant, p,
            force=args.force_null, sample_run=args.run_id,
            null_runs=[args.null_run_id] if args.null_run_id else None)
        if args.save_null_panel and null_rates is not None:
            from . import null_model
            null_model.save_panel(null_rates, args.save_null_panel, preset=p.name,
                                  controls=[os.path.basename(f) for f in args.null],
                                  run=args.null_run_id or null_model.run_token(args.null[0]),
                                  shmtrace_version=__version__)
            print(f"[null] saved panel to {args.save_null_panel}")
    elif args.null_panel:
        from . import null_model
        null_rates = null_model.load_panel(args.null_panel)
        s_null = {"applied": True, "panel": os.path.basename(args.null_panel),
                  "panel_meta": null_rates.meta, "overdispersion_rho": null_rates.rho}
        print(f"[null] using saved panel {args.null_panel} "
              f"({len(null_rates)} sites, rho {null_rates.rho:.2g})")
    if null_rates is not None and p.error_model == "geneious":
        null_rates = None                # the Geneious rule has no empirical null

    result = run(args.fastq, args.reference, args.contaminant, p, null_rates, known)
    if s_null:
        result["summary"]["null_model"] = s_null
    write_outputs(result, args.out)
    if args.rarefaction:
        depths = [int(x) for x in args.rarefaction.split(",") if x.strip()]
        rows = rarefaction(result["_aligned_all"], args.reference, p, depths,
                           null_rates=null_rates)
        write_rarefaction(rows, args.out)
        print(f"Wrote {args.out}/rarefaction.csv")
    if args.plots:
        from . import viz
        viz.render_all(result, args.out)
        print(f"Wrote {args.out}/spectrum.png, spatial.png, mutation_load.png")
    write_per_read(result, args.out)
    s = result["summary"]
    r = s["reads"]
    for b in r.get("blocks", []):
        print(f"[phasing] {b['action']}: {b['n_sites']}-site block in {b['reads']} reads "
              f"({b['why']})")
    if not r.get("blocks") and p.phase_mode != "off" and p.phase_reads:
        print(f"[phasing] {r.get('reason')}")
    for w in s["warnings"]:
        print(f"[warning] {w}")
    print(f"[{p.name}] reads_in={r['reads_in']} "
          f"passed={r['passed_prefilter']} "
          f"aligned={r['aligned_to_target']} "
          f"tested={s['calling']['candidates_tested']} "
          f"unique_variants={s['calling']['unique_variants']} "
          f"density={s['unique_variant_density']} "
          f"GC={s['gc_bias_unique']['fraction_GC']} "
          f"[{s['calling']['error_model']}, {s['calling']['multiple_testing']}]")
    ml = s["mutation_load"]
    print(f"[{p.platform}] reads={ml['reads']} mutated={ml['mutated_reads']} "
          f"mean_mutations_per_read={ml['mean_mutations_per_read']} "
          f"frequency_per_bp={ml['mutation_frequency_per_bp']}")
    for note in s["calling"].get("notes", []):
        print(f"[note] {note}")
    print(f"Wrote {args.out}/summary.json, variants.csv and per_read.csv")


if __name__ == "__main__":
    main()
