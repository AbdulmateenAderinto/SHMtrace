"""Parameter model and presets for the SHM analysis engine.

Every knob the user can set in the web UI (sliders / text boxes) maps to one
field here. Presets are just serialised Params with sensible defaults for a
given assay (JH4 intron, N114 V region, ...). Loading a preset then overlaying
user overrides is the whole configuration story.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib import resources
from pathlib import Path

import yaml


@dataclass
class Params:
    # ---- assay identity -----------------------------------------------------
    name: str = "custom"
    description: str = ""
    # Sequencing platform the reads come from: "ont" (Oxford Nanopore),
    # "pacbio" (HiFi/CCS), "illumina" (merged or single reads) or "sanger"
    # (one read per clone). The field records the platform; the matching
    # settings come from pipeline.platform_overrides() / pipeline.build_params().
    platform: str = "ont"

    # ---- read pre-processing (Step 1-2 of the current pipeline) -------------
    minimap2_preset: str = "map-ont"      # long noisy reads (Nanopore)
    min_align_len: int = 800              # read must align at least this long to target
    max_total_len: int = 2000             # discard over-length concatemers/chimeras
    max_ee_rate: float = 0.01             # expected errors / base  (Q20-equivalent)
    drop_supplementary: bool = True       # remove structural chimeras
    downsample_to: int | None = 500    # equalise depth across samples; None = keep all
    downsample_seed: int = 0
    # Per-read quality trimming before the length/EE gate:
    #   "off"  none (long-read and Illumina default)
    #   "mott" modified-Mott trim to the best segment above trim_cutoff
    #          (Sanger default; 0.01 = Q20, 0.05 = Q13)
    quality_trim: str = "off"
    trim_cutoff: float = 0.01
    # Bases below this Phred are left out of the pileup (depth, calls and
    # per-read counts). 0 = use every base; the Sanger profile uses 20.
    min_base_quality: int = 0

    # ---- competitive contaminant handling ----------------------------------
    # A read is assigned to whichever indexed sequence it aligns best to. Reads
    # whose best hit is a contaminant are dropped. Then an allele-level check:
    # a read carrying contaminant_min_alleles or more of a contaminant's bases at
    # sites where it differs from the reference is dropped too (a paralogue read
    # longer than the contaminant record, or a PCR chimera). Sites where at least
    # contaminant_germline_frac of the sample's reads carry the contaminant base
    # are shared germline, not evidence. 0 turns the check off. On the J꜀4
    # deposit, real reads carry 0-7 such bases and paralogue reads and chimeras
    # 11-62 (docs/calibration/QUANT_ACCURACY.md).
    contaminant_allele_frac: float = 0.5
    contaminant_min_alleles: int = 10
    contaminant_germline_frac: float = 0.2

    # ---- UMI consensus (umi.py) --------------------------------------------
    # "off", "flank" (UMI next to an anchor at a read end) or "header" (in the
    # read name). With UMIs on, each molecule becomes one consensus read and
    # every count downstream is of molecules.
    umi: str = "off"
    umi_ends: str = "both"           # flank: "left", "right" or "both" (dual UMIs)
    umi_length: int = 18
    umi_offset: int = 0              # bases between the UMI and its anchor
    umi_anchor_length: int = 20      # default anchors: this many reference bases at each end
    umi_left_anchor: str = ""        # or give the anchors outright (reference-forward)
    umi_right_anchor: str = ""
    umi_anchor_edits: int = 6
    umi_header_regex: str = r"(?:_|UMI[=:])([ACGTN]{6,})"
    umi_max_edits: int = 3           # per UMI end, when grouping
    umi_min_reads: int = 3           # smaller groups are dropped

    # ---- variant calling ----------------------------------------------------
    analysis_window: tuple[int, int] | None = None  # 1-based inclusive; None = whole ref
    germline_freq_threshold: float = 0.50  # positions with an alt >= this are germline -> removed
    # Opt-in, for experiments whose reference is the construct that went into the cells (an engineered
    # receptor, a base-editor substrate). There, a change acquired during diversification can be driven
    # high or to fixation by selection, and inferring germline from frequency would delete exactly the
    # variant the experiment is about. With this on, the frequency rule is off for calling and phasing.
    # It does not change the error model, which still keeps fixed differences out of the background, and
    # it does not change what a baseline difference shared with the controls does: the control-calibrated
    # test already declines to call that. Leave it off for physiological SHM.
    ancestral_reference: bool = False
    max_pvalue: float = 0.05               # keep variants with p <= this
    per_base_error: float = 0.001          # null error rate for the binomial variant test
    min_variant_reads: int = 1             # minimum supporting reads
    exclude_indels: bool = True            # SNVs only
    end_trim_bp: int = 0                   # drop calls within N bp of the reference ends (N114 = 20)

    # How each candidate is tested against sequencing error:
    #   "auto"     per-column read-quality Poisson, or the empirical null when
    #              control reads are supplied (recommended)
    #   "geneious" Geneious approximate p-value (binomial at the mean quality of
    #              the variant bases); use with
    #              multiple_testing="none" to reproduce the manuscript's rule
    #   "strand"   strand-aware beta-binomial null fitted from control reads:
    #              each strand tested against its own error rate, with
    #              overdispersion between libraries; falls back to "auto"
    #              without controls
    #   "qcal"     quality-calibrated: like Geneious, a high-quality variant base
    #              counts for more, but how often bases of that quality are
    #              miscalls is measured in the control reads (per site, base
    #              change and strand) instead of taken from the quality score
    #   "pbinom"   exact Poisson-binomial over every read's quality (LoFreq's
    #              test) instead of the Poisson approximation
    #   "sanger"   one read per clone: every substitution with base quality
    #              >= min_base_quality is called, no test (see below)
    error_model: str = "auto"
    # error_model "qcal": read each base's quality through a table fitted on the controls, keyed
    # by the reported quality and the lowest quality among its two neighbours on each side
    # (recal.py). Errors cluster, so a base next to a low-quality one is more often wrong than
    # its own quality says, and one in a clean stretch less often.
    recalibrate: bool = False
    # Give every base the same quality (FASTA_DEFAULT_Q) after the read filters. Set by
    # run_project when a control has no base qualities (FASTA, .seq, GenBank) while samples do:
    # judging quality-bearing samples against a control whose qualities are invented calls
    # sample errors at high quality as mutations, so both are judged on the same footing.
    ignore_qualities: bool = False
    # How the two strands' p-values are combined under error_model="strand":
    #   "fisher"  Fisher's method (evidence adds up across strands)
    #   "both"    the larger of the two (a call needs evidence on each strand)
    strand_combine: str = "fisher"
    # pseudo-reads pulling each site's strand rate toward the pooled rate for
    # its base change and strand (light, so a site the controls show erring
    # often keeps its own high rate)
    strand_prior_trials: float = 50.0
    # Post-test strand filter, any error model:
    #   "off"   none
    #   "bias"  drop a call whose supporting reads are split between strands
    #           differently from the column (binomial p < strand_bias_alpha);
    #           needs at least strand_bias_min_reads supporting reads to judge
    #   "both"  drop a call unless at least one supporting read is on each strand
    strand_filter: str = "off"
    strand_bias_alpha: float = 0.01
    strand_bias_min_reads: int = 3
    # Homopolymers are where nanopore basecalling errs most. Calls whose
    # position touches a run of >= homopolymer_flag identical bases are flagged
    # in variants.csv; mask_homopolymers drops them.
    homopolymer_flag: int = 4
    mask_homopolymers: bool = False
    # error_model="sanger" (one read = one clone): every substitution in the
    # trimmed read with base quality >= min_base_quality (20 if unset) is a
    # call, no statistical test. The germline rule only applies at positions
    # covered by at least this many clones; with fewer clones in the whole
    # sample it is switched off and a note says so.
    sanger_germline_min_clones: int = 4
    # The first and last stretch of a Sanger read is messy (ambiguity codes,
    # spurious indels and mismatches), and a GenBank or FASTA export carries no
    # qualities to trim it by. After alignment each clone is cut back, from
    # each end, to the first run of this many reference bases with no ambiguity
    # code, no indel and at most sanger_clean_max_mismatch mismatches; what lies
    # outside is dropped. 0 turns it off.
    sanger_clean_window: int = 20
    sanger_clean_max_mismatch: int = 1

    # ---- multiple testing ---------------------------------------------------
    # The caller tests (analysed positions x 3 alt bases) hypotheses per sample,
    # ~1,587 for the 529 bp JH4 window. Uncorrected p <= 0.05 therefore expects
    # ~79 false calls by chance alone, which dominates the apparent background.
    # See docs/calibration/FINDINGS.md.
    # "bh" (Benjamini-Hochberg), "none", or "control": the threshold is set where
    # the false calls predicted from the held-out controls (each control called
    # against the others at the samples' depth) are at most fdr_alpha of the
    # calls. BH over candidates assumes their p-values are uniform under the
    # null, which a test conditioned on seeing the allele does not give; the
    # controls measure the actual null distribution instead. Falls back to BH
    # with fewer than two controls. The recommended mode uses it.
    multiple_testing: str = "bh"
    # how the false calls expected at a threshold are taken from the decoy controls:
    # "mean" over controls, or "max" (the noisiest control; guards against a sample
    # that errs more than the average control)
    control_decoy_estimate: str = "mean"
    # a control far noisier than the others (project.noisy_decoys) is left out of the decoys,
    # so one contaminated or poor library cannot raise every sample's threshold; it still
    # contributes to the error model, and the project notes name it
    exclude_noisy_decoys: bool = True
    # processes for per-library work (parallel.py): 1 = this process only (default);
    # 0 = one fewer than the cores, up to 8. Results do not depend on it. Worker
    # processes are started by "spawn", so a script that sets jobs != 1 must run its
    # analysis under `if __name__ == "__main__":` (the shmtrace commands do).
    jobs: int = 1
    fdr_alpha: float = 0.05                # target false discovery rate
    # Which hypotheses Benjamini-Hochberg counts:
    #   "candidates" only (position, alt) pairs seen in at least one read (the
    #                calibrated nanopore behaviour; on noisy reads nearly every
    #                pair is seen, so this is close to "all")
    #   "all"        every analysed position x 3 alt bases. Needed for accurate
    #                reads (Illumina, PacBio HiFi): there only a few error
    #                singletons are seen, and counting just those lets them all
    #                through. The illumina and pacbio platform profiles set it.
    bh_family: str = "candidates"

    # ---- empirical error model (optional) -----------------------------------
    # Nanopore error is context-driven and correlated across reads, so the
    # quality-derived null underestimates it by ~3x. Given AID-deficient (or
    # otherwise mutation-free) control reads, the per-position, per-alt error
    # rate can be measured directly and used as the test null instead.
    null_prior_alpha: float = 0.5          # beta prior, pseudo-successes
    null_prior_beta: float = 500.0         # beta prior, pseudo-trials (shrinkage strength)
    # An unusually clean control can drive the fitted rate below what the
    # sequencer's own quality scores justify, which makes the test *more*
    # permissive and adds calls. Never trust the empirical rate below the
    # read-quality expectation for that column. See docs/calibration/FINDINGS.md.
    null_rate_floor: bool = True

    # ---- read phasing -------------------------------------------------------
    # A second haplotype, paralogue, or surviving contaminant appears as a block
    # of positions that co-occur in the same subset of reads. Averaging those
    # reads into one pileup inflates the background several-fold. Real SHM does
    # not co-occur this way: each read carries its own mutations.
    phase_reads: bool = True               # False is the same as phase_mode="off"
    # "auto": drop a block only if it matches one seen in control reads or in
    #         another sample; flag the rest (an expanded clone looks the same)
    # "drop": drop every block (pre-October-2026 behaviour)
    # "flag": report blocks, never drop
    # "off":  skip block detection
    phase_mode: str = "auto"
    phase_min_freq: float = 0.01           # a site must vary this much to be tested for co-occurrence
    phase_min_sites: int = 4               # need at least this many candidate sites to phase
    phase_min_cover: int = 3               # a read must cover this many sites to be classified
    phase_score_threshold: float = 0.6     # read carrying >= this fraction of alts = minor haplotype
    phase_max_drop: float = 0.45           # refuse to drop more than this share of reads

    # ---- single-read mutations ----------------------------------------------
    # Mutations seen in one read are estimated in aggregate, as the excess of
    # single-read events over the controls' rate, rather than called one by one
    # (singletons.py, docs/calibration/SINGLE_READ.md). Only events whose base
    # has at least this quality are counted; the estimate is corrected for the
    # share of real bases the cut removes. 35 suits current nanopore (Dorado
    # sup, Q capped at 50); 0 or None turns the estimate off.
    single_read_min_quality: int | None = 35

    # ---- indels ---------------------------------------------------------------
    # exclude_indels=False calls insertions and deletions (indels.py). Only
    # indels whose quality (weakest inserted base, or weaker base flanking a
    # deletion) reaches this are counted, in controls and samples alike:
    # nanopore indel errors mostly come with low-quality flanks.
    indel_min_quality: int = 30
    # A knockout library holds a few real indel molecules (PCR-made deletions,
    # for instance) in two or three reads; no error model removes those, and a
    # floor on support does. See docs/calibration/INDELS.md.
    indel_min_reads: int = 3

    # ---- clones from the V(D)J junction --------------------------------------
    # When an amplicon starts in the V region, the bases just 5' of the J
    # segment are the V(D)J junction, unique to each B-cell clone. Reads are
    # grouped by it (clones.py) so a clone counts once however many reads it
    # contributes. junction_ref_start is the 1-based reference position where
    # the J coding segment begins (J꜀4: 53); None turns this off.
    junction_ref_start: int | None = None
    junction_length: int = 20
    junction_max_edits: int = 3

    # ---- reporting ----------------------------------------------------------
    density_scale: float = 1e5             # unique-variant density is reported x this


def preset_path(name_or_path: str) -> Path:
    """Resolve a preset by bare name (packaged) or by filesystem path.

    >>> preset_path("jh4_intron").name
    'jh4_intron.yaml'
    """
    p = Path(name_or_path)
    if p.suffix in (".yaml", ".yml") and p.exists():
        return p
    stem = p.stem or str(name_or_path)
    packaged = resources.files("shmtrace.presets") / f"{stem}.yaml"
    if packaged.is_file():
        return Path(str(packaged))
    raise FileNotFoundError(
        f"No preset {name_or_path!r}. Bundled presets: {', '.join(available_presets())}"
    )


def available_presets() -> list[str]:
    return sorted(
        f.name[:-5] for f in resources.files("shmtrace.presets").iterdir()
        if f.name.endswith(".yaml")
    )


def load_preset(path: str, **overrides) -> Params:
    with open(preset_path(path)) as fh:
        data = yaml.safe_load(fh) or {}
    if data.get("analysis_window"):
        data["analysis_window"] = tuple(data["analysis_window"])
    data.update({k: v for k, v in overrides.items() if v is not None})
    return Params(**data)


def to_dict(p: Params) -> dict:
    return asdict(p)
