"""What SHMtrace did with every delivered read of a library, step by step.

    from read_accounting import account
    account(fastq, reference, contaminants, region, depth=500)

Runs SHMtrace's own read preparation (pipeline.prepare_reads) with the project's parameters and
tallies the fate it records for each read:

  delivered            reads in the FASTQ
  too_long             longer than the preset's max_total_len (no minimum length is applied)
  low_quality          mean expected error rate above max_ee_rate
  not_target           unmapped, best match to a paralogue/contaminant reference, or off target
  chimeric             supplementary alignments (two parts of the read map apart)
  short_alignment      aligned length below min_align_len
  paralogue_alleles    carries >= contaminant_min_alleles paralogue-specific bases
  minor_haplotype      a read-phased minor haplotype block (phasing)
  usable               reads left for analysis
  analysed             reads analysed (a random draw of `depth` when more are usable)
"""
from shmtrace import pipeline, project

GROUPS = {"too long": "too_long", "low quality": "low_quality", "unmapped": "not_target",
          "contaminant": "not_target", "off target": "not_target", "chimeric": "chimeric",
          "alignment too short": "short_alignment", "contaminant alleles": "paralogue_alleles",
          "minor haplotype": "minor_haplotype"}
ORDER = ["delivered", "too_long", "low_quality", "not_target", "chimeric", "short_alignment",
         "paralogue_alleles", "minor_haplotype", "usable", "analysed"]


def account(fastq, reference, contaminants, region, depth=500, mode="recommended", overrides=None):
    p = project.build_params({"region": region, "mode": mode, "overrides": overrides or {}})
    fates = {}
    aligned, stats = pipeline.prepare_reads(str(fastq), str(reference), [str(c) for c in contaminants], p,
                                            fates=fates)
    out = dict.fromkeys(ORDER, 0)
    out["delivered"] = stats.get("reads_in", 0)
    for f in fates.values():
        key = GROUPS.get(f.get("fate"))
        if key:
            out[key] += 1
    out["usable"] = len(aligned)
    out["analysed"] = min(len(aligned), depth) if depth else len(aligned)
    return out
