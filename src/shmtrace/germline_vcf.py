"""Germline variants from a VCF (Clair3, DeepVariant, NanoSNP or any caller), removed from the calls.

SHMtrace treats a change carried by at least `germline_freq_threshold` (50%) of a sample's reads
as germline. A heterozygous site in a mouse of mixed background sits near 50% and can fall just
below it, and is then called as a mutation in every read that carries it. A germline caller run
on the same animal (on a tail or non-B-cell sample, or on the reads themselves) knows better:
give its VCF, project-wide or per sample, and every substitution it lists is removed from that
sample's calls.

The VCF must be in the amplicon's coordinates: its CHROM must be the reference record's name,
as it is when the germline caller was run against the same FASTA. Records on another contig are
skipped and counted, and a REF base that disagrees with the reference is reported, since either
means the coordinates do not match. Only substitutions are read (REF and ALT one base each; a
multi-allelic ALT gives one substitution per base). A record whose sample genotype is 0/0 or
missing is not a variant in that animal and is skipped; FILTER other than PASS or "." is skipped.
"""
from __future__ import annotations

import gzip


def read_sites(path: str, contig: str, ref: str) -> tuple[set, dict]:
    """{(0-based position, alt)} of the germline substitutions on `contig`, and what was read."""
    opener = gzip.open if path.endswith(".gz") else open
    sites: set = set()
    st = {"records": 0, "substitutions": 0, "other_contig": 0, "not_snv": 0, "filtered": 0,
          "not_carried": 0, "ref_mismatch": 0, "outside_reference": 0}
    ref_u = ref.upper()
    # aligners and germline callers name a contig by the first word of its FASTA header
    names = {contig, contig.split()[0]} if contig.split() else {contig}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            st["records"] += 1
            chrom, pos, r, alts = f[0], f[1], f[3].upper(), f[4].upper()
            if chrom not in names:
                st["other_contig"] += 1
                continue
            filt = f[6] if len(f) > 6 else "."
            if filt not in ("PASS", "."):
                st["filtered"] += 1
                continue
            if len(f) > 9 and f[8]:
                keys = f[8].split(":")
                vals = f[9].split(":")
                gt = dict(zip(keys, vals, strict=False)).get("GT", "")   # trailing fields may be dropped
                alleles = [a for a in gt.replace("|", "/").split("/") if a not in ("", ".")]
                if not alleles or all(a == "0" for a in alleles):
                    st["not_carried"] += 1
                    continue
            if len(r) != 1:
                st["not_snv"] += 1
                continue
            pos0 = int(pos) - 1
            if not 0 <= pos0 < len(ref_u):
                st["outside_reference"] += 1
                continue
            if ref_u[pos0] != r:
                st["ref_mismatch"] += 1
            for alt in alts.split(","):
                if len(alt) == 1 and alt in "ACGT" and alt != r:
                    sites.add((pos0, alt))
                    st["substitutions"] += 1
                else:
                    st["not_snv"] += 1
    return sites, st


def describe(path: str, st: dict) -> str | None:
    """A warning if the VCF looks like it is not in the amplicon's coordinates."""
    name = path.rsplit("/", 1)[-1]
    if st["records"] and not st["substitutions"] and st["other_contig"]:
        return (f"germline VCF {name}: none of its {st['records']} records is on the amplicon reference; "
                "it must be called against the same FASTA, so CHROM matches the reference name")
    if st["ref_mismatch"]:
        return (f"germline VCF {name}: {st['ref_mismatch']} record(s) have a REF base that differs from the "
                "reference at that position; check that the VCF is in the amplicon's coordinates")
    return None
