"""Analysed reads as an AIRR Community Rearrangement table, for Immcantation and other AIRR tools.

One row per analysed read (per clone for Sanger). The AIRR-required fields are all present;
those SHMtrace does not determine (V/D/J calls, CIGARs, the CDR3 junction) are left empty, which
the standard allows.

`sequence_alignment` is the read with sequencing error removed: the germline base at every
aligned reference position, except where the read carries a mutation SHMtrace called in that
sample. Tools that count mutations from `sequence_alignment` against `germline_alignment`
(SHazaM's observedMutations, BASELINe, dowser) then count SHMtrace's calls, not every
mismatch of a noisy read. The raw read, aligned the same way, is in `sequence_alignment_raw`.
Both are in reference coordinates: a deleted reference base is "-", inserted bases are left
out, so the two alignments and the germline have the same length.

Fields beyond the standard carry the SHMtrace prefix where AIRR has no equivalent:
`repertoire_id` is the sample, `clone_id` the junction clone (clones.py) when the region has
one, and `shmtrace_group`, `shmtrace_mutations` and `shmtrace_positions` the sample's group and
the read's called mutations (1-based position, reference base, read base).
"""
from __future__ import annotations

import csv
import io

REQUIRED = ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call",
            "sequence_alignment", "germline_alignment", "junction", "junction_aa",
            "v_cigar", "d_cigar", "j_cigar"]
EXTRA = ["sequence_alignment_raw", "germline_start", "germline_end", "repertoire_id", "clone_id",
         "shmtrace_group", "shmtrace_mutations", "shmtrace_positions"]
FIELDS = REQUIRED + EXTRA


def read_rows(sample: str, group: str, aligned, variants, ref: str,
              clone_of: dict | None = None) -> list[dict]:
    """AIRR rows for one sample's analysed reads."""
    ref_u = ref.upper()
    called = {(v.pos - 1, v.alt) for v in variants}
    rows = []
    for r in aligned:
        cols = {rpos: base for rpos, base, _q in r.pairs}
        if not cols:
            continue
        lo, hi = min(cols), max(cols)
        raw, clean, muts = [], [], []
        for rpos in range(lo, hi + 1):
            base = cols.get(rpos)
            if base is None:                       # deleted in the read
                raw.append("-")
                clean.append("-")
                continue
            raw.append(base)
            if base != ref_u[rpos] and (rpos, base) in called:
                clean.append(base)
                muts.append(f"{rpos + 1}{ref_u[rpos]}>{base}")
            else:
                clean.append(ref_u[rpos])
        c = clone_of.get(r.qname) if clone_of else None
        rows.append({
            "sequence_id": r.qname, "sequence": r.seq, "rev_comp": "T" if getattr(r, "strand", 1) == -1 else "F",
            "productive": "", "v_call": "", "d_call": "", "j_call": "",
            "sequence_alignment": "".join(clean), "germline_alignment": ref_u[lo:hi + 1],
            "junction": "", "junction_aa": "", "v_cigar": "", "d_cigar": "", "j_cigar": "",
            "sequence_alignment_raw": "".join(raw), "germline_start": lo + 1, "germline_end": hi + 1,
            "repertoire_id": sample, "clone_id": f"{sample}_{c}" if c is not None else "",
            "shmtrace_group": group, "shmtrace_mutations": len(muts), "shmtrace_positions": ";".join(muts)})
    return rows


def write_tsv(rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def project_rows(res: dict) -> list[dict]:
    """Rows for every sample of a project result (run_project keeps the reads in `_details`)."""
    details = res.get("_details") or {}
    groups = {s["sample"]: s.get("group") or s["sample"] for s in res.get("samples", [])}
    ref = res.get("reference")
    if isinstance(ref, dict):
        ref = ref.get("sequence")
    rows = []
    for name, d in details.items():
        aligned = d.get("_aligned")
        if not aligned or not ref:
            continue
        rows += read_rows(name, groups.get(name, name), aligned, d.get("variants", []), ref, d.get("_clone_of"))
    return rows
