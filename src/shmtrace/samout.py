"""The analysed reads of a sample as a SAM alignment (and BAM, where pysam is installed).

So the same reads can be opened in IGV, Geneious or samtools. Each read the analysis kept
is taken from the sample's input file and aligned again to the run's reference with the
aligner and preset the run used, so its CIGAR, bases and qualities are those of the read
as delivered; trimming the analysis applied before aligning (Sanger's quality trim) is not
re-applied, and the ends it would have trimmed appear soft-clipped instead. Reads the
input no longer has under their analysed names (UMI consensus molecules) are left out.
"""
from __future__ import annotations

import io
import os
import tempfile
import zipfile

from .align import _revcomp, build_index
from .io_utils import read_fasta, read_fastq

PRESET = {"ont": "map-ont", "pacbio": "map-hifi", "illumina": "sr", "sanger": "map-ont"}


def _cigar(hit, length: int) -> str:
    """minimap2's CIGAR with the unaligned read ends as soft clips, reference-forward."""
    lead, trail = hit.q_st, length - hit.q_en
    if hit.strand == -1:
        lead, trail = trail, lead
    ops = "MIDNSHP=X"
    body = "".join(f"{n}{ops[op]}" for n, op in hit.cigar)
    return f"{lead}S" * bool(lead) + body + f"{trail}S" * bool(trail)


def sam(reference: str, reads_input, names: set[str], platform: str, sample: str) -> tuple[str, int]:
    """(SAM text sorted by position, reads written)."""
    ref_name, ref_seq = read_fasta(reference)[0]
    aligner, _t, _c = build_index([reference], [], PRESET.get(platform, "map-ont"))
    rows = []
    for qname, seq, qual in read_fastq(reads_input):
        if qname not in names:
            continue
        hit = max((h for h in aligner.map(seq) if h.is_primary), key=lambda h: h.mlen, default=None)
        if hit is None:
            continue
        rev = hit.strand == -1
        out_seq, out_q = (_revcomp(seq), qual[::-1]) if rev else (seq, qual)
        rows.append((hit.r_st, "\t".join([
            qname, "16" if rev else "0", ref_name, str(hit.r_st + 1), str(min(60, hit.mapq)),
            _cigar(hit, len(seq)), "*", "0", "0", out_seq, out_q or "*",
            f"NM:i:{hit.NM}", f"RG:Z:{sample}"])))
    rows.sort(key=lambda r: r[0])
    head = [f"@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:{ref_name}\tLN:{len(ref_seq)}",
            f"@RG\tID:{sample}\tSM:{sample}", "@PG\tID:shmtrace\tPN:SHMtrace"]
    return "\n".join(head + [r[1] for r in rows]) + "\n", len(rows)


def bam_bytes(sam_text: str) -> tuple[bytes, bytes] | None:
    """(BAM, its .bai index) from SAM text, or None without pysam."""
    try:
        import pysam
    except ImportError:
        return None
    with tempfile.TemporaryDirectory() as d:
        s, b = os.path.join(d, "x.sam"), os.path.join(d, "x.bam")
        with open(s, "w") as fh:
            fh.write(sam_text)
        with pysam.AlignmentFile(s, "r") as src, pysam.AlignmentFile(b, "wb", template=src) as dst:
            for rec in src:
                dst.write(rec)
        pysam.index(b)
        with open(b, "rb") as fh, open(b + ".bai", "rb") as fi:
            return fh.read(), fi.read()


def zip_bytes(per_sample: dict[str, str]) -> bytes:
    """Every sample's SAM, and BAM with its index where pysam is available."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for sample, text in per_sample.items():
            safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in sample)
            z.writestr(f"{safe}.sam", text)
            bb = bam_bytes(text)
            if bb:
                z.writestr(f"{safe}.bam", bb[0])
                z.writestr(f"{safe}.bam.bai", bb[1])
    return buf.getvalue()
