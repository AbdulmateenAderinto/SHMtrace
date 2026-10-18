"""Reconstruct the engineered TCR reference and extract reads from the Bassan 2021 deposit.

    python crm/tcr_prepare.py        -> crm/out/tcr/ (reference, provenance, QC) and FASTQ beside the SAMs

The figshare deposit (10.6084/m9.figshare.12803501) holds twelve mapped SAM files from an AID-diversified
MAGE-A1-specific TCR. Their headers name the reference `hT27`, 1,785 bp, and record
`bwa mem -x pacbio ... 4_CCS.fastq`: these are CCS reads aligned with bwa, not subreads, and the deposit is
read here without modification.

The reference sequence itself is not in the deposit. It is reconstructed exactly, not by consensus: every
alignment carries an MD tag, which names the reference base wherever a read disagrees, so CIGAR and MD
together give the reference base under every aligned position. Reads must agree; disagreements are counted
and reported rather than voted on, and a position is only written when every read that covers it agrees.

Read names are `movie/zmw/ccs`, so they identify a ZMW, not a molecular UMI. Nothing here treats them as UMIs
or as independent cells.
"""
import collections
import json
import re
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
DATA = Path.home() / "Nanoshm/public_data/TCR_figshare_12803501"
OUT = CRM / "out/tcr"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MD_RE = re.compile(r"(\d+)|(\^[A-Za-z]+)|([A-Za-z])")


def ref_bases(cigar, md, seq, pos0):
    """{0-based reference position: base} implied by this alignment's CIGAR and MD tag."""
    # the reference positions covered by match-type operations, in order, with their query bases
    aligned, qp, rp = [], 0, pos0
    for n, op in ((int(n), op) for n, op in CIGAR_RE.findall(cigar)):
        if op in "M=X":
            for k in range(n):
                aligned.append((rp + k, seq[qp + k]))
            qp += n
            rp += n
        elif op in "IS":
            qp += n
        elif op in "DN":
            rp += n
        # H and P consume neither
    out, i = {}, 0
    for num, dele, mism in MD_RE.findall(md):
        if num:
            for _ in range(int(num)):
                if i < len(aligned):
                    r, b = aligned[i]
                    out[r] = b                      # a match: the reference base is the read base
                    i += 1
        elif dele:
            pass                                    # deleted reference bases are not under an aligned base
        elif mism:
            if i < len(aligned):
                r, _b = aligned[i]
                out[r] = mism.upper()               # MD names the reference base here
                i += 1
    return out


def scan(path, votes, disagree, limit=4000):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or f[5] == "*":
                continue
            flag = int(f[1])
            if flag & 0x900:                        # secondary or supplementary
                continue
            md = next((x[5:] for x in f[11:] if x.startswith("MD:Z:")), None)
            if md is None:
                continue
            for r, b in ref_bases(f[5], md, f[9].upper(), int(f[3]) - 1).items():
                if b not in "ACGT":
                    continue
                prev = votes.get(r)
                if prev is None:
                    votes[r] = b
                elif prev != b:
                    disagree[r] += 1
            n += 1
            if n >= limit:
                return n
    return n


def main():
    sams = sorted(DATA.glob("*_mapped_new.sam"), key=lambda p: int(p.name.split("_")[0]))
    if not sams:
        sys.exit(f"no SAM files in {DATA}")
    OUT.mkdir(parents=True, exist_ok=True)
    votes, disagree, scanned = {}, collections.Counter(), {}
    for s in sams:
        scanned[s.name] = scan(s, votes, disagree)
    length = max(votes) + 1 if votes else 0
    missing = [i for i in range(length) if i not in votes]
    clean = {r: b for r, b in votes.items() if not disagree[r]}
    seq = "".join(clean.get(i, "N") for i in range(length))
    (OUT / "hT27.fasta").write_text(">hT27 reconstructed from the deposit's CIGAR+MD tags\n"
                                    + "\n".join(seq[i:i + 70] for i in range(0, length, 70)) + "\n")
    # per-sample read counts and the usable fraction
    per_sample = {}
    for s in sams:
        total = mapped = 0
        with open(s) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                f = line.split("\t", 6)
                total += 1
                if len(f) > 5 and f[5] != "*" and not (int(f[1]) & 0x900):
                    mapped += 1
        per_sample[s.name] = {"records": total, "primary_mapped": mapped}
    # reads as FASTQ, taken from the deposit unchanged (SEQ/QUAL as stored, reference-forward)
    fqdir = DATA / "fastq"          # reads stay beside the deposit, outside the repository
    fqdir.mkdir(exist_ok=True)
    written = {}
    for s_ in sams:
        stem = s_.name.replace("_mapped_new.sam", "")
        n = 0
        with open(s_) as fh, open(fqdir / f"{stem}.fastq", "w") as out:
            for line in fh:
                if line.startswith("@"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 11 or f[5] == "*" or (int(f[1]) & 0x900):
                    continue
                if f[9] == "*" or f[10] == "*":
                    continue
                out.write(f"@{f[0]}\n{f[9]}\n+\n{f[10]}\n")
                n += 1
        written[stem] = n
    prov = {"fastq_written": written, "source": "figshare 10.6084/m9.figshare.12803501 (Bassan et al. 2021, Eur J Immunol)",
            "reference_name": "hT27", "header_length": 1785, "reconstructed_length": length,
            "reconstructed_from": "CIGAR + MD tags of the deposited alignments; no consensus vote",
            "aligner_in_header": "bwa mem -x pacbio (reads are CCS, not subreads)",
            "positions_with_disagreement": len(disagree), "positions_not_covered": len(missing),
            "ambiguous_written_as_N": seq.count("N"),
            "reads_scanned_per_file": scanned, "per_sample_reads": per_sample}
    json.dump(prov, open(OUT / "reference_provenance.json", "w"), indent=1)
    print(f"reference hT27: header says {prov['header_length']} bp, reconstructed {length} bp, "
          f"{seq.count('N')} ambiguous")
    print(f"positions where reads disagreed on the reference base: {len(disagree)}")
    for name, v in per_sample.items():
        print(f"   {name:34s} records {v['records']:>6}  primary mapped {v['primary_mapped']:>6}")


if __name__ == "__main__":
    main()
