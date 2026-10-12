"""Class-switch junctions from long reads.

Class-switch recombination joins a break in the donor switch region (Sμ) to a break
in an acceptor (Sγ1, Sγ3, Sε, Sα, ...), deleting what lies between. A long read
across the join aligns in two parts: its 5' part to the donor, its 3' part to the
acceptor. From those two alignments, for each read:

  breakpoints       the last donor base and the first acceptor base, in each
                    reference's coordinates
  microhomology     bases at the join that match both donor and acceptor (the read
                    cannot tell which side they came from). 0 = blunt; 1-3 short;
                    4 or more is typical of alternative end-joining (e.g. when
                    classical NHEJ factors are missing)
  insertion         bases at the join that match neither side
  mutations         base differences (Q >= 20) within `flank` bp of the join on
                    each side, per base read: AID-induced mutations cluster there
  ambiguous         the acceptor part also aligns as well somewhere else in the
                    acceptor (switch regions are repetitive), so its breakpoint is
                    not unique
  internal deletion both parts on the donor, more than `min_deletion` bp apart: an
                    intra-switch deletion rather than a switch

Reads are first oriented to the donor. Junctions whose breakpoints agree within 3 bp
on both sides are counted once, as unique junctions: reads of one junction are PCR
copies of one recombined cell (sequencing errors move a breakpoint a base or two).

Accuracy on simulated junction reads (docs/calibration/SWITCH.md): every switched
read is found; breakpoints fall within the microhomology + 3 bp for 91-95%;
microhomology of 4 nt or more is classed correctly 97% (2% read error) and 82% (6%);
blunt and short joins are sometimes read as 1-2 nt insertions, because a sequencing
indel at the join looks like one.

The references: one FASTA record per switch region, named as given to `donor`
(default: the first record) and the acceptors (all others).
"""
from __future__ import annotations

import csv
from collections import Counter
from dataclasses import asdict, dataclass

import mappy as mp

from .io_utils import phred, read_fasta, read_fastq


@dataclass
class Junction:
    read: str
    donor: str
    donor_break: int            # 1-based, last donor base
    acceptor: str
    acceptor_break: int         # 1-based, first acceptor base
    kind: str                   # "switch" or "internal deletion"
    microhomology: int
    insertion: str
    ambiguous: bool
    mutations_near: int
    bases_near: int
    deletion_bp: int | None = None

    @property
    def key(self):
        return (self.donor, self.donor_break, self.acceptor, self.acceptor_break, self.insertion)


def _mismatches(seq, quals, h, ref, q_lo, q_hi, min_q=20):
    """(mutations, bases) of hit h within read positions [q_lo, q_hi)."""
    mut = bases = 0
    rp, qp = h.r_st, h.q_st
    for length, op in h.cigar:
        if op in (0, 7, 8):
            for k in range(length):
                q = qp + k
                if q_lo <= q < q_hi and quals[q] >= min_q:
                    bases += 1
                    mut += seq[q] != ref[rp + k]
            rp += length
            qp += length
        elif op == 1:
            qp += length
        elif op in (2, 3):
            rp += length
    return mut, bases


def call_read(name, seq, qual, aligner, refs, donor, flank=50, min_deletion=200,
              min_part=60) -> Junction | None:
    hs = list(aligner.map(seq))
    don = [h for h in hs if h.ctg == donor and h.q_en - h.q_st >= min_part]
    if not don:
        return None
    d0 = max(don, key=lambda h: h.mlen)
    if d0.strand == -1:                          # orient the read to the donor
        seq = mp.revcomp(seq)
        qual = qual[::-1]
        hs = list(aligner.map(seq))
        don = [h for h in hs if h.ctg == donor and h.q_en - h.q_st >= min_part and h.strand == 1]
        if not don:
            return None
        d0 = max(don, key=lambda h: h.mlen)
    quals = phred(qual)
    # the acceptor part: a forward hit starting downstream of the donor part's start,
    # on an acceptor, or on the donor far downstream (an internal deletion)
    cands = [h for h in hs if h.strand == 1 and h is not d0 and h.q_st > d0.q_st + min_part // 2
             and h.q_en - h.q_st >= min_part
             and (h.ctg != donor or h.r_st - d0.r_en >= min_deletion)]
    if not cands:
        return None
    a0 = max(cands, key=lambda h: h.mlen)
    dref, aref = refs[donor], refs[a0.ctg]
    # the join: the donor part ends at read d0.q_en (ref d0.r_en); the acceptor part
    # starts at read a0.q_st (ref a0.r_st). Overlap = bases both alignments claim.
    d_end_q, d_end_r = d0.q_en, d0.r_en
    a_st_q, a_st_r = a0.q_st, a0.r_st
    if a_st_q < d_end_q:                      # overlapping claims: move the split to the middle
        ov = d_end_q - a_st_q
        d_end_q -= ov - ov // 2
        d_end_r -= ov - ov // 2
        a_st_q = d_end_q
        a_st_r += ov // 2
    insertion = seq[d_end_q:a_st_q] if a_st_q > d_end_q else ""
    # microhomology: extend the donor match forward and the acceptor match backward
    # across the join; bases explained by both sides are shared
    mh = 0
    if not insertion:
        k = 0
        while (d_end_q - 1 - k >= 0 and a_st_r - 1 - k >= 0
               and seq[d_end_q - 1 - k] == aref[a_st_r - 1 - k] and k < 30):
            k += 1                               # read's last donor bases also match acceptor
        j = 0
        while (a_st_q + j < len(seq) and d_end_r + j < len(dref)
               and seq[a_st_q + j] == dref[d_end_r + j] and j < 30):
            j += 1                               # read's first acceptor bases also match donor
        mh = k + j
    ambiguous = sum(1 for h in hs if h.ctg == a0.ctg and h is not a0 and h.strand == 1
                    and abs(h.q_st - a0.q_st) < 20 and h.mlen >= 0.95 * a0.mlen) > 0
    m1, b1 = _mismatches(seq, quals, d0, dref, d_end_q - flank, d_end_q)
    m2, b2 = _mismatches(seq, quals, a0, aref, a_st_q, a_st_q + flank)
    internal = a0.ctg == donor
    return Junction(name, donor, d_end_r, a0.ctg, a_st_r + 1,
                    "internal deletion" if internal else "switch", mh, insertion, ambiguous,
                    m1 + m2, b1 + b2, (a_st_r - d_end_r) if internal else None)


def analyse(path, reference: str, donor: str | None = None, preset: str = "map-ont",
            flank: int = 50) -> dict:
    recs = read_fasta(reference)
    refs = {n: s.upper() for n, s in recs}
    donor = donor or recs[0][0]
    aligner = mp.Aligner(reference, preset=preset, best_n=10)
    if not aligner:
        raise RuntimeError("could not index the switch references")
    out, n = [], 0
    for name, seq, qual in read_fastq(path):
        n += 1
        j = call_read(name, seq, qual, aligner, refs, donor, flank=flank)
        if j:
            out.append(j)
    return {"reads": n, "junctions": out, "summary": summarise(out)}


def mh_class(m: int, ins: str) -> str:
    if ins:
        return "insertion"
    return "blunt" if m == 0 else ("1-3 nt" if m <= 3 else ">=4 nt")


def unique_junctions(js: list[Junction], slack: int = 3) -> list[Junction]:
    """One junction per recombination event: reads whose breakpoints agree within
    `slack` bp on both sides (sequencing errors move a breakpoint a base or two) are
    copies of one join. The commonest member represents the event."""
    groups: list[list[Junction]] = []
    for j in sorted(js, key=lambda x: (x.acceptor, x.donor_break, x.acceptor_break)):
        for g in groups:
            h = g[0]
            if (h.acceptor == j.acceptor and abs(h.donor_break - j.donor_break) <= slack
                    and abs(h.acceptor_break - j.acceptor_break) <= slack):
                g.append(j)
                break
        else:
            groups.append([j])
    out = []
    for g in groups:
        common = Counter((j.donor_break, j.acceptor_break, j.microhomology, j.insertion)
                         for j in g).most_common(1)[0][0]
        out.append(next(j for j in g if (j.donor_break, j.acceptor_break, j.microhomology,
                                          j.insertion) == common))
    return out


def summarise(js: list[Junction]) -> dict:
    sw = [j for j in js if j.kind == "switch"]
    uniq = {j.key: j for j in unique_junctions(sw)}
    cls = Counter(mh_class(j.microhomology, j.insertion) for j in uniq.values())
    nu = len(uniq) or 1
    mut = sum(j.mutations_near for j in sw)
    bp = sum(j.bases_near for j in sw)
    return {"junction_reads": len(sw), "unique_junctions": len(uniq),
            "internal_deletions": sum(1 for j in js if j.kind == "internal deletion"),
            "acceptors": dict(Counter(j.acceptor for j in uniq.values())),
            "microhomology_classes_pct": {k: round(100 * v / nu, 1) for k, v in cls.items()},
            "mean_microhomology": round(sum(j.microhomology for j in uniq.values()) / nu, 2),
            "ambiguous_pct": round(100 * sum(j.ambiguous for j in uniq.values()) / nu, 1),
            "mutations_near_join_per_1e3": round(1e3 * mut / bp, 3) if bp else None}


def write_csv(js: list[Junction], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(js[0])) if js else ["read"])
        w.writeheader()
        for j in js:
            w.writerow(asdict(j))


def compare_groups(per_sample: dict[str, dict], groups: dict[str, str], case: str,
                   control: str) -> list[dict]:
    """Animal-level comparison of junction features (replicates.permutation_p, Welch)."""
    from scipy import stats

    from .replicates import permutation_p
    rows = []
    feats = {"pct_mh_4plus": lambda s: s["microhomology_classes_pct"].get(">=4 nt", 0.0),
             "pct_blunt": lambda s: s["microhomology_classes_pct"].get("blunt", 0.0),
             "pct_insertion": lambda s: s["microhomology_classes_pct"].get("insertion", 0.0),
             "mean_microhomology": lambda s: s["mean_microhomology"],
             "mutations_near_join_per_1e3": lambda s: s["mutations_near_join_per_1e3"]}
    for f, get in feats.items():
        a = [get(s) for n, s in per_sample.items() if groups.get(n) == case and get(s) is not None]
        b = [get(s) for n, s in per_sample.items() if groups.get(n) == control and get(s) is not None]
        row = {"feature": f, "group": case, "control": control,
               "mean_group": round(sum(a) / len(a), 3) if a else None,
               "mean_control": round(sum(b) / len(b), 3) if b else None}
        if len(a) >= 2 and len(b) >= 2:
            row["p_welch"] = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
            row["p_permutation"] = permutation_p(a, b)
        rows.append(row)
    return rows


def main(argv=None):
    import argparse
    import json
    import os
    ap = argparse.ArgumentParser(prog="nanoshm-switch", description=__doc__.split("\n\n")[0])
    ap.add_argument("--reference", required=True, help="FASTA, one record per switch region")
    ap.add_argument("--donor", help="the donor record's name (default: the first)")
    ap.add_argument("--sample", action="append", required=True,
                    help="NAME=READS[:GROUP], repeat for each sample")
    ap.add_argument("--compare", help="CASE,CONTROL groups to compare animal by animal")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    per, groups = {}, {}
    for spec in a.sample:
        name, rest = spec.split("=", 1)
        path, _, grp = rest.partition(":")
        r = analyse(path, a.reference, a.donor)
        write_csv(r["junctions"], os.path.join(a.out, f"{name}_junctions.csv"))
        per[name] = r["summary"]
        groups[name] = grp or name
        print(name, json.dumps(r["summary"]))
    json.dump({"samples": per, "groups": groups}, open(os.path.join(a.out, "switch.json"), "w"),
              indent=1)
    if a.compare:
        case, control = a.compare.split(",")
        for row in compare_groups(per, groups, case, control):
            print(row)


if __name__ == "__main__":
    main()
