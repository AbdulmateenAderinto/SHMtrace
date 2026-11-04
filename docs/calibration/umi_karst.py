"""SHMtrace's UMI consensus on real nanopore UMI data with known truth (Karst et al. 2021).

    python docs/calibration/umi_karst.py [DATA_DIR]   (default ~/Nanoshm/public_data/karst_umi)

Data (downloaded from github.com/SorenKarst/longread_umi, GPL-3; not redistributed here):
  ont_r10_zymo_rrna.fq.gz         4,259 MinION R10 reads (guppy 3.4.4 hac) of ~4.4-kb rRNA
                                  operons from the ZymoBIOMICS mock community, each molecule
                                  tagged with an 18-nt UMI at both ends (NNNYRNNNYRNNNYRNNN)
                                  just outside the primers 8F (AGRGTTYGATYMTGGCTCAG) and
                                  2490R (CGACATCGAGGTGCCAAAC)
  zymo-ref-uniq_2019-10-28.fa     the 43 distinct operon sequences of the community: the truth

Every read is aligned to the 43 operons (minimap2), its UMIs are read next to the primers
(shmtrace.umi.extract, IUPAC-aware anchors), reads are grouped by UMI (umi.group) and
each group of 3 or more becomes one consensus (umi._consensus) on the operon most of its
reads align to. Error rates are measured against the operon each sequence matches best:
mismatches (and deletions) per aligned base, for raw reads and for consensus molecules.
"""
import gzip
import json
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path

import mappy as mp

from shmtrace import umi
from shmtrace.align import AlignedRead, _aligned_pairs
from shmtrace.config import Params
from shmtrace.io_utils import phred

DATA = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home() / "Nanoshm/public_data/karst_umi")
OUT = Path(__file__).parent / "umi_karst"
F_PRIMER = "AGRGTTYGATYMTGGCTCAG"
R_PRIMER_RC = mp.revcomp("CGACATCGAGGTGCCAAAC")


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split(";")[0]
            out[name] = []
        elif name:
            out[name].append(line.upper())
    return {k: "".join(v) for k, v in out.items()}


def read_fq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, _p, q = fh.readline().strip(), fh.readline(), fh.readline().strip()
            yield h[1:].split()[0], s, q


def to_aligned(name, seq, qual, hit):
    if hit.strand == -1:
        wseq, wq, q0 = mp.revcomp(seq), phred(qual)[::-1], len(seq) - hit.q_en
    else:
        wseq, wq, q0 = seq, phred(qual), hit.q_st
    pairs = _aligned_pairs(hit, wseq, wq, q0)
    return AlignedRead(name, hit.ctg, hit.r_st, hit.r_en, pairs, wseq, hit.strand, [], q0,
                       q0 + (hit.q_en - hit.q_st))


def err_of(seq, aligner):
    """(edit distance, aligned length) against the best-matching operon."""
    h = next((x for x in aligner.map(seq) if x.is_primary), None)
    return (h.NM, h.blen, h.ctg) if h else (None, None, None)


def grouping_stats(gid, aligned, min_reads):
    """Groups of `min_reads` or more and how pure they are: the share of a group's
    reads on its commonest operon, and the share of groups mixing two operon
    sequences that differ by more than near-identical copies of one genome do."""
    groups = defaultdict(list)
    for name, g in gid.items():
        groups[g].append(aligned[name][0].ctg)
    big = [ctgs for ctgs in groups.values() if len(ctgs) >= min_reads]
    pur = [Counter(c).most_common(1)[0][1] / len(c) for c in big]
    species = [{x.split("_")[0] for x in c} for c in big]
    return {"umis": len(groups), "groups_3plus_reads": len(big),
            "reads_in_groups_3plus": sum(len(c) for c in big),
            "group_purity_mean": round(statistics.mean(pur), 4) if pur else None,
            "groups_mixing_species_pct": round(100 * sum(len(s) > 1 for s in species)
                                               / len(big), 2) if big else None}


def main():
    refs = read_fasta(DATA / "zymo-ref-uniq_2019-10-28.fa")
    fa = DATA / "refs_clean_names.fa"           # kept with the data, not in the repository
    OUT.mkdir(parents=True, exist_ok=True)
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in refs.items()))
    al = mp.Aligner(str(fa), preset="map-ont")
    one = {k: mp.Aligner(seq=v, preset="map-ont") for k, v in refs.items()}
    p = replace(Params(), umi="flank", umi_ends="both", umi_length=18,
                umi_left_anchor=F_PRIMER, umi_right_anchor=R_PRIMER_RC, umi_anchor_edits=6,
                umi_max_edits=3, umi_min_reads=3)
    raw = {}
    aligned = {}
    for name, seq, qual in read_fq(DATA / "ont_r10_zymo_rrna.fq.gz"):
        raw[name] = (seq, qual)
        h = next((x for x in al.map(seq) if x.is_primary), None)
        if h is None or h.r_en - h.r_st < 3000:
            continue
        aligned[name] = (to_aligned(name, seq, qual, h), h)
    keys = {}
    for name, (r, _h) in aligned.items():
        u = umi.extract(r, refs[r.ctg], p)
        if u and "N" not in u:
            keys[name] = u
    raw_err = [(h.NM, h.blen) for _r, h in aligned.values()]
    raw_rate = sum(x[0] for x in raw_err) / sum(x[1] for x in raw_err)

    def molecules(gid):
        groups = defaultdict(list)
        for name, g in gid.items():
            groups[g].append(name)
        cons_rows = []
        for g, names in groups.items():
            if len(names) < p.umi_min_reads:
                continue
            ctg = Counter(aligned[n][0].ctg for n in names).most_common(1)[0][0]
            purity = sum(aligned[n][0].ctg == ctg for n in names) / len(names)
            members = []
            for n in names:
                r, h = aligned[n]
                if r.ctg != ctg:             # realign to the group's operon
                    seq, qual = raw[n]
                    h2 = next((x for x in one[ctg].map(seq) if x.is_primary), None)
                    if h2 is None:
                        continue
                    r = to_aligned(n, seq, qual, h2)
                    r.ctg = ctg
                members.append(r)
            c = umi._consensus(members, refs[ctg], f"umi{g}")
            if c is None:
                continue
            # consensus as a sequence (deleted columns left out), measured like a raw read
            cseq = "".join(b for _rp, b, _q in c.pairs)
            nm, blen, best = err_of(cseq, al)
            subs = sum(1 for rp, b, _q in c.pairs if b != refs[ctg][rp])
            cons_rows.append({"reads": len(names), "operon": ctg, "best_operon": best,
                              "purity": purity, "nm": nm, "blen": blen,
                              "subs_vs_group_operon": subs, "span": c.r_en - c.r_st})
        return groups, cons_rows

    # the directional step at +2 (the default for dual UMIs) and +1 edits beyond
    # umi_max_edits, and off
    directional = {}
    for d in (2, 1, 0):
        gid = umi.group(keys, 2 * p.umi_max_edits, d)
        _g, rows = molecules(gid)
        st = grouping_stats(gid, aligned, p.umi_min_reads)
        ok = [r for r in rows if r["nm"] is not None]
        st["consensus_error_pct"] = round(100 * sum(r["nm"] for r in ok)
                                          / sum(r["blen"] for r in ok), 4)
        st["consensus_substitution_error_pct"] = round(
            100 * sum(r["subs_vs_group_operon"] for r in ok) / sum(r["span"] for r in ok), 4)
        st["molecules_3_4_reads"] = sum(r["reads"] <= 4 for r in ok)
        directional[f"+{d}" if d else "off"] = st
    groups, cons_rows = molecules(umi.group(keys, 2 * p.umi_max_edits, umi.directional(p)))
    cons = [r for r in cons_rows if r["nm"] is not None]
    cons_rate = sum(r["nm"] for r in cons) / sum(r["blen"] for r in cons) if cons else None
    by_depth = {}
    for lo, hi in ((3, 4), (5, 9), (10, 19), (20, 10 ** 6)):
        rs = [r for r in cons if lo <= r["reads"] <= hi]
        if rs:
            by_depth[f"{lo}-{hi if hi < 10 ** 6 else '+'}"] = {
                "molecules": len(rs),
                "error_pct": round(100 * sum(r["nm"] for r in rs) / sum(r["blen"] for r in rs), 4),
                "perfect_pct": round(100 * sum(r["nm"] == 0 for r in rs) / len(rs), 1),
                "substitution_error_pct": round(100 * sum(r["subs_vs_group_operon"] for r in rs)
                                                / sum(r["span"] for r in rs), 4),
                "no_substitution_pct": round(100 * sum(r["subs_vs_group_operon"] == 0 for r in rs)
                                             / len(rs), 1)}
    summary = {
        "reads": len(raw), "aligned": len(aligned), "with_umi": len(keys),
        "umi_groups": len(groups), "molecules_3plus_reads": len(cons),
        "median_reads_per_molecule": statistics.median([r["reads"] for r in cons]) if cons else None,
        "group_purity_mean": round(statistics.mean(r["purity"] for r in cons), 4) if cons else None,
        "raw_read_error_pct": round(100 * raw_rate, 3),
        "consensus_error_pct": round(100 * cons_rate, 4) if cons_rate is not None else None,
        "error_reduction_fold": round(raw_rate / cons_rate, 1) if cons_rate else None,
        "by_reads_per_molecule": by_depth,
        "consensus_matching_group_operon_pct": round(
            100 * sum(r["best_operon"] == r["operon"] for r in cons) / len(cons), 1) if cons else None,
        "directional": directional,
    }
    json.dump({"summary": summary, "molecules": cons_rows}, open(OUT / "results.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
