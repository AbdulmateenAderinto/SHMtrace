"""Major issue 5: does acceptor-breakpoint ambiguity change the switch endpoints?

    python crm/switch_ambiguity_sensitivity.py   -> crm/out/switch_pacbio/ambiguity_sensitivity.json

`switch.py`'s own ambiguous flag reads 0% on these libraries because minimap2 returns no secondaries here;
re-aligning the acceptor segment directly measures 19-47% (Figure S16D). That is a confirmed defect in the
flag. The question the referee asked is different and is answered here: if every junction whose acceptor
segment has more than one near-equal placement is dropped, do the comparisons in Figure S16C move?

Each sample is re-analysed, its ambiguous unique junctions identified by the same direct probe, and the
per-sample summary recomputed from the junctions that remain. The group comparisons are then re-run on that
reduced set and placed beside the published ones.
"""
import json
import sys
from pathlib import Path

import mappy as mp

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(CRM))
from shmtrace import switch  # noqa: E402
from shmtrace.io_utils import read_fasta  # noqa: E402

import switch_pacbio as SP  # noqa: E402


def ambiguous_keys(junctions, seqs, acc_ref, seg=60, slack=5):
    """The junction keys whose acceptor segment has more than one near-equal placement."""
    al = mp.Aligner(seq=acc_ref, preset="sr", best_n=50)
    bad = set()
    for j in junctions:
        s = seqs.get(j.read)
        if not s:
            continue
        seg_seq = s[0][s[1]:s[1] + seg]
        if len(seg_seq) < seg:
            continue
        hits = sorted((h.mlen for h in al.map(seg_seq)), reverse=True)
        if hits and sum(1 for m in hits if m >= hits[0] - slack) > 1:
            bad.add(j.key)
    return bad


def main():
    ref_path = SP.REFDIR / "mouse_Smu_Sg1.fasta"
    refs = {n.split()[0]: s.upper() for n, s in read_fasta(str(ref_path))}
    donor = "Smu"
    kept, groups, dropped = {}, {}, {}
    for run, geno, rep in SP.SAMPLES:
        name = f"{geno}_{rep}"
        res, junctions = SP.analyse_sample(run, ref_path, refs, donor)
        # rebuild the orientation map the probe needs, exactly as analyse_sample does
        al = mp.Aligner(str(ref_path), preset=SP.PRESET, best_n=10)
        seqs = {}
        for j in junctions:
            pass
        # analyse_sample already measured the share; re-derive the keys from the same inputs
        import statistics  # noqa: F401
        from shmtrace.io_utils import read_fastq
        for nm, seq, _q in read_fastq(str(SP.DATA / f"{run}_subreads.fastq.gz")):
            d0 = max((h for h in al.map(seq) if h.ctg == donor), key=lambda h: h.mlen, default=None)
            oriented = mp.revcomp(seq) if (d0 and d0.strand == -1) else seq
            acc = [h for h in al.map(oriented) if h.ctg != donor and h.strand == 1]
            if acc:
                seqs[nm] = (oriented, max(acc, key=lambda h: h.mlen).q_st)
        sw = [j for j in junctions if j.kind == "switch"]
        bad = ambiguous_keys(list(switch.unique_junctions(sw).values()) if isinstance(
            switch.unique_junctions(sw), dict) else switch.unique_junctions(sw), seqs, refs["Sg1"])
        clean = [j for j in junctions if getattr(j, "key", None) not in bad]
        kept[name] = switch.summarise(clean)
        groups[name] = geno
        dropped[name] = {"junctions_before": len(junctions), "junctions_after": len(clean),
                         "ambiguous_keys_dropped": len(bad),
                         "published_ambiguous_pct": res["acceptor_ambiguity_direct"]["ambiguous_pct"]}
        print(f"  {name:12s} {len(junctions):>5} -> {len(clean):>5} junctions "
              f"({len(bad)} ambiguous keys dropped)", flush=True)
    comps = []
    for case in ("Shld1", "Xlf", "Shld1_Xlf"):
        comps += switch.compare_groups(kept, groups, case, "WT")
    pub = {(c["feature"], c["group"]): c
           for c in json.load(open(SP.OUT / "switch_pacbio.json"))["comparisons"]}
    rows = []
    for c in comps:
        p = pub.get((c["feature"], c["group"]))
        rows.append({**c, "published_p_welch": p["p_welch"] if p else None,
                     "published_mean_group": p["mean_group"] if p else None,
                     "published_mean_control": p["mean_control"] if p else None})
    out = {"what": "group comparisons after dropping every junction whose acceptor segment has more than "
                   "one near-equal placement", "per_sample": dropped, "comparisons_clean_vs_published": rows}
    (SP.OUT / "ambiguity_sensitivity.json").write_text(json.dumps(out, indent=1, default=str))
    print(f"\n{'feature':28s} {'group':10s} {'p clean':>9s} {'p published':>12s}  same side of 0.05?")
    for r in rows:
        a, b = r["p_welch"], r["published_p_welch"]
        print(f"  {r['feature']:26s} {r['group']:10s} {a:>9.4f} {b:>12.4f}  "
              f"{'yes' if (a < 0.05) == (b < 0.05) else 'NO'}")


if __name__ == "__main__":
    main()
