"""Does the mutation footprint follow the guide? (Ma et al. 2016, SRP075898, ABL exon 6)

    python crm/tam_footprint.py        -> crm/out/tam/footprint.json

Five libraries target the same ABL exon-6 amplicon with five different sgRNAs (#4, #6, #7, #9, #15). That
makes the arm self-controlled for localisation in a way the reporter arm is not: the question is not whether
an editor mutates, but whether the position it mutates moves when the guide moves. Each library is the
negative control for the other four at their own peaks, so no no-editor library is needed and none exists
here.

The reference is the consensus of all five libraries pooled. Pooling is deliberate: a reference taken from
one library would carry that library's own edits, and pooling five differently-targeted libraries means no
single guide's footprint can dominate any position.

These are per-position mismatch rates, not mutation calls. They contain sequencing error, which is why only
the *shape* along the amplicon is interpreted, and why the comparison is between libraries at the same
position rather than between positions within a library.
"""
import collections
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
REPO = CRM.parents[2]
sys.path.insert(0, str(REPO / "src"))
from shmtrace import io_utils  # noqa: E402

DATA = Path.home() / "Nanoshm/public_data/SRP075898"
OUT = CRM / "out/tam"
GUIDES = [("SRR4254122", "sgRNA_4"), ("SRR4254123", "sgRNA_6"), ("SRR4254124", "sgRNA_7"),
          ("SRR4254125", "sgRNA_9"), ("SRR4254127", "sgRNA_15")]
MIN_Q, LIMIT, MIN_LEN = 30, 20_000, 120
EDGE = 15                       # positions trimmed from each end before looking for a peak


def merged(run, limit=LIMIT):
    """Merged pairs of one library, as the pipeline would merge them."""
    n = 0
    for _nm, seq, qual in io_utils.merge_pairs(str(DATA / f"{run}_1.fastq.gz"),
                                               str(DATA / f"{run}_2.fastq.gz")):
        yield seq, qual
        n += 1
        if n >= limit:
            return


def anchors(top=4, anchor_len=25, min_share=0.05):
    """The distinct amplicon starts in this arm.

    These libraries are multiplexed, like the reporter arm: several amplicons share each file. Anchoring on
    the single commonest start analyses one amplicon and silently drops the guides that target the others,
    which is exactly the mistake this function exists to prevent.
    """
    c = collections.Counter()
    for run, _g in GUIDES:
        for i, (seq, _q) in enumerate(merged(run, 4000)):
            if len(seq) >= MIN_LEN:
                c[seq[:anchor_len]] += 1
            if i > 3000:
                break
    total = sum(c.values())
    return [a for a, n in c.most_common(top) if n >= min_share * total]


def build_reference(anchor, per_library=6000):
    """Consensus of this amplicon, pooled over all five libraries so no one guide's edits dominate."""
    cols = collections.defaultdict(collections.Counter)
    used = 0
    for run, _g in GUIDES:
        n = 0
        for seq, _q in merged(run, 30_000):
            if not seq.startswith(anchor):
                continue
            for i, b in enumerate(seq):
                cols[i][b] += 1
            n += 1
            used += 1
            if n >= per_library:
                break
    if not cols:
        return None, 0, 0.0
    length = max(i for i in cols if cols[i].total() >= 0.5 * cols[0].total()) + 1
    ref = "".join(cols[i].most_common(1)[0][0] for i in range(length))
    weakest = min(cols[i].most_common(1)[0][1] / cols[i].total() for i in range(length))
    return ref, used, weakest


def profile(run, anchor, ref):
    """Mismatches and covered bases per reference position, over reads sharing the anchor."""
    mism, cov, n = collections.Counter(), collections.Counter(), 0
    for seq, qual in merged(run):
        if not seq.startswith(anchor):
            continue
        n += 1
        for i in range(min(len(seq), len(ref))):
            if ord(qual[i]) - 33 < MIN_Q:
                continue
            cov[i] += 1
            if seq[i] != ref[i]:
                mism[i] += 1
    return mism, cov, n


def main():
    missing = [r for r, _g in GUIDES if not (DATA / f"{r}_1.fastq.gz").exists()]
    if missing:
        sys.exit(f"missing reads for {missing}")
    OUT.mkdir(parents=True, exist_ok=True)
    found = anchors()
    print(f"{len(found)} amplicons in this arm", flush=True)
    out = {"dataset": "SRP075898 (Ma et al. 2016), ABL exon-6 footprint arm",
           "note": "per-position mismatch rates, not mutation calls; each library is the control for the "
                   "others at their own peaks. The arm is multiplexed, so every amplicon is analysed.",
           "amplicons": {}}
    best = {}
    for ai, anchor in enumerate(found):
        ref, used, weakest = build_reference(anchor)
        if ref is None or len(ref) < 100:
            continue
        name = f"amp{ai}_{anchor[:10]}"
        print(f"\n{name}: {len(ref)} bp, pooled from {used} reads, weakest consensus "
              f"{100 * weakest:.1f}%", flush=True)
        profiles, reads = {}, {}
        for run, guide in GUIDES:
            mism, cov, n = profile(run, anchor, ref)
            profiles[guide] = [1000 * mism[i] / cov[i] if cov[i] >= 500 else None for i in range(len(ref))]
            reads[guide] = n
        rows = {}
        for guide, rate in profiles.items():
            usable = [i for i, r in enumerate(rate)
                      if r is not None and EDGE <= i < len(ref) - EDGE]
            if not usable:
                continue
            peak = max(usable, key=lambda i: rate[i])
            others = [profiles[g][peak] for g in profiles if g != guide and profiles[g][peak] is not None]
            med = sorted(others)[len(others) // 2] if others else None
            rows[guide] = {"reads_on_anchor": reads[guide], "peak_position": peak,
                           "peak_rate_per_1e3": round(rate[peak], 2),
                           "others_at_this_peak_median": round(med, 2) if med else None,
                           "fold_over_others": round(rate[peak] / med, 1) if med else None}
            f = rows[guide]["fold_over_others"] or 0
            if f > best.get(guide, (0,))[0]:
                best[guide] = (f, name, peak, rows[guide]["peak_rate_per_1e3"])
        out["amplicons"][name] = {"bp": len(ref), "weakest_consensus_pct": round(100 * weakest, 1),
                                  "libraries": rows,
                                  "profiles_per_1e3": {g: [None if r is None else round(r, 3) for r in p]
                                                       for g, p in profiles.items()}}
        print(f"   {'guide':10s} {'reads':>7} {'peak':>6} {'rate/1e3':>9} {'others':>8} {'fold':>6}",
              flush=True)
        for guide, v in rows.items():
            print(f"   {guide:10s} {v['reads_on_anchor']:>7} {v['peak_position']:>6} "
                  f"{v['peak_rate_per_1e3']:>9} {str(v['others_at_this_peak_median']):>8} "
                  f"{str(v['fold_over_others']):>6}", flush=True)
    out["best_per_guide"] = {g: {"fold": f, "amplicon": a_, "position": p, "rate_per_1e3": r}
                             for g, (f, a_, p, r) in best.items()}
    json.dump(out, open(OUT / "footprint.json", "w"), indent=1)
    print("\nbest localised excess per guide, across amplicons:", flush=True)
    for g, v in sorted(best.items()):
        print(f"   {g:10s} {v[0]:>6}x at {v[1]} position {v[2]} ({v[3]} per 1,000)", flush=True)


if __name__ == "__main__":
    main()
