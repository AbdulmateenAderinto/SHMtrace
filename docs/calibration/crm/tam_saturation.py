"""Does the number of mutated sites follow the number of guides? (Ma et al. 2016, saturation arm)

    python crm/tam_saturation.py        -> crm/out/tam/saturation.json

Four libraries target the same locus with 3, 5, 7 and 10 sgRNAs. If the calls track the guides rather than
the editor, the count of sharply localised sites should rise with the number of guides. The arm is
self-controlled in the same way as the footprint arm: no no-editor library exists, and none is needed,
because at any one position the other three libraries say what that position looks like without that guide.

A site is counted as localised when a library's rate there is at least `FOLD` times the lowest of the
*other three* libraries at that position, and clears an absolute floor. Leaving the library itself out of its
own baseline matters: scored against the lowest of all four, the library with the fewest guides is usually
the minimum, so it would be compared against itself and score zero by construction.

These are per-position mismatch rates, not mutation calls: they still contain sequencing error, which is why
only the counts of localised sites are compared, and only between libraries at the same positions.

Run this in the foreground; backgrounded it dies part-way (see NEXT_STEPS.md).
"""
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
sys.path.insert(0, str(CRM))
import tam_footprint as F  # noqa: E402

OUT = F.OUT
LIBS = [("SRR4254105", 3), ("SRR4254104", 5), ("SRR4254103", 7), ("SRR4254102", 10)]
FOLD, FLOOR = 5.0, 2.0      # times the lowest library at that position, and per 1,000 bases


def main():
    missing = [r for r, _n in LIBS if not (F.DATA / f"{r}_1.fastq.gz").exists()]
    if missing:
        sys.exit(f"missing reads for {missing}")
    OUT.mkdir(parents=True, exist_ok=True)
    F.GUIDES = [(r, f"{n}_sgRNAs") for r, n in LIBS]        # reuse the footprint machinery on this arm
    found = F.anchors()
    print(f"{len(found)} amplicons in the saturation arm", flush=True)
    out = {"dataset": "SRP075898 (Ma et al. 2016), saturation arm",
           "rule": f"a site counts when a library is >= {FOLD}x the lowest of the OTHER libraries at that "
                   f"position and >= {FLOOR} per 1,000 (leave-one-out, so no library is its own baseline)",
           "note": "per-position mismatch rates, not mutation calls",
           "amplicons": {}}
    totals = {f"{n}_sgRNAs": 0 for _r, n in LIBS}
    for ai, anchor in enumerate(found):
        ref, used, weakest = F.build_reference(anchor)
        if ref is None or len(ref) < 100:
            continue
        name = f"amp{ai}_{anchor[:10]}"
        rates = {}
        for run, label in F.GUIDES:
            mism, cov, n = F.profile(run, anchor, ref)
            rates[label] = [1000 * mism[i] / cov[i] if cov[i] >= 500 else None for i in range(len(ref))]
        hot = {label: [] for label in rates}
        for i in range(F.EDGE, len(ref) - F.EDGE):
            vals = {label: rates[label][i] for label in rates}
            if any(v is None for v in vals.values()):
                continue
            for label in rates:
                others = [v for k, v in vals.items() if k != label]
                lo = max(min(others), 0.05)                # a floor, so a zero cannot make the fold infinite
                if vals[label] >= FLOOR and vals[label] >= FOLD * lo:
                    hot[label].append(i)
        out["amplicons"][name] = {"bp": len(ref), "weakest_consensus_pct": round(100 * weakest, 1),
                                  "localised_sites": {k: v for k, v in hot.items()}}
        for label, sites in hot.items():
            totals[label] += len(sites)
        print(f"\n{name}: {len(ref)} bp, weakest consensus {100 * weakest:.1f}%", flush=True)
        for label in (f"{n}_sgRNAs" for _r, n in LIBS):
            print(f"   {label:10s} localised sites: {len(hot[label])}", flush=True)
    out["totals"] = totals
    json.dump(out, open(OUT / "saturation.json", "w"), indent=1)
    print("\nlocalised sites in total, by how many guides the library carried:", flush=True)
    for _r, n in LIBS:
        print(f"   {n:>2} sgRNAs: {totals[f'{n}_sgRNAs']}", flush=True)


if __name__ == "__main__":
    main()
