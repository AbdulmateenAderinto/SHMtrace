"""Data for the figures that show how SHMtrace makes a call and what the data look like.

    python docs/calibration/crm/intuition_data.py      -> out/intuition_<dataset>.json

Each dataset is prepared and called exactly as in benchmark_real.py (recommended mode, the dataset's depth,
every AID-active library against all its AID-null controls, every control against the other controls).
Saved per dataset:
  calls           every call in each library: position, ref, alt, supporting reads, depth, allele fraction
  candidates      each AID-active library's tested candidates (p value, kept) and the FDR threshold
  decoys          the controls' candidate p values and window bases, as the control FDR uses them
  profile         per-position mismatch rate (base quality >= 20, germline sites excluded) in every library
  quality         per reported base quality: control bases and mismatches (germline sites excluded)
J_H4 only:
  anatomy         a called hotspot and the control site with the most recurrent error: alt reads by strand,
                  depth, alt base qualities in the wild-type library and each control, its p value and the
                  threshold
  raster          for one wild-type and one Aicda-/- library: each read's clone and its mismatch positions
                  (quality >= 20, germline sites excluded), called sites marked
"""
import json
import sys
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path

from shmtrace import batch, clones, null_model, pipeline
from shmtrace.config import load_preset
from shmtrace.io_utils import read_fasta

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402

OUT = Path(__file__).parent / "out"
QMIN = 20
RASTER = ("WT3", "KO1")
# the library the anatomy is shown in: Sμ WT3 carries the recurrent site at the controls' own level
ANATOMY_CASE = {"su_lps": "WT3"}


def mismatches(aligned, ref, lo, hi, qmin=0):
    """[(read index, pos0, base, quality, strand)] for every mismatch inside the window."""
    out = []
    for i, a in enumerate(aligned):
        for rpos, base, q in a.pairs:
            if base and lo - 1 <= rpos <= hi - 1 and base != ref[rpos] and q >= qmin:
                out.append((i, rpos, base, q, a.strand))
    return out


def depth_by_pos(aligned, lo, hi, qmin=0):
    d = Counter()
    for a in aligned:
        for rpos, base, q in a.pairs:
            if base and lo - 1 <= rpos <= hi - 1 and q >= qmin:
                d[rpos] += 1
    return d


def site_view(aligned, ref, pos0, alt):
    """Alt reads by strand, depth by strand and the alt bases' qualities at one site."""
    v = {"alt_plus": 0, "alt_minus": 0, "depth_plus": 0, "depth_minus": 0, "alt_quals": []}
    for a in aligned:
        for rpos, base, q in a.pairs:
            if rpos != pos0 or not base:
                continue
            s = "plus" if a.strand > 0 else "minus"
            v[f"depth_{s}"] += 1
            if base == alt:
                v[f"alt_{s}"] += 1
                v["alt_quals"].append(int(q))
    return v


def run(key):
    ds = DATASETS[key]
    p = replace(load_preset(ds["preset"]), **ds.get("overrides", {}))
    p = replace(p, **pipeline.mode_overrides("recommended"))
    ref_name, ref = read_fasta(ds["reference"])[0]
    lo, hi = pipeline._window(p, ref)
    if p.end_trim_bp:
        lo, hi = max(lo, p.end_trim_bp + 1), min(hi, len(ref) - p.end_trim_bp)
    samples = {n: s["fastq"] for n, s in ds["samples"].items()}
    controls = sorted(n for n, s in ds["samples"].items() if s["control"])
    b = batch.prepare_batch(samples, ds["reference"], ds.get("contaminants", []), p, depth=ds["depth"],
                            controls=set(controls))
    nulls, res = {}, {}
    for name in samples:
        peers = tuple(c for c in controls if c != name)
        if peers not in nulls:
            nulls[peers] = null_model.fit_aligned([b["reads"][c] for c in peers], ref, p, lo, hi,
                                                  calibrate=p.error_model == "qcal", ref_name=ref_name)[0]
        res[name] = pipeline.analyse(b["reads"][name], ds["reference"], p, b["stats"][name],
                                     null_rates=nulls[peers], depth=b["depth"])
        print(key, name, len(res[name]["variants"]), flush=True)
    aligned = {n: r["_aligned"] for n, r in res.items()}
    # germline differences: >= 50% of a library's reads, removed as SHMtrace does
    germ = set()
    for al in aligned.values():
        d, mm = depth_by_pos(al, lo, hi), Counter((x[1], x[2]) for x in mismatches(al, ref, lo, hi))
        germ |= {pa for pa, c in mm.items() if d[pa[0]] and c / d[pa[0]] >= 0.5}
    germ_pos = {pa[0] for pa in germ}
    out = {"dataset": key, "window": [lo, hi], "depth": b["depth"], "controls": controls,
           "cases": [n for n in samples if n not in controls], "calls": {}, "candidates": {}, "profile": {}}
    for n, r in res.items():
        out["calls"][n] = [[v.pos, v.ref, v.alt, v.supporting_reads, v.depth, round(v.frequency, 6)]
                           for v in r["variants"]]
        d = depth_by_pos(aligned[n], lo, hi, QMIN)
        mm = Counter(x[1] for x in mismatches(aligned[n], ref, lo, hi, QMIN) if x[1] not in germ_pos)
        out["profile"][n] = {str(pos + 1): [mm[pos], d[pos]] for pos in range(lo - 1, hi) if d[pos]}
        if n not in controls:
            s = r["summary"]["calling"]
            out["candidates"][n] = {"p": [c[1] for c in s.get("_candidates", [])],
                                    "kept": [c[2] for c in s.get("_candidates", [])],
                                    "frequency": [c[0] for c in s.get("_candidates", [])],
                                    "threshold": s.get("control_fdr_threshold"),
                                    "bases": null_model.window_bases(aligned[n], lo, hi)}
    full = nulls[tuple(controls)]
    out["decoys"] = {c: [list(pv), bases] for c, (pv, bases) in full.decoys.items()}
    # reported quality against observed error in the controls
    qb, qm = Counter(), Counter()
    for c in controls:
        for a in aligned[c]:
            for rpos, base, q in a.pairs:
                if base and lo - 1 <= rpos <= hi - 1 and rpos not in germ_pos:
                    qb[int(q)] += 1
                    qm[int(q)] += base != ref[rpos]
    out["quality"] = {str(q): [qm[q], qb[q]] for q in sorted(qb)}
    out["anatomy"] = anatomy(key, res, aligned, ref, controls, lo, hi, germ_pos)
    if key == "jh4":
        out["raster"] = raster(res, aligned, ref, p, lo, hi, germ_pos)
    OUT.mkdir(exist_ok=True)
    json.dump(out, open(OUT / f"intuition_{key}.json", "w"))
    print(OUT / f"intuition_{key}.json")


def geneious_recurrent(key, controls):
    """Sites the Geneious rule called in AID-null libraries and SHMtrace called nowhere (benchmark_real.py),
    most libraries first: {(pos0, alt): (controls called, cases called)}."""
    f = OUT / f"benchmark_real_{key}.json"
    if not f.exists():
        return {}
    rows = json.load(open(f))["rows"]
    rec = {tuple(s) for r in rows if r["mode"] == "recommended" for s in r["sites"]}
    out = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["mode"] == "geneious":
            for pos, rb, alt in (tuple(s) for s in r["sites"]):
                if (pos, rb, alt) not in rec:
                    out[(pos - 1, alt)][0 if r["sample"] in controls else 1] += 1
    return {k: tuple(v) for k, v in out.items() if v[0]}


def anatomy(key, res, aligned, ref, controls, lo, hi, germ_pos):
    """A true call and a recurrent error, both tested in one AID-active library."""
    rec = geneious_recurrent(key, controls)
    cases = [n for n in res if n not in controls]
    # the error: the site the Geneious rule called in the most controls and cases, SHMtrace in none
    err = max(rec, key=lambda k: (rec[k][0], rec[k][1])) if rec else None
    # the case: the library with the most quality >= 20 alt reads at that site, else the most mutated
    hq = {n: Counter((pos, base) for _i, pos, base, _q, _s in mismatches(aligned[n], ref, lo, hi, QMIN))
          for n in res}
    case = ANATOMY_CASE.get(key) or (max(cases, key=lambda n: hq[n][err]) if err else
                                     max(cases, key=lambda n: len(res[n]["variants"])))
    cand = {(pos - 1, alt): pv for pos, alt, pv in res[case]["summary"]["calling"]["_candidate_keys"]}
    called = {(v.pos - 1, v.alt): v for v in res[case]["variants"]}
    pooled_alt = sum((hq[c] for c in controls), Counter())
    # the true call: a WRC/GYW hotspot carried by 5-20 reads with no control support at all
    def wrc(k):      # C in a WRC motif on the top strand
        return ref[k[0]] == "C" and ref[k[0] - 2:k[0]] in ("AA", "AG", "TA", "TG")
    hot = ([k for k, v in called.items() if 5 <= v.supporting_reads <= 20 and pooled_alt[k] == 0 and wrc(k)]
           or [k for k, v in called.items() if pooled_alt[k] == 0] or list(called))
    true = max(hot, key=lambda k: called[k].supporting_reads) if hot else None
    sites = {}
    for label, site in (("true", true), ("error", err)):
        if site is None:
            continue
        pos0, alt = site
        sites[label] = {"pos": pos0 + 1, "ref": ref[pos0], "alt": alt, "context": ref[pos0 - 2:pos0 + 3],
                        "p": cand.get((pos0, alt), 1.0), "called": (pos0, alt) in called,
                        "geneious_calls": rec.get((pos0, alt)),
                        "libraries": {n: site_view(aligned[n], ref, pos0, alt) for n in [case, *controls]}}
    s = res[case]["summary"]["calling"]
    return {"case": case, "threshold": s.get("control_fdr_threshold"), "sites": sites}


def raster(res, aligned, ref, p, lo, hi, germ_pos):
    out = {}
    for n in RASTER:
        called = {(v.pos - 1, v.alt) for v in res[n]["variants"]}
        cl = clones.call_clones(aligned[n], res[n]["variants"], p, lo, hi) or {}
        clone_of = cl.get("clone_of", {})
        reads = defaultdict(list)
        for i, pos, base, _q, _s in mismatches(aligned[n], ref, lo, hi, QMIN):
            if pos in germ_pos:          # germline differences, in nearly every read
                continue
            reads[i].append([pos + 1, base, int((pos, base) in called)])
        out[n] = [{"clone": clone_of.get(a.qname), "muts": reads.get(i, []), "start": a.r_st + 1,
                   "end": a.r_en} for i, a in enumerate(aligned[n])]
    return out


if __name__ == "__main__":
    for k in sys.argv[1:] or ("jh4", "n114", "su_lps", "pim1"):
        run(k)
