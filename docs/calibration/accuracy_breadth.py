"""How accurate is SHMtrace beyond J꜀4 and N114, and how much does the simulation flatter it?

    python docs/calibration/accuracy_breadth.py cache           # merged, subsampled controls
    python docs/calibration/accuracy_breadth.py qmaps           # planted-base quality maps
    python docs/calibration/accuracy_breadth.py run breadth     # (or realism, crossctx; --jobs N)
    python docs/calibration/accuracy_breadth.py report

QUANT_ACCURACY.md measures accuracy on two nanopore amplicons. This runs the same design (lineage_sim
on each set's own mutation-free controls, one held out at a time, the others as the null) on every
public set with three or more of them:

  ldha    Illumina J꜀4 70 bp, PRJNA973112, 4 AID-null mice
  ch12    Illumina CH12F3 VDJ, PRJNA918596, 4 empty-vector replicates
  hmces   Illumina Ramos IGH VDJ, PRJNA820976, 3 AID-KO clones
  su_lps  ONT Sμ, LPS blasts, PRJNA1481264, 3 Aicda-/-
  su_ch12 ONT Sμ, CH12, PRJNA1481264, 3 Aicda-/-
  pim1    ONT Pim1 (an AID off-target), PRJNA1481264, 3 Aicda-/-
  ot1/ot2 ONT CH12 off-targets OT1 and OT2, PRJNA1481264, 3 Aicda-/- each

Arms:
  breadth   planted as in QUANT_ACCURACY (donor quality, one read per molecule), two depths, two
            draws per held-out control
  realism   at one depth, the planted side made less ideal one step at a time (lineage_sim docstring):
            qmap (mutant-base quality), pcr3 / pcr10 (3 or 10 reads per molecule, log-normal skew 1,
            first-cycle PCR errors at 1e-4 per base, not in the truth), deep (bushy trees, 4 levels),
            and all (qmap + pcr3 + deep). J꜀4 and N114 join this arm.
  crossctx  Sμ libraries from LPS blasts called against CH12 controls and the reverse: the null from
            another cell context, as when a lab's controls were not prepared beside its samples.

Reads stay outside the repository: the cache is ~/Nanoshm/public_data/accuracy_cache/. Results go to
docs/calibration/accuracy_breadth/ (runs.jsonl, appended as runs finish, so a run can be resumed;
qmaps.json; report.md).
"""
import argparse
import json
import os
import random
import statistics
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "crm"))
from datasets import DATASETS  # noqa: E402
from public_illumina import ch12, hmces, ldha  # noqa: E402

from shmtrace import lineage_sim, pipeline, project  # noqa: E402
from shmtrace.io_utils import read_fasta, read_fastq  # noqa: E402

CACHE = Path.home() / "Nanoshm/public_data/accuracy_cache"
OUT = HERE / "accuracy_breadth"
PCR = {"pcr3": {"copies": 3, "skew": 1.0, "early_error": 1e-4},
       "pcr10": {"copies": 10, "skew": 1.0, "early_error": 1e-4}}
SCEN = dict(lineage_sim.SCENARIOS, polyclonal_deep={**lineage_sim.SCENARIOS["polyclonal"], "levels": 4},
            subclonal_deep={**lineage_sim.SCENARIOS["subclonal"], "levels": 4})


def _from_crm(key, scen, depths, main):
    d = DATASETS[key]
    s = d["samples"]
    return {"region": d["preset"], "reference": d["reference"], "contaminants": d.get("contaminants", []),
            "extra": {"overrides": d.get("overrides", {})} if d.get("overrides") else {},
            "controls": {k: v["fastq"] for k, v in s.items() if v["control"]},
            "cases": {k: v["fastq"] for k, v in s.items() if not v["control"]},
            "scenarios": scen, "depths": depths, "main": main}


def _from_illumina(fn, scen, depths, main):
    pr = fn()
    return {"region": pr["region"], "reference": pr["reference"], "contaminants": pr.get("contaminants", []),
            "extra": {"platform": "illumina", "overrides": pr.get("overrides", {})},
            "controls": {s["name"]: s["path"] for s in pr["samples"] if s["control"]},
            "cases": {s["name"]: s["path"] for s in pr["samples"] if not s["control"]},
            "scenarios": scen, "depths": depths, "main": main}


SETS = {
    "ldha": _from_illumina(ldha, ["polyclonal", "low", "null"], (1000, 5000), 5000),
    "ch12": _from_illumina(ch12, ["subclonal", "low", "null"], (1000, 5000), 5000),
    "hmces": _from_illumina(hmces, ["subclonal", "low", "null"], (1000, 5000), 5000),
    "su_lps": _from_crm("su_lps", ["polyclonal", "low", "null"], (500, 2000), 2000),
    "su_ch12": _from_crm("su_ch12", ["subclonal", "low", "null"], (500, 2000), 2000),
    "pim1": _from_crm("pim1", ["low", "null"], (500, 2000), 2000),
    "ot1": _from_crm("ot1", ["low", "null"], (500, 2000), 2000),
    "ot2": _from_crm("ot2", ["low", "null"], (500, 2000), 2000),
    # the paper's two, for the realism arm only (QUANT_ACCURACY.md has their breadth)
    "jh4": _from_crm("jh4", ["polyclonal", "null"], (500,), 500),
    "n114": _from_crm("n114", ["subclonal", "null"], (500,), 500),
}
PAPER = ("jh4", "n114")
BREADTH = [k for k in SETS if k not in PAPER]
STUDY = {"ldha": "PRJNA973112", "ch12": "PRJNA918596", "hmces": "PRJNA820976", "jh4": "jh4", "n114": "n114",
         **{k: "PRJNA1481264" for k in ("su_lps", "su_ch12", "pim1", "ot1", "ot2")}}


def cached(key, name):
    return CACHE / key / f"{name}.fastq"


def cmd_cache(a):
    """Each library's reads, Illumina pairs merged, at most 2.5 x the largest depth, by a fixed seed; all of
    J꜀4's and N114's, whose raw files lose most reads to the paralogue and length filters."""
    for key, s in SETS.items():
        cap = None if key in PAPER else int(2.5 * max(s["depths"]))
        (CACHE / key).mkdir(parents=True, exist_ok=True)
        for name, path in {**s["controls"], **s["cases"]}.items():
            out = cached(key, name)
            if out.exists():
                continue
            reads = list(read_fastq(path))
            if cap and len(reads) > cap:
                reads = random.Random(17).sample(reads, cap)
            tmp = out.with_suffix(".tmp")
            with open(tmp, "w") as fh:
                for n, sq, q in reads:
                    fh.write(f"@{n}\n{sq}\n+\n{q}\n")
            tmp.rename(out)
            print(key, name, len(reads), flush=True)


def _params(s):
    return project.build_params({"region": s["region"], "mode": "recommended", **s["extra"]})


def _quals(key, n=3000):
    """Reference- and mutant-base qualities at sites real in the AID-positive libraries: an alt seen in
    5+ reads at 0.5%+ of them, and at a twentieth of that or less in the controls."""
    s = SETS[key]
    p = _params(s)
    ref = read_fasta(s["reference"])[0][1]
    libs = {}
    for name in [*s["controls"], *s["cases"]]:
        al, _ = pipeline.prepare_reads(str(cached(key, name)), s["reference"], s["contaminants"], p)
        random.Random(1).shuffle(al)
        libs[name] = (name in s["controls"], al[:n])
    # the real load: substitutions per 10^3 aligned bases in the window, each AID-positive library's excess
    # over the controls' mean; it needs no calls, so it does not lean on SHMtrace
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)

    def mm(al):
        n = m = 0
        for r in al:
            for rpos, base, _q in r.pairs:
                if base and lo <= rpos <= hi:
                    n += 1
                    m += base != ref[rpos]
        return 1e3 * m / n if n else 0.0
    ctl_mm = statistics.mean(mm(al) for c, al in libs.values() if c)
    excess = sorted(mm(al) - ctl_mm for c, al in libs.values() if not c)
    cnt, dep = defaultdict(Counter), defaultdict(Counter)
    for ctl, al in libs.values():
        for r in al:
            for rpos, base, _q in r.pairs:
                if base:
                    dep[ctl][rpos] += 1
                    cnt[ctl][(rpos, base)] += 1
    real = set()
    for (pos, b), k in cnt[False].items():
        if b == ref[pos]:
            continue
        fs = k / dep[False][pos]
        fc = cnt[True].get((pos, b), 0) / max(1, dep[True][pos])
        if k >= 5 and fs >= 0.005 and fc <= fs / 20:
            real.add((pos, b))
    sites = {pos for pos, _b in real}
    alt_q, ref_q = [], []
    for ctl, al in libs.values():
        if ctl:
            continue
        for r in al:
            for rpos, base, q in r.pairs:
                if rpos in sites and base:
                    if (rpos, base) in real:
                        alt_q.append(q)
                    elif base == ref[rpos]:
                        ref_q.append(q)
    return {"control_mismatch_x1e3": round(ctl_mm, 4), "case_excess_x1e3": [round(x, 4) for x in excess],
            "sites": len(real), "alt_n": len(alt_q), "ref_n": len(ref_q),
            "alt_mean": round(statistics.mean(alt_q), 2) if alt_q else None,
            "ref_mean": round(statistics.mean(ref_q), 2) if ref_q else None,
            "alt_lt30_pct": round(100 * sum(q < 30 for q in alt_q) / len(alt_q), 2) if alt_q else None,
            "ref_lt30_pct": round(100 * sum(q < 30 for q in ref_q) / len(ref_q), 2) if ref_q else None,
            "distinct_q": len({q for _c, al in libs.values() for r in al for _p, _b, q in r.pairs}),
            "map": lineage_sim.quality_map(ref_q, alt_q)}


def cmd_qmaps(a):
    OUT.mkdir(exist_ok=True)
    with ProcessPoolExecutor(a.jobs) as ex:
        res = dict(zip(SETS, ex.map(_quals, SETS), strict=True))
    # Files with one quality on every base (PRJNA973112 and PRJNA918596 as distributed) keep it: a planted
    # base of any other quality would stand out. A set with too few real sites borrows the pooled map of
    # the other sets of its study, whose basecaller shares its quality scale.
    for k, r in res.items():
        r["map_source"] = k
        if r["distinct_q"] <= 1:
            r["map_source"], r["map"] = "flat qualities: none", {}
        elif r["alt_n"] < 1000:
            donors = [d for d in res if STUDY[d] == STUDY[k] and d != k and res[d]["alt_n"] >= 1000
                      and res[d]["distinct_q"] > 1]
            r["map_source"] = ",".join(donors) or "none"
            r["map"] = _pool([res[d]["map"] for d in donors]) if donors else {}
    json.dump(res, open(OUT / "qmaps.json", "w"), indent=1)
    for k, r in res.items():
        print(k, {x: r[x] for x in r if x != "map"})


def _pool(maps):
    keys = sorted({q for m in maps for q in m})
    return {q: round(statistics.median([m.get(q, q) for m in maps])) for q in keys}


def _qmap(key):
    m = json.load(open(OUT / "qmaps.json"))[key]["map"]
    return {int(q): int(v) for q, v in m.items()}


LOAD_FLOOR = 0.2      # per 10^3 bp: below this the cases show no excess the reads can measure


def _calib(key):
    return json.load(open(OUT / "qmaps.json"))[key]


def _freq(sc, s, n=1500, seeds=(1, 2, 3)):
    """A scenario's planted mutations per 10^3 bp of window, read-weighted, from its design alone."""
    p = _params(s)
    ref = read_fasta(s["reference"])[0][1]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)
    v = []
    for seed in seeds:
        t = lineage_sim.design(n, ref, lo, hi, sc, random.Random(seed))
        v.append(1e3 * sum(map(len, t["read_muts"].values())) / (n * (hi - lo + 1)))
    return statistics.mean(v)


def scenario(key, scen, deep=False):
    """The scenario, its loads scaled so the main scenario plants the cases' median real excess and "low" a
    quarter of it, or all of it where "low" is the only one (the off-targets). J꜀4 and N114 keep
    QUANT_ACCURACY's scenarios."""
    sc = SCEN[f"{scen}_deep"] if deep and scen != "null" else SCEN[scen]
    if key in PAPER or scen == "null":
        return sc
    ex = _calib(key)["case_excess_x1e3"]
    main = [x for x in SETS[key]["scenarios"] if x not in ("low", "null")]
    target = max(LOAD_FLOOR, statistics.median(ex)) / (4 if scen == "low" and main else 1)
    f = target / _freq(sc, SETS[key])
    return {**sc, **{x: sc[x] * f for x in ("trunk", "branch", "private", "dominant_load") if x in sc}}


def task(t):
    key, arm, variant, scen, depth, k, rep, peers_from = t
    s = SETS[key]
    names = sorted(s["controls"])
    ctl = names[k]
    if peers_from:
        peers = [str(cached(peers_from, n)) for n in sorted(SETS[peers_from]["controls"])]
    else:
        peers = [str(cached(key, n)) for n in names if n != ctl]
    seed = 7 + 31 * k + 1009 * rep + depth
    kw = {}
    if variant in ("qmap", "all"):
        kw["qmap"] = _qmap(key)
    if variant in PCR:
        kw["pcr"] = PCR[variant]
    if variant == "all":
        kw["pcr"] = PCR["pcr3"]
    sc = scenario(key, scen, deep=variant in ("deep", "all"))
    r = lineage_sim.run_one(str(cached(key, ctl)), peers, s["reference"], s["contaminants"], s["region"],
                            sc, seed, "recommended", depth=depth, project_extra=s["extra"], **kw)
    r.pop("roc", None)
    r.pop("lineages", None)
    r.update(set=key, arm=arm, variant=variant, scenario=scen, control=ctl, rep=rep, peers_from=peers_from)
    return r


def tasks(arm):
    out = []
    if arm == "breadth":
        for key in BREADTH:
            s = SETS[key]
            for scen in s["scenarios"]:
                for depth in s["depths"]:
                    for k in range(len(s["controls"])):
                        for rep in (0, 1):
                            out.append((key, arm, "base", scen, depth, k, rep, None))
    elif arm == "realism":
        for key, s in SETS.items():
            scen = [x for x in s["scenarios"] if x != "low"]
            for variant in ("base", "qmap", "pcr3", "pcr10", "deep", "all"):
                for sc in scen:
                    if sc == "null" and variant in ("qmap", "deep"):
                        continue        # nothing planted: these change nothing
                    for k in range(len(s["controls"])):
                        out.append((key, arm, variant, sc, s["main"], k, 0, None))
    elif arm == "crossctx":
        for key, other in (("su_lps", "su_ch12"), ("su_ch12", "su_lps")):
            for scen in ("polyclonal" if key == "su_lps" else "subclonal", "null"):
                for k in range(3):
                    for rep in (0, 1):
                        out.append((key, arm, "base", scen, SETS[key]["main"], k, rep, other))
    return out


def _done():
    f = OUT / "runs.jsonl"
    if not f.exists():
        return set()
    return {tuple(json.loads(x)["_task"]) for x in open(f)}


def cmd_run(a):
    OUT.mkdir(exist_ok=True)
    todo = [t for t in tasks(a.arm) if tuple(t) not in _done()]
    if a.limit:
        todo = todo[:a.limit]
    print(f"{len(todo)} runs", flush=True)
    with ProcessPoolExecutor(a.jobs) as ex, open(OUT / "runs.jsonl", "a") as fh:
        futs = {ex.submit(task, t): t for t in todo}
        for i, f in enumerate(as_completed(futs), 1):
            t = futs[f]
            try:
                r = f.result()
            except Exception as e:      # one failed library should not lose the rest
                print("FAILED", t, repr(e), flush=True)
                continue
            r["_task"] = list(t)
            fh.write(json.dumps(r, default=lambda o: sorted(o) if isinstance(o, set) else str(o)) + "\n")
            fh.flush()
            print(f"{i}/{len(todo)}", *t[:6], flush=True)


# ---------------------------------------------------------------- report

PAIRS = {"frequency": ("est_read_weighted_x1e3", "true_read_weighted_x1e3"),
         "unique": ("est_unique_density_x1e5", "true_unique_density_x1e5"),
         "combined": ("est_combined_density_x1e5", "true_unique_density_x1e5")}


def _ratios(rs, m):
    e, t = PAIRS[m]
    return [r[e] / r[t] for r in rs if r.get(e) is not None and r.get(t)]


def _cell(v):
    if not v:
        return "–"
    if len(v) == 1:
        return f"{v[0]:.2f}"
    return f"{statistics.mean(v):.2f} ({min(v):.2f}–{max(v):.2f})"


def _diff(rs):
    v = [r["est_read_weighted_x1e3"] - r["true_read_weighted_x1e3"] for r in rs
         if r.get("est_read_weighted_x1e3") is not None and r.get("true_read_weighted_x1e3") is not None]
    t = [r["true_read_weighted_x1e3"] for r in rs if r.get("true_read_weighted_x1e3") is not None]
    return f"{statistics.mean(t):.2f} / {statistics.mean(v):+.2f}" if v else "–"


def _fc(rs):
    v = [r["false_calls"] for r in rs]
    return f"{statistics.mean(v):.1f}" + (f" (max {max(v)})" if v else "")


def _recall(rs):
    v = [r["recall_2plus"] for r in rs if r.get("recall_2plus") is not None]
    return f"{statistics.mean(v):.2f}" if v else "–"


def _null_density(rs):
    v = [r.get("est_unique_density_x1e5") for r in rs if r.get("est_unique_density_x1e5") is not None]
    return f"{statistics.mean(v):.2f}" if v else "–"


LABEL = {"ldha": "Illumina J꜀4 (PRJNA973112)", "ch12": "Illumina CH12 VDJ (PRJNA918596)",
         "hmces": "Illumina Ramos VDJ (PRJNA820976)", "su_lps": "ONT Sμ, LPS blasts",
         "su_ch12": "ONT Sμ, CH12", "pim1": "ONT Pim1", "ot1": "ONT OT1 (CH12)", "ot2": "ONT OT2 (CH12)",
         "jh4": "ONT J꜀4 (this study)", "n114": "ONT N114 (this study)"}


def cmd_report(a):
    runs = [json.loads(x) for x in open(OUT / "runs.jsonl")]
    g = defaultdict(list)
    for r in runs:
        g[(r["arm"], r["set"], r["variant"], r["scenario"], r["depth"])].append(r)
    out = ["# Accuracy beyond J꜀4 and N114: report", "",
           "Estimate / truth: mean over held-out controls and draws, with the range in brackets "
           "(1.00 = exact). Null rows give the false calls and the unique density estimated where the truth "
           "is zero. Recall is of mutations carried by two or more analysed reads. Truth / error: the planted "
           "mutation frequency and the mean estimate minus it, per 10³ bp; at the loads of these studies a "
           "library holds a few dozen planted mutations, and the ratios swing with them.", ""]
    hdr = ["| Set | Scenario | Depth | Runs | Truth / error | Mutation frequency | Unique density | Combined density | "
           "False calls | Recall, 2+ reads |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for arm, title in (("breadth", "Breadth: every public set with three or more mutation-free controls"),
                       ("crossctx", "Sμ with controls from the other cell context")):
        keys = sorted((k for k in g if k[0] == arm), key=lambda k: (list(SETS).index(k[1]), k[3], k[4]))
        if not keys:
            continue
        out += [f"## {title}", ""] + hdr
        for k in keys:
            rs = g[k]
            row = [LABEL[k[1]], k[3], str(k[4]), str(len(rs))]
            if k[3] == "null":
                row += ["–", "–", f"{_null_density(rs)} (est)", "–"]
            else:
                row += [_diff(rs)] + [_cell(_ratios(rs, m)) for m in PAIRS]
            row += [_fc(rs), _recall(rs) if k[3] != "null" else "–"]
            out.append("| " + " | ".join(row) + " |")
        out.append("")
    keys = [k for k in g if k[0] == "realism"]
    if keys:
        out += ["## Realism: the planted side made less ideal, one step at a time", "",
                "| Set | Scenario | Planted | Runs | Truth / error | Mutation frequency | Unique density | "
                "Combined density | False calls | Recall, 2+ reads |", "|---|---|---|---:|---:|---:|---:|---:|---:|---:|"]
        order = ("base", "qmap", "pcr3", "pcr10", "deep", "all")
        for k in sorted(keys, key=lambda k: (list(SETS).index(k[1]), k[3], order.index(k[2]))):
            rs = g[k]
            row = [LABEL[k[1]], k[3], k[2], str(len(rs))]
            if k[3] == "null":
                row += ["–", "–", f"{_null_density(rs)} (est)", "–"]
            else:
                row += [_diff(rs)] + [_cell(_ratios(rs, m)) for m in PAIRS]
            row += [_fc(rs), _recall(rs) if k[3] != "null" else "–"]
            out.append("| " + " | ".join(row) + " |")
        out.append("")
    md = "\n".join(out)
    open(OUT / "report.md", "w").write(md)
    print(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("cache", "qmaps", "run", "report"))
    ap.add_argument("arm", nargs="?", choices=("breadth", "realism", "crossctx"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--limit", type=int, help="run only the first N outstanding tasks (a timing check)")
    a = ap.parse_args()
    {"cache": cmd_cache, "qmaps": cmd_qmaps, "run": cmd_run, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
