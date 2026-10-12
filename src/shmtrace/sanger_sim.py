"""Validate the Sanger pipeline on simulated clones whose truth is known.

A real Sanger set has no answer key, so this module writes one: a project of
simulated animals, each a set of B-cell lineages, sequenced as .ab1 traces
with realistic quality (noisy starts, fading ends, rare confident errors),
plus the traces that go wrong in practice (failed reads, paralogue clones,
short reads). It then plants the bookkeeping errors this pipeline's QC is
meant to catch, writes a sample sheet as a person would, runs the pipeline
exactly as on real data, and scores it against the truth:

  calling   per animal, the distinct mutations called against those truly
            present (inside the window, at bases a clone read at Q >= 20):
            recall, precision (in animals with 20 or more true mutations), the
            false calls per 10^5 bases read (the single-pass Sanger error floor
            at Q20, also what the controls measure), and the animal's planted
            germline difference, which must not be called.
  exclusion failed, paralogue and other unusable traces must not be analysed.
  QC        every planted error must raise its flag; clones with no planted
            error should raise none ("false flags").

Planted errors (each in a known clone):
  misfiled_labelled    a mutated clone filed under another animal, its tube
                       label still naming its own (label_mismatch, and
                       shared_mutations with its lineage mates)
  misfiled_relabelled  the same, but the label was also changed, so only the
                       shared mutations can reveal it (shared_mutations)
  retyped_copy         a clone's sequence entered again under another animal
                       as a FASTA record, like a re-typed GenBank entry
                       (duplicate_clone)
  file_copy            a trace file copied under another animal (the sample
                       sheet's byte-for-byte warning, and duplicate_clone)
  control_contaminant  a mutated clone filed under a control (control_mutations)
  lineage_misfiled     a whole lineage filed under another animal, labels
                       intact (label_mismatch; sharing cannot see it, since no
                       mate remains in the true animal)

    shmtrace-validate-sanger --reps 10 -o sanger_validation/
"""
from __future__ import annotations

import csv
import json
import math
import os
import random
import shutil
from collections import defaultdict
from pathlib import Path

from .io_utils import read_fasta, write_abif
from .simulate import _alt_for, _weight

HERE = Path(__file__).parent
REFERENCE = str(HERE / "references/mouse_jh4_intron.fasta")
CONTAMINANT = str(HERE / "references/mouse_jh4_contaminant.fasta")
EXPECTED = {"misfiled_labelled": {"label_mismatch", "shared_mutations"},
            "misfiled_relabelled": {"shared_mutations"},
            "retyped_copy": {"duplicate_clone"},
            "file_copy": {"duplicate_file", "duplicate_clone"},
            "control_contaminant": {"control_mutations"},
            "lineage_misfiled": {"label_mismatch"}}


def _poisson(lam: float, rng: random.Random) -> int:
    k, t, L = 0, rng.random(), math.exp(-lam)
    while t > L:
        k += 1
        t *= rng.random()
    return k


def _read(seq: str, rng: random.Random, start=None, length=None) -> tuple[str, list[int]]:
    """A Sanger read of `seq`: low quality for the first bases, a clean middle with
    rare confident errors, and a fading end. Errors follow the qualities."""
    s = rng.randint(0, 20) if start is None else start
    n = length or rng.randint(620, 760)
    frag = seq[s:s + n]
    head, tail = rng.randint(12, 30), rng.randint(60, 160)
    out, quals = [], []
    for i, b in enumerate(frag):
        if i < head:
            q = int(5 + 30 * i / head)
        elif i >= len(frag) - tail:
            q = int(max(3, 38 * (len(frag) - i) / tail))
        else:
            q = rng.randint(32, 55)
        err = 10 ** (-q / 10)
        if rng.random() < err:
            b = rng.choice([x for x in "ACGT" if x != b]) if q > 8 or rng.random() < .7 else "N"
        out.append(b)
        quals.append(q)
    return "".join(out), quals


def _apply(ref: str, muts: dict) -> str:
    s = list(ref)
    for p0, alt in muts.items():
        s[p0] = alt
    return "".join(s)


def simulate(out: str, seed: int = 1, n_wt: int = 3, n_ko: int = 3, clones: int = 30,
             plant: bool = True) -> dict:
    """Write traces, a sample sheet and the truth for one simulated project."""
    rng = random.Random(seed)
    out = Path(out)
    if out.exists():
        shutil.rmtree(out)
    (out / "traces").mkdir(parents=True)
    ref = read_fasta(REFERENCE)[0][1].upper()
    para = read_fasta(CONTAMINANT)[1][1].upper()
    lo, hi = 100, 628                                  # 0-based JH4 window (101-629)
    # AID targets the JH4 intron next to the J segment: mutations fall in and
    # around the analysed window, thinning beyond it
    pool = [i for i in range(lo - 60, hi + 120)]
    wts = [_weight(ref, i) for i in pool]

    def draw(n: int, avoid=()) -> dict:
        got: dict = {}
        while len(got) < n:
            i = rng.choices(pool, wts)[0]
            if i not in got and i not in avoid:
                got[i] = _alt_for(ref[i], rng)
        return got

    animals = [f"WT{i + 1}" for i in range(n_wt)] + [f"KO{i + 1}" for i in range(n_ko)]
    clone_rows: dict = {}
    lineages: dict = defaultdict(list)                 # (animal, lineage) -> clone ids
    germline: dict = {}
    k = 0
    strain = draw(2)                                   # differences every animal carries
    private: dict = {}
    for ai, a in enumerate(animals):
        wt = a.startswith("WT")
        g = dict(strain)
        if rng.random() < 0.2:                         # now and then, one of its own
            private[a] = draw(1, avoid=set(strain))
            g.update(private[a])
        germline[a] = [(p + 1, alt) for p, alt in g.items()]
        n_lin = rng.randint(4, 8) if wt else 1
        sizes = [rng.random() + .2 for _ in range(n_lin)]
        load = (3, 7, 11)[ai % 3] if wt else 0
        founders = [draw(_poisson(load, rng)) if wt and j else {} for j in range(n_lin)]
        for c in range(clones):
            j = rng.choices(range(n_lin), sizes)[0]
            muts = dict(founders[j])
            if wt and j:
                muts.update(draw(_poisson(2, rng), avoid=set(muts)))
            k += 1
            cid = f"P{seed}-{k:04d}"
            clone_rows[cid] = {"true": a, "filed": a, "lineage": j, "kind": "clone",
                               "muts": {p: alt for p, alt in muts.items()},
                               "label": f"{a.lower()}-{c + 1}", "planted": ""}
            lineages[(a, j)].append(cid)
        for kind in ("failed", "paralogue", "short"):
            for _ in range(rng.randint(1, 2)):
                k += 1
                cid = f"P{seed}-{k:04d}"
                clone_rows[cid] = {"true": a, "filed": a, "lineage": -1, "kind": kind, "muts": {},
                                   "label": f"{a.lower()}-x{k}", "planted": ""}

    extra_files: list = []
    if plant:
        wts_ = [a for a in animals if a.startswith("WT")]
        kos = [a for a in animals if a.startswith("KO")]

        def in_window(c):
            return sum(lo <= p0 <= hi for p0 in clone_rows[c]["muts"])

        def mutated_in(a, min_load=5, min_mates=3):
            """Clones of animal a with min_load mutations in the window and at least
            min_mates clones in their lineage (so mates remain after a move)."""
            out_ = []
            for (an, j), ids in lineages.items():
                if an == a and j and len(ids) >= min_mates:
                    out_ += [c for c in ids if in_window(c) >= min_load
                             and not clone_rows[c]["planted"]]
            return out_ or [c for c in clone_rows if clone_rows[c]["true"] == a
                            and clone_rows[c]["kind"] == "clone" and not clone_rows[c]["planted"]
                            and in_window(c) >= 3]
        x, y, z = wts_[0], wts_[1], wts_[-1]
        c1 = rng.choice(mutated_in(z))
        clone_rows[c1].update(filed=x, planted="misfiled_labelled")
        c2 = rng.choice([c for c in mutated_in(y) if c != c1])
        clone_rows[c2].update(filed=x, planted="misfiled_relabelled",
                              label=f"{x.lower()}-{clones + 1}")
        c3 = rng.choice(mutated_in(z))
        clone_rows[c3]["planted"] = "retyped_copy_source"
        extra_files.append(("retyped_copy", c3, y))
        c4 = rng.choice(mutated_in(y))
        clone_rows[c4]["planted"] = "file_copy_source"
        extra_files.append(("file_copy", c4, z))
        c5 = rng.choice(mutated_in(z))
        clone_rows[c5].update(filed=kos[0], planted="control_contaminant")
        # a whole small lineage of x, filed under z
        cand = [(a_, j) for (a_, j), ids in lineages.items() if a_ == x and j and 2 <= len(ids) <= 5
                and all(not clone_rows[c]["planted"] for c in ids)
                and any(len(clone_rows[c]["muts"]) >= 3 for c in ids)]
        if cand:
            for c in lineages[rng.choice(cand)]:
                clone_rows[c].update(filed=z, planted="lineage_misfiled")

    # sequence everything
    alias = {a: a.lower() for a in animals}
    sheet = []
    traces: dict = {}
    for cid, c in clone_rows.items():
        a = c["true"]
        seq = _apply(ref, {**{p - 1: alt for p, alt in germline[a]}, **c["muts"]})
        if c["kind"] == "failed":
            s = "".join(rng.choice("ACGTN") for _ in range(rng.randint(150, 400)))
            q = [rng.randint(2, 9) for _ in s]
        elif c["kind"] == "paralogue":
            s, q = _read(para, rng)
        elif c["kind"] == "short":
            s, q = _read(seq, rng, start=rng.randint(150, 350), length=rng.randint(150, 200))
        else:
            s, q = _read(seq, rng)
        f = out / "traces" / f"{cid}_{c['label']}.ab1"
        write_abif(str(f), s, q, cid)
        traces[cid] = (f, s, q)
        sheet.append({"file": f"traces/{f.name}", "sample": c["filed"],
                      "group": "WT" if c["filed"].startswith("WT") else "KO",
                      "control": "no" if c["filed"].startswith("WT") else "yes",
                      "clone_id": cid, "label": c["label"], "exclude": "",
                      "aliases": alias[c["filed"]]})
    for kind, src, to in extra_files:
        cid = f"{src}-{kind}"
        f, s, q = traces[src]
        if kind == "file_copy":
            dest = out / "traces" / f"{cid}.ab1"
            shutil.copyfile(f, dest)
        else:                                          # the clean part, re-typed as FASTA
            good = [i for i, x in enumerate(q) if x >= 25]
            dest = out / "traces" / f"{cid}.fasta"
            dest.write_text(f">{cid}\n{s[good[0]:good[-1] + 1]}\n")
        clone_rows[cid] = {**clone_rows[src], "filed": to, "planted": kind, "copy_of": src,
                           "label": f"{to.lower()}-{len(sheet)}"}
        sheet.append({"file": f"traces/{dest.name}", "sample": to, "group": "WT",
                      "control": "no", "clone_id": cid, "label": clone_rows[cid]["label"],
                      "exclude": "", "aliases": alias[to]})
    with open(out / "samplesheet.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sheet[0]))
        w.writeheader()
        w.writerows(sheet)
    truth = {"seed": seed, "window": [lo + 1, hi + 1], "germline": germline,
             "private_germline": {a: [[p + 1, alt] for p, alt in g.items()] for a, g in private.items()},
             "clones": {cid: {**c, "muts": [[p + 1, alt] for p, alt in sorted(c["muts"].items())]}
                        for cid, c in clone_rows.items()}}
    json.dump(truth, open(out / "truth.json", "w"), indent=1)
    return truth


def score(res: dict, truth: dict) -> dict:
    lo, hi = truth["window"]
    clones = truth["clones"]
    details = res["_details"]
    analysed = {a.qname for d in details.values() for a in d["_aligned"]}
    # calling, per filed animal
    calling = []
    for animal, d in details.items():
        called = {(v.pos, v.alt) for v in d["variants"] if lo <= v.pos <= hi}
        true_set = set()
        for a in d["_aligned"]:
            c = clones.get(a.qname)
            if c is None:
                continue
            read = {r + 1 for r, b, q in a.pairs if q >= 20 and b in "ACGT"}
            true_set |= {(p, alt) for p, alt in c["muts"] if lo <= p <= hi and p in read}
        germ = {tuple(x) for x in truth["germline"].get(animal, [])}
        tp = len(called & true_set)
        bases = sum(r["bases"] for r in d["per_read"])
        calling.append({"animal": animal, "true": len(true_set), "called": len(called), "tp": tp,
                        "recall": tp / len(true_set) if true_set else None,
                        "precision": tp / len(called) if called else None,
                        "germline_called": len(called & germ), "bases": bases,
                        "false_per_1e5": round(1e5 * (len(called) - tp) / bases, 1) if bases else None})
    # exclusion
    bad = [cid for cid, c in clones.items() if c["kind"] in ("failed", "paralogue")]
    wrongly_kept = [cid for cid in bad if cid in analysed]
    # QC
    flags = defaultdict(set)
    for f in (res.get("qc") or {}).get("flags", []):
        if f["severity"] != "info":
            flags[f["clone"]].add(f["check"])
    planted = {cid: c["planted"] for cid, c in clones.items() if c["planted"] in EXPECTED}
    detected = {cid: bool(flags.get(cid, set()) & EXPECTED[kind]) for cid, kind in planted.items()}
    involved = set(planted) | {c.get("copy_of") for c in clones.values() if c.get("copy_of")}
    # a clone whose lineage mate was misfiled shares mutations with it, legitimately flagged
    mates = defaultdict(set)
    for cid, c in clones.items():
        if c["kind"] == "clone" and c["lineage"] > 0:
            mates[(c["true"], c["lineage"])].add(cid)
    for cid in list(planted):
        c = clones[cid]
        involved |= mates.get((c["true"], c["lineage"]), set())
    # true mutations in half or more of an animal's clones (an expanded lineage, or a
    # hotspot hit in several lineages) are removed by the germline rule; each such
    # animal must be flagged
    lineage_germline = []
    for animal, d in details.items():
        germ = {tuple(x) for x in truth["germline"].get(animal, [])}
        lost = [(x[0], x[1]) for x in d["summary"]["calling"].get("germline_sites") or []
                if lo <= x[0] <= hi and (x[0], x[1]) not in germ]
        if lost:
            lineage_germline.append({"animal": animal, "lost": lost, "flagged": any(
                f["check"] == "germline_or_lineage" and f["sample"] == animal
                for f in (res.get("qc") or {}).get("flags", []))})
    false_flags = sorted((cid, sorted(v - {"load_outlier"})) for cid, v in flags.items()
                         if cid not in involved and cid in clones and v - {"load_outlier"})
    def shared_with_mates(cid):
        """Most window mutations the clone shares with a mate still in its true animal."""
        c = clones[cid]
        mine = {tuple(m) for m in c["muts"] if lo <= m[0] <= hi}
        return max([len(mine & {tuple(m) for m in clones[o]["muts"]}) for o in
                    mates.get((c["true"], c["lineage"]), set())
                    if o != cid and clones[o]["filed"] == c["true"]] or [0])

    return {"calling": calling, "excluded_ok": not wrongly_kept, "wrongly_kept": wrongly_kept,
            "planted": [{"clone": cid, "error": planted[cid], "detected": detected[cid],
                         "flags": sorted(flags.get(cid, set())),
                         "shared_with_mates": shared_with_mates(cid)} for cid in planted],
            "false_flags": false_flags, "lineage_germline": lineage_germline,
            "germline_flags_without_lineage": [f["sample"] for f in (res.get("qc") or {}).get("flags", [])
                                               if f["check"] == "germline_or_lineage"
                                               and f["sample"] not in {x["animal"] for x in lineage_germline}
                                               and f["sample"] not in truth.get("private_germline", {})],
            "clean_clones": sum(1 for cid in clones if cid not in involved and cid in analysed)}


def run_once(out: str, seed: int = 1, plant: bool = True) -> dict:
    from . import project, samplesheet
    truth = simulate(out, seed=seed, plant=plant)
    sh = samplesheet.load(os.path.join(out, "samplesheet.csv"))
    proj = sh.project(REFERENCE, "jh4_intron", contaminants=[CONTAMINANT])
    res = project.run_project(proj)
    s = score(res, truth)
    s["sheet_warnings"] = sh.warnings
    return s


def summarise(runs: list[dict]) -> dict:
    rec = [c["recall"] for r in runs for c in r["calling"] if c["recall"] is not None]
    pre = [c["precision"] for r in runs for c in r["calling"]
           if c["precision"] is not None and c["true"] >= 20]
    fp = sum(c["called"] - c["tp"] for r in runs for c in r["calling"])
    bp = sum(c["bases"] for r in runs for c in r["calling"])
    by_err = defaultdict(lambda: [0, 0])
    for r in runs:
        for p in r["planted"]:
            by_err[p["error"]][0] += p["detected"]
            by_err[p["error"]][1] += 1
    clean = sum(r["clean_clones"] for r in runs)
    # misfiled clones whose label was changed too: caught by shared mutations alone,
    # by how many they share with their lineage mates
    by_shared = defaultdict(lambda: [0, 0])
    for r in runs:
        for p_ in r["planted"]:
            if p_["error"] == "misfiled_relabelled":
                k = min(p_["shared_with_mates"], 8)
                by_shared[k][0] += p_["detected"]
                by_shared[k][1] += 1
    ff = sum(len(r["false_flags"]) for r in runs)
    return {"runs": len(runs),
            "recall_min": min(rec), "recall_mean": sum(rec) / len(rec),
            "precision_min_mutated": min(pre), "precision_mean_mutated": sum(pre) / len(pre),
            "false_calls_per_1e5_bases": round(1e5 * fp / bp, 1),
            "germline_called": sum(c["germline_called"] for r in runs for c in r["calling"]),
            "unusable_traces_analysed": sum(len(r["wrongly_kept"]) for r in runs),
            "detection": {k: f"{d}/{n}" for k, (d, n) in sorted(by_err.items())},
            "false_flags": f"{ff} of {clean} clean clones",
            "relabelled_misfiling_by_shared_mutations": {
                (f"{k}+" if k == 8 else str(k)): f"{d}/{n}" for k, (d, n) in sorted(by_shared.items())},
            "lost_to_germline_rule_flagged": f"{sum(x['flagged'] for r in runs for x in r['lineage_germline'])}"
                                        f"/{sum(len(r['lineage_germline']) for r in runs)}",
            "germline_flags_without_cause": sum(len(r["germline_flags_without_lineage"])
                                                           for r in runs)}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog="shmtrace-validate-sanger", description=__doc__.split("\n\n")[0])
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    runs = []
    for i in range(a.reps):
        r = run_once(os.path.join(a.out, f"rep{i + 1}"), seed=a.seed + i)
        runs.append(r)
        print(f"rep {i + 1}: " + json.dumps(summarise([r])))
    s = summarise(runs)
    json.dump({"summary": s, "runs": runs}, open(os.path.join(a.out, "validation.json"), "w"),
              indent=1, default=list)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
