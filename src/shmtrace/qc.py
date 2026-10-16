"""Checks on Sanger clones that flag, never fix.

The errors that most change a Sanger result are not calling errors but
bookkeeping ones: a clone filed under the wrong animal, one trace counted
twice, a contaminating clone in a knockout control. None of these can be
corrected automatically with confidence, so each is reported, with the
clone's provenance, for a person to decide.

Checks, each a row in `flags` with a severity of "problem" (likely wrong),
"check" (worth a look) or "info":

  shared_mutations   a clone shares more somatic mutations with a clone of
                     another animal than independent animals share by chance.
                     Mutations are made independently in each animal, so a
                     clone sharing many with another animal's clone most likely
                     came from that animal (misfiled, or a copied trace).
  duplicate_file     two clones are the same file, byte for byte (a copied trace).
  duplicate_clone    the same, and the two clones also agree at every base they
                     both read: one trace counted twice, or two clones of one
                     lineage, which can only have come from one animal.
  label_mismatch     the clone's tube label names another sample (needs the
                     sample sheet's aliases column).
  control_error_site a call at a position where the control clones mismatch too
  control_mutations  a clone of a control (e.g. AID-knockout) animal carries
                     more mutations than its own base qualities predict
                     (calibrated on the other control clones; Poisson,
                     Bonferroni over control clones): contamination, or a
                     mislabelled clone.
  load_outlier       a clone carries far more mutations than the rest of its
                     animal's clones.
  short_clone        a clone reads less than a third of the analysis window.
  germline_or_lineage  differences the germline rule removed (>= 50% of an
                     animal's clones) that no other animal shows. A strain
                     difference shows in every animal of the strain; a lineage
                     filling half an animal's clones looks the same to the
                     rule, and its mutations are then lost. Per animal.

`excluded` lists every input clone that was not analysed, and why: by the
sample sheet, or by the pipeline (no quality segment, unmapped, contaminant or
paralogue, no clean stretch, minor haplotype, ...).

The chance-sharing test: each mutation (position, base) is weighted by the
number of animals it appears in, counting an animal once however many of its
clones carry it, so a lineage expanded in one animal does not count as a
hotspot. For clones A and B, the expected number of shared mutations is
sum over A's mutations of min(1, load_B * w(m)), w being the weights scaled to
sum to 1. The tail probability of the observed sharing under a Poisson with
that mean is Bonferroni-corrected over every pair of clones from different
animals. Simulation (tests/test_qc.py) checks that it fires on planted
misfilings and stays quiet on correctly filed clones.
"""
from __future__ import annotations

import re
from collections import Counter, defaultdict

from scipy.stats import poisson

from . import motifs

SEVERITY_ORDER = {"problem": 0, "check": 1, "info": 2}


def _clone_table(details: dict, lo: int, hi: int, minq: int) -> dict:
    """clone -> {sample, muts: set[(pos1, alt)], load, bases, calls: {pos0: base},
    exp_err: errors its base qualities predict among the bases it read}"""
    out = {}
    for sample, res in details.items():
        sites = {(v.pos, v.alt) for v in res["variants"]}
        for a in res["_aligned"]:
            calls, exp_err = {}, 0.0
            for r, b, q in a.pairs:
                if lo <= r <= hi and q >= minq and b in "ACGT":
                    calls[r] = b
                    exp_err += 10 ** (-q / 10)
            muts = {(r + 1, b) for r, b in calls.items() if (r + 1, b) in sites}
            out[a.qname] = {"sample": sample, "muts": muts, "load": len(muts),
                            "bases": len(calls), "calls": calls, "exp_err": exp_err}
    return out


def _weights(clones: dict, skip: tuple = ()) -> dict:
    """mutation -> share of (animal, mutation) occurrences, each animal once."""
    by_animal = defaultdict(set)
    for name, c in clones.items():
        if name in skip:
            continue
        by_animal[c["sample"]] |= c["muts"]
    n = Counter(m for ms in by_animal.values() for m in ms)
    tot = sum(n.values()) or 1
    return {m: k / tot for m, k in n.items()}


def shared_mutation_pairs(clones: dict, controls: set, alpha: float = 0.01,
                          min_shared: int = 3) -> list[dict]:
    names = sorted(clones)
    w_all = _weights(clones)
    cand = []
    for i, a in enumerate(names):
        A = clones[a]
        if A["load"] < min_shared:
            continue
        for b in names[i + 1:]:
            B = clones[b]
            if B["sample"] == A["sample"] or B["load"] < min_shared:
                continue
            shared = A["muts"] & B["muts"]
            if len(shared) >= min_shared:
                cand.append((a, b, shared))
    # every pair of clones from different animals that could share min_shared is a
    # test, whether or not it was a candidate
    per_sample = Counter(c["sample"] for c in clones.values() if c["load"] >= min_shared)
    n_tests = max(1, (sum(per_sample.values()) ** 2 - sum(v * v for v in per_sample.values())) // 2)
    rows = []
    for a, b, shared in cand:
        A, B = clones[a], clones[b]
        # weights without these two clones, so their own sharing is not the null
        w = _weights(clones, skip=(a, b)) if len(cand) < 5000 else w_all
        floor = min(w.values(), default=1e-3) / 2
        lam = sum(min(1.0, B["load"] * w.get(m, floor)) for m in A["muts"])
        lam = max(lam, 1e-9)
        p = float(poisson.sf(len(shared) - 1, lam))
        padj = min(1.0, p * n_tests)
        if padj >= alpha:
            continue
        common = [r for r in A["calls"] if r in B["calls"]]
        differ = sum(A["calls"][r] != B["calls"][r] for r in common)
        rows.append({"a": a, "b": b, "shared": len(shared), "expected": round(lam, 2),
                     "p_adj": padj, "overlap_bp": len(common), "differences": differ,
                     "identical": differ == 0 and len(common) >= 150})
    return rows


def label_tokens(label: str) -> set[str]:
    return {t for t in re.split(r"[^0-9a-z]+", (label or "").lower()) if t}


def run(details: dict, samples: list[dict], ref: str, lo: int, hi: int, minq: int = 20,
        provenance: dict | None = None, fates: dict | None = None,
        sheet_excluded: list | None = None, alpha: float = 0.01,
        aliases: dict | None = None, control_sites: dict | None = None) -> dict:
    provenance = provenance or {}
    controls = {s["name"] for s in samples if s.get("control")}
    clones = _clone_table(details, lo, hi, minq)
    flags: list[dict] = []

    def prov(name):
        p = provenance.get(name) or provenance.get(name.split("#")[0]) or {}
        return {k: p.get(k, "") for k in ("order", "well", "label", "file")}

    def flag(check, severity, clone, detail, other="", other_sample=""):
        flags.append({"check": check, "severity": severity, "sample": clones[clone]["sample"],
                      "clone": clone, "mutations": clones[clone]["load"],
                      "other_sample": other_sample, "other_clone": other, "detail": detail,
                      **prov(clone)})

    # 1-2. clones sharing mutations across animals
    pairs = shared_mutation_pairs(clones, controls, alpha)
    # identical clones in two animals: two animals cannot make the same mutated clone,
    # so identity is evidence enough without the sharing test
    seen_pairs = {(r["a"], r["b"]) for r in pairs}
    names = sorted(clones)
    for i, a in enumerate(names):
        A = clones[a]
        if A["load"] < 3:
            continue
        for b in names[i + 1:]:
            B = clones[b]
            if B["sample"] == A["sample"] or (a, b) in seen_pairs or len(A["muts"] & B["muts"]) < 3:
                continue
            common = [r for r in A["calls"] if r in B["calls"]]
            if len(common) >= 150 and all(A["calls"][r] == B["calls"][r] for r in common):
                pairs.append({"a": a, "b": b, "shared": len(A["muts"] & B["muts"]),
                              "expected": None, "p_adj": None, "overlap_bp": len(common),
                              "differences": 0, "identical": True})
    best_other = defaultdict(lambda: defaultdict(int))   # clone -> other sample -> pairs
    for r in pairs:
        A, B = clones[r["a"]], clones[r["b"]]
        best_other[r["a"]][B["sample"]] += 1
        best_other[r["b"]][A["sample"]] += 1
        if r["identical"]:
            for x, y in ((r["a"], r["b"]), (r["b"], r["a"])):
                flag("duplicate_clone", "problem", x,
                     f"identical to {y} ({clones[y]['sample']}) at all {r['overlap_bp']} bases both "
                     f"read, sharing {r['shared']} mutations; two animals cannot make the same "
                     "clone, so these come from one animal (the same trace twice, or one lineage)",
                     y, clones[y]["sample"])
    own_links = defaultdict(int)    # clone -> significant links within its own animal
    for a in clones:
        if clones[a]["load"] < 3:
            continue
        for b in clones:
            if b != a and clones[b]["sample"] == clones[a]["sample"] \
                    and len(clones[a]["muts"] & clones[b]["muts"]) >= 3:
                own_links[a] += 1
    for name, others in best_other.items():
        if any(f["clone"] == name and f["check"] == "duplicate_clone" for f in flags):
            continue
        top, n = max(others.items(), key=lambda kv: kv[1])
        mine = own_links.get(name, 0)
        sev = "problem" if n > mine else "check"
        flag("shared_mutations", sev, name,
             f"shares more mutations than chance with {n} clone(s) of {top}"
             + (f" and with {mine} of its own animal's" if mine else
                ", and with none of its own animal's") + "; it may belong to " + top,
             other_sample=top,
             other=", ".join(sorted(r["b"] if r["a"] == name else r["a"] for r in pairs
                                    if name in (r["a"], r["b"]) and
                                    clones[r["b"] if r["a"] == name else r["a"]]["sample"] == top))[:200])

    # 2a. the same file under two names
    import hashlib
    by_hash: dict = {}
    for name in sorted(clones):
        f = prov(name)["file"]
        if not f:
            continue
        try:
            h = hashlib.sha1(open(f, "rb").read()).hexdigest()
        except OSError:
            continue
        if h in by_hash:
            o = by_hash[h]
            for x, y in ((name, o), (o, name)):
                flag("duplicate_file", "problem", x,
                     f"the same file, byte for byte, as {y} ({clones[y]['sample']}): one trace "
                     "counted twice", y, clones[y]["sample"])
        else:
            by_hash[h] = name

    # 2b. the tube label names another sample
    if aliases:
        for name, c in clones.items():
            toks = label_tokens(prov(name)["label"])
            mine = set(aliases.get(c["sample"], []))
            named = sorted(s for s, al in aliases.items() if s != c["sample"] and toks & set(al))
            if named and not toks & mine:
                flag("label_mismatch", "check", name,
                     f"labelled '{prov(name)['label']}', which names {', '.join(named)}",
                     other_sample=named[0])

    # 3. mutations in a control animal's clones, beyond what the clone's own base
    # qualities predict (scaled by how far the controls' calls exceed their predicted
    # errors, leaving the clone out); Bonferroni over the control clones
    ctrl = [n for n, c in clones.items() if c["sample"] in controls]
    tot_calls = sum(clones[n]["load"] for n in ctrl)
    tot_exp = sum(clones[n]["exp_err"] for n in ctrl)
    for name in ctrl:
        c = clones[name]
        if c["load"] < 3:
            continue
        scale = max(1.0, (tot_calls - c["load"]) / max(tot_exp - c["exp_err"], 1e-9))
        lam = max(scale * c["exp_err"], 1e-9)
        if float(poisson.sf(c["load"] - 1, lam)) * len(ctrl) >= alpha:
            continue
        hot = sum(1 for pos, _alt in c["muts"]
                  if motifs.is_site("WRC", ref, pos - 1) or motifs.is_site("GYW", ref, pos - 1)
                  or motifs.is_site("WA", ref, pos - 1) or motifs.is_site("TW", ref, pos - 1))
        flag("control_mutations", "problem", name,
             f"{c['load']} mutations in a control animal's clone where its base qualities "
             f"predict {lam:.1f} errors ({hot} at AID/Pol-eta hotspots): contamination or a "
             "mislabelled clone")

    # 3b. a call where the control clones mismatch too: a position the chemistry gets wrong
    for pos_s, (n_ctl, depth_ctl) in sorted((control_sites or {}).items()):
        pos = int(pos_s)
        for name, c in sorted(clones.items()):
            if c["sample"] in controls:
                continue
            alts = [a for pp, a in c["muts"] if pp == pos]
            if not alts:
                continue
            flag("control_error_site", "info", name,
                 f"the call at {pos}{alts[0]} is at a position where {n_ctl} of {depth_ctl} control "
                 "clones also mismatch: the chemistry misreads it, so weigh this call accordingly "
                 "(nothing was removed)")

    # 4. a clone far more mutated than the rest of its animal
    by_sample = defaultdict(list)
    for c in clones.values():
        by_sample[c["sample"]].append(c["load"])
    for name, c in clones.items():
        loads = sorted(by_sample[c["sample"]])
        if c["sample"] in controls or len(loads) < 8:
            continue
        q1, q3 = loads[len(loads) // 4], loads[(3 * len(loads)) // 4]
        cut = q3 + 3 * max(q3 - q1, 1)
        if c["load"] >= max(6, cut):
            flag("load_outlier", "check", name,
                 f"{c['load']} mutations where this animal's clones have a middle half of "
                 f"{q1}-{q3}")

    # 5. short clones
    width = hi - lo + 1
    for name, c in clones.items():
        if c["bases"] < width / 3:
            flag("short_clone", "info", name, f"reads {c['bases']} of the window's {width} bases")

    # 6. differences removed as germline in only one animal
    removed = {sm: {(x[0], x[1]): x for x in (res["summary"]["calling"].get("germline_sites") or [])
                    if lo + 1 <= x[0] <= hi + 1}
               for sm, res in details.items()}
    for sm, sites in removed.items():
        others = set().union(*[set(v) for k, v in removed.items() if k != sm]) if len(removed) > 1 else set()
        if len(removed) < 2:
            continue
        own = sorted(k for k in sites if k not in others)
        if own:
            flags.append({"check": "germline_or_lineage", "severity": "check", "sample": sm,
                          "clone": "", "mutations": len(own), "other_sample": "",
                          "other_clone": "",
                          "detail": f"{len(own)} difference(s) removed as germline that no other "
                          "animal shows: " + ", ".join(f"{p_}{a_} ({sites[(p_, a_)][2]}/"
                                                      f"{sites[(p_, a_)][3]} clones)"
                                                      for p_, a_ in own[:10])
                          + ". A germline difference would usually show in the other animals of "
                          "the strain; an expanded lineage in half this animal's clones would "
                          "not. If they are a lineage's mutations they are missing from this "
                          "animal's counts.",
                          "order": "", "well": "", "label": "", "file": ""})

    excluded = []
    for e in sheet_excluded or []:
        excluded.append({"clone": e["clone_id"], "sample": e["sample"], "by": "sample sheet",
                         "reason": e["reason"], **{k: e.get(k, "") for k in
                                                   ("order", "well", "label", "file")}})
    for sample, fs in (fates or {}).items():
        for q, rec in sorted(fs.items()):
            excluded.append({"clone": q, "sample": sample, "by": "pipeline",
                             "reason": rec.get("fate", ""), **prov(q)})

    flags.sort(key=lambda f: (SEVERITY_ORDER[f["severity"]], f["check"], f["sample"], f["clone"]))
    # what ties each pair together, as [1-based position, base], so a person (or a check against
    # another library of the same mice) can see the lineage rather than only its size
    for r in pairs:
        if r["a"] in clones and r["b"] in clones:
            r["shared_mutations"] = sorted([p, b] for p, b in clones[r["a"]]["muts"] & clones[r["b"]]["muts"])
    summary = Counter((f["check"], f["severity"]) for f in flags)
    return {"flags": flags, "excluded": excluded, "pairs": pairs,
            "duplicate_groups": duplicate_groups(flags, clones),
            "summary": [{"check": c, "severity": s, "clones": n} for (c, s), n in sorted(summary.items())],
            "clones_analysed": len(clones)}


def duplicate_groups(flags: list[dict], clones: dict) -> list[dict]:
    """One row per set of clones that are the same sequence, instead of a row per pair.

    `duplicate_clone` and `duplicate_file` are flagged pairwise, so a dozen copies of one
    sequence become dozens of flags (134 of them in one real project). The pairs are the
    edges of a graph whose connected components are the sequences themselves; this returns
    those components, each with the samples it spans, which is what a person has to decide
    about.
    """
    parent: dict = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    kinds: dict = defaultdict(set)
    for f in flags:
        if f["check"] not in ("duplicate_clone", "duplicate_file") or not f.get("other_clone"):
            continue
        a, b = f["clone"], f["other_clone"]
        kinds[find(a)].add(f["check"])
        kinds[find(b)].add(f["check"])
        parent[find(a)] = find(b)
    members: dict = defaultdict(list)
    for c in parent:
        members[find(c)].append(c)
    out = []
    for cs in members.values():
        cs = sorted(cs)
        samples = sorted({clones[c]["sample"] for c in cs if c in clones})
        kind = set()
        for c in cs:
            kind |= kinds.get(find(c), set())
        sets = [clones[c]["muts"] for c in cs if c in clones]
        out.append({"clones": cs, "n_clones": len(cs), "samples": samples, "n_samples": len(samples),
                    "checks": sorted(kind),
                    "shared_mutations": sorted([p, b] for p, b in set.intersection(*sets)) if sets else [],
                    "mutations": max((clones[c]["load"] for c in cs if c in clones), default=0),
                    "severity": "problem" if len(samples) > 1 else "check"})
    out.sort(key=lambda g: (-g["n_samples"], -g["n_clones"], g["clones"][0]))
    return out
