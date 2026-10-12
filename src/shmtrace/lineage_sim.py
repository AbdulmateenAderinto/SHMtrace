"""End-to-end simulation with lineage truth: clones, trees and expansion.

`simulate.py` plants mutations at chosen read counts and scores the caller.
What it cannot test is anything that depends on clone structure: whether the
junction finds the clones, whether one expanded clone biases the read-weighted
frequency, and whether the clone-weighted number recovers the truth. Real
libraries have no known clone structure, so this module makes one.

A library is built in the spirit of SHMulate (Yaari et al. 2013), on real
reads:

1. Reads of a mutation-free control are dealt into clones. Clone sizes fall
   off like a power law; optionally one clone takes a fixed share of the
   library (the expanded, or PCR-jackpotted, clone Wt2 turned out to carry).
2. Each clone gets a random 20-bp junction, written over the bases 5' of J in
   every read that has one, with substitutions at the read's own error rate
   and occasional single-base indels, so clustering faces nanopore-like noise.
3. Each clone gets a lineage: trunk mutations in all its reads, a branch in a
   random part of them, a sub-branch nested in that, and private mutations in
   single reads. Positions are drawn by a 5-mer targeting model when one is
   given (otherwise by WRC/GYW and WA/TW weights), substitutions by the model's
   substitution probabilities (otherwise G:C transitions favoured).
4. The library runs through `project.run_project` like any other, with the
   other controls as its null, and is scored against the truth.

The background errors are the control's own, so the false-call side is real;
the planted side is idealised in the way `simulate.py` describes. Hold the
planted control out of the null.

Scores:

* **Calls.** Recall by how many analysed reads carry the mutation, false calls
  (called but never planted), and an ROC: recall against false calls as the
  FDR threshold sweeps, from every candidate the test saw. The sweep is the
  test's own decision; post-test filters (strand, homopolymer) are applied
  only at the operating point.
* **Clones.** Reads with a junction, clones found against clones present, and
  the adjusted Rand index between true and recovered read-to-clone assignment.
* **Frequency.** The read-weighted and clone-weighted estimates against their
  true values, so the bias an expanded clone introduces can be measured, not
  just argued.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import tempfile
from collections import Counter, defaultdict

from . import clones, pipeline, project, simulate, targeting
from .align import build_index
from .io_utils import read_fasta, read_fastq

BASES = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")

# A scenario is a library design. Loads are mean mutations per clone (trunk),
# per branch, and per read (private); clone loads vary as a gamma of shape 2.
SCENARIOS = {
    # many clones, none dominant: read- and clone-weighting should agree
    "polyclonal": {"clones": 40, "zipf": 1.0, "dominant": 0.0,
                   "trunk": 6.0, "branch": 2.0, "private": 0.3},
    # Wt2's situation: one unmutated clone holds 60% of the reads
    "expanded": {"clones": 40, "zipf": 1.0, "dominant": 0.6, "dominant_load": 0.0,
                 "trunk": 6.0, "branch": 2.0, "private": 0.3},
    # lightly mutated, as in a weak or early response
    "low": {"clones": 40, "zipf": 1.0, "dominant": 0.0,
            "trunk": 1.0, "branch": 0.5, "private": 0.1},
    # a cell line expressing AID (the N114 hybridoma subclones): one parental
    # sequence, small subclones sharing a mutation or two, most mutations in one read
    "subclonal": {"clones": 60, "zipf": 0.6, "dominant": 0.0,
                  "trunk": 0.6, "branch": 0.4, "private": 0.15},
    # nothing planted except clone structure: every call is false
    "null": {"clones": 40, "zipf": 1.0, "dominant": 0.0,
             "trunk": 0.0, "branch": 0.0, "private": 0.0},
}


# ---------------------------------------------------------------- the library

def clone_sizes(n_reads: int, n_clones: int, zipf: float, dominant: float,
                rng: random.Random) -> list[int]:
    """Reads per clone: an optional dominant clone first, the rest power-law."""
    sizes = []
    rest = n_reads
    if dominant > 0:
        sizes.append(max(1, round(dominant * n_reads)))
        rest -= sizes[0]
        n_clones -= 1
    n_clones = max(1, min(n_clones, rest))
    w = [1.0 / (i + 1) ** zipf for i in range(n_clones)]
    tot = sum(w)
    tail = [1] * n_clones
    left = rest - n_clones
    # largest-remainder apportionment of what is left, by weight
    shares = [left * x / tot for x in w]
    for i, s in enumerate(shares):
        tail[i] += int(s)
    order = sorted(range(n_clones), key=lambda i: -(shares[i] - int(shares[i])))
    for i in order[:rest - sum(tail)]:
        tail[i] += 1
    rng.shuffle(tail)
    return sizes + sorted(tail, reverse=True)


def random_junctions(n: int, length: int, rng: random.Random, min_dist: int = 8) -> list[str]:
    """n random junctions, each at least `min_dist` edits from every other."""
    out: list[str] = []
    while len(out) < n:
        j = "".join(rng.choice(BASES) for _ in range(length))
        if all(clones._edit(j, o, min_dist - 1) >= min_dist for o in out):
            out.append(j)
    return out


def _site_sampler(ref: str, lo: int, hi: int, model):
    """Positions (0-based) and weights for AID-like placement."""
    if model is not None:
        prof = targeting.window_profile(model, ref, lo, hi)["mutability"]
        pos = sorted(prof)
        return pos, [prof[i] for i in pos]
    pos = [i for i in range(lo, hi + 1) if ref[i] in BASES]
    return pos, [simulate._weight(ref, i) for i in pos]


def _alt(ref: str, i: int, model, rng: random.Random) -> str:
    r = ref[i]
    if model is not None:
        probs = model.substitution.get(ref[i - 2:i + 3]) or {}
        choices = [(b, probs.get(b, 0.0)) for b in BASES if b != r]
        if sum(w for _b, w in choices) > 0:
            return rng.choices([b for b, _w in choices], weights=[w for _b, w in choices])[0]
    return simulate._alt_for(r, rng)


def _poisson(lam: float, rng: random.Random) -> int:
    if lam <= 0:
        return 0
    L, k, prod = math.exp(-lam), 0, rng.random()
    while prod > L:
        k += 1
        prod *= rng.random()
    return k


def design(n_reads: int, ref: str, lo: int, hi: int, scenario: dict, rng: random.Random,
           model=None, junction_length: int = 20, j_start0: int = 52) -> dict:
    """The truth: each read's clone, and each read's mutations (pos0, alt)."""
    sizes = clone_sizes(n_reads, scenario["clones"], scenario["zipf"],
                        scenario.get("dominant", 0.0), rng)
    juncs = random_junctions(len(sizes), junction_length, rng)
    # where each clone's junction meets J: J segments are trimmed by a few
    # bases at recombination (in the J꜀4 controls, alignments begin -3 to +16
    # bp from J), and the trimming is the clone's
    jstarts = [j_start0 + rng.randint(-3, 12) for _ in sizes]
    pos, w = _site_sampler(ref, lo, hi, model)

    def draw(k: int, taken: set) -> list[tuple[int, str]]:
        out = []
        while len(out) < k and len(taken) < len(pos):
            i = rng.choices(pos, weights=w)[0]
            if i not in taken:
                taken.add(i)
                out.append((i, _alt(ref, i, model, rng)))
        return out

    order = list(range(n_reads))
    rng.shuffle(order)
    read_clone, read_muts, lineages = {}, defaultdict(set), []
    start = 0
    for c, size in enumerate(sizes):
        members = order[start:start + size]
        start += size
        dominant = c == 0 and scenario.get("dominant", 0.0) > 0
        load = scenario.get("dominant_load", scenario["trunk"]) if dominant else scenario["trunk"]
        # clone loads vary: gamma, shape 2, with this mean
        lam = rng.gammavariate(2.0, load / 2.0) if load > 0 else 0.0
        taken: set = set()
        trunk = draw(_poisson(lam, rng), taken)
        branch, sub, bset, sset = [], [], [], []
        if len(members) >= 4 and not dominant:
            bset = rng.sample(members, max(1, int(len(members) * rng.uniform(0.3, 0.6))))
            branch = draw(_poisson(scenario["branch"], rng), taken)
            if len(bset) >= 2:
                sset = rng.sample(bset, max(1, int(len(bset) * rng.uniform(0.3, 0.6))))
                sub = draw(_poisson(scenario["branch"], rng), taken)
        for r in members:
            read_clone[r] = c
            read_muts[r].update(trunk)
            if r in bset:
                read_muts[r].update(branch)
            if r in sset:
                read_muts[r].update(sub)
            priv = draw(_poisson(scenario["private"], rng), set(taken))
            read_muts[r].update(priv)
        lineages.append({"clone": c, "reads": size, "junction": juncs[c], "j_start": jstarts[c] + 1,
                         "trunk": len(trunk), "branch": len(branch), "subbranch": len(sub)})
    return {"read_clone": read_clone, "read_muts": dict(read_muts),
            "junctions": juncs, "jstarts": jstarts, "lineages": lineages}


def _noisy(junction: str, quals: str, rng: random.Random, indel_rate: float) -> tuple[str, str]:
    """The junction as a nanopore read would show it: errors at the read's own rate."""
    out_s, out_q = [], []
    for k, b in enumerate(junction):
        q = quals[k] if k < len(quals) else "+"
        e = 10 ** (-(ord(q) - 33) / 10)
        u = rng.random()
        if u < indel_rate / 2:
            continue                                   # deleted
        if u < indel_rate:
            out_s.append(rng.choice(BASES))            # inserted before it
            out_q.append(q)
        out_s.append(rng.choice([x for x in BASES if x != b]) if rng.random() < e else b)
        out_q.append(q)
    return "".join(out_s), "".join(out_q)


def build_library(reads: list[tuple[str, str, str]], reference: str, truth: dict, p,
                  rng: random.Random, indel_rate: float = 0.04) -> tuple[list, dict]:
    """Write the designed clones into the reads.

    Mutations are written at reference positions the read aligns (not deleted);
    the junction replaces the bases 5' of the alignment start when that start
    is at J. Returns the new reads and, per read name, what it actually carries.
    """
    ref = read_fasta(reference)[0][1]
    aligner, _t, _c = build_index([reference], [], p.minimap2_preset)
    j0 = (p.junction_ref_start or 0) - 1
    L = p.junction_length
    out, carried = [], {}
    for idx, (name, seq, qual) in enumerate(reads):
        muts = truth["read_muts"].get(idx, set())
        c = truth["read_clone"][idx]
        hit = next((h for h in aligner.map(seq) if h.is_primary), None)
        if hit is None:
            out.append((name, seq, qual))
            carried[name] = {"clone": c, "muts": set(), "junction": False}
            continue
        r2q = simulate._ref_to_query(hit, len(seq))
        fs = list(seq if hit.strand == 1 else seq.translate(_COMP)[::-1])
        fq = qual if hit.strand == 1 else qual[::-1]
        got = set()
        for pos0, alt in muts:
            if pos0 in r2q and ref[pos0] != alt:
                fs[r2q[pos0]] = alt
                got.add((pos0, alt))
        fs = "".join(fs)
        has_j = False
        q0 = r2q.get(hit.r_st)
        t = truth["jstarts"][c]
        if (p.junction_ref_start and q0 is not None and q0 >= L
                and abs(hit.r_st - j0) <= clones.SLACK):
            # rebuild the read so the clone's junction abuts J at the clone's
            # own boundary t: every read of a clone shares its J trimming, as
            # in a real clone, whatever the donor read's own geometry was
            js, jq = _noisy(truth["junctions"][c], fq[q0 - L:q0], rng, indel_rate)
            if t >= hit.r_st:
                qt = next((r2q[x] for x in range(t, t + 10) if x in r2q), None)
                if qt is not None:
                    fs = fs[:q0 - L] + js + fs[qt:]
                    fq = fq[:q0 - L] + jq + fq[qt:]
                    got = {m for m in got if m[0] >= t}
                    has_j = True
            else:
                germ = ref[t:hit.r_st]
                fs = fs[:q0 - L] + js + germ + fs[q0:]
                fq = fq[:q0 - L] + jq + fq[q0] * len(germ) + fq[q0:]
                has_j = True
        if hit.strand == -1:
            fs, fq = fs.translate(_COMP)[::-1], fq[::-1]
        out.append((name, fs, fq))
        carried[name] = {"clone": c, "muts": got, "junction": has_j}
    return out, carried


# ---------------------------------------------------------------- scoring

def adjusted_rand(a: list, b: list) -> float | None:
    """Adjusted Rand index of two labellings of the same items."""
    n = len(a)
    if n < 2:
        return None
    c2 = lambda x: x * (x - 1) / 2   # noqa: E731
    pairs = Counter(zip(a, b, strict=True))
    sum_ij = sum(c2(v) for v in pairs.values())
    sum_a = sum(c2(v) for v in Counter(a).values())
    sum_b = sum(c2(v) for v in Counter(b).values())
    expected = sum_a * sum_b / c2(n)
    top = (sum_a + sum_b) / 2
    return 1.0 if top == expected else round((sum_ij - expected) / (top - expected), 4)


def qvalues(pvals: list[float], m: int | None = None) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, with m hypotheses (>= len(pvals))."""
    n = len(pvals)
    m = max(m or n, n)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    best = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        best = min(best, pvals[i] * m / rank)
        q[i] = best
    return q


ALPHAS = (1e-6, 1e-5, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)


def roc(calling: dict, truth_support: dict, min_reads: int = 2) -> list[dict]:
    """Recall (of mutations carried by >= min_reads analysed reads) and false
    calls as the FDR threshold sweeps over every candidate the test saw."""
    keys = calling.get("_candidate_keys") or []
    if not keys:
        return []
    pv = [k[2] for k in keys]
    bh = calling.get("multiple_testing") == "bh"
    q = qvalues(pv, calling.get("hypotheses")) if bh else pv
    target = {k for k, n in truth_support.items() if n >= min_reads}
    out = []
    for a in ALPHAS:
        called = {(pos - 1, alt) for (pos, alt, _p), qq in zip(keys, q, strict=True) if qq <= a}
        tp = len(called & target)
        fp = len(called - set(truth_support))
        out.append({"alpha": a, "recall": round(tp / len(target), 4) if target else None,
                    "true_calls": tp, "false_calls": fp})
    return out


SUPPORT_BINS = ((1, 1), (2, 2), (3, 5), (6, 10), (11, 10 ** 9))


def _bin(k: int) -> str:
    for a, b in SUPPORT_BINS:
        if a <= k <= b:
            return str(a) if a == b else (f"{a}-{b}" if b < 10 ** 9 else f">{a - 1}")
    return "0"


def score(res: dict, row: dict, cl: dict | None, carried: dict, p, lo: int, hi: int) -> dict:
    """Everything the simulation knows against what the caller reported."""
    analysed = [r.qname for r in res["_aligned"]]
    bp = hi - lo + 1
    support = Counter(m for n in analysed for m in carried.get(n, {}).get("muts", ())
                      if lo <= m[0] <= hi)
    called = {(v.pos - 1, v.alt) for v in res["variants"]}
    by_bin: dict = defaultdict(lambda: [0, 0])
    for m, k in support.items():
        by_bin[_bin(k)][1] += 1
        by_bin[_bin(k)][0] += m in called
    truth_set = set(support)
    tp, fp = len(called & truth_set), len(called - truth_set)

    # frequencies: the truth computed the way each estimate is defined
    n = len(analysed)
    true_read_weighted = 1e3 * sum(support.values()) / (n * bp) if n else None
    members = defaultdict(list)
    for name in analysed:
        if name in carried:
            members[carried[name]["clone"]].append(name)
    clone_true = []
    for _c, names in members.items():
        cnt = Counter(m for nm in names for m in carried[nm]["muts"] if lo <= m[0] <= hi)
        need = 1 if len(names) == 1 else (len(names) + 1) // 2
        clone_true.append(sum(1 for k in cnt.values() if k >= need))
    true_clone_weighted = 1e3 * sum(clone_true) / (len(clone_true) * bp) if clone_true else None

    out = {
        "reads_analysed": n,
        "planted_mutations": len(truth_set),
        "recall_by_support": {k: {"called": v[0], "planted": v[1],
                                  "recall": round(v[0] / v[1], 3) if v[1] else None}
                              for k, v in sorted(by_bin.items(), key=lambda kv: _bin_order(kv[0]))},
        "true_calls": tp, "false_calls": fp,
        "recall_2plus": _recall(called, support, 2),
        "fdr": round(fp / (tp + fp), 4) if tp + fp else 0.0,
        "roc": roc(res["summary"]["calling"], support),
        "true_read_weighted_x1e3": _r(true_read_weighted),
        "est_read_weighted_x1e3": row.get("mutation_frequency_x1e3"),
        "true_clone_weighted_x1e3": _r(true_clone_weighted),
        "est_clone_weighted_x1e3": row.get("clone_weighted_frequency_x1e3"),
        "true_clones_in_reads": len(members),
    }
    # densities: distinct mutations per base sequenced, the manuscript's measure, on
    # the estimator's own denominator (analysed positions x mean depth)
    s_call = res["summary"]["calling"]
    denom = (s_call.get("analysed_positions") or 0) * max(1.0, s_call.get("mean_depth") or 0.0)
    if denom:
        out.update({
            "true_unique_density_x1e5": _r(1e5 * len(truth_set) / denom),
            "est_unique_density_x1e5": row.get("unique_density_x1e5"),
            "est_combined_density_x1e5": row.get("combined_density_x1e5"),
            "true_single_read_mutations": sum(1 for k in support.values() if k == 1),
            "est_single_read_mutations": row.get("single_read_estimate"),
        })
    # clone recovery among reads with a junction the clusterer could read
    j0 = (p.junction_ref_start or 0) - 1
    keys = {}
    for r in res["_aligned"]:
        k = clones.junction(r, j0, p.junction_length) if p.junction_ref_start else None
        if k and r.qname in carried:
            keys[r.qname] = k
    if keys:
        found = clones.cluster(keys, p.junction_max_edits)
        names = sorted(keys)
        truth_lab = [carried[x]["clone"] for x in names]
        out.update({"reads_with_junction": len(keys),
                    "true_clones_with_junction": len(set(truth_lab)),
                    "clones_found": len(set(found.values())),
                    "adjusted_rand": adjusted_rand(truth_lab, [found[x] for x in names]),
                    "true_largest_pct": round(100 * max(Counter(truth_lab).values()) / len(keys), 1),
                    "est_largest_pct": cl["largest_clone_pct"] if cl else None})
    return out


def _bin_order(label: str) -> int:
    return int(label.lstrip(">").split("-")[0]) + (1 if label.startswith(">") else 0)


def _recall(called: set, support: Counter, k: int) -> float | None:
    t = {m for m, n in support.items() if n >= k}
    return round(len(called & t) / len(t), 4) if t else None


def _r(x):
    return None if x is None else round(x, 4)


# ---------------------------------------------------------------- running

def run_one(control: str, peers: list[str], reference: str, contaminants: list[str],
            region: str, scenario: dict, seed: int, mode: str = "recommended",
            model_path: str | None = None, place_by_model: bool = True,
            workdir: str | None = None, depth: int | None = None) -> dict:
    """Build one library from `control` and score it, with `peers` as the null."""
    p = project.build_params({"region": region, "mode": mode})
    ref = read_fasta(reference)[0][1]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)
    model = targeting.load(model_path) if model_path else None
    rng = random.Random(seed)
    reads = list(read_fastq(control))
    truth = design(len(reads), ref, lo, hi, scenario, rng,
                   model if place_by_model else None, p.junction_length,
                   (p.junction_ref_start or 53) - 1)
    lib, carried = build_library(reads, reference, truth, p, rng)
    workdir = workdir or tempfile.mkdtemp(prefix="shmtrace_lineage_")
    path = os.path.join(workdir, f"sim_{seed}.fastq")
    with open(path, "w") as fh:
        for name, seq, qual in lib:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
    proj = {"name": "sim", "region": region, "mode": mode, "reference": reference,
            "contaminants": contaminants, "depth": depth,
            "samples": [{"name": "SIM", "path": path, "group": "SIM"}]
            + [{"name": f"C{i}", "path": c, "group": "CTRL", "control": True}
               for i, c in enumerate(peers)]}
    if model_path:
        proj["targeting_model"] = model_path
    try:
        out = project.run_project(proj)
    finally:
        os.remove(path)
    row = next(s for s in out["samples"] if s["sample"] == "SIM")
    res = out["_details"]["SIM"]
    sc = score(res, row, (out.get("clones") or {}).get("SIM"), carried, p, lo, hi)
    sc.update({"seed": seed, "control": os.path.basename(control), "mode": mode, "depth": depth,
               "s5f_index": row.get("s5f_index"), "s5f_p": row.get("s5f_p"),
               "lineages": truth["lineages"]})
    return sc


def run_suite(controls: list[str], reference: str, contaminants: list[str], region: str,
              scenarios: list[str], modes: list[str], reps: int = 1,
              model_path: str | None = None, seed: int = 1,
              progress=lambda msg: None, place_by_model: bool = True) -> dict:
    """Every scenario x mode, each control held out in turn, `reps` draws each."""
    runs = []
    for name in scenarios:
        for mode in modes:
            for k, ctl in enumerate(controls):
                peers = [c for c in controls if c != ctl]
                for rep in range(reps):
                    s = seed + 1000 * rep + 10 * k
                    progress(f"{name} / {mode} / {os.path.basename(ctl)} / rep {rep}")
                    r = run_one(ctl, peers, reference, contaminants, region,
                                SCENARIOS[name], s, mode, model_path, place_by_model)
                    r["scenario"] = name
                    runs.append(r)
    return {"runs": runs, "summary": summarise(runs)}


def summarise(runs: list[dict]) -> list[dict]:
    rows = []
    for key in dict.fromkeys((r["scenario"], r["mode"]) for r in runs):
        rs = [r for r in runs if (r["scenario"], r["mode"]) == key]
        mean = lambda f: _mean([r.get(f) for r in rs])   # noqa: E731
        bins: dict = defaultdict(lambda: [0, 0])
        for r in rs:
            for b, v in r["recall_by_support"].items():
                bins[b][0] += v["called"]
                bins[b][1] += v["planted"]
        tp = sum(r["true_calls"] for r in rs)
        fp = sum(r["false_calls"] for r in rs)
        roc_rows = []
        for i, a in enumerate(ALPHAS):
            pts = [r["roc"][i] for r in rs if len(r["roc"]) > i]
            if pts:
                roc_rows.append({"alpha": a,
                                 "recall": _mean([x["recall"] for x in pts]),
                                 "false_calls": _mean([x["false_calls"] for x in pts])})
        rows.append({
            "scenario": key[0], "mode": key[1], "runs": len(rs),
            "recall_by_support": {b: round(v[0] / v[1], 3) if v[1] else None
                                  for b, v in sorted(bins.items(), key=lambda kv: _bin_order(kv[0]))},
            "recall_2plus": mean("recall_2plus"),
            "false_calls_per_library": round(fp / len(rs), 2),
            "fdr": round(fp / (tp + fp), 4) if tp + fp else 0.0,
            "adjusted_rand": mean("adjusted_rand"),
            "clones_found": mean("clones_found"),
            "true_clones_with_junction": mean("true_clones_with_junction"),
            "true_largest_pct": mean("true_largest_pct"),
            "est_largest_pct": mean("est_largest_pct"),
            "true_read_weighted_x1e3": mean("true_read_weighted_x1e3"),
            "est_read_weighted_x1e3": mean("est_read_weighted_x1e3"),
            "true_clone_weighted_x1e3": mean("true_clone_weighted_x1e3"),
            "est_clone_weighted_x1e3": mean("est_clone_weighted_x1e3"),
            "true_unique_density_x1e5": mean("true_unique_density_x1e5"),
            "est_unique_density_x1e5": mean("est_unique_density_x1e5"),
            "est_combined_density_x1e5": mean("est_combined_density_x1e5"),
            "s5f_index": mean("s5f_index"),
            "roc": roc_rows,
        })
    return rows


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 4) if xs else None


def report_markdown(res: dict) -> str:
    rows = res["summary"]
    bins = list(dict.fromkeys(b for r in rows for b in r["recall_by_support"]))
    f = lambda x, d=2: "–" if x is None else f"{x:.{d}f}"   # noqa: E731
    out = ["# Lineage simulation", "",
           "Recall by the number of analysed reads carrying a planted mutation; false "
           "calls are calls never planted.", "",
           "| Scenario | Mode | Runs | " + " | ".join(f"{b} reads" for b in bins)
           + " | False calls | FDR |", "|---|---|---:|" + "---:|" * (len(bins) + 2)]
    for r in rows:
        out.append(f"| {r['scenario']} | {r['mode']} | {r['runs']} | "
                   + " | ".join(f(r["recall_by_support"].get(b)) for b in bins)
                   + f" | {r['false_calls_per_library']} | {f(r['fdr'])} |")
    out += ["", "Clones and frequencies (mutations per 10³ bp), true against estimated.", "",
            "| Scenario | Mode | Adjusted Rand | Clones found / present | Largest clone % "
            "| Read-weighted | Clone-weighted | S5F index |",
            "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        out.append(
            f"| {r['scenario']} | {r['mode']} | {f(r['adjusted_rand'])} | "
            f"{f(r['clones_found'], 1)} / {f(r['true_clones_with_junction'], 1)} | "
            f"{f(r['est_largest_pct'], 0)} / {f(r['true_largest_pct'], 0)} | "
            f"{f(r['est_read_weighted_x1e3'])} / {f(r['true_read_weighted_x1e3'])} | "
            f"{f(r['est_clone_weighted_x1e3'])} / {f(r['true_clone_weighted_x1e3'])} | "
            f"{f(r['s5f_index'])} |")
    return "\n".join(out)


def plot_roc(res: dict, path: str) -> None:
    """Recall of mutations in two or more reads against false calls per library."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [r for r in res["summary"] if r["scenario"] != "null" and r["roc"]]
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    # one colour per scenario, one line style per mode, so each pair compares
    scen = list(dict.fromkeys(r["scenario"] for r in rows))
    modes = list(dict.fromkeys(r["mode"] for r in rows))
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    styles = ("-", "--", ":", "-.")
    for r in rows:
        xs = [max(0.1, x["false_calls"]) for x in r["roc"]]
        ys = [x["recall"] for x in r["roc"]]
        line, = ax.plot(xs, ys, marker="o", ms=3, label=f"{r['scenario']} · {r['mode']}",
                        color=colours[scen.index(r["scenario"]) % len(colours)],
                        ls=styles[modes.index(r["mode"]) % len(styles)])
        k = ALPHAS.index(0.05)
        if k < len(r["roc"]):
            ax.plot(xs[k], ys[k], marker="s", ms=8, mfc="none", color=line.get_color())
    ax.set_xscale("log")
    ax.set_xlabel("False calls per library (log)")
    ax.set_ylabel("Recall, mutations in ≥2 reads")
    ax.set_ylim(0, 1)
    ax.set_title("Threshold sweep; open square = FDR 0.05", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Simulate libraries with known clones and lineages from mutation-free "
                    "controls, run them end to end, and score calls, clones and frequencies")
    ap.add_argument("--control", action="append", required=True,
                    help="mutation-free reads (repeat; each is held out in turn, the rest "
                         "form the null; at least two)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--contaminant", action="append", default=[])
    ap.add_argument("--region", required=True, help="preset, e.g. jh4_intron")
    ap.add_argument("--scenarios", default="polyclonal,expanded,low,null",
                    help=f"comma-separated, from {', '.join(SCENARIOS)}")
    ap.add_argument("--modes", default="recommended,geneious")
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--targeting-model", help="5-mer model CSV: places mutations and scores them")
    ap.add_argument("--place", choices=("model", "motif"), default="model",
                    help="draw mutation positions from the targeting model (default, when "
                         "one is given) or from WRC/GYW and WA/TW weights; 'motif' with a "
                         "model tests the targeting index on mutations it did not place")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="lineage_sim")
    a = ap.parse_args()
    if len(a.control) < 2:
        ap.error("give at least two --control files: one is planted, the others are the null")
    res = run_suite(a.control, a.reference, a.contaminant, a.region,
                    a.scenarios.split(","), a.modes.split(","), a.reps,
                    a.targeting_model, a.seed, progress=lambda m: print(m, flush=True),
                    place_by_model=a.place == "model")
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "lineage_sim.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: sorted(o) if isinstance(o, set) else str(o))
    md = report_markdown(res)
    with open(os.path.join(a.out, "report.md"), "w") as fh:
        fh.write(md)
    plot_roc(res, os.path.join(a.out, "roc.png"))
    print(md)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------- UMI libraries

def _random_umi(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def build_umi_library(reads: list, reference: str, truth: dict, p, rng: random.Random,
                      copies_mean: float = 6.0, umi_len: int = 18, ends: str = "both",
                      umi_q: int = 15) -> tuple[list, dict]:
    """Molecules with UMIs, each sequenced several times.

    `truth` is a design() over molecules. Each molecule's copies are *different*
    donor reads with the molecule's mutations (and junction) written in, so every
    copy carries its own real sequencing errors. The molecule's UMI is written
    just outside the reference ends with errors at Phred `umi_q` and 2% indels,
    as a nanopore read of a UMI would show. A copy whose read does not reach a
    reference end gets no UMI there, as in a real library. Returns the reads and,
    per read, its molecule and the mutations it carries.
    """
    ref = read_fasta(reference)[0][1]
    n_mol = len(truth["read_clone"])
    order = list(range(len(reads)))
    rng.shuffle(order)
    assign, k = {}, 0
    for m in range(n_mol):
        copies = 1 + _poisson(max(copies_mean - 1, 0), rng)
        for _ in range(copies):
            if k >= len(order):
                break
            assign[order[k]] = m
            k += 1
    donors = sorted(assign)
    sub = [reads[i] for i in donors]
    read_truth = {"read_clone": {j: truth["read_clone"][assign[i]] for j, i in enumerate(donors)},
                  "read_muts": {j: truth["read_muts"].get(assign[i], set())
                                for j, i in enumerate(donors)},
                  "junctions": truth["junctions"], "jstarts": truth["jstarts"]}
    lib, carried = build_library(sub, reference, read_truth, p, rng)
    umis = {m: (_random_umi(umi_len, rng), _random_umi(umi_len, rng)) for m in range(n_mol)}
    aligner, _t, _c = build_index([reference], [], p.minimap2_preset)
    qch = chr(umi_q + 33)
    out = []
    for j, (name, seq, qual) in enumerate(lib):
        m = assign[donors[j]]
        carried[name]["molecule"] = m
        carried[name]["umi"] = umis[m]
        hit = next((h for h in aligner.map(seq) if h.is_primary), None)
        if hit is None:
            out.append((name, seq, qual))
            continue
        fs = seq if hit.strand == 1 else seq.translate(_COMP)[::-1]
        fq = qual if hit.strand == 1 else qual[::-1]
        q_st = hit.q_st if hit.strand == 1 else len(seq) - hit.q_en
        q_en = hit.q_en if hit.strand == 1 else len(seq) - hit.q_st
        left, right = umis[m]
        if ends in ("right", "both") and hit.r_en >= len(ref) - 3:
            at = q_en + (len(ref) - hit.r_en)
            u, uq = _noisy(right, qch * umi_len, rng, 0.02)
            fs, fq = fs[:at] + u + fs[at:], fq[:at] + uq + fq[at:]
        if ends in ("left", "both") and hit.r_st <= 3:
            at = max(0, q_st - hit.r_st)
            u, uq = _noisy(left, qch * umi_len, rng, 0.02)
            fs, fq = fs[:at] + u + fs[at:], fq[:at] + uq + fq[at:]
        out.append((name, fs, fq))
    return out, carried


def score_umi(res: dict, row: dict, carried: dict, lo: int, hi: int) -> dict:
    """Score a UMI run: each consensus read stands for the molecule most of its reads
    came from; truth is that molecule's mutations at positions the consensus covers."""
    bp = hi - lo + 1
    mol_muts: dict = defaultdict(set)
    for c in carried.values():
        if "molecule" in c:
            mol_muts[c["molecule"]] |= {m for m in c["muts"] if lo <= m[0] <= hi}
    support, purity, per_mol = Counter(), [], []
    for r in res["_aligned"]:
        mem = [carried[n]["molecule"] for n in getattr(r, "members", []) if n in carried]
        if not mem:
            continue
        m, k = Counter(mem).most_common(1)[0]
        purity.append(k / len(mem))
        covered = {rp for rp, _b, _q in r.pairs}
        muts = {x for x in mol_muts[m] if x[0] in covered}
        per_mol.append(len(muts))
        support.update(muts)
    called = {(v.pos - 1, v.alt) for v in res["variants"]}
    truth_set = set(support)
    tp, fp = len(called & truth_set), len(called - truth_set)
    n = len(per_mol)
    s_call = res["summary"]["calling"]
    denom = (s_call.get("analysed_positions") or 0) * max(1.0, s_call.get("mean_depth") or 0.0)
    return {"molecules_analysed": n, "planted_mutations": len(truth_set),
            "true_calls": tp, "false_calls": fp,
            "recall": round(tp / len(truth_set), 4) if truth_set else None,
            "recall_single_molecule": _recall_exact(called, support, 1),
            "group_purity": round(sum(purity) / len(purity), 4) if purity else None,
            "true_read_weighted_x1e3": _r(1e3 * sum(per_mol) / (n * bp)) if n else None,
            "est_read_weighted_x1e3": row.get("mutation_frequency_x1e3"),
            "true_unique_density_x1e5": _r(1e5 * len(truth_set) / denom) if denom else None,
            "est_unique_density_x1e5": row.get("unique_density_x1e5"),
            "est_combined_density_x1e5": row.get("combined_density_x1e5")}


def _recall_exact(called: set, support: Counter, k: int) -> float | None:
    t = {m for m, n in support.items() if n == k}
    return round(len(called & t) / len(t), 4) if t else None


def score_umi_reads(res: dict, row: dict, carried: dict, lo: int, hi: int) -> dict:
    """Score a UMI library analysed read by read, against the molecules' truth:
    calls against the molecules' mutations; frequency per read (each copy counts),
    and unique density on the read denominator."""
    bp = hi - lo + 1
    support, per = Counter(), []
    for r in res["_aligned"]:
        c = carried.get(r.qname)
        if not c:
            continue
        covered = {rp for rp, _b, _q in r.pairs}
        muts = {m for m in c["muts"] if lo <= m[0] <= hi and m[0] in covered}
        per.append(len(muts))
        support.update(muts)
    mol_support = Counter()
    seen = defaultdict(set)
    for r in res["_aligned"]:
        c = carried.get(r.qname)
        if c and "molecule" in c:
            for m in c["muts"]:
                if lo <= m[0] <= hi:
                    seen[m].add(c["molecule"])
    for m, mols in seen.items():
        mol_support[m] = len(mols)
    called = {(v.pos - 1, v.alt) for v in res["variants"]}
    truth_set = set(support)
    tp, fp = len(called & truth_set), len(called - truth_set)
    n = len(per)
    s_call = res["summary"]["calling"]
    denom = (s_call.get("analysed_positions") or 0) * max(1.0, s_call.get("mean_depth") or 0.0)
    return {"reads_analysed": n, "planted_mutations": len(truth_set),
            "true_calls": tp, "false_calls": fp,
            "recall": round(tp / len(truth_set), 4) if truth_set else None,
            "recall_single_molecule": _recall_exact(called, mol_support, 1),
            "true_read_weighted_x1e3": _r(1e3 * sum(per) / (n * bp)) if n else None,
            "est_read_weighted_x1e3": row.get("mutation_frequency_x1e3"),
            "true_unique_density_x1e5": _r(1e5 * len(truth_set) / denom) if denom else None,
            "est_unique_density_x1e5": row.get("unique_density_x1e5"),
            "est_combined_density_x1e5": row.get("combined_density_x1e5")}
