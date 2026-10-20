"""CASE vs CONTROL comparative SHM analysis — the SHMTool-superset layer.

Unlike SHMTool, we start from raw FASTQ (each file = one independent source, e.g.
one mouse): every file is run through the single-sample engine, then mutations
are pooled per category and compared with base-composition-corrected frequencies
and proportion tests, exactly in the spirit of MacCarthy et al. (2009).

Reported for CONTROL vs CASE (unique mutations, base-composition corrected):
  * overall mutation frequency
  * A:T vs G:C mutation frequencies (and per-base)
  * 12-class from->to spectrum
  * AID / Pol-eta hotspot and coldspot frequencies (+ WRC-vs-GYW strand bias)
  * optional subregion S vs its complement S'
Each comparison gets a two-proportion test (R prop.test equivalent) with stars.
Non-unique ("total") counts are also reported, with the no-UMI caveat.

LEGACY, DESCRIPTIVE ROUTE. This command pools every file's mutations and tests
sites or reads as if independent, and it calls each file without a
control-fitted error model (no empirical null, no decoy thresholding). Its p
values and stars describe pooled counts; they are not evidence at the level of
biological replicates, and the calling is not the calibrated analysis the paper
describes. For that, name the background controls in a sample sheet and use
`shmtrace-sheet run` (or a project in `shmtrace-app`): leave-one-out control
fitting, control-based thresholding, and animal_tests with one value per animal.
Every report from this command says so, and lists the error model each file used.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from scipy.stats import chi2_contingency

from . import batch, motifs, pipeline
from .config import Params, load_preset
from .io_utils import read_fasta

BASES = ("A", "C", "G", "T")


# ---- statistics -------------------------------------------------------------
def prop_test(x1: int, L1: int, x2: int, L2: int) -> float:
    """Two-proportion chi-square with continuity correction (R prop.test)."""
    if min(L1, L2) == 0 or (x1 + x2) == 0:
        return 1.0
    table = [[x1, max(0, L1 - x1)], [x2, max(0, L2 - x2)]]
    try:
        _, p, _, _ = chi2_contingency(table, correction=True)
        return float(p)
    except ValueError:
        return 1.0


def stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


# ---- per-category aggregation ----------------------------------------------
@dataclass
class Category:
    name: str
    n_files: int
    unique: set          # {(pos0, ref, alt)} pooled across files
    nonunique: int       # summed supporting reads
    total_reads: int


class Unique(list):
    """Unique mutations per file, as SHMTool counts them.

    SHMTool counts an identical mutation once *within* a file (one independent
    source, e.g. one mouse) but once per file across files, with L = sites x
    files. So the counts are kept per file and summed, not pooled into one set.
    """

    def __len__(self):                         # total unique mutations, summed over files
        return sum(len(s) for s in list.__iter__(self))

    def __iter__(self):                        # every (file's) unique mutation
        for s in list.__iter__(self):
            yield from s

    def files(self):
        return list.__iter__(self)


def _collect(results: list[dict], restrict: set | None):
    uniq, nonuniq, reads = Unique(), {}, 0
    for res in results:
        reads += res["summary"]["reads"]["aligned_to_target"]
        mine = set()
        for v in res["variants"]:
            if restrict is not None and (v.pos - 1) not in restrict:
                continue
            mine.add((v.pos - 1, v.ref, v.alt))
            key = (v.pos - 1, v.ref, v.alt)
            nonuniq[key] = nonuniq.get(key, 0) + v.supporting_reads
        list.append(uniq, mine)
    return uniq, nonuniq, reads


def _run_groups(control_files, case_files, reference, contaminants, p):
    """Prepare every file together (shared haplotype evidence, common depth)."""
    samples = {f"control_{i}": f for i, f in enumerate(control_files)}
    samples.update({f"case_{i}": f for i, f in enumerate(case_files)})
    b = batch.prepare_batch(samples, reference, contaminants, p)
    out = {name: pipeline.analyse(b["reads"][name], reference, p, b["stats"][name],
                                  depth=b["depth"])
           for name in samples}
    ctrl = [out[n] for n in samples if n.startswith("control_")]
    case = [out[n] for n in samples if n.startswith("case_")]
    b["calling"] = _calling_provenance(out)
    return ctrl, case, b


def _site_positions(ref: str, window, restrict: set | None) -> list[int]:
    lo = (window[0] - 1) if window else 0
    hi = (window[1] - 1) if window else len(ref) - 1
    return [i for i in range(lo, hi + 1)
            if (restrict is None or i in restrict)]


def _base_site_counts(ref: str, positions: list[int]) -> dict[str, int]:
    c = {b: 0 for b in BASES}
    for i in positions:
        b = ref[i].upper()
        if b in c:
            c[b] += 1
    return c


ROUTE = "legacy pooled descriptive comparison (SHMTool-style)"
LEGACY_NOTE = ("Legacy pooled comparison: mutations are pooled across files and every p value "
               "is a two-proportion test that treats sites or reads as independent. These p "
               "values and stars are descriptive; they are not evidence at the level of "
               "biological replicates. Files were called without a control-fitted error model. "
               "For the calibrated analysis with one value per animal, name the controls in a "
               "sample sheet and run `shmtrace-sheet run` (or use a project in shmtrace-app).")


def _calling_provenance(named: dict) -> dict:
    """Which error model and multiple-testing rule each file was actually called with,
    and every calling note (a requested model that fell back, for example)."""
    per_file = {}
    for name, res in named.items():
        c = res["summary"]["calling"]
        per_file[name] = {"error_model": c.get("error_model"),
                          "multiple_testing": c.get("multiple_testing"),
                          "notes": list(c.get("notes") or [])}
    models = Counter(v["error_model"] for v in per_file.values())
    return {"route": ROUTE, "error_models_used": dict(models), "per_file": per_file}


def _freq_row(label, ctrl_x, ctrl_L, case_x, case_L):
    fc = ctrl_x / ctrl_L if ctrl_L else 0.0
    fa = case_x / case_L if case_L else 0.0
    p = prop_test(ctrl_x, ctrl_L, case_x, case_L)
    direction = "↓" if fa < fc else ("↑" if fa > fc else "=")
    return {"metric": label, "control_freq": fc, "case_freq": fa,
            "control_n": ctrl_x, "case_n": case_x,
            "p_value": p, "sig": stars(p), "change": direction}


def compare(control_files: list[str], case_files: list[str],
            reference: str, contaminants: list[str], p: Params,
            subregion: list[int] | None = None) -> dict:
    ref = read_fasta(reference)[0][1]
    restrict = set(x - 1 for x in subregion) if subregion else None

    ctrl_res, case_res, b = _run_groups(control_files, case_files, reference, contaminants, p)
    cu, cnu_sites, cr = _collect(ctrl_res, restrict)
    au, anu_sites, ar = _collect(case_res, restrict)
    cnu, anu = sum(cnu_sites.values()), sum(anu_sites.values())
    # reads analysed per group = depth x files (each file analysed at b["depth"])
    cQ, aQ = b["depth"] * len(control_files), b["depth"] * len(case_files)
    Gc, Ga = len(control_files), len(case_files)

    positions = _site_positions(ref, p.analysis_window, restrict)
    base_sites = _base_site_counts(ref, positions)
    N_total = sum(base_sites.values())

    def uc(uniq, froms=None, tos=None):     # unique count filtered
        return sum(1 for (_p, r, a) in uniq
                   if (froms is None or r in froms) and (tos is None or a in tos))

    rows = []
    # overall (each site has 3 possible substitutions -> x3)
    rows.append(_freq_row("Overall mutation frequency",
                          len(cu), N_total * Gc * 3, len(au), N_total * Ga * 3))
    # A:T vs G:C consolidated
    rows.append(_freq_row("A:T mutations",
                          uc(cu, froms={"A", "T"}), (base_sites["A"] + base_sites["T"]) * Gc * 3,
                          uc(au, froms={"A", "T"}), (base_sites["A"] + base_sites["T"]) * Ga * 3))
    rows.append(_freq_row("G:C mutations",
                          uc(cu, froms={"G", "C"}), (base_sites["G"] + base_sites["C"]) * Gc * 3,
                          uc(au, froms={"G", "C"}), (base_sites["G"] + base_sites["C"]) * Ga * 3))
    # G:C transitions (C->T, G->A) — the AID replication-bypass signature
    rows.append(_freq_row("G:C transitions (C→T, G→A)",
                          uc(cu, froms={"C"}, tos={"T"}) + uc(cu, froms={"G"}, tos={"A"}),
                          (base_sites["C"] + base_sites["G"]) * Gc,
                          uc(au, froms={"C"}, tos={"T"}) + uc(au, froms={"G"}, tos={"A"}),
                          (base_sites["C"] + base_sites["G"]) * Ga))

    # Non-unique ("total") mutations: every read counts, L = sites x reads. The
    # three substitutions at a site are mutually exclusive within one read, so
    # L is not multiplied by 3 (MacCarthy et al. 2009). Without UMIs, PCR
    # duplicates inflate these counts, so read them as an upper bound.
    def nu(sites, froms=None, tos=None):
        return sum(n for (_p, r, a), n in sites.items()
                   if (froms is None or r in froms) and (tos is None or a in tos))
    nonunique_rows = [
        _freq_row("Overall (non-unique)", cnu, N_total * cQ, anu, N_total * aQ),
        _freq_row("A:T (non-unique)", nu(cnu_sites, {"A", "T"}),
                  (base_sites["A"] + base_sites["T"]) * cQ,
                  nu(anu_sites, {"A", "T"}), (base_sites["A"] + base_sites["T"]) * aQ),
        _freq_row("G:C (non-unique)", nu(cnu_sites, {"G", "C"}),
                  (base_sites["G"] + base_sites["C"]) * cQ,
                  nu(anu_sites, {"G", "C"}), (base_sites["G"] + base_sites["C"]) * aQ),
    ]

    # motif hot/cold-spots
    motif_sites = motifs.site_counts(ref, positions)
    ctrl_pos = [pp for (pp, _r, _a) in cu]
    case_pos = [pp for (pp, _r, _a) in au]
    ctrl_hits = motifs.score_variant_positions(ref, ctrl_pos)
    case_hits = motifs.score_variant_positions(ref, case_pos)
    motif_rows = []
    for cat, members in motifs.CATEGORIES.items():
        cx = sum(ctrl_hits[m] for m in members)
        ax = sum(case_hits[m] for m in members)
        cL = sum(motif_sites[m] for m in members) * Gc
        aL = sum(motif_sites[m] for m in members) * Ga
        motif_rows.append({**_freq_row(cat, cx, cL, ax, aL), "members": members})
    # strand bias: WRC (top) vs GYW (bottom)
    top, bot = motifs.STRAND_PAIR
    strand = {
        "control_WRC": ctrl_hits[top], "control_GYW": ctrl_hits[bot],
        "case_WRC": case_hits[top], "case_GYW": case_hits[bot],
        "note": "similar WRC and GYW mutation counts imply no strand bias",
    }

    # 12-class spectrum (unique counts, both categories), one 2x2 test per
    # from->to class with L = sites of the from-base x files (SHMTool)
    spec = {}
    for r in BASES:
        spec[r] = {}
        for a in BASES:
            if a == r:
                continue
            cx, ax = uc(cu, {r}, {a}), uc(au, {r}, {a})
            row = _freq_row(f"{r}>{a}", cx, base_sites[r] * Gc, ax, base_sites[r] * Ga)
            spec[r][a] = {"control": cx, "case": ax, "p_value": row["p_value"],
                          "sig": row["sig"], "control_freq": row["control_freq"],
                          "case_freq": row["case_freq"]}

    calling = b["calling"]
    warnings = [LEGACY_NOTE]
    fallbacks = sorted({n for v in calling["per_file"].values() for n in v["notes"]
                        if "used the read-quality model instead" in n})
    if fallbacks:
        warnings.append("Requested error model unavailable without control reads; files were "
                        "called with: " + ", ".join(f"{k} ({v} file(s))" for k, v in
                                                     sorted(calling["error_models_used"].items()))
                        + ". " + " ".join(fallbacks))
    return {
        "analysis_route": ROUTE,
        "warnings": warnings,
        "calling": calling,
        "params": {"preset": p.name, "analysis_window": p.analysis_window,
                   "subregion_sites": len(restrict) if restrict else None,
                   "common_depth": b["depth"],
                   "haplotype_blocks_dropped": len(b["known_blocks"])},
        "groups": {"control_files": len(control_files), "case_files": len(case_files),
                   "control_reads": cr, "case_reads": ar,
                   "control_unique": len(cu), "case_unique": len(au),
                   "control_nonunique": cnu, "case_nonunique": anu},
        "base_composition": base_sites,
        "frequency_rows": rows,
        "nonunique_rows": nonunique_rows,
        "motif_rows": motif_rows,
        "strand_bias": strand,
        "spectrum": spec,
        "_control_pos": ctrl_pos, "_case_pos": case_pos,     # for plots
    }


def report_markdown(res: dict) -> str:
    g = res["groups"]
    out = ["# SHM comparative analysis — CONTROL vs CASE (legacy, pooled, descriptive)", ""]
    out += [f"> **Warning:** {w}" + "\n" for w in res.get("warnings", [])]
    models = (res.get("calling") or {}).get("error_models_used")
    if models:
        out += ["Error model used per file: "
                + ", ".join(f"{k} × {v}" for k, v in sorted(models.items())) + ".", ""]
    out += [
           f"Preset **{res['params']['preset']}**, window "
           f"{res['params']['analysis_window']}. "
           f"Control = {g['control_files']} file(s), {g['control_reads']} reads; "
           f"Case = {g['case_files']} file(s), {g['case_reads']} reads. "
           f"Every file analysed at {res['params'].get('common_depth')} reads.", "",
           "Frequencies are base-composition corrected (unique mutations). "
           "Pooled two-proportion tests (sites or reads treated as independent; "
           "not replicate-level): * p<0.05, ** p<0.01, *** p<0.001.", "",
           "## Mutation profile", "",
           "| Metric | Control | Case | Δ | pooled p |",
           "|---|---:|---:|:--:|---|"]
    for r in res["frequency_rows"]:
        out.append(f"| {r['metric']} | {r['control_freq']:.2e} | "
                   f"{r['case_freq']:.2e} | {r['change']} | {r['sig']} |")
    out += ["", "## Non-unique mutations (every read counted; no-UMI upper bound)", "",
            "| Metric | Control | Case | Δ | pooled p |", "|---|---:|---:|:--:|---|"]
    for r in res.get("nonunique_rows", []):
        out.append(f"| {r['metric']} | {r['control_freq']:.2e} | "
                   f"{r['case_freq']:.2e} | {r['change']} | {r['sig']} |")
    out += ["", "## Substitution classes (unique, per from-base site)", "",
            "| Class | Control | Case | Δ | pooled p |", "|---|---:|---:|:--:|---|"]
    for r, row in res.get("spectrum", {}).items():
        for a, c in row.items():
            if "p_value" not in c:
                continue
            ch = "↑" if c["case_freq"] > c["control_freq"] else (
                "↓" if c["case_freq"] < c["control_freq"] else "=")
            out.append(f"| {r}→{a} | {c['control_freq']:.2e} | {c['case_freq']:.2e} | "
                       f"{ch} | {c['sig']} |")
    out += ["", "## Motif targeting (hot/cold-spots)", "",
            "| Motif class | Control | Case | Δ | pooled p |",
            "|---|---:|---:|:--:|---|"]
    for r in res["motif_rows"]:
        out.append(f"| {r['metric']} ({'/'.join(r['members'])}) | "
                   f"{r['control_freq']:.2e} | {r['case_freq']:.2e} | "
                   f"{r['change']} | {r['sig']} |")
    s = res["strand_bias"]
    out += ["", "## Strand bias (AID hotspot top vs bottom strand)", "",
            f"Control WRC/GYW = {s['control_WRC']}/{s['control_GYW']}; "
            f"Case WRC/GYW = {s['case_WRC']}/{s['case_GYW']}. {s['note']}."]
    return "\n".join(out)


def main():
    import argparse
    import json
    import os
    ap = argparse.ArgumentParser(
        description="Legacy CASE vs CONTROL comparison: pooled, SHMTool-style, descriptive. "
                    "Files are called without a control-fitted error model and p values pool "
                    "files. For calibrated, replicate-level inference use `shmtrace-sheet run`.")
    ap.add_argument("--control", nargs="+", required=True, help="control FASTQ(s)")
    ap.add_argument("--case", nargs="+", required=True, help="case FASTQ(s)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--contaminant", action="append", default=[])
    ap.add_argument("--preset", required=True)
    ap.add_argument("--subregion", help="e.g. 200-400 (1-based, within window)")
    ap.add_argument("--out", default="compare_results")
    ap.add_argument("--plots", action="store_true")
    ap.add_argument("--mode", default="recommended",
                    help="calling profile, as for shmtrace --mode")
    args = ap.parse_args()

    p = load_preset(args.preset, **pipeline.mode_overrides(args.mode))
    sub = None
    if args.subregion:
        a, b = args.subregion.split("-")
        sub = list(range(int(a), int(b) + 1))
    import sys
    print("shmtrace-compare: " + LEGACY_NOTE, file=sys.stderr)
    res = compare(args.control, args.case, args.reference, args.contaminant, p, sub)

    os.makedirs(args.out, exist_ok=True)
    clean = {k: v for k, v in res.items() if not k.startswith("_")}
    json.dump(clean, open(os.path.join(args.out, "comparison.json"), "w"), indent=2)
    open(os.path.join(args.out, "report.md"), "w").write(report_markdown(res))
    if args.plots:
        from . import viz
        viz.comparison_plots(res, args.out)
    print(report_markdown(res))
    print(f"\nWrote {args.out}/comparison.json, report.md"
          + (", plots" if args.plots else ""))


if __name__ == "__main__":
    main()
