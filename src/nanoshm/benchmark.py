"""Benchmark nanoshm against Geneious on a set of samples.

    nanoshm-benchmark --sheet samples.csv --reference ref.fasta \
        --contaminant contam.fasta --preset jh4_intron --out bench/

The sample sheet is a CSV with one row per library:

    sample  fastq            null  run     geneious_calls  geneious_export
    Wt2     WT2.fastq        no    DHJGKB  161             exports/Wt2.csv
    KO1     KO1.fastq        yes   DHJGKB  2               exports/KO1.csv

Only `sample` and `fastq` are required. `null` marks a mutation-free control
(e.g. Aicda-/-): those rows supply the empirical null for the others, fitted
leave-one-out so a control is never tested against itself, and only from the
same `run`. `geneious_calls` is a per-sample count to compare against;
`geneious_export` is a Geneious variant table, which enables site-level
precision, recall and Jaccard.

Every mode is run at one common depth (see batch.py), so counts and densities
are comparable across samples and modes.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys

from . import batch, null_model, pipeline
from .config import Params, load_preset
from .io_utils import read_fasta

MODES = ("recommended", "geneious")
_CHANGE = re.compile(r"^\s*([ACGT])\s*(?:->|→|>)\s*([ACGT])\s*$", re.I)


# ---- reading inputs ---------------------------------------------------------
def read_sheet(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        rows = list(csv.DictReader(fh, dialect=dialect))
    base = os.path.dirname(os.path.abspath(path))
    out = []
    for r in rows:
        r = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
        if not r.get("sample") or not r.get("fastq"):
            continue
        for key in ("fastq", "geneious_export"):
            if r.get(key) and not os.path.isabs(r[key]):
                r[key] = os.path.join(base, r[key])
        r["null"] = r.get("null", "").lower() in ("yes", "y", "true", "1", "null", "control")
        r["geneious_calls"] = int(r["geneious_calls"]) if r.get("geneious_calls") else None
        out.append(r)
    return out


def _num(x: str) -> float | None:
    x = (x or "").strip().replace(",", "")
    if not x:
        return None
    try:
        return float(x[:-1]) / 100 if x.endswith("%") else float(x)
    except ValueError:
        return None


def read_geneious_export(path: str, offset: int = 0, max_p: float | None = 0.05,
                         max_freq: float | None = 0.5,
                         window: tuple[int, int] | None = None) -> set[tuple[int, str, str]]:
    """Parse a Geneious variant table into {(1-based pos, ref, alt)}.

    Geneious "Find Variations/SNPs" exports usually list every candidate, so the
    manuscript's filters are re-applied here: `Variant P-Value (approximate)`
    <= max_p and `Variant Frequency` < max_freq (pass None to skip either), and
    an optional 1-based inclusive `window`.

    Accepts the CSV/TSV that Geneious writes when the annotations of a "Find
    Variations/SNPs" track are exported (columns such as Minimum, Change,
    Polymorphism Type), and also nanoshm's own variants.csv (position, ref,
    alt). `offset` is added to positions when the Geneious reference frame
    starts elsewhere.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        rows = list(csv.DictReader(fh, dialect=dialect))
    calls = set()
    for r in rows:
        r = {(k or "").strip().lower(): (v or "").strip() for k, v in r.items()}
        ptype = r.get("polymorphism type", r.get("type", "")).lower()
        if ptype and not any(t in ptype for t in ("snp", "substitution", "transition",
                                                  "transversion", "polymorphism")):
            continue
        if r.get("ref") and r.get("alt") and r.get("position"):
            pos, ref, alt = int(float(r["position"])), r["ref"].upper(), r["alt"].upper()
        else:
            m = _CHANGE.match(r.get("change", ""))
            start = r.get("minimum") or r.get("position") or r.get("min")
            if not m or not start:
                continue
            pos, ref, alt = int(float(start.replace(",", ""))), m.group(1).upper(), m.group(2).upper()
        pv = _num(r.get("variant p-value (approximate)", ""))
        fr = _num(r.get("variant frequency", ""))
        if max_p is not None and pv is not None and pv > max_p:
            continue
        if max_freq is not None and fr is not None and fr >= max_freq:
            continue
        pos += offset
        if window and not (window[0] <= pos <= window[1]):
            continue
        if len(ref) == 1 and len(alt) == 1 and ref != alt:
            calls.add((pos, ref, alt))
    return calls


# ---- scoring ----------------------------------------------------------------
def site_concordance(ours: set, theirs: set) -> dict:
    shared = len(ours & theirs)
    return {"shared": shared, "nanoshm_only": len(ours - theirs),
            "geneious_only": len(theirs - ours),
            "precision_vs_geneious": round(shared / len(ours), 3) if ours else None,
            "recall_vs_geneious": round(shared / len(theirs), 3) if theirs else None,
            "jaccard": round(shared / len(ours | theirs), 3) if (ours | theirs) else None}


def _mode_params(base: Params, mode: str) -> Params:
    from dataclasses import replace
    return replace(base, **pipeline.mode_overrides(mode))


# ---- running ----------------------------------------------------------------
def run_benchmark(sheet: list[dict], reference: str, contaminants: list[str], base: Params,
                  modes=MODES, depth: int | None = None, offset: int = 0) -> dict:
    ref_seq = read_fasta(reference)[0][1]
    lo, hi = pipeline._window(base, ref_seq)
    samples = {r["sample"]: r["fastq"] for r in sheet}
    controls = {r["sample"] for r in sheet if r["null"]}
    runs = {r["sample"]: r.get("run") or null_model.run_token(r["fastq"]) for r in sheet}
    win = base.analysis_window or (1, len(ref_seq))
    if base.end_trim_bp:
        win = (max(win[0], base.end_trim_bp + 1), min(win[1], len(ref_seq) - base.end_trim_bp))
    exports = {r["sample"]: read_geneious_export(r["geneious_export"], offset, window=win)
               for r in sheet if r.get("geneious_export")}

    rows, details = [], {}
    ref_name = read_fasta(reference)[0][0]
    # one null per (control set, read preparation); only qcal needs the slow
    # leave-one-control-out calibration
    null_cache: dict = {}
    for mode in modes:
        p = _mode_params(base, mode)
        b = batch.prepare_batch(samples, reference, contaminants, p,
                                depth=depth, controls=controls)
        for r in sheet:
            name = r["sample"]
            null_rates, null_note = None, "none"
            if p.error_model != "geneious":
                peers = [c for c in controls if c != name and runs.get(c) == runs.get(name)]
                if peers:
                    calibrate = p.error_model == "qcal"
                    key = (tuple(sorted(peers)), p.phase_mode, p.strand_prior_trials, calibrate)
                    if key not in null_cache:
                        null_cache[key] = null_model.fit_aligned(
                            [b["reads"][c] for c in peers], ref_seq, p, lo, hi,
                            calibrate=calibrate, ref_name=ref_name)[0]
                    null_rates = null_cache[key]
                    null_note = "+".join(sorted(peers))
            print(f"[{mode}] {name}", file=sys.stderr, flush=True)
            res = pipeline.analyse(b["reads"][name], reference, p, b["stats"][name],
                                   null_rates=null_rates, depth=b["depth"])
            s = res["summary"]
            ours = {(v.pos, v.ref, v.alt) for v in res["variants"]}
            row = {"sample": name, "mode": mode, "null_control": r["null"],
                   "depth": s["reads"]["aligned_to_target"],
                   "candidates": s["calling"]["candidates_tested"],
                   "calls": s["calling"]["unique_variants"],
                   "geneious_calls": r["geneious_calls"],
                   "density_x1e5": s["unique_variant_density"],
                   "fraction_GC": s["gc_bias_unique"]["fraction_GC"],
                   "empirical_null_from": null_note,
                   "reads_dropped_by_phasing": s["reads"].get("reads_dropped_minor", 0),
                   "blocks_flagged": sum(1 for x in s["reads"].get("blocks", [])
                                         if x["action"] == "flag")}
            if r["geneious_calls"] is not None:
                row["delta_vs_geneious"] = row["calls"] - r["geneious_calls"]
            if name in exports:
                row.update(site_concordance(ours, exports[name]))
            rows.append(row)
            details[f"{mode}/{name}"] = res
    return {"rows": rows, "details": details, "modes": list(modes),
            "depth": depth, "reference_frame_offset": offset}


def separation(rows: list[dict], mode: str) -> dict:
    sig = [r["density_x1e5"] for r in rows if r["mode"] == mode and not r["null_control"]]
    bg = [r["density_x1e5"] for r in rows if r["mode"] == mode and r["null_control"]]
    if not sig or not bg:
        return {}
    ms, mb = sum(sig) / len(sig), sum(bg) / len(bg)
    return {"signal_mean": round(ms, 2), "background_mean": round(mb, 2),
            "fold": round(ms / mb, 1) if mb else None}


def report_markdown(res: dict) -> str:
    rows = res["rows"]
    has_sites = any("jaccard" in r for r in rows)
    out = ["# nanoshm versus Geneious", "",
           f"Depth: every sample analysed at {rows[0]['depth'] if rows else '?'} reads. "
           "Calls are distinct position-substitution pairs.", ""]
    head = ["Sample", "Mode", "Calls", "Geneious", "Δ", "Density ×10⁻⁵", "Null from",
            "Phasing dropped"]
    if has_sites:
        head += ["Shared", "Precision", "Recall", "Jaccard"]
    out += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows:
        cells = [r["sample"] + (" (control)" if r["null_control"] else ""), r["mode"],
                 r["calls"], r["geneious_calls"] if r["geneious_calls"] is not None else "",
                 r.get("delta_vs_geneious", ""), r["density_x1e5"], r["empirical_null_from"],
                 r["reads_dropped_by_phasing"]]
        if has_sites:
            cells += [r.get("shared", ""), r.get("precision_vs_geneious", ""),
                      r.get("recall_vs_geneious", ""), r.get("jaccard", "")]
        out.append("| " + " | ".join(str(c) for c in cells) + " |")
    out += ["", "## Signal versus background", ""]
    for mode in res["modes"]:
        s = separation(rows, mode)
        if s:
            out.append(f"- **{mode}**: samples {s['signal_mean']} vs controls "
                       f"{s['background_mean']} ×10⁻⁵ ({s['fold']}-fold)")
    gene = [r for r in rows if r["mode"] == res["modes"][0] and r["geneious_calls"] is not None]
    if gene:
        gs = [r["geneious_calls"] for r in gene if not r["null_control"]]
        gb = [r["geneious_calls"] for r in gene if r["null_control"]]
        if gs and gb and sum(gb):
            out.append(f"- **Geneious (calls)**: samples {sum(gs) / len(gs):.1f} vs controls "
                       f"{sum(gb) / len(gb):.1f} ({(sum(gs) / len(gs)) / (sum(gb) / len(gb)):.0f}-fold)")
    if not has_sites:
        out += ["", "Site-level concordance needs a `geneious_export` column in the sample sheet."]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description="Benchmark nanoshm against Geneious")
    ap.add_argument("--sheet", required=True, help="sample sheet CSV/TSV (see module docs)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--contaminant", action="append", default=[])
    ap.add_argument("--preset", required=True)
    ap.add_argument("--modes", default=",".join(MODES),
                    help="comma-separated: recommended,geneious")
    ap.add_argument("--depth", type=int, help="analysis depth (default: common depth)")
    ap.add_argument("--geneious-offset", type=int, default=0,
                    help="add this to Geneious positions to reach the nanoshm reference frame")
    ap.add_argument("--out", default="benchmark")
    args = ap.parse_args()

    sheet = read_sheet(args.sheet)
    p = load_preset(args.preset)
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    res = run_benchmark(sheet, args.reference, args.contaminant, p, modes,
                        args.depth, args.geneious_offset)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "benchmark.csv"), "w", newline="") as fh:
        keys = []
        for r in res["rows"]:
            keys += [k for k in r if k not in keys]
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(res["rows"])
    for key, detail in res["details"].items():
        d = os.path.join(args.out, key)
        pipeline.write_outputs(detail, d)
    with open(os.path.join(args.out, "benchmark.json"), "w") as fh:
        json.dump({k: v for k, v in res.items() if k != "details"}, fh, indent=2)
    md = report_markdown(res)
    with open(os.path.join(args.out, "report.md"), "w") as fh:
        fh.write(md)
    print(md)
    print(f"\nWrote {args.out}/benchmark.csv, report.md and per-sample outputs")


if __name__ == "__main__":
    main()
