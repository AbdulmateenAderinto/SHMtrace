"""Run a whole project: many samples, groups, controls, one reference.

This is the engine behind the desktop app. A project is

    {"name": "...", "region": "<preset name or .yaml path>",
     "reference": "<fasta>", "contaminants": [...],
     "platform": "ont" | "pacbio" | "illumina" | "sanger",
     "mode": "recommended" | "geneious" | ...,
     "depth": null | int,
     "samples": [{"name": "WT1", "path": "...", "group": "WT", "control": false}, ...]}

Every sample is prepared together (shared phasing knowledge, one common depth),
each sample is called against a null fitted from the project's control samples
(leaving the sample itself out when it is a control), and per-sample and
per-group tables are computed. Everything returned is JSON-serialisable.
"""
from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from dataclasses import replace

from . import batch, metrics, motifs, null_model, pipeline, singletons
from .config import load_preset
from .io_utils import read_fasta

BASES = "ACGT"
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def build_params(project: dict):
    """Preset, then platform, then mode, then the user's own overrides.

    The platform has to reach `mode_overrides`: on Sanger traces "recommended"
    means the Sanger caller, and without the platform it would silently become
    the quality-calibrated test, which needs more than one read per clone.
    """
    plat = project.get("platform") or "ont"
    mode = project.get("mode") or "recommended"
    if hasattr(pipeline, "build_params"):
        over = {k: (tuple(v) if k == "analysis_window" else v)
                for k, v in (project.get("overrides") or {}).items() if v not in (None, "")}
        return pipeline.build_params(project["region"], platform=plat, mode=mode, **over)
    over = dict(pipeline.mode_overrides(mode))
    for k, v in (project.get("overrides") or {}).items():
        if v not in (None, ""):
            over[k] = tuple(v) if k == "analysis_window" else v
    return load_preset(project["region"], **over)


def _welch(a: list[float], b: list[float]) -> float | None:
    if len(a) < 2 or len(b) < 2:
        return None
    from scipy import stats
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def _welch_log(a: list[float], b: list[float]) -> float | None:
    """Welch on log10 densities, which is what SHM densities are distributed on.

    Mutation densities span orders of magnitude and are roughly log-normal, so a
    test on the raw numbers is dominated by the one most-mutated animal. A sample
    with no mutations has no logarithm; it is left-censored to a third of the
    smallest observed value, which pulls the difference *towards* no difference.
    """
    vals = [v for v in a + b if v > 0]
    if not vals or len(a) < 2 or len(b) < 2:
        return None
    floor = min(vals) / 3
    return _welch([math.log10(max(v, floor)) for v in a],
                  [math.log10(max(v, floor)) for v in b])


def _mwu(a: list[float], b: list[float]) -> float | None:
    if not a or not b:
        return None
    from scipy import stats
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def _clean(x):
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    return x


def sample_metrics(name: str, res: dict, ref: str, lo: int, hi: int, group: str,
                   control: bool, scale: float = 1e5) -> dict:
    s = res["summary"]
    denom = (s["calling"].get("analysed_positions") or 0) * max(
        1.0, s["calling"].get("mean_depth") or 0.0)
    one_event = round(scale / denom, 4) if denom else None
    vs = res["variants"]
    per_read = res.get("per_read_mutations") or []
    nreads = len(per_read)
    bp = hi - lo + 1
    total_mut = sum(per_read)
    gc = sum(1 for v in vs if v.ref in "GC")
    ts = sum(1 for v in vs if (v.ref, v.alt) in TRANSITIONS)
    hot = sum(1 for v in vs if v.ref in "GC" and (
        motifs.is_site("WRC", ref, v.pos - 1) or motifs.is_site("GYW", ref, v.pos - 1)))
    wa = sum(1 for v in vs if v.ref in "AT" and (
        motifs.is_site("WA", ref, v.pos - 1) or motifs.is_site("TW", ref, v.pos - 1)))
    at = len(vs) - gc
    return {
        "sample": name, "group": group, "control": control,
        "reads": nreads,
        "unique_mutations": len(vs),
        "unique_density_x1e5": s["unique_variant_density"],
        # the density one single mutation would give in this sample: the scale a
        # sample with no calls is censored to on a log plot, as the paper does
        "one_event_density_x1e5": one_event,
        "mutation_frequency_x1e3": round(1e3 * total_mut / (nreads * bp), 4) if nreads else None,
        "mutated_reads_pct": round(100 * sum(1 for m in per_read if m) / nreads, 2) if nreads else None,
        "mean_mutations_per_read": round(total_mut / nreads, 4) if nreads else None,
        "pct_GC": round(100 * gc / len(vs), 1) if vs else None,
        "pct_transitions": round(100 * ts / len(vs), 1) if vs else None,
        "pct_GC_in_WRC_GYW": round(100 * hot / gc, 1) if gc else None,
        "pct_AT_in_WA_TW": round(100 * wa / at, 1) if at else None,
        "analysed_bp": bp,
        "model_used": s["calling"].get("error_model", ""),
        "warnings": "; ".join(list(s.get("warnings") or [])
                              + list(s["calling"].get("notes") or [])),
    }


def add_single_read(m: dict, est: dict | None, res: dict, scale: float) -> None:
    """Add the single-read estimate, and the density it completes, to a sample row.

    combined density = (every per-site call + estimated single-read mutations
    the test did not call) per base sequenced. The test is reliable for
    mutations in several reads and cannot see most of those in one; the
    estimate leaves out what the test called, so nothing is counted twice.
    """
    if not est:
        return
    s = res["summary"]["calling"]
    denom = (s.get("analysed_positions") or 0) * max(1.0, s.get("mean_depth") or 0.0)
    called = len(res["variants"])
    multi = sum(1 for v in res["variants"] if v.supporting_reads >= 2)
    f = scale / denom if denom else 0.0
    m.update({
        "multi_read_mutations": multi,
        "single_read_observed": est["observed"],
        "single_read_background": est["expected"],
        "single_read_estimate": est["estimate"],
        "single_read_density_x1e5": round(est["estimate"] * f, 3),
        "single_read_density_ci_low": round(est["ci_low"] * f, 3),
        "single_read_density_ci_high": round(est["ci_high"] * f, 3),
        # a count cannot be negative: below zero the estimate says "none seen"
        "combined_density_x1e5": round((called + max(0.0, est["estimate"])) * f, 3),
        "combined_density_ci_low": round((called + max(0.0, est["ci_low"])) * f, 3),
        "combined_density_ci_high": round((called + max(0.0, est["ci_high"])) * f, 3),
    })


def _slug(s: str) -> str:
    """A group name safe to embed in a column header."""
    return "".join(c if (c.isalnum() or c in "-") else "_" for c in str(s)).strip("_") or "group"


CENSORED = ("unique_density_x1e5", "combined_density_x1e5", "single_read_density_x1e5")


def _censored(rows: list[dict], k: str) -> list[float]:
    """Metric values with a zero density raised to that sample's one-event scale.

    A sample that found nothing has no logarithm. The paper's convention is to
    give it the density one single mutation would have produced in that sample,
    which is the smallest difference the assay could have seen.
    """
    out = []
    for r in rows:
        v = r.get(k)
        if v is None:
            continue
        # the single-read estimate is a difference: it can fall below zero, or
        # just above it, and anything under one event is indistinguishable from none
        one = r.get("one_event_density_x1e5")
        if k in CENSORED and one and v < one:
            v = one
        out.append(v)
    return out


def group_table(samples: list[dict], metric_keys: list[str]) -> list[dict]:
    groups = defaultdict(list)
    for r in samples:
        groups[r["group"]].append(r)
    ctrl_groups = [g for g, rs in groups.items() if all(r["control"] for r in rs)]
    ctrl = ctrl_groups[0] if ctrl_groups else None
    out = []
    for g, rs in groups.items():
        row = {"group": g, "n": len(rs), "control": g == ctrl}
        for k in metric_keys:
            vals = [r[k] for r in rs if r.get(k) is not None]
            row[f"{k}_mean"] = round(statistics.mean(vals), 4) if vals else None
            row[f"{k}_sd"] = round(statistics.stdev(vals), 4) if len(vals) > 1 else None
            if ctrl and g != ctrl:
                cv = [r[k] for r in groups[ctrl] if r.get(k) is not None]
                row[f"{k}_p_welch_vs_{_slug(ctrl)}"] = _clean(_welch(vals, cv))
                row[f"{k}_p_welch_log_vs_{_slug(ctrl)}"] = _clean(
                    _welch_log(_censored(rs, k), _censored(groups[ctrl], k)))
                row[f"{k}_p_mwu_vs_{_slug(ctrl)}"] = _clean(_mwu(vals, cv))
        row["compared_with"] = "" if (not ctrl or g == ctrl) else ctrl
        out.append(row)
    return out


def run_project(project: dict, progress=lambda frac, msg: None) -> dict:
    p = build_params(project)
    reference = project["reference"]
    contaminants = project.get("contaminants") or []
    ref_name, ref = read_fasta(reference)[0]
    lo, hi = pipeline._window(p, ref)
    lo, hi = max(lo, p.end_trim_bp), min(hi, len(ref) - 1 - p.end_trim_bp)
    if hi < lo:
        raise ValueError(
            f"The reference {ref_name} is {len(ref)} bp, which leaves nothing to analyse "
            f"after the {project['region']} preset's window and end trims. Pick a preset "
            "that matches this amplicon, or set the analysis window yourself.")
    samples = project["samples"]
    names = [s["name"] for s in samples]
    if len(set(names)) != len(names):
        raise ValueError("sample names must be unique")
    controls = {s["name"] for s in samples if s.get("control")}
    sanger = (project.get("platform") == "sanger")

    progress(0.05, "aligning and filtering reads")
    b = batch.prepare_batch({s["name"]: s["path"] for s in samples}, reference,
                            contaminants, p, depth=project.get("depth"), controls=controls)
    depth = None if sanger else b["depth"]

    # The quality-calibrated model needs one leave-one-control-out calibration.
    # It is measured once, from every control, and reused for each control's
    # own leave-itself-out null, so a project of N samples calibrates once.
    qscale = None
    if p.error_model == "qcal" and not sanger and controls:
        if len(controls) >= 2:
            progress(0.1, "calibrating the background model")
            qscale = null_model.calibrate_qscale(
                [b["reads"][c] for c in sorted(controls)], ref, ref_name, p, lo, hi)[0]
        else:
            # one control cannot be left out of itself, so the inflation factor
            # cannot be measured; the default is the one calibrated on JH4
            qscale = null_model.DEFAULT_QSCALE

    nulls: dict = {}
    per_sample, details = [], {}
    # single-read mutations need controls and one shared depth
    single_q = p.single_read_min_quality if (not sanger and depth is not None) else None
    tallies: dict = {}
    single_read: dict = {}
    for i, s in enumerate(samples):
        name = s["name"]
        progress(0.15 + 0.75 * i / len(samples), f"calling {name}")
        null = None
        peers = tuple(sorted(c for c in controls if c != name))
        if p.error_model != "geneious" and not sanger and peers:
            if peers not in nulls:
                fit = null_model.fit_aligned([b["reads"][c] for c in peers], ref, p,
                                             lo, hi, calibrate=False, ref_name=ref_name)[0]
                if qscale is not None:
                    fit.qscale = qscale
                nulls[peers] = fit
            null = nulls[peers]
        res = pipeline.analyse(b["reads"][name], reference, p, b["stats"][name],
                               null_rates=null, depth=depth)
        m = sample_metrics(name, res, ref, lo, hi, s.get("group") or name, name in controls,
                           scale=p.density_scale)
        if single_q and peers:
            for n in (name, *peers):
                if n not in tallies:
                    tallies[n] = singletons.tally(b["reads"][n], ref, lo, hi)
            single_read[name] = singletons.estimate(
                tallies[name], [tallies[c] for c in peers], single_q,
                p.germline_freq_threshold, ref,
                called={(v.pos - 1, v.alt) for v in res["variants"]})
            add_single_read(m, single_read[name], res, p.density_scale)
        m["null_from"] = "+".join(peers) if null is not None else "none"
        per_sample.append(m)
        details[name] = res

    progress(0.92, "summarising")
    notes = []
    if sanger:
        notes += sorted({n for s in per_sample for n in ([s["warnings"]] if s["warnings"] else [])
                         if "germline" in n})
    elif not controls:
        notes.append("No control samples, so the sequencing background could not be "
                     "measured. Every mutation was judged from base qualities alone, "
                     "which calls more background. Mark your mutation-free samples as "
                     "controls, or sequence one on the same run.")
    elif len(controls) < 2 and p.error_model == "qcal":
        notes.append(f"Only one control, so the background model could not be calibrated "
                     f"by leaving a control out; the default factor "
                     f"({null_model.DEFAULT_QSCALE:g}, measured on JH4) was used instead. "
                     f"Two or more controls per run is better.")
    if depth is not None and any(s["reads"] < depth for s in per_sample):
        notes.append("Some samples have fewer reads than the common depth, and mutation "
                     "density rises with depth, so compare those with care.")

    keys = ["unique_density_x1e5", "combined_density_x1e5", "single_read_density_x1e5",
            "mutation_frequency_x1e3", "mutated_reads_pct",
            "mean_mutations_per_read", "pct_GC", "pct_transitions", "pct_GC_in_WRC_GYW"]
    groups = group_table(per_sample, keys)

    # spectrum per group (unique mutations, ref>alt counts) and per sample
    spec_group = defaultdict(Counter)
    spec_sample = {}
    positional = defaultdict(lambda: defaultdict(int))   # group -> pos -> reads mutated
    pos_by_sample = defaultdict(lambda: defaultdict(int))  # sample -> pos -> reads mutated
    group_reads = Counter()
    sample_reads = Counter()
    variants_rows = []
    load_hist = {}
    gmap = {s["name"]: (s.get("group") or s["name"]) for s in samples}
    for name, res in details.items():
        g = gmap[name]
        c = Counter(f"{v.ref}>{v.alt}" for v in res["variants"])
        spec_sample[name] = dict(c)
        spec_group[g].update(c)
        sample_reads[name] = len(res.get("per_read_mutations") or [])
        group_reads[g] += sample_reads[name]
        for v in res["variants"]:
            positional[g][v.pos] += v.supporting_reads
            pos_by_sample[name][v.pos] += v.supporting_reads
            variants_rows.append({"sample": name, "group": g, "pos": v.pos, "ref": v.ref,
                                  "alt": v.alt, "reads": v.supporting_reads, "depth": v.depth,
                                  "frequency": round(v.frequency, 6),
                                  "p_value": _clean(v.p_value),
                                  "plus_reads": v.plus_reads, "minus_reads": v.minus_reads,
                                  "context": v.context, "homopolymer": v.homopolymer,
                                  "WRC_GYW": motifs.is_site("WRC", ref, v.pos - 1)
                                  or motifs.is_site("GYW", ref, v.pos - 1)})
        load_hist[name] = dict(Counter(res.get("per_read_mutations") or []))

    pos_rows = []
    for g, d in positional.items():
        n = group_reads[g] or 1
        for pos in sorted(d):
            pos_rows.append({"group": g, "pos": pos, "ref": ref[pos - 1],
                             "mutated_reads": d[pos], "frequency_pct": round(100 * d[pos] / n, 4)})

    pos_sample_rows = []
    for name, d in pos_by_sample.items():
        n = sample_reads[name] or 1
        for pos in sorted(d):
            pos_sample_rows.append({"sample": name, "group": gmap[name], "pos": pos,
                                    "mutated_reads": d[pos],
                                    "frequency_pct": round(100 * d[pos] / n, 4)})

    hotspots = [i + 1 for i in range(lo, hi + 1)
                if motifs.is_site("WRC", ref, i) or motifs.is_site("GYW", ref, i)]
    progress(1.0, "done")
    return {
        "name": project.get("name") or "nanoshm project",
        "params": {"region": project["region"], "platform": project.get("platform") or "ont",
                   "mode": project.get("mode") or "recommended", "depth": depth,
                   "error_model": p.error_model, "analysis_window": [lo + 1, hi + 1],
                   "reference_name": ref_name, "reference_length": len(ref)},
        "samples": per_sample,
        # per-sample detail of the single-read estimate (singletons.estimate)
        "single_read": single_read,
        "groups": groups,
        "notes": notes,
        "spectrum_group": {g: dict(c) for g, c in spec_group.items()},
        "spectrum_sample": spec_sample,
        "positional": pos_rows,
        "positional_sample": pos_sample_rows,
        "variants": variants_rows,
        "mutation_load": load_hist,
        "reference": ref,
        "hotspots": hotspots,
        "group_order": list(dict.fromkeys(gmap[s["name"]] for s in samples)),
        "_details": details,
    }


def public(res: dict) -> dict:
    return {k: v for k, v in res.items() if not k.startswith("_")}


__all__ = ["run_project", "public", "build_params", "metrics", "replace"]
