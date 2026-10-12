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

from . import (
    batch,
    browse,
    clones,
    indels,
    metrics,
    motifs,
    null_model,
    pipeline,
    qc,
    singletons,
    targeting,
)
from .config import load_preset
from .io_utils import read_fasta
from .reads_filter import downsample

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
    # calls packed into a few bases are a localised artefact, not SHM: real
    # co-occurring blocks share reads and are removed by phasing instead
    clusters = metrics.positional_clusters(vs, bp)
    cluster_notes = [
        f"{c['calls']} calls fall in {c['end'] - c['start'] + 1} bp at "
        f"{c['start']}-{c['end']} (p = {c['p_adjusted']:.1e}); check the spectrum "
        "there before counting them as mutations" for c in clusters]
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
        # C>T and G>A as a share of all calls: what AID gives when the uracil
        # is left unrepaired, and a quick check that calls look like AID's
        "pct_CT_GA": round(100 * sum(1 for v in vs if (v.ref, v.alt) in (("C", "T"), ("G", "A"))) / len(vs), 1) if vs else None,
        "pct_GC_in_WRC_GYW": round(100 * hot / gc, 1) if gc else None,
        "pct_AT_in_WA_TW": round(100 * wa / at, 1) if at else None,
        "analysed_bp": bp,
        "model_used": s["calling"].get("error_model", ""),
        "clustered_calls": sum(c["calls"] for c in clusters),
        "warnings": "; ".join(list(s.get("warnings") or [])
                              + list(s["calling"].get("notes") or [])
                              + cluster_notes),
    }


def corrected_spectrum(counts, base_sites: dict, n_files: int) -> dict:
    """Unique mutations of each class per site of its from-base, per file.

    The denominator SHMTool uses (L = sites x files). A raw count of C>T and
    one of A>G are not comparable when the window holds twice as many C as A;
    divided by the sites that could have mutated, they are.
    """
    out = {}
    for r in "ACGT":
        for a in "ACGT":
            if a == r:
                continue
            L = base_sites.get(r, 0) * max(1, n_files)
            out[f"{r}>{a}"] = round(counts.get(f"{r}>{a}", 0) / L, 6) if L else None
    return out


def sanger_germline(reads_by_sample: dict, ref: str, lo: int, hi: int, p) -> set:
    """Germline differences from the reference, measured over a Sanger project's clones.

    A sample of fewer than `sanger_germline_min_clones` clones cannot tell a
    germline difference (carried by every clone of the mouse) from a mutation,
    so its own germline rule is off. The other mice of the same strain can: a
    change carried by at least `germline_freq_threshold` of all the project's
    clones covering the position, among at least that many clones, is germline.
    Without this, one clone from an Aicda-/- mouse reported the strain's five
    germline differences in J꜀4 as five mutations.
    """
    from collections import Counter
    depth, alts = Counter(), Counter()
    for aligned in reads_by_sample.values():
        for a in aligned:
            for pos, base, _q in a.pairs:
                if lo <= pos <= hi and base in "ACGT":
                    depth[pos] += 1
                    if base != ref[pos].upper():
                        alts[(pos, base)] += 1
    return {(pos, b) for (pos, b), k in alts.items()
            if depth[pos] >= p.sanger_germline_min_clones
            and k / depth[pos] >= p.germline_freq_threshold}


def coverage_dips(depth: list, start: int, frac: float = 0.5, min_run: int = 10) -> str:
    """A warning naming stretches covered by under `frac` of the median depth.

    Density divides by analysed positions times mean depth, so it assumes even
    coverage. A stretch well below the median is where that assumption fails.
    """
    if not depth:
        return ""
    med = sorted(depth)[len(depth) // 2]
    if med <= 0:
        return ""
    runs, run_start = [], None
    for i, d in enumerate(list(depth) + [med]):
        low = d < frac * med
        if low and run_start is None:
            run_start = i
        elif not low and run_start is not None:
            if i - run_start >= min_run:
                runs.append((start + run_start, start + i - 1))
            run_start = None
    if not runs:
        return ""
    where = ", ".join(f"{a}-{b}" for a, b in runs[:4])
    return (f"coverage falls below half the median depth ({med}) at {where}; "
            "densities assume even coverage there")


def add_indels(m: dict, called: list, res: dict, scale: float) -> None:
    """Indel counts for a sample row, kept apart from the substitution density."""
    s = res["summary"]["calling"]
    denom = (s.get("analysed_positions") or 0) * max(1.0, s.get("mean_depth") or 0.0)
    m.update({
        "indels_called": len(called),
        "insertions": sum(1 for v in called if v.kind == "ins"),
        "deletions": sum(1 for v in called if v.kind == "del"),
        "frameshift_indels": sum(1 for v in called if v.frameshift),
        "indel_density_x1e5": round(len(called) * scale / denom, 3) if denom else None,
    })


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

    # project["browse"]: keep the alignment of every read, kept or dropped, so
    # the result can be inspected read by read. Off by default: it re-runs the
    # read preparation to record what became of each one.
    want_browse = bool(project.get("browse"))

    progress(0.05, "aligning and filtering reads")
    fates: dict | None = {} if sanger else None
    b = batch.prepare_batch({s["name"]: s["path"] for s in samples}, reference,
                            contaminants, p, depth=project.get("depth"), controls=controls,
                            fates=fates)
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
    project_germline = None           # Sanger: germline seen across all clones
    # single-read mutations need controls and one shared depth
    single_q = p.single_read_min_quality if (not sanger and depth is not None) else None
    tallies: dict = {}
    single_read: dict = {}
    control_background: dict = {}
    alignments: dict = {}
    coverage: dict = {}
    clone_results: dict = {}
    # project["saturation"]: True for a sensible ladder up to the analysed
    # depth, or an explicit list of depths. Off by default.
    sat = project.get("saturation")
    if sat is True and depth:
        saturation_depths = sorted({max(10, int(depth * f))
                                    for f in (0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0)})
    elif isinstance(sat, (list, tuple)):
        saturation_depths = sorted(int(d) for d in sat)
    else:
        saturation_depths = []
    # indels: tested against the peers' indel rate, or thresholded for Sanger
    want_indels = not p.exclude_indels
    indel_piles: dict = {}
    indel_nulls: dict = {}
    indel_rows: list = []
    indel_qscale = None
    indel_minq = (p.min_base_quality or 20) if sanger else p.indel_min_quality
    if want_indels and not sanger and len(controls) >= 1:
        for c in sorted(controls):
            indel_piles[c] = indels.Pileup(b["reads"][c], ref, lo, hi, indel_minq)
        indel_qscale = indels.calibrate([indel_piles[c] for c in sorted(controls)], ref, p)[0]
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
        few = sanger and len(b["reads"][name]) < p.sanger_germline_min_clones
        if few and project_germline is None:
            project_germline = sanger_germline(b["reads"], ref, lo, hi, p)
        res = pipeline.analyse(b["reads"][name], reference, p, b["stats"][name],
                               null_rates=null, depth=depth,
                               germline=project_germline if few else None)
        if want_browse:
            progress(0.15 + 0.75 * i / len(samples), f"recording the alignment of {name}")
            fates: dict = {}
            kept, _st = pipeline.prepare_reads(s["path"], reference, contaminants, p,
                                               known_blocks=b.get("known_blocks"), fates=fates)
            analysed = {a.qname for a in res["_aligned"]}
            for a in kept:
                if a.qname not in analysed:
                    fates.setdefault(a.qname, {"fate": "not sampled"})
            alignments[name] = browse.alignment(res["_aligned"], fates, reference,
                                                res["variants"])
        if saturation_depths:
            # the whole depth curve, not just the one analysed point. Opt-in:
            # it re-calls variants at every depth.
            res["rarefaction"] = pipeline.rarefaction(
                b["reads"][name], reference, p,
                [d for d in saturation_depths if d <= len(b["reads"][name])],
                draws=project.get("saturation_draws", 3), null_rates=null)
        if null is not None and not control_background:
            # the fitted background a non-control sample is judged against:
            # the highest per-alt error rate at each position (viz_project)
            per_pos: dict = {}
            for (pos0, _alt), rate in null.items():
                if lo <= pos0 <= hi:
                    per_pos[pos0 + 1] = max(per_pos.get(pos0 + 1, 0.0), rate)
            control_background.update(per_pos)
        m = sample_metrics(name, res, ref, lo, hi, s.get("group") or name, name in controls,
                           scale=p.density_scale)
        cl = None if sanger else clones.call_clones(res["_aligned"], res["variants"], p, lo, hi)
        if cl:
            clone_results[name] = cl
            m.update({"reads_with_junction": cl["reads_with_junction"],
                      "clones": cl["clones"],
                      "largest_clone_pct": cl["largest_clone_pct"],
                      "clone_weighted_frequency_x1e3": cl["clone_weighted_frequency_x1e3"]})
            if cl["largest_clone_pct"] >= 50:
                note = (f"one clone holds {cl['largest_clone_pct']:.0f}% of the reads with a "
                        "junction; read-weighted numbers are dominated by it, so compare "
                        "clone_weighted_frequency_x1e3")
                m["warnings"] = "; ".join(x for x in [m["warnings"], note] if x)
        dep = res["summary"]["calling"].get("_depth") or []
        if dep:
            coverage[name] = {"start": lo + 1, "depth": dep}
            low = coverage_dips(dep, lo + 1)
            if low:
                m["warnings"] = "; ".join(x for x in [m["warnings"], low] if x)
        if single_q and peers:
            # every library at the analysed depth, drawn as analyse() draws it: the
            # estimate is divided by the analysed reads' bases, so counting events in
            # a deeper library would inflate it by the ratio of the two
            pd = replace(p, downsample_to=depth)
            for n in (name, *peers):
                if n not in tallies:
                    tallies[n] = singletons.tally(downsample(b["reads"][n], pd)[0], ref, lo, hi)
            single_read[name] = singletons.estimate(
                tallies[name], [tallies[c] for c in peers], single_q,
                p.germline_freq_threshold, ref,
                called={(v.pos - 1, v.alt) for v in res["variants"]})
            add_single_read(m, single_read[name], res, p.density_scale)
        if want_indels and (sanger or peers):
            pile = indels.Pileup(res["_aligned"], ref, lo, hi, indel_minq)
            null_i = None
            if not sanger:
                if peers not in indel_nulls:
                    indel_nulls[peers] = indels.IndelNull([indel_piles[c] for c in peers], ref)
                    indel_nulls[peers].qscale = indel_qscale
                null_i = indel_nulls[peers]
            called_i, _st = indels.call(pile, null_i, p, sanger=sanger)
            add_indels(m, called_i, res, p.density_scale)
            for v in called_i:
                indel_rows.append({"sample": name, "group": s.get("group") or name,
                                   "pos": v.pos, "kind": v.kind, "length": v.length,
                                   "bases": v.seq, "reads": v.supporting_reads,
                                   "depth": v.depth, "frequency": round(v.frequency, 6),
                                   "p_value": _clean(v.p_value), "plus_reads": v.plus_reads,
                                   "minus_reads": v.minus_reads, "homopolymer": v.homopolymer,
                                   "frameshift": v.frameshift})
        m["null_from"] = "+".join(peers) if null is not None else "none"
        per_sample.append(m)
        details[name] = res

    progress(0.92, "summarising")
    notes = []
    qc_out = None
    if sanger:
        # bookkeeping checks: misfiled or duplicated clones, contaminated controls (qc.py)
        sheet = project.get("sample_sheet") or {}
        qc_out = qc.run(details, samples, ref, lo, hi, minq=p.min_base_quality or 20,
                        provenance=sheet.get("provenance"), fates=fates,
                        sheet_excluded=sheet.get("excluded"), aliases=sheet.get("aliases"))
        notes += sheet.get("warnings") or []
        n_prob = len({f["clone"] for f in qc_out["flags"] if f["severity"] == "problem"})
        if n_prob:
            notes.append(f"{n_prob} clone(s) flagged as likely misfiled, duplicated or "
                         "contaminating; see the QC table. Nothing was changed: decide each one "
                         "and record the decision in the sample sheet.")
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
    if want_indels and not sanger and len(controls) < 2:
        notes.append("Indels need a control other than the sample itself, so they were "
                     + ("not called." if not controls else
                        "not called for the control. Two or more controls per run is better."))
    if depth is not None and any(s["reads"] < depth for s in per_sample):
        notes.append("Some samples have fewer reads than the common depth, and mutation "
                     "density rises with depth, so compare those with care.")

    keys = ["unique_density_x1e5", "combined_density_x1e5", "single_read_density_x1e5",
            "clone_weighted_frequency_x1e3",
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

    # AID targeting against a 5-mer model, when one is given (targeting.py)
    targeting_out = None
    if project.get("targeting_model"):
        tm = targeting.load(project["targeting_model"])
        prof = targeting.window_profile(tm, ref, lo, hi)
        per_s, per_g = {}, {}
        for name, res in details.items():
            per_s[name] = targeting.targeting_index([v.pos - 1 for v in res["variants"]], prof)
        for g in set(gmap.values()):
            pos = [v.pos - 1 for n, res in details.items() if gmap[n] == g
                   for v in res["variants"]]
            per_g[g] = targeting.targeting_index(pos, prof)
        for row in per_sample:
            t = per_s.get(row["sample"]) or {}
            row["s5f_index"], row["s5f_p"] = t.get("index"), t.get("p")
        notes_by = {(v["sample"], v["pos"], v["alt"]): None for v in variants_rows}
        for name, res in details.items():
            for a in targeting.annotate(res["variants"], tm, ref, prof):
                notes_by[(name, a["pos"], a["alt"])] = a
        for v in variants_rows:
            a = notes_by.get((v["sample"], v["pos"], v["alt"])) or {}
            v["s5f_mutability_rel"] = a.get("s5f_mutability_rel")
            v["s5f_substitution"] = a.get("s5f_substitution")
        targeting_out = {
            "model": tm.name, "window_mean": prof["mean"],
            "profile": {str(i + 1): round(m / prof["mean"], 4)
                        for i, m in sorted(prof["mutability"].items())} if prof["mean"] else {},
            "samples": per_s, "groups": per_g,
            "expected_spectrum": targeting.expected_spectrum(tm, ref, prof),
            # the model itself, for the hedgehog plot (1,024 values)
            "model_mutability": {k: round(v, 8) for k, v in tm.mutability.items()}}

    base_sites = {b: sum(1 for i in range(lo, hi + 1) if ref[i].upper() == b) for b in "ACGT"}
    n_in_group = Counter(gmap[n] for n in details)
    spec_corr_group = {g: corrected_spectrum(c, base_sites, n_in_group[g])
                       for g, c in spec_group.items()}
    spec_corr_sample = {n: corrected_spectrum(Counter(c), base_sites, 1)
                        for n, c in spec_sample.items()}

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
        # Sanger: flags for a person to review, and every clone not analysed (qc.py)
        "qc": qc_out,
        # per-sample detail of the single-read estimate (singletons.estimate)
        "single_read": single_read,
        "indels": indel_rows,
        # position -> fitted control error rate, for the background diagnostic
        "control_background": {str(k): round(v, 8) for k, v in sorted(control_background.items())},
        "groups": groups,
        "notes": notes,
        "spectrum_group": {g: dict(c) for g, c in spec_group.items()},
        # unique mutations per site of the from-base, per animal: the SHMTool
        # convention, so a base the window is short of is not under-counted
        "spectrum_corrected_group": spec_corr_group,
        "spectrum_corrected_sample": spec_corr_sample,
        "base_composition": base_sites,
        "targeting": targeting_out,
        "coverage": coverage,
        # per sample: clones from the V(D)J junction, each counted once (clones.py)
        "clones": clone_results,
        "spectrum_sample": spec_sample,
        "positional": pos_rows,
        "positional_sample": pos_sample_rows,
        "variants": variants_rows,
        "mutation_load": load_hist,
        "reference": ref,
        "hotspots": hotspots,
        "group_order": list(dict.fromkeys(gmap[s["name"]] for s in samples)),
        "_details": details,
        "_alignments": alignments,
    }


def public(res: dict) -> dict:
    """The result without the bulk per-sample detail, plus which figures it supports."""
    from . import viz_diagnostics
    out = {k: v for k, v in res.items() if not k.startswith("_")}
    base = ["density_by_group", "mutation_frequency_by_group", "mutated_reads_by_group",
            "substitution_matrix", "spectrum_by_group", "positional_map", "mutation_load"]
    if any(s.get("pct_GC_in_WRC_GYW") is not None for s in res.get("samples", [])):
        base.append("hotspot_targeting")
    return {**out, "available_figures": base + viz_diagnostics.drawable(res),
            # the alignments themselves stay on the server; the Reads tab asks
            # for the window it is showing
            "has_alignments": bool(res.get("_alignments"))}


__all__ = ["run_project", "public", "build_params", "metrics", "replace"]
