"""Every number the Cell Reports Methods article and the STAR protocol report, computed from the
analysis outputs, so the text never carries a number by hand.

    from crm_numbers import compute; N = compute()      # dict of formatted strings
    python docs/calibration/crm/manuscript/crm_numbers.py   # print them all

Inputs: docs/calibration/figures/* (the manuscript runs), docs/calibration/crm/out/* (benchmarks,
spike-ins, simulations, public data), the Sanger calls, and the simulation outputs of
REPLICATES_AND_POWER and UMI (data-independent).
"""
import csv
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import lab_details  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import indel_spike  # noqa: E402

CRM = Path(__file__).resolve().parents[1]
SANGER_ARCHIVE_TRACES = 447    # .ab1 files in the eight 2024 orders; see SANGER_VS_NANOPORE.md
OUT = CRM / "out"
CAL = CRM.parent
FIG = CAL / "figures"


def _load(p):
    return json.load(open(p))


def geo(x):
    return 10 ** np.mean(np.log10(x))


def f0(x):
    return f"{x:,.0f}"


def f1(x):
    return f"{x:.1f}"


def f2(x):
    return f"{x:.2f}"


def pct(x):
    return f"{100 * x:.0f}"


def series(xs):
    """'221, 344, and 73' (a list in running text, with the serial comma)."""
    xs = [str(x) for x in xs]
    return xs[0] if len(xs) == 1 else f"{xs[0]} and {xs[1]}" if len(xs) == 2 else ", ".join(xs[:-1]) + f", and {xs[-1]}"


def word(n):
    return ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"][n] if n < 10 else f"{n:,}"


def sci(p):
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10^{{−{abs(int(e))}}}"


def pfmt(p):
    return sci(p) if p < 0.001 else (f"{p:.4f}" if p < 0.01 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}")


def spike(key, caller, tag=""):
    rows = [r for r in _load(OUT / f"spikein_{key}{tag}.json")["rows"] if r["caller"] == caller]
    rec = defaultdict(lambda: [0, 0])
    for r in rows:
        for k, (g, n) in r["recall_by_reads"].items():
            rec[int(k)][0] += g
            rec[int(k)][1] += n
    tp, fp = sum(r["true_calls"] for r in rows), sum(r["false_calls"] for r in rows)
    by = {k: v[0] / v[1] for k, v in rec.items() if v[1]}
    per_ctl = {}
    for c in sorted({r["control"] for r in rows}):
        rr = [r for r in rows if r["control"] == c]
        t_, f_ = sum(r["true_calls"] for r in rr), sum(r["false_calls"] for r in rr)
        per_ctl[c] = f_ / (t_ + f_) if t_ + f_ else 0.0
    return {"by": by, "recall": sum(v[0] for v in rec.values()) / sum(v[1] for v in rec.values()),
            "fp": fp / len(rows), "fdr": fp / (tp + fp) if tp + fp else 0.0, "per_ctl": per_ctl,
            "ge3": sum(rec[k][0] for k in rec if k >= 3) / max(1, sum(rec[k][1] for k in rec if k >= 3))}


def bench(key):
    rows = _load(OUT / f"benchmark_real_{key}.json")["rows"]
    out = {}
    for mode in {r["mode"] for r in rows}:
        rs = [r for r in rows if r["mode"] == mode]
        dens = lambda r: r["density_x1e5"] if r["calls"] else 1e5 / r["bases"]   # noqa: E731
        case = [r for r in rs if not r["control"]]
        ctl = [r for r in rs if r["control"]]
        out[mode] = {"fold": geo([dens(r) for r in case]) / geo([dens(r) for r in ctl]),
                     "ctl_calls": sum(r["calls"] for r in ctl),
                     "case_calls": {r["sample"]: r["calls"] for r in case},
                     "ctl_each": {r["sample"]: r["calls"] for r in ctl},
                     "reads": [r["reads"] for r in rs]}
    return out


def compute():
    N = {}
    from validation_summary import numbers as validation_numbers
    N.update(validation_numbers(OUT))
    # ---------------- spike-ins
    s = {k: {c: spike(k, c) for c in ("recommended", "geneious", "lofreq", "lofreq_nobaq", "lofreq_permissive")}
         for k in ("jh4", "n114", "su_lps", "pim1")}
    r = s["jh4"]["recommended"]
    N.update(sp_jh4_recall=pct(r["recall"]), sp_jh4_fp=f1(r["fp"]), sp_jh4_fdr=f"{r['fdr']:.3f}",
             sp_jh4_ge3=pct(r["ge3"]), sp_jh4_1=pct(r["by"][1]))
    g = s["jh4"]["geneious"]
    N.update(sp_jh4_gen_1=pct(g["by"][1]), sp_jh4_gen_fp=f1(g["fp"]), sp_jh4_gen_fdr=f"{g['fdr']:.3f}")
    N.update(sp_jh4_lfnb=pct(s["jh4"]["lofreq_nobaq"]["recall"]), sp_jh4_lfp=pct(s["jh4"]["lofreq_permissive"]["recall"]),
             sp_jh4_lfp_fp=f1(s["jh4"]["lofreq_permissive"]["fp"]))
    r = s["n114"]["recommended"]
    worst = max(r["per_ctl"], key=r["per_ctl"].get)
    others = [v for k, v in r["per_ctl"].items() if k != worst]
    N.update(sp_n114_recall=pct(r["recall"]), sp_n114_1=pct(r["by"][1]), sp_n114_fp=f1(r["fp"]),
             sp_n114_fdr=f"{r['fdr']:.3f}", sp_n114_worst=worst, sp_n114_worst_fdr=f"{r['per_ctl'][worst]:.2f}",
             sp_n114_other_fdr=f"{min(others):.3f}–{max(others):.3f}",
             sp_n114_lfnb=pct(s["n114"]["lofreq_nobaq"]["recall"]),
             sp_n114_gen_recall=pct(s["n114"]["geneious"]["recall"]), sp_n114_gen_fp=f1(s["n114"]["geneious"]["fp"]),
             sp_n114_gen_fdr=f"{s['n114']['geneious']['fdr']:.3f}",
             # The two Benjamini-Hochberg figures come from different places, and the file names invite a
             # mistake: spikein_jh4_bh.json IS the J_H4 BH run, but spikein_n114_bh.json is NOT the N114 one.
             # It is the run that uses the first, noisier C3 library as a control (FDR 0.105, reported below as
             # c31_bh_fdr). N114's own BH result is the recommended_bh caller inside spikein_n114.json.
             sp_n114_bh_fdr=f"{spike('n114', 'recommended_bh')['fdr']:.3f}",
             sp_jh4_bh_fdr=f"{spike('jh4_bh', 'recommended')['fdr']:.3f}")
    for k in ("su_lps", "pim1"):
        r, g = s[k]["recommended"], s[k]["geneious"]
        lf = [s[k][c]["recall"] for c in ("lofreq", "lofreq_nobaq", "lofreq_permissive")]
        N.update({f"sp_{k}_recall": pct(r["recall"]), f"sp_{k}_fdr": f"{r['fdr']:.3f}",
                  f"sp_{k}_gen_recall": pct(g["recall"]), f"sp_{k}_gen_fdr": f"{g['fdr']:.3f}",
                  f"sp_{k}_lf_lo": pct(min(lf)), f"sp_{k}_lf_hi": pct(max(lf)),
                  f"sp_{k}_lfp_fdr": f"{s[k]['lofreq_permissive']['fdr']:.2f}"})
    N["sp_ext_lfp_fdr"] = f"{min(s['su_lps']['lofreq_permissive']['fdr'], s['pim1']['lofreq_permissive']['fdr']):.2f}–" \
                          f"{max(s['su_lps']['lofreq_permissive']['fdr'], s['pim1']['lofreq_permissive']['fdr']):.2f}"
    N["sp_ext_lf_lo"] = pct(min(min(s[k][c]["recall"] for c in ("lofreq", "lofreq_nobaq", "lofreq_permissive"))
                                for k in ("su_lps", "pim1")))
    N["sp_ext_lf_hi"] = pct(max(max(s[k][c]["recall"] for c in ("lofreq", "lofreq_nobaq", "lofreq_permissive"))
                                for k in ("su_lps", "pim1")))
    fdrs = [s[k]["recommended"]["fdr"] for k in ("jh4", "su_lps", "pim1")]
    N["sp_fdr_range_3"] = f"{min(fdrs):.3f}–{max(fdrs):.3f}"
    # ---------------- N114 sensitivity: the first library of empty-vector subclone C3 as a control
    sen = OUT / "sensitivity_c3_1"
    rows = [r for r in _load(sen / "spikein_n114.json")["rows"] if r["caller"] == "recommended"]
    tp, fp_ = sum(r["true_calls"] for r in rows), sum(r["false_calls"] for r in rows)
    held = [r for r in rows if r["control"] == "C3"]
    htp, hfp = sum(r["true_calls"] for r in held), sum(r["false_calls"] for r in held)
    N["c31_bh_fdr"] = f"{spike('n114_bh', 'recommended')['fdr']:.3f}"      # Benjamini-Hochberg, first C3 library
    N["c31_fdr"] = f"{fp_ / (tp + fp_):.3f}"
    N["c31_held_fdr"] = f"{hfp / (htp + hfp):.2f}"
    N["c31_fp"] = f1(fp_ / len(rows))
    bs = _load(sen / "benchmark_real_n114.json")["rows"]
    rs = [r for r in bs if r["mode"] == "recommended"]
    dens = lambda r: r["density_x1e5"] if r["calls"] else 1e5 / r["bases"]   # noqa: E731
    N["c31_sep"] = f0(geo([dens(r) for r in rs if not r["control"]]) / geo([dens(r) for r in rs if r["control"]]))
    N["c31_ctl_calls"] = str(sum(r["calls"] for r in rs if r["control"]))
    # ---------------- lineage simulation
    lin = {(x["scenario"], x["mode"]): x for x in _load(OUT / "lineage_sim_jh4/lineage_sim.json")["summary"]}
    N.update(lin_1_rec=pct(lin[("polyclonal", "recommended")]["recall_by_support"]["1"]),
             lin_1_gen=pct(lin[("polyclonal", "geneious")]["recall_by_support"]["1"]),
             lin_null_rec=f1(lin[("null", "recommended")]["false_calls_per_library"]),
             lin_null_gen=f1(lin[("null", "geneious")]["false_calls_per_library"]))
    pr, pg = lin[("polyclonal", "recommended")], lin[("polyclonal", "geneious")]
    N.update(lin_fdr_rec=f"{pr['fdr']:.3f}", lin_fdr_gen=f"{pg['fdr']:.3f}",
             lin_fp_rec=f1(pr["false_calls_per_library"]), lin_fp_gen=f1(pg["false_calls_per_library"]))
    ari = [x["adjusted_rand"] for x in lin.values() if x["adjusted_rand"] is not None]
    N["lin_ari"] = f"{min(ari):.2f}–{max(ari):.2f}"
    e = lin[("expanded", "recommended")]
    N.update(lin_rw_est=f1(e["est_read_weighted_x1e3"]), lin_rw_true=f1(e["true_read_weighted_x1e3"]),
             lin_cw_est=f1(e["est_clone_weighted_x1e3"]), lin_cw_true=f1(e["true_clone_weighted_x1e3"]),
             lin_cw_pct=pct(e["est_clone_weighted_x1e3"] / e["true_clone_weighted_x1e3"]))
    # ---------------- the high-allele-fraction rule on a mutated dominant clone (out/mask_sim_jh4)
    mf = OUT / "mask_sim_jh4/lineage_sim.json"
    if mf.exists():
        runs = defaultdict(list)
        for r in _load(mf)["runs"]:
            runs[r["scenario"]].append(r["score"] if "score" in r else r)
        err = lambda rs, k: [100 * (x[f"est_{k}_x1e3"] - x[f"true_{k}_x1e3"]) / x[f"true_{k}_x1e3"] for x in rs]  # noqa: E731
        for tag, sc in (("40", "expanded_mutated_40"), ("60", "expanded_mutated_60"), ("80", "expanded_mutated_80")):
            rs = runs.get(sc) or []
            if not rs:
                continue
            lost = [x["planted_at_mask"] for x in rs]
            N[f"mask{tag}_lost"] = f"{min(lost)}–{max(lost)}" if min(lost) != max(lost) else str(lost[0])
            N[f"mask{tag}_kept"] = str(sum(x["planted_at_mask_called"] for x in rs))
            # the worst run of the three, reported as the size of the shortfall; the prose says "below the
            # truth", so a key that is not a shortfall must not be written as one (minor 24)
            for key, metric in (("rw", "read_weighted"), ("cw", "clone_weighted")):
                worst = min(err(rs, metric))
                N[f"mask{tag}_{key}"] = f"{abs(worst):.0f}" if worst < 0 else f"+{worst:.0f}"
        N["mask_runs"] = f0(len(runs.get("expanded_mutated_60") or []))
    # ---------------- quantification
    q = _load(OUT / "quant_accuracy_jh4/runs.json") + _load(OUT / "quant_accuracy_n114/runs.json")

    def ratio(amp, scen, d, e_, t_):
        rs = [x for x in q if x["amplicon"] == amp and x["scenario"] == scen and x["depth"] == d
              and x["mode"] == "recommended" and x.get(t_)]
        return np.mean([x[e_] / x[t_] for x in rs])
    rw = [ratio(a, sc, 500, "est_read_weighted_x1e3", "true_read_weighted_x1e3")
          for a, scs in (("JH4", ("polyclonal", "expanded", "low")), ("N114", ("subclonal", "low"))) for sc in scs]
    N["q_rw_dev"] = f"{100 * max(abs(1 - x) for x in rw):.0f}"
    N["q_ud_500"] = f2(ratio("JH4", "polyclonal", 500, "est_unique_density_x1e5", "true_unique_density_x1e5"))
    N["q_ud_1000"] = f2(ratio("JH4", "polyclonal", 1000, "est_unique_density_x1e5", "true_unique_density_x1e5"))
    N["q_ud_n114_1000"] = f2(ratio("N114", "subclonal", 1000, "est_unique_density_x1e5", "true_unique_density_x1e5"))
    cd = [ratio(a, sc, 500, "est_combined_density_x1e5", "true_unique_density_x1e5")
          for a, scs in (("JH4", ("polyclonal", "expanded", "low")), ("N114", ("subclonal", "low"))) for sc in scs]
    N["q_cd"] = f"{min(cd):.2f}–{max(cd):.2f}"
    rows = _load(OUT / "single_read_spike.json")
    N["sr_rec"] = pct(np.mean([x["delta"] / x["planted"] for x in rows if x["planted"]]))
    # ---------------- real libraries
    b = {k: bench(k) for k in ("jh4", "n114", "su_lps", "su_ch12", "pim1", "pim1_deep")}
    sep = {k: b[k]["recommended"]["fold"] for k in ("jh4", "n114", "su_lps", "su_ch12")}
    N.update({f"sep_{k}": f0(v) for k, v in sep.items()})
    N["sep_lo"], N["sep_hi"] = f0(min(sep.values())), f0(max(sep.values()))
    gen = [b[k]["geneious"]["fold"] for k in sep]
    lfs = [b[k][c]["fold"] for k in sep for c in ("lofreq", "lofreq_nobaq", "lofreq_permissive")]
    N["gen_lo"], N["gen_hi"] = f0(min(gen)), f0(max(gen))
    N["lf_lo"], N["lf_hi"] = f0(min(lfs)) if min(lfs) >= 10 else f"{min(lfs):.0f}", f0(max(lfs))
    ctl_rec = [b[k]["recommended"]["ctl_calls"] for k in ("jh4", "n114", "su_lps", "su_ch12", "pim1")]
    ctl_gen = [b[k]["geneious"]["ctl_calls"] for k in ("jh4", "n114", "su_lps", "su_ch12", "pim1")]
    N["ctl_rec"], N["ctl_gen"] = f"{min(ctl_rec)}–{max(ctl_rec)}", f"{min(ctl_gen)}–{max(ctl_gen)}"
    # the Summary quotes the drop in AID-null calls per assay against the uncorrected rule, so it is
    # taken per assay from the same two counts and reported as a range, not pooled across assays
    red = [100 * (1 - r / g) for r, g in zip(ctl_rec, ctl_gen) if g]
    N["ctl_red_lo"], N["ctl_red_hi"] = f"{min(red):.0f}", f"{max(red):.0f}"
    j = b["jh4"]
    N["jh4_wt_rec"] = series(j["recommended"]["case_calls"][x] for x in ("WT2", "WT3", "WT4"))
    N["jh4_ko_rec"] = series(j["recommended"]["ctl_each"][x] for x in ("KO1", "KO2", "KO3"))
    N["n114_aid_rec"] = series(b["n114"]["recommended"]["case_calls"][x] for x in ("C1", "D3", "D5"))
    N["jh4_wt_gen"] = series(j["geneious"]["case_calls"][x] for x in ("WT2", "WT3", "WT4"))
    N["jh4_ko_gen"] = series(j["geneious"]["ctl_each"][x] for x in ("KO1", "KO2", "KO3"))
    N["jh4_wt_lf"] = series(j["lofreq"]["case_calls"][x] for x in ("WT2", "WT3", "WT4"))
    N["jh4_wt2_lfnb"] = str(j["lofreq_nobaq"]["case_calls"]["WT2"])
    N["jh4_wt2_rec"] = str(j["recommended"]["case_calls"]["WT2"])
    su = list(b["su_lps"]["recommended"]["case_calls"].values()) + list(b["su_ch12"]["recommended"]["case_calls"].values())
    suc = list(b["su_lps"]["recommended"]["ctl_each"].values()) + list(b["su_ch12"]["recommended"]["ctl_each"].values())
    sg = list(b["su_lps"]["geneious"]["case_calls"].values()) + list(b["su_ch12"]["geneious"]["case_calls"].values())
    sgc = list(b["su_lps"]["geneious"]["ctl_each"].values()) + list(b["su_ch12"]["geneious"]["ctl_each"].values())
    N["su_gen_wt_calls"], N["su_gen_ctl_calls"] = f"{min(sg)}–{max(sg)}", f"{min(sgc)}–{max(sgc)}"
    N["su_wt_calls"], N["su_ctl_calls"] = f"{min(su)}–{max(su)}", f"{min(suc)}–{max(suc)}"
    pm = b["pim1"]
    N["pim1_rec_fold"], N["pim1_rec_ctl"] = f1(pm["recommended"]["fold"]), str(pm["recommended"]["ctl_calls"])
    N["pim1_gen_fold"], N["pim1_gen_ctl"] = f1(pm["geneious"]["fold"]), str(pm["geneious"]["ctl_calls"])
    pd_ = b["pim1_deep"]
    N["pim1d_reads"] = f"{min(pd_['recommended']['reads']):,}–{max(pd_['recommended']['reads']):,}"
    rc = list(pd_["recommended"]["case_calls"].values())
    rk = list(pd_["recommended"]["ctl_each"].values())
    gc = list(pd_["geneious"]["case_calls"].values())
    gk = list(pd_["geneious"]["ctl_each"].values())
    N["pim1d_rec"] = f"{min(rc)}–{max(rc)}"
    N["pim1d_rec_ctl"] = "none" if max(rk) == 0 else f"{min(rk)}–{max(rk)}"
    N["pim1d_gen"], N["pim1d_gen_ctl"] = f"{min(gc)}–{max(gc)}", f"{min(gk)}–{max(gk)}"
    # ---------------- calibration (leave-one-control-out)
    for k, run in (("jh4", "jh4_shmtrace"), ("n114", "n114_shmtrace"), ("full", "jh4_shmtrace_full")):
        cal = _load(FIG / run / "project.json")["calibration"]
        tr = cal["qscale_trace"]
        N[f"loco_{k}_face"] = f1(100 * tr["1.0"][0] / tr["1.0"][1])
        N[f"loco_{k}_face_n"] = f"{tr['1.0'][0]} of {tr['1.0'][1]:,}"
        N[f"loco_{k}_margin"] = f"{cal['qscale']:.0f}"
        last = tr[str(float(cal["qscale"]))]
        N[f"loco_{k}_final"] = f1(100 * last[0] / last[1])
    margins = {N["loco_jh4_margin"], N["loco_n114_margin"], N["loco_full_margin"]}
    N["loco_margins"] = (f"{N['loco_jh4_margin']}×" if len(margins) == 1 else
                         f"{N['loco_jh4_margin']}× on J_{{H}}4, {N['loco_full_margin']}× on the full length, and "
                         f"{N['loco_n114_margin']}× on N114 V_{{H}}")
    # ---------------- J_H4 manuscript run, Sanger comparison
    jp = _load(FIG / "jh4_shmtrace/project.json")
    S = {x["sample"]: x for x in jp["samples"]}
    N["m_wt_calls"] = ", ".join(str(S[x]["unique_mutations"]) for x in ("WT2", "WT3", "WT4"))
    wtc = [S[x]["unique_mutations"] for x in ("WT2", "WT3", "WT4")]
    koc = [S[x]["unique_mutations"] for x in ("KO1", "KO2", "KO3")]
    N["m_wt_range"], N["m_ko_range"] = f"{min(wtc)}–{max(wtc)}", f"{min(koc)}–{max(koc)}"
    N["loco_margins_x"] = N["loco_margins"]
    N["m_ko_calls"] = ", ".join(str(S[x]["unique_mutations"]) for x in ("KO1", "KO2", "KO3"))
    st = _load(FIG / "sanger_vs_nanopore/recommended_manuscript/stats.json")
    # The Sanger accounting the Results quote, derived so the three figures cannot drift apart:
    # the archive, what the sample sheet files to these six mice, and what the analysis returned.
    # The 447 .ab1 files of the eight 2024 orders (SANGER_VS_NANOPORE.md) sit outside the
    # repository, so that one is declared and checked only when the directory is to hand. None of
    # this touches the earlier Feng orders, which hold the Fam72a material and are a different
    # campaign: no trace from them is in this archive or this sample sheet.
    N["sanger_clones"] = str(sum(st["density_sanger"]["sequences"].values()))
    with open(CAL / "samplesheets/jh4_manuscript.csv") as fh:
        N["sanger_sheet"] = str(sum(1 for _ in csv.DictReader(fh)))
    archive = Path("~/Documents/WT JH4 SANGER AB1").expanduser()
    if archive.is_dir():
        seen = len(list(archive.rglob("*.ab1")))
        assert seen == SANGER_ARCHIVE_TRACES, (
            f"the Sanger archive holds {seen} chromatograms, not the {SANGER_ARCHIVE_TRACES} the "
            "Results quote; update SANGER_ARCHIVE_TRACES and SANGER_VS_NANOPORE.md together")
    N["sanger_archive"] = str(SANGER_ARCHIVE_TRACES)
    # the Figure S12 accounting: what the sheet submitted after its own exclusions, and the clones
    # each mouse returned. "103 total" is a count of clones, not a power of ten, which is why it is
    # built here rather than written into the legend where a stray 10^{3} could pass unnoticed.
    with open(CAL / "samplesheets/jh4_manuscript.csv") as fh:
        flagged = sum(1 for r in csv.DictReader(fh) if r["exclude"].strip())
    N["sanger_excluded"] = str(flagged)
    N["sanger_submitted"] = str(int(N["sanger_sheet"]) - flagged)
    seq = st["density_sanger"]["sequences"]
    wt, ko = [seq[m] for m in ("WT2", "WT3", "WT4")], [seq[m] for m in ("A1", "A2", "A3")]
    N["sanger_wt_clones"] = ", ".join(str(x) for x in wt[:-1]) + f" and {wt[-1]}"
    N["sanger_ko_clones"] = ", ".join(str(x) for x in ko[:-1]) + f" and {ko[-1]}"
    N["sanger_wt_total"], N["sanger_ko_total"] = str(sum(wt)), str(sum(ko))
    N["sanger_failed"] = str(int(N["sanger_submitted"]) - int(N["sanger_clones"]))
    do, ds = st["density_ont"], st["density_sanger"]
    N.update(ont_fold=f0(do["fold"]), ont_p=pfmt(do["p_welch_log"]), sanger_fold=f0(ds["fold"]),
             sanger_p=pfmt(ds["p_welch_log"]), ont_wt_geo=f0(do["wt_geomean"]), ont_ko_geo=f2(do["ko_geomean"]),
             sanger_wt_geo=f0(ds["wt_geomean"]), sanger_ko_geo=f0(ds["ko_geomean"]))
    sp = st["spectrum"]
    N.update(r_spec=f2(sp["pearson"]["r"]), r_spec_ci=f"{sp['pearson']['ci'][0]:.2f}–{sp['pearson']['ci'][1]:.2f}",
             slope=f2(sp["pearson"]["slope"]), interaction_p=f2(sp["rm_anova_interaction"]["p"]),
             sidak_min=f2(sp["paired_sidak_min"]), n_sanger=f0(sp["sanger_unique"]), n_ont=f0(sp["ont_unique"]))
    pm = _load(CRM / "out/fig5_per_mouse.json")["summary"]["gc_ts_diff"]   # crm/fig5_per_mouse.py
    N["pm_gcts"] = f"{min(-x for x in pm):.0f}–{max(-x for x in pm):.0f}"
    be = st["beyond"]
    N["r_pos"] = f2(be["positional_r"])
    hs, ho = np.mean([m[0] for m in be["hotspot"]["sanger"]]), np.mean([m[0] for m in be["hotspot"]["ont"]])
    as_, ao = np.mean([m[1] for m in be["hotspot"]["sanger"]]), np.mean([m[1] for m in be["hotspot"]["ont"]])
    N.update(hot_gc_s=f0(hs), hot_gc_o=f0(ho), hot_at_s=f0(as_), hot_at_o=f0(ao),
             hot_gc_exp=f0(be["hotspot_expected"][0]), hot_at_exp=f0(be["hotspot_expected"][1]))
    # ---------------- full length
    fp_ = _load(FIG / "jh4_shmtrace_full/project.json")
    FS = {x["sample"]: x for x in fp_["samples"]}
    dens = lambda x: x["unique_density_x1e5"] or x["one_event_density_x1e5"]   # noqa: E731
    wtf = geo([dens(FS[x]) for x in ("WT2", "WT3", "WT4")])
    kof = geo([dens(FS[x]) for x in ("KO1", "KO2", "KO3")])
    N.update(full_wt=f0(wtf), full_ko=f2(kof), full_fold=f0(wtf / kof))
    pos = [v["pos"] for v in fp_["variants"] if v["sample"].startswith("WT")]
    N["full_450"] = pct(np.mean(np.array(pos) <= 90 + 450))
    first = sum(1 for x in pos if x < 90 + 250)
    last = sum(1 for x in pos if x > 1230 - 250)
    N["full_decline"] = f0(first / max(1, last))

    def spec(vs):
        c = defaultdict(int)
        for v in vs:
            if v["sample"].startswith("WT"):
                c[v["ref"] + v["alt"]] += 1
        t = sum(c.values())
        return c, t
    c1, t1 = spec(jp["variants"])
    c2, t2 = spec(fp_["variants"])
    keys = sorted(set(c1) | set(c2))
    x = [100 * c1[k] / t1 for k in keys]
    y = [100 * c2[k] / t2 for k in keys]
    N.update(s1_r=f"{stats.pearsonr(x, y)[0]:.3f}", s1_n529=f0(t1), s1_nfull=f0(t2),
             s1_maxdiff=f1(max(abs(a - b_) for a, b_ in zip(x, y, strict=True))))
    # ---------------- targeting, clones, trees
    t = _load(OUT / "targeting_jh4.json")
    wi = [t["wt"][x]["index"] for x in ("WT2", "WT3", "WT4")]
    wp = [t["wt"][x]["p"] for x in ("WT2", "WT3", "WT4")]
    qu = t["quintiles"]
    N.update(s5f_wt=f"{min(wi):.2f}–{max(wi):.2f}", s5f_wt_p=f"{max(wp):.3f}",
             s5f_ko=f2(t["ko_all_positions"]["index"]), s5f_ko_p=f2(t["ko_all_positions"]["p"]),
             s5f_dose=f1(qu[-1]["wt_calls_per_position"] / qu[0]["wt_calls_per_position"]))
    ex = {x["sample"]: x for x in _load(OUT / "extras_jh4/project.json")["samples"]}
    N.update(wt2_largest=f0(ex["WT2"]["largest_clone_pct"]),
             wt2_rw=f1(ex["WT2"]["mutation_frequency_x1e3"]), wt2_cw=f1(ex["WT2"]["clone_weighted_frequency_x1e3"]),
             wt3_rw=f1(ex["WT3"]["mutation_frequency_x1e3"]), wt3_cw=f1(ex["WT3"]["clone_weighted_frequency_x1e3"]),
             wt4_rw=f1(ex["WT4"]["mutation_frequency_x1e3"]), wt4_cw=f1(ex["WT4"]["clone_weighted_frequency_x1e3"]))
    tr = _load(FIG / "mutation_trees/trees.json")
    w3 = tr["wt3_figure_clone"]
    nodes = w3["nodes"]
    roots = [n for n in nodes if n["parent"] is None]
    root = max(roots, key=lambda n: n["reads"])
    kids = [n for n in nodes if n["parent"] == root["id"]]
    row = next(x for x in tr["clones"]["WT3"] if x["clone"] == w3["clone"])
    fg = tr["four_gamete"]["jh4"][f"WT3:{w3['clone']}"]
    N.update(tree_reads=str(w3["reads"]), tree_fit=pct(w3["fit"]), tree_founders=word(len(root["mutations"])),
             tree_sub=word(len(kids)), tree_depth=str(row["max_depth"]),
             tree_incompat=f"{fg['incompatible']} of {fg['pairs']:,}", tree_incompat_pct=f1(fg["incompatible_pct"]))
    wt2 = tr["clones"]["WT2"][0]
    N["wt2_clone_reads"] = str(wt2["reads"])
    # ---------------- depth series
    for k in ("jh4", "n114"):
        case = ("WT2", "WT3", "WT4") if k == "jh4" else ("C1", "D3", "D5")
        dr = OUT / f"depth_draws_{k}.json"
        if dr.exists():
            # repeated draws: the share each library had reached at 500 reads, averaged over draws, and the
            # spread over draws, so the figure's and the text's numbers come from the same runs. The denominator
            # is the most that library ever yielded over every depth and draw, as in Figure S9E.
            runs = _load(dr)["runs"]
            per = {}
            for r in runs:
                per.setdefault((r["sample"], r["depth"]), []).append(r["calls"])
            shares = []
            for x in case:
                mx = max(c for (sx, _d), v in per.items() if sx == x for c in v)
                shares += [c / mx for c in per[(x, 500)]]
            N[f"depth_{k}_500"] = pct(np.mean(shares))
            N[f"depth_{k}_500_range"] = f"{pct(min(shares))}–{pct(max(shares))}"
        else:
            d = _load(FIG / f"depth_shmtrace/{k}.json")
            mx = {x: max(d[dep][x]["calls"] for dep in d) for x in case}
            N[f"depth_{k}_500"] = pct(np.mean([d["500"][x]["calls"] / mx[x] for x in case]))
            N[f"depth_{k}_500_range"] = "one draw per depth"
    # ---------------- N114
    np_ = _load(FIG / "n114_shmtrace/project.json")
    NS = {x["sample"]: x for x in np_["samples"]}
    N["n114_calls"] = ", ".join(str(NS[x]["unique_mutations"]) for x in ("C1", "D3", "D5"))
    N["n114_ctl"] = ", ".join(str(NS[x]["unique_mutations"]) for x in ("C3", "C4", "D1"))
    # ---------------- public CH12F3 and Ramos details
    ch = _load(OUT / "public/ch12/project.json")
    ev = [x for x in ch["samples"] if x["group"] == "EV"]
    nev = sum(x["unique_mutations"] for x in ev)
    N["ch12_ev"] = (f"{word(nev)} call{'s' if nev != 1 else ''} in {sum(x['reads'] for x in ev):,} empty-vector reads")
    nr = sum(x["unique_mutations"] for x in ch["samples"] if x["group"] == "R50E")
    N["ch12_r50e"] = "made no calls" if nr == 0 else f"made {word(nr)} call{'s' if nr != 1 else ''}"
    mut = sorted(v["pos"] for v in ch["variants"] if v["group"] == "R174ER178E")
    best = max(((lo, sum(lo <= x < lo + 9 for x in mut)) for lo in range(min(mut), max(mut) + 1)), key=lambda z: z[1])
    inw = [x for x in mut if best[0] <= x < best[0] + 9]
    N["ch12_cluster"] = f"{min(inw)}–{max(inw)}"
    N["ch12_cluster_pct"] = pct(len(inw) / len(mut))
    N["ch12_cluster_bp"] = str(max(inw) - min(inw) + 1)
    hm = _load(OUT / "public/hmces/project.json")
    gc = [x["pct_GC"] for x in hm["samples"] if x["group"] != "AID_KO" and x.get("pct_GC") is not None]
    N["hmces_ko3_gc"] = f"{next(x for x in hm['samples'] if x['sample'] == 'AIDKO_3')['pct_GC']:.0f}"
    hmc = [x["unique_mutations"] for x in hm["samples"] if x["group"] in ("WT", "HMCES_KO")]
    N["hm_calls"] = f"{min(hmc)}–{max(hmc)}"
    N["hmces_ko3_calls"] = str(next(x for x in hm["samples"] if x["sample"] == "AIDKO_3")["unique_mutations"])
    k3 = _load(OUT / "public/hmces_k3/project.json")
    mut = [x for x in k3["samples"] if x["group"] in ("WT", "HMCES_KO")]
    calls = [x["unique_mutations"] for x in mut]
    gcs = [x["pct_GC"] for x in mut if x.get("pct_GC") is not None]
    N["hmk3_calls"] = f"{min(calls)}–{max(calls)}"
    N["hmk3_gc"] = f"{min(gcs):.0f}–{max(gcs):.0f}" if gcs else "–"
    N["hmces_gc_range"] = f"{min(gc):.0f}–{max(gc):.0f}" if gc else "[no substitution calls; see note]"
    # sensitive mode (recommended at FDR 0.20), same spike-in design (spikein.py --tag _sensitive)
    for k in ("jh4", "n114"):
        r = spike(k, "recommended_fdr20", "_sensitive")
        N.update({f"sens_{k}_1": pct(r["by"][1]), f"sens_{k}_fdr": f"{r['fdr']:.2f}",
                  f"sens_{k}_fp": f1(r["fp"]), f"sens_{k}_all": pct(r["recall"])})
        assert r["by"][1] > spike(k, "geneious")["by"][1], "the text says sensitive beats the Geneious rule"
    # ClairS (somatic, with a matched normal) and Clair3 (germline) on the same planted mutations
    # (spikein.py --clairs / --clair3), ClairS on clean high-fraction plants (clairs_highaf.py) and on the
    # real wild-type libraries (clairs_real.py)
    for k in ("jh4", "n114"):
        cs, csp = spike(k, "clairs", "_clairs"), spike(k, "clairs_permissive", "_clairs")
        c3, c3p = spike(k, "clair3", "_clair3"), spike(k, "clair3_permissive", "_clair3")
        N.update({f"cs_{k}": pct(max(cs["recall"], csp["recall"])), f"cs_{k}_fp": f1(max(cs["fp"], csp["fp"])),
                  f"c3_{k}": pct(max(c3["recall"], c3p["recall"]))})
        rows = [r for r in _load(OUT / f"spikein_{k}_clairs.json")["rows"] if r["caller"] == "clairs"]
        N[f"cs_{k}_n"] = f"{sum(r['planted'] for r in rows):,}"
    hi = _load(OUT / "clairs_highaf_jh4.json")["rows"]
    band = defaultdict(lambda: [0, 0])
    for r in hi:
        if r["caller"] != "clairs":
            continue
        for kk, (g, n) in r["recall_by_reads"].items():
            b = 5 if int(kk) < 38 else (10 if int(kk) < 76 else 20)
            band[b][0] += g
            band[b][1] += n
    N.update({f"cs_hi_{b}": pct(band[b][0] / band[b][1]) for b in (5, 10, 20)})
    N["cs_hi_n"] = f0(band[20][1])
    if (OUT / "deepsomatic_highaf_jh4.json").exists():      # DeepSomatic on the same plants (clairs_highaf.py)
        band = defaultdict(lambda: [0, 0])
        for r in _load(OUT / "deepsomatic_highaf_jh4.json")["rows"]:
            if r["caller"] == "deepsomatic_permissive":            # as drawn in Figure 3C
                for kk, (g, n) in r["recall_by_reads"].items():
                    band[5 if int(kk) < 38 else (10 if int(kk) < 76 else 20)][0] += g
                    band[5 if int(kk) < 38 else (10 if int(kk) < 76 else 20)][1] += n
        N.update({f"ds_hi_{b}": pct(band[b][0] / band[b][1]) for b in (5, 10, 20)})
    # phrases that read naturally when a caller recalled nothing, over every dataset it was run on
    names = {"jh4": "J_{H}4", "n114": "N114", "su_lps": "Sμ", "pim1": "*Pim1*"}
    for tag, ftag, settings in (("cs", "_clairs", ("clairs", "clairs_permissive")),
                                ("c3", "_clair3", ("clair3", "clair3_permissive")),
                                ("dv", "_deepvariant", ("deepvariant", "deepvariant_permissive")),
                                ("ds", "_deepsomatic", ("deepsomatic", "deepsomatic_permissive"))):
        done = [k for k in names if (OUT / f"spikein_{k}{ftag}.json").exists()]
        if not done:
            continue
        rec = {k: max(spike(k, s, ftag)["recall"] for s in settings) for k in done}
        fp = max(spike(k, s, ftag)["fp"] for k in done for s in settings)
        n_pl = sum(r["planted"] for k in done for r in _load(OUT / f"spikein_{k}{ftag}.json")["rows"]
                   if r["caller"] == settings[0])
        where = ", ".join(names[k] for k in done[:-1]) + f", and {names[done[-1]]}" if len(done) > 2 else \
            " and ".join(names[k] for k in done)
        parts = [f"{pct(rec[k])}% on {names[k]}" for k in done]
        parts[0] = parts[0].replace("% on", "% of planted mutations on", 1)
        N[f"{tag}_recall"] = (f"none of the {n_pl:,} planted mutations on {where}" if not any(rec.values()) else
                              ", ".join(parts[:-1]) + f", and {parts[-1]}")
        # recall at the defaults alone, and for mutations in 20 reads (the most any spike-in had)
        N[f"{tag}_default"] = "none" if not any(spike(k, settings[0], ftag)["recall"] for k in done) else "some"
        r20 = [max(spike(k, s, ftag)["by"].get(20, 0) for s in settings) for k in done]
        N[f"{tag}_20"] = f"{pct(min(r20))}–{pct(max(r20))}"
        N[f"{tag}_fp"] = "no" if fp == 0 else f1(fp)
        N[f"{tag}_sets"] = word(len(done))
        N[f"{tag}_max"] = pct(max(rec.values()))
    # two callers in one sentence: "Clair3 and DeepVariant each recalled none of ..." when both recalled nothing
    for tag, (x, y), (lx, ly) in (("som", ("cs", "ds"), ("ClairS", "DeepSomatic")),
                                  ("germ", ("c3", "dv"), ("Clair3", "DeepVariant"))):
        if f"{x}_recall" not in N or f"{y}_recall" not in N:
            continue
        if N[f"{x}_recall"] == N[f"{y}_recall"] and N[f"{x}_recall"].startswith("none"):
            N[f"{tag}_recall"] = f"{lx} and {ly} each recalled {N[f'{x}_recall']}"
        else:
            N[f"{tag}_recall"] = f"{lx} recalled {N[f'{x}_recall']}, and {ly} {N[f'{y}_recall']}"
    # the somatic callers on the real AID-active libraries: the most passed in total, over their two settings
    for caller, k, tag in (("deepsomatic", "jh4", "ds_real_pass"), ("deepsomatic", "n114", "ds_real_n114_pass"),
                           ("clairs", "n114", "cs_real_n114_pass")):
        if (OUT / f"{caller}_real_{k}.json").exists():
            rr = _load(OUT / f"{caller}_real_{k}.json")["libraries"]
            N[tag] = word(max(sum(rr[w][s]["pass_somatic"] for w in rr) for s in (caller, f"{caller}_permissive")))
    real = _load(OUT / "clairs_real_jh4.json")["libraries"]
    cands = real["WT3"]["clairs_permissive"]["candidates"]
    N.update(cs_real_pass=word(sum(real[w]["clairs_permissive"]["pass_somatic"] for w in real)),
             cs_real_empty=word(sum(1 for w in real if not real[w]["clairs_permissive"]["candidates"])),
             cs_real_cand=f0(len(cands)), cs_real_rej=f0(sum(c["filter"] != "PASS" for c in cands)),
             cs_real_hi=word(sum(c["af"] >= 0.2 for c in cands)),
             cs_real_hi_pass=(lambda k, n: "neither" if (k, n) == (0, 2) else ("none" if k == 0 else word(k)))(
                 sum(c["af"] >= 0.2 and c["filter"] == "PASS" for c in cands), sum(c["af"] >= 0.2 for c in cands)),
             cs_real_reads=f"{min(c['alt_reads'] for c in cands)}–{max(c['alt_reads'] for c in cands)}",
             cs_real_naf=f"{100 * max(c['normal_af'] for c in cands):.1f}")
    # the CH12 AID off-targets OT1 and OT2 (PRJNA1481264, second pass; run_ircm_offtargets.py)
    OFF = Path.home() / "Nanoshm/public_data/PRJNA1481264"
    for k in ("ot1", "ot2"):
        f = OFF / f"offtarget_{k}_2000.json"
        if not f.exists():
            continue
        o = _load(f)
        by = defaultdict(list)
        for s in o["samples"]:
            by[s["group"]].append(s)
        # controls left out of the decoys by the noisy-control rule, named so the text can report them
        drop = [m.group(1) for n in o.get("notes", [])
                if (m := re.search(r"Control (\S+) is far noisier", n))]
        N[f"{k}_ctl_drop"] = series(drop) if drop else ""
        N[f"{k}_wt"] = series(s["unique_mutations"] for s in by["vehicle"])
        N[f"{k}_dko"] = series(s["unique_mutations"] for s in by["dKO"])
        N[f"{k}_ko"] = series(s["unique_mutations"] for s in by["AID_KO"])
        # No fold is quoted for these libraries. Controls with no calls would have to be floored before a ratio
        # could be taken, and the floor, not the data, would then set its size; the groups are compared by their
        # densities and a test on them instead (the reviewer's point; see REVIEW_RESPONSE.md, major issue 1).
        grp = {g["group"]: g for g in o["groups"]}
        for g, tag in (("vehicle", "wt"), ("dKO", "dko"), ("AID_KO", "ko")):
            N[f"{k}_{tag}_dens"] = f2(grp[g]["unique_density_x1e5_mean"])
            N[f"{k}_{tag}_sd"] = f2(grp[g]["unique_density_x1e5_sd"])
        # No test is quoted for the double knockout against the AID-competent libraries. Both groups would
        # have to be floored first, every OT2 dKO library having seen nothing, and the floor then decides the
        # answer: on log frequency floored at 1e-4 OT1 gave p = 0.034, but 0.0085 at 1e-3 and 0.106 at 1e-6,
        # crossing the significance threshold on the constant alone, while OT2's floored dKO group has zero
        # variance and makes the t statistic meaningless. The counts are unambiguous and are what the text
        # reports: 1, 1 and 0 calls at OT1 and 0, 0 and 0 at OT2, against hundreds in the AID-competent
        # libraries. The dKO-against-control comparison below needs no such floor.
        N[f"{k}_p_dko_ko"] = pfmt(grp["dKO"]["unique_density_x1e5_p_welch_log_vs_AID_KO"])
        vs = [v for s in by["vehicle"] for v in o["variants"] if v["sample"] == s["sample"]]
        N[f"{k}_gc"] = pct(sum(v["ref"] in "CG" for v in vs) / len(vs))
        N[f"{k}_n"] = f0(len(vs))
        g = next(x for x in o["groups"] if x["group"] == "vehicle")
        N[f"{k}_p"] = pfmt(g["mutation_frequency_x1e3_p_welch_log_vs_AID_KO"])
        N[f"{k}_p_dens"] = pfmt(g["unique_density_x1e5_p_welch_log_vs_AID_KO"])
    # Figure 2 and Figures S2 and S10: anatomy of a call, the control FDR, where real SHM sits, the background along
    # the amplicon and reported quality against measured error (intuition_data.py, fig2_call.py)
    from bisect import bisect_right
    su = _load(OUT / "intuition_su_lps.json")
    an = su["anatomy"]
    case = an["case"]
    tr, er = an["sites"]["true"], an["sites"]["error"]

    def reads(site, n):
        v = site["libraries"][n]
        return v["alt_plus"] + v["alt_minus"], v["depth_plus"] + v["depth_minus"]

    def med(qs):
        return sorted(qs)[len(qs) // 2]
    k_t, d_t = reads(tr, case)
    k_e, d_e = reads(er, case)
    ko = max((n for n in er["libraries"] if n != case), key=lambda n: reads(er, n)[0])
    k_k, d_k = reads(er, ko)
    ctl_t = [reads(tr, n)[0] for n in tr["libraries"] if n != case]
    hq_t = [q for n in tr["libraries"] if n != case for q in tr["libraries"][n]["alt_quals"] if q >= 20]
    N.update(ana_case=case, ana_t_pos=f0(tr["pos"]), ana_t_ctx=tr["context"], ana_t_reads=f"{k_t} of {d_t:,}",
             ana_t_q=f0(med(tr["libraries"][case]["alt_quals"])), ana_t_ctl=f"{min(ctl_t)}–{max(ctl_t)}",
             ana_t_ctl_hq=word(len(hq_t)), ana_t_p=f"{tr['p']:.0e}".replace("e-", " × 10^{−") + "}",
             ana_e_pos=f0(er["pos"]), ana_e_ctx=er["context"], ana_e_reads=f"{k_e} of {d_e:,}",
             ana_e_q=f0(med(er["libraries"][case]["alt_quals"])), ana_e_ko=ko, ana_e_ko_reads=f"{k_k} of {d_k:,}",
             ana_e_ko_q=f0(med(er["libraries"][ko]["alt_quals"])),
             ana_e_gen_ctl=word(er["geneious_calls"][0]), ana_e_gen_case=word(er["geneious_calls"][1]),
             # the recurrent candidate's p: the test gives exactly 1 here, which reads better as "1"
             # than as pfmt's "1.00", and still moves if the data ever stop saturating the test
             ana_e_p="1" if er["p"] >= 0.995 else pfmt(er["p"]))
    c = su["candidates"][case]
    N.update(ana_thr=f"{c['threshold']:.3f}", ana_kept=f0(bisect_right(sorted(c["p"]), c["threshold"])),
             ana_cand=f"{len(c['p']):,}")
    # Figure S7C: UMI consensus on the simulated libraries, four mutation-bearing
    # conditions plus the two AID-null ones the panels leave out
    umi = _load(CAL / "umi_sim/runs.json")
    cells = [("JH4", "polyclonal"), ("JH4", "low"), ("N114", "subclonal"), ("N114", "low")]
    pick = lambda c, kind, f: [r[kind][f] for r in umi if (r["amplicon"], r["scenario"]) == c]
    mean = lambda v: sum(v) / len(v)
    N["umi_fp_poly"] = f2(mean(pick(cells[0], "umi", "false_calls")))
    rr = [100 * mean(pick(c, "reads", "recall")) for c in cells]
    ru = [100 * mean(pick(c, "umi", "recall")) for c in cells]
    N["umi_rec_reads_lo"], N["umi_rec_reads_hi"] = f0(min(rr)), f0(max(rr))
    N["umi_rec_umi_lo"], N["umi_rec_umi_hi"] = f0(min(ru)), f0(max(ru))
    N["umi_null_jh4"] = f2(mean(pick(("JH4", "null"), "reads", "false_calls")))
    N["umi_null_n114"] = f2(mean(pick(("N114", "null"), "reads", "false_calls")))
    N["umi_molecules"] = f0(mean([x for c in cells for x in pick(c, "umi", "molecules_analysed")]))
    # Figure S20: how far one PacBio read reaches, for the read-span comparison
    sw_pb = CRM / "out/switch_pacbio/switch_pacbio.json"
    if sw_pb.exists():
        sm = _load(sw_pb)["samples"]
        jl = sorted(v["read_len_median_junction"] for v in sm.values())
        N["pacbio_len"] = f"{round(statistics.median(jl)):,}"
        N["pacbio_n"] = f0(len(sm))
        rpj = sorted(v["reads_per_unique_junction"] for v in sm.values())
        N["pacbio_rpj"] = f"{rpj[0]:.1f}\u2013{rpj[-1]:.1f}"
    fig2n = _load(OUT / "fig2_numbers.json")
    N.update(af_n=f"{fig2n['af_calls']:,}", af_lt1=pct(fig2n["af_share_below"]["1"] / 100),
             af_lt5=pct(fig2n["af_share_below"]["5"] / 100), af_lt20=f"{fig2n['af_share_below']['20']:.1f}")
    jr = _load(OUT / "intuition_jh4.json")["raster"]
    N.update(rast_wt_called=f"{sum(m[2] for r in jr['WT3'] for m in r['muts']):,}",
             rast_wt_other=f0(sum(1 - m[2] for r in jr["WT3"] for m in r["muts"])),
             rast_ko_other=f0(sum(1 - m[2] for r in jr["KO1"] for m in r["muts"])))
    hiq = []          # measured error at reported quality >= 40 in each run's controls
    for k in ("jh4", "n114", "su_lps", "pim1"):
        q = _load(OUT / f"intuition_{k}.json")["quality"]
        m, b = (sum(v[i] for qq, v in q.items() if int(qq) >= 40) for i in (0, 1))
        hiq.append(m / b)
    N.update(q40_lo=f"{min(hiq) * 1e5:.0f}", q40_hi=f"{max(hiq) * 1e5:.0f}")
    # ectopic AID in non-B cells (PRJNA170048, Perez-Duran et al.): the dead-AID control reveals the
    # UNG-dependent spectrum that every-mismatch counting blurs (out/public/morange)
    mo = _load(OUT / "public/morange/project.json")
    au = _load(OUT / "public/morange/authors_counts.json")["libraries"]
    mv = defaultdict(list)
    for v in mo["variants"]:
        mv[v["sample"]].append(v)
    ms = {s["sample"]: s for s in mo["samples"]}
    for arm in ("Control", "Ugi", "UNG"):
        vs = [v for e in (1, 2, 3) for v in mv[f"Exp{e}_{arm}_AID"]]
        gc = [v for v in vs if v["ref"] in "CG"]
        N[f"mo_{arm.lower()}_ts"] = pct(sum((v["ref"], v["alt"]) in (("C", "T"), ("G", "A"))
                                            for v in gc) / len(gc))
        N[f"mo_{arm.lower()}_ts_raw"] = pct(statistics.mean(
            au[f"Exp{e}_{arm}_AID"]["pct_transition_at_GC"] for e in (1, 2, 3)) / 100)
        N[f"mo_{arm.lower()}_calls"] = f0(sum(ms[f"Exp{e}_{arm}_AID"]["unique_mutations"] for e in (1, 2, 3)))
        N[f"mo_{arm.lower()}_dead"] = f0(sum(ms[f"Exp{e}_{arm}_AIDE58Q"]["unique_mutations"] for e in (1, 2, 3)))
        # Two traps here, both of which inflate the fold.
        # First, mutation_frequency_x1e3 is rounded to four decimals, so a dead-AID library with a few
        # mutations in 1.24e8 bases stores as exactly 0.0 and looks mutation-free; the frequency is therefore
        # recomputed from the counts. Second, a library that genuinely saw nothing sits at the frequency one
        # event would give on its own denominator, the convention of Figure 4F and Table S6, not at a fixed
        # constant: a constant sets the fold by itself (1e-6 gave Ugi 8,102-fold, 1e-5 would give 810).
        def g(k, arm=arm):
            def one(e):
                x = ms[f"Exp{e}_{arm}_{k}"]
                t, b = x.get("mutations_total") or 0, x["frequency_bases"]
                return 1000.0 * t / b if t else 1000.0 / b
            return statistics.geometric_mean([one(e) for e in (1, 2, 3)])
        N[f"mo_{arm.lower()}_fold"] = f0(g("AID") / g("AIDE58Q"))
        # how much of each denominator is a floor rather than a measurement, so the text can say so
        N[f"mo_{arm.lower()}_dead_zero"] = f0(sum(
            not (ms[f"Exp{e}_{arm}_AIDE58Q"].get("mutations_total") or 0) for e in (1, 2, 3)))
        gr = lambda k, arm=arm: statistics.geometric_mean(  # noqa: E731
            [au[f"Exp{e}_{arm}_{k}"]["freq_x1e3"] for e in (1, 2, 3)])
        N[f"mo_{arm.lower()}_fold_raw"] = f1(gr("AID") / gr("AIDE58Q"))
    # the repair-context swing as an effect with an interval, paired within experiment (minor 7)
    try:
        _sh = np.array(json.load(open(OUT / "fig7_morange.json"))["shmtrace"], dtype=float)
        _d = _sh[:, 0] - _sh[:, 2]                      # Ugi minus UNG, one value per experiment
        _se = _d.std(ddof=1) / len(_d) ** 0.5
        _tc = stats.t.ppf(0.975, len(_d) - 1)
        N["mo_arm_diff"] = f"{_d.mean():.0f}"
        N["mo_arm_ci"] = f"{_d.mean() - _tc * _se:.0f}–{_d.mean() + _tc * _se:.0f}"
        N["mo_arm_p"] = "paired t, p = " + pfmt(float(stats.ttest_rel(_sh[:, 0], _sh[:, 2]).pvalue))
    except (OSError, ValueError, KeyError):
        N["mo_arm_diff"] = N["mo_arm_ci"] = N["mo_arm_p"] = "[rebuild Figure 7]"
    N["mo_dead_ts"] = (f'{min(au[f"Exp{e}_{a}_AIDE58Q"]["pct_transition_at_GC"] for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG")):.0f}'
                       f'–{max(au[f"Exp{e}_{a}_AIDE58Q"]["pct_transition_at_GC"] for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG")):.0f}')
    allv = [v for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG") for v in mv[f"Exp{e}_{a}_AID"]]
    N["mo_calls"] = f0(len(allv))
    N["mo_gc"] = pct(sum(v["ref"] in "CG" for v in allv) / len(allv))
    N["mo_wrc"] = (f'{min(ms[f"Exp{e}_{a}_AID"]["pct_GC_in_WRC_GYW"] for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG")):.0f}'
                   f'–{max(ms[f"Exp{e}_{a}_AID"]["pct_GC_in_WRC_GYW"] for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG")):.0f}')
    N["mo_raw_lo"] = f2(min(au[f"Exp{e}_{a}_AIDE58Q"]["freq_x1e3"] for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG")))
    N["mo_raw_hi"] = f2(max(au[f"Exp{e}_{a}_AIDE58Q"]["freq_x1e3"] for e in (1, 2, 3) for a in ("Control", "Ugi", "UNG")))
    N["mo_noisy_calls"] = f0(ms["Exp3_Control_AIDE58Q"]["unique_mutations"])
    N["mo_noisy_gc"] = f1(ms["Exp3_Control_AIDE58Q"]["pct_GC"])
    N["mo_noisy_gc_raw"] = f1(au["Exp3_Control_AIDE58Q"]["pct_GC"])
    N["mo_clean_gc_raw"] = (f'{au["Exp1_Control_AIDE58Q"]["pct_GC"]:.0f} and {au["Exp2_Control_AIDE58Q"]["pct_GC"]:.0f}')
    # indels planted into held-out nanopore controls (indel_spike.py)
    for k in ("jh4", "n114"):
        r = indel_spike.summary(k)
        N.update({f"ind_{k}_ge3": pct(r["ge3"]), f"ind_{k}_fp": f2(r["fp"]) if r["fp"] else "no"})
    r = indel_spike.summary("jh4")
    N.update(ind_jh4_hp=pct(r["by_kind_ge3"]["del1_hp"]), ind_jh4_ins2=pct(r["by_kind_ge3"]["ins2"]),
             ind_jh4_long=pct(r["by_kind_ge3"]["dellong"]),
             ind_draws=word(_load(OUT / "indel_spike_jh4.json")["draws"]))
    # a published Sanger study re-analysed automatically (Barbulescu et al. 2024, cohort 2; sanger_cohort2.py)
    f = OUT / "sanger_cohort2.json"
    if f.exists():
        c2 = _load(f)
        pub = _load(CRM / "published/barbulescu2024_jh4.json")
        per = c2["per_mouse"]
        B4 = "GCAT"
        def pspec(m):  # noqa: E306
            mat = pub["mice"][m]["matrix"]
            return {f"{B4[i]}>{B4[j]}": mat[i][j] for i in range(4) for j in range(4) if i != j}
        def gcs(sp):  # noqa: E306
            tot = sum(sp.values())
            return (100 * (sp.get("C>T", 0) + sp.get("G>A", 0)) / tot,
                    100 * sum(sp.get(k, 0) for k in ("C>A", "C>G", "G>C", "G>T")) / tot)
        mice = list(pub["mice"])
        g_of = lambda m: per[m]["group"]  # noqa: E731
        sh_f = {m: per[m]["freq_x1e4"] for m in per}
        # the paper's own plotted values (its source data) where they exist, else its printed tables
        src = pub.get("source_data", {})
        pb_f = {m: src.get("fig2h_x1e4", {}).get(m, 1e4 * pub["mice"][m]["unique_mutations"] / pub["mice"][m]["nt"])
                for m in mice}
        wt = [m for m in mice if g_of(m) == "WT"]
        mk = [m for m in mice if g_of(m) == "Maea_KO"]
        mean = lambda d, ms: float(np.mean([d[m] for m in ms]))  # noqa: E731
        N["c2_traces"] = f"{c2['traces']}"
        N["c2_clones"] = f"{c2['clones_analysed']}"
        N["c2_seconds"] = f"{c2['seconds']:.0f}"
        N["c2_fold"] = f1(mean(sh_f, wt) / mean(sh_f, mk))
        N["c2_pub_fold"] = f1(mean(pb_f, wt) / mean(pb_f, mk))
        N["c2_p"] = pfmt(float(stats.ttest_ind([sh_f[m] for m in wt], [sh_f[m] for m in mk], alternative="greater")[1]))
        N["c2_pub_p"] = f"{pub['fig2h_p']['WT_vs_Maea_KO']}"
        ts = {m: gcs(per[m]["spectrum"]) for m in wt + mk}
        tsp = {m: tuple(src["fig2i_ts_tv_pct"][m]) if m in src.get("fig2i_ts_tv_pct", {}) else gcs(pspec(m))
               for m in wt + mk}
        for k, i in (("ts", 0), ("tv", 1)):
            N[f"c2_{k}_wt"] = f"{np.mean([ts[m][i] for m in wt]):.0f}"
            N[f"c2_{k}_mk"] = f"{np.mean([ts[m][i] for m in mk]):.0f}"
            N[f"c2_{k}_p"] = pfmt(float(stats.ttest_ind([ts[m][i] for m in wt], [ts[m][i] for m in mk])[1]))
            N[f"c2_{k}_pub_wt"] = f"{np.mean([tsp[m][i] for m in wt]):.0f}"
            N[f"c2_{k}_pub_mk"] = f"{np.mean([tsp[m][i] for m in mk]):.0f}"
        N["c2_ts_pub_p"], N["c2_tv_pub_p"] = f"{pub['fig2i_p']['Ts']}", f"{pub['fig2i_p']['Tv']}"
        N["c2_rho"] = f2(float(stats.spearmanr([pb_f[m] for m in mice], [sh_f[m] for m in mice])[0]))
        N["c2_ratio"] = f1(10 ** float(np.mean([math.log10(sh_f[m] / pb_f[m]) for m in mice])))
        types = sorted(pspec(mice[0]))
        xs, ys = [], []
        for grp in (wt, mk):
            sh = {t_: sum(per[m]["spectrum"].get(t_, 0) for m in grp) for t_ in types}
            pb = {t_: sum(pspec(m)[t_] for m in grp) for t_ in types}
            xs += [pb[t_] / sum(pb.values()) for t_ in types]
            ys += [sh[t_] / sum(sh.values()) for t_ in types]
        N["c2_spec_r"] = f2(float(stats.pearsonr(xs, ys)[0]))
        grp = {x["group"]: x for x in c2["project"]["groups"]}
        ko = [sh_f[m] for m in per if per[m]["group"] == "Aicda_KO"]
        N["c2_ko_freq"] = f1(float(np.mean(ko)))
        N["c2_ko_n"] = f"{sum(per[m]['clones'] for m in per if per[m]['group'] == 'Aicda_KO')}"
        N["c2_ko_calls"] = word(sum(s["unique_mutations"] for s in c2["project"]["samples"] if s["group"] == "Aicda_KO"))
        N["c2_wt_ko_p"] = pfmt(grp["WT"]["mutation_frequency_x1e3_p_welch_log_vs_Aicda_KO"])
        N["c2_mk_ko_p"] = pfmt(grp["Maea_KO"]["mutation_frequency_x1e3_p_welch_log_vs_Aicda_KO"])
        N["c2_dup_flags"] = f"{sum(1 for x in c2['project']['qc']['flags'] if x['check'] == 'duplicate_clone')}"
        N["c2_filed_twice"] = word(sum(1 for x in c2["project"]["qc"]["flags"] if x["check"] == "duplicate_file") // 2)
    # the second published Sanger study: Barbulescu et al. (2025), Mkln1 (sanger_cohort2.py mkln)
    f = OUT / "sanger_mkln.json"
    if f.exists():
        mk = _load(f)
        pubm = _load(CRM / "published/barbulescu2025_jh4_mkln.json")
        perm = mk["per_mouse"]
        B5 = "ATCG"
        def pspec_m(m):  # noqa: E306  (rows = To, columns = From, as printed)
            mat = pubm["mice"][m]["matrix"]
            return {f"{B5[j]}>{B5[i]}": mat[i][j] for i in range(4) for j in range(4) if i != j}
        def tv_m(sp):  # noqa: E306
            tot = sum(sp.values())
            return 100 * sum(sp.get(k, 0) for k in ("C>A", "C>G", "G>C", "G>T")) / tot
        wtm = [m for m in perm if perm[m]["group"] == "WT"]
        kom = [m for m in perm if perm[m]["group"] == "Mkln1_KO"]
        pw = [1e4 * pubm["mice"][m]["f_mut"] for m in pubm["mice"] if pubm["mice"][m]["group"] == "WT"]
        pk = [1e4 * pubm["mice"][m]["f_mut"] for m in pubm["mice"] if pubm["mice"][m]["group"] == "Mkln1_KO"]
        shw, shk = [perm[m]["freq_x1e4"] for m in wtm], [perm[m]["freq_x1e4"] for m in kom]
        N["mk_traces"] = f"{mk['traces']}"
        N["mk_fold"] = f1(float(np.mean(shw) / np.mean(shk)))
        N["mk_pub_fold"] = f1(float(np.mean(pw) / np.mean(pk)))
        N["mk_p"] = pfmt(float(stats.ttest_ind(shw, shk, alternative="greater")[1]))
        N["mk_pub_p"] = pfmt(float(stats.ttest_ind(pw, pk, alternative="greater")[1]))
        tvw, tvk = [tv_m(perm[m]["spectrum"]) for m in wtm], [tv_m(perm[m]["spectrum"]) for m in kom]
        ptw = [tv_m(pspec_m(m)) for m in pubm["mice"] if pubm["mice"][m]["group"] == "WT"]
        ptk = [tv_m(pspec_m(m)) for m in pubm["mice"] if pubm["mice"][m]["group"] == "Mkln1_KO"]
        N["mk_tv_wt"], N["mk_tv_ko"] = f"{np.mean(tvw):.0f}", f"{np.mean(tvk):.0f}"
        N["mk_tv_p"] = pfmt(float(stats.ttest_ind(tvw, tvk)[1]))
        N["mk_tv_pub_wt"], N["mk_tv_pub_ko"] = f"{np.mean(ptw):.0f}", f"{np.mean(ptk):.0f}"
        N["mk_tv_pub_p"] = pfmt(float(stats.ttest_ind(ptw, ptk)[1]))
        TY = sorted(pspec_m(next(iter(pubm["mice"]))))
        xs, ys = [], []
        for grp_s, grp_p, key in ((wtm, "WT", "mk_spec_r_wt"), (kom, "Mkln1_KO", "mk_spec_r_ko")):
            sh = {k: sum(perm[m]["spectrum"].get(k, 0) for m in grp_s) for k in TY}
            pb = {k: sum(pspec_m(m)[k] for m in pubm["mice"] if pubm["mice"][m]["group"] == grp_p) for k in TY}
            gx, gy = [pb[k] / sum(pb.values()) for k in TY], [sh[k] / sum(sh.values()) for k in TY]
            N[key] = f2(float(stats.pearsonr(gx, gy)[0]))
            xs += gx
            ys += gy
        N["mk_spec_r"] = f2(float(stats.pearsonr(xs, ys)[0]))
        N["mk_paralog"] = f"{mk['excluded'].get('contaminant', 0)}"
    # a published Sanger study from another laboratory (Thientosapol et al. 2018; external/swhel_samhd1.py)
    sw_f = CRM.parent / "external/out/swhel_samhd1.json"
    if sw_f.exists():
        sw = _load(sw_f)
        swp = _load(CRM / "published/thientosapol2018_swhel.json")
        per = sw["per_mouse"]
        arms = ["Samhd1_WT", "Samhd1_KO", "Samhd1_KO_UGI", "Msh2_KO", "Msh2_Samhd1_KO"]
        mice = {g: [m for m in per if per[m]["group"] == g] for g in arms}
        tot: dict = {g: Counter() for g in arms}
        for v in sw["project"]["variants"]:
            tot[per[v["sample"]]["group"]][f"{v['ref']}>{v['alt']}"] += v["reads"]
        types = sorted(swp["groups"]["Samhd1_WT"]["m"])
        xs = [swp["groups"][g]["m"][ty] / sum(swp["groups"][g]["m"].values()) for g in arms for ty in types]
        ys = [tot[g][ty] / sum(tot[g].values()) for g in arms for ty in types]
        N["sw_clones"] = f0(sum(per[m]["clones"] for g in arms for m in mice[g]))
        N["sw_bp"] = "558"
        N["sw_muts"] = f0(sum(sum(tot[g].values()) for g in arms))
        N["sw_pub_muts"] = f0(sum(sum(swp["groups"][g]["m"].values()) for g in arms))
        N["sw_r"] = f"{float(stats.pearsonr(xs, ys)[0]):.4f}"

        def _cls(m, keys):
            sp = per[m]["spectrum"]
            s = sum(sp.values())
            return 100 * sum(sp.get(k, 0) for k in keys) / s if s else 0.0
        CL = {"at_tv": ("A>T", "A>C", "T>A", "T>G"), "at_ts": ("A>G", "T>C"),
              "gc_tv": ("C>A", "C>G", "G>C", "G>T"), "gc_ts": ("C>T", "G>A")}
        for k, keys in CL.items():
            a_ = [_cls(m, keys) for m in mice["Samhd1_WT"]]
            b_ = [_cls(m, keys) for m in mice["Samhd1_KO"]]
            pv = float(stats.ttest_ind(a_, b_)[1])
            # all five classes carry a bare p value, as the sentences that use them supply "p ="
            N[f"sw_{k}_p"] = pfmt(pv)
        # total mutation per animal, Welch on log frequency (the study tested it by ANOVA with Holm-Sidak, p = 0.07)
        lf = lambda ms: [math.log10(per[m]["freq_x1e4"]) for m in ms]  # noqa: E731
        N["sw_tot_p"] = pfmt(float(stats.ttest_ind(lf(mice["Samhd1_WT"]), lf(mice["Samhd1_KO"]),
                                                   equal_var=False)[1]))
        # One deposited wild-type host name carries about twice the clones of any other, so it likely stands for
        # two animals; the animal-level tests are repeated without it rather than splitting it without evidence.
        big = max(mice["Samhd1_WT"], key=lambda m: per[m]["clones"])
        rest = [m for m in mice["Samhd1_WT"] if m != big]
        N["sw_big_host"], N["sw_big_n"] = big.split("HOST-")[-1], f0(per[big]["clones"])
        N["sw_wt_hosts"] = f0(len(mice["Samhd1_WT"]))
        N["sw_drop_p"] = series(pfmt(float(stats.ttest_ind([_cls(m, CL[k]) for m in rest],
                                                            [_cls(m, CL[k]) for m in mice["Samhd1_KO"]])[1]))
                                for k in ("at_tv", "gc_tv", "gc_ts", "at_ts"))
        N["sw_drop_tot_p"] = pfmt(float(stats.ttest_ind(lf(rest), lf(mice["Samhd1_KO"]), equal_var=False)[1]))
        AT = ("A>G", "T>C", "A>T", "A>C", "T>A", "T>G")
        N["sw_at_wt"] = f"{np.mean([_cls(m, AT) for m in mice['Samhd1_WT']]):.0f}"
        N["sw_at_msh2"] = f"{np.mean([_cls(m, AT) for m in mice['Msh2_KO']]):.0f}"
        # the evidence-based clone set (sanger_evidence.py): lineages filed under two mice resolved by nanopore
        ev_path = OUT / "sanger_mkln_evidence.json"
        if ev_path.exists():
            ev = _load(ev_path)["per_mouse"]
            ew = [ev[m]["freq_x1e4"] for m in ev if ev[m]["group"] == "WT"]
            ek = [ev[m]["freq_x1e4"] for m in ev if ev[m]["group"] == "Mkln1_KO"]
            N["mk_ev_fold"] = f1(float(np.mean(ew) / np.mean(ek)))
            N["mk_ev_p"] = pfmt(float(stats.ttest_ind(np.log10(ew), np.log10(ek), equal_var=False)[1]))
    # ---------------- PacBio class-switch junctions (crm/switch_pacbio.py)
    sw_path = CRM / "out/switch_pacbio/switch_pacbio.json"
    if sw_path.exists():
        sw = _load(sw_path)
        S_, G_ = sw["samples"], sw["groups"]
        by = {g: [v for k, v in S_.items() if G_[k] == g] for g in set(G_.values())}
        comp = {(r["group"], r["feature"]): r for r in sw["comparisons"]}
        uj = [v["unique_junctions"] for v in S_.values()]
        rpj = [v["reads_per_unique_junction"] for v in S_.values()]
        amb = [v["acceptor_ambiguity_direct"]["ambiguous_pct"] for v in S_.values()]
        mh = {g: [v["mean_microhomology"] for v in vs] for g, vs in by.items()}
        p4 = {g: [v["microhomology_classes_pct"].get(">=4 nt", 0.0) for v in vs] for g, vs in by.items()}
        bl = {g: [v["microhomology_classes_pct"].get("blunt", 0.0) for v in vs] for g, vs in by.items()}
        N.update(pb_uj_lo=f0(min(uj)), pb_uj_hi=f0(max(uj)),
                 pb_rpj_lo=f1(min(rpj)), pb_rpj_hi=f0(max(rpj)),
                 pb_amb_lo=f0(min(amb)), pb_amb_hi=f0(max(amb)),
                 pb_mh_wt=f2(np.mean(mh["WT"])), pb_mh_xlf=f2(np.mean(mh["Xlf"])),
                 pb_mh_shld=f2(np.mean(mh["Shld1"])), pb_mh_dko=f2(np.mean(mh["Shld1_Xlf"])),
                 pb_p_xlf=pfmt(comp[("Xlf", "mean_microhomology")]["p_welch"]),
                 pb_p_shld=pfmt(comp[("Shld1", "mean_microhomology")]["p_welch"]),
                 pb_p_dko=pfmt(comp[("Shld1_Xlf", "mean_microhomology")]["p_welch"]),
                 pb_4_wt=f0(np.mean(p4["WT"])), pb_4_xlf=f0(np.mean(p4["Xlf"])), pb_4_dko=f0(np.mean(p4["Shld1_Xlf"])),
                 pb_bl_wt=f0(np.mean(bl["WT"])), pb_bl_dko=f0(np.mean(bl["Shld1_Xlf"])),
                 pb_p4_dko=pfmt(comp[("Shld1_Xlf", "pct_mh_4plus")]["p_welch"]),
                 pb_pbl_dko=pfmt(comp[("Shld1_Xlf", "pct_blunt")]["p_welch"]))
    # ---------------- public Illumina Smu, PI3K inhibitors (crm/public_smu_pi3k.py)
    smu_path = CRM / "out/public_smu/public_smu.json"
    if smu_path.exists():
        sm = _load(smu_path)
        tw, tk = {}, {}
        for t, v in sm["treatments"].items():
            tw[t] = [x["unique_mutations"] for x in v["samples"] if x["group"] == "WT"]
            tk[t] = [x["unique_mutations"] for x in v["samples"] if x["group"] == "AID_KO"]
        allw = [x for t in tw for x in tw[t]]
        allk = [x for t in tk for x in tk[t]]
        N.update(smu_wt_dmso=f1(np.mean(tw["DMSO"])), smu_wt_idela=f1(np.mean(tw["Idelalisib"])),
                 smu_wt_duve=f1(np.mean(tw["Duvelisib"])), smu_ko_hi=f0(max(allk)),
                 smu_wt_lo=f0(min(allw)), smu_wt_hi=f0(max(allw)),
                 smu_gt_p=pfmt(float(stats.mannwhitneyu(allw, allk, alternative="greater").pvalue)),
                 smu_idela_p=pfmt(float(stats.ttest_ind(tw["Idelalisib"], tw["DMSO"], equal_var=False).pvalue)),
                 smu_duve_p=pfmt(float(stats.ttest_ind(tw["Duvelisib"], tw["DMSO"], equal_var=False).pvalue)))
        rp = sm.get("repeats_paired")
        if rp:
            rw = [x["unique_mutations"] for x in rp["samples"] if x["group"] == "WT"]
            rk = [x["unique_mutations"] for x in rp["samples"] if x["group"] == "AID_KO"]
            fw = [x["mutation_frequency_x1e3"] for x in rp["samples"] if x["group"] == "WT"]
            fk = [x["mutation_frequency_x1e3"] for x in rp["samples"] if x["group"] == "AID_KO"]
            N.update(smur_wt=f1(np.mean(rw)), smur_ko=f1(np.mean(rk)), smur_ko_tot=word(sum(rk)),
                     smur_wt_lo=f0(min(rw)), smur_wt_hi=f0(max(rw)), smur_ko_hi=f0(max(rk)),
                     smur_wt_freq=f2(np.mean(fw)), smur_ko_freq=f2(np.mean(fk)),
                     smur_p=pfmt(float(stats.mannwhitneyu(rw, rk, alternative="greater").pvalue)))
            s5 = [x["s5f_index"] for x in rp["samples"] if x["group"] == "WT" and x["s5f_index"]]
            ctga = [x["pct_CT_GA"] for x in rp["samples"] if x["group"] == "WT" and x["pct_CT_GA"]]
            N.update(smur_s5f_lo=f2(min(s5)), smur_s5f_hi=f2(max(s5)),
                     smur_ctga_lo=f0(min(ctga)), smur_ctga_hi=f0(max(ctga)))
    # ---------------- engineered AID diversification (crm/liu_mutators.py, tam_*.py, tcr_bassan.py,
    # devilder_clones.py). crm/out is not under version control, so each block is optional.
    liu_p = CRM / "out/liu/liu_mutators.json"
    if liu_p.exists():
        ls = _load(liu_p)["samples"]
        get = lambda pfx, k: [v[k] for n_, v in ls.items() if n_.startswith(pfx) and v[k] is not None]  # noqa: E731
        for tag, aid, apo in (("wt", "AID_wt293T", "A3A_wt293T"), ("ko", "AID_ungmsh2", "A3A_ungmsh2")):
            for key, short in (("pct_GC_in_WRC_GYW", "wrc"), ("pct_CT_GA", "ctga")):
                x, y = get(aid, key), get(apo, key)
                N[f"eng_{tag}_{short}_aid"] = f0(np.mean(x))
                N[f"eng_{tag}_{short}_apo"] = f0(np.mean(y))
                N[f"eng_{tag}_{short}_p"] = pfmt(float(stats.ttest_ind(x, y, equal_var=False).pvalue))
        nc = [v["unique_mutations"] for n_, v in ls.items() if n_.startswith("NC")]
        N["eng_nc_hi"] = f0(max(nc))
    foot_p = CRM / "out/tam/footprint.json"
    if foot_p.exists():
        fb = _load(foot_p)["best_per_guide"]
        clean = [v for v in fb.values() if v["rate_per_1e3"] > 50]
        N["eng_guides_clean"] = word(len(clean))
        N["eng_guide_lo"] = f0(min(v["rate_per_1e3"] for v in clean))
        N["eng_guide_hi"] = f0(max(v["rate_per_1e3"] for v in clean))
        N["eng_guide_fold_lo"] = f0(min(v["fold"] for v in clean))
        N["eng_guide_fold_hi"] = f0(max(v["fold"] for v in clean))
    sat_p = CRM / "out/tam/saturation.json"
    if sat_p.exists():
        t_ = _load(sat_p)["totals"]
        N["eng_sat"] = ", ".join(str(t_[f"{k}_sgRNAs"]) for k in (3, 5, 7, 10))
    tam_p = CRM / "out/tam/tam_reporter.json"
    if tam_p.exists():
        cal = _load(tam_p)["calibrated_run"]["samples"]
        N["eng_tam_ctrl"] = f0(cal["mock"]["unique_mutations"])
        N["eng_tam_case"] = word(cal["dCas9_AIDX_sgRNA"]["unique_mutations"])
        raw = _load(tam_p)["raw_mismatch_per_1e3"]
        amp = "amp_CCTGACCTACGGCGTGCAGT"
        N["eng_tam_raw_case"] = f0(raw["dCas9_AIDX_sgRNA"][amp])
        N["eng_tam_raw_ctrl"] = f0(raw["mock"][amp])
    dev_p = CRM / "out/devilder/devilder_clones.json"
    if dev_p.exists():
        rnds = _load(dev_p)["rounds"]
        ann = [c for v in rnds.values() for c in v if c["published_mutations"]]
        ok = [c for c in ann if not c["published_not_recovered"] and not c["recovered_not_published"]]
        N["eng_dev_n"] = f0(len(ann))
        N["eng_dev_ok"] = f0(len(ok))
        N["eng_dev_mut"] = f0(sum(len(c["published_mutations"]) for c in ann))
        # four clusters recur between the two rounds, so the entries outnumber the clusters
        N["eng_dev_distinct"] = f0(len({c["cluster"] for c in ann}))
    tcr_p = CRM / "out/tcr/tcr_bassan.json"
    if tcr_p.exists():
        ts = _load(tcr_p)["samples"]
        base = ts["h5_wt"]["mutation_frequency_x1e3"]
        oth = [v["mutation_frequency_x1e3"] for k, v in ts.items() if k != "h5_wt"]
        N["eng_tcr_base"] = f2(base)
        N["eng_tcr_lo"] = f2(min(oth))
        N["eng_tcr_hi"] = f2(max(oth))
    import shmtrace
    N["shmtrace_version"] = shmtrace.__version__
    N.update(lab_details.placeholders())
    # Table S7 and its STAR Methods paragraph: the default against two stricter profiles
    # (crm/strict_sensitivity.py)
    st_path = OUT / "strict_sensitivity.json"
    if st_path.exists():
        st = _load(st_path)
        lin, sp, pj = st["lineage"]["polyclonal"], st["spikein"], st["projects"]
        N.update(st_lin_fp_def=f1(lin["default"]["false_per_lib"]), st_lin_fp=f1(lin["strict"]["false_per_lib"]),
                 st_lin_fp_max=f1(lin["strict_max"]["false_per_lib"]),
                 st_lin_1_def=pct(lin["default"]["recall_1read"]), st_lin_1=pct(lin["strict"]["recall_1read"]),
                 st_lin_fdr=f"{lin['strict']['fdr']:.3f}",
                 st_jh4_rec_def=pct(sp["jh4"]["default"]["recall"]), st_jh4_rec=pct(sp["jh4"]["strict"]["recall"]),
                 st_jh4_fdr=f"{sp['jh4']['strict']['fdr']:.3f}", st_n114_rec=pct(sp["n114"]["strict"]["recall"]),
                 st_n114_fdr=f"{sp['n114']['strict']['fdr']:.3f}",
                 st_jh4_rec_max=pct(sp["jh4"]["strict_max"]["recall"]))
        strict = [v for ds in pj.values() for k, v in ds.items() if k != "default"]
        N.update(st_p_hi=pfmt(max(v["p_log_welch"] for v in strict)),
                 st_fold_lo=f0(min(v["fold"] for v in strict)), st_fold_hi=f0(max(v["fold"] for v in strict)))
    # Figure 4I: the usable reads each library yielded, against the 500 drawn for analysis. This was the
    # point of Figure S1E, which was dropped on 21 October because it repeated Figure 4I.
    usable = [v["usable"] for run in ("jh4_shmtrace", "n114_shmtrace")
              for v in _load(FIG / f"{run}/table_s1.json").values()]
    if usable:
        N["reads_usable_lo"] = f"{min(usable):,}"
        N["reads_usable_hi"] = f"{max(usable):,}"
        N["reads_usable_n"] = str(len(usable))
    # Figure S3: how many figures the J_H4 project run actually drew, of the 26 a run can produce
    jh4_run = FIG / "jh4_shmtrace/project.json"
    if jh4_run.exists():
        N["s3_drawn"] = word(len(_load(jh4_run).get("available_figures", [])))
    # Figure S4: the false-discovery target against planted truth at six targets (figs/figS4_fdr_calibration.py)
    fc_path = OUT / "fdr_calibration.json"
    if fc_path.exists():
        fc = _load(fc_path)
        at = lambda t: {k: v[t] for k, v in fc.items()}  # noqa: E731
        d, one = at("0.05"), at("0.01")
        pc1 = lambda x: f"{100 * x:.1f}"  # noqa: E731
        loose = [(k, float(t), v["fdp_pooled"]) for k, per in fc.items() for t, v in per.items() if float(t) > 0.05]
        over = [x for x in loose if x[2] > x[1]]
        three = [k for k in fc if k != "n114"]
        LBL = {"jh4": "J_{H}4", "n114": "N114", "su_lps": "Sμ", "pim1": "*Pim1*"}
        N.update(fc_def_lo=pc1(min(v["fdp_pooled"] for v in d.values())),
                 fc_def_hi=pc1(max(v["fdp_pooled"] for v in d.values())),
                 fc_def_lib_max=pc1(max(max(v["fdp_libraries"]) for v in d.values())),
                 fc_libs=f0(d["jh4"]["libraries"]),
                 fc_loose_over=series([f"{LBL[k]} at {100 * t:g}% ({pc1(f)}%)" for k, t, f in over]),
                 fc_strict_lo=pc1(min(one[k]["fdp_pooled"] for k in three)),
                 fc_strict_hi=pc1(max(one[k]["fdp_pooled"] for k in three)),
                 fc_n114_1=pc1(one["n114"]["fdp_pooled"]),
                 fc_floor_lo=f1(min(v["false_per_library"] for v in one.values())),
                 fc_floor_hi=f1(max(v["false_per_library"] for v in one.values())),
                 fc_def_fp_lo=f1(min(v["false_per_library"] for v in d.values())),
                 fc_def_fp_hi=f1(max(v["false_per_library"] for v in d.values())),
                 fc_jh4_rec_def=pct(d["jh4"]["recall"]), fc_jh4_rec_1=pct(one["jh4"]["recall"]),
                 fc_su_rec_def=pct(d["su_lps"]["recall"]), fc_su_rec_1=pct(one["su_lps"]["recall"]))
    # Figure S5: how many controls, and which (crm/controls_design.py, figs/figS5_controls.py)
    cd_path = OUT / "controls_design.json"
    if cd_path.exists():
        cd = _load(cd_path)
        J, S = cd["jh4"], cd["su"]
        mean = lambda xs: sum(xs) / len(xs)  # noqa: E731
        lin = {n: J[f"polyclonal|{n}|none|on"] for n in (1, 2, 3, 4)}
        nulls = [J[f"null|{n}|none|on"] for n in (1, 2, 3, 4)]
        noisy = [J[f"polyclonal|{n}|one_noisy|{r}"] for n in (3, 4) for r in ("on", "off")]
        other = [S["other experiment (2)"], S["other experiment (3)"]]
        pc1 = lambda x: f"{100 * x:.1f}"  # noqa: E731
        N.update(cd_runs=f0(lin[1]["runs"]),
                 cd_rec_1=pct(mean(lin[1]["recall"])), cd_rec_2=pct(mean(lin[2]["recall"])),
                 cd_rec_34_lo=pct(min(mean(lin[n]["recall"]) for n in (3, 4))),
                 cd_rec_34_hi=pct(max(mean(lin[n]["recall"]) for n in (3, 4))),
                 cd_fp_1=f1(mean(lin[1]["false"])),
                 cd_fp_lo=f1(min(mean(lin[n]["false"]) for n in (2, 3, 4))),
                 cd_fp_hi=f1(max(mean(lin[n]["false"]) for n in (2, 3, 4))),
                 cd_fdp_lo=pc1(min(lin[n]["fdp_pooled"] for n in (2, 3, 4))),
                 cd_fdp_hi=pc1(max(lin[n]["fdp_pooled"] for n in (2, 3, 4))),
                 cd_null_any=word(sum(c["any_call"] for c in nulls)), cd_null_n=f0(sum(c["runs"] for c in nulls)),
                 cd_noisy_fp_lo=f1(min(mean(c["false"]) for c in noisy)),
                 cd_noisy_fp_hi=f1(max(mean(c["false"]) for c in noisy)),
                 cd_noisy_flagged=f0(min(J[f"polyclonal|{n}|one_noisy|on"]["set_aside"] for n in (3, 4))),
                 cd_clean_flag_lo=f0(min(J[f"polyclonal|{n}|none|on"]["set_aside"] for n in (3, 4))),
                 cd_clean_flag_hi=f0(max(J[f"polyclonal|{n}|none|on"]["set_aside"] for n in (3, 4))),
                 cd_noisy_null_any=word(J["null|4|one_noisy|on"]["any_call"]),
                 cd_su_rec_matched=pct(S["matched (2)"]["recall"]),
                 cd_su_rec_other_lo=pct(min(a["recall"] for a in other)),
                 cd_su_rec_other_hi=pct(max(a["recall"] for a in other)),
                 cd_su_fdp_matched=pc1(S["matched (2)"]["fdp_pooled"]),
                 cd_su_fdp_other_lo=pc1(min(a["fdp_pooled"] for a in other)),
                 cd_su_fdp_other_hi=pc1(max(a["fdp_pooled"] for a in other)),
                 cd_su_fp_matched=f1(S["matched (2)"]["false_per_library"]),
                 cd_su_fp_pooled=f1(S["pooled (5)"]["false_per_library"]),
                 cd_su_rec_pooled=pct(S["pooled (5)"]["recall"]),
                 cd_su_fdp_pooled=pc1(S["pooled (5)"]["fdp_pooled"]))
    # The dilution series: real hypermutation diluted into AID-null reads (crm/titration.py). Its figure was removed at the author’s request; the numbers stay in the Results.
    ti_path = OUT / "titration_summary.json"
    if ti_path.exists():
        ti = _load(ti_path)
        don = {k: v for k, v in ti.items() if isinstance(v, dict) and "points" in v}
        det = lambda loc: [100 * v["detected_from"] for v in don.values() if v["locus"] == loc]  # noqa: E731
        rat = [x for v in don.values() if v["locus"] == "n114" for x in v["ratio_25_50"]]
        jh4 = {v["donor"]: v for v in don.values() if v["locus"] == "jh4"}
        most = max(jh4.values(), key=lambda v: v["unique_at_full"])
        least = min(jh4.values(), key=lambda v: v["detected_from"] is None or -v["detected_from"])
        N.update(ti_reps=f0(len(next(iter(don.values()))["points"][0]["observed"])),
                 ti_n114_det_lo=f"{min(det('n114')):g}", ti_n114_det_hi=f"{max(det('n114')):g}",
                 ti_n114_ratio_lo=f2(min(rat)), ti_n114_ratio_hi=f2(max(rat)),
                 ti_jh4_det_lo=f"{min(det('jh4')):g}", ti_jh4_det_hi=f"{max(det('jh4')):g}",
                 ti_jh4_most=most["donor"].replace("WT", "Wt"), ti_jh4_least=least["donor"].replace("WT", "Wt"))
        dflt = ti.get("n114_default", {})
        if dflt:
            N.update(ti_d3_def_50=f2(dflt["D3|0.5"]), ti_d3_def_100=f2(dflt["D3|1.0"]))
    # Figure S11: the Immcantation route on the same N114 reads (crm/immcantation_compare.py, figs/figS11_immcantation.py)
    im_path = OUT / "immcantation_summary.json"
    if im_path.exists():
        im = _load(im_path)
        R = im["routes"]
        N.update(im_vcall=im["v_call"], im_vid=f1(im["v_identity"]), im_vlen=f0(im["v_length"]),
                 im_fold_imgt=f0(R["imgt"]["fold"]), im_fold_par=f0(R["parental"]["fold"]),
                 im_fold_q20=f0(R["parental_q20"]["fold"]), im_fold_sh=f0(R["shmtrace"]["fold"]),
                 im_ctl_imgt=f2(1e3 * R["imgt"]["empty_geo"]), im_ctl_q20=f2(1e3 * R["parental_q20"]["empty_geo"]),
                 im_aid_q20=f1(1e3 * R["parental_q20"]["aid_geo"]), im_aid_sh=f1(1e3 * R["shmtrace"]["aid_geo"]))
        # a range whose ends round alike is printed once
        for k in ("cd_rec_34", "cd_su_rec_other", "cd_fdp", "cd_su_fdp_other", "ti_n114_det", "ti_n114_ratio"):
            lo, hi = N[f"{k}_lo"], N[f"{k}_hi"]
            N[k] = lo if lo == hi else f"{lo}–{hi}"
    return N


if __name__ == "__main__":
    for k, v in compute().items():
        print(f"{k:22} {v}")
