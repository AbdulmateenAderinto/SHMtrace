"""Every number the Cell Reports Methods article and the STAR protocol report, computed from the
analysis outputs, so the text never carries a number by hand.

    from crm_numbers import compute; N = compute()      # dict of formatted strings
    python docs/calibration/crm/manuscript/crm_numbers.py   # print them all

Inputs: docs/calibration/figures/* (the manuscript runs), docs/calibration/crm/out/* (benchmarks,
spike-ins, simulations, public data), the Sanger calls, and the simulation outputs of
REPLICATES_AND_POWER and UMI (data-independent).
"""
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import lab_details  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import indel_spike  # noqa: E402

CRM = Path(__file__).resolve().parents[1]
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
                         f"{N['loco_n114_margin']}× on N114")
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
    do, ds = st["density_ont"], st["density_sanger"]
    N.update(ont_fold=f0(do["fold"]), ont_p=pfmt(do["p_welch_log"]), sanger_fold=f0(ds["fold"]),
             sanger_p=pfmt(ds["p_welch_log"]), ont_wt_geo=f0(do["wt_geomean"]), ont_ko_geo=f2(do["ko_geomean"]),
             sanger_wt_geo=f0(ds["wt_geomean"]), sanger_ko_geo=f0(ds["ko_geomean"]))
    sp = st["spectrum"]
    N.update(r_spec=f2(sp["pearson"]["r"]), r_spec_ci=f"{sp['pearson']['ci'][0]:.2f}–{sp['pearson']['ci'][1]:.2f}",
             slope=f2(sp["pearson"]["slope"]), interaction_p=f2(sp["rm_anova_interaction"]["p"]),
             sidak_min=f2(sp["paired_sidak_min"]), n_sanger=f0(sp["sanger_unique"]), n_ont=f0(sp["ont_unique"]))
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
             s1_maxdiff=f1(max(abs(a - b_) for a, b_ in zip(x, y))))
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
        d = _load(FIG / f"depth_shmtrace/{k}.json")
        case = ("WT2", "WT3", "WT4") if k == "jh4" else ("C1", "D3", "D5")
        mx = {x: max(d[dep][x]["calls"] for dep in d) for x in case}
        N[f"depth_{k}_500"] = pct(np.mean([d["500"][x]["calls"] / mx[x] for x in case]))
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
                  f"sens_{k}_fp": f1(r["fp"])})
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
            if r["caller"] == "deepsomatic_permissive":            # as drawn in Figure 3H
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
    # Figure 2 and Figures S2 and S6: anatomy of a call, the control FDR, where real SHM sits, the background along
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
             ana_e_gen_ctl=word(er["geneious_calls"][0]), ana_e_gen_case=word(er["geneious_calls"][1]))
    c = su["candidates"][case]
    N.update(ana_thr=f"{c['threshold']:.3f}", ana_kept=f0(bisect_right(sorted(c["p"]), c["threshold"])),
             ana_cand=f"{len(c['p']):,}")
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
        g = lambda k, arm=arm: statistics.geometric_mean(   # noqa: E731
            [max(ms[f"Exp{e}_{arm}_{k}"]["mutation_frequency_x1e3"] or 0, 1e-6) for e in (1, 2, 3)])
        N[f"mo_{arm.lower()}_fold"] = f0(g("AID") / g("AIDE58Q"))
        gr = lambda k, arm=arm: statistics.geometric_mean(  # noqa: E731
            [au[f"Exp{e}_{arm}_{k}"]["freq_x1e3"] for e in (1, 2, 3)])
        N[f"mo_{arm.lower()}_fold_raw"] = f1(gr("AID") / gr("AIDE58Q"))
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
    import shmtrace
    N["shmtrace_version"] = shmtrace.__version__
    N.update(lab_details.placeholders())
    return N


if __name__ == "__main__":
    for k, v in compute().items():
        print(f"{k:22} {v}")
