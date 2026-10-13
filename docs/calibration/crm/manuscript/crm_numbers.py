"""Every number the Cell Reports Methods article and the STAR protocol report, computed from the
analysis outputs, so the text never carries a number by hand.

    from crm_numbers import compute; N = compute()      # dict of formatted strings
    python docs/calibration/crm/manuscript/crm_numbers.py   # print them all

Inputs: docs/calibration/figures/* (the manuscript runs), docs/calibration/crm/out/* (benchmarks,
spike-ins, simulations, public data), the Sanger calls, and the simulation outputs of
REPLICATES_AND_POWER and UMI (data-independent).
"""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

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


def word(n):
    return ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"][n] if n < 10 else f"{n:,}"


def sci(p):
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10^{{−{abs(int(e))}}}"


def pfmt(p):
    return sci(p) if p < 0.001 else (f"{p:.4f}" if p < 0.01 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}")


def spike(key, caller):
    rows = [r for r in _load(OUT / f"spikein_{key}.json")["rows"] if r["caller"] == caller]
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
             sp_n114_bh_fdr=f"{spike('n114_bh', 'recommended')['fdr']:.3f}",
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
    N["jh4_wt_rec"] = ", ".join(str(j["recommended"]["case_calls"][x]) for x in ("WT2", "WT3", "WT4"))
    N["jh4_ko_rec"] = ", ".join(str(j["recommended"]["ctl_each"][x]) for x in ("KO1", "KO2", "KO3"))
    N["jh4_wt_gen"] = ", ".join(str(j["geneious"]["case_calls"][x]) for x in ("WT2", "WT3", "WT4"))
    N["jh4_ko_gen"] = ", ".join(str(j["geneious"]["ctl_each"][x]) for x in ("KO1", "KO2", "KO3"))
    N["jh4_wt_lf"] = ", ".join(str(j["lofreq"]["case_calls"][x]) for x in ("WT2", "WT3", "WT4"))
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
                         f"{N['loco_jh4_margin']}× on J_{{H}}4, {N['loco_full_margin']}× on the full length and "
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
    N["hmces_ko3_calls"] = str(next(x for x in hm["samples"] if x["sample"] == "AIDKO_3")["unique_mutations"])
    k3 = _load(OUT / "public/hmces_k3/project.json")
    mut = [x for x in k3["samples"] if x["group"] in ("WT", "HMCES_KO")]
    calls = [x["unique_mutations"] for x in mut]
    gcs = [x["pct_GC"] for x in mut if x.get("pct_GC") is not None]
    N["hmk3_calls"] = f"{min(calls)}–{max(calls)}"
    N["hmk3_gc"] = f"{min(gcs):.0f}–{max(gcs):.0f}" if gcs else "–"
    N["hmces_gc_range"] = f"{min(gc):.0f}–{max(gc):.0f}" if gc else "[no substitution calls; see note]"
    return N


if __name__ == "__main__":
    for k, v in compute().items():
        print(f"{k:22} {v}")
