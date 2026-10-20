"""Data for the SciCanvas Figure S10 (figS15_swhel.js).

    python docs/calibration/crm/scicanvas/swhel_data.py

Reads external/out/swhel_samhd1.json (SHMtrace on the deposited clones) and
published/thientosapol2018_swhel.json (the paper's Table 1), reduces each panel to the numbers it
plots, and writes them into figS15_swhel.js between the DATA markers.
"""
import collections
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from inject import inject  # noqa: E402

CRM = HERE.parent
RES = json.load(open(CRM.parent / "external/out/swhel_samhd1.json"))
PUB = json.load(open(CRM / "published/thientosapol2018_swhel.json"))
BUILDER = HERE / "figS15_swhel.js"
TYPES = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]
TS_GC, TV_GC = ("C>T", "G>A"), ("C>A", "C>G", "G>C", "G>T")
TS_AT, TV_AT = ("A>G", "T>C"), ("A>T", "A>C", "T>A", "T>G")
ORDER = ["Samhd1_WT", "Samhd1_KO", "Samhd1_KO_UGI", "Msh2_KO", "Msh2_Samhd1_KO"]
SHORT = {"Samhd1_WT": "wt/het", "Samhd1_KO": "Samhd1", "Samhd1_KO_UGI": "Samhd1+UGI",
         "Msh2_KO": "Msh2", "Msh2_Samhd1_KO": "Msh2 Samhd1"}


def main():
    per = RES["per_mouse"]
    mice = {g: [m for m in per if per[m]["group"] == g] for g in ORDER}
    # total mutations, as the paper's Table 1 counts them: every supporting clone counted
    tot = collections.defaultdict(collections.Counter)
    for v in RES["project"]["variants"]:
        tot[per[v["sample"]]["group"]][f"{v['ref']}>{v['alt']}"] += v["reads"]

    # A: mutations per sequence, per group
    A = []
    for g in ORDER:
        p = PUB["groups"][g]
        n = sum(per[m]["clones"] for m in mice[g])
        A.append({"group": g, "label": SHORT[g], "pub": p["per_sequence"],
                  "sh": round(sum(tot[g].values()) / n, 3),
                  "pub_seqs": p["sequences"], "sh_seqs": n,
                  "pub_muts": sum(p["m"].values()), "sh_muts": sum(tot[g].values())})

    # B: the 12-type spectrum, published against SHMtrace, as a share of each group's mutations
    B, xs, ys = [], [], []
    for g in ORDER:
        p = PUB["groups"][g]["m"]
        tp, ts_ = sum(p.values()), sum(tot[g].values())
        for t in TYPES:
            x, y = 100 * p[t] / tp, 100 * tot[g][t] / ts_
            B.append([round(x, 3), round(y, 3), g, t])
            xs.append(x)
            ys.append(y)
    r = float(stats.pearsonr(xs, ys)[0])

    # C: the four classes the paper's Fig. 2B tests, wild type against Samhd1-null
    def pct(m, keys):
        sp = per[m]["spectrum"]
        s = sum(sp.values())
        return 100 * sum(sp.get(k, 0) for k in keys) / s if s else None

    def pubpct(g, keys):
        p = PUB["groups"][g]["m"]
        return 100 * sum(p[k] for k in keys) / sum(p.values())

    C = []
    for lab, keys in (("A:T Tv", TV_AT), ("A:T Ts", TS_AT), ("C:G Tv", TV_GC), ("C:G Ts", TS_GC)):
        a = [pct(m, keys) for m in mice["Samhd1_WT"]]
        b = [pct(m, keys) for m in mice["Samhd1_KO"]]
        C.append({"label": lab,
                  "wt": [round(v, 2) for v in a], "ko": [round(v, 2) for v in b],
                  "pub_wt": round(pubpct("Samhd1_WT", keys), 2), "pub_ko": round(pubpct("Samhd1_KO", keys), 2),
                  "p": float(stats.ttest_ind(a, b).pvalue)})

    # D: the published claims as compact verdicts, and a footnote of totals
    def stars(p):
        return "ns" if p >= 0.05 else ("*" if p >= 0.01 else ("**" if p >= 0.001 else "****"))
    fw = [per[m]["freq_x1e4"] for m in mice["Samhd1_WT"]]
    fk = [per[m]["freq_x1e4"] for m in mice["Samhd1_KO"]]
    pf = float(stats.ttest_ind(np.log10(fw), np.log10(fk), equal_var=False).pvalue)
    D = [["A:T transversions", "down ****", "down " + stars(C[0]["p"])],
         ["A:T transitions", "unchanged", "unchanged" if C[1]["p"] >= 0.05 else "changed " + stars(C[1]["p"])],
         ["C:G transversions", "down ****", "down " + stars(C[2]["p"])],
         ["C:G transitions", "up **", "up " + stars(C[3]["p"])],
         # the paper tested total mutation by one-way ANOVA, SHMtrace by Welch on per-animal log frequency: the
         # same direction under two different tests, which the panel marks rather than presenting as a match
         ["Total mutationᵃ", "lower, ns", "lower " + stars(pf)]]
    pub_m = sum(sum(PUB["groups"][g]["m"].values()) for g in ORDER)
    sh_m = sum(sum(tot[g].values()) for g in ORDER)
    foot = (f"{sh_m:,} of the {pub_m:,} mutations the paper counts; r = {r:.4f} over the 60 spectrum cells of "
            f"(B); {sum(per[m]['clones'] for g in ORDER for m in mice[g]):,} clones, "
            f"{sum(len(mice[g]) for g in ORDER)} animals, {RES['seconds']:.0f} s.")
    data = {"A": A, "B": B, "r": round(r, 4), "C": C, "D": D, "foot": foot}
    n = inject(BUILDER, data)
    print(f"{BUILDER.name}: data {n:+d} bytes; r = {r:.4f}; "
          f"{sum(sum(tot[g].values()) for g in ORDER)} mutations")


if __name__ == "__main__":
    main()
