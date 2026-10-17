"""Per-mouse substitution classes, Sanger against nanopore, for Figure 5E.

    python crm/fig5_per_mouse.py        -> crm/out/fig5_per_mouse.json

Figure 5B pools the three wild-type mice, and the first review asked for the paired per-mouse view: a pooled
correlation can hide mice that disagree. This recomputes each mouse's calls on both platforms from the same
inputs Figure 5 uses (the manuscript Sanger sample sheet; the J_H4 nanopore project at 500 reads; the 529-bp
window) and, before writing anything, checks that pooling them reproduces Figure 5B's pooled spectrum exactly.
If it does not, it stops, because then the panel would not be showing the same data as the rest of the figure.

The four classes are descriptive groupings of the 12 substitution types, not an attribution to repair
pathways.
"""
import json
import sys
from pathlib import Path

CRM = Path(__file__).resolve().parent
CAL = CRM.parent
sys.path.insert(0, str(CAL))
sys.path.insert(0, str(CRM))
import sanger_vs_nanopore as sv  # noqa: E402

from shmtrace import project, samplesheet  # noqa: E402

REF = CAL.parents[1] / "src/shmtrace/references"
STATS = CAL / "figures/sanger_vs_nanopore/recommended_manuscript/stats.json"
CLASSES = {"G:C transitions": ("C>T", "G>A"), "G:C transversions": ("C>A", "C>G", "G>C", "G>T"),
           "A:T transitions": ("A>G", "T>C"), "A:T transversions": ("A>T", "A>C", "T>A", "T>G")}


def calls():
    sheet = samplesheet.load(str(CAL / "samplesheets/jh4_manuscript.csv"))
    res = project.run_project(sheet.project(str(REF / "mouse_jh4_intron.fasta"), "jh4_intron",
                                            contaminants=[str(REF / "mouse_jh4_contaminant.fasta")]))
    sanger = {m: [f"{v.ref}>{v.alt}" for v in res["_details"][m]["variants"] if sv.LO <= v.pos <= sv.HI]
              for m in sv.WT}
    ont_proj = json.load(open(CAL / "figures/jh4_shmtrace/project.json"))
    ont = {m: [f"{v['ref']}>{v['alt']}" for v in ont_proj["variants"]
               if v["sample"] == m and sv.LO <= v["pos"] <= sv.HI] for m in sv.WT}
    return {"sanger": sanger, "ont": ont}


def main():
    c = calls()
    target = json.load(open(STATS))["spectrum"]["pooled_props"]
    for plat, d in c.items():
        allt = [t for m in sv.WT for t in d[m]]
        pooled = {t: round(100 * allt.count(t) / len(allt), 2) for t in sv.MUT}
        if pooled != target[plat]:
            sys.exit(f"{plat}: pooled spectrum does not reproduce Figure 5B; not writing the panel data")
    out = {"mice": sv.WT, "classes": list(CLASSES), "calls": {}, "share": {}}
    for plat, d in c.items():
        out["calls"][plat] = {m: len(d[m]) for m in sv.WT}
        out["share"][plat] = {m: [100 * sum(t in CLASSES[k] for t in d[m]) / len(d[m]) for k in CLASSES]
                              for m in sv.WT}
    diffs = [out["share"]["ont"][m][i] - out["share"]["sanger"][m][i] for m in sv.WT for i in range(len(CLASSES))]
    gcts = [out["share"]["ont"][m][0] - out["share"]["sanger"][m][0] for m in sv.WT]
    out["summary"] = {"max_abs_diff": max(abs(x) for x in diffs), "gc_ts_diff": gcts,
                      "rm_anova_interaction_p": json.load(open(STATS))["spectrum"]["rm_anova_interaction"]["p"]}
    (CRM / "out").mkdir(exist_ok=True)
    json.dump(out, open(CRM / "out/fig5_per_mouse.json", "w"), indent=1)
    print("reproduces Figure 5B pooled spectrum; max per-mouse class difference "
          f"{out['summary']['max_abs_diff']:.1f} points; G:C transitions {[round(x, 1) for x in gcts]}")


if __name__ == "__main__":
    main()
