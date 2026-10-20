"""The default caller against two stricter profiles: Table S7 and its STAR Methods paragraph.

    python crm/spikein.py jh4 n114 --no-lofreq --tag _strict \\
        --profiles recommended,recommended_fdr01,recommended_max
    shmtrace-simulate-lineages ... --modes recommended,strict,strict_max --out crm/out/lineage_sim_jh4_strict
    python crm/strict_sensitivity.py [--rerun-projects]      -> crm/out/strict_sensitivity.json

The default's control-based threshold targets a false-discovery proportion of 0.05, so a library with many
true mutations admits more false calls in absolute terms. "strict" lowers the target to 0.01; "strict_max"
keeps 0.05 but counts the noisiest control's decoys instead of their mean. Both are compared with the
default on the planted-mutation spike-ins, the simulated lineages, and the J_H4 and N114 projects, where the
question is whether the group conclusions depend on the default's permissiveness.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

CRM = Path(__file__).resolve().parent
OUT = CRM / "out"
PROFILES = {"default": ("recommended", "recommended"), "strict": ("recommended_fdr01", "strict"),
            "strict_max": ("recommended_max", "strict_max")}     # name: (spike-in profile, project/lineage mode)


def projects():
    """The J_H4 and N114 projects at each profile (as jh4_shmtrace.py and n114_shmtrace.py, 500 reads)."""
    from shmtrace import project
    sys.path.insert(0, str(CRM.parent))
    import jh4_shmtrace as j  # noqa: E402
    import n114_shmtrace as n  # noqa: E402
    sets = {"jh4": dict(region="jh4_intron", reference=str(j.REF / "mouse_jh4_intron.fasta"),
                        contaminants=[str(j.REF / "mouse_jh4_contaminant.fasta")],
                        samples=[{"name": k, "path": str(j.fastq(f, "delivered")), "group": g,
                                  "control": g == "Aicda-/-"} for k, (f, g) in j.LIBS.items()]),
            "n114": dict(region="n114_vregion", reference=str(n.REF), contaminants=[],
                         samples=[{"name": k, "path": str(n.fastq(f, "delivered")), "group": g,
                                   "control": g == "pMX-Empty"} for k, (f, g) in n.LIBS.items()])}
    out = {}
    for key, s in sets.items():
        for _name, (_sp, mode) in PROFILES.items():
            res = project.public(project.run_project({"name": f"{key} {mode}", "mode": mode, "depth": 500, **s}))
            out[f"{key}/{mode}"] = {"samples": res["samples"], "groups": res["groups"],
                                    "animal_tests": res.get("animal_tests")}
    json.dump(out, open(OUT / "strict_projects.json", "w"), default=str)
    return out


def main(rerun=False):
    res = {"spikein": {}, "lineage": {}, "projects": {}}
    for k in ("jh4", "n114"):
        agg = defaultdict(lambda: {"tp": 0, "pl": 0, "fp": 0, "n": 0, "one": [0, 0]})
        for r in json.load(open(OUT / f"spikein_{k}_strict.json"))["rows"]:
            a = agg[r["caller"]]
            a["tp"] += r["true_calls"]
            a["pl"] += r["planted"]
            a["fp"] += r["false_calls"]
            a["n"] += 1
            a["one"][0] += r["recall_by_reads"]["1"][0]
            a["one"][1] += r["recall_by_reads"]["1"][1]
        res["spikein"][k] = {name: {"false_per_lib": agg[sp]["fp"] / agg[sp]["n"],
                                    "fdr": agg[sp]["fp"] / max(1, agg[sp]["tp"] + agg[sp]["fp"]),
                                    "recall": agg[sp]["tp"] / agg[sp]["pl"],
                                    "recall_1read": agg[sp]["one"][0] / agg[sp]["one"][1]}
                             for name, (sp, _m) in PROFILES.items()}
    mode_name = {m: name for name, (_sp, m) in PROFILES.items()}
    for row in json.load(open(OUT / "lineage_sim_jh4_strict/lineage_sim.json"))["summary"]:
        res["lineage"].setdefault(row["scenario"], {})[mode_name[row["mode"]]] = {
            "false_per_lib": row["false_calls_per_library"], "fdr": row.get("fdr"),
            "recall_1read": (row.get("recall_by_support") or {}).get("1"), "recall_2plus": row.get("recall_2plus")}
    P = projects() if rerun or not (OUT / "strict_projects.json").exists() else json.load(open(OUT / "strict_projects.json"))
    for key, v in P.items():
        ds, mode = key.split("/")
        g = next(x for x in v["groups"] if not x.get("control"))
        t = {r["metric"]: r for r in v["animal_tests"] or []}["unique_density"]
        case = [s["unique_mutations"] for s in v["samples"] if not s["control"]]
        ctl = [s["unique_mutations"] for s in v["samples"] if s["control"]]
        res["projects"].setdefault(ds, {})[mode_name[mode]] = {
            "case_calls": case, "control_calls": ctl,
            "p_log_welch": next(val for k_, val in g.items() if k_.startswith("unique_density_x1e5_p_welch_log_vs_")),
            "p_negbin": t["p_negbinomial"], "fold": t["fold_geometric"]}
    json.dump(res, open(OUT / "strict_sensitivity.json", "w"), indent=1)
    print(json.dumps(res, indent=1)[:3000])


if __name__ == "__main__":
    main("--rerun-projects" in sys.argv)
