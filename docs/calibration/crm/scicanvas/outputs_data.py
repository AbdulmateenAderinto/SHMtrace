"""Numbers for the SHMtrace outputs overview (figS21_outputs.js), from the J_H4 project run.

    python crm/scicanvas/outputs_data.py

Reads docs/calibration/figures/jh4_shmtrace (project.json, table_s1.json) and writes them into the builder's DATA,
one small example per kind of output: read accounting, the control background, density per animal, the spectrum,
position along the amplicon and clones.
"""
import json
from pathlib import Path

from inject import inject

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[1] / "figures" / "jh4_shmtrace"
CLASSES = ["A>T", "A>C", "A>G", "T>A", "T>C", "T>G", "C>A", "C>T", "C>G", "G>A", "G>T", "G>C"]


def main():
    res = json.loads((RUN / "project.json").read_text())
    reads = json.loads((RUN / "table_s1.json").read_text())
    samples = res["samples"]
    group = {s["sample"]: s["group"] for s in samples}
    lo, hi = res["params"]["analysis_window"]

    # read accounting: % of delivered reads left after each gate (as qc_waterfall draws it)
    qc = {}
    for name, r in reads.items():
        n = r["delivered"]
        gated = n - r["too_long"] - r["low_quality"]
        on = gated - r["not_target"]
        usable = r["usable"]
        qc[name] = [100.0, 100 * gated / n, 100 * on / n, 100 * usable / n, 100 * r["analysed"] / n]

    bg = res["control_background"]
    xs = sorted(int(k) for k in bg)
    background = {"pos": xs, "rate": [bg[str(x)] for x in xs]}

    tests = {t["metric"]: t for t in res["animal_tests"]}
    dens = {"rows": [[s["sample"], s["group"], s["unique_density_x1e5"], s["one_event_density_x1e5"]]
                     for s in samples],
            "p": tests["unique_density"]["p_welch_log"], "fold": tests["unique_density"]["fold_geometric"]}

    wt = [s["sample"] for s in samples if not s["control"]]
    spec = {"classes": CLASSES,
            "group": [100 * res["spectrum_group"]["Wt"].get(c, 0) / sum(res["spectrum_group"]["Wt"].values())
                      for c in CLASSES],
            "animals": [[100 * res["spectrum_sample"][a].get(c, 0) / max(1, sum(res["spectrum_sample"][a].values()))
                         for c in CLASSES] for a in wt],
            "gc_hot": [s["pct_GC_in_WRC_GYW"] for s in samples if not s["control"]]}

    pos = {"wt": [[p["pos"], p["frequency_pct"]] for p in res["positional"] if p["group"] == "Wt"],
           "ko": [[p["pos"], p["frequency_pct"]] for p in res["positional"] if p["group"] != "Wt"],
           "hot": [h for h in res["hotspots"] if lo <= h <= hi], "window": [lo, hi]}

    clones = [[name, group[name], row["reads"], row["mutations"]]
              for name, c in res["clones"].items() for row in c["rows"]]
    largest = {name: c["largest_clone_pct"] for name, c in res["clones"].items()}

    data = {"qc": qc, "bg": background, "dens": dens, "spec": spec, "pos": pos,
            "clones": clones, "largest": largest, "drawn": len(res["available_figures"])}
    print(f"{inject(HERE / 'figS21_outputs.js', data):+d} characters")


if __name__ == "__main__":
    main()
