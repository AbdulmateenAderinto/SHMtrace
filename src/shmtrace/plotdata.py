"""The numbers behind every figure, as tidy tables to graph however you like.

The app's figures are a quick look. These tables are for drawing your own: one row per
observation, one column per variable, plain column names with the unit in the name, and a
data dictionary that says what each column is. They load as they are into R (ggplot2),
Python (pandas, seaborn), Excel or Prism, so the choice of chart, scale, colours and
statistics stays with the person making the figure.

- `<figure>.csv`: exactly what that figure in the app draws.
- `per_sample_long.csv`: every per-sample measure stacked (sample, group, measure, value,
  unit), for plotting any measure against any grouping.
- `group_statistics_long.csv`: the groups table, one statistic per row.
- `prism/<measure>.csv`: one column per group, for pasting into a Prism column table.
- `data_dictionary.csv` and `README.txt`.

Nothing here is computed anew: every number is one the run already reported.
"""
from __future__ import annotations

import csv
import io
import re
import zipfile

# per-sample measures: key in the samples table -> (column name here, unit, what it is)
MEASURES = {
    "unique_density_x1e5": ("unique_mutations_per_1e5_bp", "per 10^5 bp",
                            "unique mutations (each position and base once) per 100,000 analysed bases"),
    "mutation_frequency_x1e3": ("mutations_per_1e3_bp", "per 10^3 bp",
                                "every mutated base in every read, per 1,000 bases sequenced"),
    "mutated_reads_pct": ("reads_mutated_pct", "%", "reads carrying at least one mutation"),
    "mean_mutations_per_read": ("mutations_per_read_mean", "per read", "mean mutations per read"),
    "pct_GC_in_WRC_GYW": ("gc_mutations_in_wrc_gyw_pct", "%", "G:C mutations at WRC/GYW hotspots"),
    "pct_AT_in_WA_TW": ("at_mutations_in_wa_tw_pct", "%", "A:T mutations at WA/TW hotspots"),
    "pct_GC": ("mutations_at_gc_pct", "%", "unique mutations at a G or C"),
    "pct_transitions": ("transitions_pct", "%", "transitions among unique mutations"),
    "clone_weighted_frequency_x1e3": ("clone_weighted_mutations_per_1e3_bp", "per 10^3 bp",
                                      "mutation frequency with each clone counted once"),
    "single_read_density_x1e5": ("one_read_mutations_per_1e5_bp", "per 10^5 bp",
                                 "one-read mutations net of the controls' background"),
    "combined_density_x1e5": ("all_mutations_per_1e5_bp", "per 10^5 bp",
                              "called mutations plus the one-read estimate"),
    "unique_mutations": ("unique_mutations", "count", "unique mutations called"),
    "reads": ("reads", "count", "reads analysed"),
}

# column -> (unit, meaning), for the data dictionary; columns named as in MEASURES are added
COLUMNS = {
    "sample": ("", "sample (library) name"), "group": ("", "group the sample belongs to"),
    "control": ("yes/no", "the sample is a control (AID-null, empty vector …)"),
    "measure": ("", "which per-sample measure the row holds"), "value": ("see unit", "the measure's value"),
    "unit": ("", "unit of value"),
    "zero_drawn_at": ("per 10^5 bp", "where a sample with no mutations is drawn on a log axis: the density one mutation would give"),
    "analysed_bases": ("bases", "denominator of the unique-mutation density"),
    "bases_sequenced": ("bases", "denominator of the mutation frequency"),
    "mutations_total": ("count", "every mutated base in every read"),
    "gc_mutations": ("count", "G:C mutations the hotspot percentage rests on"),
    "from_base": ("", "reference base"), "to_base": ("", "base read instead"),
    "substitution": ("", "from>to"),
    "count": ("count", "unique mutations"),
    "percent_of_group": ("%", "share of the group's unique mutations"),
    "percent_of_sample": ("%", "share of the sample's unique mutations"),
    "mean_percent": ("%", "mean over the group's samples of each sample's share"),
    "sd_percent": ("%", "standard deviation of that share over the samples"),
    "samples": ("count", "samples with at least one call"),
    "group_calls": ("count", "the group's unique mutations; the app leaves out a group with fewer than 20"),
    "mutations_per_read": ("count", "mutations in one read"),
    "load_class": ("", "mutations per read, binned as the app's donut chart bins them"),
    "share_pct": ("%", "share of the group's reads"),
    "position": ("1-based", "position on the analysis reference"),
    "reference_base": ("", "reference base at the position"),
    "mutated_reads": ("count", "reads mutated at the position"),
    "wrc_gyw_hotspot": ("yes/no", "the position is in a WRC/GYW motif"),
    "compared_with": ("", "control group the P values test against"),
    "statistic": ("", "mean, sd, or a test: welch, welch_log (Welch on log10), mann_whitney, and sensitivity variants"),
}

LOAD_CLASSES = [("0", 0, 0), ("1-3", 1, 3), ("4-10", 4, 10), ("11-20", 11, 20), ("20+", 21, 10 ** 9)]
BASES = "ATCG"


def _ctrl(s) -> str:
    return "yes" if s.get("control") else "no"


def per_sample(res: dict, key: str, extra: dict[str, str] | None = None) -> list[dict]:
    col = MEASURES[key][0]
    out = []
    for s in res["samples"]:
        row = {"sample": s["sample"], "group": s["group"], "control": _ctrl(s), col: s.get(key)}
        for src, name in (extra or {}).items():
            row[name] = s.get(src)
        out.append(row)
    return out


def per_sample_long(res: dict) -> list[dict]:
    return [{"sample": s["sample"], "group": s["group"], "control": _ctrl(s),
             "measure": name, "value": s[key], "unit": unit}
            for s in res["samples"] for key, (name, unit, _m) in MEASURES.items() if s.get(key) is not None]


def substitution_matrix(res: dict) -> list[dict]:
    out = []
    for g in res["group_order"]:
        c = res["spectrum_group"].get(g) or {}
        tot = sum(c.values())
        for a in BASES:
            for b in BASES:
                if a != b:
                    n = c.get(f"{a}>{b}", 0)
                    out.append({"group": g, "from_base": a, "to_base": b, "count": n,
                                "percent_of_group": round(100 * n / tot, 3) if tot else None})
    return out


def spectrum_by_sample(res: dict) -> list[dict]:
    gmap = {s["sample"]: s["group"] for s in res["samples"]}
    out = []
    for s in res["samples"]:
        c = (res.get("spectrum_sample") or {}).get(s["sample"]) or {}
        tot = sum(c.values())
        for a in BASES:
            for b in BASES:
                if a != b:
                    n = c.get(f"{a}>{b}", 0)
                    out.append({"sample": s["sample"], "group": gmap[s["sample"]], "control": _ctrl(s),
                                "substitution": f"{a}>{b}", "from_base": a, "to_base": b, "count": n,
                                "percent_of_sample": round(100 * n / tot, 3) if tot else None})
    return out


def spectrum_by_group(res: dict) -> list[dict]:
    """The spectrum chart: per group, the mean and SD of each sample's share."""
    per = spectrum_by_sample(res)
    out = []
    for g in res["group_order"]:
        rows = [r for r in per if r["group"] == g and r["percent_of_sample"] is not None]
        calls = sum((res["spectrum_group"].get(g) or {}).values())
        for a in BASES:
            for b in BASES:
                if a == b:
                    continue
                xs = [r["percent_of_sample"] for r in rows if r["substitution"] == f"{a}>{b}"]
                m = sum(xs) / len(xs) if xs else None
                sd = (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5 if len(xs) > 1 else None
                out.append({"group": g, "substitution": f"{a}>{b}", "from_base": a, "to_base": b,
                            "mean_percent": None if m is None else round(m, 3),
                            "sd_percent": None if sd is None else round(sd, 3),
                            "samples": len(xs), "count": (res["spectrum_group"].get(g) or {}).get(f"{a}>{b}", 0),
                            "group_calls": calls})
    return out


def mutations_per_read(res: dict) -> list[dict]:
    gmap = {s["sample"]: s["group"] for s in res["samples"]}
    return [{"sample": s, "group": gmap.get(s, ""), "mutations_per_read": int(k), "reads": n}
            for s, hist in res["mutation_load"].items() for k, n in sorted(hist.items(), key=lambda x: int(x[0]))]


def mutations_per_read_classes(res: dict) -> list[dict]:
    """The donut chart: reads per load class, per group."""
    per = mutations_per_read(res)
    out = []
    for g in res["group_order"]:
        rows = [r for r in per if r["group"] == g]
        tot = sum(r["reads"] for r in rows)
        for name, lo, hi in LOAD_CLASSES:
            n = sum(r["reads"] for r in rows if lo <= r["mutations_per_read"] <= hi)
            out.append({"group": g, "load_class": name, "reads": n,
                        "share_pct": round(100 * n / tot, 3) if tot else None})
    return out


def positional_by_group(res: dict) -> list[dict]:
    hot = set(res.get("hotspots") or [])
    return [{"group": x["group"], "position": x["pos"], "reference_base": x.get("ref", ""),
             "mutated_reads": x["mutated_reads"], "reads_mutated_pct": x["frequency_pct"],
             "wrc_gyw_hotspot": "yes" if x["pos"] in hot else "no"} for x in res["positional"]]


def positional_by_sample(res: dict) -> list[dict]:
    ref, hot = res.get("reference") or "", set(res.get("hotspots") or [])
    return [{"sample": x["sample"], "group": x["group"], "position": x["pos"],
             "reference_base": ref[x["pos"] - 1] if 0 < x["pos"] <= len(ref) else "",
             "mutated_reads": x["mutated_reads"], "reads_mutated_pct": x["frequency_pct"],
             "wrc_gyw_hotspot": "yes" if x["pos"] in hot else "no"} for x in res.get("positional_sample") or []]


_STAT = re.compile(r"^(?P<m>.+?)_(?P<s>mean|sd|p_welch_log|p_welch|p_mwu|sensitivity_p_welch_log_half_event"
                   r"|sensitivity_p_welch_log_third_of_smallest)(?:_vs_.*)?$")
_STAT_NAMES = {"p_welch": "welch", "p_welch_log": "welch_log", "p_mwu": "mann_whitney",
               "sensitivity_p_welch_log_half_event": "welch_log_zeros_at_half_an_event",
               "sensitivity_p_welch_log_third_of_smallest": "welch_log_zeros_at_a_third_of_the_smallest"}


def group_statistics_long(res: dict) -> list[dict]:
    out = []
    for g in res.get("groups") or []:
        for k, v in g.items():
            m = _STAT.match(k)
            if not m or v is None:
                continue
            key = m["m"]
            out.append({"group": g["group"], "control": _ctrl(g), "compared_with": g.get("compared_with") or "",
                        "measure": MEASURES.get(key, (key,))[0], "statistic": _STAT_NAMES.get(m["s"], m["s"]),
                        "value": v})
    return out


def prism(res: dict) -> dict[str, list[list]]:
    """Per measure, one column per group (Prism's column table); rows are the samples."""
    out = {}
    for key, (name, _u, _m) in MEASURES.items():
        cols = {g: [s.get(key) for s in res["samples"] if s["group"] == g and s.get(key) is not None]
                for g in res["group_order"]}
        if not any(cols.values()):
            continue
        n = max(len(v) for v in cols.values())
        out[f"prism/{name}.csv"] = [list(cols)] + [[(cols[g][i] if i < len(cols[g]) else "") for g in cols]
                                                   for i in range(n)]
    return out


def tables(res: dict) -> dict[str, list[dict]]:
    """Every tidy table, keyed by file name; a figure's table is named after the figure."""
    t = {
        "density_by_group.csv": per_sample(res, "unique_density_x1e5", {
            "one_event_density_x1e5": "zero_drawn_at", "unique_mutations": "unique_mutations",
            "density_bases": "analysed_bases"}),
        "mutation_frequency_by_group.csv": per_sample(res, "mutation_frequency_x1e3", {
            "mutations_total": "mutations_total", "frequency_bases": "bases_sequenced"}),
        "mutated_reads_by_group.csv": per_sample(res, "mutated_reads_pct", {"reads": "reads"}),
        "hotspot_targeting.csv": per_sample(res, "pct_GC_in_WRC_GYW", {"gc_mutations": "gc_mutations"}),
        "substitution_matrix.csv": substitution_matrix(res),
        "spectrum_by_group.csv": spectrum_by_group(res),
        "spectrum_by_sample.csv": spectrum_by_sample(res),
        "mutation_load.csv": mutations_per_read_classes(res),
        "mutations_per_read.csv": mutations_per_read(res),
        "positional_map.csv": positional_by_group(res),
        "positional_by_sample.csv": positional_by_sample(res),
        "per_sample_long.csv": per_sample_long(res),
        "group_statistics_long.csv": group_statistics_long(res),
    }
    return {k: v for k, v in t.items() if v}


def dictionary(tabs: dict[str, list[dict]]) -> list[dict]:
    named = {name: (unit, meaning) for name, unit, meaning in MEASURES.values()}
    out = []
    for f, rows in tabs.items():
        for c in dict.fromkeys(k for r in rows for k in r):
            unit, meaning = COLUMNS.get(c) or named.get(c) or ("", "")
            out.append({"file": f, "column": c, "unit": unit, "meaning": meaning})
    return out


README = """These tables hold the numbers behind SHMtrace's figures, for drawing your own.

One row per observation, one column per variable; units are in the column names or in
data_dictionary.csv. Each figure's table is named after the figure (density_by_group.csv is
the mutation density chart, and so on). per_sample_long.csv stacks every per-sample measure,
so any measure can be plotted against any grouping. The prism/ folder has one file per
measure with one column per group, for pasting into a Prism column table.

Every number is one the run already reported; nothing is recalculated. A sample with no
mutations has a density of 0; on a log axis SHMtrace draws it at zero_drawn_at, the
density one mutation would have given.

In R:      d <- read.csv("density_by_group.csv")
           ggplot(d, aes(group, unique_mutations_per_1e5_bp)) + geom_point() + scale_y_log10()
In Python: d = pandas.read_csv("density_by_group.csv")
           seaborn.stripplot(data=d, x="group", y="unique_mutations_per_1e5_bp")
"""


def _csv(rows) -> str:
    buf = io.StringIO()
    if rows and isinstance(rows[0], dict):
        w = csv.DictWriter(buf, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)
    else:
        csv.writer(buf).writerows(rows)
    return buf.getvalue()


def files(res: dict) -> dict[str, str]:
    """Every file of the plot-data folder, keyed by its path inside the folder."""
    tabs = tables(res)
    out = {name: _csv(rows) for name, rows in tabs.items()}
    out.update({name: _csv(rows) for name, rows in prism(res).items()})
    out["data_dictionary.csv"] = _csv(dictionary(tabs))
    out["README.txt"] = README
    return out


def zip_bytes(res: dict) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in files(res).items():
            z.writestr(name, text)
    return buf.getvalue()
