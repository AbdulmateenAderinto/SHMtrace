"""Indels in real libraries: the paper's nanopore projects and three Sanger clone collections, re-run with
indel calling on and nothing else changed.

    python crm/real_indels.py        -> crm/out/real_indels.json

CALIBRATION EVIDENCE, NOT A MANUSCRIPT FIGURE. These libraries have no known indel phenotype and no group
comparison is tested here, so nothing in this file is a biological result and none of it is plotted in the
paper. It is kept for two things it did establish:

  - The 28-bp deletion at position 576, with 5-bp microhomology at its ends, is the whole of the J_H4
    planted-indel benchmark's false-call rate (Figure 3I), and it recurs in AID-null libraries of two
    different mice and on both platforms. It is a PCR product, not a calling error.
  - The 2-bp insertion/deletion pair at 278/280 appeared in all four wild-type mice of the Maea cohort. It is
    carried only by clones SHMtrace's own QC already flags as one molecule filed under several animals, so it
    is one event, not four. Indels carried only by such clones are set aside here and listed in `set_aside`.

The planted-indel benchmark (indel_spike.py, Figure 3) measures recall on synthetic events. This asks what the
caller finds in real data, and whether the platforms agree:

  jh4_ont       J_H4 nanopore, the paper's six libraries at 500 reads (datasets.py, as jh4_extras.py)
  jh4_sanger    Sanger clones of the same six mice (samplesheets/jh4_manuscript.csv, as fig5_per_mouse.py)
  n114_ont      N114 nanopore, the paper's six subclones at 500 reads
  maea_sanger   Barbulescu 2024 cohort 2, Wt / Maea-deficient / Aicda-/- (sanger_cohort2.py cohort2)
  mkln_sanger   Barbulescu 2025, Wt / Mkln1-/- / Aicda-/- (sanger_cohort2.py mkln)

Every run is the published configuration with `exclude_indels: false` added, so the substitution calls are
unchanged and are reported beside the indels as a check. Nanopore indels are tested against the AID-null
peers' own indel rate and need three supporting reads at Q30 (indels.py). Sanger indels are every indel with
its bases at Q20, as substitutions are: one read per clone, so there is no test, and the AID-null clones'
count is the background to read them against.

Two bookkeeping rules, both applied here and both reported:
  - A Sanger indel carried only by clones that SHMtrace's QC flags as one molecule filed under two or more
    animals (duplicate_clone / duplicate_file) is set aside, since it cannot be assigned to an animal; it is
    listed under `set_aside`. Within one animal, several clones carrying one indel are one lineage and the
    indel counts once, as for substitutions.
  - The Mkln1 sheet re-uses the Wt and Aicda-/- mice of jh4_manuscript.csv (SANGER.md), so only its Mkln1-/-
    mice are new; summaries count the shared mice once.

Two cross-platform checks for the same six J_H4 mice (`checks`): whether nanopore events at >= 2% of reads
turn up in that mouse's Sanger clones, for substitutions and for indels alike, so an indel discordance can be
told from a sampling one; and, for each nanopore indel, the lineage its reads share (substitutions carried by
>= 80% of them), how many of that lineage's reads carry it, and whether that lineage is among the Sanger clones.
"""
import json
import sys
from collections import Counter
from pathlib import Path

CRM = Path(__file__).resolve().parent
CAL = CRM.parent
REPO = CAL.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(CRM))
from shmtrace import project, samplesheet  # noqa: E402

from datasets import DATASETS  # noqa: E402

REF = REPO / "src/shmtrace/references"
OUT = CRM / "out/real_indels.json"
INDELS_ON = {"exclude_indels": False}


def nanopore(key, region):
    ds = DATASETS[key]
    return {"name": f"{key}_indels", "region": region, "reference": ds["reference"],
            "contaminants": ds.get("contaminants", []), "mode": "recommended", "depth": ds["depth"],
            "overrides": dict(INDELS_ON),
            "samples": [{"name": n, "path": s["fastq"], "group": "AID-null" if s["control"] else "AID",
                         "control": s["control"]} for n, s in ds["samples"].items()]}


def sanger(sheet):
    sh = samplesheet.load(str(CAL / "samplesheets" / sheet))
    proj = sh.project(str(REF / "mouse_jh4_intron.fasta"), "jh4_intron",
                      contaminants=[str(REF / "mouse_jh4_contaminant.fasta")])
    proj["overrides"] = {**(proj.get("overrides") or {}), **INDELS_ON}
    return proj


RUNS = {"jh4_ont": lambda: nanopore("jh4", "jh4_intron"),
        "jh4_sanger": lambda: sanger("jh4_manuscript.csv"),
        "n114_ont": lambda: nanopore("n114", "n114_vregion"),
        "maea_sanger": lambda: sanger("cohort2_barbulescu2024.csv"),
        "mkln_sanger": lambda: sanger("jh4_mkln_barbulescu2025.csv")}


def _ref():
    return "".join(x.strip() for x in open(REF / "mouse_jh4_intron.fasta") if not x.startswith(">")).upper()


def _carries(a, row):
    want = "D" if row["kind"] == "del" else "I"
    for r, k, b, *_ in a.indels:
        length = b if isinstance(b, int) else len(b)
        if r + 1 == row["pos"] and k == want and length == row["length"]:
            return True
    return False


def summarise(pub, res=None, sanger_run=False):
    rows = pub.get("indels") or []
    dup = {f["clone"] for f in pub.get("qc", {}).get("flags", [])
           if f["check"] in ("duplicate_clone", "duplicate_file")} if sanger_run else set()
    set_aside = []
    if sanger_run and res is not None:
        kept = []
        for r in rows:
            carriers = [a.qname for a in res["_details"][r["sample"]]["_aligned"] if _carries(a, r)]
            if carriers and all(c in dup for c in carriers):
                set_aside.append({**{k: r[k] for k in ("sample", "pos", "kind", "length")}, "clones": carriers})
            else:
                kept.append(r)
        rows = kept
    samples = {}
    for s in pub["samples"]:
        mine = [r for r in rows if r["sample"] == s["sample"]]
        samples[s["sample"]] = {
            "group": s.get("group"), "control": bool(s.get("control")), "reads": s.get("reads"),
            "unique_mutations": s.get("unique_mutations"),
            "indels_as_analysed": s.get("indels_called", 0),          # before setting aside cross-animal copies
            "indels_called": len(mine), "insertions": sum(r["kind"] == "ins" for r in mine),
            "deletions": sum(r["kind"] == "del" for r in mine),
            "frameshift_indels": s.get("frameshift_indels", 0),
            "indel_density_x1e5": s.get("indel_density_x1e5"),
            "unique_density_x1e5": s.get("unique_density_x1e5"),
            "indels": [{k: r[k] for k in ("pos", "kind", "length", "bases", "reads", "depth", "frequency",
                                          "homopolymer", "frameshift")} for r in mine]}
    return samples, set_aside


def main(only=None):
    out = json.loads(OUT.read_text()) if OUT.exists() else {}
    for key, build in RUNS.items():
        if only and key not in only:
            continue
        print(f"== {key}", flush=True)
        res = project.run_project(build())
        pub = project.public(res)
        samples, set_aside = summarise(pub, res, sanger_run=key.endswith("_sanger"))
        groups = sorted({v["group"] for v in samples.values()})
        out[key] = {"notes": [n for n in pub.get("notes", []) if "indel" in n.lower()],
                    "samples": samples, "set_aside": set_aside,
                    "by_group": {g: {"samples": sum(v["group"] == g for v in samples.values()),
                                     "indels": sum(v["indels_called"] for v in samples.values() if v["group"] == g),
                                     "deletions": sum(v["deletions"] for v in samples.values() if v["group"] == g),
                                     "insertions": sum(v["insertions"] for v in samples.values() if v["group"] == g),
                                     "substitutions": sum(v["unique_mutations"] or 0 for v in samples.values()
                                                          if v["group"] == g)}
                                 for g in groups}}
        if key in ("jh4_ont", "jh4_sanger"):
            RESULTS[key] = res
        OUT.write_text(json.dumps(out, indent=1, default=str))
        for e in set_aside:
            print(f"   set aside: {e['kind']}{e['length']}@{e['pos']} in {e['sample']}, clones {e['clones']} "
                  "(one molecule filed under two or more animals)")
        for name, v in samples.items():
            lens = Counter(i["length"] for i in v["indels"])
            print(f"   {name:16s} {str(v['group']):10s} reads {v['reads']!s:>5}  subs {v['unique_mutations']!s:>4}  "
                  f"indels {v['indels_called']:>3} (del {v['deletions']}, ins {v['insertions']})  "
                  f"lengths {dict(sorted(lens.items()))}", flush=True)


RESULTS = {}


def checks():
    """The two cross-platform checks on the same six J_H4 mice (see the module docstring)."""
    out = json.loads(OUT.read_text())
    ont, san = RESULTS.get("jh4_ont"), RESULTS.get("jh4_sanger")
    if ont is None or san is None:
        sys.exit("checks need jh4_ont and jh4_sanger run in this call")
    ref = _ref()
    called = {(v["pos"], v["alt"]) for v in project.public(ont)["variants"]}
    pub_ont = project.public(ont)

    def subs(a):
        return {(r + 1, b) for r, b, q in a.pairs
                if r is not None and b in "ACGT" and b != ref[r] and q >= 20 and 101 <= r + 1 <= 629}

    co, link = {}, []
    for m in ("WT2", "WT3", "WT4"):
        reads, clones = ont["_details"][m]["_aligned"], san["_details"][m]["_aligned"]
        n = len(clones)
        csubs = [subs(c) for c in clones]
        hi_s = [v for v in pub_ont["variants"] if v["sample"] == m and (v.get("frequency") or 0) >= 0.02
                and 101 <= v["pos"] <= 629]
        hi_i = [i for i in out["jh4_ont"]["samples"][m]["indels"] if i["frequency"] >= 0.02]
        co[m] = {"sanger_clones": n,
                 "substitutions_ge_2pct": len(hi_s),
                 "substitutions_expected_in_sanger": round(sum(1 - (1 - v["frequency"]) ** n for v in hi_s), 1),
                 "substitutions_seen_in_sanger": sum(any((v["pos"], v["alt"]) in c for c in csubs) for v in hi_s),
                 "indels_ge_2pct": len(hi_i),
                 "indels_expected_in_sanger": round(sum(1 - (1 - i["frequency"]) ** n for i in hi_i), 1),
                 "indels_seen_in_sanger": sum(any(_carries(c, i) for c in clones) for i in hi_i)}
        for i in out["jh4_ont"]["samples"][m]["indels"]:
            carriers = [a for a in reads if _carries(a, i)]
            cnt = Counter(x for a in carriers for x in subs(a) if x in called)
            lineage = {x for x, c in cnt.items() if c >= 0.8 * len(carriers)}
            row = {"mouse": m, "indel": f"{i['kind']}{i['length']}@{i['pos']}", "reads": len(carriers),
                   "lineage_substitutions": len(lineage)}
            if len(lineage) >= 3:
                lin_r = [a for a in reads if len(lineage & subs(a)) >= 0.8 * len(lineage)]
                lin_c = [c for c, cs in zip(clones, csubs) if len(lineage & cs) >= 0.8 * len(lineage)]
                row.update(lineage_reads=len(lin_r), lineage_reads_with_indel=sum(_carries(a, i) for a in lin_r),
                           sanger_clones_of_lineage=len(lin_c),
                           sanger_clones_of_lineage_with_indel=sum(_carries(c, i) for c in lin_c))
            link.append(row)
    out["checks"] = {"cooccurrence": co, "linkage": link}
    OUT.write_text(json.dumps(out, indent=1, default=str))
    for m, v in co.items():
        print(f"   {m}: substitutions seen in Sanger {v['substitutions_seen_in_sanger']} of "
              f"{v['substitutions_expected_in_sanger']} expected; indels {v['indels_seen_in_sanger']} of "
              f"{v['indels_expected_in_sanger']}")
    for r in link:
        print("  ", r)


if __name__ == "__main__":
    args = set(sys.argv[1:])
    main(args or None)
    if not args or {"jh4_ont", "jh4_sanger"} <= args:
        checks()
