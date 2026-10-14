"""The three public Illumina datasets, re-run with the current SHMtrace for the methods paper.

    python docs/calibration/crm/public_illumina.py ldha ch12 hmces

  ldha   PRJNA973112 (GSE232658): splenic GC B cells, JH4 (positions 53-122); Wt and Ldha-KO against
         Aicda-null (AID-Cre/Cre) mice; 20,000 reads per mouse
  ch12   PRJNA918596: CH12F3 Ung-/-Msh2-/-Aicda-/- clone 120, VDJ amplicon; ectopic AID, the R50E and
         R174E/R178E mutants and empty vector (controls); 5,000 reads per replicate
  hmces  PRJNA820976: Ramos IGH VDJ; Wt and HMCES-KO clones against AID-KO clones; indels on;
         20,000 reads per clone
Settings are those of the original drivers (~/Nanoshm/public_data/run_prjna*.py), with the 5-mer
targeting model. Writes out/public/<name>/ (project.json and the report).
"""
import json
import sys
from pathlib import Path

from shmtrace import project, report

PUB = Path.home() / "Nanoshm/public_data"
REPO = Path(__file__).resolve().parents[3]
REF = REPO / "src/shmtrace/references"
S5F = str(Path.home() / "Nanoshm/models/HH_S5F.csv")
OUT = Path(__file__).parent / "out/public"


def ldha():
    s = PUB / "PRJNA973112_sub"
    group = lambda m: "AID_KO" if "AIDKO" in m else ("WT" if "_WT" in m else "LDHA_KO")   # noqa: E731
    samples = [{"name": m.replace("JH4_", "").replace("JH_", ""), "path": f"{s}/{r}_1.fastq.gz+{s}/{r}_2.fastq.gz",
                "group": group(m), "control": "AIDKO" in m}
               for m, r in (line.split() for line in open(s / "mouse_runs.txt"))]
    return {"name": "PRJNA973112 JH4 (LDHA)", "region": "jh4_intron", "platform": "illumina",
            "reference": str(REF / "mouse_jh4_intron.fasta"), "contaminants": [str(REF / "mouse_jh4_contaminant.fasta")],
            # 5,000 reads per mouse: the 70-bp window is resequenced far beyond one read per molecule, and at
            # 20,000 the run no longer finished in reasonable time with the current code
            "depth": 5000, "targeting_model": S5F, "samples": samples,
            # clones off: the 70-bp window gives no use for them, and junction clustering of 20,000 reads per
            # mouse is quadratic
            "overrides": {"analysis_window": [53, 122], "end_trim_bp": 0, "min_align_len": 60, "exclude_indels": False,
                          "junction_ref_start": None}}


def ch12():
    s = PUB / "PRJNA918596_ch12/pooled"
    conds = ["AID", "R50E", "R174ER178E", "EV"]
    samples = [{"name": f"{c}_r{r}", "path": f"{s}/{c}_r{r}_1.fastq.gz+{s}/{c}_r{r}_2.fastq.gz",
                "group": c, "control": c == "EV"} for c in conds for r in "1234"]
    return {"name": "PRJNA918596 CH12F3", "region": "generic_amplicon", "platform": "illumina",
            "reference": str(PUB / "PRJNA918596_ch12_ref.fasta"), "depth": 5000, "targeting_model": S5F,
            "samples": samples, "overrides": {"analysis_window": [21, 249], "end_trim_bp": 0, "min_align_len": 200,
                                              "max_total_len": 400, "exclude_indels": False}}


def hmces():
    s = PUB / "PRJNA820976_sub"
    libs = {"AIDKO_1": "SRR18513363", "AIDKO_2": "SRR18513362", "AIDKO_3": "SRR18513361",
            "HMCES_1": "SRR18513314", "HMCES_2": "SRR18513313", "HMCES_3": "SRR18513312",
            "HMCES_4": "SRR18513310", "HMCES_5": "SRR18513365", "HMCES_6": "SRR18513364",
            "WT_1": "SRR18513320", "WT_2": "SRR18513319", "WT_3": "SRR18513318",
            "WT_4": "SRR18513317", "WT_5": "SRR18513316", "WT_6": "SRR18513315"}
    grp = {"AIDKO": "AID_KO", "HMCES": "HMCES_KO", "WT": "WT"}
    samples = [{"name": n, "path": f"{s}/{r}_1.fastq.gz+{s}/{r}_2.fastq.gz", "group": grp[n.split("_")[0]],
                "control": n.startswith("AIDKO")} for n, r in libs.items()]
    return {"name": "PRJNA820976 Ramos (HMCES)", "region": "generic_amplicon", "platform": "illumina",
            "reference": str(PUB / "PRJNA820976_ramos_ref.fasta"), "depth": 20000, "targeting_model": S5F,
            "samples": samples, "overrides": {"analysis_window": [26, 446], "end_trim_bp": 0, "min_align_len": 400,
                                              "max_total_len": 600, "exclude_indels": False}}


def hmces_k3():
    """Ramos with AID-KO clone 3 analyzed as a sample: it carries non-AID-like mutations, and as a
    decoy it raises the control-based threshold of every clone (SHMtrace names it in the notes)."""
    proj = hmces()
    proj["name"] += ", AIDKO_3 as a sample"
    for s in proj["samples"]:
        if s["name"] == "AIDKO_3":
            s["control"] = False
    return proj


if __name__ == "__main__":
    for key in sys.argv[1:]:
        proj = {"ldha": ldha, "ch12": ch12, "hmces": hmces, "hmces_k3": hmces_k3}[key]()
        proj["overrides"] = {**proj.get("overrides", {}), "jobs": 0}     # all but one core; same results
        res = project.run_project(proj)
        out = OUT / key
        out.mkdir(parents=True, exist_ok=True)
        report.write_all(res, str(out))
        json.dump(project.public(res), open(out / "project.json", "w"), default=str)
        for r in sorted(res["samples"], key=lambda r: (r["group"], r["sample"])):
            print(key, r["sample"], r["group"], r["reads"], r["unique_mutations"], r["mutation_frequency_x1e3"],
                  r.get("pct_CT_GA"), r.get("s5f_index"), r.get("indels_called"), flush=True)
