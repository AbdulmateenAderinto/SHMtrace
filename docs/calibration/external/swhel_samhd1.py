"""A published Sanger study from another laboratory, re-analysed automatically.

    python docs/calibration/external/swhel_samhd1.py <genbank-clone-dir>

Thientosapol et al., PNAS 115:4921-4926 (2018), "SAMHD1 enhances immunoglobulin hypermutation by
promoting transversion mutation" (Jolly laboratory, Centenary Institute, Sydney). The study
hypermutates a SW_HEL VH10 knock-in V region in germinal-centre B cells and reports that
inactivating Samhd1 raises transitions at C:G and lowers transversions at A:T and at C:G.

Its 1,795 Sanger clones are deposited in GenBank (accessions MH198*-MH200*), one record per clone,
each naming its genotype, its host animal and its well. Nothing else about the study is used here:
the germline is the iterated majority consensus of the clones themselves, each clone is assigned to
the animal its own record names, and the groups are the genotypes its own records state.

Unlike the two studies of Figures S8 and S9 this one deposits sequences rather than chromatograms,
so it exercises the calling, the spectrum and the group tests but not the chromatogram reading. It
has no Aicda-/- group, so there is no control-calibrated background; Sanger error is about 2 in 10^5
bases (measured in the AID-null clones of both other studies), some 500-fold below the mutation
frequency here, so it cannot account for the differences tested.

Writes out/swhel_samhd1.json and a sample sheet beside the clones.
"""
import collections
import csv
import json
import os
import re
import sys
import time
from pathlib import Path

CAL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CAL.parents[1] / "src"))
from shmtrace import project, samplesheet  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
# the genotype each record states, mapped to the comparison it belongs to; the UGI arm is a
# separate experiment (uracil glycosylase inhibited), and Msh2 is a different gene
GROUP = {("WT-SWHEL+GFP", "Samhd1-wt or -het"): "Samhd1_WT",
         ("SAMHD1wt/25+GFP-SWHEL", "Samhd1-wt or -het"): "Samhd1_WT",
         ("SAMHD1WT/25+GFP-SWHEL", "Samhd1-wt or -het"): "Samhd1_WT",
         ("Samhd125/25-SWhel+GFP", "Samhd1-ko/ko"): "Samhd1_KO",
         ("Samhd141/41-SWhel+GFP", "Samhd1-ko/ko"): "Samhd1_KO",
         ("SAMHD1ko-SWhel+ugi", "SAMHD1-ko/ko"): "Samhd1_KO_UGI",
         ("MSH2ko-SWhel+GFP", "Msh2-ko/ko"): "Msh2_KO",
         ("MkoSko-SWhel+GFP", "Msh2-ko/ko, Samhd1-ko/ko"): "Msh2_Samhd1_KO"}
TS_GC, TV_GC = ("C>T", "G>A"), ("C>A", "C>G", "G>C", "G>T")
TS_AT, TV_AT = ("A>G", "T>C"), ("A>T", "A>C", "T>A", "T>G")


def records(clone_dir):
    """Every deposited clone as {accession, clone name, genotype, group, animal, well}."""
    out = []
    for f in sorted(os.listdir(clone_dir)):
        if not f.endswith(".gb"):
            continue
        t = open(os.path.join(clone_dir, f)).read()
        clone = re.search(r'/clone="([^"]+)"', t)
        note = re.search(r'/note="([^"]*)"', t, re.S)
        if not clone:
            continue
        name = clone.group(1)
        geno = ""
        if note:
            m = re.search(r"genotype:\s*([^\";\n]+)", note.group(1))
            geno = m.group(1).strip() if m else ""
        parts = name.split("_")
        animal = next((p for p in parts if p.startswith("HOST")), "")
        well = next((p for p in parts if p.startswith("SEQ")), "")
        grp = GROUP.get((parts[0], geno))
        out.append({"file": os.path.join(clone_dir, f), "accession": f[:-3], "clone": name,
                    "arm": parts[0], "genotype": geno, "group": grp,
                    # an animal id repeats between arms, so the sample is the arm's animal
                    "sample": f"{parts[0]}:{animal}", "well": well})
    return out


def write_sheet(recs, path):
    cols = ["file", "sample", "group", "control", "clone_id", "order", "well", "label", "exclude",
            "aliases", "note"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        for r in recs:
            w.writerow({"file": r["file"], "sample": r["sample"], "group": r["group"], "control": "",
                        "clone_id": r["accession"], "order": "PNAS2018", "well": r["well"],
                        "label": r["clone"], "exclude": "" if r["group"] else "not part of a tested arm",
                        "aliases": "", "note": f"genotype {r['genotype']}"})
    return path


def share(sp, keys):
    tot = sum(sp.values())
    return 100 * sum(sp.get(k, 0) for k in keys) / tot if tot else None


def main(clone_dir):
    recs = records(clone_dir)
    print(f"{len(recs)} deposited clones; arms " + ", ".join(
        f"{k}={v}" for k, v in sorted(collections.Counter(r["group"] for r in recs).items()) if k))
    sheet = write_sheet(recs, os.path.join(os.path.dirname(clone_dir), "swhel_samplesheet.csv"))
    ref = os.path.join(os.path.dirname(clone_dir), "swhel_vh10_consensus.fasta")
    t0 = time.time()
    sh = samplesheet.load(sheet)
    proj = sh.project(ref, str(HERE / "swhel_vh.yaml"))
    res = project.run_project(proj)
    seconds = time.time() - t0
    pub = project.public(res)
    per = {}
    for s in pub["samples"]:
        m = s["sample"]
        sp = pub["spectrum_sample"].get(m, {})
        per[m] = {"group": s["group"], "clones": s["reads"], "unique_mutations": s["unique_mutations"],
                  "bases": s["frequency_bases"],
                  "freq_x1e4": 1e4 * s["unique_mutations"] / s["frequency_bases"] if s["frequency_bases"] else 0,
                  "spectrum": sp,
                  "ts_gc": share(sp, TS_GC), "tv_gc": share(sp, TV_GC),
                  "ts_at": share(sp, TS_AT), "tv_at": share(sp, TV_AT),
                  "pct_GC_in_WRC_GYW": s.get("pct_GC_in_WRC_GYW"), "pct_AT_in_WA_TW": s.get("pct_AT_in_WA_TW")}
    OUT.mkdir(exist_ok=True)
    json.dump({"seconds": seconds, "clones_analysed": pub["qc"]["clones_analysed"],
               "excluded": collections.Counter(e["reason"] for e in pub["qc"]["excluded"]),
               "qc_flags": collections.Counter(f["check"] for f in pub["qc"]["flags"]),
               "per_mouse": per, "project": pub}, (OUT / "swhel_samhd1.json").open("w"), default=str)
    print(f"{pub['qc']['clones_analysed']} clones analysed in {seconds:.0f} s")
    for g in ("Samhd1_WT", "Samhd1_KO", "Samhd1_KO_UGI", "Msh2_KO", "Msh2_Samhd1_KO"):
        ms = [m for m in per if per[m]["group"] == g]
        if not ms:
            continue
        f = [per[m]["freq_x1e4"] for m in ms]
        print(f"  {g:16} {len(ms)} mice, {sum(per[m]['clones'] for m in ms):4} clones, "
              f"freq {sum(f)/len(f):6.1f} x1e-4")


if __name__ == "__main__":
    main(sys.argv[1])
