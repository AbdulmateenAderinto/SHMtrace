"""Export the Sanger calls the manuscript figures read.

    python export_sanger_calls.py [SANGER_DIR]      (default ~/Nanoshm/sanger_ab1_original)

SANGER_DIR is a set organised by organize_jh4_ab1.py, one folder per mouse. The
manuscript uses the original classification (--original): clones where they were filed.
Writes sanger_calls/<mouse> corrected Sanger Annotations.csv (Geneious-shaped, positions
on the 529-bp reference) and sanger_calls/corrected_sanger.json.
"""
import csv
import json
import sys
from pathlib import Path

from nanoshm import project

ROOT = Path(__file__).resolve().parents[3]
REF = str(ROOT / "src/nanoshm/references/mouse_jh4_intron.fasta")
CON = [str(ROOT / "src/nanoshm/references/mouse_jh4_contaminant.fasta")]
MICE = {"Wt 2": "WT2", "Wt 3": "WT3", "Wt 4": "WT4",
        "Aicda−/− 1": "Aicda_A1", "Aicda−/− 2": "Aicda_A2", "Aicda−/− 3": "Aicda_A3"}


def main(src: Path):
    out_dir = Path(__file__).resolve().parent / "sanger_calls"
    out_dir.mkdir(exist_ok=True)
    r = project.run_project({
        "name": "sanger", "region": "jh4_intron", "platform": "sanger", "reference": REF,
        "contaminants": CON,
        "samples": [{"name": f, "path": str(src / f), "group": "WT" if f.startswith("WT") else "KO",
                     "control": not f.startswith("WT")} for f in MICE.values()]})
    out = {"source": str(src)}
    for lab, f in MICE.items():
        d = r["_details"][f]
        vs = [v for v in d["variants"] if 101 <= v.pos <= 629]
        out[lab] = {"clones": len(d["_aligned"]), "unique": len(vs),
                    "bases": sum(x["bases"] for x in d["per_read"]),
                    "loads": [x["mutations"] for x in d["per_read"]],
                    "variants": [(v.pos, v.ref, v.alt, v.supporting_reads) for v in vs]}
        with open(out_dir / f"{f} corrected Sanger Annotations.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Name", "Minimum", "Polymorphism Type", "Reference Nucleotide(s)",
                        "Variant Nucleotide(s)", "Variant Frequency", "Variant P-Value (approximate)",
                        "Coverage"])
            for v in vs:
                kind = "transition" if (v.ref + v.alt) in ("AG", "GA", "CT", "TC") else "transversion"
                w.writerow([f"{v.ref}->{v.alt}", v.pos - 100, f"SNP ({kind})", v.ref, v.alt,
                            f"{100 * v.supporting_reads / max(v.depth, 1):.1f}%", "0", v.depth])
        print(f"{lab:12} clones {out[lab]['clones']:2} unique {out[lab]['unique']:3} "
              f"bases {out[lab]['bases']} (= {out[lab]['bases'] / 529:.1f} full clones)")
    json.dump(out, open(out_dir / "corrected_sanger.json", "w"), indent=1)


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 else Path.home() / "Nanoshm/sanger_ab1_original")
