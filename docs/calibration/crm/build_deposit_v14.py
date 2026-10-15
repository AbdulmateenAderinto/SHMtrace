"""Data deposit Ver14: Ver13 plus the other-caller runs added after Ver13.

    python docs/calibration/crm/build_deposit_v14.py

Copies ONT_SHM_DataDeposit_Ver13_2026-10-15 to ONT_SHM_DataDeposit_Ver14_<date> (Ver13 is left untouched) and
adds, under Raw_Data/05_SHMtrace_Analysis/CRM_and_STAR/Other_callers/:
  ClairS and Clair3 on the Smu and Pim1 spike-ins (Ver13 had J_H4 and N114 only)
  DeepVariant and DeepSomatic on the spike-ins of all four loci
  DeepSomatic on clean mutations planted at 5, 10 and 20% of reads, and on the real wild-type J_H4 libraries
and replaces Benchmark_real_libraries/Table_S2_benchmark.xlsx (now with all four callers), then regenerates
FILE_MANIFEST.csv and SHA256SUMS.txt and writes VER14_NOTES.md and VER14_ADDITIONS_MANIFEST.csv.
"""
import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_deposit_v12 import SUB, sha256, write_manifests  # noqa: E402

SRC = SUB / "ONT_SHM_DataDeposit_Ver13_2026-10-15"
DST = SUB / f"ONT_SHM_DataDeposit_Ver14_{date.today().isoformat()}"
OUT = Path(__file__).parent / "out"
ADD = {
    "Other_callers": [OUT / f"spikein_{k}_{c}.json" for k in ("su_lps", "pim1") for c in ("clairs", "clair3")]
                     + [OUT / f"spikein_{k}_{c}.json" for k in ("jh4", "n114", "su_lps", "pim1")
                        for c in ("deepvariant", "deepsomatic")]
                     + [OUT / "deepsomatic_highaf_jh4.json", OUT / "deepsomatic_real_jh4.json"],
    "Benchmark_real_libraries": [OUT / "Table_S2_benchmark.xlsx"],
}
README_ROWS = """
## Added in Ver14

| Folder / file | Contents | Used for |
|---|---|---|
| `Other_callers/spikein_{su_lps,pim1}_{clairs,clair3}.json` | ClairS and Clair3 on the Sμ and *Pim1* planted mutations (Ver13 had J_H4 and N114) | Results; Table S2 |
| `Other_callers/spikein_*_deepvariant.json` | DeepVariant 1.9.0 (ONT_R104; defaults and candidates from 1% of reads) on the planted mutations of all four loci | Results; Table S2 |
| `Other_callers/spikein_*_deepsomatic.json` | DeepSomatic 1.9.0 (ONT; defaults and candidates from 0.5% of reads), the other controls as its matched normal | Results; Table S2 |
| `Other_callers/deepsomatic_highaf_jh4.json` | DeepSomatic on clean mutations planted at 5, 10 and 20% of reads | Results; Figure 2H |
| `Other_callers/deepsomatic_real_jh4.json` | DeepSomatic on the real wild-type J_H4 libraries: every candidate, its FILTER, QUAL, allele fraction and support | Results |
| `Benchmark_real_libraries/Table_S2_benchmark.xlsx` | regenerated, now with ClairS, Clair3, DeepVariant and DeepSomatic on all four loci | Table S2 |
"""


def main():
    if DST.exists():
        raise SystemExit(f"{DST} exists; remove it or change the date")
    missing = [p for ps in ADD.values() for p in ps if not p.exists()]
    if missing:
        raise SystemExit("missing sources:\n  " + "\n  ".join(map(str, missing)))
    shutil.copytree(SRC, DST, ignore=shutil.ignore_patterns(".DS_Store", "* 2.md", "* 2"))
    raw = DST / "Raw_Data"
    base = raw / "05_SHMtrace_Analysis" / "CRM_and_STAR"
    added = []
    for folder, sources in ADD.items():
        d = base / folder
        d.mkdir(exist_ok=True)
        for src in sources:
            shutil.copy2(src, d / src.name)
            added.append(d / src.name)
    readme = base / "README.md"
    readme.write_text(readme.read_text() + README_ROWS)
    man = raw / "04_Manifests_and_Integrity"
    notes = man / "VER14_NOTES.md"
    notes.write_text(f"""# Ver14 notes ({date.today().isoformat()})

- Adds the other-caller runs made for the Cell Reports Methods article after Ver13, under
  `05_SHMtrace_Analysis/CRM_and_STAR/Other_callers/` (see its README, "Added in Ver14"): ClairS and Clair3 on
  the Sμ and *Pim1* spike-ins, DeepVariant and DeepSomatic on the spike-ins of all four loci, and DeepSomatic on
  clean high-fraction plants and on the real wild-type J_H4 libraries. Table S2 is regenerated.
- No read file changes; `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over `Raw_Data/`.
""")
    with open(man / "VER14_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in [*added, notes, readme]:
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
    top = DST / "README.md"
    t = top.read_text()
    t = t.replace("## What changed from Ver12", f"""## What changed from Ver13 ({date.today().isoformat()})

- Adds ClairS and Clair3 on Sμ and *Pim1*, and DeepVariant and DeepSomatic on all four loci, under
  `05_SHMtrace_Analysis/CRM_and_STAR/Other_callers/`; see `04_Manifests_and_Integrity/VER14_NOTES.md`. No read
  file changes.

## What changed from Ver12""", 1)
    top.write_text(t)
    write_manifests(raw)
    print(DST, f"{len(added)} files added")


if __name__ == "__main__":
    main()
