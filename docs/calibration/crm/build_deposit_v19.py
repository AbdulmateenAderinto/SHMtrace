"""Data deposit Ver19: Ver18 plus the ectopic-AID validation set (PRJNA170048, Figure 7F).

    python docs/calibration/crm/build_deposit_v15.py

Copies ONT_SHM_DataDeposit_Ver18_2026-10-15 to ONT_SHM_DataDeposit_Ver19_<date> (Ver18 is left untouched) and
adds Raw_Data/05_SHMtrace_Analysis/CRM_and_STAR/Ectopic_AID/ (the SHMtrace project and the original study's
own per-library counts), then regenerates
FILE_MANIFEST.csv and SHA256SUMS.txt and writes VER19_NOTES.md and VER19_ADDITIONS_MANIFEST.csv.
"""
import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_deposit_v12 import SUB, sha256, write_manifests  # noqa: E402

SRC = SUB / "ONT_SHM_DataDeposit_Ver18_2026-10-15"
DST = SUB / f"ONT_SHM_DataDeposit_Ver19_{date.today().isoformat()}"
OUT = Path(__file__).parent / "out"
ADD = {"Ectopic_AID": [OUT / "public/morange/project.json", OUT / "public/morange/authors_counts.json"]}
README_ROWS = """
## Added in Ver19

| Folder / file | Contents | Used for |
|---|---|---|
| `Ectopic_AID/project.json` | SHMtrace on PRJNA170048 (Perez-Duran et al. 2012): AID-ER vs catalytically dead AID-E58Q-ER in NIH-3T3 cells on an mOrange reporter, in three repair contexts, 18 libraries | Results; Figure 7F |
| `Ectopic_AID/authors_counts.json` | the same libraries counted as the original study did (every mismatch at quality > 20, no control subtraction), from GSE39108 Supplementary Table I | Results; Figure 7F |
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
    notes = man / "VER19_NOTES.md"
    notes.write_text(f"""# Ver19 notes ({date.today().isoformat()})

- Adds the ectopic-AID validation set under `05_SHMtrace_Analysis/CRM_and_STAR/Ectopic_AID/` (see its README,
  "Added in Ver19").
- No read file changes; `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over `Raw_Data/`.
""")
    with open(man / "VER19_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in [*added, notes, readme]:
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
    top = DST / "README.md"
    t = top.read_text()
    t = t.replace("## What changed from Ver17", f"""## What changed from Ver18 ({date.today().isoformat()})

- Adds the ectopic-AID validation set under `05_SHMtrace_Analysis/CRM_and_STAR/Ectopic_AID/`; see
  `04_Manifests_and_Integrity/VER19_NOTES.md`. No read file changes.

## What changed from Ver17""", 1)
    top.write_text(t)
    write_manifests(raw)
    print(DST, f"{len(added)} files added")


if __name__ == "__main__":
    main()
