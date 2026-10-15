"""Data deposit Ver15: Ver14 plus DeepSomatic on the real N114 libraries.

    python docs/calibration/crm/build_deposit_v15.py

Copies ONT_SHM_DataDeposit_Ver14_2026-10-15 to ONT_SHM_DataDeposit_Ver15_<date> (Ver14 is left untouched) and
adds Raw_Data/05_SHMtrace_Analysis/CRM_and_STAR/Other_callers/deepsomatic_real_n114.json, then regenerates
FILE_MANIFEST.csv and SHA256SUMS.txt and writes VER15_NOTES.md and VER15_ADDITIONS_MANIFEST.csv.
"""
import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_deposit_v12 import SUB, sha256, write_manifests  # noqa: E402

SRC = SUB / "ONT_SHM_DataDeposit_Ver14_2026-10-15"
DST = SUB / f"ONT_SHM_DataDeposit_Ver15_{date.today().isoformat()}"
OUT = Path(__file__).parent / "out"
ADD = {"Other_callers": [OUT / "deepsomatic_real_n114.json"]}
README_ROWS = """
## Added in Ver15

| Folder / file | Contents | Used for |
|---|---|---|
| `Other_callers/deepsomatic_real_n114.json` | DeepSomatic on the real pMX-*Aicda* N114 subclones, the empty-vector libraries pooled as its normal: every candidate, its FILTER, QUAL, allele fraction and support | Results |
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
    notes = man / "VER15_NOTES.md"
    notes.write_text(f"""# Ver15 notes ({date.today().isoformat()})

- Adds DeepSomatic on the real pMX-*Aicda* N114 subclones, under `05_SHMtrace_Analysis/CRM_and_STAR/Other_callers/`
  (see its README, "Added in Ver15").
- No read file changes; `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over `Raw_Data/`.
""")
    with open(man / "VER15_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in [*added, notes, readme]:
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
    top = DST / "README.md"
    t = top.read_text()
    t = t.replace("## What changed from Ver13", f"""## What changed from Ver14 ({date.today().isoformat()})

- Adds DeepSomatic on the real N114 libraries under `05_SHMtrace_Analysis/CRM_and_STAR/Other_callers/`; see
  `04_Manifests_and_Integrity/VER15_NOTES.md`. No read file changes.

## What changed from Ver13""", 1)
    top.write_text(t)
    write_manifests(raw)
    print(DST, f"{len(added)} files added")


if __name__ == "__main__":
    main()
