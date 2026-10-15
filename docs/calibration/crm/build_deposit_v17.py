"""Data deposit Ver17: Ver16 plus the data behind Figures 2, S7 and S8.

    python docs/calibration/crm/build_deposit_v15.py

Copies ONT_SHM_DataDeposit_Ver16_2026-10-15 to ONT_SHM_DataDeposit_Ver17_<date> (Ver16 is left untouched) and
adds Raw_Data/05_SHMtrace_Analysis/CRM_and_STAR/Call_anatomy/ (intuition_data.py outputs), then regenerates
FILE_MANIFEST.csv and SHA256SUMS.txt and writes VER17_NOTES.md and VER17_ADDITIONS_MANIFEST.csv.
"""
import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_deposit_v12 import SUB, sha256, write_manifests  # noqa: E402

SRC = SUB / "ONT_SHM_DataDeposit_Ver16_2026-10-15"
DST = SUB / f"ONT_SHM_DataDeposit_Ver17_{date.today().isoformat()}"
OUT = Path(__file__).parent / "out"
ADD = {"Call_anatomy": [OUT / f"intuition_{k}.json" for k in ("jh4", "n114", "su_lps", "pim1")]
                        + [OUT / "fig2_numbers.json", OUT / "figS8_numbers.json"]}
README_ROWS = """
## Added in Ver17

| Folder / file | Contents | Used for |
|---|---|---|
| `Call_anatomy/intuition_*.json` | per locus: every call with its allele fraction, each AID-active library's tested candidates and FDR threshold, the controls' decoy p values, per-position mismatch rates, reported quality against measured error, two annotated sites, and (J_H4) per-read mismatch rasters | Figures 2, S7 and S8 |
| `Call_anatomy/fig2_numbers.json`, `figS8_numbers.json` | the summary numbers quoted in the text | Results, Discussion |
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
    notes = man / "VER17_NOTES.md"
    notes.write_text(f"""# Ver17 notes ({date.today().isoformat()})

- Adds the data behind Figures 2, S7 and S8 under `05_SHMtrace_Analysis/CRM_and_STAR/Call_anatomy/` (see its
  README, "Added in Ver17").
- No read file changes; `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over `Raw_Data/`.
""")
    with open(man / "VER17_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in [*added, notes, readme]:
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
    top = DST / "README.md"
    t = top.read_text()
    t = t.replace("## What changed from Ver15", f"""## What changed from Ver16 ({date.today().isoformat()})

- Adds the data behind Figures 2, S7 and S8 under `05_SHMtrace_Analysis/CRM_and_STAR/Call_anatomy/`; see
  `04_Manifests_and_Integrity/VER17_NOTES.md`. No read file changes.

## What changed from Ver15""", 1)
    top.write_text(t)
    write_manifests(raw)
    print(DST, f"{len(added)} files added")


if __name__ == "__main__":
    main()
