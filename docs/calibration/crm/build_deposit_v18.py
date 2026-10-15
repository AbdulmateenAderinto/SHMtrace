"""Data deposit Ver18: Ver17 with Call_anatomy/figS8_numbers.json renamed figS2_numbers.json.

    python docs/calibration/crm/build_deposit_v18.py

The supplementary figures were renumbered in order of first citation, and reported quality against measured
error became Figure S2. Copies ONT_SHM_DataDeposit_Ver17_2026-10-15 to ONT_SHM_DataDeposit_Ver18_<date>
(Ver17 is left untouched), renames the file (its content is unchanged), updates the CRM_and_STAR README's
figure numbers, regenerates FILE_MANIFEST.csv and SHA256SUMS.txt, and writes VER18_NOTES.md and
VER18_ADDITIONS_MANIFEST.csv.
"""
import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_deposit_v12 import SUB, sha256, write_manifests  # noqa: E402

SRC = SUB / "ONT_SHM_DataDeposit_Ver17_2026-10-15"
DST = SUB / f"ONT_SHM_DataDeposit_Ver18_{date.today().isoformat()}"


def main():
    if DST.exists():
        raise SystemExit(f"{DST} exists; remove it or change the date")
    shutil.copytree(SRC, DST, ignore=shutil.ignore_patterns(".DS_Store", "* 2.md", "* 2"))
    raw = DST / "Raw_Data"
    base = raw / "05_SHMtrace_Analysis" / "CRM_and_STAR"
    old, new = base / "Call_anatomy/figS8_numbers.json", base / "Call_anatomy/figS2_numbers.json"
    old.rename(new)
    readme = base / "README.md"
    t = readme.read_text()
    a = ("| `Call_anatomy/intuition_*.json` | per locus: every call with its allele fraction, each AID-active "
         "library's tested candidates and FDR threshold, the controls' decoy p values, per-position mismatch rates, "
         "reported quality against measured error, two annotated sites, and (J_H4) per-read mismatch rasters | "
         "Figures 2, S7 and S8 |")
    b = "| `Call_anatomy/fig2_numbers.json`, `figS8_numbers.json` | the summary numbers quoted in the text | Results, Discussion |"
    if t.count(a) != 1 or t.count(b) != 1:
        raise SystemExit("the Ver17 README rows were not found as expected")
    t = t.replace(a, a.replace("Figures 2, S7 and S8", "Figures 2, S2 and S6"))
    t = t.replace(b, b.replace("figS8_numbers.json", "figS2_numbers.json"))
    t += """
## Changed in Ver18

The supplementary figures were renumbered in order of first citation in the article: reported quality against
measured error is now Figure S2 (was S8) and the background along the amplicon is Figure S6 (was S7).
`Call_anatomy/figS8_numbers.json` is renamed `figS2_numbers.json`; its content is unchanged.
"""
    readme.write_text(t)
    man = raw / "04_Manifests_and_Integrity"
    notes = man / "VER18_NOTES.md"
    notes.write_text(f"""# Ver18 notes ({date.today().isoformat()})

- Renames `05_SHMtrace_Analysis/CRM_and_STAR/Call_anatomy/figS8_numbers.json` to `figS2_numbers.json` (content
  unchanged), after the article's supplementary figures were renumbered in order of first citation, and updates
  the figure numbers in that folder's README.
- No read file changes; `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over `Raw_Data/`.
""")
    with open(man / "VER18_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in (new, notes, readme):
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
    top = DST / "README.md"
    t = top.read_text()
    t = t.replace("## What changed from Ver16", f"""## What changed from Ver17 ({date.today().isoformat()})

- Renames `Call_anatomy/figS8_numbers.json` to `figS2_numbers.json` after the supplementary figures were
  renumbered; see `04_Manifests_and_Integrity/VER18_NOTES.md`. No read file changes.

## What changed from Ver16""", 1)
    top.write_text(t)
    write_manifests(raw)
    print(DST, "renamed", old.name, "->", new.name)


if __name__ == "__main__":
    main()
