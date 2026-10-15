"""Data deposit Ver13: Ver12 plus the analyses added to the Cell Reports Methods article after Ver12.

    python docs/calibration/crm/build_deposit_v13.py

Copies ONT_SHM_DataDeposit_Ver12_2026-10-13 to ONT_SHM_DataDeposit_Ver13_<date> (Ver12 is left untouched) and
adds, under Raw_Data/05_SHMtrace_Analysis/CRM_and_STAR/:
  Other_callers/        ClairS and Clair3 on the planted mutations, ClairS on clean high-fraction plants and on the
                        real wild-type J_H4 libraries
  Spike_ins/            sensitive mode (FDR 0.20), indels planted into nanopore controls, base-quality
                        recalibration (tested, off by default)
  Benchmark_real_libraries/Table_S2_benchmark.xlsx, regenerated with ClairS and Clair3
then regenerates FILE_MANIFEST.csv and SHA256SUMS.txt and writes VER13_NOTES.md and VER13_ADDITIONS_MANIFEST.csv.
"""
import csv
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_deposit_v12 import SUB, sha256, write_manifests  # noqa: E402

SRC = SUB / "ONT_SHM_DataDeposit_Ver12_2026-10-13"
DST = SUB / f"ONT_SHM_DataDeposit_Ver13_{date.today().isoformat()}"
OUT = Path(__file__).parent / "out"
ADD = {
    "Other_callers": [OUT / f"spikein_{k}_{c}.json" for k in ("jh4", "n114") for c in ("clairs", "clair3")]
                     + [OUT / "clairs_real_jh4.json", OUT / "clairs_highaf_jh4.json"],
    "Spike_ins": [OUT / f"spikein_{k}_{t}.json" for k in ("jh4", "n114") for t in ("sensitive", "recal", "recal2")]
                 + [OUT / f"indel_spike_{k}.json" for k in ("jh4", "n114")],
    "Benchmark_real_libraries": [OUT / "Table_S2_benchmark.xlsx"],
}
README_ROWS = """
## Added in Ver13

| Folder / file | Contents | Used for |
|---|---|---|
| `Other_callers/spikein_*_clairs.json` | ClairS (defaults and snv_min_af 0.005), the other controls as its matched normal, on the planted mutations of `Spike_ins/` | Results (ClairS and Clair3); Table S2 |
| `Other_callers/spikein_*_clair3.json` | Clair3 (defaults and snp_min_af 0.01) on the same planted mutations | Results; Table S2 |
| `Other_callers/clairs_highaf_jh4.json` | ClairS on clean mutations planted at 5, 10 and 20% of reads | Results |
| `Other_callers/clairs_real_jh4.json` | ClairS on the real wild-type J_H4 libraries: every candidate, its FILTER, QUAL, allele fraction and support | Results |
| `Spike_ins/spikein_*_sensitive.json` | the recommended mode at FDR 0.10 and 0.20, and the calibrated test uncorrected | Discussion (sensitive mode) |
| `Spike_ins/indel_spike_*.json` | insertions and deletions planted into held-out nanopore controls | Results (indels) |
| `Spike_ins/spikein_*_recal*.json` | base-quality recalibration from the controls: a remap (`_recal`) and the shift kept as an option (`_recal2`) | SHMtrace docs (RECAL.md) |
| `Benchmark_real_libraries/Table_S2_benchmark.xlsx` | regenerated, now with ClairS and Clair3 | Table S2 |
"""


def main():
    if DST.exists():
        raise SystemExit(f"{DST} exists; remove it or change the date")
    shutil.copytree(SRC, DST, ignore=shutil.ignore_patterns(".DS_Store", "* 2.md", "* 2"))
    raw = DST / "Raw_Data"
    base = raw / "05_SHMtrace_Analysis" / "CRM_and_STAR"
    added = []
    for folder, sources in ADD.items():
        d = base / folder
        d.mkdir(exist_ok=True)
        for src in sources:
            if not src.exists():
                raise SystemExit(f"missing source: {src}")
            shutil.copy2(src, d / src.name)
            added.append(d / src.name)
    readme = base / "README.md"
    readme.write_text(readme.read_text() + README_ROWS)
    man = raw / "04_Manifests_and_Integrity"
    notes = man / "VER13_NOTES.md"
    notes.write_text(f"""# Ver13 notes ({date.today().isoformat()})

- Adds the analyses made for the Cell Reports Methods article after Ver12, under
  `05_SHMtrace_Analysis/CRM_and_STAR/` (see its README, "Added in Ver13"): ClairS and Clair3 on the planted
  mutations and on the real J_H4 libraries, the sensitive-mode, indel and base-quality-recalibration spike-ins,
  and Table S2 regenerated with ClairS and Clair3.
- No read file changes; `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over `Raw_Data/`.
""")
    with open(man / "VER13_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in [*added, notes, readme]:
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
    top = DST / "README.md"
    t = top.read_text()
    t = t.replace("## What changed from Ver11", f"""## What changed from Ver12 ({date.today().isoformat()})

- Adds the later Cell Reports Methods analyses (other callers, sensitive mode, indels, recalibration) under
  `05_SHMtrace_Analysis/CRM_and_STAR/`; see `04_Manifests_and_Integrity/VER13_NOTES.md`. No read file changes.

## What changed from Ver11""", 1)
    top.write_text(t)
    write_manifests(raw)
    print(DST, f"{len(added)} files added")


if __name__ == "__main__":
    main()
