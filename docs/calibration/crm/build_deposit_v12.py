"""Data deposit Ver12: Ver11 plus the second library of N114 empty-vector subclone C3.

    python docs/calibration/crm/build_deposit_v12.py

Copies ONT_SHM_DataDeposit_Ver11_2026-10-13 to ONT_SHM_DataDeposit_Ver12_<date> (Ver11 is left
untouched; its sync-conflict duplicates "README 2.md" and "Raw_Data 2" are not copied), adds
75K3B3_2_C3_Empty_March_2024_2.fastq to Raw_Data/03_N114_ONT/Original_FASTQ/ from the hybridoma
manuscript folder, adds it to SAMPLE_MANIFEST.csv and FASTQ_VALIDATION.csv, regenerates FILE_MANIFEST.csv
and SHA256SUMS.txt over Raw_Data/, and writes VER12_NOTES.md and VER12_ADDITIONS_MANIFEST.csv.
"""
import csv
import hashlib
import json
import shutil
from datetime import date
from pathlib import Path

SUB = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4")
SRC = SUB / "ONT_SHM_DataDeposit_Ver11_2026-10-13"
DST = SUB / f"ONT_SHM_DataDeposit_Ver12_{date.today().isoformat()}"
HYB = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Hybridoma manuscript/75K3B3_fastq")
NEW = "75K3B3_2_C3_Empty_March_2024_2"
REL = f"03_N114_ONT/Original_FASTQ/{NEW}.fastq"
CAL = Path(__file__).resolve().parents[1]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fastq_check(p):
    lines = Path(p).read_text().splitlines()
    ok = all(len(lines[i + 1]) == len(lines[i + 3]) for i in range(0, len(lines) - 3, 4))
    return len(lines) // 4, len(lines) % 4 == 0, ok


def main():
    if DST.exists():
        raise SystemExit(f"{DST} exists; remove it or change the date")
    shutil.copytree(SRC, DST, ignore=lambda d, names: [n for n in names
                                                      if Path(d) == SRC and n in ("README 2.md", "Raw_Data 2")])
    raw = DST / "Raw_Data"
    man = raw / "04_Manifests_and_Integrity"
    shutil.copy2(HYB / f"{NEW}.fastq", raw / REL)
    n, div4, lens = fastq_check(raw / REL)

    # SAMPLE_MANIFEST: the new library, and what each C3 library is used for
    acct = json.load(open(CAL / "figures/n114_shmtrace/table_s1.json"))["C3"]
    rows = list(csv.reader(open(man / "SAMPLE_MANIFEST.csv")))
    out = []
    for r in rows:
        if len(r) > 3 and r[0] == "ONT N114" and r[3] == "75K3B3_1_C3_Empty_March_2024_1":
            r[-1] = ("First library of subclone C3: analysed in the JIM manuscript (Ver21); in the Cell Reports "
                     "Methods article, a sensitivity analysis only")
            out.append(r)
            out.append(["ONT N114", "pMX-empty", "C3 (second library)", NEW, REL, str(n), str(acct["usable"]),
                        str(acct["analysed"]), "",
                        "Second library of subclone C3: analysed in place of the first in the Cell Reports Methods "
                        "article (the first was noisier than the other controls)"])
        else:
            out.append(r)
    csv.writer(open(man / "SAMPLE_MANIFEST.csv", "w", newline="")).writerows(out)
    with open(man / "FASTQ_VALIDATION.csv", "a") as fh:      # the file's own format: path quoted
        fh.write(f'"{REL}",{n},{"yes" if div4 else "no"},{"yes" if lens else "no"}\n')

    notes = man / "VER12_NOTES.md"
    notes.write_text(f"""# Ver12 notes ({date.today().isoformat()})

- Adds `{REL}`, the second library of N114 empty-vector subclone C3 ({n:,} reads, copied unchanged
  from the hybridoma manuscript's `75K3B3_fastq/` folder). Both C3 libraries were prepared in March 2024.
- The Cell Reports Methods article and its STAR Protocols companion analyse this library as the C3 control;
  the first library (`75K3B3_1_C3_Empty_March_2024_1.fastq`) was noisier than the other controls and is
  reported there as a sensitivity analysis. The JIM manuscript (Ver21) analyses the first library, and the
  contents of `05_SHMtrace_Analysis/` are still that manuscript's analyses.
- `SAMPLE_MANIFEST.csv` gains a row for the new library; `FASTQ_VALIDATION.csv` gains its check.
- `FILE_MANIFEST.csv` and `SHA256SUMS.txt` are regenerated over every file under `Raw_Data/`.
- Ver11's `README 2.md` and `Raw_Data 2/` (sync-conflict duplicates) are not carried over.
""")

    files = sorted(p for p in raw.rglob("*") if p.is_file() and p.name != ".DS_Store"
                   and p.name not in ("SHA256SUMS.txt", "FILE_MANIFEST.csv"))
    digests = {p: sha256(p) for p in files}
    with open(man / "FILE_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["relative_path", "size_bytes", "sha256", "top_level_category"])
        for p in files:
            rel = p.relative_to(raw).as_posix()
            w.writerow([rel, p.stat().st_size, digests[p], rel.split("/")[0]])
    digests[man / "FILE_MANIFEST.csv"] = sha256(man / "FILE_MANIFEST.csv")
    with open(man / "SHA256SUMS.txt", "w") as fh:
        for p in sorted(digests):
            fh.write(f"{digests[p]}  {p.relative_to(raw).as_posix()}\n")
    with open(man / "VER12_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in (raw / REL, notes):
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])

    readme = DST / "README.md"
    t = readme.read_text().replace("Ver11 (2026-10-13)", f"Ver12 ({date.today().isoformat()})")
    t = t.replace("`VER11_NOTES.md`, `VER11_ADDITIONS_MANIFEST.csv`",
                  "`VER11_NOTES.md`, `VER12_NOTES.md`, `VER12_ADDITIONS_MANIFEST.csv`")
    t = t.replace("## What changed from Ver10", f"""## What changed from Ver11

- Adds the second library of N114 empty-vector subclone C3 (`{REL}`), which the Cell Reports
  Methods article analyses as the C3 control; see `04_Manifests_and_Integrity/VER12_NOTES.md`.
- Manifests and checksums regenerated.

## What changed from Ver10""")
    readme.write_text(t)
    print(DST)
    print(f"added {REL}: {n:,} reads, lines divisible by 4: {div4}, sequence/quality lengths match: {lens}")
    print(f"{len(files)} files hashed")


if __name__ == "__main__":
    main()
