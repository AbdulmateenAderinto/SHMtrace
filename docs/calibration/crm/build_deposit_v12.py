"""Data deposit Ver12: Ver11 plus the second library of N114 empty-vector subclone C3.

    python docs/calibration/crm/build_deposit_v12.py                  # Ver11 + the second C3 library
    python docs/calibration/crm/build_deposit_v12.py --crm-analysis   # add the CRM/STAR analyses

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

    print(DST)
    print(f"added {REL}: {n:,} reads, lines divisible by 4: {div4}, sequence/quality lengths match: {lens}")
    _readme(DST, n)
    write_additions(raw)
    write_manifests(raw)


def write_additions(raw):
    """VER12_ADDITIONS_MANIFEST.csv, written in full before the files are hashed."""
    base = raw / "05_SHMtrace_Analysis" / "CRM_and_STAR"
    with open(raw / "04_Manifests_and_Integrity/VER12_ADDITIONS_MANIFEST.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for p in (raw / REL, raw / "04_Manifests_and_Integrity/VER12_NOTES.md"):
            w.writerow([p.relative_to(raw).as_posix(), p.stat().st_size, sha256(p)])
        if base.exists():
            w.writerow([base.relative_to(raw).as_posix() + "/ (folder)",
                        sum(p.stat().st_size for p in base.rglob("*") if p.is_file()), "see SHA256SUMS.txt"])


def write_manifests(raw):
    man = raw / "04_Manifests_and_Integrity"
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
    print(f"{len(files)} files hashed")


def _readme(DST, n):
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


OUT = Path(__file__).parent / "out"
FIGS = CAL / "figures"
IRCM_REF = Path.home() / "Nanoshm/public_data/PRJNA1481264/references"
# deposit folder -> sources (files or folders), for the Cell Reports Methods article and STAR protocol
CRM = {
    "JH4_529bp_window": [FIGS / "jh4_shmtrace"],
    "JH4_full_length_90-1230": [FIGS / "jh4_shmtrace_full"],
    "N114_V_region": [FIGS / "n114_shmtrace"],
    "N114_sensitivity_first_C3_library": [OUT / "sensitivity_c3_1"],
    "Depth_series": [FIGS / "depth_shmtrace"],
    "Sanger_vs_ONT": [FIGS / "sanger_vs_nanopore/recommended_manuscript"],
    "Mutation_trees": [FIGS / "mutation_trees"],
    "Clones_and_targeting": [OUT / "extras_jh4", OUT / "extras_jh4_full", OUT / "extras_n114",
                             OUT / "targeting_jh4.json", OUT / "targeting_ch12.json"],
    "Benchmark_real_libraries": [*sorted(OUT.glob("benchmark_real_*.json")), OUT / "Table_S2_benchmark.xlsx"],
    "Spike_ins": [OUT / f"spikein_{k}.json" for k in ("jh4", "n114", "su_lps", "pim1")]
                 + [OUT / "spikein_jh4_bh.json", OUT / "spikein_n114_bh.json", OUT / "spikein_n114_4ctl.json",
                    OUT / "single_read_spike.json"],
    "Simulations": [OUT / "lineage_sim_jh4", OUT / "quant_accuracy_jh4", OUT / "quant_accuracy_n114"],
    "Public_Illumina": [OUT / "public" / k for k in ("ldha", "ch12", "hmces", "hmces_k3")],
    "Public_nanopore_IRCM_references": [IRCM_REF / "Igh_Su.fasta", IRCM_REF / "Pim1.fasta"],
    "UMI_consensus_Karst": [CAL / "umi_karst/results.json"],
}

CRM_README = """# SHMtrace analyses for the Cell Reports Methods article and the STAR Protocols protocol

Every analysis behind those two manuscripts, as SHMtrace wrote it (recommended mode: control-error
model with a leave-one-control-out margin and the control-based, target-decoy FDR). The folders above
this one hold the JIM manuscript (Ver21) analyses, made with the earlier Benjamini-Hochberg threshold.
Reads: `../../02_JH4_ONT/Original_FASTQ/` and `../../03_N114_ONT/Original_FASTQ/`; for N114, empty-vector
subclone C3 is its second library (`75K3B3_2_C3_Empty_March_2024_2.fastq`).

| Folder | Contents | Used for |
|---|---|---|
| `JH4_529bp_window/` | J_H4 positions 101-629: `project.json`, report, tables, `table_s1.json` | Figures 3-5; Tables S1, S3 |
| `JH4_full_length_90-1230/` | J_H4 positions 90-1,230 | Figure 5; Figure S6 |
| `N114_V_region/` | N114 positions 21-635 (C3 = second library) | Figure 3; Figure S3; Table S1 |
| `N114_sensitivity_first_C3_library/` | the same with the first C3 library: main run, benchmark, spike-ins, depth, quantification | Results (sensitivity); Table S3 |
| `Depth_series/` | 30-1,000 reads per library, controls calibrated at each depth | Figure S2 |
| `Sanger_vs_ONT/` | Sanger and nanopore statistics over the 529-bp window | Figure 4; Table S3 |
| `Mutation_trees/` | junction clones and trees, four-gamete tests | Figure 5D; Figure S6 |
| `Clones_and_targeting/` | clone-weighted frequencies (`extras_*`), S5F targeting indices | Figures 5B-5C, 6C |
| `Benchmark_real_libraries/` | SHMtrace, the Geneious rule and LoFreq (three settings) on the same reads, six datasets; Table S2 | Figure 3; Table S2 |
| `Spike_ins/` | planted mutations per held-out control and caller; `_bh`: Benjamini-Hochberg with the earlier controls (J_H4; N114 with the first C3 library), and caller `recommended_bh` in `spikein_n114.json` with the second; `n114_4ctl`: all four empty-vector libraries; `single_read_spike.json` | Figure 2A-2C; Figure S5A |
| `Simulations/` | B cell lineage simulations and quantification accuracy by depth | Figure 2D-2F |
| `Public_Illumina/` | PRJNA973112 (Ldha), PRJNA918596 (CH12F3), PRJNA820976 (Ramos HMCES; `hmces_k3`: AID-KO clone 3 analysed as a sample) | Figure 6A-6D |
| `Public_nanopore_IRCM_references/` | consensus references built from PRJNA1481264 controls (Smu, Pim1) | Figures 2-3 |
| `UMI_consensus_Karst/` | UMI consensus error by reads per molecule (results only; the reads are Karst et al.'s) | Figure 6E |

Reproduce with `docs/calibration/crm/` of the SHMtrace repository (`rerun_A.sh`-`rerun_D.sh`,
`rerun_N114.sh`, then the figure scripts and `manuscript/`); every number in the text is computed from
these files by `manuscript/crm_numbers.py`.
"""


def add_crm(dst=None):
    """Add the CRM/STAR analyses to the deposit as 05_SHMtrace_Analysis/CRM_and_STAR/ and re-hash."""
    dst = Path(dst) if dst else DST
    raw = dst / "Raw_Data"
    base = raw / "05_SHMtrace_Analysis" / "CRM_and_STAR"
    if base.exists():
        shutil.rmtree(base)            # this builder's own folder, rebuilt from the sources
    base.mkdir(parents=True)
    added = []
    for folder, sources in CRM.items():
        d = base / folder
        d.mkdir()
        for src in sources:
            if not src.exists():
                raise SystemExit(f"missing source: {src}")
            target = d / src.name
            if src.is_dir():
                shutil.copytree(src, target, ignore=shutil.ignore_patterns(".DS_Store", "*.bam", "*.bai"))
            else:
                shutil.copy2(src, target)
            added.append(target)
    (base / "README.md").write_text(CRM_README)
    t = (raw / "05_SHMtrace_Analysis/README.md").read_text()
    if "CRM_and_STAR" not in t:
        t += ("\n`CRM_and_STAR/` holds the analyses of the Cell Reports Methods article and the STAR Protocols "
              "protocol (control-based FDR; N114 C3 = its second library); see its README.\n")
        (raw / "05_SHMtrace_Analysis/README.md").write_text(t)
    notes = raw / "04_Manifests_and_Integrity/VER12_NOTES.md"
    t = notes.read_text()
    if "CRM_and_STAR" not in t:
        t = t.replace("contents of `05_SHMtrace_Analysis/` are still that manuscript's analyses.",
                      "the other folders of `05_SHMtrace_Analysis/` are that manuscript's analyses.")
        t += ("- Adds `05_SHMtrace_Analysis/CRM_and_STAR/`: every analysis behind the Cell Reports Methods "
              "article and the STAR Protocols protocol (see its README).\n")
        notes.write_text(t)
    write_additions(raw)
    write_manifests(raw)
    print(f"{base}: {sum(1 for p in base.rglob('*') if p.is_file())} files")


if __name__ == "__main__":
    import sys
    if "--crm-analysis" in sys.argv:
        add_crm()
    else:
        main()
