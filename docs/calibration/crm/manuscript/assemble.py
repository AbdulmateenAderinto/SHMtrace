"""Assemble the Cell Reports Methods and STAR Protocols submission folders beside the earlier packages.

    python docs/calibration/crm/manuscript/assemble.py
Run build_crm.py, build_star.py and build_letters.py (into crm/build) first. Writes, in submission_v4:
  SHMtrace_CellReportsMethods_<date>/   manuscript (docx, PDF), highlights/eTOC/motivation, cover letter,
                                         figures 1-6 (PDF, 600/1,000-dpi TIFF), graphical abstract,
                                         Document S1, Table S2, README
  SHMtrace_STARProtocols_<date>/        protocol (docx, PDF), cover letter, figures 1-4, graphical abstract
"""
import shutil
import subprocess
from datetime import date
from pathlib import Path

CRM = Path(__file__).resolve().parents[1]
BUILD = CRM / "build"
FIG = CRM / "figures"
SFIG = CRM / "star/figures"
SUB = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4")
SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
STAMP = date.today().isoformat()


def pdf(docx, outdir):
    subprocess.run([SOFFICE, "--headless", "--convert-to", "pdf", str(docx), "--outdir", str(outdir)],
                   check=True, capture_output=True)


def main():
    a = SUB / f"SHMtrace_CellReportsMethods_{STAMP}"
    s = SUB / f"SHMtrace_STARProtocols_{STAMP}"
    for d in (a, s):
        if d.exists():
            shutil.rmtree(d)
    for sub in ("1_Manuscript", "2_Figures", "3_Graphical_Abstract", "4_Supplemental", "5_Cover_Letter"):
        (a / sub).mkdir(parents=True)
    for f in ("SHMtrace_CellReportsMethods.docx", "Highlights_eTOC_Motivation.docx"):
        shutil.copy(BUILD / f, a / "1_Manuscript")
        pdf(a / "1_Manuscript" / f, a / "1_Manuscript")
    shutil.copy(BUILD / "CoverLetter_CellReportsMethods.docx", a / "5_Cover_Letter")
    pdf(a / "5_Cover_Letter/CoverLetter_CellReportsMethods.docx", a / "5_Cover_Letter")
    shutil.copy(FIG / "Figure1.pdf", a / "2_Figures/Figure1.pdf")
    shutil.copy(FIG / "Figure1_1000dpi.tif", a / "2_Figures/Figure1.tif")
    for i in range(2, 7):
        for ext in ("pdf", "tif"):
            shutil.copy(FIG / f"Figure{i}.{ext}", a / f"2_Figures/Figure{i}.{ext}")
    for ext, src in (("pdf", "GraphicalAbstract.pdf"), ("tif", "GraphicalAbstract_1000dpi.tif"),
                     ("png", "GraphicalAbstract_preview.png")):
        shutil.copy(FIG / src, a / f"3_Graphical_Abstract/GraphicalAbstract.{ext}")
    shutil.copy(BUILD / "Document_S1.docx", a / "4_Supplemental")
    pdf(a / "4_Supplemental/Document_S1.docx", a / "4_Supplemental")
    shutil.copy(CRM / "out/Table_S2_benchmark.xlsx", a / "4_Supplemental/Table_S2.xlsx")
    (a / "README.md").write_text(README_A)

    for sub in ("1_Manuscript", "2_Figures", "3_Graphical_Abstract", "4_Cover_Letter"):
        (s / sub).mkdir(parents=True)
    shutil.copy(BUILD / "SHMtrace_STARProtocols.docx", s / "1_Manuscript")
    pdf(s / "1_Manuscript/SHMtrace_STARProtocols.docx", s / "1_Manuscript")
    shutil.copy(BUILD / "CoverLetter_STARProtocols.docx", s / "4_Cover_Letter")
    pdf(s / "4_Cover_Letter/CoverLetter_STARProtocols.docx", s / "4_Cover_Letter")
    shutil.copy(SFIG / "STAR_Figure1.pdf", s / "2_Figures/Figure1.pdf")
    shutil.copy(SFIG / "STAR_Figure1_1000dpi.tif", s / "2_Figures/Figure1.tif")
    for i in range(2, 5):
        for ext in ("pdf", "tif"):
            shutil.copy(SFIG / f"STAR_Figure{i}.{ext}", s / f"2_Figures/Figure{i}.{ext}")
    shutil.copy(SFIG / "STAR_Figure1.pdf", s / "3_Graphical_Abstract/GraphicalAbstract.pdf")
    shutil.copy(SFIG / "STAR_Figure1_preview.png", s / "3_Graphical_Abstract/GraphicalAbstract.png")
    (s / "README.md").write_text(README_S)
    print(a)
    print(s)


README_A = f"""# Cell Reports Methods submission ({STAMP})

| Folder | Contents |
|---|---|
| `1_Manuscript/` | Article (docx, PDF) in journal order, figures embedded for review; highlights, eTOC blurb and motivation as the separate file the journal asks for |
| `2_Figures/` | Figures 1–6 as vector PDF and TIFF (Figure 1: SciCanvas, 1,000 dpi; Figures 2–6: 600 dpi) |
| `3_Graphical_Abstract/` | Graphical abstract (PDF, 1,000-dpi TIFF, PNG preview) |
| `4_Supplemental/` | Document S1 (Figures S1–S6, Tables S1 and S3), Table S2 (Excel) |
| `5_Cover_Letter/` | Cover letter draft |

Before submission: fill the bracketed items ([supplier], [Cat#], [RRID], accessions, DOIs, primer
sequences), confirm the author list and contacts, and attach the declaration of interests form.
Add the second library of N114 empty-vector subclone C3 (75K3B3_2_C3_Empty_March_2024_2.fastq, from
Hybridoma manuscript/75K3B3_fastq) to the data deposit: the N114 analyses use it in place of the first.
Built by docs/calibration/crm/manuscript/ in the SHMtrace repository (README in docs/calibration/crm/).
"""

README_S = f"""# STAR Protocols submission ({STAMP})

| Folder | Contents |
|---|---|
| `1_Manuscript/` | Protocol (docx, PDF): summary, highlights, before you begin, key resources table, materials and equipment, 26 numbered steps with timing, expected outcomes, quantification, limitations, troubleshooting, resource availability |
| `2_Figures/` | Figures 1–4 (PDF, TIFF) |
| `3_Graphical_Abstract/` | Graphical abstract (the protocol overview) |
| `4_Cover_Letter/` | Cover letter draft |

Before submission: fill the bracketed laboratory details ([to confirm]: buffers, antibody dilutions,
primer sequences, catalog numbers, provider requirements), the technical contact's email, and the
reference to the companion article once accepted.
"""

if __name__ == "__main__":
    main()
