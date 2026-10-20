"""Bundle the manuscript, supplement and figures for an outside reviewer.

    python docs/calibration/crm/make_review_package.py   -> build/SHMtrace_review_package(.zip)

Run after manuscript/build_crm.py. The Word files are converted to PDF here (LibreOffice), Table S4's image is
rendered from the supplement PDF itself, and every standalone figure image is copied from crm/figures, the same
files the supplement embeds, so the images a reviewer sees cannot drift from the documents (the first review
found a stale Table S4 image and a stale Figure S10 image). The brief is review_brief.md.
"""
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import docx

CRM = Path(__file__).resolve().parent
BUILD, FIGS = CRM / "build", CRM / "figures"
PKG = BUILD / "SHMtrace_review_package"
SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
FIGURES = ["GraphicalAbstract_preview", "Figure1_preview", "Figure2", "Figure3", "Figure4", "Figure5", "Figure6",
           "Figure7", "FigureS1_workflow_preview", "FigureS2", "FigureS3", "FigureS5", "FigureS6",
           "FigureS7_bookkeeping", "FigureS8_sanger_cohort2", "FigureS9_sanger_mkln", "FigureS10_swhel",
           "FigureS12_platforms", "FigureS13_overview", "FigureS14_engine", "FigureS15",
           "FigureS16_switch_pacbio", "FigureS17_public_smu", "FigureS18_engineered"]


def run(*cmd, cwd=None):
    subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    docs = ["SHMtrace_CellReportsMethods.docx", "Document_S1.docx"]
    for d in docs:
        if not (BUILD / d).exists():
            sys.exit(f"{d} missing: run manuscript/build_crm.py first")
    run(SOFFICE, "--headless", "--convert-to", "pdf", *docs, "--outdir", str(BUILD), cwd=BUILD)
    if PKG.exists():
        shutil.rmtree(PKG)
    (PKG / "figures").mkdir(parents=True)
    shutil.copy(CRM / "review_brief.md", PKG / "00_REVIEW_BRIEF.md")
    for src, name in (("SHMtrace_CellReportsMethods", "01_Manuscript"), ("Document_S1", "02_Supplement_Document_S1")):
        shutil.copy(BUILD / f"{src}.pdf", PKG / f"{name}.pdf")
        run("pdftotext", "-layout", str(PKG / f"{name}.pdf"), str(PKG / f"{name}.txt"))
    front = docx.Document(BUILD / "Highlights_eTOC_Motivation.docx")
    (PKG / "03_Highlights_eTOC_Motivation.txt").write_text("\n".join(p.text for p in front.paragraphs))
    shutil.copy(CRM / "out/Table_S1_benchmark.xlsx", PKG / "04_Table_S1_benchmark.xlsx")
    for f in FIGURES:
        shutil.copy(FIGS / f"{f}.png", PKG / "figures" / f"{f.replace('_preview', '')}.png")
    # Table S4 exists only as a table in the supplement: render its page rather than keep a separate image
    page = next(n for n in range(1, 60)
                if subprocess.run(["pdftotext", "-f", str(n), "-l", str(n), str(PKG / "02_Supplement_Document_S1.pdf"),
                                   "-"], capture_output=True, text=True).stdout.lstrip().startswith("Table S4"))
    run("pdftoppm", "-f", str(page), "-l", str(page), "-r", "200", "-png", "-singlefile",
        str(PKG / "02_Supplement_Document_S1.pdf"), str(PKG / "figures" / "TableS4"))
    shutil.copy(PKG / "figures" / "TableS4.png", BUILD / "TableS4_preview.png")
    zp = BUILD / "SHMtrace_review_package.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(PKG.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(BUILD))
    print(f"{zp} ({zp.stat().st_size / 1e6:.1f} MB); Table S4 from supplement page {page}")


if __name__ == "__main__":
    main()
