"""Build the STAR Protocols article.

    python docs/calibration/crm/manuscript/build_star.py OUTDIR

Writes OUTDIR/SHMtrace_STARProtocols.docx in the journal's order: title, authors, summary, highlights,
graphical abstract note, before you begin, key resources table, materials and equipment, step-by-step
method details (numbered, with timing), expected outcomes, quantification and statistical analysis,
limitations, troubleshooting, resource availability, acknowledgments, author contributions,
declarations, figure legends with figures, references.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import crm_methods as M  # noqa: E402
import crm_numbers  # noqa: E402
import star_text as S  # noqa: E402
from build_crm import fill_all  # noqa: E402
from docbuild import Manuscript, word_count  # noqa: E402

CRM = Path(__file__).resolve().parents[1]
FIG = CRM / "star/figures"
SPECIAL = ("CRITICAL:", "Note:", "Pause point:", "Optional:")


def steps(m, items, counter):
    for st in items:
        tag = next((t for t in SPECIAL if st.startswith(t)), None)
        if tag:
            m.para(f"**{tag}**{st[len(tag):]}", indent=8, size=10.5)
        else:
            counter[0] += 1
            m.para(f"{counter[0]}.\t{st}", indent=0)


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    n = crm_numbers.compute()
    for k in dir(S):
        if k.isupper():
            setattr(S, k, fill_all(getattr(S, k), n))
    assert word_count([("p", S.SUMMARY)]) <= 80, f"summary {word_count([('p', S.SUMMARY)])} words"
    for h in S.HIGHLIGHTS:
        assert len(h.replace("_{H}", "H")) <= 85, f"highlight too long: {h}"
    m = Manuscript()
    m.para(S.TITLE, size=16, bold=True)
    m.para(S.AUTHORS)
    for a in S.AFFIL:
        m.para(a, size=10, space_after=0)
    m.heading("Summary", 1)
    m.para(S.SUMMARY)
    m.para(S.SUMMARY_TAIL)
    m.heading("Highlights", 1)
    m.bullets(S.HIGHLIGHTS)
    m.heading("Graphical abstract", 1)
    m.para(S.GA_NOTE)
    # before you begin (its numbered steps are counted separately, as in the template)
    m.heading("Before you begin", 1)
    pre = [0]
    for kind, payload in S.BEFORE:
        if kind == "h":
            m.heading(payload, 2)
        elif kind == "p":
            m.para(payload)
        else:
            steps(m, payload, pre)
    m.heading("Key resources table", 1)
    # the other callers benchmarked in the article are not used in the protocol
    krt = [r for r in M.KRT if not r[0].startswith(("ClairS", "DeepVariant"))]
    m.table(M.KRT_HEADER, fill_all(krt, n), widths=[75, 50, 45], size=8)
    m.heading("Materials and equipment", 1)
    for title, header, rows in S.MATERIALS_TABLES:
        m.table(header, rows, title=f"**{title}**", size=9)
        m.para("")
    m.para(S.MATERIALS_NOTE)
    m.heading("Step-by-step method details", 1)
    n = [0]
    for title, timing, items in S.STEPS:
        m.heading(title, 2)
        m.para(f"**{timing}**", size=10.5)
        steps(m, items, n)
    m.heading("Expected outcomes", 1)
    for p in S.EXPECTED:
        m.para(p)
    m.heading("Quantification and statistical analysis", 1)
    for p in S.QUANT:
        m.para(p)
    m.heading("Limitations", 1)
    m.para(S.LIMITATIONS)
    m.heading("Troubleshooting", 1)
    for prob, issue, sol in S.TROUBLE:
        m.heading(prob, 2)
        m.para(issue)
        m.heading("Potential solution", 3)
        m.para(sol)
    m.heading("Resource availability", 1)
    for h, t in S.RESOURCE:
        m.heading(h, 2)
        m.para(t)
    m.heading("Acknowledgments", 1)
    m.para(S.ACK)
    m.heading("Author contributions", 1)
    m.para(S.CONTRIB)
    m.heading("Declaration of interests", 1)
    m.para(S.INTERESTS)
    m.heading("Declaration of generative AI and AI-assisted technologies in the manuscript preparation process", 1)
    m.para(S.AI_DECL)
    m.page_break()
    m.heading("Figure titles and legends", 1)
    files = ["STAR_Figure1_preview.png", "STAR_Figure2.png", "STAR_Figure3.png", "STAR_Figure4.png"]
    for (title, legend), f in zip(S.LEGENDS, files):
        m.figure(FIG / f, width_mm=165)
        m.para(f"**{title}**")
        m.para(legend, size=10)
    m.references()
    path = out / "SHMtrace_STARProtocols.docx"
    m.save(path)
    print(path, f"| steps {n[0]} | summary {word_count([('p', S.SUMMARY)])} words | references {len(m.order)}")


if __name__ == "__main__":
    main(sys.argv[1])
