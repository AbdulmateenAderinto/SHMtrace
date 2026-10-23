"""V2's separate submission file: highlights, eTOC blurb and Motivation, as Cell Press asks.

Starts from the canonical build (docs/calibration/crm/build/Highlights_eTOC_Motivation.docx, 21 Oct), whose numbers
build_crm.py filled from the analysis, and changes only what V2 changed: with Philip Barbulescu added the paper has
three authors, so the eTOC names "Aderinto et al." (Cell Press convention). The journal limits build_crm.py checks are
checked again: four highlights of at most 85 characters, an eTOC blurb of at most 50 words.

    ../../../../../.venv/bin/python highlights_v2.py
"""
import sys
from pathlib import Path

import docx

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CANON = BASE.parents[3] / "docs/calibration/crm/build/Highlights_eTOC_Motivation.docx"
sys.path.insert(0, str(BASE.parent / "codex_version/sources"))
from docbuild import Manuscript  # noqa: E402

sections, current = {}, None
for p in docx.Document(str(CANON)).paragraphs:
    t = p.text.strip()
    if t in ("Highlights", "eTOC blurb", "Motivation"):
        current = sections.setdefault(t, [])
    elif t and current is not None:
        current.append(t)

highlights, (etoc,), (motivation,) = sections["Highlights"], sections["eTOC blurb"], sections["Motivation"]
assert etoc.startswith("Aderinto and Martin present"), etoc[:40]
etoc = etoc.replace("Aderinto and Martin present", "Aderinto et al. present", 1)

# "1–4 calls" says nothing without the background it is out of (author, 23 Oct). The highlight now gives the
# proportion the Summary reports: the drop in AID-null calls against the uncorrected rule on the same libraries,
# per assay. Filled from crm_numbers.compute(), the numbers the Summary's {ctl_red_lo}–{ctl_red_hi}% come from
# (needs SHMtrace, so run with .venv).
REPO = BASE.parents[3]
for p in ("docs/calibration/crm/manuscript", "docs/calibration/crm", "src"):
    sys.path.insert(0, str(REPO / p))
import crm_numbers  # noqa: E402

N = crm_numbers.compute()
assert highlights[1].startswith("AID-null controls retain"), highlights[1]
highlights[1] = f"AID-null control calls fall {N['ctl_red_lo']}–{N['ctl_red_hi']}% vs an uncorrected rule; true mutations are kept"
# the graphical abstract types the same numbers into its builder; they must match the analysis too
ga = (BASE / "scicanvas/crm_ga.js").read_text()
assert f"{N['ctl_red_lo']}–{N['ctl_red_hi']}% fewer control calls" in ga, "graphical abstract is stale"
assert f"r = {N['r_spec']}" in ga, "graphical abstract is stale"

# highlights 3 and 4 (author, 23 Oct): SHMtrace is SHM-specific as well as error-aware, and its readouts transfer
# across platforms. PacBio is left out of 4 because its dataset is class-switch junctions, not AID signatures, and
# the Sanger route reports background rather than suppressing it, hence "low background"
assert highlights[2].startswith("Nanopore gives spectra similar to Sanger"), highlights[2]
assert highlights[3].startswith("Re-running three published Sanger studies"), highlights[3]
highlights[2] = "Calls feed SHM readouts: AID hotspot targeting, spectra, clones and mutation trees"
highlights[3] = "AID signatures are recovered from Sanger, nanopore and Illumina with low background"

# Motivation (author, 23 Oct): the gap is not only cloning but the shape of the software built since. The SHM
# toolkits follow B-cell receptor repertoire sequencing and start from assembled repertoires, so a defined reporter
# locus read in bulk has had no modern route; nanopore reads plus control-calibrated calling give it one. Written
# from what is verifiable: Sible et al. 2021 is the current J_H4 protocol, and Immcantation, IgBLAST, SHazaM and
# MiXCR are repertoire tools (the Discussion and Figure S13 make the same point). It says "has followed", not
# "is only for", since MiXCR documents a nanopore preset; nothing is claimed about any tool's availability.
assert motivation.startswith("At genomic reporter loci"), motivation[:40]
assert motivation.endswith("discrepancies explicit."), motivation[-40:]
motivation = (
    "At genomic reporter loci such as the J_{H}4 intron, somatic hypermutation is still measured mostly by cloning "
    "PCR products and sequencing colonies one at a time, which is slow and samples few molecules. The software "
    "built for somatic hypermutation since has followed B-cell receptor repertoire sequencing, starting from "
    "assembled repertoires and summarising mutations already identified; a defined reporter locus, read in bulk on "
    "whatever platform a study uses and judged against that study's own background, has had no equivalent modern "
    "route. Long nanopore reads could supply one, replacing cloning while preserving linkage between mutations "
    "across the amplicon, but their per-read error rate is close to the mutation rate of interest, and "
    "general-purpose variant callers judge each base by its quality score rather than by the errors a given library "
    "actually makes. We therefore built SHMtrace to calibrate each call against matched AID-null libraries "
    "sequenced in the same run, to treat animals rather than reads as replicates, and to make the whole workflow, "
    "from reads or chromatograms to group statistics, reproducible.")

assert len(highlights) == 4 and all(len(h) <= 85 for h in highlights), [len(h) for h in highlights]
assert len(etoc.split()) <= 50, len(etoc.split())

h = Manuscript(line_numbers=False)
h.para("SHMtrace V2 – for Alberto review: highlights, eTOC blurb and Motivation (a separate file at submission)",
       size=10)
h.heading("Highlights", 1)
h.bullets(highlights)
h.heading("eTOC blurb", 1)
h.para(etoc)
h.heading("Motivation", 1)
h.para(motivation)
out = BASE / "SHMtrace_V2_Highlights_eTOC_Motivation.docx"
h.save(out)
print(out.name, "| highlights (chars):", [len(x) for x in highlights], "| eTOC words:", len(etoc.split()),
      "| Motivation words:", len(motivation.split()))
