"""Cover letters for the Cell Reports Methods article and the STAR Protocols protocol (drafts).

    python docs/calibration/crm/manuscript/build_letters.py OUTDIR
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import crm_text as T  # noqa: E402
import star_text as S  # noqa: E402
from docbuild import Manuscript  # noqa: E402

CRM_LETTER = [
    "Dear Editor,",
    f"We submit our manuscript, \"{T.TITLE},\" for consideration as an Article in *Cell Reports Methods*.",
    "Somatic hypermutation (SHM) is still measured mostly by cloning PCR products and Sanger sequencing colonies "
    "one at a time. Long reads could replace cloning, but a single read errs at rates close to the mutation rates "
    "of interest, and general-purpose variant callers judge each base by its own quality score rather than by the "
    "errors a given library actually makes. We describe SHMtrace, which calibrates every candidate mutation "
    "against matched AID-null libraries sequenced in the same run, tests groups with the animal as the unit, and "
    "runs the whole analysis from delivered reads or Sanger chromatograms, from the command line or a desktop "
    "application.",
    "The paper provides the three things the journal asks of a method: the method, a comparison with related "
    "methods, and validation. We compared SHMtrace with LoFreq and with the uncorrected test used in commercial "
    "software on the same reads, using planted mutations, simulated B cell lineages and nanopore libraries of "
    "four loci from two laboratories. SHMtrace kept the AID-null background at 0–3 calls per assay while keeping the "
    "wild-type signal, and separated AID-active from AID-null libraries 92- to 137-fold in four of five assays, where the "
    "alternatives reached 2- to 62-fold. Nanopore and Sanger sequencing of the same mice gave concordant spectra, "
    "and full-length reads added whole-amplicon maps, clones and mutation trees. The same engine reproduced "
    "published findings on three public Illumina datasets, including a negative result and an indel phenotype.",
    "We believe the work will interest immunologists who measure SHM and anyone calling rare variants against "
    "biological controls. A companion protocol is planned for *STAR Protocols*. All data, code and the software "
    "will be public at publication. The manuscript has not been published or submitted elsewhere, and all authors "
    "approved it. We declare no competing interests.",
    "Sincerely,",
    "Alberto Martin, on behalf of the authors",
]
STAR_LETTER = [
    "Dear Editor,",
    f"We submit \"{S.TITLE}\" as a Protocol for *STAR Protocols*.",
    "The protocol takes a laboratory from Peyer's patch germinal-center B cells to calibrated SHM measurements "
    "without cloning: sorting, J_{H}4 amplification, full-length nanopore sequencing through a service provider, "
    "and analysis with SHMtrace, open-source software with a desktop application that calls mutations against "
    "matched AID-null controls and compares groups animal by animal. It includes a key resources table, timing for "
    "every step, expected outcomes from real data, and troubleshooting drawn from problems we and others met, "
    "including co-amplified paralogs, dominant clones, error clusters and a contaminated control.",
    "The protocol accompanies our Cell Reports Methods submission describing and validating SHMtrace. It has not "
    "been published or submitted elsewhere, and we declare no competing interests.",
    "Sincerely,",
    "Abdulmateen Aderinto and Alberto Martin",
]


def main(outdir):
    out = Path(outdir)
    for name, paras in (("CoverLetter_CellReportsMethods.docx", CRM_LETTER),
                        ("CoverLetter_STARProtocols.docx", STAR_LETTER)):
        m = Manuscript(line_numbers=False)
        for p in paras:
            m.para(p)
        m.save(out / name)
        print(out / name)


if __name__ == "__main__":
    main(sys.argv[1])
