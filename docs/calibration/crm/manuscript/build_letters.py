"""Cover letters for the Cell Reports Methods article and the STAR Protocols protocol (drafts).

    python docs/calibration/crm/manuscript/build_letters.py OUTDIR
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import crm_numbers  # noqa: E402
import crm_text as T  # noqa: E402
from build_crm import fill  # noqa: E402
import star_text as S  # noqa: E402
from docbuild import Manuscript  # noqa: E402

CRM_LETTER = [
    "Dear Editor,",
    f"We are pleased to submit our manuscript, \"{T.TITLE},\" for consideration as an Article in *Cell Reports "
    f"Methods*.",
    "Somatic hypermutation (SHM) is still measured mostly by cloning PCR products and sequencing colonies one at a "
    "time, which is slow and samples few molecules. Long nanopore reads could replace cloning, but a single read "
    "makes errors at rates close to the mutation rates being measured, and general-purpose variant callers judge "
    "each base by its quality score rather than by the errors that a given library actually makes. In this "
    "manuscript, we describe SHMtrace, which tests each candidate mutation against matched AID-deficient libraries "
    "sequenced in the same run, compares groups with the animal as the unit, and runs the complete analysis from "
    "delivered reads or Sanger chromatograms, from the command line or a desktop application.",
    "The manuscript provides the three elements that the journal asks of a method: the method itself, a "
    "comparison with related methods, and validation. We compared SHMtrace with LoFreq and with the uncorrected "
    "test used in commercial software on the same reads, using planted mutations, simulated B cell lineages, and "
    "nanopore libraries of four loci from two laboratories. SHMtrace kept the AID-deficient background at "
    "{ctl_rec} calls per assay while retaining the wild-type signal, and it separated AID-active from AID-deficient "
    "libraries {sep_lo}- to {sep_hi}-fold in four of five assays, compared with {lf_lo}- to {gen_hi}-fold for the "
    "alternatives. Nanopore and Sanger sequencing of the same mice gave concordant mutation spectra, and "
    "full-length reads additionally resolved whole-amplicon mutation maps, B cell clones, and mutation trees. The "
    "same analysis reproduced published findings in three public Illumina datasets, including a negative result "
    "and an indel phenotype.",
    "We believe that this work will interest immunologists who measure SHM and researchers who call rare variants "
    "against biological controls. A companion protocol is being submitted to *STAR Protocols*. All data, code, and "
    "software will be publicly available at publication. This manuscript has not been published or submitted "
    "elsewhere, all authors have approved it, and we declare no competing interests.",
    "Sincerely,",
    "Alberto Martin, on behalf of the authors",
]
STAR_LETTER = [
    "Dear Editor,",
    f"We are pleased to submit \"{S.TITLE}\" for consideration as a Protocol in *STAR Protocols*.",
    "This protocol takes a laboratory from Peyer's patch germinal-center B cells to calibrated measurements of "
    "somatic hypermutation without cloning. It covers cell sorting, J_{H}4 amplification, full-length nanopore "
    "sequencing through a service provider, and analysis with SHMtrace, open-source software with a desktop "
    "application that tests mutations against matched AID-deficient controls and compares groups animal by "
    "animal. The protocol includes a key resources table, timing for every step, expected outcomes from real "
    "data, and troubleshooting based on problems that we and others have encountered, including co-amplified "
    "paralogs, dominant clones, error clusters, and a contaminated control.",
    "The protocol accompanies our Cell Reports Methods submission, which describes and validates SHMtrace. It has "
    "not been published or submitted elsewhere, and we declare no competing interests.",
    "Sincerely,",
    "Abdulmateen Aderinto and Alberto Martin",
]


def main(outdir):
    out = Path(outdir)
    n = crm_numbers.compute()
    for name, paras in (("CoverLetter_CellReportsMethods.docx", CRM_LETTER),
                        ("CoverLetter_STARProtocols.docx", STAR_LETTER)):
        m = Manuscript(line_numbers=False)
        for p in paras:
            m.para(fill(p, n))
        m.save(out / name)
        print(out / name)


if __name__ == "__main__":
    main(sys.argv[1])
