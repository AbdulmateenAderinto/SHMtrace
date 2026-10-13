"""Ver18: every ONT result from SHMtrace (JH4 and N114), from Ver17.

    python revise_manuscript_v18.py VER17.docx VER18.docx      (needs python-docx)

Ver17 (revise_manuscript_sanger.py) re-analysed the Sanger clones and N114. Ver18 replaces
the remaining Geneious ONT results with SHMtrace's, all on the deposit's realised 500-read
subsets (the inputs Methods 2.7 names):
  JH4  docs/calibration/jh4_shmtrace.py (529-bp window) and --full-length (positions 90-1,230)
  N114 docs/calibration/n114_shmtrace.py (now the sub500 files, as for JH4)
  depth series docs/calibration/depth_shmtrace.py; Sanger-ONT agreement sanger_vs_nanopore.py
The calling Methods become one SHMtrace section for both assays. Figures 1-4 and
Supplementary S1-S7 and S9 are replaced. Every replaced passage is highlighted yellow;
gene names, J_H4 and 10^-5 in it are re-formatted (italic, subscript, superscript).
SOURCE is never modified, and the script stops if any expected text is not found.
"""
import copy
import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_COLOR_INDEX
from docx.oxml.ns import qn
from docx.text.run import Run

from revise_manuscript_sanger import set_cell, swap_image

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
SV = HERE.parent / "figures/sanger_vs_nanopore/recommended_manuscript"

# (paragraph starts with, old text, new text); old "" with new = whole paragraph is REWRITE
EDITS = [
    # ---- Abstract
    ("Somatic hypermutation (SHM) diversifies",
     "(Geneious for JH4; for N114, a quality-calibrated test against the empty-vector controls)",
     "(SHMtrace, a quality-calibrated test of each candidate against the assay's negative controls at a "
     "false-discovery rate of 0.05)"),
    ("Somatic hypermutation (SHM) diversifies", "(Pearson r = 0.95)", "(Pearson r = 0.91)"),
    ("Somatic hypermutation (SHM) diversifies",
     "ONT detected at least a 91-fold group difference, with no calls in any empty-vector control,",
     "ONT detected at least an 89-fold group difference, with one call across the three empty-vector controls,"),
    ("Somatic hypermutation (SHM) diversifies",
     "the single realised subsampling series recovered 95–97% of the maximum number of filtered positions "
     "observed at 1,000 reads; AID-deficient controls retained a post-filter apparent-variant background of "
     "~1.2 × 10−5 per analysed base.",
     "the single realised subsampling series recovered a mean of 91% of the filtered positions observed at the "
     "deepest analysed depth; AID-deficient controls retained a post-filter apparent-variant background of "
     "~0.5 × 10−5 per analysed base."),
    ("Somatic hypermutation (SHM) diversifies",
     "lacks unique molecular identifiers, does not correct candidate-level p-values for multiplicity, and did not",
     "lacks unique molecular identifiers and did not"),
    ("Abbreviations:", "CI, confidence interval.", "CI, confidence interval; FDR, false-discovery rate."),
    # ---- Methods 2.8: one SHMtrace calling section
    ("Quality-control-passing reads were aligned independently", "",
     "Quality-control-passing reads were analysed for each sample with SHMtrace v0.2, written for this study "
     "(Data and code availability). Reads were aligned with minimap2 (map-ont) to the 1,230-nt JH4 analysis "
     "reference or the 655-nt N114 analysis reference, which were distinct from the shorter mapping references "
     "used during preprocessing. Every single-nucleotide difference from the reference within the analysed "
     "window was a candidate position–substitution; a single supporting read was sufficient for nomination, "
     "corresponding to a nominal ~0.2% allele fraction at the 500-read working depth. No minimum "
     "allele-frequency criterion was imposed at nomination."),
    ("The Geneious approximate variant p-value was used", "",
     "Each candidate was tested against the sequencing background measured in the assay's negative-control "
     "libraries (Aicda−/− for JH4; empty vector for N114). The reads supporting it at or above each base "
     "quality were compared with the number expected from the controls' error rate for that substitution, "
     "strand and quality, blended with the position's own control rate. Because quality scores overstate "
     "certainty, the control rates were multiplied by a factor chosen by leave-one-control-out: the smallest "
     "at which no more than 1% of tested positions in a held-out control reached p < 0.01. This factor was 5 "
     "for both assays; at face value, 56 of 727 such JH4 tests and 72 of 580 N114 tests reached p < 0.01. "
     "Benjamini–Hochberg correction was applied across all positions and alternate alleles at a "
     "false-discovery rate of 0.05, and each control library was itself called against a background fitted "
     "from the other two, never from itself. A passing call is a filtered apparent variant, not an "
     "individually confirmed mutation; the controls provide empirical assay-level estimates of residual "
     "apparent calls after the complete workflow. The same reads were originally analysed in Geneious Prime "
     "2026.1.2 (Kearse et al., 2012) with an unadjusted approximate per-candidate p ≤ 0.05; those calls are "
     "retained in the data deposit and are compared in Section 3.1."),
    ("Germline and strain-background positions were handled",
     "Competitive allele selection and the recurrent-contaminant screen reduce,",
     "SHMtrace additionally removed reads carrying the recurrent-contaminant allele at ≥ 50% of its diagnostic "
     "sites. Competitive allele selection and the recurrent-contaminant screens reduce,"),
    ("For N114, the 651-nt reference used for preprocessing",
     "the 655-nt parental sequence used for Geneious variant analysis",
     "the 655-nt parental sequence used for variant analysis"),
    ("For N114, the 651-nt reference used for preprocessing",
     "refer to the 655-nt Geneious reference.", "refer to the 655-nt analysis reference."),
    ("For N114, the 651-nt reference used for preprocessing",
     " Removing high-frequency positions affected at most two calls per pMX-Aicda subclone and none of the "
     "empty-vector subclones.", ""),
    ("Sanger clones were re-analysed from their original", "with shmtrace v0.2", "with SHMtrace v0.2"),
    ("For the separate full-length ONT analysis",
     "per-position read depth was computed from the alignments.",
     "per-position read depth was computed from the alignments, and calls were made in a separate SHMtrace "
     "run over the full interval with otherwise identical settings (leave-one-control-out factor 10)."),
    ("Depth-standardised unique-variant densities were log10-transformed",
     "After re-analysis every sample contained at least one call, so no left-censoring was required.",
     "A library with no calls (one Aicda−/− ONT library and two empty-vector subclones) was placed at its "
     "one-event density, one call divided by its reads and analysed bases; the corresponding fold differences "
     "are therefore lower bounds. Every Sanger sample contained at least one call."),
    # ---- Results 3.1
    ("In these samples, Wt mice had markedly higher",
     "p = 0.0037 for ONT (~47-fold)",
     "p = 0.0018 for ONT (≥ 108-fold, the one Aicda−/− library without calls placed at its one-event density)"),
    ("In these samples, Wt mice had markedly higher", "~23 versus ~1.2 × 10−5", "~23 versus ~0.55 × 10−5"),
    ("After all filters, Aicda−/− libraries retained", "",
     "After all filters, the three Aicda−/− libraries retained 0, 3 and 1 calls, a geometric-mean apparent "
     "background of ~0.55 × 10−5 distinct calls per analysed base (pooled, 4 calls in 787,576 analysed "
     "read-bases, about one per 200,000). Expressed as Q = −10 × log10(rate), this is a Q53-equivalent "
     "assay-level indicator, not raw ONT accuracy or an individual-call false-positive probability. Residual "
     "calls may include technical error, PCR error, low-frequency haplotypes, contamination, or genuine "
     "AID-independent variation. In the three Wt libraries, 422, 578 and 323 candidate position–substitutions "
     "were tested and 155, 347 and 69 retained; in the Aicda−/− libraries, 243, 240 and 244 were tested "
     "(Supplementary Table S1). The original Geneious analysis of the same reads, with an unadjusted "
     "per-candidate p ≤ 0.05, retained a similar number of Wt calls (161, 342 and 68) but more Aicda−/− calls "
     "(2, 6 and 3; ~1.2 × 10−5) and a smaller group difference (~47-fold; p = 0.0037). Calling against the "
     "controls with false-discovery control therefore mainly reduced the background."),
    # ---- Results 3.2
    ("Because ONT reads span the complete amplicon",
     "(36 versus 1.7 × 10−5; ~21-fold)", "(32 versus 0.35 × 10−5; ~92-fold)"),
    ("Because ONT reads span the complete amplicon",
     "density declined roughly ten-fold toward the 3′ end",
     "density declined more than ten-fold toward the 3′ end"),
    ("Because ONT reads span the complete amplicon",
     "(Pearson r = 0.989; all differences ≤1.5 percentage points)",
     "(Pearson r = 0.974; all differences ≤ 2.2 percentage points)"),
    ("Because ONT reads span the complete amplicon", "(739 versus 571;", "(680 versus 571;"),
    # ---- Results 3.3
    ("The matched ONT and Sanger JH4 datasets",
     "Pearson r = 0.95 and slope = 1.05", "Pearson r = 0.91 and slope = 1.02"),
    ("The matched ONT and Sanger JH4 datasets", "(p = 0.78)", "(p = 0.50)"),
    ("The matched ONT and Sanger JH4 datasets", "(smallest adjusted p = 0.12)", "(smallest adjusted p = 0.46)"),
    ("The matched ONT and Sanger JH4 datasets", "(25-bp bins; Pearson r = 0.87)", "(25-bp bins; Pearson r = 0.88)"),
    ("The matched ONT and Sanger JH4 datasets",
     "Sanger 58%, ONT 55%; A:T calls at WA/TW: 70% and 69%",
     "Sanger 58%, ONT 56%; A:T calls at WA/TW: 70% on both"),
    # ---- Results 3.4
    ("We next tested the workflow in the fixed N114",
     "of the 655-nt Geneious reference", "of the 655-nt analysis reference"),
    ("At matched 500-read depth, pMX-Aicda subclones",
     "density of 30 × 10−5, whereas no empty-vector control had any call (Figure 4A). With each control placed "
     "at its one-event density (0.33 × 10−5), the group difference was at least 91-fold (exploratory Welch "
     "p = 0.0045).",
     "density of 29 × 10−5, whereas the three empty-vector controls had one call between them (0, 0 and 1; "
     "Figure 4A). With the two controls without calls placed at their one-event density (0.33 × 10−5), the "
     "group difference was at least 89-fold (exploratory Welch p = 0.0034)."),
    ("At matched 500-read depth, pMX-Aicda subclones", "whereas controls had none", "whereas the controls had one call in total"),
    # ---- Results 3.5
    ("We evaluated one independently drawn, non-nested reservoir", "",
     "We evaluated one independently drawn, non-nested reservoir subsample per library at 30, 60, 125, 250, "
     "500, and 1,000 reads (nominal seed 42) and re-called variants at each depth, with the controls "
     "calibrated at the same depth (Supplementary Figures S6 and S7). In AID-active samples, the number of "
     "filtered position–substitution calls rose with depth; at 500 reads, mean counts were 91% of each "
     "sample's maximum for both JH4 and N114. Controls retained 0–3 calls from 60 reads upward. At 30 reads, "
     "however, the JH4 Aicda−/− libraries themselves received 12–29 calls each: with 30 reads per control the "
     "background model cannot be calibrated, and calls at that depth are not interpretable. These percentages "
     "describe this finite dataset and single subsampling realisation; they do not estimate recovery of the "
     "unknown true mutation repertoire or establish a universal saturation depth."),
    ("Depth-standardised unique-variant density declines", "",
     "Below ~250 reads, a mutation seen in a single read rarely passes the false-discovery threshold against "
     "the controls, so libraries whose mutations are mostly confined to single reads are under-called at low "
     "depth. Wt2, whose calls are mostly single-read, had 3 and 2 calls at 60 and 125 reads against 110 at "
     "250, whereas the clonally shared mutations of Wt3 were called at every depth. Depth-standardised density "
     "therefore did not show the simple decline expected of a saturating numerator: it rose with depth in N114 "
     "(from ~9 to 29 × 10−5 between 30 and 500 reads) and was irregular in JH4 (Supplementary Figure S6). "
     "Greater depth does not improve an individual read's accuracy; it samples more amplified products and "
     "gives each candidate more opportunity to reach significance. Cross-platform support therefore comes from "
     "the shared group-level AID dependence and descriptively similar aggregate spectrum, not equality of "
     "absolute density or confirmation of each position."),
    # ---- Discussion
    ("The present data do not establish the accuracy of every",
     "The unadjusted Geneious per-candidate threshold used for JH4 also does not control false discovery across "
     "the many positions and alternate alleles examined; the N114 calls were made at a false-discovery rate "
     "of 0.05.",
     "Calls were made at a false-discovery rate of 0.05 against a background fitted from the controls, but "
     "that background is itself estimated from three libraries per assay, and it cannot be calibrated at very "
     "low depth (Supplementary Figure S6)."),
    ("The present implementation is intentionally substitution-focused",
     "multiplicity-controlled and externally calibrated calling", "externally calibrated calling"),
    # ---- Conclusions
    ("This proof-of-concept study shows that", "(Pearson r = 0.95)", "(Pearson r = 0.91)"),
    ("This proof-of-concept study shows that",
     "at least a 91-fold difference between AID-expressing and empty-vector N114 subclones, with no calls in "
     "any empty-vector control.",
     "at least an 89-fold difference between AID-expressing and empty-vector N114 subclones, with one call "
     "across the three empty-vector controls."),
    ("This proof-of-concept study shows that", "~1.2 × 10−5 apparent variants", "~0.5 × 10−5 apparent variants"),
    # ---- Data and code
    ("Data and code availability:", "reference sequences, Geneious variant-call exports,",
     "reference sequences, SHMtrace variant calls and project files, the original Geneious variant-call exports,"),
    ("Data and code availability:", "All analysis and figure-generation code will be made publicly available "
     "upon publication.", "SHMtrace (v0.2) and all analysis and figure-generation code will be made publicly "
     "available upon publication."),
    # ---- Supplementary Table legends
    ("Supplementary Table S1. Per-sample sequencing", "",
     "Supplementary Table S1. Per-sample sequencing and candidate-variant summary. Reads delivered by the "
     "sequencing provider and reads satisfying assay-specific target selection, amplicon-length, and "
     "expected-error criteria are shown for ONT libraries; Sanger entries report independent clones. Reads or "
     "clones analysed per sample are followed by the distinct position–substitution candidates tested, the "
     "number retained by SHMtrace's quality-calibrated test against the assay's negative controls at a "
     "false-discovery rate of 0.05 (Methods), and the resulting depth-standardised unique-variant density "
     "(× 10−5). For Sanger, every substitution with base quality ≥ 20 is a call, so candidates and calls are "
     "the same. JH4 densities use 529 analysed bp (for Sanger, the bases each clone read within them); Sanger "
     "clones were re-analysed from the original chromatograms and kept with the mouse they were originally "
     "filed under. N114 densities use positions 21–635 of the 655-nt analysis reference (615 bp). Passing "
     "calls are filtered apparent variants, not individually confirmed mutations."),
    ("Group geometric means with descriptive 95%", "Wt ONT 59 [26–129]", "Wt ONT 59 [26–134]"),
    ("Group geometric means with descriptive 95%", "Aicda−/− ONT 1.2 [0.8–2.3]",
     "Aicda−/− ONT 0.55 [0.38–1.14] (one library with no calls, at its one-event density)"),
    ("Group geometric means with descriptive 95%",
     "empty vector: no calls (below 0.33, the one-event density); pMX-Aicda 30 [19–53]",
     "empty vector 0.33 (two subclones with no calls, at their one-event density, and one with one call); "
     "pMX-Aicda 29 [22–50]"),
    # ---- Figure legends
    ("Figure 1. Experimental design",
     "(5) Geneious candidate substitutions with approximate p ≤ 0.05 and frequency < 50% were retained",
     "(5) Candidate substitutions passing SHMtrace's quality-calibrated test against the negative controls "
     "(false-discovery rate 0.05) and with frequency < 50% were retained"),
    ("Figure 2. Depth-standardised", "after the stated fixed-position and Geneious p-value filters",
     "after the stated fixed-position filter (ONT calls: SHMtrace, false-discovery rate 0.05 against the "
     "Aicda−/− controls)"),
    ("Figure 2. Depth-standardised", "every mouse had at least one call.",
     "every Sanger mouse had at least one call. One Aicda−/− ONT library had none and is drawn (open circle) "
     "and tested at its one-event density (0.38 × 10−5), so the ONT fold difference is a lower bound."),
    ("Figure 2. Depth-standardised", "~47-fold (ONT; p = 0.0037)", "≥ 108-fold (ONT; p = 0.0018)"),
    ("Figure 2. Depth-standardised", "Wt ONT 59 [26–129], Aicda−/− ONT 1.2 [0.8–2.3]",
     "Wt ONT 59 [26–134], Aicda−/− ONT 0.55 [0.38–1.14]"),
    ("Figure 2. Depth-standardised", "shading, overlap with the shared 529-bp window.",
     "shading, the shared 529-bp window (positions 101–629). Panel C is a separate SHMtrace run over "
     "positions 90–1,230 with otherwise identical settings."),
    ("Figure 3. ONT yields",
     "For ONT, positions with frequency ≥ 50% relative to the C57BL/6 reference and Geneious approximate "
     "p > 0.05 were removed.",
     "ONT calls were made with SHMtrace at a false-discovery rate of 0.05 against the Aicda−/− controls, after "
     "removing positions with frequency ≥ 50% relative to the C57BL/6 reference."),
    ("Figure 3. ONT yields", "interaction was p = 0.78", "interaction was p = 0.50"),
    ("Figure 3. ONT yields", "(Pearson r = 0.95; slope = 1.05;", "(Pearson r = 0.91; slope = 1.02;"),
    ("Figure 3. ONT yields", "adjusted paired p = 0.43 and 1.00", "adjusted paired p = 0.10 and 0.97"),
    ("Figure 4. Filtered apparent variants", "of the 655-nt Geneious reference", "of the 655-nt analysis reference"),
    ("Figure 4. Filtered apparent variants",
     "No empty-vector subclone had any call; each is drawn (open circle) and tested at its one-event density "
     "(0.33 × 10−5), so the group difference is a lower bound: the exploratory Welch test on log10-transformed "
     "densities gave p = 0.0045 and at least a 91-fold difference (30 versus < 0.33 × 10−5).",
     "Two empty-vector subclones had no calls and the third had one; those without calls are drawn (open "
     "circles) and tested at their one-event density (0.33 × 10−5), so the group difference is a lower bound: "
     "the exploratory Welch test on log10-transformed densities gave p = 0.0034 and at least an 89-fold "
     "difference (29 versus 0.33 × 10−5)."),
    ("Figure 4. Filtered apparent variants", "(n = 298) and empty vector (n = 0)",
     "(n = 287) and empty vector (n = 1)"),
    ("Supplementary Figure S1. Pooled substitution",
     "after the fixed-position and Geneious p-value filters",
     "after the fixed-position filter and SHMtrace calling (false-discovery rate 0.05; the full-length "
     "interval is a separate run with otherwise identical settings)"),
    ("Supplementary Figure S1. Pooled substitution", "(1,141 bp; n = 739)", "(1,141 bp; n = 680)"),
    ("Supplementary Figure S1. Pooled substitution",
     "(Pearson r = 0.989; all differences ≤1.5 percentage points)",
     "(Pearson r = 0.974; all differences ≤ 2.2 percentage points)"),
    ("Supplementary Figure S3. Individual JH4", "(n = 161, 342, and 68 unique variants",
     "(n = 155, 347, and 69 unique variants"),
    ("Supplementary Figure S3. Individual JH4", "(n = 2, 6, and 3)", "(n = 0, 3, and 1)"),
    ("Supplementary Figure S3. Individual JH4",
     "Variants were restricted to single-nucleotide substitutions with variant frequency < 50% and Geneious "
     "approximate variant p ≤ 0.05.",
     "Calls are single-nucleotide substitutions retained by SHMtrace's quality-calibrated test against the "
     "Aicda−/− controls at a false-discovery rate of 0.05, after removing positions at which ≥ 50% of reads "
     "differ from the reference (Methods)."),
    ("Supplementary Figure S4. Individual N114", "(n = 79, 58, and 161 unique variants",
     "(n = 69, 68, and 150 unique variants"),
    ("Supplementary Figure S4. Individual N114", "(n = 0, 0, and 0)", "(n = 0, 0, and 1)"),
    ("Supplementary Figure S5. Recurrence-weighted", "(n = 1,441 read-weighted observations); (B) three "
     "empty-vector controls (n = 0)", "(n = 1,455 read-weighted observations); (B) three empty-vector "
     "controls (n = 4)"),
    ("Supplementary Figure S6. Sequencing depth", "", None),          # rewritten below
    ("Supplementary Figure S7. Filtered apparent-variant",
     "At 500 reads, mean counts were 97% for JH4 and 95% for N114 relative to 1,000 reads.",
     "At 500 reads, mean counts were 91% of each sample's maximum for both JH4 and N114. The shaded 30-read "
     "JH4 point is not interpretable: the background model could not be calibrated at that depth "
     "(Supplementary Figure S6)."),
    ("Supplementary Figure S7. Filtered apparent-variant", "estimate 95–97% sensitivity", "estimate 91% sensitivity"),
    ("Supplementary Figure S9. Sanger and ONT agree",
     "ONT calls use the manuscript's Geneious rule (reproduced in shmtrace);",
     "ONT calls are the SHMtrace calls of Figure 3;"),
    ("Supplementary Figure S9. Sanger and ONT agree", "(Pearson r = 0.87)", "(Pearson r = 0.88)"),
    ("Supplementary Figure S2. Individual JH4", "(shmtrace v0.2, Sanger mode)", "(SHMtrace v0.2, Sanger mode)"),
]
S6_LEGEND_BODY = (
    "Sequencing depth versus filtered position–substitution discovery and depth-standardised unique-variant "
    "density. JH4 and N114 libraries were sampled once at 30, 60, 125, 250, 500, and 1,000 reads with "
    "independent, non-nested reservoir samples and nominal seed 42; variants were re-called at each depth with "
    "SHMtrace, the controls calibrated at the same depth. (A, C) Distinct calls versus depth; open circles, "
    "geometric mean; whiskers, geometric SD (n = 3 per group). AID-active samples rose toward the 1,000-read "
    "count, and controls retained 0–3 calls from 60 reads upward. Shading: at 30 reads the JH4 background "
    "model could not be calibrated, and each Aicda−/− library received 12–29 calls, so that point is not "
    "interpretable. (B, D) Distinct calls divided by depth and analysed bases (× 10−5; geometric mean ± "
    "geometric SD; zeros set to one only for the geometric calculation). The metric is not monotonic: at low "
    "depth, mutations seen in a single read rarely pass the false-discovery threshold, so the numerator is "
    "under-counted. Dashed line, 500-read working depth. One realised subset was analysed at each "
    "sample/depth combination; these operational curves do not estimate sensitivity against a known truth set.")

# Supplementary Table S1 ONT rows: reads analysed, candidates tested, calls, density (x1e-5)
S1_ONT = {"Wt2": ("500", "422", "155", "58.8"), "Wt3": ("500", "578", "347", "133.5"),
          "Wt4": ("500", "323", "69", "26.2"), "KO1": ("497", "243", "0", "0.00"),
          "KO2": ("500", "240", "3", "1.14"), "KO3": ("498", "244", "1", "0.38"),
          "1_C3": ("500", "206", "0", "0.00"), "3_C4": ("494", "206", "0", "0.00"),
          "4_D1": ("500", "168", "1", "0.33"), "6_C1": ("500", "307", "69", "22.8"),
          "7_D3": ("500", "271", "68", "22.5"), "10_D5": ("499", "393", "150", "49.7")}
S1_HEADER = {7: ("Candidates after ≥50% filter", "Candidates tested"),
             8: ("Filtered calls (p ≤ 0.05)", "Filtered calls (FDR ≤ 0.05)")}
# Supplementary Table S2: (old, new), searched across the table's cells
S2 = [
    ("p = 0.0037 (~47-fold)", "p = 0.0018 (≥ 108-fold; the control with no calls at one event)"),
    ("p = 0.0045 (≥ 91-fold; controls with no calls at one event)",
     "p = 0.0034 (≥ 89-fold; controls with no calls at one event)"),
    ("r = 0.95 [0.82–0.99]; slope = 1.05", "r = 0.91 [0.69–0.97]; slope = 1.02"),
    ("interaction p = 0.78; none significant (minimum adjusted p = 0.12)",
     "interaction p = 0.50; none significant (minimum adjusted p = 0.46)"),
    ("Ts adjusted p = 0.43; Tv adjusted p = 1.00", "Ts adjusted p = 0.10; Tv adjusted p = 0.97"),
    ("≈1.2 × 10−5 per analysed base (Q49-equivalent", "≈0.55 × 10−5 per analysed base (Q53-equivalent"),
]
# legend start -> regenerated figure
IMAGES = [
    ("Figure 1. Experimental design", FIG / "Figure1_study_design_pipeline.png"),
    ("Figure 2. Depth-standardised", FIG / "Figure2_JH4_combined.png"),
    ("Figure 3. ONT yields", FIG / "Figure3_concordance.png"),
    ("Figure 4. Filtered apparent variants", FIG / "Figure4_hybridoma_combined.png"),
    ("Supplementary Figure S1. Pooled substitution", FIG / "FigureS_spectrum_529_vs_fulllength.png"),
    ("Supplementary Figure S3. Individual JH4", FIG / "FigureS6_individual_JH4_ONT_matrices.png"),
    ("Supplementary Figure S4. Individual N114", FIG / "FigureS7_individual_N114_matrices.png"),
    ("Supplementary Figure S5. Recurrence-weighted", FIG / "FigureS7_recurrence_weighted.png"),
    ("Supplementary Figure S6. Sequencing depth", FIG / "FigureS3_depth_discovery_density.png"),
    ("Supplementary Figure S7. Filtered apparent-variant", FIG / "FigureS4_saturation.png"),
    ("Supplementary Figure S9. Sanger and ONT agree", SV / "fig3_agreement_position_hotspot.png"),
]
DROP_PARAS = ["N114 variant calling. N114 calls were made with SHMtrace"]   # merged into 2.8

# typography restored in replaced text: italic gene names, J_H4, superscripts
TOKENS = re.compile(r"(Aicda−/−|Aicda|JH4|\bWt(?=\d|\b)|10−5)")


def replace_text(paragraph, old, new):
    """Replace `old` (within the paragraph's direct runs) by `new`, highlighted yellow.

    Only the runs `old` touches are rewritten, in place, so hyperlinks and bookmarks
    elsewhere in the paragraph keep their position."""
    runs = paragraph.runs
    full = "".join(r.text for r in runs)
    i = full.find(old)
    if i < 0:
        raise RuntimeError(f"not found: {old[:70]!r} in {paragraph.text[:60]!r}")
    j = i + len(old)
    pos, hit = 0, []
    for r in runs:
        a, b = pos, pos + len(r.text)
        pos = b
        if b > i and a < j or (a == b == i):
            hit.append((a, b, r))
    if not hit:                                        # empty old at a run boundary
        hit = [next((a, b, r) for a, b, r in
                    ((sum(len(x.text) for x in runs[:k]), sum(len(x.text) for x in runs[:k + 1]), runs[k])
                     for k in range(len(runs))) if a <= i <= b)]
    first, last = hit[0], hit[-1]
    pieces = []
    if first[0] < i:
        pieces.append((first[2].text[: i - first[0]], first[2], False))
    pieces.append((new, first[2], True))
    if last[1] > j:
        pieces.append((last[2].text[j - last[0]:], last[2], False))
    anchor = first[2]._r
    for text, src, mark in pieces:
        if not text:
            continue
        el = copy.deepcopy(src._r)
        anchor.addprevious(el)
        r = Run(el, paragraph)
        r.text = text
        if mark:
            r.font.highlight_color = WD_COLOR_INDEX.YELLOW
            # new text takes the paragraph's plain style; restyle() re-applies
            # italics and sub/superscripts to gene names and exponents
            r.italic, r.font.subscript, r.font.superscript = None, None, None
    for _a, _b, r in hit:
        r._r.getparent().remove(r._r)


def _clone(run, text):
    el = copy.deepcopy(run._r)
    run._r.addprevious(el)
    r = Run(el, run._parent)
    r.text = text
    return r


def restyle(paragraph):
    """Split each highlighted run around gene names and exponents and format the pieces.

    A highlighted run that is italic as a whole (Ver17 inserted some text in the style
    of an adjacent gene name) is first returned to roman."""
    for run in paragraph.runs:
        if (run.font.highlight_color == WD_COLOR_INDEX.YELLOW and run.italic
                and not TOKENS.fullmatch(run.text.strip())):
            run.italic = None
    for run in list(paragraph.runs):
        if run.font.highlight_color != WD_COLOR_INDEX.YELLOW or not TOKENS.search(run.text):
            continue
        for piece in [x for x in TOKENS.split(run.text) if x]:
            if piece == "Aicda−/−":
                _clone(run, "Aicda").italic = True
                _clone(run, "−/−").font.superscript = True
            elif piece == "Aicda" or piece == "Wt":
                _clone(run, piece).italic = True
            elif piece == "JH4":
                _clone(run, "J").italic = True
                h = _clone(run, "H")
                h.italic, h.font.subscript = True, True
                _clone(run, "4").italic = True
            elif piece == "10−5":
                _clone(run, "10")
                _clone(run, "−5").font.superscript = True
            else:
                _clone(run, piece)
        run._r.getparent().remove(run._r)


def main(src, out):
    doc = Document(src)
    paras = doc.paragraphs

    def find(start):
        hits = [p for p in paras if p.text.startswith(start)]
        if len(hits) != 1:
            raise RuntimeError(f"{len(hits)} paragraphs start with {start!r}")
        return hits[0]

    touched = set()
    for start, old, new in EDITS:
        p = find(start)
        if new is None:
            continue
        replace_text(p, old or "".join(r.text for r in p.runs), new)
        touched.add(p._p)
    # S6 legend: keep the bold label run, replace the body
    s6 = find("Supplementary Figure S6. Sequencing depth")
    replace_text(s6, s6.text[len("Supplementary Figure S6. "):], S6_LEGEND_BODY)
    touched.add(s6._p)
    for start in DROP_PARAS:
        p = find(start)
        p._p.getparent().remove(p._p)

    s1 = doc.tables[1]
    for col, (old, new) in S1_HEADER.items():
        set_cell(s1.rows[0].cells[col], old, new)
    seen = set()
    for row in s1.rows:
        key = row.cells[0].text.strip()
        if key in S1_ONT and row.cells[2].text.strip() == "ONT":
            seen.add(key)
            for col, val in zip((5, 7, 8, 9), S1_ONT[key], strict=True):
                cell = row.cells[col]
                if cell.text.strip() != val:
                    set_cell(cell, cell.text.strip(), val)
    if seen != set(S1_ONT):
        raise RuntimeError(f"S1 rows missing: {set(S1_ONT) - seen}")
    t1 = [c for row in doc.tables[0].rows for c in row.cells
          if "observed post-filter Aicda−/− background ~1.2 × 10−5" in c.text]
    if not t1:
        raise RuntimeError("Table 1 background cell not found")
    set_cell(t1[0], "background ~1.2 × 10−5", "background ~0.55 × 10−5")
    s2 = doc.tables[2]
    for old, new in S2:
        cells = [c for row in s2.rows for c in row.cells if old in c.text]
        if not cells:
            raise RuntimeError(f"S2 text not found: {old!r}")
        set_cell(cells[0], old, new)

    for start, path in IMAGES:
        legend = find(start)
        idx = [p._p for p in doc.paragraphs].index(legend._p)
        img_p = next(doc.paragraphs[k] for k in range(idx - 1, idx - 4, -1)
                     if list(doc.paragraphs[k]._p.iter(qn("a:blip"))))
        swap_image(doc, img_p, path)

    for p in doc.paragraphs:
        restyle(p)
    for t in doc.tables:
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    restyle(p)
    doc.core_properties.comments = ("Ver18: all ONT results from SHMtrace (JH4 and N114, deposit sub500 "
                                    "subsets); Sanger as Ver17 (shmtrace docs/calibration/manuscript_corrected)")
    doc.save(out)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
