"""Claude figure revision 2 (v1 + Codex corrections + independent-review fixes): rewrite the Codex manuscript's recorded document operations for the revised figures.

Reads codex_version/document_operations.json (never modified), applies the figure moves, citation remapping and
legend rewrites listed below, and writes claude_figures_v2/manuscript/document_operations.json. Numbers in the
legends are the values already reported in the Codex text or computed by the figure scripts from crm/out; none is
changed.
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
REVIEW = BASE.parent
SRC = REVIEW / "codex_version/document_operations.json"
FIG = BASE / "figures"
OLD_FIG = str(REVIEW / "codex_version/figures")

ops = json.loads(SRC.read_text())


def walk(x, fn):
    if isinstance(x, str):
        return fn(x)
    if isinstance(x, list):
        return [walk(v, fn) for v in x]
    if isinstance(x, dict):
        return {k: walk(v, fn) for k, v in x.items()}
    return x


# 1. every image now comes from this revision's figures folder (unchanged figures were copied there)
ops = [[n, walk(a, lambda s: s.replace(OLD_FIG, str(FIG))), k] for n, a, k in ops]

# 2. citations of panels that moved. Placeholders first, so no remapping is applied twice.
MOVES = [("Figure 3E", "Figure S8A"), ("Figure 3F", "Figure S8B"), ("Figure 3G", "Figure S8C"),
         ("Figure S8A", "Figure S8D"), ("Figure 3H", "Figure S8E"), ("Figure S8B", "Figure S8F"),
         ("Figure S8C", "Figure S8G"), ("Figure S8D", "Figure S8H"), ("Figure 4I", "Figure S1E")]
log = []


def remap(s):
    out = s
    for i, (a, _b) in enumerate(MOVES):
        out = re.sub(re.escape(a) + r"(?![0-9A-Za-z])", f"@@{i}@@", out)
    for i, (_a, b) in enumerate(MOVES):
        out = out.replace(f"@@{i}@@", b)
    if out != s:
        log.append((s[:60], [m for m in MOVES if m[0] in s]))
    return out


# the old S8 legend is rewritten whole below, so remap everything except figure-block legends
for o in ops:
    if o[0] in ("para", "caption", "table", "heading"):
        o[1] = walk(o[1], remap)


def find(kind, start):
    idx = [i for i, o in enumerate(ops) if o[0] == kind and isinstance(o[1][0 if kind == "caption" else 1], str)
           and o[1][0 if kind == "caption" else 1].startswith(start)]
    assert len(idx) == 1, (kind, start, idx)
    return idx[0]


def replace_in_para(start, before, after):
    i = next(i for i, o in enumerate(ops) if o[0] == "para" and o[1][0].startswith(start))
    assert ops[i][1][0].count(before) == 1, before
    ops[i][1][0] = ops[i][1][0].replace(before, after)


# 3. small text edits that follow from the figures
replace_in_para("ClairS[@zheng2026]", "occurred in <1% of reads (Figure 2E)",
                "occurred in <1% of reads (Figures 2E and 3C)")
replace_in_para("An independent laboratory's SW_{HEL}", "(r = 1.0000 across 60 spectrum cells)",
                "(r > 0.9999 across 60 spectrum cells)")
ops[0][1][0] = "SHMtrace Claude figure revision 2"
ops[1][1][0] = "Version 1.3 Codex 1 + Claude figures 2 · 22 October 2026"

# 4. main-figure legends
LEG = {
"Figure 2.": ("Figure 2. Matched controls distinguish candidate substitutions from recurrent background",
 "(A, B) C-to-T support by strand at a hotspot (A) and a recurrent-background position (B) in Sμ WT3 and its "
 "three AID-null controls (PRJNA1481264; 2,000 reads/library), on one scale. Counts include all qualities; "
 "parentheses give Q≥20 subsets; p, SHMtrace's test. (C) Base qualities of the supporting reads, aligned under each "
 "site; bars, medians; dotted line, Q20. (D) Candidates kept and depth-scaled decoy-estimated false calls across "
 "p-value thresholds; the threshold is the largest at which the decoy estimate stays within 5% of the candidates "
 "kept. This empirical estimate assumes representative controls; it is not a guaranteed FDR. (E) Allele fractions "
 "of 1,197 retained calls in 12 AID-active libraries (500 reads for J_{H}4/N114 V_{H}; 2,000 for Sμ/*Pim1*); "
 "dotted line, default ≥50% exclusion. See Figures 3C and S2."),
"Figure 3.": ("Figure 3. Planted mutations quantify recall and false calls",
 "(A, B) Planted-substitution recall by supporting reads: 60 sites × 3 draws × 3 held-out controls/locus, 90 "
 "events/support point; vertical lines, range over the three controls. Observed FDP, pooled false/(true+false) "
 "calls. (C) Permissive somatic callers on J_{H}4 with peer controls as matched normal: solid, the spike-ins of A "
 "(1–20 of 500 reads); dashed, clean plants at 5–20% of reads, 30 sites/fraction/control. Grey, cumulative share of "
 "the real SHMtrace calls of Figure 2E by allele fraction. (D) Recall against false calls/library, every caller at "
 "four loci; the x axis is linear from 0 to 0.1 and logarithmic above, so zero stays zero. BAQ, base-alignment "
 "quality. See Figures S4–S8 and S24; Tables "
 "S1 and S6."),
"Figure 4.": ("Figure 4. Matched controls reduce apparent substitution background in real libraries",
 "(A–E) Unique variants/10^{5} analyzed bases per caller: blue, AID-competent and red, AID-null libraries; "
 "shaded strip below the axis break, libraries with no call. Axis ranges differ by panel. J_{H}4: "
 "three mice/genotype; N114: three subclones/arm, 500 reads/library. External Sμ/*Pim1*: 3–4 case and three null "
 "libraries, 2,000 reads each. Welch tests on log density, Holm-adjusted across five callers/panel; zero-call "
 "libraries enter at one event on their own exposure. ns, not significant; \\*p<0.05; \\*\\*p<0.01; "
 "\\*\\*\\*p<0.001; with three libraries per group, ns does not show a caller fails. (F) Geometric-mean case ÷ control density, same zero convention; CH12 OT1/OT2, AID off-target loci 1/2 in CH12 cells, "
 "with Wt and double-knockout groups separate (n=3 each). Colour is descriptive, not a test. (G) Calls summed over "
 "three controls/assay. (H) Held-out candidate tests with p<0.01 by error margin; dotted line, 1% expected if "
 "calibrated; rings, selected margin. Panels use peer-control caller-comparison calibration; "
 "off-target Results use production-project calibration (STAR Methods). See Figures "
 "S9–S11; Tables S1 and S3."),
"Figure 5.": ("Figure 5. Matched Sanger and nanopore libraries yield similar mutation profiles",
 "Matched 529-bp data from three Wt and three AID-null mice[@barbulescu2025]. (A) Within-platform unique-variant "
 "densities: distinct Q≥20 changes among 26–40 clones/Wt mouse (Sanger) or control-calibrated calls among 500 reads "
 "(nanopore), so each platform is read on its own. Bars, geometric means; p, log-Welch; nanopore zeros use one event. "
 "(B) Pooled twelve-type shares (471 Sanger/638 nanopore calls); dashed line, identity; r describes similarity of "
 "pattern, not interchangeable values. (C) Four classes by mouse: bars, mean±SD; shapes, mice. "
 "Platform×type testing uses all twelve types, matched by mouse. (D) "
 "Positional shares in 25-bp bins; the shaded final bin is 4 bp. (E) WRC/GYW use among G:C calls and WA/TW among "
 "A:T calls; bars, mouse means; points, mice; dashed lines, share of eligible positions in each motif. See Figures "
 "S13–S16; Tables S4 and S5."),
"Figure 6.": ("Figure 6. Long reads retain positional profiles and linked mutation patterns",
 "Libraries as Figure 5. (A) J_{H}4 variants/10^{5} bases in 50-bp bins, positions 90–1,230; thin lines, Wt mice; "
 "thick lines, genotype means, n=3 each. Shading, Sanger window; arrow, interval read only by nanopore. (B) Calls per "
 "position by S5F mutability quintile (~106 positions each), divided by each series' mean: pooled Wt calls versus "
 "germline-excluded null mismatches. The five points are bins of one window, not replicates. Targeting index, mean "
 "mutability at called positions ÷ window mean; permutation tests. (C) Junction-defined clones among 500 "
 "reads/mouse: dark, largest; light, others; hatching, junction unreadable. Right, read-weighted (open) and "
 "majority-rule clone-weighted (filled) frequency. (D) Largest Wt3 clone: square, inferred ancestor; ticks, "
 "mutations by class; numbers, reads carrying every mutation to that node; distance, accumulated mutations. 48 of 69 "
 "reads fit one path; 21 carry mutations of two branches and count at each. This is a read reconstruction, not "
 "proven cellular ancestry. See Figure S17."),
"Figure 7.": ("Figure 7. Public Illumina datasets test analysis across biological contexts",
 "(A) Illumina J_{H}4 "
 "frequencies (PRJNA973112); points, mice; strip below the axis break, no call; bars, geometric means with zeros at "
 "one event on each library's exposure. Tests: log-Welch versus null; Mann–Whitney for *Ldha*/Wt. (B) "
 "*Aicda*^{−/−}*Ung*^{−/−}*Msh2*^{−/−} CH12F3 expressing AID variants or empty vector (PRJNA918596), n=4/arm; "
 "display as A. (C) Substitution classes of each arm's pooled calls; numbers in bars are calls; S5F indices above. "
 "(D) Ramos indels called per clone at 20,000 reads each (PRJNA820976); n=6 Wt/HMCES-KO and 3 AID-KO; bars, means. "
 "(E) mOrange repair arms (PRJNA170048): circles, SHMtrace calls; squares, active-AID raw Q>20 mismatches; crosses, "
 "inactive-AID raw mismatches; lines join one experiment's three arms. Empty-arm inactive-AID libraries calibrate "
 "calling; 4–40 calls/library give coarse shares. See Table S1."),
}
for key, (title, legend) in LEG.items():
    i = find("caption", key)
    ops[i][1] = [title, legend]

# 5. Figure S1 gains the read accounting (was Figure 4I) as a second page
i = find("figure_block", "Figure S1.")
paths, title, legend = ops[i][1]
legend = legend.replace(" What became of every delivered read is in Figure 4I.", "")
legend = legend.replace(" Experimental details:", " (E) Fate of every delivered read of the J_{H}4 and N114 V_{H} "
    "libraries through SHMtrace's read preparation, as a share of delivered reads; numbers above bars, reads "
    "delivered. Experimental details:", 1)
ops[i][1] = [paths + [str(FIG / "FigureS1_read_fates.png")], title, legend]

# 6. Figure S8: downstream simulations (was Figure 3E-3H and S8A-B) on page 1, UMI consensus on page 2
i = find("figure_block", "Figure S8.")
old = ops[i][1][2]
b_txt = old[old.index("(B) ") + 4:old.index(" (C) ")]
c_txt = old[old.index("(C) ") + 4:old.index(" (D) ")].replace("Figure 3E", "Figure S8A")
d_txt = old[old.index("(D) ") + 4:]
page1 = ("What the read-library route's estimates and group comparisons deliver once calls are made. (A) Recall by read support in simulated J_{H}4 libraries of "
    "approximately 40 clones (power-law clone sizes; S5F-weighted mutations shared by a clone, a branch or one read), "
    "three libraries per scenario; vertical lines, range over the three libraries; text, false calls per library in "
    "all-null libraries. (B) Read- and clone-weighted estimates against their own truth; bars, mean estimates; "
    "black ticks, mean truth; stems join each library's truth to its estimate, because the truth differs between "
    "libraries. In the expanded scenario an unmutated clone contributes 60% of reads. (C) Estimate ÷ truth against "
    "depth, means and ranges over three held-out controls; solid, J_{H}4 polyclonal; dashed, N114 V_{H} subclonal. "
    "Combined density adds estimated uncalled single-read events. (D) Change in the aggregate single-read estimate "
    "against one-read mutations planted into held-out J_{H}4 control libraries (three controls; three replicates "
    "per amount); dashed line, identity. (E) Share of 300 simulated experiments per condition with p < 0.05, three "
    "or five animals per group, between-animal coefficient of variation 1. The two blocks answer different "
    "questions: with no true difference the rates are false-positive rates (dotted line, the nominal 5%); with a "
    "two- or four-fold difference they are power. Pooled-read tests treat reads rather than animals as the unit. "
    "Monte Carlo 95% half-width ≤ 5.6 points; the exact two-sided permutation test cannot fall below p = 0.10 at "
    "n = 3. (F) " + b_txt)
page2 = ("(G) " + c_txt + " (H) " + d_txt).replace(
    "which is why it sits beside the simulation in C", "which is why it sits beside the simulation in G")
blk1 = ["figure_block", [[str(FIG / "FigureS8_downstream.png")],
        "Figure S8. Downstream estimates and group tests in simulation, related to Figure 3 (page 1 of 2)", page1], {}]
blk2 = ["figure_block", [[str(FIG / "FigureS8_umi.png")],
        "Figure S8. UMI consensus in simulation and on real reads, related to Figure 3 (page 2 of 2)", page2], {}]
ops[i:i + 1] = [blk1, ["page_break", [], {}], blk2]

# 7. Figure S16: precision and the table's name follow the rebuilt figure
i = find("figure_block", "Figure S16.")
ops[i][1][2] = ops[i][1][2].replace("(D) Each published claim against what SHMtrace finds",
                                    "(D) Published and reanalysed effects, side by side")
assert "Published and reanalysed" in ops[i][1][2]

# a pre-existing markup defect: the allele name's asterisk opened an italic run over the rest of the legend
i = find("figure_block", "Figure S12.")
assert ops[i][1][2].count("IGHV1-47-36*01") == 1
ops[i][1][2] = ops[i][1][2].replace("IGHV1-47-36*01", "IGHV1-47-36\\*01")

# 7b. supplementary figures ported from the Codex v2 revision (S9, S11, S13, S17 page 1, S20B): legends follow the
# new displays; the experimental detail of the original legends is kept
LEGEND_EDITS = {
 "Figure S9. N114 V-region results": [
  ("(A) Unique-variant density in each subclone (n = 3 per vector). An empty-vector subclone with no call is drawn at "
   "one event on its own exposure, because a density of exactly zero cannot be placed on a logarithmic axis.",
   "(A) Unique-variant density in each subclone (n = 3 per vector), on a linear axis so that a subclone with no call "
   "sits at zero."),
  ("Each cell gives the unique-variant count and its share of that group's calls.",
   "Each cell gives the unique-variant count; both matrices share one colour scale. The single empty-vector call is "
   "shown, but a spectrum of one call is not interpreted.")],
 "Figure S9. Per-library spectra": [
  ("and each cell gives the unique-variant count and its share of that library's calls.",
   "and each cell gives the unique-variant count, on one colour scale per page. A library with no call has no "
   "spectrum, and totals of one or two calls are descriptive only.")],
 "Figure S11.": [
  ("in each AID-null library (red) and averaged over the AID-active libraries (blue); J_{H}4, 500 reads per library; "
   "Sμ (PRJNA1481264), 2,000. Insets, the controls alone on their own scale, on which the small recurrent artifacts "
   "are visible.",
   "averaged over the AID-active libraries (blue, top track, its own scale) and in each AID-null library (red, one "
   "track each, one shared scale per locus, on which the small recurrent artifacts are visible); J_{H}4, 500 reads "
   "per library; Sμ (PRJNA1481264), 2,000.")],
 "Figure S13.": [
  ("Top, per-mouse Sanger substitution matrices at Q≥20.",
   "Top, per-mouse counts of distinct Sanger substitutions at Q≥20, From (rows) → To (columns), on one colour "
   "scale; the AID-null totals of one or two are descriptive only.")],
 "Figure S17. Full-length spectra": [
  ("with each cell giving the mutation count and its share of all mutations; the two sets of proportions are highly "
   "correlated (Pearson r = 0.989).",
   "with each cell giving the mutation count (A, B), and the two sets of shares plotted against each other (C; "
   "Pearson r = 0.989, r² = 0.977). The window lies inside the interval, so the comparison is descriptive and not an independent "
   "validation.")],
 "Figure S20.": [
  ("(B) Unique-junction microhomology classes, averaged across three mice: blunt, 1–3 nt, ≥4 nt and junctions with "
   "inserted bases.",
   "(B) Unique-junction microhomology classes in each mouse: blunt, 1–3 nt, ≥4 nt and junctions with inserted "
   "bases.")],
}
for key, reps in LEGEND_EDITS.items():
    i = find("figure_block", key)
    for a, b in reps:
        assert ops[i][1][2].count(a) == 1, (key, a[:60])
        ops[i][1][2] = ops[i][1][2].replace(a, b)

# 7c. text corrections from the independent review (CLAUDE_INDEPENDENT_REVIEW_OF_CODEX_WORK.md, H2-H8)
replace_in_para("Cloned-amplicon Sanger sequencing", "recalled at measured false-discovery rates of 0.050 and 0.045",
                "recalled at observed false-discovery proportions of 0.050 and 0.045")
replace_in_para("We planted 60 AID-weighted", "(observed FDR 0.050)", "(observed false-discovery proportion, FDP, 0.050)")
replace_in_para("We planted 60 AID-weighted", "(FDR 0.083)", "(FDP 0.083)")
for before in ("(FDR 0.045;", "(FDR 0.050);", "at FDR 0.055 and 0.035", "at FDR 0.114 and 0.050", "permissive FDR 0.27"):
    replace_in_para("On N114 V_{H}, SHMtrace recalled", before, before.replace("FDR", "FDP"))
replace_in_para("We next simulated approximately 40", "with FDR 0.037 versus 0.012", "with FDP 0.037 versus 0.012")
replace_in_para("Seen this way, comparator performance", "at observed FDR 0.18 and 0.16", "at observed FDP 0.18 and 0.16")
replace_in_para("SHMtrace retained 1–4 calls", "reached 90% and 59% of discoveries at the largest tested depth",
                "reached 90% and 59% of the most each library yielded at any depth")
replace_in_para("SHMtrace processed 182 Sanger",
                "Absolute unique-variant densities are not directly comparable between workflows because the Sanger "
                "analysis counts mutations in individual cloned products, whereas the nanopore analysis counts distinct "
                "variants among 500 reads.",
                "Both densities count distinct variants but are not directly comparable: Sanger counts every Q≥20 "
                "change in 26–40 single clones, nanopore counts control-calibrated calls among 500 reads that "
                "resample expanded clones.")
replace_in_para("A shared framework accommodates", "and the same AID-null contrast",
                "and separated the genotypes in the same direction")
replace_in_para("The gap this fills is specific",
                "counting every mismatch to the germline separated AID from empty-vector N114 libraries by 2-fold where "
                "SHMtrace separated them by 188-fold",
                "counting every mismatch on the same N114 reads separated AID from empty-vector libraries by 2- to "
                "32-fold, depending on germline and masking, where SHMtrace separated them by 188-fold")
replace_in_para("Reanalysis of public data tested", "r = 1.0000", "r > 0.9999")
replace_in_para("What this enables in practice",
                "would require approximately twenty, with considerable uncertainty in that estimate",
                "would require about 24, and more if three mice underestimate the variation")
# STAR Methods (not counted): restore two qualifications the Codex compression dropped
replace_in_para("The margin c is chosen by leaving one control out at a time.",
                "The margin c is chosen by leaving one control out at a time.",
                "The margin c is chosen by leaving one control out at a time. It is an empirical calibration factor: it "
                "does not imply that sequencing errors are five or ten times the fitted rate, and it does not by itself "
                "set the false-discovery level, which the decoy threshold addresses separately. Because the same few "
                "controls choose the margin and supply the decoys, the false-discovery estimate is conditional on how "
                "well they represent the samples' background; it is not a guarantee.")

# 7d. the two ablations of out/spikein_*_ablation.json (Table S1, Spike-ins sheet), run on the same planted libraries
replace_in_para("On N114 V_{H}, SHMtrace recalled",
                "LoFreq recalled 0–26%, with permissive FDP 0.27–0.52 (Figure 3D).",
                "LoFreq recalled 0–26%, with permissive FDP 0.27–0.52 (Figure 3D). Two ablations on the same libraries "
                "isolated what the controls contribute: without them, a read-quality test made 5.7–46.6 false calls "
                "per library; testing each site against its control error rate made 1.3–2.8 but recalled 27–64%, "
                "against SHMtrace's 63–90% (Table S1).")
i = next(i for i, o in enumerate(ops) if o[0] == "para" and o[1][0].startswith("Table S1. Caller comparison"))
assert ops[i][1][0].count("for all tested callers in Figure 3 (Excel)") == 1
ops[i][1][0] = ops[i][1][0].replace("for all tested callers in Figure 3 (Excel)",
                                    "for all tested callers and two SHMtrace ablations (Excel)")
replace_in_para("For each library, the reads analyzed by SHMtrace were written out",
                "Three replicates were performed for each held-out control.",
                "Three replicates were performed for each held-out control. Two ablations used the same planted "
                "libraries (`crm/spikein.py --profiles ctlnull_bh,qualonly_pbinom --tag _ablation`; Table S1): "
                "ctlnull_bh tests each site against its per-position error rate in the other controls with "
                "Benjamini–Hochberg correction, without quality calibration or the decoy threshold, and qualonly_pbinom "
                "applies the exact Poisson-binomial test to read qualities with no control null. In these runs the "
                "recommended mode reproduced the original spike-in results exactly.")

# 7e. SciCanvas figures revised (Figure 1, S1 page 2, S3, S21, S23; scicanvas/ in this revision)
i = find("caption", "Figure 1.")
a = "undergo a control-fitted test and empirical false-discovery threshold;"
assert ops[i][1][1].count(a) == 1
ops[i][1][1] = ops[i][1][1].replace(a, "undergo a control-fitted test and an empirical, decoy-based threshold;")
SC_EDITS = {
 "Figure S1.": [
  ("numbers above bars, reads delivered. Experimental details:",
   "numbers above bars, reads delivered. (F) Identity to the J_{H}4 target and to its paralog after competitive "
   "alignment, for the first 400 reads of the bundled Wt3 example library (every second target-assigned read "
   "drawn); each read is assigned to whichever it matches better. Experimental details:")],
 "Figure S3.": [
  (" with a note on how to read it and what to watch for; the other figures of each group are listed by file name.",
   " with a note on how to read it and what to watch for; the further figures of each group are listed in the "
   "software documentation.")],
 "Figure S21.": [
  ("Bottom, measurements obtained from the two workflows in this study and their agreement in the six matched mice.",
   "Bottom, a table setting the two in-house platforms side by side: the unit each counts, its span, its calling "
   "rule, what each gives alone, the Wt/AID-null separation and the main limits, then their agreement in the six "
   "matched mice.")],
}
for key, reps in SC_EDITS.items():
    i = find("figure_block", key)
    for a, b in reps:
        assert ops[i][1][2].count(a) == 1, (key, a[:50])
        ops[i][1][2] = ops[i][1][2].replace(a, b)
i = find("figure_block", "Figure S23.")
ops[i][1][1] = "Figure S23. Rules underlying variant calling, related to STAR Methods"
ops[i][1][2] = ("The four statistical stages of the read-library route, each with the rule and the equation the software "
    "implements, and each shown working on the J_{H}4 libraries rather than stated in the abstract. Read preparation "
    "and competitive alignment are routine preprocessing and are shown in Figure S1. (A) Reported base quality "
    "against the mismatch rate measured in *Aicda*^{−/−} libraries; the dashed line is the error probability the "
    "quality score nominally implies, and the gap between them is what the margin corrects. (B) Two candidates at "
    "equal depth in the same library, with their supporting reads in the wild-type and control libraries, and why "
    "one is called and the other is not. (C) Candidate counts and held-out-decoy counts across p-value cutoffs, "
    "with the threshold that selection returns; the threshold is an empirical estimate that depends on the controls "
    "representing the samples' background. (D) The per-animal unique-variant density that group comparisons are "
    "made on. Experimental details: panels are drawn from `intuition_jh4.json` and the bundled example library.")

# 8. supplementary legends for figures whose display changed
EDITS = {
 "Figure S4.": [("dashed line, agreement; shading, over target; a pooled proportion of zero is drawn on the 0 row.",
                 "dashed line, realised = target; shading, over target; rings, the default 5%. The y axis is "
                 "logarithmic; a pooled proportion of zero has no logarithm and is drawn on its own 0 row below an "
                 "axis break.")],
 "Figure S5.": [("(A) Recall of planted mutations by number of controls.",
                 "(A) Recall of planted mutations by number of controls; the shaded one-control condition uses a "
                 "different procedure (below) and is set apart from the decoy-calibrated designs."),
                ("(C) One control carrying recurrent artifacts",
                 "(C) All controls clean, or one control carrying recurrent artifacts")],
 "Figure S7.": [("Recall below the three-read floor is zero by definition rather than by measurement (`indel_min_reads`), "
                 "so those counts are drawn open, on the shading, and the line starts at three reads instead of "
                 "running through them.",
                 "Below three reads the calling rule (`indel_min_reads`) calls no indel, so that region is shaded as "
                 "outside the rule and no recall is drawn there: it is not a measured failure.")],
 "Figure S12.": [("(A) Mutation frequency per library; filled, pMX-*Aicda*; open, empty vector; triangles, libraries "
                  "with no SHMtrace call, drawn at one event; bars, geometric means.",
                  "(A) Mutation frequency per library, the primary comparison; filled, pMX-*Aicda*; open, empty "
                  "vector; libraries with no SHMtrace call sit in the 0 strip below the axis break; bars, geometric "
                  "means with those zeros at one event. Beneath each route, its reference, the bases it keeps and "
                  "what it counts. (B) Secondary: geometric-mean separation; a near-zero denominator inflates it."),
                 ("(B) Geometric-mean separation of pMX-*Aicda* from empty-vector libraries. ", "")],
 "Figure S24.": [("(A) All-null resamples with any call and 95% Clopper–Pearson intervals;",
                  "Two outcomes are distinguished: any discovery in an all-null resample (A), and the false-discovery "
                  "proportion of each mutation-bearing resample, its false calls over all its calls (C). The shaded "
                  "case-only artifact, absent from every control, marks the boundary of control-based calibration. "
                  "(A) All-null resamples with any call and 95% Clopper–Pearson intervals;")],
}
for key, reps in EDITS.items():
    i = find("figure_block", key)
    for a, b in reps:
        assert ops[i][1][2].count(a) == 1, (key, a[:50])
        ops[i][1][2] = ops[i][1][2].replace(a, b)

# 7f. restore the flow of the main text and tighten the main captions (flow_edits.py)
import importlib.util as _ilu  # noqa: E402
_spec = _ilu.spec_from_file_location("flow_edits", HERE / "flow_edits.py")
FE = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(FE)
_NUM = re.compile(r"\d+(?:[.,]\d+)*")
_canon = json.loads((REVIEW / "filled_blocks.json").read_text())
_known = set(_NUM.findall(json.dumps(_canon, ensure_ascii=False)))
for o in ops:
    _known |= set(_NUM.findall(json.dumps(o[1], ensure_ascii=False)))
for start, text in FE.PARAS.items():
    idx = [i for i, o in enumerate(ops) if o[0] == "para" and o[1][0].startswith(start)]
    assert len(idx) == 1, (start, idx)
    new_nums = set(_NUM.findall(text)) - _known
    assert not new_nums, (start, new_nums)
    ops[idx[0]][1][0] = text
for key, text in FE.CAPTIONS.items():
    i = find("caption", key)
    new_nums = set(_NUM.findall(text)) - _known
    assert not new_nums, (key, new_nums)
    ops[i][1][1] = text

# 7g. methods detail moved verbatim from Results to STAR Methods (STAR Methods does not count toward the limit).
# Each entry: (Results paragraph opening, exact text removed, STAR Methods paragraph opening or None if STAR Methods
# already states it, text appended there).
MOVES_TO_METHODS = [
 ("The central step is the variant test",
  " We preferred it to Benjamini–Hochberg correction[@benjamini1995] because candidates arise only where a read already "
  "differs, so the tested set is itself chosen by the data.", None, None),
 ("SHMtrace starts from the files",
  " Reads that align better to a paralog, are chimeric or carry paralog-specific bases are removed, and the same number "
  "of usable reads is then drawn from every library.", "SHMtrace reads FASTQ files",
  " Reads that align better to a paralog, are chimeric or carry paralog-specific bases are removed, and the same number "
  "of usable reads is then drawn from every library."),
 ("Read preparation mattered as well", "WHOLE", "SHMtrace reads FASTQ files",
  " Read preparation removed 67% of delivered J_{H}4 reads overall (Figure S1E), and competitive alignment alone was "
  "not enough: without an additional paralog-base filter, simulated paralog reads produced 35–54 false variants per "
  "library."),
 ("We next simulated approximately 40",
  " Mutations were introduced at positions sampled according to S5F mutability, including mutations shared across an "
  "entire clone, restricted to individual branches, or unique to individual sequences (Figure S8A).", None, None),
 ("We next simulated approximately 40",
  " Junction clustering recovered clone assignments with adjusted Rand indices[@hubert1985] of 0.94–0.98.",
  "Lineage simulations partitioned one null library",
  " Junction clustering recovered clone assignments with adjusted Rand indices[@hubert1985] of 0.94–0.98."),
 ("At 500 reads per library",
  ", and in a separate test the estimate rose by 95% of the single-read mutations planted (Figure S8D)", None, None),
 ("Because both of those studies come from",
  " One host that may pool several animals was retained, and excluding it did not change the conclusions (STAR "
  "Methods).", None, None),
 ("Finally, PacBio",
  " Because reads per unique junction vary between samples, counts reflect PCR recovery rather than switching "
  "efficiency, and 19–47% of junctions place equally well at more than one acceptor position (STAR Methods).",
  "PRJNA831666[@vincendeau2022] holds",
  " Because reads per unique junction vary between samples, counts reflect PCR recovery rather than switching "
  "efficiency, and 19–47% of unique junctions place equally well at more than one acceptor position."),
 ("Long reads also keep each molecule's junction",
  " Overall, 70% of mutated reads fitted one tree, and 7/1,128 mutation pairs failed the four-gamete test.",
  "Junction grouping uses the twenty bases",
  " Overall, 70% of mutated reads fitted one tree, and 7/1,128 mutation pairs failed the four-gamete test."),
 ("On dual-UMI R10 rRNA-operon reads", "WHOLE", "PRJNA973112 (J_{H}4, positions 53–122",
  " On these dual-UMI R10 rRNA-operon reads[@karst2021], consensus reduced error 57-fold (95% CI 35–99) to 0.17%, "
  "with residual errors concentrated in homopolymer indels; at ≥10 reads per family, substitution error approached one "
  "per 10,000 bases (Figure S8H). This is a technical consensus benchmark rather than an SHM experiment."),
 ("Two CH12 AID off-target loci", "; 0.024 on frequency", "Groups were compared with one value per library",
  " For the two CH12 off-target loci, Welch tests on mutation frequency gave p = 0.024 (off-target 1) and 0.030 "
  "(off-target 2)."),
 ("Two CH12 AID off-target loci", " and 0.030", None, None),
]
for start, cut, dest, text in MOVES_TO_METHODS:
    idx = [i for i, o in enumerate(ops) if o[0] == "para" and o[1][0].startswith(start)]
    assert len(idx) == 1, (start, idx)
    if cut == "WHOLE":
        ops[idx[0]][0] = "_moved"
    else:
        assert ops[idx[0]][1][0].count(cut) == 1, (start, cut[:40])
        ops[idx[0]][1][0] = ops[idx[0]][1][0].replace(cut, "")
    if dest:
        d = [i for i, o in enumerate(ops) if o[0] == "para" and o[1][0].startswith(dest)]
        assert len(d) == 1, (dest, d)
        ops[d[0]][1][0] = ops[d[0]][1][0].rstrip() + text
ops = [o for o in ops if o[0] != "_moved"]
replace_in_para("We next simulated approximately 40", "following a power-law distribution.",
                "following a power-law distribution (Figure S8A).")
# the single-read sentence keeps its figure citations after the clause moved
replace_in_para("At 500 reads per library", "restored 96–100% of the true density.",
                "restored 96–100% of the true density (Figures S8C and S8D).")

# 7h. author-chosen cuts (22 Oct): B, Discussion summary sentences; D1, the Clair3/DeepVariant sentence; C, caption
# phrases the panels now state themselves (direct labels, reading guides) or the Results already give
B_D1 = [
 ("Measuring somatic hypermutation has long meant",
  " Together, these results support a practical route from deeply sequenced amplicons to interpretable SHM "
  "measurements across biological replicates.", ""),
 ("Seen this way, comparator performance",
  " Together, these outputs make the sensitivity–background tradeoff explicit and help investigators choose a calling "
  "regime appropriate to their experiment.", ""),
 ("Seen this way, comparator performance",
  " Clair3 and DeepVariant found none of the planted mutations, as expected of callers built to genotype a sample "
  "against a reference rather than to find a variant carried by a few reads.", ""),
 ("Applying the same engine to other laboratories",
  " By retaining assignments, counting conventions and quality flags, SHMtrace makes legacy sequence collections "
  "easier to revisit, compare and reuse.", ""),
 ("What this enables in practice",
  " More broadly, by bringing control-calibrated mutation detection, whole-amplicon profiles and comparisons across "
  "biological replicates into one reproducible workflow, SHMtrace makes long-read SHM data more informative and "
  "existing sequence collections easier to reanalyse transparently, providing a practical foundation for studying how "
  "AID targeting and DNA repair shape antibody diversification.",
  " More broadly, SHMtrace offers a practical foundation for studying how AID targeting and DNA repair shape antibody "
  "diversification."),
]
for start, a, b in B_D1:
    replace_in_para(start, a, b)
C_CUTS = {
 "Figure 1.": [(" SHMtrace estimates background by substitution, strand, quality and position from matched AID-null "
                "controls.", ""),
               ("(B) Inputs and read preparation; Sanger traces are trimmed at Q≥20 instead of the expected-error gate "
                "and are not downsampled.", "(B) Inputs and read preparation.")],
 "Figure 2.": [("; dotted line, ≥50% exclusion", "")],
 "Figure 3.": [("; the x axis is linear from 0 to 0.1 and logarithmic above", "")],
 "Figure 4.": [(" Axis ranges differ by panel.", ""),
               (" J_{H}4 and N114: three mice or subclones per group at 500 reads; Sμ and *Pim1*: 3–4 case and three "
                "null libraries at 2,000 reads.", ""),
               (" Panels use the caller-comparison calibration (STAR Methods).", "")],
 "Figure 5.": [(" (471 Sanger and 638 nanopore calls)", ""), ("; bars, mean±SD; shapes, mice", "")],
 "Figure 6.": [("; thin lines, Wt mice; thick lines, genotype means. Shading, the Sanger window.", "."),
               ("; the five points are bins of one window, not replicates", ""),
               (" (dark, largest; hatched, junction unreadable)", ""),
               (": square, inferred ancestor; ticks, mutations; numbers, reads carrying every mutation up to that "
                "node.", "."),
               (" The tree reconstructs reads, not proven cellular ancestry.", "")],
 "Figure 7.": [("; numbers in bars are calls", ""),
               (": circles, SHMtrace calls; squares and crosses, raw Q>20 mismatches for active and inactive AID; lines "
                "join one experiment's three arms", "")],
}
for key, reps in C_CUTS.items():
    i = find("caption", key)
    for a, b in reps:
        assert ops[i][1][1].count(a) == 1, (key, a[:50])
        ops[i][1][1] = ops[i][1][1].replace(a, b)

# 7i. V2 for Alberto's review (22 Oct): the author's own edits to the draft .docx (diffed against a fresh build), the
# answer to the author's comment on the Mkln1 paragraph, Philip Barbulescu as an author, and the version label
replace_in_para("HMCES-deficient Ramos clones",
                "A rule developed after inspecting this dataset therefore excluded it from threshold selection, while "
                "retaining it in model fitting and reporting the exclusion.",
                "After inspecting this dataset, we developed a rule that excluded it from threshold selection while "
                "retaining it for model fitting and reporting the exclusion.")
replace_in_para("Applying the same engine to other laboratories",
                "turned out to test it harder than any internal benchmark",
                "turned out to be a harder test than any internal benchmark")
replace_in_para("During the preparation of this work", "After using this tool,", "After using these tools,")
# author comment "is this correct?": the numbers match the analysis outputs, but "similar agreement" overstated the
# knockout match (spectrum r = 0.62; Discussion lists the Mkln1 discrepancies)
replace_in_para("A second study from the same laboratory", "A second study from the same laboratory gave similar "
                "agreement.", "A second study from the same laboratory agreed on frequency but less closely on the "
                "knockout spectrum.")
# word cap: drop a minor clause (candidate D2)
replace_in_para("Two C-to-T candidates in an external", "; the remainder could not be distinguished confidently from "
                "background.", ".")
# Philip Barbulescu performed the laboratory's Sanger sequencing
i = next(i for i, o in enumerate(ops) if o[0] == "para" and o[1][0].startswith("Abdulmateen Aderinto^{1,3}"))
assert ops[i][1][0] == "Abdulmateen Aderinto^{1,3} and Alberto Martin^{1,2,*}"
ops[i][1][0] = "Abdulmateen Aderinto^{1,3}, Philip Barbulescu^{1,4} and Alberto Martin^{1,2,*}"
# Philip has since moved to the SGC (UHN-affiliated); Cell Press convention: affiliation where the work was done,
# plus a present-address note
j = next(j for j, o in enumerate(ops) if o[0] == "para" and o[1][0] == "^{3}Technical contact")
ops.insert(j + 1, ["para", ["^{4}Present address: Structural Genomics Consortium, Toronto, ON M5G 1L7, Canada"],
                   {"size": 10, "space_after": 0}])
for key, cut in (("Figure 4.", ", with Wt and double-knockout groups shown separately"),
                 ("Figure 2.", "; dotted line, Q20")):
    k = find("caption", key)
    assert ops[k][1][1].count(cut) == 1, (key, cut)
    ops[k][1][1] = ops[k][1][1].replace(cut, "")
replace_in_para("Conceptualization, A.A. and A.M.", "investigation, A.A.;", "investigation, A.A. and P.B.;")
i = find("caption", "Figure 1.")
assert ops[i][1][1].count(" See also Figure S1.") == 1
ops[i][1][1] = ops[i][1][1].replace(" See also Figure S1.", "")
ops[0][1][0] = "SHMtrace V2 – for Alberto review"
ops[1][1][0] = "22 October 2026"

# 7j. figure consistency pass (22 Oct): the J_H4 mice are Wt2-Wt4 and KO1-KO3 and their genotype Aicda-/- in every
# figure, as in the text (AID-null stays the term for controls in general); Figure 5A now reads Aicda-/-, so its
# legend does too. Figure S17 page 2 is redrawn at page width (src/figS17_trees.py), same trees.
k = find("caption", "Figure 5.")
assert ops[k][1][1].count("three Wt and three AID-null mice") == 1
ops[k][1][1] = ops[k][1][1].replace("three Wt and three AID-null mice", "three Wt and three *Aicda*^{−/−} mice")
# Figure S9A carries the separation again (85-fold, p = 0.0071, as Figure 4 and Table S1 report it); its legend says
# how the test treats the zero-call subclones the linear axis shows at zero
k = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and o[1][1].startswith("Figure S9. N114 V-region"))
old = "on a linear axis so that a subclone with no call sits at zero."
assert ops[k][1][2].count(old) == 1
ops[k][1][2] = ops[k][1][2].replace(old, old + " Fold change (geometric means) and p (Welch test on log density) "
                                    "enter zero-call subclones at one event, as in Figure 4.")
# page 1 redrawn in the style of Figures 5 and 7 (23 Oct): the two matrices became one spectrum panel
old = ("(B, C) Pooled From (rows) → To (columns) substitution matrices for the pMX-*Aicda* subclones (B) and the "
       "empty-vector subclones (C). Each cell gives the unique-variant count; both matrices share one colour scale. "
       "The single empty-vector call is shown, but a spectrum of one call is not interpreted.")
new = ("(B) Substitution spectrum of the pMX-*Aicda* subclones: bars, shares of the pooled calls, with the count of "
       "each type beneath it; points, each subclone. The single empty-vector call is named but not plotted, since "
       "one call is not a spectrum.")
assert ops[k][1][2].count(old) == 1
ops[k][1][2] = ops[k][1][2].replace(old, new)


def legend(image, old, new, title=None):
    """Swap one passage of the legend of the figure block that shows image (and its title, if given)."""
    i = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and any(image in f for f in o[1][0]))
    assert ops[i][1][2].count(old) == 1, (image, old)
    ops[i][1][2] = ops[i][1][2].replace(old, new)
    if title:
        assert ops[i][1][1].count(title[0]) == 1
        ops[i][1][1] = ops[i][1][1].replace(*title)
    return i


# S9 page 2 removed (23 Oct): its bars were the pooled Wt shares already shown in Figure S17A and Figure 5B, and the
# per-mouse pattern is Figure 5C's; S9 is now one page, and S13 and S17 page 1 are restyled as spectra
i = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and o[1][1].startswith("Figure S9. Per-library spectra"))
del ops[i]
i = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and o[1][1].startswith("Figure S9. N114 V-region"))
for part, old in ((1, " (page 1 of 2)"), (2, " Spectra of the individual libraries are shown on the following page.")):
    assert ops[i][1][part].count(old) == 1, old
    ops[i][1][part] = ops[i][1][part].replace(old, "")
# S13's top panel removed (23 Oct): its bars were Figure 5B's Sanger shares, type for type, and the per-mouse pattern
# is Figure 5C's; page 1 keeps the readable-window panel
legend("FigureS13_sanger_matrices.png",
       "Top, per-mouse counts of distinct Sanger substitutions at Q≥20, From (rows) → To (columns), on one colour "
       "scale; the AID-null totals of one or two are descriptive only. Bottom, fraction of Wt clones",
       "Fraction of Wt clones")
i = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and any("FigureS13_sanger_matrices" in f for f in o[1][0]))
ops[i][1][0] = [f for f in ops[i][1][0] if not f.endswith("FigureS13_sanger_matrices.png")]
# "(page 1 of 2)" counted the two images across two pages; one image is left, on one page
assert ops[i][1][1].count(" (page 1 of 2)") == 1
ops[i][1][1] = ops[i][1][1].replace(" (page 1 of 2)", "")
legend("FigureS17_spectra_full_length.png",
       "Pooled wild-type (*Wt*) substitution spectra across the shared 529-bp window and across the reliably covered "
       "1,141-bp interval, with each cell giving the mutation count (A, B), and the two sets of shares plotted "
       "against each other (C; Pearson r = 0.989, r² = 0.977).",
       "(A) Pooled wild-type (*Wt*) substitution spectra across the shared 529-bp window (lighter bars) and across "
       "the reliably covered 1,141-bp interval (darker bars), with both counts beneath each type. (B) The two sets "
       "of shares plotted against each other (Pearson r = 0.989, r² = 0.977).")

# S2 redrawn at full width with its two panels side by side (23 Oct): it was drawn 85 mm wide and placed at page
# width, nearly twice its size
legend("FigureS2.png", "Top, the mismatch rate", "(A) The mismatch rate")
legend("FigureS2.png", "Bottom, the number of control bases in each bin.", "(B) The number of control bases in each bin.")
# S12 gives SHMtrace 188-fold on N114 where Figures 4F and S9A give 85-fold: both hold, for different measures
legend("FigureS12_immcantation.png",
       "(B) Secondary: geometric-mean separation; a near-zero denominator inflates it.",
       "(B) Secondary: geometric-mean separation; a near-zero denominator inflates it. SHMtrace's 188-fold here is "
       "a ratio of mutation frequencies over the V segment; the 85-fold of Figures 4F and S9A is a ratio of "
       "unique-variant densities in the same libraries, which count each distinct substitution once.")

# 7k. the read accounting (S1E-F) becomes its own Figure S2 (23 Oct). It is first cited with S1, so it takes the
# next number and S2-S24 become S3-S25. Renumbering runs last, over every string, and touches only supplementary
# figure numbers in a list governed by "Figure"/"Figures": tables, the S5F model and other papers' "Supplementary
# Figure" panels are left alone. Every change is printed for checking.
CITE = re.compile(r"\b(Supplementary Figures?|Supplementary Fig\.|Figures?|Tables?|Document|Data)(\s+)"
                  r"((?:S?\d+[A-Za-z]?(?:–S?\d+[A-Za-z]?)?)(?:(?:,\s*|\s+and\s+|–)(?:S?\d+[A-Za-z]?"
                  r"(?:–S?\d+[A-Za-z]?)?))*)")


def bump(m):
    if m.group(1) not in ("Figure", "Figures"):
        return m.group(0)

    def one(t):
        n, panel = int(t.group(1)), t.group(2)
        if n == 1 and panel in ("E", "F"):
            return "S2" + {"E": "A", "F": "B"}[panel]
        return f"S{n + 1}{panel}" if n >= 2 else t.group(0)
    return m.group(1) + m.group(2) + re.sub(r"S(\d+)([A-Za-z]?)", one, m.group(3))


renumbered = []


def renumber(s):
    out = CITE.sub(bump, s)
    if out != s:
        renumbered.append((s, out))
    return out


ops = [[n, walk(a, renumber), k] for n, a, k in ops]
for before, after in renumbered:
    for x, y in zip(CITE.finditer(before), CITE.finditer(after), strict=True):
        if x.group(0) != y.group(0):
            print(f"  renumbered: {x.group(0)!r} -> {y.group(0)!r}")
# the S1 block keeps A-D (the workflow); E-F go to a new Figure S2 right after it
i = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and o[1][1].startswith("Figure S1."))
paths, title, text = ops[i][1]
cut = text.index(" (E) Fate of every delivered read")
details = text.index(" Experimental details:")
fates = [p for p in paths if p.endswith("FigureS1_read_fates.png")]
ops[i][1] = [[p for p in paths if p not in fates], title, text[:cut]]
s2 = text[cut + 1:details].replace("(E) Fate", "(A) Fate").replace("(F) Identity", "(B) Identity")
ops.insert(i + 1, ["figure_block", [fates, "Figure S2. Read accounting, related to Figure 1", s2 + text[details:]],
                   ops[i][2]])
# as part of S1 the read accounting was cited early; on its own it was first cited in STAR Methods, out of order.
# The Results sentence that opens read preparation now cites it (one word)
replace_in_para("SHMtrace starts from the files", "(Figures 1 and S1)", "(Figures 1, S1 and S2)")
# Figure S24 (was S23) points at read preparation and competitive alignment, now in S1 and S2
legend("FigureS23_engine.png", "are routine preprocessing and are shown in Figure S1.",
       "are routine preprocessing and are shown in Figures S1 and S2.")

# 7l. first citations in order at the end too (23 Oct): the stress tests (S25) are cited in the Figure 3 caption,
# before the overview (S23) and the engine (S24), which only STAR Methods cites. They become S23, S24 and S25,
# and the stress-test block moves ahead of the other two. A range that ends on one of them (the supplemental index,
# "Figures S1–S25") counts figures, so it is left alone.
SWAP = {23: 24, 24: 25, 25: 23}
swapped = []


def swap(m):
    if m.group(1) not in ("Figure", "Figures"):
        return m.group(0)

    def one(t):
        if m.group(3)[max(0, t.start() - 1):t.start()] == "–":
            return t.group(0)
        n = int(t.group(1))
        return f"S{SWAP.get(n, n)}{t.group(2)}"
    out = m.group(1) + m.group(2) + re.sub(r"S(\d+)([A-Za-z]?)", one, m.group(3))
    if out != m.group(0):
        swapped.append((m.group(0), out))
    return out


ops = [[n, walk(a, lambda s: CITE.sub(swap, s)), k] for n, a, k in ops]
for before, after in swapped:
    print(f"  reordered: {before!r} -> {after!r}")
blocks = {re.match(r"Figure S(\d+)\.", o[1][1]).group(1): i for i, o in enumerate(ops) if o[0] == "figure_block"}
stress = ops.pop(blocks["23"])
ops.insert(min(blocks["24"], blocks["25"]), stress)

# 7m. the engine figure (S25) shows all six stages again (23 Oct, author: the engine-specific figure is the place
# for the full detail). Read preparation and competitive alignment return as A and B; the statistical stages
# become C-F. The panels are the same as Figure S2's, drawn here as stages of the engine.
legend("FigureS23_engine.png",
       "The four statistical stages of the read-library route, each with the rule and the equation the software "
       "implements, and each shown working on the J_{H}4 libraries rather than stated in the abstract. Read "
       "preparation and competitive alignment are routine preprocessing and are shown in Figures S1 and S2. "
       "(A) Reported base quality",
       "The six stages of the read-library route, each with the rule and the equation the software implements, and "
       "each shown working on the J_{H}4 libraries rather than stated in the abstract. (A) Read accounting in each "
       "library: reads removed by the length and expected-error gates, reads that do not align best to the target, "
       "and the reads left usable (the individual filters are in Figure S2A). (B) Alignment identity to the target "
       "against identity to its paralog in the bundled example library; each read is assigned to whichever it "
       "matches better (Figure S2B). (C) Reported base quality",
       ("Rules underlying variant calling", "Rules underlying read preparation and variant calling"))
for old, new in (("(D) The per-animal", "(F) The per-animal"), ("(C) Candidate counts", "(E) Candidate counts"),
                 ("(B) Two candidates", "(D) Two candidates")):
    legend("FigureS23_engine.png", old, new)
legend("FigureS23_engine.png", "panels are drawn from `intuition_jh4.json` and the bundled example library.",
       "panels are drawn from `intuition_jh4.json`, the per-library read-accounting records and the bundled example "
       "library.")
# with six panels the legend ran onto a page of its own; each card states its rule, so the legend is tightened
# (same content, panel by panel) until figure and legend share the page
i = next(i for i, o in enumerate(ops) if o[0] == "figure_block" and any("S23_engine" in f for f in o[1][0]))
ops[i][1][2] = (
    "The six stages of the read-library route; each card gives the rule and equation the software implements, and "
    "its panel shows the stage on the J_{H}4 libraries. (A) Reads gated, not assigned to the target, "
    "and usable, per library; the same data as Figure S2A, which shows each filter. (B) Identity to the target "
    "against identity to its paralog in the bundled example library; the same data as Figure S2B. (C) Reported "
    "base "
    "quality against the mismatch rate in *Aicda*^{−/−} libraries; dashed line, the error the quality score "
    "implies; the gap is what the margin corrects. (D) Two candidates at equal depth, one called and one not. "
    "(E) Candidates and held-out-decoy counts across "
    "p-value cutoffs, with the selected threshold, an empirical estimate that depends on the controls "
    "representing the samples' background. (F) Per-animal unique-variant density, on which groups are compared. "
    "Experimental details: panels are drawn from `intuition_jh4.json`, the per-library read-accounting records and "
    "the bundled example library.")

# 7n. the Motivation's framing carried into the main text (author, 23 Oct): the gap is not only that cloning is
# slow, it is that the software built for SHM since has followed B-cell receptor repertoire sequencing, so a
# reporter locus read in bulk has had no equivalent route. Two sentences carry it, in the Summary and at the end of
# the Introduction's "two common analysis practices" paragraph, where the Discussion's gap paragraph then picks it
# up. "Has instead followed", not "is only for": MiXCR documents a nanopore preset. The trims that follow pay for
# the added words; none changes a number or a claim, and the Summary stays within its own 150-word limit.
MAIN_TEXT = [
    ("Existing SHM software summarises previously identified mutations[@maccarthy2009,@gupta2015,@yaari2013] "
     "but does not use matched AID-null controls to distinguish candidate variants from technical background.",
     "SHM software has instead followed B-cell receptor repertoire sequencing, summarising mutations already "
     "identified[@maccarthy2009,@gupta2015,@yaari2013] rather than using matched AID-null controls to separate "
     "candidate variants from technical background."),
    ("Cloned-amplicon Sanger sequencing provides accurate but low-throughput measurements of somatic hypermutation "
     "(SHM).",
     "Somatic hypermutation (SHM) is measured accurately but at low throughput by cloned-amplicon Sanger "
     "sequencing, while existing software targets B-cell receptor repertoires."),
    # paying for the two additions, in the same two sections
    ("Nanopore reads preserved mutation linkage across the amplicon and yielded substitution spectra similar to "
     "Sanger sequencing of the same mice.",
     "Nanopore reads preserved linkage and gave spectra similar to Sanger in the same mice."),
    ("twelve loci from eleven laboratories across four sequencing platforms. Across five nanopore assays",
     "twelve loci from eleven laboratories across four platforms. Across five nanopore assays"),
    ("libraries showed 66–141-fold separation in unique-variant density",
     "libraries separated 66–141-fold in unique-variant density"),
    ("Further datasets extended analysis to engineered diversification and class-switch junctions.",
     "Datasets covered engineered diversification and class-switch junctions."),
    ("mutations, recalled at observed false-discovery proportions of 0.050 and 0.045, simulated",
     "mutations (false-discovery proportions 0.050 and 0.045), simulated"),
    # the Summary describes the comparator instead of naming it (author, 23 Oct); the Results and Discussion
    # still name Geneious. "General-purpose ... not tuned for SHM" is what its own cited title describes: an
    # integrated desktop platform for the organization and analysis of sequence data[@kearse2012]
    ("by 64\u201390% relative to the uncorrected Geneious rule.",
     "by 64–90% relative to an uncorrected rule from a general-purpose sequence-analysis platform not "
     "tuned for SHM."),
    # "GUI" was dropped: SHMtrace ships a desktop app of its own, so the contrast is general-purpose versus
    # SHM-specific, not interface. Paid for by a redundant word in the Introduction.
    ("In the configurations evaluated here, neither measures",
     "In the configurations evaluated, neither measures"),
    ("under the stated zero-call convention", "under the zero-call convention"),
    ("and compare mutation profiles across biological replicates",
     "and compare profiles across biological replicates"),
    ("so mismatches that recur at one position", "so mismatches recurring at one position"),
    ("What all platforms, however, share is the need to separate AID-driven mutation from their own error.",
     "All platforms, however, share the need to separate AID-driven mutation from their own error."),
    ("biological replicates, starting from delivered reads or Sanger chromatograms",
     "biological replicates, from delivered reads or Sanger chromatograms"),
]
for _old, _new in MAIN_TEXT:
    _hit = 0
    for _o in ops:
        for _j, _s in enumerate(_o[1]):
            if isinstance(_s, str) and _old in _s:
                _o[1][_j] = _s.replace(_old, _new)
                _hit += 1
    assert _hit == 1, (_hit, _old[:60])

(HERE / "document_operations.json").write_text(json.dumps(ops, indent=2, ensure_ascii=False))
(HERE / "citation_remap_log.json").write_text(json.dumps(log, indent=1, ensure_ascii=False))
print("operations:", len(ops), "; paragraphs remapped:", len(log))
