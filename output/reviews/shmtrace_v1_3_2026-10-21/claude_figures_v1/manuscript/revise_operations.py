"""Claude figure revision: rewrite the Codex manuscript's recorded document operations for the revised figures.

Reads codex_version/document_operations.json (never modified), applies the figure moves, citation remapping and
legend rewrites listed below, and writes claude_figures_v1/manuscript/document_operations.json. Numbers in the
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
ops[0][1][0] = "SHMtrace Claude figure revision"
ops[1][1][0] = "Version 1.3 Codex 1 + Claude figures 1 · 22 October 2026"

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
 "events/support point; vertical lines, range over the three controls. FDR labels are pooled false/(true+false) "
 "calls. (C) Permissive somatic callers on J_{H}4 with peer controls as matched normal: solid, the spike-ins of A "
 "(1–20 of 500 reads); dashed, clean plants at 5–20% of reads, 30 sites/fraction/control. Grey, cumulative share of "
 "the real SHMtrace calls of Figure 2E by allele fraction. (D) Recall against false calls/library, every caller at "
 "four loci; the x axis is linear from 0 to 0.1 and logarithmic above, so zero stays zero. BAQ, base-alignment "
 "quality. Lineage, quantification and group-test simulations are in Figure S8. See Figures S4–S8 and S24; Tables "
 "S1 and S6."),
"Figure 4.": ("Figure 4. Matched controls reduce apparent substitution background in real libraries",
 "(A–E) Unique variants/10^{5} analyzed bases per caller (rows): AID-competent libraries above each row's centre, "
 "AID-null below; shaded strip left of the axis break, libraries with no call. Axis ranges differ by panel. J_{H}4: "
 "three mice/genotype; N114: three subclones/arm, 500 reads/library. External Sμ/*Pim1*: 3–4 case and three null "
 "libraries, 2,000 reads each. Welch tests on log density, Holm-adjusted across five callers/panel; zero-call "
 "libraries enter at one event on their own exposure. ns, not significant; \\*p<0.05; \\*\\*p<0.01; "
 "\\*\\*\\*p<0.001. (F) Geometric-mean case ÷ control density, same zero convention; OT1/OT2, off-targets 1/2, "
 "with Wt and double-knockout groups separate (n=3 each). Colour is descriptive, not a test. (G) Calls summed over "
 "three controls/assay. (H) Observed held-out candidate tests with p<0.01 against the error margin; dotted line, 1% "
 "expected under calibration; rings, selected margin. Panels use peer-control caller-comparison calibration; "
 "off-target Results use production-project calibration (STAR Methods). Read accounting, Figure S1E. See Figures "
 "S9–S11; Tables S1 and S3."),
"Figure 5.": ("Figure 5. Matched Sanger and nanopore libraries yield similar mutation profiles",
 "Matched 529-bp data from three Wt and three null mice[@barbulescu2025]. (A) Within-platform densities: Sanger "
 "counts Q≥20 substitutions across 26–40 clones/Wt mouse; nanopore, called variants among 500 reads. The units "
 "differ, so each platform is read on its own. Bars, geometric means; p, log-Welch; nanopore zeros use one event. "
 "(B) Pooled twelve-type shares (471 Sanger/638 nanopore calls); dashed line, identity; r describes similarity of "
 "pattern, not interchangeable values. (C) Four classes by mouse: lines join each mouse's Sanger and nanopore "
 "shares; shapes, mice; shading, mean±SD. Platform×type testing uses all twelve types, matched by mouse. (D) "
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
 "Panels are grouped by question: burden (A, B), signature (C), indels (D), repair (E). (A) Illumina J_{H}4 "
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

(HERE / "document_operations.json").write_text(json.dumps(ops, indent=2, ensure_ascii=False))
(HERE / "citation_remap_log.json").write_text(json.dumps(log, indent=1, ensure_ascii=False))
print("operations:", len(ops), "; paragraphs remapped:", len(log))
