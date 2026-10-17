"""Build the Cell Reports Methods article and its companion files.

    python docs/calibration/crm/manuscript/build_crm.py OUTDIR

Writes OUTDIR/SHMtrace_CellReportsMethods.docx (manuscript: title page, summary, keywords, introduction,
results, discussion, limitations, resource availability, declarations, supplemental index, figure legends
with figures, Table 1, STAR Methods with the key resources table, references), OUTDIR/Highlights_eTOC.docx
(highlights, eTOC blurb, motivation, as the journal asks: a separate file) and OUTDIR/Document_S1.docx.
Numbers for the public Illumina sets are read from out/public/*/project.json; the single-read
recovery from out/single_read_spike.json. Journal limits are checked and the build stops if any is
exceeded.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import crm_methods as M  # noqa: E402
import crm_numbers  # noqa: E402
import crm_text as T  # noqa: E402
from docbuild import Manuscript, word_count  # noqa: E402

CRM = Path(__file__).resolve().parents[1]
FIG = CRM / "figures"
OUT = CRM / "out"
CAL = CRM.parent
MC = CAL / "manuscript_corrected/figures"


def geo(x):
    return 10 ** np.mean(np.log10(x))


def public_numbers():
    n = {}
    rows = json.load(open(OUT / "single_read_spike.json"))
    n["sr_rec"] = f"{100 * np.mean([r['delta'] / r['planted'] for r in rows if r['planted']]):.0f}"
    pdeep = [r for r in json.load(open(OUT / "benchmark_real_pim1_deep.json"))["rows"]]
    rec = [r for r in pdeep if r["mode"] == "recommended"]
    gen = [r for r in pdeep if r["mode"] == "geneious"]
    rng = lambda rs, c: (min(r["calls"] for r in rs if r["control"] == c), max(r["calls"] for r in rs if r["control"] == c))  # noqa: E731
    reads = [r["reads"] for r in rec]
    a, b = rng(rec, False), rng(rec, True)
    ga, gb = rng(gen, False), rng(gen, True)
    n["pim1_deep"] = (f"With all usable reads ({min(reads):,}–{max(reads):,} per library), SHMtrace called "
                      f"{a[0]}–{a[1]} mutations per AID-competent library and "
                      f"{'none' if b[1] == 0 else f'{b[0]}–{b[1]}'} per control, and the Geneious rule called "
                      f"{ga[0]}–{ga[1]} and {gb[0]}–{gb[1]}, respectively. At a weak AID target, the signal therefore "
                      "remains close to the background even at the full depth available (Table S2).")
    # Illumina J_H4 (Ldha)
    p = json.load(open(OUT / "public/ldha/project.json"))
    f = {g: [s["mutation_frequency_x1e3"] for s in p["samples"] if s["group"] == g] for g in ("WT", "LDHA_KO", "AID_KO")}
    g = {x["group"]: x for x in p["groups"]}
    n["ldha_fold"] = f"{np.mean(f['WT']) / np.mean(f['AID_KO']):.0f}"
    pw = g['WT']['mutation_frequency_x1e3_p_welch_log_vs_AID_KO']
    n["ldha_p"] = crm_numbers.sci(pw) if pw < 0.001 else f"{pw:.3f}"
    n["ldha_nwt"], n["ldha_nko"] = len(f["WT"]), len(f["AID_KO"])
    from scipy import stats
    n["ldha_p_ko"] = f"{stats.mannwhitneyu(f['WT'], f['LDHA_KO']).pvalue:.2f}"
    # CH12F3
    p = json.load(open(OUT / "public/ch12/project.json"))
    t = {(x["group"], x["metric"]): x for x in p["animal_tests"]}
    n["ch12_fold"] = f"{t[('AID', 'mutation_frequency')]['fold']:.0f}"
    g = {x["group"]: x for x in p["groups"]}
    n["ch12_p"] = f"{g['AID']['mutation_frequency_x1e3_p_welch_log_vs_EV']:.5f}".rstrip("0")
    calls = [v for v in p["variants"]]
    aid = [v for v in calls if v["group"] == "AID"]
    mut = [v for v in calls if v["group"] == "R174ER178E"]
    ctga = lambda vs: 100 * sum(v["ref"] + v["alt"] in ("CT", "GA") for v in vs) / len(vs)   # noqa: E731
    n["ch12_ctga"] = f"{ctga(aid):.0f}"
    n["ch12_mut_ctga"] = f"{ctga(mut):.0f}"
    s5 = json.load(open(OUT / "targeting_ch12.json"))
    n["ch12_aid_s5f"] = f"{s5['AID']['index']:.2f}"
    n["ch12_mut_s5f"] = f"{s5['R174ER178E']['index']:.2f}"
    # HMCES
    p = json.load(open(OUT / "public/hmces/project.json"))
    s = {x["sample"]: x for x in p["samples"]}
    dels = {gname: [x["deletions"] for x in p["samples"] if x["group"] == gname] for gname in ("WT", "HMCES_KO")}
    ins = {gname: [x["insertions"] for x in p["samples"] if x["group"] == gname] for gname in ("WT", "HMCES_KO")}
    n["hmces_del_fold"] = f"{np.mean(dels['HMCES_KO']) / np.mean(dels['WT']):.1f}"
    n["hmces_ins_fold"] = f"{np.mean(ins['HMCES_KO']) / np.mean(ins['WT']):.1f}"
    n["hmces_del_p"] = crm_numbers.sci(stats.ttest_ind(np.log10(np.array(dels['HMCES_KO']) + 1), np.log10(np.array(dels['WT']) + 1), equal_var=False).pvalue)
    n["hmces_ko3_gc"] = f"{s['AIDKO_3']['pct_GC']:.0f}"
    return n


def fill(text, n):
    for k, v in n.items():
        text = text.replace("{" + k + "}", str(v))
    if re.search(r"(?<![\^_])\{[a-z0-9_]+\}", text):
        raise RuntimeError(f"unfilled placeholder in: {text[:80]}")
    return text


def fill_all(x, n):
    """Fill placeholders in a string or any nesting of lists and tuples of strings."""
    if isinstance(x, str):
        return fill(x, n)
    if isinstance(x, (list, tuple)):
        return type(x)(fill_all(y, n) for y in x)
    return x


def fill_modules(n):
    for mod in (T, M):
        for k in dir(mod):
            if k.isupper():
                setattr(mod, k, fill_all(getattr(mod, k), n))
    global STATS, SUPP_FIGS
    STATS, SUPP_FIGS = fill_all(STATS, n), fill_all(SUPP_FIGS, n)


def check_figure_numbers(n):
    """The graphical abstract is drawn in SciCanvas with its numbers typed in: they must match the analysis."""
    ga = (CRM / "scicanvas/crm_ga.js").read_text()
    for text in (f"{n['ctl_rec']} calls per control set", f"r = {n['r_spec']}"):
        if text not in ga:
            raise RuntimeError(f"graphical abstract is stale: expected '{text}' in scicanvas/crm_ga.js")


def check_limits():
    words = lambda s: len(s.split())   # noqa: E731
    sw = word_count([("p", T.SUMMARY)])
    assert sw <= 150, f"summary {sw} words"
    for h in T.HIGHLIGHTS:
        assert len(h) <= 85, f"highlight {len(h)} chars: {h}"
    assert words(T.ETOC) <= 50, f"eTOC {words(T.ETOC)} words"
    assert len(T.TITLE) <= 145, "title too long"
    return sw


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    n = {**crm_numbers.compute(), **public_numbers()}
    fill_modules(n)
    check_figure_numbers(n)
    sw = check_limits()
    m = Manuscript()
    # ---- title page
    m.para(T.TITLE, size=16, bold=True)
    m.para(T.AUTHORS)
    for a in T.AFFIL:
        m.para(a, size=10, space_after=0)
    m.heading("Summary", 1)
    m.para(T.SUMMARY)
    m.heading("Keywords", 1)
    m.para(T.KEYWORDS)
    # ---- main text
    main_blocks = []
    m.heading("Introduction", 1)
    for p in T.INTRO:
        m.para(p)
        main_blocks.append(("p", p))
    m.heading("Results", 1)
    for kind, text in T.RESULTS:
        text = fill(text, n)
        main_blocks.append((kind, text))
        if kind == "h":
            m.heading(text, 2)
        else:
            m.para(text)
    m.heading("Discussion", 1)
    for p in T.DISCUSSION:
        m.para(p)
        main_blocks.append(("p", p))
    m.heading("Limitations of the study", 2)
    m.para(T.LIMITATIONS)
    main_blocks.append(("p", T.LIMITATIONS))
    m.heading("Resource availability", 1)
    for h, t in T.RESOURCE:
        m.heading(h, 2)
        m.para(t)
    m.heading("Acknowledgments", 1)
    m.para(T.ACK)
    m.heading("Author contributions", 1)
    m.para(T.CONTRIB)
    m.heading("Declaration of interests", 1)
    m.para(T.INTERESTS)
    m.heading("Declaration of generative AI and AI-assisted technologies in the manuscript preparation process", 1)
    m.para(T.AI_DECL)
    m.heading("Supplemental information", 1)
    for s in M.SUPPLEMENT:
        m.para(s)
    # ---- figures
    m.page_break()
    m.heading("Figure titles and legends", 1)
    for i, (title, legend) in enumerate(M.LEGENDS, 1):
        m.figure(FIG / f"Figure{i}{'_preview' if i == 1 else ''}.png", width_mm=165)
        m.para(f"**{title}**")
        m.para(fill(legend, n), size=10)
    # ---- table
    m.page_break()
    m.heading("Tables", 1)
    m.table(M.TABLE1_HEADER, M.TABLE1, title=f"**{M.TABLE1_TITLE}**", footnote=M.TABLE1_NOTE,
            widths=[30, 22, 30, 36, 30, 22])
    # ---- STAR Methods
    m.page_break()
    m.heading("STAR★Methods", 1)
    m.heading("Key resources table", 2)
    m.table(M.KRT_HEADER, M.KRT, widths=[75, 50, 45], size=8)
    m.heading("Experimental model and study participant details", 2)
    for h, t in M.EXPERIMENTAL:
        m.heading(h, 3)
        m.para(t)
    m.heading("Method details", 2)
    for kind, t in M.METHODS:
        if kind == "h2":
            m.heading(t, 3)
        else:
            m.para(t)
    m.heading("Quantification and statistical analysis", 2)
    for p in M.QUANT:
        m.para(p)
    m.references()
    path = out / "SHMtrace_CellReportsMethods.docx"
    m.save(path)
    mw = word_count(main_blocks)
    print(f"{path}\nsummary {sw} words; main text (introduction to limitations) {mw} words; references {len(m.order)}")
    assert mw < 7000, "main text over 7,000 words"

    # ---- highlights, eTOC, motivation (separate file)
    h = Manuscript(line_numbers=False)
    h.heading("Highlights", 1)
    h.bullets(T.HIGHLIGHTS)
    h.heading("eTOC blurb", 1)
    h.para(T.ETOC)
    h.heading("Motivation", 1)
    h.para(T.MOTIVATION)
    h.save(out / "Highlights_eTOC_Motivation.docx")
    print("highlights:", [len(x) for x in T.HIGHLIGHTS], "eTOC words:", len(T.ETOC.split()))

    # ---- Document S1
    s = Manuscript(line_numbers=False)
    s.para("Supplemental information", size=14, bold=True)
    s.para(T.TITLE, size=11)
    for title, img, legend in SUPP_FIGS:
        s.page_break()
        for i in img if isinstance(img, list) else [img]:
            s.figure(i, width_mm=165)
        s.para(f"**{title}**")
        s.para(legend, size=10)
    # Table S1: what SHMtrace did with every delivered read, and the calls (J_H4 and N114)
    s.page_break()
    rows = []
    for run, label in (("jh4", "J_{H}4"), ("n114", "N114")):
        acct = json.load(open(CAL / f"figures/{run}_shmtrace/table_s1.json"))
        for lib, a in acct.items():
            rows.append([label, lib, f"{a['delivered']:,}", f"{a['too_long']:,}", f"{a['low_quality']:,}",
                         f"{a['not_target']:,}", f"{a['chimeric'] + a['short_alignment']:,}",
                         f"{a['paralogue_alleles']:,}", f"{a['usable']:,}", a["analysed"], a["candidates_tested"],
                         a["calls"]])
    s.table(["Assay", "Library", "Delivered", "Over-length", "Expected error > 0.01", "Unmapped or paralog",
             "Chimeric or short alignment", "Paralog alleles", "Usable", "Analyzed", "Candidates tested",
             "Calls (FDR ≤ 0.05)"], rows, size=7.5,
            title="**Table S1. Read preparation and calls per library, related to Figure 4**",
            footnote="Counts of reads removed at each step of SHMtrace's read preparation, from the delivered FASTQ "
                     "files; 500 usable reads were drawn at random (seed 0) for analysis. WT, wild type; KO, "
                     "*Aicda*^{−/−}; C1, D3, D5, pMX-*Aicda* subclones; C3, C4, D1, empty-vector subclones.")
    # Table S3: statistics
    s.page_break()
    s.table(["Comparison", "Test", "Result"], STATS, size=8.5,
            title="**Table S3. Statistical summary, related to Figures 3–7**",
            footnote="With three animals per group, p values and intervals are descriptive. Folds with a control "
                     "library without calls placed at its one-event density are lower bounds.")
    # Table S4: SHMtrace against the two published Sanger analyses, metric by metric (sanger_compare.py)
    sys.path.insert(0, str(CRM))
    import sanger_compare
    s.page_break()
    s.table(sanger_compare.HEADER, sanger_compare.table_rows(), size=7.5,
            title="**Table S4. SHMtrace and two published Sanger analyses, metric by metric, related to Figures S8 "
                  "and S9**",
            footnote="Each cell gives wild type \u2192 knockout, the arrow the direction of a significant change, and "
                     "the stars a two-tailed test across mice (Welch on log_{10} for frequency, Student t on "
                     "per-mouse percentages otherwise): \\*\\*\\*, p < 0.001; \\*\\*, p < 0.01; \\*, p < 0.05; ns, not "
                     "significant. Frequencies count unique mutations over the bases sequenced, as both papers do: "
                     "SHMtrace counts the bases each clone read at base quality \u2265 20 and sets paralog clones "
                     "aside, whereas both papers counted every filed clone at one fixed length (577 bp for Maea, "
                     "556 bp for Mkln1, recovered from their tables). \"Counted as published\" applies that "
                     "convention to SHMtrace's mutations. The evidence-based clone set (Mkln1 only) moves or sets "
                     "aside the clones whose filing the sequences contradict, by rules fixed in advance (STAR "
                     "Methods). "
                     "Published values are the Maea paper's Source Data (frequency and G:C shares) or computed from the papers' "
                     "per-mouse tables (the other rows, and every Mkln1 row), the *Mkln1* tables read as "
                     "printed. Maea, *Maea*^{fl/fl}*Mb1*^{cre} (Barbulescu et al., 2024); Mkln1, *Mkln1*^{\u2212/\u2212} "
                     "(Barbulescu et al., 2025). N.R., not reported.")
    s.save(out / "Document_S1.docx")
    print("Document S1 written")


STATS = [
    ("J_{H}4, Wt vs *Aicda*^{−/−} (nanopore, unique-variant density)", "Welch t on log_{10}; n = 3 vs 3",
     "≥ {ont_fold}-fold; p = {ont_p}"),
    ("J_{H}4, Wt vs *Aicda*^{−/−} (Sanger, unique-variant density)", "Welch t on log_{10}; n = 3 vs 3",
     "{sanger_fold}-fold; p = {sanger_p}"),
    ("N114, pMX-*Aicda* vs empty vector", "geometric means; n = 3 vs 3", "≥ {sep_n114}-fold"),
    ("J_{H}4, full length (positions 90–1,230)", "geometric means", "{full_wt} vs {full_ko} × 10^{−5}; ~{full_fold}-fold"),
    ("Spectrum, Sanger vs nanopore (12 types, 3 Wt mice)", "Pearson r [95% CI]; slope",
     "r = {r_spec} [{r_spec_ci}]; slope {slope}"),
    ("Spectrum, per type", "repeated-measures interaction; paired t with Šídák",
     "interaction p = {interaction_p}; smallest adjusted p = {sidak_min}"),
    ("Position along the window, Sanger vs nanopore", "Pearson r, 25-bp bins", "r = {r_pos}"),
    ("Targeting index, Wt calls (J_{H}4)", "permutation, 2,000 position sets", "{s5f_wt}; p ≤ {s5f_wt_p}"),
    ("Targeting index, *Aicda*^{−/−} mismatch positions", "permutation", "{s5f_ko}; p = {s5f_ko_p}"),
    ("Group tests without a true difference (simulated)", "share with p < 0.05, 300 experiments",
     "pooled reads 0.91–0.94; negative binomial 0.050–0.057; Welch 0.037–0.043"),
    ("Leave-one-control-out margin", "held-out control tests with p < 0.01 at the samples' depth",
     "J_{H}4 {loco_jh4_face}% at 1×, {loco_jh4_final}% at {loco_jh4_margin}×; N114 {loco_n114_face}% at 1×, "
     "{loco_n114_final}% at {loco_n114_margin}×"),
    ("Spike-in FDR, control-based vs Benjamini–Hochberg", "observed false/(true + false) calls",
     "J_{H}4 {sp_jh4_fdr} vs {sp_jh4_bh_fdr}; N114 {sp_n114_fdr} vs {sp_n114_bh_fdr} (other C3 library: "
     "{c31_fdr} vs {c31_bh_fdr})"),
]

SUPP_FIGS = [
    ("Figure S1. Detailed workflow and read accounting, related to Figure 1",
     CAL / "manuscript_corrected/figure1_scicanvas/out/Figure_workflow_detailed_preview.png",
     "(A) Experimental design for the J_{H}4 and N114 assays. (B and C) Read preparation and variant calling as run by "
     "SHMtrace from the delivered reads. (D) The Sanger route. (E) The fate of every delivered read per library; line, "
     "the 500 reads analyzed."),
    ("Figure S2. Reported base quality against the error measured in AID-null controls, related to Figure 2",
     FIG / "FigureS2.png",
     "Mismatch rate among the AID-null control bases of each reported quality (bins of two quality values with "
     "at least 2,000 bases; germline sites excluded) for the four nanopore loci, sequenced in three runs; dashed "
     "line, the rate a quality score promises, 10^{−Q/10}."),
    ("Figure S3. Single-read estimate, power and UMI consensus on simulated libraries, related to Figure 3",
     FIG / "FigureS3.png",
     "(A) Change in the aggregate single-read estimate against one-read mutations planted into held-out J_{H}4 "
     "control libraries (three controls; three replicates per amount); dashed line, identity. (B) Power at α = 0.05 "
     "to detect a two-fold change with the animal-level test, by animals per group, from the J_{H}4 wild-type pilot "
     "(1,000 simulated experiments per point). (C) False calls per library with and without UMI consensus in "
     "simulated UMI libraries built from real control reads (three held-out controls)."),
    ("Figure S4. N114 V-region results and per-library spectra, related to Figure 4",
     [MC / "Figure4_hybridoma_combined.png", MC / "FigureS6_individual_JH4_ONT_matrices.png",
      MC / "FigureS7_individual_N114_matrices.png"],
     "Top, N114 unique-variant density per subclone (n = 3 per vector; empty-vector subclones without calls at one "
     "event) and pooled spectra. Middle and bottom, From → To substitution matrices for each J_{H}4 mouse and N114 "
     "subclone; each cell gives unique variants and their share of the sample's calls."),
    ("Figure S5. Sequencing depth and unique-variant discovery, related to Figure 4",
     [MC / "FigureS3_depth_discovery_density.png", MC / "FigureS4_saturation.png"],
     "Each library was analyzed at 30–1,000 reads (one random draw of usable reads per depth), with controls "
     "calibrated at the same depth. Top, unique variants and density against depth (geometric mean ± geometric SD, "
     "n = 3 per group). Bottom, calls at each depth relative to each AID-active sample's maximum (mean ± SD); at 500 "
     "reads, {depth_jh4_500}% (J_{H}4) and {depth_n114_500}% (N114)."),
    ("Figure S6. The background along the amplicon and in each read, related to Figure 4",
     FIG / "FigureS6.png",
     "(A and B) Mismatches at base quality ≥ 20 at every position of the analysis window, germline sites "
     "excluded, in each AID-null library (red) and averaged over the AID-active libraries (blue); J_{H}4, 500 "
     "reads per library; Sμ (PRJNA1481264), 2,000. r, Spearman correlation between pairs of controls. (C and D) "
     "Every mismatch at quality ≥ 20 (germline sites excluded) in each of the 500 reads of a wild-type (C) and an "
     "*Aicda*^{−/−} (D) J_{H}4 library; colored, at sites SHMtrace called; gray, other mismatches. Reads are sorted "
     "by clone (dotted lines between the largest clones) and by mutation load."),
    ("Figure S7. Sanger re-analysis, related to Figure 5",
     [MC / "FigureS5_individual_Sanger_matrices.png", CAL / "figures/sanger_vs_nanopore/figS_sanger_window_original.png"],
     "Top, substitution matrices for each mouse's Sanger clones (re-analyzed from 447 chromatograms; base quality "
     "≥ 20). Bottom, the share of wild-type clones with a high-quality base at each J_{H}4 position and calls per 10^{3} "
     "bases read; shading, the 529-bp window used for comparison."),
    ("Figure S8. A published Sanger study re-analyzed automatically, related to Figure 5",
     [FIG / "FigureS8_sanger_cohort2.png"],
     "J_{H}4 clones of germinal-center B cells from four wild-type (blue), four *Maea*^{fl/fl}*Mb1*^{cre} (green) and "
     "*Aicda*^{−/−} (red) mice[@barbulescu2024], each clone assigned to the mouse it was filed under in the original "
     "analysis. Filled symbols and blue lines, SHMtrace from the chromatograms; open squares and brown lines, the "
     "published values, replotted from the paper's per-mouse tables (Supplementary Fig. 9b) and clone distributions "
     "(Fig. 2g). (A) Mutations per 10^{4} bp in each mouse; lines, means. (B) The same per mouse, SHMtrace against "
     "published; dashed line, identity. (C) Clones by mutations per clone, pooled by genotype, binned as published. "
     "(D) G:C transitions and transversions as a share of each mouse's mutations. (E) Pooled 12-type substitution "
     "spectra, SHMtrace against published. (F) The published comparisons and SHMtrace's, with the one-tailed "
     "t-test used in the original study; \"background explains\" is the share of each group's calls that the "
     "*Aicda*^{\u2212/\u2212} clones' own error rate accounts for (STAR Methods), and the analysis time is for "
     "the whole project from one sample sheet."),
    ("Figure S9. A second published Sanger study re-analyzed automatically, related to Figure 5",
     [FIG / "FigureS9_sanger_mkln.png"],
     "J_{H}4 clones of germinal-center B cells from three wild-type (blue), three *Mkln1*^{\u2212/\u2212} "
     "(purple) and three *Aicda*^{\u2212/\u2212} (red) mice[@barbulescu2025], each clone assigned to the mouse "
     "folder it was filed under in the original analysis. Filled symbols, SHMtrace from the chromatograms; open "
     "squares, the published values, replotted from the paper's per-mouse tables (Supplementary Figure S2A). "
     "(A) Unique mutations per 10^{4} bases sequenced, per mouse; lines, means. (B) G:C transitions and "
     "transversions as a share of each mouse's mutations. (C) Pooled 12-type substitution spectra, SHMtrace "
     "against published, with r for each genotype. Of the three published *Mkln1*^{\u2212/\u2212} tables, only "
     "mouse #1's spectrum matches the clones of a folder (K2); tables #2 and #3 match no folder or set of clones in "
     "the original project, which limits the knockout's agreement. (D) The comparisons, published and SHMtrace's; the published P values are computed "
     "from the per-mouse tables as printed, and the values the paper's summary bars report (Supplementary "
     "Figure S2C) are listed separately. Paralog clones are J_{H}4-paralog amplicons, which SHMtrace analyzes "
     "apart from J_{H}4; \"background explains\" is as in Figure S8."),
    ("Figure S10. A published Sanger study from another laboratory re-analyzed automatically, related to Figure 5",
     [FIG / "FigureS10_swhel.png"],
     "The 1,795 Sanger clones of Thientosapol et al. (2018), who hypermutated a SW_{HEL} V_{H}10 knock-in V region "
     "in germinal-center B cells and reported that inactivating *Samhd1* raises transitions at C:G and lowers "
     "transversions at A:T and at C:G. SHMtrace was given only the deposited GenBank records: the germline is the "
     "iterated majority consensus of the clones themselves, each clone is assigned to the animal its own record "
     "names, and the arms are the genotypes its records state. (A) Mutations per sequence in each arm, SHMtrace "
     "against the paper's Table 1. (B) The 12 substitution types in each of the five arms, as a share of that "
     "arm's mutations. (C) The four classes the paper tests, one point per animal for SHMtrace and the published "
     "value as an open square; stars, two-tailed t-test across animals. (D) Each published claim against what "
     "SHMtrace finds. Unlike Figures S8 and S9 this study deposits sequences rather than chromatograms, so the "
     "chromatogram reading and the bookkeeping checks are not exercised; it has no *Aicda*^{\u2212/\u2212} arm, "
     "and Sanger error (about 2 in 10^{5} bases, measured in the AID-null clones of those two studies) is some "
     "500-fold below the mutation frequency here."),
    ("Figure S11. Full-length spectra and mutation trees in each mouse, related to Figure 6",
     [MC / "FigureS_spectrum_529_vs_fulllength.png", CAL / "figures/mutation_trees/mutation_trees_by_mouse.png"],
     "Top, pooled wild-type spectra over the 529-bp window and the full 1,141-bp interval (Pearson r = {s1_r}). "
     "Bottom, the most branched clone of each wild-type mouse drawn as a mutation tree."),
    ("Figure S12. What Sanger clones and nanopore reads each measure, related to Figures 5 and 6",
     [FIG / "FigureS12_platforms.png"],
     "Top, one unit of data from each platform on one base-pair scale. A Sanger clone is one molecule read once "
     "over the 529-bp window that the clones read cleanly, so each clone's mutations are read directly and "
     "sequencing errors are rare, though not absent. A nanopore read spans the whole 1,230-bp amplicon; the eight "
     "reads drawn are real reads of wild-type library Wt3, from four clones, with their called mutations in place, "
     "and errors are drawn at the rate measured in the *Aicda*^{\u2212/\u2212} controls over the window "
     "(0.14% of bases at any quality, 85% of them below Q20, open circles; 0.02% at Q \u2265 20, filled), at "
     "random positions. A mutation is called where more reads carry it than the controls "
     "explain; single-read mutations that cannot be told from error one at a time "
     "are estimated in aggregate. Bottom, what the Sanger clones and the nanopore reads of this study each gave, "
     "and what both gave in the six mice of Figure 5, with the agreement measured there. The split describes the "
     "two workflows evaluated here, not the platforms in general: Sanger clones of an amplicon that spans the "
     "V(D)J junction can also be grouped into lineages, UMI consensus also gives one value per molecule, and how "
     "much of a locus a long read spans depends on the assay. The nanopore fold counts controls without calls at "
     "one event. Sanger error rate: mismatches at base quality \u2265 20 in the AID-null clones of the two "
     "published studies of Figures S8 and S9."),
    ("Figure S13. An overview of SHMtrace, related to STAR Methods",
     [FIG / "FigureS13_overview.png"],
     "From libraries to group tests. Reads or chromatograms are assigned to animals through a sample sheet, "
     "gated and aligned competitively against the target, its paralog and any contaminant; the AID-null controls "
     "measure the background, each candidate in a read is tested against it, and the controls, called against one "
     "another as decoys, set the threshold. Sanger clones are counted instead, with the background reported beside "
     "them (Figure 1C). The read-length limit is set per locus (STAR Methods), and changes carried by half the "
     "reads or more are set aside as likely germline, a heuristic that can also remove a real mutation shared by "
     "a dominant clone. Results are summarised per animal and compared between groups with one value per "
     "animal. Bookkeeping checks are reported beside the results and never applied. Figure S14 gives the rule and "
     "equation behind each step."),
    ("Figure S14. How SHMtrace turns reads into calls, related to STAR Methods",
     [FIG / "FigureS14_engine.png"],
     "Each of the six stages of the analysis, with the rule it applies and that rule working on the J_{H}4 "
     "libraries. (A) The fate of every delivered read, per library: gated by length or expected error, sent "
     "elsewhere by competitive alignment, or usable (Figure S1 separates the gates). (B) Each read's identity "
     "to the J_{H}4 target against its identity to the "
     "paralogue; a read is analyzed at whichever it matches better. (C) Reported base quality against the error "
     "measured in the *Aicda*^{\u2212/\u2212} controls, with the rate the quality score promises (dashed). "
     "(D) Two candidates at the same depth in the same library, one called and one not, with the reads carrying "
     "the change in the wild-type and in the control libraries. (E) Candidates kept against the p-value "
     "threshold for a wild-type library and for the held-out controls used as decoys, and the threshold their "
     "crossing sets. (F) Unique-variant density per animal, the value the group tests use. Plotted from "
     "`intuition_jh4.json`, the per-library read accounting and the bundled example library; the equations are "
     "the ones the code applies."),
]

if __name__ == "__main__":
    main(sys.argv[1])
