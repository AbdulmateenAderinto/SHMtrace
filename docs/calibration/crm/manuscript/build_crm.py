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
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import crm_methods as M  # noqa: E402
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
    n["pim1_deep"] = (f"but with all usable reads ({min(reads):,}–{max(reads):,} per library) SHMtrace called {a[0]}–{a[1]} "
                      f"mutations in every AID-competent library and {'none' if b[1] == 0 else f'{b[0]}–{b[1]}'} in the "
                      f"controls, whereas the Geneious rule called {ga[0]}–{ga[1]} and {gb[0]}–{gb[1]}. A weak AID target "
                      "therefore needs depth, and the control calibration is what keeps the extra depth from adding "
                      "background (Table S2).")
    # Illumina J_H4 (Ldha)
    p = json.load(open(OUT / "public/ldha/project.json"))
    f = {g: [s["mutation_frequency_x1e3"] for s in p["samples"] if s["group"] == g] for g in ("WT", "LDHA_KO", "AID_KO")}
    g = {x["group"]: x for x in p["groups"]}
    n["ldha_fold"] = f"{np.mean(f['WT']) / np.mean(f['AID_KO']):.0f}"
    n["ldha_p"] = f"{g['WT']['mutation_frequency_x1e3_p_welch_log_vs_AID_KO']:.3f}"
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
    n["hmces_del_p"] = f"{stats.ttest_ind(np.log10(np.array(dels['HMCES_KO']) + 1), np.log10(np.array(dels['WT']) + 1), equal_var=False).pvalue:.1e}"
    n["hmces_ko3_gc"] = f"{s['AIDKO_3']['pct_GC']:.0f}"
    return n


def fill(text, n):
    for k, v in n.items():
        text = text.replace("{" + k + "}", str(v))
    if "{" in text.replace("^{", "").replace("_{", ""):
        raise RuntimeError(f"unfilled placeholder in: {text[:80]}")
    return text


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
    n = public_numbers()
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
            title="**Table S1. Read preparation and calls per library, related to Figure 3**",
            footnote="Counts of reads removed at each step of SHMtrace's read preparation, from the delivered FASTQ "
                     "files; 500 usable reads were drawn at random (seed 0) for analysis. WT, wild type; KO, "
                     "*Aicda*^{−/−}; C1, D3, D5, pMX-*Aicda* subclones; C3, C4, D1, empty-vector subclones.")
    # Table S3: statistics
    s.page_break()
    s.table(["Comparison", "Test", "Result"], STATS, size=8.5,
            title="**Table S3. Statistical summary, related to Figures 2–6**",
            footnote="With three animals per group, p values and intervals are descriptive. Folds with a control "
                     "library without calls placed at its one-event density are lower bounds.")
    s.save(out / "Document_S1.docx")
    print("Document S1 written")


STATS = [
    ("J_{H}4, Wt vs *Aicda*^{−/−} (nanopore, distinct-call density)", "Welch t on log_{10}; n = 3 vs 3",
     "≥ 137-fold; p = 0.0025"),
    ("J_{H}4, Wt vs *Aicda*^{−/−} (Sanger, distinct-call density)", "Welch t on log_{10}; n = 3 vs 3", "38-fold; p = 0.0009"),
    ("N114, pMX-*Aicda* vs empty vector", "Welch t on log_{10}; n = 3 vs 3", "≥ 92-fold; p = 0.0048"),
    ("J_{H}4, full length (positions 90–1,230)", "geometric means", "39 vs 0.51 × 10^{−5}; ~78-fold"),
    ("Spectrum, Sanger vs nanopore (12 types, 3 Wt mice)", "Pearson r [95% CI]; slope", "r = 0.93 [0.77–0.98]; slope 1.00"),
    ("Spectrum, per type", "repeated-measures interaction; paired t with Šídák", "interaction p = 0.61; smallest adjusted p = 0.68"),
    ("Position along the window, Sanger vs nanopore", "Pearson r, 25-bp bins", "r = 0.92"),
    ("Targeting index, Wt calls (J_{H}4)", "permutation, 2,000 position sets", "1.22–1.36; p ≤ 0.015"),
    ("Targeting index, *Aicda*^{−/−} mismatch positions", "permutation", "1.00; p = 0.50"),
    ("Group tests without a true difference (simulated)", "share with p < 0.05, 300 experiments",
     "pooled reads 0.91–0.94; negative binomial 0.050–0.057; Welch 0.037–0.043"),
    ("Leave-one-control-out margin", "held-out control tests with p < 0.01", "J_{H}4 4.0% at 1×, 0.9% at 5×; N114 5.6% at 1×, 0.3% at 5×"),
]

SUPP_FIGS = [
    ("Figure S1. Detailed workflow and read accounting, related to Figure 1",
     CAL / "manuscript_corrected/figure1_scicanvas/out/Figure_workflow_detailed_preview.png",
     "(A) Experimental design for the J_{H}4 and N114 assays. (B and C) Read preparation and variant calling as run by "
     "SHMtrace from the delivered reads. (D) The Sanger route. (E) The fate of every delivered read per library; line, "
     "the 500 reads analyzed."),
    ("Figure S2. Sequencing depth and distinct-call discovery, related to Figure 3",
     [MC / "FigureS3_depth_discovery_density.png", MC / "FigureS4_saturation.png"],
     "Each library was analyzed at 30–1,000 reads (one random draw of usable reads per depth), with controls "
     "calibrated at the same depth. Top, distinct calls and density against depth (geometric mean ± geometric SD, "
     "n = 3 per group). Bottom, calls at each depth relative to each AID-active sample's maximum (mean ± SD); at 500 "
     "reads, 94% (J_{H}4) and 91% (N114)."),
    ("Figure S3. N114 V-region results and per-library spectra, related to Figure 3",
     [MC / "Figure4_hybridoma_combined.png", MC / "FigureS6_individual_JH4_ONT_matrices.png",
      MC / "FigureS7_individual_N114_matrices.png"],
     "Top, N114 distinct-call density per subclone (n = 3 per vector; empty-vector subclones without calls at one "
     "event) and pooled spectra. Middle and bottom, From → To substitution matrices for each J_{H}4 mouse and N114 "
     "subclone; each cell gives distinct calls and their share of the sample's calls."),
    ("Figure S4. Sanger re-analysis, related to Figure 4",
     [MC / "FigureS5_individual_Sanger_matrices.png", CAL / "figures/sanger_vs_nanopore/figS_sanger_window_original.png"],
     "Top, substitution matrices for each mouse's Sanger clones (re-analyzed from 447 chromatograms; base quality "
     "≥ 20). Bottom, the share of wild-type clones with a high-quality base at each J_{H}4 position and calls per 10^{3} "
     "bases read; shading, the 529-bp window used for comparison."),
    ("Figure S5. Single-read estimate, power and UMI consensus on simulated libraries, related to Figure 2",
     FIG / "FigureS5.png",
     "(A) Change in the aggregate single-read estimate against one-read mutations planted into held-out J_{H}4 "
     "control libraries (three controls; three replicates per amount); dashed line, identity. (B) Power at α = 0.05 "
     "to detect a two-fold change with the animal-level test, by animals per group, from the J_{H}4 wild-type pilot "
     "(1,000 simulated experiments per point). (C) False calls per library with and without UMI consensus in "
     "simulated UMI libraries built from real control reads (three held-out controls)."),
    ("Figure S6. Full-length spectra and mutation trees in each mouse, related to Figure 5",
     [MC / "FigureS_spectrum_529_vs_fulllength.png", CAL / "figures/mutation_trees/mutation_trees_by_mouse.png"],
     "Top, pooled wild-type spectra over the 529-bp window and the full 1,141-bp interval (Pearson r = 0.986). "
     "Bottom, the most branched clone of each wild-type mouse drawn as a mutation tree."),
]

if __name__ == "__main__":
    main(sys.argv[1])
