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
    ("Figure S2. Sequencing depth and unique-variant discovery, related to Figure 3",
     [MC / "FigureS3_depth_discovery_density.png", MC / "FigureS4_saturation.png"],
     "Each library was analyzed at 30–1,000 reads (one random draw of usable reads per depth), with controls "
     "calibrated at the same depth. Top, unique variants and density against depth (geometric mean ± geometric SD, "
     "n = 3 per group). Bottom, calls at each depth relative to each AID-active sample's maximum (mean ± SD); at 500 "
     "reads, {depth_jh4_500}% (J_{H}4) and {depth_n114_500}% (N114)."),
    ("Figure S3. N114 V-region results and per-library spectra, related to Figure 3",
     [MC / "Figure4_hybridoma_combined.png", MC / "FigureS6_individual_JH4_ONT_matrices.png",
      MC / "FigureS7_individual_N114_matrices.png"],
     "Top, N114 unique-variant density per subclone (n = 3 per vector; empty-vector subclones without calls at one "
     "event) and pooled spectra. Middle and bottom, From → To substitution matrices for each J_{H}4 mouse and N114 "
     "subclone; each cell gives unique variants and their share of the sample's calls."),
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
     "Top, pooled wild-type spectra over the 529-bp window and the full 1,141-bp interval (Pearson r = {s1_r}). "
     "Bottom, the most branched clone of each wild-type mouse drawn as a mutation tree."),
]

if __name__ == "__main__":
    main(sys.argv[1])
