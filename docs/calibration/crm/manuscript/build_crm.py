"""Build the Cell Reports Methods article and its companion files.

    python docs/calibration/crm/manuscript/build_crm.py OUTDIR

Writes OUTDIR/SHMtrace_CellReportsMethods.docx (manuscript: title page, summary, keywords, introduction,
results, discussion, limitations, resource availability, declarations, supplemental index, figure legends
with figures, STAR Methods with the key resources table, references), OUTDIR/Highlights_eTOC.docx
(highlights, eTOC blurb, motivation, as the journal asks: a separate file) and OUTDIR/Document_S1.docx.
Numbers for the public Illumina sets are read from out/public/*/project.json; the single-read
recovery from out/single_read_spike.json. Journal limits are checked and the build stops if any is
exceeded.
"""
import collections
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
    # the effect and its interval: failing to reach significance is not evidence of no difference
    lw, lk = np.log10([x for x in f["WT"] if x > 0]), np.log10([x for x in f["LDHA_KO"] if x > 0])
    se = (lw.var(ddof=1) / len(lw) + lk.var(ddof=1) / len(lk)) ** 0.5
    dof = ((lw.var(ddof=1) / len(lw) + lk.var(ddof=1) / len(lk)) ** 2
           / ((lw.var(ddof=1) / len(lw)) ** 2 / (len(lw) - 1)
              + (lk.var(ddof=1) / len(lk)) ** 2 / (len(lk) - 1)))
    d, tcrit = lk.mean() - lw.mean(), stats.t.ppf(0.975, dof)
    n["ldha_fold_ko"] = f"{10 ** d:.2f}"
    n["ldha_ci_ko"] = f"{10 ** (d - tcrit * se):.2f}–{10 ** (d + tcrit * se):.2f}"
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


NUMWORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight",
           9: "nine", 10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen",
           15: "fifteen", 16: "sixteen", 17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}


def dataset_counts():
    """The loci and laboratory counts the Summary, Discussion and Table S6 footnote all quote, taken
    from the table itself through M.TABLE_S6_LOCI rather than written out in three places.

    Every row of the table must be classified and every classification must name a row, so adding a
    dataset fails the build instead of leaving the counts stale. A locus is counted once however many
    rows share it, and a row whose locus is not an AID target counts towards neither total."""
    rows = {r[0] for r in M.TABLE_S6}
    declared = set(M.TABLE_S6_LOCI)
    assert rows == declared, (f"Table S6 and its locus classification disagree; rows not classified: "
                              f"{sorted(rows - declared)}; classified but absent: {sorted(declared - rows)}")
    cls, labs = {}, set()
    for dataset, (lab, loci) in M.TABLE_S6_LOCI.items():
        for locus, c in loci:
            if c == "none":
                continue
            assert cls.setdefault(locus, c) == c, f"{locus} is classified two ways"
            if lab:
                labs.add(lab)
    counts = collections.Counter(cls.values())
    n_loci = len(cls)
    assert counts["ig"] + counts["off"] + counts["construct"] == n_loci
    # the prose spells these out; past twenty it would read badly as a word, so use the digits
    w = lambda i: NUMWORD.get(i, str(i))   # noqa: E731
    return {"ds_loci": w(n_loci), "ds_labs": w(len(labs)), "ds_shm": w(counts["ig"]),
            "ds_offtarget": w(counts["off"]), "ds_engineered": w(counts["construct"])}


def check_figure_numbers(n):
    """The graphical abstract is drawn in SciCanvas with its numbers typed in: they must match the analysis."""
    ga = (CRM / "scicanvas/crm_ga.js").read_text()
    for text in (f"{n['ctl_rec']} calls per control set", f"r = {n['r_spec']}"):
        if text not in ga:
            raise RuntimeError(f"graphical abstract is stale: expected '{text}' in scicanvas/crm_ga.js")


def check_paired_settings(n):
    """STAR Methods reports the shared-artifact and case-only resamples once, "in both settings".
    That wording is only honest while exclusion on and off agree, so the build checks it rather
    than trusting the prose: if they ever diverge, the sentence has to name both numbers again."""
    for cond in ("one_noisy", "shared", "case_only"):
        for field in ("any", "n"):
            on, off = n[f"fdv_null_{cond}_on_{field}"], n[f"fdv_null_{cond}_off_{field}"]
            assert on == off, (f"STAR Methods says the {cond} resamples gave the same result with "
                               f"exclusion on and off, but {field} is {on} vs {off}; the text must "
                               "report both settings separately")
    # the masking passage quotes one figure for the 60% and 80% dominant-clone libraries
    assert n["mask60_lost"] == n["mask80_lost"], (
        f"STAR Methods quotes one loss figure for 60% and 80% dominance, but they are "
        f"{n['mask60_lost']} and {n['mask80_lost']}; the text must give both")


def journal_word_count():
    """What Cell Reports Methods counts against its 7,000-word limit.

    "Research articles should be under 7,000 words, not including the reference list or the STAR
    Methods" (cell.com/cell-reports-methods/information-for-authors/article-types, read 2026-10-18).
    Those two are the only exclusions named, so everything else in the main-text manuscript counts:
    the summary, keywords, resource availability, the declarations, the supplemental-information
    index and, by far the largest of them, the figure titles and legends. Counting only
    introduction to limitations, as this build did until now, understated the total by about a
    quarter and let the manuscript sit well over the limit while the build reported it as compliant.

    Document S1 is a separate supplemental file, so supplemental figure legends are not counted."""
    W = lambda *t: word_count([("p", x) for x in t if x])   # noqa: E731
    parts = [
        ("title, authors, affiliations", W(T.TITLE, T.AUTHORS, *T.AFFIL)),
        ("summary", W(T.SUMMARY)),
        ("keywords", W(T.KEYWORDS)),
        ("introduction", W(*T.INTRO)),
        ("results", W(*[x for _k, x in T.RESULTS])),
        ("discussion", W(*T.DISCUSSION)),
        ("limitations of the study", W(*T.LIMITATIONS)),
        ("resource availability", W(*[h for h, _x in T.RESOURCE], *[x for _h, x in T.RESOURCE])),
        ("acknowledgments", W(T.ACK)),
        ("author contributions", W(T.CONTRIB)),
        ("declaration of interests", W(T.INTERESTS)),
        ("declaration of generative AI", W(T.AI_DECL)),
        ("supplemental information index", W(*M.SUPPLEMENT)),
        ("figure titles and legends", W(*[x for x, _l in M.LEGENDS], *[l for _x, l in M.LEGENDS])),
    ]
    return sum(v for _k, v in parts), parts


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
    n = {**crm_numbers.compute(), **public_numbers(), **dataset_counts()}
    fill_modules(n)
    check_figure_numbers(n)
    check_paired_settings(n)
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
    for p in T.LIMITATIONS:      # the author splits the limitations across four paragraphs
        m.para(p)
        main_blocks.append(("p", p))
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
    jw, parts = journal_word_count()
    print(f"{path}\nsummary {sw} words; introduction to limitations {mw} words; references {len(m.order)}")
    print(f"journal count (all but references and STAR Methods) {jw} words of 7,000")
    if jw >= 7000:
        for k, v in sorted(parts, key=lambda kv: -kv[1]):
            print(f"    {k:34s} {v:>5}")

    # legends and figures live in different files and have drifted apart before, so check them against
    # each other here rather than trusting a later read-through (legend_check.py)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import legend_check
    assert legend_check.main() == 0, "legends disagree with the figures they describe"
    # details only the authors can supply are left as "[description]": name them rather than let a placeholder
    # reach a reviewer unnoticed (lab_details.py holds them)
    import lab_details
    still = lab_details.missing()
    if still:
        print(f"submission blockers: {len(still)} laboratory detail(s) still to supply")
        for d in still:
            print("  -", d)

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
        imgs = img if isinstance(img, list) else [img]
        s.figure_block(imgs, title, legend)
    # Table S1: what SHMtrace did with every delivered read, and the calls (J_H4 and N114)
    s.page_break()
    rows = []
    for run, label in (("jh4", "J_{H}4"), ("n114", "N114")):
        acct = json.load(open(CAL / f"figures/{run}_shmtrace/table_s1.json"))
        proj = json.load(open(CAL / f"figures/{run}_shmtrace/project.json"))
        by = {x["sample"]: x for x in proj["samples"]}
        num = lambda v: "–" if v in (None, "") else f"{v:g}"                        # noqa: E731
        for lib, a in acct.items():
            m = by.get(lib, {})
            rows.append([label, lib, f"{a['delivered']:,}", f"{a['too_long']:,}", f"{a['low_quality']:,}",
                         f"{a['not_target']:,}", f"{a['chimeric'] + a['short_alignment']:,}",
                         f"{a['paralogue_alleles']:,}", f"{a['usable']:,}", a["analysed"], a["candidates_tested"],
                         a["calls"], num(m.get("unique_density_x1e5")), num(m.get("mutation_frequency_x1e3")),
                         num(m.get("clone_weighted_frequency_x1e3"))])
    s.table(["Assay", "Library", "Delivered", "Over-length", "Expected error > 0.01", "Unmapped or paralog",
             "Chimeric or short alignment", "Paralog alleles", "Usable", "Analyzed", "Candidates tested",
             "Calls (FDR ≤ 0.05)", "Unique-variant density (×10^{−5})", "Read-weighted frequency (×10^{−3})",
             "Clone-weighted frequency (×10^{−3})"], rows, size=7.5,
            title="**Table S1. Read preparation and calls per library, related to Figure 4**",
            footnote="Per-library read removal at each preparation step, starting from delivered FASTQ files. Five "
                     "hundred usable reads were randomly selected (seed 0). All three mutation endpoints are reported "
                     "(STAR Methods): unique-variant density counts each distinct call once per library over "
                     "analyzed base observations; read-weighted frequency counts its support across reads; "
                     "clone-weighted frequency gives each junction-defined clone equal weight. A dash indicates that "
                     "no clone had an assignable junction. WT, wild type; KO, *Aicda*^{−/−}; C1, D3 and D5, "
                     "pMX-*Aicda* subclones; C3, C4 and D1, empty-vector subclones.")
    # Table S3: statistics
    s.page_break()
    s.table(["Comparison", "Test", "Result"], STATS, size=8.5,
            title="**Table S3. Statistical summary, related to Figures 3–7**",
            footnote="With three animals per group, p values and intervals should be interpreted descriptively. Fold "
                     "separation depends on the one-event convention used for zero-call controls and does not "
                     "represent a statistical lower confidence bound.")
    # Table S4: SHMtrace against the two published Sanger analyses, metric by metric (sanger_compare.py)
    sys.path.insert(0, str(CRM))
    import sanger_compare
    s.page_break()
    s.table(sanger_compare.HEADER, sanger_compare.table_rows(), size=7.5,
            title="**Table S4. SHMtrace and two published Sanger analyses, metric by metric, related to Figures S8 "
                  "and S9**",
            footnote="Cells report Wt → knockout values; arrows indicate the direction of significant changes. Stars "
                     "denote two-sided tests across mice: Welch tests on log_{10} frequency and Student t tests on "
                     "per-mouse percentages; \\*\\*\\*, p < 0.001; \\*\\*, p < 0.01; \\*, p < 0.05; ns, not "
                     "significant. Frequencies count unique mutations over clone number × fixed sequence length. "
                     "SHMtrace uses analyzed clones and the shared 529-bp window, excluding paralog amplicons. The "
                     "publications count all filed clones across their full consensus lengths, recovered from the "
                     "tables as 577 bp for Maea and 556 bp for *Mkln1*. “Counted as published” applies those "
                     "denominators to SHMtrace mutation counts; the paired rows therefore differ only in denominator. "
                     "The evidence-based *Mkln1* clone set reassigns or sets aside clones with sequence evidence "
                     "contradicting their filing, using prespecified rules (STAR Methods). Published Maea frequency "
                     "and G:C values are from Source Data; other Maea values and all *Mkln1* values were computed from "
                     "the printed per-mouse tables. Maea, *Maea*^{fl/fl}*Mb1*^{cre} (Barbulescu et al., 2024); "
                     "*Mkln1*, *Mkln1*^{−/−} (Barbulescu et al., 2025). N.R., not reported.")
    # Table S5: the same five questions of all three re-analyses (agreement_table.py)
    import agreement_table
    s.page_break()
    s.table(agreement_table.HEADER, agreement_table.table_rows(), widths=[20, 30, 22, 30, 28, 22, 30], size=6.5,
            title="**Table S5. What \"reproduced\" means for each re-analyzed study, related to Figures S8\u2013S10**",
            footnote=agreement_table.NOTE)
    # Table S6: the dataset inventory. It was Table 1 while the main text had room for a table; with
    # seven figures, figures and tables together reach the journal's limit, so it travels here.
    s.page_break()
    s.table(M.TABLE_S6_HEADER, M.TABLE_S6, title=f"**{M.TABLE_S6_TITLE}**", footnote=M.TABLE_S6_NOTE,
            widths=[28, 20, 28, 32, 28, 20, 18], size=7.5)
    s.save(out / "Document_S1.docx")
    print("Document S1 written")
    # last, so that both documents are on disk to work from even while the manuscript is over length
    assert jw < 7000, (f"manuscript is {jw} words by the journal's count, {jw - 6999} over; the limit "
                       "excludes only the reference list and the STAR Methods")


STATS = [
    ("J_{H}4, Wt vs *Aicda*^{−/−} (nanopore, unique-variant density)", "Welch t on log_{10}; n = 3 vs 3",
     "{ont_fold}-fold; p = {ont_p}"),
    ("J_{H}4, Wt vs *Aicda*^{−/−} (Sanger, unique-variant density)", "Welch t on log_{10}; n = 3 vs 3",
     "{sanger_fold}-fold; p = {sanger_p}"),
    ("N114, pMX-*Aicda* vs empty vector", "geometric means; n = 3 vs 3", "{sep_n114}-fold"),
    ("J_{H}4, full length (positions 90–1,230)", "geometric means", "{full_wt} vs {full_ko} × 10^{−5}; ~{full_fold}-fold"),
    # descriptive, not a population interval: 12 compositional categories pooled from 3 mice, so a Fisher
    # interval here describes the types and not the animals. Results drops it for the same reason.
    ("Spectrum, Sanger vs nanopore (12 types, 3 Wt mice)", "Pearson r (descriptive); slope",
     "r = {r_spec}; slope {slope}"),
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
    # minor 14: the sensitive mode's claims, with the default beside them so the trade is visible
    ("Sensitive mode (FDR 0.20) vs default (0.05), J_{H}4 spike-ins",
     "single-read recall; overall recall; observed FDR; false calls per library",
     "{sens_jh4_1}% vs {sp_jh4_1}%; {sens_jh4_all}% vs {sp_jh4_recall}%; {sens_jh4_fdr} vs {sp_jh4_fdr}; "
     "{sens_jh4_fp} vs {sp_jh4_fp}"),
    ("Sensitive mode (FDR 0.20) vs default (0.05), N114 spike-ins",
     "single-read recall; overall recall; observed FDR; false calls per library",
     "{sens_n114_1}% vs {sp_n114_1}%; {sens_n114_all}% vs {sp_n114_recall}%; {sens_n114_fdr} vs "
     "{sp_n114_fdr}; {sens_n114_fp} vs {sp_n114_fp}"),
    # the two comparisons that did not reach significance, reported as an effect with an interval:
    # not reaching significance is not evidence of no effect, and the interval says what is still possible
    ("*Ldha*-deleted vs Wt (Illumina J_{H}4, mutation frequency)", "Welch interval on log_{10}; Mann–Whitney p; n = 7 vs 7",
     "{ldha_fold_ko}-fold [95% CI {ldha_ci_ko}]; p = {ldha_p_ko}"),
    ("Repair context, transitions at C:G (NIH-3T3, SHMtrace calls)",
     "one value per experiment; n = 3 per arm", "{mo_ugi_ts}% (Ugi) vs {mo_ung_ts}% (UNG); {mo_arm_p}"),
    ("Spike-in FDR, control-based vs Benjamini–Hochberg", "observed false/(true + false) calls",
     "J_{H}4 {sp_jh4_fdr} vs {sp_jh4_bh_fdr}; N114 {sp_n114_fdr} vs {sp_n114_bh_fdr} (other C3 library: "
     "{c31_fdr} vs {c31_bh_fdr})"),
]

SUPP_FIGS = [
    ("Figure S1. Detailed workflow and read accounting, related to Figure 1",
     FIG / "FigureS1_workflow_preview.png",
     "(A) Experimental design for J_{H}4 and N114. (B, C) Read preparation and variant calling from delivered sequence "
     "files, corresponding to the read-library route in Figure 1C. By default, alleles carried by ≥50% of reads are "
     "excluded as likely germline differences; this heuristic can also exclude acquired mutations (STAR Methods). (D) "
     "Sanger analysis counts substitutions at base quality ≥20 across analyzed clones and the 529-bp window. The "
     "mismatch rate measured in AID-null clones is reported separately without subtraction. Each sequencing clone "
     "represents one PCR-derived molecule sequenced once. (E) Accounting for every delivered read in each library; "
     "line, the 500 reads selected for analysis."),
    ("Figure S2. Reported base qualities and observed control mismatch rates, related to Figure 2",
     FIG / "FigureS2.png",
     "Top, mismatch rates in AID-null controls for the four nanopore loci, grouped into two-quality-score bins "
     "containing ≥2,000 bases after germline exclusion. Shading, 95% Clopper–Pearson intervals; downward triangles, "
     "zero-mismatch bins plotted at the upper confidence bound; dashed line, the nominal error probability, "
     "10^{−Q/10}. Bottom, number of control bases per bin. Sparse observations above Q40 increase uncertainty. Each "
     "locus is calibrated using its own controls. J_{H}4 libraries shared sequencing order DHJGKB, and N114 libraries "
     "shared order 75K3B3. For Sμ and *Pim1* (PRJNA1481264), the deposited records do not establish which libraries "
     "shared a flow cell."),
    ("Figure S3. Single-read estimate, power, UMI consensus and planted indels, related to Figure 3",
     FIG / "FigureS3.png",
     "(A) Change in the aggregate single-read estimate against one-read mutations planted into held-out J_{H}4 control "
     "libraries (three controls; three replicates per amount); dashed line, identity. (B) Power to detect a two-fold "
     "difference between two AID-proficient groups at α = 0.05, two-sided, with the Welch test on log rates, by "
     "animals per group, simulated from the J_{H}4 wild-type pilot (3 mice, 500 reads each; 1,000 experiments per "
     "point; negative-binomial counts with the pilot's clumping, animals log-normal with the pilot's between-animal "
     "coefficient of variation, CV). Dotted line, 80% power; for unique density the estimate reaches it at about 24 "
     "animals per group, so the plotted range does not by itself identify the number required. Shading, the same curve "
     "with the CV at the lower end of its 95% interval from 3 mice (0.44) and at twice the estimate (1.9), a stated "
     "pessimistic case: with 2 degrees of freedom the interval's upper end is effectively unbounded, so the CV, not "
     "the simulation, is what limits this calculation. (C and E) Simulated UMI libraries built from real control "
     "reads, three held-out controls each: false calls per library, and the recall of planted mutations that goes "
     "with them. Consensus removes every false call in three of the four conditions and leaves {umi_fp_poly} in "
     "the fourth, annotated on the panel because a zero-height bar cannot be told from a missing one, while recall "
     "rises from {umi_rec_reads_lo}\u2013{umi_rec_reads_hi}% to {umi_rec_umi_lo}\u2013{umi_rec_umi_hi}%. The two "
     "AID-null conditions are not plotted; their false calls fall from {umi_null_jh4} and {umi_null_n114} per "
     "library to zero. Grouping reads by tag leaves {umi_molecules} molecules from the 500 reads analysed, so the "
     "gain is bought with depth. (D) Recall of indels planted "
     "into reads of a held-out AID-null library, by the number of reads carrying each: 10 insertions or deletions of "
     "1–40 bp per draw, 8 draws per control, 240 planted indels per amplicon; vertical lines, the range of the three "
     "held-out controls; shading, the one- and two-read region the three-read floor excludes (STAR Methods). Legend, "
     "recall of the indels carried by three or more reads; above the panel, the false calls each cost. On J_{H}4 these "
     "are not eight events but one: the same 28-bp deletion, in the same four reads of one control, in every draw. It "
     "carries 5-bp microhomology at its ends, so it is most likely a molecule the material already held rather than a "
     "calling error."),
    ("Figure S4. N114 V-region results, related to Figure 4 (page 1 of 2)",
     [MC / "Figure4_hybridoma_combined.png"],
     "(A) N114 unique-variant density per subclone, n = 3 per vector. Empty-vector subclones without calls are plotted "
     "at one event. (B, C) Pooled substitution matrices for pMX-*Aicda* and empty-vector subclones, respectively. "
     "Cells report unique-variant counts and their shares of group calls. Panel A measures unique-variant density. "
     "Individual-library spectra are shown on the following page."),
    ("Figure S4. Per-library spectra, related to Figure 4 (page 2 of 2)",
     [MC / "FigureS6_individual_JH4_ONT_matrices.png", MC / "FigureS7_individual_N114_matrices.png"],
     "Substitution matrices for individual J_{H}4 Wt and *Aicda*^{−/−} libraries (top) and N114 pMX-*Aicda* and "
     "empty-vector subclones (bottom). For each component, A denotes AID-active libraries and B denotes AID-null "
     "libraries; panel letters restart in the lower component. Cells report unique-variant counts and their shares of "
     "the corresponding library's calls."),
    ("Figure S5. Sequencing depth and unique-variant discovery, related to Figure 4",
     FIG / "FigureS5.png",
     "Every library was analyzed at 30–1,000 reads, and the whole analysis was repeated five times at each depth with "
     "a different random draw of usable reads, the AID-null controls redrawn and recalibrated at that depth each time, "
     "so the spread shown is the spread a single run would have had. (A and B) J_{H}4 unique variants and "
     "unique-variant density against depth; (C and D) the same for N114. Points, the five draws of each library; "
     "lines, their mean; blue, AID-active libraries, red, AID-null. A library with no calls has a density of zero, "
     "which has no logarithm: those draws are drawn on the \"0\" row below the dotted break rather than floored to a "
     "small number, and they are common (J_{H}4 KO1 and KO3 called nothing at any depth in any draw). Density falls "
     "with depth because its numerator, distinct calls, saturates while its denominator keeps growing, so a density is "
     "only comparable at a stated depth. (E) What 500 reads, the depth used throughout, had reached in each AID-active "
     "library: calls at 500 reads as a percentage of that library's mean at 1,000, the largest depth tested; points, "
     "draws; bars, means; {depth_jh4_500}% (J_{H}4) and {depth_n114_500}% (N114) on average. This is a comparison with "
     "the largest depth tested here, not evidence that discovery has reached an asymptote. At 30 reads one draw is not "
     "informative: J_{H}4 Wt2 gave between 7 and 60 calls across its five draws, against 212 to 244 at 500 reads."),
    ("Figure S6. The background along the amplicon and in each read, related to Figure 4",
     FIG / "FigureS6.png",
     "(A and B) Mismatches at base quality ≥ 20 at every position of the analysis window, germline sites excluded, in "
     "each AID-null library (red) and averaged over the AID-active libraries (blue); J_{H}4, 500 reads per library; Sμ "
     "(PRJNA1481264), 2,000. Insets, the controls alone on their own scale, on which the small recurrent artifacts are "
     "visible. r, Spearman correlation between pairs of controls: the controls agree only weakly position by position "
     "(r = 0.12–0.20 over three pairs in J_{H}4 and 0.29–0.45 in Sμ). Recurrent artifacts are therefore a tendency of "
     "the chemistry and sequence context rather than a fixed per-position fingerprint, which is why the error model "
     "pools over substitution type, strand and base quality and shrinks each position towards that pooled rate instead "
     "of trusting its own count (STAR Methods). (C and D) Every mismatch at quality ≥ 20 (germline sites excluded) in "
     "each of the 500 reads of a wild-type (C) and an *Aicda*^{−/−} (D) J_{H}4 library; colored, at sites SHMtrace "
     "called; gray, other mismatches. Reads are sorted by clone (dotted lines between the largest clones) and by "
     "mutation load."),
    # The failure split (22 paralog amplicons, 12 short of a Q>=20 segment, 11 unmapped, three too
    # short) is still written out: the pipeline keeps only the traces it analysed, so the per-trace
    # gating reasons cannot be recomputed from its output. The total and the paralog count are the
    # two the archive manifest confirms independently.
    ("Figure S7. Sanger re-analysis, related to Figure 5 (page 1 of 2)",
     [MC / "FigureS5_individual_Sanger_matrices.png", CAL / "figures/sanger_vs_nanopore/figS_sanger_window_original.png"],
     "Top, per-mouse Sanger substitution matrices at Q≥20. Bottom, fraction of Wt clones with a high-quality base at "
     "each J_{H}4 position and calls per 10^{3} observed bases; shading, the shared 529-bp window. Of "
     "{sanger_archive} chromatograms in the eight archived sequencing orders, which cover that study’s wild-type, "
     "*Mkln1*-deficient and AID-null mice and also hold six traces of an unrelated amplicon, {sanger_sheet} were "
     "filed to the six mice used for nanopore sequencing; {sanger_excluded} of those were excluded in the sample "
     "sheet and not processed. Of the remaining {sanger_submitted}, {sanger_failed} failed analysis: 22 paralog amplicons, 12 "
     "without a contiguous "
     "20-base Q≥20 segment, 11 unmapped and three too short. The {sanger_clones} analyzed clones comprised "
     "{sanger_wt_clones} Wt clones ({sanger_wt_total} total) and {sanger_ko_clones} AID-null clones "
     "({sanger_ko_total} total). Bookkeeping checks flag "
     "potential problems without automatically removing clones."),
    ("Figure S7. Sanger bookkeeping checks, related to Figure 5 (page 2 of 2)",
     [FIG / "FigureS7_bookkeeping.png"],
     "Forty simulated projects each contained one instance of each bookkeeping error. Detection is shown against the "
     "planted errors; the shared-lineage test used a lineage represented by 125 clones. Label changes were detectable "
     "only when affected clones shared mutations with another animal's clones. No false flags occurred among 6,018 "
     "correctly assigned clones."),
    ("Figure S8. Reanalysis of published Maea Sanger data, related to Figure 5",
     [FIG / "FigureS8_sanger_cohort2.png"],
     "J_{H}4 clones of germinal-center B cells from four wild-type (blue), four *Maea*^{fl/fl}*Mb1*^{cre} (purple) and "
     "*Aicda*^{−/−} (red) mice[@barbulescu2024], each clone assigned to the mouse it was filed under in the original "
     "analysis. Filled circles and teal mean lines, SHMtrace from the chromatograms; open squares and brown mean "
     "lines, the published values, replotted from the paper's Source Data (Fig. 2g–i) and per-mouse tables "
     "(Supplementary Fig. 9b). The study reports three *Aicda*^{−/−} mice; four are filed in the original project "
     "and are all analyzed here (F). (A) Mutations per 10^{4} bp in each mouse; lines, means. (B) The same per "
     "mouse, SHMtrace against published; dashed line, identity. (C) Clones by mutations per clone, pooled by "
     "genotype, binned as published; the published bins come from the Fig. 2g Source Data, which count 96 and 183 "
     "clones, whereas the per-mouse tables behind Table S4 count 97 and 187. "
     "(D) G:C transitions and transversions as a share of each mouse's mutations. (E) Pooled 12-type substitution "
     "spectra, SHMtrace against published. (F) The published comparisons and SHMtrace's, with the one-tailed "
     "t-test used in the original study; \"background explains\" is the share of each group's calls that the "
     "*Aicda*^{\u2212/\u2212} clones' own error rate accounts for (STAR Methods), and the analysis time is for "
     "the whole project from one sample sheet, once that sheet had been built."),
    ("Figure S9. Reanalysis of published *Mkln1* Sanger data, related to Figure 5",
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
    ("Figure S10. Reanalysis of published SW_{HEL} Sanger sequences, related to Figure 5",
     [FIG / "FigureS10_swhel.png"],
     "The 1,795 Sanger clones of Thientosapol et al. (2018), who hypermutated a SW_{HEL} V_{H}10 knock-in V region "
     "in germinal-center B cells and reported that inactivating *Samhd1* raises transitions at C:G and lowers "
     "transversions at A:T and at C:G. SHMtrace was given only the deposited GenBank records: the germline is the "
     "iterated majority consensus of the clones themselves, each clone is assigned to the animal its own record "
     "names, and the arms are the genotypes its records state. (A) Mutations per sequence in each arm, SHMtrace "
     "against the paper's Table 1. (B) The 12 substitution types in each of the five arms, as a share of that "
     "arm's mutations. (C) The four classes the paper tests, one point per animal for SHMtrace and the published "
     "value as an open square; stars, two-tailed t-test across animals. (D) Each published claim against what "
     "SHMtrace finds; for total mutation (a) the two columns use different tests, the paper's one-way ANOVA and "
     "SHMtrace's Welch test on per-animal log frequency, so they agree in direction but not in inference. Unlike Figures S8 and S9 this study deposits sequences rather than chromatograms, so the "
     "chromatogram reading and the bookkeeping checks are not exercised; it has no *Aicda*^{\u2212/\u2212} arm, "
     "so the Sanger error rate has to be carried over from the AID-null clones of those two studies "
     "(about 20 in 10^{5} bases at Q \u2265 20 in each), roughly 10-fold below the mutation frequency "
     "measured here."),
    ("Figure S11. Full-length spectra, related to Figure 6 (page 1 of 2)",
     [MC / "FigureS_spectrum_529_vs_fulllength.png"],
     "Pooled Wt substitution spectra across the 529-bp shared window and the full 1,141-bp interval (Pearson r = "
     "{s1_r}). Mutation trees are shown on the following page."),
    ("Figure S11. Mutation trees in each mouse, related to Figure 6 (page 2 of 2)",
     [CAL / "figures/mutation_trees/mutation_trees_by_mouse.png"],
     "Mutation-tree reconstructions for the most branched clone in each Wt mouse. Trees summarize relationships among "
     "sequenced reads (Figure 6); the corresponding cellular lineages were not independently established."),
    ("Figure S12. Sequencing units, coverage and Illumina read preparation, related to Figures 5–7",
     [FIG / "FigureS12_platforms.png"],
     "Top, Sanger and nanopore sequencing units drawn on a common base-pair scale. A Sanger clone represents one "
     "PCR-derived molecule sequenced over the shared 529-bp window. The eight nanopore reads are experimental Wt3 "
     "reads from four clones spanning the 1,230-bp amplicon, with called mutations shown over positions 90–1,230. "
     "Error marks are simulated because uncalled mismatches cannot individually be classified as technical errors or "
     "true mutations. They were placed randomly at the AID-null mismatch rate after germline exclusion (0.14% of "
     "bases); 85% fall below Q20 (pale orange), with the remainder at Q≥20 (solid orange). Both are smaller than "
     "mutation markers. Individually unresolved single-read mutations are estimated in aggregate. Middle, overlapping "
     "Illumina pairs are merged before read preparation. Examples show median pairs among the first 20,000 pairs from "
     "a Wt J_{H}4 mouse (PRJNA973112; 95.5% merged) and 14,679 pairs from an AID-expressing CH12F3 replicate "
     "(PRJNA918596; 96.7% merged). Agreeing overlap bases have summed qualities capped at Q60; disagreements retain "
     "the higher-quality call or become N at tied qualities. Adapter read-through is removed, and pairs with <20-bp "
     "overlap are discarded. Overlap disagreements occurred at 0.13% and 0.55% of bases, respectively. Right, coverage "
     "of each study locus by a single processed read. Merged reads enter the Figure 1 workflow using the minimap2 sr "
     "preset, an expected-error limit of 0.02 per base and a minimum alignment requirement of no more than 100 bp. "
     "PRJNA973112 and PRJNA918596 contain uniformly Q30 bases, yielding one quality level in the error model; "
     "PRJNA820976 retains binned qualities (`crm/scicanvas/s12_illumina.py`). Bottom, measurements obtained from the "
     "two workflows in this study and their agreement in the six matched mice. These comparisons are "
     "assay-specific: Sanger amplicons spanning V(D)J junctions can support lineage grouping, UMI consensus represents "
     "one molecule per family, and long-read coverage depends on amplicon design. Nanopore fold separation assigns "
     "zero-call controls one event. Sanger background is the Q≥20 mismatch rate in AID-null clones from Figures S8 and "
     "S9."),
    ("Figure S13. SHMtrace analysis workflow, related to STAR Methods",
     [FIG / "FigureS13_overview.png"],
     "Sequence files are assigned to biological replicates through a sample sheet, filtered and aligned competitively "
     "against target, paralog and contaminant references. Matched AID-null libraries provide the background model for "
     "read-library substitution tests; held-out controls serve as decoys for threshold selection. Sanger substitutions "
     "are counted separately, with measured background reported alongside the counts (Figure 1C). Length limits are "
     "locus-specific (STAR Methods). The default ≥50% allele filter excludes likely germline differences but can also "
     "remove acquired mutations in dominant clones. Group comparisons use one value per animal. Illustrative panels "
     "show the error-rate grid by substitution and quality, expected versus observed candidate support, and calls "
     "versus decoys across thresholds. Outputs include density, mutation frequency, base-composition-corrected 12-type "
     "spectra, WRC/GYW and WA/TW targeting, positional profiles, clones, mutation trees and indels; Sanger outputs "
     "also report expected background counts. Bookkeeping flags do not automatically remove traces. A trace "
     "deliberately assigned to two animals is retained in both, and suspected accidental copies are flagged for "
     "review. Only rows explicitly marked as excluded in the sample sheet are omitted, with reasons recorded. "
     "Sample-sheet rows are illustrative. Figure S14 details the analytical rules and equations."),
    ("Figure S14. Rules underlying read preparation and variant calling, related to STAR Methods",
     [FIG / "FigureS14_engine.png"],
     "Six analytical stages are illustrated with J_{H}4 libraries. (A) Per-library read accounting: length and "
     "expected-error filtering, competitive alignment and usable reads (individual filters in Figure S1). (B) Target "
     "versus paralog alignment identity; assignment follows the better match. (C) Reported base quality versus "
     "*Aicda*^{−/−} mismatch rate; dashed line, nominal quality-derived error probability. (D) Two candidates at equal "
     "depth in the same library, showing supporting reads in Wt and control libraries and their different calling "
     "outcomes. (E) Candidate and held-out-decoy counts across p-value cutoffs, with the selected threshold. (F) "
     "Per-animal unique-variant density used for group comparisons. Panels use `intuition_jh4.json`, per-library "
     "read-accounting records and the bundled example library; equations describe the implemented calculations."),
    ("Figure S15. End-to-end stress tests using resampled control reads, related to STAR Methods",
     FIG / "FigureS15.png",
     "AID-null J_{H}4 reads were resampled and analyzed with `run_project`, including margin selection, decoy "
     "thresholding and noisy-control exclusion. Cases used KO1. Two-control analyses used one partition each from KO2 "
     "and KO3; three-control analyses used two KO2 partitions and one KO3 partition. Case and control read records "
     "did not overlap. These are conditional resampling experiments using existing libraries, rather than independent "
     "preparations or a new sequencing run. Unless stated otherwise, analyses used three controls and 500 reads. Teal, "
     "exclusion enabled; amber, exclusion disabled, paired on identical inputs. (A) All-null resamples with any call "
     "and 95% Clopper–Pearson intervals; n = {fdv_null_none_on_n} per setting for unmodified background and one noisy control, and n = {fdv_null_case_only_on_n} "
     "for case-only and shared artifacts. Any discovery gives an all-null false-discovery proportion (FDP) of 1. (B) "
     "Mean false calls ± SEM. (C) Individual FDPs in mutation-bearing resamples, five per condition and setting; black "
     "bars, means; dashed line, nominal 5% target. (D) Exclusion rate under a heterogeneous null; n = 5 per setting "
     "with two controls at 250 or 1,000 reads, and n = 25 with three controls at 500 reads. Artifacts occupied six "
     "positions with nominal 6% support at each. Shared artifacts used identical alternate alleles and avoided planted "
     "biological mutations. Injection preserved read orientation and quality order. Shared all-null artifacts produced "
     "calls in {fdv_null_shared_on_any}/{fdv_null_shared_on_n} resamples with exclusion enabled and "
     "{fdv_null_shared_off_any}/{fdv_null_shared_off_n} disabled, versus "
     "{fdv_null_case_only_on_any}/{fdv_null_case_only_on_n} in both case-only settings. Small resample counts limit "
     "conclusions about error control (STAR Methods)."),
    ("Figure S16. Class-switch junctions in a public PacBio dataset, related to Figure 7",
     FIG / "FigureS16_switch_pacbio.png",
     "Sμ–Sγ1 products from Wt, *Shld1*^{−/−}, *Xlf*^{−/−} and *Shld1*^{−/−}*Xlf*^{−/−} primary B cells, three mice per "
     "genotype and one PacBio run per mouse (PRJNA831666). (A) Delivered reads, junction-bearing reads and unique "
     "junctions per library. Unique junctions group breakpoints within 3 bp on both sides. Reads per unique junction "
     "vary several-fold between samples; counts therefore measure long-range PCR-product recovery rather than "
     "switching efficiency. (B) Unique-junction microhomology classes, averaged across three mice: blunt, 1–3 nt, ≥4 "
     "nt and junctions with inserted bases. (C) Mean microhomology per mouse; bars, genotype means; p, Welch tests "
     "against Wt, n = 3 per genotype. Permutation p values cannot fall below 0.1 at this sample size and are provided "
     "in `crm/out/switch_pacbio/comparisons.csv`. (D) Fraction of unique junctions with an alternative "
     "acceptor-segment placement within five matched bases of the best alignment, assessed by independent segment "
     "realignment. The module's internal ambiguity flag is 0% throughout because minimap2 reports no secondary "
     "alignments; breakpoints within repeats therefore remain approximate. Removing ambiguous junctions excludes 2–21% "
     "per sample and preserves the classification relative to p = 0.05 for 13 of 15 genotype–feature comparisons. "
     "*Shld1*^{−/−} ≥4-nt microhomology becomes stronger and *Xlf*^{−/−} insertions weaker (STAR Methods); not all 15 "
     "comparisons appear in panel C. The dataset lacks AID-null controls and does not evaluate substitution-test "
     "calibration."),
    ("Figure S17. Public S\u03bc mutation libraries under PI3K\u03b4 blockade, related to Figure 7",
     FIG / "FigureS17_public_smu.png",
     "Sμ amplicons from Wt and *Aicda*^{−/−} primary B cells treated with DMSO, idelalisib or duvelisib, three "
     "replicates per genotype and treatment (GSE77788). (A) Library coverage of the 651-bp reference, whose first half "
     "contains switch repeats. Six of nine Wt libraries are 300-bp single-end, whereas all *Aicda*^{−/−} libraries are "
     "paired-end, confounding genotype with library layout. Mates do not overlap. Window B is covered by one read from "
     "every library; window C covers repeats accessible only through the second mate. (B) Unique mutations per library "
     "in window B at 20,000 reads; bars, group means. *Aicda*^{−/−} libraries remain at background across treatments. "
     "Wt means increase with both inhibitors, but within-treatment variation exceeds between-treatment differences, "
     "leaving the drug effect unresolved in this window. (C) Descriptive genotype comparison in window C using the "
     "twelve paired libraries (three Wt and nine *Aicda*^{−/−}). Each Wt treatment has only one replicate, so pooled "
     "treatment groups do not provide a balanced treatment comparison. Targeting index and C→T/G→A share refer to Wt "
     "calls."),
    ("Figure S18. Engineered AID diversification in four public datasets, related to Figure 7",
     FIG / "FigureS18_engineered.png",
     "(A) One GFP substrate edited by a panel of deaminases delivered alike (SRP159221), three replicates per "
     "arm: the share of calls at WRC/GYW (circles) and the share that are C\u2192T or G\u2192A (squares), in "
     "wild-type and in UNG/MSH2-null 293T. Bars, group means. The motif share separates AID from APOBEC3A in "
     "both backgrounds; the transition share does so only in wild-type cells, because losing uracil excision "
     "drives transitions whichever enzyme made the uracil. (B) Mismatches per 1,000 bases at Q \u2265 30 along "
     "one ABL exon-6 amplicon, for five libraries differing only in which sgRNA they carry (SRP075898). Each "
     "track is scaled to its own peak, whose height is given, so a weak guide stays visible; the arm is "
     "multiplexed and all three of its amplicons were analysed. Shading marks the + strand protospacers, "
     "located by taking each PAM drawn in Ma et al.'s Supplementary Figure 9 and reading the 20 bases 5' "
     "of it; that figure does not number its guides, so the bands are unlabelled. Every library's highest "
     "peak falls inside one. #6, #9 and #15 peak at distinct positions at "
     "116\u2013137 per 1,000 against a background near 1. #4 passes a fold test on a low local baseline but "
     "reaches only 6 per 1,000, and #7's peak lies where two other guides are also high, so its fold "
     "understates it. These are per-position rates, not mutation calls. (C) Positions where a library sits at "
     "least five-fold above the lowest of the other three at that base, against the number of sgRNAs the "
     "library carries. The baseline leaves each library out of its own comparison. All four pools target the "
     "same 100-bp window of ABL exon 6 (positions 201–300), and only seven guides fit inside it, so the "
     "ten-guide pool adds three that flank it[@ma2016]; the zero at three guides is therefore too little "
     "signal to clear the five-fold threshold, not guides aimed elsewhere. (D) An AID-diversified MAGE-A1 TCR (figshare "
     "12803501) called against the 1,785-bp ancestral construct reconstructed from the deposit: the one "
     "available baseline library against ten diversified and selected populations. With a single baseline, and "
     "that baseline belonging to one founder, no false-discovery rate is claimed. Analyses in "
     "`crm/liu_mutators.py`, `crm/tam_footprint.py`, `crm/tam_saturation.py` and `crm/tcr_bassan.py`."),
]

if __name__ == "__main__":
    main(sys.argv[1])
