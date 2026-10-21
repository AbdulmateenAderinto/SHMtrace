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
# the combined review document's version (write_complete): raise it whenever a version goes out to readers
VERSION = "1.1"


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
                      "remains close to the background even at the full depth available (Table S1).")
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
    # zeros at one event on the mouse's own exposure, the paper's one zero-count rule (replicates.ZERO_EVENTS);
    # dropping them, as this did, removed the Ldha-deleted mouse with no calls from the estimate
    def logf(grp):
        return np.log10([max(s_["mutation_frequency_x1e3"], 1e3 / s_["frequency_bases"])
                         for s_ in p["samples"] if s_["group"] == grp])
    lw, lk = logf("WT"), logf("LDHA_KO")
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
    # clones at the same depth: a clone with no deletions counts one (the one zero-count rule), not k + 1
    n["hmces_del_p"] = crm_numbers.sci(stats.ttest_ind(np.log10(np.maximum(dels['HMCES_KO'], 1)),
                                                       np.log10(np.maximum(dels['WT'], 1)), equal_var=False).pvalue)
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
    """The loci and laboratory counts the Summary, Discussion and Table S2 footnote all quote, taken
    from the table itself through M.TABLE_S6_LOCI rather than written out in three places.

    Every row of the table must be classified and every classification must name a row, so adding a
    dataset fails the build instead of leaving the counts stale. A locus is counted once however many
    rows share it, and a row whose locus is not an AID target counts towards neither total."""
    rows = {r[0] for r in M.TABLE_S6}
    declared = set(M.TABLE_S6_LOCI)
    assert rows == declared, (f"Table S2 and its locus classification disagree; rows not classified: "
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
    # Figure S20's read-span bars are drawn from lengths typed into the builder, so the PacBio one
    # is checked against the switch data the legend quotes
    s12 = (CRM / "scicanvas/figS20_platforms.js").read_text()
    pb = n["pacbio_len"].replace(",", "")
    if f", {pb}, PACB]" not in s12:
        raise RuntimeError(f"Figure S20 is stale: expected the PacBio bar at {pb} bp in "
                           "scicanvas/figS20_platforms.js")


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


def write_article(m, n):
    """The article from its title page to the end of the STAR Methods, without the reference list.

    Returns the main-text blocks the word count reads. The submitted article and the combined review
    document both call this, so they cannot drift apart."""
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
        m.caption(title, fill(legend, n))
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
    return main_blocks


def write_supplement(s):
    """Every supplementary figure and table, in number order (Document S1, and the combined document)."""
    # in number order, which renumber_supp.py keeps equal to the order the article first cites them
    for title, img, legend in sorted(SUPP_FIGS, key=lambda f: int(re.match(r"Figure S(\d+)", f[0])[1])):
        s.page_break()
        imgs = img if isinstance(img, list) else [img]
        s.figure_block(imgs, title, legend)
    # Table S2: the dataset inventory. It was Table 1 while the main text had room for a table; with
    # seven figures, figures and tables together reach the journal's limit, so it travels here.
    s.page_break()
    s.table(M.TABLE_S6_HEADER, M.TABLE_S6, title=f"**{M.TABLE_S6_TITLE}**", footnote=M.TABLE_S6_NOTE,
            widths=[28, 20, 28, 32, 28, 20, 18], size=7.5)
    # Table S3: what SHMtrace did with every delivered read, and the calls (J_H4 and N114)
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
            title="**Table S3. Read preparation and calls per library, related to Figure 4**",
            footnote="Per-library read removal at each preparation step, starting from delivered FASTQ files. Five "
                     "hundred usable reads were randomly selected (seed 0). All three mutation endpoints are reported "
                     "(STAR Methods): unique-variant density counts each distinct call once per library over "
                     "analyzed base observations; read-weighted frequency counts its support across reads; "
                     "clone-weighted frequency gives each junction-defined clone equal weight. A dash indicates that "
                     "no clone had an assignable junction. WT, wild type; KO, *Aicda*^{−/−}; C1, D3 and D5, "
                     "pMX-*Aicda* subclones; C3, C4 and D1, empty-vector subclones.")
    # Table S4: SHMtrace against the two published Sanger analyses, metric by metric (sanger_compare.py)
    sys.path.insert(0, str(CRM))
    import sanger_compare
    s.page_break()
    s.table(sanger_compare.HEADER, sanger_compare.table_rows(), size=7.5,
            title="**Table S4. SHMtrace and two published Sanger analyses, metric by metric, related to Figures S13 "
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
            title="**Table S5. What \"reproduced\" means for each re-analyzed study, related to Figures S13–S15**",
            footnote=agreement_table.NOTE)
    # Table S6: statistics
    s.page_break()
    s.table(["Comparison", "Test", "Result"], STATS, size=8.5,
            title="**Table S6. Statistical summary, related to Figures 3–7**",
            footnote="With three animals per group, p values and intervals should be interpreted descriptively. Fold "
                     "separation depends on the one-event convention used for zero-call controls and does not "
                     "represent a statistical lower confidence bound.")
    # Table S7: the default threshold against two stricter settings (strict_sensitivity.py)
    st_path = OUT / "strict_sensitivity.json"
    if st_path.exists():
        st = json.load(open(st_path))
        prof = (("default", "Default"), ("strict", "α 0.01"), ("strict_max", "Noisiest decoy"))
        rows = []
        for ds, lab in (("jh4", "J_{H}4"), ("n114", "N114")):
            for k, name in prof:
                r = st["spikein"][ds][k]
                rows.append((f"Spike-ins, {lab}", name, f"{r['false_per_lib']:.1f}", f"{r['fdr']:.3f}",
                             f"{100 * r['recall']:.0f}%", f"{100 * r['recall_1read']:.0f}%"))
        for sc, lab in (("polyclonal", "polyclonal"), ("expanded", "expanded clone"), ("low", "low load"),
                        ("null", "null")):
            for k, name in prof:
                r = st["lineage"][sc][k]
                one = "–" if r["recall_1read"] is None else f"{100 * r['recall_1read']:.0f}%"
                fdr = "1 (all false)" if sc == "null" else f"{r['fdr']:.3f}"
                rows.append((f"Lineages, {lab}", name, f"{r['false_per_lib']:.1f}", fdr, "–", one))
        s.page_break()
        s.table(["Benchmark", "Setting", "False calls per library", "FDR", "Recall", "Single-read recall"], rows,
                widths=[34, 26, 24, 22, 18, 24], size=7.5, title=f"**{M.TABLE_S7_TITLE}**")
        prow = []
        for ds, lab in (("jh4", "J_{H}4"), ("n114", "N114")):
            for k, name in prof:
                r = st["projects"][ds][k]
                prow.append((lab, name, ", ".join(map(str, r["case_calls"])), ", ".join(map(str, r["control_calls"])),
                             crm_numbers.pfmt(r["p_log_welch"]), f"{r['fold']:.0f}"))
        s.table(["Project", "Setting", "Calls, mutation-bearing", "Calls, controls", "p (log-Welch)", "Fold"], prow,
                widths=[22, 26, 34, 26, 22, 18], size=7.5, footnote=M.TABLE_S7_NOTE)


def version_label():
    """The combined document's version, with the date and the commit it was built from."""
    import datetime
    import subprocess
    git = ["git", "-C", str(Path(__file__).resolve().parent)]
    short = subprocess.run([*git, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run([*git, "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                           text=True).stdout.strip()
    when = datetime.date.today()
    return (f"{VERSION}, {when.day} {when:%B %Y}, commit {short or 'unknown'}"
            + (" with uncommitted changes" if dirty else ""))


def write_complete(out, n):
    """One review document: the article, its figures and STAR Methods, every supplementary figure and
    table, and one reference list for both. Figure, table and reference citations are internal links."""
    label = version_label()
    c = Manuscript(line_numbers=False, links=True, footer_note=f"SHMtrace, complete manuscript, version {VERSION}")
    c.para("SHMtrace: complete manuscript with supplemental information", size=14, bold=True)
    c.para(f"Version {label}", size=11)
    c.para("Every figure, table and reference citation in this document is a link to its target. The "
           "article is followed by Document S1 (Figures S1–S23 and Tables S2–S7) and one reference list "
           "for both. Table S1 and Data S1 are separate files.", size=10)
    c.page_break()
    write_article(c, n)
    c.page_break()
    c.heading("Document S1. Supplemental figures and tables", 1)
    write_supplement(c)
    c.references()
    path = out / f"SHMtrace_CellReportsMethods_complete_v{VERSION}.docx"
    c.save(path)
    print(f"{path.name} written (version {label})")


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    n = {**crm_numbers.compute(), **public_numbers(), **dataset_counts()}
    fill_modules(n)
    check_figure_numbers(n)
    check_paired_settings(n)
    sw = check_limits()
    m = Manuscript()
    main_blocks = write_article(m, n)
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
    # Document S1 carries no reference list of its own and its footnotes say the superscripts cite
    # the manuscript's. Seed its citation order from the manuscript so they actually do: numbered
    # independently, every supplemental superscript pointed at the wrong paper.
    s.order = list(m.order)
    s.para("Supplemental information", size=14, bold=True)
    s.para(T.TITLE, size=11)
    write_supplement(s)
    s.save(out / "Document_S1.docx")
    print("Document S1 written")
    write_complete(out, n)
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
     "The workflow of Figure 1 in full, from the design of the two experiments to the fate of every delivered read. "
     "(A) Experimental design for the J_{H}4 and N114 libraries. (B, C) Read preparation and variant calling from "
     "the delivered sequence files, which is the read-library route of Figure 1C. By default, an allele carried by "
     "50% or more of the reads is excluded as a likely germline difference; because a mutation driven to that "
     "frequency within a dominant clone is removed with it, this is a heuristic and not a test (STAR Methods). (D) "
     "The Sanger route counts every substitution at base quality ≥20 across the analysed clones and the 529-bp "
     "window. The mismatch rate measured in the AID-null clones is reported beside those counts rather than "
     "subtracted from them, because it estimates the background of the whole set and not of any one call. A "
     "sequencing clone is one PCR-derived molecule sequenced once, so it says nothing about which B cell the "
     "molecule came from. (E) What happened to every delivered read in each library; the line marks the 500 reads "
     "drawn for analysis. Experimental details: as in Figure 5 for J_{H}4 and Figure 4 for N114."),
    ("Figure S2. Reported base qualities and observed control mismatch rates, related to Figure 2",
     FIG / "FigureS2.png",
     "A reported base quality is a prediction, and this is how far it holds at each of the four nanopore loci. Top, "
     "the mismatch rate measured in the AID-null controls, grouped into two-quality-score bins holding at least 2,000 "
     "bases after germline positions are excluded. Shading, 95% Clopper–Pearson intervals; downward triangles, bins "
     "with no mismatch, drawn at the upper confidence bound because a rate of zero cannot be placed on a logarithmic "
     "axis; dashed line, the error probability the quality score nominally implies, 10^{−Q/10}. Bottom, the number of "
     "control bases in each bin. Observations thin out above Q40, so the uncertainty grows there. This is why each "
     "locus is calibrated against its own controls rather than against the reported qualities. Experimental details: "
     "the J_{H}4 libraries were sequenced under order DHJGKB and the N114 libraries under order 75K3B3; for Sμ and "
     "*Pim1* (PRJNA1481264) the deposited records do not establish which libraries shared a flow cell."),
    ("Figure S7. Single-read estimate, power, UMI consensus and planted indels, related to Figure 3",
     FIG / "FigureS7.png",
     "Four checks on what the read-library route can and cannot resolve. (A) Change in the aggregate single-read "
     "estimate against one-read mutations planted into held-out J_{H}4 control libraries (three controls; three "
     "replicates per amount); dashed line, identity. (B) Power at α = 0.05, two-sided, from the Welch test on log "
     "rates, by animals per group and by the size of the difference sought, simulated from the J_{H}4 wild-type "
     "pilot (3 mice, 500 reads each; 1,000 experiments per point; negative-binomial counts with the pilot's "
     "clumping, animals log-normal with the pilot's between-animal coefficient of variation, CV 0.95). Dotted "
     "line, 80% power. How many animals are needed depends on the difference being sought, not on the method "
     "alone: a 2-fold difference between two AID-proficient groups needs about 24 per group at this CV, whereas at "
     "10-fold three animals give 63% power and five give 97%, and at 50-fold three give 97%. The comparisons this "
     "paper makes are AID-competent against AID-null and run 38–141-fold, which is why they separate at three "
     "animals per group. The CV is itself estimated from three mice and so is the loosest term in the "
     "calculation. (C) Simulated UMI libraries built from real control reads, three held-out controls each: false "
     "calls per library when every read is used and when reads are first collapsed to a consensus by tag. The four "
     "conditions are the two amplicons under two planted mutation loads: J_{H}4 polyclonal (about 40 clones) and "
     "J_{H}4 low load, and N114 V_{H} subclonal (one dominant clone) and N114 V_{H} low load, the same scenarios "
     "as Figure 3D. Consensus removes every false call in three of the four and leaves {umi_fp_poly} in the "
     "fourth, annotated on the panel because a zero-height bar cannot be told from a missing one, while recall of "
     "the planted mutations rises from {umi_rec_reads_lo}–{umi_rec_reads_hi}% to {umi_rec_umi_lo}–"
     "{umi_rec_umi_hi}%. The two AID-null conditions are not plotted; their false calls fall from {umi_null_jh4} "
     "and {umi_null_n114} per library to zero. Grouping reads by tag leaves {umi_molecules} molecules from the 500 "
     "reads analysed, so the gain is bought with depth. (D) Recall of indels planted into reads of a held-out "
     "AID-null library, by the number of reads carrying each: 10 insertions or deletions of 1–40 bp per draw, 8 "
     "draws per control, 240 planted indels per amplicon; vertical lines, the range of the three held-out "
     "controls; shading, the one- and two-read region the three-read floor excludes (STAR Methods). Legend, recall "
     "of the indels carried by three or more reads; above the panel, the false calls each cost. On J_{H}4 these "
     "are not eight events but one: the same 28-bp deletion, in the same four reads of one control, in every "
     "draw. It carries 5-bp microhomology at its ends, so it is most likely a molecule the material already held "
     "rather than a calling error."),
    ("Figure S8. N114 V-region results, related to Figure 4 (page 1 of 2)",
     [FIG / "FigureS8_n114.png"],
     "AID expression drives mutation of the fixed immunoglobulin V-region substrate carried by N114 hybridoma cells. "
     "(A) Unique-variant density in each subclone (n = 3 per vector). An empty-vector subclone with no call is drawn "
     "at one event on its own exposure, because a density of exactly zero cannot be placed on a logarithmic axis. "
     "(B, C) Pooled From (rows) → To (columns) substitution matrices for the pMX-*Aicda* subclones (B) and the "
     "empty-vector subclones (C). Each cell gives the unique-variant count and its share of that group's calls. "
     "Because N114 cells carry a fixed V-region substrate, do not undergo ongoing V(D)J recombination and are not "
     "selected here for antigen binding, the excess in pMX-*Aicda* subclones reflects AID-dependent mutation rather "
     "than rearrangement or selection. That substrate carries an opal stop in CDR2, which silences antibody "
     "production and whose reversion the original study reads out as hypermutation[@martin2002]; the calls "
     "here are made against the line's own sequence, stop included, and count mutation across the whole V "
     "region rather than reversion at that codon (STAR Methods). Experimental details: N114 cells were transduced with pMX-*Aicda* or matched "
     "empty vector, both encoding GFP and puromycin resistance; GFP-high cells were subcloned by limiting dilution "
     "and three puromycin-resistant subclones per vector were maintained for 1 month before DNA collection. The "
     "~0.65 kb V region was amplified with Platinum SuperFi II high-fidelity polymerase and sequenced as "
     "individually barcoded R10.4.1/Kit 14 libraries, analysed at 500 reads per subclone. A unique variant is one "
     "substitution at one reference position, counted once per subclone however many reads carry it. Spectra of the "
     "individual libraries are shown on the following page."),
    ("Figure S8. Per-library spectra, related to Figure 4 (page 2 of 2)",
     [FIG / "FigureS8_spectra_jh4.png", FIG / "FigureS8_spectra_n114.png"],
     "The same spectra resolved to one matrix per library, so that the pooled patterns of Figure 4 and the preceding "
     "page can be seen to hold in each animal and subclone rather than in the pool alone. Top, J_{H}4 wild-type "
     "(*Wt*) and *Aicda*^{−/−} libraries; bottom, N114 pMX-*Aicda* and empty-vector subclones. Within each component "
     "(A) denotes the AID-active libraries and (B) the AID-null libraries, and the panel letters restart in the lower "
     "component. Matrices are oriented From (rows) → To (columns), and each cell gives the unique-variant count and "
     "its share of that library's calls. Experimental details: the libraries are those of Figures 4 and 5, analysed "
     "at 500 reads each."),
    ("Figure S9. Sequencing depth and unique-variant discovery, related to Figure 4",
     FIG / "FigureS9.png",
     "Every library was analyzed at 30–1,000 reads, and the whole analysis was repeated five times at each depth with "
     "a different random draw of usable reads, the AID-null controls redrawn and recalibrated at that depth each time, "
     "so the spread shown is the spread a single run would have had. (A and B) J_{H}4 unique variants and "
     "unique-variant density against depth; (C and D) the same for N114. Points, the five draws of each library; "
     "lines, their mean; blue, AID-active libraries, red, AID-null. A library with no calls has a density of zero, "
     "which has no logarithm: those draws are drawn on the \"0\" row below the dotted break rather than floored to a "
     "small number, and they are common (J_{H}4 KO1 and KO3 called nothing at any depth in any draw). Density falls "
     "with depth because its numerator, distinct calls, saturates while its denominator keeps growing, so a density is "
     "only comparable at a stated depth. (E) How much of each AID-active library's own variant set each depth "
     "recovered. The denominator is the most that library ever yielded, over every depth and every draw, so a library "
     "is read against its own discovery rather than against its 1,000-read value alone; lines, the range of five "
     "draws; vertical line, the 500 reads used throughout. Discovery rises steeply and then flattens, as further "
     "reads re-sample positions already seen rather than revealing new ones, which is what saturation looks like for "
     "a finite set of truly mutated positions. The J_{H}4 libraries are near their own total by 500 reads "
     "({depth_jh4_500}% on average), whereas the N114 V_{H} libraries, which carry fewer mutations, are still "
     "climbing ({depth_n114_500}%). At 30 reads one draw is not informative: J_{H}4 Wt2 gave between 7 and 60 calls "
     "across its five draws, against 212 to 244 at 500 reads."),
    ("Figure S10. The background along the amplicon and in each read, related to Figure 4",
     FIG / "FigureS10.png",
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
    ("Figure S12. Sanger re-analysis, related to Figure 5 (page 1 of 2)",
     [FIG / "FigureS12_sanger_matrices.png", FIG / "FigureS12_sanger_window.png"],
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
    ("Figure S12. Sanger bookkeeping checks, related to Figure 5 (page 2 of 2)",
     [FIG / "FigureS12_bookkeeping.png"],
     "Archived Sanger projects accumulate clerical errors, so SHMtrace flags them; this measures how often it "
     "succeeds. Forty simulated projects were built, each carrying one instance of every bookkeeping error, and "
     "detection is shown against the errors planted. The shared-lineage test used a lineage represented by 125 "
     "clones. A changed label was detectable only when the affected clones shared mutations with another animal's "
     "clones, because a label on its own leaves no trace in the sequence. No correctly assigned clone was flagged in "
     "error among the 6,018 tested, so the flags are specific as well as sensitive. The flags report a clone for "
     "review and never remove it automatically (STAR Methods)."),
    ("Figure S13. Reanalysis of published Maea Sanger data, related to Figure 5",
     [FIG / "FigureS13_sanger_cohort2.png"],
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
    ("Figure S14. Reanalysis of published *Mkln1* Sanger data, related to Figure 5",
     [FIG / "FigureS14_sanger_mkln.png"],
     "J_{H}4 clones of germinal-center B cells from three wild-type (blue), three *Mkln1*^{\u2212/\u2212} "
     "(purple) and three *Aicda*^{\u2212/\u2212} (red) mice[@barbulescu2025], each clone assigned to the mouse "
     "folder it was filed under in the original analysis. Filled symbols, SHMtrace from the chromatograms; open "
     "squares, the published values, replotted from the paper's per-mouse tables (Supplementary Figure S2A). "
     "(A) Unique mutations per 10^{4} bases sequenced, per mouse; lines, means. (B) G:C transitions and "
     "transversions as a share of each mouse's mutations. (C) Pooled 12-type substitution spectra, SHMtrace "
     "against published, with r for each genotype. Of the three published *Mkln1*^{\u2212/\u2212} tables, only "
     "mouse #1's spectrum matches the clones of a folder (K2); tables #2 and #3 match no folder or set of clones in "
     "the original project, which limits the knockout's agreement. (D) The comparisons, published and SHMtrace's; "
     "the published P values are computed "
     "from the per-mouse tables as printed, and the values the paper's summary bars report (Supplementary "
     "Figure S2C) are listed separately. Paralog clones are J_{H}4-paralog amplicons, which SHMtrace analyzes "
     "apart from J_{H}4; \"background explains\" is as in Figure S13."),
    ("Figure S15. Reanalysis of published SW_{HEL} Sanger sequences, related to Figure 5",
     [FIG / "FigureS15_swhel.png"],
     "The 1,795 Sanger clones of Thientosapol et al. (2018), who hypermutated a SW_{HEL} V_{H}10 knock-in V region "
     "in germinal-center B cells and reported that inactivating *Samhd1* raises transitions at C:G and lowers "
     "transversions at A:T and at C:G. SHMtrace was given only the deposited GenBank records: the germline is the "
     "iterated majority consensus of the clones themselves, each clone is assigned to the animal its own record "
     "names, and the arms are the genotypes its records state. (A) Mutations per sequence in each arm, SHMtrace "
     "against the paper's Table 1. (B) The 12 substitution types in each of the five arms, as a share of that "
     "arm's mutations. (C) The four classes the paper tests, one point per animal for SHMtrace and the published "
     "value as an open square; stars, two-tailed t-test across animals. (D) Each published claim against what "
     "SHMtrace finds; for total mutation (a) the two columns use different tests, the paper's one-way ANOVA and "
     "SHMtrace's Welch test on per-animal log frequency, so they agree in direction but not in inference. Unlike "
     "Figures S13 and S14 this study deposits sequences rather than chromatograms, so the "
     "chromatogram reading and the bookkeeping checks are not exercised; it has no *Aicda*^{\u2212/\u2212} arm, "
     "so the Sanger error rate has to be carried over from the AID-null clones of those two studies "
     "(about 20 in 10^{5} bases at Q \u2265 20 in each), roughly 10-fold below the mutation frequency "
     "measured here."),
    ("Figure S16. Full-length spectra, related to Figure 6 (page 1 of 2)",
     [FIG / "FigureS16_spectra_full_length.png"],
     "Extending the analysis from the Sanger window to the full amplicon increases the number of mutations sampled "
     "without altering the signature they carry. Pooled wild-type (*Wt*) substitution spectra across the shared "
     "529-bp window and across the reliably covered 1,141-bp interval, with each cell giving the mutation count and "
     "its share of all mutations; the two sets of proportions are highly correlated (Pearson r = {s1_r}). "
     "Experimental details: the three *Wt* mice of Figure 5, analysed at 500 reads per mouse. Mutation trees from "
     "the same libraries are shown on the following page."),
    ("Figure S16. Mutation trees in each mouse, related to Figure 6 (page 2 of 2)",
     [FIG / "FigureS16_mutation_trees.png"],
     "The read-level reconstruction of Figure 6D, repeated in every mouse. Each tree is the most branched clone of "
     "one wild-type (*Wt*) mouse, built from the mutations its reads share. A tree summarises relationships among "
     "sequenced reads and not among cells: the reads of one clone may derive from the same molecule amplified more "
     "than once, so the corresponding cellular lineages were not independently established. Experimental details: "
     "the three *Wt* mice of Figure 5, analysed at 500 reads per mouse."),
    ("Figure S20. What a read means on each of the four platforms, related to Figures 5–7",
     [FIG / "FigureS20_platforms.png"],
     "Top row, the two platforms applied to the same six mice, drawn on a common base-pair scale. A Sanger clone represents one "
     "PCR-derived molecule sequenced over the shared 529-bp window. The eight nanopore reads are experimental Wt3 "
     "reads from four clones spanning the 1,230-bp amplicon, with called mutations shown over positions 90–1,230. "
     "Error marks are simulated because uncalled mismatches cannot individually be classified as technical errors or "
     "true mutations. They were placed randomly at the AID-null mismatch rate after germline exclusion (0.14% of "
     "bases); 85% fall below Q20 (pale orange), with the remainder at Q≥20 (solid orange). Both are smaller than "
     "mutation markers. Individually unresolved single-read mutations are estimated in aggregate. Second row, the "
     "two platforms read from public deposits: overlapping "
     "Illumina pairs are merged before read preparation. Examples show median pairs among the first 20,000 pairs from "
     "a Wt J_{H}4 mouse (PRJNA973112; 95.5% merged) and 14,679 pairs from an AID-expressing CH12F3 replicate "
     "(PRJNA918596; 96.7% merged). Agreeing overlap bases have summed qualities capped at Q60; disagreements retain "
     "the higher-quality call or become N at tied qualities. Adapter read-through is removed, and pairs with <20-bp "
     "overlap are discarded. Overlap disagreements occurred at 0.13% and 0.55% of bases, respectively. Merged reads "
     "enter the Figure 1 workflow using the minimap2 sr "
     "preset, an expected-error limit of 0.02 per base and a minimum alignment requirement of no more than 100 bp. "
     "PRJNA973112 and PRJNA918596 contain uniformly Q30 bases, yielding one quality level in the error model; "
     "PRJNA820976 retains binned qualities (`crm/scicanvas/s12_illumina.py`). Beside it, one PacBio read across a "
     "class-switch junction: a HiFi read carries the Sμ donor, the join and the Sγ1 acceptor on one molecule, "
     "so the microhomology at the join and any inserted bases are read directly rather than inferred. The bar is "
     "the median junction-bearing read over that study’s {pacbio_n} samples, {pacbio_len} bases, and where the break "
     "falls varies between reads, so the donor and acceptor split drawn is schematic. Microhomology classes are "
     "that study’s own unique junctions, averaged over its samples. Reads per unique junction run {pacbio_rpj} "
     "between samples, so junction counts measure long-range PCR recovery and not switching efficiency, and "
     "that deposit has no AID-null arm. Below both rows, coverage of each study locus by a single processed read, "
     "the PacBio switch junctions among them, whose bar is clipped because at scale it would reduce every other "
     "platform to a sliver. Bottom, measurements obtained from the "
     "two workflows in this study and their agreement in the six matched mice. These comparisons are "
     "assay-specific: Sanger amplicons spanning V(D)J junctions can support lineage grouping, UMI consensus represents "
     "one molecule per family, and long-read coverage depends on amplicon design. Nanopore fold separation assigns "
     "zero-call controls one event. Sanger background is the Q≥20 mismatch rate in AID-null clones from Figures S13 and S14."),
    ("Figure S21. SHMtrace analysis workflow, related to STAR Methods",
     [FIG / "FigureS21_overview.png"],
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
     "Sample-sheet rows are illustrative. Figure S22 details the analytical rules and equations."),
    ("Figure S22. Rules underlying read preparation and variant calling, related to STAR Methods",
     [FIG / "FigureS22_engine.png"],
     "The rule and the equation behind each stage of the analysis, each shown working on the J_{H}4 libraries rather "
     "than stated in the abstract. (A) Read accounting in one library: the length and expected-error filters, "
     "competitive alignment and the reads left usable (the individual filters are in Figure S1). (B) Alignment "
     "identity to the target against identity to its paralog; each read is assigned to whichever it matches better. "
     "(C) Reported base quality against the mismatch rate measured in *Aicda*^{−/−} libraries; the dashed line is the "
     "error probability the quality score nominally implies, and the gap between them is what the margin corrects. "
     "(D) Two candidates at equal depth in the same library, with their supporting reads in the wild-type and control "
     "libraries, and why one is called and the other is not. (E) Candidate counts and held-out-decoy counts across "
     "p-value cutoffs, with the threshold that selection returns. (F) The per-animal unique-variant density that "
     "group comparisons are made on. Experimental details: panels are drawn from `intuition_jh4.json`, the "
     "per-library read-accounting records and the bundled example library; the equations shown are those the software "
     "implements."),
    ("Figure S23. End-to-end stress tests using resampled control reads, related to STAR Methods",
     FIG / "FigureS23.png",
     "AID-null J_{H}4 reads were resampled and analyzed with `run_project`, including margin selection, decoy "
     "thresholding and noisy-control exclusion. Cases used KO1. Two-control analyses used one partition each from KO2 "
     "and KO3; three-control analyses used two KO2 partitions and one KO3 partition. Case and control read records "
     "did not overlap. These are conditional resampling experiments using existing libraries, rather than independent "
     "preparations or a new sequencing run. Unless stated otherwise, analyses used three controls and 500 reads. Teal, "
     "exclusion enabled; amber, exclusion disabled, paired on identical inputs. (A) All-null resamples with any call "
     "and 95% Clopper–Pearson intervals; n = {fdv_null_none_on_n} per setting for unmodified background and one "
     "noisy control, and n = {fdv_null_case_only_on_n} "
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
    ("Figure S19. Class-switch junctions in a public PacBio dataset, related to Figure 7",
     FIG / "FigureS19_switch_pacbio.png",
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
    ("Figure S18. Public S\u03bc mutation libraries under PI3K\u03b4 blockade, related to Figure 7",
     FIG / "FigureS18_public_smu.png",
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
    ("Figure S17. Engineered AID diversification, related to Figure 7",
     FIG / "FigureS17_engineered.png",
     "(A) One GFP substrate edited by a panel of deaminases delivered alike (SRP159221), three replicates per "
     "arm: the share of calls at WRC/GYW (circles) and the share that are C\u2192T or G\u2192A (squares), in "
     "wild-type and in UNG/MSH2-null 293T. Bars, group means. The motif share separates AID from APOBEC3A in "
     "both backgrounds; the transition share does so only in wild-type cells, because losing uracil excision "
     "drives transitions whichever enzyme made the uracil. This is the one engineered dataset with replicated "
     "no-editor controls, so its calls are control-calibrated. (B) An AID-diversified antibody V_{H} against "
     "the mutations its authors list for the same clusters[@devilder2019]. Both columns come from that study's "
     "Additional file 1, which publishes, for each selected clone cluster, its nucleotide sequence and the "
     "authors' own list of the mutations they found in it: SHMtrace was given those sequences with the "
     "parental construct as reference, and its calls are compared with that list position by position, each "
     "cluster kept separate and nothing collapsed to a consensus. Each marker is one cluster entry, "
     "its area the number of entries sharing that pair of counts; dashed line, equality. SHMtrace recovers "
     "{eng_dev_mut} of the {eng_dev_mut} listed mutations and calls none the authors do not, over {eng_dev_n} "
     "entries across the study's two selection rounds. Those entries cover {eng_dev_distinct} distinct "
     "clusters, four of which recur between rounds, so the comparison is made against the deposit's own "
     "per-round listing rather than a deduplicated set. This is the only engineered analysis whose comparator "
     "is an independently curated list rather than another run of this software. Analyses in "
     "`crm/liu_mutators.py` and `crm/devilder_clones.py`."),
    ("Figure S4. The control-based false-discovery threshold against planted truth, related to Figure 3",
     FIG / "FigureS4_fdr_calibration.png",
     "The spike-in libraries of Figure 3 (60 substitutions planted into 1–20 reads of each held-out AID-null "
     "control, three draws per control; {fc_libs} libraries per locus) were called at six target false-discovery "
     "rates, every other setting at its default and only the other controls serving as null and decoys. (A) "
     "Realised false-discovery proportion, false over all calls pooled across a locus's libraries, against the "
     "target; dashed line, agreement; shading, over target; a pooled proportion of zero is drawn on the 0 row. "
     "(B) False calls per library. (C) Recall of planted substitutions. Vertical lines, the default target of "
     "5%. At the default, pooled proportions were {fc_def_lo}–{fc_def_hi}% ({fc_def_fp_lo}–{fc_def_fp_hi} false "
     "calls per library), and single libraries reached {fc_def_lib_max}%. Below 5%, false calls levelled off at "
     "{fc_floor_lo}–{fc_floor_hi} per library while recall fell, consistent with library-specific errors that the "
     "other controls do not share and with decoy counts too coarse to resolve small targets (STAR Methods). The "
     "realised proportion depends on the number of true mutations per library (here 60) and is not a per-library "
     "guarantee. Analysis in `crm/spikein.py` (`--tag _fdrcal`) and `crm/figs/figS4_fdr_calibration.py`."),
    ("Figure S5. How many AID-null controls, and which, related to Figure 3",
     FIG / "FigureS5_controls.png",
     "Everything the read-library route does rests on the controls, so their number and their provenance were "
     "varied through the whole `run_project` procedure. J_{H}4: the case is KO1 carrying planted polyclonal "
     "lineages, or none; the controls are halves of KO2 and KO3, so the third and fourth are technical partitions "
     "of the same two animals rather than new mice; {cd_runs} repeats per setting; points, repeats; bars, means. "
     "(A) Recall of planted mutations by number of controls. With one control there are no decoys and the "
     "threshold falls back to Benjamini–Hochberg (BH), which is why recall drops: the empirical threshold needs a "
     "control to hold out. (B) False calls per library in lineage-bearing (filled) and all-null (open) libraries; "
     "above, the pooled false-discovery proportion of the lineage-bearing libraries. No all-null library produced "
     "a call at any number of controls. (C) One control carrying recurrent artifacts at six positions in 6% of its "
     "reads, with the noisy-control rule off and on, beside the clean set. The rule set the noisy control aside in "
     "{cd_noisy_flagged} of {cd_runs} runs and a clean control in {cd_clean_flag_lo}–{cd_clean_flag_hi}; with four "
     "controls and the rule on, {cd_noisy_null_any} of {cd_runs} all-null libraries produced a call. Two controls "
     "are the point at which the procedure works as intended, and a third adds little. Controls from a different "
     "experiment were also tested, on Sμ: each LPS AID-null library, carrying 60 planted substitutions in three "
     "draws, was called against its own experiment's other two AID-null controls, against two or three AID-null "
     "controls from the CH12 experiment of the same BioProject, or against all five. Recall was "
     "{cd_su_rec_matched}% with matched controls and {cd_su_rec_other}% with CH12 controls, at false-discovery "
     "proportions of {cd_su_fdp_matched}% and {cd_su_fdp_other}%; pooling all five lowered false calls from "
     "{cd_su_fp_matched} to {cd_su_fp_pooled} per library. The SRA reads carry no run identifiers, so that "
     "comparison is between experiments rather than between verified sequencing runs. Analysis in "
     "`crm/controls_design.py` and `crm/figs/figS5_controls.py`."),
    ("Figure S3. What a project run writes, and how to read it, related to Figure 1",
     FIG / "FigureS3_outputs_preview.png",
     "Top, a project run: reads per animal and two or more AID-null controls in; a report, figures in PNG, SVG or "
     "PDF, Excel and CSV tables and the full result (with an AIRR rearrangement file) out. Below, SHMtrace's figures "
     "grouped by the question each answers, one per group redrawn from the J_{H}4 nanopore run (three Wt and three "
     "*Aicda*^{−/−} mice, 500 reads each), with a note on how to read it and what to watch for; the other figures of "
     "each group are listed by file name. (A) Read accounting through each filter. (B) The error rate the AID-null "
     "controls showed at each position; dashed line, median. (C) Unique-variant density per animal; open circles, "
     "animals with no call, drawn at one event; bars, geometric means; p, Welch test on log density. (D) The twelve "
     "substitution types; bars, Wt pooled; points, mice. (E) Share of reads mutated at each position, Wt above the "
     "line and *Aicda*^{−/−} below; ticks, WRC/GYW hotspots. (F) Clones defined by the V(D)J junction, by size and "
     "number of mutations. Seven figures are drawn for every project and up to nineteen more when the run has the "
     "data, so a run of this project writes {s3_drawn} of them. Each panel is an output of the run itself, not a schematic. Built with `crm/scicanvas/figS3_outputs.js`."),
    ("Figure S6. Real hypermutation diluted into AID-null reads, related to Figure 3",
     FIG / "FigureS6_titration.png",
     "Each AID-active library (donor) was mixed with one AID-null library (diluent) at 1–100% donor reads in 500 "
     "usable reads, {ti_reps} draws per fraction, and called against the other two AID-null libraries (J_{H}4: "
     "donors Wt2–Wt4, diluent KO1, controls KO2 and KO3; N114: donors C1, D3 and D5, diluent D1, controls C3 and C4). "
     "(A, B) Mean observed mutation frequency against the frequency expected from donor and diluent alone; lines, "
     "the range of draws; dashed line, equality; a frequency of zero is drawn on the 0 row. N114 used "
     "`ancestral_reference` (STAR Methods). (C) Unique mutations called, mean of draws, against the share of donor "
     "reads. On N114, frequency was {ti_n114_ratio} of the expected value at 25–50% donor reads and calls began at "
     "{ti_n114_det}%. On J_{H}4, calling was all-or-nothing and began at {ti_jh4_det_lo}–{ti_jh4_det_hi}% donor "
     "reads depending on the donor, consistent with a control-based threshold that admits no call until enough "
     "candidates exceed the controls' background. Analysis in `crm/titration.py` and `crm/figs/figS6_titration.py`."),
    ("Figure S11. SHMtrace and the IgBLAST–SHazaM route on the same N114 reads, related to Figure 4",
     FIG / "FigureS11_immcantation.png",
     "The 500 reads SHMtrace analyzed from each N114 library were assigned by IgBLAST ({im_vcall}; the parental "
     "line is {im_vid}% identical to it over {im_vlen} bp) and their substitutions over the V segment counted by "
     "SHazaM, with the IMGT germline, with the parental V segment as germline, and with the parental germline after "
     "masking bases below Q20; SHMtrace's calls over the same segment are shown for comparison. (A) Mutation "
     "frequency per library; filled, pMX-*Aicda*; open, empty vector; triangles, libraries with no SHMtrace call, "
     "drawn at one event; bars, geometric means. (B) Geometric-mean separation of pMX-*Aicda* from empty-vector "
     "libraries. The SHazaM routes count every substitution in every read, including real single-read mutations "
     "that SHMtrace does not call, so their AID-library values are higher ({im_aid_q20} against {im_aid_sh} per "
     "10^{3} bp with masking); their empty-vector values ({im_ctl_imgt} and {im_ctl_q20} per 10^{3} bp without and "
     "with the parental germline and masking) are the background they cannot remove. MiXCR was not run. Analysis in "
     "`crm/immcantation_compare.py` and `crm/figs/figS11_immcantation.py`."),
]

if __name__ == "__main__":
    main(sys.argv[1])
