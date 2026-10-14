"""Ver22: the JIM manuscript on SHMtrace v0.2's control-based FDR and the second C3 library, from Ver21.

    python revise_manuscript_v22.py VER21.docx VER22.docx      (needs python-docx)

Ver21 called variants with Benjamini-Hochberg at FDR 0.05 and used the first library of N114
empty-vector subclone C3. Ver22 reports the analyses as rerun on 13 October 2026:
  - the FDR is estimated from the controls themselves (held-out controls as decoys), as in
    SHMtrace v0.2's recommended mode and the Cell Reports Methods article;
  - N114 subclone C3 is its second library (75K3B3_2_C3_Empty_March_2024_2), the first being
    noisier than the other controls.
Every number below is computed from the analysis outputs (docs/calibration/figures, the CRM
numbers, the Sanger calls), never typed in. Ver21's highlights are cleared first, so the yellow
text in Ver22 is what changed since Ver21. SOURCE is never modified, and the script stops if any
expected text is not found.
"""
import copy
import json
import math
import random
import statistics
import sys
from pathlib import Path

import numpy as np
from docx import Document
from docx.enum.text import WD_COLOR_INDEX
from docx.oxml.ns import qn
from scipy import stats

from revise_manuscript_sanger import set_cell, swap_image
from revise_manuscript_v18 import replace_text, restyle

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
CAL = HERE.parent / "figures"
SV = CAL / "sanger_vs_nanopore/recommended_manuscript"
S1F = HERE / "figure1_scicanvas/out"
sys.path.insert(0, str(HERE.parent / "crm/manuscript"))
import crm_numbers  # noqa: E402

ELIAS = ("Elias, J.E., Gygi, S.P., 2007. Target-decoy search strategy for increased confidence in large-scale "
         "protein identifications by mass spectrometry. Nat. Methods 4, 207–214. https://doi.org/10.1038/nmeth1019")


def geo(x):
    return 10 ** np.mean(np.log10(x))


def boot(x, rng, n=20000):
    v = sorted(geo([rng.choice(x) for _ in x]) for _ in range(n))
    return v[int(.025 * n)], v[int(.975 * n) - 1]


def a_n(n):
    """'a' or 'an' before a number read aloud (an 8-, 11-, 18-, 80-...)."""
    s = str(n)
    return "an" if s.startswith("8") or s in ("11", "18") or (len(s) in (5, 6) and s[:2] in ("11", "18")) else "a"


def word(n):
    return ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"][n] \
        if n <= 10 else f"{n:,}"


def values():
    V = dict(crm_numbers.compute())
    jh4 = {s["sample"]: s for s in json.load(open(CAL / "jh4_shmtrace/project.json"))["samples"]}
    n114 = {s["sample"]: s for s in json.load(open(CAL / "n114_shmtrace/project.json"))["samples"]}
    rng = random.Random(1)
    dens = lambda S, k: S[k]["unique_density_x1e5"] or S[k]["one_event_density_x1e5"]   # noqa: E731
    owt, oko = [dens(jh4, k) for k in ("WT2", "WT3", "WT4")], [dens(jh4, k) for k in ("KO1", "KO2", "KO3")]
    pa, pe = [dens(n114, k) for k in ("C1", "D3", "D5")], [dens(n114, k) for k in ("C3", "C4", "D1")]
    V["ont_wt_ci"] = "{:.0f}–{:.0f}".format(*boot(owt, rng))
    V["ont_ko_zero"] = sum(1 for k in ("KO1", "KO2", "KO3") if not jh4[k]["unique_mutations"])
    V["ko_calls"] = [jh4[k]["unique_mutations"] for k in ("KO1", "KO2", "KO3")]
    V["wt_calls"] = [jh4[k]["unique_mutations"] for k in ("WT2", "WT3", "WT4")]
    V["n114_a"] = [n114[k]["unique_mutations"] for k in ("C1", "D3", "D5")]
    V["n114_e"] = [n114[k]["unique_mutations"] for k in ("C3", "C4", "D1")]
    V["n114_e_zero"] = sum(1 for c in V["n114_e"] if not c)
    V["n114_geo"] = f"{geo(pa):.0f}"
    V["n114_ctl_geo"] = f"{geo(pe):.2f}"
    V["n114_ci"] = "{:.0f}–{:.0f}".format(*boot(pa, rng))
    V["n114_fold"] = f"{geo(pa) / geo(pe):.0f}"
    V["n114_p"] = f"{stats.ttest_ind(np.log10(pa), np.log10(pe), equal_var=False).pvalue:.4f}"
    gc = sum(n114[k]["pct_GC"] * n114[k]["unique_mutations"] / 100 for k in ("C1", "D3", "D5"))
    V["n114_gc"] = f"{100 * gc / sum(V['n114_a']):.0f}"
    # the pooled background and its Q
    bases = 500 * 529 * 3
    ko_tot = sum(V["ko_calls"])
    V["ko_tot"], V["ko_bases"] = ko_tot, bases
    V["ko_one_per"] = f"{round(bases / ko_tot, -5):,.0f}"
    V["q_bg"] = f"{-10 * math.log10(float(V['ont_ko_geo']) * 1e-5):.0f}"
    # read accounting (Supplementary Table S1)
    acc = {run: json.load(open(CAL / f"{run}_shmtrace/table_s1.json")) for run in ("jh4", "n114")}
    every = [v for t in acc.values() for v in t.values()]
    usable = [v["usable"] for v in every]
    V["delivered_lo"] = f"{min(v['delivered'] for v in every):,}"
    V["delivered_hi"] = f"{max(v['delivered'] for v in every):,}"
    V["usable_mean"] = f"{statistics.mean(usable):,.0f}"
    V["usable_mean_n"] = statistics.mean(usable)
    V["usable_median"] = f"{math.floor(statistics.median(usable) + 0.5):,}"
    for run, t in acc.items():
        D = sum(v["delivered"] for v in t.values())
        part = lambda *ks: round(100 * sum(v[k] for v in t.values() for k in ks) / D)   # noqa: E731
        V[f"{run}_rm"] = (part("too_long"), part("low_quality"),
                          part("not_target", "chimeric", "short_alignment", "paralogue_alleles", "minor_haplotype"),
                          part("usable"))
    # cost per QC-passing read (Table 1): US$35 per library over the mean usable reads
    per_read = 35 / V["usable_mean_n"]
    V["cost_read"] = f"{per_read:.3f}"
    V["cost_fold"] = f"{(315 / 37) / per_read:.0f}"
    # depth series
    for a, cases in (("jh4", ("WT2", "WT3", "WT4")), ("n114", ("C1", "D3", "D5"))):
        d = json.load(open(CAL / f"depth_shmtrace/{a}.json"))
        ctl = [d[k][s]["calls"] for k in d for s in d[k] if s not in cases]
        V[f"{a}_ctl_range"] = f"{min(ctl)}–{max(ctl)}"
        V[f"{a}_dgeo"] = {k: geo([d[k][s]["density"] for s in cases]) for k in d}
        V[f"{a}_depth"] = d
    d = V["jh4_depth"]
    V["wt3_250"], V["wt3_max"] = d["250"]["WT3"]["calls"], max(d[k]["WT3"]["calls"] for k in d)
    V["wt2_500"], V["wt2_1000"] = d["500"]["WT2"]["calls"], d["1000"]["WT2"]["calls"]
    # read-weighted N114 observations (Supplementary Figure S5)
    sr = {s["sample"]: s for s in json.load(open(CAL / "n114_shmtrace/project.json"))["samples"]}
    V["n114_rw_a"] = f"{sum(sr[k]['mutations_total'] for k in ('C1', 'D3', 'D5')):,}"
    V["n114_rw_e"] = f"{sum(sr[k]['mutations_total'] for k in ('C3', 'C4', 'D1')):,}"
    # Sanger-ONT spectrum statistics not in the CRM numbers
    st = json.load(open(SV / "stats.json"))["spectrum"]
    V["tstv"] = tuple(f"{x:.2f}" for x in st["gc_ts_tv"]["sidak_p"])
    return V


def three(xs):
    return f"{xs[0]}, {xs[1]} and {xs[2]}"


def edits(V):
    nz = V["ont_ko_zero"]
    ko_zero = {1: "the one Aicda−/− library without calls placed at its", 2: "the two Aicda−/− libraries without calls "
               "placed at their", 3: "the three Aicda−/− libraries without calls placed at their"}[nz]
    e_tot = sum(V["n114_e"])
    e_phrase = "no calls" if e_tot == 0 else ("a single call" if e_tot == 1 else f"{e_tot} calls")
    nf, of = V["n114_fold"], V["ont_fold"]
    jd, nd = V["jh4_depth"], V["n114_depth"]
    pct500 = {a: round(100 * statistics.mean(d["500"][s]["calls"] / max(d[k][s]["calls"] for k in d) for s in cases))
              for a, d, cases in (("jh4", jd, ("WT2", "WT3", "WT4")), ("n114", nd, ("C1", "D3", "D5")))}
    lo, hi = sorted(pct500.values())
    jg, ng = V["jh4_dgeo"], V["n114_dgeo"]
    rm, rn = V["jh4_rm"], V["n114_rm"]
    return [
        # ---- Abstract
        ("Somatic hypermutation (SHM) diversifies", "(Pearson r = 0.93)", f"(Pearson r = {V['r_spec']})"),
        ("Somatic hypermutation (SHM) diversifies",
         "ONT detected at least a 92-fold group difference, with no calls in any empty-vector control,",
         f"ONT detected at least {a_n(nf)} {nf}-fold group difference, with {e_phrase} among the three empty-vector "
         f"controls,"),
        ("Somatic hypermutation (SHM) diversifies", "a mean of 91–94%", f"a mean of {lo}–{hi}%"),
        ("Somatic hypermutation (SHM) diversifies", "AID-deficient controls", "AID-null controls"),
        # ---- Methods 2.7: read accounting
        ("Each barcoded library yielded", "2,336–5,000 delivered reads",
         f"{V['delivered_lo']}–{V['delivered_hi']} delivered reads"),
        ("Each barcoded library yielded", "(mean 1,620; median 1,456)",
         f"(mean {V['usable_mean']}; median {V['usable_median']})"),
        ("Each barcoded library yielded", "28% at alignment", f"{rm[2]}% at alignment"),
        ("Each barcoded library yielded", "leaving 33% usable", f"leaving {rm[3]}% usable"),
        ("Each barcoded library yielded", "25%, 22%, 9% and 43%", f"{rn[0]}%, {rn[1]}%, {rn[2]}% and {rn[3]}%"),
        # ---- Methods 2.8: the control-based FDR
        ("Each candidate was tested against the sequencing background",
         "This factor was 5 for both assays; at face value, 55 of 1,391 such JH4 tests and 69 of 1,225 N114 tests "
         "reached p < 0.01. Benjamini–Hochberg correction was applied across all positions and alternate alleles at "
         "a false-discovery rate of 0.05, and each control library was itself called against a background fitted "
         "from the other two, never from itself.",
         f"This factor was {V['loco_jh4_margin']} for JH4 and {V['loco_n114_margin']} for N114; at face value, "
         f"{V['loco_jh4_face_n']} such JH4 tests and {V['loco_n114_face_n']} N114 tests reached p < 0.01. The same "
         "held-out tests set the significance threshold. Each held-out control served as a decoy, and for each sample "
         "SHMtrace kept the largest p-value threshold at which the number of decoy candidates expected to pass, "
         "scaled to the sample's tested bases, was at most 5% of the sample's passing candidates, a target–decoy "
         "estimate of the false-discovery rate (Elias and Gygi, 2007). Each control library was itself called "
         "against a background fitted from the other two, never from itself, with only the other controls as "
         "decoys."),
        ("For the separate full-length ONT analysis", "(leave-one-control-out factor 5)",
         f"(leave-one-control-out factor {V['loco_full_margin']})"),
        ("Depth-standardised unique-variant densities were log10-transformed",
         "(one Aicda−/− ONT library and all three empty-vector subclones)",
         f"({word(nz)} Aicda−/− ONT {'library' if nz == 1 else 'libraries'} and {word(V['n114_e_zero'])} "
         f"empty-vector subclones)"),
        # ---- Results 3.1
        ("In these samples, Wt mice had markedly higher",
         "p = 0.0025 for ONT (≥ 137-fold, the one Aicda−/− library without calls placed at its one-event density)",
         f"p = {V['ont_p']} for ONT (≥ {of}-fold, {ko_zero} one-event density)"),
        ("In these samples, Wt mice had markedly higher", "~66 × 10−5 by ONT", f"~{V['ont_wt_geo']} × 10−5 by ONT"),
        ("After all filters, the three Aicda−/− libraries retained", "0, 2 and 1 calls", f"{three(V['ko_calls'])} calls"),
        ("After all filters, the three Aicda−/− libraries retained",
         "(pooled, 3 calls in ~790,000 analysed read-bases, about one per 260,000)",
         f"(pooled, {V['ko_tot']} calls in ~{round(V['ko_bases'], -4):,} analysed read-bases, about one per "
         f"{V['ko_one_per']})"),
        ("After all filters, the three Aicda−/− libraries retained", "a Q53-equivalent", f"a Q{V['q_bg']}-equivalent"),
        ("After all filters, the three Aicda−/− libraries retained", "205, 340 and 73 retained",
         f"{three(V['wt_calls'])} retained"),
        # ---- 3.2
        ("Because ONT reads span the complete amplicon", "(39 versus 0.51 × 10−5; ~78-fold)",
         f"({V['full_wt']} versus {V['full_ko']} × 10−5; ~{V['full_fold']}-fold)"),
        ("Because ONT reads span the complete amplicon", "Approximately 70% of Wt calls",
         f"Approximately {V['full_450']}% of Wt calls"),
        ("Because ONT reads span the complete amplicon", "about eight-fold",
         f"about {word(int(V['full_decline']))}-fold"),
        ("Because ONT reads span the complete amplicon", "(Pearson r = 0.986; all differences ≤ 1.8 percentage points)",
         f"(Pearson r = {V['s1_r']}; all differences ≤ {V['s1_maxdiff']} percentage points)"),
        ("Because ONT reads span the complete amplicon", "(795 versus 618;", f"({V['s1_nfull']} versus {V['s1_n529']};"),
        ("Because each read spans the amplicon, the reads of one B-cell clone", "(7 of 1,176)",
         f"({V['tree_incompat']})"),
        # ---- 3.3
        ("The matched ONT and Sanger JH4 datasets produced", "Pearson r = 0.93 and slope = 1.00",
         f"Pearson r = {V['r_spec']} and slope = {V['slope']}"),
        ("The matched ONT and Sanger JH4 datasets produced", "interaction (p = 0.61)",
         f"interaction (p = {V['interaction_p']})"),
        ("The matched ONT and Sanger JH4 datasets produced", "(smallest adjusted p = 0.68)",
         f"(smallest adjusted p = {V['sidak_min']})"),
        ("The matched ONT and Sanger JH4 datasets produced", "(25-bp bins; Pearson r = 0.92)",
         f"(25-bp bins; Pearson r = {V['r_pos']})"),
        # ---- 3.4 N114
        ("At matched 500-read depth, pMX-Aicda subclones",
         "30 × 10−5, whereas none of the three empty-vector controls had a call (Figure 4A). With the controls placed "
         "at their one-event density (0.33 × 10−5), the group difference was at least 92-fold (exploratory Welch "
         "p = 0.0048).",
         f"{V['n114_geo']} × 10−5, whereas the three empty-vector controls had {e_phrase} between them (Figure 4A). "
         f"With the controls without calls placed at their one-event density ({V['n114_ctl_geo']} × 10−5, the "
         f"density of one call), the group difference was at least {nf}-fold (exploratory Welch p = {V['n114_p']})."),
        ("At matched 500-read depth, pMX-Aicda subclones", "(96% of calls at G:C)", f"({V['n114_gc']}% of calls at G:C)"),
        ("At matched 500-read depth, pMX-Aicda subclones", "whereas the controls had no calls",
         f"whereas the controls had {e_phrase} between them"),
        # ---- 3.5 depth
        ("We analysed each library at 30, 60, 125, 250, 500 and 1,000 reads", "were 94% (JH4) and 91% (N114)",
         f"were {pct500['jh4']}% (JH4) and {pct500['n114']}% (N114)"),
        ("We analysed each library at 30, 60, 125, 250, 500 and 1,000 reads",
         "0–3 per JH4 Aicda−/− library, and 0–5 per N114 empty-vector subclone at 30 and 60 reads and none from 125 "
         "reads upward",
         f"{V['jh4_ctl_range']} per JH4 Aicda−/− library and {V['n114_ctl_range']} per N114 empty-vector subclone"),
        ("Depth-standardised density therefore declined with depth",
         "~460 × 10−5 at 30 reads to 66 × 10−5 at 500 and 36 × 10−5 at 1,000 reads, and the pMX-Aicda mean from ~110 "
         "to 30 and 16 × 10−5",
         f"~{round(jg['30'], -1):.0f} × 10−5 at 30 reads to {jg['500']:.0f} × 10−5 at 500 and {jg['1000']:.0f} × 10−5 "
         f"at 1,000 reads, and the pMX-Aicda mean from ~{ng['30']:.0f} to {ng['500']:.0f} and {ng['1000']:.0f} × 10−5"),
        ("Depth-standardised density therefore declined with depth", "(332 of a maximum of 340)",
         f"({V['wt3_250']} of a maximum of {V['wt3_max']})"),
        ("Depth-standardised density therefore declined with depth", "(205 and 247)",
         f"({V['wt2_500']} and {V['wt2_1000']})"),
        # ---- Discussion: cost per QC-passing read
        ("The workflow also scales favourably", "(a mean of 1,620 usable reads per sample), the corresponding value is "
         "~US$0.022 per read, an approximately 394-fold difference. Thus, the 394-fold estimate",
         f"(a mean of {V['usable_mean']} usable reads per sample), the corresponding value is ~US${V['cost_read']} per "
         f"read, an approximately {V['cost_fold']}-fold difference. Thus, the {V['cost_fold']}-fold estimate"),
        # ---- Conclusions
        ("This proof-of-concept study shows", "(Pearson r = 0.93)", f"(Pearson r = {V['r_spec']})"),
        ("This proof-of-concept study shows",
         "at least a 92-fold difference between AID-expressing and empty-vector N114 subclones, with no calls in any "
         "of the three empty-vector controls",
         f"at least {a_n(nf)} {nf}-fold difference between AID-expressing and empty-vector N114 subclones, with "
         f"{e_phrase} among the three empty-vector controls"),
        # ---- Table 1 legend
        ("Table 1. Illustrative per-sample comparison", "mean of 1,620 usable reads",
         f"mean of {V['usable_mean']} usable reads"),
        # ---- Supplementary Table S2 note
        ("Group geometric means with descriptive", "Wt ONT 66 [28–132]",
         f"Wt ONT {V['ont_wt_geo']} [{V['ont_wt_ci']}]"),
        ("Group geometric means with descriptive", "(one library with no calls, at its one-event density)",
         f"({word(nz)} {'library' if nz == 1 else 'libraries'} with no calls, at the one-event density)"),
        ("Group geometric means with descriptive",
         "empty vector 0.33 (no subclone had a call; all three at their one-event density); pMX-Aicda 30 [20–55]",
         f"empty vector {V['n114_ctl_geo']} ({three(V['n114_e'])} calls; a subclone without calls at the one-event "
         f"density, which equals one call); pMX-Aicda {V['n114_geo']} [{V['n114_ci']}]"),
    ]


def legends(V):
    nz, nf = V["ont_ko_zero"], V["n114_fold"]
    e_tot = sum(V["n114_e"])
    jd, nd = V["jh4_depth"], V["n114_depth"]
    pct500 = {a: round(100 * statistics.mean(d["500"][s]["calls"] / max(d[k][s]["calls"] for k in d) for s in cases))
              for a, d, cases in (("jh4", jd, ("WT2", "WT3", "WT4")), ("n114", nd, ("C1", "D3", "D5")))}
    return [
        ("Figure 1. Study overview", "(Benjamini–Hochberg false-discovery rate ≤ 0.05)",
         "(false-discovery rate ≤ 0.05, estimated with the controls as decoys)"),
        ("Figure 2. Depth-standardised",
         "One Aicda−/− ONT library had none and is drawn (open circle) and tested at its one-event density",
         f"{word(nz).capitalize()} Aicda−/− ONT libraries had none and are drawn (open circles) and tested at their "
         "one-event density" if nz > 1 else
         "One Aicda−/− ONT library had none and is drawn (open circle) and tested at its one-event density"),
        ("Figure 2. Depth-standardised", "≥ 137-fold (ONT; p = 0.0025)", f"≥ {V['ont_fold']}-fold (ONT; p = {V['ont_p']})"),
        ("Figure 2. Depth-standardised", "Wt ONT 66 [28–132]", f"Wt ONT {V['ont_wt_geo']} [{V['ont_wt_ci']}]"),
        ("Figure 3. ONT yields", "(n = 471 Sanger and 618 ONT", f"(n = {V['n_sanger']} Sanger and {V['n_ont']} ONT"),
        ("Figure 3. ONT yields", "interaction was p = 0.61", f"interaction was p = {V['interaction_p']}"),
        ("Figure 3. ONT yields", "(Pearson r = 0.93; slope = 1.00;", f"(Pearson r = {V['r_spec']}; slope = {V['slope']};"),
        ("Figure 3. ONT yields", "adjusted paired p = 0.15 and 0.98", f"adjusted paired p = {V['tstv'][0]} and {V['tstv'][1]}"),
        ("Figure 4. Filtered apparent variants",
         "None of the three empty-vector subclones had a call; they are drawn (open circles) and tested at their "
         "one-event density (0.33 × 10−5), so the group difference is a lower bound: the exploratory Welch test on "
         "log10-transformed densities gave p = 0.0048 and at least a 92-fold difference (30 versus 0.33 × 10−5).",
         f"The three empty-vector subclones had {three(V['n114_e'])} calls; those without a call are drawn (open "
         f"circles) and tested at their one-event density ({V['n114_ctl_geo']} × 10−5, the density of one call), so "
         f"the group difference is a lower bound: the exploratory Welch test on log10-transformed densities gave "
         f"p = {V['n114_p']} and at least {a_n(nf)} {nf}-fold difference ({V['n114_geo']} versus "
         f"{V['n114_ctl_geo']} × 10−5)."),
        ("Figure 4. Filtered apparent variants", "for pMX-Aicda (n = 300) and empty vector (n = 0)",
         f"for pMX-Aicda (n = {sum(V['n114_a'])}) and empty vector (n = {e_tot})"),
        ("Supplementary Figure S1. Pooled substitution", "(n = 618 distinct calls)", f"(n = {V['s1_n529']} distinct calls)"),
        ("Supplementary Figure S1. Pooled substitution", "n = 795)", f"n = {V['s1_nfull']})"),
        ("Supplementary Figure S1. Pooled substitution", "(Pearson r = 0.986; all differences ≤ 1.8 percentage points)",
         f"(Pearson r = {V['s1_r']}; all differences ≤ {V['s1_maxdiff']} percentage points)"),
        ("Supplementary Figure S3. Individual JH4", "(n = 205, 340, and 73 unique variants",
         f"(n = {V['wt_calls'][0]}, {V['wt_calls'][1]}, and {V['wt_calls'][2]} unique variants"),
        ("Supplementary Figure S3. Individual JH4", "(n = 0, 2, and 1)",
         f"(n = {V['ko_calls'][0]}, {V['ko_calls'][1]}, and {V['ko_calls'][2]})"),
        ("Supplementary Figure S4. Individual N114", "(n = 75, 59, and 166 unique variants",
         f"(n = {V['n114_a'][0]}, {V['n114_a'][1]}, and {V['n114_a'][2]} unique variants"),
        ("Supplementary Figure S4. Individual N114", "(n = 0, 0, and 0)",
         f"(n = {V['n114_e'][0]}, {V['n114_e'][1]}, and {V['n114_e'][2]})"),
        ("Supplementary Figure S5. Recurrence-weighted", "(n = 1,478 read-weighted observations)",
         f"(n = {V['n114_rw_a']} read-weighted observations)"),
        ("Supplementary Figure S5. Recurrence-weighted", "(n = 0; no calls)",
         f"(n = {V['n114_rw_e']}; {e_tot} {'call' if e_tot == 1 else 'calls'})"),
        ("Supplementary Figure S6. ",
         "JH4 controls retained 0–3 calls at every depth, and N114 controls 0–5 calls at 30 and 60 reads and none "
         "from 125 reads upward",
         f"JH4 controls retained {V['jh4_ctl_range']} calls and N114 controls {V['n114_ctl_range']} calls at every depth"),
        ("Supplementary Figure S7. ", "approached the 1,000-read count by 500 reads in this dataset",
         "relative to the 1,000-read count in this dataset"),
        ("Supplementary Figure S7. ", "were 94% (JH4) and 91% (N114)",
         f"were {pct500['jh4']}% (JH4) and {pct500['n114']}% (N114)"),
        ("Supplementary Figure S9. Sanger and ONT agree", "(Pearson r = 0.92)", f"(Pearson r = {V['r_pos']})"),
        ("Supplementary Figure S10. Mutation tree", "7 of 1,176 mutation pairs",
         f"{V['tree_incompat']} mutation pairs"),
        ("Supplementary Figure S11. ", "multiplied by a five-fold margin chosen by leave-one-control-out",
         f"multiplied by a margin ({word(int(V['loco_jh4_margin']))}-fold for JH4, "
         f"{word(int(V['loco_n114_margin']))}-fold for N114) chosen by leave-one-control-out"),
        ("Supplementary Figure S11. ", "were kept at a Benjamini–Hochberg false-discovery rate of 0.05",
         "were kept at a false-discovery rate of 0.05 estimated with the controls as decoys"),
    ]


def tables(V):
    nz = V["ont_ko_zero"]
    return [  # (table index, old cell text, new)
        (0, "~US$0.022 per QC-passing read available", f"~US${V['cost_read']} per QC-passing read available"),
        (2, "p = 0.0025 (≥ 137-fold; the control with no calls at one event)",
         f"p = {V['ont_p']} (≥ {V['ont_fold']}-fold; " +
         ("the control" if nz == 1 else f"the {word(nz)} controls") + " with no calls at one event)"),
        (2, "p = 0.0048 (≥ 92-fold; all three controls without calls, at one event)",
         f"p = {V['n114_p']} (≥ {V['n114_fold']}-fold; controls without calls at one event)"),
        (2, "r = 0.93 [0.77–0.98]; slope = 1.00", f"r = {V['r_spec']} [{V['r_spec_ci']}]; slope = {V['slope']}"),
        (2, "interaction p = 0.61; none significant (minimum adjusted p = 0.68)",
         f"interaction p = {V['interaction_p']}; none significant (minimum adjusted p = {V['sidak_min']})"),
        (2, "Ts adjusted p = 0.15; Tv adjusted p = 0.98", f"Ts adjusted p = {V['tstv'][0]}; Tv adjusted p = {V['tstv'][1]}"),
        (2, "Q53-equivalent", f"Q{V['q_bg']}-equivalent"),
    ]


IMAGES = [
    ("Figure 1. ", S1F / "Figure1_overview_preview.png"),
    ("Figure 2. Depth-standardised", FIG / "Figure2_JH4_combined.png"),
    ("Figure 3. ONT yields", FIG / "Figure3_concordance.png"),
    ("Figure 4. Filtered apparent variants", FIG / "Figure4_hybridoma_combined.png"),
    ("Supplementary Figure S1. Pooled substitution", FIG / "FigureS_spectrum_529_vs_fulllength.png"),
    ("Supplementary Figure S3. Individual JH4", FIG / "FigureS6_individual_JH4_ONT_matrices.png"),
    ("Supplementary Figure S4. Individual N114", FIG / "FigureS7_individual_N114_matrices.png"),
    ("Supplementary Figure S5. Recurrence-weighted", FIG / "FigureS7_recurrence_weighted.png"),
    ("Supplementary Figure S6. ", FIG / "FigureS3_depth_discovery_density.png"),
    ("Supplementary Figure S7. ", FIG / "FigureS4_saturation.png"),
    ("Supplementary Figure S9. Sanger and ONT agree", SV / "fig3_agreement_position_hotspot.png"),
    ("Supplementary Figure S10. Mutation tree", CAL / "mutation_trees/mutation_trees_wt3.png"),
    ("Supplementary Figure S11. ", S1F / "Figure_workflow_detailed_preview.png"),
]

# Supplementary Table S1 ONT rows: table key -> (run, library, new key)
S1_ROWS = {"Wt2": ("jh4", "WT2"), "Wt3": ("jh4", "WT3"), "Wt4": ("jh4", "WT4"), "KO1": ("jh4", "KO1"),
           "KO2": ("jh4", "KO2"), "KO3": ("jh4", "KO3"), "1_C3": ("n114", "C3"), "3_C4": ("n114", "C4"),
           "4_D1": ("n114", "D1"), "6_C1": ("n114", "C1"), "7_D3": ("n114", "D3"), "10_D5": ("n114", "D5")}
S1_RENAME = {"1_C3": "2_C3"}     # the second library of subclone C3


def s1_values():
    out = {}
    for key, (run, lib) in S1_ROWS.items():
        acct = json.load(open(CAL / f"{run}_shmtrace/table_s1.json"))[lib]
        proj = {s["sample"]: s for s in json.load(open(CAL / f"{run}_shmtrace/project.json"))["samples"]}[lib]
        d = proj["unique_density_x1e5"] or 0.0
        dens = f"{d:.1f}" if d >= 10 else f"{d:.2f}"
        out[key] = (f"{acct['delivered']:,}", f"{acct['usable']:,}", f"{acct['analysed']}",
                    f"{acct['candidates_tested']}", f"{acct['calls']}", dens)
    return out


def clear_highlights(doc):
    def runs():
        for p in doc.paragraphs:
            yield from p.runs
        for t in doc.tables:
            for row in t.rows:
                for c in row.cells:
                    for p in c.paragraphs:
                        yield from p.runs
    for r in runs():
        if r.font.highlight_color is not None:
            r.font.highlight_color = None


def add_reference(doc, after_start, text):
    """A new reference paragraph after the one starting `after_start`, in its style, highlighted."""
    ref = next(p for p in doc.paragraphs if p.text.startswith(after_start))
    new = copy.deepcopy(ref._p)
    ref._p.addnext(new)
    from docx.text.paragraph import Paragraph
    par = Paragraph(new, ref._parent)
    for r in par.runs[1:]:
        r._r.getparent().remove(r._r)
    for h in new.iter(qn("w:hyperlink")):
        h.getparent().remove(h)
    par.runs[0].text = text
    par.runs[0].font.highlight_color = WD_COLOR_INDEX.YELLOW


def main(src, out):
    V = values()
    doc = Document(src)
    clear_highlights(doc)

    def find(start):
        hits = [p for p in doc.paragraphs if p.text.startswith(start)]
        if len(hits) != 1:
            raise RuntimeError(f"{len(hits)} paragraphs start with {start!r}")
        return hits[0]

    work = [(find(s), old, new) for s, old, new in edits(V) + legends(V)]
    for p, old, new in work:
        if old != new:
            replace_text(p, old, new)
    add_reference(doc, "Edgar, R.C.", ELIAS)

    for ti, old, new in tables(V):
        cells = [c for row in doc.tables[ti].rows for c in row.cells if old in c.text]
        if not cells:
            raise RuntimeError(f"table {ti}: not found {old!r}")
        if old != new:
            set_cell(cells[0], old, new)
    s1 = doc.tables[1]
    vals, seen = s1_values(), set()
    for row in s1.rows:
        key = row.cells[0].text.strip()
        if key in vals and row.cells[2].text.strip() == "ONT":
            seen.add(key)
            if key in S1_RENAME:
                set_cell(row.cells[0], key, S1_RENAME[key])
            for col, val in zip((3, 4, 5, 7, 8, 9), vals[key], strict=True):
                cell = row.cells[col]
                if cell.text.strip() != val:
                    set_cell(cell, cell.text.strip(), val)
    if seen != set(vals):
        raise RuntimeError(f"S1 rows missing: {set(vals) - seen}")

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
    doc.core_properties.comments = ("Ver22: control-based FDR (SHMtrace v0.2) and the second library of N114 "
                                    "subclone C3 (shmtrace docs/calibration/manuscript_corrected/revise_manuscript_v22.py)")
    doc.save(out)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
