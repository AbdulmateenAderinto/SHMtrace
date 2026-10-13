"""Text of the Cell Reports Methods article (SHMtrace). Markup: see docbuild.py.

Numbers for the public Illumina reanalyses are filled from out/public/*/project.json by
build_crm.py (keys in braces, e.g. {ldha_fold}); everything else is written out and traced to
the analysis outputs listed in STAR Methods.
"""

TITLE = ("Control-calibrated profiling of somatic hypermutation from cloning-free long reads, Sanger traces "
         "and short reads")
RUNNING = "Control-calibrated SHM profiling"
AUTHORS = "Abdulmateen Aderinto^{1,3} and Alberto Martin^{1,2,*}"
AFFIL = ["^{1}Department of Immunology, University of Toronto, Toronto, ON M5S 1A8, Canada",
         "^{2}Lead contact", "^{3}Technical contact",
         "*Correspondence: alberto.martin@utoronto.ca"]

SUMMARY = (
    "Somatic hypermutation (SHM) is usually measured by cloning amplicons for Sanger sequencing, because a single "
    "long read errs at rates close to the mutation rates being measured. We describe SHMtrace, which tests every "
    "candidate substitution against the errors measured in matched activation-induced cytidine deaminase (AID)-null "
    "libraries from the same sequencing run and compares groups one animal at a time. In planted mutations, simulated "
    "B cell lineages and libraries from six loci, five laboratories and three sequencing platforms, SHMtrace kept the "
    "AID-null background at {ctl_rec} calls per assay while retaining wild-type signal, separating AID-active from "
    "AID-null libraries {sep_lo}- to {sep_hi}-fold in four of five nanopore assays, where a quality-based caller "
    "reached {lf_lo}- to {lf_hi}-fold and an uncorrected test {gen_lo}- to {gen_hi}-fold. Nanopore and Sanger sequencing of the same mice gave concordant spectra. Full-length reads "
    "added whole-amplicon maps, B cell clones and mutation trees.")
KEYWORDS = ("somatic hypermutation, activation-induced cytidine deaminase, nanopore sequencing, variant calling, "
            "negative controls, false-discovery rate, pseudoreplication, immunoglobulin, mutation trees, Sanger "
            "sequencing")
MOTIVATION = (
    "Measuring somatic hypermutation still relies on cloning PCR products and sequencing colonies one at a time, "
    "which is slow and samples few molecules. Long reads could replace cloning, but their per-read error rate is "
    "close to the mutation rate of interest, and general-purpose variant callers judge each base by its quality "
    "score rather than by the errors a given library actually makes. We built SHMtrace to calibrate every call "
    "against matched AID-null libraries sequenced in the same run, to count animals rather than reads as "
    "replicates, and to make the whole workflow, from reads or chromatograms to group statistics, reproducible.")
HIGHLIGHTS = [
    "Each candidate mutation is tested against errors of matched AID-null libraries",
    "Background in AID-null controls falls to {ctl_rec} calls per assay while signal is kept",
    "Nanopore and Sanger sequencing of the same mice give concordant mutation spectra",
    "Long reads add whole-amplicon maps, B cell clones and mutation trees",
]
ETOC = ("Aderinto and Martin describe SHMtrace, which measures somatic hypermutation from nanopore, Sanger or "
        "Illumina reads by testing each candidate mutation against errors measured in matched AID-null libraries. "
        "Across six loci and three platforms it suppresses background while keeping signal, and long reads add "
        "clones and mutation trees.")

INTRO = [
    "Somatic hypermutation (SHM) diversifies immunoglobulin genes in germinal-center B cells. It is initiated by "
    "activation-induced cytidine deaminase (AID), which deaminates cytosines in immunoglobulin DNA; the resulting "
    "uracils are copied by replication or processed by uracil excision and mismatch repair into the full spectrum "
    "of substitutions at G:C and A:T base pairs[@muramatsu2000,@dinoia2007,@methot2017]. Because each repair route "
    "leaves its own footprint, the rate, position and spectrum of mutations are the standard readouts for studies "
    "of AID, its cofactors and the repair pathways that act on its lesions[@dinoia2002,@rada2004,@rada2002,@zeng2001]. "
    "In mice the most widely used reporter is the intron downstream of J_{H}4, a non-coding region that mutates "
    "like a V gene but is not subject to antigen selection[@jolly1997,@lebecque1990].",
    "The assay has changed little in decades. Amplicons are cloned into bacteria and sequenced colony by colony, "
    "which takes days of hands-on work and yields tens of sequences per animal[@sible2021]. High-throughput "
    "alternatives exist, but short reads fragment the amplicon and lose the linkage between mutations on one "
    "molecule[@martin2023]. Long reads from Oxford Nanopore Technologies (ONT) now reach modal accuracies above 99% "
    "with R10.4.1 flow cells[@zhang2023,@sereika2022,@sanderson2024] and have been used to sequence immunoglobulin "
    "genes where consensus can absorb error[@minervini2021,@subassatish2022]. Calling rare substitutions across a "
    "polyclonal population from single reads is harder: a wild-type J_{H}4 intron carries roughly one mutation per "
    "100 bases, while the residual error of a single read, after filtering, is within one or two orders of magnitude "
    "of that, and it is concentrated at particular positions and contexts[@delahaye2021].",
    "Two habits of analysis make this worse. First, general-purpose low-frequency variant callers judge each base by "
    "its own quality score[@wilm2012], and simple rules applied in commercial software test candidates one by one "
    "without correcting for the number of tests[@kearse2012]. Neither knows the errors a particular library makes: "
    "errors that recur at a position, from PCR, homopolymers or a co-amplified paralog, look like mutations. "
    "Second, reads from one animal are often pooled and treated as independent replicates, which overstates "
    "certainty because reads of one animal share clones and PCR history[@lazic2010]. SHM studies already include the "
    "control that addresses both problems: libraries from AID-deficient animals or cells, prepared and sequenced "
    "alongside the samples. Existing SHM software summarizes mutations once they are called "
    "[@maccarthy2009,@gupta2015,@yaari2013], but does not use these controls to decide which differences are "
    "mutations.",
    "Here we describe SHMtrace, which calibrates every call against matched AID-null libraries from the same run, "
    "tests groups with the animal as the unit, and runs the whole workflow from delivered reads or Sanger "
    "chromatograms to group statistics. We validated it against planted mutations, simulated B cell lineages, "
    "six loci from five laboratories, three sequencing platforms and Sanger sequencing of the same mice, and "
    "compared it with LoFreq[@wilm2012] and with the uncorrected test used in Geneious[@kearse2012].",
]

RESULTS = [
    ("h", "Design: calls are calibrated against the controls, not against quality scores"),
    ("p", "SHMtrace works from the files a sequencing provider or core facility returns (Figure 1). Nanopore and "
          "Illumina reads are gated on length and expected error[@edgar2015], aligned with minimap2[@li2018] against "
          "the target together with competing references, and reads that align better to a paralog, are chimeric "
          "or carry ten or more paralog-specific bases are removed (Figure 1B). A random draw of the same number "
          "of usable reads is analyzed from every library. Sanger chromatograms enter through a sample sheet and are "
          "trimmed to their clean core; UMI-tagged reads are collapsed to one consensus per molecule."),
    ("p", "The central step is the test (Figure 1A and 1C). For every substitution, strand and base quality, "
          "SHMtrace measures how often the AID-null libraries of the same run show that change, and asks whether a "
          "candidate is supported by more reads, at the qualities observed, than this background explains. Because "
          "quality scores overstate certainty, the control error rates are multiplied by a margin chosen by "
          "leave-one-control-out: each control is called against the others, and the smallest margin at which no "
          "more than 1% of its tests reach p < 0.01 is kept. The same held-out calls, made at the samples' depth, set the "
          "threshold: each control is a decoy library, and SHMtrace keeps the most permissive p-value cutoff at which "
          "the decoys, scaled to the bases tested, predict no more than 5% of a library's calls to be false, a "
          "target–decoy estimate of the false-discovery rate (FDR)[@elias2007] that, unlike the Benjamini–Hochberg "
          "procedure[@benjamini1995], does not assume uniform p-values for candidates selected because they were seen. "
          "Sites at which half or more of the reads differ from the reference are removed as germline or parental. "
          "Each control library is itself called against the other controls, never against itself or its own decoy, "
          "so control calls are an honest estimate of residual background."),
    ("p", "Groups are compared with one value per animal by a negative-binomial test with a parametric bootstrap, "
          "a Welch test on log rates and a permutation test (Figure 1D). Every project also reports the substitution "
          "spectrum, positional clustering of calls, a 5-mer targeting index against the SHazaM S5F model[@yaari2013], "
          "an aggregate estimate of mutations carried by single reads, clones identified by their V(D)J junction, "
          "mutation trees, repair-route indices and a power calculator. The same engine runs from the command line "
          "and from a desktop application."),

    ("h", "Planted mutations and simulated lineages calibrate recall and false calls"),
    ("p", "We first measured performance where the truth is known (Figure 2). Into reads from one AID-null library, "
          "with the other controls as background, we planted 60 substitutions at AID-weighted sites, each in 1–20 "
          "reads, keeping every read's own base qualities, and repeated this three times for each held-out control. "
          "On J_{H}4, SHMtrace recalled {sp_jh4_recall}% of planted mutations with {sp_jh4_fp} false calls per library (observed "
          "FDR {sp_jh4_fdr}), and {sp_jh4_ge3}% of those carried by three or more reads (Figure 2A). The uncorrected Geneious "
          "rule recalled more single-read mutations ({sp_jh4_gen_1}% against {sp_jh4_1}%) but made more false calls "
          "({sp_jh4_gen_fp} per library, FDR {sp_jh4_gen_fdr}). LoFreq at its defaults called none of the planted mutations: its base-alignment "
          "quality correction, designed for short reads, suppresses nanopore calls almost entirely. With that "
          "correction switched off, LoFreq recalled {sp_jh4_lfnb}%, all in mutations carried by five or more reads, and a "
          "permissive setting raised recall to {sp_jh4_lfp}% at {sp_jh4_lfp_fp} false calls per library (Figures 2A and 2C)."),
    ("p", "N114 V-region reads behaved similarly (Figure 2B): SHMtrace recalled {sp_n114_recall}% of planted "
          "mutations, including {sp_n114_1}% of those in a single read, against {sp_n114_lfnb}% for LoFreq without "
          "base-alignment quality, with {sp_n114_fp} false calls per library (FDR {sp_n114_fdr}); the Geneious rule "
          "recalled {sp_n114_gen_recall}% with {sp_n114_gen_fp} false calls (FDR {sp_n114_gen_fdr}). Empty-vector "
          "subclone C3 had been sequenced as two libraries, and these results use the second. With the first in its "
          "place the FDR was {c31_fdr}, and {c31_held_fdr} when that library was the one held out: it was noisier "
          "than the peers that serve as its decoys, which no control-based estimate can anticipate. The two "
          "libraries shared no calls, so the noise was the library's, not the subclone's. Benjamini–Hochberg correction, which assumes uniform p values "
          "among candidates that were selected because they were seen, gave the same FDR with these controls "
          "({sp_n114_bh_fdr}) but {c31_bh_fdr} with the first C3 library, against {c31_fdr} for the control-based "
          "threshold. On the independent laboratory's Sμ and *Pim1* control libraries (2,000 reads), "
          "SHMtrace recalled {sp_su_lps_recall}% and {sp_pim1_recall}% of planted mutations at an observed FDR of "
          "{sp_su_lps_fdr} and {sp_pim1_fdr}, the Geneious rule {sp_su_lps_gen_recall}% and {sp_pim1_gen_recall}% "
          "(FDR {sp_su_lps_gen_fdr} and {sp_pim1_gen_fdr}), and LoFreq {sp_ext_lf_lo}–{sp_ext_lf_hi}% across its "
          "three settings, with an FDR of {sp_ext_lfp_fdr} when made permissive (Figure 2C)."),
    ("p", "Planted single sites cannot test what depends on clones, so we simulated B cell lineages on the J_{H}4 "
          "control reads: about 40 clones per library with power-law sizes, each with its own V(D)J junction, trunk, "
          "branch and private mutations placed by the S5F model (Figure 2D). Mutations in three or more reads were "
          "recalled by both SHMtrace and the Geneious rule; SHMtrace recalled {lin_1_rec}% of single-read mutations against {lin_1_gen}%, at an FDR of "
          "{lin_fdr_rec} against {lin_fdr_gen} ({lin_fp_rec} and {lin_fp_gen} false calls per library: a threshold "
          "set as a share of calls admits more errors where there are many real mutations), and in null libraries "
          "with nothing planted it made {lin_null_rec} false calls per library against {lin_null_gen}. Clones "
          "were recovered from the junction with an adjusted Rand index[@hubert1985] of {lin_ari}. When one unmutated "
          "clone held 60% of the reads, the read-weighted mutation frequency was correctly low ({lin_rw_est} against a truth of "
          "{lin_rw_true} per 10^{3} bp), while the clone-weighted frequency recovered {lin_cw_pct}% of the clones' true "
          "load ({lin_cw_est} against {lin_cw_true}; Figure 2E)."),
    ("p", "On the same simulated libraries the read-weighted mutation frequency was within {q_rw_dev}% of the truth at 500 "
          "reads on J_{H}4 and N114 (Figure 2F). The number of distinct mutations per base sequenced, the measure used "
          "for SHM density, was {q_ud_500} of the truth at 500 reads and {q_ud_1000} at 1,000 on J_{H}4 ({q_ud_n114_1000} on N114), "
          "falling with depth because a per-site test confirms a smaller share of single-read mutations as columns deepen; with the "
          "aggregate single-read estimate added it was {q_cd} across scenarios. Densities are therefore compared at matched depth. In "
          "planted-mutation tests of that estimate, 25, 50 and 100 one-read mutations raised it by {sr_rec}% of the "
          "number planted on average (Figure S5)."),
    ("p", "Finally, we asked how often group tests declare a difference when there is none (Figure 2G). In 300 "
          "simulated experiments per condition with a between-animal coefficient of variation of 1, as measured in "
          "wild-type J_{H}4, pooling reads across animals produced p < 0.05 in 91% (3 vs 3) and 94% (5 vs 5) of "
          "experiments without a true difference. The negative-binomial test held the nominal 5% (5.0% and 5.7%) "
          "and had the most power of the calibrated tests; with three animals per group, a two-fold change was "
          "detected in 12% of experiments and a four-fold change in 30%."),

    ("h", "On real libraries, matched controls set the background"),
    ("p", "We then compared callers on real libraries with matched AID-null controls: our J_{H}4 intron libraries "
          "from germinal-center B cells of three wild-type and three *Aicda*^{−/−} mice, N114 hybridoma subclones "
          "transduced with AID or empty vector, and three loci from an independent laboratory's public nanopore "
          "deposit, the Sμ switch region in LPS-activated B cells and in CH12 cells, and the AID off-target *Pim1* "
          "(Table 1; Figure 3). Every caller analyzed exactly the same reads."),
    ("p", "SHMtrace made {ctl_rec} calls in the three controls of each assay, against {ctl_gen} for the Geneious rule, while "
          "keeping a wild-type signal similar to that rule's (Figure 3A–3E and 3G). The separation between AID-active and "
          "AID-null libraries (geometric means, controls without calls placed at one event) was {sep_jh4}-fold on J_{H}4, "
          "{sep_n114}-fold on N114, {sep_su_lps}-fold on Sμ in LPS-activated B cells and {sep_su_ch12}-fold in CH12 cells, "
          "against {gen_lo}- to {gen_hi}-fold for the Geneious rule and {lf_lo}- to {lf_hi}-fold for the three LoFreq settings (Figure 3F). On J_{H}4, "
          "SHMtrace called {jh4_wt_rec} mutations in the wild-type libraries and {jh4_ko_rec} in the controls; the "
          "Geneious rule called {jh4_wt_gen} and {jh4_ko_gen}. LoFreq with base-alignment quality called {jh4_wt_lf} "
          "wild-type mutations; without it, Wt2, a library dominated by one clone with many single-read mutations, "
          "fell from {jh4_wt2_rec} calls to {jh4_wt2_lfnb}. On the external Sμ sets SHMtrace called {su_wt_calls} wild-type mutations per "
          "library, against {su_gen_wt_calls} for the Geneious rule, with {su_ctl_calls} per control against "
          "{su_gen_ctl_calls}."),
    ("p", "*Pim1*, where AID acts weakly, marked the limit. At 2,000 reads no caller separated AID-competent from "
          "*Aicda*^{−/−} libraries convincingly (SHMtrace {pim1_rec_fold}-fold, {pim1_rec_ctl} control calls; Geneious rule {pim1_gen_fold}-fold, "
          "{pim1_gen_ctl} control calls), {pim1_deep}"),
    ("p", "The leave-one-control-out margin did what it was designed to do (Figure 3H). Taken at face value, the "
          "control error rates let {loco_jh4_face}% of held-out J_{H}4 tests and {loco_n114_face}% of N114 tests reach "
          "p < 0.01 at the samples' depth; the chosen margins ({loco_margins}) brought each to 1% or below "
          "({loco_jh4_final}%, {loco_full_final}% and {loco_n114_final}%). Read preparation removed 44–79% of the "
          "reads delivered for each J_{H}4 library (67% overall), chiefly at the expected-error gate (29%) and at "
          "competitive alignment (25% unmapped or aligning best to the paralog or contaminant; Figure 3I). Without the paralog-allele "
          "filter, full-length paralog reads that align better to the 1,230-bp target than to the shorter paralog "
          "record contributed 35–54 false calls per simulated library."),

    ("h", "Nanopore and Sanger sequencing agree in the same mice"),
    ("p", "The same six mice had been analyzed by cloned-amplicon Sanger sequencing[@barbulescu2025]. We re-analyzed "
          "their 447 chromatograms with SHMtrace's Sanger mode and compared both platforms over the 529-bp window "
          "the Sanger clones read cleanly (Figure 4). Both separated wild-type from *Aicda*^{−/−} mice (Sanger "
          "{sanger_fold}-fold, p = {sanger_p}; nanopore at least {ont_fold}-fold, p = {ont_p}; Welch tests on log densities, n = 3 per "
          "group; Figure 4A). Pooled substitution spectra agreed (Pearson r = {r_spec}, 95% CI {r_spec_ci}, slope {slope}; "
          "Figure 4B), as did the distribution of calls along the intron (r = {r_pos} in 25-bp bins; Figure 4C) and "
          "hotspot use (G:C calls at WRC/GYW, {hot_gc_s}% Sanger and {hot_gc_o}% nanopore; A:T calls at WA/TW, {hot_at_s}% and "
          "{hot_at_o}%, against {hot_gc_exp}% and {hot_at_exp}% of such sites in the window; Figure 4D). A repeated-measures test found no platform-by-type "
          "interaction (p = {interaction_p}). The absolute densities differ, by design: a Sanger clone is one molecule read "
          "once, and nanopore densities count distinct mutations among 500 reads."),
    ("p", "Re-analyzing the chromatograms also exposed bookkeeping errors in the original clone collection, such as "
          "clones filed under the wrong mouse and renamed copies of other experiments' traces. SHMtrace's Sanger "
          "mode therefore works from an explicit sample sheet and flags, without correcting, copied traces, clones "
          "shared between animals beyond chance, mutated clones under controls and mislabeled tubes. In 40 simulated "
          "projects with one of each error planted, every copied or re-typed trace, contaminated control and "
          "misfiled clone with its label intact was flagged, as were all 125 clones of misfiled lineages, with no "
          "flag on 6,018 correctly filed clones (Figure 4E). A misfiled clone whose label had also been changed was "
          "flagged in 22 of 40 projects, depending on how many mutations it shared with its lineage."),

    ("h", "Long reads add whole-amplicon maps, clones and mutation trees"),
    ("p", "Each nanopore read spans the 1,230-bp J_{H}4 amplicon, so the whole product can be analyzed (Figure 5A). "
          "Over the reliably covered positions 90–1,230, wild-type libraries carried {full_fold} times the density of the "
          "controls; {full_450}% of wild-type calls fell in the first 450 bp, and the number of calls fell about "
          "{full_decline}-fold "
          "between the first and last 250 bp, consistent with the decline of SHM about 1 kb from the V(D)J "
          "junction[@lebecque1990]. Spectra over the full interval matched those over the Sanger window (r = {s1_r})."),
    ("p", "Calls in wild-type libraries favored the positions the S5F model considers most mutable: the targeting "
          "index, the model's mutability at mutated positions relative to the window, was {s5f_wt} (p ≤ {s5f_wt_p}, "
          "permutation), and calls per position rose {s5f_dose}-fold from the least to the most mutable fifth of positions "
          "(Figure 5B). Every position at which an *Aicda*^{−/−} read mismatched the reference, which is sequencing "
          "error, scored {s5f_ko} (p = {s5f_ko_p}). This check is independent of the error model and identifies calls that do "
          "not look like AID's."),
    ("p", "Because J_{H}4 amplicons start in the rearranged V region, the bases upstream of J_{H}4 are the V(D)J "
          "junction, a barcode of the B cell clone (Figure 5C). In Wt2, one clone held {wt2_largest}% of the reads with a "
          "junction and was largely unmutated, which lowered the read-weighted mutation frequency to {wt2_rw} per "
          "10^{3} bp; weighting each clone once gave {wt2_cw}. Wt3 and Wt4, without a dominant clone, barely moved "
          "({wt3_rw} to {wt3_cw}; {wt4_rw} to {wt4_cw}). Within a clone, the reads can be ordered into a mutation tree (Figure 5D). The largest "
          "Wt3 clone ({tree_reads} reads) carried {tree_founders} founder mutations in every mutated read and then "
          "divided into {tree_sub} sub-lineages up to {tree_depth} mutations deep; {tree_fit}% of its mutated reads "
          "fitted a single tree, and {tree_incompat} mutation pairs failed the four-gamete test. Without unique molecular identifiers (UMIs), PCR recombination "
          "can create such incompatible pairs, so these trees are read-level reconstructions."),

    ("h", "The same engine works across platforms, loci and biological questions"),
    ("p", "Three public Illumina datasets tested SHMtrace where the platform, laboratory and question differed "
          "(Figure 6). In splenic germinal-center B cells sequenced across J_{H}4 (positions 53–122), wild-type mice "
          "carried on average {ldha_fold} times the mutation frequency of AID-null mice (Welch test on log rates, "
          "p = {ldha_p}; n = {ldha_nwt} and {ldha_nko}), and mice lacking LDHA in activated B cells did not differ "
          "from wild type (Mann-Whitney p = {ldha_p_ko}), "
          "the conclusion of the original study[@sharma2023] (Figure 6A)."),
    ("p", "In CH12F3 cells lacking AID, UNG and MSH2, ectopic AID raised the mutation frequency {ch12_fold}-fold over "
          "empty vector (p = {ch12_p}; four replicates each), with {ch12_ev}[@wang2023]. "
          "Without uracil excision and mismatch repair, AID lesions can only be copied as C→T or G→A, and "
          "{ch12_ctga}% of AID calls were (Figure 6B). The catalytically inactive R50E mutant {ch12_r50e}. A "
          "second mutant appeared to mutate more than AID, but its calls had the wrong spectrum (C→T/G→A "
          "{ch12_mut_ctga}%), {ch12_cluster_pct}% of them within {ch12_cluster_bp} bp (positions {ch12_cluster}) across replicates, and scored a targeting index of "
          "{ch12_mut_s5f} against {ch12_aid_s5f} for AID: an error cluster specific to those libraries, which a "
          "control-fitted null cannot see because the controls do not share it. SHMtrace flagged it by positional "
          "clustering and the spectrum (Figure 6C)."),
    ("p", "In Ramos B cells, loss of HMCES is reported to increase SHM-associated deletions specifically[@wu2022]. "
          "With indel calling on, HMCES-deficient clones carried {hmces_del_fold}-fold more deletions than wild type "
          "(p = {hmces_del_p}, n = 6 clones each) against {hmces_ins_fold}-fold for insertions (Figure 6D). One "
          "AID-knockout control clone carried {hmces_ko3_calls} substitution calls that were not AID-like "
          "({hmces_ko3_gc}% at G:C); its two siblings had none. As a decoy, that clone raised the control-based "
          "threshold so far that no substitution was called in any other clone, and SHMtrace now names such a "
          "control in its notes. Analyzed as a sample instead, with the two clean controls, wild-type and "
          "HMCES-deficient clones carried {hmk3_calls} substitution calls each, {hmk3_gc}% of them at G:C, as "
          "expected of AID."),
    ("p", "UMIs remove the PCR and sequencing error of single reads. On public R10 nanopore reads of rRNA operons "
          "with dual UMIs[@karst2021], SHMtrace's consensus reduced the error per base 57-fold, to 0.17% with "
          "residual errors mostly in homopolymer indels; substitutions fell to about 1 in 10,000 bases with 10 or "
          "more reads per molecule (Figure 6E), 100-fold below a wild-type J_{H}4 mutation frequency."),
]

DISCUSSION = [
    "SHMtrace replaces the two assumptions that make low-frequency calling from single reads unreliable for SHM. "
    "Instead of trusting quality scores, it measures the background from the AID-null libraries that SHM studies "
    "already include, at the resolution of substitution, strand and quality, and it checks that calibration by "
    "holding each control out in turn. Instead of counting reads as replicates, it tests animals. Both choices "
    "matter in practice: a quality-based caller either missed most real mutations on nanopore reads or, made "
    "permissive, called as much background as the uncorrected rule, and pooling reads produced false positives in "
    "more than 90% of experiments without a true difference.",
    "The comparison with LoFreq should be read for what it is. LoFreq was designed for deep short-read data, where "
    "base qualities are well calibrated and the base-alignment quality correction removes alignment artifacts; on "
    "nanopore reads that correction suppressed almost every call, and without it LoFreq still needed five or more "
    "supporting reads. This is not a defect of LoFreq on its intended data but an illustration that a per-base error "
    "model cannot know the library-specific background of an SHM assay, which only a matched control measures. The "
    "Geneious rule recalled more single-read planted mutations on J_{H}4 and *Pim1* and fewer on N114 and Sμ, "
    "and made more false calls on every amplicon; in polyclonal simulations, where single-read mutations sit among many real calls, SHMtrace recalled "
    "{lin_1_rec}% of them and the Geneious rule {lin_1_gen}%, at more false calls but within the nominal FDR. For counting mutations confined to one read, SHMtrace reports an aggregate background-subtracted "
    "estimate rather than individual calls.",
    "Long reads change what an SHM experiment can report. The whole amplicon is read at once, so the 3′ decline of "
    "mutation is visible without tiling; the V(D)J junction identifies clones, so an expanded or jackpotted clone no "
    "longer dominates the mutation frequency; and the mutations of a clone can be ordered into a tree. Sanger "
    "sequencing of the same mice, re-analyzed from the chromatograms, agreed on spectrum, position and hotspot use. "
    "Re-analysis also showed that the weakest link in the traditional assay is bookkeeping, which SHMtrace now "
    "checks explicitly.",
    "Several features are designed for routine use: the analysis starts from delivered reads, every step's read "
    "counts are reported, the random draw is reproducible from a fixed seed, and a desktop application exposes the "
    "same engine with reports a non-specialist can read. The power calculator makes a practical point plain: "
    "between-animal variation, not sequencing depth, limits SHM experiments, and with wild-type variation as in "
    "J_{H}4, a two-fold change between AID-proficient groups needs more than ten animals per group.",
]
LIMITATIONS = (
    "SHMtrace requires AID-null control libraries from the same sequencing run; across runs, errors and basecallers "
    "change, and a cross-run background ranks samples but does not calibrate them. Its calls are filtered apparent "
    "variants: without UMIs, early PCR errors, PCR recombination and duplicates are indistinguishable from mutations "
    "and from expanded clones, and the trees are read-level reconstructions. Single-read mutations are recovered in "
    "aggregate, not individually, and the distinct-mutation density depends on depth, so libraries must be compared "
    "at matched depth. The FDR is estimated from the controls themselves, so a control library noisier than its peers, "
    "as the first library of N114 subclone C3, lets the FDR exceed the nominal level, "
    "and a control far noisier than its peers suppresses calls in every sample until it is set aside. An artifact present in one library and absent from the controls is not removed by "
    "calibration; spectrum, clustering and targeting checks flag such cases but need a person to judge them. At a "
    "weak AID target such as *Pim1*, signal approaches the background at the depths tested. Our own validation used "
    "three animals per group, one nanopore run per assay and one locus per cell system; the thresholds were chosen on "
    "our libraries and transferred to the public sets without retuning, but a held-out sequencing run of new animals "
    "would test them further. Indel calling was validated on Illumina reads only, and UMI consensus on rRNA operons, "
    "not on an SHM library. Where UMIs or duplex reads are available, molecule-level consensus is preferable to "
    "calibration of single reads, and SHMtrace supports both.")

RESOURCE = [
    ("Lead contact", "Requests for further information and resources should be directed to and will be fulfilled by "
                     "the lead contact, Alberto Martin (alberto.martin@utoronto.ca)."),
    ("Materials availability", "This study did not generate new unique reagents."),
    ("Data and code availability", "Nanopore reads of the J_{H}4 and N114 libraries, the Sanger chromatograms and "
     "clone records, SHMtrace project files, per-library read accounting and all benchmark outputs have been "
     "deposited at [repository] and are publicly available as of the date of publication; accession numbers are "
     "listed in the key resources table. The public datasets analyzed are listed in Table 1 and the key resources "
     "table. SHMtrace (v0.2), including the desktop application, and all analysis and figure code have been "
     "deposited at [Zenodo DOI] and are publicly available as of the date of publication. Any additional "
     "information required to reanalyze the data reported in this paper is available from the lead contact upon "
     "request."),
]
ACK = ("We thank the University of Toronto flow cytometry facility and members of the Martin laboratory for "
       "discussions. Long-read sequencing was performed by Plasmidsaurus, and cloned-amplicon Sanger sequencing by "
       "The Centre for Applied Genomics (SickKids, Toronto). This work was supported by the Canadian Institutes of "
       "Health Research (PJT-180269).")
CONTRIB = ("Conceptualization, A.A. and A.M.; methodology, A.A.; software, A.A.; validation, A.A.; formal analysis, "
           "A.A.; investigation, A.A.; data curation, A.A.; writing – original draft, A.A.; writing – review & "
           "editing, A.A. and A.M.; visualization, A.A.; supervision, A.M.; funding acquisition, A.M.")
INTERESTS = "The authors declare no competing interests."
AI_DECL = ("During the preparation of this work, the authors used Anthropic Claude to assist with code drafting, "
           "language editing and figure layout. After using this tool, the authors reviewed and edited the content "
           "as needed and take full responsibility for the content of the published article.")
