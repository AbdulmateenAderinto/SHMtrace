"""Editorial proposals only. Importing this file does not edit manuscript sources.

Indices refer to source_blocks.json, captured from v1.3. Numerical placeholders stay live.
"""

import json
from pathlib import Path

RESULTS = {
    1: "SHMtrace starts from delivered reads or Sanger chromatograms (Figures 1 and S1). Read libraries undergo length and expected-error filtering[@edgar2015], competitive alignment with minimap2[@li2018], paralog checks and matched-depth sampling. Sanger traces are assigned through a sample sheet and quality-trimmed; UMI-tagged reads can be collapsed to molecular consensuses (STAR Methods).",
    2: "For each substitution, SHMtrace tests read support against control error rates stratified by substitution type, strand and quality, incorporating position-specific background (Figures 1A and 1C). A leave-one-control-out margin adjusts underestimated error (Figure S2); held-out controls then set an empirical target–decoy threshold targeting 5% false discoveries[@elias2007]. This estimate depends on representative controls and is not a general FDR guarantee. Data-dependent candidate selection motivates empirical calibration rather than assuming valid Benjamini–Hochberg p-values[@benjamini1995]. Alleles in ≥50% of reads are excluded as likely germline or parental differences; genuine dominant-lineage mutations can also cross this heuristic threshold, and exclusions are reported.",
    3: "Substitution calls can have one supporting read. A separate Sanger route counts Q≥20 substitutions and reports expected control background without subtraction (Figure 1C). Windows, germline handling, quality checks and comparisons across biological replicates are shared. Sanger reanalyses therefore evaluate counting and provenance, rather than the calibrated read-library test.",
    4: "Two C-to-T candidates in Sμ {ana_case} had similar qualities and support on both strands (Figure 2). Position {ana_t_pos}, in a WRC hotspot, lacked Q≥20 control support and was retained (p = {ana_t_p}); position {ana_e_pos} recurred in controls and was rejected (p = {ana_e_p}). Of {ana_cand} candidates, {ana_kept} met the threshold (p ≤ {ana_thr}); the remainder were unresolved against background.",
    5: "Group comparisons use one value per biological replicate, with log-rate Welch, bootstrapped negative-binomial and permutation tests (Figure 1D). Outputs include spectra, hotspot use, S5F targeting[@yaari2013], clustering, aggregate single-read estimates, junction-defined clones, read-level trees and power estimates (Figure S3). Desktop and command-line interfaces share the analysis engine.",
    7: "We planted 60 substitutions per held-out control, preserving qualities and varying support from 1–20 reads (Figure 3). On J_{H}4, SHMtrace recalled {sp_jh4_recall}% overall and {sp_jh4_ge3}% of events with ≥3 supporting reads, with {sp_jh4_fp} false calls per library (observed proportion {sp_jh4_fdr}). In contrast, Geneious recalled more single-read events ({sp_jh4_gen_1}% versus {sp_jh4_1}%) with {sp_jh4_gen_fp} false calls. Default LoFreq recalled none; BAQ-off and permissive settings recalled {sp_jh4_lfnb}% and {sp_jh4_lfp}% (Figures 3A and 3D; Table S1).",
    8: "ClairS[@zheng2026] and DeepSomatic[@park2026] received peer controls as matched normals. Neither recovered planted events at defaults. Across 2,160 substitutions, permissive ClairS recovered none; DeepSomatic recalled 19%, 18%, 5% and 3% on J_{H}4, N114 V_{H}, Sμ and *Pim1*, respectively. Recall improved at higher allele fractions (Figure 3C), whereas {af_lt1}% of {af_n} retained SHMtrace calls occurred in <1% of reads (Figure 2E). Germline callers Clair3 and DeepVariant recovered none, consistent with their intended use (Table S1).",
    9: "On N114 V_{H}, SHMtrace recalled {sp_n114_recall}% overall and {sp_n114_1}% of single-read substitutions, with {sp_n114_fp} false calls per library (observed proportion {sp_n114_fdr}; Figure 3B). On external Sμ and *Pim1*, recall was {sp_su_lps_recall}% and {sp_pim1_recall}% at proportions {sp_su_lps_fdr} and {sp_pim1_fdr}, respectively (Figure 3D; Table S1).",
    10: "At the default 5% target, pooled false-discovery proportions were {fc_def_lo}–{fc_def_hi}% across loci (Figure S4). A 1% target yielded {fc_strict_lo}–{fc_strict_hi}% on J_{H}4, Sμ and *Pim1* and {fc_n114_1}% on N114 V_{H}; J_{H}4 recall fell from {fc_jh4_rec_def}% to {fc_jh4_rec_1}%. Stricter targets therefore did not consistently improve calibration.",
    11: "J_{H}4 lineage recall rose from {cd_rec_1}% with one control to {cd_rec_2}% with two and {cd_rec_34}% with three or four (Figure S5). One control requires a Benjamini–Hochberg fallback without decoy calibration. All-null simulations produced no calls ({cd_null_n} runs). On Sμ, controls from another experiment gave similar recall, and pooling controls reduced false calls; sequencing-run identity was unavailable.",
    12: "Diluting real hypermutated reads into null reads tested mutations with biological clonality and clustering (Figure S6). N114 V_{H} frequencies were {ti_n114_ratio} of expectation at 25–50% donor reads, with detection from {ti_n114_det}%. J_{H}4 detection depended strongly on the donor, beginning at {ti_jh4_det_lo}% for {ti_jh4_most} and {ti_jh4_det_hi}% for {ti_jh4_least}.",
    13: "For planted 1–40-bp indels with ≥3 supporting reads, recall was {ind_jh4_ge3}% on J_{H}4 and {ind_n114_ge3}% on N114 V_{H}, with {ind_jh4_fp} and zero false calls per library (Figure S7). One- and two-read events were excluded by design. Single-base homopolymer deletions were least reliably recovered on J_{H}4 ({ind_jh4_hp}%).",
    14: "Simulations of approximately 40 junction-defined clones included shared, branch and private mutations (Figure 3E). SHMtrace recalled {lin_1_rec}% of single-read events versus {lin_1_gen}% for Geneious, at false-discovery proportions {lin_fdr_rec} versus {lin_fdr_gen}; null-library false calls averaged {lin_null_rec} versus {lin_null_gen}. Clone assignments had adjusted Rand indices[@hubert1985] of {lin_ari}. With an unmutated clone contributing 60% of reads, majority-rule clone weighting recovered {lin_cw_pct}% of its corresponding true burden (Figure 3F).",
    15: "A mutated dominant clone exposed the ≥50% filter's limitation: mutations crossing the threshold at 60–80% clone abundance were removed, reducing read- and clone-weighted estimates by up to {mask80_rw}% and {mask80_cw}%. No loss occurred at 40%; excluded sites remain available for inspection (STAR Methods).",
    16: "Read-weighted frequencies were within {q_rw_dev}% of simulated truth at 500 reads (Figure 3G). Unique-variant density varied with depth because repeated variants increase exposure without adding discoveries, while rare-event detection also changes. Adding the aggregate single-read estimate recovered 96–100% of truth; independently planted single-read events were recovered at {sr_rec}% in aggregate (Figure S8A). Unique-density comparisons require matched depth.",
    17: "Across 300 simulations per condition, pooled-read tests falsely rejected in 91–94% of null experiments, versus 5.0–5.7% for the negative-binomial test (Figure 3H). With three animals per group, negative-binomial power was only 12% for two-fold and 30% for four-fold differences. Read depth therefore cannot replace biological replication.",
    19: "Caller comparisons used identical reads from J_{H}4 in three Wt and three *Aicda*^{−/−} mice, AID-expressing and empty-vector N114 subclones, and external Sμ and *Pim1* libraries (Figure 4; Tables S2 and S3).",
    20: "SHMtrace retained {ctl_rec} calls across three controls per assay versus {ctl_gen} for Geneious (Figures 4A–4E and 4G). Geometric-mean unique-density separation was {sep_jh4}-, {sep_n114}-, {sep_su_lps}- and {sep_su_ch12}-fold for J_{H}4, N114 V_{H}, primary-cell Sμ and CH12 Sμ, respectively, using one-event values for zero-call controls (Figure 4F). At 500 reads, J_{H}4 and N114 V_{H} reached {depth_jh4_500}% and {depth_n114_500}% of discoveries at the largest tested depth, without establishing saturation (Figures S9 and S10). Background mismatches were position-dependent and separable from abundant clone-shared calls (Figure S11).",
    21: "On the same N114 reads, IgBLAST[@ye2013]/SHazaM[@gupta2015] mismatch counting gave {im_fold_imgt}-fold AID/empty-vector frequency separation using IMGT germline, {im_fold_par}-fold using the parental sequence and {im_fold_q20}-fold after Q20 masking, versus {im_fold_sh}-fold for SHMtrace (Figure S12). These are frequency contrasts, distinct from Figure 4's unique-density endpoint.",
    22: "At weakly mutated *Pim1*, SHMtrace gave {pim1_rec_fold}-fold separation at 2,000 reads, with four control calls; Geneious gave {pim1_gen_fold}-fold with eleven. Full available depth still yielded overlapping case/control counts, leaving the signal close to background (Table S1).",
    23: "At two CH12 off-targets[@seija2026], Wt unique densities were {ot1_wt_dens} and {ot2_wt_dens} per 10^{5} bases versus {ot1_ko_dens} and {ot2_ko_dens} in *Aicda*^{−/−} controls (log-density Welch p = {ot1_p_dens} and {ot2_p_dens}). MLLT1/MLLT3 double knockout reduced densities to {ot1_dko_dens} and {ot2_dko_dens}, with no detected difference from null controls (p = {ot1_p_dko_ko} and {ot2_p_dko_ko}). Zero counts precluded a convention-independent fold estimate.",
    24: "Leave-one-control-out error margins of five for the shared J_{H}4 window and ten for full-length J_{H}4 and N114 V_{H} reduced held-out candidate tests with p<0.01 to {loco_jh4_final}%, {loco_full_final}% and {loco_n114_final}%, respectively (Figure 4H). This calibrates the error model; the separate decoy threshold estimates the false-discovery proportion.",
    25: "Read preparation removed 67% of delivered J_{H}4 reads overall (Figure 4I). Competitive alignment required an additional paralog-base filter: without it, simulated paralog reads produced 35–54 false variants per library (STAR Methods).",
    27: "Of {sanger_submitted} processed Sanger chromatograms from the six nanopore-sequenced mice, {sanger_clones} passed analysis over the shared 529-bp window (Figures 5 and S13). Wt/null separation was {sanger_fold}-fold for Sanger and {ont_fold}-fold for nanopore (p = {sanger_p} and {ont_p}; n=3 mice per group; nanopore zeros assigned one event). Pooled spectra correlated at r = {r_spec} and positional profiles at r = {r_pos}, with similar hotspot enrichment. No platform-by-type interaction was detected (p = {interaction_p}), without establishing equivalence. Absolute densities remain platform-specific: Sanger counts clone-level substitutions; nanopore counts distinct variants among 500 reads.",
    29: "Reanalysis of published *Maea* J_{H}4 chromatograms[@maccarthy2009,@barbulescu2024] recovered a {c2_fold}-fold frequency reduction versus {c2_pub_fold}-fold published (one-tailed p = {c2_p} and {c2_pub_p}; Figure S14). Spectra correlated at r = {c2_spec_r}, but denominator-harmonized frequencies averaged {c2_ratio}-fold published values. Shared trace assignments and identical mutated-clone pairs were retained as filed and flagged; detailed discrepancies are in Table S4 and STAR Methods.",
    30: "Reanalysis of {mk_traces} *Mkln1* chromatograms[@barbulescu2025] found an {mk_fold}-fold reduction versus {mk_pub_fold}-fold published (p = {mk_p} and {mk_pub_p}). The G:C transversion increase persisted, while spectrum agreement varied between Wt and knockout (r = {mk_spec_r_wt} and {mk_spec_r_ko}; Figure S15; Table S4).",
    31: "An independent laboratory's SW_{HEL} study[@thientosapol2018] provided {sw_clones} sequences across five arms (Figure S16). SHMtrace recovered {sw_muts} of {sw_pub_muts} tabulated mutations (r = {sw_r} across 60 spectrum cells), matching mutations per sequence to two decimals. The four substitution-class comparisons agreed with the publication. Total mutation decreased in both analyses, but significance differed under different tests. A potentially pooled host was retained; excluding it did not change the reanalysis conclusions (STAR Methods). The *Msh2*^{−/−} arm showed the expected A:T collapse. Deposited sequences evaluate counting, spectra and group statistics, without testing chromatogram reading.",
    33: "Long reads extended J_{H}4 analysis from 529 to 1,141 bp (positions 90–1,230; Figure 6A). Wt density exceeded controls {full_fold}-fold; {full_450}% of mutations lay in the first 450 bp, with a seven-fold first-to-last-250-bp decline, consistent with earlier work[@lebecque1990]. Full-length and shared-window spectra remained similar (r = {s1_r}; Figure S17).",
    34: "Wt calls were enriched at S5F-predicted mutable positions (index 1.22–1.36; permutation p≤{s5f_wt_p}), whereas null mismatches were not (index {s5f_ko}; p = {s5f_ko_p}; Figure 6B). S5F is absent from the calling test, providing a separate check of SHM-associated targeting.",
    35: "Junction-defined clone weighting changed Wt2 frequency from {wt2_rw} to {wt2_cw} per 10^{3} bp because one largely unmutated lineage dominated the library (Figure 6C). In the largest Wt3 lineage, {tree_fit}% of mutated reads fitted one tree (Figure 6D). Without UMIs, residual error and PCR recombination remain possible; these reconstruct read relationships rather than proven cellular ancestry.",
    37: "In public Illumina J_{H}4 data[@sharma2023], Wt mutation frequency exceeded AID-null frequency {ldha_fold}-fold (log-Welch p = {ldha_p}; n={ldha_nwt} and {ldha_nko}). *Ldha* deletion showed no detected difference from Wt (Mann–Whitney p = {ldha_p_ko}; Figure 7A).",
    38: "In repair-deficient CH12F3 cells[@wang2023], ectopic AID increased frequency {ch12_fold}-fold over empty vector (p = {ch12_p}; n=4 each), entirely through C→T/G→A; R50E yielded no calls (Figure 7B). R174E/R178E calls instead showed an atypical spectrum, positional cluster and low targeting index, warranting provenance checks before biological interpretation (Figure 7C).",
    39: "HMCES-deficient Ramos clones[@wu2022] had {hmces_del_fold}-fold more deletions and {hmces_ins_fold}-fold more insertions (deletion p = {hmces_del_p}; n=6 per group; Figure 7D). Substitution calls depended strongly on a noisy AID-null library. A rule developed after inspecting this dataset excluded it from decoy thresholding, retaining it in model fitting and reporting the exclusion. Adjusted calls were {hm_calls} per Wt/HMCES-deficient clone; treating that control as a sample instead yielded {hmk3_calls}. This sensitivity limits substitution-burden inference (STAR Methods); the rule also acted in mOrange.",
    40: "Dual-UMI rRNA reads[@karst2021] provided a technical consensus benchmark: total error fell 57-fold (95% CI 35–99) to 0.17%, with residual homopolymer indels; substitution error approached 10^{−4} at ≥10 reads per family (Figure S8D). This does not validate UMI-tagged SHM analysis.",
    41: "In the mOrange study[@perezduran2012], inactive-AID controls calibrated active AID across empty-vector, Ugi and UNG contexts (Figure 7E). SHMtrace retained {mo_calls} substitutions, all at G:C, with {mo_wrc}% at WRC/GYW hotspots. Active/inactive separation was {mo_control_fold}-, {mo_ugi_fold}- and {mo_ung_fold}-fold, respectively; Ugi zeros used one-event values. Ugi-to-UNG transition shares fell from {mo_ugi_ts}% to {mo_ung_ts}%. The paired mean difference was {mo_arm_diff} percentage points (95% CI {mo_arm_ci}; {mo_arm_p}); few calls per library limit precision.",
    42: "On a common GFP substrate[@liu2018], AID showed greater WRC/GYW targeting than APOBEC3A in Wt and UNG/MSH2-deficient cells, while repair-deficient transition shares were similar (p = {eng_ko_ctga_p}; Figure S18A). No-editor libraries yielded at most one call in 20,000 reads. Transition enrichment alone therefore did not identify the enzyme.",
    43: "The second engineered dataset, a CRISPR-X antibody-evolution study[@devilder2019], supplied an independently curated mutation list. SHMtrace recovered {eng_dev_mut} of {eng_dev_mut} published mutations across {eng_dev_n} selected-cluster entries, adding none (Figure S18B). This tests variant identity in selected sequences, rather than population frequency or technical false-discovery rates.",
    44: "Public PI3Kδ-inhibitor Sμ libraries[@compagno2017] confounded genotype with read layout (Figure S19). In the repeat interval reached by paired libraries, nine null libraries yielded {smur_ko_tot} call versus {smur_wt_lo}–{smur_wt_hi} per Wt library (descriptive pooled p = {smur_p}); retained calls showed SHM-associated targeting. Only one Wt replicate per treatment reached these repeats, preventing a replicated drug comparison.",
    45: "PacBio Sμ–Sγ1 junctions[@vincendeau2022] tested junction structure rather than calibrated mutation calls (Figure S20). Mean microhomology was {pb_mh_wt} nt in Wt and {pb_mh_shld} nt in *Shld1*^{−/−}, rising to {pb_mh_xlf} nt without XLF and {pb_mh_dko} nt in double mutants. Double mutants showed more ≥4-nt microhomology and fewer blunt junctions; *Shld1*^{−/−} showed more insertions. These patterns support recovery of a known repair phenotype. Counts reflect PCR recovery rather than switching efficiency; {pb_amb_lo}–{pb_amb_hi}% of junctions had alternative equally good acceptor placements (STAR Methods).",
}

LEGENDS = [
    "(A) Quality-based versus control-calibrated decisions, schematically. SHMtrace estimates background by substitution, strand, quality and position from matched null controls. Check marks/crosses denote correct/incorrect schematic decisions. (B) Accepted inputs and preparation. Sanger Q≥20 trimming replaces the expected-error gate; Sanger clones are not downsampled. (C) Read-library substitutions undergo a control-fitted test and empirical false-discovery threshold; single-read support is allowed. Sanger Q≥20 substitutions are counted with expected background reported separately, without subtraction. Both share windows, quality checks and the heuristic ≥50% allele exclusion, with an opt-in known-ancestor exception. (D) Outputs and comparisons across biological replicates. Reads need not represent independent original molecules. FDR, false-discovery rate. See also Figure S1.",
    "(A, B) C-to-T support by strand in Sμ WT3 and controls (PRJNA1481264; 2,000 reads/library). Counts include all qualities; parentheses give Q≥20 subsets. Positions 165 and 97 show hotspot-associated and recurrent-background candidates. p denotes SHMtrace's test. (C) Supporting-base qualities; bars, medians; dotted line, Q20. (D) Retained candidates and depth-scaled decoy estimates across thresholds; the empirical 5% criterion retained 88 calls. (E) Allele fractions of retained calls at 500 reads for J_{H}4/N114 V_{H} and 2,000 for Sμ/*Pim1*. Shading marks 5–20% support tested in Figure 3C, rather than a universal caller boundary; dotted line, default 50% exclusion. See also Figure S2.",
    "(A, B) Planted-substitution recall by support: 60 sites × 3 draws × 3 held-out controls/locus, 90 events/support point; lines span controls. FDR labels are pooled false/(true+false) calls. (C) Permissive somatic callers on J_{H}4: solid, 1–20/500 reads; dashed, 5–20% support, 30 sites/fraction/control with peer controls as normal. (D) Recall versus false calls/library across four loci; zero is displayed at 0.05 on the logarithmic axis. (E) Simulated approximately 40-clone libraries, three/scenario; lines span libraries. (F) Read- and clone-weighted estimates versus corresponding truth; stems join each library's pair. An unmutated expanded clone contributes 60% of reads. (G) Estimate/truth by depth, means and ranges over three controls; solid, J_{H}4; dashed, N114. Combined density adds estimated uncalled single-read events. (H) Rejection rates across 300 experiments/condition, n=3 or 5/group; dotted line, 5%. Monte Carlo 95% half-width ≤5.6 points. Exact two-sided permutation p≥0.10 at n=3. BAQ, base-alignment quality. See Figures S4–S8 and S24; Tables S1 and S6.",
    "(A–E) Unique variants/10^{5} analyzed bases; circles, libraries; triangles, no calls. J_{H}4: three mice/genotype; N114: three subclones/arm, 500 reads/library. External Sμ/*Pim1*: 3–4 case and three null libraries, 2,000 reads each. Welch tests on log density compare case/null per caller, Holm-adjusted across five callers/panel. Zero-call tests use one event/library exposure. ns, not significant; \\*p<0.05; \\*\\*p<0.01; \\*\\*\\*p<0.001. (F) Geometric-mean case/null separation with the same zero convention; off-target Wt and double-knockout groups are separate (n=3 each). (G) Total calls in three controls/assay. (H) Held-out candidate tests with p<0.01 versus error margin; dotted line, 1%. (I) Delivered-read disposition. Panels use peer-control caller-comparison calibration; off-target Results use production-project calibration, so counts can differ (STAR Methods). See Figures S9–S11; Tables S1 and S3.",
    "Matched 529-bp data from three Wt and three null mice[@barbulescu2025]. (A) Within-platform densities; 26–40 Sanger clones/Wt mouse and 500 nanopore reads. Bars, geometric means; p, log-Welch. Nanopore zeros use one event. Densities are distinct platform-specific endpoints. (B) Pooled twelve-type shares ({n_sanger} Sanger/{n_ont} nanopore calls), without composition adjustment; dashed line, identity; r, descriptive correlation. (C) Four classes by mouse; bars, mean±SD; shapes identify mice. Platform×type testing uses all twelve types, matched by mouse. (D) Positional shares in 25-bp bins. (E) WRC/GYW use among G:C calls and WA/TW among A:T calls; bars, mouse means; points, mice; dashed lines, sequence opportunities. See Figures S13–S16; Tables S4 and S5.",
    "Libraries as Figure 5. (A) J_{H}4 variants/10^{5} bases in 50-bp bins, positions 90–1,230; thin lines, Wt mice; thick lines, genotype means, n=3 each. Shading, Sanger window. (B) Calls by S5F mutability quintile, normalized to each series' mean: Wt calls versus germline-excluded null mismatches. Targeting indices use permutation tests. (C) Junction-defined lineages among 500 reads/mouse: dark, largest; light, others; hatching, unassigned. Frequencies are read- and majority-rule clone-weighted. (D) Largest Wt3 lineage: square, reference; ticks, substitution classes; node sizes/counts, nested read support; distance, accumulated mutations. Forty-eight of 69 reads fit one tree; 21 are incompatible. This is a read reconstruction, not proven cellular ancestry. See Figure S17.",
    "(A) Illumina J_{H}4 frequencies (PRJNA973112); points, mice; triangles, no calls; bars, geometric means with one-event zero floors. Tests: log-Welch versus null; Mann–Whitney for *Ldha*/Wt. (B) Repair-deficient CH12F3 expressing AID variants or empty vector (PRJNA918596), n=4/arm; display as A. (C) Substitution classes and S5F indices; n counts calls. (D) Ramos indels at 20,000 reads/clone (PRJNA820976); n=6 Wt/HMCES-KO and 3 AID-KO; bars, means. (E) Paired mOrange repair arms across three experiments (PRJNA170048): circles, SHMtrace; squares, active-AID raw Q>20 mismatches; crosses, inactive-AID raw mismatches. Empty-arm inactive-AID libraries calibrate calling; 4–40 calls/library give coarse shares. See Table S1.",
]

# Preserve all six paragraphs and most of their wording; target the inaccurate openings
# and replace technical endings with specific affirmative take-home statements.
DISCUSSION_REPLACE = {
    0: [
        (
            "This study shows that the choice is no longer necessary, provided the background is measured rather than assumed.",
            "This study supports a cloning-free SHM workflow by estimating technical background from controls.",
        ),
        (
            "SHMtrace takes whatever a platform delivers, reads or chromatograms, and replaces the two assumptions that make low-frequency calling from single reads unreliable for SHM.",
            "SHMtrace connects delivered reads or chromatograms to mutation profiles and comparisons across biological replicates.",
        ),
        (
            "Instead of treating reads as replicates, it tests animals.",
            "Instead of treating reads as replicates, it tests independent biological samples.",
        ),
        (
            "were recovered at the false-discovery proportions the threshold targeted",
            "were recovered at false-discovery proportions near the default target",
        ),
    ],
    1: [
        (
            "learn from tumour–normal pairs where the variant sits in a large share of reads",
            "were developed for tumour–normal sequencing, with different training data and calling assumptions",
        ),
        ("already part of every well-designed SHM experiment", "where available"),
        (
            "The held-out controls then set the threshold empirically, in the spirit of target–decoy analysis in proteomics[@elias2007], rather than by assuming the distribution the p-values would need for a correction such as Benjamini–Hochberg[@benjamini1995].",
            "Its distinctive contribution is to connect control-informed mutation detection with interpretable profiles and comparisons across biological replicates in one reproducible workflow.",
        ),
    ],
    2: [
        (
            "a per-base error model fitted without matched controls has no way to reach this assay's library-specific background",
            "the configurations tested here did not capture this assay's background as effectively as matched-control calibration",
        ),
        ("varied between runs", "differed across the evaluated datasets"),
        (
            "and events still unresolved individually are estimated in aggregate rather than discarded.",
            "while events unresolved individually can be estimated in aggregate. Together, these outputs make the sensitivity–background tradeoff explicit and help investigators choose a calling regime appropriate to their experiment.",
        ),
    ],
    3: [
        (
            "Because the background is measured rather than modelled, the same engine reads all four platforms, and what each contributes becomes a question about the experiment rather than the software.",
            "A shared framework accommodates four platforms through distinct analyses: control-calibrated read-library substitutions, direct Sanger counting and PacBio junction characterization.",
        ),
        (
            "there is no read support to weigh against a background, and so nothing for a per-call test to work on, and the AID-null clones measure a background reported beside the counts, with how much of each animal's total it explains.",
            "the workflow therefore counts quality-filtered substitutions directly and reports the control-derived background separately.",
        ),
        (
            "so the cloning-free route reports the biology the classical assay reports.",
            "supporting the use of cloning-free long reads to recover established SHM patterns while adding whole-amplicon linkage.",
        ),
    ],
    4: [
        (
            "Applying the same engine to other laboratories' data turned out to test it harder than any internal benchmark.",
            "Reanalysis of public data tested SHMtrace against independent biological and counting expectations.",
        ),
        (
            "Deposits that include explicit clone-to-animal assignments and the original counting conventions, alongside sequences and chromatograms, would make such comparisons routine.",
            "By retaining assignments, counting conventions and quality flags, SHMtrace makes legacy sequence collections easier to revisit, compare and reuse.",
        ),
    ],
    5: [
        ("which is why it transferred without retuning", "supporting its application"),
        (
            "Molecule-level consensus is the natural next addition: UMI collapsing of long reads already reduced technical error 57-fold[@karst2021], and combining it with control-calibrated calling should push the measurable range below what single reads can resolve.",
            "Further validation on independently prepared runs and UMI-tagged SHM libraries can establish how far this framework extends. SHMtrace brings control-calibrated mutation detection, whole-amplicon profiles and comparisons across biological replicates into one reproducible workflow. By preserving mutation linkage and documenting analytical choices, it makes long-read SHM data more informative and supports transparent reanalysis of existing sequence collections. This provides a practical foundation for studying how AID targeting and DNA repair shape antibody diversification.",
        ),
    ],
}

RESOURCE_TEXT = "Nanopore reads, Sanger traces, clone records, project files and benchmark outputs will be deposited at {lab_data_doi}; public datasets and accessions are listed in Table S2 and the key resources table. SHMtrace v{shmtrace_version}, its desktop application, and analysis and figure code will be deposited at {lab_code_doi}. Both deposits will be public by publication. The archive will include checksums, sample sheets, read accounting, figure input tables, environment and tool versions, and commands linking analyses to reported outputs. Additional reanalysis information is available from the lead contact."

# These additions retain particular quantitative details outside the counted narrative.
# A migration ledger generated by evaluate_proposal.py identifies any other placeholders
# omitted from Results so they can be reviewed before implementation.
METHODS_ADDITIONS = [
    (
        "Sanger counting and provenance",
        "For the matched-platform comparison, {sanger_sheet} traces were assigned to the six mice, {sanger_excluded} excluded in the sample sheet and {sanger_submitted} processed; {sanger_failed} failed analysis and {sanger_clones} were retained (Figure S13). This distinguishes the complete {sanger_archive}-trace archive from the matched comparison.",
    ),
    (
        "Caller operating ranges",
        "On real J_{H}4 libraries, ClairS reported {cs_real_pass} mutation and DeepSomatic at most {ds_real_pass}, compared with 221, 344 and 73 SHMtrace calls; N114 totals were {cs_real_n114_pass}, at most {ds_real_n114_pass}, and 70, 49 and 173, respectively. ClairS rejected both Wt3 candidates supported by ≥20% of reads. Among {af_n} retained calls, {af_lt1}% were in <1%, {af_lt5}% in <5% and {af_lt20}% in <20% of reads. These results describe the tested settings and datasets rather than universal caller limits.",
    ),
    (
        "Sanger comparison details",
        "In matched mice, G:C hotspot shares were {hot_gc_s}% for Sanger and {hot_gc_o}% for nanopore versus {hot_gc_exp}% of eligible positions; A:T shares were {hot_at_s}% and {hot_at_o}% versus {hot_at_exp}%. No equivalence conclusion follows from the non-significant interaction at n=3. In the independent SW_{HEL} comparison, the analysis window was {sw_bp} bp; per-class p-values were {sw_at_tv_p}, {sw_gc_tv_p}, {sw_gc_ts_p} and {sw_at_ts_p}. The published total-mutation p=0.07 used one-way ANOVA; the reanalysis p={sw_tot_p} used per-animal log-rate Welch. The deposit's {sw_wt_hosts} hosts differ from the reported 13; the potentially pooled host was retained and excluded in sensitivity analysis. A:T mutation shares fell from {sw_at_wt}% in Wt to {sw_at_msh2}% in *Msh2* mutants.",
    ),
    (
        "PacBio junction comparison details",
        "Mean-microhomology comparisons against Wt gave p={pb_p_shld} for *Shld1*, {pb_p_xlf} for *Xlf* and {pb_p_dko} for the double mutant. The ≥4-nt share rose from {pb_4_wt}% to {pb_4_dko}% (p={pb_p4_dko}), and blunt junctions fell from {pb_bl_wt}% to {pb_bl_dko}% (p={pb_pbl_dko}). *Shld1* insertion-bearing junctions were {pb_ins_shld}% versus {pb_ins_wt}% in Wt (p={pb_p_ins_shld}). Reads per unique junction ranged {pacbio_rpj}; inferred structures do not measure switching efficiency.",
    ),
]

# Select a moderate reduction: preserve the author's dilution paragraph, the full
# matched-mouse comparison, most benchmark reporting and the main biological results.
# The alternatives above remain available for later editorial choices.
COMPACT_RESULTS = RESULTS.copy()
_baseline = json.loads((Path(__file__).parent / "source_blocks.json").read_text())
_shorten = {1, 2, 4, 5, 8, 16, 20, 24, 25, 29, 31, 39, 41, 45}
RESULTS = {i: t for i, t in COMPACT_RESULTS.items() if i in _shorten}
RESULTS[7] = _baseline["RESULTS"][7][1].replace(
    "and {sp_jh4_ge3}% of mutations were supported by at least three reads",
    "{sp_jh4_ge3}% among events supported by at least three reads",
)
# Include "and" after the preceding comma so the two recall denominators are explicit.
RESULTS[7] = RESULTS[7].replace("overall, {sp_jh4_ge3}%", "overall and {sp_jh4_ge3}%")
RESULTS[9] = _baseline["RESULTS"][9][1].replace(
    "{sp_su_lps_gen_fdr} and {sp_n114_gen_fdr}", "{sp_su_lps_gen_fdr} and {sp_pim1_gen_fdr}"
)
RESULTS[19] = _baseline["RESULTS"][19][1].replace(
    "and three loci from an ", "and two loci in three assay contexts from an "
)
RESULTS[23] = _baseline["RESULTS"][23][1].replace(
    "statistically indistinguishable from ", "with no detected difference from "
)
RESULTS[30] = _baseline["RESULTS"][30][1].replace(
    "chromatograms[@barbulescu2024]", "chromatograms[@barbulescu2025]"
)
RESULTS[27] = _baseline["RESULTS"][27][1].replace(
    "processed 182 Sanger chromatograms", "processed {sanger_submitted} Sanger chromatograms"
)
LEGENDS[0] += (
    " The PacBio evaluation concerns junction structure, without validating control-calibrated mutation calls."
)
LEGENDS[6] = LEGENDS[6].replace(
    "Repair-deficient CH12F3", "*Aicda*^{−/−}*Ung*^{−/−}*Msh2*^{−/−} CH12F3"
)
DISCUSSION_REPLACE[0].append(
    (
        "Instead of trusting quality scores, it measures the background",
        "Rather than relying on quality scores alone, it estimates the background",
    )
)
DISCUSSION_REPLACE[0].append(
    (
        "of simulated experiments without a true difference.",
        "of simulated experiments without a true difference. Together, these results support a practical route from deeply sequenced amplicons to interpretable SHM measurements across biological replicates.",
    )
)
DISCUSSION_REPLACE[1].append(
    (
        "Neither class measures what a particular SHM library does wrong.",
        "The evaluated configurations did not use AID-null calibration in the manner tested here.",
    )
)

METHODS_ADDITIONS += [
    (
        "Two-candidate example",
        "At Sμ {ana_case}, position {ana_t_pos} ({ana_t_ctx}) had {ana_t_reads} supporting reads, median Q{ana_t_q}; position {ana_e_pos} ({ana_e_ctx}) had {ana_e_reads}, median Q{ana_e_q}. Both had support on both strands. {ana_e_ko} carried the latter in {ana_e_ko_reads} reads at median Q{ana_e_q}. Geneious called it in {ana_e_gen_ctl} controls and all {ana_e_gen_case} cases. These examples distinguish control-supported calling decisions, without independent truth validation of each native candidate.",
    ),
    (
        "Real-library caller details",
        "Figure 4's unique-density separation ranged {gen_lo}–{gen_hi}-fold for Geneious and {lf_lo}–{lf_hi}-fold across LoFreq settings. J_{H}4 Wt calls were {jh4_wt2_rec}, 344 and 73 versus 0, 2 and 0 in nulls; Geneious gave 203, 337 and 73 versus 2, 6 and 5. Default LoFreq gave 0, 11 and 0; BAQ-off gave {jh4_wt2_lfnb} in Wt2. External Sμ Wt counts were {su_wt_calls} for SHMtrace and {su_gen_wt_calls} for Geneious, with {su_ctl_calls} versus {su_gen_ctl_calls} per null. The displayed Wt/null raster libraries had {rast_wt_other}/{rast_ko_other} other Q≥20 mismatches after germline exclusion; Wt added {rast_wt_called} mismatches at called sites. These are individual mismatch observations, distinct from unique variant counts.",
    ),
    (
        "Maea comparison details",
        "SHMtrace processed {c2_traces} chromatograms in {c2_seconds} s after manual assignment curation. G:C transitions fell from {c2_ts_wt}% to {c2_ts_mk}% (same published shares; p={c2_ts_p} versus {c2_ts_pub_p} published), while transversions rose from {c2_tv_wt}% to {c2_tv_mk}% versus {c2_tv_pub_wt}% to {c2_tv_pub_mk}% published (p={c2_tv_p} versus {c2_tv_pub_p}). Mouse rankings correlated at ρ={c2_rho}. Controls had {c2_ko_calls} mutations in {c2_ko_n} clones (p={c2_wt_ko_p} versus Wt and {c2_mk_ko_p} versus knockout). Shared trace assignments and {c2_dup_flags} identical mutated-clone pairs remained as filed and were flagged.",
    ),
    (
        "Ramos control sensitivity",
        "The noisy null library carried {hmces_ko3_calls} substitutions, {hmces_ko3_gc}% at G:C, while the other two had none. Retaining it as a decoy suppressed all sample substitution calls. Excluding it from decoys, while retaining it in model fitting, yielded {hmces_gc_range}% G:C calls; treating it as a sample instead yielded {hmk3_gc}%. The rule was developed from this dataset, so this is a sensitivity analysis rather than independent validation of that rule.",
    ),
    (
        "mOrange raw-background comparison",
        "The source study's Q>20 every-mismatch counts gave inactive-AID background of {mo_raw_lo}–{mo_raw_hi} per 10^{3} bases, approximately two-fold active/inactive separation in the empty arm, and {mo_dead_ts}% G:C transitions in inactive backgrounds. SHMtrace's inactive Ugi and UNG arms retained {mo_ugi_dead} and {mo_ung_dead} calls, respectively. Calling across contexts used inactive-AID empty-arm controls. The paired effect and interval use one transition-share value per experiment; the pooled {mo_ugi_ts}%/{mo_ung_ts}% shares have a different weighting and need not subtract to the paired mean difference.",
    ),
    (
        "Somatic caller spike-in false calls",
        "Across the four planted datasets, permissive DeepSomatic produced at most {ds_fp} false calls per library under the evaluated settings.",
    ),
]

METHODS_REPLACEMENTS = [
    (
        "Four datasets in which a deaminase is directed deliberately",
        "Two datasets in which a deaminase is directed deliberately",
    ),
]
