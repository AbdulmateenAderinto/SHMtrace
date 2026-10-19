"""Text of the Cell Reports Methods article (SHMtrace). Markup: see docbuild.py.

Numbers for the public Illumina reanalyses are filled from out/public/*/project.json by
build_crm.py (keys in braces, e.g. {ldha_fold}); everything else is written out and traced to
the analysis outputs listed in STAR Methods.
"""

TITLE = ("SHMtrace calls somatic hypermutation and other AID-driven mutation with a measured "
         "false-discovery rate")
RUNNING = "Control-calibrated SHM profiling"
AUTHORS = "Abdulmateen Aderinto^{1,3} and Alberto Martin^{1,2,*}"
AFFIL = ["^{1}Department of Immunology, University of Toronto, Toronto, ON M5S 1A8, Canada",
         "^{2}Lead contact", "^{3}Technical contact",
         "\\*Correspondence: alberto.martin@utoronto.ca"]

SUMMARY = (
    "Somatic hypermutation (SHM) at genomic reporter loci is often measured by cloning amplicons and sequencing "
    "colonies, because single reads err at rates close to the mutation rate. We describe SHMtrace, which tests "
    "candidate substitutions in read libraries against errors measured in matched activation-induced cytidine "
    "deaminase (AID)-deficient controls from the same run and compares biological replicates. "
    "Sanger mismatches are counted at a quality threshold, with background reported separately. We evaluated "
    "SHMtrace using planted mutations, simulated lineages, and fourteen loci from thirteen laboratories "
    "across four platforms. It retained {ctl_rec} AID-null calls per assay and separated AID-active from AID-null libraries "
    "{sep_lo}–{sep_hi}-fold in four of five nanopore assays, compared with {lf_lo}–{lf_hi}-fold and "
    "{gen_lo}–{gen_hi}-fold for quality-based and uncorrected tests. Nanopore and Sanger sequencing of the same "
    "mice gave similar pooled spectra. Re-analysis of three published Sanger studies recovered principal "
    "comparisons while identifying unresolved differences.")

KEYWORDS = ("somatic hypermutation, activation-induced cytidine deaminase, nanopore sequencing, variant calling, "
            "negative controls, false-discovery rate, pseudoreplication, immunoglobulin, mutation trees, Sanger "
            "sequencing")
MOTIVATION = (
    "At genomic reporter loci such as the J_{H}4 intron, somatic hypermutation is still measured mostly by "
    "cloning PCR products and sequencing colonies one at a time, which is slow and samples few molecules. "
    "Our aim was to analyze a defined genomic reporter while preserving linkage between mutations across the "
    "amplicon. Long nanopore reads could replace cloning for this assay, but their "
    "per-read error rate is close to the mutation rate of interest, and general-purpose variant callers judge "
    "each base by its quality score rather than by the errors a given library actually makes. We therefore "
    "built SHMtrace to calibrate each call against matched AID-null libraries sequenced in the same run, to "
    "treat animals rather than reads as replicates, and to make the whole workflow, from reads or chromatograms "
    "to group statistics, reproducible. The Sanger route also supports re-analysis of deposited traces once clone identities and sample assignments "
    "have been curated, making counting choices and remaining discrepancies explicit.")
HIGHLIGHTS = [
    # four, as Cell Press guidance asks; each at most 85 characters
    "Candidate mutations in reads are tested against errors of matched AID-null libraries",
    "AID-null controls retain {ctl_rec} calls per assay while true mutations are kept",
    "Nanopore gives spectra similar to Sanger, plus whole-amplicon maps and mutation trees",
    "Re-running three published Sanger studies recovers their principal comparisons",
]
ETOC = ("Aderinto and Martin present SHMtrace for somatic hypermutation analysis of amplicon reads and Sanger traces. Read-library candidates are tested against matched control errors; Sanger mismatches are counted with background reported separately. Benchmarks assess recall, false calls, biological comparisons and reproducibility across four platforms.")

INTRO = [
    "Somatic hypermutation (SHM) diversifies immunoglobulin (Ig) genes in germinal-center B cells. Activation-induced "
    "cytidine deaminase (AID) initiates this process by deaminating cytosines in Ig DNA. Replication over the "
    "resulting uracils, together with uracil excision and mismatch repair, generates substitutions at G:C and A:T base "
    "pairs[@muramatsu2000,@dinoia2007,@methot2017]. Mutation frequency, position and spectrum therefore provide "
    "complementary readouts of AID activity and lesion processing4–7. The mouse J_{H}4 downstream intron is widely "
    "used as a genomic reporter because it undergoes SHM without directly contributing to antigen binding8,9.",
    "This assay has changed little in decades. Amplicons are cloned into bacteria and sequenced colony "
    "by colony, which takes days of hands-on work and yields tens of sequences per animal[@sible2021]. "
    "Our goal is to analyze this genomic reporter while retaining linkage between mutations across the "
    "amplicon; short reads may not span the complete region[@martin2023]. Long reads from Oxford Nanopore Technologies (ONT) now reach "
    "modal accuracies above 99% with R10.4.1 flow cells[@zhang2023,@sereika2022,@sanderson2024] and have been "
    "used to sequence Ig genes where consensus sequencing can absorb errors[@minervini2021,@subassatish2022]. "
    "Calling rare substitutions from single reads across a polyclonal population is harder. A wild-type J_{H}4 "
    "intron carries roughly one mutation per 100 bases, whereas the residual error of a single filtered read is "
    "within one or two orders of magnitude of that rate and is concentrated at particular positions and "
    "sequence contexts[@delahaye2021].",
    "Two common analysis practices make this problem worse. First, general-purpose low-frequency variant callers "
    "judge each base by its own quality score[@wilm2012], and the simple rules used in commercial software test "
    "candidates one at a time without correcting for the number of tests[@kearse2012]. In the configurations "
    "evaluated here, neither measures the errors a particular library makes, so mismatches that recur at one "
    "position, from PCR, sequence context or a co-amplified paralog, can be counted as SHM mutations. "
    "Second, reads from one animal are often "
    "pooled and treated as independent replicates, which overstates certainty because reads from the same "
    "animal share clones and PCR history[@lazic2010]. SHM studies already include the control that addresses "
    "the first: libraries from AID-null animals or cells, sequenced alongside the samples. The second is a "
    "matter of design, one value per animal. "
    "Existing SHM software summarizes mutations once called[@maccarthy2009,@gupta2015,@yaari2013] "
    "but does not use these controls to decide which differences are mutations.",
    "Here we describe SHMtrace, which calibrates each call against matched AID-null libraries from the "
    "same run, compares groups with the animal as the unit, and runs from "
    "delivered reads or chromatograms to group statistics. We validated it with planted mutations, "
    "simulated lineages, fourteen loci from thirteen laboratories across four platforms, Sanger "
    "sequencing of the same mice, and deposited clones from three published Sanger studies, and compared it against LoFreq[@wilm2012], the long-read somatic callers "
    "ClairS[@zheng2026] and DeepSomatic[@park2026], the germline callers Clair3[@zheng2022] and "
    "DeepVariant[@poplin2018], and the uncorrected Geneious test[@kearse2012].",
]

RESULTS = [
    ("h", "SHMtrace calibrates substitutions against AID-null controls"),
    ("p", "SHMtrace starts from the files that a sequencing provider or core facility returns (Figures 1 and S1). Nanopore "
          "and Illumina reads are first filtered by length and expected error[@edgar2015] and aligned with "
          "minimap2[@li2018] against the target together with competing references. Reads that align better to a "
          "paralog, are chimeric, or carry ten or more paralog-specific bases are removed (Figure 1B), and the "
          "same number of usable reads is then drawn at random from every library. Sanger chromatograms are "
          "assigned to animals through a sample sheet and trimmed to their high-quality core, and reads tagged "
          "with unique molecular identifiers (UMIs) are collapsed to one consensus per molecule."),
    ("p", "The central step is the variant test (Figures 1A and 1C). For each substitution type, strand, and base "
          "quality, SHMtrace measures how often the AID-null libraries of the same run show that change, and asks "
          "whether a candidate is supported by more reads, at the qualities observed, than this background "
          "explains. Because quality scores overstate certainty at the highest qualities (Figure S2), the control "
          "rates are multiplied by a margin chosen by leaving one control out at a time (STAR Methods). The "
          "held-out controls then serve as decoys, in an approach inspired by target–decoy analysis[@elias2007]: "
          "SHMtrace keeps the most permissive cutoff at which estimated false calls are at most 5% of calls. This is "
          "an empirical estimate targeting 5%, not a guarantee: it relies on a few controls representing the "
          "sample's background, and the same controls also choose the margin. Candidates arise only where a read "
          "already differs, so the set tested is data-dependent and its null p values need not be superuniform, "
          "as Benjamini–Hochberg[@benjamini1995] requires; SHMtrace counts false calls in controls "
          "selected the same way. Alternative alleles carried by at least half the reads at a position are "
          "set aside as likely germline or parental differences. This is a heuristic, not a test: a mutation "
          "shared by a lineage that dominates a library crosses the same threshold, so every site it removes "
          "is reported with its read count and frequency."),
    ("p", "The read-library test can retain a candidate supported by a single read; repeated sequencing of the "
          "same original molecule is not required. SHMtrace uses a separate Sanger counting route (Figure 1C): "
          "every substitution at base quality ≥ 20 is counted, and the AID-null clones' mismatch rate gives "
          "the expected background count, reported separately without subtraction. The germline "
          "filter, window, automatic checks and animal-level comparisons are shared. The Sanger results below "
          "therefore test trace reading, counting, bookkeeping and group comparisons, not the "
          "control-calibrated test."),
    ("p", "Two C-to-T candidates in an external Sμ wild-type library illustrate the test ({ana_case}; Figure 2). "
          "Position {ana_t_pos}, within a WRC hotspot ({ana_t_ctx}), was supported by 20/2,000 reads; position "
          "{ana_e_pos} ({ana_e_ctx}) by 15/1,941. Both occurred on both strands at similar median qualities ({ana_t_q} "
          "and {ana_e_q}). No control supported position {ana_t_pos} at Q20, whereas {ana_e_ko} supported position "
          "{ana_e_pos} in 14/1,879 reads at median Q33. SHMtrace called the first (p = {ana_t_p}) but rejected the "
          "second (p = 1). The Geneious rule called the latter in {ana_e_gen_ctl} controls and all {ana_e_gen_case} "
          "wild-type libraries. Decoy thresholding retained {ana_kept} of {ana_cand} candidates at p ≤ {ana_thr} "
          "(Figure 2D)."),
    ("p", "Group comparisons use one value per biological replicate, with negative-binomial parametric-bootstrap, "
          "log-rate Welch and permutation tests (Figure 1D). Outputs include substitution spectra, positional "
          "clustering, motif shares, an S5F-based 5-mer targeting index23, an aggregate estimate of uncalled "
          "single-read mutations, junction-defined clones, read-level mutation trees and power estimates. Command-line "
          "and desktop interfaces use the same engine."),

    ("h", "Planted mutations and simulated lineages quantify recall and false calls"),
    ("p", "We planted 60 AID-weighted substitutions into each held-out control library, preserving base qualities and "
          "assigning 1–20 supporting reads per site, with three repetitions per control (Figure 3). On J_{H}4, "
          "SHMtrace recalled {sp_jh4_recall}% overall and {sp_jh4_ge3}% of mutations supported by at least three "
          "reads, with {sp_jh4_fp} false calls per library (observed FDR {sp_jh4_fdr}). Geneious recalled more "
          "single-read mutations ({sp_jh4_gen_1}% versus {sp_jh4_1}%) but produced {sp_jh4_gen_fp} false calls per "
          "library (FDR {sp_jh4_gen_fdr}). Default LoFreq recalled none; disabling its short-read-oriented "
          "base-alignment quality correction yielded {sp_jh4_lfnb}% recall, restricted to mutations with at least five "
          "supporting reads. Permissive settings increased recall to {sp_jh4_lfp}% with {sp_jh4_lfp_fp} false calls "
          "per library (Figures 3A and 3C)."),
    ("p", "ClairS[@zheng2026] and DeepSomatic[@park2026] were supplied with the other AID-null libraries as matched "
          "normals. Neither recovered planted mutations at default settings. Across 2,160 substitutions on J_{H}4, "
          "N114, Sμ and *Pim1*, permissive ClairS recovered none; permissive DeepSomatic recalled 19%, 18%, 5% and 3%, "
          "respectively, with at most {ds_fp} false calls per library. Recall increased at higher allele fractions, "
          "approaching complete recovery at 20% on J_{H}4 (Figure 3H). On real J_{H}4 libraries, ClairS reported "
          "{cs_real_pass} mutation and DeepSomatic at most {ds_real_pass} across three mice, compared with 221, 344 "
          "and 73 SHMtrace calls. Corresponding N114 totals were {cs_real_n114_pass} and at most {ds_real_n114_pass}, "
          "compared with 70, 49 and 173. ClairS rejected both Wt3 candidates supported by at least 20% of reads. Most "
          "SHMtrace-retained calls had much lower fractions: among {af_n} calls across four loci, {af_lt1}% occurred "
          "in <1% of reads, {af_lt5}% in <5% and {af_lt20}% in <20% (Figure 2E). Germline callers Clair3 and "
          "DeepVariant recovered none of the 2,160 planted substitutions, consistent with their intended operating "
          "range (Table S2)."),
    ("p", "On N114, SHMtrace recalled {sp_n114_recall}% of planted substitutions, including {sp_n114_1}% supported by "
          "one read, with {sp_n114_fp} false calls per library (FDR {sp_n114_fdr}; Figure 3B). Geneious recalled "
          "{sp_n114_gen_recall}% with {sp_n114_gen_fp} false calls (FDR {sp_n114_gen_fdr}); LoFreq without "
          "base-alignment quality recalled {sp_n114_lfnb}%. Of **two** C3 control libraries, the one with more usable "
          "reads was used primarily. Substituting the noisier library raised SHMtrace FDR to {c31_fdr}, or "
          "{c31_held_fdr} when that library was tested as the case. The **two** libraries shared no calls, consistent "
          "with library-specific noise, although sampling also limits overlap. Benjamini–Hochberg gave FDR "
          "{sp_n114_fdr} with the primary controls and {c31_bh_fdr} with the alternative C3 library, compared with "
          "SHMtrace's {c31_fdr}. On external Sμ and *Pim1* controls at 2,000 reads, SHMtrace recalled "
          "{sp_su_lps_recall}% and {sp_pim1_recall}% at FDR {sp_su_lps_fdr} and {sp_pim1_fdr}. Geneious recalled "
          "{sp_su_lps_gen_recall}% and {sp_pim1_gen_recall}% at FDR {sp_su_lps_gen_fdr} and {sp_n114_gen_fdr}; LoFreq "
          "recalled {sp_ext_lf_lo}–{sp_ext_lf_hi}%, with permissive FDR {sp_ext_lfp_fdr} (Figure 3C)."),
    ("p", "Indel validation used ten planted insertions or deletions of 1–40 bp per library, supported by 1–20 reads, "
          "with {ind_draws} draws per control and error models fitted to the other controls. Among events supported by "
          "at least three reads, recall was {ind_jh4_ge3}% on J_{H}4 and {ind_n114_ge3}% on N114, with {ind_jh4_fp} "
          "and zero false calls per library. Events supported by one or **two** reads were excluded by the three-read "
          "requirement. J_{H}4 single-base homopolymer deletions were least reliably recovered ({ind_jh4_hp}%); recall "
          "was {ind_jh4_ins2}% for 2-bp insertions and {ind_jh4_long}% for deletions ≥10 bp (Figure S3D)."),
    ("p", "We next simulated approximately 40 power-law-sized B-cell clones per J_{H}4 library, with S5F-weighted "
          "trunk, branch and private mutations (Figure 3D). SHMtrace recalled {lin_1_rec}% of single-read mutations "
          "versus {lin_1_gen}% for Geneious, with FDR {lin_fdr_rec} versus {lin_fdr_gen} and {lin_fp_rec} versus "
          "{lin_fp_gen} false calls per library. In null libraries, false calls averaged {lin_null_rec} versus "
          "{lin_null_gen}. Junction clustering recovered clone assignments with adjusted Rand indices32 of {lin_ari}. "
          "When an unmutated clone contributed 60% of reads, read-weighted frequency was {lin_rw_est} versus a true "
          "{lin_rw_true} per 10^{3} bp, while the majority-rule clone-weighted measure recovered {lin_cw_pct}% of its "
          "corresponding true burden ({lin_cw_est} versus {lin_cw_true}; Figure 3E)."),
    ("p", "The dominant unmutated clone could not trigger the ≥50% allele filter. Repeating the simulation with a "
          "mutated dominant clone showed no loss at 40% abundance, but all mutations crossing the threshold at 60% or "
          "80% were removed. Read-weighted and clone-weighted estimates fell by as much as {mask80_rw}% and "
          "{mask80_cw}%, respectively (STAR Methods). Excluded sites are reported for inspection."),
    ("p", "Read-weighted frequency was within {q_rw_dev}% of truth at 500 reads for J_{H}4 and N114 (Figure 3F). "
          "J_{H}4 unique-variant density reached {q_ud_500} of truth at 500 reads and {q_ud_1000} at 1,000 reads; N114 "
          "reached {q_ud_n114_1000} at 1,000. This measure depends on depth because distinct calls saturate while "
          "analyzed base observations increase, and rare-event detection also changes. Adding the aggregate "
          "single-read estimate brought density to {q_cd} of truth across scenarios. When 25, 50 or 100 single-read "
          "mutations were planted, the aggregate estimate increased by {sr_rec}% of the planted number on average "
          "(Figure S3). Unique densities should therefore be compared at matched depth."),
    ("p", "Group-test calibration was assessed in 300 simulations per condition, with between-animal coefficient of "
          "variation 1 (Figure 3G). Under no true difference, pooled-read tests returned p < 0.05 in 91% of 3-versus-3 "
          "and 94% of 5-versus-5 experiments. Negative-binomial rejection rates were 5.0% and 5.7%. With three animals "
          "per group, that test detected two-fold and four-fold changes in 12% and 30% of simulations, illustrating "
          "the limited power at this sample size."),

    ("h", "Matched controls set the background in real libraries"),
    ("p", "Real-library comparisons used identical prepared reads for every caller (Figure 4; Table S6). Assays "
          "comprised J_{H}4 from three wild-type (Wt) and three *Aicda*^{−/−} mice, AID-transduced and empty-vector "
          "N114 subclones, and external Sμ libraries from LPS-activated primary B cells and CH12 cells, together with "
          "the AID off-target *Pim1*."),
    ("p", "SHMtrace retained {ctl_rec} calls across the three controls per assay, versus {ctl_gen} for Geneious "
          "(Figures 4A–4E and 4G). Geometric-mean unique-density separation, with zero-call controls assigned "
          "one-event values, was {sep_jh4}-fold for J_{H}4, {sep_n114}-fold for N114, {sep_su_lps}-fold for "
          "primary-cell Sμ and {sep_su_ch12}-fold for CH12 Sμ, compared with 16–62-fold for Geneious and 2–40-fold for "
          "LoFreq configurations (Figure 4F). J_{H}4 Wt counts were {jh4_wt2_rec}, 344 and 73 versus 0, 2 and 0 in "
          "controls; Geneious yielded 203, 337 and 73 versus 2, 6 and 5. Default LoFreq called 0, 11 and 0 Wt "
          "variants. Disabling base-alignment quality yielded {jh4_wt2_lfnb} in Wt2 versus {jh4_wt2_rec} for SHMtrace. "
          "At 500 reads, J_{H}4 and N114 had reached {depth_jh4_500}% and {depth_n114_500}% of calls observed at the "
          "largest tested depth (Figure S5), without establishing saturation. External Sμ Wt counts were {su_wt_calls} "
          "for SHMtrace and {su_gen_wt_calls} for Geneious, with {su_ctl_calls} versus {su_gen_ctl_calls} per control. "
          "Control mismatches were position-dependent (Figures S6A and S6B). After germline exclusion, representative "
          "Wt and null J_{H}4 libraries contained {rast_wt_other} and {rast_ko_other} other high-quality mismatches, "
          "while the Wt library added {rast_wt_called} at called sites, often shared within clones (Figures S6C and "
          "S6D)."),
    ("p", "*Pim1* approached the detection limit. At 2,000 reads, separation was {pim1_rec_fold}-fold for SHMtrace "
          "with four control calls and {pim1_gen_fold}-fold for Geneious with eleven. Using all available reads "
          "(3,238–3,447 per library) yielded 3–10 versus 0–3 SHMtrace calls in AID-competent and null libraries, and "
          "5–7 versus 0–4 for Geneious. Signal remained close to background (Table S2)."),
    ("p", "CH12 off-targets OT1 and OT2 tested a genetic perturbation expected to reduce AID targeting33. OT1 Wt "
          "libraries carried 79, 56 and 57 calls versus 0, 2 and 1 in *Aicda*^{−/−} controls, with densities of "
          "{ot1_wt_dens} versus {ot1_ko_dens} per 10^{5} bases (Welch p = {ot1_p_dens} on log density; {ot1_p} on "
          "frequency). OT2 counts were 28, 23 and 27 versus 0, 0 and 1, with densities of {ot2_wt_dens} versus "
          "{ot2_ko_dens} (p = {ot2_p_dens} and {ot2_p}). G:C substitutions accounted for {ot1_gc}% of {ot1_n} OT1 "
          "calls and {ot2_gc}% of {ot2_n} OT2 calls. MLLT1/MLLT3 double knockout reduced counts to 1, 1 and 0 at OT1 "
          "and zero throughout OT2, statistically indistinguishable from *Aicda*^{−/−} controls (p = {ot1_p_dko_ko} "
          "and {ot2_p_dko_ko}). No fold is reported because a zero-count convention would determine its magnitude. At "
          "OT1, {ot1_ctl_drop} was excluded from threshold-setting decoys, leaving **two** controls (Limitations)."),
    ("p", "Without margin inflation, {loco_jh4_face}% of held-out J_{H}4 and {loco_n114_face}% of N114 candidate tests "
          "reached p < 0.01. Selected margins of 5× for J_{H}4, 10× for full-length J_{H}4 and 10× for N114 reduced "
          "these proportions to {loco_jh4_final}%, {loco_full_final}% and {loco_n114_final}% (Figure 4H). Read "
          "preparation removed 44–79% of delivered J_{H}4 reads, 67% overall, principally through expected-error "
          "filtering (29%) and competitive alignment (25%; Figure 4I). Without the paralog-allele filter, full-length "
          "paralog reads contributed 35–54 false calls per simulated library because they could align better to the "
          "longer target than to the shorter competing reference."),

    ("h", "Nanopore and Sanger profiles agree in matched mice"),
    ("p", "Of 447 archived chromatograms34, 184 were assigned to the six nanopore-sequenced mice and 134 yielded "
          "analyzed clones (Figures 5 and S7). Across the shared 529-bp window, both platforms separated Wt from "
          "*Aicda*^{−/−} mice: Sanger, {sanger_fold}-fold (p = 8.8 × 10−4); nanopore, {ont_fold}-fold with zero-call "
          "controls assigned one event (p = {ont_p}; log-density Welch tests, n = 3 per group). Pooled spectra "
          "correlated across twelve substitution types (descriptive r = {r_spec}; slope = {slope}), as did positional "
          "profiles in 25-bp bins (r = {r_pos}). WRC/GYW shares among G:C mutations were {hot_gc_s}% for Sanger and "
          "{hot_gc_o}% for nanopore; WA/TW shares among A:T mutations were {hot_at_s}% and {hot_at_o}%, versus window "
          "opportunities of {hot_gc_exp}% and {hot_at_exp}%. The platform-by-type interaction was not detected (p = "
          "{interaction_p}), which does not establish equivalence with three mice. Absolute unique densities are not "
          "directly comparable because Sanger samples individual cloned products whereas nanopore counts distinct "
          "variants among 500 reads."),
    ("p", "Chromatogram reanalysis identified misfiled clones and renamed copies of traces. SHMtrace therefore flags "
          "copied files, unexpectedly shared mutated clones, mutated controls and label inconsistencies without "
          "automatically correcting assignments. In 40 simulated projects, all copied or retyped traces, contaminated "
          "controls and misfiled clones retaining their labels were flagged, including all 125 clones in misfiled "
          "lineages; no flags occurred among 6,018 correctly filed clones (Figure S7). When labels were also changed, "
          "detection fell to 22/40 projects and depended on shared mutations."),

    ("h", "Sanger reanalysis recovers published comparisons and identifies discrepancies"),
    ("p", "A published J_{H}4 study compared four Wt and four *Maea*^{fl/fl}*Mb1*^{cre} mice with *Aicda*^{−/−} "
          "controls using SHMTool[@barbulescu2024,@maccarthy2009]. After reconstruction of clone assignments, SHMtrace "
          "processed {c2_traces} chromatograms in {c2_seconds} s; curation took substantially longer (Figure S8; STAR "
          "Methods). Mutation frequency was {c2_fold}-fold lower in *Maea*-deficient mice versus {c2_pub_fold}-fold in "
          "the publication (one-tailed t tests, p = {c2_p} and {c2_pub_p}). G:C transitions fell from {c2_ts_wt}% to "
          "{c2_ts_mk}% (published, {c2_ts_pub_wt}% to {c2_ts_pub_mk}%; p = {c2_ts_p} and {c2_ts_pub_p}), and "
          "transversions rose from {c2_tv_wt}% to {c2_tv_mk}% ({c2_tv_pub_wt}% to {c2_tv_pub_mk}%; p = {c2_tv_p} and "
          "{c2_tv_pub_p}). Spectra correlated at r = {c2_spec_r}; mouse rankings at ρ = {c2_rho}. Under the published "
          "clone-count × 577-bp convention, frequencies averaged {c2_ratio}-fold the published values. Controls "
          "carried {c2_ko_calls} mutations in {c2_ko_n} clones (p = {c2_wt_ko_p} and {c2_mk_ko_p} versus Wt and "
          "knockout). Four traces assigned to {c2_ko_calls} mice and {c2_dup_flags} identical mutated-clone pairs were "
          "retained as filed and flagged."),

    ("p", "Reanalysis of {mk_traces} Mkln1-study chromatograms34 found an {mk_fold}-fold frequency reduction versus "
          "the published {mk_pub_fold}-fold (p = {mk_p} and {mk_pub_p}). G:C transversions increased from {mk_tv_wt}% "
          "to {mk_tv_ko}% versus {mk_tv_pub_wt}% to {mk_tv_pub_ko}% in the published tables (p = {mk_tv_p} and "
          "{mk_tv_pub_p}). Spectrum correlations were r = {mk_spec_r_wt} for Wt and {mk_spec_r_ko} for knockout; only "
          "one published knockout table could be matched to its folder (Figure S9; Table S4). Excluding {mk_paralog} "
          "mostly unmutated paralog clones changed denominators. Applying the original counting convention brought "
          "knockout frequencies close to the published estimates, although Wt discrepancies remained. Three lineages "
          "were filed under **two** mice each; matched nanopore data supported reassignment of **two** and could not "
          "place the third. Applying prespecified resolution rules preserved the result ({mk_ev_fold}-fold; p = "
          "{mk_ev_p})."),

    ("p", "An independent SW_{HEL} V_{H}10 study tested *Samhd1*, *Msh2* and UGI perturbations36. Analysis of its "
          "{sw_clones} GenBank clones used the majority-consensus reference and a {sw_bp}-bp shared window, close to "
          "the published 560 bp (Figure S10). SHMtrace recovered {sw_muts} of {sw_pub_muts} tabulated mutations across "
          "five arms and twelve substitution types (descriptive r = {sw_r}); mutations per sequence agreed to **two** "
          "decimal places. Samhd1-associated changes agreed in direction and significance for all four substitution "
          "classes: A:T and G:C transversions decreased (p = 9.7 × 10−11 and {sw_gc_tv_p}), G:C transitions increased "
          "(p = {sw_gc_ts_p}), and A:T transitions showed no detectable change ({sw_at_ts_p}). Total mutation "
          "decreased significantly under a per-animal log-frequency Welch test (p = {sw_tot_p}), unlike the published "
          "five-arm ANOVA (p = 0.07); these are different comparisons. The deposit identifies twelve Wt hosts rather "
          "than thirteen, with one unusually large clone collection; excluding that host preserved the conclusions. "
          "*Msh2* loss reduced the A:T share from {sw_at_wt}% to {sw_at_msh2}%. This sequence-only dataset tests "
          "counting and group analysis, not chromatogram processing or AID-null calibration."),
    ("p", "Counting conventions explain part of the J_{H}4 discrepancies: SHMtrace uses analyzed clones and the shared "
          "529-bp window, whereas the publications include all filed clones at a fixed consensus length. Applying the "
          "published denominators improved knockout agreement but did not resolve Wt differences (Tables S4 and S5)."),

    ("h", "Long reads retain whole-amplicon profiles and within-lineage linkage"),
    ("p", "The J_{H}4 amplicon extends analysis to positions 90–1,230 (Figure 6A). Wt unique density was 229-fold "
          "above controls; approximately {full_450}% of mutations occurred within the first 450 bp, with a seven-fold "
          "decline between the first and last 250 bp, consistent with the reported decline downstream of V(D)J9. "
          "Full-interval and Sanger-window spectra correlated at r = {s1_r} (Figure S11)."),
    ("p", "Wt calls favored S5F-mutable positions: targeting indices were {s5f_wt} (permutation p ≤ {s5f_wt_p}), and "
          "calls per position increased {s5f_dose}-fold between the lowest and highest mutability quintiles (Figure "
          "6B). *Aicda*^{−/−} mismatch positions gave an index of {s5f_ko} (p = {s5f_ko_p}). This context-based "
          "assessment is separate from error-model fitting and can flag patterns inconsistent with expected AID "
          "targeting."),
    ("p", "Junction sequences upstream of J_{H}4 group reads into B-cell lineages (Figure 6C). In Wt2, one largely "
          "unmutated lineage accounted for {wt2_largest}% of junction-assigned reads. Read-weighted frequency was "
          "{wt2_rw} per 10^{3} bp, compared with {wt2_cw} when each lineage's majority-rule mutation profile was "
          "counted once. Wt3 and Wt4 changed less ({wt3_rw} to {wt3_cw} and {wt4_rw} to {wt4_cw}). The largest Wt3 "
          "lineage comprised {tree_reads} reads and {tree_founders} founder mutations, with {tree_sub} inferred "
          "sublineages extending to seventeen mutations (Figure 6D). Seventy percent of mutated reads fitted one tree; "
          "7/1,128 mutation pairs failed the four-gamete test. Without UMIs, PCR recombination and residual errors "
          "remain possible, so these are read-level reconstructions rather than proven cellular lineages."),

    ("h", "Public datasets test transfer across platforms and experimental designs"),
    ("p", "Three public Illumina studies tested different biological contexts (Figure 7). At J_{H}4 positions 53–122, "
          "Wt mice showed 56-fold greater mutation frequency than AID-null mice (log-rate Welch p = {ldha_p}; n = "
          "{ldha_nwt} and {ldha_nko}). *Ldha* deletion showed no detectable difference from Wt (Mann–Whitney p = "
          "{ldha_p_ko}), consistent with the original study[@sharma2023] (Figure 7A)."),
    ("p", "In Aicda−/−Ung−/−Msh2−/− CH12F3 cells, ectopic AID increased frequency {ch12_fold}-fold over empty vector "
          "(p = {ch12_p}; n = 4 each), with one control `call` in 20,000 reads38. All AID calls were C→T/G→A, "
          "consistent with loss of uracil excision and mismatch repair; R50E yielded none (Figure 7B). R174E/R178E "
          "appeared more active but gave only {ch12_mut_ctga}% transitions, a {ch12_cluster_pct}% cluster within "
          "positions {ch12_cluster} and targeting index {ch12_mut_s5f} versus {ch12_aid_s5f} for AID. These features "
          "suggest a library-specific artifact absent from controls (Figure 7C). Establishing its source would require "
          "additional sequencing provenance or independent material."),
    ("p", "HMCES-deficient Ramos clones showed {hmces_del_fold}-fold more deletions than Wt (p = {hmces_del_p}; n = 6 "
          "each), versus {hmces_ins_fold}-fold more insertions39 (Figure 7D). One AID-knockout control carried "
          "{hmces_ko3_calls} substitutions, {hmces_ko3_gc}% at G:C, while the other **two** had none. Including it as "
          "a decoy suppressed all substitution calls. The noisy-control rule excluded it from threshold-setting decoys "
          "while retaining it in model fitting and reporting the exclusion; the rule also activated in OT1 and mOrange "
          "analyses. Retained Wt/HMCES-knockout counts were {hm_calls} per clone, with {hmces_gc_range}% at G:C. "
          "Treating the noisy clone as a sample instead yielded {hmk3_calls} calls and {hmk3_gc}% at G:C, "
          "demonstrating sensitivity to control assignment."),
    ("p", "On dual-UMI R10 rRNA-operon reads40, consensus reduced error 57-fold (95% CI 35–99) to 0.17%, with residual "
          "errors concentrated in homopolymer indels. At ≥10 reads per family, substitution error approached one per "
          "10,000 bases (Figure 7E). This is a technical consensus benchmark rather than an SHM experiment."),
    ("p", "An mOrange reporter study used AID-ER or catalytically dead AID-E58Q-ER in NIH-3T3 cells, with empty, Ugi "
          "or UNG repair vectors and 36-bp sequencing41. Raw Q>20 mismatch counting gave dead-AID backgrounds of "
          "{mo_raw_lo}–{mo_raw_hi} per 10^{3} bases, leaving only two-fold AID separation in the empty arm. "
          "Calibration against that arm's dead-AID libraries retained {mo_calls} calls across nine AID libraries, all "
          "at G:C and {mo_wrc}% in WRC/GYW motifs. Dead-AID Ugi and UNG samples retained zero and **two** calls. "
          "Separation was {mo_control_fold}-, {mo_ugi_fold}- and {mo_ung_fold}-fold; the zero-call Ugi denominator "
          "used one event. Transition shares fell from {mo_ugi_ts}% with Ugi to {mo_ung_ts}% with UNG, whereas raw "
          "backgrounds contained {mo_dead_ts}% transitions (Figure 7F). The paired difference was {mo_arm_diff} "
          "percentage points (95% CI {mo_arm_ci}; p = 0.0021), based on few calls per library."),
    ("p", "Four engineered-diversification datasets tested analysis against experimental starting constructs (Figure "
          "S18). On a common GFP substrate42, no-editor libraries yielded at most one `call` in 20,000 reads. AID and "
          "APOBEC3A differed in WRC/GYW share in Wt cells ({eng_wt_wrc_aid}% versus {eng_wt_wrc_apo}%; p = "
          "{eng_wt_wrc_p}) and UNG/MSH2-null cells ({eng_ko_wrc_aid}% versus {eng_ko_wrc_apo}%; p = 9.0 × 10−4), but "
          "transition shares converged in the repair-null background ({eng_ko_ctga_aid}% versus {eng_ko_ctga_apo}%; p "
          "= {eng_ko_ctga_p}). Transition enrichment alone therefore did not distinguish these mutators. TAM "
          "profiles43 showed guide-associated peaks and increasing relative-excess sites with guide-pool size. "
          "However, its single-control reporter calibration returned {eng_tam_case} calls at maximal raw excess and "
          "{eng_tam_ctrl} in the control; that arm is reported as mismatch rates. Selected antibody and TCR "
          "sequences44,45 provided known-reference demonstrations without replicated AID-null controls. Additional "
          "analyses examined PI3Kδ-inhibitor Sμ libraries46 and PacBio CSR junctions[@vincendeau2022], with "
          "assay-specific limitations (Figures S16 and S17)."),
]

DISCUSSION = [
    "SHMtrace combines empirical background modeling with biological-replicate comparisons. Control mismatch rates "
    "inform substitution tests by type, strand and quality, while held-out libraries assess calibration. This "
    "addresses recurrent background inadequately captured by quality scores alone. Independent-sample inference "
    "addresses a separate problem: pooling reads produced false positives in more than 90% of null simulations. The "
    "contribution is a reproducible framework linking these decisions to mutation profiles and transparent read "
    "accounting.",
    "Comparator performance depends on the assay and operating range. LoFreq was designed for deep short-read data "
    "with calibrated qualities; its base-alignment correction suppressed nearly all nanopore calls here, and disabling "
    "it still limited recovery of low-support variants. These results do not establish poor performance on its "
    "intended data. Above Q40, measured control error stopped decreasing and varied between runs ({q40_lo}–{q40_hi} × "
    "10^{−5}; Figure S2), supporting empirical assessment of background. ClairS and DeepSomatic recovered clean "
    "variants at higher allele fractions but few real SHM calls; chemistry/model mismatches further qualify "
    "external-data comparisons. Germline callers were expected to miss rare variants. Geneious provided greater "
    "single-read recall on J_{H}4 and *Pim1* but more false calls on each amplicon. Unresolved single-read events can "
    "instead be estimated in aggregate. Sensitive mode increases the empirical target to 0.20: single-read recall "
    "reached {sens_jh4_1}% on J_{H}4 and {sens_n114_1}% on N114, with observed FDR {sens_jh4_fdr} and {sens_n114_fdr}. "
    "These trade-offs should guide settings according to the experimental question.",
    "Long reads preserve the positional decline of SHM and linkage across the amplicon. Junction grouping additionally "
    "permits lineage-based summaries, although the majority-rule clone-weighted measure differs from the mean burden "
    "of all cells and excludes unassigned reads. Mutation trees provide read-level descriptions of shared variants, "
    "subject to PCR and sequencing limitations. Matched Sanger analyses supported agreement in pooled spectra, "
    "position and hotspot use, while also exposing sample-tracking errors (Figure S12).",
    "Published sequence collections can support reproducible reanalysis when clones are linked to animals and "
    "genotypes. SHMtrace recovered the principal J_{H}4 frequency and G:C comparisons and closely matched SW_{HEL} "
    "counts across five arms and twelve substitution types. Agreement was incomplete: Wt frequencies differed after "
    "denominator harmonization, the published *Mkln1* A:T decrease was not significant in the reanalysis, and the "
    "SW_{HEL} total-mutation inference differed under a different test (Tables S4 and S5; Figures S8–S10). Deposits "
    "should therefore include explicit clone-to-animal assignments and original counting conventions, alongside "
    "sequences and chromatograms.",
    "Read accounting, fixed-seed sampling and a shared desktop/command-line engine support routine use. The power "
    "analysis also illustrates the importance of biological replication: under variation estimated from three Wt "
    "J_{H}4 mice, detecting a two-fold difference requires approximately twenty animals per group, with considerable "
    "uncertainty in that estimate (Figure S3B). Increasing sequencing depth cannot substitute for independent "
    "biological samples.",
]
LIMITATIONS = (
    "Control calibration depends on representative background libraries prepared and sequenced with the samples. "
    "Transfer across runs or basecallers was not established, and one core nanopore run per assay limits "
    "generalization. Without UMIs, early PCR errors, recombination and duplicates cannot be separated reliably from "
    "mutations or expanded lineages; trees remain read-level reconstructions. The ≥50% heuristic can remove genuine "
    "lineage-trunk mutations, although excluded sites are reported and known-reference experiments can use an opt-in "
    "exception. Uncalled single-read mutations are estimated only in aggregate, and unique density requires matched "
    "depth. Noisy controls can raise observed FDR above its target; outlier exclusion requires at least three controls "
    "and must be inspected. Sample-specific artifacts absent from controls may persist despite spectrum and clustering "
    "flags. Weak targets such as *Pim1* remained close to background. Indels were tested with nanopore spike-ins and "
    "one Illumina phenotype, not a nanopore biological indel standard; UMI validation used rRNA rather than SHM. "
    "Molecule-level consensus is preferable when available. Two Sanger studies originated from the authors' "
    "laboratory, and one required reconstructed assignments; the independent study supplied sequences without "
    "chromatograms or AID-null controls. These examples validate different workflow components and do not collectively "
    "establish universal error control.")

RESOURCE = [
    ("Lead contact", "Requests for further information and resources should be directed to and will be fulfilled by "
                     "the lead contact, Alberto Martin (alberto.martin@utoronto.ca)."),
    ("Materials availability", "This study did not generate new unique reagents."),
    ("Data and code availability", "J_{H}4 and N114 nanopore reads, Sanger chromatograms, clone records, project "
                                   "files, read accounting and benchmark outputs will be deposited at {lab_data_doi} "
                                   "and released by publication. Public-data accessions are listed in Table S6 and the "
                                   "key resources table. SHMtrace v{shmtrace_version}, its desktop application, and "
                                   "analysis/figure code will be deposited at [SHMtrace Zenodo DOI]. The planned "
                                   "deposit includes SHA-256 manifests, sample sheets, figure-source tables, "
                                   "environment and tool versions, the manuscript-build commit, and commands linking "
                                   "reported values to inputs and outputs. SHMtrace alignment uses mappy; LoFreq "
                                   "preparation uses standalone minimap2, with both versions listed. Additional "
                                   "information is available from the lead contact."),
]
ACK = ("We thank the University of Toronto flow cytometry facility and members of the Martin laboratory for "
       "discussions. Long-read sequencing was performed by Plasmidsaurus, and cloned-amplicon Sanger sequencing by "
       "The Centre for Applied Genomics (SickKids, Toronto). This work was supported by the Canadian Institutes of "
       "Health Research (PJT-180269).")
CONTRIB = ("Conceptualization, A.A. and A.M.; methodology, A.A.; software, A.A.; validation, A.A.; formal analysis, "
           "A.A.; investigation, A.A.; data curation, A.A.; writing – original draft, A.A.; writing – review & "
           "editing, A.A. and A.M.; visualization, A.A.; supervision, A.M.; funding acquisition, A.M.")
INTERESTS = "The authors declare no competing interests."
AI_DECL = ("During the preparation of this work, the authors used Anthropic Claude and OpenAI Codex to assist with drafting the "
           "analysis and figure code, language editing and figure layout. "
           "After using this tool, the authors reviewed and edited the content "
           "as needed and take full responsibility for the content of the published article.")
