"""STAR Methods, key resources table, figure legends, the dataset table (Table S6) and the
supplemental index for the
Cell Reports Methods article. Markup: see docbuild.py. Bracketed items ([...]) are details only the
authors can supply (catalog numbers, accessions) and are listed for completion in the cover note."""

KRT_HEADER = ["REAGENT or RESOURCE", "SOURCE", "IDENTIFIER"]
KRT = [
    ("**Antibodies**", "", ""),
    ("Anti-mouse B220 (CD45R)–FITC, clone RA3-6B2", "{lab_b220_source}", "{lab_b220_id}"),
    ("Anti-mouse GL7–eFluor 660, clone GL7", "{lab_gl7_source}", "{lab_gl7_id}"),
    ("Anti-mouse Fas (CD95)–BV421, clone Jo2", "{lab_fas_source}", "{lab_fas_id}"),
    ("**Chemicals, peptides, and recombinant proteins**", "", ""),
    ("Platinum SuperFi II PCR Master Mix", "Thermo Fisher Scientific", "{lab_superfi_id}"),
    ("Proteinase K", "{lab_protk_source}", "{lab_protk_id}"),
    ("**Critical commercial assays**", "", ""),
    ("Monarch PCR & DNA Cleanup Kit", "New England Biolabs", "{lab_monarch_id}"),
    ("Premium PCR (full-length amplicon) sequencing, R10.4.1, Kit 14", "Plasmidsaurus", "https://www.plasmidsaurus.com"),
    ("**Deposited data**", "", ""),
    ("J_{H}4 and N114 nanopore reads, Sanger chromatograms and SHMtrace outputs", "This paper", "{lab_sra}; {lab_data_doi}"),
    ("J_{H}4 Sanger chromatograms, *Maea*^{fl/fl}*Mb1*^{cre} study (cohort 2)", "Barbulescu et al.[@barbulescu2024]; "
     "this paper", "{lab_data_doi}"),
    ("Nanopore amplicons, Sμ and *Pim1*, AID-competent and *Aicda*^{−/−}", "NCBI SRA", "BioProject PRJNA1481264"),
    ("Illumina J_{H}4, *Ldha* study", "Sharma et al.[@sharma2023]", "BioProject PRJNA973112; GEO GSE232658"),
    ("Illumina CH12F3 VDJ, ectopic AID", "Wang et al.[@wang2023]", "BioProject PRJNA918596"),
    ("Illumina Ramos IGH VDJ, HMCES", "Wu et al.[@wu2022]", "BioProject PRJNA820976"),
    ("Illumina mOrange reporter, ectopic AID in NIH-3T3 cells", "Pérez-Durán et al.[@perezduran2012]",
     "BioProject PRJNA170048; GEO GSE39108"),
    ("R10 nanopore reads with dual UMIs, ZymoBIOMICS operons", "Karst et al.[@karst2021]", "github.com/SorenKarst/longread_umi"),
    ("PacBio Sμ–Sγ1 class-switch products, *Shld1* and *Xlf* mice", "Vincendeau et al.[@vincendeau2022]",
     "BioProject PRJNA831666"),
    ("Illumina Sμ mutation libraries, PI3Kδ inhibitors", "Compagno et al.[@compagno2017]",
     "BioProject PRJNA311552; GEO GSE77788"),
    ("Illumina GFP reporter, diversifying base editors", "Liu et al.[@liu2018]", "SRA SRP159221"),
    ("Illumina reporter and ABL exon 6, targeted AID mutagenesis", "Ma et al.[@ma2016]", "SRA SRP075898"),
    ("PacBio MAGE-A1 TCR, AID diversification", "Bassan et al.[@bassan2021]",
     "figshare 10.6084/m9.figshare.12803501"),
    ("Anti-HLA-A2 antibody clone sequences", "Devilder et al.[@devilder2019]",
     "the study's Additional file 1"),
    ("**Experimental models: Cell lines**", "", ""),
    ("N114 hybridoma", "Martin et al.[@martin2002]", "N/A"),
    ("**Experimental models: Organisms/strains**", "", ""),
    ("Mouse: C57BL/6N wild type", "{lab_wt_source}", "{lab_wt_rrid}"),
    ("Mouse: *Aicda*^{−/−} (C57BL/6 background)", "Muramatsu et al.[@muramatsu2000]", "{lab_aicda_rrid}"),
    ("**Oligonucleotides**", "", ""),
    ("J_{H}4 intron primers", "Jolly et al.[@jolly1997]", "{lab_jh4_primers}"),
    ("N114 V-region primers", "Martin et al.[@martin2002]", "{lab_n114_primers}"),
    ("**Recombinant DNA**", "", ""),
    ("pMX-*Aicda* and pMX empty vector (GFP, puromycin resistance)", "Martin laboratory[@feng2021]", "N/A"),
    ("**Software and algorithms**", "", ""),
    ("SHMtrace v{shmtrace_version} (command line and desktop application)", "This paper", "{lab_code_doi}; {lab_github}"),
    ("minimap2 2.30 (mappy 2.31)", "Li[@li2018]", "https://github.com/lh3/minimap2"),
    ("LoFreq 2.1.5", "Wilm et al.[@wilm2012]", "https://github.com/CSB5/lofreq"),
    ("ClairS 0.5.1 and Clair3 2.0.3 (Docker image hkubal/clairs)", "Zheng et al.[@zheng2026,@zheng2022]",
     "https://github.com/HKU-BAL/ClairS"),
    ("DeepVariant 1.9.0 and DeepSomatic 1.9.0 (Docker images google/deepvariant, google/deepsomatic)",
     "Poplin et al.[@poplin2018]; Park et al.[@park2026]", "https://github.com/google/deepvariant"),
    ("SAMtools 1.24", "Danecek et al.[@danecek2021]", "https://github.com/samtools/samtools"),
    ("SHazaM 1.2.0 (HH_S5F model)", "Gupta et al.[@gupta2015]; Yaari et al.[@yaari2013]", "https://immcantation.readthedocs.io"),
    ("Python 3.13, NumPy 2.5, SciPy 1.18", "Virtanen et al.[@virtanen2020]", "https://scipy.org"),
    ("SRA Toolkit 3.4.1", "NCBI", "https://github.com/ncbi/sra-tools"),
    ("Geneious Prime 2026.1.2 (comparison only)", "Kearse et al.[@kearse2012]", "https://www.geneious.com"),
    ("**Other**", "", ""),
    ("R10.4.1 flow cells, Ligation Sequencing Kit V14 (via provider)", "Oxford Nanopore Technologies", "N/A"),
]

EXPERIMENTAL = [
    ("Mice", "Four-month-old Wt and *Aicda*^{−/−} mice were maintained on a predominantly C57BL/6 background. Both "
             "sexes were included (n = 3 per genotype; {lab_mouse_sex}), "
             "but sex-specific effects were not evaluated. The J_{H}4 comparison used the same mice previously "
             "analyzed by cloned-amplicon Sanger sequencing[@barbulescu2025] to permit paired platform assessment. "
             "Procedures were approved by the University of Toronto Animal Care Committee (protocol 20011472)."),
    ("Cell lines", "N114 authentication: {lab_n114_auth}. "
                   "Mycoplasma status: {lab_n114_myco}. Cells48,52 were retrovirally transduced with pMX-*Aicda* or "
                   "matched empty vector, both encoding GFP and puromycin resistance. GFP-high cells were selected and "
                   "subcloned by limiting dilution; three puromycin-resistant subclones per vector were maintained for "
                   "1 month before DNA collection. One library per subclone was analyzed. When repeat sequencing was "
                   "available (C1 and D3 from stocks frozen on different dates; **two** C3 libraries), the library "
                   "with the most usable reads was selected. The alternative C3 library was retained for sensitivity "
                   "analysis."),
]

METHODS = [
    ("h2", "Germinal-center B cell isolation and DNA extraction"),
    ("p", "Peyer's patches were harvested post-mortem into cold FACS buffer ({lab_facs_buffer}) and dissociated "
          "into single-cell suspensions by pressing them through a 70-µm strainer, then washed once. Cells were "
          "stained for 20–30 min on ice with {lab_stain} and washed twice. Germinal-center B cells "
          "(B220^{+}GL7^{+}Fas^{+}) were isolated by fluorescence-activated cell sorting on a BD FACSymphony S6 SE "
          "or FACSAria II, "
          "collecting approximately 10,000 cells per mouse; wild-type and *Aicda*^{−/−} mice were processed on "
          "the same day with the same reagents, so that genotype is not confounded with handling. Genomic DNA "
          "was released from the sorted cells, or from 5 × 10^{6} hybridoma cells, by proteinase K digestion "
          "({lab_lysis}), and lysate quality was checked by A260/280 and A260/230 ratios."),
    ("h2", "Amplification and sequencing"),
    ("p", "The 1,230-bp J_{H}4 intron amplicon[@jolly1997,@barbulescu2024] and the 655-bp N114 V-region "
          "amplicon[@martin2002] were amplified with Platinum SuperFi II high-fidelity DNA polymerase. Each 50-µL "
          "reaction held 1× SuperFi II master mix, 0.5 µM of each primer and {lab_lysate_volume} of lysate, with "
          "a no-template control in every setup; wild-type and AID-null samples were amplified together. The "
          "J_{H}4 primers were {lab_jh4_primers}, and the N114 primers {lab_n114_primers}. Cycling was 98 °C for "
          "5 min; 35 cycles of 98 °C for 30 s, 68.8 °C for 30 s, and 72 °C for 1 min; and a final 72 °C for 5 min. "
          "Products were checked on a 1% agarose gel for a single band of the expected size, then column-purified "
          "and submitted as separate, individually barcoded samples for full-length "
          "amplicon sequencing (Plasmidsaurus Premium PCR service), which uses ligation library preparation, "
          "R10.4.1 flow cells, Kit 14 chemistry, and super-accuracy basecalling "
          "(dna_r10.4.1_e8.2_400bps_sup@v5.2.0). Between 2,336 and 5,000 reads were delivered per library. For the "
          "Sanger benchmark, J_{H}4 amplicons from separate PCRs of the same DNA had been cloned and "
          "sequenced[@barbulescu2025], and the 447 resulting chromatograms were re-analyzed here."),
    ("h2", "Read preparation"),
    ("p", "Inputs are FASTQ files, optionally gzipped. Overlapping Illumina pairs are merged where appropriate; "
          "dataset-specific single-end handling is described below. Reads exceeding the locus-specific length limit "
          "(2,000 bp for J_{H}4; 900 bp for N114) or mean expected error of 0.01 per base, calculated from Phred "
          "scores28, are excluded. Competitive minimap2 alignment uses `map-ont` for nanopore or sr for Illumina, "
          "assigning each read to its best-matching reference. Unmapped reads, competing-reference assignments, "
          "supplementary alignments and alignments shorter than 800 bp for J_{H}4 or 400 bp for N114 are removed. "
          "J_{H}4 competitors comprise an approximately 12%-divergent Agouti-type paralog consensus and a recurrent "
          "contaminant haplotype. To remove full-length paralogs that align preferentially to the longer target, an "
          "additional filter excludes reads carrying ≥10 competitor-specific bases at positions where those alleles "
          "occur in <20% of sample reads. Genuine target reads carry 0–7 such bases, compared with 11–62 in "
          "paralogs/chimeras. Every exclusion is recorded. Retained reads are sampled with Python random, seed 0: 500 "
          "per JH4/N114 library and 2,000 per IRCM library. Windows are J_{H}4 101–629 (529 bp) or 90–1,230, N114 "
          "21–635, Sμ 21–631 (611 bp), *Pim1* 21–882 (862 bp), OT1 21–988 (968 bp) and OT2 21–962 (942 bp). Without a "
          "preset, the default excludes twenty bases at each reference end."),
    ("h2", "Control-calibrated variant calling"),
    ("p", "Each observed single-nucleotide difference within the window is a candidate, including substitutions "
          "supported by one read. Control error rates are estimated by substitution type, strand and quality "
          "threshold, combining a pooled rate with position-specific counts. For a candidate with k supporting reads "
          "ranked q_{1} ≥ q_{2} ≥ …, the test evaluates, for each j ≤ k, the probability of at least j errors with "
          "quality ≥qj, using strand-specific depth and control rates multiplied by margin c. The reported p value is "
          "the smallest tail probability multiplied by k."),
    ("p", "The margin c is chosen by leaving one control out at a time. Each control is tested against an error "
          "model built from the other controls, and c is the smallest value in the grid 1, 2, 3, 5, 10, 20, and "
          "50 at which no more than 1% of candidate tests pooled across held-out controls reach p < 0.01. In the analyses reported "
          "here, c was {loco_margins}. The same held-out tests, made at the samples' analysis depth, are then used "
          "to set the significance threshold. Each held-out control is kept as a decoy, defined by its sorted p "
          "values and the number of bases it tested. For a sample and a candidate threshold t, the expected number "
          "of false calls is the mean, across decoys, of the number of decoy p values ≤ t, scaled by the ratio of "
          "the sample's tested bases to the decoy's when that ratio exceeds 1. Decoys are never scaled down to a "
          "shallower library, whose coarser p values they would otherwise underestimate. SHMtrace keeps the "
          "largest t at which the expected number of false calls is at most 0.05 times the number of sample p "
          "values ≤ t, an empirical false-discovery estimate inspired by target–decoy analysis[@elias2007]. A control is never used as its own decoy: "
          "when a control is tested, only the other controls serve as decoys. A control whose held-out candidates "
          "at p ≤ 0.01 per base exceed five times the median of the other controls, and number at least five, is "
          "left out of the decoys (unless no decoy would remain) and named in the report; it still contributes to "
          "the error model, and is never removed from the analysis. The comparison needs a majority to judge "
          "against, so it is made only when at least two other controls remain, that is with three or more "
          "controls; with exactly two, the more variable library is named in the report but not left out, "
          "because which of the two is unusual cannot be established. The rule was written for one public "
          "dataset, the Ramos clones of Figure 7D, in which a single *Aicda*^{−/−} clone suppressed every call; "
          "it was then applied to the other datasets and was already fixed when the OT1 and OT2 libraries were "
          "analyzed; at OT1, control KO2 was also excluded from the decoys. It is therefore a rule derived "
          "from data rather than a pre-registered one, and it is reported here so that its effect can be judged. "
          "With fewer than two controls, or for indels, the Benjamini\u2013Hochberg "
          "procedure[@benjamini1995] across positions and alternative alleles is used instead. "
          "The number of controls therefore selects between four modes, and the output records which was "
          "used in its `multiple_testing` field, so a reader never has to infer it. With **no control** the "
          "background is not measured at all: every candidate is judged from base qualities, the correction "
          "is Benjamini\u2013Hochberg, indels are not called, and the run carries a note saying so. With "
          "**one control** that library supplies the error model, but the margin cannot be fitted by leaving "
          "a control out and the default of 10, measured on J_{H}4, is used; a control cannot be its own "
          "decoy, so there is no decoy set and the correction is again Benjamini\u2013Hochberg; the "
          "noisy-control rule cannot act, since it needs two others to compare against. With **two** each is "
          "the other's decoy, the margin is fitted but can reach the top of its grid, and a disagreement "
          "between the pair is reported without acting on it, because with two there is no majority to say "
          "which is wrong. With **three or more** the margin is fitted by leave-one-out, the decoys are the "
          "peers, and a far noisier control is set aside automatically. The TAM reporter arm of Figure S18 is "
          "the one-control case, and its calibration inverted (Results); a single control is enough to run "
          "but not to calibrate, and results from it should be reported as mismatch rates without an "
          "FDR claim. "

          "Benjamini–Hochberg is not avoided because the p values are discrete: it controls the false-discovery "
          "rate for valid discrete p values, conservatively. It is avoided because of how candidates arise. A "
          "position enters the test only when at least one read already differs from the reference there, so "
          "the tested set is chosen by the same data that supply the evidence, and the null distribution of a "
          "candidate's p value is conditioned on that selection rather than superuniform on (0, 1]. The "
          "control-based threshold applies the same candidate-selection rule to samples and held-out controls, "
          "so it measures the corresponding background empirically. This does not establish formal FDR control: "
          "the decoys are few, their fitted nulls differ from the sample null, and their exchangeability with "
          "samples is an assumption. Figure S15 tests selected consequences of these limitations. Alternative "
          "alleles carried by at least 50% of covering reads are removed as germline or parental differences; "
          "other alleles at the same position remain eligible. Without controls, recommended mode uses a "
          "Poisson approximation to read-quality error probabilities; the separate pbinom mode implements "
          "the exact Poisson-binomial calculation."),
    ("p", "Per-sample summaries report G:C share, transition share, G:C substitutions at WRC/GYW and A:T substitutions "
          "at WA/TW, with the corresponding motif opportunities in the analysis window. These describe spectrum and "
          "context, not uniquely identifiable repair-pathway contributions: uracil excision, mismatch repair and "
          "replication over uracil produce overlapping outcomes. Repair interpretations in *Ung*^{−/−}*Msh2*^{−/−} "
          "cells and Ugi/UNG fibroblasts are therefore grounded in the experimental perturbations rather than inferred "
          "from spectra alone."),
    ("h2", "Testing the calling procedure as a whole"),
    ("p", "End-to-end assessment used partitions and planted mutations from three existing *Aicda*^{−/−} J_{H}4 "
          "libraries (`crm/fdr_validation.py`). Cases used KO1 reads. Two-control analyses used one partition each "
          "from KO2/KO3; three-control analyses used **two** KO2 partitions and one KO3 partition. Case/control sets "
          "shared no read records, but without UMIs this does not establish original-molecule independence. Partitions "
          "are technical resamples rather than independent biological controls. Each `run_project` analysis reselected "
          "its margin, threshold and exclusions. Conditions were unmodified background, one noisy control, case-only "
          "artifact and shared case/control artifact. Artifacts occupied six positions at nominal 6% support; shared "
          "artifacts used identical alternate alleles, and case artifacts avoided planted biological sites. Control "
          "number (two/three), depth (250/500/1,000) and mutation burden (none/approximately 400 true recovered calls) "
          "were varied. Exclusion-on/off analyses used identical inputs (Figure S15). These stress tests assess "
          "behavior conditional on existing libraries, not transfer to a new run."),
    ("p", "Under unmodified all-null background, discoveries occurred in {fdv_null_none_on_any}/{fdv_null_none_on_n} "
          "exclusion-on and {fdv_null_none_off_any}/{fdv_null_none_off_n} exclusion-off resamples; with one noisy "
          "control, both were {fdv_null_one_noisy_off_any}/{fdv_null_one_noisy_off_n}. Clopper–Pearson 95% upper "
          "bounds of {fdv_null_none_on_upper}% and {fdv_null_none_off_upper}% for the unmodified condition do not "
          "establish a 5% bound. Shared-artifact resamples yielded {fdv_null_shared_on_any}/{fdv_null_shared_on_n} "
          "discoveries in both settings, each with a {fdv_null_shared_on_upper}% upper bound. Case-only artifacts "
          "yielded discoveries in {fdv_null_case_only_on_any}/{fdv_null_case_only_on_n} under both settings, averaging "
          "{fdv_null_case_only_on_false} false calls. With biological mutations and exclusion enabled, mean "
          "per-library FDP was {fdv_polyclonal_none_on_fdp}%, {fdv_polyclonal_one_noisy_on_fdp}%, "
          "{fdv_polyclonal_shared_on_fdp}% and {fdv_polyclonal_case_only_on_fdp}% for unmodified, noisy-control, "
          "shared-artifact and case-only conditions, respectively (five resamples each). Per-library FDP, pooled "
          "false/(true+false), and the all-null probability of any discovery are reported separately. These small "
          "conditional experiments do not establish general error control."),
    ("p", "The noisy control was excluded in {fdv_null_one_noisy_on_excluded}/{fdv_null_one_noisy_on_n} all-null "
          "rule-on resamples. With biological mutations, exclusion-on/off analyses averaged {fdv_polyclonal_one_noisy_on_true}/{fdv_polyclonal_one_noisy_off_true} true calls and "
          "{fdv_polyclonal_one_noisy_on_false}/{fdv_polyclonal_one_noisy_off_false} false calls across five paired resamples. At least three controls are required to assess an "
          "outlier against peers. When the peer median is zero, ≥5 held-out candidates at p≤0.01 trigger exclusion; an "
          "unmodified control was consequently excluded in {fdv_null_none_on_excluded}/{fdv_null_none_on_n} all-null "
          "resamples. Every exclusion is recorded and warrants review, including sensitivity to retaining that "
          "control."),
    ("p", "The default mean-decoy estimator was compared with `control_decoy_estimate=\"max\"`, which requires the "
          "threshold to satisfy the noisiest retained decoy (`crm/fdr_estimator.py`). On unmodified mutation-bearing "
          "inputs, mean FDP was {fdv_polyclonal_none_on_fdp}% versus {fdv_max_polyclonal_none_on_fdp}%, with "
          "{fdv_polyclonal_none_on_true} versus {fdv_max_polyclonal_none_on_true} true calls. For case-only artifacts, "
          "FDP was {fdv_polyclonal_case_only_on_fdp}% versus {fdv_max_polyclonal_case_only_on_fdp}%. Five paired "
          "resamples per condition used identical inputs. This compares implemented estimators; it is neither "
          "independent validation nor proof that the maximum guarantees FDR control."),
    ("p", "Reported runtimes are wall-clock seconds for one `run_project` `call` on an Apple M2 Pro (12 cores, 32 GB; "
          "macOS), from prepared input to written results. They include file reading, mappy alignment, calling and "
          "sample summaries, but exclude reference construction and external benchmark containers. Sanger projects "
          "contain a few hundred chromatograms; nanopore projects contain six libraries analyzed at 500–2,000 reads. "
          "Input recovery and curation, including reconstruction of clone assignments from historical projects, took "
          "substantially longer and are excluded from compute times."),
    ("p", "The Geneious rule was reproduced as an uncorrected approximate "
          "p ≤ 0.05 per candidate, which matched Geneious exports site by site (Jaccard index 0.97–0.98 on Wt "
          "J_{H}4)."),
    ("h2", "What each mutation measure counts"),
    ("p", "Three measures are reported per library, and they answer different questions; group tests then use one "
          "value per animal. A unique variant is a position and substituted base counted once per library, "
          "however many reads carry it; unique means distinct within that library and does not imply an "
          "independent mutational origin, since one lineage's trunk mutation is one unique variant. "
          "**Unique-variant density** is D / B × 10^{5}, where D is the number of unique variants and B the base "
          "observations inside the analysis window: the covered window positions times their mean read depth, not "
          "the number of reference positions. Because B grows with depth while D saturates, this density falls as "
          "depth rises even at a constant mutation rate, so it is compared only between libraries analyzed at the "
          "same depth (Figure S5). **Read-weighted mutation frequency** is Σ m_{i} / (N × L) × 10^{3}, where "
          "m_{i} is the number of called sites read i carries at or above the configured minimum base quality (0 for these read-library analyses), N the reads analyzed and L the "
          "window length; a mutation shared by many reads is counted once in each of them, so an expanded clone "
          "raises it. **Clone-weighted frequency** groups reads by their V(D)J junction, keeps a mutation when at "
          "least half the reads of its clone carry it (any read of a single-read clone), and is Σ c_{k} / (K × L) "
          "× 10^{3} over the K clones, so each clone counts once whatever its size. For wild-type J_{H}4 library "
          "Wt2 at 500 reads over the 529-bp window these give 221 unique variants over 263,865 base observations "
          "(83.8 × 10^{−5}), 747 read-level mutations over 264,500 bases (2.82 × 10^{−3}), and 207 clonal "
          "mutations over 24 clones (16.3 × 10^{−3}); the clone-weighted value is the largest because one clone "
          "holds 81% of the reads. For Sanger clones, where each clone is one molecule read once, the frequency "
          "denominator is the clones analyzed times L, with L the 529-bp window; the published studies used their "
          "own full consensus length and counted every filed clone, which is the convention Table S4 reports "
          "separately."),
    ("p", "Position-specific substitution rates are shrunk toward pooled type/strand/quality rates using m = 500 "
          "pseudo-observations, with the pooled rate as a floor. Recurrently erroneous sites retain higher rates, "
          "while zero observed errors do not imply error-free sequence. An unobserved substitution class receives half "
          "a pseudo-count over observed control bases to keep quality-tail probabilities above zero. Indels use a "
          "separate control model because their errors are strongly systematic. Each exact event is tested against a "
          "class-rate floor defined by insertion/deletion, length up to 4 and homopolymer run up to 5, with positional "
          "shrinkage and a fitted inflation margin. Benjamini–Hochberg corrects candidate tests. At least three "
          "supporting reads and Q≥30 at inserted or deletion-flanking bases are required. UMI grouping permits up to "
          "three edits per tag end. Abundance-ranked greedy merging admits a neighboring tag only when it is at most "
          "half as abundant. Families with fewer than three reads are discarded; remaining families yield "
          "quality-weighted consensus sequences. The aggregate single-read estimator scales Q≥35 excess mismatches by "
          "the fraction of correctly read sample bases reaching that quality. Stability across Q30–Q40 was assessed "
          "for J_{H}4."),
    ("h2", "Single-read estimate, clones, mutation trees, and quality checks"),
    ("p", "Single-read mutations not called individually are estimated in aggregate. SHMtrace counts Q≥35 single-read "
          "substitutions outside the accepted `call` set, subtracts the control expectation scaled to sample bases, "
          "and divides by the fraction of correctly read bases reaching that quality. The 95% interval incorporates "
          "Poisson variation, between-control overdispersion and control-rate uncertainty."),
    ("p", "Two different things are called a clone in this field, and this paper keeps them apart. A *Sanger "
          "clone* is one bacterial colony, that is one PCR product carried by one transformant and sequenced "
          "once: a unit of sequencing, which says nothing about which B cell the molecule came from, and two "
          "Sanger clones of one animal may well come from the same B cell lineage. A *B cell clone* is a "
          "lineage, defined here by a shared V(D)J junction, which is what the clone-weighted frequency and the "
          "mutation trees use. The two coincide only by accident. Where the text says a measure counts each "
          "clone once, it means the B cell clone; the Sanger denominator is clones analyzed × 529 bp, where "
          "those clones are the sequencing units. Sanger amplicons of the J_{H}4 intron do not span the "
          "junction, so B cell clones cannot be identified in them, and clone-weighted frequencies are reported "
          "only for the nanopore libraries."),
    ("p", "Junction grouping uses the twenty bases immediately upstream of alignment starts occurring within 18 bp of "
          "J_{H}4. Sequences are clustered greedily around abundance-ranked centroids with edit distance ≤3 and a free "
          "5′ end. Each lineage's majority profile retains calls supported by at least half its reads. Tree nodes "
          "merge mutations supported by ≥2 lineage reads with Jaccard overlap ≥0.9; each node's parent is the smallest "
          "node containing ≥90% of its supporting reads. A read fits when its nodes occupy one ancestral path. The "
          "four-gamete test requires each of the four presence/absence combinations for a mutation pair in at least "
          "**two** reads."),
    ("p", "The targeting index is mean S5F mutability[@yaari2013] at mutated positions divided by the window mean. A "
          "one-sided permutation test uses 2,000 random position sets; the model was exported from SHazaM 1.2.0. "
          "Positional-clustering flags identify 10-base windows with more calls than expected under a Poisson "
          "distribution, with Bonferroni correction across windows (p < 10^{−3}). Figure S13 summarizes the workflow; "
          "Figure S14 illustrates the analytical rules using J_{H}4 libraries. Sanger traces are assigned through a "
          "sample sheet and trimmed to their longest reliably aligned Q≥20 segment. All Q≥20 substitutions are counted "
          "after germline filtering. Quality-control flags identify copied files, unexpectedly shared clones between "
          "animals (Poisson test weighted by mutation frequency across animals), excess mutations in controls, label "
          "mismatches and mutation-load outliers. Flags are reported for review without automatic exclusion."),
    ("p", "The Sanger route counts Q≥20 substitutions without per-call significance testing. AID-null clones provide "
          "the observed Q≥20 mismatch rate after their own germline differences are excluded, stratified by quality "
          "and position. Applying these rates to each sample's bases gives the expected background count, reported "
          "beside observations without subtraction. Sites mismatched in more than one control clone are annotated on "
          "corresponding sample calls. This background estimate contextualizes Sanger counts but does not determine "
          "acceptance of individual substitutions."),
    ("p", "For the Maea Sanger study[@barbulescu2024], chromatograms from four 96-well plates and the pilot plate of "
          "cohort **two** were downloaded from the sequencing facility portal. Mouse assignments were reconstructed "
          "from the original project folders using screen recordings. Four traces filed under **two** mice were "
          "retained under both assignments. Published per-mouse values came from Source Data: Figure 2h for mutation "
          "frequencies, including *Aicda*^{−/−} mice; Figure 2i for G:C transitions and transversions; and Figure 2g "
          "for clone distributions. Twelve-type counts were transcribed from Supplementary Figure 9b, using printed "
          "cell values where they conflicted with printed totals. Wt and *Maea*-deficient groups were compared using "
          "the original one-tailed unpaired t test; comparisons with *Aicda*^{−/−} used Welch tests on log "
          "frequencies. Data S1 contains the sample sheet and transcribed tables."),
    ("p", "For the *Mkln1* study[@barbulescu2025], mouse assignments were reconstructed from original project folders "
          "using a screen recording. A separate folder contained 58 later-run *Mkln1*^{−/−} traces named for the same "
          "three mice. These were included under their named mice because the published substitution counts exceed the "
          "earlier preprint by eight mutations, matching the number and types of SHMtrace calls in this batch. The "
          "files do not establish which later clones entered the published tables or how the three knockout folders "
          "map to published mice #1–#3; comparisons therefore use groups rather than matched mice. Each later clone "
          "carries at most **two** mutations and cannot be placed using the matched nanopore libraries, because "
          "assignment requires a read carrying at least three of its mutations."),
    ("p", "A sensitivity analysis used an evidence-based *Mkln1* clone set with rules fixed before analysis and "
          "independent of published values. Paralog amplicons and trace copies were set aside as in the primary "
          "analysis. Lineages filed under **two** mice, identified by duplicate-clone and shared-mutation checks, were "
          "sought in matched J_{H}4 nanopore libraries after excluding positions mismatched in >0.5% of *Aicda*^{−/−} "
          "reads. A read supported a lineage when it carried more than half its shared mutations and at least three. "
          "Lineages supported by at least three reads in exactly one mouse library were reassigned to that mouse; "
          "other multiply filed lineages were set aside. AID-null clones flagged for more mutations than predicted by "
          "their base qualities would also have been excluded; none met this criterion. Clones without an original "
          "mouse assignment remained excluded, as in the primary analysis."),
    ("p", "For the SW_{HEL} study[@thientosapol2018], all 1,795 GenBank clones (MH198782–MH200576) were downloaded and "
          "treated as Sanger sequences; the sample sheet records each accession, host and well. GenBank provides no "
          "base qualities, so Q30 was assigned as a processing convention rather than measured quality. The "
          "majority-consensus reference converged in **two** iterations. The analysis window covered the interval "
          "represented in ≥99% of clones (558 bp of the 576-bp consensus). Host and genotype assignments followed "
          "deposited records, with the UGI group kept separate. No *Aicda*^{−/−} group was available for background "
          "estimation. Published per-group counts were transcribed from the source study's Table 1, which counts each "
          "sequence independently. SHMtrace totals were compared using this convention; unique-mutation counts are "
          "lower, as expected. The publication lists 1,796 sequences versus 1,795 deposited, owing to one fewer "
          "*Msh2*^{−/−} record (445 versus 446); no deposited record was excluded. The Wt deposit names twelve hosts "
          "versus thirteen reported, and AT40 contains 88 clones, approximately twice the other collections, "
          "suggesting that **two** animals may share this label. Pooled counts are unaffected. Animal-level tests use "
          "deposited host names and were repeated without AT40; clones were not arbitrarily split between inferred "
          "animals. Excluding AT40 gave p = 1.6 × 10−10, 8.5 × 10−6, 5.0 × 10−10 and 0.18 for the four substitution "
          "classes, and 2.6 × 10−4 for total mutation, versus {sw_at_tv_p}, {sw_gc_tv_p}, {sw_gc_ts_p}, 0.14 and 3.8 × "
          "10−4 with all hosts."),
    ("p", "The default high-frequency-allele rule excludes changes present in ≥50% of covering reads as potential "
          "germline or parental differences. Its effect on acquired mutations was evaluated by simulating mutated "
          "dominant lineages on J_{H}4 AID-null reads, as in Figure 3D. Dominant clones occupied 40%, 60% or 80% of "
          "reads and received trunk mutations with mean six, as did other clones. Three libraries per abundance were "
          "tested against the other **two** controls. At {mask60_rw}%, no planted mutation crossed the exclusion "
          "threshold. At 60% and 80%, {mask60_lost} mutations per library crossed it, and all were removed before "
          "testing. Read-weighted frequency fell by as much as {mask60_rw}% and {mask80_rw}%, respectively; "
          "clone-weighted frequency fell by as much as {mask60_cw}% and {mask80_cw}%. Libraries without "
          "threshold-crossing mutations were unaffected. Reports list each excluded site with supporting-read count, "
          "depth and frequency to permit inspection of dominant mutated lineages."),
    ("h2", "Benchmark against other callers"),
    ("p", "Figure 4 and Table S2 use `crm/benchmark_real.py`, which invokes calling directly to ensure identical input "
          "reads across callers. Production analyses use `run_project`. A library-level comparison "
          "(`crm/route_check.py`) found identical calls across all twelve JH4/N114 libraries. External nanopore counts "
          "differed because margins were fitted from each control's **two** peers in the benchmark and across all "
          "three controls in the production project. Examples were 88 versus 101 Sμ Wt calls and 6 versus 1 *Pim1* "
          "calls. Benchmark/production separation estimates were J_{H}4 141/141, N114 85/84, Sμ LPS 66/84, Sμ CH12 "
          "94/171, OT1 53/46 and OT2 28/24; *Pim1* remained weak at 1.8/0.8. Figure 4 consistently reports the "
          "benchmark route, while Table S1 reports production outputs for the JH4/N114 libraries with matching calls."),
    ("p", "For each library, the reads analyzed by SHMtrace were written out, aligned with minimap2 using the same "
          "preset with secondary alignments turned off, and sorted and indexed with SAMtools. They were then "
          "called with LoFreq 2.1.5 `call`, restricted to the analysis window, at three settings: default settings "
          "(base-alignment quality on, Bonferroni correction, default filters); `-B` (base-alignment quality off); "
          "and `-B -a 0.05 -b 1 --no-default-filter` (uncorrected p ≤ 0.05, filters off). Substitutions at an "
          "allele frequency below 0.5 were kept, matching SHMtrace's germline filter, and all callers used the same "
          "density denominator for each library. In the spike-in tests, 60 sites per replicate were chosen with "
          "WRC/GYW motifs weighted 4-fold and WA/TW motifs 2-fold, and each substitution was written into 1, 2, 3, "
          "5, 10, or 20 reads covering that site while keeping each read's own base quality. Three replicates were "
          "performed for each held-out control. The IRCM runs (PRJNA1481264) were downloaded with the SRA Toolkit, "
          "validated with vdb-validate, and subsampled to 6,000 reads (seed 0). The reference for each locus was "
          "the majority consensus of its first *Aicda*^{−/−} library, built with minimap2 and SAMtools consensus "
          "and iterated to convergence. The OT1 and OT2 off-targets of the same BioProject were added later and "
          "analyzed the same way, at 2,000 reads, with references built by the same consensus procedure "
          "(1,008 and 982 bp, converged in two iterations). Other loci of that BioProject were set aside: "
          "*Cd83* and *Ly6e* yielded no calls in any library, including at 30,000 reads, and the remaining "
          "cell contexts provide two libraries per group or a single *Aicda*^{−/−} library, which is below "
          "what the control-based test and the animal-level comparison need."),
    ("p", "ClairS 0.5.1[@zheng2026] and DeepSomatic 1.9.0[@park2026] (somatic callers), and Clair3 2.0.3[@zheng2022] "
          "and DeepVariant 1.9.0[@poplin2018] (germline callers), were run from official Docker images on the same "
          "planted mutations and window-restricted alignments used for LoFreq. Somatic normals pooled usable reads "
          "from the other **two** AID-null libraries, which also fitted SHMtrace background models. Default nanopore "
          "models were ClairS `ont_r10_dorado_sup_5khz`, Clair3 `r1041_e82_400bps_sup_v520`, DeepVariant `ONT_R104` "
          "and DeepSomatic `ONT`. Permissive analyses used ClairS `--snv_min_af 0.005 --min_coverage 4`; Clair3 "
          "`--snp_min_af 0.01 --min_coverage 4`; and candidate thresholds of 1% for DeepVariant and 0.5% for "
          "DeepSomatic, with **two** supporting reads. Only FILTER PASS substitutions below 0.5 allele fraction were "
          "retained, matching other callers. PRJNA1481264 reports PromethION sequencing without chemistry details. The "
          "same R10.4.1 models were used for Sμ and *Pim1*, introducing a possible chemistry/basecaller mismatch that "
          "could reduce comparator performance independently of SHM biology. These results evaluate the specified "
          "configurations and do not establish general caller superiority. Higher-fraction experiments planted thirty "
          "substitutions into 25, 50 or 100 of 500 reads (5%, 10% or 20%) from each held-out J_{H}4 control. ClairS "
          "and DeepSomatic were also applied to real J_{H}4 Wt libraries against pooled *Aicda*^{−/−} reads and to "
          "N114 pMX-*Aicda* subclones against pooled empty-vector reads, using 500 reads per library (seed 0)."),
    ("h2", "Simulations"),
    ("p", "Lineage simulations partitioned one null library into approximately forty power-law-sized clones. Each "
          "received a random 20-bp junction, J trimming from −3 to +12 bp, and S5F-weighted trunk, branch, subbranch "
          "and private mutations with gamma-distributed loads. In the expanded-clone scenario, one unmutated clone "
          "supplied 60% of reads. Each control was held out against the other **two**. Quantification was assessed at "
          "125, 250, 500 and 1,000 reads on JH4/N114 (168 libraries). Group simulations used log-normal animal rates "
          "(coefficient of variation 1), negative-binomial counts and 262,913 base observations per animal, across 300 "
          "experiments per condition. This represents approximately 526 covered bases per read at 500 reads, rather "
          "than complete coverage of the 529-bp window; observed per-read coverage ranged from 516 to 528 bases. Forty "
          "Sanger projects included realistic trace qualities, unusable traces and one of each planted tracking error."),
    ("h2", "Public Illumina and UMI datasets"),
    ("p", "Public Illumina projects used AID-null controls, Illumina settings and the S5F model: PRJNA973112, J_{H}4 "
          "positions 53–122, 5,000 merged pairs per mouse; PRJNA918596, CH12F3 positions 21–249, 5,000 reads per "
          "replicate; PRJNA820976, Ramos positions 26–446, 20,000 reads per clone with indel calling. CH12F3/Ramos "
          "references were control-read majority consensuses. Mutation frequency, rather than unique density, was "
          "compared. Overlapping pairs required ≥20 bp overlap with ≤10% mismatches (Figure S12). Distributed "
          "PRJNA973112/PRJNA918596 files have flat Q30 qualities, giving a single quality level; Ramos retains binned "
          "qualities. For mOrange PRJNA170048, the analyzed interval is original positions 75–693, equivalent to "
          "31–649 in the reconstructed 679-bp reference, whose first base corresponds to original position 45 (offset "
          "44). Analyses used 200,000 of the first 600,000 reads per library, approximately 10,600 per position. "
          "Empty-repair-arm AID-E58Q-ER libraries fitted the background; dead-AID Ugi/UNG libraries were analyzed as "
          "samples. Reference reconstruction and raw-mismatch comparisons used the original supplementary table's "
          "position, reference-base and count columns. UMI evaluation used 4,259 R10 rRNA-operon reads with dual 18-nt "
          "tags40. Grouping allowed up to three edits per end and directional rare-tag merging. Families of at least "
          "three reads produced quality-weighted consensuses, scored against 43 known mock-community operon sequences."),
    ("h2", "Public PacBio class-switch junctions"),
    ("p", "PRJNA831666[@vincendeau2022] contains Sμ–Sγ1 products from Wt, *Shld1*^{−/−}, *Xlf*^{−/−} and "
          "*Shld1*^{−/−}*Xlf*^{−/−} primary B cells, three mice per genotype with one deposited PacBio run each. "
          "Although filenames indicate subreads, reads aligning over ≥600 bp of the bundled Sμ fragment had median "
          "identity 99.8%; `map-hifi` was therefore used instead of `map-ont`. This identity measure does not "
          "independently establish CCS processing. References were extracted from GRCm39 chromosome 12 in minus-strand "
          "orientation: donor Sμ with forward-primer site, 113,386,000–113,392,500; acceptor Sγ1–Cγ1 through the "
          "reverse-primer site, 113,290,000–113,304,000. Read endpoints defined these windows. "
          "`crm/build_switch_reference.py` requires both project J_{H}4 and bundled Sμ to align within the fetched "
          "locus before writing. `shmtrace.switch` oriented reads to the donor and selected the **two** best forward "
          "donor/acceptor alignments (`crm/switch_pacbio.py`). Breakpoints within 3 bp on both sides were grouped as "
          "one unique junction. Short off-target reads were common (median approximately 350 bp versus 1.3 kb for "
          "junction reads), so delivered and analyzed counts are reported separately. Genotype comparisons use one "
          "value per mouse, with Welch and permutation tests; exact two-sided permutation p cannot fall below 0.1 at "
          "n=3 per group. Fifteen comparisons comprise five features against three genotypes. Unadjusted p values and "
          "means are reported; nine comparisons pass 0.05 both before and after joint Benjamini–Hochberg correction. "
          "Acceptor ambiguity was assessed by realigning each junction's first sixty acceptor bases with sr and "
          "identifying placements within five matched bases of the best. The module flag is uninformative here because "
          "secondary alignments are absent. Excluding ambiguous junctions in `crm/switch_ambiguity_sensitivity.py` "
          "removed 2–21% per sample; 13/15 comparisons retained their significance classification. *Shld1*^{−/−} "
          "microhomology ≥4 nt changed from p=0.065 to 0.005, while *Xlf*^{−/−} insertions changed from p=0.016 to "
          "0.060. These **two** comparisons remain provisional. The engineered-diversification analysis treats AID "
          "versus APOBEC3A as its prespecified contrast within each genotype (Figure S18A), without additional "
          "multiplicity adjustment."),
    ("h2", "Public engineered-diversification datasets"),
    ("p", "Engineered-diversification analyses used `ancestral_reference` against the experimental starting sequence. "
          "This opt-in mode retains high-frequency reference differences that could arise through selection and leaves "
          "physiological defaults unchanged (`crm/germline_policy_demo.py`). SRP159221[@liu2018] contains a 238-bp "
          "EGFP reporter edited by several deaminases, with three replicates per arm. Pairs overlap by 97.8%; the "
          "no-editor consensus supplies the reference. The prespecified comparison was DBE-AIDmono versus DBE-humanA3a "
          "in Wt and UNG/MSH2-null 293T at 20,000 reads per library (`crm/liu_mutators.py`). Motif/spectrum shares "
          "summarize accepted calls and do not determine inclusion. SRP075898[@ma2016] multiplexes amplicons within "
          "files. Each amplicon was analyzed against the others as competing references. `SRR3608504` exports paired "
          "and single-end blocks despite its single-end label; only the single-end block was used to harmonize read "
          "representation. Its reporter has one no-editor control and is reported as mismatch rates. ABL exon-6 guide "
          "and saturation libraries were compared with one another as descriptive local baselines "
          "(`crm/tam_reporter.py`, `crm/tam_footprint.py`, `crm/tam_saturation.py`). Fragmented genome-wide off-target "
          "material was not analyzed through amplicon gates. The MAGE-A1 TCR deposit[@bassan2021] supplies mapped CCS "
          "reads without a reference file. CIGAR/MD reconstruction yielded 1,785 bp, matching the SAM header and "
          "agreeing across all covered bases. Identity to the physical ancestral construct was not independently "
          "verified. Sample `11_h12_HA` contained nineteen reads and was excluded; the single baseline, `1_h5_wt`, "
          "belongs to h5 and is not treated as a matched AID-null control for other founders (`crm/tcr_prepare.py`, "
          "`crm/tcr_bassan.py`). Its selected-population read frequencies are endpoints of selection, not mutation "
          "rates. Devilder anti-HLA-A2 cluster sequences44 and published mutation lists were extracted from Additional "
          "file 1, preserving each cluster (`crm/devilder_clones.py`). SHMtrace recovered all 150 listed mutations "
          "across 42 clusters without additional calls: seventeen clusters/thirty-one mutations in round one and "
          "twenty-five clusters/119 mutations in round **two**. This assesses agreement with independently curated "
          "variant identities. No equivalent curated list was available for the TCR libraries. Selected Devilder "
          "sequences support variant comparisons but not population frequencies, unselected distributions or technical "
          "error estimates."),
    ("h2", "Public Illumina Sμ mutation libraries"),
    ("p", "GSE77788[@compagno2017] comprises Sμ libraries from Wt/Aicda−/− primary B cells under DMSO, idelalisib or "
          "duvelisib, three replicates per genotype/treatment. All nine null libraries are paired, whereas six Wt "
          "libraries are 300-bp single-end. Against the bundled 651-bp reference, R1 covers 401–650, single reads "
          "351–650 and R2 1–250. Most pairs do not overlap; merging retains only 9–21%, preferentially short products. "
          "Analyses therefore use R1 alone or the single-end file over shared positions 410–625, excluding primer/end "
          "sequence. Each treatment is analyzed separately against its own null controls at 20,000 reads per library. "
          "The repeat interval, positions 1–327, is reached by R2; a second analysis uses positions 20–240 in the "
          "twelve paired libraries. Its three Wt libraries represent one per treatment, so the pooled genotype "
          "contrast is descriptive rather than a balanced treatment-adjusted design (`crm/public_smu_pi3k.py`)."),
]

QUANT = [
    "Group comparisons use one value per independent animal, subclone, clone, culture or experiment, as specified in "
    "Table S6. Reads are not biological replicates. The prespecified primary test is Welch on log_{10} rates; "
    "negative-binomial and permutation analyses provide sensitivity checks (Figure 3G). Negative-binomial counts k "
    "have mean n·r and variance n·r + φ(n·r)^{2}, where n is observed bases and r the group rate. A common dispersion "
    "φ is fitted by maximum likelihood under the full model, without cross-feature moderation. The common-rate "
    "likelihood-ratio statistic is evaluated using 499 parametric-bootstrap samples with fitted φ; zero-count groups "
    "retain rate zero. Exact permutation inference is used where feasible. At three libraries per group, twenty "
    "allocations yield minimum two-sided p=0.10, explaining zero permutation power at 0.05. This resolution limit does "
    "not establish parametric assumptions. Each Figure 3G condition contains 300 experiments; Wilson 95% half-widths "
    "are at most 5.6 percentage points. At n=3 under the null, rejection rates were Welch 3.7% (95% CI 2.1–6.4), "
    "negative-binomial 5.0% (3.1–8.1), and pooled reads 91.3% (87.6–94.0). Log-rate tests assign zero-call libraries "
    "their one-event value. Geometric-mean folds against these controls depend on that convention and are not "
    "confidence bounds. Spectrum analyses use descriptive Pearson correlation across twelve pooled categories, a "
    "category-based Fisher interval, repeated-measures platform×type testing with mouse as subject, and paired "
    "Šídák-adjusted t tests. The Fisher interval does not represent between-mouse uncertainty. Sample sizes and "
    "dispersion measures are specified in legends. Analyses used Python 3.13, NumPy and SciPy[@virtanen2020].",
]

LEGENDS = [
    ("Figure 1. SHMtrace tests read-library variants against controls and counts Sanger substitutions separately",
     "(A) Quality-based tests can retain recurrent technical mismatches when reported qualities do not reflect their "
     "source. SHMtrace compares observed support with mismatch rates measured by substitution, strand and quality in "
     "matched AID-null libraries. Check marks and crosses indicate correct and incorrect decisions in the schematic. "
     "(B) Read preparation. Sanger traces use Q≥20 trimming instead of step 1; they undergo competitive alignment and "
     "relevant paralog checks but are not downsampled. (C) Nanopore, Illumina and UMI-consensus libraries use a "
     "control-fitted test and empirical false-discovery threshold; a `call` can have one supporting read. Sanger "
     "substitutions at Q≥20 are counted directly, with control-derived expected background reported without "
     "subtraction. The workflows share default high-frequency-allele handling, windows and quality checks. The ≥50% "
     "filter is a heuristic; known-reference experiments have an opt-in exception. (D) Shared outputs include "
     "sample-level comparisons, aggregate single-read estimates, mutation profiles and read-level lineage summaries. "
     "Reads need not identify independent original molecules. FDR, false-discovery rate. See also Figure S1."),
    ("Figure 2. Matched controls distinguish candidate substitutions from recurrent background",
     "(A, B) Two C-to-T candidates in Sμ WT3 (PRJNA1481264; 2,000 reads per library). Bars show support on + and − "
     "strands in WT3 and *Aicda*^{−/−} controls. Counts include all qualities; parentheses show the descriptive Q≥20 "
     "subset. The test uses all qualities. (A) Position 165 in a WRC hotspot. (B) Position 97, a recurrent mismatch "
     "called by the uncorrected Geneious rule in WT3 and **two** controls. p denotes SHMtrace's test. (C) "
     "Supporting-base qualities; bars, medians; dotted line, Q20 for display. (D) Candidates retained at p≤t and false "
     "calls estimated from held-out controls, scaled to tested bases. The largest threshold satisfying the empirical "
     "5% criterion retained 88 calls. (E) Allele fractions among retained SHMtrace calls across JH4/N114 at 500 reads "
     "and Sμ/Pim1 at 2,000. Shading marks the 5–20% range tested with permissive ClairS/DeepSomatic on J_{H}4 (Figure "
     "3H), not a universal caller limit. The dotted 50% line marks a heterozygous germline fraction and the default "
     "exclusion threshold. See also Figure S2."),
    ("Figure 3. Planted mutations and simulated lineages quantify recall and false calls",
     "(A, B) Recall by supporting-read count for substitutions planted into held-out J_{H}4 and N114 controls. Per "
     "locus: 60 sites × 3 draws × 3 controls; each support-count point pools 90 mutations. Lines span held-out "
     "controls. Labels give pooled false/(true+false) calls across nine libraries; complete results are in Table S2. "
     "Comparators are SHMtrace, Geneious and LoFreq with default, BAQ-off or permissive settings. (C) Overall recall "
     "versus false calls per library on JH4/N114 (500 reads) and Sμ/Pim1 (2,000); zero false calls are displayed at "
     "0.05. (D) J_{H}4 lineage simulations: approximately 40 clones, with mean trunk/branch/private loads of 6/2/0.3 "
     "for polyclonal and 1/0.5/0.1 for low-load scenarios. Three 500-read libraries per scenario; lines span "
     "libraries. Text gives null-library false calls. (E) Read- and clone-weighted frequencies versus corresponding "
     "truth, means of three libraries. The expanded clone is unmutated and contributes 60% of reads. Stems connect "
     "each library's truth and estimate. (F) Estimate/truth versus depth, means and ranges across three held-out "
     "controls; J_{H}4 polyclonal, solid; N114 subclonal, dashed. Truth uses planted mutations present in analyzed "
     "reads and each estimator's denominator. Combined density adds the aggregate estimate of uncalled single-read "
     "events. (G) Fraction of 300 experiments yielding two-sided p<0.05, with three or five animals per group. Dotted "
     "line, nominal 5%; Monte Carlo 95% half-width ≤5.6 percentage points. At n=3 per group, twenty allocations give a "
     "minimum exact two-sided permutation p of 0.10. (H) Permissive somatic-caller recall on J_{H}4 at 1–20/500-read "
     "support (solid) and 5%, 10% or 20% support (dashed): thirty sites per fraction per held-out control, with the "
     "other **two** controls pooled as normal. SHMtrace spike-in recall is shown for comparison. BAQ, base-alignment "
     "quality. See also Figure S3 and Tables S2 and S3."),
    ("Figure 4. Matched controls reduce apparent substitution background in real libraries",
     "(A–E) Unique variants per 10^{5} analyzed base observations for callers applied to identical reads; triangles "
     "indicate zero calls. J_{H}4: three Wt/three *Aicda*^{−/−} mice; N114: three AID/three empty-vector subclones, "
     "500 reads each. External Sμ/Pim1: 3–4 AID-competent and three null libraries, 2,000 reads each. (F) "
     "Geometric-mean case/control separation with zero-call libraries assigned one-event densities. OT1/OT2 Wt and "
     "MLLT1/MLLT3 double-knockout groups (three each) are compared separately with their shared null controls. (G) "
     "Calls summed across three controls per assay. (H) Held-out candidate tests with p<0.01 versus error-rate margin; "
     "dotted line, 1%. (I) Delivered-read accounting. Panels use the caller-comparison route: each library's error "
     "model, margin and threshold are calibrated from peer controls. OT1/OT2 Results use the production project's "
     "shared margin/decoys, self-decoy exclusion and noisy-control rule; counts can consequently differ. Both use "
     "empirical control-based thresholding. See also Figures S4–S6 and Tables S1 and S2."),
    ("Figure 5. Matched Sanger and nanopore libraries yield similar mutation profiles",
     "(A) Unique variants per 10^{5} bases per mouse, n=3 per genotype; Sanger, 26–40 clones per Wt mouse; nanopore, "
     "500 reads. Bars, geometric means; p, Welch tests on log densities. Zero-call nanopore controls use one-event "
     "values. Absolute densities are platform-specific endpoints. (B) Pooled twelve-type substitution shares "
     "({n_sanger} Sanger and {n_ont} nanopore calls across three Wt mice), without base-composition adjustment; dashed "
     "line, identity; r, descriptive Pearson correlation. The displayed Fisher interval describes pooled categories "
     "rather than between-mouse uncertainty. (C) Call shares in 25-bp bins across the shared 529-bp interval. (D) G:C "
     "calls at WRC/GYW and A:T calls at WA/TW; bars, three-mouse means; points, mice; dashed lines, window "
     "opportunities. (E) Four displayed substitution classes per mouse and platform; lines join matched mice. The "
     "platform×type test uses all twelve types with mouse as subject. Pooled values generate panel B. See also Figures "
     "S7–S10 and Tables S4 and S5."),
    ("Figure 6. Long reads retain positional profiles and linked mutation patterns",
     "(A) Unique variants per 10^{5} base observations in 50-bp bins over J_{H}4 positions 90–1,230; thin lines, Wt "
     "mice; thick lines, genotype means, n=3 each. Shading marks the Sanger window. (B) Calls per position by "
     "S5F-mutability quintile, relative to each series' mean; Wt calls versus null mismatches, with germline sites "
     "excluded. Targeting indices use permutation tests. (C) Junction-defined lineages among each mouse's 500 reads: "
     "dark, largest; light, others; hatched, unassigned (21–45%). Frequencies are read-weighted and majority-rule "
     "clone-weighted measures. (D) Read-level reconstruction of the largest Wt3 lineage. The square denotes the "
     "unmutated reference; ticks, substitution classes; nodes, linked mutation sets, sized/labeled by supporting "
     "reads; distance, accumulated mutations. Node counts are nested and do not sum. Forty-eight of 69 reads follow "
     "one path; twenty-one contain incompatible patterns consistent with residual error, PCR recombination or "
     "recurrent mutation. These are not proven cellular lineages. See also Figure S11."),
    ("Figure 7. Public datasets test analysis across platforms and biological contexts",
     "(A) Illumina J_{H}4 frequencies in Wt, *Ldha*-deleted and AID-null mice (PRJNA973112). Points, mice; triangles, "
     "zero-call samples; bars, geometric means with one-event floors. Tests: two-sided log-rate Welch against null and "
     "Mann–Whitney for *Ldha* versus Wt. (B) Aicda−/−Ung−/−Msh2−/− CH12F3 replicates expressing AID, R50E, R174E/R178E "
     "or empty vector (PRJNA918596; n=4 each); display as A. (C) Substitution classes and S5F targeting for AID and "
     "R174E/R178E; n denotes calls, not replicates. Mutant calls cluster at positions {ch12_cluster}. (D) Ramos "
     "deletion/insertion events in Wt and HMCES-knockout clones (n=6 each) and AID-knockout controls (n=3), 20,000 "
     "reads per clone (PRJNA820976); bars, means. (E) UMI-consensus error in rRNA operons[@karst2021] by family size; "
     "n, families. Dashed line, raw-read error; diamond, 103 pooled families with a bootstrap 95% interval from 2,000 "
     "family resamples, holding raw error fixed. (F) C:G transition shares in three mOrange experiments across empty, "
     "Ugi and UNG repair contexts (PRJNA170048): circles, SHMtrace; squares, raw Q>20 mismatches; crosses, dead-AID "
     "raw mismatches. Lines join repair arms within experiments. Calling uses empty-arm dead-AID controls; 4–40 called "
     "C:G substitutions per library give coarse shares. See also Table S2."),
]

TABLE_S6_TITLE = ("Table S6. Datasets used to develop and validate SHMtrace, related to STAR Methods")
TABLE_S6_HEADER = ["Dataset", "Platform", "Locus (window)", "Design", "AID-null controls",
                 "Unit of analysis", "Used for"]
TABLE_S6 = [
    ("J_{H}4, this study", "ONT R10.4.1", "J_{H}4 intron (529 bp; 1,141 bp)", "GC B cells, 3 Wt mice",
     "3 *Aicda*^{−/−} mice, same run", "mouse", "Figures 2–6"),
    ("N114, this study", "ONT R10.4.1", "N114 V region (615 bp)", "3 pMX-*Aicda* subclones",
     "3 empty-vector subclones", "subclone\u1d43", "Figures 2–4"),
    ("J_{H}4 Sanger, same mice", "Sanger", "J_{H}4 intron (529 bp)", "26–40 clones per Wt mouse",
     "10–11 clones per *Aicda*^{−/−} mouse", "mouse\u1d47", "Figure 5"),
    ("J_{H}4 Sanger, published study[@barbulescu2024]", "Sanger", "J_{H}4 intron (529 bp)",
     "4 Wt, 4 *Maea*^{fl/fl}*Mb1*^{cre} mice; 20–57 clones each", "3 *Aicda*^{−/−} mice reported, 4 filed", "mouse", "Figure S8"),
    ("J_{H}4 Sanger, published study[@barbulescu2025]", "Sanger", "J_{H}4 intron (529 bp)",
     "3 Wt, 3 *Mkln1*^{−/−} mice; 24–93 clones each", "3 *Aicda*^{−/−} mice", "mouse", "Figure S9"),
    ("PRJNA1481264, Sμ", "ONT", "Sμ (611 bp)", "3 LPS blasts; 4 CH12", "3 *Aicda*^{−/−} each", "culture", "Figures 2 and 4"),
    ("PRJNA1481264, *Pim1*", "ONT", "*Pim1* (862 bp)", "3 LPS blasts", "3 *Aicda*^{−/−}", "culture", "Figures 2 and 4"),
    ("PRJNA1481264, OT1 and OT2", "ONT", "AID off-targets (968 and 942 bp)",
     "3 CH12 wild type, 3 *Mllt1*/*Mllt3* double knockout", "3 *Aicda*^{−/−} each", "culture", "Figure 4"),
    ("PRJNA973112[@sharma2023]", "Illumina", "J_{H}4 (70 bp)", "7 Wt, 7 *Ldha*-deleted mice", "4 AID-null mice",
     "mouse", "Figure 7"),
    ("PRJNA918596[@wang2023]", "Illumina", "VDJ (229 bp)", "AID, R50E, R174E/R178E (4 each)",
     "4 empty-vector replicates", "replicate\u1d9c", "Figure 7"),
    ("PRJNA820976[@wu2022]", "Illumina", "IGH VDJ (421 bp)", "6 Wt, 6 HMCES-KO clones", "3 AID-KO clones",
     "clone", "Figure 7"),
    ("Karst et al.[@karst2021]", "ONT R10, dual UMI", "rRNA operons (4.4 kb)", "103 molecules", "known sequences",
     "molecule", "Figure 7"),
    ("SW_{HEL} Sanger, published study[@thientosapol2018]", "Sanger (GenBank)", "SW_{HEL} V_{H}10 (558 bp)",
     "*Samhd1*, *Msh2* and double mutants, UGI arm; 1,795 clones", "none", "host animal", "Figure S10"),
    ("PRJNA170048[@perezduran2012]", "Illumina (36 bp)", "mOrange reporter (619 bp)",
     "AID-ER in NIH-3T3 cells, 3 experiments x 3 repair contexts", "catalytically dead AID-E58Q-ER",
     "experiment", "Figure 7"),
    ("PRJNA831666[@vincendeau2022]", "PacBio", "Sμ–Sγ1 switch junctions", "3 Wt, 3 *Shld1*, 3 *Xlf*, 3 double",
     "none (repair genotypes, no AID-null arm)", "mouse", "Figure S16"),
    ("GSE77788[@compagno2017]", "Illumina", "Sμ (216 bp shared; 221 bp repeats)",
     "9 Wt, 3 treatments", "9 *Aicda*^{−/−}, same treatments", "culture", "Figure S17"),
    ("SRP159221[@liu2018]", "Illumina", "EGFP reporter (238 bp)", "AID and APOBEC editors, 3 each",
     "3 no-editor libraries", "culture", "Figure S18A"),
    ("SRP075898[@ma2016]", "Illumina", "EGFP reporter and ABL exon 6 (multiplexed)",
     "tethered and untethered editors; 5 sgRNAs; 3–10 sgRNA pools", "1 no-editor library (insufficient)",
     "library", "Figure S18B, S18C"),
    ("figshare 12803501[@bassan2021]", "PacBio", "engineered TCR construct (1,785 bp)",
     "10 diversified and selected populations", "1 parental library of one founder", "population", "Figure S18D"),
    ("Devilder et al.[@devilder2019]", "Sanger (published sequences)", "antibody V_{H} (423 bp)",
     "42 selected clusters, 2 rounds", "parental construct as reference", "cluster", "STAR Methods"),
    ("Simulations", "–", "J_{H}4, N114", "planted sites, lineages, group tests, Sanger projects", "held-out controls",
     "simulated library", "Figures 3 and 5"),
]
TABLE_S6_NOTE = ("GC, germinal-center; Wt, wild type; KO, knockout; ONT, Oxford Nanopore Technologies; UMI, unique molecular identifier. a, independent subclones from one transduction; b, the same mice as the nanopore comparison; c, deposit replicates whose experimental origin is unspecified. The table covers fourteen AID-targeted loci from thirteen laboratories, counting J_{H}4 once across the rows that share it: six somatic hypermutation targets, four AID-targeted loci that are not SHM (Sμ, *Pim1*, OT1 and OT2), and four engineered substrates. The rRNA operons are excluded; the Sμ–Sγ1 junctions are class switching at Sμ and are not counted again. Lengths are analysis windows; the mOrange amplicon is 679 bp. Superscript numbers cite the main manuscript's reference list.")

SUPPLEMENT = [
    "Document S1. Figures S1–S18 and Tables S1 and S3–S6",
    "Table S2. Caller comparison: the five calling settings on real libraries in Figure 4, and planted-mutation results for all tested callers in Figure 3 (Excel)",
    "Data S1. Benchmark outputs (JSON) and SHMtrace project files for every analysis",
]
