"""STAR Methods, key resources table, figure legends, Table 1 and supplemental index for the
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
    ("Mice", "Wt and *Aicda*^{−/−} mice, 4 months of age, were maintained on a predominantly C57BL/6 background. "
             "Both male and female mice were included (n = 3 per genotype). Animals were not stratified by sex "
             "because the study was designed to validate a measurement method rather than to test sex-dependent "
             "variation. The J_{H}4 analysis intentionally used the same mice previously analyzed by "
             "cloned-amplicon Sanger sequencing[@barbulescu2025], enabling a mouse-matched platform comparison. All "
             "animal procedures were approved by the University of Toronto Animal Care Committee (protocol "
             "20011472)."),
    ("Cell lines", "N114 hybridoma cells[@martin2002,@martin2002b] were retrovirally transduced with pMX-*Aicda* or "
                   "the matched empty vector; both constructs encoded GFP and puromycin resistance. GFP-high cells "
                   "were selected and subcloned by limiting dilution, and three puromycin-resistant subclones per "
                   "vector were maintained for 1 month before DNA collection. Each subclone is represented by one "
                   "library. Where a subclone had been sequenced more than once (C1 and D3, from stocks frozen on "
                   "different dates; C3, as two libraries), the library with the most usable reads was analyzed, "
                   "and for C3 the other library is reported as a sensitivity analysis (Results)."),
]

METHODS = [
    ("h2", "Germinal-center B cell isolation and DNA extraction"),
    ("p", "Peyer's patches were harvested post-mortem, dissociated into single-cell suspensions, and stained with "
          "anti-B220, GL7, and anti-Fas. Germinal-center B cells (B220^{+}GL7^{+}Fas^{+}) were isolated by "
          "fluorescence-activated cell sorting (BD FACSymphony S6 SE). Genomic DNA was extracted from "
          "approximately 10,000 sorted germinal-center B cells or 5 × 10^{6} hybridoma cells by proteinase K "
          "digestion."),
    ("h2", "Amplification and sequencing"),
    ("p", "The 1,230-bp J_{H}4 intron amplicon[@jolly1997,@barbulescu2024] and the 655-bp N114 V-region "
          "amplicon[@martin2002] were amplified with SuperFi II high-fidelity DNA polymerase. Cycling was 98 °C for "
          "5 min; 35 cycles of 98 °C for 30 s, 68.8 °C for 30 s, and 72 °C for 1 min; and a final 72 °C for 5 min. "
          "Amplicons were column-purified and submitted as separate, individually barcoded samples for full-length "
          "amplicon sequencing (Plasmidsaurus Premium PCR service), which uses ligation library preparation, "
          "R10.4.1 flow cells, Kit 14 chemistry, and super-accuracy basecalling "
          "(dna_r10.4.1_e8.2_400bps_sup@v5.2.0). Between 2,336 and 5,000 reads were delivered per library. For the "
          "Sanger benchmark, J_{H}4 amplicons from separate PCRs of the same DNA had been cloned and "
          "sequenced[@barbulescu2025], and the 447 resulting chromatograms were re-analyzed here."),
    ("h2", "Read preparation"),
    ("p", "SHMtrace reads FASTQ files, which may be gzipped or, for Illumina data, paired files that are merged on "
          "their overlap. Reads are then processed in the following order. First, reads longer than a locus-specific "
          "limit (2,000 bp for J_{H}4 and 900 bp for N114) are removed, as are reads with a mean expected error "
          "above 0.01 per base, calculated from Phred qualities[@edgar2015]. Second, reads are aligned with "
          "minimap2 (map-ont preset, or sr for Illumina) against an index containing the target and any competing "
          "references, and each read is assigned to the sequence it aligns to best. Unmapped reads, reads assigned "
          "to a competing reference, reads with a supplementary alignment (chimeras), and reads aligned over less "
          "than a minimum length (800 bp for J_{H}4 and 400 bp for N114) are removed. For J_{H}4, the competing "
          "references were a J_{H}4-like paralog (the consensus of Agouti-type clones, approximately 12% "
          "divergent) and a recurrent contaminant haplotype. Because a full-length paralog read can align better "
          "to the longer target, reads carrying 10 or more bases specific to a competing reference, at sites where "
          "fewer than 20% of the sample's reads carry them, are also removed. Genuine J_{H}4 reads carry 0–7 such "
          "bases, whereas paralog reads and chimeras carry 11–62. The fate of every read is recorded. Usable reads "
          "were then downsampled at random with a fixed seed (Python random, seed 0) to the analysis depth: 500 "
          "reads for J_{H}4 and N114 and 2,000 reads for the IRCM loci. The analysis windows were positions "
          "101–629 (529 bp, the window shared with Sanger) or 90–1,230 (full length) of the J_{H}4 reference and "
          "positions 21–635 of the N114 reference."),
    ("h2", "Control-calibrated variant calling"),
    ("p", "Every single-nucleotide difference from the reference within the analysis window is a candidate, and a "
          "single supporting read is sufficient. For each of the 12 substitution types, each strand, and each base "
          "quality, SHMtrace estimates from the AID-null libraries of the same run how often control reads "
          "show that change at or above that quality, combining the pooled rate with each position's own control "
          "count. For a candidate supported by k reads, ranked by quality q_{1} ≥ q_{2} ≥ …, SHMtrace calculates "
          "for each j ≤ k the probability of observing j or more such miscalls at quality ≥ q_{j}, given the read "
          "depth on each strand and the control error rates multiplied by a margin c. The candidate's p value is "
          "the smallest of these probabilities multiplied by k."),
    ("p", "The margin c is chosen by leaving one control out at a time. Each control is tested against an error "
          "model built from the other controls, and c is the smallest value in the grid 1, 2, 3, 5, 10, 20, and "
          "50 at which no more than 1% of a held-out control's tests reach p < 0.01. In the analyses reported "
          "here, c was {loco_margins}. The same held-out tests, made at the samples' analysis depth, are then used "
          "to set the significance threshold. Each held-out control is kept as a decoy, defined by its sorted p "
          "values and the number of bases it tested. For a sample and a candidate threshold t, the expected number "
          "of false calls is the mean, across decoys, of the number of decoy p values ≤ t, scaled by the ratio of "
          "the sample's tested bases to the decoy's when that ratio exceeds 1. Decoys are never scaled down to a "
          "shallower library, whose coarser p values they would otherwise underestimate. SHMtrace keeps the "
          "largest t at which the expected number of false calls is at most 0.05 times the number of sample p "
          "values ≤ t, a target–decoy estimate of the FDR[@elias2007]. A control is never used as its own decoy: "
          "when a control is tested, only the other controls serve as decoys. A control whose held-out candidates "
          "at p ≤ 0.01 per base exceed five times the median of the other controls, and number at least five, is "
          "left out of the decoys (unless no decoy would remain) and named in the report; it still contributes to "
          "the error model. With fewer than two controls, or for indels, the Benjamini–Hochberg "
          "procedure[@benjamini1995] across positions and alternative alleles is used instead. Sites at which 50% "
          "or more of reads differ from the reference are removed as germline or parental differences. Without "
          "controls, SHMtrace falls back to a Poisson-binomial test based on read qualities."),
    ("p", "A unique variant is a position–substitution pair counted once per sample. Unique-variant density is the "
          "number of unique variants per analyzed base × 10^{5}, and mutation frequency is the number of "
          "supporting reads per base sequenced. The Geneious rule was reproduced as an uncorrected approximate "
          "p ≤ 0.05 per candidate, which matched Geneious exports site by site (Jaccard index 0.97–0.98 on Wt "
          "J_{H}4)."),
    ("h2", "Single-read estimate, clones, mutation trees, and quality checks"),
    ("p", "Mutations carried by a single read are estimated in aggregate. SHMtrace counts single-read events at "
          "base quality ≥ 35 that were not otherwise called, subtracts the control rate of such events scaled to "
          "the sample's bases, and divides by the proportion of correctly read bases at that quality. The 95% "
          "interval combines Poisson noise, between-control overdispersion, and the uncertainty of the control "
          "rate."),
    ("p", "Clones are defined by the 20 bases immediately 5′ of the alignment start in reads that begin within "
          "18 bp of J_{H}4. These junction sequences are clustered greedily against abundance-ranked centroids at "
          "an edit distance of ≤ 3, allowing a free 5′ end. A clone's mutations are the calls carried by at least "
          "half of its reads, and the clone-weighted frequency counts each clone once. To build mutation trees, "
          "calls present in two or more reads of a clone and carried by nearly the same reads (Jaccard index "
          "≥ 0.9) are merged into one node, and each node's parent is the smallest node whose reads include at "
          "least 90% of its own. A read fits the tree when its nodes lie on a single path from the root. Two calls "
          "fail the four-gamete test when reads carrying both, either one alone, and neither are all present, each "
          "in at least two reads."),
    ("p", "The targeting index is the mean S5F mutability[@yaari2013] at mutated positions divided by the mean over "
          "the window, with a one-sided permutation p value from 2,000 random sets of positions; the model was "
          "exported from SHazaM 1.2.0. Positional clustering flags any 10-base window containing more calls than a "
          "Poisson distribution allows (Bonferroni correction across windows, p < 10^{−3}). Each stage of the "
          "analysis, the rule it applies and that rule working on the J_{H}4 libraries are shown together in "
          "Figure S10. Sanger chromatograms "
          "are assigned to animals through a sample sheet and trimmed to their longest stretch of base quality "
          "≥ 20 that aligns cleanly, and every substitution at quality ≥ 20 is counted after the germline filter. "
          "Quality-control flags are reported but never applied automatically. They cover copied files, clones "
          "shared between animals beyond chance (Poisson test weighted by mutation frequency across animals), "
          "control clones with excess mutations, label mismatches, and outliers in mutation load."),
    ("p", "A Sanger clone is one molecule read once, so a call cannot be tested by counting reads. What the "
          "AID-null clones can measure is the error itself: the share of their bases at quality \u2265 20 that "
          "mismatch the reference, after their own germline differences are removed, per reported quality and "
          "per position. Each library is then charged its own bases at those rates, giving the calls the "
          "background alone would produce, reported per animal beside its calls; a position at which more than "
          "one control clone mismatches is noted on any call made there. Nothing is removed by this: it is the "
          "scale against which the calls are read, as the control-fitted null is for the long-read path."),
    ("p", "For the published Sanger study[@barbulescu2024], the J_{H}4 chromatograms of the four 96-well plates and "
          "the pilot plate of its second cohort were downloaded from the sequencing facility's portal. Each clone "
          "was assigned to the mouse folder it was filed under in the original analysis project, read from "
          "screen recordings of that project; four traces filed under two mice were analyzed as filed. The "
          "published per-mouse values (sequences, nucleotides, unique mutations and 12-type substitution counts) "
          "were transcribed from the paper's Supplementary Figure 9b and the clone distributions from its Figure "
          "2g. Where a printed table's cells and totals differ, frequencies use the printed totals, which "
          "reproduce the published Figure 2h, and spectra use the cells. The wild-type and *Maea*-deficient "
          "groups were compared with the one-tailed unpaired t-test used in the original study, and with Welch "
          "tests on log frequencies against the *Aicda*^{−/−} mice. The sample sheet and the transcribed "
          "tables are part of Data S1."),
    ("h2", "Benchmark against other callers"),
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
    ("p", "Two somatic callers that compare a case library with a matched normal, ClairS 0.5.1[@zheng2026] and "
          "DeepSomatic 1.9.0[@park2026], and two germline callers, Clair3 2.0.3[@zheng2022] and DeepVariant "
          "1.9.0[@poplin2018], were run from their official Docker images on the same planted mutations, using "
          "the same alignments as LoFreq and restricted to the analysis window. For the somatic callers, the "
          "normal was the pooled usable reads of the other two AID-null libraries, the libraries from which "
          "SHMtrace's null is fitted. Each caller was run at its defaults for R10.4.1 nanopore reads (ClairS "
          "`ont_r10_dorado_sup_5khz`; Clair3 `r1041_e82_400bps_sup_v520`; DeepVariant `ONT_R104`; DeepSomatic "
          "`ONT`) and at the lowest allele fraction it accepts (ClairS `--snv_min_af 0.005 --min_coverage 4`; "
          "Clair3 `--snp_min_af 0.01 --min_coverage 4`; DeepVariant and DeepSomatic candidate thresholds of 1% "
          "and 0.5% of reads, with two supporting reads). PRJNA1481264 reports PromethION sequencing without the "
          "chemistry, so the same R10.4.1 models were used for Sμ and *Pim1*. Substitutions with FILTER PASS "
          "below 0.5 allele fraction were kept, as for the other callers. To find the allele fraction at which "
          "the somatic callers begin to work, 30 substitutions were also planted into 25, 50, and 100 of 500 "
          "reads (5%, 10%, and 20%) of each held-out J_{H}4 control, which otherwise carries few mutations. ClairS "
          "and DeepSomatic were also run on the real AID-active libraries: the wild-type J_{H}4 libraries against "
          "the pooled *Aicda*^{−/−} reads, and the pMX-*Aicda* N114 subclones against the pooled empty-vector reads "
          "(500 reads each, seed 0)."),
    ("h2", "Simulations"),
    ("p", "For lineage simulations, the reads of one AID-null library were divided into approximately 40 "
          "clones with power-law sizes. Each clone received a random 20-bp junction and J trimming (−3 to +12 bp), "
          "and trunk, branch, sub-branch, and private mutations were placed according to the S5F model with "
          "gamma-distributed mutation loads. In the expanded-clone scenario, one unmutated clone made up 60% of "
          "reads. Each control was held out in turn, with the other two as background. Quantification accuracy "
          "was measured with the same design at 125, 250, 500, and 1,000 reads for J_{H}4 and N114 (168 "
          "libraries). Group tests were simulated with log-normal animal rates (coefficient of variation of 1) and "
          "negative-binomial counts on 500 reads × 526 bp per animal, with 300 experiments per condition. Forty "
          "Sanger projects were simulated with realistic trace qualities, unusable traces, and one of each "
          "bookkeeping error planted."),
    ("h2", "Public Illumina and UMI datasets"),
    ("p", "PRJNA973112 (J_{H}4, positions 53–122, 5,000 merged read pairs per mouse), PRJNA918596 (CH12F3 VDJ, "
          "positions 21–249, 5,000 reads per replicate), and PRJNA820976 (Ramos VDJ, positions 26–446, 20,000 "
          "reads per clone, with indel calling) were analyzed with the same analysis engine and Illumina settings, "
          "using the AID-null groups as controls and the S5F model. PRJNA170048 (mOrange reporter, positions "
          "75–693 of the original numbering, 200,000 of the first 600,000 reads per library, about 10,600 reads "
          "per position) was analyzed the same way, with the catalytically dead AID-E58Q-ER libraries of the "
          "empty-vector arm as controls and the dead-AID libraries of the other two arms left as samples; its "
          "reference was rebuilt from the position and reference-nucleotide columns of the study's own "
          "supplementary table, and the same table's per-position counts were used for the comparison with "
          "every-mismatch counting. The CH12F3 and Ramos references were the "
          "majority consensus of the control reads. Because these libraries are sequenced far beyond one read per "
          "molecule, mutation frequency rather than unique-variant density was compared. UMI consensus was "
          "evaluated on 4,259 R10 reads of rRNA operons carrying 18-nt UMIs at both ends[@karst2021]. Reads were "
          "grouped by UMI (allowing up to three edits per end, with directional merging of rare UMIs), and each "
          "group of three or more reads was collapsed into one quality-weighted consensus, which was scored against "
          "the 43 known operon sequences of the mock community."),
]

QUANT = [
    "Groups were compared with one value per animal (or per independent subclone or clone). The negative-binomial "
    "test fits each group's rate with one shared dispersion, as in DESeq2[@love2014], and obtains p values by "
    "parametric bootstrap. The Welch test uses log_{10} rates, with libraries without calls placed at their "
    "one-event density, so that folds are lower bounds. The permutation test is exact. Folds are ratios of "
    "geometric means. Spectra were compared by pooled Pearson correlation across the 12 substitution types with "
    "a Fisher interval, by a repeated-measures test of the platform-by-substitution-type interaction with mouse "
    "as the matching factor, and by paired Šídák-adjusted t tests. Sample sizes, statistical tests, and "
    "dispersion measures are given in each figure legend. Analyses used Python 3.13 with NumPy and "
    "SciPy[@virtanen2020].",
]

LEGENDS = [
    ("Figure 1. SHMtrace calibrates every call against matched AID-null libraries",
     "(A) A read-level test judges each base only by its own quality score, so sequencing errors that recur at one "
     "position (from PCR, homopolymers, or paralog reads) pass the test and are wrongly counted as SHM mutations. "
     "SHMtrace instead measures the error rate of each substitution, strand, and base quality in AID-null libraries "
     "of the same run and calls a mutation only when its support exceeds that background. ✓, correct decision; "
     "✕, wrong decision. "
     "(B) Inputs and read preparation. (C) Calling and automatic checks; numbered steps run in order. (D) Outputs. "
     "FDR, false-discovery rate. See also Figure S1."),
    ("Figure 2. Matched controls separate mutations from recurrent errors that look the same",
     "(A and B) Two C-to-T changes in the same Sμ wild-type library (WT3; PRJNA1481264, 2,000 reads per library): "
     "share of reads carrying each, on the + strand (dark) and − strand (light), in WT3 (blue) and in each "
     "*Aicda*^{−/−} control (red); numbers, reads carrying the change of reads covering the site. (A) A WRC hotspot "
     "(underlined C); (B) a site of recurrent background, which the uncorrected Geneious rule called in WT3 and "
     "in two controls. p, SHMtrace's test against the controls. (C) Base qualities of the reads carrying each "
     "change; bars, medians; dotted line, Q20. (D) The control FDR for the same library: candidates with "
     "p ≤ t (blue) against the false calls that the held-out controls predict at t (red; mean over controls, "
     "scaled to the library's bases); SHMtrace keeps calls up to the largest t at which the predicted false calls "
     "stay below 5% of the calls (dotted line), here the dashed threshold. (E) Share of reads carrying each call "
     "made by SHMtrace in the AID-active libraries of four loci (J_{H}4 and N114, 500 reads; Sμ and *Pim1*, "
     "2,000 reads); shaded, where ClairS and DeepSomatic recall most clean planted mutations (Figure 3H); dotted "
     "line, a heterozygous germline genotype. See also Figure S2."),
    ("Figure 3. Planted mutations and simulated lineages calibrate recall and false calls",
     "(A and B) Recall of mutations planted into reads of a held-out AID-null library (J_{H}4, A; N114, B), by the "
     "number of reads carrying each, for SHMtrace, the uncorrected Geneious rule and LoFreq at three settings; 60 "
     "sites × 3 replicates × 3 held-out controls per amplicon. (C) Overall recall against false calls per library for "
     "every caller and amplicon (colors as in A; shapes, amplicon: J_{H}4 and N114 at 500 reads, Sμ and *Pim1* "
     "from PRJNA1481264 at 2,000 reads); callers without false calls plotted at 0.05. "
     "(D) Recall in simulated B cell lineages by read support (polyclonal, solid; low load, dashed), and false calls "
     "in null libraries; 3 libraries per scenario. (E) Read-weighted and clone-weighted mutation frequency against "
     "the truth (black bars) in simulated libraries; expanded clone, one unmutated clone with 60% of reads. "
     "(F) Estimate over truth against analysis depth for three measures (J_{H}4 polyclonal, solid; N114 subclonal, "
     "dashed); mean of 3 held-out controls. (G) Share of 300 simulated experiments per condition with p < 0.05 for "
     "four group tests (3 or 5 animals per group); dotted line, 0.05. (H) Recall of the somatic callers ClairS "
     "and DeepSomatic on J_{H}4, at their permissive settings, against the share of reads carrying each mutation: "
     "the spike-ins of A (shaded; 1–20 of 500 reads; solid) and 30 mutations planted into 5%, 10%, and 20% of 500 "
     "reads of each held-out AID-null library (dashed, open symbols), with the other two libraries as the matched "
     "normal; SHMtrace on the spike-ins for comparison. "
     "See also Figure S3 and Table S2."),
    ("Figure 4. On real libraries, matched controls set the background",
     "(A–E) Unique variants per 10^{5} analyzed bases in AID-competent (blue) and AID-null (red) libraries for each "
     "caller on the same reads; triangles, libraries without calls. J_{H}4, 3 wild-type and 3 *Aicda*^{−/−} mice; "
     "N114, 3 pMX-*Aicda* and 3 empty-vector subclones (500 reads each); Sμ and *Pim1* (PRJNA1481264), 3–4 "
     "AID-competent and 3 *Aicda*^{−/−} libraries (2,000 reads each). (F) Separation: geometric mean of AID-competent "
     "over AID-null densities, libraries without calls at their one-event density. (G) Calls summed over the three "
     "AID-null libraries. (H) Share of held-out control tests reaching p < 0.01 against the margin on control error "
     "rates; dotted line, the 1% target. (I) Fate of every delivered read in each library. See also Figures S4 "
     "and S5 and Table S2."),
    ("Figure 5. Nanopore and Sanger sequencing agree in the same mice",
     "(A) Unique variants per 10^{5} bases per mouse (n = 3 per genotype) from Sanger clones (26–40 per wild-type "
     "mouse, re-analyzed from chromatograms) and from nanopore reads (500 per mouse); bars, geometric means; p, Welch "
     "test on log densities; nanopore fold is a lower bound (one control without calls at one event). (B) Pooled "
     "share of each of 12 substitution types ({n_sanger} Sanger, {n_ont} nanopore calls, 3 wild-type mice); dashed line, "
     "identity; r, Pearson with 95% CI. (C) Share of calls in 25-bp bins along the 529-bp window. (D) Share of G:C "
     "calls at WRC/GYW and of A:T calls at WA/TW; bars, mean of 3 mice; points, mice; dashed lines, the window's "
     "share of such sites. (E) Planted bookkeeping errors flagged by the Sanger pipeline in 40 simulated projects. "
     "See also Figure S7."),
    ("Figure 6. Long reads add whole-amplicon maps, clones and mutation trees",
     "(A) Unique variants per 10^{5} bases in 50-bp bins over J_{H}4 positions 90–1,230; thin lines, wild-type mice; "
     "thick lines, means (n = 3 per genotype); shading, the window shared with Sanger. (B) Calls per position by "
     "quintile of S5F mutability for wild-type calls and for every *Aicda*^{−/−} mismatch (germline sites excluded), "
     "relative to each series' mean; targeting index, permutation test. (C) Reads with a V(D)J junction by clone "
     "(dark, largest clone) and the read- and clone-weighted mutation frequencies. (D) Mutation tree of the largest "
     "Wt3 clone: square, unmutated ancestor; ticks, mutations colored by class; circles, nodes, sized and labeled "
     "by the reads carrying them; horizontal distance, mutations since the ancestor. See also Figure S9."),
    ("Figure 7. The same engine works across platforms, loci and biological questions",
     "(A) Mutation frequency per mouse in J_{H}4 (Illumina, PRJNA973112): wild type, *Ldha* deleted in activated "
     "B cells, and AID-null; bars, geometric means; p, Welch test on log rates. (B) Mutation frequency per replicate "
     "in CH12F3 *Ung*^{−/−}*Msh2*^{−/−}*Aicda*^{−/−} cells with ectopic AID, the R50E and R174E/R178E mutants and empty "
     "vector (PRJNA918596; n = 4). (C) Spectrum and targeting of AID and R174E/R178E calls; the mutant's calls cluster "
     "at positions {ch12_cluster}. (D) Deletions and insertions per clone in Ramos wild-type and HMCES-deficient clones "
     "(PRJNA820976; n = 6). (E) Error per base of UMI consensus against raw R10 reads by reads per molecule "
     "(rRNA operons[@karst2021]; 103 molecules). (F) Share of transitions among mutations at C:G in NIH-3T3 "
     "cells expressing AID-ER and an mOrange reporter (PRJNA170048), by the repair context set by a third "
     "vector: SHMtrace's calls against every mismatch at quality > 20, as counted in the original study; "
     "shaded, the same share in the catalytically dead AID-E58Q-ER libraries (mean of three experiments each). "
     "See also Table S2."),
]

TABLE1_TITLE = "Table 1. Datasets used to develop and validate SHMtrace"
TABLE1_HEADER = ["Dataset", "Platform", "Locus (window)", "Design", "AID-null controls", "Used for"]
TABLE1 = [
    ("J_{H}4, this study", "ONT R10.4.1", "J_{H}4 intron (529 bp; 1,141 bp)", "GC B cells, 3 Wt mice",
     "3 *Aicda*^{−/−} mice, same run", "Figures 2–6"),
    ("N114, this study", "ONT R10.4.1", "N114 V region (615 bp)", "3 pMX-*Aicda* subclones",
     "3 empty-vector subclones", "Figures 2–4"),
    ("J_{H}4 Sanger, same mice", "Sanger", "J_{H}4 intron (529 bp)", "26–40 clones per Wt mouse",
     "10–11 clones per *Aicda*^{−/−} mouse", "Figure 5"),
    ("J_{H}4 Sanger, published study[@barbulescu2024]", "Sanger", "J_{H}4 intron (529 bp)",
     "4 Wt, 4 *Maea*^{fl/fl}*Mb1*^{cre} mice; 20–57 clones each", "4 *Aicda*^{−/−} mice", "Figure S8"),
    ("PRJNA1481264, Sμ", "ONT", "Sμ (611 bp)", "3 LPS blasts; 4 CH12", "3 *Aicda*^{−/−} each", "Figures 2 and 4"),
    ("PRJNA1481264, *Pim1*", "ONT", "*Pim1* (862 bp)", "3 LPS blasts", "3 *Aicda*^{−/−}", "Figures 2 and 4"),
    ("PRJNA1481264, OT1 and OT2", "ONT", "AID off-targets (1,008 and 982 bp)",
     "3 CH12 wild type, 3 *Mllt1*/*Mllt3* double knockout", "3 *Aicda*^{−/−} each", "Figure 4"),
    ("PRJNA973112[@sharma2023]", "Illumina", "J_{H}4 (70 bp)", "7 Wt, 7 *Ldha*-deleted mice", "4 AID-null mice",
     "Figure 7"),
    ("PRJNA918596[@wang2023]", "Illumina", "VDJ (229 bp)", "AID, R50E, R174E/R178E (4 each)",
     "4 empty-vector replicates", "Figure 7"),
    ("PRJNA820976[@wu2022]", "Illumina", "IGH VDJ (421 bp)", "6 Wt, 6 HMCES-KO clones", "3 AID-KO clones",
     "Figure 7"),
    ("Karst et al.[@karst2021]", "ONT R10, dual UMI", "rRNA operons (4.4 kb)", "103 molecules", "known sequences",
     "Figure 7"),
    ("PRJNA170048[@perezduran2012]", "Illumina (36 bp)", "mOrange reporter (679 bp)",
     "AID-ER in NIH-3T3 cells, 3 experiments x 3 repair contexts", "catalytically dead AID-E58Q-ER",
     "Figure 7"),
    ("Simulations", "–", "J_{H}4, N114", "planted sites, lineages, group tests, Sanger projects", "held-out controls",
     "Figures 3 and 5"),
]
TABLE1_NOTE = ("GC, germinal-center; Wt, wild type; KO, knockout; ONT, Oxford Nanopore Technologies; UMI, unique "
               "molecular identifier.")

SUPPLEMENT = [
    "Document S1. Figures S1–S10 and Tables S1 and S3",
    "Table S2. Per-library calls and densities for every caller and dataset (Excel)",
    "Data S1. Benchmark outputs (JSON) and SHMtrace project files for every analysis",
]
