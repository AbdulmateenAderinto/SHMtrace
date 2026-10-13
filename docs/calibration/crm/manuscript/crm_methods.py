"""STAR Methods, key resources table, figure legends, Table 1 and supplemental index for the
Cell Reports Methods article. Markup: see docbuild.py. Bracketed items ([...]) are details only the
authors can supply (catalog numbers, accessions) and are listed for completion in the cover note."""

KRT_HEADER = ["REAGENT or RESOURCE", "SOURCE", "IDENTIFIER"]
KRT = [
    ("**Antibodies**", "", ""),
    ("Anti-mouse B220 (CD45R)", "[supplier]", "[Cat# ; RRID]"),
    ("Anti-mouse GL7", "[supplier]", "[Cat# ; RRID]"),
    ("Anti-mouse Fas (CD95)", "[supplier]", "[Cat# ; RRID]"),
    ("**Chemicals, peptides, and recombinant proteins**", "", ""),
    ("Platinum SuperFi II DNA polymerase", "Thermo Fisher Scientific", "[Cat#]"),
    ("Proteinase K", "[supplier]", "[Cat#]"),
    ("**Critical commercial assays**", "", ""),
    ("Monarch PCR & DNA Cleanup Kit", "New England Biolabs", "[Cat#]"),
    ("Premium PCR (full-length amplicon) sequencing, R10.4.1, Kit 14", "Plasmidsaurus", "https://www.plasmidsaurus.com"),
    ("**Deposited data**", "", ""),
    ("J_{H}4 and N114 nanopore reads, Sanger chromatograms and SHMtrace outputs", "This paper", "[SRA BioProject; Zenodo DOI]"),
    ("Nanopore amplicons, Sμ and *Pim1*, AID-competent and *Aicda*^{−/−}", "NCBI SRA", "BioProject PRJNA1481264"),
    ("Illumina J_{H}4, *Ldha* study", "Sharma et al.[@sharma2023]", "BioProject PRJNA973112; GEO GSE232658"),
    ("Illumina CH12F3 VDJ, ectopic AID", "Wang et al.[@wang2023]", "BioProject PRJNA918596"),
    ("Illumina Ramos IGH VDJ, HMCES", "Wu et al.[@wu2022]", "BioProject PRJNA820976"),
    ("R10 nanopore reads with dual UMIs, ZymoBIOMICS operons", "Karst et al.[@karst2021]", "github.com/SorenKarst/longread_umi"),
    ("**Experimental models: Cell lines**", "", ""),
    ("N114 hybridoma", "Martin et al.[@martin2002]", "N/A"),
    ("**Experimental models: Organisms/strains**", "", ""),
    ("Mouse: C57BL/6 wild type", "[source]", "[RRID]"),
    ("Mouse: *Aicda*^{−/−} (C57BL/6 background)", "Muramatsu et al.[@muramatsu2000]", "[RRID: MGI]"),
    ("**Oligonucleotides**", "", ""),
    ("J_{H}4 intron primers", "Jolly et al.[@jolly1997]", "[sequences in Table S1]"),
    ("N114 V-region primers", "Martin et al.[@martin2002]", "[sequences in Table S1]"),
    ("**Recombinant DNA**", "", ""),
    ("pMX-*Aicda* and pMX empty vector (GFP, puromycin resistance)", "Martin laboratory[@feng2021]", "N/A"),
    ("**Software and algorithms**", "", ""),
    ("SHMtrace v0.2 (command line and desktop application)", "This paper", "[Zenodo DOI]; [GitHub URL]"),
    ("minimap2 2.30 (mappy 2.31)", "Li[@li2018]", "https://github.com/lh3/minimap2"),
    ("LoFreq 2.1.5", "Wilm et al.[@wilm2012]", "https://github.com/CSB5/lofreq"),
    ("SAMtools 1.24", "Danecek et al.[@danecek2021]", "https://github.com/samtools/samtools"),
    ("SHazaM 1.2.0 (HH_S5F model)", "Gupta et al.[@gupta2015]; Yaari et al.[@yaari2013]", "https://immcantation.readthedocs.io"),
    ("Python 3.13, NumPy 2.5, SciPy 1.18", "Virtanen et al.[@virtanen2020]", "https://scipy.org"),
    ("SRA Toolkit 3.4.1", "NCBI", "https://github.com/ncbi/sra-tools"),
    ("Geneious Prime 2026.1.2 (comparison only)", "Kearse et al.[@kearse2012]", "https://www.geneious.com"),
    ("**Other**", "", ""),
    ("R10.4.1 flow cells, Ligation Sequencing Kit V14 (via provider)", "Oxford Nanopore Technologies", "N/A"),
]

EXPERIMENTAL = [
    ("Mice", "Wild-type and *Aicda*^{−/−} mice on a predominantly C57BL/6 background, aged 4 months, of both sexes, "
             "were used (n = 3 per genotype); sex-specific effects were not evaluated because the study validates a "
             "measurement method. These were the same mice analyzed previously by cloned-amplicon Sanger "
             "sequencing[@barbulescu2025]. Procedures were approved by the University of Toronto Animal Care "
             "Committee (protocol 20011472)."),
    ("Cell lines", "N114 hybridoma cells[@martin2002,@martin2002b] were transduced with pMX-*Aicda* or the empty "
                   "vector (both encoding GFP and puromycin resistance), sorted for high GFP and subcloned by limiting "
                   "dilution; three puromycin-resistant subclones per vector were kept for 1 month before DNA "
                   "collection. Each subclone is represented by one library: where a subclone had been sequenced "
                   "more than once (C1 and D3, from stocks frozen on different dates; C3, as two libraries), the "
                   "library with the most usable reads was analyzed, and for C3 the other library is reported as a "
                   "sensitivity analysis (Results)."),
]

METHODS = [
    ("h2", "Germinal-center B cells and DNA"),
    ("p", "Peyer's patches were dissociated and stained with anti-B220, anti-GL7 and anti-Fas, and B220^{+}GL7^{+}Fas^{+} "
          "germinal-center B cells were sorted (BD FACSymphony S6 SE). Genomic DNA was extracted from about 10,000 "
          "sorted cells, or 5 × 10^{6} hybridoma cells, by proteinase K digestion."),
    ("h2", "Amplification and sequencing"),
    ("p", "The 1,230-bp J_{H}4 intron amplicon[@jolly1997,@barbulescu2024] and the 655-bp N114 V region[@martin2002] "
          "were amplified with SuperFi II polymerase (98 °C for 5 min; 35 cycles of 98 °C for 30 s, 68.8 °C for 30 s "
          "and 72 °C for 1 min; 72 °C for 5 min) and column-purified. Each amplicon was submitted as a separate "
          "barcoded sample for full-length amplicon sequencing (Plasmidsaurus Premium PCR: ligation library, "
          "R10.4.1 flow cells, Kit 14, super-accuracy basecalling with dna_r10.4.1_e8.2_400bps_sup@v5.2.0). "
          "2,336–5,000 reads were delivered per library. For the Sanger benchmark, J_{H}4 amplicons from separate PCRs "
          "of the same DNA were cloned and sequenced[@barbulescu2025]; 447 chromatograms were re-analyzed here."),
    ("h2", "Read preparation"),
    ("p", "SHMtrace reads FASTQ (optionally gzipped, or Illumina pairs, merged on overlap) and applies, in order: an "
          "upper read-length limit (2,000 bp for J_{H}4, 900 bp for N114, set per locus), a mean expected error of at "
          "most 0.01 per base computed from Phred qualities[@edgar2015], alignment with minimap2 (map-ont, or sr for "
          "Illumina) against an index holding the target and any competing references, assignment of each read to the "
          "sequence it aligns to best, and removal of unmapped reads, reads assigned to a competing reference, reads "
          "with a supplementary alignment (chimeras) and reads aligned over less than a minimum length (800 bp J_{H}4, "
          "400 bp N114). For J_{H}4, the competing references were a J_{H}4-like paralog (the consensus of "
          "Agouti-type clones, about 12% divergent) and a recurrent contaminant haplotype. Because a full-length "
          "paralog read can align better to the longer target, reads carrying 10 or more bases specific to a "
          "competing reference (at sites where fewer than 20% of the sample's reads carry them) were then removed; real "
          "J_{H}4 reads carry 0–7 such bases and paralog reads and chimeras 11–62. The fate of every read is recorded. "
          "Usable reads were downsampled at random with a fixed seed (Python random, seed 0) to the analysis depth "
          "(500 reads for J_{H}4 and N114; 2,000 for the IRCM loci). Analysis windows were positions 101–629 (529 bp, "
          "the window shared with Sanger) or 90–1,230 (full length) of the J_{H}4 reference and positions 21–635 of "
          "the N114 reference."),
    ("h2", "Control-calibrated calling"),
    ("p", "Every single-nucleotide difference from the reference in the window is a candidate, one supporting read "
          "sufficing. For each substitution class (12), strand (2) and base quality, the rate at which control reads "
          "show that change at or above that quality is estimated from the AID-null libraries of the same run, pooled "
          "with each position's own control count. For a candidate with k supporting reads sorted by quality "
          "q_{1} ≥ q_{2} ≥ …, SHMtrace computes for each j ≤ k the probability of observing j or more such miscalls "
          "at quality ≥ q_{j}, given the column depth on each strand and the control rates multiplied by a margin c; "
          "the candidate's p value is the smallest over j, multiplied by k. The margin is chosen by "
          "leave-one-control-out from the grid 1, 2, 3, 5, 10, 20 and 50 as the smallest value at which at most 1% of "
          "a held-out control's tests reach p < 0.01; it was {loco_margins} in the analyses reported. Each held-out "
          "control, downsampled to the samples' analysis depth and called against the other controls, is kept as a "
          "decoy: its sorted p values and the bases it tested. For a sample with p values P and a threshold t, the "
          "expected number of false calls is the mean over decoys of the decoy p values ≤ t, scaled by the ratio of "
          "the sample's tested bases to the decoy's when that exceeds 1 (decoys are never scaled down to a shallower "
          "library, whose coarser p values they would otherwise under-represent); SHMtrace keeps the largest t at which this is at most 0.05 times "
          "the number of sample p values ≤ t (target–decoy FDR[@elias2007]). A control is never its own decoy: when a control is called, only "
          "its peers' decoys are used. With fewer than two controls, or for indels, the Benjamini–Hochberg "
          "procedure[@benjamini1995] across positions and alternative alleles is used instead. Sites at which 50% or "
          "more of reads differ from the reference are removed. Without controls, SHMtrace falls back to a Poisson-binomial test on "
          "read qualities. A distinct call is a position–substitution pair counted once per sample; the "
          "distinct-call density is calls per analyzed base × 10^{5}, and the mutation frequency is supporting "
          "reads per base sequenced. The Geneious rule was reproduced as an uncorrected approximate p ≤ 0.05 per "
          "candidate, which matched Geneious exports site by site (Jaccard index 0.97–0.98 on wild-type J_{H}4)."),
    ("h2", "Single-read estimate, clones, trees and checks"),
    ("p", "Mutations carried by one read are estimated in aggregate: one-read events at base quality ≥ 35 not "
          "otherwise called, minus the control rate of such events scaled to the sample's bases, divided by the share "
          "of correctly read bases at that quality, with a 95% interval combining Poisson noise, between-control "
          "overdispersion and the uncertainty of the control rate. Clones are defined by the 20 bases 5′ of the "
          "alignment start of reads beginning within 18 bp of J_{H}4, clustered greedily against abundance-ranked "
          "centroids at edit distance ≤ 3 with a free 5′ edge; a clone's mutations are calls carried by at least half "
          "of its reads, and the clone-weighted frequency counts each clone once. Within a clone, calls in two or more "
          "reads carried by nearly the same reads (Jaccard ≥ 0.9) form one node, each node's parent is the smallest "
          "node whose reads include at least 90% of its own, a read fits the tree when its nodes lie on one path from "
          "the root, and two calls fail the four-gamete test when reads carrying both, either alone and neither are "
          "all present (each in at least two reads). The targeting index is the mean S5F mutability[@yaari2013] at "
          "mutated positions over the window mean, with a one-sided permutation p value (2,000 random position sets); "
          "the model was exported from SHazaM 1.2.0. Positional clustering flags any 10-base window holding more calls "
          "than a Poisson scatter allows (Bonferroni over windows, p < 10^{−3}). Sanger chromatograms are assigned to "
          "animals by a sample sheet, trimmed to their longest stretch of base quality ≥ 20 that aligns cleanly, and "
          "every substitution at quality ≥ 20 is counted after the germline rule; quality-control flags (copied files, "
          "clones shared between animals beyond chance by a Poisson test weighted by mutation frequency across animals, "
          "control clones with excess mutations, label mismatches, load outliers) are reported, never applied."),
    ("h2", "Benchmark against LoFreq and the Geneious rule"),
    ("p", "For each library, the reads SHMtrace analyzed were written out and aligned with minimap2 (same preset, "
          "secondary alignments off), sorted and indexed with SAMtools, and called with LoFreq 2.1.5 `call` restricted "
          "to the analysis window at three settings: defaults (base-alignment quality on, Bonferroni correction, "
          "default filters); `-B` (base-alignment quality off); and `-B -a 0.05 -b 1 --no-default-filter` "
          "(uncorrected p ≤ 0.05, filters off). Substitutions at allele frequency < 0.5 were kept, matching "
          "SHMtrace's germline rule. All callers share one density denominator per library. In spike-in tests, 60 "
          "sites per replicate were chosen with WRC/GYW weighted 4× and WA/TW 2×, and each substitution was written "
          "into 1, 2, 3, 5, 10 or 20 reads covering it, keeping the read's own base quality; three replicates per "
          "held-out control. The IRCM runs (PRJNA1481264) were downloaded with the SRA Toolkit, validated with "
          "vdb-validate and subsampled to 6,000 reads (seed 0); each locus' reference was the majority consensus of "
          "its first *Aicda*^{−/−} library (minimap2 and SAMtools consensus, iterated to convergence)."),
    ("h2", "Simulations"),
    ("p", "Lineage simulations dealt the reads of one AID-null library into about 40 clones with power-law sizes, "
          "gave each a random 20-bp junction and J trimming (−3 to +12 bp), and planted trunk, branch, sub-branch and "
          "private mutations placed by the S5F model, with gamma-distributed loads; in the expanded scenario one "
          "unmutated clone held 60% of reads. Each control was held out in turn with the other two as background. "
          "Quantification accuracy was measured on the same design at 125, 250, 500 and 1,000 reads for J_{H}4 and "
          "N114 (168 libraries). Group tests were simulated as log-normal animal rates (coefficient of variation 1) "
          "with negative-binomial counts on 500 reads × 526 bp per animal, 300 experiments per condition. Sanger "
          "projects (40) were simulated with realistic trace qualities, unusable traces and one of each bookkeeping "
          "error planted."),
    ("h2", "Public Illumina and UMI datasets"),
    ("p", "PRJNA973112 (J_{H}4, positions 53–122, 5,000 merged pairs per mouse), PRJNA918596 (CH12F3 VDJ, positions "
          "21–249, 5,000 reads per replicate) and PRJNA820976 (Ramos VDJ, positions 26–446, 20,000 reads per clone, "
          "indels on) were analyzed with the same engine and Illumina settings, the AID-null groups as controls, and "
          "the S5F model; the CH12F3 and Ramos references were majority consensuses of the control reads. These "
          "libraries are resequenced far beyond one read per molecule, so mutation frequency, not distinct-call "
          "density, was compared. UMI consensus was evaluated on 4,259 R10 reads of rRNA operons with 18-nt UMIs at "
          "both ends[@karst2021]: reads were grouped by UMI (up to three edits per end, directional merging of rare "
          "UMIs) and each group of three or more reads became one quality-weighted consensus, scored against the "
          "community's 43 known operon sequences."),
]

QUANT = [
    "Groups were compared with one value per animal (or independent subclone or clone). The negative-binomial test "
    "fits each group's rate with one shared dispersion, as in DESeq2[@love2014], and obtains p by parametric bootstrap; "
    "the Welch test uses log_{10} rates, with libraries without calls placed at their one-event density so that folds "
    "are lower bounds; the permutation test is exact. Folds are ratios of geometric means. Spectra were compared by "
    "pooled Pearson correlation across the 12 substitution types with a Fisher interval, a repeated-measures "
    "platform-by-type interaction with mouse as the matching factor and paired Šídák-adjusted t tests. Sample sizes, "
    "tests and dispersion measures are given in each figure legend. Analyses used Python 3.13 with NumPy and "
    "SciPy[@virtanen2020].",
]

LEGENDS = [
    ("Figure 1. SHMtrace calibrates every call against matched AID-null libraries",
     "(A) A read-level test judges each base by its own quality score, so errors that recur at one position pass as "
     "mutations; SHMtrace measures the error rate of each substitution, strand and quality in AID-null libraries of "
     "the same run and calls only support beyond it. ✓, correct decision; ✕, wrong decision. "
     "(B) Inputs and read preparation. (C) Calling and automatic checks; numbered steps run in order. (D) Outputs. "
     "FDR, false-discovery rate. See also Figure S1."),
    ("Figure 2. Planted mutations and simulated lineages calibrate recall and false calls",
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
     "four group tests; dotted line, 0.05. See also Figure S5."),
    ("Figure 3. On real libraries, matched controls set the background",
     "(A–E) Distinct calls per 10^{5} analyzed bases in AID-competent (blue) and AID-null (red) libraries for each "
     "caller on the same reads; triangles, libraries without calls. J_{H}4, 3 wild-type and 3 *Aicda*^{−/−} mice; "
     "N114, 3 pMX-*Aicda* and 3 empty-vector subclones (500 reads each); Sμ and *Pim1* (PRJNA1481264), 3–4 "
     "AID-competent and 3 *Aicda*^{−/−} libraries (2,000 reads each). (F) Separation: geometric mean of AID-competent "
     "over AID-null densities, libraries without calls at their one-event density. (G) Calls summed over the three "
     "AID-null libraries. (H) Share of held-out control tests reaching p < 0.01 against the margin on control error "
     "rates; dotted line, the 1% target. (I) Fate of every delivered read in each library. See also Figures S2 "
     "and S3 and Table S2."),
    ("Figure 4. Nanopore and Sanger sequencing agree in the same mice",
     "(A) Distinct calls per 10^{5} bases per mouse (n = 3 per genotype) from Sanger clones (26–40 per wild-type "
     "mouse, re-analyzed from chromatograms) and from nanopore reads (500 per mouse); bars, geometric means; p, Welch "
     "test on log densities; nanopore fold is a lower bound (one control without calls at one event). (B) Pooled "
     "share of each of 12 substitution types ({n_sanger} Sanger, {n_ont} nanopore calls, 3 wild-type mice); dashed line, "
     "identity; r, Pearson with 95% CI. (C) Share of calls in 25-bp bins along the 529-bp window. (D) Share of G:C "
     "calls at WRC/GYW and of A:T calls at WA/TW; bars, mean of 3 mice; points, mice; dashed lines, the window's "
     "share of such sites. (E) Planted bookkeeping errors flagged by the Sanger pipeline in 40 simulated projects. "
     "See also Figure S4."),
    ("Figure 5. Long reads add whole-amplicon maps, clones and mutation trees",
     "(A) Distinct calls per 10^{5} bases in 50-bp bins over J_{H}4 positions 90–1,230; thin lines, wild-type mice; "
     "thick lines, means (n = 3 per genotype); shading, the window shared with Sanger. (B) Calls per position by "
     "quintile of S5F mutability for wild-type calls and for every *Aicda*^{−/−} mismatch (germline sites excluded), "
     "relative to each series' mean; targeting index, permutation test. (C) Reads with a V(D)J junction by clone "
     "(dark, largest clone) and the read- and clone-weighted mutation frequencies. (D) Mutation tree of the largest "
     "Wt3 clone: square, unmutated ancestor; ticks, mutations colored by class; circles, nodes, sized and labeled "
     "by the reads carrying them; horizontal distance, mutations since the ancestor. See also Figure S6."),
    ("Figure 6. The same engine works across platforms, loci and biological questions",
     "(A) Mutation frequency per mouse in J_{H}4 (Illumina, PRJNA973112): wild type, *Ldha* deleted in activated "
     "B cells, and AID-null; bars, geometric means; p, Welch test on log rates. (B) Mutation frequency per replicate "
     "in CH12F3 *Ung*^{−/−}*Msh2*^{−/−}*Aicda*^{−/−} cells with ectopic AID, the R50E and R174E/R178E mutants and empty "
     "vector (PRJNA918596; n = 4). (C) Spectrum and targeting of AID and R174E/R178E calls; the mutant's calls cluster "
     "at positions {ch12_cluster}. (D) Deletions and insertions per clone in Ramos wild-type and HMCES-deficient clones "
     "(PRJNA820976; n = 6). (E) Error per base of UMI consensus against raw R10 reads by reads per molecule "
     "(rRNA operons[@karst2021]; 103 molecules). See also Table S2."),
]

TABLE1_TITLE = "Table 1. Datasets used to develop and validate SHMtrace"
TABLE1_HEADER = ["Dataset", "Platform", "Locus (window)", "Design", "AID-null controls", "Used for"]
TABLE1 = [
    ("J_{H}4, this study", "ONT R10.4.1", "J_{H}4 intron (529 bp; 1,141 bp)", "GC B cells, 3 Wt mice",
     "3 *Aicda*^{−/−} mice, same run", "Figures 2–5"),
    ("N114, this study", "ONT R10.4.1", "N114 V region (615 bp)", "3 pMX-*Aicda* subclones",
     "3 empty-vector subclones", "Figures 2 and 3"),
    ("J_{H}4 Sanger, same mice", "Sanger", "J_{H}4 intron (529 bp)", "26–40 clones per Wt mouse",
     "10–11 clones per *Aicda*^{−/−} mouse", "Figure 4"),
    ("PRJNA1481264, Sμ", "ONT", "Sμ (611 bp)", "3 LPS blasts; 4 CH12", "3 *Aicda*^{−/−} each", "Figure 3"),
    ("PRJNA1481264, *Pim1*", "ONT", "*Pim1* (862 bp)", "3 LPS blasts", "3 *Aicda*^{−/−}", "Figure 3"),
    ("PRJNA973112[@sharma2023]", "Illumina", "J_{H}4 (70 bp)", "7 Wt, 7 *Ldha*-deleted mice", "4 AID-null mice",
     "Figure 6"),
    ("PRJNA918596[@wang2023]", "Illumina", "VDJ (229 bp)", "AID, R50E, R174E/R178E (4 each)",
     "4 empty-vector replicates", "Figure 6"),
    ("PRJNA820976[@wu2022]", "Illumina", "IGH VDJ (421 bp)", "6 Wt, 6 HMCES-KO clones", "3 AID-KO clones",
     "Figure 6"),
    ("Karst et al.[@karst2021]", "ONT R10, dual UMI", "rRNA operons (4.4 kb)", "103 molecules", "known sequences",
     "Figure 6"),
    ("Simulations", "–", "J_{H}4, N114", "planted sites, lineages, group tests, Sanger projects", "held-out controls",
     "Figures 2 and 4"),
]
TABLE1_NOTE = ("GC, germinal-center; Wt, wild type; KO, knockout; ONT, Oxford Nanopore Technologies; UMI, unique "
               "molecular identifier.")

SUPPLEMENT = [
    "Document S1. Figures S1–S6 and Tables S1 and S3",
    "Table S2. Per-library calls and densities for every caller and dataset (Excel)",
    "Data S1. Benchmark outputs (JSON) and SHMtrace project files for every analysis",
]
