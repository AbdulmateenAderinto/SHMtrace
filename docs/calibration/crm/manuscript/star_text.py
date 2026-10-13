"""Text of the STAR Protocols article: cloning-free SHM profiling with nanopore sequencing and SHMtrace.

Markup: see docbuild.py. Steps are lists of strings; a step starting with "CRITICAL:", "Note:",
"Pause point:" or "Optional:" is set as such. [to confirm] marks a laboratory detail the authors
must supply (dilutions, catalog numbers, primer sequences); none is invented here.
"""

TITLE = ("Protocol for cloning-free measurement of somatic hypermutation in mouse germinal-center B cells by "
         "nanopore sequencing and SHMtrace")
AUTHORS = "Abdulmateen Aderinto^{1,3,*} and Alberto Martin^{1,2,4,**}"
AFFIL = ["^{1}Department of Immunology, University of Toronto, Toronto, ON M5S 1A8, Canada",
         "^{2}Lead contact", "^{3}Technical contact",
         "^{4}Senior author",
         "*Correspondence: [technical contact email, to confirm]",
         "**Correspondence: alberto.martin@utoronto.ca"]
SUMMARY = (
    "Somatic hypermutation is usually measured by cloning amplicons and Sanger sequencing colonies. Here, we "
    "present a protocol to measure it without cloning. We describe sorting germinal-center B cells from mouse "
    "Peyer's patches, amplifying the J_{H}4 intron, sequencing full-length amplicons on nanopore flow cells, and "
    "analyzing the reads with SHMtrace, which calls mutations against matched AID-null controls and compares groups "
    "animal by animal. The analysis also accepts Sanger chromatograms and Illumina reads.")
SUMMARY_TAIL = ("For complete details on the use and execution of this protocol, please refer to Aderinto and Martin "
                "[companion Cell Reports Methods article, to be cited once accepted].")
HIGHLIGHTS = [
    "Sort germinal-center B cells and amplify the full 1.2-kb J_{H}4 intron without cloning",
    "Sequence full-length amplicons on nanopore flow cells through a service provider",
    "Call mutations against matched AID-null controls with SHMtrace, app or command line",
    "Compare groups one animal at a time and check every step's read counts",
]
GA_NOTE = "Graphical abstract: Figure 1 of this protocol, as a separate file."

BEFORE = [
    ("p", "This protocol measures somatic hypermutation (SHM) in the intron downstream of J_{H}4 in germinal-center "
          "(GC) B cells from mouse Peyer's patches[@jolly1997,@sible2021]. The bench part replaces cloning with "
          "full-length amplicon sequencing on nanopore flow cells; the analysis part uses SHMtrace, open-source "
          "software that tests every candidate mutation against the errors measured in AID-null libraries from the "
          "same sequencing run. We have also used it for a hybridoma V region[@martin2002], and the analysis applies "
          "unchanged to other amplicons (a V gene, a switch region, an AID off-target), to Sanger chromatograms and "
          "to Illumina reads (see the companion article)."),
    ("h", "Institutional permissions"),
    ("p", "Experiments with mice must be approved by the relevant institutional animal care committee. Our "
          "procedures were approved by the University of Toronto Animal Care Committee (protocol 20011472)."),
    ("h", "Design the experiment"),
    ("steps", [
        "Include AID-null controls in every sequencing run: at least three *Aicda*^{−/−} mice (or other AID-null "
        "samples) prepared and submitted alongside the samples.",
        "CRITICAL: the controls measure the background of the run. Errors and basecalling differ between runs, so "
        "controls from another run can rank samples but cannot calibrate them.",
        "Choose the number of animals before sequencing. Mice, not reads, limit SHM experiments: with wild-type "
        "variation in J_{H}4, three mice per group detect differences of tens-fold (AID dependence), but a "
        "two-fold change needs more than ten mice per group. SHMtrace's power calculator (step 24) estimates this "
        "from a pilot.",
        "Plan for 500 usable reads per library. With about 5,000 reads delivered per library, 20–56% are usable "
        "after read preparation (Figure 3B), which leaves a margin.",
    ]),
    ("h", "Prepare the software"),
    ("steps", [
        "Install SHMtrace (Python 3.10 or later; macOS, Linux or Windows) from the release archive: "
        "`python -m venv shmtrace-env`, then `shmtrace-env/bin/pip install shmtrace-<version>.tar.gz`. Alternatively, "
        "install the desktop application from the release page.",
        "Check the installation: `shmtrace --help` prints the options; the test suite (`pytest` in the unpacked "
        "source) passes in a few minutes.",
        "Note: SHMtrace bundles the J_{H}4 reference (1,230 bp) and its paralog and contaminant references, and "
        "settings for J_{H}4 (`jh4_intron`), a hybridoma V region (`n114_vregion`) and any other amplicon "
        "(`generic_amplicon`).",
    ]),
    ("h", "Check the primers for a co-amplified paralog"),
    ("steps", [
        "Optional: for a new amplicon, sequence one AID-null library first and run SHMtrace on it alone. A large "
        "share of reads removed at competitive alignment, or calls at ten or more sites in the same reads, "
        "indicates a co-amplified paralog; add its consensus as a contaminant reference (troubleshooting, "
        "problem 1).",
    ]),
]

MATERIALS_TABLES = [
    ("PCR reaction (per sample)", ["Reagent", "Final concentration", "Amount"], [
        ("Platinum SuperFi II PCR master mix (2×)", "1×", "25 µL"),
        ("Forward primer, 10 µM", "0.5 µM", "2.5 µL"),
        ("Reverse primer, 10 µM", "0.5 µM", "2.5 µL"),
        ("Genomic DNA (lysate of ~10,000 cells)", "–", "[volume, to confirm]"),
        ("Nuclease-free water", "–", "to 50 µL"),
        ("**Total**", "–", "**50 µL**"),
    ]),
    ("PCR cycling conditions", ["Steps", "Temperature", "Time", "Cycles"], [
        ("Initial denaturation", "98 °C", "5 min", "1"),
        ("Denaturation", "98 °C", "30 s", "35 cycles"),
        ("Annealing", "68.8 °C", "30 s", ""),
        ("Extension", "72 °C", "1 min", ""),
        ("Final extension", "72 °C", "5 min", "1"),
        ("Hold", "4 °C", "forever", ""),
    ]),
]
MATERIALS_NOTE = (
    "Computer: any recent laptop (8 GB of memory). A J_{H}4 project of six libraries (about 5,000 reads each) runs "
    "in about 1.5 min from the command line on one core of a laptop (Apple M2 Pro). Alternatives: other high-fidelity "
    "polymerases can be used, but the controls must be amplified with the same enzyme; any provider of full-length "
    "nanopore amplicon sequencing (or an in-house ligation kit with native barcoding) can replace the service used "
    "here, provided reads are delivered as FASTQ with base qualities.")

STEPS = [
    ("Isolate germinal-center B cells from Peyer's patches", "Timing: 4–5 h", [
        "Harvest Peyer's patches post-mortem into cold FACS buffer (PBS with 2% fetal bovine serum) [to confirm].",
        "Dissociate them into a single-cell suspension by pressing through a 70-µm strainer; wash once.",
        "Stain with anti-B220, anti-GL7 and anti-Fas for 20–30 min on ice [dilutions and fluorophores, to confirm]; "
        "wash twice.",
        "Sort B220^{+}GL7^{+}Fas^{+} GC B cells (we used a BD FACSymphony S6 SE) into a tube; collect about 10,000 "
        "cells per mouse.",
        "CRITICAL: process wild-type and AID-null mice on the same day with the same reagents.",
        "Pause point: pelleted cells can be stored at −80 °C.",
    ]),
    ("Extract genomic DNA", "Timing: 2 h, or overnight", [
        "Lyse the sorted cells by proteinase K digestion [buffer, volume and incubation, to confirm], then "
        "inactivate the proteinase K by heating.",
        "Check DNA quality and quantity by NanoDrop (A260/280 and A260/230).",
        "Pause point: lysates can be stored at −20 °C.",
    ]),
    ("Amplify the J_{H}4 intron", "Timing: 3 h", [
        "Set up the PCR as in the PCR reaction table with the J_{H}4 primers[@jolly1997,@barbulescu2024] "
        "[sequences, to confirm], including a no-template control.",
        "Run the cycling program in the PCR cycling table.",
        "Check 5 µL on a 1% agarose gel: a single band at about 1.2 kb.",
        "Purify the product with a column PCR clean-up kit and elute in 15–20 µL water.",
        "CRITICAL: amplify wild-type and AID-null samples in the same PCR setup.",
    ]),
    ("Sequence full-length amplicons", "Timing: 15 min hands-on; results overnight to 2 days", [
        "Quantify the purified amplicons and submit each as a separate sample for full-length amplicon nanopore "
        "sequencing (we used Plasmidsaurus Premium PCR: R10.4.1 flow cells, Kit 14, super-accuracy basecalling), "
        "following the provider's concentration and volume requirements [to confirm for your provider].",
        "Download the FASTQ files (one per library) and keep them unmodified; SHMtrace works from them directly.",
        "Note: about 2,300–5,000 reads per library were delivered in our runs.",
    ]),
    ("Set up the SHMtrace project", "Timing: 15 min", [
        "Write a sample sheet (CSV) with one row per library and the columns `file`, `sample`, `group` and "
        "`control` (yes for AID-null libraries), for example: `reads/DHJGKB_5_AID_KO_1.fastq,KO1,Aicda-/-,yes`.",
        "Check it: `shmtrace-sheet check samplesheet.csv` reports every problem (missing files, duplicate names, a "
        "sample given two groups) with its line number.",
        "In the application (Figure 2A): choose the region (`jh4_intron`), the platform (`ont`) and the calling mode "
        "(`recommended`); choose the built-in J_{H}4 reference, which adds the paralog reference automatically; type "
        "the sheet's path and click Load sheet; check that the controls are ticked.",
    ]),
    ("Run the analysis", "Timing: 2–10 min", [
        "In the application, name the project and click Analyse. From the command line: "
        "`shmtrace-sheet run samplesheet.csv --platform ont --region jh4_intron --reference mouse_jh4_intron.fasta "
        "--contaminant mouse_jh4_contaminant.fasta -o results/`.",
        "CRITICAL: check the platform line `shmtrace-sheet run` prints first. Without `--platform`, the platform "
        "is inferred from the files (reads, nanopore; .ab1 chromatograms, Sanger) and a sheet mixing both is "
        "refused; give `--platform illumina` for Illumina reads. Reads analyzed as Sanger clones are counted "
        "without a test (troubleshooting, problem 6).",
        "SHMtrace prepares the reads of each library (length and expected-error gates, competitive alignment "
        "against the paralog, the paralog-allele filter), draws 500 usable reads at random with a fixed seed, "
        "calibrates the background from the controls by leave-one-control-out, calls each library, and compares "
        "the groups.",
    ]),
    ("Read the report", "Timing: 30 min", [
        "Results tab (Figure 2B): the group comparison with fold and p value, densities per animal, spectra, "
        "positional maps and hotspot use. `report.html` in the output folder holds the same.",
        "Check the controls first: each control's calls are made against the other controls; expect {m_ko_range} calls "
        "per control at 500 reads (Figure 3A).",
        "Check the calibration margin (`calibration` in `project.json`): the smallest margin at which at most 1% of "
        "held-out control tests reach p < 0.01; it was {loco_margins_x} in our analyses. The same held-out "
        "calls set the false-discovery threshold (`control_fdr_threshold` in each sample's calling summary).",
        "Check the read accounting (`qc_waterfall.png`; Figure 3B): reads removed at each step per library.",
        "Read the warnings: a dominant clone, positional clustering of calls, too few reads for the target depth, "
        "or a control with excess mutations each produce a message (troubleshooting).",
        "Tables tab (Figure 2C) and `tables.xlsx`: per-sample calls, densities and frequencies, every call with its "
        "support, and the clones.",
    ]),
    ("Plan the next experiment", "Timing: 5 min", [
        "Power tab (Figure 2D): with this analysis as the pilot, choose the measure, the fold change to detect and "
        "the numbers of animals and reads; the table gives the power of the animal-level test. From the command "
        "line: `shmtrace-power --pilot results/project.json --group WT --metric density --fold 0.5 --animals 3,5,8,10`.",
    ]),
    ("Optional: analyze Sanger chromatograms", "Timing: 30 min", [
        "List the .ab1 files in a sample sheet (one row per trace, with the animal as `sample`), check it, and run "
        "`shmtrace-sheet run sheet.csv --platform sanger --region jh4_intron --reference mouse_jh4_intron.fasta "
        "-o sanger/`.",
        "Read `qc_flags.csv`: copied traces, clones shared between animals, mutated control clones and label "
        "mismatches are flagged, never corrected; decide each one before interpreting the results.",
    ]),
]

EXPECTED = [
    "With three wild-type and three *Aicda*^{−/−} mice at 500 reads per library, expect about {m_wt_range} distinct calls "
    "per wild-type library and {m_ko_range} per control, a group difference of at least {ont_fold}-fold (Figure 3A), a spectrum with both "
    "G:C and A:T mutations and transitions favored (Figure 3C), and calls concentrated in the first 450 bp of the "
    "intron (Figure 3D). Libraries dominated by one clone, which is common in Peyer's patches, carry a low "
    "read-weighted frequency; the clone-weighted frequency corrects it (Figure 4A).",
    "Over the 529 bp shared with Sanger sequencing, nanopore and Sanger data from the same mice agree on spectrum "
    "(Pearson r = 0.93), position (r = 0.92) and hotspot use, while absolute densities differ by design (see the "
    "companion article).",
]

QUANT = [
    "SHMtrace reports, per library, the distinct-call density (distinct position–substitution calls per analyzed "
    "base × 10^{5}), the mutation frequency (supporting reads per base sequenced), the share of mutated reads, the "
    "spectrum, hotspot use, an aggregate single-read estimate and, for J_{H}4, clone-weighted frequencies. Groups "
    "are compared with one value per animal: a negative-binomial test with a parametric bootstrap, a Welch test on "
    "log rates and an exact permutation test. Pooled-read tests are shown only as a contrast; they declare a "
    "difference in more than 90% of simulated experiments without one. Compare densities only between libraries "
    "analyzed at the same depth, and controls without calls are placed at their one-event density, so folds are "
    "lower bounds. Calls are filtered apparent variants; exclusion criteria are the read-preparation steps and "
    "the germline rule (sites at ≥ 50% of reads), both reported.",
]

LIMITATIONS = (
    "The protocol needs AID-null controls in the same sequencing run; without them SHMtrace falls back to a "
    "quality-score test with more background. Without unique molecular identifiers (UMIs), PCR duplicates, early PCR "
    "errors and recombination cannot be separated from mutations and clonal expansion; supporting reads are not "
    "independent molecules, and mutation trees are read-level reconstructions. Mutations carried by a single read "
    "are estimated in aggregate, not identified one by one. The distinct-call density depends on depth. An artifact "
    "present in one library and absent from the controls is flagged by the spectrum, positional clustering and "
    "targeting checks but not removed. At weakly targeted loci, more reads are needed (thousands per library at "
    "*Pim1*). The thresholds were set on R10.4.1/Kit 14 super-accuracy reads; other chemistries or basecallers need "
    "the same controls, which SHMtrace re-measures in every run.")

TROUBLE = [
    ("Problem 1", "Few usable reads, many removed at competitive alignment (step 21).",
     "The primers co-amplify a related sequence. In J_{H}4, a paralog about 12% divergent made up to a quarter of "
     "delivered reads. Use the bundled paralog reference; for a new amplicon, build a consensus of the off-target "
     "reads from an AID-null library and pass it with `--contaminant`. If most reads fail the expected-error gate "
     "instead, check the basecalling model and read quality with the provider."),
    ("Problem 2", "An AID-null control carries calls (step 19).",
     "Look at the control's spectrum and calls. Calls that look like AID's (G:C bias, WRC/GYW hotspots) suggest "
     "contamination or a mislabeled sample; calls with another spectrum suggest a library-specific artifact. A control far "
     "noisier than the others also suppresses calls in every sample, because each control serves as a decoy for "
     "the false-discovery threshold; SHMtrace names such a control in the report's notes. In a public Ramos "
     "dataset, one AID-knockout clone carried {hmces_ko3_calls} calls, {hmces_ko3_gc}% at G:C, against none in its "
     "two siblings, and no substitution was called in any wild-type clone until it was analyzed as a sample "
     "({hmk3_calls} calls per clone, {hmk3_gc}% at G:C). Report the analysis with and without such a control, "
     "and resequence if possible."),
    ("Problem 3", "One clone dominates a library (step 22 warning).",
     "Common in Peyer's patches: in one wild-type mouse, a nearly unmutated clone held {wt2_largest}% of reads with a V(D)J "
     "junction (Figure 4A). Report the clone-weighted frequency alongside the read-weighted one, and use the "
     "distinct-call density for group comparisons."),
    ("Problem 4", "Calls cluster in a few bases (step 22 warning).",
     "A cluster of calls recurring across replicates in about 10 bases, not co-occurring within reads, with a "
     "spectrum unlike AID's and a targeting index below 1, is an error cluster (Figure 4B). Exclude the window or the "
     "libraries and say so; do not count it as mutation."),
    ("Problem 5", "Too few reads per library, or controls with calls at low depth (step 16).",
     "Below about 250 reads, single-read mutations rarely pass the test and densities are under-counted; with very "
     "few control reads the background is poorly calibrated. Analyze all libraries at one depth of at least 500 "
     "usable reads; SHMtrace warns when a library has fewer reads than the target depth (Figure 4C)."),
    ("Problem 6", "Controls show tens of calls each and no calibration margin is reported (step 16).",
     "The run was analyzed in Sanger mode, by `--platform sanger` or by an earlier development version of SHMtrace, which "
     "assumed Sanger when no platform was given. Run it again with `--platform ont` (or without `--platform`, to "
     "infer it); the report's Data line names the platform used."),
]

RESOURCE = [
    ("Lead contact", "Further information and requests for resources and reagents should be directed to and will "
                     "be fulfilled by the lead contact, Alberto Martin (alberto.martin@utoronto.ca)."),
    ("Technical contact", "Technical questions on executing this protocol should be directed to and will be "
                          "answered by the technical contact, Abdulmateen Aderinto ([email, to confirm])."),
    ("Materials availability", "This study did not generate new unique reagents."),
    ("Data and code availability", "The nanopore reads, Sanger chromatograms and SHMtrace outputs shown here are "
     "deposited with the companion article ([accessions, to confirm]). SHMtrace v0.2, its desktop application and "
     "documentation are available at [Zenodo DOI and GitHub URL, to confirm]."),
]
ACK = ("We thank the University of Toronto flow cytometry facility and members of the Martin laboratory. "
       "Sequencing was performed by Plasmidsaurus and The Centre for Applied Genomics (SickKids, Toronto). This work "
       "was supported by the Canadian Institutes of Health Research (PJT-180269).")
CONTRIB = ("Conceptualization, A.A. and A.M.; methodology, A.A.; software, A.A.; investigation, A.A.; writing – "
           "original draft, A.A.; writing – review & editing, A.A. and A.M.; supervision, A.M.; funding "
           "acquisition, A.M.")
INTERESTS = "The authors declare no competing interests."
AI_DECL = ("During the preparation of this work, the authors used Anthropic Claude to assist with code drafting, "
           "language editing and figure layout. After using this tool, the authors reviewed and edited the content "
           "as needed and take full responsibility for the content of the published article.")

LEGENDS = [
    ("Figure 1. Protocol overview",
     "Bench steps (top) and analysis steps (bottom), with timing. GC, germinal-center; ONT, Oxford Nanopore "
     "Technologies."),
    ("Figure 2. Set up and run the SHMtrace project",
     "(A) Set up: region, platform, mode, reference (the paralog reference is added automatically) and the sample "
     "sheet with controls ticked. (B) Results: the group comparison and densities per animal (three wild-type and "
     "three *Aicda*^{−/−} mice, 500 reads each; bars, geometric means). (C) Tables: per-sample summary. (D) Power: "
     "this analysis as a pilot for the next experiment."),
    ("Figure 3. Expected outcomes",
     "(A) Distinct calls per 10^{5} bases per mouse (n = 3 per genotype; 500 reads each); bars, geometric means; "
     "triangle, a control without calls at its one-event density; p, Welch test on log densities. (B) The fate of "
     "every delivered read per library. (C) Pooled share of each substitution type in wild-type calls. (D) Distinct "
     "calls along the full amplicon (50-bp bins); thin lines, mice; thick lines, means."),
    ("Figure 4. Troubleshooting examples",
     "(A) A dominant clone: reads with a V(D)J junction by clone in each wild-type mouse, and the read- against "
     "clone-weighted frequency (problem 3). (B) An error cluster: positions of calls in CH12F3 libraries with "
     "ectopic AID and with an AID mutant whose calls cluster at {ch12_cluster} (problem 4; public data, PRJNA918596). "
     "(C) Depth: distinct calls in wild-type and AID-null libraries at 30–1,000 reads (geometric mean ± geometric SD, "
     "n = 3; problem 5)."),
]
