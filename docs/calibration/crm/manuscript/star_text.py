"""Text of the STAR Protocols article: cloning-free SHM profiling with nanopore sequencing and SHMtrace.

Markup: see docbuild.py. Steps are lists of strings; a step starting with "CRITICAL:", "Note:",
"Pause point:" or "Optional:" is set as such. {lab_...} marks a laboratory detail the authors
must supply (lab_details.py); none is invented here.
"""

TITLE = ("Protocol for cloning-free measurement of somatic hypermutation in mouse germinal-center B cells by "
         "nanopore sequencing and SHMtrace")
AUTHORS = "Abdulmateen Aderinto^{1,3,*} and Alberto Martin^{1,2,4,**}"
AFFIL = ["^{1}Department of Immunology, University of Toronto, Toronto, ON M5S 1A8, Canada",
         "^{2}Lead contact", "^{3}Technical contact",
         "^{4}Senior author",
         "\\*Correspondence: {lab_tech_email}",
         "\\*\\*Correspondence: alberto.martin@utoronto.ca"]
SUMMARY = (
    "Somatic hypermutation is usually measured by cloning amplicons and sequencing individual colonies. Here, "
    "we present a protocol for measuring it without cloning. We describe steps for sorting germinal-center B "
    "cells from mouse Peyer's patches, amplifying the J_{H}4 intron, sequencing full-length amplicons on "
    "nanopore flow cells, and analyzing the reads with SHMtrace, which tests mutations against matched "
    "AID-null controls and compares groups animal by animal.")
SUMMARY_TAIL = ("For complete details on the use and execution of this protocol, please refer to Aderinto and Martin "
                "{lab_crm_citation}.")
HIGHLIGHTS = [
    "Sort germinal-center B cells and amplify the full 1.2-kb J_{H}4 intron without cloning",
    "Sequence full-length amplicons on nanopore flow cells through a service provider",
    "Test mutations against matched AID-null controls with the SHMtrace app or CLI",
    "Compare groups one animal at a time and check the read counts at every step",
]
GA_NOTE = "Graphical abstract: Figure 1 of this protocol, as a separate file."

BEFORE = [
    ("p", "This protocol measures somatic hypermutation (SHM) in the intron downstream of J_{H}4 in germinal-center "
          "(GC) B cells from mouse Peyer's patches[@jolly1997,@sible2021]. At the bench, cloning is replaced by "
          "full-length amplicon sequencing on nanopore flow cells. The analysis uses SHMtrace, open-source software "
          "that tests each candidate mutation against the errors measured in AID-null libraries from the same "
          "sequencing run. We have also used this protocol for a hybridoma V region[@martin2002], and the analysis "
          "applies without changes to other amplicons (for example, a V gene, a switch region, or an AID off-target "
          "gene), to Sanger chromatograms, and to Illumina reads (see the companion article)."),
    ("h", "Institutional permissions"),
    ("p", "Experiments with mice must be approved by the relevant institutional animal care committee. Our "
          "procedures were approved by the University of Toronto Animal Care Committee (protocol 20011472)."),
    ("h", "Design the experiment"),
    ("steps", [
        "Include AID-null controls in every sequencing run. Prepare at least three *Aicda*^{−/−} mice (or other "
        "AID-null samples) and submit them alongside the samples.",
        "CRITICAL: The controls measure the background of the sequencing run. Because errors and basecalling differ "
        "between runs, controls from another run can rank samples but cannot calibrate them.",
        "Choose the number of animals before sequencing. Animals, not reads, limit the sensitivity of SHM "
        "experiments. With the between-animal variation observed in wild-type J_{H}4, three mice per group detect "
        "large differences such as AID dependence, but a two-fold change requires more than ten mice per group. "
        "SHMtrace's power calculation (step 24) estimates this from a pilot experiment.",
        "Plan for 500 usable reads per library. With approximately 5,000 reads delivered per library, 20–56% "
        "remain usable after read preparation (Figure 3B), which leaves a margin.",
    ]),
    ("h", "Prepare the software"),
    ("steps", [
        "Install SHMtrace (Python 3.10 or later; macOS, Linux, or Windows) from the release archive: "
        "`python -m venv shmtrace-env`, then `shmtrace-env/bin/pip install shmtrace-<version>.tar.gz`. "
        "Alternatively, install the desktop application from the release page.",
        "Check the installation. `shmtrace --help` prints the available options, and the test suite (`pytest` in "
        "the unpacked source) passes in a few minutes.",
        "Note: SHMtrace includes the J_{H}4 reference (1,230 bp), its paralog and contaminant references, and "
        "settings for J_{H}4 (`jh4_intron`), a hybridoma V region (`n114_vregion`), and any other amplicon "
        "(`generic_amplicon`).",
    ]),
    ("h", "Check the primers for a co-amplified paralog"),
    ("steps", [
        "Optional: For a new amplicon, sequence one AID-null library first and analyze it with SHMtrace on its "
        "own. If a large proportion of reads is removed at competitive alignment, or if calls appear at ten or "
        "more sites in the same reads, the primers are likely co-amplifying a paralog. Add the paralog's consensus "
        "sequence as a contaminant reference (troubleshooting, problem 1).",
    ]),
]

MATERIALS_TABLES = [
    ("PCR reaction (per sample)", ["Reagent", "Final concentration", "Amount"], [
        ("Platinum SuperFi II PCR master mix (2×)", "1×", "25 µL"),
        ("Forward primer, 10 µM", "0.5 µM", "2.5 µL"),
        ("Reverse primer, 10 µM", "0.5 µM", "2.5 µL"),
        ("Genomic DNA (lysate of ~10,000 cells)", "–", "{lab_lysate_volume}"),
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
    "Computer: any recent laptop with 8 GB of memory is sufficient. A J_{H}4 project of six libraries "
    "(approximately 5,000 reads each) runs in about 1.5 min from the command line on one core of a laptop "
    "(Apple M2 Pro). "
    "Alternatives: other high-fidelity polymerases can be used, provided that the controls are amplified with "
    "the same enzyme. Any provider of full-length nanopore amplicon sequencing, or an in-house ligation kit with "
    "native barcoding, can replace the service used here, provided that reads are delivered as FASTQ files with "
    "base qualities.")

STEPS = [
    ("Isolate germinal-center B cells from Peyer's patches", "Timing: 4–5 h", [
        "Harvest Peyer's patches post-mortem into cold FACS buffer ({lab_facs_buffer}).",
        "Dissociate the Peyer's patches into a single-cell suspension by pressing them through a 70-µm strainer, "
        "and wash the cells once.",
        "Stain the cells with anti-B220, GL7, and anti-Fas for 20–30 min on ice ({lab_stain}), and wash them twice.",
        "Sort B220^{+}GL7^{+}Fas^{+} GC B cells into a tube (we used a BD FACSymphony S6 SE), collecting "
        "approximately 10,000 cells per mouse.",
        "CRITICAL: Process wild-type and AID-null mice on the same day with the same reagents.",
        "Pause point: Pelleted cells can be stored at −80 °C.",
    ]),
    ("Extract genomic DNA", "Timing: overnight", [
        "Lyse the sorted cells by proteinase K digestion, and then inactivate the proteinase K by heating "
        "({lab_lysis}).",
        "Check DNA quality and quantity by NanoDrop (A260/280 and A260/230 ratios).",
        "Pause point: Lysates can be stored at −20 °C.",
    ]),
    ("Amplify the J_{H}4 intron", "Timing: 3 h", [
        "Set up the PCR as described in the PCR reaction table, using the J_{H}4 "
        "primers[@jolly1997,@barbulescu2024] ({lab_jh4_primers}), and include a no-template control.",
        "Run the cycling program described in the PCR cycling table.",
        "Check 5 µL of the product on a 1% agarose gel; a single band of approximately 1.2 kb is expected.",
        "Purify the product with a column PCR clean-up kit and elute it in 15–20 µL of water.",
        "CRITICAL: Amplify wild-type and AID-null samples in the same PCR setup.",
    ]),
    ("Sequence full-length amplicons", "Timing: 15 min hands-on; results in 1–2 days", [
        "Quantify the purified amplicons and submit each one as a separate sample for full-length amplicon "
        "nanopore sequencing, following the provider's concentration and volume requirements. "
        "We used the Plasmidsaurus Premium PCR service (R10.4.1 flow cells, Kit 14, and "
        "super-accuracy basecalling).",
        "Download the FASTQ files (one per library) and keep them unmodified, because SHMtrace works from them "
        "directly.",
        "Note: Approximately 2,300–5,000 reads per library were delivered in our runs.",
    ]),
    ("Set up the SHMtrace project", "Timing: 15 min", [
        "Write a sample sheet (CSV file) with one row per library and the columns `file`, `sample`, `group`, and "
        "`control` (yes for AID-null libraries). For example: "
        "`reads/DHJGKB_5_AID_KO_1.fastq,KO1,Aicda-/-,yes`.",
        "Note: The `control` column may name the kind of control instead of yes: `aid_null`, `empty_vector`, "
        "`plasmid` (a cloned template), or `germline_tissue` (for example, tail DNA), each sequenced on the same run. "
        "A `baseline` library (the cells before the experiment) is used as the control only when no mutation-free "
        "control is given, because it may carry earlier mutations. Do not use a control from another run.",
        "Check the sample sheet with `shmtrace-sheet check samplesheet.csv`, which reports every problem (for "
        "example, missing files, duplicate names, or a sample assigned to two groups) with its line number.",
        "In the application (Figure 2A), choose the region (`jh4_intron`), the platform (`ont`), and the calling "
        "mode (`recommended`). Choose the built-in J_{H}4 reference, which adds the paralog reference "
        "automatically. Then type the path to the sample sheet, click Load sheet, and check that the controls "
        "are ticked.",
        "Optional: To recover more mutations carried by a single read, choose the `sensitive` mode "
        "(`--mode sensitive`), which tests at a false-discovery rate of 0.20 instead of 0.05. About one call in "
        "five is then expected to be background, so state the mode when reporting the results.",
    ]),
    ("Run the analysis", "Timing: 2–10 min", [
        "In the application, name the project and click Analyse. From the command line, run "
        "`shmtrace-sheet run samplesheet.csv --platform ont --region jh4_intron --reference mouse_jh4_intron.fasta "
        "--contaminant mouse_jh4_contaminant.fasta -o results/`.",
        "CRITICAL: Check the platform line that `shmtrace-sheet run` prints first. If `--platform` is not given, "
        "the platform is inferred from the files (nanopore for reads and Sanger for .ab1 chromatograms), and a "
        "sample sheet that mixes both is refused. Give `--platform illumina` for Illumina reads. Reads analyzed "
        "as Sanger clones are counted without a statistical test (troubleshooting, problem 6).",
        "SHMtrace prepares the reads of each library (length and expected-error filters, competitive alignment "
        "against the paralog, and the paralog-allele filter) and draws 500 usable reads at random with a fixed "
        "seed. It then calibrates the background from the controls by leaving one control out at a time, tests "
        "each library, and compares the groups.",
    ]),
    ("Read the report", "Timing: 30 min", [
        "Open the Results tab (Figure 2B), which shows the group comparison with its fold change and p value, "
        "densities per animal, spectra, positional maps, and hotspot use. The file `report.html` in the output "
        "folder contains the same information.",
        "Check the controls first. Each control is tested against the other controls, and approximately "
        "{m_ko_range} calls per control are expected at 500 reads (Figure 3A).",
        "Check the calibration margin (`calibration` in `project.json`), which is the smallest margin at which "
        "no more than 1% of held-out control tests reach p < 0.01. It was {loco_margins_x} in our analyses. The "
        "same held-out tests set the false-discovery threshold (`control_fdr_threshold` in each sample's "
        "calling summary).",
        "Check the read accounting (`qc_waterfall.png`; Figure 3B), which shows how many reads were removed at "
        "each step in each library.",
        "Read the warnings. A dominant clone, positional clustering of calls, too few reads for the target depth, "
        "or a control with excess mutations each produces a message (see troubleshooting).",
        "Open the Tables tab (Figure 2C) or `tables.xlsx`, which list per-sample calls, densities, and "
        "frequencies, every call with its supporting reads, and the clones.",
        "Optional: To continue in Immcantation or other AIRR Community tools, use `airr_rearrangements.tsv`, which "
        "lists every analyzed read with sequencing error removed: each read keeps only the mutations SHMtrace called.",
    ]),
    ("Plan the next experiment", "Timing: 5 min", [
        "Open the Power tab (Figure 2D). Using this analysis as the pilot, choose the measure, the fold change to "
        "detect, and the numbers of animals and reads; the resulting table gives the power of the animal-level "
        "test. From the command line, run "
        "`shmtrace-power --pilot results/project.json --group WT --metric density --fold 0.5 --animals 3,5,8,10`.",
    ]),
    ("Optional: analyze Sanger chromatograms", "Timing: 30 min", [
        "List the .ab1 files in a sample sheet (one row per trace, with the animal as `sample`), check the sheet, "
        "and run `shmtrace-sheet run sheet.csv --platform sanger --region jh4_intron --reference "
        "mouse_jh4_intron.fasta -o sanger/`.",
        "Read `qc_flags.csv`, which flags, but never corrects, copied traces, clones shared between animals, "
        "mutated control clones, and label mismatches. Decide on each flag before interpreting the results.",
    ]),
]

EXPECTED = [
    "With three wild-type and three *Aicda*^{−/−} mice at 500 reads per library, expect approximately "
    "{m_wt_range} unique variants per wild-type library and {m_ko_range} per control, which corresponds to a "
    "group difference of at least {ont_fold}-fold (Figure 3A). The spectrum should include both G:C and A:T "
    "mutations with a bias toward transitions (Figure 3C), and mutations should be concentrated in the first "
    "450 bp of the intron (Figure 3D). Libraries dominated by one clone, which is common in Peyer's patches, have "
    "a low read-weighted mutation frequency; the clone-weighted frequency corrects for this (Figure 4A).",
    "Over the 529 bp shared with Sanger sequencing, nanopore and Sanger data from the same mice agree on the "
    "mutation spectrum (Pearson r = {r_spec}), on the position of mutations (r = {r_pos}), and on hotspot use. "
    "Absolute densities differ by design and should be compared within, not between, platforms (see the "
    "companion article).",
]

QUANT = [
    "For each library, SHMtrace reports the unique-variant density (unique position–substitution variants per "
    "analyzed base × 10^{5}), the mutation frequency (supporting reads per base sequenced), the proportion of "
    "mutated reads, the substitution spectrum, hotspot use, an aggregate single-read estimate, and, for J_{H}4, "
    "clone-weighted frequencies. Groups are compared with one value per animal using a negative-binomial test "
    "with a parametric bootstrap, a Welch test on log rates, and an exact permutation test. Tests that pool reads "
    "across animals are shown only for contrast, because they report a difference in more than 90% of "
    "simulated experiments without one. Compare densities only between libraries analyzed at the same depth. "
    "Controls without calls are placed at their one-event density, so fold changes are lower bounds. Calls are "
    "filtered apparent variants; the exclusion criteria are the read-preparation steps and the germline filter "
    "(sites at which ≥ 50% of reads differ from the reference), and both are reported.",
]

LIMITATIONS = (
    "This protocol requires AID-null controls in the same sequencing run; without them, SHMtrace falls back "
    "to a test based on quality scores, which calls more background. Without unique molecular identifiers "
    "(UMIs), PCR duplicates, early PCR errors, and recombination cannot be separated from true mutations and "
    "clonal expansion. Supporting reads are therefore not independent molecules, and mutation trees are "
    "read-level reconstructions. Mutations carried by a single read are estimated in aggregate rather than "
    "identified individually, and unique-variant density depends on depth. An artifact present in one library "
    "but absent from the controls is flagged by the spectrum, positional clustering, and targeting checks but "
    "is not removed. Weakly targeted loci need more reads (thousands per library for *Pim1*). The thresholds "
    "were established on R10.4.1/Kit 14 super-accuracy reads; other chemistries or basecallers require the same "
    "controls, which SHMtrace re-measures in every run.")

TROUBLE = [
    ("Problem 1", "Few reads are usable, and many are removed at competitive alignment (step 21).",
     "The primers co-amplify a related sequence. In J_{H}4, a paralog approximately 12% divergent made up as much "
     "as a quarter of the delivered reads. Use the included paralog reference. For a new amplicon, build a "
     "consensus of the off-target reads from an AID-null library and provide it with `--contaminant`. If "
     "most reads instead fail the expected-error filter, check the basecalling model and read quality with the "
     "provider."),
    ("Problem 2", "An AID-null control carries calls (step 19).",
     "Examine the control's spectrum and calls. Calls that resemble AID activity (a G:C bias and WRC/GYW "
     "hotspots) suggest contamination or a mislabeled sample, whereas calls with a different spectrum suggest a "
     "library-specific artifact. Because each control serves as a decoy for the false-discovery threshold, a "
     "control far noisier than the others can suppress calls in every sample. With three or more controls, "
     "SHMtrace therefore sets aside such a control as a decoy and names it in the report's notes "
     "(`decoys_set_aside` and `decoy_noise_x1e5` in `calibration`); with only two controls, it cannot tell which "
     "one is unusual. In a public Ramos dataset, one AID-knockout clone carried {hmces_ko3_calls} calls, "
     "{hmces_ko3_gc}% at G:C, whereas the other two control clones had none. With that clone set aside as a "
     "decoy, wild-type and HMCES-deficient clones carried {hm_calls} calls each. Report the analysis with and "
     "without such a control, and resequence it if possible."),
    ("Problem 3", "One clone dominates a library (warning in step 22).",
     "This is common in Peyer's patches. In one wild-type mouse, a nearly unmutated clone made up {wt2_largest}% "
     "of the reads with a V(D)J junction (Figure 4A). Report the clone-weighted frequency alongside the "
     "read-weighted frequency, and use unique-variant density for group comparisons."),
    ("Problem 4", "Calls cluster within a few bases (warning in step 22).",
     "A cluster of calls that recurs across replicates within about 10 bases, does not co-occur within reads, "
     "has a spectrum unlike that of AID, and has a targeting index below 1 is an error cluster (Figure 4B). "
     "Exclude the affected window or libraries and report the exclusion; do not count these calls as mutations."),
    ("Problem 5", "Libraries have too few reads, or controls show calls at low depth (step 16).",
     "Below approximately 250 reads, mutations carried by a single read rarely pass the test and densities are "
     "underestimated, and with very few control reads the background is poorly calibrated. Analyze all "
     "libraries at a common depth of at least 500 usable reads. SHMtrace warns when a library has fewer reads "
     "than the target depth (Figure 4C)."),
    ("Problem 6", "Controls show tens of calls each, and no calibration margin is reported (step 16).",
     "The run was analyzed in Sanger mode, either because `--platform sanger` was given or because an earlier "
     "development version of SHMtrace assumed Sanger data when no platform was given. Run the analysis again "
     "with `--platform ont`, or without `--platform` so that the platform is inferred; the Data line of the "
     "report names the platform used."),
]

RESOURCE = [
    ("Lead contact", "Further information and requests for resources and reagents should be directed to and will "
                     "be fulfilled by the lead contact, Alberto Martin (alberto.martin@utoronto.ca)."),
    ("Technical contact", "Technical questions on executing this protocol should be directed to and will be "
                          "answered by the technical contact, Abdulmateen Aderinto ({lab_tech_email})."),
    ("Materials availability", "This study did not generate new unique reagents."),
    ("Data and code availability", "The nanopore reads, Sanger chromatograms and SHMtrace outputs shown here are "
     "deposited with the companion article ({lab_sra}; {lab_data_doi}). SHMtrace v{shmtrace_version}, its desktop application and "
     "documentation are available at {lab_code_doi} and {lab_github}."),
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
     "(A) Setup: region, platform, calling mode, reference (the paralog reference is added automatically), and "
     "the sample sheet with controls ticked. (B) Results: the group comparison and densities per animal (three "
     "wild-type and three *Aicda*^{−/−} mice, 500 reads each; bars, geometric means). (C) Tables: per-sample "
     "summary. (D) Power: this analysis used as the pilot for the next experiment."),
    ("Figure 3. Expected outcomes",
     "(A) Unique variants per 10^{5} bases per mouse (n = 3 per genotype; 500 reads each); bars, geometric means; "
     "triangle, a control without calls shown at its one-event density; p, Welch test on log densities. (B) Fate "
     "of every delivered read in each library. (C) Pooled proportion of each substitution type among wild-type "
     "mutations. (D) Unique variants along the full amplicon (50-bp bins); thin lines, individual mice; thick "
     "lines, means."),
    ("Figure 4. Troubleshooting examples",
     "(A) A dominant clone: reads with a V(D)J junction by clone in each wild-type mouse, and the read-weighted "
     "versus clone-weighted mutation frequency (problem 3). (B) An error cluster: positions of calls in CH12F3 "
     "libraries with ectopic AID and with an AID mutant whose calls cluster at positions {ch12_cluster} (problem "
     "4; public data, PRJNA918596). (C) Depth: unique variants in wild-type and AID-null libraries at "
     "30–1,000 reads (geometric mean ± geometric SD, n = 3; problem 5)."),
]
