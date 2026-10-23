"""Step 7f of revise_operations.py: restore the flow of the main text after the Codex compression.

Each entry replaces a whole paragraph, found by its opening words. The aim is connected prose: transitions,
reasons and the order of the argument restored, mostly from the authors' v1.3 wording, while keeping every later
correction (citation keys, FDP terminology, review fixes H2-H8, the ablation sentence). No number is changed; the
check in revise_operations.py fails if a number appears that is not already in v1.3 or the current text. Figure
captions are tightened to pay for the restored words; their detail remains in STAR Methods and Document S1.
"""

PARAS = {
# ───────── Results: calibration ─────────
"SHMtrace starts from delivered reads":
"SHMtrace starts from the files that a sequencing provider or core facility returns (Figures 1 and S1). Nanopore and "
"Illumina reads are first filtered by length and expected error[@edgar2015] and then aligned with minimap2[@li2018] "
"against the target together with competing references. Reads that align better to a paralog, are chimeric or carry "
"paralog-specific bases are removed, and the same number of usable reads is then drawn from every library. Sanger "
"chromatograms, by contrast, are assigned to animals through a sample sheet and trimmed to their high-quality core, and "
"UMI-tagged reads can be collapsed to one consensus per molecule (STAR Methods).",

"For each substitution, SHMtrace tests read support":
"The central step is the variant test (Figures 1A and 1C). For each substitution type, strand and base quality, "
"SHMtrace measures how often the AID-null libraries of the same run show that change, and then asks whether a candidate "
"is supported by more reads than this background explains. Because quality scores overstate certainty at the highest "
"qualities (Figure S2), the control rates are multiplied by a margin chosen by leaving one control out at a time. The "
"held-out controls then serve as decoys, in an approach inspired by target–decoy analysis[@elias2007], and SHMtrace "
"keeps the most permissive cutoff at which estimated false calls are at most 5% of calls. This is an empirical estimate "
"rather than a guarantee, since it relies on a few controls representing the samples' background. We preferred it to "
"Benjamini–Hochberg correction[@benjamini1995] because candidates arise only where a read already differs, so the "
"tested set is itself chosen by the data. Finally, alleles carried by at least half the reads are set aside as likely "
"germline or parental differences; because a mutation shared by a dominant lineage can cross the same threshold, this "
"is a heuristic, and every excluded site is reported.",

"The read-library test can retain a candidate":
"Because it weighs read support against measured background, the read-library test can retain a candidate supported by "
"a single read; repeated sequencing of the same original molecule is not required. Sanger clones instead take a "
"separate counting route (Figure 1C): every substitution at base quality ≥ 20 is counted, and the AID-null clones' "
"mismatch rate gives the expected background, reported beside the counts rather than subtracted. The germline filter, "
"window, automatic checks and animal-level comparisons are shared. The Sanger results below therefore test trace "
"reading, counting, bookkeeping and group comparisons, not the control-calibrated test.",

"Two C-to-T candidates in Sμ WT3":
"Two C-to-T candidates in an external wild-type Sμ library (WT3) illustrate how matched controls separate mutations "
"from recurrent background (Figure 2). Both had support on both strands and similar base qualities. However, no "
"control supported the substitution at position 165, within a WRC hotspot, at Q≥20, whereas the substitution at "
"position 97 recurred in the controls. SHMtrace therefore retained position 165 (p = 5 × 10^{−19}) and rejected "
"position 97 (p = 1). Across WT3, 88 of 743 candidates met the threshold (p ≤ 0.009); the remainder could not be "
"distinguished confidently from background.",

"Group comparisons use one value per biological replicate":
"Downstream, SHMtrace compares groups using one summary value per biological replicate, with Welch tests on "
"log-transformed rates, bootstrapped negative-binomial tests and permutation tests (Figure 1D). It also describes "
"spectra, hotspot use, S5F targeting[@yaari2013] and clustering, estimates uncalled single-read mutations in aggregate, "
"groups reads into junction-defined clones and reconstructs read-level trees (Figure S3). The desktop and command-line "
"interfaces share one analysis engine.",

# ───────── Results: planted mutations ─────────
"ClairS[@zheng2026] and DeepSomatic[@park2026] received":
"We also gave the long-read somatic callers ClairS[@zheng2026] and DeepSomatic[@park2026] the other AID-null libraries "
"as matched normals. Neither recovered planted mutations at default settings, and even with permissive settings ClairS "
"recovered none of 2,160 substitutions, while DeepSomatic recalled 19%, 18%, 5% and 3% on J_{H}4, N114 V_{H}, Sμ and "
"*Pim1*, respectively. Their recall improved at higher allele fractions (Figure 3C), but real SHM mostly lies below "
"that range: 53% of 1,197 retained SHMtrace calls occurred in <1% of reads (Figures 2E and 3C). The germline callers "
"Clair3 and DeepVariant recovered none, consistent with their intended use (Table S1).",

"On N114 V_{H}, SHMtrace recalled 90%":
"On N114 V_{H}, SHMtrace recalled 90% of planted substitutions, including 58% supported by one read, with 2.6 false "
"calls per library (FDP 0.045; Figure 3B). Here Geneious performed similarly, recalling 88% with 2.8 false calls (FDP "
"0.050), whereas LoFreq without base-alignment quality recalled 34%. On external Sμ and *Pim1* controls at 2,000 reads, "
"SHMtrace recalled 63% and 72% at FDP 0.055 and 0.035, compared with 59% and 64% for Geneious at FDP 0.114 and 0.050; "
"LoFreq recalled 0–26%, with permissive FDP 0.27–0.52 (Figure 3D). To isolate what the controls contribute, two "
"ablations used the same libraries. Without controls, a read-quality test made 5.7–46.6 false calls per library; "
"testing each site against its control error rate instead made 1.3–2.8, but recalled only 27–64%, against SHMtrace's "
"63–90% (Table S1). Measured background therefore removes most false calls, and the calibrated test recovers the "
"sensitivity a simple per-site comparison loses.",

"Calling the same libraries at targets":
"We then called the same libraries at targets from 1% to 30% to test whether the threshold delivers the rate it "
"targets. At the default 5%, the pooled false-discovery proportion was 3.5–5.5% across the four loci, and looser "
"targets stayed at or below their stated rate except J_{H}4 at 10% (13.8%) (Figure S4A). Targets below 5%, however, "
"did not deliver their stated rate on J_{H}4, Sμ or *Pim1* (2.0–6.1% at a 1% target; N114 V_{H}, 0.5%): false calls "
"per library levelled off while recall fell, from 81% to 17% on J_{H}4 (Figures S4B and S4C).",

"Varying the controls through the whole procedure":
"Because everything rests on the controls, we next varied them through the whole procedure on J_{H}4 (Figure S5). One "
"control leaves no decoys, so the threshold falls back to Benjamini–Hochberg; it recalled 68% of planted lineage "
"mutations, against 93% with two controls and 95% with three or four, with false-discovery proportions of 2.8–3.1% "
"for the two- to four-control designs. All-null libraries yielded no calls with any number of controls (40 runs), and "
"one control carrying recurrent artifacts changed false calls little, with or without the noisy-control rule. On Sμ, "
"AID-null controls from a second experiment performed like the matched controls (66–68% versus 63% recall), and "
"pooling all five lowered false calls from 2.2 to 0.8 per library.",

"Diluting real hypermutated reads":
"Planted mutations lack the clonality and clustering of real SHM, so we also diluted real hypermutated reads into "
"AID-null reads (Figure S6). On N114 V_{H}, mutation frequency was 0.87–0.97 of the expected value at 25–50% donor "
"reads, and calls were made from 1–5% donor reads. On J_{H}4, in contrast, detection was all-or-nothing and depended "
"on the donor: calls began at 2% donor reads for Wt3 and at 50% for Wt2, with no call below that.",

"Indel validation used ten planted":
"Indels were validated in the same way, with ten planted insertions or deletions of 1–40 bp per library, supported by "
"1–20 reads, eight draws per control and error models fitted to the other controls. Among events supported by at "
"least three reads, recall was 72% on J_{H}4 and 91% on N114 V_{H}, with 0.33 and zero false calls per library; events "
"supported by one or two reads fall below the three-read requirement. Single-base homopolymer deletions in J_{H}4 were "
"the least reliably recovered (19%), whereas recall was 100% for 2-bp insertions and 89% for deletions ≥10 bp "
"(Figure S7).",

"The dominant unmutated clone could not trigger":
"Because that dominant clone was unmutated, it could not trigger the ≥50% allele filter. When we repeated the "
"simulation with a mutated dominant clone, nothing was lost at 40% abundance, but all mutations crossing the threshold "
"at 60% or 80% were removed, and read-weighted and clone-weighted estimates fell by as much as 62% and 34%, "
"respectively (STAR Methods). Excluded sites are therefore reported for inspection.",

"Read-weighted frequencies were within 2%":
"At 500 reads per library, read-weighted mutation frequencies were within 2% of the simulated values (Figure S8C). "
"Unique-variant density was less consistently recovered, because it counts each distinct substitution once: as depth "
"increases, additional reads carry mutations already counted, enlarging the denominator without adding variants. "
"Adding the aggregate estimate of uncalled single-read mutations restored 96–100% of the true density, and in a "
"separate test the estimate rose by 95% of the single-read mutations planted (Figure S8D). Comparisons of "
"unique-variant density therefore require matched sequencing depth.",

"Group-test calibration was assessed":
"Finally, we asked whether the group tests hold their stated error rate. In 300 simulations per condition, with a "
"between-animal coefficient of variation of 1, pooled-read tests returned p < 0.05 in 91% of 3-versus-3 and 94% of "
"5-versus-5 experiments without a true difference (Figure S8E), whereas negative-binomial rejection rates were 5.0% "
"and 5.7%. With three animals per group, however, that test detected two-fold and four-fold changes in only 12% and "
"30% of simulations, illustrating the limited power at this sample size.",

# ───────── Results: real libraries ─────────
"SHMtrace retained 1–4 calls across three controls":
"Across these assays, SHMtrace retained 1–4 calls in the three controls, compared with 5–20 for Geneious (Figures "
"4A–4E and 4G). As a result, geometric-mean unique-density separation, with zero-call controls assigned one-event "
"values, reached 141-, 85-, 66- and 94-fold for J_{H}4, N114 V_{H}, primary-cell Sμ and CH12 Sμ, respectively, "
"compared with 16–62-fold for Geneious and 2–40-fold for LoFreq configurations (Figure 4F). At 500 reads, J_{H}4 and "
"N114 V_{H} had reached 90% and 59% of the most each library yielded at any depth, without establishing saturation "
"(Figures S9 and S10). Background mismatches were position-dependent, yet separable from the abundant clone-shared "
"calls of Wt libraries (Figure S11).",

"The same N114 V_{H} reads were analyzed by IgBLAST":
"Repertoire tools meet the same background from the other side. When the same N114 V_{H} reads were analyzed by "
"IgBLAST[@ye2013] and SHazaM[@gupta2015], which count every mismatch to the germline (Figure S12), the AID libraries "
"were only 2-fold above the empty-vector libraries with the IMGT germline, because the parental line differs from it "
"(99.3% identity) and sequencing errors count as mutations. Using the parental V as germline raised this to 6-fold, "
"and masking bases below Q20 to 32-fold, compared with 188-fold for SHMtrace.",

"At *Pim1*, mutation calls in AID-competent":
"Not every locus separated so clearly. At *Pim1*, mutation calls in AID-competent libraries remained close to the "
"AID-null background: at 2,000 reads per library, SHMtrace yielded 1.8-fold separation, with four calls across the "
"controls, compared with 1.6-fold and eleven control calls for Geneious. Increasing depth to all available reads "
"(3,238–3,447 per library) did not clearly resolve the signal, with SHMtrace detecting 3–10 variants per AID-competent "
"library and 0–3 per AID-null library, and Geneious 5–7 and 0–4, respectively (Table S1).",

"The two CH12 AID off-target loci":
"Two CH12 AID off-target loci, off-target 1 and off-target 2, then tested a genetic perturbation expected to reduce "
"AID targeting[@seija2026]. Off-target 1 Wt libraries carried 69, 50, and 55 calls versus 0, 2, and 0 in "
"*Aicda*^{−/−} controls, with densities of 3.07 versus 0.04 per 10^{5} bases (Welch p = 9.7 × 10^{−4} on log density; "
"0.024 on frequency). Off-target 2 counts were 42, 28, and 23 versus 2, 0, and 0, with densities of 1.68 versus 0.04 "
"(p = 5.5 × 10^{−4} and 0.030), and G:C substitutions accounted for 97% of calls at both loci. MLLT1/MLLT3 double "
"knockout reduced counts to 2, 0, and 2 at off-target 1 (0.07 per 10^{5} bases) and 1, 0, and 0 at off-target 2 "
"(0.02), with no detected difference from *Aicda*^{−/−} controls (p = 0.51 and 0.43). We report no fold here, because "
"a zero-count convention would determine its magnitude.",

"Leave-one-control-out error margins":
"These results depend on the error margin chosen from the controls. Without it, 5.4% of J_{H}4 and 9.6% of N114 "
"V_{H} held-out control tests fell below p < 0.01; margins of five for the shared J_{H}4 window and ten for "
"full-length J_{H}4 and N114 V_{H} reduced these to 1.0%, 0.7% and 0.5%, respectively (Figure 4H). The margin "
"calibrates the error model, whereas the separate decoy threshold estimates the false-discovery proportion.",

"Read preparation removed 67%":
"Read preparation mattered as well. It removed 67% of delivered J_{H}4 reads overall (Figure S1E), and competitive "
"alignment alone was not enough: without an additional paralog-base filter, simulated paralog reads produced 35–54 "
"false variants per library (STAR Methods).",

# ───────── Results: Sanger reanalysis ─────────
"Reanalysis of published *Maea* J_{H}4":
"We next re-analysed published Sanger data. A J_{H}4 study had compared four Wt and four "
"*Maea*^{fl/fl}*Mb1*^{cre} mice with *Aicda*^{−/−} controls using SHMTool[@maccarthy2009,@barbulescu2024]. From its "
"chromatograms, SHMtrace recovered a 4.6-fold frequency reduction, against 4.4-fold in the publication (one-tailed "
"p = 2.3 × 10^{−4} and 0.0143; Figure S14). Spectra correlated at r = 0.99, and denominator-harmonized frequencies "
"averaged 1.0-fold the published values. Shared trace assignments and identical mutated-clone pairs were retained as "
"filed and flagged; the detailed discrepancies are in Table S4 and STAR Methods.",

"Reanalysis of 402 *Mkln1*-study":
"A second study from the same laboratory gave similar agreement. Reanalysis of its 402 *Mkln1*-study "
"chromatograms[@barbulescu2025] found an 8.0-fold frequency reduction versus the published 7.8-fold (p = 0.0038 and "
"0.014), and G:C transversions increased from 15% to 29%, against 16% to 27% in the published tables (p = 0.0078 and "
"0.0040). Spectrum correlations were r = 0.97 for Wt and 0.62 for knockout (Figure S15; Table S4).",

"An independent laboratory's SW_{HEL} study":
"Because both of those studies come from the laboratory that produced our libraries, we also re-analysed a Sanger "
"study from an unrelated laboratory at a different locus. Thientosapol et al. hypermutated a SW_{HEL} V_{H}10 knock-in "
"V region in germinal-centre B cells and deposited 1,795 sequences across five arms[@thientosapol2018] (Figure S16). "
"From these records alone, SHMtrace counted 2,804 of the 2,809 mutations the study tabulates (r > 0.9999 across 60 "
"spectrum cells), and mutations per sequence agreed to two decimal places in every arm. The four substitution classes "
"the study tested moved as reported, and total mutation decreased in both analyses, although significance differed "
"under the different tests used. One host that may pool several animals was retained, and excluding it did not change "
"the conclusions (STAR Methods). Reassuringly, the *Msh2*^{−/−} arm, which SHMtrace was not told was a repair mutant, "
"showed the expected collapse of A:T mutations. Because this study deposits sequences rather than chromatograms, it "
"tests counting, spectra and group statistics but not trace reading.",

# ───────── Results: long reads ─────────
"Junction sequences upstream of J_{H}4":
"Long reads also keep each molecule's junction, which groups reads into B-cell lineages (Figure 6C). In Wt2, one "
"largely unmutated lineage accounted for 81% of junction-assigned reads, so the read-weighted frequency of 2.8 per "
"10^{3} bp rose to 16.3 when each lineage's majority-rule mutation profile was counted once. Wt3 and Wt4 changed less "
"(22.5 to 24.9 and 4.9 to 5.3). The largest Wt3 lineage comprised 69 reads and four founder mutations, with five "
"inferred sublineages extending to 17 mutations (Figure 6D). Overall, 70% of mutated reads fitted one tree, and "
"7/1,128 mutation pairs failed the four-gamete test. Without UMIs, however, PCR recombination and residual errors "
"remain possible, so these are read-level reconstructions rather than proven cellular lineages.",

# ───────── Results: public datasets ─────────
"Three public Illumina studies tested":
"To test transfer beyond nanopore, we analysed three public Illumina studies of different biological contexts "
"(Figure 7). At J_{H}4 positions 53–122, Wt mice showed 56-fold greater mutation frequency than AID-null mice "
"(log-rate Welch p = 0.005; n = 7 and 4), whereas *Ldha* deletion showed no detectable difference from Wt "
"(Mann–Whitney p = 0.75), consistent with the original study[@sharma2023] (Figure 7A).",

"HMCES-deficient Ramos clones[@wu2022]":
"HMCES-deficient Ramos clones[@wu2022] carried 3.8-fold more deletions and 2.5-fold more insertions than Wt clones "
"(deletions, p = 1.0 × 10^{−4}; n = 6 per group; Figure 7D). Substitution calling in this deposit, however, depended "
"strongly on one noisy AID-null library. A rule developed after inspecting this dataset therefore excluded it from "
"threshold selection, while retaining it in model fitting and reporting the exclusion. With this adjustment, Wt and "
"HMCES-deficient clones yielded 52–99 substitution calls each, whereas treating the noisy library as a sample instead "
"yielded 277–413. Control assignment thus materially affects the inferred substitution burden here (STAR Methods), "
"and the same rule also acted in the mOrange analysis.",

"In the mOrange study[@perezduran2012]":
"A published mOrange reporter study compared active and catalytically inactive AID in cells carrying an empty "
"vector, Ugi or UNG[@perezduran2012]. Using the inactive-AID libraries to calibrate calling, SHMtrace retained 115 "
"substitutions, all at G:C and 75–100% at WRC/GYW hotspots (Figure 7E). Separation between active and inactive AID "
"was 20-fold with empty vector, 1,002-fold with Ugi and 103-fold with UNG, with Ugi zeros assigned one event. Among "
"retained substitutions, the transition share fell from 95% with Ugi to 4% with UNG, whereas uncorrected inactive-AID "
"backgrounds contained 43–66% transitions. The paired mean difference was 92 percentage points (95% CI 74–110; "
"paired t, p = 0.0021), although few calls per library limit precision.",

"PacBio Sμ–Sγ1 junctions[@vincendeau2022]":
"Finally, PacBio Sμ–Sγ1 reads[@vincendeau2022] carry the donor, the join and the acceptor on one molecule, so we "
"applied the junction module to a class-switch phenotype it was not tuned for (Figure S20). Mean microhomology was "
"1.01 nt in wild type and unchanged in *Shld1*^{−/−} (1.02 nt), but rose without XLF (1.69 nt) and in double mutants "
"(1.59 nt), which also showed more ≥4-nt microhomology and fewer blunt junctions; *Shld1*^{−/−} junctions instead "
"carried more insertions. This is the shift toward microhomology-mediated end joining expected when classical "
"non-homologous end joining is impaired. Because reads per unique junction vary between samples, counts reflect PCR "
"recovery rather than switching efficiency, and 19–47% of junctions place equally well at more than one acceptor "
"position (STAR Methods).",

# ───────── Discussion ─────────
"Measuring somatic hypermutation has long meant":
"Measuring somatic hypermutation has long meant choosing between reading molecules accurately and reading enough of "
"them. Our results suggest that this choice is no longer necessary, provided the background is measured rather than "
"assumed. SHMtrace takes whatever a platform delivers, reads or chromatograms, and replaces the two assumptions that "
"make low-frequency calling from single reads unreliable for SHM. First, rather than relying on quality scores alone, "
"it estimates the background in the AID-null libraries that SHM studies already include, at the level of substitution "
"type, strand and base quality, and checks that calibration by holding out each control in turn. Second, instead of "
"treating reads as replicates, it tests independent biological samples. Both replacements held: 81% of planted J_{H}4 "
"mutations and 90% on N114 V_{H} were recovered at false-discovery proportions near the default target, while a "
"quality-based caller either missed most true mutations or, made permissive, called as much background as the "
"uncorrected rule, and pooling reads produced false positives in more than 90% of simulated experiments without a "
"true difference. Together, these results support a practical route from deeply sequenced amplicons to interpretable "
"SHM measurements across biological replicates.",

"The gap this fills is specific":
"The gap this fills is specific, and it sits between two mature bodies of software. Repertoire toolkits identify and "
"summarise mutations well: IgBLAST assigns germlines[@ye2013], Change-O and SHazaM quantify load[@gupta2015], SHMTool "
"compares datasets[@maccarthy2009], and S5F models targeting[@yaari2013]. They begin, however, from mismatches already "
"taken to be real, which is why counting every mismatch on the same N114 reads separated AID from empty-vector "
"libraries by 2- to 32-fold, depending on germline and masking, where SHMtrace separated them by 188-fold. "
"General-purpose variant callers do the opposite: they decide what is real, but from error models fitted to data of a "
"different kind. LoFreq reads base qualities[@wilm2012], while ClairS and DeepSomatic were developed for tumour–normal "
"sequencing, with different training data and calling assumptions[@zheng2026,@park2026]. As configured here, neither "
"class measures what a particular SHM library does wrong. SHMtrace closes that circuit by making the AID-null control, "
"where available, the source of the error model, and by carrying the measurement through to the animal-level "
"comparison that the biology requires[@lazic2010]. Its distinctive contribution is therefore to connect "
"control-informed mutation detection with interpretable profiles and comparisons across biological replicates in one "
"reproducible workflow.",

"Seen this way, comparator performance":
"Seen this way, comparator performance is a matter of operating range rather than quality. LoFreq was designed for "
"deep short-read data with calibrated qualities; its base-alignment correction suppressed nearly all nanopore calls "
"here, and disabling it still limited recovery of low-support variants. That is not a defect on its intended data; "
"rather, as configured here, it did not capture this assay's library-specific background as well as matched-control "
"calibration did. The measurements suggest why. Above Q40, control error stopped decreasing and differed across the "
"evaluated datasets (1–9 × 10^{−5}; Figure S2), so the quality score runs out of information exactly where SHM calling "
"needs it. ClairS and DeepSomatic, similarly, recovered clean variants at higher allele fractions but few real SHM "
"calls, and chemistry and model mismatches further qualify the external-data comparisons. Clair3 and DeepVariant found "
"none of the planted mutations, as expected of callers built to genotype a sample against a reference rather than to "
"find a variant carried by a few reads. Geneious, in turn, provided greater single-read recall on J_{H}4 and *Pim1* "
"but more false calls on each amplicon. These are choices a study can make deliberately: sensitive mode raises the "
"empirical target to 0.20 and recovers 72% of single-read mutations on J_{H}4 and 87% on N114 V_{H} at observed FDP "
"0.18 and 0.16, while events unresolved individually can still be estimated in aggregate. Together, these outputs make "
"the sensitivity–background tradeoff explicit and help investigators choose a calling regime appropriate to their "
"experiment.",

"A shared framework accommodates four platforms":
"Because the background is measured rather than modelled, one framework serves all four platforms, each through the "
"analysis its data support. A Sanger clone is one molecule read once, so there is no read support to weigh against a "
"background; the workflow therefore counts quality-filtered substitutions directly and reports the control-derived "
"background beside them. Nanopore and Illumina libraries, by contrast, carry matched AID-null arms, so each candidate "
"is tested against the errors those controls make. Nanopore spans the whole amplicon, which made the 3′ decline "
"visible without tiling, identified clones from the V(D)J junction and ordered their mutations into trees; Illumina "
"gives depth without that linkage. PacBio combines length with per-base accuracy, and reading the donor, the join and "
"the acceptor on one molecule recovered a known end-joining phenotype from a deposit with no AID-null arm, using "
"junction shape in place of calibrated calls. Figure S21 sets out what each measured. Agreement between platforms is "
"itself evidence: nanopore and Sanger reads from the same six mice gave substitution spectra correlated at r = 0.92 "
"and separated the genotypes in the same direction, supporting the use of cloning-free long reads to recover "
"established SHM patterns while adding whole-amplicon linkage.",

"Reanalysis of public data tested SHMtrace":
"Applying the same engine to other laboratories' data turned out to test it harder than any internal benchmark. "
"Published sequence collections support reproducible reanalysis when clones are linked to animals and genotypes: "
"SHMtrace recovered the principal J_{H}4 frequency and G:C comparisons, matched SW_{HEL} counts across five arms and "
"twelve substitution types (2,804 of 2,809 tabulated mutations; r > 0.9999), and recovered 150 of 150 mutations an "
"antibody-evolution study curated by hand, calling none its authors did not. Agreement was not complete, and the "
"exceptions are as informative: Wt frequencies differed after denominator harmonisation, the published *Mkln1* A:T "
"decrease was not significant in the reanalysis, and the SW_{HEL} total-mutation inference differed under a different "
"test (Tables S4 and S5; Figures S14–S16). The Sanger re-analysis also surfaced bookkeeping problems, which SHMtrace "
"now checks and reports rather than silently resolving. By retaining assignments, counting conventions and quality "
"flags, SHMtrace makes legacy sequence collections easier to revisit, compare and reuse.",

"What this enables in practice":
"What this enables in practice is a short experiment-to-answer path that keeps its own evidence: read accounting, "
"fixed-seed sampling and one engine behind the desktop and command-line interfaces, with every exclusion, threshold and "
"flagged site written out for inspection. The power analysis sets the honest expectation. How many animals a study "
"needs depends on the difference it seeks: under variation estimated from three Wt J_{H}4 mice, the AID-competent "
"against AID-null contrasts made here separate at three animals per group, whereas a two-fold difference between two "
"AID-proficient groups would require about 24, and more if three mice underestimate the variation (Figure S8F). "
"Sequencing depth cannot substitute for independent biological samples, and no analysis can recover a comparison the "
"design did not support. Within that boundary, however, the framework extends readily. It needs only a target, "
"matched AID-null libraries and a sample sheet, which is why it could be applied to off-target loci, switch regions, "
"engineered diversification and base-editor substrates across twelve loci from eleven laboratories. Further "
"validation on independently prepared runs and UMI-tagged SHM libraries will establish how far it extends. More "
"broadly, by bringing control-calibrated mutation detection, whole-amplicon profiles and comparisons across biological "
"replicates into one reproducible workflow, SHMtrace makes long-read SHM data more informative and existing sequence "
"collections easier to reanalyse transparently, providing a practical foundation for studying how AID targeting and "
"DNA repair shape antibody diversification.",
}

CAPTIONS = {
"Figure 1.":
"(A) Quality-based and control-calibrated decisions, drawn schematically. SHMtrace estimates background by "
"substitution, strand, quality and position from matched AID-null controls. (B) Inputs and read preparation; Sanger "
"traces are trimmed at Q≥20 instead of the expected-error gate and are not downsampled. (C) Read libraries are tested "
"against a control-fitted error model with an empirical, decoy-based threshold, and single-read support is allowed. "
"Sanger Q≥20 substitutions are counted, with expected background reported beside them rather than subtracted. Both "
"routes share the window, quality checks and the heuristic ≥50% allele exclusion. (D) Outputs, with comparisons "
"across biological replicates. PacBio reads are analysed for junction structure only. See also Figure S1.",
"Figure 2.":
"(A, B) Support for two C-to-T candidates, at a hotspot (A) and a recurrent-background position (B), in Sμ WT3 and its "
"three AID-null controls (PRJNA1481264; 2,000 reads per library). Counts include all qualities, with Q≥20 subsets in "
"parentheses; p, SHMtrace's test. (C) Base qualities of the supporting reads; bars, medians; dotted line, Q20. (D) "
"Candidates kept and decoy-estimated false calls across p-value thresholds; the threshold is the largest at which the "
"estimate stays within 5% of the candidates kept, an empirical estimate rather than a guaranteed FDR. (E) Allele "
"fractions of 1,197 retained calls in 12 AID-active libraries; dotted line, ≥50% exclusion. See Figures 3C and S2.",
"Figure 3.":
"(A, B) Recall of planted substitutions by supporting reads: 60 sites × 3 draws × 3 held-out controls per locus; "
"vertical lines, range over controls; FDP, pooled false over all calls. (C) Permissive somatic callers on J_{H}4, "
"given peer controls as matched normals: solid, the spike-ins of A; dashed, clean plants at 5–20% of reads. Grey, "
"cumulative share of the real SHMtrace calls in Figure 2E. (D) Recall against false calls per library for every caller "
"at four loci; the x axis is linear from 0 to 0.1 and logarithmic above. BAQ, base-alignment quality. See Figures "
"S4–S8 and S24; Tables S1 and S6.",
"Figure 4.":
"(A–E) Unique variants per 10^{5} analyzed bases for each caller: blue, AID-competent libraries; red, AID-null; shaded "
"strip below the axis break, libraries with no call. Axis ranges differ by panel. J_{H}4 and N114: three mice or "
"subclones per group at 500 reads; Sμ and *Pim1*: 3–4 case and three null libraries at 2,000 reads. Welch tests on "
"log density, Holm-adjusted across callers, with zero-call libraries at one event; ns, not significant, which with "
"three libraries per group does not show that a caller fails; \\*p<0.05; \\*\\*p<0.01; \\*\\*\\*p<0.001. (F) "
"Geometric-mean case/control separation, same zero convention; CH12 OT1/OT2, AID off-target loci in CH12 cells, with "
"Wt and double-knockout groups shown separately. (G) Calls summed over three controls. (H) Held-out control tests "
"with p<0.01 by error margin; dotted line, 1% expected if calibrated; rings, selected margin. Panels use the "
"caller-comparison calibration (STAR Methods). See Figures S9–S11; Tables S1 and S3.",
"Figure 5.":
"Matched 529-bp data from three Wt and three AID-null mice[@barbulescu2025]. (A) Unique-variant densities within "
"each platform: distinct Q≥20 changes among 26–40 clones per Wt mouse (Sanger) and control-calibrated calls among 500 "
"reads (nanopore). Bars, geometric means; p, log-Welch; nanopore zeros at one event. (B) Pooled shares of the twelve "
"substitution types (471 Sanger and 638 nanopore calls); dashed line, identity. The correlation describes similarity "
"of pattern, not interchangeable values. (C) Four classes by mouse; bars, mean±SD; shapes, mice. (D) Positional "
"shares in 25-bp bins; the final shaded bin is 4 bp. (E) Hotspot use; bars, mouse means; dashed lines, share of "
"eligible positions in each motif. See Figures S13–S16; Tables S4 and S5.",
"Figure 6.":
"Libraries as in Figure 5. (A) J_{H}4 variants per 10^{5} bases in 50-bp bins across positions 90–1,230; thin lines, "
"Wt mice; thick lines, genotype means. Shading, the Sanger window. (B) Calls per position by S5F mutability quintile, "
"relative to each series' mean, for Wt calls and germline-excluded null mismatches; the five points are bins of one "
"window, not replicates. (C) Junction-defined clones among 500 reads per mouse (dark, largest; hatched, junction "
"unreadable), with read- and clone-weighted frequency alongside. (D) The largest Wt3 clone: square, inferred "
"ancestor; ticks, mutations; numbers, reads carrying every mutation up to that node. Of 69 reads, 48 fit one path and "
"21 carry mutations of two branches. The tree reconstructs reads, not proven cellular ancestry. See Figure S17.",
"Figure 7.":
"(A) J_{H}4 mutation frequency in Illumina data (PRJNA973112); points, mice; strip below the break, no call; bars, "
"geometric means with zeros at one event. Tests: log-Welch versus AID-null; Mann–Whitney for *Ldha* versus Wt. (B) "
"*Aicda*^{−/−}*Ung*^{−/−}*Msh2*^{−/−} CH12F3 cells expressing AID variants or empty vector (PRJNA918596; n = 4 per "
"arm). (C) Substitution classes of each arm's pooled calls; numbers in bars are calls. (D) Ramos indels per clone at "
"20,000 reads (PRJNA820976); bars, means. (E) mOrange repair arms (PRJNA170048): circles, SHMtrace calls; squares and "
"crosses, raw Q>20 mismatches for active and inactive AID; lines join one experiment's three arms. See Table S1.",
}
