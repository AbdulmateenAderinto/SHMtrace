# Review brief: SHMtrace manuscript

Please review this manuscript as a careful, critical peer reviewer for *Cell Reports Methods* would, and also as
a co-author checking it before submission. Be specific and blunt. Praise is not needed.

## What it is

**Title:** Control-calibrated profiling of somatic hypermutation from cloning-free long reads, Sanger traces and
short reads. Aderinto and Martin, University of Toronto. Unpublished; not yet submitted.

**Target journal:** *Cell Reports Methods* (Cell Press). A shorter companion version for *STAR Protocols* exists and
is not included.

**In one paragraph:** SHMtrace measures somatic hypermutation (SHM) from amplicon sequencing without cloning. Each
candidate mutation is tested against the errors measured in matched AID-null (*Aicda*⁻/⁻) libraries sequenced in the
same run, and groups are compared with the animal, not the read, as the unit. The paper validates it with:
- planted mutations and simulated B cell lineages;
- ten loci from seven laboratories on three platforms (nanopore, Sanger, Illumina);
- nanopore against Sanger sequencing of the same mice;
- automatic re-analysis of three published Sanger studies from their deposited clones: two from the authors'
  laboratory and one external (Thientosapol et al. 2018, PNAS).

## Files

| File | Contents |
|---|---|
| `01_Manuscript.pdf` | Main text, figures 1-7 with legends, STAR Methods, references |
| `01_Manuscript.txt` | The same as plain text (layout preserved), for easier quoting |
| `02_Supplement_Document_S1.pdf` | Document S1: supplementary figures S1-S18 and tables S1, S3, S4, S5 |
| `02_Supplement_Document_S1.txt` | The same as plain text |
| `03_Highlights_eTOC_Motivation.txt` | Highlights, eTOC blurb, Motivation statement |
| `04_Table_S2_benchmark.xlsx` | Table S2, the benchmark spreadsheet |
| `figures/` | Full-resolution images of the main figures, the graphical abstract, most supplementary figures and Table S4, for close inspection |

Figures S4, S7 and S11 appear only inside Document S1.

## What to check

1. **Claims against evidence.** For each claim in the summary, highlights and results, is it supported by a figure or
   table, and is it stated no more strongly than the data allow? Flag overclaiming, including words like
   "unaided", "reproduces", "concordant" and "validated".
2. **Statistics.**
   - Are the tests appropriate? Animals are meant to be the unit throughout; flag any place where reads or clones are
     treated as replicates.
   - Are the multiple-testing and false-discovery claims justified?
   - Are P values, fold changes and correlations reported with enough context (n, test, sidedness)?
3. **Internal consistency.** Do numbers agree between the summary, results text, figure panels, legends and Table S4?
   List every mismatch you find with both locations.
4. **Methods.** Could someone reproduce each analysis from the STAR Methods? What is missing (parameters, versions,
   thresholds, data and code availability)?
5. **Controls and confounders.** Batch effects, sequencing-run effects, PCR duplicates and lineage expansion,
   germline differences, paralogous amplicons, and the choice of AID-null libraries as the error model.
6. **The re-analyses of published studies** (Figures S8, S9 and S10; Table S4). Are the comparisons fair? Are the
   remaining differences explained adequately? Is the language about published analyses neutral and factual?
7. **Figures.** Readability at print size, axis labels and units, colour use (colour-blind safety), legends that
   fully define symbols and error bars, and panel ordering that matches the text.
8. **Journal requirements.** The build checks summary ≤ 150 words, highlights ≤ 85 characters each and eTOC
   ≤ 50 words. Also check the STAR Methods structure, the key resources table if present, the "Limitations of the
   study" section and the lead contact statement.
9. **Writing.** Clarity for a non-specialist immunologist, jargon, sentences that are too long, and anything
   ambiguous.

## Things the authors already know

These are disclosed in the manuscript. Comment on whether the disclosure is adequate, but do not treat them as
discoveries:
- In the *Mkln1* re-analysis (Figure S9), only one of the three published knockout per-mouse tables could be matched to
  the clones of a mouse folder in the original project, which limits knockout agreement (r = 0.62 against 0.97 for
  wild type).
- The *Mkln1* and *Maea* studies counted paralogous amplicons as sequences, and the published frequencies divide by
  a fixed amplicon length. Table S4 reports SHMtrace's numbers both ways.
- One nanopore assay out of five did not reach the AID-null separation of the others (see the summary).

## How to report

1. A three-to-five-sentence overall assessment, including whether you would recommend major revision, minor
   revision or acceptance, and why.
2. **Major issues**, numbered. For each: where it is (page, figure or panel, or a short quote), what is wrong, why it
   matters, and a concrete fix.
3. **Minor issues**, numbered, in the same format.
4. **Inconsistencies**: a table of location A, location B, and what disagrees.
5. Any reference you think is missing. Give it only if you are certain it exists; otherwise describe the kind of
   work that should be cited. Do not invent citations.

When unsure, say so rather than guessing. Quote the manuscript exactly when pointing at text.
