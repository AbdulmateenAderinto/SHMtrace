# Prompt for an external review of the SHMtrace manuscript

Paste the block below into ChatGPT and attach, or paste after it,
`crm/build/SHMtrace_manuscript.txt` (the article, ~20,800 words including STAR Methods) and, if you want the
figures reviewed too, `crm/build/SHMtrace_Document_S1.txt` plus the figure images from
`crm/figures/`. The reviewer cannot run the code or see the data, and the prompt says so; it is written to
get judgement on what is written, not a re-analysis.

---

## The prompt

You are reviewing a methods manuscript submitted to Cell Reports Methods. Act as a critical, numerate referee
with expertise in somatic hypermutation, variant calling from sequencing reads, and applied statistics. Your
job is to find what is wrong, overstated or unclear, not to summarise or praise.

**What the paper claims.** SHMtrace calls somatic hypermutation from sequencing reads by measuring the error
background in AID-null control libraries sequenced in the same run, instead of trusting per-base quality
scores. It reports an empirically estimated false-discovery rate targeting 5%, set by using held-out controls
as decoys. The paper validates this on the authors' own nanopore and Sanger data and on several public
datasets across platforms.

**What you can and cannot do.** You have only the text, and the figures if supplied. You cannot run the
software, inspect the code, or see the underlying data. Do not claim to have verified anything you have not.
Where a claim can only be checked against code or data, say so and say exactly what you would need. Do not
invent references, accession numbers, statistics or quotations; if you think a citation is missing, describe
what is needed rather than supplying one.

**Review these in order.**

1. **Internal consistency of the numbers.** The manuscript is dense with figures, percentages, p-values and
   sample sizes. Check that numbers in the abstract, Results, figure legends and STAR Methods agree with each
   other; that every n is stated and consistent; that percentages are of the denominator they claim; and that
   no number appears in two places with two values. List every discrepancy with its location.

2. **Whether the central claim is earned.** The FDR is *estimated from the controls themselves*, and the same
   controls also choose the calibration margin. Is the paper honest and precise about the difference between
   an empirical target and a guarantee? Is the circularity acknowledged where it matters? Would a reader come
   away believing something stronger than the evidence supports?

3. **Statistical soundness.** Much of the validation uses three animals or libraries per group. Assess whether
   the tests are appropriate at that size, whether effect sizes and uncertainty are reported alongside
   p-values, whether any comparison is presented as significant when it is underpowered, and whether anything
   is pooled that should not be. Note any place where replicates, reads or clones are conflated.

4. **The comparisons with other tools.** The paper compares against the Geneious rule, LoFreq, ClairS and
   DeepSomatic. Is each tool used in a setting it was designed for, or at least is the mismatch stated? Is the
   comparison fair, and does the paper say where another tool is better?

5. **Overreach in the public re-analyses.** Several datasets have known weaknesses the paper should be
   handling explicitly: deposits with a single control library; deposits where read layout is confounded with
   the biological groups; selected-clone collections that carry no population frequencies; datasets with no
   AID-null arm at all. For each public dataset, judge whether the text states plainly what that dataset can
   and cannot establish, or whether it borrows credibility from the stronger datasets.

6. **The germline heuristic.** Changes carried by at least half of a library's reads are set aside as
   germline. Is the consequence — that a mutation shared by a dominant clone is removed — stated clearly
   enough, and is the opt-in exception for experiments with a known ancestral reference adequately justified
   and bounded?

7. **Reproducibility and completeness.** Could a competent group repeat this from the STAR Methods alone?
   Identify missing parameters, undefined terms, unstated software versions, and any analysis mentioned in
   Results but absent from Methods.

8. **Writing.** Flag sentences that are ambiguous, that overclaim, that bury a caveat, or that a non-specialist
   editor would stumble over. Flag any hedging so heavy it obscures what was actually found.

**Output.** Give me, in this order:

- A referee's summary verdict in under 150 words: what this paper establishes, and the single biggest
  weakness.
- **Major issues**, numbered, each with: the location, what is wrong, why it matters, and what would fix it.
- **Minor issues and corrections**, as a list with locations.
- **Numerical discrepancies**, as a table: location, the two conflicting values, which looks right.
- **Questions you would ask the authors**, where you cannot decide from the text alone.

Rank issues by how much they threaten the conclusions, not by how easy they are to fix. If you believe a
claim is wrong, say so directly. If you are uncertain, say that instead of hedging into vagueness.
