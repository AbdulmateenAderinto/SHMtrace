# Figure grouping audit, 21 October 2026

Every main and supplementary figure, read panel by panel, asking one question: **does each figure hold
one argument, and do panels that belong together sit together?** Checked against the v1.2 build
(commit `ac002b1`), 7 main figures and 22 supplementary figures over 24 pages.

**Verdict: the grouping is sound in 24 of 29 figures.** (Items 1 and 2 were fixed on 21 October, so 26 now.) Five have a panel that sits apart from the
material it belongs with, and two of those are duplications across the main/supplementary boundary.
None is wrong, and the paper reads correctly as it stands; what follows is where a reader has to carry
something across the figure to use it.

Fixed while auditing (see also item 2): Figure S1D said "447 traces (.ab1), 8 orders", which was comment 10 of the
20 October round. The text was corrected then, the figure was not. It now reads "184 traces (.ab1),
these six mice". The same figure's bare `N114` labels became `N114 V_H`.

---

## What is grouped well

| Figure | Why it holds together |
|---|---|
| **1** | Problem (A) → inputs and preparation (B) → the two calling routes (C) → outputs (D). A reader can follow it top to bottom without going back. |
| **2** | A–D are all the *same two candidates* in one Sμ library: their support, their qualities, the threshold they are judged by. One worked example, carried through. |
| **4 A–G** | A–E is one panel per dataset; F and G are the two caller × dataset matrices that summarise them. The summary sits immediately after the thing it summarises. |
| **5** | Six mice, two platforms, five views of the same comparison. |
| **6** | Two clean pairs: A, B on *where* mutations fall; C, D on *clonal structure*. |
| **S9** | A, B (J_H4) then C, D (N114 V_H) — the same two views per amplicon, in the same order. |
| **S10** | A, B along the amplicon; C, D within the read. Two scales, two pairs. |
| **S13 / S14 / S15** | Three reanalyses built on one spine: density → spectrum → a table of published claims against SHMtrace's. Parallel structure makes them readable as a set. |
| **S23** | Four views of one stress test, A→D from false discovery to exclusion. |

---

## Five places where related material sits apart

### 1. ~~Figure 3H is separated from Figure 3A–C~~ — FIXED 21 October

Figure 3 is a 2 × 4 grid. A, B are spike-in recall by supporting reads; C is recall against false
calls; **H is the same benchmark run against permissive somatic callers** (ClairS, DeepSomatic). H
asks A's question with a different caller family, but sits in the far corner with D, E, F, G between.

Between them lie three *different* questions: D (lineage simulations), E, F (quantification against
truth), G (which statistical test to use for group comparisons).

**Done.** The somatic-caller panel moved up beside the spike-ins and the panels were relettered so
they read in position order: row 1 is now A, B, C (somatic callers), D (recall against false calls) —
every caller comparison, closing with the summary; row 2 is E (lineages), F (clone estimates),
G (quantification), H (group tests). The legend opens by saying which row answers which question.

The old letters C–H all moved (C→D, D→E, E→F, F→G, G→H, H→C), so every citation was remapped in one
pass across the text, the methods, the legends and two figure scripts. `star_text.py` was left alone:
its "Figure 3" is the STAR Protocols article's own.

### 2. ~~Figure 4I duplicates Figure S1E~~ — FIXED 21 October

A–H compare five callers over six datasets. **I is the fate of every delivered read** for this study's
twelve libraries — no callers in it, a different scope, and a bookkeeping question rather than a
calling one.

It also overlaps **Figure S1E**, "Read accounting in SHMtrace", which shows the same twelve libraries'
delivered reads with the 500-read line. S1E gives three categories and totals; 4I gives seven
categories as percentages. Same analysis, two figures.

**Done.** Panel E was removed from Figure S1, whose seven-category counterpart in Figure 4I is the
more informative of the two. Figure S1 is now "Detailed workflow", A–D, and 153 mm tall instead of 188;
its legend points at Figure 4I for the read accounting.

S1E carried one thing 4I did not: absolute counts with the 500-read line, which showed at a glance that
every library had reads to spare. Figure 4I is percentages only, so that point now travels as a number
in its legend — each library yielded 962–2,905 usable reads against the 500 analysed — computed in
`crm_numbers.py` from the same accounting, not written in by hand.

### 3. Figure 7E is a method validation among biological datasets

Figure 7's panels are public AID datasets: A (J_H4, Illumina), B and C (CH12F3 AID mutants — the same
dataset, correctly adjacent), D (Ramos indels), F (ectopic AID across repair backgrounds).

**E is rRNA operon reads**, measuring error remaining after UMI consensus. It is a sequencing-accuracy
measurement, not an AID experiment, and it is the one panel in the figure with no AID in it.

It is also the second half of a split topic: **Figure S7C is UMI consensus on simulated libraries**,
and 7E is UMI consensus on real data. The evidence for one claim sits either side of the
main/supplementary boundary.

**Suggested:** move 7E next to S7C so the UMI consensus case is made in one place. Figure 7 then holds
five panels, all of them public AID datasets.

### 4. Figure S7 is four unrelated checks in one figure

Its own title says so: "Single-read estimate, power, UMI consensus and planted indels". A (single-read
estimation), B (statistical power), C (UMI consensus) and D (indel recall) share no data, no axis and
no question — only the fact that each supports Figure 3 somewhere.

This is the weakest grouping in the paper. It is a defensible way to hold down the supplementary count,
and splitting it would add figure numbers, so it may be worth keeping as is. If 7E moves here (item 3),
C and the new panel at least form a pair inside it.

### 5. Figure 5B and 5E are the same comparison, two panels apart

B is the pooled twelve-type substitution spectrum, Sanger against nanopore. E is that same comparison
resolved per mouse. Between them sit C (positional) and D (hotspot motifs).

Minor, and the present order has its own logic — pooled result first, then the two other views, then
back to the per-mouse detail. But a reader checking whether the pooled agreement holds per mouse has to
skip two panels to do it. B and E adjacent would read more directly.

---

## Smaller notes

- **Figure S16's two pages are two different analyses** — full-length spectra (page 1) and mutation
  trees (page 2) — sharing a number because both support Figure 6. Defensible, since page 1 supports
  6A/6B and page 2 supports 6D, but they are not one figure in the way S8's two pages are.
- **S13 orders its panels differently from S14.** S13 puts the G:C panel after the clone-size
  distribution (A, B, C, D = density, agreement, clones, G:C); S14 puts it directly after density
  (A, B = density, G:C). Aligning S13 to S14's order would make the three reanalysis figures scan
  identically.
- **Figure 2E widens from the worked example to every dataset.** It is the one panel in Figure 2 not
  about the two Sμ candidates. It sits last, which is the right place for it, so this is noted rather
  than flagged.

---

## What is left

Items 1 and 2 are fixed. Items 3, 4 and 5 remain:

- **3** — Figure 7E (UMI consensus on rRNA reads) is a method validation among public AID datasets, and
  its simulated counterpart is Figure S7C. Moving 7E beside S7C would put that evidence in one place
  and leave Figure 7 as five panels of AID data.
- **4** — Figure S7 is four unrelated checks under one number. Defensible as a way of holding down the
  supplementary count; item 3 would at least give it a matching pair.
- **5** — Figure 5B and 5E are the pooled and per-mouse views of one comparison, two panels apart.

None is a duplication, so the paper is correct either way.
