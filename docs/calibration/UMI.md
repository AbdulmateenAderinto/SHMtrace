# UMI consensus: one read per molecule

*12 October 2026. `shmtrace/umi.py`; simulation `docs/calibration/umi_sim.py`, output in `umi_sim/`.*


**In the app:** Set up → Advanced → *UMI consensus*. Choose where the UMI sits and how long it is; the anchors default to the reference ends. The results report molecules per sample.

## What it does

A unique molecular identifier (UMI) is a random tag put on each template molecule before PCR, so all reads sharing a UMI are copies of one molecule. With `umi` on, shmtrace:

1. **Reads each read's UMI.**
   - `flank`: the UMI is next to an anchor at a read end. By default the anchor is the first or last 20 reference bases (the primer), found with up to 6 edits. Use `umi_ends` = left, right or both for dual UMIs.
   - `header`: the UMI is in the read name, as umi_tools and similar tools write it.
2. **Groups reads by UMI**, allowing `umi_max_edits` (3) edits per UMI end, among reads aligned to the same reference. Groups are built greedily by abundance against centroids. A UMI that lost or gained a base, or was cut short at the read end, still finds its centroid. With UMIs at both ends, a UMI seen once or twice may also join a group at least three times its size within two more edits, as UMI-tools' "directional" method allows, because it is likelier an error copy than a molecule of its own. With one UMI this step is off (see *Grouping* below).
3. **Builds a consensus** for each group of `umi_min_reads` (3) or more reads. At each reference position the base or deletion with the highest quality-weighted likelihood wins, and its quality is the posterior against the runner-up, capped at Q60. A disputed base therefore gets a low quality and is not trusted as a mutation.

Everything downstream (calling, controls, densities, clones, QC) then counts molecules. Controls must be UMI libraries prepared the same way.

### Settings

| Setting | Default | |
|---|---|---|
| `umi` | off | `flank` or `header` |
| `umi_ends` | both | `left`, `right` or `both` |
| `umi_length` | 18 | |
| `umi_offset` | 0 | bases between the UMI and its anchor |
| `umi_anchor_length` | 20 | or give `umi_left_anchor` / `umi_right_anchor` |
| `umi_anchor_edits` | 6 | |
| `umi_max_edits` | 3 | per UMI end |
| `umi_directional` | auto | extra edits for the directional step: 2 with UMIs at both ends, off with one |
| `umi_min_reads` | 3 | |
| `umi_header_regex` | `(?:_\|UMI[=:])([ACGTN]{6,})` | |

## Validation

There is no UMI data from the lab yet, so it was validated by simulation on real reads. For each amplicon, each mutation-free control was held out in turn:
- Its reads became the copies of simulated molecules: about 6 copies per molecule, and each copy a *different* real read with the molecule's mutations written in. Every copy therefore carries its own real nanopore errors.
- An 18-nt UMI was written just outside the reference end(s), with errors at Q15 plus 2% indels. Copies that don't reach a reference end get no UMI there.
- The other two controls became UMI libraries with nothing planted, used as the null.

J꜀4 reads begin at their V(D)J junction rather than at the reference start, so J꜀4 used a 3′ UMI only; N114 used both ends. The same simulated reads were analysed with UMI consensus and read by read, and both were scored against the molecules' true mutations. Means over three held-out controls:

| Amplicon · scenario | Analysis | Units | False calls | Recall | Recall of single-molecule mutations | Frequency est/true | Unique density est/true |
|---|---|---:|---:|---:|---:|---:|---:|
| J꜀4 · polyclonal | every read | 500 | 7.0 | 0.94 | 0.80 | 0.99 | 0.97 |
| J꜀4 · polyclonal | **UMI consensus** | 155 | **0.7** | **0.99** | **0.99** | 0.98 | 0.99 |
| J꜀4 · low | every read | 500 | 2.3 | 0.91 | 0.77 | 1.02 | 0.94 |
| J꜀4 · low | **UMI consensus** | 156 | **0.0** | **1.00** | **1.00** | 1.00 | 1.00 |
| N114 · subclonal | every read | 500 | 6.3 | 0.97 | 0.94 | 1.03 | 1.08 |
| N114 · subclonal | **UMI consensus** | 127 | **0.0** | **1.00** | **1.00** | 1.00 | 1.00 |
| N114 · low | every read | 500 | 5.7 | 0.97 | 0.94 | 1.03 | 1.07 |
| N114 · low | **UMI consensus** | 129 | **0.0** | **1.00** | **1.00** | 1.00 | 1.00 |
| Null (both amplicons) | UMI consensus | 126–153 | 0.0 | | | | |

Every consensus came from a single molecule's reads (group purity 1.00).

In this simulation, read-by-read analysis already benefits from about six copies of each molecule, so its numbers are better than for an ordinary library. The gain from UMIs is largest where it matters: mutations on a single molecule, and the background.

### What getting UMIs to work showed

- **Reads lost to filtering.** About two-thirds of raw reads fail the pipeline's filters. A design that counts on 3 or more copies per molecule needs about 3× more raw reads than molecules × 3.
- **Errors in the UMI itself.** The UMI is read with nanopore errors, so copies of one molecule differ from each other by 0–6 edits. Plain centroid grouping split molecules; the directional rule fixed this. A UMI at the very end of a read had often lost a base and was rejected; now a UMI short by up to `umi_max_edits` bases is kept.

## Grouping: indels and the directional step

*24 October 2026.* Both problems were found while porting this module to VDJtrace.

**Candidate lookup.** `umi.group` finds candidate centroids by exact pieces of the UMI. With at most *e* edits, one of *e* + 1 pieces is untouched. Pieces used to be cut from each UMI at its own length and compared position by position. A UMI one base short (cut at the read end, which extraction allows) or with an indel then shared no piece with its centroid and became a molecule of its own, even one edit away (`ACGTACGTACG` against `ACGTACGTACGT`). Only the directional step's full scan sometimes hid this. Now every centroid is cut into pieces of one length, set by the commonest UMI length, and every substring of the query is looked up. Centroids too short to cut are always checked. `tests/test_umi.py` covers truncated, deleted and inserted UMIs that the old lookup split off.

**Per reference.** Reads are grouped only with reads aligned to the same reference: one molecule has one. VDJtrace does the same within one locus, V gene and J gene. On real Takara 12-nt UMIs, this cut the share of groups mixing rearrangements from 40% to 3.5% at 2 edits.

**The directional step** lets a UMI seen once or twice join a group at least 3× its size up to `umi_directional` edits beyond `umi_max_edits`. VDJtrace found that at +2 edits it mixed different rearrangements in 95% of groups on real 12-nt Takara UMIs, which are T-rich. So it was measured here at +2, +1 and off.

*Karst et al. R10 reads* (2 × 18-nt UMIs, 6 edits; `umi_karst.py`, `umi_karst/results.json` → `directional`):

| Directional | UMI groups | Molecules (≥ 3 reads) | Reads in them | Molecules of 3–4 reads | Consensus error | Substitution error | Groups mixing species |
|---|---:|---:|---:|---:|---:|---:|---:|
| **+2** (default, dual UMIs) | 213 | 105 | 3,962 | 6 | 0.21% | 0.027% | 0% |
| +1 | 272 | 107 | 3,901 | 8 | 0.26% | 0.048% | 0% |
| off | 360 | 110 | 3,812 | 11 | 0.38% | 0.088% | 0% |

Without the step, about 150 reads fall out of molecules and five more small molecules appear. These are most likely pieces split off real molecules. With 3–4 reads each they carry high error, so substitution error triples. The step merged no reads of different species at any setting. For comparison, the old lookup at +2 gave 216 UMI groups and 103 molecules.

*Simulated UMIs at scale* (`umi_sim.py --grouping`, `umi_sim/grouping.md`). N molecules, each with 1 + Poisson(5) copies. The UMI is read at Q15 with 2% indels, and 5% of copies are cut 1–2 bases short. UMIs are uniform, or T-rich (35% T, as in the Takara UMIs).

| Molecules | UMI | Composition | Directional | Groups mixing molecules | Purity | Molecules split | Reads in groups ≥ 3 |
|---:|---|---|---|---:|---:|---:|---:|
| 1,000 | 18 nt | uniform | +2 | 6.2% | 0.989 | 0.5% | 95.5% |
| 1,000 | 18 nt | uniform | +1 | 0.5% | 0.999 | 0.6% | 94.9% |
| 1,000 | 18 nt | uniform | **off** | 0.2% | 1.000 | 0.7% | 93.1% |
| 1,000 | 18 nt | T-rich | +2 | 8.6% | 0.986 | 0.4% | 95.4% |
| 1,000 | 18 nt | T-rich | **off** | 0.2% | 1.000 | 0.7% | 92.9% |
| 3,000 | 18 nt | uniform | +2 | 15.9% | 0.972 | 0.3% | 95.6% |
| 3,000 | 18 nt | uniform | +1 | 1.8% | 0.997 | 0.4% | 94.9% |
| 3,000 | 18 nt | uniform | **off** | 0.5% | 0.999 | 0.5% | 93.1% |
| 3,000 | 18 nt | T-rich | +2 | 22.4% | 0.958 | 0.3% | 95.9% |
| 3,000 | 18 nt | T-rich | +1 | 3.0% | 0.995 | 0.4% | 94.8% |
| 3,000 | 18 nt | T-rich | **off** | 1.0% | 0.998 | 0.5% | 93.0% |
| 1,000 | 2 × 18 nt | T-rich | **+2** | 0.0% | 1.000 | 1.6% | 94.8% |
| 1,000 | 2 × 18 nt | T-rich | off | 0.0% | 1.000 | 1.9% | 92.5% |
| 3,000 | 2 × 18 nt | uniform | **+2** | 0.0% | 1.000 | 1.4% | 95.0% |
| 3,000 | 2 × 18 nt | uniform | off | 0.0% | 1.000 | 1.8% | 92.7% |
| 3,000 | 2 × 18 nt | T-rich | **+2** | 0.0% | 1.000 | 1.6% | 94.8% |
| 3,000 | 2 × 18 nt | T-rich | off | 0.0% | 1.000 | 1.9% | 92.5% |

Every row is a mean of three replicates; the full table, 150 molecules included, is in `umi_sim/grouping.md`. At 150 molecules a single UMI at +2 mixed molecules in under 1% of groups, and the other settings in none. In the pipeline simulation on the J꜀4 example controls (`examples_compare/KO1–3`, about 68 molecules), all three settings gave the same molecules, calls and purity (1.00).

**Decision.**
- With two 18-nt UMIs (36 nt, up to 8 edits at +2), the step recovers split molecules on real reads and mixed none in simulation, so it stays at +2.
- With one 18-nt UMI (up to 5 edits at +2), it mixed molecules in 6–9% of groups at 1,000 molecules and 16–22% at 3,000. That is the failure VDJtrace saw on 12-nt UMIs, while the step gained only about 2% of reads, so it is off.

`umi_directional` overrides either choice.

## On real nanopore UMI data

`docs/calibration/umi_karst.py`, output in `umi_karst/results.json`.

**Data:** the Karst et al. (2021) test set, from github.com/SorenKarst/longread_umi. It is used here and not redistributed.
- 4,259 MinION R10 reads (guppy 3.4.4, high accuracy) of 4.4-kb rRNA operons from the ZymoBIOMICS mock community.
- Each molecule carries an 18-nt UMI at both ends, just outside the degenerate primers 8F and 2490R. SHMtrace's anchor search now accepts IUPAC codes for primers like these.
- The truth: the 43 distinct operon sequences of the community.

**Method:** every read was aligned to the 43 operons. Its UMIs were read next to the primers with `umi.extract`, reads were grouped with `umi.group`, and each group of 3 or more became one consensus with `umi._consensus`, on the operon most of its reads align to. Error rates are measured against the best-matching operon.

| | Raw reads | UMI consensus |
|---|---:|---:|
| Reads with both UMIs read | 4,082 of 4,259 | |
| Molecules (3 or more reads; median 33) | | 105 |
| Error per aligned base, all kinds | 9.5% | 0.21% (46-fold lower) |
| Consensus matching its own operon best | | 100% |

**By reads per molecule:**

| Reads per molecule | Molecules | All errors | Substitutions | Molecules with no substitution |
|---|---:|---:|---:|---:|
| 3–4 | 6 | 2.2% | 0.28% | 0% |
| 5–9 | 7 | 0.45% | 0.023% | 43% |
| 10–19 | 15 | 0.10% | 0.011% | 53% |
| 20 or more | 77 | 0.053% | 0.010% | 56% |

**What this shows:**
- **Substitutions fall to about 1 in 10,000 bases from 10 reads per molecule.** Substitutions are the only error that matters for SHM calls. That is 100-fold below a wild-type J꜀4 mutation frequency of about 10⁻², and it holds with no control-based correction at all.
- **What's left is mostly indels in homopolymers.** They don't affect substitution calling.
- **The authors' pipeline** (two rounds of racon and medaka polishing) reaches a lower total error, mainly by fixing those indels. SHMtrace's column-wise consensus doesn't attempt that.
- **Grouping is reliable on real UMIs:** 95% of a group's reads align to the same operon copy, and the consensus always matches it best. The remaining 5% are reads that the near-identical operon copies make ambiguous.
- **Design guidance:** aim for 10 or more reads per molecule. Molecules with 3–4 reads still carry about 0.3% substitution error.

`tests/test_umi_karst.py` checks these numbers wherever the data are present.

## Other public data checked

- **Ramos IGHV4-34 SHM** (PLoS Comput Biol 2021, 10.1371/journal.pcbi.1009323): induced V-region SHM in Ramos cells with UMIs (MiSeq 2×300), plus a background control. The supplement (S4, S7) gives only the UMI consensus sequences as FASTA, with duplicate counts in S10. There are no raw reads or UMIs, so it can't test this step. It could serve as an error-corrected SHM case-against-control dataset for a V region.
- **R2C2+UMI human IGH** (PNAS Nexus 2024, PRJNA1005293): only processed consensus reads were deposited.

## Limits

- **No PCR errors were simulated.** An error in an early cycle is shared by many copies and survives consensus. Its rate with a high-fidelity polymerase (about 10⁻⁶ per base per doubling) is far below the nanopore error rate.
- **Real UMI data comes from 16S rRNA operons, not SHM.** No public dataset has raw UMI reads of a fixed SHM locus. The lab's own UMI JH4 libraries would be the decisive test.
- **Molecules with fewer than 3 reads are dropped.** That's the right trade for accuracy, but it costs depth.
