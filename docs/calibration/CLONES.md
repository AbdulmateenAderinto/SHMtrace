# Clones from the V(D)J junction

*12 October 2026. Data: the six bundled J꜀4 libraries (497 reads each).*

## The problem

Without UMIs every read counts once. A clone that expanded in the germinal centre, or a template PCR copied many times, counts as often as it was read. The project review called template counting its main unresolved problem.

## The idea, and where it came from

The LDHA study validated in [`ILLUMINA.md`](ILLUMINA.md) collapsed its reads "utilizing the VDJ recombination as a cellular barcode". The J꜀4 amplicon carries the same barcode. It starts in the V region, so the bases just 5' of the J꜀4 coding segment are the V(D)J junction, made by recombination and unique to each B-cell clone.

On the germline reference those bases are germline sequence, which a rearranged read does not carry. So the aligner starts a rearranged read's alignment at J꜀4 (reference 53), and the unaligned bases before it are the junction. In Wt3, 281 of 497 reads begin aligning at 50–59, with a median 76 bp of unaligned sequence ahead.

`clones.py` takes the 20 bases before the alignment start of every read that begins within 18 bp of J꜀4. It clusters them greedily, most abundant first, joining a junction to a clone centroid within edit distance 3 (centroids, not chains, so near neighbours cannot merge two clones). Each clone's mutations are the called mutations carried by at least half of its reads that cover the position. A clone counts once.

This gives counts per clone, not per molecule. Two cells of one clone, or two PCR copies of one cell's template, look the same. Nothing short of a UMI separates those. For the question that matters here — is one lineage dominating the library? — the clone is the right unit, and it needs no change to the protocol.

## Is the junction really a barcode?

If it is, reads sharing a junction should share their mutations: they inherited them. Reads with different junctions should share almost none.

| Library | Reads with a junction | Clones | Mutation sharing (Jaccard) within a clone | Between clones | Ratio |
|---|---:|---:|---:|---:|---:|
| Wt3 | 352 | 41 | 0.510 | 0.024 | **22×** |
| Wt4 | 258 | 18 | 0.873 | 0.003 | **265×** |
| Wt2 | 391 | 28 | 0.004 | 0.002 | 2× (see below) |

For Wt3 and Wt4 the junction separates lineages cleanly. The knockouts carry almost no mutations, so there is nothing to share.

## What it found: Wt2 is one clone, mostly

Wt2's junctions are dominated by one sequence: 322 of 329 reads in its largest cluster carry the identical 20-mer (`TTACTACGGTAGTAGCTAGG`, a D segment joined to J꜀4 through two N bases). Identity at 98% rules out a clustering artefact; sequencing errors would scatter it.

**That clone is almost unmutated.** 6% of its reads carry any called mutation, against 56% of Wt2's other reads; 0.07 mutations per read against 5.03. It makes up 84% of Wt2's reads with a junction, and it is why Wt2's within-clone sharing looks low: an unmutated clone has nothing to share.

Whether it is a genuinely expanded clone or a PCR jackpot of one cell's template, it is one lineage, and read-weighted numbers count it 330 times.

## Read-weighted against clone-weighted

| Library | Mutations per 10³ bp, per read | Per clone | Largest clone |
|---|---:|---:|---:|
| Wt2 | 1.74 | **10.67** | 84% |
| Wt3 | 21.34 | 24.25 | 21% |
| Wt4 | 4.78 | 4.62 | 13% |
| KO1 | 0.00 | 0.00 | 26% |
| KO2 | 0.13 | 0.17 | 17% |
| KO3 | 0.02 | 0.06 | 14% |

Wt2 moves sixfold; the libraries without a dominant clone barely move. Wild type against *Aicda*-/-, Welch on log values: **P = 0.0033** clone-weighted against P = 0.0086 read-weighted, because Wt2 is no longer an outlier dragging the wild-type group down.

## What is reported

For presets with `junction_ref_start` (J꜀4: 53), every project now reports `reads_with_junction`, `clones`, `largest_clone_pct` and `clone_weighted_frequency_x1e3` per sample, with group statistics. It also writes `clones.csv` (every clone: reads, junction, its mutations) and a figure with each clone as one dot. A sample whose largest clone holds half its junction reads or more gets a warning to compare the clone-weighted frequency.

Nothing here changes the calls. The read-weighted numbers stay as they were; the clone-weighted one sits alongside.

## Limits

- **Per clone, not per molecule.** A UMI is still the only way to count templates.
- **Only amplicons that start in the V region.** The N114 hybridoma is monoclonal by design and has one junction; the Sμ and generic presets have none, so the setting is off for them.
- **Reads that do not reach the junction are left out of the clone count.** In these libraries 22–79% of analysed reads carry one; the per-read numbers still use every read.
- **The 50% consensus rule** treats a mutation in fewer than half of a clone's reads as subclonal. For a large clone that is under-counting real subclonal mutations; for a small one, a single read's error could pass it. A lineage tree (IgTreeZ) would place each mutation at its node instead; this is the trunk only.
