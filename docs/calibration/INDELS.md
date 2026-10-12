# Insertions and deletions

*11 October 2026. Data: the six bundled J꜀4 libraries (`examples_compare/`, run DHJGKB, 497 reads each). Off by default: turn on with `exclude_indels: false`, or "Call insertions and deletions" under Advanced in the app.*

## Why Geneious left them out, and what changes

Nanopore reads carry about one indel error each. On J꜀4 the knockouts show 1.05–1.10 per read, 41–46% of them in homopolymers, and about 225 distinct one-read indels per library. Read quality cannot separate a real indel from that.

But the errors are systematic. The same deletion in the same homopolymer recurs at a similar rate in every library. A mutation-free control therefore measures, for each exact indel, how often it appears by error (`nanoshm/indels.py`):

1. Every read's indels are left-aligned. An aligner may place a deletion anywhere along a repeat, and left-aligning gives one representation per event.
2. Only indels of quality ≥ Q30 are counted, in controls and samples alike. Quality here means the weakest inserted base, or the weaker base flanking a deletion. Nanopore indel errors mostly come with low-quality flanks, while a real indel is read cleanly.
3. Each indel's error rate is the controls' rate for that exact event. Where the controls cover it thinly, the rate is pulled toward the rate for its class (insertion or deletion, length up to 4, homopolymer run up to 5), and never goes below the class rate.
4. The rates are inflated by a factor fitted by testing each control against the others, as for the quality-calibrated substitution model. The calls are then corrected by Benjamini-Hochberg.
5. An indel needs at least 3 supporting reads.

Indels are reported in their own table (`indels.csv`) and in their own columns: `indels_called`, `insertions`, `deletions`, `frameshift_indels` and `indel_density_x1e5`. They are not added to the substitution density, so densities stay comparable with the paper's.

Under the Sanger platform (one read per clone, low indel error), every indel with quality ≥ Q20 is called, as substitutions are.

## Spike-in

`indel_spike.py` plants 10 indels per draw into one knockout. Kind and number of reads are drawn independently, and planted reads keep their own base qualities. The library is called against a null fitted and calibrated on the other two knockouts only. 8 draws per knockout, 240 planted indels.

Recall by supporting reads, at the defaults:

| | 1 | 2 | 3 | 5 | 10 | 20 |
|---|---:|---:|---:|---:|---:|---:|
| all kinds | 0.00 | 0.00 | 0.43 | 0.72 | 0.94 | 0.97 |

| Kind | Planted at ≥3 reads | Recall |
|---|---:|---:|
| 1-bp deletion, outside a homopolymer | 21 | 0.76 |
| 1-bp deletion, in a homopolymer (run ≥ 3) | 16 | 0.19 |
| 1-bp insertion | 18 | 0.72 |
| 2-bp deletion | 35 | 0.80 |
| 2-bp insertion | 24 | 0.96 |
| 3-bp deletion | 32 | 0.84 |
| 10–40-bp deletion | 19 | 0.95 |

False calls: 0.71 per library. They are three distinct events, repeated across draws:
- a TT insertion in a TTT homopolymer at 246 in KO3 (3 reads);
- a TT insertion at 433 in KO2 (3 reads), at one of the germline-difference sites, where reads with the germline haplotype may align awkwardly;
- a 1-bp insertion in 5 reads in KO3.

Each is present in one knockout only, so no background fitted from the other libraries can predict it.

How the settings were chosen:

| Indel quality ≥ | Minimum reads | Recall at 3 / 5 / 10 reads | False calls per library |
|---|---|---|---|
| 0 | 1 | 0.20 / 0.52 / 0.74 | 0.58 |
| 20 | 1 | 0.43 / 0.61 / 0.91 | 0.71 |
| 30 | 1 | 0.60 / 0.67 / 0.89, plus some 2-read events | 1.08 |
| **30** | **3** | **0.43 / 0.72 / 0.94** | **0.71** |

The quality floor is what makes small indels detectable. At Q0, a 1-bp deletion outside a homopolymer needed 20 reads; at Q30 it needs 3 to 5. The 3-read floor removes two-read events that are real molecules in a knockout, such as a 28-bp deletion in 2 reads, probably made in PCR. It counts reads at ≥ Q30, so a planted indel in 3 reads is only called when all three pass the quality floor. That is why recall at 3 reads falls from 0.60 to 0.43.

## Real libraries

| Library | Indels called |
|---|---|
| Wt2 | 0 |
| Wt3 | 10: deletions of 1 bp (109, 212, 441, 443), 14 bp (156), 33 bp (440), 156 bp (205) and 195 bp (450); insertions of G (446) and AC (475). Each in 5–22 reads, none seen in any knockout |
| Wt4 | 0 |
| KO1 | 0 |
| KO2 | 1 (the TT insertion at 433, 3 reads) |
| KO3 | 1 (the TT insertion at 246, 3 reads) |

Wt3 is the most mutated and most clonal library. Its large deletions fit AID-initiated deletions carried by an expanded clone. The 1-bp events at 441 and 443 sit beside the 33-bp deletion that starts at 440, so the same clone may be carrying them.

## Limits

- **Controls are required**, at least one besides the sample itself. With fewer than two controls, the inflation factor cannot be fitted properly. With exactly two, at Q30, it reaches the top of its grid (50).
- **1-bp deletions in homopolymers** are the commonest nanopore error and stay mostly undetectable.
- **Strand.** The bundled example reads are all one strand, so strand balance cannot be checked here. The `plus_reads` and `minus_reads` columns allow it on mixed-strand data.
- **One run, one locus, and a sequencing-only check.** The background is calibrated per run. The planted indels are idealised. Real indels should be checked against Sanger clones of the same samples. That is the next calibration this needs.

## Reproduce

```bash
python docs/calibration/indel_spike.py 8        # defaults: Q30, 3 reads
python docs/calibration/indel_spike.py 8 0      # no quality floor
```
