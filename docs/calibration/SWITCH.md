# Class-switch junctions from long reads

*12 October 2026. `nanoshm/switch.py`, `nanoshm-switch`; tests in `tests/test_switch.py`.*

## What it measures

Class-switch recombination joins a break in Sμ to a break in an acceptor switch region and deletes what lies between. A long read across the join aligns in two parts: its 5′ part to Sμ and its 3′ part to the acceptor. For each read the module reports:

| | |
|---|---|
| **Breakpoints** | The last Sμ base and the first acceptor base, in each region's own coordinates. |
| **Microhomology** | Bases at the join that match both sides. Blunt joins have none; 1–3 nt is short; 4 nt or more is characteristic of alternative end-joining, as when classical NHEJ (Ku, XRCC4, Lig4) is impaired. |
| **Insertion** | Bases at the join that match neither side. |
| **Mutations near the join** | Q ≥ 20 differences within 50 bp of the join, per base read. AID-induced mutations cluster there. |
| **Ambiguity** | The acceptor part aligns equally well elsewhere in the repetitive acceptor, so its breakpoint is not unique. |
| **Internal deletions** | Both parts land in Sμ, more than 200 bp apart. |

Reads with breakpoints within 3 bp on both sides are one event: PCR copies of one recombined cell. They are counted once as a unique junction. Per sample, the module reports acceptor use, microhomology classes, mean microhomology, the ambiguous share and mutations near the join. Groups are compared animal by animal (Welch and permutation).

```bash
nanoshm-switch --reference switch_regions.fa --donor Smu \
    --sample WT1=wt1.fastq:WT --sample WT2=wt2.fastq:WT --sample KO1=ko1.fastq:KO ... \
    --compare KO,WT -o switch_out
```

The reference is one FASTA record per switch region (Sμ and the acceptors), taken from the genome assembly the lab uses. No switch-region references are bundled beyond the Sμ fragment.

## Accuracy on simulated junction reads

The simulation uses:
- switch-like regions of 3 kb, random sequence studded with each region's repeat motif;
- junction reads with known breakpoints, of three kinds in equal numbers: random joins (mostly blunt or 1–3 nt), joins chosen to have 4–8 nt microhomology, and joins with 2–5 nt insertions;
- substitution and indel errors throughout, and half of the reads reverse-complemented.

| Read error (sub + indel) | Switched reads found | Breakpoints within microhomology + 3 bp |
|---|---|---|
| 2% (R10.4.1-like) | 460 / 460 | 95% |
| 6% (older chemistry) | 469 / 469 | 91% |

Microhomology class, true (rows) against called, at 2% read error:

| True \ called | Blunt | 1–3 nt | ≥ 4 nt | Insertion |
|---|---:|---:|---:|---:|
| Blunt | **83%** | 2% | 7% | 9% |
| 1–3 nt | 3% | **80%** | 8% | 9% |
| ≥ 4 nt | 0% | 1% | **97%** | 1% |
| Insertion | 5% | 7% | 3% | **86%** |

At 6% read error, the share called correctly falls to 69% (blunt), 68% (1–3 nt), 82% (≥ 4 nt) and 90% (insertion).

- **What's reliable:** long microhomology, the alternative end-joining signature, is classed correctly 97% of the time at current nanopore accuracy, and breakpoints are placed within a few bases.
- **The main error:** blunt and short joins are sometimes called 1–2 nt insertions, because a sequencing indel right at the join looks exactly like one. Comparisons between groups sequenced the same way are fair, since the error affects them equally, but absolute insertion rates are overestimated. A consensus per molecule (UMIs) would remove these errors.

## Limits

- **Simulated only.** No real switch-junction reads have been tested yet. A published dataset of Sμ–Sγ1 junctions, or the lab's own CSR PCRs, would be the real test.
- **Repetitive acceptors:** a breakpoint inside a long repeat array may be placed on the wrong repeat copy. Such reads are flagged as ambiguous, and the flag is reported per sample.
- **UMI consensus is not yet available for junction reads.** It works for amplicons against one reference.
