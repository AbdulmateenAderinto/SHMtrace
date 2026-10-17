# Simulated libraries with known clones

*12 October 2026. `shmtrace-simulate-lineages` (`shmtrace/lineage_sim.py`). Backgrounds: the three J꜀4 Aicda-/- libraries, each planted in turn with the other two as the null. One draw per control per scenario; model HH_S5F (see [`TARGETING.md`](TARGETING.md)).*

## Why

The spike-in test (`shmtrace-validate`) plants mutations at chosen read counts. It cannot test anything that depends on clones: whether the junction finds them, how much an expanded clone biases the read-weighted frequency, and whether the clone-weighted number gets the truth back. No real library has a known clone structure, so this simulation builds one.

## How a library is built

In the spirit of SHMulate (Yaari et al. 2013), but on real reads, so the sequencing errors are the sequencer's own:

1. The control's reads are dealt into about 40 clones with power-law sizes. In the *expanded* scenario one clone, left unmutated, takes 60% of the reads, as in Wt2 ([`CLONES.md`](CLONES.md)).
2. Each clone gets a random 20-bp junction and its own J trimming (−3 to +12 bp, the range the controls' alignments show). Every read with a junction is rebuilt so that the clone's junction meets J at that point, then given substitutions at the read's own error rate and indels at 0.2% per base, the rate measured in these controls.
3. Each clone gets a lineage. Trunk mutations go into all its reads, with the load varying between clones (gamma, shape 2). A branch goes into 30–60% of the reads, a sub-branch into part of that, and private mutations into single reads. Positions come from the 5-mer model's mutability, and substitutions from its substitution probabilities.
4. The library runs through `run_project` like any other and is scored against what was planted.

| Scenario | Clones | Trunk / branch / private (mean mutations) | Dominant clone |
|---|---:|---|---|
| polyclonal | 40 | 6 / 2 / 0.3 | none |
| expanded | 40 | 6 / 2 / 0.3 | 60% of reads, unmutated |
| expanded_mutated_40 | 40 | 6 / 2 / 0.3 | 40% of reads, mutated (trunk 6) |
| expanded_mutated_60 | 40 | 6 / 2 / 0.3 | 60% of reads, mutated (trunk 6) |
| expanded_mutated_80 | 40 | 6 / 2 / 0.3 | 80% of reads, mutated (trunk 6) |
| low | 40 | 1 / 0.5 / 0.1 | none |
| null | 40 | nothing planted | none |

## Calls

Recall by the number of analysed reads carrying a planted mutation (mean of three libraries), and calls that were never planted:

| Scenario | Mode | 1 read | 2 | 3–5 | 6–10 | >10 | False calls per library |
|---|---|---:|---:|---:|---:|---:|---:|
| polyclonal | recommended | 0.38 | 0.75 | 0.99 | 0.98 | 0.99 | 3.3 |
| polyclonal | geneious | 0.56 | 0.84 | 0.99 | 0.99 | 0.99 | 2.3 |
| expanded | recommended | 0.35 | 0.68 | 0.95 | 0.96 | 0.99 | 4.0 |
| expanded | geneious | 0.59 | 0.83 | 0.97 | 0.97 | 0.99 | 3.3 |
| low | recommended | 0.01 | 0.58 | 0.99 | 0.98 | 1.00 | 1.7 |
| low | geneious | 0.61 | 0.85 | 0.97 | 0.98 | 1.00 | 3.3 |
| **null** | **recommended** | – | – | – | – | – | **1.3** |
| **null** | **geneious** | – | – | – | – | – | **4.0** |

Anything in three or more reads is found in both modes. The modes differ on one- and two-read mutations, and they split the trade the way the spike-in test found: the recommended model makes a third of Geneious's false calls in a library with nothing in it, and recalls fewer rare mutations.

**One result goes against the recommended mode, and it is reported here as found.** Sweeping the threshold (`roc.png`; recall of mutations in two or more reads against false calls) puts Geneious's curve slightly *above* the recommended one in the polyclonal and expanded libraries. At about one false call per library, Geneious's ordering recalls 0.96 of the polyclonal mutations against 0.92. The recommended test discounts positions where the controls themselves err, and a real mutation that lands on such a position pays for it. Geneious does not discount them, and it pays instead in the null, where its false calls are concentrated. Which costs more depends on how mutated the library is, so this is an open question about the recommended model's quality inflation (`qscale`), not a reason to change it on one simulation. The planted bases also carry the donor base's quality, which flatters every test equally.

The expanded clone costs the recommended mode recall (0.88 of mutations in two or more reads, against 0.96 polyclonal). The 60% unmutated clone adds depth at every position, so the same number of mutant reads is a smaller fraction against a larger expected error count.

## Clones

| Scenario | Adjusted Rand | Clones found / present | Largest clone, found / true |
|---|---:|---:|---:|
| polyclonal | 0.97 | 43 / 35 | 21% / 21% |
| expanded | 0.96 | 37 / 29 | 58% / 59% |
| low | 0.96 | 46 / 35 | 20% / 21% |
| null | 0.95 | 45 / 35 | 20% / 21% |

The junction recovers the clones: an adjusted Rand index of 0.95–0.97, against 1 for perfect agreement and 0 for chance. The largest clone is measured to within a point. Clone counts run about 8 high. The extras are small splinters made by junction errors, and they barely touch the clone-weighted frequency.

The simulation also improved the clustering. With error-free junctions the clusterer was already near perfect (0.97), and the residual came from a flaw in the distance. The 20-bp window's 5' edge is an arbitrary cut, so an indel inside it slides one extra base in at that edge, and plain edit distance charged that as a second error. Leaving the 5' edge free changed nothing at the real indel rate and raised recovery from 0.96 to 0.99 at ten times that rate. It is now the clusterer's distance. On the real libraries it merged a few splinters: Wt3 41 → 37 clones, with sharpened barcode statistics (`CLONES.md` has the re-run numbers).

## Frequencies: does clone-weighting recover the truth?

True against estimated, mutations per 10³ bp (true values computed the way each estimate is defined):

| Scenario | Read-weighted, est / true | Clone-weighted, est / true |
|---|---:|---:|
| polyclonal | 16.0 / 16.4 | 13.0 / 11.7 |
| **expanded** | **5.4 / 5.8** | **10.6 / 12.4** |
| low | 2.1 / 2.2 | 2.4 / 2.1 |
| null | 0.05 / 0 | 0.06 / 0 |

Both estimators are accurate for what they measure. What matters is the expanded row: the read-weighted number is *correctly* 5.8, because 60% of the reads really are unmutated, while the clones carry 12.4. Clone-weighting recovers 85% of that, against 46% for the read-weighted number. That is the Wt2 situation, now with known truth: 1.74 per read against 12.00 per clone there is the same distortion, and the per-clone number is the one to believe about the lineages.

The clone-weighted estimate runs 10–15% off either way: low where the 50% rule drops branch mutations that missed being called, high where splinter clones keep a single read's calls.

## Targeting

| Mutations placed by | S5F targeting index |
|---|---:|
| the 5-mer model | 1.48–1.79 |
| WRC/GYW and WA/TW weights only | 1.23 |
| real wild-type calls (`TARGETING.md`) | 1.19–1.34 |

An index near 1.3, as in the real libraries, is what motif-level AID preference gives under this model. The model's own placement reaches about 1.5.

## Running it

```bash
shmtrace-simulate-lineages --region jh4_intron \
  --control KO1.fastq --control KO2.fastq --control KO3.fastq \
  --reference mouse_jh4_intron.fasta --contaminant mouse_jh4_contaminant.fasta \
  --targeting-model HH_S5F.csv --out lineage_sim
```

`--scenarios`, `--modes`, `--reps` and `--place motif` vary the design. It writes `report.md`, `lineage_sim.json` (every run, every lineage) and `roc.png`. The 24-run suite above takes about 8 minutes.

## What the high-allele-fraction rule costs

A change carried by at least 50% of the reads covering it is set aside as a germline or parental difference.
A mutation shared by a dominant lineage reaches the same frequency, so the `expanded_mutated_*` scenarios
measure what that rule removes. Each planted mutation is counted as exposed when the reads carrying it reach
the threshold at its position; `planted_at_mask` and `planted_at_mask_called` are reported per run.

| Dominant clone | Planted mutations exposed, per library | Recovered | Read-weighted frequency | Clone-weighted |
|---|---:|---:|---|---|
| 40% of reads, mutated | 0 | – | within 5% of truth | within 14% |
| 60% of reads, mutated | 0–4 | none | down to 40% below truth | down to 26% below |
| 80% of reads, mutated | 0–4 | none | down to 62% below truth | down to 34% below |

Three libraries per abundance, each of the three J꜀4 AID-null controls planted in turn. Nothing that crosses
the threshold survives: the rule runs before any test. The loss is bounded by the dominant clone's own trunk,
not by the size of the library, and it falls hardest on the read-weighted frequency, which those trunk
mutations dominate. Libraries where nothing crossed the threshold were unaffected.

Every site the rule removes is reported per sample as `high_fraction_sites` (position, alternate base, reads,
depth and frequency), with `high_fraction_sites_removed` as the count. A true germline difference sits near
100% (the five J꜀4 lab-haplotype sites are at 99.7%), whereas a dominant lineage's trunk sits near that
clone's share of the reads, so the two can usually be told apart by inspection.

## Limits

- **The planted side is idealised.** A planted base has the donor read's quality and no alignment artefact of its own. False calls are real; recall is an upper bound.
- **The lineages are simple:** a trunk, one branch and one nested sub-branch, which is shallower than a real germinal-centre tree.
- **No PCR.** There is no amplification bias or PCR error, and no chimeras. A clone's reads are independent draws.
- **Three backgrounds, one draw each.** The differences between modes are small next to the spread between controls.
