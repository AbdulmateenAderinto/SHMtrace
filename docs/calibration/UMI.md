# UMI consensus: one read per molecule

*12 October 2026. `nanoshm/umi.py`; simulation `docs/calibration/umi_sim.py`, output in `umi_sim/`.*

## What it does

A unique molecular identifier (UMI) is a random tag put on each template molecule before PCR, so all reads sharing a UMI are copies of one molecule. With `umi` on, nanoshm:

1. **Reads each read's UMI.**
   - `flank`: the UMI is next to an anchor at a read end. By default the anchor is the first or last 20 reference bases (the primer), found with up to 6 edits. Use `umi_ends` = left, right or both for dual UMIs.
   - `header`: the UMI is in the read name, as umi_tools and similar tools write it.
2. **Groups reads by UMI**, allowing `umi_max_edits` (3) edits per UMI end. Groups are built greedily by abundance against centroids. A UMI seen once or twice may join a group at least three times its size within two more edits, as UMI-tools' "directional" method allows, because it is far likelier an error copy than a molecule of its own.
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

## Limits

- **No PCR errors were simulated.** An error in an early cycle is shared by many copies and survives consensus. Its rate with a high-fidelity polymerase (about 10⁻⁶ per base per doubling) is far below the nanopore error rate.
- **No real UMI data has been tested yet.** Public nanopore UMI amplicon data (for example the dual-UMI method of Karst et al., 2021) would be the next test.
- **Molecules with fewer than 3 reads are dropped.** That's the right trade for accuracy, but it costs depth.
