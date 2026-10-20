# Next steps for the all-usable-datasets pass

Branch `public-validation-2026-10`. Written 17 October 2026 so this can be picked up cold. The ledger with
every candidate and its status is `DATASET_LEDGER.md`; the two completed datasets are written up in
`PUBLIC_VALIDATION.md`.

## State at the last checkpoint (17 Oct, late)

Done and committed, all on branch `public-validation-2026-10`:

| | dataset | outcome |
|---|---|---|
| policy | — | `ancestral_reference`, tests first; 305 tests pass |
| E1 | Bassan TCR | reference reconstructed exactly; 11 libraries called; policy sensitivity measured |
| E5 | Devilder clones | 42/42 annotated clusters reproduce the published mutation lists |
| E2 | Liu base editors | motif separates AID from A3A where the transition share fails |
| E6 | Ma TAM, three arms | position follows the guide (8.6-369x at distinct sites); sites scale with guide number (0/32/37/66); single-control calibration shown to invert |
| E3 | Chen TRACE | **blocked**: archive resolves assay but not condition; evidence in the ledger |

Still open:

- **Run the TAM scripts in the foreground.** `crm/tam_reporter.py` and `crm/tam_footprint.py` both die
  silently part-way when started with `nohup ... &`, always at the point where merging starts, and both
  complete normally in the foreground. Their outputs now match the ledger. The cause was not chased down;
  if you background them, check the output rather than trusting the exit.
- The TAM **saturation arm** is done (`crm/tam_saturation.py`): 0, 32, 37 and 66 localised sites for 3, 5, 7
  and 10 guides. The `reporter-footprint` pair was inspected and rejected: UGI is confounded
  with read layout, read length and amplicon composition, with one library per condition (see the ledger).
  Nothing in SRP075898 remains analysable without the fragmented-read route.
- **7.7 GB of partly-fetched off-target runs sit in `~/Nanoshm/public_data/SRP075898/`** from before the
  download was restricted to the amplicon arms. They are not used by anything here and are re-fetchable from
  `manifest.tsv`; delete them if you want the space back. The run accessions are in `off_runs.txt` beside
  them. They were left in place rather than deleted unprompted.
- The TAM amplicon arms have finished downloading into `~/Nanoshm/public_data/SRP075898/` (3.5 GB of amplicon
  runs; the 17 GB off-target arm is deliberately excluded). The footprint and saturation arms are not yet
  analysed and would extend the localisation result with a guide titration.
- Nothing from this pass is in the manuscript yet. The main text is at 6,999 of 7,000 words, so any Results
  sentence needs words cut elsewhere; see the word-budget note in `PUBLIC_VALIDATION.md`.

## Completed in this pass## Completed in this pass## Completed in this pass

```bash
PY=~/Nanoshm/nanoshm_v0.2/.venv/bin/python
cd ~/Nanoshm/nanoshm_v0.2/docs/calibration
$PY crm/build_switch_reference.py && $PY -W ignore crm/switch_pacbio.py   # PRJNA831666 -> Figure S19
$PY -W ignore crm/public_smu_pi3k.py                                      # GSE77788    -> Figure S18
$PY -W ignore crm/germline_policy_demo.py                                 # synthetic, no download
```

Reads for the two public datasets are in `~/Nanoshm/public_data/PRJNA831666/` and
`~/Nanoshm/public_data/GSE77788/`. Nothing is pushed.

## The decision that blocks E1–E5

`DATASET_LEDGER.md` ends with a demonstrated two-part problem in the germline policy, with the exact code
locations. **No code has been changed.** Before any engineered-diversification dataset is analysed, the
authors need to choose between:

1. **Do nothing.** Then E1, E2, E5 (and the selected arms of E3, E4, E6) cannot report the variants selection
   enriched, and fixed variants are silently absent. Any analysis of those datasets would have to state that.
2. **Add an opt-in `ancestral_reference` policy**, as sketched at the end of the ledger, with regression
   coverage written first. This is a change to core calling code and needs the authors' agreement.

Option 2 is the one the prompt favours ("the smallest justified solution"), but it is not a change to make
unilaterally late in a revision.

## Downloads not yet requested

Nothing below has been downloaded. Approximate compressed sizes, in the prompt's execution order:

| Step | Dataset | Size |
|---|---|---|
| E1 | Bassan TCR (figshare) | ~448 MB |
| E2 | Liu 2018, SRP159221 | ~1.02 GB |
| E5 | Devilder supplement DOCX | small |
| E3 | Chen TRACE, PRJNA555784 | ~2.36 GB |
| E6 | Ma TAM, SRP075898 | ~19.3 GB |
| E4 | Roadmap, PRJNA649296 | ~35.4 GB |
| — | Lescale PRJNA1235282 | 3.83 GB |
| — | Marton PRJNA1304899 | 6.42 GB |
| — | Álvarez-Prado GSE102944 | ~22 GB |
| — | Wu PRJNA1337487 (amplicon subset) | ~50 GB |

That is roughly 140 GB before any blocked candidate is reconsidered. Check free disk before starting: the two
datasets already fetched take 0.6 GB.

## Order to resume in

1. **Settle the germline policy** (above). Everything engineered depends on it.
2. **E1, Bassan TCR** — smallest, and the clearest new scientific role. Required first:
   inspect the SAM headers for upstream filtering; reconstruct the recoded construct reference from the
   parental file and the paper rather than a germline TCR; resolve the founder/h5 relationships; check whether
   primer UMIs survive. Do not call h5 a matched AID-off control for all founders.
3. **E2, Liu 2018** — map labels to true biological replicates before anything else; confirm the "NC"
   construct from the paper. Keep APOBEC arms labelled separately from AID.
4. **E5, Devilder** — download the supplement, extract the nucleotide strings with provenance, analyse as
   clones against the parental construct. No population frequencies from a selected-clone set.
5. **E3/E6/E4** in that order, resolving library roles before downloading the large ones. For E4 in
   particular, find out which libraries sequence the target and which count guides; guide counts are not
   mutation observations and would waste most of the 35 GB.
6. **The remaining CSR sets** only if a second PacBio CSR example is wanted. Lescale and Marton are the same
   research group as PRJNA831666, so they do not add an independent laboratory.

## Standing constraints for whoever continues

- The main text is at 6,999 of 7,000 words. Any new Results sentence needs words cut elsewhere; see the word
  budget note in `PUBLIC_VALIDATION.md`.
- `switch.py`'s ambiguity flag is a confirmed defect, documented and deliberately unpatched
  (`PUBLIC_VALIDATION.md`, `HANDOVER.md` item 5). Do not quote `ambiguous_pct` from the module.
- Keep reads outside the repository, in `~/Nanoshm/public_data/<accession>/`, and do not redistribute them.
- Do not treat reads, repeat runs, amplicons or selected descendants as independent biological replicates.
- Do not contact authors, push, or buy anything.
