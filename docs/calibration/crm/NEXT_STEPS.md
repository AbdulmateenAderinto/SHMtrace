# Next steps for the all-usable-datasets pass

Branch `public-validation-2026-10`. Written 17 October 2026 so this can be picked up cold. The ledger with
every candidate and its status is `DATASET_LEDGER.md`; the two completed datasets are written up in
`PUBLIC_VALIDATION.md`.

## State at the last checkpoint (17 Oct, evening)

Done and committed:

- **Germline policy.** `ancestral_reference` implemented, tests first
  (`tests/test_ancestral_reference.py`, six tests); 305 tests pass, `ruff check src` clean. The patch table
  and the scope deliberately left out are in `DATASET_LEDGER.md`.
- **E1 Bassan TCR.** Reference reconstructed exactly from CIGAR+MD (1,785 bp, zero disagreements); eleven
  libraries called; policy sensitivity measured. `crm/out/tcr/`.
- **E5 Devilder clones.** All 42 annotated clusters reproduce the authors' nucleotide mutation lists
  exactly. `crm/out/devilder/`.

In flight:

- **E2 Liu SRP159221.** All 386 files downloaded (1.02 GB). `crm/liu_mutators.py` was running at the
  checkpoint; rerun it if `crm/out/liu/liu_mutators.json` is absent. The prespecified contrast is fixed in
  the script's `ARMS`: DBE-AIDmono against DBE-humanA3a, in wild-type and UNG/MSH2-null 293T, with the NC
  arm as control. Pairs merge cleanly here (97.8%, a 238-bp EGFP amplicon), unlike GSE77788.

Not started: **E3** (Chen TRACE, PRJNA555784, 2.36 GB, 55 paired runs grouped as On-target / Targeted /
MEK1 screen / Off-target / T7 variants), **E6** (Ma TAM, SRP075898, 19.3 GB).

Two traps already hit, worth not repeating:

- `curl` inside a `while read` loop eats stdin and silently truncates a download. Both manifests carry
  expected sizes; the loops redirect `curl < /dev/null` and verify size. Downloading with
  `xargs -P 8` is roughly ten times faster than the sequential loop.
- The SRP159221 labels are inconsistent about separators (`DBE_AID_mono` against `DBE-AIDmono `, some with a
  trailing space). `crm/liu_mutators.py` lists them exactly as deposited; do not normalise them blindly.

## Completed in this pass## Completed in this pass

```bash
PY=~/Nanoshm/nanoshm_v0.2/.venv/bin/python
cd ~/Nanoshm/nanoshm_v0.2/docs/calibration
$PY crm/build_switch_reference.py && $PY -W ignore crm/switch_pacbio.py   # PRJNA831666 -> Figure S16
$PY -W ignore crm/public_smu_pi3k.py                                      # GSE77788    -> Figure S17
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
