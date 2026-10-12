# nanoshm quickstart

From read files to a figure. The full account is in [`USER_GUIDE.md`](USER_GUIDE.md).

1. Open the app (or run `nanoshm-app`). It serves on your own machine only.
2. **Set up, box 1.** Pick the **Region** (`jh4_intron`, `n114_vregion`,
   `igh_smu`, `generic_amplicon`, `sanger_clones`), the **platform** (`ont`,
   `pacbio`, `illumina`, `sanger`) and leave the mode on `recommended`, which
   becomes the Sanger caller when the platform is `sanger`.
3. **Box 2.** Choose a **built-in reference**, or drop your own FASTA and your
   contaminant FASTA. Picking the built-in JH4 intron adds the JH4 contaminant
   for you.
4. **Box 3.** Drop the read files in. A folder becomes one sample, R1/R2 files pair up, and anything else is one sample per file. Press
   **Group by name**, then **Controls by name**, then read down the list and fix
   what the two buttons got wrong. Check the count: it should name the right
   number of controls.
5. **Box 4.** Name the project and press **Analyse**. A few minutes for six
   nanopore libraries.
6. **Results.** Check the headline line: the window and reference length, the
   reads per sample, and the control fold difference. Seven figures below.
7. **Tables.** Samples table: check `reads` is the same for every sample,
   `null_from` names controls rather than `none`, and `warnings` is empty.
   Groups table: the means, the SDs and the Welch and Mann-Whitney p-values
   against the control group.
8. **Export.** **Everything (.zip)** for the record. For a figure you will edit,
   download it as **SVG** or **PDF** from the figure list; PNG is for slides.

Things to check before believing a number:

- every sample analysed at the same depth, and no depth warnings;
- at least two controls ticked, three if you have them, all from the same
  sequencing run;
- `pct_GC` well above half and `pct_GC_in_WRC_GYW` high in the mutating group;
- the control group's density near zero;
- `clustered_calls` zero. Anything above zero means calls packed into a few
  bases, which is usually an artefact in that one library rather than
  hypermutation, and no control can correct for it;
- the spectrum looks like AID's. If the sample is repair-deficient (UNG or
  MSH2 knocked out), `pct_CT_GA` should be near 100.

With three samples per group the Mann-Whitney p cannot go below 0.1, so report
the Welch p and show every sample as a dot.

Two densities are reported. `unique_density_x1e5` counts only what the per-site
test called, and is the conservative, Geneious-comparable number.
`combined_density_x1e5` adds an estimate of the mutations carried by a single
read, which no per-site test can confirm at these depths, and comes with a 95%
interval. Use the first to compare with published work and the second when you
want the fullest estimate of mutation load.
