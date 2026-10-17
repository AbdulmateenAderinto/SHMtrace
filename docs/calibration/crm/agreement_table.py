"""Table S5: what "reproduced" means for each re-analysed study, on one set of questions.

    python crm/agreement_table.py        -> prints the rows (the manuscript builder imports them)

Figures S8-S10 each compare SHMtrace with a published analysis, but "agrees" means something different in
each: the same clones or not, the same endpoint or not, the same denominator or not, the same test or not.
This table asks the same five questions of all three, so a reader can see which claims are a like-for-like
reproduction and which only agree in direction.

Nothing here is computed from prose: the answers are the facts recorded in the study definitions
(sanger_compare.STUDIES, published/*.json) and in the runs themselves.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent

HEADER = ["Study (figure)", "Same clones?", "Same endpoint?", "Same denominator?", "Same test?",
          "Published values are", "Agreement"]

# One row per study. Each answer states the fact, not a verdict.
ROWS = [
    ("*Maea*^{fl/fl}*Mb1*^{cre}, Barbulescu et al. 2024 (Figure S8)",
     "Yes, the deposited chromatograms, each clone kept with the mouse it was filed under. SHMtrace sets "
     "paralog amplicons aside, so it analyzes 92/179 of the 97/187 clones the per-mouse tables count",
     "Yes for frequency and the G:C shares; the WRC/GYW and background columns are SHMtrace's only",
     "No. The paper divides by every filed clone × 577 bp; SHMtrace by the clones it analyzes × the 529-bp "
     "window they read cleanly. Table S4's \"counted as published\" row applies the paper's denominator to "
     "SHMtrace's mutations",
     "In the figure, yes: Figure S8F uses the one-tailed t the original study used. Table S4 reports "
     "two-tailed tests throughout, so its stars are the stricter of the two",
     "Source Data for frequency and the G:C shares (Fig. 2g–i); per-mouse tables (Supplementary Fig. 9b) for "
     "the rest",
     "Same direction and significance for frequency and both G:C classes. The Fig. 2g Source Data count 96/183 "
     "clones and the per-mouse tables 97/187, a discrepancy in the published record, not in the re-analysis"),
    ("*Mkln1*^{−/−}, Barbulescu et al. 2025 (Figure S9)",
     "Partly. Clones are kept with the mouse folder they were filed under, but of the three published "
     "knockout tables only mouse #1's spectrum matches a folder; #2 and #3 match no folder or set of clones "
     "in the project",
     "Yes for frequency and the G:C shares; WRC/GYW and background are SHMtrace's only",
     "No. The paper uses every filed clone × 556 bp; SHMtrace the clones it analyzes × 529 bp. An "
     "evidence-based clone set (Table S4) is reported as a sensitivity analysis",
     "In the figure, yes: Figure S9D uses the one-tailed t the original study used. Table S4 reports "
     "two-tailed tests throughout",
     "Computed from the per-mouse tables as printed (Supplementary Figure S2A). The paper's summary bars "
     "(Supplementary Figure S2C) give different values and are listed separately",
     "Same direction and significance for frequency and G:C transversions. The knockout spectrum agrees less "
     "well, which the untraceable tables #2 and #3 are enough to explain"),
    ("SW_{HEL} *Samhd1*, Thientosapol et al. 2018 (Figure S10)",
     "Yes in number, 1,795 deposited GenBank records, but not in provenance: the deposit names 12 wild-type "
     "hosts where the paper reports 13, and one record carries about twice every other host's clones",
     "Yes for the four substitution classes the paper tests; SHMtrace adds per-animal values the paper does "
     "not report",
     "No. The paper collates in a 560-bp window and counts every sequence independently; SHMtrace uses the "
     "consensus of the clones themselves and counts unique mutations per animal",
     "No for total mutation: the paper uses one-way ANOVA across its five arms, SHMtrace a Welch test on "
     "per-animal log frequency between two. The four class comparisons use a two-tailed t across animals",
     "The paper's Table 1, transcribed from the PMC copy",
     "All four substitution classes agree in direction and significance. Total mutation agrees in direction "
     "only, under two different tests, so it is not a like-for-like reproduction"),
]

NOTE = ("No deposit supplies raw traces for the SW_{HEL} study, so its chromatogram reading and bookkeeping "
        "checks are not exercised; it also has no *Aicda*^{−/−} arm. \"Same clones?\" asks whether SHMtrace saw "
        "the same molecules, \"same endpoint?\" whether the quantity is the one the paper reported, \"same "
        "denominator?\" whether it is divided the same way, and \"same test?\" whether the comparison is the "
        "one the paper made. Where any answer is no, agreement in direction is not a like-for-like "
        "reproduction, and the figures say which is which.")


def table_rows():
    return [list(r) for r in ROWS]


if __name__ == "__main__":
    for r in ROWS:
        print(r[0])
        for h, v in zip(HEADER[1:], r[1:], strict=True):
            print(f"   {h:22s} {v}")
