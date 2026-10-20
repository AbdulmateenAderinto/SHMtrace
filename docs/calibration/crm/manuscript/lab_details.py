"""Laboratory details only the authors can supply, in one place for both papers and the cover letters.

Each entry is (description, value). Leave the value as None until it is known: the text then shows
"[description]" and `python lab_details.py` lists what is still missing (a value holding "[...]" is
listed too, as only partly known). Values taken from the lab's published methods (Feng et al., 2021;
Barbulescu et al., 2024, 2025) are marked "# paper". The texts refer to an entry
as {lab_<key>}.
"""

DETAILS = {
    # contacts
    "tech_email": ("technical contact email", "a.aderinto@mail.utoronto.ca"),
    # antibodies (key resources table)
    # paper: Barbulescu 2024 Nat Commun, 2025 J Immunol (Peyer's patch panel)
    "b220_source": ("B220 supplier", "SouthernBiotech"),
    # vendor catalogue; RRID from a distributor listing, check at rrid.site
    "b220_id": ("B220 Cat#; RRID", "Cat# 1665-02; RRID: AB_2795114"),
    # paper: Feng 2021 (clone GL7, eBioscience); the 2024/2025 papers print 'RA3-6B2; SouthernBiotech', a typo
    "gl7_source": ("GL7 supplier", "Thermo Fisher Scientific (eBioscience)"),
    # vendor catalogue (100 ug size)
    "gl7_id": ("GL7 Cat#; RRID", "Cat# 50-5902-82; RRID: [to confirm]"),
    # paper: Barbulescu 2024, 2025
    "fas_source": ("Fas supplier", "BD Biosciences"),
    # BD product page
    "fas_id": ("Fas Cat#; RRID", "Cat# 562633; RRID: AB_2737690"),
    # reagents and kits
    # Platinum SuperFi II PCR Master Mix, 100 reactions (the STAR table uses the 2x master mix)
    "superfi_id": ("SuperFi II Cat#", "12368010"),
    # paper: Feng 2021; Barbulescu 2024, 2025
    "protk_source": ("proteinase K supplier", "Invitrogen (Thermo Fisher Scientific)"),
    "protk_id": ("proteinase K Cat#", None),
    "monarch_id": ("Monarch kit Cat#", None),
    # mice
    # paper: Barbulescu 2025 (the JIM mice are that study's cohort)
    "wt_source": ("C57BL/6 source", "Bred in house (wild-type littermates)"),
    "wt_rrid": ("C57BL/6 RRID", "N/A (bred in house)"),
    "aicda_rrid": ("Aicda-/- RRID (MGI)", None),
    # per-mouse sex: the Methods give "both sexes"; reporting guidelines ask for the sex of each animal
    "mouse_sex": ("sex of each of the 6 J_{H}4 mice (e.g. Wt2 M, Wt3 F, ...)", None),
    # cell lines: authentication and mycoplasma status, or an explicit statement that they were not done
    "n114_auth": ("how the N114 hybridoma was authenticated, or 'not authenticated'", None),
    "n114_myco": ("N114 mycoplasma testing result and date, or 'not tested'", None),
    # primers
    # Barbulescu 2025 Table S3 (Jh4-F/Jh4-R, with EcoRI/SpeI tails); 85%/81% of Wt3 nanopore reads carry them
    "jh4_primers": ("J_{H}4 primer sequences",
                    "forward 5′-GGAATTCGCCTGACATCTGAGGACTCTGC-3′, reverse 5′-GACTAGTCCTCTCCAGTTTCGGCTGAATCC-3′"),
    # EPI2ME folder vprimers.fasta; they bound the 655-nt N114 reference exactly
    "n114_primers": ("N114 primer sequences", "forward 5′-TTACCTGGGTCTATGGCAGT-3′, reverse 5′-TGAAGGCTCAGAATCCCCC-3′"),
    # bench details (STAR protocol)
    # paper: Barbulescu 2025 Peyer's patch panel; its gating text says FVD eF780, its panel eF506
    "stain": ("antibody dilutions and fluorophores", "anti-B220–FITC (RA3-6B2, 1:400), GL7–eFluor 660 (GL7, 1:100), and anti-Fas–BV421 (Jo2, 1:200), with anti-CD16/32 Fc block (clone 93, 1:100) and Fixable Viability Dye eFluor 506 (1:1,000)"),
    "facs_buffer": ("FACS buffer composition, e.g. PBS with 2% fetal bovine serum", None),
    # paper: Barbulescu 2024, 2025
    "lysis": ("lysis buffer, volume, and incubation", "55 °C overnight, then 95 °C for 15 min; [lysis buffer and volume to confirm]"),
    "lysate_volume": ("lysate volume per PCR", None),
    # deposits
    "sra": ("SRA BioProject", None),
    "data_doi": ("data deposit DOI", None),
    "code_doi": ("SHMtrace Zenodo DOI", None),
    "github": ("GitHub URL", "https://github.com/AbdulmateenAderinto/SHMtrace"),
    "crm_citation": ("companion Cell Reports Methods article, to be cited once accepted", None),
}


def placeholders():
    """{lab_<key>: value or "[description]"} for crm_numbers.compute()."""
    return {f"lab_{k}": v if v is not None else f"[{d}]" for k, (d, v) in DETAILS.items()}


def missing():
    """Entries still unknown, and entries filled only in part (their value still holds a "[...]")."""
    return [d for d, v in DETAILS.values() if v is None or "[" in v]


if __name__ == "__main__":
    left = missing()
    print(f"{len(DETAILS) - len(left)} of {len(DETAILS)} details filled; still missing:")
    for d in left:
        print(" -", d)
    if any("DOI" in d or "SRA" in d or "accession" in d.lower() for d in left):
        print("\nWhile a deposit identifier is unresolved, Data and code availability says the files "
              "\"will be deposited\". Switch those two sentences to \"have been deposited\" only once "
              "the identifiers resolve.")
