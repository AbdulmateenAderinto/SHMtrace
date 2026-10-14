"""Laboratory details only the authors can supply, in one place for both papers and the cover letters.

Each entry is (description, value). Leave the value as None until it is known: the text then shows
"[description]" and `python lab_details.py` lists what is still missing. The texts refer to an entry
as {lab_<key>}.
"""

DETAILS = {
    # contacts
    "tech_email": ("technical contact email", "a.aderinto@mail.utoronto.ca"),
    # antibodies (key resources table)
    "b220_source": ("B220 supplier", None),
    "b220_id": ("B220 Cat#; RRID", None),
    "gl7_source": ("GL7 supplier", None),
    "gl7_id": ("GL7 Cat#; RRID", None),
    "fas_source": ("Fas supplier", None),
    "fas_id": ("Fas Cat#; RRID", None),
    # reagents and kits
    "superfi_id": ("SuperFi II Cat#", None),
    "protk_source": ("proteinase K supplier", None),
    "protk_id": ("proteinase K Cat#", None),
    "monarch_id": ("Monarch kit Cat#", None),
    # mice
    "wt_source": ("C57BL/6 source", None),
    "wt_rrid": ("C57BL/6 RRID", None),
    "aicda_rrid": ("Aicda-/- RRID (MGI)", None),
    # primers
    "jh4_primers": ("J_{H}4 primer sequences", None),
    "n114_primers": ("N114 primer sequences", None),
    # bench details (STAR protocol)
    "stain": ("antibody dilutions and fluorophores", None),
    "facs_buffer": ("FACS buffer composition, e.g. PBS with 2% fetal bovine serum", None),
    "lysis": ("lysis buffer, volume, and incubation", None),
    "lysate_volume": ("lysate volume per PCR", None),
    # deposits
    "sra": ("SRA BioProject", None),
    "data_doi": ("data deposit DOI", None),
    "code_doi": ("SHMtrace Zenodo DOI", None),
    "github": ("GitHub URL", None),
    "crm_citation": ("companion Cell Reports Methods article, to be cited once accepted", None),
}


def placeholders():
    """{lab_<key>: value or "[description]"} for crm_numbers.compute()."""
    return {f"lab_{k}": v if v is not None else f"[{d}]" for k, (d, v) in DETAILS.items()}


def missing():
    return [d for d, v in DETAILS.values() if v is None]


if __name__ == "__main__":
    left = missing()
    print(f"{len(DETAILS) - len(left)} of {len(DETAILS)} details filled; still missing:")
    for d in left:
        print(" -", d)
