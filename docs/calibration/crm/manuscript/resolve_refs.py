"""Resolve the manuscripts' references against PubMed (title search) and format them in Cell style.

    python docs/calibration/crm/manuscript/resolve_refs.py
Writes refs.json: key -> {pmid, cell (formatted string), title, ...}. Paced at < 3 requests/s.
Nothing is typed by hand: every author list, journal, volume, pages and DOI comes from PubMed.
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERIES = {   # key: PubMed query (a title, optionally with an author)
    "muramatsu2000": "Class switch recombination and hypermutation require activation-induced cytidine deaminase[Title]",
    "dinoia2007": "Molecular mechanisms of antibody somatic hypermutation[Title] AND Di Noia[Author]",
    "methot2017": "Molecular Mechanisms of Somatic Hypermutation and Class Switch Recombination[Title] AND Methot[Author]",
    "jolly1997": "Rapid methods for the analysis of immunoglobulin gene hypermutation[Title]",
    "lebecque1990": "Boundaries of somatic mutation in rearranged immunoglobulin genes[Title]",
    "sible2021": "somatic hypermutation in the JH4 intron of germinal center B cells from mouse Peyer's patches[Title]",
    "dinoia2002": "Altering the pathway of immunoglobulin hypermutation by inhibiting uracil-DNA glycosylase[Title]",
    "rada2004": "Mismatch recognition and uracil excision provide complementary paths[Title]",
    "rogozin1992": "Somatic hypermutagenesis in immunoglobulin genes[Title] AND Rogozin[Author]",
    "wilm2012": "LoFreq: a sequence-quality aware, ultra-sensitive variant caller[Title]",
    "li2018": "Minimap2: pairwise alignment for nucleotide sequences[Title]",
    "edgar2015": "Error filtering, pair assembly and error correction for next-generation sequencing reads[Title]",
    "yaari2013": "Models of somatic hypermutation targeting and substitution based on synonymous mutations[Title]",
    "gupta2015": "Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data[Title]",
    "maccarthy2009": "SHMTool: a webserver for comparative analysis of somatic hypermutation datasets[Title]",
    "kearse2012": "Geneious Basic: an integrated and extendable desktop software platform[Title]",
    "karst2021": "High-accuracy long-read amplicon sequences using unique molecular identifiers[Title]",
    "lazic2010": "The problem of pseudoreplication in neuroscientific studies[Title]",
    "love2014": "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2[Title]",
    "smith2017": "UMI-tools: modeling sequencing errors in Unique Molecular Identifiers[Title]",
    "zhang2023": "The newest Oxford Nanopore R10.4.1 full-length 16S rRNA sequencing[Title]",
    "sereika2022": "Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes[Title]",
    "sanderson2024": "Evaluation of the accuracy of bacterial genome reconstruction with Oxford Nanopore R10.4.1[Title]",
    "martin2002": "Activation-induced cytidine deaminase turns on somatic hypermutation in hybridomas[Title]",
    "martin2002b": "Martin A[Author] AND Scharff MD[Author] AND 2002[dp] AND Proc Natl Acad Sci U S A[Journal] AND transgene",
    "barbulescu2025": "MKLN1-dependent GID4/CTLH E3 ubiquitin ligase complex assemblies[Title]",
    "barbulescu2024": "FAM72A degrades UNG2 through the GID/CTLH complex[Title]",
    "feng2021": "FAM72A antagonizes UNG2 to promote mutagenic repair during antibody maturation[Title]",
    "sharma2023": "Distinct metabolic requirements regulate B cell activation and germinal center responses[Title]",
    "wang2023": "Mesoscale DNA feature in antibody-coding sequence facilitates somatic hypermutation[Title]",
    "wu2022": "HMCES protects immunoglobulin genes specifically from deletions during somatic hypermutation[Title]",
    "yaari2013sim": "Quantifying selection in high-throughput Immunoglobulin sequencing data sets[Title]",
    "minervini2021": "Nanopore sequencing approach for immunoglobulin gene analysis in chronic lymphocytic leukemia[Title]",
    "subassatish2022": "NAb-seq: an accurate, rapid, and cost-effective method for antibody long-read sequencing[Title]",
    "deng2024": "R2C2 + UMI: Combining concatemeric and unique molecular identifier-based consensus sequencing[Title]",
    "martin2023": "The IgH Eµ-MAR regions promote UNG-dependent error-prone repair[Title]",
    "delahaye2021": "Sequencing DNA with nanopores: Troubles and biases[Title]",
    "zeng2001": "DNA polymerase eta is an A-T mutator in somatic hypermutation of immunoglobulin variable genes[Title]",
    "rada2002": "Immunoglobulin isotype switching is inhibited and somatic hypermutation perturbed in UNG-deficient mice[Title]",
    "muramatsu2000": "11007474[uid]", "rogozin1992": "1420357[uid]", "karst2021": "33432244[uid]",
    "sanderson2024": "38713194[uid]", "virtanen2020": "32015543[uid]",
    "martin2002b": "Martin A[Author] AND Scharff MD[Author] AND 2002[dp] AND Proc Natl Acad Sci U S A[Journal] AND transgene",
    "zeng2001": "Zeng X[Author] AND Gearhart PJ[Author] AND 2001[dp] AND Nat Immunol[Journal]",
    "danecek2021": "Twelve years of SAMtools and BCFtools[Title]",
}


def get(url):
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read())
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(url)


def authors_cell(names):
    """Cell style: 'Last, F.M., Last, F., and Last, F.' ; 'et al.' only after 10 names."""
    fmt = []
    for n in names:
        parts = n.split(" ")
        last, ini = " ".join(parts[:-1]), parts[-1]
        fmt.append(f"{last}, {'.'.join(ini)}.")
    if len(fmt) > 10:
        return ", ".join(fmt[:10]) + ", et al."
    if len(fmt) == 1:
        return fmt[0]
    return ", ".join(fmt[:-1]) + ", and " + fmt[-1]


def cell_journal(abbrev, full):
    """PubMed's abbreviation with Cell's periods: a word is abbreviated (gets a period) unless it
    appears whole in the full journal name (Cell, Nature, Bioinformatics)."""
    words = {w.lower().strip(".,:") for w in full.split()}
    return " ".join(w if w.lower() in words or not w.isalpha() else w + "." for w in abbrev.split())


def cell_pages(p):
    """'553-63' -> '553–563' (full end page, en dash)."""
    if "-" not in p:
        return p
    a, b = p.split("-", 1)
    if a.isdigit() and b.isdigit() and len(b) < len(a):
        b = a[: len(a) - len(b)] + b
    elif a[:1].isalpha() and b.isdigit():          # e2193-e2207 style left alone
        pass
    return f"{a}–{b}"


def main():
    out = json.load(open(HERE / "refs.json")) if (HERE / "refs.json").exists() else {}
    for key, q in QUERIES.items():
        if key in out:
            continue
        ids = get(E + "esearch.fcgi?db=pubmed&retmode=json&term=" + urllib.parse.quote(q))["esearchresult"]["idlist"]
        time.sleep(0.45)
        if len(ids) != 1:
            print("CHECK", key, ids)
            out[key] = {"query": q, "ids": ids}
            continue
        s = get(E + f"esummary.fcgi?db=pubmed&retmode=json&id={ids[0]}")["result"][ids[0]]
        time.sleep(0.45)
        names = [a["name"] for a in s["authors"] if a.get("authtype") == "Author"]
        doi = next((x["value"] for x in s.get("articleids", []) if x["idtype"] == "doi"), "")
        year = s["pubdate"][:4]
        title = s["title"].rstrip(".")
        journal = cell_journal(s["source"], s.get("fulljournalname", ""))
        cell = (f"{authors_cell(names)} ({year}). {title}. {journal} {s.get('volume', '')}, "
                f"{cell_pages(s.get('pages', ''))}. {doi}.").replace(" ,", ",")
        out[key] = {"pmid": ids[0], "cell": cell, "title": title, "year": year, "first": names[0] if names else "",
                    "authors": names, "journal": journal, "volume": s.get("volume", ""),
                    "pages": cell_pages(s.get("pages", "")), "doi": doi}
        print(key, "->", cell[:110])
        json.dump(out, open(HERE / "refs.json", "w"), indent=1, ensure_ascii=False)
    json.dump(out, open(HERE / "refs.json", "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
