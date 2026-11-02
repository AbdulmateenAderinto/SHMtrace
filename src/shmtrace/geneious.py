"""Called mutations as Geneious-style annotation tables, one CSV per sample.

Geneious's "Find Variations/SNPs" step is where most labs met per-site calls, and its
exported annotation table is what their spreadsheets, scripts and lab notebooks read.
This writes SHMtrace's calls in that shape, so a sample can be dropped into the same
downstream steps. The first columns are Geneious's own (Name, Minimum, Maximum, Length,
Change, Polymorphism Type, the nucleotides, Coverage, Variant Frequency, Variant Raw
Frequency, Variant P-Value (approximate)); SHMtrace's evidence and the review mark
follow, under their own names.

Positions are 1-based on the analysis reference, as in variants.csv. Only
substitutions are written; insertions and deletions have their own table.
benchmark.read_geneious_export reads these files back.
"""
from __future__ import annotations

import csv
import io
import re
import zipfile

COLUMNS = ["Name", "Minimum", "Maximum", "Length", "Change", "Polymorphism Type",
           "Reference Nucleotide(s)", "Variant Nucleotide(s)", "Coverage",
           "Variant Frequency", "Variant Raw Frequency", "Variant P-Value (approximate)",
           "Forward Reads", "Reverse Reads", "Context", "Homopolymer Run", "WRC/GYW Hotspot",
           "Sample", "Group", "Review"]
TRANSITIONS = {"AG", "GA", "CT", "TC"}


def review_key(row: dict) -> str:
    """How a call is named in a review: sample, position and change."""
    return f"{row['sample']}|{row['pos']}|{row['ref']}>{row['alt']}"


def _p(v) -> str:
    if v is None or v == "":
        return ""
    v = float(v)
    return "0" if v == 0 else f"{v:.3g}"


def rows(res: dict, sample: str, review: dict | None = None, drop_rejected: bool = False) -> list[dict]:
    """One sample's substitution calls as Geneious annotation rows, by position."""
    review = review or {}
    out = []
    for v in sorted((x for x in res.get("variants") or [] if x["sample"] == sample),
                    key=lambda x: (x["pos"], x["alt"])):
        mark = (review.get(review_key(v)) or {}).get("mark", "")
        if drop_rejected and mark == "rejected":
            continue
        kind = "transition" if v["ref"] + v["alt"] in TRANSITIONS else "transversion"
        out.append({
            "Name": f"{v['ref']}->{v['alt']}", "Minimum": v["pos"], "Maximum": v["pos"], "Length": 1,
            "Change": f"{v['ref']} -> {v['alt']}", "Polymorphism Type": f"SNP ({kind})",
            "Reference Nucleotide(s)": v["ref"], "Variant Nucleotide(s)": v["alt"],
            "Coverage": v["depth"], "Variant Frequency": f"{100 * (v.get('frequency') or 0):.2f}%",
            "Variant Raw Frequency": v["reads"], "Variant P-Value (approximate)": _p(v.get("p_value")),
            "Forward Reads": v.get("plus_reads", ""), "Reverse Reads": v.get("minus_reads", ""),
            "Context": v.get("context", ""), "Homopolymer Run": v.get("homopolymer", ""),
            "WRC/GYW Hotspot": "yes" if v.get("WRC_GYW") else "",
            "Sample": sample, "Group": v.get("group", ""), "Review": mark,
        })
    return out


def csv_text(rows_: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS)
    w.writeheader()
    w.writerows(rows_)
    return buf.getvalue()


def file_name(sample: str) -> str:
    """'<sample> SHMtrace Annotations.csv', with characters a file system refuses replaced."""
    return re.sub(r'[\\/:*?"<>|]', "_", sample) + " SHMtrace Annotations.csv"


def files(res: dict, drop_rejected: bool = False) -> dict[str, str]:
    """Every sample's table, keyed by file name, in the samples' order."""
    review = res.get("review") or {}
    return {file_name(s["sample"]): csv_text(rows(res, s["sample"], review, drop_rejected))
            for s in res.get("samples") or []}


def zip_bytes(res: dict, drop_rejected: bool = False) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in files(res, drop_rejected).items():
            z.writestr(name, text)
    return buf.getvalue()
