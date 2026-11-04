"""Codons, amino acids and named regions for amplicons that code, and R/S counts by region.

SHMtrace knows nothing of reading frames by itself: a J_H4 intron has none, and a V region's
frame and CDR boundaries depend on the clone. When the user gives the coding stretch of the
reference (its first base is the first base of a codon) and, optionally, named regions such
as FR1, CDR1 …, each called substitution is placed in its codon and classed as silent,
replacement or stop, and the calls are counted per region, the replacement-to-silent ratio
of the classical CDR/FR comparison.

Each call is judged against the reference codon on its own, as Geneious's Find Variations
does; two mutations in one codon of one read are not combined. Calls outside the coding
stretch are "non-coding" and take no part in the R/S counts.
"""
from __future__ import annotations

import re

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def parse(annotation: dict | None, ref_len: int) -> dict | None:
    """Check and normalise {"cds": "start-end" or [s, e], "regions": [...] or "NAME s-e" lines}."""
    if not annotation:
        return None
    out = {"cds": None, "regions": []}
    cds = annotation.get("cds")
    if cds:
        s, e = _span(cds, "the coding region")
        if not (1 <= s < e <= ref_len):
            raise ValueError(f"the coding region must lie within 1-{ref_len}")
        out["cds"] = [s, e]
    regions = annotation.get("regions") or []
    if isinstance(regions, str):
        regions = [line for line in regions.splitlines() if line.strip()]
    for item in regions:
        if isinstance(item, str):
            m = re.match(r"^\s*(\S+)\s+(\d+)\s*[-–]\s*(\d+)\s*$", item)
            if not m:
                raise ValueError(f"a region is a name and a span, like 'CDR1 76-99': {item!r}")
            item = {"name": m[1], "start": int(m[2]), "end": int(m[3])}
        s, e = int(item["start"]), int(item["end"])
        if not (1 <= s <= e <= ref_len):
            raise ValueError(f"region {item['name']} must lie within 1-{ref_len}")
        out["regions"].append({"name": str(item["name"])[:30], "start": s, "end": e})
    out["regions"].sort(key=lambda r: r["start"])
    return out if out["cds"] or out["regions"] else None


def _span(v, what) -> tuple[int, int]:
    if isinstance(v, (list, tuple)) and len(v) == 2:
        return int(v[0]), int(v[1])
    m = re.match(r"^\s*(\d+)\s*[-–]\s*(\d+)\s*$", str(v))
    if not m:
        raise ValueError(f"give {what} as start-end, for example 1-366")
    return int(m[1]), int(m[2])


def region_of(regions: list[dict], pos: int) -> str:
    return next((r["name"] for r in regions if r["start"] <= pos <= r["end"]), "")


def classify(ref: str, cds: list[int] | None, pos: int, alt: str) -> dict:
    """The codon a substitution falls in and what it does to the amino acid."""
    if not cds or not (cds[0] <= pos <= cds[1]):
        return {"effect": "non-coding"}
    k = (pos - cds[0]) // 3
    start = cds[0] + 3 * k                      # 1-based first base of the codon
    if start + 2 > cds[1]:
        return {"effect": "non-coding"}         # a partial codon at the end
    codon = ref[start - 1:start + 2].upper()
    i = pos - start
    mutant = codon[:i] + alt.upper() + codon[i + 1:]
    a, b = CODON.get(codon, "X"), CODON.get(mutant, "X")
    effect = ("silent" if a == b else "stop" if b == "*" else "replacement")
    return {"codon_number": k + 1, "codon_ref": codon, "codon_alt": mutant,
            "aa_ref": a, "aa_alt": b, "effect": effect}


def apply(res: dict) -> None:
    """Fill res["coding"] (per call) and res["rs_by_region"] from res["annotation"]."""
    ann = res.get("annotation")
    if not ann:
        res.pop("coding", None)
        res.pop("rs_by_region", None)
        return
    ref, cds, regions = res.get("reference") or "", ann.get("cds"), ann.get("regions") or []
    coding, counts = {}, {}
    gmap = {s["sample"]: s.get("group", "") for s in res.get("samples") or []}
    for v in res.get("variants") or []:
        c = classify(ref, cds, v["pos"], v["alt"])
        c["region"] = region_of(regions, v["pos"])
        coding[f"{v['sample']}|{v['pos']}|{v['ref']}>{v['alt']}"] = c
        for reg in ({c["region"], "all"} if c["region"] else {"all"}):
            row = counts.setdefault((v["sample"], reg), {"silent": 0, "replacement": 0, "stop": 0, "non-coding": 0})
            row[c["effect"]] += 1
    rows = []
    for (sample, reg), n in sorted(counts.items(), key=lambda x: (x[0][0], x[0][1] != "all", x[0][1])):
        rows.append({"sample": sample, "group": gmap.get(sample, ""), "region": reg,
                     "replacement": n["replacement"], "silent": n["silent"], "stop": n["stop"],
                     "non_coding": n["non-coding"],
                     "rs_ratio": round(n["replacement"] / n["silent"], 3) if n["silent"] else None})
    res["coding"], res["rs_by_region"] = coding, rows
