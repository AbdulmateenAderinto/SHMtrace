"""Read alignments for visual inspection, before and after processing.

The pipeline decides which reads to keep by rules that are each defensible on
their own, and the numbers it reports stand or fall on those decisions. This
module hands back the alignment itself so the decisions can be looked at: every
input read in its aligned place against the reference, labelled with what became
of it, and the called mutations marked.

It is the view the Geneious "assemble, then find variants" step gave, except
that a dropped read is still drawn, in the position it would have occupied, with
the reason it went.

A read is reduced to its differences from the reference, not its sequence: the
span it covers, the positions where it disagrees and what it reads there, and
the quality of those bases. For a 1.2 kb amplicon at a 1% error rate that is a
dozen entries per read instead of 1,200, which is what makes a few thousand
reads cheap enough to send to a browser.
"""
from __future__ import annotations

from .io_utils import read_fasta

# how a read left the pipeline, in the order the stages run. The browser
# colours by this and offers it as a filter.
FATES = [
    ("kept", "analysed"),
    ("not sampled", "dropped at depth normalisation"),
    ("minor haplotype", "removed as a co-occurring block"),
    ("alignment too short", "aligned span below the minimum"),
    ("chimeric", "more than one primary alignment"),
    ("contaminant", "best hit was a contaminant reference"),
    ("off target", "best hit was neither target nor contaminant"),
    ("low quality", "expected-error rate above the gate"),
    ("too long", "longer than the amplicon window allows"),
    ("no quality segment", "nothing left after quality trimming"),
    ("unmapped", "did not align"),
]
FATE_ORDER = {name: i for i, (name, _d) in enumerate(FATES)}


def _diffs(pairs, ref: str, minq: int = 0) -> list:
    """(position, base, quality) wherever a read disagrees with the reference."""
    out = []
    for rpos, base, q in pairs:
        if base != ref[rpos].upper() and base in "ACGT" and q >= minq:
            out.append([rpos + 1, base, int(q)])
    return out


def alignment(aligned, fates: dict, reference: str, variants=None,
              max_reads: int = 4000) -> dict:
    """One sample's reads, kept and dropped, ready to draw.

    `aligned` is what the pipeline kept; `fates` is what it recorded about
    everything else. Reads are ordered by start position so the browser can
    draw them as a pileup, and capped at `max_reads` with the kept ones first,
    because a dropped read is only ever evidence about the sample while a kept
    read is evidence about the result.
    """
    ref_name, ref = read_fasta(reference)[0]
    rows = []
    for a in aligned:
        rows.append({"name": a.qname, "fate": "kept", "start": a.r_st + 1,
                     "end": a.r_en, "strand": getattr(a, "strand", 1),
                     "diffs": _diffs(a.pairs, ref),
                     "indels": [[pos + 1, op, (x if isinstance(x, str) else int(x))]
                                for pos, op, x, _q in getattr(a, "indels", ())]})
    for name, rec in fates.items():
        row = {"name": name, "fate": rec.get("fate", "dropped"),
               "start": (rec.get("r_st", 0) or 0) + 1, "end": rec.get("r_en") or 0,
               "strand": rec.get("strand", 1), "diffs": [], "indels": []}
        if rec.get("pairs"):
            row["diffs"] = _diffs(rec["pairs"], ref)
        rows.append(row)
    # Ordered by fate, then by position within a fate. Grouping this way is
    # what makes the view answer "before and after" at a glance: the analysed
    # reads form one band and each reason for dropping forms its own.
    rows.sort(key=lambda r: (FATE_ORDER.get(r["fate"], 99), r["start"], r["name"]))
    truncated = max(0, len(rows) - max_reads)
    rows = rows[:max_reads]
    counts: dict = {}
    for r in rows:
        counts[r["fate"]] = counts.get(r["fate"], 0) + 1
    called = [{"pos": v.pos, "ref": v.ref, "alt": v.alt, "reads": v.supporting_reads}
              for v in (variants or [])]
    return {"reference_name": ref_name, "reference": ref, "reads": rows,
            "counts": counts, "truncated": truncated, "called": called}


def window(data: dict, start: int = 1, end: int | None = None,
           fates: set | None = None, offset: int = 0, limit: int = 400) -> dict:
    """A slice of an alignment: the reads overlapping a region, paged.

    Kept separate from `alignment` so the browser can scroll a large pileup
    without the server rebuilding it each time.
    """
    end = end or len(data["reference"])
    rows = [r for r in data["reads"]
            if r["end"] >= start and (r["start"] <= end or not r["end"])
            and (fates is None or r["fate"] in fates)]
    total = len(rows)
    page = rows[offset:offset + limit]
    return {
        "reference_name": data["reference_name"],
        "sequence": data["reference"][start - 1:end],
        "start": start, "end": end,
        "reads": page, "total": total, "offset": offset,
        "counts": data["counts"], "truncated": data["truncated"],
        "called": [c for c in data["called"] if start <= c["pos"] <= end],
    }
