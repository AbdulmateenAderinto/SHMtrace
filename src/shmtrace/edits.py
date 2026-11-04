"""Hand edits to Sanger base calls, recorded and applied by re-running.

Geneious lets a person correct a miscalled base in a trace. Here an edit is recorded against
a clone at a reference position (with the base the trace called there, the base it should
be, a note and when), kept with the project, and applied by running the project again: each
edited clone is written as a FASTQ copy of its trace's base calls with the edits made, under
the clone's own name, and the new run reads that copy in place of the .ab1. Qualities are
left as the trace gave them. The new run lists every edit in its notes and keeps them, so
the result always says which calls rest on a corrected base.
"""
from __future__ import annotations

import time
from pathlib import Path

from .chromatogram import mapping
from .io_utils import _list_item, read_abif

BASES = set("ACGTN")
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def locate(trace_file: str, reference: str, pos: int) -> tuple[int, bool, str]:
    """(index in the trace's base calls, read reverse?, the called base along the reference)."""
    _name, seq, _qual = read_abif(trace_file)
    r2i, rev = mapping(seq, reference)
    if pos not in r2i:
        raise ValueError(f"the clone does not cover position {pos}")
    b = seq[r2i[pos]]
    return r2i[pos], rev, (_COMP.get(b, b) if rev else b)


def record(edits: list[dict], sample: str, clone: str, trace_file: str, reference: str,
           ref_seq: str, pos: int, to: str, note: str = "") -> list[dict]:
    """Add, change or (to="") remove the edit at (sample, clone, pos); returns the new list."""
    keep = [e for e in edits if (e["sample"], e["clone"], e["pos"]) != (sample, clone, pos)]
    if not to:
        return keep
    to = to.upper()
    if to not in BASES:
        raise ValueError("a base is one of A, C, G, T or N")
    _i, _rev, called = locate(trace_file, reference, pos)
    if called == to:
        raise ValueError(f"the trace already calls {to} at {pos}")
    return keep + [{"sample": sample, "clone": clone, "pos": pos, "reference": ref_seq[pos - 1].upper(),
                    "called": called, "to": to, "note": note[:300], "when": time.strftime("%Y-%m-%dT%H:%M:%S")}]


def write_edited(trace_file: str, clone: str, reference: str, edits: list[dict], out: Path) -> Path:
    """The clone's base calls with its edits made, as one FASTQ record named after the clone."""
    _name, seq, qual = read_abif(trace_file)
    bases = list(seq)
    r2i, rev = mapping(seq, reference)
    for e in edits:
        i = r2i.get(e["pos"])
        if i is None:
            raise ValueError(f"{clone} no longer covers position {e['pos']}")
        bases[i] = _COMP[e["to"]] if rev else e["to"]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(f"@{clone}\n{''.join(bases)}\n+\n{qual}\n")
    return out


def edited_project(proj: dict, edits: list[dict], folder: Path) -> dict:
    """A copy of the project whose edited clones read their edited copies."""
    by_clone: dict = {}
    for e in edits:
        by_clone.setdefault((e["sample"], e["clone"]), []).append(e)
    reference = proj["reference"]
    # always from the samples as first given, so a re-run of a re-run edits the traces again
    original = proj.get("original_samples") or proj["samples"]
    samples = []
    for s in original:
        items = s.get("path")
        if isinstance(items, (list, tuple)) and any(k[0] == s["name"] for k in by_clone):
            new = []
            for item in items:
                f, rid = _list_item(item)
                es = by_clone.get((s["name"], rid)) or by_clone.get((s["name"], Path(f).stem))
                if es:
                    fq = write_edited(f, rid or Path(f).stem, reference, es,
                                      folder / s["name"] / f"{rid or Path(f).stem}.edited.fastq")
                    new.append([str(fq), rid or Path(f).stem])
                else:
                    new.append(item)
            s = {**s, "path": new}
        samples.append(s)
    return {**proj, "samples": samples, "base_edits": edits, "original_samples": original}


def notes(edits: list[dict]) -> list[str]:
    if not edits:
        return []
    listed = "; ".join(f"{e['sample']} {e['clone']} {e['pos']} {e['called']}→{e['to']}"
                       + (f" ({e['note']})" if e["note"] else "") for e in edits)
    return [f"{len(edits)} Sanger base call{'s' if len(edits) != 1 else ''} edited by hand before this run: {listed}."]
