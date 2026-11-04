"""A Sanger trace around one reference position, drawn along the reference.

Check calls shows the peaks under a call so a person can tell a clean peak from a shoulder
or a mixed one. This reads the four trace channels and the base calls' peak positions from
the .ab1 file, aligns the clone's own base calls to the reference (so the window is found
whatever trimming the analysis applied), and returns the stretch of trace that covers a
few bases either side of the position, in reference orientation: a clone read against the
reference is reversed and its channels complemented, as Geneious shows it.
"""
from __future__ import annotations

import struct
from functools import lru_cache

from .align import _revcomp, build_index
from .io_utils import parse_abif

_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def read_trace(path: str) -> dict:
    """{"channels": {base: [intensity, ...]}, "peaks": [...], "bases": str, "quals": [...]}."""
    with open(path, "rb") as fh:
        tags = parse_abif(fh.read())
    order = (tags.get(("FWO_", 1)) or b"GATC").decode("latin-1")[:4]
    channels = {}
    for i, base in enumerate(order):
        raw = tags.get(("DATA", 9 + i))
        if raw is None:
            raise ValueError(f"{path}: no trace data (DATA {9 + i})")
        channels[base] = list(struct.unpack(f">{len(raw) // 2}h", raw))
    bases_b = tags.get(("PBAS", 2)) or tags.get(("PBAS", 1)) or b""
    peaks_b = tags.get(("PLOC", 2)) if ("PBAS", 2) in tags else tags.get(("PLOC", 1))
    quals_b = tags.get(("PCON", 2)) if ("PBAS", 2) in tags else tags.get(("PCON", 1))
    bases = bases_b.decode("latin-1").upper().replace("\x00", "")
    peaks = list(struct.unpack(f">{len(peaks_b) // 2}h", peaks_b)) if peaks_b else []
    quals = list(quals_b[:len(bases)]) if quals_b else [0] * len(bases)
    n = min(len(bases), len(peaks))
    return {"channels": channels, "peaks": peaks[:n], "bases": bases[:n], "quals": quals[:n]}


@lru_cache(maxsize=8)
def _aligner(reference: str):
    aligner, _t, _c = build_index([reference], [], "map-ont")
    return aligner


def _ref_to_query(hit) -> dict[int, int]:
    """0-based reference position -> 0-based index in the reference-forward query."""
    out, q, r = {}, hit.q_st, hit.r_st
    for length, op in hit.cigar:
        if op in (0, 7, 8):
            for k in range(length):
                out[r + k] = q + k
            q += length
            r += length
        elif op == 1:
            q += length
        elif op in (2, 3):
            r += length
    return out


def mapping(seq: str, reference: str) -> tuple[dict[int, int], bool]:
    """(1-based reference position -> 0-based index in `seq` as read, read reverse?)."""
    hit = max((h for h in _aligner(reference).map(seq) if h.is_primary), key=lambda h: h.mlen, default=None)
    if hit is None:
        raise ValueError("the clone's base calls do not align to the reference")
    rev, n = hit.strand == -1, len(seq)
    # minimap2 reports query coordinates on the read as given; for a reverse hit the
    # reference-forward query is the reverse complement, so remap onto it

    class _H:                                   # the hit, with its start on the forward query
        q_st, r_st, cigar = (n - hit.q_en if rev else hit.q_st), hit.r_st, hit.cigar
    return {r + 1: (n - 1 - q if rev else q) for r, q in _ref_to_query(_H).items()}, rev


def window(path: str, reference: str, pos: int, flank: int = 10) -> dict:
    """The trace over reference positions pos-flank..pos+flank (1-based), in reference orientation."""
    tr = read_trace(path)
    seq = tr["bases"]
    r2i, rev = mapping(seq, reference)
    lo, hi = pos - flank, pos + flank
    cols = []
    for p in range(lo, hi + 1):
        i = r2i.get(p)                          # index into the read's own base calls
        if i is None:
            cols.append({"pos": p, "base": None, "quality": None, "peak": None})
            continue
        b = seq[i]
        cols.append({"pos": p, "base": _COMPLEMENT.get(b, b) if rev else b,
                     "quality": tr["quals"][i], "peak": tr["peaks"][i]})
    peaks = [c["peak"] for c in cols if c["peak"] is not None]
    if not peaks:
        raise ValueError(f"the clone does not cover position {pos}")
    spacing = max(4, (max(peaks) - min(peaks)) // max(1, len(peaks) - 1))
    a, b = max(0, min(peaks) - spacing), max(peaks) + spacing
    out_ch = {}
    for base, ys in tr["channels"].items():
        seg = ys[a:b + 1]
        out_ch[_COMPLEMENT.get(base, base) if rev else base] = seg[::-1] if rev else seg
    for c in cols:                              # peak positions within the returned stretch
        if c["peak"] is not None:
            c["peak"] = (b - c["peak"]) if rev else (c["peak"] - a)
    return {"pos": pos, "start": lo, "end": hi, "reverse": rev, "channels": out_ch,
            "columns": cols, "length": b - a + 1}
