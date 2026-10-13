"""Inline barcodes: split a pooled run into samples, and report a screen.

Many conditions in one run (a CRISPR screen of SHM regulators, a panel of cell
lines) are often pooled with an inline barcode in the PCR primer: the read carries
<barcode><primer>...<primer'><barcode'>. Native ONT barcodes are demultiplexed by
the basecaller; inline ones are split here.

Each read is oriented by finding the left primer anchor near its start (in the read
or its reverse complement), within `anchor_edits`. The `length` bases just outside
the anchor are the barcode, and are matched to the sample sheet's barcodes: the
best match must be within `max_edits` and at least `margin` edits better than the
second best, or the read is left unassigned rather than risk putting it in the wrong
sample. With dual barcodes (`barcode` and `barcode2` columns), both ends must agree
on the sample.

A sample sheet with a `barcode` column (and the same pooled FASTQ in `file` for
every row) is demultiplexed into one file per sample before the analysis
(samplesheet.Sheet.project); `screen_table` then ranks every condition against the
control with the animal-level tests (replicates.py) and Benjamini-Hochberg q-values.
"""
from __future__ import annotations

import os
from collections import Counter

from .io_utils import read_fastq
from .umi import _edit, _semiglobal

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _flank(seq: str, anchor: str, length: int, side: str, offset: int, window: int,
           anchor_edits: int) -> str | None:
    """The `length` bases outside `anchor` at one end of `seq`, or None."""
    if side == "left":
        region = seq[:window + len(anchor)]
        ed, s, _e = _semiglobal(region, anchor)
        if ed > anchor_edits:
            return None
        end = s - offset
        return seq[max(0, end - length):end] if end - length >= -2 else None
    region_start = max(0, len(seq) - window - len(anchor))
    ed, _s, e = _semiglobal(seq[region_start:], anchor)
    if ed > anchor_edits:
        return None
    start = region_start + e + offset
    out = seq[start:start + length]
    return out if len(out) >= length - 2 else None


def _best(code: str, barcodes: dict[str, str], max_edits: int, margin: int) -> str | None:
    if code is None:
        return None
    scored = sorted((_edit(code, bc, max_edits + margin), name) for name, bc in barcodes.items())
    if not scored or scored[0][0] > max_edits:
        return None
    if len(scored) > 1 and scored[1][0] - scored[0][0] < margin:
        return None
    return scored[0][1]


def assign(seq: str, left_anchor: str, right_anchor: str, bc1: dict, bc2: dict | None,
           length: int, offset: int = 0, window: int = 120, anchor_edits: int = 4,
           max_edits: int = 2, margin: int = 2) -> str | None:
    """The sample a read belongs to, or None (no anchor, no clear barcode, or ends disagree)."""
    for s in (seq, _rc(seq)):
        left = _flank(s, left_anchor, length, "left", offset, window, anchor_edits)
        if left is None:
            continue
        a = _best(left, bc1, max_edits, margin)
        if bc2:
            right = _flank(s, right_anchor, length, "right", offset, window, anchor_edits)
            b = _best(right, bc2, max_edits, margin)
            return a if (a and a == b) else None       # both ends must name the sample
        return a
    return None


def demultiplex(fastq: str, samples: dict[str, tuple[str, str | None]], left_anchor: str,
                right_anchor: str, out_dir: str, length: int | None = None, **kw) -> dict:
    """Split `fastq` into one FASTQ per sample. samples: name -> (barcode, barcode2 or None).

    Returns {"paths": {sample: path}, "stats": {...}}. Anchors are given
    reference-forward: the left primer as it reads at the reference start, the right
    primer as it reads at the reference end (so the right barcode follows it).
    """
    os.makedirs(out_dir, exist_ok=True)
    bc1 = {n: b[0].upper() for n, b in samples.items()}
    bc2 = {n: b[1].upper() for n, b in samples.items() if b[1]} or None
    length = length or len(next(iter(bc1.values())))
    handles = {n: open(os.path.join(out_dir, f"{n}.fastq"), "w") for n in samples}
    counts = Counter()
    try:
        for name, seq, qual in read_fastq(fastq):
            who = assign(seq.upper(), left_anchor.upper(), right_anchor.upper(), bc1, bc2,
                         length, **kw)
            counts[who or "_unassigned"] += 1
            if who:
                handles[who].write(f"@{name}\n{seq}\n+\n{qual}\n")
    finally:
        for h in handles.values():
            h.close()
    total = sum(counts.values())
    return {"paths": {n: os.path.join(out_dir, f"{n}.fastq") for n in samples},
            "stats": {"reads": total, "unassigned": counts["_unassigned"],
                      "per_sample": {n: counts[n] for n in samples}}}


def bh(pvals: list[float | None]) -> list[float | None]:
    """Benjamini-Hochberg q-values (None stays None)."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    q: list[float | None] = [None] * len(pvals)
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


def screen_table(res: dict, metric: str = "unique_density", test: str = "p_negbinomial") -> list[dict]:
    """Every condition against the control, ranked, with BH q-values across conditions."""
    rows = [r for r in (res.get("animal_tests") or []) if r["metric"] == metric]
    qs = bh([r.get(test) for r in rows])
    out = []
    for r, q in zip(rows, qs, strict=True):
        out.append({"condition": r["group"], "control": r["control"],
                    "animals": r["animals_case"], "fold": r.get("fold"),
                    "fold_is_lower_bound": r.get("fold_is_lower_bound", False),
                    "fold_ci": r.get("fold_ci"), "p": r.get(test), "q": q})
    out.sort(key=lambda r: (r["q"] is None, r["q"] if r["q"] is not None else 1, -(r["fold"] or 0)))
    return out


def screen_plot(rows: list[dict], path: str, q_cut: float = 0.1) -> None:
    """Volcano: log2 fold against -log10 q, conditions past `q_cut` labelled."""
    import math

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs = [math.log2(r["fold"]) if r["fold"] else 0 for r in rows]
    ys = [-math.log10(max(r["q"], 1e-6)) if r["q"] is not None else 0 for r in rows]
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    hit = [r["q"] is not None and r["q"] < q_cut for r in rows]
    ax.scatter(xs, ys, c=["#C2548A" if h else "#999" for h in hit], s=28, zorder=3)
    for x, y, r, h in zip(xs, ys, rows, hit, strict=True):
        if h:
            ax.annotate(r["condition"], (x, y), xytext=(4, 3), textcoords="offset points", fontsize=8)
    ax.axhline(-math.log10(q_cut), ls="--", color="#bbb", lw=1)
    ax.axvline(0, color="#ddd", lw=1)
    ax.set_xlabel("log2 fold vs control (animal-level)")
    ax.set_ylabel("-log10 q (Benjamini-Hochberg)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main(argv=None):
    import argparse
    import csv
    import json
    ap = argparse.ArgumentParser(prog="shmtrace-screen",
                                 description="Rank every condition of a project against its "
                                             "control (animal-level tests, BH q-values)")
    ap.add_argument("result", help="project.json")
    ap.add_argument("--metric", default="unique_density",
                    choices=("unique_density", "mutation_frequency"))
    ap.add_argument("--q", type=float, default=0.1, help="q-value to label on the plot")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    rows = screen_table(json.load(open(a.result)), a.metric)
    with open(os.path.join(a.out, "screen.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["condition"])
        w.writeheader()
        w.writerows(rows)
    screen_plot(rows, os.path.join(a.out, "screen.png"), a.q)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
