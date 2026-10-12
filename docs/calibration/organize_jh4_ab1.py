"""Organise the J꜀4 Sanger traces by mouse, correcting the labels the sequence contradicts.

    python docs/calibration/organize_jh4_ab1.py SOURCE_DIR OUT_DIR [GENBANK_DIR]

SOURCE_DIR holds the sequencing provider's order folders (berrumar_<order>_chromats/*.ab1).
Every trace is copied, never moved, to OUT_DIR/<mouse>/<mouse>_<order>-<well>_<label>.ab1,
and OUT_DIR/manifest.csv records where each came from and why it went where it did.

The mouse is the one the file's label names (w2/wt2 -> WT2, a1 -> Aicda_A1, k3 -> MKLN_K3,
...), except where the sequence contradicts it (CORRECTIONS, evidence in
docs/calibration/SANGER_VS_NANOPORE.md). Traces that are not usable J꜀4 go to
OUT_DIR/_excluded/<reason>/, decided by nanoshm's own read preparation: the J꜀4 paralogue
or other contaminant, no alignment, or too short an alignment.

GENBANK_DIR (optional) holds the hand-collated GenBank clone exports. A record whose name
matches no trace is a clone sequenced in an order whose traces are not in SOURCE_DIR; it is
written to OUT_DIR/<mouse>/genbank_only_<mouse>.fasta (no qualities: read at a flat Q30) so a
mouse's folder holds all of its clones. nanoshm reads a folder of .ab1 and .fasta files as
one sample.
"""
from __future__ import annotations

import csv
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

from nanoshm import pipeline, project
from nanoshm.io_utils import parse_abif, read_abif, read_fasta

ROOT = Path(__file__).resolve().parents[2]
REF = str(ROOT / "src/nanoshm/references/mouse_jh4_intron.fasta")
CON = [str(ROOT / "src/nanoshm/references/mouse_jh4_contaminant.fasta")]

# GenBank records with no trace whose label is contradicted by the sequence: the K4
# traces 88-k4-28 and 96-k4-36 renamed as WT2 clones (identical sequences).
GENBANK_EXCLUDE = {"88-w2-28-user_added": "renamed copy of trace 434355-H11 88-k4-28 (MKLN K4)",
                   "96-w2-36-user_added": "renamed copy of trace 434355-H12 96-k4-36 (MKLN K4)"}

# (order, well) -> (mouse, reason). Each is a contradiction between label and sequence.
CORRECTIONS = {
    ("434141", "E8"): ("WT3", "labelled w2-1, but 18 mutations where WT2's other 34 traces carry "
                              "at most 3; sits where w3-1 belongs in the w3 series (G8 w3-3, H8 "
                              "w3-4, A9 w3-5) while w2-1 is already at A7 on this plate"),
    ("434141", "F8"): ("WT3", "labelled w2-2, but 13 mutations (WT2 traces: at most 3); sits where "
                              "w3-2 belongs, while w2-2 is already at B7 on this plate"),
    ("434141", "A11"): ("MKLN_K2", "labelled w4-9, but carries the same 7-mutation lineage as K2 "
                                   "clones k2-1 (E11, same plate) and k2-24 (order 434355); not in "
                                   "any WT4 nanopore read"),
    ("432625", "E6"): ("_excluded/unassigned", "labelled a1-1 (Aicda-/-), but 37 AID-type mutations "
                                               "(64% of G:C changes at WRC/GYW); an AID-null mouse "
                                               "cannot make this, and it shares no lineage with any "
                                               "clone, so its mouse is unknown"),
}


def mouse_from_label(label: str) -> str:
    s = label.lower()
    for pat, m in ((r"jh6", "_excluded/JH6_amplicon"), (r"(mkln_ko_)(\d)", "MKLN_K{}"),
                   (r"(jh4k)(\d)", "MKLN_K{}"), (r"(^|-)k(\d)", "MKLN_K{}"),
                   (r"(^|-)a(\d)", "Aicda_A{}"), (r"(^|-)wt?2-", "WT2"), (r"(^|-)wt?3-", "WT3"),
                   (r"(^|-)w4", "WT4"), (r"w1", "_excluded/unassigned")):
        hit = re.search(pat, s)
        if hit:
            return m.format(hit.group(2)) if "{}" in m else m
    return "_excluded/unassigned"


def well_of(path: str) -> tuple[str, str]:
    t = parse_abif(open(path, "rb").read())

    def text(tag):
        v = t.get((tag, 1))
        return v[1:1 + v[0]].decode(errors="replace") if v else ""
    import struct
    d = t.get(("RUND", 1))
    date = "{:04d}-{:02d}-{:02d}".format(*struct.unpack(">hbb", d[:4])) if d else ""
    return text("TUBE"), date


def main(src: str, out: str, genbank: str | None = None) -> None:
    src_p, out_p = Path(src), Path(out)
    traces = []
    for f in sorted(src_p.glob("*/*.ab1")):
        order = f.parent.name.replace("berrumar_", "").replace("_chromats", "")
        well, date = well_of(str(f))
        label = re.sub(r"-user_added$", "", f.stem)
        traces.append({"path": f, "order": order, "well": well, "date": date, "label": label})

    # nanoshm's own read preparation decides which traces are usable J꜀4
    p = project.build_params({"region": "jh4_intron", "platform": "sanger"})
    ref = read_fasta(REF)[0][1].upper()
    with tempfile.TemporaryDirectory() as tmp:
        fq = Path(tmp) / "all.fastq"
        with open(fq, "w") as fh:
            for t in traces:
                _n, s, q = read_abif(str(t["path"]))
                fh.write(f"@{t['order']}|{t['well']}\n{s}\n+\n{q}\n")
        fates: dict = {}
        aligned, _ = pipeline.prepare_reads(str(fq), REF, CON, p, fates=fates)
    germ = project.sanger_germline({"all": aligned}, ref, 100, 628, p)
    muts = {a.qname: sum(1 for pos, b, q in a.pairs if 100 <= pos <= 628 and q >= 20 and b in "ACGT"
                         and b != ref[pos] and (pos, b) not in germ) for a in aligned}
    fate_dir = {"contaminant": "_excluded/paralogue_or_contaminant", "unmapped": "_excluded/failed_read",
                "alignment too short": "_excluded/failed_read", "chimeric": "_excluded/failed_read"}

    if out_p.exists():
        shutil.rmtree(out_p)
    rows = []
    for t in traces:
        key = f"{t['order']}|{t['well']}"
        by_label = mouse_from_label(t["label"])
        mouse, status, reason = by_label, "kept", ""
        if (t["order"], t["well"]) in CORRECTIONS:
            mouse, reason = CORRECTIONS[(t["order"], t["well"])]
            status = "excluded" if mouse.startswith("_excluded") else "relabelled"
        elif key in fates:
            f = fates[key]["fate"]
            mouse, status, reason = fate_dir.get(f, "_excluded/failed_read"), "excluded", f"nanoshm: {f}"
        elif by_label.startswith("_excluded"):
            status, reason = "excluded", "label names no mouse in this study"
        tag = by_label.split("/")[-1]
        name = f"{mouse.split('/')[-1] if not mouse.startswith('_excluded') else tag}_{t['order']}-{t['well']}_{t['label']}.ab1"
        dest = out_p / mouse / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(t["path"], dest)
        rows.append({"mouse": mouse, "file": str(dest.relative_to(out_p)), "status": status,
                     "original_label": t["label"], "mouse_by_label": by_label, "order": t["order"],
                     "well": t["well"], "run_date": t["date"],
                     "mutations_q20_101_629": muts.get(key, ""), "reason": reason,
                     "original_path": str(t["path"])})
    if genbank:
        from nanoshm.io_utils import read_fastq
        stems = {t["path"].stem for t in traces}
        recs: dict = {}
        for f in sorted(Path(genbank).glob("*/*.gb")):
            for n, s, _q in read_fastq(str(f)):
                if n not in stems and (n not in recs or len(s) > len(recs[n])):
                    recs[n] = s
        per_mouse: dict = {}
        for n, s in sorted(recs.items()):
            label = re.sub(r"-user_a\w*$", "", n)
            m = mouse_from_label(label)
            reason = GENBANK_EXCLUDE.get(n, "")
            if reason:
                m = "_excluded/unassigned"
            if m.startswith("_excluded") or not re.match(r"^(WT|Aicda)", m):
                status = "excluded"
                reason = reason or ("label names no mouse in this study" if m.startswith("_excluded")
                                    else "MKLN clone with no trace (not needed here)")
                if not m.startswith("_excluded"):
                    m = "_excluded/genbank_only_mkln"
            else:
                status = "genbank_only"
                reason = "no trace in SOURCE_DIR; GenBank sequence read at flat Q30"
                per_mouse.setdefault(m, []).append((label, s))
            rows.append({"mouse": m, "file": (f"{m}/genbank_only_{m}.fasta" if status == "genbank_only"
                                              else ""), "status": status, "original_label": label,
                         "mouse_by_label": mouse_from_label(label), "order": "", "well": "",
                         "run_date": "", "mutations_q20_101_629": "", "reason": reason,
                         "original_path": f"GenBank record {n}"})
        for m, recs_m in per_mouse.items():
            (out_p / m).mkdir(parents=True, exist_ok=True)
            with open(out_p / m / f"genbank_only_{m}.fasta", "w") as fh:
                fh.write("".join(f">{label}\n{s}\n" for label, s in recs_m))
    with open(out_p / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["mouse"], r["order"], r["well"])))
    from collections import Counter
    for (m, s), n in sorted(Counter((r["mouse"], r["status"]) for r in rows).items()):
        print(f"{m:40} {s:10} {n}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
