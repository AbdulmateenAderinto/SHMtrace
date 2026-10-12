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

FOLDER_MOUSE = {"W2": "WT2", "W3": "WT3", "W4": "WT4", "A1": "Aicda_A1", "A2": "Aicda_A2",
                "A3": "Aicda_A3", "K2": "MKLN_K2", "K3": "MKLN_K3", "K4": "MKLN_K4"}

# (GenBank folder, mouse by label) -> (mouse, reason): whole batches the folder contradicts,
# with evidence at the batch level (docs/calibration/SANGER_VS_NANOPORE.md)
GROUP_MOVES = {
    ("W2", "WT3"): ("WT3", "filed under W2, labelled w3/wt3: five of this batch share 8-13 "
                           "mutations with single WT3 nanopore reads and none with WT2's; its loads "
                           "match WT3's clones (p = 0.28), not WT2's own (p = 1.5e-6)"),
    ("W3", "WT2"): ("WT2", "filed under W3, labelled w2: near-unmutated (11 of 12 carry none), "
                           "which 22 of WT3's 23 clones are not; WT2's own clones are"),
}

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
    ("414309", "F7"): ("_excluded/unassigned", "labelled jh6-a1-e6 (Aicda-/-), but 5 confident AID-type "
                                               "mutations (Q26-56, four at WRC/GYW or WA/TW, mostly A:T "
                                               "transitions); a contaminant, removed at the user's "
                                               "decision"),
    ("432625", "E6"): ("_excluded/unassigned", "labelled a1-1 (Aicda-/-), but 37 AID-type mutations "
                                               "(64% of G:C changes at WRC/GYW); an AID-null mouse "
                                               "cannot make this, and it shares no lineage with any "
                                               "clone, so its mouse is unknown"),
}


def mouse_from_label(label: str) -> str:
    s = label.lower()
    jh6 = re.search(r"jh6-(\w+?)-", s)
    if jh6:        # 'jh6' clones are J꜀4-intron clones from another primer set: their mouse is
        s = jh6.group(1) + "-"           # the tag after it (jh6-a3-g9 is an A3 clone)
    for pat, m in ((r"(mkln_ko_)(\d)", "MKLN_K{}"),
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
    # share of the 529-bp window (101-629) each trace reads at Q>=20. Below MIN_USABLE a
    # trace is excluded: mutations / (sequences x 529) assumes every sequence covers it
    usable = {a.qname: sum(1 for pos, b, q in a.pairs if 100 <= pos <= 628 and q >= 20 and b in "ACGT")
              / 529 for a in aligned}
    muts = {a.qname: sum(1 for pos, b, q in a.pairs if 100 <= pos <= 628 and q >= 20 and b in "ACGT"
                         and b != ref[pos] and (pos, b) not in germ) for a in aligned}
    # only a trace nanoshm cannot read is excluded: no stretch of good quality (a noisy
    # chromatogram), no or too short an alignment, or a chimera. A clean but short trace
    # stays, and counts the bases it reads.
    fate_dir = {"contaminant": "_excluded/paralogue_or_contaminant"}

    # the original organisation: the GenBank folder each record was filed in (the "24"
    # folder is the union of all and says nothing). It is the default; the label is used
    # only for traces with no GenBank record.
    gb_folder: dict = {}
    if genbank:
        for f in Path(genbank).glob("*/*.gb"):
            fo = f.parent.name.split("_-_")[0]
            if fo == "24":
                continue
            name = open(f).readline().split()[1]
            gb_folder.setdefault(name, set()).add(fo)

    def folder_of(stem):
        for n in (stem, stem[:23]):        # Geneious cuts names at 23 characters
            if n in gb_folder and len(gb_folder[n]) == 1:
                return next(iter(gb_folder[n]))
        return None

    if out_p.exists():
        shutil.rmtree(out_p)
    rows = []
    for t in traces:
        key = f"{t['order']}|{t['well']}"
        by_label = mouse_from_label(t["label"])
        fo = folder_of(t["path"].stem)
        default = FOLDER_MOUSE.get(fo, by_label) if fo else by_label
        mouse, status, reason = default, "kept", (f"GenBank folder {fo}" if fo else "label (no GenBank record)")
        group = GROUP_MOVES.get((fo, by_label))
        if (t["order"], t["well"]) in CORRECTIONS:
            mouse, reason = CORRECTIONS[(t["order"], t["well"])]
            status = "excluded" if mouse.startswith("_excluded") else "relabelled"
        elif key in fates:
            f = fates[key]["fate"]
            mouse, status, reason = fate_dir.get(f, "_excluded/noisy_or_failed"), "excluded", f"nanoshm: {f}"
        elif group:
            mouse, status, reason = group[0], "relabelled", group[1]
        elif default.startswith("_excluded"):
            status, reason = "excluded", "no GenBank folder, and the label names no mouse in this study"
        tag = by_label.split("/")[-1]
        name = f"{mouse.split('/')[-1] if not mouse.startswith('_excluded') else tag}_{t['order']}-{t['well']}_{t['label']}.ab1"
        dest = out_p / mouse / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(t["path"], dest)
        rows.append({"mouse": mouse, "file": str(dest.relative_to(out_p)), "status": status,
                     "original_label": t["label"], "mouse_by_label": by_label, "order": t["order"],
                     "well": t["well"], "run_date": t["date"],
                     "mutations_q20_101_629": muts.get(key, ""),
                     "usable_q20_101_629": round(usable[key], 3) if key in usable else "",
                     "genbank_folder": fo or "",
                     "reason": reason,
                     "original_path": str(t["path"])})
    if genbank:
        from nanoshm.io_utils import read_fastq
        # Geneious cuts GenBank LOCUS names at 23 characters
        # ("23-jhextra-w3-2-user_ad" is the trace "23-jhextra-w3-2-user_added")
        stems = {t["path"].stem for t in traces} | {t["path"].stem[:23] for t in traces}
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
