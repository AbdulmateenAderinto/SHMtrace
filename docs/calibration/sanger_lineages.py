"""Which mouse does each Sanger lineage belong to? Founder mutations, Sanger and nanopore.

    python docs/calibration/sanger_lineages.py [SHEET]   (default: samplesheets/jh4_manuscript.csv)

Clones of one B-cell lineage share the lineage's founder (trunk) mutations; mice make
their mutations independently, so a lineage belongs to one mouse. For the wild-type
clones (as filed in the sheet):

1. lineages: clones linked when they share 3 or more mutations that are also at least
   half of the lighter clone's mutations (connected components);
   founders: mutations carried by at least 80% of a lineage's clones that read the site.
2. nanopore: each lineage's founders are looked up in every wild-type mouse's nanopore
   library (all reads after filtering); a read carrying 4 or more counts, and only
   lineages with 4+ founders are tested. The owner is the mouse with 3+ such reads and
   at least five times any other mouse's.
3. a tree of the mutated clones (average linkage on the share of mutations not in
   common), coloured by the mouse each clone is filed under, labelled with its tube label.

Writes figures/sanger_lineages/{lineages.csv, tree.png, tree.pdf}.
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402

from shmtrace import pipeline, project, samplesheet  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "src/shmtrace/references"
OUT = Path(__file__).parent / "figures/sanger_lineages"
DEP = Path.home() / ("Documents/Projects/2026/PhD - Fundamental Immunology/Martin Lab/"
                     "Paper #1 - Nanopore sequencing/June 2nd analyses/submission_v4/"
                     "ONT_SHM_DataDeposit_Ver8_2026-07-30/Raw_Data/02_JH4_ONT/Original_FASTQ")
ONT = {"WT2": DEP / "DHJGKB_2_Wildtype_2.fastq", "WT3": DEP / "DHJGKB_3_Wildtype_3.fastq",
       "WT4": DEP / "DHJGKB_4_Wildtype_4.fastq"}
COL = {"WT2": "#C2548A", "WT3": "#2B6CB0", "WT4": "#2F855A"}
LO, HI = 101, 629


def label_mouse(label):
    toks = set(re.split(r"[^0-9a-z]+", label.lower()))
    for m, keys in (("WT2", {"w2", "wt2"}), ("WT3", {"w3", "wt3"}), ("WT4", {"w4", "wt4"})):
        if toks & keys:
            return m
    return ""


def main(sheet_path):
    OUT.mkdir(parents=True, exist_ok=True)
    sh = samplesheet.load(sheet_path)
    proj = sh.project(str(REF / "mouse_jh4_intron.fasta"), "jh4_intron",
                      contaminants=[str(REF / "mouse_jh4_contaminant.fasta")])
    res = project.run_project(proj)
    prov = sh.provenance()
    clones = {}
    for mouse, d in res["_details"].items():
        if not mouse.startswith("WT"):
            continue
        sites = {(v.pos, v.alt) for v in d["variants"] if LO <= v.pos <= HI}
        for a in d["_aligned"]:
            read = {rp + 1 for rp, b, q in a.pairs if q >= 20 and b in "ACGT" and LO <= rp + 1 <= HI}
            muts = {(rp + 1, b) for rp, b, q in a.pairs if q >= 20 and (rp + 1, b) in sites}
            lab = prov.get(a.qname, {}).get("label", "")
            clones[a.qname] = {"filed": mouse, "label": lab, "label_mouse": label_mouse(lab),
                               "muts": muts, "read": read}

    # lineages: share >= 3 mutations
    names = sorted(clones)
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            both = clones[a]["read"] & clones[b]["read"]
            A = {m for m in clones[a]["muts"] if m[0] in both}
            B = {m for m in clones[b]["muts"] if m[0] in both}
            sh = len(A & B)
            # 3+ shared, and most of the lighter clone's mutations: hotspots alone
            # (a few shared among many) would chain unrelated clones together
            if sh >= 3 and sh >= 0.5 * min(len(A), len(B)):
                parent[find(a)] = find(b)
    groups = defaultdict(list)
    for n in names:
        if clones[n]["muts"]:
            groups[find(n)].append(n)

    # nanopore libraries, all reads after filtering
    p = project.build_params({"region": "jh4_intron", "mode": "recommended"})
    ont_reads = {}
    for m, f in ONT.items():
        aligned, _ = pipeline.prepare_reads(str(f), str(REF / "mouse_jh4_intron.fasta"),
                                            [str(REF / "mouse_jh4_contaminant.fasta")], p)
        ont_reads[m] = [{(rp + 1, b) for rp, b, q in r.pairs if q >= 20} for r in aligned]

    rows = []
    for k, (_root, members) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
        cnt = Counter(m for n in members for m in clones[n]["muts"])
        founders = []
        for mut, c in cnt.items():
            readers = sum(1 for n in members if mut[0] in clones[n]["read"])
            if readers and c / readers >= 0.8 and (len(members) == 1 or c >= 2):
                founders.append(mut)
        founders = sorted(founders)
        # a read must carry 4 of the founders; lineages with fewer cannot be told from
        # hotspot coincidence. The owner is the mouse with 3+ such reads and at least
        # five times any other's.
        testable = len(founders) >= 4
        ont = {m: (sum(1 for r in reads if len(r & set(founders)) >= 4) if testable else 0)
               for m, reads in ont_reads.items()}
        top = sorted(ont.values(), reverse=True)
        tot = sum(ont.values())
        owner = (max(ont, key=ont.get) if testable and top[0] >= 3 and top[0] >= 5 * top[1]
                 else "")
        rows.append({
            "lineage": f"L{k + 1}", "clones": len(members),
            "filed": "; ".join(f"{m} {c}" for m, c in Counter(clones[n]["filed"] for n in members).most_common()),
            "labels": "; ".join(f"{m or '?'} {c}" for m, c in Counter(clones[n]["label_mouse"] for n in members).most_common()),
            "founders": len(founders),
            "ont_WT2": ont["WT2"], "ont_WT3": ont["WT3"], "ont_WT4": ont["WT4"],
            "nanopore_says": (owner or ("unclear" if tot else "not seen")) if testable
            else "too few founders",
            "members": " ".join(f"{clones[n]['label']}({clones[n]['filed']})" for n in members),
        })
        for n in members:
            clones[n]["lineage"] = f"L{k + 1}"
            clones[n]["owner"] = owner
    with open(OUT / "lineages.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # tree of the mutated clones
    mut = [n for n in names if clones[n]["muts"]]
    D = np.zeros((len(mut), len(mut)))
    for i, a in enumerate(mut):
        for j in range(i + 1, len(mut)):
            b = mut[j]
            both = clones[a]["read"] & clones[b]["read"]
            A = {m for m in clones[a]["muts"] if m[0] in both}
            B = {m for m in clones[b]["muts"] if m[0] in both}
            u = A | B
            D[i, j] = D[j, i] = 1 - len(A & B) / len(u) if u else 1.0
    Z = linkage(squareform(D), "average")
    fig, ax = plt.subplots(figsize=(7.5, max(6, 0.2 * len(mut))))
    dn = dendrogram(Z, orientation="left", labels=mut, ax=ax, color_threshold=0,
                    above_threshold_color="#888", leaf_font_size=7)
    ax.set_yticklabels([])
    for ypos, leaf in zip(ax.get_yticks(), dn["ivl"], strict=True):
        c = clones[leaf]
        txt = f"{c['label']}  [{c.get('lineage', '')}{'→' + c['owner'] if c.get('owner') and c['owner'] != c['filed'] else ''}]"
        ax.text(-0.02, ypos, txt, ha="right", va="center", fontsize=6.5, color=COL[c["filed"]],
                transform=ax.get_yaxis_transform())
        ax.text(1.005, ypos, f"{len(c['muts'])}", ha="left", va="center", fontsize=6,
                color="#666", transform=ax.get_yaxis_transform())
    for m, col in COL.items():
        ax.plot([], [], color=col, lw=6, label=f"filed under {m}")
    ax.legend(loc="lower left", fontsize=8, frameon=False, bbox_to_anchor=(0, 1.0), ncol=3)
    ax.set_xlabel("share of mutations not in common (average linkage)")
    ax.set_title("Wild-type Sanger clones. [lineage→mouse whose nanopore reads carry its "
                 "founders, where it differs from the filing]; right: mutations", fontsize=8, pad=22)
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.subplots_adjust(left=0.32)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"tree.{ext}", dpi=170, bbox_inches="tight")
    for r in rows:
        print({k: v for k, v in r.items() if k != "members"})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).parent / "samplesheets/jh4_manuscript.csv"))
