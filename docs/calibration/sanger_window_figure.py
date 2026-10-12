"""Is the 529-bp Sanger window (101-629) defensible? Coverage, mutation density and the J꜀4
exon/intron boundary along the amplicon, wild-type Sanger clones.

    python docs/calibration/sanger_window_figure.py [--sanger-dir ~/Nanoshm/sanger_ab1_by_mouse]

Writes figures/sanger_vs_nanopore/figS_sanger_window.{png,pdf,svg} and window.json.
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nanoshm import motifs, pipeline, project  # noqa: E402
from nanoshm.io_utils import read_abif, read_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REF = str(ROOT / "src/nanoshm/references/mouse_jh4_intron.fasta")
CON = [str(ROOT / "src/nanoshm/references/mouse_jh4_contaminant.fasta")]
OUT = ROOT / "docs/calibration/figures/sanger_vs_nanopore"
WT = ("WT2", "WT3", "WT4")
KO = ("Aicda_A1", "Aicda_A2", "Aicda_A3")


def clones(folder: Path, mice, tmp: Path):
    fq = tmp / "clones.fastq"
    with open(fq, "w") as fh:
        for m in mice:
            for f in sorted((folder / m).glob("*.ab1")):
                _n, s, q = read_abif(str(f))
                fh.write(f"@{m}|{f.name}\n{s}\n+\n{q}\n")
    p = project.build_params({"region": "jh4_intron", "platform": "sanger"})
    return pipeline.prepare_reads(str(fq), REF, CON, p)[0], p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sanger-dir", default=str(Path.home() / "Nanoshm/sanger_ab1_by_mouse"))
    ap.add_argument("--name", default="figS_sanger_window", help="output file name, without extension")
    a = ap.parse_args()
    import tempfile
    ref = read_fasta(REF)[0][1].upper()
    L = len(ref)
    with tempfile.TemporaryDirectory() as t:
        wt, p = clones(Path(a.sanger_dir), WT, Path(t))
        ko, _ = clones(Path(a.sanger_dir), KO, Path(t))
    germ = project.sanger_germline({"wt": wt}, ref, 0, L - 1, p)
    n = len(wt)
    cov = [0] * L
    muts = set()
    for c in wt:
        for pos, b, q in c.pairs:
            if q >= 20 and b in "ACGT":
                cov[pos] += 1
                if b != ref[pos] and (pos, b) not in germ:
                    muts.add((c.qname.split("|")[0], pos, b))
    ko_bases = ko_mm = 0
    for c in ko:
        for pos, b, q in c.pairs:
            if 100 <= pos <= 628 and q >= 20 and b in "ACGT":
                ko_bases += 1
                ko_mm += b != ref[pos] and (pos, b) not in germ
    exon_end = ref.find("CTCCTCAGGTAAG") + 8           # last coding base, 1-based

    # unique mutations per mouse, per 25-bp bin, over the clones covering each bin
    bins = list(range(51, 676, 25))
    dens = []
    for lo in bins:
        hi = lo + 24
        k = sum(1 for _m, pos, _b in muts if lo - 1 <= pos <= hi - 1)
        bp = sum(cov[i] for i in range(lo - 1, min(hi, L)))
        dens.append(1e3 * k / bp if bp else 0)

    def window(lo, hi):
        spans = sum(1 for c in wt if min(x for x, _b, _q in c.pairs) <= lo - 1
                    and max(x for x, _b, _q in c.pairs) >= hi - 1)
        k = sum(1 for _m, pos, _b in muts if lo - 1 <= pos <= hi - 1)
        gc = [(pos) for _m, pos, b in muts if lo - 1 <= pos <= hi - 1 and ref[pos] in "GC"]
        hot = sum(1 for pos in gc if motifs.is_site("WRC", ref, pos) or motifs.is_site("GYW", ref, pos))
        return {"clones_spanning": spans, "of": n, "unique_mutations": k,
                "gc_at_wrc_gyw": f"{hot}/{len(gc)}"}

    stats = {"clones": n, "exon_end": exon_end, "intron_start": exon_end + 1,
             "ko_mismatches_in_window": f"{ko_mm}/{ko_bases}",
             "windows": {f"{lo}-{hi}": window(lo, hi) for lo, hi in ((101, 629), (74, 629), (101, 656))},
             "coverage_at": {x: cov[x - 1] for x in (53, 60, 74, 101, 300, 500, 600, 629, 640, 656, 670)}}
    json.dump(stats, open(OUT / f"{a.name}.json", "w"), indent=1)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.spines.top": False,
                         "savefig.dpi": 300})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1], "hspace": .12})
    for ax in (ax1, ax2):
        ax.axvspan(101, 629, color="#2F6BA6", alpha=.10, lw=0)
        ax.axvspan(74, 101, color="#E69F00", alpha=.12, lw=0)
        ax.axvspan(629, 656, color="#E69F00", alpha=.12, lw=0)
        ax.axvline(exon_end + .5, color="#9A3439", lw=1.2, ls="--")
    xs = list(range(1, L + 1))
    ax1.plot(xs, [100 * c / n for c in cov], color="#1F4E79", lw=1.6)
    ax1.set_ylabel(f"Clones with a Q≥20 base (%)\n(wild type, n = {n})")
    ax1.set_ylim(0, 105)
    ax1.set_xlim(40, 690)
    ax1.text(exon_end - 3, 75, "← JH4\ncoding", fontsize=9, color="#9A3439", va="top", ha="right")
    ax1.text(exon_end + 4, 75, "intron →", fontsize=9, color="#9A3439", va="top")
    ax1.text(365, 8, "529-bp window, 101–629", ha="center", color="#1F4E79", fontsize=10, fontweight="bold")
    ax1.text(87, 8, "+27\n(74–100)", ha="center", color="#9A6000", fontsize=8.5)
    ax1.text(642, 8, "+27\n(630–656)", ha="center", color="#9A6000", fontsize=8.5)
    ax2.bar([b + 12 for b in bins], dens, width=23, color="#5A9BD5", edgecolor="#1F4E79", lw=.6)
    ax2.set_ylabel("Unique mutations per 10³\nbases sequenced, 25-bp bins")
    ax2.set_xlabel("Position in the J꜀4 reference (bp)".replace("J꜀4", "JH4"))
    for ax in (ax1, ax2):
        ax.yaxis.grid(True, ls=(0, (5, 4)), color="#C9C9C9")
        ax.set_axisbelow(True)
    ax1.set_title("Where the Sanger clones are readable, and where the mutations are", loc="left",
                  fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{a.name}.{ext}", bbox_inches="tight")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
