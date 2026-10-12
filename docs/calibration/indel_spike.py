"""Spike-in check of the indel caller on the J꜀4 knockouts.

Plants insertions and deletions of known size into one knockout at a known
number of reads, calls it against a null fitted (and calibrated) on the other
two knockouts only, and scores recall by read support and false calls.

    python docs/calibration/indel_spike.py [draws per knockout] [indel min quality]
"""
import random
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from nanoshm import batch, indels, simulate
from nanoshm.io_utils import read_fasta, read_fastq
from nanoshm.project import build_params
from nanoshm.variants import homopolymer_run

ROOT = Path(__file__).resolve().parents[2]
REF = str(ROOT / "src/nanoshm/references/mouse_jh4_intron.fasta")
CON = str(ROOT / "src/nanoshm/references/mouse_jh4_contaminant.fasta")
EX = ROOT / "examples_compare"
KOS = ["KO1", "KO2", "KO3"]
DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 5
MINQ = int(sys.argv[2]) if len(sys.argv) > 2 else None   # None: the default
SUPPORT = (1, 2, 3, 5, 10, 20)

p = build_params({"region": "jh4_intron"})
ref_name, ref = read_fasta(REF)[0]
lo, hi = p.analysis_window[0] - 1, p.analysis_window[1] - 1
tmp = Path(tempfile.mkdtemp())


def pick_sites(rng, n=10):
    """Spaced sites of mixed kinds; returns (pos0, kind, bases, support, label)."""
    kinds = ["del1_hp", "del1", "ins1", "del2", "ins2", "del3", "dellong"]
    sites, used = [], []
    i = tries = 0
    while len(sites) < n and tries < 10000:
        tries += 1
        kind = rng.choice(kinds)
        L = {"del1_hp": 1, "del1": 1, "ins1": 1, "del2": 2, "ins2": 2, "del3": 3,
             "dellong": rng.randint(10, 40)}[kind]
        pos = rng.randint(lo + 5, hi - L - 5)
        hp = homopolymer_run(ref, pos)
        if (kind == "del1_hp") != (hp >= 3):
            continue
        if any(abs(pos - u) < 10 + max(L, ul) for u, ul in used):
            continue
        if kind.startswith("ins"):
            bases = "".join(rng.choice("ACGT") for _ in range(L))
            site = (pos, "ins", bases)
        else:
            site = (pos, "del", ref[pos:pos + L])
        sites.append((*site, rng.choice(SUPPORT), kind))
        used.append((pos, L))
        i += 1
    return sites


def pileup_of(path, name):
    b = batch.prepare_batch({name: str(path)}, REF, [CON], p)
    return indels.Pileup(b["reads"][name], ref, lo, hi,
                         p.indel_min_quality if MINQ is None else MINQ)


recall = defaultdict(lambda: [0, 0])          # support -> [found, planted]
recall_kind = defaultdict(lambda: [0, 0])     # kind, at support >= 3
grid = defaultdict(lambda: [0, 0])            # (kind, support) -> [found, planted]
false_calls = []
false_detail = []                             # (kind, length, reads, held-out library)
for held in KOS:
    peers = [k for k in KOS if k != held]
    peer_piles = [pileup_of(EX / f"{c}.fastq", c) for c in peers]
    qscale, _ = indels.calibrate(peer_piles, ref, p)
    null = indels.IndelNull(peer_piles, ref)
    null.qscale = qscale
    reads = list(read_fastq(str(EX / f"{held}.fastq")))
    for draw in range(DRAWS):
        rng = random.Random(100 * draw + KOS.index(held))
        sites = pick_sites(rng)
        spiked, planted = simulate.spike_indels(reads, REF, [s[:4] for s in sites], p, rng)
        f = tmp / f"{held}_{draw}.fastq"
        with open(f, "w") as fh:
            for name, seq, qual in spiked:
                fh.write(f"@{name}\n{seq}\n+\n{qual}\n")
        called, _st = indels.call(pileup_of(f, "S"), null, p)
        got = {(v.pos - 1, v.kind, v.seq) for v in called}
        kind_of = {}
        for (pos, kind, bases, k, label) in sites:
            n = indels.normalise(ref, pos, "D" if kind == "del" else "I",
                                 len(bases) if kind == "del" else bases)
            kind_of[n] = label
        for key, k in planted.items():
            if k == 0:
                continue
            recall[k][1] += 1
            recall[k][0] += key in got
            grid[(kind_of[key], k)][1] += 1
            grid[(kind_of[key], k)][0] += key in got
            if k >= 3:
                recall_kind[kind_of[key]][1] += 1
                recall_kind[kind_of[key]][0] += key in got
        false_calls.append(len(got - set(planted)))
        for v in called:
            if (v.pos - 1, v.kind, v.seq) not in planted:
                false_detail.append((v.kind, v.length, v.supporting_reads, held))
        print(f"{held} draw {draw}: planted {sum(1 for v in planted.values() if v)}, "
              f"called {len(got)}, false {false_calls[-1]}, qscale {qscale}", flush=True)

print("\nrecall by supporting reads")
for k in sorted(recall):
    f, n = recall[k]
    print(f"  {k:>2} reads: {f}/{n} = {f / n:.2f}")
print("recall by kind (support >= 3)")
for kind, (f, n) in sorted(recall_kind.items()):
    print(f"  {kind:>8}: {f}/{n} = {f / n:.2f}")
print("recall, kind x supporting reads")
print("          " + "".join(f"{k:>8}" for k in SUPPORT))
for kind in sorted({k for k, _ in grid}):
    cells = []
    for k in SUPPORT:
        f, n = grid.get((kind, k), (0, 0))
        cells.append(f"{f}/{n}".rjust(8) if n else "       -")
    print(f"{kind:>10}" + "".join(cells))
print(f"false calls per library: {sum(false_calls) / len(false_calls):.2f} "
      f"({sum(false_calls)} in {len(false_calls)})")
print("distinct false events (kind, length, reads, held-out library):")
for ev in sorted(set(false_detail)):
    print("  ", ev, "x", false_detail.count(ev))
