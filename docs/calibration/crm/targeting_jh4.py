"""AID targeting against sequencing error on J_H4, with the SHazaM HH_S5F 5-mer model.

    python docs/calibration/crm/targeting_jh4.py
For the 529-bp window: the targeting index of each Wt library's calls; the index of every position
at which any analysed Aicda-/- read mismatches the reference (sequencing error, no AID); and the
dose-response: distinct Wt calls and Aicda-/- mismatch events per position, by quintile of model
mutability. Same reads as the manuscript (delivered, 500 drawn per library, seed 0).
Writes out/targeting_jh4.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

from shmtrace import pipeline, targeting
from shmtrace.config import load_preset
from shmtrace.io_utils import read_fasta
from shmtrace.reads_filter import downsample

sys.path.insert(0, str(Path(__file__).parent))
from datasets import DATASETS  # noqa: E402

DS = DATASETS["jh4"]
OUT = Path(__file__).parent / "out"
model = targeting.load(str(Path.home() / "Nanoshm/models/HH_S5F.csv"))
ref = read_fasta(DS["reference"])[0][1].upper()
p = load_preset("jh4_intron")
lo, hi = p.analysis_window[0] - 1, p.analysis_window[1] - 1
prof = targeting.window_profile(model, ref, lo, hi)
proj = json.load(open(Path(__file__).parents[1] / "figures/jh4_shmtrace/project.json"))
res = {"wt": {}, "ko_mismatch": {}}
wt_pos = {}
for s in ("WT2", "WT3", "WT4"):
    pos = [v["pos"] - 1 for v in proj["variants"] if v["sample"] == s]
    wt_pos[s] = pos
    res["wt"][s] = targeting.targeting_index(pos, prof)
ko_events = {}
for s in ("KO1", "KO2", "KO3"):
    al, _ = pipeline.prepare_reads(DS["samples"][s]["fastq"], DS["reference"], DS["contaminants"], p)
    al, _ = downsample(al, p)
    positions, events = set(), 0
    for a in al:
        for rpos, base, _q in a.pairs:
            if base and lo <= rpos <= hi and base.upper() != ref[rpos]:
                positions.add(rpos)
                ko_events[rpos] = ko_events.get(rpos, 0) + 1
                events += 1
    res["ko_mismatch"][s] = {**targeting.targeting_index(sorted(positions), prof), "events": events}
# germline differences and shared polymorphisms (>= 20% of the knockout reads) are not sequencing error
n_ko_reads = 3 * p.downsample_to
germ = {k for k, v in ko_events.items() if v >= 0.2 * n_ko_reads}
ko_events = {k: v for k, v in ko_events.items() if k not in germ}
res["germline_positions_excluded"] = sorted(int(x) + 1 for x in germ)
allko = sorted(ko_events)
res["ko_all_positions"] = targeting.targeting_index(allko, prof)
# dose-response by mutability quintile
mut = prof["mutability"]
order = sorted(mut, key=mut.get)
quint = np.array_split(order, 5)
wt_all = [q for s in wt_pos for q in wt_pos[s]]
res["quintiles"] = [{"q": i + 1, "positions": len(qs), "mean_mutability": float(np.mean([mut[x] for x in qs])),
                     "wt_calls_per_position": sum(1 for x in wt_all if x in set(qs)) / len(qs),
                     "ko_mismatch_events_per_position": sum(ko_events.get(int(x), 0) for x in qs) / len(qs)}
                    for i, qs in enumerate(quint)]
OUT.mkdir(exist_ok=True)
json.dump(res, open(OUT / "targeting_jh4.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float)[:3000])
