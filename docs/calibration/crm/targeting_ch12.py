"""S5F targeting index of pooled calls per CH12F3 condition (PRJNA918596 rerun), window 21-249.
    python docs/calibration/crm/targeting_ch12.py   -> out/targeting_ch12.json"""
import json
from pathlib import Path

from shmtrace import targeting
from shmtrace.io_utils import read_fasta

OUT = Path(__file__).parent / "out"
p = json.load(open(OUT / "public/ch12/project.json"))
ref = read_fasta(str(Path.home() / "Nanoshm/public_data/PRJNA918596_ch12_ref.fasta"))[0][1].upper()
model = targeting.load(str(Path.home() / "Nanoshm/models/HH_S5F.csv"))
prof = targeting.window_profile(model, ref, 20, 248)
res = {}
for g in ("AID", "R174ER178E", "R50E", "EV"):
    pos = [v["pos"] - 1 for v in p["variants"] if v["group"] == g]
    res[g] = {**targeting.targeting_index(pos, prof), "calls": len(pos)}
json.dump(res, open(OUT / "targeting_ch12.json", "w"), indent=1)
print(res)
