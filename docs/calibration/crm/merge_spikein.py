"""Replace SHMtrace's (recommended) rows in out/spikein_<key>.json with those of a rerun (_rec).
    python docs/calibration/crm/merge_spikein.py jh4 n114 su_lps pim1"""
import json
import sys
from pathlib import Path

OUT = Path(__file__).parent / "out"
for k in sys.argv[1:]:
    main = json.load(open(OUT / f"spikein_{k}.json"))
    new = json.load(open(OUT / f"spikein_{k}_rec.json"))
    rows = [r for r in main["rows"] if r["caller"] != "recommended"] + \
           [r for r in new["rows"] if r["caller"] == "recommended"]
    main["rows"] = rows
    main["recommended_from"] = "control-based FDR rerun (spikein_<key>_rec.json)"
    json.dump(main, open(OUT / f"spikein_{k}.json", "w"), indent=1)
    print(k, sum(r["caller"] == "recommended" for r in rows), "recommended rows merged")
