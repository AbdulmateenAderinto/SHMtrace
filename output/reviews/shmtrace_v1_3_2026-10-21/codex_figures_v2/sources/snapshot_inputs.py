"""Copy only resources actually read by the figure builds (not the full 1.9-GB analysis folder)."""
import hashlib,json,shutil
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
ROOT=next(p for p in BASE.parents if (p/'docs/calibration/crm').is_dir())
files=list(json.loads((BASE/'input_hashes.json').read_text()))
for ext in ('svg','pdf'):
 files.append(str(ROOT/f'docs/calibration/figures/sanger_vs_nanopore/figS_sanger_window_original.{ext}'))
mapping={};hashes={}
for path in files:
 p=Path(path);q=BASE/'inputs'/p.relative_to(ROOT);q.parent.mkdir(parents=True,exist_ok=True)
 shutil.copy2(p,q);mapping[path]=str(q.relative_to(BASE));hashes[str(q.relative_to(BASE))]=hashlib.sha256(q.read_bytes()).hexdigest()
(BASE/'input_snapshot_map.json').write_text(json.dumps(mapping,indent=2))
(BASE/'input_snapshot_hashes.json').write_text(json.dumps(hashes,indent=2))
print(len(mapping),'resources copied;',round(sum((BASE/v).stat().st_size for v in mapping.values())/1024**2,1),'MB')
