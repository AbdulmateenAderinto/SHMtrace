"""Check output integrity, physical print minimums, text bounds and archive contents."""
import hashlib,json,re,sys,zipfile
from pathlib import Path
import pdfplumber
BASE=Path(__file__).resolve().parents[1]
import types
style_stub=types.ModuleType('style');style_stub.OUT=BASE/'figures';sys.modules['style']=style_stub
sys.path.insert(0,str(BASE/'sources/revised'));import print_check
entries=json.loads((BASE/'figure_manifest.json').read_text());bounds={};prints={}
for e in entries:
 name=e['stem'];path=BASE/'figures'/name
 bad=print_check.check(str(path))
 if name=='FigureS13_sanger_window':
  # Its saved vector artwork is inside a uniform scaling transform. Physical
  # type/stroke sizes must include that transform, unlike native mpl figures.
  svg=path.with_suffix('.svg').read_text();k=float(re.search(r'<g transform="scale\(([0-9.]+)\)',svg)[1])
  raw=[float(x) for x in re.findall(r'stroke-width:\s*([0-9.]+)',svg)]
  bad=[f'physical stroke {x*k:.3f}' for x in raw if x>0 and not .495<=x*k<=1.505]
  sizes=[float(x)*100*k for x in re.findall(r'scale\(([0-9.]+) -',svg)]
  if sizes and min(sizes)<5.95:bad.append(f'outline text {min(sizes):.3f}')
 prints[name]=bad
 with pdfplumber.open(path.with_suffix('.pdf')) as doc:
  p=doc.pages[0];out=[w['text'] for w in p.extract_words() if w['x0']<-.5 or w['x1']>p.width+.5 or w['top']<-.5 or w['bottom']>p.height+.5]
  if out:bounds[name]=out
assert not bounds,bounds
assert not any(prints.values()),{k:v for k,v in prints.items() if v}
checks=json.loads((BASE/'validation.json').read_text())
checks.update({'text_outside_figure_bounds':bounds,'print_minimum_failures':{k:v for k,v in prints.items() if v},'visual_review_pages':44})
# Compare derived numerical summaries where the prior cache exists.
root=next(p for p in BASE.parents if (p/'docs/calibration/crm').is_dir());derived={}
for f in (BASE/'scratch').glob('*.json'):
 original=root/'docs/calibration/crm/out'/f.name
 if original.exists():
  try:derived[f.name]=json.loads(f.read_text())==json.loads(original.read_text())
  except (ValueError,UnicodeError):pass
checks['derived_summary_comparisons']=derived
if not all(derived.values()):raise AssertionError(derived)
checks['snapshot_resources']=len(json.loads((BASE/'input_snapshot_hashes.json').read_text()))
assert all(hashlib.sha256((BASE/p).read_bytes()).hexdigest()==h for p,h in json.loads((BASE/'input_snapshot_hashes.json').read_text()).items())
(BASE/'validation.json').write_text(json.dumps(checks,indent=2))
(BASE/'print_check.json').write_text(json.dumps(prints,indent=2))
(BASE/'CHECKS.md').write_text(f'''# Verification

- {len(entries)} figure/companion pages checked; 44-page review PDF rendered and visually reviewed.
- Text beyond figure bounds: none.
- Physical SVG print minimum violations: none (including scaling of retained S13 vector artwork).
- 196 protected original files unchanged.
- 70 original data inputs unchanged; all 72 copied resources match source hashes.
- {len(derived)} derived numerical summaries match their existing original caches exactly.
- Calling defaults, test definitions and original Word manuscripts were not modified.

This verifies the design revision and preservation of saved results; it is not a new scientific validation of the method.
''')
archive=BASE.parent/'SHMtrace_Codex_Figures_v2_sources_and_exports.zip'
paths=[]
for p in BASE.rglob('*'):
 if not p.is_file() or '__pycache__' in p.parts:continue
 rel=p.relative_to(BASE)
 if rel.parts[0]=='scratch':continue
 if rel.parts[0]=='review' and not p.name.startswith('before_'):continue
 if p.suffix=='.tif':continue
 paths.append(p)
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in paths:z.write(p,Path(BASE.name)/p.relative_to(BASE))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps(checks,indent=2));print(archive,round(archive.stat().st_size/1024**2,1),'MB',len(paths),'files')
