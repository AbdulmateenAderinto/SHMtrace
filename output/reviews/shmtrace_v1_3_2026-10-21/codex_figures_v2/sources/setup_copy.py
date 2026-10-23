"""Freeze the reviewed figure sources and protect originals before design changes."""
import hashlib,json,shutil
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
ROOT=next(p for p in BASE.parents if (p/'docs/calibration/crm').is_dir())
CRM=ROOT/'docs/calibration/crm'
OLD=BASE.parent/'codex_version'
for name in ('figures','sources/originals','sources/revised','scratch','review','manuscript'):
 (BASE/name).mkdir(exist_ok=True,parents=True)
files=list((CRM/'figs').glob('*.py'))+list((CRM/'scicanvas').glob('*.js'))
for p in files: shutil.copy2(p,BASE/'sources/originals'/p.name)
for p in (CRM/'figs').glob('*.py'):
 src=OLD/'sources'/p.name
 shutil.copy2(src if src.exists() and p.name.startswith(('fig2_','fig3_','fig4_','fig5_')) else p,BASE/'sources/revised'/p.name)
# Keep independent copies of the original images and original manuscript operation list.
for p in (OLD/'figures').glob('*.png'):
 if not p.name.startswith('review_'): shutil.copy2(p,BASE/'review'/('before_'+p.name))
shutil.copy2(OLD/'document_operations.json',BASE/'manuscript/original_operations.json')
# Freeze sources and exports in their current state, including existing uncommitted work.
protected=files+list((CRM/'manuscript').glob('*.py'))+list((CRM/'figures').glob('*'))
protected += [OLD/'SHMtrace_CellReportsMethods_complete_v1.3_Codex_1.docx']
manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected if p.is_file()}
(BASE/'original_hashes.json').write_text(json.dumps(manifest,indent=2))
style=(BASE/'sources/revised/style.py').read_text()
style=style.replace('OUT = Path(__file__).resolve().parents[1] / "figures"',f'OUT = Path({str(BASE/"figures")!r})')
style=style.replace('DATA = Path(__file__).resolve().parents[1] / "out"',f'DATA = Path({str(CRM/"out")!r})')
style=style.replace('"font.size": 6.5','"font.size": 7').replace('"axes.labelsize": 6.5','"axes.labelsize": 7')
style=style.replace('"xtick.labelsize": 6','"xtick.labelsize": 6.5').replace('"ytick.labelsize": 6','"ytick.labelsize": 6.5')
style=style.replace('"legend.fontsize": 6','"legend.fontsize": 6.5')
(BASE/'sources/revised/style.py').write_text(style)
# Redirect explicit derived-output writes; inputs remain the saved results.
for p in (BASE/'sources/revised').glob('fig*.py'):
 s=p.read_text()
 import re
 s=re.sub(r'open\(DATA / (["\'][^"\']+["\']), ["\']w["\']\)',r'open(__import__("pathlib").Path(__file__).resolve().parents[2] / "scratch" / \1, "w")',s)
 s=s.replace('Path(__file__).resolve().parents[1] / "published/',f'Path({str(CRM)!r}) / "published/')
 s=s.replace('fig.savefig(DATA / "FigureS14_sanger_cohort2_matplotlib.png", dpi=300)','__import__("style").save(fig, "FigureS14_sanger_cohort2")')
 s=s.replace('fig.savefig(DATA / "FigureS15_sanger_mkln_matplotlib.png", dpi=300)','__import__("style").save(fig, "FigureS15_sanger_mkln")')
 p.write_text(s)
print(BASE)
