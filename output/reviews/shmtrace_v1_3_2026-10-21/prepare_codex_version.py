"""Prepare an isolated manuscript revision and figure sources; do not edit v1.3."""
# ruff: noqa: E402
import hashlib
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parents[2]
CRM = ROOT / 'docs/calibration/crm'
DEST = HERE / 'codex_version'
ASSETS = DEST / 'figures'
SOURCES = DEST / 'sources'
for p in (DEST, ASSETS, SOURCES, DEST/'scratch'):
    p.mkdir(exist_ok=True)
sys.path.insert(0, str(CRM/'manuscript'))
sys.path.insert(0, str(CRM))
import build_crm as B
import proposed_edits as P
import agreement_table

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

protected = list((CRM/'manuscript').glob('*.py')) + list((CRM/'figures').glob('*.png'))
protected += [CRM/'build/SHMtrace_CellReportsMethods_complete_v1.3.docx', CRM/'build/SHMtrace_CellReportsMethods.docx']
manifest = {str(p):sha(p) for p in protected}
(DEST/'original_file_hashes.json').write_text(json.dumps(manifest,indent=2))
raw = json.loads((HERE/'source_blocks.json').read_text())
for key, value in raw.items():
    assert json.loads(json.dumps(getattr(B.M if key=='LEGENDS' else B.T,key))) == value, key
numbers = {**B.crm_numbers.compute(), **B.public_numbers(), **B.dataset_counts()}
assert numbers == json.loads((HERE/'numbers.json').read_text()), 'Analysis values changed since review'
(DEST/'numbers.json').write_text(json.dumps(numbers,indent=2,ensure_ascii=False))
new = json.loads((HERE/'proposal_source_blocks.json').read_text())
for key,value in new.items():
    setattr(B.M if key=='LEGENDS' else B.T,key,value)

# Small language and interpretation corrections identified in the review.
B.T.INTRO[0] = B.T.INTRO[0].replace('Activation-induced cytidine deaminase (AID) which',
    'activation-induced cytidine deaminase (AID), which').replace('are therefore key readouts','are key readouts')
B.T.RESULTS[11] = (B.T.RESULTS[11][0], B.T.RESULTS[11][1].replace(
    'at false-discovery proportions of {cd_fdp}%.',
    'with false-discovery proportions of {cd_fdp}% for the two- to four-control designs.'))
B.M.LEGENDS[2] = (B.M.LEGENDS[2][0], B.M.LEGENDS[2][1].replace(
    'zero is displayed at 0.05 on the logarithmic axis.',
    'zero remains zero on a symmetric-log axis, linear from zero to 0.1 false calls/library.'))
B.M.LEGENDS[4] = (B.M.LEGENDS[4][0], B.M.LEGENDS[4][1].replace(
    'Positional shares in 25-bp bins.', 'Positional shares in 25-bp bins, with the final four-base interval marked separately.'))
B.M.LEGENDS[3] = (B.M.LEGENDS[3][0], B.M.LEGENDS[3][1].replace(
    'off-target Wt and double-knockout groups', 'OT1/OT2 denote off-targets 1/2; their Wt and double-knockout groups'))
B.M.LEGENDS[6] = ('Figure 7. Public Illumina datasets test analysis across biological contexts',B.M.LEGENDS[6][1])

# Integrate the detail bank into the relevant existing Methods subsections.
destinations = {
    'Sanger counting and provenance':'Single-read estimate, clones, mutation trees, and quality checks',
    'Sanger comparison details':'Single-read estimate, clones, mutation trees, and quality checks',
    'Maea comparison details':'Single-read estimate, clones, mutation trees, and quality checks',
    'Two-candidate example':'Control-calibrated variant calling',
    'Caller operating ranges':'Benchmark against other callers',
    'Real-library caller details':'Benchmark against other callers',
    'Somatic caller spike-in false calls':'Benchmark against other callers',
    'PacBio junction comparison details':'Public PacBio class-switch junctions',
    'Ramos control sensitivity':'Public Illumina and UMI datasets',
    'mOrange raw-background comparison':'Public Illumina and UMI datasets',
}
methods=[]
for i,(kind,text) in enumerate(B.M.METHODS):
    for before,after in P.METHODS_REPLACEMENTS:
        text=text.replace(before,after)
    methods.append((kind,text))
    if i+1==len(B.M.METHODS) or B.M.METHODS[i+1][0]=='h2':
        section=next(t for k,t in reversed(methods) if k=='h2')
        for label,detail in P.METHODS_ADDITIONS:
            if destinations[label]==section:
                methods.append(('p',detail))
B.M.METHODS=methods
B.M.TABLE_S6 = [tuple('Figure S8D' if i==6 and 'Karst et al.' in row[0] else
    'Figure 3; STAR Methods' if i==6 and row[0]=='Simulations' else x
    for i,x in enumerate(row)) for row in B.M.TABLE_S6]
# The agreement table retained pre-renumbering destinations outside crm_methods.py.
mapping={'8':'14','9':'15','10':'16'}
agreement_table.ROWS=[tuple(re.sub(r'Figure S(8|9|10)(?=[A-Z\s)])',
    lambda m:'Figure S'+mapping[m[1]],x) for x in row) for row in agreement_table.ROWS]

# Copy every embedded image into the isolated version and remap only this build.
image_paths=[CRM/'figures'/f'Figure{i}{"_preview" if i==1 else ""}.png' for i in range(1,8)]
image_paths += [Path(p) for _,paths,_ in B.SUPP_FIGS for p in (paths if isinstance(paths,list) else [paths])]
for p in set(image_paths):
    shutil.copy2(p,ASSETS/p.name)
B.SUPP_FIGS=[(title,[str(ASSETS/Path(p).name) for p in paths] if isinstance(paths,list)
    else str(ASSETS/Path(paths).name),legend) for title,paths,legend in B.SUPP_FIGS]
B.FIG=ASSETS

def replace_once(s,before,after):
    assert s.count(before)==1, before[:100]
    return s.replace(before,after)

# Isolated matplotlib builders; all data reads still use the original saved outputs.
fig2=(CRM/'figs/fig2_call.py').read_text()
fig2=fig2.replace('(Fig. 3H)','(Fig. 3C)')
fig2=fig2.replace('open(DATA / "fig2_numbers.json", "w")',f'open({str(DEST/"scratch/fig2_numbers.json")!r}, "w")')
fig3=(CRM/'figs/fig3_calibration.py').read_text()
fig3=replace_once(fig3,'max(s[c]["fp"], 0.05)','s[c]["fp"]')
fig3=replace_once(fig3,'ax.set_xscale("log")\n    ax.set_xlim(0.04, 20)',
    'ax.set_xscale("symlog", linthresh=0.1, linscale=0.55)\n    ax.set_xlim(0.04, 20)')
fig3=replace_once(fig3,'ax.set_xlim(0.04, 20)','ax.set_xlim(-0.03, 20)\n    ax.set_xticks([0, 0.1, 1, 10], ["0", "0.1", "1", "10"])\n    ax.minorticks_off()')
fig3=fig3.replace('False calls per library (log)','False calls per library (symlog)')
fig3=replace_once(fig3,'    ax.text(0.045, 0.13, "0 plotted\\nat 0.05", fontsize=6, color=SOFT, transform=ax.transData)\n','')
fig4=(CRM/'figs/fig4_controls.py').read_text()
fig4=fig4.replace('168 * MM','176 * MM').replace('bottom=0.06,','bottom=0.11,')
fig4=fig4.replace('labs = [lbl for _k, lbl, _p in avail]',
    'labs = [lbl.replace("Off-target 1", "OT1").replace("Off-target 2", "OT2") for _k, lbl, _p in avail]')
fig4=fig4.replace('letter(ax, "FG"[p], dx=-0.3)','letter(ax, "FG"[p], dx=-0.22)')
fig4=fig4.replace('open(DATA / "fig4_tests.json", "w")',f'open({str(DEST/"scratch/fig4_tests.json")!r}, "w")')
fig5=(CRM/'figs/fig5_sanger.py').read_text()
fig5=replace_once(fig5,'mids = (edges[:-1] + edges[1:]) / 2',
    'mids = (edges[:-1] + edges[1:]) / 2\n    mids[-1] = (626 + 629) / 2\n    ax.axvspan(626, 630, color=SOFT, alpha=0.15, lw=0)')
fig5=replace_once(fig5,'Where calls fall, 25-bp bins','Where calls fall (terminal bin: 4 bp)')
for name,src in [('fig2_call',fig2),('fig3_calibration',fig3),('fig4_controls',fig4),('fig5_sanger',fig5)]:
    (SOURCES/f'{name}.py').write_text(src)
sys.path.insert(0,str(CRM/'figs'))
import style
style.OUT=ASSETS
for name in ('fig2_call','fig3_calibration','fig4_controls','fig5_sanger'):
    spec=importlib.util.spec_from_file_location(name,SOURCES/f'{name}.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.main()

# Updated four-platform graphic. Current source already has the new S13/S17 refs;
# the native preview had been stale. Keep both source and exported image in this version.
s21=(CRM/'scicanvas/figS21_platforms.js').read_text()
s21=replace_once(s21,'447 chromatograms from the six mice of Figure 5',
    f'{numbers["sanger_sheet"]} assigned traces; {numbers["sanger_clones"]} analyzed clones (Figure 5)')
s21=s21.replace('r = 0.92 (95% CI 0.74–0.98), slope 0.97','r = 0.92 (pooled categories), slope 0.97')
s21=s21.replace('tubes (Figure S13)','tubes (STAR Methods)')
(SOURCES/'figS21_platforms.js').write_text(s21)

B.fill_modules(numbers)
B.check_paired_settings(numbers)
B.check_limits()
total,parts=B.journal_word_count()
assert total<7000,total

class Recorder:
    def __init__(self):self.ops=[]
    def __getattr__(self,name):
        def record(*args,**kwargs):self.ops.append([name,args,kwargs])
        return record

r=Recorder()
r.para('SHMtrace Codex revision',size=16,bold=True)
r.para('Version 1.3 Codex 1 · 21 October 2026',size=11)
r.page_break()
B.write_article(r,numbers)
r.page_break();r.heading('Document S1. Supplemental figures and tables',1)
B.write_supplement(r)
r.references()
def encode(x):
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x))
(DEST/'document_operations.json').write_text(json.dumps(r.ops,indent=2,ensure_ascii=False,default=encode))
(DEST/'counts.json').write_text(json.dumps({'version':'1.3 Codex 1','journal_count':total,
    'parts':dict(parts),'original':9139},indent=2))
# Archive the original builder and typography helpers alongside the revision.
for name in ('crm_text.py','crm_methods.py','build_crm.py','crm_numbers.py','docbuild.py','refs.json'):
    shutil.copy2(CRM/'manuscript'/name,SOURCES/name)
refs=json.loads((SOURCES/'refs.json').read_text())
refs['sible2021']['cell']=refs['sible2021']['cell'].replace('J. Vis. Exp., .','J. Vis. Exp. 170, e61551.')
refs['danecek2021']['cell']=refs['danecek2021']['cell'].replace('Gigascience 10, .','Gigascience 10, giab008.')
(SOURCES/'refs.json').write_text(json.dumps(refs,indent=1,ensure_ascii=False))
# Guard that preparation has not touched any canonical source or asset.
assert all(sha(p)==h for p,h in manifest.items())
print('Isolated Codex revision prepared:',DEST)
print('Journal words:',total,'; original files unchanged:',len(manifest))
