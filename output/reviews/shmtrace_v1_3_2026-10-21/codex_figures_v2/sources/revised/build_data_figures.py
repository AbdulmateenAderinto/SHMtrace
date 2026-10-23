"""Rebuild copied data figures, preserving plotted source values and recording input hashes.

Runs the copied builders on existing saved results. The guard redirects derived JSON
writes into this package and rejects any other attempted external write.
"""
import builtins,hashlib,importlib.util,json,re,sys,traceback
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import style
from style import CASE,CASEF,CTRL,CTRLF,INK,MID,SOFT,MM,FULL,CALLERS

HERE=Path(__file__).parent
BASE=HERE.parents[1]
ROOT=next(p for p in BASE.parents if (p/'docs/calibration/crm').is_dir())
CAL=ROOT/'docs/calibration'
realopen=builtins.open
inputs={}
SNAPSHOT=json.loads((BASE/'input_snapshot_map.json').read_text()) if (BASE/'input_snapshot_map.json').exists() else {}

def guarded(file,mode='r',*a,**kw):
 if isinstance(file,(str,Path)):
  p=Path(file).resolve()
  if any(c in mode for c in 'wax+') and not p.is_relative_to(BASE):
   if p.suffix=='.json' and (p.is_relative_to(CAL) or (p.parent.name=='scratch' and p.parent.parent.name=='codex_version')):
    p=BASE/'scratch'/p.name;file=p
   else: raise RuntimeError(f'Blocked external write: {p}')
  elif p.is_relative_to(CAL) and p.suffix in ('.json','.csv','.tsv') and p.exists():
   inputs[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
  if not any(c in mode for c in 'wax+') and str(p) in SNAPSHOT:
   file=BASE/SNAPSHOT[str(p)]
 return realopen(file,mode,*a,**kw)
builtins.open=guarded

BASE.joinpath('figures').mkdir(exist_ok=True)
rawsave=style.save
# Original copied builders all call this dynamically assigned wrapper.
current=''

def clean(ax,label=None,title=None):
 for t in list(ax.texts):
  if re.fullmatch('[A-I]',t.get_text()):t.remove()
 if label:ax.text(-min(.14,.027/ax.get_position().width),1.06,label,transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
 if title:ax.set_title(title,fontsize=7.5,pad=9)

def size(fig,w=174,h=140):
 fig.set_size_inches(w*MM,h*MM)
 for t in fig.findobj(Text):
  if t.get_text() and t.get_fontsize()<6.5:t.set_fontsize(6.5)
 for ax in fig.axes:
  ax.tick_params(pad=2)

def remove_figure_legends(fig):
 for l in list(fig.legends):l.remove()

def place(ax,rect,label=None,title=None):
 ax.set_position(rect);clean(ax,label,title)

def note(fig,s,y=0.015):fig.text(0.08,y,s,fontsize=6.5,color=MID,va='bottom')

def finish(fig,name):
 size(fig,*fig.get_size_inches()/MM)
 rawsave(fig,name)
 plt.close(fig)

def f2(fig):
 size(fig,h=151)
 # A and B form the worked comparison; C spans below them.
 for a,r,l in zip(fig.axes,([.09,.66,.37,.24],[.59,.66,.37,.24],[.09,.38,.87,.16],
                          [.09,.08,.38,.21],[.59,.08,.37,.21]),'ABCDE'):
  place(a,r,l)
 a,b,c,d,e=fig.axes
 c.set_title('Base quality overlaps at both sites',fontsize=7.5,pad=8)
 # Keep the four quality groups compact in a wide plot.
 c.set_xlim(-.5,3.5)
 d.set_title('Control-derived empirical threshold',fontsize=7.5,pad=8)
 e.set_title('Allele fractions of retained calls',fontsize=7.5,pad=8)
 # The comparison explanation belongs in the figure legend, below the shared key.
 finish(fig,'Figure2')

def group_test(ax,rep,power):
 tests=[('p_pooled_reads','Pooled reads','#9a3439'),('p_welch_log','Welch, log','#274b7e'),
        ('p_negbinomial','Negative binomial','#1f6f8b'),('p_permutation','Permutation','#9aa5b1')]
 cells=[(3,2.),(3,4.),(5,2.),(5,4.)] if power else [(3,1.),(5,1.)]
 for j,(k,lab,col) in enumerate(tests):
  v=[np.mean([r[k]<.05 for r in rep if r['n']==n and r['fold']==f]) for n,f in cells]
  ax.bar(np.arange(len(cells))+(j-1.5)*.19,v,.17,color=col,label=lab)
 ax.set_ylim(0,1.3);ax.set_yticks([0,.2,.4,.6,.8,1]);ax.set_ylabel('Experiments with p < 0.05')
 ax.set_xticks(range(len(cells)),[f'{n} vs {n}\n{f:g}× difference' for n,f in cells])
 if not power:ax.axhline(.05,color=INK,lw=.6,ls=':')
 ax.legend(fontsize=6.5,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.04),columnspacing=1,handlelength=1)


def f3(fig,mode):
 size(fig,h=148 if mode=='main' else 164)
 aa=list(fig.axes);remove_figure_legends(fig)
 if mode=='main':
  # Creation order in the preserved source is A,B,D,E,F,G,H,C.
  keep=[0,1,7,2]
  for i,a in enumerate(aa):
   if i not in keep:a.remove()
  rects=[[.09,.60,.36,.29],[.59,.60,.36,.29],[.09,.13,.36,.29],[.59,.13,.36,.29]]
  titles=['J$_H$4 spike-ins','N114 V$_H$ spike-ins','Low- and high-fraction spike-ins','Sensitivity–false-call trade-off']
  for i,r,l,t in zip(keep,rects,'ABCD',titles):place(aa[i],r,l,t)
  hs=[Line2D([],[],color=CALLERS[c][1],marker=CALLERS[c][2],ms=3,lw=1,label=CALLERS[c][0])
      for c in ('recommended','geneious','lofreq','lofreq_nobaq','lofreq_permissive')]
  fig.legend(handles=hs,loc='upper center',ncol=3,bbox_to_anchor=(.52,1.0),fontsize=6.5)
  # Replace the already-corrected symlog axes label with a plain explanation.
  aa[2].set_xlabel('False calls per library')
  aa[2].text(.02,-.29,'Linear: 0–0.1; logarithmic above 0.1',transform=aa[2].transAxes,fontsize=6.5,color=MID)
  for ax in (aa[0],aa[1]):
   for t in ax.texts:
    if t.get_text().startswith('FDR:'):t.set_text(t.get_text().replace('FDR:','Observed FDP:'))
  # Existing direct labels in C are retained; move the SHMtrace label inside the axes.
  for t in aa[7].texts:
   if t.get_text()=='SHMtrace':t.set_position((2.2,.97));t.set_va('top')
  note(fig,'Simulation-based quantification and group testing are shown separately in Figure S25.')
  finish(fig,'Figure3')
 else:
  for i,a in enumerate(aa):
   if i not in (3,4,5):a.remove()
  rects=[[.09,.63,.24,.25],[.42,.63,.24,.25],[.75,.63,.23,.25]]
  for i,r,l in zip((3,4,5),rects,'ABC'):place(aa[i],r,l)
  rep=json.load(open(CAL/'replicates_calibration/runs.json'))
  ax=fig.add_axes([.09,.16,.38,.28]);group_test(ax,rep,False);clean(ax,'D','False positives: no group difference')
  ax=fig.add_axes([.59,.16,.38,.28]);group_test(ax,rep,True);clean(ax,'E','Power: true group differences')
  fig.legend(handles=[Line2D([],[],color=CALLERS[c][1],marker=CALLERS[c][2],lw=1,ms=3,label=CALLERS[c][0]) for c in ('recommended','geneious')],loc='upper center',ncol=2,bbox_to_anchor=(.5,.99),fontsize=6.5)
  # No extra error bars are invented. Existing ranges remain intact.
  note(fig,'A–C: original simulation summaries and ranges. D–E: 300 experiments per condition; animals define groups.')
  finish(fig,'FigureS25_downstream')

def zero_strip(ax):
 # Existing zero points have the same plotted floor within each dataset.
 ys=[]
 for coll in ax.collections:
  try:
   p=coll.get_paths()[0].vertices
   # a downward-triangle collection: 4 vertices including closure
   if p.shape[0]==4:
    ys.extend(coll.get_offsets()[:,1].tolist())
  except (IndexError,AttributeError):pass
 if ys:
  z=min(ys);ax.axhspan(z*.65,z*1.55,color='#edf0f3',lw=0,zorder=0)
  ax.text(.02,z,'0 calls',transform=ax.get_yaxis_transform(),fontsize=6.5,va='center',ha='left',color=MID)

def f4(fig,mode):
 aa=list(fig.axes);remove_figure_legends(fig)
 if mode=='qc':
  for a in aa[:-1]:a.remove()
  size(fig,h=95);ax=aa[-1]
  place(ax,[.09,.30,.67,.57],'A','Delivered reads: retained and removed')
  ax.legend(loc='center left',bbox_to_anchor=(1.02,.5),fontsize=6.5)
  note(fig,'One stacked bar per library; each segment is a recorded read-preparation outcome.')
  finish(fig,'FigureS26_read_accounting');return
 aa[-1].remove();size(fig,h=188)
 pos=[[.085,.755,.24,.17],[.41,.755,.24,.17],[.735,.755,.24,.17],
      [.085,.46,.24,.17],[.41,.46,.24,.17],
      [.11,.13,.22,.22],[.43,.13,.22,.22],[.75,.13,.225,.22]]
 short=['SHM','Gen.','LF','LF\n−BAQ','LF\nperm.']
 for i,a in enumerate(aa[:-1]):
  place(a,pos[i],'ABCDEFGH'[i])
  if i<5:
   a.set_xticks(range(5),short,rotation=0,ha='center');zero_strip(a)
   a.set_title(a.get_title(loc='left').replace(' (ONT)',''),fontsize=7.2,pad=8)
  if i in (5,6):
   a.set_xticks(range(5),short,rotation=0,ha='center')
   labs=[x.get_text().replace('Off-target 1','OT1').replace('Off-target 2','OT2') for x in a.get_yticklabels()]
   a.set_yticks(a.get_yticks(),labs)
   img=a.images[0]
   cbax=fig.add_axes([pos[i][0],.067,pos[i][2],.009])
   cb=fig.colorbar(img,cax=cbax,orientation='horizontal');cb.ax.tick_params(labelsize=6.5,length=2)
   cb.set_label('log$_{10}$ ratio' if i==5 else 'Control calls',fontsize=6.5,labelpad=1)
 # The original ratio heatmap colours encode log10 ratios while cells show raw ratios.
 hs=[Line2D([],[],marker='o',ls='',ms=4,mfc=CASEF,mec=CASE,label='AID-competent'),
     Line2D([],[],marker='o',ls='',ms=4,mfc=CTRLF,mec=CTRL,label='AID-null')]
 fig.legend(handles=hs,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.0),fontsize=6.5)
 fig.text(.735,.62,'Method key',fontsize=7.5,fontweight='bold',va='top')
 fig.text(.735,.56,'SHM: SHMtrace   Gen.: Geneious\nLF: LoFreq defaults\n−BAQ: BAQ disabled\nperm.: permissive settings',fontsize=6.5,color=MID,linespacing=1.6,va='top')
 fig.text(.735,.435,'Dots: individual libraries\nStars: Holm-adjusted tests\nOT1/OT2: off-targets 1/2',fontsize=6.5,color=MID,linespacing=1.6,va='top')
 note(fig,'F: geometric-mean case/control density ratio; zero-call libraries enter the ratio at one event/exposure.',.012)
 finish(fig,'Figure4')

def f5(fig):
 size(fig,h=147)
 aa=fig.axes
 rects=[[.085,.62,.24,.27],[.42,.62,.22,.27],[.085,.13,.59,.28],[.78,.13,.19,.28],[.73,.62,.24,.27]]
 # original creation order A,B,D,E,C
 for ax,r,l in zip(aa,rects,'ABDEC'):place(ax,r,l)
 aa[0].set_title('Within-platform genotype contrasts',fontsize=7.2,pad=9)
 aa[0].set_xticks([0,1,2.4,3.4],['Wt','AID-null','Wt','AID-null'])
 # C gets its own wider middle strip to make paired mouse observations readable.
 aa[4].set_position([.68,.62,.30,.27]);aa[1].set_position([.375,.62,.21,.27]);aa[0].set_position([.085,.62,.19,.27])
 aa[4].set_title('Paired spectra, three mice',fontsize=7.2,pad=9)
 aa[4].set_xticklabels(['G:C\nTs','G:C\nTv','A:T\nTs','A:T\nTv'])
 aa[0].set_xticklabels(['Wt','Null','Wt','Null'])
 for t in aa[0].texts:
  t.set_text(t.get_text().replace('Sanger clones','Sanger').replace('Nanopore reads','Nanopore'))
 aa[1].set_title('Pooled spectrum',fontsize=7.2,pad=9)
 aa[0].set_title('Genotype contrasts',fontsize=7.2,pad=9)
 aa[2].set_title('Position along the shared window',fontsize=7.2,pad=9)
 aa[3].legend(handles=[Patch(color='#a3621b',label='Sanger'),Patch(color='#1f6f8b',label='Nanopore'),
      Line2D([],[],ls='--',color=INK,lw=.7,label='Motif availability')],fontsize=6.5,loc='upper left',bbox_to_anchor=(-.1,1.05))
 note(fig,'Pattern comparisons use within-platform shares; absolute densities count different units. Lines in C pair the same mouse.')
 finish(fig,'Figure5')

def f6(fig):
 size(fig,h=150);aa=fig.axes
 for a,r,l in zip(aa,[[.085,.61,.57,.26],[.76,.61,.21,.26],[.085,.18,.22,.25],[.60,.18,.37,.25]],'ABCD'):place(a,r,l)
 a=aa[2]
 for t in list(a.texts):
  if '→' in t.get_text() or 'per 10' in t.get_text():t.remove()
 a.set_xlim(0,100);a.set_title('Junction-defined clones',fontsize=7.2,pad=9)
 b=fig.add_axes([.405,.18,.12,.25])
 ex={s['sample']:s for s in json.load(open(style.DATA/'extras_jh4/project.json'))['samples']}
 for i,s in enumerate(('WT2','WT3','WT4')):
  vals=[ex[s]['mutation_frequency_x1e3'],ex[s]['clone_weighted_frequency_x1e3']]
  b.plot([0,1],vals,color=SOFT,lw=.6)
  b.scatter([0,1],vals,s=18,marker=['o','s','^'][i],color=CASE,zorder=3,label=f'Wt {s[-1]}')
 b.set_xticks([0,1],['Read','Clone']);b.set_ylabel('Mutations per 10$^3$ bp')
 b.set_xlim(-.3,1.3);b.set_title('Weighting',fontsize=7.2,pad=9)
 b.legend(fontsize=6.5,loc='upper left',handletextpad=.2)
 aa[3].set_title('Reconstruction of one Wt3 clone',fontsize=7.2,pad=20)
 for t in aa[3].texts:
  if 'fit one path' in t.get_text():
   t.set_text('69 reads: 48 fit one path; 21 do not\nNode labels: nested read counts')
   t.set_position((.5,1.025));t.set_fontsize(6.5)
 from fig6_longreads import CLASS
 aa[3].get_legend().remove() if aa[3].get_legend() else None
 fig.legend(handles=[Line2D([],[],color=c,lw=1.4,label=lab) for lab,c in CLASS.values()],fontsize=6.5,loc='lower center',bbox_to_anchor=(.53,.078),ncol=3,handlelength=1)
 note(fig,'C: hatched reads have no assigned junction. D: ticks mark mutation classes; circle area increases with read count.',.035)
 note(fig,'Branching is a reconstruction from observed sequences, not a directly measured genealogy.',.012)
 finish(fig,'Figure6')

def f7(fig):
 size(fig,h=158);aa=list(fig.axes)
 rects=[[.09,.64,.24,.25],[.415,.64,.24,.25],[.74,.64,.23,.25],[.09,.21,.32,.24]]
 titles=['Burden: J$_H$4 mice (Illumina)','Burden: AID variants in CH12F3','Signature of the retained calls','Indels in Ramos (Illumina)']
 for a,r,l,t in zip(aa[:4],rects,'ABCD',titles):
  place(a,r,l,t)
  for tx in list(a.texts):
   if tx.get_fontstyle()=='italic':tx.remove()
 for a in aa[:2]:zero_strip(a)
 # Bars in D no longer obscure the biological observations.
 for p in aa[3].patches:p.set_alpha(.18)
 for c in aa[3].collections:
  if len(c.get_offsets()):c.set_sizes([17])
 # Rebuild E with explicit assay-specific tracks and experiment marker identity.
 aa[4].remove()
 mo=json.load(open(style.DATA/'public/morange/project.json'))
 au=json.load(open(style.DATA/'public/morange/authors_counts.json'))['libraries']
 def share(sample):
  vs=[v for v in mo['variants'] if v['sample']==sample and v['ref'] in 'CG']
  return 100*sum(v['ref']+v['alt'] in ('CT','GA') for v in vs)/len(vs) if vs else np.nan
 arms=('Ugi','Control','UNG');xx=np.arange(3)
 for k,(rect,title) in enumerate((([.53,.21,.19,.24],'E1  Control-calibrated calls'),([.80,.21,.18,.24],'E2  Q > 20 mismatches'))):
  ax=fig.add_axes(rect)
  for e,mk in zip((1,2,3),('o','s','^')):
   if k==0:
    vals=[share(f'Exp{e}_{arm}_AID') for arm in arms]
    ax.plot(xx,vals,color=CASE,lw=.6,marker=mk,ms=4,label=f'Exp {e}')
   else:
    for suffix,col,ls in (('AID',MID,'-'),('AIDE58Q',CTRL,'--')):
     vals=[au[f'Exp{e}_{arm}_{suffix}']['pct_transition_at_GC'] for arm in arms]
     ax.plot(xx,vals,color=col,lw=.6,ls=ls,marker=mk,ms=3.2)
  ax.set_ylim(-4,104);ax.set_xticks(xx,['Ugi','Empty','UNG']);ax.set_title(title,fontsize=7.2,pad=9)
  ax.set_ylabel('Transitions among C:G events (%)' if k==0 else '')
  if k==0:ax.legend(fontsize=6.5,loc='lower left',bbox_to_anchor=(-.08,-.52),ncol=3,handlelength=.7,columnspacing=.7)
  else:
   ax.legend(handles=[Line2D([],[],color=MID,lw=.6,label='Active AID'),Line2D([],[],color=CTRL,ls='--',lw=.6,label='E58Q')],
      fontsize=6.5,loc='lower left',bbox_to_anchor=(-.1,-.52),handlelength=1.2)
 note(fig,'E: Pérez-Durán non-B-cell dataset. Each shape is one experiment; lines connect its three repair conditions.',.012)
 finish(fig,'Figure7')


def supp(fig,name):
 # Generous physical size without adding decorative containers.
 w,h=fig.get_size_inches()/MM
 # All package figures are delivered at the same 174-mm width.
 # Enlarge the half-width quality plot proportionally.
 if w<100:fig.set_size_inches(FULL,h*MM)
 size(fig,w=174,h=max(h,90) if name!='FigureS8_estimates' else 130)
 if name=='FigureS2':
  fig.subplots_adjust(left=.12,right=.97,top=.91,bottom=.14)
  note(fig,'Shading: 95% binomial intervals. Triangles: upper bounds with zero mismatches. Bins require ≥2,000 bases.')
 elif name=='FigureS4_fdr_calibration':
  fig.subplots_adjust(top=.88,bottom=.25,wspace=.55)
  fig.axes[0].set_title('Observed versus target',fontsize=7.2,pad=9)
  fig.axes[1].set_title('False calls per library',fontsize=7.2,pad=9)
  fig.axes[2].set_title('Recall',fontsize=7.2,pad=9)
  for t in list(fig.axes[0].texts):t.remove()
  note(fig,'A: pooled false-discovery proportion; linear to 1%, logarithmic above. Vertical line: 5% default target.')
 elif name=='FigureS5_controls':
  fig.subplots_adjust(bottom=.25,wspace=.5)
  for a in fig.axes:
   for t in a.texts:
    if 'BH, no' in t.get_text():t.set_text('1 control:\nBH fallback')
  fig.axes[-1].set_xticklabels(['Clean','Rule\noff','Rule\non']*2)
  fig.axes[-1].set_xlabel('One noisy control: exclusion rule')
  note(fig,'One-control designs use a BH fallback; two or more controls support empirical decoy calibration.')
 elif name=='FigureS6_titration':
  fig.subplots_adjust(top=.9,bottom=.17,hspace=.75,wspace=.4)
  note(fig,'Dashed diagonal: expected = observed. Lower strips mark zero calls; curves retain the original library identities.')
 elif name=='FigureS7_indels':
  fig.subplots_adjust(bottom=.2,wspace=.42)
  note(fig,'Shaded support levels are excluded by the ≥3-read rule. Numbers beside class bars count planted indels.')
 elif name=='FigureS8_estimates':
  # A/B and C/D stay under one supplementary number, but use separate pages.
  for a in fig.axes:a.set_position([.10 if a in fig.axes[::2] else .60,.28,.35,.50])
  keep=(0,1) if current.endswith(':estimates') else (2,3)
  for i,a in enumerate(list(fig.axes)):
   if i not in keep:a.remove()
  fig.set_size_inches(FULL,91*MM)
  note(fig,'A–B: estimation and pilot-based power.' if keep==(0,1) else 'C–D: simulated and measured UMI consensus performance; n denotes molecules.')
  name='FigureS8_estimates' if keep==(0,1) else 'FigureS8_umi'
 elif name=='FigureS10':
  fig.subplots_adjust(bottom=.17,hspace=.6)
  note(fig,'Recovery in E uses each library’s maximum discovery across all depths/draws; lines retain library identity.')
 elif name=='FigureS11':
  fig.subplots_adjust(bottom=.16,hspace=.65)
  note(fig,'Top profiles share positional coordinates with the read maps; rows in C–D are ordered by clone.')
 elif name=='FigureS12_immcantation':
  fig.subplots_adjust(bottom=.35,wspace=.45)
  fig.axes[0].set_position([.12,.25,.52,.65])
  fig.axes[1].set_position([.78,.25,.19,.65])
  fig.axes[0].set_title('Absolute library-level frequencies',fontsize=7.2,pad=9)
  fig.axes[1].set_title('Secondary ratio summary',fontsize=7.2,pad=9)
  fig.axes[1].set_xticklabels(['IMGT','Parent','Q20\nmask','SHM'])
  note(fig,'Routes differ in germline reference and quality masking. Ratios use a one-event floor for zero-call libraries.')
 elif name in ('FigureS14_sanger_cohort2','FigureS15_sanger_mkln'):
  # The numerical reproduction table moves into its own editable CSV/Markdown.
  ax=fig.axes[-1]
  rows=[[t.get_text() for t in ax.texts]]
  (BASE/'scratch'/f'{name}_table_text.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
  ax.set_title('Published and reanalysed results',fontsize=7.2,pad=9)
  # Preserve table here; separate readable table is produced in custom_figures.py.
  size(fig,h=150 if name.startswith('FigureS14') else 110)
 elif name=='FigureS18_engineered':
  fig.subplots_adjust(bottom=.25,wspace=.45)
  note(fig,'A: distinct deaminases on one substrate. B: recovery of the external cluster list; marker area encodes entries.')
 elif name=='FigureS19_public_smu':
  fig.subplots_adjust(bottom=.2,hspace=.75)
  for tx in fig.axes[2].texts:
   if 'every library here is paired' in tx.get_text():
    tx.set_text('Paired reads in every library.\nWt: one library per treatment;\npooled descriptive contrast.')
    tx.set_position((.5,-.30));tx.set_ha('center')
  note(fig,'B uses the interval shared by all libraries; C uses the repeat interval in paired-read libraries only.')
 elif name=='FigureS20_switch_pacbio':
  size(fig,h=148)
  # Replace group-averaged bars by one composition bar per biological sample.
  import figS20_switch_pacbio as m
  ax=fig.axes[1];ax.clear()
  d=json.load(open(style.DATA/'switch_pacbio/switch_pacbio.json'))
  x=0;ticks=[]
  for g in m.GEN:
   ss=[v for k,v in d['samples'].items() if d['groups'][k]==g];ticks.append(x+(len(ss)-1)/2)
   for v in ss:
    bottom=0
    for cl in m.CLS:
     val=v['microhomology_classes_pct'].get(cl,0)
     ax.bar(x,val,.8,bottom=bottom,color=m.CLSC[cl],edgecolor='white',lw=.5)
     bottom+=val
    x+=1
   x+=1
  ax.set_ylim(0,100);ax.set_xticks(ticks,[m.LAB[g] for g in m.GEN]);ax.set_ylabel('Each sample’s unique junctions (%)')
  clean(ax,'B','Junction composition, sample by sample')
  ax.legend(handles=[Patch(color=m.CLSC[c],label=c.replace('>=','≥')) for c in m.CLS],fontsize=6.5,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.3))
  fig.subplots_adjust(bottom=.2,hspace=.8,wspace=.45)
  note(fig,'Each dot/bar is a sample. D measures placement ambiguity directly, separately from microhomology length.')
 elif name=='FigureS24':
  size(fig,h=151);fig.subplots_adjust(bottom=.21,hspace=.8,wspace=.45)
  note(fig,'A: probability of ≥1 call in an all-null resample (95% exact binomial intervals). B: mean false calls (upper SEM).',.045)
  note(fig,'C: per-library false / total calls with planted signal. Case-only artifacts illustrate a limit of matched-control calibration.',.012)
 if name not in ('FigureS8_estimates','FigureS8_umi'):
  for i,ax in enumerate(fig.axes):
   if i<9:clean(ax,'ABCDEFGHI'[i])
 if name in ('FigureS6_titration','FigureS11','FigureS19_public_smu','FigureS10'):
  bottom=.26 if name=='FigureS19_public_smu' else .22
  boxes=[a.get_position() for a in fig.axes];lo=min(b.y0 for b in boxes);hi=max(b.y1 for b in boxes)
  for a,b in zip(fig.axes,boxes):
   y=bottom+(b.y0-lo)*(.90-bottom)/(hi-lo);h=b.height*(.90-bottom)/(hi-lo)
   a.set_position([b.x0,y,b.width,h])
 finish(fig,name)


def dispatch(fig,name):
 if name=='Figure2':f2(fig)
 elif name=='Figure3':f3(fig,'main' if current.endswith(':main') else 'supp')
 elif name=='Figure4':f4(fig,'main' if current.endswith(':main') else 'qc')
 elif name=='Figure5':f5(fig)
 elif name=='Figure6':f6(fig)
 elif name=='Figure7':f7(fig)
 else:supp(fig,name)
style.save=dispatch

JOBS=[('fig2_call.py',''),('fig3_calibration.py','main'),('fig3_calibration.py','supp'),
      ('fig4_controls.py','main'),('fig4_controls.py','qc'),('fig5_sanger.py',''),('fig6_longreads.py',''),('fig7_platforms.py','')]
JOBS += [(p,'') for p in ('figS2_quality.py','figS4_fdr_calibration.py','figS5_controls.py','figS6_titration.py','figS7_indels.py')]
JOBS += [('figS8_estimates.py','estimates'),('figS8_estimates.py','umi')]
JOBS += [(p,'') for p in ('figS10_depth.py','figS11_background.py','figS12_immcantation.py','figS14_sanger_cohort2.py',
 'figS15_sanger_mkln.py','figS18_engineered.py','figS19_public_smu.py','figS20_switch_pacbio.py','figS24_fdr_validation.py')]

if __name__=='__main__':
 style.setup();errors=[]
 selected=sys.argv[1:]
 for script,variant in JOBS:
  if selected and script not in selected:continue
  current=f'{script}:{variant}'
  try:
   spec=importlib.util.spec_from_file_location(Path(script).stem,HERE/script)
   mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
   mod.main()
  except Exception as e:
   traceback.print_exc();errors.append([current,str(e)]);plt.close('all')
 ip=BASE/'input_hashes.json'
 previous=json.loads(ip.read_text()) if ip.exists() else {}
 previous.update(inputs);ip.write_text(json.dumps(previous,indent=2))
 (BASE/'build_errors.json').write_text(json.dumps(errors,indent=2))
 if errors:raise SystemExit(1)
