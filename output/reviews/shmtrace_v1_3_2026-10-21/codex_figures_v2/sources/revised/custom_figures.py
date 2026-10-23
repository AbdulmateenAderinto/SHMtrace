"""Open-layout schematics and saved-result validation figures for the isolated revision."""
import csv,json,re,sys
from collections import Counter
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle,Patch
from matplotlib.lines import Line2D
import style
import build_data_figures as B
from style import CASE,CASEF,CTRL,CTRLF,INK,MID,SOFT,MM,FULL
BASE=B.BASE;CAL=B.CAL;CRM=CAL/'crm';OUT=BASE/'figures'
TEAL='#1f6f8b';SANG='#a3621b';MUT='#6a4c93';ERR='#d0632a'
style.setup()

def fig(h=120):return plt.figure(figsize=(FULL,h*MM))
def blank(h=120):
 f=fig(h);a=f.add_axes([.035,.035,.93,.93]);a.set(xlim=(0,1),ylim=(0,1));a.axis('off');return f,a

def text(a,x,y,s,sz=7,col=INK,bold=False,ha='left',va='top',**kw):
 return a.text(x,y,s,fontsize=sz,color=col,fontweight='bold' if bold else 'normal',ha=ha,va=va,**kw)
def arrow(a,x,y,u,v,col=MID):a.annotate('',xy=(u,v),xytext=(x,y),arrowprops={'arrowstyle':'->','lw':.8,'color':col,'mutation_scale':8})
def pile(a,x,y,w,hits,n=6):
 for r in range(n):a.plot([x,x+w],[y-r*.025]*2,color='#bdc8d1',lw=1.4)
 for f,rows,col in hits:
  for r in rows:a.plot(x+w*f,y-r*.025,'o',ms=3,color=col)

def data(name):
 p=BASE/'sources/originals'/f'{name}.js'
 d=json.loads(re.search(r'/\*DATA\*/(.*?)/\*END\*/',p.read_text(),re.S)[1])
 (BASE/'sources'/f'{name}_data.json').write_text(json.dumps(d,indent=2,ensure_ascii=False))
 return d

def save(f,name):B.finish(f,name)

def figure1():
 f,a=blank(133)
 text(a,.0,.99,'A  Matched inputs',8,bold=True)
 text(a,.34,.99,'B  Control-informed calibration',8,bold=True,col=TEAL)
 text(a,.80,.99,'C  Interpretable calls',8,bold=True)
 text(a,.0,.91,'Sample reads',7,bold=True,col=CASE)
 pile(a,.01,.85,.23,[(.3,[0,2,4],MUT),(.78,[1,3,5],ERR)])
 text(a,.0,.65,'AID-null control reads',7,bold=True,col=CTRL)
 pile(a,.01,.59,.23,[(.78,[0,2,4],ERR)])
 text(a,.0,.39,'Matched preparation and run\nDepth-matched usable reads',6.5,col=MID)
 arrow(a,.26,.72,.32,.72)
 text(a,.35,.90,'Observed change versus\ncontrol-derived background',7.5,bold=True)
 text(a,.35,.78,'Candidate absent from controls',7,col=MUT)
 text(a,.35,.73,'Evidence can exceed background',6.5,col=MID)
 text(a,.35,.64,'Candidate recurrent in controls',7,col=ERR)
 text(a,.35,.59,'Evidence can match background',6.5,col=MID)
 text(a,.35,.48,'Candidate tests → empirical decoys',7,bold=True,col=TEAL)
 text(a,.35,.42,'Controls set the calling threshold\nEmpirical false-discovery target: 5%',6.5,col=MID)
 arrow(a,.72,.72,.78,.72)
 text(a,.80,.89,'Retained variants',7.5,bold=True)
 text(a,.80,.81,'Mutation burden\nSubstitution spectrum\nPosition and hotspots\nClones and reconstruction',7,col=MID,linespacing=1.8)
 text(a,.0,.28,'D  Platform-specific measurement routes',8,bold=True)
 text(a,.0,.20,'Nanopore / Illumina',7,bold=True,col=TEAL)
 text(a,.0,.145,'Control-calibrated variants\nOptional UMI consensus',6.5,col=MID)
 text(a,.38,.20,'Sanger clones',7,bold=True,col=SANG)
 text(a,.38,.145,'Count changes at Q ≥ 20\nReport background alongside',6.5,col=MID)
 text(a,.76,.20,'PacBio switch reads',7,bold=True,col=CTRL)
 text(a,.76,.145,'Junction structure\nMicrohomology and insertions',6.5,col=MID)
 text(a,.0,.025,'Schematic pileups illustrate the decision logic; colours distinguish candidate changes and recurrent background.',6.5,col=MID)
 save(f,'Figure1')

def s1():
 f,a=blank(111)
 text(a,0,.99,'A  J$_H$4: germinal-centre B cells from mice',8,bold=True)
 text(a,.01,.86,'Wild type / AID-null\n3 mice per group',7,col=CASE)
 arrow(a,.24,.83,.30,.83)
 text(a,.32,.86,'Peyer’s patch B cells\nB220$^+$ GL7$^+$ Fas$^+$',7)
 arrow(a,.59,.83,.65,.83)
 text(a,.67,.86,'High-fidelity PCR\nJ$_H$4 intron amplicon',7)
 a.plot([.05,.94],[.67,.67],color=TEAL,lw=1.4)
 text(a,.05,.63,'Nanopore: reliable interval 90–1,230 (1,141 bp)',7,col=TEAL)
 a.plot([.06,.48],[.54,.54],color=SANG,lw=1.4)
 text(a,.05,.50,'Sanger comparison: 101–629 (529 bp)',7,col=SANG)
 text(a,0,.38,'B  N114: fixed V-region substrate in hybridoma subclones',8,bold=True)
 text(a,.01,.27,'pMX-Aicda / empty vector\nGFP and puromycin',7,col=CASE)
 arrow(a,.25,.23,.31,.23)
 text(a,.33,.27,'GFP-high selection\nLimiting dilution',7)
 arrow(a,.57,.23,.63,.23)
 text(a,.65,.27,'3 subclones per vector\n1 month before DNA collection',7)
 text(a,.01,.10,'Nanopore V-region amplicon: ~655 bp; analysed interval: 21–635 (615 bp)',7,col=TEAL)
 text(a,.01,.025,'Computational processing is shown in Figure 1; quantitative read accounting is in Figure S26.',6.5,col=MID)
 save(f,'FigureS1_workflow')

def s3():
 d=data('figS3_outputs');f=fig(119)
 a=f.add_axes([.09,.33,.24,.49]);B.clean(a,'A','Read preparation')
 for n,v in d['qc'].items():a.plot(range(len(v)),v,color=CASE if n.startswith('WT') else CTRL,lw=.7,alpha=.7)
 a.set_xticks(range(5),['Input','Quality','Target','Usable','Draw'],rotation=30,ha='right');a.set_ylabel('Delivered reads retained (%)')
 b=f.add_axes([.42,.33,.23,.49]);B.clean(b,'B','Per-animal mutation burden')
 for name,grp,val,one in d['dens']['rows']:
  x=0 if name.startswith('WT') else 1
  idx=int(name[-1])-2 if x==0 else int(name[-1])-1
  b.scatter(x+(idx-1)*.08,val,s=20,color=CASE if x==0 else CTRL,zorder=3)
 b.set_xticks([0,1],['Wt','AID-null']);b.set_ylabel('Unique variants per 10$^5$ bases')
 c=f.add_axes([.75,.33,.23,.49]);B.clean(c,'C','Substitution spectrum')
 c.bar(range(12),d['spec']['group'],color=[SANG]*6+[CASE]*6,width=.8)
 c.set_xticks(range(12),[s.replace('>','→') for s in d['spec']['classes']],rotation=90);c.set_ylabel('Share of calls (%)')
 B.note(f,'Example outputs from the same J$_H$4 project: one line per library in A, one point per mouse in B.',.10)
 B.note(f,'The report also includes positional profiles, hotspot use, clone summaries, and downloadable tables.',.055)
 save(f,'FigureS3_outputs')

BASES='ATCG'
def matrix(ax,counts,title,vmax=None):
 total=sum(counts.values());arr=np.full((4,4),np.nan)
 for i,x in enumerate(BASES):
  for j,y in enumerate(BASES):
   if i!=j:arr[i,j]=counts.get(x+'>'+y,0)
 masked=np.ma.masked_invalid(arr)
 ax.imshow(masked,cmap='Blues',vmin=0,vmax=vmax or max(1,np.nanmax(arr)))
 for i in range(4):
  for j in range(4):
   if i!=j:
    v=arr[i,j];ax.text(j,i,str(int(v)),ha='center',va='center',fontsize=7,
     color='white' if v>(vmax or np.nanmax(arr))*.55 and v else INK)
 ax.set_xticks(range(4),list(BASES));ax.set_yticks(range(4),list(BASES));ax.tick_params(length=0)
 ax.set_xlabel('To');ax.set_ylabel('From');ax.set_title(f'{title}\n{total:,} calls',fontsize=7.2,pad=9)
 for s in ax.spines.values():s.set_visible(False)
 if total==0:ax.text(.5,-.36,'No calls; spectrum undefined',transform=ax.transAxes,fontsize=6.5,ha='center',color=MID)
 elif total<5:ax.text(.5,-.36,'Descriptive only: very small count',transform=ax.transAxes,fontsize=6.5,ha='center',color=MID)

def counts(vs):return Counter(v['ref']+'>'+v['alt'] for v in vs)

def matrices():
 for key,names,title in [('jh4',['WT2','WT3','WT4','KO1','KO2','KO3'],'J$_H$4'),
                         ('n114',['C1','D3','D5','C3','C4','D1'],'N114 V$_H$')]:
  p=json.load(open(CAL/f'figures/{key}_shmtrace/project.json'))
  f=fig(139);cc=[counts([v for v in p['variants'] if v['sample']==n]) for n in names]
  # Common absolute-count scale within this figure; no separate 0–100% control scale.
  vmax=max(max(c.values(),default=0) for c in cc)
  for i,(n,c) in enumerate(zip(names,cc)):
   ax=f.add_axes([.09+(i%3)*.32,.63-(i//3)*.44,.20,.24]);matrix(ax,c,n,vmax)
  B.note(f,f'{title}: absolute counts, one matrix per library. All panels use the same colour scale; diagonal cells are not substitutions.')
  save(f,f'FigureS9_spectra_{key}')
 # Pooled N114 and its burden comparison.
 p=json.load(open(CAL/'figures/n114_shmtrace/project.json'));f=fig(102)
 a=f.add_axes([.10,.30,.25,.51]);B.clean(a,'A','Density per subclone')
 for g,col in [('pMX-Aicda',CASE),('pMX-Empty',CTRL)]:
  ss=[s for s in p['samples'] if s['group']==g];x=0 if g=='pMX-Aicda' else 1
  a.scatter(x+np.linspace(-.1,.1,len(ss)),[s['unique_density_x1e5'] for s in ss],s=22,color=col)
 a.set_xticks([0,1],['pMX-Aicda','Empty vector']);a.set_ylabel('Unique variants per 10$^5$ bases');a.set_ylim(bottom=-1)
 for k,(g,l) in enumerate([('pMX-Aicda','B'),('pMX-Empty','C')]):
  ax=f.add_axes([.47+k*.29,.30,.19,.48]);matrix(ax,counts([v for v in p['variants'] if v['group']==g]),'pMX-Aicda' if g=='pMX-Aicda' else 'Empty vector')
  B.clean(ax,l)
 B.note(f,'Matrices show absolute unique-variant counts. Empty-vector calls are too few for a stable percentage spectrum.')
 save(f,'FigureS9_n114')
 # Sanger per-mouse matrices, original corrected clone calls.
 d=json.load(open(CAL/'manuscript_corrected/sanger_calls/corrected_sanger.json'))
 names=[n for n in d if n!='source'];cc=[Counter(v[1]+'>'+v[2] for v in d[n]['variants']) for n in names]
 vmax=max(max(c.values(),default=0) for c in cc);f=fig(140)
 for i,(n,c) in enumerate(zip(names,cc)):
  ax=f.add_axes([.09+(i%3)*.32,.63-(i//3)*.44,.20,.24]);matrix(ax,c,n,vmax)
 B.note(f,'Sanger calls at Q ≥ 20. Absolute counts share one colour scale; tiny control totals are not interpreted as spectra.')
 save(f,'FigureS13_sanger_matrices')


def s17():
 p=json.load(open(CAL/'figures/jh4_shmtrace/project.json'))
 q=json.load(open(CAL/'figures/jh4_shmtrace_full/project.json'))
 c1=counts([v for v in p['variants'] if v['sample'].startswith('WT')]);c2=counts([v for v in q['variants'] if v['sample'].startswith('WT')])
 f=fig(108);vmax=max(max(c1.values()),max(c2.values()))
 for i,(c,title) in enumerate([(c1,'Shared 529-bp window'),(c2,'Reliable 1,141-bp interval')]):
  a=f.add_axes([.09+i*.31,.34,.20,.43]);matrix(a,c,title,vmax);B.clean(a,'AB'[i])
 a=f.add_axes([.73,.32,.24,.49]);B.clean(a,'C','Spectrum across nested windows')
 ts=[x+'>'+y for x in BASES for y in BASES if x!=y]
 x=[100*c1[t]/sum(c1.values()) for t in ts];y=[100*c2[t]/sum(c2.values()) for t in ts]
 a.scatter(x,y,s=20,color=MUT);a.plot([0,22],[0,22],ls='--',lw=.6,color=SOFT)
 a.set(xlim=(0,22),ylim=(0,22),xlabel='529-bp calls (%)',ylabel='1,141-bp calls (%)')
 B.note(f,'The shared window is contained within the full interval: this is a descriptive comparison of nested data, not independent validation.',.08)
 save(f,'FigureS17_spectra_full_length')
 # Separate, full-width tree for each mouse (retains the original selected clone).
 from fig6_longreads import draw_tree,CLASS
 t=json.load(open(CAL/'figures/mutation_trees/trees.json'))
 for j,(n,cl) in enumerate((('WT2',0),('WT3',0),('WT3',14))):
  nodes=t['wt3_figure_clone']['nodes'] if (n,cl)==('WT3',0) else t['trees'][f'{n}:{cl}'];meta=next(v for v in t['clones'][n] if v['clone']==cl)
  f=fig(102);a=f.add_axes([.08,.28,.88,.51]);draw_tree(a,nodes)
  B.clean(a,'ABC'[j],f'{n}: clone {cl}, {meta["reads"]} reads; {100*meta["fit"]:.0f}% of mutation-bearing reads fit')
  hs=[Line2D([],[],color=c,lw=1.4,label=lab) for lab,c in CLASS.values()]
  f.legend(handles=hs,loc='lower center',bbox_to_anchor=(.5,.12),ncol=3,fontsize=6.5)
  B.note(f,'Node labels are nested read counts; ticks are mutations. Branches represent reconstructed sequence relationships, not observed ancestry.')
  save(f,f'FigureS17_tree_{n}_clone{cl}')


def external_sanger():
 # Reuse exact injected plot values; no reanalysis of external chromatograms.
 for name,n,ko in [('figS14_cohort2',14,'Maea_KO'),('figS15_mkln',15,'Mkln1_KO')]:
  d=data(name);f=fig(137 if n==14 else 105)
  groups=['WT',ko,'Aicda_KO'];cols={groups[0]:CASE,ko:MUT,'Aicda_KO':CTRL}
  a=f.add_axes([.10,.60 if n==14 else .27,.23,.29 if n==14 else .52]);B.clean(a,'A','Frequency per mouse')
  for i,g in enumerate(groups):
   v=d['A'][g]
   for key,dx,mk,face in [('sh',-.15,'o',cols[g]),('pub',.15,'s','white')]:
    vs=v[key]
    a.scatter(i+dx+np.linspace(-.05,.05,len(vs)),vs,s=22,marker=mk,facecolor=face,edgecolor=cols[g],lw=.6)
  a.set_xticks(range(3),['Wt','Maea KO' if n==14 else 'Mkln1 KO','AID-null']);a.set_ylabel('Mutations per 10$^4$ bp');a.set_ylim(bottom=-3)
  a.legend(handles=[Line2D([],[],marker='o',ls='',color=CASE,label='SHMtrace'),Line2D([],[],marker='s',ls='',mfc='white',mec=INK,label='Published')],fontsize=6.5,loc='upper right')
  a=f.add_axes([.43,.60 if n==14 else .27,.23,.29 if n==14 else .52]);B.clean(a,'B','Spectrum comparison')
  specs=d['E'] if n==14 else d['C']
  for g,rows in specs.items():a.scatter([r[0] for r in rows],[r[1] for r in rows],s=19,color=cols[g],label='Wt' if g=='WT' else 'KO')
  a.plot([0,40],[0,40],ls='--',lw=.6,color=SOFT);a.set(xlabel='Published (%)',ylabel='SHMtrace (%)',xlim=(0,40),ylim=(0,40));a.legend(fontsize=6.5)
  a=f.add_axes([.76,.60 if n==14 else .27,.21,.29 if n==14 else .52]);B.clean(a,'C','G:C mutation classes')
  rows=d['D'] if n==14 else d['B']
  for i,r in enumerate(rows):
   for key,dx,mk,face in [('sh',-.14,'o',cols[r['group']]),('pub',.14,'s','white')]:
    a.scatter(i+dx+np.linspace(-.04,.04,len(r[key])),r[key],s=16,marker=mk,facecolor=face,edgecolor=cols[r['group']],lw=.6)
  a.set_xticks(range(len(rows)),[r['kind']+'\n'+('Wt' if r['group']=='WT' else 'KO') for r in rows]);a.set_ylabel('Share of all mutations (%)')
  if n==14:
   a=f.add_axes([.10,.18,.23,.24]);B.clean(a,'D','Agreement by mouse')
   for x,y,g in d['B']['points']:a.scatter(x,y,s=20,color=cols[g])
   a.plot([8,400],[8,400],ls='--',lw=.6,color=SOFT);a.set(xscale='log',yscale='log',xlabel='Published /10$^4$ bp',ylabel='SHMtrace /10$^4$ bp')
   a=f.add_axes([.60,.18,.37,.24]);B.clean(a,'E','Mutation load per clone')
   shades=plt.get_cmap('Blues')(np.linspace(.1,.9,6))
   for i,r in enumerate(d['C']['rows']):
    left=0
    for val,col in zip(r['pct'],shades):a.barh(i,val,left=left,height=.65,color=col,edgecolor='white',lw=.5);left+=val
   a.set_yticks(range(4),[('Wt' if r['group']=='WT' else 'Maea KO')+' '+r['who']+f' (n={r["n"]})' for r in d['C']['rows']]);a.set_xlabel('Clones (%)')
   a.legend(handles=[Patch(color=c,label=b) for b,c in zip(d['C']['bins'],shades)],fontsize=6.5,ncol=6,loc='upper center',bbox_to_anchor=(.5,1.32),handlelength=.8,columnspacing=.6)
  B.note(f,'Filled circles: reanalysed; open squares: published. Detailed test reproduction and provenance are in the companion table.')
  save(f,f'FigureS{n}_sanger_'+('cohort2' if n==14 else 'mkln'))
  table=d['F'] if n==14 else d['D']
  foot=('Published and reanalysed tests can differ in included sequences or definitions; see STAR Methods and the source-study notes.' if n==14 else 'a: Ts/Tv as drawn in the paper’s Fig. S2C. b: published column computed from its tables. The KO spectrum comparison uses the published KO table traceable to folder K2.')
  table_figure(table,f'FigureS{n}_comparison_table',f'Figure S{n}: published and reanalysed results',foot)
 # SWHEL
 d=data('figS16_swhel');f=fig(121)
 cols=dict(zip([r['group'] for r in d['A']],[CASE,MUT,TEAL,SANG,MID]))
 a=f.add_axes([.16,.33,.23,.48]);B.clean(a,'A','Mutation load')
 for i,r in enumerate(d['A']):
  a.scatter(r['sh'],i+.12,s=23,color=cols[r['group']]);a.scatter(r['pub'],i-.12,s=23,marker='s',facecolor='white',edgecolor=cols[r['group']],lw=.7)
 a.set_yticks(range(5),[r['label'].replace('+',' + ') for r in d['A']]);a.invert_yaxis();a.set_xlabel('Mutations per clone')
 a=f.add_axes([.50,.33,.20,.48]);B.clean(a,'B','Substitution spectrum')
 for x,y,g,t in d['B']:a.scatter(x,y,s=14,color=cols[g])
 a.plot([0,50],[0,50],ls='--',lw=.6,color=SOFT);a.set(xlim=(0,50),ylim=(0,50),xlabel='Published (%)',ylabel='SHMtrace (%)')
 a=f.add_axes([.80,.33,.18,.48]);B.clean(a,'C','Mutation classes')
 for i,r in enumerate(d['C']):
  for vals,off,col in [(r['wt'],-.15,CASE),(r['ko'],.15,MUT)]:a.scatter(i+off+np.linspace(-.04,.04,len(vals)),vals,s=13,color=col)
 a.set_xticks(range(4),[r['label'].replace(' ','\n') for r in d['C']]);a.set_ylabel('Each mouse’s mutations (%)')
 B.note(f,'2,804 of 2,809 published mutations recovered across 1,795 clones and 40 animals; percentages are descriptive summaries.',.13)
 B.note(f,'Filled circles: SHMtrace; open squares: published. Detailed effects and differing tests are in the companion table.',.075)
 save(f,'FigureS16_swhel')
 table_figure(d['D'],'FigureS16_comparison_table','Figure S16: published and reanalysed effects',
 'The published total-mutation comparison used one-way ANOVA; SHMtrace uses Welch on log frequency. Other tests retain their original definitions.')


def table_figure(rows,name,title,foot='Published and reanalysed tests can differ in included sequences or definitions; see STAR Methods and the source-study notes.'):
 rows=[r[:3] for r in rows]
 with open(BASE/'sources'/f'{name}.csv','w',newline='') as fh:
  w=csv.writer(fh);w.writerow(['Measure','Published','SHMtrace']);w.writerows(rows)
 f=fig(80 if len(rows)<7 else 108);a=f.add_axes([.05,.22,.90,.63]);a.axis('off')
 table=a.table(cellText=rows,colLabels=['Measure','Published','SHMtrace'],colWidths=[.54,.23,.23],cellLoc='left',loc='center')
 table.auto_set_font_size(False);table.set_fontsize(7);table.scale(1,1.65)
 for (r,c),cell in table.get_celld().items():
  cell.set_linewidth(.5);cell.set_edgecolor('#bdc8d1')
  if r==0:cell.set_facecolor('#edf0f3');cell.set_text_props(weight='bold')
 f.text(.05,.94,title,fontsize=8,fontweight='bold')
 import textwrap
 f.text(.05,.10,'\n'.join(textwrap.wrap(foot,115)),fontsize=6.5,color=MID)
 save(f,name)


def s21():
 f,a=blank(125);text(a,0,.99,'What each platform measures',8,bold=True)
 rows=[('Sanger','Cloned amplicon','Q ≥ 20 changes counted;\nbackground reported alongside','Per-clone load and spectrum;\nshared readable window'),
       ('Nanopore','Amplicon read','Control-calibrated variants;\noptional UMI consensus','Long positional span;\njunction-defined clones'),
       ('Illumina','Merged paired read','Control-calibrated variants;\nreference and layout stated','Mutation burden and signature;\nspan limited by read layout'),
       ('PacBio','Switch-junction read','Junction alignment and\nunique-junction summaries','Microhomology, insertions;\nplacement ambiguity checked')]
 xs=[0,.18,.39,.71];heads=['Platform','Measurement unit','Analysis route','Output / qualification']
 for x,s in zip(xs,heads):text(a,x,.87,s,7,bold=True)
 for i,r in enumerate(rows):
  y=.75-i*.17
  for j,(x,s) in enumerate(zip(xs,r)):text(a,x,y,s,6.7,col=TEAL if j==0 else INK,bold=j==0,linespacing=1.7)
  a.plot([0,.99],[y-.125]*2,color='#bdc8d1',lw=.5)
 text(a,0,.035,'Measurement units and counting rules differ across routes; absolute densities should not be treated as interchangeable.',6.5,col=MID)
 save(f,'FigureS21_platforms')

def s22():
 f,a=blank(114);text(a,0,.99,'Distinct analysis routes, matched to the measurement',8,bold=True)
 for y,title,route,out,col in [(.78,'Amplicon reads','Control-derived background → candidate tests → empirical threshold','Retained variants; per-animal summaries',TEAL),
  (.50,'Sanger clones','Quality-trimmed bases → Q ≥ 20 changes counted','Per-clone changes; control background alongside',SANG),
  (.22,'PacBio switch reads','Junction alignment → unique junctions → structural summaries','Microhomology and insertion phenotypes',CTRL)]:
  text(a,0,y,title,7.5,col=col,bold=True)
  text(a,.27,y,route,7)
  text(a,.27,y-.09,out,6.5,col=MID)
 text(a,0,.035,'Control subtraction is not applied to Sanger clone counts. Switch-junction structure is a separate endpoint from substitution calling.',6.5,col=MID)
 save(f,'FigureS22_routes')

def s23():
 d=data('figS23_engine');f=fig(141)
 a=f.add_axes([.09,.54,.23,.30]);B.clean(a,'A','Measured control background')
 a.plot(d['p3']['x'],d['p3']['y'],color=CTRL,lw=1,marker='o',ms=2)
 a.plot([0,50],[0,50],ls='--',color=SOFT,lw=.6);a.set(xlabel='Reported quality',ylabel='Empirical quality')
 a=f.add_axes([.43,.54,.22,.30]);B.clean(a,'B','Candidate evidence')
 for i,r in enumerate(d['p4']):
  a.bar(i-.18,r['case'],.32,color=CASE);a.bar(i+.18,r['ctl'],.32,color=CTRL)
 a.set_xticks(range(2),['\n'.join(r['label']) for r in d['p4']]);a.set_ylabel('Reads with candidate change')
 a.legend(handles=[Patch(color=CASE,label='Sample'),Patch(color=CTRL,label='Controls')],fontsize=6.5,loc='upper right')
 a=f.add_axes([.75,.54,.23,.30]);B.clean(a,'C','Empirical threshold')
 x=10**np.array(d['p5']['grid']);a.plot(x,d['p5']['case'],color=CASE,lw=1,label='Sample')
 for v in d['p5']['decoys']:a.plot(x,v,color=CTRL,lw=.6,alpha=.5)
 a.axvline(d['p5']['t'],ls='--',color=INK,lw=.6);a.set(xscale='log',yscale='symlog',xlabel='Candidate p-value threshold',ylabel='Calls retained',ylim=(0,None))
 a.legend(handles=[Line2D([],[],color=CASE,lw=1,label='Sample'),Line2D([],[],color=CTRL,lw=.6,label='Decoys')],fontsize=6.5,loc='upper left')
 f.text(.085,.31,'Estimate background',fontsize=7.5,fontweight='bold',color=TEAL)
 f.text(.085,.25,'Pool by substitution, strand, and quality;\nshrink position-level estimates toward pooled rates.',fontsize=6.7,color=MID,linespacing=1.8)
 f.text(.43,.31,'Test candidates',fontsize=7.5,fontweight='bold',color=TEAL)
 f.text(.43,.25,'Compare observed support with\ncontrol-derived expected support.',fontsize=6.7,color=MID,linespacing=1.8)
 f.text(.75,.31,'Calibrate selection',fontsize=7.5,fontweight='bold',color=TEAL)
 f.text(.75,.25,'Held-out control calls supply decoys;\nselect an empirical threshold.',fontsize=6.7,color=MID,linespacing=1.8)
 B.note(f,'J$_H$4 example from the saved statistical-engine figure data. Preprocessing is documented separately in Figures 1 and S26.')
 save(f,'FigureS23_engine')


def background():
 from figS11_background import rates,raster_panel
 f=fig(187)
 for j,(key,title) in enumerate((('jh4','J$_H$4'),('su_lps','Sμ'))):
  d=json.load(open(style.DATA/f'intuition_{key}.json'));lo,hi=d['window'];xs=np.arange(lo,hi+1)
  arrays=[100*np.nanmean([rates(d,n,lo,hi) for n in d['cases']],axis=0)]
  arrays += [100*rates(d,n,lo,hi) for n in d['controls']]
  names=['AID-active mean']+list(d['controls'])
  ctrlmax=max(np.nanmax(v) for v in arrays[1:])*1.08
  for k,(vals,n) in enumerate(zip(arrays,names)):
   a=f.add_axes([.17+j*.49,.825-k*.082,.30,.065])
   a.plot(xs,vals,color=CASE if k==0 else CTRL,lw=.6)
   a.set_xlim(lo,hi);a.set_ylim(0,np.nanmax(vals)*1.08 if k==0 else ctrlmax)
   a.set_ylabel(n,fontsize=6.5,rotation=0,ha='right',va='center',labelpad=7)
   a.set_yticks([0,round((np.nanmax(vals) if k==0 else ctrlmax)*.8,1)])
   a.tick_params(labelsize=6.5)
   if k<3:a.set_xticks([])
   else:a.set_xlabel(f'{title} position (bp)',fontsize=6.5,labelpad=2)
   if k==0:B.clean(a,'AB'[j],f'{title}: mismatch rate at Q ≥ 20 (%)')
 d=json.load(open(style.DATA/'intuition_jh4.json'))
 for j,(n,label) in enumerate((('WT3','wild type'),('KO1','Aicda$^{-/-}$'))):
  a=f.add_axes([.10+j*.49,.16,.37,.26]);raster_panel(a,d,n,label);B.clean(a,'CD'[j])
 B.note(f,'A–B: controls have separate tracks with a shared scale within each locus; the AID-active mean uses its own scale.',.05)
 B.note(f,'C–D: individual J$_H$4 reads, ordered by clone and mutation load. Colours mark called sites; grey marks other mismatches.',.025)
 save(f,'FigureS11')

if __name__=='__main__':
 jobs=[figure1,s1,s3,matrices,s17,external_sanger,s21,s22,s23,background]
 selected=sys.argv[1:]
 for fn in jobs:
  if not selected or fn.__name__ in selected:fn()
 # Merge provenance rather than replacing records from data-figure builders.
 p=BASE/'input_hashes.json';old=json.loads(p.read_text()) if p.exists() else {};old.update(B.inputs);p.write_text(json.dumps(old,indent=2))
