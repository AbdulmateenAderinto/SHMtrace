"""Package vector figures, a review PDF, captions and a manuscript integration map."""
import hashlib,html,io,json,re,shutil,textwrap,zipfile
from pathlib import Path
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
from reportlab.lib.utils import simpleSplit
from reportlab.lib.pagesizes import letter
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
FIG=BASE/'figures'
ROOT=next(p for p in BASE.parents if (p/'docs/calibration/crm').is_dir())
OLD=BASE.parent/'codex_version'

# Each item represents one review page, not a newly numbered supplementary figure.
entries=[]
def add(stem,title,change,old=None):entries.append(dict(stem=stem,title=title,change=change,old=old or stem))
for n,t,c in [
 (1,'Control-informed mutation calling','Open three-stage schematic; separate Sanger and PacBio measurement routes.'),
 (2,'Matched controls distinguish candidates from background','Enlarged paired examples, full-width quality comparison, and uncluttered allele-fraction distribution.'),
 (3,'Sensitivity and false calls in spike-in benchmarks','Calling performance stays in the main figure; downstream simulations move to S25.'),
 (4,'Matched controls in real libraries','Wider small multiples, zero-call strips and heatmap keys; read accounting moves to S26.'),
 (5,'Matched Sanger and Nanopore mutation profiles','Paired mouse observations replace dominant bars; motif availability is labelled explicitly.'),
 (6,'Full-length information and linked mutation patterns','Separated clone composition and weighting; clarified reconstruction and nested node counts.'),
 (7,'External datasets across biological contexts','Burden, signature and indels are clearer; E is split into calibrated-call and raw-mismatch tracks.')]:
 add(f'Figure{n}',f'Figure {n}. {t}',c,'Figure1_preview' if n==1 else None)
ss=[
 ('FigureS1_workflow','S1. Experimental design','Experimental design separated from computational processing.','FigureS1_workflow_preview'),
 ('FigureS2','S2. Reported quality and observed control error','Aligned support panel; explicit interval and zero-mismatch definitions.',None),
 ('FigureS3_outputs','S3. Representative project outputs','Three readable outputs replace six small instruction cards.','FigureS3_outputs_preview'),
 ('FigureS4_fdr_calibration','S4. Empirical calibration against planted truth','True zero FDP retained; linear-to-log scale explained; default target marked.',None),
 ('FigureS5_controls','S5. Control number and noisy-control handling','BH fallback distinguished from decoy calibration; conditions named explicitly.',None),
 ('FigureS6_titration','S6. Dilution into AID-null reads','More room for axes; zero-call and identity-line explanations retained.',None),
 ('FigureS7_indels','S7. Indel spike-ins','Exclusion rule and planted-event counts stated explicitly.',None),
 ('FigureS8_estimates','S8. Estimation and pilot-based power (A–B)','Estimation/power separated from UMI panels under the same figure number.',None),
 ('FigureS8_umi','S8. UMI consensus (C–D)','A separate page preserves readable labels and molecule counts.','FigureS8_estimates'),
 ('FigureS9_n114','S9. N114 results: pooled view','Absolute counts replace unstable control percentages; density includes all six subclones.',None),
 ('FigureS9_spectra_jh4','S9. Per-library J_H4 spectra','Common absolute-count colour scale; zero-event spectra marked undefined.',None),
 ('FigureS9_spectra_n114','S9. Per-library N114 spectra','Common absolute-count colour scale; tiny totals explicitly qualified.',None),
 ('FigureS10','S10. Sequencing depth and discovery','Clarifies the library-specific maximum-discovery denominator.',None),
 ('FigureS11','S11. Positional background and individual reads','Separate control tracks replace overlaid profiles and small insets.',None),
 ('FigureS12_immcantation','S12. Immcantation comparison','Absolute values lead; ratios remain a secondary summary with an explicit zero floor.',None),
 ('FigureS13_sanger_matrices','S13. Sanger calls by mouse','Counts share one scale; small control totals are descriptive only.',None),
 ('FigureS13_sanger_window','S13. Readable Sanger window','Saved vector curves retained; denominator clarification added without reprocessing traces.',None),
 ('FigureS14_sanger_cohort2','S14. Maea Sanger reanalysis','Larger data panels; numerical reproduction table moved to a companion page.',None),
 ('FigureS14_comparison_table','S14. Published and reanalysed results','Same numerical table, separately readable and exported as CSV.','FigureS14_sanger_cohort2'),
 ('FigureS15_sanger_mkln','S15. Mkln1 Sanger reanalysis','Larger burden, spectrum and class plots; table separated.',None),
 ('FigureS15_comparison_table','S15. Published and reanalysed results','Same numerical table with source-study footnotes and CSV.','FigureS15_sanger_mkln'),
 ('FigureS16_swhel','S16. SWHEL sequence reanalysis','Near-unity correlation precision removed; recovery counts made explicit.',None),
 ('FigureS16_comparison_table','S16. Published and reanalysed effects','Effect comparison and differing test definitions moved to a readable table.','FigureS16_swhel'),
 ('FigureS17_spectra_full_length','S17. Spectra across nested windows','Absolute counts and an explicit qualification that windows overlap.',None),
 ('FigureS17_tree_WT2_clone0','S17. WT2 clone 0 reconstruction','Full-width tree; fit is scoped to mutation-bearing reads.','FigureS17_mutation_trees'),
 ('FigureS17_tree_WT3_clone0','S17. WT3 clone 0 reconstruction','Full-width tree; read counts and reconstructed relationships explained.','FigureS17_mutation_trees'),
 ('FigureS17_tree_WT3_clone14','S17. WT3 clone 14 reconstruction','Preserves the actual third example rather than incorrectly describing it as WT4.','FigureS17_mutation_trees'),
 ('FigureS18_engineered','S18. Engineered diversification','Local validation questions and external-list recovery remain explicit.',None),
 ('FigureS19_public_smu','S19. Public Sμ libraries','More separation between shared-window and repeat-window comparisons.',None),
 ('FigureS20_switch_pacbio','S20. PacBio switch junctions','One composition bar per biological sample replaces group-mean stacks.',None),
 ('FigureS21_platforms','S21. Platform measurements','Compact comparison table replaces dense diagrams and prose cards.',None),
 ('FigureS22_routes','S22. Distinct measurement routes','A route-comparison plate replaces the redundant workflow overview.','FigureS22_overview'),
 ('FigureS23_engine','S23. Background, testing and empirical calibration','Three readable statistical steps replace six formula cards.',None),
 ('FigureS24','S24. Artifact stress tests','Any-discovery probability distinguished from per-library FDP; interval definitions explicit.',None),
 ('FigureS25_downstream','S25. Downstream inference simulations','Original Figure 3E–G retained as A–C; original H separated into false positives (D) and power (E).','Figure3'),
 ('FigureS26_read_accounting','S26. Read accounting','Original Figure 4I moved intact to a dedicated companion figure.','Figure4')]
for s,t,c,o in ss:add(s,'Figure '+t,c,o)
# Validate assets and create predictable 600-dpi TIFF exports.
for e in entries:
 for ext in ('png','pdf','svg'):
  if not (FIG/f'{e["stem"]}.{ext}').is_file():raise FileNotFoundError(e['stem']+'.'+ext)
 p=FIG/f'{e["stem"]}.png'
 Image.open(p).convert('RGB').save(p.with_suffix('.tif'),dpi=(600,600),compression='tiff_lzw')
(BASE/'figure_manifest.json').write_text(json.dumps(entries,indent=2,ensure_ascii=False))

# Review PDF retains vector figures at their intended 174-mm width, including the tall figures.
W,H=letter;pw=174/25.4*72;writer=PdfWriter()
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=letter)
c.setTitle('SHMtrace Codex figure revision 2');c.setFont('Helvetica-Bold',20);c.drawString(55,H-85,'SHMtrace: revised figures')
c.setFont('Helvetica',12);c.drawString(55,H-112,'Codex design revision 2 · 22 October 2026')
lines=[
 'Seven main figures and the supplementary figure set, rebuilt as a separate copy.',
 'Data figures use saved results; no caller defaults or statistical tests were changed.',
 'Figure 3E–H moved to S25; Figure 4I moved to S26.',
 'S8 spans two pages. S14–S16 comparison tables and S17 trees have dedicated pages.',
 'Numbers are review identifiers. See the integration map before manuscript submission.',
 'The original Word manuscript is unchanged. Updated captions and panel maps accompany this pack.',
 'Editable SVG, Python sources, PNG, PDF and TIFF exports are in the package folder.'
]
y=H-165
for line in lines:
 for s in simpleSplit(line,'Helvetica',10,490):c.setFont('Helvetica',10);c.drawString(55,y,s);y-=16
 y-=14
c.setFont('Helvetica-Bold',11);c.drawString(55,y-10,'Review order')
y-=35
for s in ['Main figures 1–7','Supplementary figures S1–S24, including companion pages','S25: downstream simulations; S26: read accounting']:
 c.setFont('Helvetica',10);c.drawString(55,y,s);y-=20
c.save();writer.add_page(PdfReader(buf).pages[0]);writer.add_outline_item('Overview',0)
for i,e in enumerate(entries):
 src=PdfReader(FIG/f'{e["stem"]}.pdf').pages[0];sw=float(src.mediabox.width);sh=float(src.mediabox.height)
 k=pw/sw;fh=sh*k
 if fh>H-135:raise ValueError(f'Figure too tall at print width: {e["stem"]}')
 buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=letter);c.setFont('Helvetica-Bold',11)
 c.drawString((W-pw)/2,H-39,e['title'].replace('_',''))
 c.setFont('Helvetica',8);c.setFillColorRGB(.3,.34,.38)
 for j,line in enumerate(simpleSplit(e['change'],'Helvetica',8,pw)):c.drawString((W-pw)/2,52-j*11,line)
 c.setFont('Helvetica',8);c.drawRightString(W-44,22,f'{i+2}');c.save()
 page=PdfReader(buf).pages[0];page.merge_transformed_page(src,Transformation().scale(k).translate((W-pw)/2,H-68-fh),over=False)
 writer.add_page(page);writer.add_outline_item(e['title'],len(writer.pages)-1)
writer.add_metadata({'/Title':'SHMtrace Codex figure revision 2','/Author':'SHMtrace figure review'})
writer.write(BASE/'SHMtrace_Codex_Figures_v2.pdf')

# Before/after gallery for local inspection. Images and downloads stay within this folder.
parts=['<!doctype html><html><head><meta charset="utf-8"><title>SHMtrace revised figures</title><style>body{font:16px Arial;margin:32px;color:#1d2733;max-width:1500px}h1{font-size:28px}section{margin:45px 0 70px}figure{margin:0}img{width:100%;height:auto} .pair{display:grid;grid-template-columns:1fr 1fr;gap:32px}small{color:#4a5562}a{color:#1f6f8b} .links a{margin-right:16px}@media(max-width:800px){.pair{grid-template-columns:1fr}}</style></head><body><h1>SHMtrace figure revision 2</h1><p>Separate copies; original manuscript unchanged. Open the PDF for a sequential review at print width.</p><p><a href="SHMtrace_Codex_Figures_v2.pdf">Figure review PDF</a> · <a href="FIGURE_INTEGRATION.md">Panel and legend integration</a></p>']
for e in entries:
 old=BASE/'review'/f'before_{e["old"]}.png'
 parts.append(f'<section id="{e["stem"]}"><h2>{html.escape(e["title"])}</h2><p>{html.escape(e["change"])}</p><div class="pair">')
 if old.exists():parts.append(f'<figure><figcaption>Before</figcaption><img loading="lazy" src="review/{old.name}"></figure>')
 parts.append(f'<figure><figcaption>Revised</figcaption><img loading="lazy" src="figures/{e["stem"]}.png"></figure></div><p class="links">'+''.join(f'<a href="figures/{e["stem"]}.{ext}">{ext.upper()}</a>' for ext in ('svg','pdf','png','tif'))+'</p></section>')
parts.append('</body></html>');(BASE/'index.html').write_text('\n'.join(parts))

# Full original legends remain available; provide explicit replacements for changed content.
ops=json.loads((BASE/'manuscript/original_operations.json').read_text())
original=[]
for name,args,kw in ops:
 if name=='caption':original.append((args[0],args[1]))
 elif name=='figure_block':original.append((args[1],args[2]))
(BASE/'manuscript/original_legends.md').write_text('\n\n'.join('## '+t+'\n\n'+s for t,s in original))
new={
 'Figure1':('SHMtrace uses controls to interpret read-library variants',
 '(A) Matched sample and AID-null read libraries, illustrated schematically. Purple and orange marks distinguish a candidate change from recurrent background. (B) Candidate support is tested against control-derived background, with empirical decoys informing threshold selection at a 5% target. Absence or presence in controls alone does not determine a call. (C) Retained calls support burden, spectrum, positional and clone-aware summaries. (D) Measurement routes differ: Nanopore and Illumina read variants use control calibration; Sanger clone substitutions are counted at Q ≥ 20 with background reported alongside; PacBio switch reads supply junction-structure endpoints.'),
 'Figure2':('Matched controls distinguish candidates from recurrent background',
 original[1][1].replace('Shading marks 5–20% support tested in Figure 3C, rather than a universal caller boundary; dotted line, default 50% exclusion.','Dotted line, default 50% exclusion; somatic-caller comparisons are in Figure 3C.')),
 'Figure3':('Planted mutations quantify calling sensitivity and false calls',
 '(A, B) Recall by supporting-read count in held-out J_H4 and N114 controls: 60 sites × 3 draws × 3 held-out controls per locus; 90 planted events per support point. Vertical ranges span controls; annotations give observed pooled false-discovery proportions. (C) J_H4 recall at low and high mutation fractions; filled solid curves show the low-fraction spike-ins, and open dashed curves the clean high-fraction plants for permissive somatic callers. Grey shading marks the low-fraction spike-in range. (D) Recall over all planted events versus mean false calls per library, by caller colour and dataset shape. The horizontal axis is linear from 0 to 0.1 and logarithmic above 0.1; zero is not replaced by a positive value. Downstream inference simulations are in Figure S25.'),
 'Figure4':('Matched controls reduce apparent substitution background in real libraries',
 original[3][1].split('(I)')[0].strip()+' Read accounting formerly in panel I is shown in Figure S26. Grey strips identify zero-call libraries. Heatmap F cells show raw case/control ratios; colour encodes log10(max(ratio, 1)). G cells show summed control calls.'),
 'Figure5':('Matched Sanger and Nanopore libraries yield similar mutation patterns',
 'Matched 529-bp data from three Wt and three AID-null mice. (A) Within-platform unique-variant densities: Sanger and Nanopore count different units and do not share an absolute measurement scale. Original fold contrasts and tests are retained. (B) Pooled 12-type spectrum; correlation and slope are descriptive comparisons across compositional categories. (C) Four substitution classes for each Wt mouse: marker shape identifies the mouse, colour identifies platform, and lines connect its paired observations; Ts, transitions; Tv, transversions. (D) Original positional shares in 25-bp bins; the final four-base interval is shaded and placed at its actual midpoint. Correlation is descriptive. (E) Per-mouse hotspot shares, with faint mean bars; dashed lines show motif availability in the reference window. Original full experimental and statistical details are retained in manuscript/original_legends.md.'),
 'Figure6':('Long reads retain positional profiles and linked mutation patterns',
 original[5][1].replace('21 are incompatible.','21 do not fit a single reconstructed path.')+' In C, clone composition and the read-versus-clone weighting comparison are now separate adjacent plots, with unchanged values. D shows reconstructed sequence relationships rather than directly observed ancestry; node labels count nested read carriers.'),
 'Figure7':('Public datasets test analysis across biological contexts',
 original[6][1]+' In D, individual sample points are emphasized over faint mean bars. E1 shows calibrated active-AID calls; E2 shows the original Q > 20 mismatch summaries for active AID and E58Q. Shapes identify experiments 1–3 and lines connect each experiment across the Ugi, empty and UNG conditions. No observations or tests were changed.'),
 'FigureS25_downstream':('Downstream inference on the original simulations',
 '(A–C) Former Figure 3E–G: lineage recall, clone-weighted versus read-weighted estimates against truth, and estimate/truth versus depth. Existing ranges and sample-level observations are retained. (D, E) Former Figure 3H separated into false-positive rates under no group difference (D) and power under two- and fourfold differences (E), with three or five animals per group and 300 experiments per condition. The original test outcomes are unchanged.'),
 'FigureS26_read_accounting':('Read-preparation outcomes per library',
 'Former Figure 4I. One stacked bar per J_H4 or N114 library; segments give the recorded fractions of delivered reads retained as usable or removed for the stated preparation reason. Counts, categories and denominators are unchanged.'),
 'FigureS1_workflow':('Experimental design for the two in-house assay contexts',
 '(A) J_H4 germinal-centre B-cell assay, with three Wt and three AID-null mice. Reliable Nanopore and shared Sanger intervals are labelled. (B) N114 hybridoma V-region assay, with three subclones per vector following GFP-high selection and limiting dilution. Computational processing and quantitative read accounting are presented separately.'),
 'FigureS3_outputs':('Representative outputs from a J_H4 project',
 '(A) Read accounting, one trajectory per library. (B) Unique-variant density, one point per mouse. (C) Pooled substitution spectrum. This plate shows three representative outputs; positional profiles, hotspot analyses and clone summaries remain available from the software and in their dedicated manuscript figures.'),
 'FigureS11':('Position-specific background and individual-read mismatch maps',
 '(A, B) Saved Q ≥ 20 mismatch profiles for J_H4 and Sμ. AID-active means use their own scale; AID-null libraries occupy separate tracks sharing a scale within each locus. This replaces the original overlaid control profiles and insets without changing their values. (C, D) Original J_H4 read maps: rows are ordered by clone and mutation load; coloured marks lie at called sites and grey marks are other mismatches. Original between-control correlations remain reported in the manuscript.'),
 'FigureS14_sanger_cohort2':('Reanalysis of published Maea Sanger data',
 '(A) Mutation frequency per mouse, reanalysed versus published. (B) Pooled spectrum comparison. (C) G:C transitions/transversions by mouse. (D) Per-mouse agreement. (E) Mutation-load distribution with published and reanalysed clone counts. The former numerical panel F is a companion comparison table and CSV. Original numerical values, test definitions and study provenance are retained.'),
 'FigureS15_sanger_mkln':('Reanalysis of published Mkln1 Sanger data',
 '(A) Frequency per mouse. (B) Pooled spectrum comparison. (C) G:C transitions and transversions. Filled circles are reanalysed values; open squares are published values. The former numerical panel D is a companion table and CSV; its source-study qualifications remain attached.'),
 'FigureS16_swhel':('Reanalysis of published SWHEL sequences',
 '(A) Mutations per clone in each study arm. (B) Pooled 12-type spectra across five arms, shown descriptively without excessive correlation precision. (C) Original per-mouse class percentages for Wt/heterozygous and Samhd1-deficient animals. The former panel D is a companion table documenting effects and differing tests. Reanalysis recovered 2,804 of the 2,809 mutations counted by the source study across 1,795 clones and 40 animals.'),
 'FigureS21_platforms':('Platform-specific measurement units and outputs',
 'The table distinguishes Sanger clone counting, control-calibrated read-library variants, and PacBio switch-junction analysis. It states the counting unit, analysis route, output and principal qualification for each platform. Detailed matched-platform validation remains in Figures 5–6 and the corresponding supplement.'),
 'FigureS22_routes':('Distinct analysis routes matched to the measurement',
 'The plate distinguishes read-library variant calling, Sanger clone counting and switch-junction structure. It replaces a redundant software overview. Sanger background is reported alongside counts rather than subtracted; junction structure is a separate endpoint from substitution calling.'),
 'FigureS23_engine':('Control background, candidate testing and empirical threshold selection',
 '(A) Original J_H4 reported-versus-measured quality example. (B) Saved support for two example candidates in the sample and controls. (C) Retained sample candidates and control decoys over the p-value threshold, with the saved selected threshold marked. The plate focuses on the three statistical stages; detailed preprocessing is documented separately.')
}
new.update({
 'FigureS9_n114':('N114 V-region burden and pooled counts','(A) Unique-variant density in three pMX-Aicda and three pMX-Empty subclones, plotted on a linear axis so zero is retained at zero. (B, C) Pooled From-to substitution counts for the active and empty-vector groups. The active group has 292 calls and the control group one call. A control percentage spectrum based on one call is not interpreted. Each subclone was analysed at 500 reads; full assay details are retained in the original legend archive.'),
 'FigureS9_spectra_jh4':('Per-library J_H4 substitution counts','Top row: WT2, WT3 and WT4. Bottom row: KO1, KO2 and KO3. Rows of each matrix are From bases and columns are To bases. Cells give absolute unique-variant counts; all six matrices share one count colour scale. Zero-event spectra are undefined; the two-call control total is descriptive only. Libraries were analysed at 500 reads each.'),
 'FigureS9_spectra_n114':('Per-library N114 substitution counts','Top row: C1, D3 and D5 pMX-Aicda subclones. Bottom row: C3, C4 and D1 pMX-Empty subclones. Matrix cells show absolute unique-variant counts at 500 reads per library. All matrices share one count colour scale. Zero-event spectra are undefined; the single control call is descriptive only.'),
 'FigureS13_sanger_matrices':('Per-mouse Sanger substitution counts','Top row: Wt mice 2, 3 and 4. Bottom row: the three AID-null mice. Matrix rows are From bases and columns are To bases; cells show unique substitution counts at Q ≥ 20, with one count colour scale across all mice. Control totals of one or two changes are descriptive only. Full clone inclusion and exclusion accounting is retained in the original legend archive.'),
 'FigureS13_sanger_window':('Readable Sanger window and mutation density','The original saved vector coverage and density curves are retained. Coverage is the fraction of Wt clones with a Q ≥ 20 base at each reference position. Density is the number of unique changes per 1,000 observed bases within each 25-bp bin. Shading marks the shared 101–629 interval and the original alternative edge intervals. No raw chromatograms were reprocessed for this design revision.')
})
# Explicit supplemental caption corrections without discarding original experimental details.
original_by_stem={}
for name,args,kw in ops:
 if name=='figure_block':
  for path in args[0]:original_by_stem[Path(path).stem]=(args[1],args[2])
for e in entries:
 if e['stem'] in new:continue
 title,cap=original_by_stem.get(e['old'],(e['title'],e['change']))
 stem=e['stem']
 if stem.startswith('FigureS9') or stem=='FigureS13_sanger_matrices':
  cap=re.sub(r'and its share of[^.]*\.', '.',cap)
  cap+=' Revised matrices show absolute counts only; all library matrices within one plate share a count colour scale. Zero-event spectra are undefined and very small totals are descriptive only.'
 if stem=='FigureS20_switch_pacbio':cap+=' Panel B now shows each biological sample separately, using the same saved per-sample composition values instead of their group means.'
 if stem=='FigureS4_fdr_calibration':cap+=' Zero is retained at zero on a symmetric-log vertical axis, linear to 1% and logarithmic above it.'
 if stem.startswith('FigureS17_tree_'):
  cap='Selected reconstruction from the original saved trees. Nodes give nested read-carrier counts; ticks mark mutation classes. Fit percentages refer to mutation-bearing reads, not all reads in a clone. The three examples are WT2 clone 0, WT3 clone 0, and WT3 clone 14, not one example from each mouse.'
 if stem.startswith('FigureS17_spectra'):cap+=' The two windows are nested; their descriptive correlation is not an independent validation statistic.'
 if 'comparison_table' in stem:cap='Companion table for the preceding reanalysis figure. Values and qualifications are retained from its original comparison panel; an editable CSV is in sources/.'
 new[stem]=(title,cap)
# Correct stale cross-references in retained supplementary legends.
for stem,(t,c) in list(new.items()):
 c=c.replace('Figure 3E','Figure S25A').replace('Figure 3F','Figure S25B').replace('Figure 3G','Figure S25C').replace('Figure 3H','Figure S25D–E').replace('Figure 4I','Figure S26A')
 new[stem]=(t,c)
(BASE/'revised_legends.json').write_text(json.dumps(new,indent=2,ensure_ascii=False))
maptext='''# Figure revision 2: integration map

This is a separate figure package. The original v1.3 and Codex 1 Word manuscripts have not been overwritten. Use the revised legends below together with `manuscript/original_legends.md`, which preserves the full original experimental detail. Do not replace full Methods detail with the short review-PDF notes.

## Panel moves

| Before | Revised destination |
|---|---|
| Figure 3A–D | Figure 3A–D, same scientific panels |
| Figure 3E | Figure S25A |
| Figure 3F | Figure S25B |
| Figure 3G | Figure S25C |
| Figure 3H | Figure S25D (false positives) and S25E (power) |
| Figure 4I | Figure S26A |
| Figure S8A–B / C–D | Separate pages, same figure number and panel letters |
| Figure S14A | Figure S14A |
| Figure S14B (mouse agreement) | Figure S14D |
| Figure S14C (clone load) | Figure S14E |
| Figure S14D (G:C classes) | Figure S14C |
| Figure S14E (spectrum) | Figure S14B |
| Figure S14F | S14 companion comparison table / CSV |
| Figure S15B / C | Figure S15C / B |
| Figure S15D | S15 companion comparison table / CSV |
| Figure S16D | S16 companion comparison table / CSV |
| Figure S17 combined trees | Three full-width pages: WT2 clone 0; WT3 clone 0; WT3 clone 14 |
| Figure 7E | Figure 7E1 (calibrated calls) and E2 (Q > 20 mismatch summaries) |

Figures S25/S26 are provisional appended review identifiers so existing S1–S24 references stay recognizable. Before submission, renumber supplements by first citation and update all references together. If comparison tables become numbered supplementary tables, assign their final table numbers at that stage.

Figure 1, S1, S3, S21, S22 and S23 have new content arrangements: replace their legends rather than retaining old panel descriptions. Figure 5C is now paired points, not bars; its inferential test result remains in the manuscript/Methods rather than being changed. Figure 6C has adjacent composition and weighting plots. Figure S20B is sample-level composition, not a group-average stack. Matrix percentages were removed from S9/S13/S17 graphics for readability; counts and denominators remain available and percentages are derivable.

## Caption integration notes

- In Figure 2, remove any obsolete shaded somatic-caller range description; benchmarking lives in Figure 3C.
- In Figure 4, remove the old panel-I sentence. The statistical tests and geometric-mean zero convention are unchanged.
- In S3, old six-output panel references no longer apply: A is QC, B burden, C spectrum.
- In S11, replace overlay/inset descriptions with separate-track descriptions.
- In S17, the selected tree examples are not one per mouse. Do not introduce a WT4 reconstruction absent from the saved selected-tree data.
- S14–S16 companion tables preserve the original numerical comparisons; no new statistical tests were run.
- S13 coverage curves are the existing saved vector curves with clearer denominator text; raw chromatograms were not reprocessed.

## Revised legends

'''
maptext+='\n\n'.join('### '+stem+' — '+title+'\n\n'+cap for stem,(title,cap) in new.items())
(BASE/'FIGURE_INTEGRATION.md').write_text(maptext)

# Verify all originals and all recorded data inputs after all rendering.
def verify(manifest):
 d=json.loads((BASE/manifest).read_text());return {p:'missing' if not Path(p).exists() else 'changed' for p,h in d.items() if not Path(p).exists() or hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h}
validation={'original_files_checked':len(json.loads((BASE/'original_hashes.json').read_text())),
 'original_changes':verify('original_hashes.json'),'input_files_checked':len(json.loads((BASE/'input_hashes.json').read_text())),
 'input_changes':verify('input_hashes.json'),'figure_pages':len(entries),'review_pdf_pages':len(writer.pages)}
(BASE/'validation.json').write_text(json.dumps(validation,indent=2))
if validation['original_changes'] or validation['input_changes']:raise RuntimeError(validation)
print(json.dumps(validation,indent=2))
