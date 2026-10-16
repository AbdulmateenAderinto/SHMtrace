() => {
  // Cell Reports Methods, Figure S12: an overview of SHMtrace, from libraries to group tests. Built with
  // manuscript_corrected/figure1_scicanvas/build.js into crm/figures/FigureS12_overview; Figure S13 gives the
  // rule and equation behind each step. Cell Press rules as in crm_figure1.js (174 mm; Arial 6-8 pt;
  // strokes 0.5-1.5 pt). Mutations purple, sequencing errors orange-red, AID-null controls red,
  // AID-competent animals blue, SHMtrace teal.
  const W = 1400, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8', FAINT = '#cfd8df';
  const TEAL = '#1f6f8b', CTRL = '#9a3439', CASE = '#274b7e', SANG = '#a3621b', ILL = '#3d5a80';
  const MUT = '#6a4c93', ERR = '#d0632a';
  const F = 'arial', TITLE = 22.5, HEAD = 20.5, SM = 17.5;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: SM, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: MID, width: 2.2, headSize: 0.7, ...e }); o.push(c); return c; };
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const lines = (arr, x, y, e = {}, dy = 22.5) => arr.map((s, k) => T(s, x, y + k * dy, { color: MID, ...e }));
  const head = (x0, x1, title, sub) => {
    TC(title, (x0 + x1) / 2, 6, { fontSize: TITLE, bold: true });
    TC(sub, (x0 + x1) / 2, 36, { color: SOFT, italic: true });
    R(x0, 66, x1 - x0, 2.6, { fill: INK });
  };
  // a pileup: reads as grey bars, mutations or errors as dots
  const pile = (x, y, w, n, dots) => {
    for (let r = 0; r < n; r++) R(x, y + r * 15, w, 7, { fill: FAINT, radius: 3.5 });
    dots.forEach(([f, r, c]) => E(x + f * w, y + r * 15 + 3.5, 10, { fill: c }));
  };
  const cardTitle = (x, y, text, colour = INK) => T(text, x + 16, y + 12, { fontSize: HEAD, bold: true, color: colour });

  // ───────── column 1: libraries ─────────
  const c1 = 20, c1w = 240;
  head(c1, c1 + c1w, 'LIBRARIES', 'nanopore, Illumina or Sanger');
  TC('AID-competent animals', c1 + c1w / 2, 86, { bold: true, color: CASE });
  pile(c1 + 20, 118, c1w - 40, 8, [[0.2, 0, MUT], [0.55, 0, MUT], [0.35, 1, MUT], [0.8, 2, MUT], [0.2, 3, MUT],
    [0.62, 4, MUT], [0.45, 5, MUT], [0.7, 5, ERR], [0.3, 6, MUT], [0.85, 7, MUT]]);
  TC('AID-null controls (two or more)', c1 + c1w / 2, 254, { bold: true, color: CTRL });
  pile(c1 + 20, 286, c1w - 40, 6, [[0.7, 1, ERR], [0.7, 4, ERR], [0.35, 3, ERR]]);
  const plat = [['s-nanopore-sequencing-strand-through-pore', 'nanopore', 56], ['s-flow-cell-clusters-illumina', 'Illumina', 70],
    ['s-sanger-chromatogram-trace-only', 'Sanger', 56]];
  plat.forEach(([id, lab, sz], k) => {
    const cx = c1 + 40 + k * 80;
    I(id, cx, 420, sz);
    TC(lab, cx, 456, { color: k === 2 ? SANG : (k === 1 ? ILL : INK) });
  });
  const sheet = card(c1, 500, c1w, 230);
  cardTitle(c1, 500, 'Sample sheet');
  lines(['Each file, or each Sanger', 'clone, is assigned to an', 'animal, a group, and', 'control or not. Order and', 'well stay with every flag.'], c1 + 16, 548);
  E(c1 + 10, 755, 10, { fill: MUT }); T('mutation', c1 + 22, 745, { color: MUT });
  E(c1 + 110, 755, 10, { fill: ERR }); T('sequencing error', c1 + 122, 745, { color: ERR });

  // ───────── column 2: read preparation ─────────
  const c2 = 290, c2w = 250;
  head(c2, c2 + c2w, 'READ PREPARATION', 'every read is gated');
  const gates = card(c2, 86, c2w, 196);
  cardTitle(c2, 86, 'Gates');
  lines(['Reads longer than the', 'amplicon, or averaging more', 'than 0.01 errors per base,', 'are dropped. Sanger: the', 'longest stretch at Q ≥ 20.'], c2 + 16, 132);
  const align = card(c2, 300, c2w, 250, { stroke: TEAL, strokeWidth: 2 });
  cardTitle(c2, 300, 'Competitive alignment', TEAL);
  lines(['Each read is aligned to the', 'target, its paralogue and any', 'contaminant at once, and', 'goes to its best match. A', 'paralogue can be analysed', 'as a target of its own.'], c2 + 16, 346);
  // three references, one read landing on the first
  [['target', CASE], ['paralogue', MUT], ['other', SOFT]].forEach(([lab, col], k) => {
    R(c2 + 18 + k * 76, 500, 62, 8, { fill: col, radius: 4 });
    TC(lab, c2 + 49 + k * 76, 512, { color: col });
  });
  const acct = card(c2, 568, c2w, 162);
  cardTitle(c2, 568, 'Read accounting');
  lines(['The fate of every delivered', 'read, with its reason, for', 'every library (Figure S1).'], c2 + 16, 614);

  // ───────── column 3: the calibrated analysis ─────────
  const c3 = 570, c3w = 530, half = (c3w - 16) / 2;
  head(c3, c3 + c3w, 'AUTOMATED ANALYSIS', 'calibrated by the AID-null controls');
  const bg = card(c3, 86, half, 214, { fill: '#fdf3f3', stroke: '#e3bcbf' });
  cardTitle(c3, 86, 'Background', CTRL);
  lines(['Measured in the controls,', 'per position, change, strand', 'and base quality, and', 'inflated by one factor', 'fitted leaving each', 'control out in turn.'], c3 + 16, 132);
  const test = card(c3 + half + 16, 86, half, 214, { fill: '#f0f4fa', stroke: '#bfd0e4' });
  cardTitle(c3 + half + 16, 86, 'Candidate test', CASE);
  lines(['A Poisson tail over every', 'depth of evidence, against', 'that background. A change', 'in half the reads or more', 'is germline (Sanger: per', 'sequencing batch).'], c3 + half + 32, 132);
  const thr = card(c3, 318, c3w, 220, { stroke: TEAL, strokeWidth: 2 });
  cardTitle(c3, 318, 'Threshold set by the controls (target–decoy)', TEAL);
  lines(['Each control is called against the others, giving decoy',
    'p-values. The threshold is the largest at which the false',
    'calls the decoys predict are at most 5% of the calls kept.', '',
    'Sanger: every substitution at Q ≥ 20 counts; the AID-null',
    'clones measure the error, and each animal is told how many',
    'of its calls that background explains.'], c3 + 16, 364);
  const qc = card(c3, 556, c3w, 174, { stroke: SANG });
  cardTitle(c3, 556, 'Bookkeeping checks (reported, never applied)', SANG);
  lines(['Copied traces and duplicate clones, as groups; clones',
    'shared between animals beyond chance; mutated control',
    'clones; label mismatches; load outliers. A trace filed under',
    'two animals on purpose is kept under both; an accidental',
    'copy is counted once.'], c3 + 16, 600);

  // ───────── column 4: per animal, then per group ─────────
  const c4 = 1130, c4w = 250;
  head(c4, c4 + c4w, 'RESULTS', 'per animal, then per group');
  const prof = card(c4, 86, c4w, 382);
  cardTitle(c4, 86, 'Mutation profiles');
  lines(['unique-variant density', 'mutation frequency', '12-type spectrum,', 'corrected for base',
    'composition', 'AID (WRC/GYW) and', 'Polη (WA/TW) hotspots', 'position along the', 'amplicon',
    'indels, clones and', 'mutation trees', 'calls the background', 'explains (Sanger)'], c4 + 16, 132, {}, 24);
  const grp = card(c4, 486, c4w, 244, { fill: '#edf5f8', stroke: '#b3d0db' });
  cardTitle(c4, 486, 'Group tests', TEAL);
  lines(['Welch t on log₁₀ density,', 'one value per animal, so', 'clones of one lineage', 'cannot stand in for', 'independent mice.'], c4 + 16, 532);
  // a dot plot: two groups with means
  // a dot plot: one dot per animal, AID-competent against AID-null, means as bars
  R(c4 + 30, 694, c4w - 60, 1.6, { fill: OUTL });
  [[c4 + 85, [656, 662, 668], '#6d9ad0', CASE, 'case'], [c4 + 165, [678, 682, 686], '#d07a7e', CTRL, 'control']].forEach(([cx, ys, f, e, lab]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 13, yy, 11, { fill: f, stroke: e, strokeWidth: 1.6 }));
    const m = ys.reduce((p, q) => p + q) / 3;
    R(cx - 26, m - 1.2, 52, 2.4, { fill: e });
    TC(lab, cx, 698, { color: e });
  });

  // ───────── flow ─────────
  LINK(sheet, 'e', acct, 'w');
  LINK(gates, 's', align, 'n');
  LINK(gates, 'e', bg, 'w');
  LINK(bg, 'e', test, 'w');
  LINK(test, 's', thr, 'n', { style: 'elbow' });
  LINK(thr, 'e', prof, 'w', { style: 'elbow' });
  LINK(prof, 's', grp, 'n');
  LINK(acct, 'e', qc, 'w');

  return { id: 'crm-figS12', name: 'Figure S12', width: W, height: 790, background: '#ffffff', objects: o,
    alt: 'Figure S12. An overview of SHMtrace: libraries and a sample sheet, read preparation, the control-calibrated analysis and bookkeeping checks, per-animal profiles and group tests.' };
}
