() => {
  // Cell Reports Methods, Figure S22: an overview of SHMtrace, from libraries to group tests. Built with
  // manuscript_corrected/figure1_scicanvas/build.js into crm/figures/FigureS22_overview; Figure S23 gives the
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
  const SEG = (x1, y1, x2, y2, w, fill) => {
    const len = Math.hypot(x2 - x1, y2 - y1), cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    return R(cx - len / 2, cy - w / 2, len, w, { fill, radius: w / 2, rot: (Math.atan2(y2 - y1, x2 - x1) * 180) / Math.PI });
  };
  const DASH = (x1, y1, x2, y2, w, fill, on = 5, off = 4) => {
    const len = Math.hypot(x2 - x1, y2 - y1), ux = (x2 - x1) / len, uy = (y2 - y1) / len;
    for (let d = 0; d < len; d += on + off) { const e = Math.min(d + on, len); SEG(x1 + ux * d, y1 + uy * d, x1 + ux * e, y1 + uy * e, w, fill); }
  };
  const X = (cx, y) => TC('✕', cx, y, { fontSize: 18, bold: true, color: '#c0392b' });
  const OK = (cx, y) => TC('✓', cx, y, { fontSize: 18, bold: true, color: CASE });
  const PINK = '#d07a7e', PALE = '#e9eef2';
  const cardTitle = (x, y, text, colour = INK) => T(text, x + 16, y + 12, { fontSize: HEAD, bold: true, color: colour });

  // ───────── column 1: libraries ─────────
  const c1 = 20, c1w = 240;
  head(c1, c1 + c1w, 'LIBRARIES', 'reads or Sanger traces');
  TC('AID-competent animals', c1 + c1w / 2, 86, { bold: true, color: CASE });
  pile(c1 + 20, 118, c1w - 40, 8, [[0.2, 0, MUT], [0.55, 0, MUT], [0.35, 1, MUT], [0.8, 2, MUT], [0.2, 3, MUT],
    [0.62, 4, MUT], [0.45, 5, MUT], [0.7, 5, ERR], [0.3, 6, MUT], [0.85, 7, MUT]]);
  TC('AID-null controls (two or more)', c1 + c1w / 2, 254, { bold: true, color: CTRL });
  pile(c1 + 20, 286, c1w - 40, 6, [[0.7, 1, ERR], [0.7, 4, ERR], [0.35, 3, ERR]]);
  // four platforms in a 2 x 2 grid, each label beside its icon
  const plat = [['s-nanopore-sequencing-strand-through-pore', 'nanopore', 38, INK],
    ['s-flow-cell-clusters-illumina', 'Illumina', 42, ILL],
    ['s-long-read-sequencer-pacbio-style', 'PacBio', 38, CTRL],
    ['s-sanger-chromatogram-trace-only', 'Sanger', 40, SANG]];
  plat.forEach(([id, lab, sz, colour], k) => {
    const x = c1 + (k % 2) * 142, cy = 404 + Math.floor(k / 2) * 50;
    I(id, x + 24, cy, sz);
    T(lab, x + 54, cy - 11, { color: colour });
  });
  const sheet = card(c1, 500, c1w, 230);
  cardTitle(c1, 500, 'Sample sheet');
  // a small sample sheet: file, animal, control
  [['file', 'animal', 'control'], ['a.fastq', 'Wt 2', 'no'], ['b.fastq', 'KO 1', 'yes'], ['c12.ab1', 'Wt 3', 'no']].forEach((row, r) => {
    const y = 552 + r * 30;
    if (r > 0) R(c1 + 12, y - 6, c1w - 24, 1.2, { fill: FAINT });
    row.forEach((cell, j) => T(cell, c1 + 16 + [0, 86, 156][j], y, { bold: r === 0, color: r === 0 ? INK : (j === 2 && cell === 'yes' ? CTRL : MID) }));
  });
  lines(['labels stay with every flag'], c1 + 16, 680, { color: SOFT });
  E(c1 + 10, 755, 10, { fill: MUT }); T('mutation', c1 + 22, 745, { color: MUT });
  E(c1 + 110, 755, 10, { fill: ERR }); T('sequencing error', c1 + 122, 745, { color: ERR });

  // ───────── column 2: read preparation ─────────
  const c2 = 290, c2w = 250;
  head(c2, c2 + c2w, 'READ PREPARATION', 'every read is gated');
  const gates = card(c2, 86, c2w, 196);
  cardTitle(c2, 86, 'Gates');
  // one read too long, one full of errors: both out
  const gx = c2 + 22;
  R(gx, 140, 150, 7, { fill: FAINT, radius: 3.5 });
  R(gx, 158, 196, 7, { fill: PALE, radius: 3.5 }); X(gx + 210, 148);
  R(gx, 176, 150, 7, { fill: FAINT, radius: 3.5 });
  R(gx, 194, 150, 7, { fill: PALE, radius: 3.5 });
  [0.1, 0.3, 0.45, 0.62, 0.8].forEach((f) => E(gx + f * 150, 197.5, 8, { fill: '#eab497' })); X(gx + 168, 184);
  lines(['too long, or > 1% errors/base', 'Sanger: Q ≥ 20 core kept'], c2 + 16, 220);
  const align = card(c2, 300, c2w, 250, { stroke: TEAL, strokeWidth: 2 });
  cardTitle(c2, 300, 'Competitive alignment', TEAL);
  // one read, three references; the read lands on its best match
  R(c2 + 24, 362, 62, 7, { fill: FAINT, radius: 3.5 });
  [0.15, 0.5, 0.85].forEach((f) => SEG(c2 + 24 + f * 62, 372, c2 + 24 + f * 62, 404, 1.4, CASE));
  [['target', CASE], ['paralogue', MUT], ['other', SOFT]].forEach(([lab, col], k) => {
    R(c2 + 18 + k * 76, 408, 62, 8, { fill: col, radius: 4 });
    TC(lab, c2 + 49 + k * 76, 420, { color: col });
  });
  OK(c2 + 49, 444);
  lines(['each read goes to its best', 'match; a paralogue can be', 'a target of its own'], c2 + 16, 476);
  const acct = card(c2, 568, c2w, 162);
  cardTitle(c2, 568, 'Read accounting');
  // the fate of every read as one stacked bar
  [[0.64, TEAL], [0.14, OUTL], [0.12, MUT], [0.06, ERR], [0.04, SOFT]].reduce((acc, [f, c]) => {
    R(c2 + 18 + acc * (c2w - 36), 616, f * (c2w - 36), 18, { fill: c }); return acc + f; }, 0);
  T('used', c2 + 18, 640, { color: TEAL });
  T('dropped', c2 + c2w - 84, 640, { color: MID });
  lines(['each drop with its reason,', 'every library (Figure S2)'], c2 + 16, 670, { color: SOFT });

  // ───────── column 3: the calibrated analysis ─────────
  const c3 = 570, c3w = 530, half = (c3w - 16) / 2;
  head(c3, c3 + c3w, 'AUTOMATED ANALYSIS', 'calibrated by the AID-null controls');
  const bg = card(c3, 86, half, 214, { fill: '#fdf3f3', stroke: '#e3bcbf' });
  cardTitle(c3, 86, 'Background', CTRL);
  // an error-rate grid: change x quality
  const tint = ['#f6dcdc', '#e9b4b6', '#fbeeee', '#d98c90', '#efc8ca', '#fdf6f6', '#c7686d', '#e2a3a6'];
  for (let r = 0; r < 3; r++) for (let q = 0; q < 9; q++) R(c3 + 22 + q * 22, 138 + r * 20, 19, 17, { fill: tint[(q * 3 + r * 5) % 8] });
  lines(['per position, change,', 'strand and quality,', '+ a leave-one-out margin'], c3 + 16, 212, {}, 21);
  const test = card(c3 + half + 16, 86, half, 214, { fill: '#f0f4fa', stroke: '#bfd0e4' });
  cardTitle(c3 + half + 16, 86, 'Candidate test', CASE);
  // expected error counts (bars) and the observed count far in the tail
  const tx = c3 + half + 34;
  [44, 66, 50, 28, 14, 6, 3, 1].forEach((h, j) => R(tx + j * 18, 200 - h, 14, h, { fill: '#bfd0e4' }));
  SEG(tx - 4, 200, tx + 214, 200, 1.4, INK);
  SEG(tx + 196, 138, tx + 196, 200, 2.6, MUT); E(tx + 196, 136, 11, { fill: MUT });
  TC('observed', tx + 196, 206, { color: MUT });
  lines(['Poisson tail against the', 'background; in ≥ 50% of', 'reads: germline, set aside'], c3 + half + 32, 230, {}, 21);
  const thr = card(c3, 318, c3w, 220, { stroke: TEAL, strokeWidth: 2 });
  cardTitle(c3, 318, 'Threshold set by the controls (target–decoy)', TEAL);
  // candidates (filled) against decoys from controls called against each other (open); the threshold
  const dx = c3 + 24, dw = 250;
  [0.0, 0.05, 0.1, 0.15, 0.2, 0.26, 0.6, 0.85].forEach((f) => E(dx + f * dw, 384, 12, { fill: MUT }));
  [0.4, 0.5, 0.58, 0.68, 0.76, 0.84, 0.92, 0.99].forEach((f) => E(dx + f * dw, 412, 12, { fill: '#ffffff', stroke: CTRL, strokeWidth: 1.5 }));
  SEG(dx - 6, 430, dx + dw + 6, 430, 1.4, INK);
  DASH(dx + 0.31 * dw, 368, dx + 0.31 * dw, 428, 1.6, TEAL);
  T('p-value →', dx + dw - 70, 436, { color: SOFT });
  TC('kept', dx + 0.14 * dw, 436, { color: TEAL });
  lines(['calls (filled) against', 'decoys: each control', 'called against the', 'others (open)'], c3 + 300, 366, {}, 21);
  lines(['threshold: predicted false calls ≤ 5% of calls kept'], c3 + 16, 466, { color: TEAL, bold: true });
  lines(['Sanger: Q ≥ 20 changes counted; background shown beside'], c3 + 16, 498, { color: SANG });
  const qc = card(c3, 556, c3w, 174, { stroke: SANG });
  cardTitle(c3, 556, 'Bookkeeping checks (reported, never applied)', SANG);
  // five checks as five small drawings
  const bw5 = (c3w - 32) / 5, bX = (j) => c3 + 16 + j * bw5, by = 604;
  I('s-sanger-chromatogram-trace-only', bX(0) + bw5 / 2 - 8, by + 22, 50); I('s-sanger-chromatogram-trace-only', bX(0) + bw5 / 2 + 8, by + 30, 50);
  R(bX(1) + 18, by + 22, 64, 7, { fill: SANG, radius: 3.5 });
  E(bX(1) + 22, by + 4, 14, { fill: '#6d9ad0' }); E(bX(1) + 78, by + 4, 14, { fill: '#6d9ad0' });
  SEG(bX(1) + 24, by + 10, bX(1) + 40, by + 21, 1.4, SOFT); SEG(bX(1) + 76, by + 10, bX(1) + 60, by + 21, 1.4, SOFT);
  R(bX(2) + 14, by + 20, 72, 10, { fill: '#ffffff', stroke: CTRL, strokeWidth: 1.5, radius: 5 });
  [0.2, 0.5, 0.8].forEach((f) => E(bX(2) + 14 + f * 72, by + 25, 8, { fill: MUT }));
  T('Wt 2', bX(3) + 12, by + 8, { color: MID }); TC('≠', bX(3) + 62, by + 8, { bold: true, color: SANG }); T('Wt 4', bX(3) + 70, by + 8, { color: MID });
  [0.1, 0.2, 0.25, 0.35, 0.4, 0.95].forEach((f, j) => E(bX(4) + 10 + f * 80, by + 22 - (j === 5 ? 0 : (j % 2) * 8), 9, { fill: j === 5 ? SANG : OUTL }));
  [['copied', 'trace'], ['clone in', 'two animals'], ['mutated', 'control'], ['label', 'mismatch'], ['load', 'outlier']].forEach(([a, b], j) => {
    TC(a, bX(j) + bw5 / 2, by + 50, { color: MID }); TC(b, bX(j) + bw5 / 2, by + 72, { color: MID });
  });

  // ───────── column 4: per animal, then per group ─────────
  const c4 = 1130, c4w = 250;
  head(c4, c4 + c4w, 'RESULTS', 'per animal, then per group');
  const prof = card(c4, 86, c4w, 382);
  cardTitle(c4, 86, 'Mutation profiles');
  // six profiles as six glyphs, two per row
  const pw2 = (c4w - 36) / 2, pX = (j) => c4 + 12 + j * (pw2 + 12), pY = (i) => 140 + i * 108;
  const plab = (j, i, text) => TC(text, pX(j) + pw2 / 2, pY(i) + 70, { color: MID });
  // density and frequency: wild type against control
  R(pX(0) + 26, pY(0) + 8, 26, 56, { fill: MUT }); R(pX(0) + 60, pY(0) + 50, 26, 14, { fill: PINK });
  SEG(pX(0) + 14, pY(0) + 64, pX(0) + pw2 - 6, pY(0) + 64, 1.4, INK); plab(0, 0, 'density');
  // spectrum
  [[0.4, '#2c5d8a'], [0.2, '#7fa7cf'], [0.25, '#d9952b'], [0.15, '#f0c47e']].reduce((acc, [f, c]) => {
    R(pX(1) + 30, pY(0) + 6 + acc * 58, 46, f * 58, { fill: c }); return acc + f; }, 0);
  plab(1, 0, 'spectrum');
  // hotspots: a motif within the sequence
  ['#e9eef2', '#f2b880', '#e07b39', '#f2b880', '#e9eef2', '#e9eef2'].forEach((c, j) => R(pX(0) + 10 + j * 15, pY(1) + 30, 13, 20, { fill: c }));
  SEG(pX(0) + 46.5, pY(1) + 6, pX(0) + 46.5, pY(1) + 26, 2.4, MUT); plab(0, 1, 'hotspots');
  // position along the amplicon
  [[0, 40], [14, 18], [28, 30], [42, 12], [56, 26], [70, 42], [84, 48]].reduce((p, q) => {
    if (p) SEG(pX(1) + 10 + p[0], pY(1) + 10 + p[1], pX(1) + 10 + q[0], pY(1) + 10 + q[1], 2, CASE); return q; }, null);
  SEG(pX(1) + 6, pY(1) + 62, pX(1) + pw2 - 6, pY(1) + 62, 1.4, INK); plab(1, 1, 'position');
  // clones and trees
  const tx0 = pX(0) + 18, ty0 = pY(2) + 32;
  R(tx0 - 5, ty0 - 5, 10, 10, { fill: INK }); SEG(tx0, ty0, tx0 + 26, ty0, 2, INK);
  [[-22, 50], [-7, 64], [8, 42], [22, 56]].forEach(([dy, len]) => { SEG(tx0 + 26, ty0, tx0 + 26, ty0 + dy, 1.4, SOFT); SEG(tx0 + 26, ty0 + dy, tx0 + 26 + len, ty0 + dy, 2, INK); });
  plab(0, 2, 'clones, trees');
  // indels: a deletion gap and an insertion
  R(pX(1) + 8, pY(2) + 30, 36, 8, { fill: FAINT, radius: 4 }); R(pX(1) + 62, pY(2) + 30, 34, 8, { fill: FAINT, radius: 4 });
  DASH(pX(1) + 46, pY(2) + 34, pX(1) + 60, pY(2) + 34, 1.4, ERR, 3, 3);
  R(pX(1) + 78, pY(2) + 14, 4, 14, { fill: MUT }); plab(1, 2, 'indels');
  const grp = card(c4, 486, c4w, 244, { fill: '#edf5f8', stroke: '#b3d0db' });
  cardTitle(c4, 486, 'Group tests', TEAL);
  lines(['Welch t on log₁₀ density,', 'one value per animal'], c4 + 16, 532);
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

  return { id: 'crm-figS13', name: 'Figure S23', width: W, height: 790, background: '#ffffff', objects: o,
    alt: 'Figure S23. An overview of SHMtrace: libraries and a sample sheet, read preparation, the control-calibrated analysis and bookkeeping checks, per-animal profiles and group tests.' };
}
