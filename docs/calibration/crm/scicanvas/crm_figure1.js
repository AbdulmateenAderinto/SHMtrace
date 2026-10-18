() => {
  // Cell Reports Methods, Figure 1: why SHMtrace calibrates against matched AID-null libraries, and
  // what it does. A, read-level against control-calibrated testing; B, inputs and read preparation;
  // C, calling and checks; D, outputs. Cell Press rules as in figure1_overview.js (174 mm; Arial
  // 6-8 pt; strokes 0.5-1.5 pt; no red with green). Mutations purple, sequencing errors orange-red,
  // AID-null controls red, SHMtrace teal.
  const W = 1400, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8', FAINT = '#cfd8df';
  const TEAL = '#1f6f8b', CTRL = '#9a3439', CASE = '#274b7e', SANG = '#a3621b', ILL = '#3d5a80';
  const MUT = '#6a4c93', ERR = '#d0632a', OKC = '#274b7e', BAD = '#c0392b';
  const F = 'arial', TITLE = 22.5, HEAD = 21, LBL = 19.5, SM = 17.5;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const SEG = (x1, y1, x2, y2, w, fill) => {
    const len = Math.hypot(x2 - x1, y2 - y1), cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    return R(cx - len / 2, cy - w / 2, len, w, { fill, radius: w / 2, rot: (Math.atan2(y2 - y1, x2 - x1) * 180) / Math.PI });
  };
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: INK, width: 2.4, headSize: 0.8, ...e }); o.push(c); return c; };
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const badge = (n, cx, cy, fill) => { const b = Make.badge(n, cx, cy, { fill, labelSize: SM, w: 32, h: 32, x: cx - 16, y: cy - 16 }); o.push(b); return b; };
  const panel = (letter, title, x, y) => { T(letter, x, y, { fontSize: TITLE, bold: true }); T(title, x + 28, y, { fontSize: TITLE, bold: true }); };
  const lines = (arr, x, y, e = {}, dy = 23) => arr.map((s, k) => T(s, x, y + k * dy, { fontSize: SM, color: MID, ...e }));
  const mark = (cx, y, ok) => TC(ok ? '✓' : '✕', cx, y, { fontSize: 20, bold: true, color: ok ? OKC : BAD });

  // a pileup: reads as grey bars; dots at columns cols[c] for reads listed in hits[c]
  const pile = (x, y, w, n, hits, label) => {
    for (let r = 0; r < n; r++) R(x, y + r * 15, w, 7, { fill: FAINT, radius: 3.5 });
    hits.forEach(([fx, reads, colour]) => reads.forEach((r) => E(x + fx * w, y + r * 15 + 3.5, 10, { fill: colour })));
    if (label) TC(label, x + w / 2, y + n * 15 + 4, { fontSize: SM, color: MID });
  };

  // dashed segment, horizontal or vertical
  const DASH = (x1, y1, x2, y2, w, fill, on = 7, off = 5) => {
    const len = Math.hypot(x2 - x1, y2 - y1), ux = (x2 - x1) / len, uy = (y2 - y1) / len;
    for (let d = 0; d < len; d += on + off) { const e = Math.min(d + on, len); SEG(x1 + ux * d, y1 + uy * d, x1 + ux * e, y1 + uy * e, w, fill); }
  };
  const PINK = '#d07a7e', PALE = '#e9eef2', WIN = '#e4f0f4';

  // ───────── A. Where the background comes from ─────────
  // Each half shows the idea as a drawing; the reasoning is in the legend.
  panel('A', 'Where the background comes from', 30, 10);
  const ax = 30, aw = 660;
  // A1: a quality-based test. Mutations and recurrent errors both carry high base quality, so quality cannot
  // tell them apart
  card(ax, 48, aw, 232);
  T('Quality-based test', ax + 16, 58, { fontSize: HEAD, bold: true });
  T('(e.g. LoFreq, Geneious)', ax + 205, 61, { fontSize: SM, color: SOFT });
  T('✓ correct  ✕ wrong', ax + aw - 172, 61, { fontSize: SM, color: SOFT });
  const cols = [[0.3, [0, 2, 4], MUT], [0.72, [1, 3, 5], ERR]];
  pile(ax + 20, 100, 240, 7, cols);
  mark(ax + 20 + 0.3 * 240, 196, true); mark(ax + 20 + 0.72 * 240, 196, false);
  [['mutation', 'called', 0.3, MUT], ['error', 'miscalled', 0.72, ERR]].forEach(([w, v, f, c]) => { TC(w, ax + 20 + f * 240, 228, { fontSize: SM, color: c }); TC(v, ax + 20 + f * 240, 251, { fontSize: SM, color: c }); });
  const qx = ax + 330, qw = 290, qa = 190, qq = (q) => qx + (qw * q) / 40;
  TC('same base quality', qx + qw / 2, 104, { fontSize: SM, bold: true, color: INK });
  R(qq(0), 138, qq(20) - qq(0), 46, { fill: PALE });
  TC('low', qq(10), 150, { fontSize: SM, color: SOFT });
  [[31, MUT], [34, MUT], [37, MUT], [39, MUT]].forEach(([q, c]) => E(qq(q), 148, 12, { fill: c }));
  [[30, ERR], [33, ERR], [36, ERR], [38, ERR]].forEach(([q, c]) => E(qq(q), 172, 12, { fill: c }));
  SEG(qx, qa, qx + qw, qa, 1.6, INK);
  [0, 20, 40].forEach((q) => { SEG(qq(q), qa, qq(q), qa + 7, 1.6, INK); TC(`Q${q}`, qq(q), qa + 9, { fontSize: SM, color: MID }); });
  TC('errors from PCR, homopolymers,', qx + qw / 2, 226, { fontSize: SM, color: MID });
  TC('paralog reads look like mutations', qx + qw / 2, 249, { fontSize: SM, color: MID });
  // A2: the control-calibrated test. The same two columns, now read against matched AID-null controls: the
  // mutation is beyond the controls' background, the error is not
  card(ax, 294, aw, 292, { stroke: TEAL, strokeWidth: 2 });
  T('Control-calibrated test (SHMtrace)', ax + 16, 304, { fontSize: HEAD, bold: true, color: TEAL });
  pile(ax + 20, 346, 200, 6, [[0.3, [0, 2, 4], MUT], [0.75, [1, 3, 5], ERR]], 'sample');
  pile(ax + 240, 346, 120, 6, [[0.75, [0, 2, 3, 5], ERR]], 'controls');
  mark(ax + 20 + 0.3 * 200, 458, true); mark(ax + 20 + 0.75 * 200, 458, true);
  [['mutation', 'called', 0.3, MUT], ['error', 'rejected', 0.75, ERR]].forEach(([w, v, f, c]) => { TC(w, ax + 20 + f * 200, 488, { fontSize: SM, color: c }); TC(v, ax + 20 + f * 200, 511, { fontSize: SM, color: c }); });
  // bars: reads carrying each change, in the sample and in the controls
  const gx0 = ax + 410, base = 500, u = 20;
  T('reads with the change', gx0, 346, { fontSize: SM, color: MID });
  R(gx0, 378, 8, 14, { fill: MUT }); R(gx0 + 8, 378, 8, 14, { fill: ERR });
  T('sample', gx0 + 22, 375, { fontSize: SM, color: MID });
  R(gx0 + 98, 378, 16, 14, { fill: PINK });
  T('controls', gx0 + 120, 375, { fontSize: SM, color: MID });
  [['mutation', MUT, 3, 0, gx0 + 6], ['error', ERR, 3, 4, gx0 + 130]].forEach(([name, col, smp, ctl, gx]) => {
    R(gx, base - smp * u, 36, smp * u, { fill: col });
    R(gx + 42, base - Math.max(ctl * u, 3), 36, Math.max(ctl * u, 3), { fill: PINK });
    TC(name, gx + 39, base + 8, { fontSize: SM, color: col });
  });
  SEG(gx0, base, gx0 + 224, base, 1.6, INK);
  T('AID-null controls from the same run; empirical FDR target 0.05', ax + 20, 552, { fontSize: SM, color: TEAL, bold: true });

  // ───────── B. Inputs and read preparation ─────────
  panel('B', 'Inputs and read preparation', 720, 10);
  const bx = 720, bw = 650;
  const inputs = [
    ['s-nanopore-sequencing-strand-through-pore', 'Nanopore', 'full-length amplicons', INK, 74],
    ['s-sanger-chromatogram-trace-only', 'Sanger', '.ab1 traces + sample sheet', SANG, 74],
    ['s-flow-cell-clusters-illumina', 'Illumina', 'paired reads, merged', ILL, 104],
    ['s-umi-tagged-pcr-products', 'UMI-tagged', 'one consensus per molecule', TEAL, 96],
  ];
  const iw = (bw - 3 * 12) / 4, inCards = [];
  inputs.forEach(([icon, title, sub, c, size], k) => {
    const x = bx + k * (iw + 12);
    inCards.push(card(x, 48, iw, 182));
    I(icon, x + iw / 2, 100, size);
    TC(title, x + iw / 2, 150, { bold: true, color: c });
    const words = sub.split(' '), half = Math.ceil(words.length / 2);
    TC(words.slice(0, half).join(' '), x + iw / 2, 176, { fontSize: SM, color: MID });
    TC(words.slice(half).join(' '), x + iw / 2, 198, { fontSize: SM, color: MID });
  });
  const prep = card(bx, 294, bw, 292, { stroke: TEAL });
  T('Read preparation', bx + 16, 304, { fontSize: HEAD, bold: true, color: TEAL });
  // four steps as four small drawings, each with a two-word name and a two-line caption
  const pw = (bw - 32) / 4, px = (k) => bx + 16 + k * pw;
  const prepSteps = [
    [['Length and', 'error gates'], ['over-length or', '> 1% errors/base']],
    [['Target or', 'paralog?'], ['aligned to both,', 'best match kept']],
    [['Paralog', 'alleles'], ['≥ 10 paralog', 'bases: removed']],
    [['Matched', 'depth'], ['random draw,', 'same per library']],
  ];
  prepSteps.forEach(([title, cap], k) => {
    badge(k + 1, px(k) + 18, 362, TEAL);
    T(title[0], px(k) + 40, 342, { bold: true, fontSize: SM });
    T(title[1], px(k) + 40, 364, { bold: true, fontSize: SM });
    TC(cap[0], px(k) + pw / 2, 492, { fontSize: SM, color: MID });
    TC(cap[1], px(k) + pw / 2, 514, { fontSize: SM, color: MID });
  });
  const X = (cx, y) => TC('✕', cx, y, { fontSize: 18, bold: true, color: BAD });
  // 1: one read too long, one full of errors; both out
  let x0 = px(0) + 14;
  R(x0, 404, 96, 7, { fill: FAINT, radius: 3.5 });
  R(x0, 422, 124, 7, { fill: PALE, radius: 3.5 }); X(x0 + 116, 412);
  R(x0, 440, 96, 7, { fill: FAINT, radius: 3.5 });
  R(x0, 458, 96, 7, { fill: PALE, radius: 3.5 });
  [0.1, 0.3, 0.45, 0.62, 0.8].forEach((f) => E(x0 + f * 96, 461.5, 8, { fill: '#eab497' })); X(x0 + 112, 448);
  // 2: a read placed against target and paralog references; it matches the target
  x0 = px(1) + 12;
  R(x0, 410, 84, 7, { fill: FAINT, radius: 3.5 });
  [0.15, 0.45, 0.8].forEach((f) => SEG(x0 + f * 84, 418, x0 + f * 84, 432, 1.6, TEAL));
  R(x0, 434, 84, 7, { fill: TEAL, radius: 3.5 }); T('target', x0 + 90, 427, { fontSize: SM, color: TEAL });
  R(x0, 458, 84, 7, { fill: OUTL, radius: 3.5 }); T('paralog', x0 + 90, 451, { fontSize: SM, color: SOFT });
  // 3: a read carrying many paralog-specific bases is removed; one with a single mutation is kept
  x0 = px(2) + 16;
  R(x0, 416, 104, 7, { fill: FAINT, radius: 3.5 });
  for (let j = 0; j < 10; j++) R(x0 + 6 + j * 9.6, 412, 3, 15, { fill: SOFT });
  X(x0 + 120, 406);
  R(x0, 452, 104, 7, { fill: FAINT, radius: 3.5 }); R(x0 + 60, 448, 3, 15, { fill: MUT });
  TC('✓', x0 + 120, 442, { fontSize: 18, bold: true, color: OKC });
  // 4: three libraries of unequal depth cut to the same number of reads
  x0 = px(3) + 16;
  [7, 4, 5].forEach((n, j) => {
    for (let r = 0; r < n; r++) R(x0 + j * 40, 470 - r * 10, 32, 6, { fill: r < 4 ? '#8d99a6' : PALE, radius: 3 });
  });
  DASH(x0 - 4, 438, x0 + 116, 438, 1.6, TEAL);
  // Sanger traces share steps 2 and 3; trimming replaces step 1, and step 4 cannot apply because a clone is one
  // molecule and every clone is kept
  T('Sanger: trimming at Q ≥ 20 replaces 1; every clone is kept, so no 4', bx + 16, 548, { fontSize: SM, color: SANG });
  inCards.forEach((k) => LINK(k, 's', prep, 'n', { color: MID, style: 'elbow', headSize: 0.6 }));

  // ───────── C. Calling and checks ─────────
  // Two implemented routes: read-library candidates are tested against a control error model;
  // repeated observations of the same original molecule are not required. The Sanger route
  // counts quality-filtered substitutions without a per-call test, and
  // the AID-null clones give the background those counts carry, reported beside them rather than subtracted.
  // The two routes rejoin at the germline filter, the window and the automatic checks.
  const cy = 620;
  panel('C', 'Calling: a test for reads, a count for Sanger clones', 30, cy);
  const cw = (W - 60 - 3 * 18) / 4, colX = (k) => 30 + k * (cw + 18);
  const R1 = cy + 58, R2 = cy + 228, CH = 132, SHH = R2 + CH - R1;
  T('Reads (nanopore, Illumina, UMI): tested against controls', 30, cy + 32, { fontSize: SM, bold: true, color: TEAL });
  T('Both routes', colX(2), cy + 32, { fontSize: SM, bold: true, color: MID });
  T('Sanger clones: counted', 30, cy + 202, { fontSize: SM, bold: true, color: SANG });
  const step = (k, y, h, title, colour, n) => {
    const x = colX(k);
    const c = card(x, y, cw, h, { stroke: colour });
    badge(n, x + 24, y + 26, colour);
    T(title, x + 46, y + 14, { bold: true });
    return [c, x];
  };
  const side = (arr, x, y) => lines(arr, x, y, { color: MID }, 22);
  // 1 (reads): an error rate for every substitution type, strand and quality bin, plus a margin
  const [c1, x1] = step(0, R1, CH, 'Control error model', TEAL, 1);
  const tint = ['#d6e8ee', '#a9cfdb', '#e8f2f5', '#7fb3c4', '#c3dde6', '#eef6f8', '#4f92a8', '#b6d6e0'];
  for (let r = 0; r < 3; r++) for (let q = 0; q < 8; q++) R(x1 + 16 + q * 16, R1 + 56 + r * 16, 14, 14, { fill: tint[(q * 3 + r * 5) % 8] });
  side(['rate per type,', 'strand, quality', '+ held-out margin'], x1 + 160, R1 + 50);
  // 2 (reads): candidates against decoys from held-out controls; the threshold keeps estimated FDR <= 0.05
  const [c2, x2] = step(1, R1, CH, 'Test and FDR', TEAL, 2);
  [[0.0], [0.11], [0.22], [0.55], [0.8]].forEach(([f]) => E(x2 + 20 + f * 120, R1 + 66, 11, { fill: MUT }));
  [[0.36], [0.5], [0.62], [0.74], [0.88], [0.97]].forEach(([f]) => E(x2 + 20 + f * 120, R1 + 90, 11, { fill: '#ffffff', stroke: CTRL, strokeWidth: 1.5 }));
  SEG(x2 + 16, R1 + 106, x2 + 144, R1 + 106, 1.6, INK);
  DASH(x2 + 20 + 0.31 * 120, R1 + 52, x2 + 20 + 0.31 * 120, R1 + 104, 1.6, TEAL, 5, 4);
  side(['decoys (open)', 'set the threshold:', 'FDR target 0.05'], x2 + 160, R1 + 50);
  // 1 (Sanger): every substitution at Q >= 20 in a trace is counted
  const [s1, xs1] = step(0, R2, CH, 'Count the mismatches', SANG, 1);
  I('s-sanger-chromatogram-trace-only', xs1 + 80, R2 + 86, 96);
  E(xs1 + 84, R2 + 60, 11, { fill: MUT });
  side(['every substitution', 'at quality ≥ 20', 'is counted'], xs1 + 160, R2 + 50);
  // 2 (Sanger): the AID-null clones' count is shown beside the calls, not subtracted from them
  const [s2, xs2] = step(1, R2, CH, 'Background reported', SANG, 2);
  R(xs2 + 28, R2 + 54, 40, 48, { fill: MUT }); R(xs2 + 84, R2 + 90, 40, 12, { fill: PINK });
  SEG(xs2 + 16, R2 + 102, xs2 + 138, R2 + 102, 1.6, INK);
  TC('Wt', xs2 + 48, R2 + 104, { fontSize: SM, color: MID }); TC('AID-null', xs2 + 104, R2 + 104, { fontSize: SM, color: MID });
  side(['shown beside', 'the calls, never', 'subtracted'], xs2 + 160, R2 + 50);
  // 3 (both): germline alleles out; calls only inside the covered window
  const [sh1, xg] = step(2, R1, SHH, 'Germline and window', MID, 3);
  const gx = xg + 20, gw = cw - 40, wl = gx + 0.18 * gw, wr = gx + 0.86 * gw;
  R(wl, R1 + 58, wr - wl, 110, { fill: WIN });
  TC('covered window', (wl + wr) / 2, R1 + 62, { fontSize: SM, color: TEAL });
  const gcol = 0.6, hits3 = [[0.3, 1], [0.42, 3], [0.74, 0], [0.5, 4], [0.07, 2], [0.95, 5]];
  for (let r = 0; r < 6; r++) {
    R(gx, R1 + 88 + r * 14, gw, 7, { fill: FAINT, radius: 3.5 });
    E(gx + gcol * gw, R1 + 91.5 + r * 14, 10, { fill: INK });
  }
  hits3.forEach(([f, r]) => {
    const inWin = gx + f * gw > wl && gx + f * gw < wr;
    E(gx + f * gw, R1 + 91.5 + r * 14, 10, inWin ? { fill: MUT } : { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.5 });
  });
  X(gx + gcol * gw, R1 + 172);
  lines(['window first, then alleles in', '≥ 50% of reads set aside and', 'reported, not dropped, before', 'any test; indels optional', '(opt-in: known ancestor)'], xg + 16, R1 + 192, { color: MID }, 20);
  // 4 (both): the automatic checks, as four glyphs
  const [sh2, xc] = step(3, R1, SHH, 'Automatic checks', MID, 4);
  const gwid = (cw - 36) / 2, gX = (j) => xc + 12 + j * (gwid + 12), gY = (i) => R1 + 62 + i * 116;
  // spectrum: transitions share of G:C calls
  [[46, '#2c5d8a'], [22, '#7fa7cf'], [30, '#d9952b'], [18, '#f0c47e']].forEach(([h, c], j) => R(gX(0) + 18 + j * 26, gY(0) + 60 - h, 20, h, { fill: c }));
  SEG(gX(0) + 12, gY(0) + 60, gX(0) + 124, gY(0) + 60, 1.6, INK);
  TC('spectrum', gX(0) + gwid / 2, gY(0) + 66, { fontSize: SM, color: MID });
  // positional clustering: calls bunch along the amplicon
  SEG(gX(1) + 8, gY(0) + 44, gX(1) + gwid - 8, gY(0) + 44, 1.6, OUTL);
  [0.12, 0.16, 0.2, 0.23, 0.27, 0.6, 0.66, 0.7, 0.9].forEach((f, j) => E(gX(1) + 8 + f * (gwid - 16), gY(0) + 32 - (j % 2) * 10, 9, { fill: MUT }));
  TC('clustering', gX(1) + gwid / 2, gY(0) + 66, { fontSize: SM, color: MID });
  // 5-mer targeting: calls fall at mutable 5-mers
  const heat = ['#fbe7cf', '#f2b880', '#e07b39', '#fbe7cf', '#f6d2a8', '#c75a1e', '#fbe7cf', '#f2b880'];
  heat.forEach((c, j) => R(gX(0) + 12 + j * 14, gY(1) + 30, 13, 18, { fill: c }));
  [2, 5].forEach((j) => SEG(gX(0) + 18.5 + j * 14, gY(1) + 8, gX(0) + 18.5 + j * 14, gY(1) + 26, 2.4, MUT));
  TC('5-mer targeting', gX(0) + gwid / 2, gY(1) + 66, { fontSize: SM, color: MID });
  // Sanger bookkeeping: a clone filed under the wrong animal is flagged
  const fx = gX(1) + 30;
  SEG(fx, gY(1) + 4, fx, gY(1) + 58, 2.4, INK);
  R(fx + 1, gY(1) + 4, 40, 24, { fill: SANG });
  R(fx + 52, gY(1) + 24, 44, 8, { fill: FAINT, radius: 4 }); R(fx + 52, gY(1) + 40, 44, 8, { fill: SANG, radius: 4 });
  TC('Sanger filing', gX(1) + gwid / 2, gY(1) + 66, { fontSize: SM, color: MID });
  LINK(c1, 'e', c2, 'w', { headSize: 0.7 });
  LINK(s1, 'e', s2, 'w', { headSize: 0.7 });
  LINK(c2, 'e', sh1, 'w', { headSize: 0.7, style: 'elbow' });
  LINK(s2, 'e', sh1, 'w', { headSize: 0.7, style: 'elbow' });
  LINK(sh1, 'e', sh2, 'w', { headSize: 0.7 });

  // ───────── D. Outputs ─────────
  const dy = cy + 396;
  panel('D', 'Outputs', 30, dy);
  const tw = (W - 60 - 5 * 14) / 6, tiles = [];
  const tile = (k, title, sub) => {
    const x = 30 + k * (tw + 14);
    tiles.push(card(x, dy + 38, tw, 228));
    T(title[0], x + 12, dy + 48, { bold: true, fontSize: SM });
    if (title[1]) T(title[1], x + 12, dy + 70, { bold: true, fontSize: SM });
    lines(sub, x + 12, dy + 196, { color: SOFT }, 21);
    return x;
  };
  // animal-level tests: two groups of dots with means
  let x = tile(0, ['Animal-level tests', 'and power'], ['NB, Welch,', 'permutation']);
  [[x + 44, [dy + 108, dy + 118, dy + 130], '#6d9ad0', CASE], [x + 128, [dy + 160, dy + 166, dy + 170], '#d07a7e', CTRL]].forEach(([cx, ys, f, e]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 13, yy, 12, { fill: f, stroke: e, strokeWidth: 1.5 }));
    const m = ys.reduce((a, b) => a + b) / 3; SEG(cx - 26, m, cx + 26, m, 2.4, e);
  });
  SEG(x + 14, dy + 184, x + tw - 14, dy + 184, 1.6, INK);
  // single-read estimate: observed vs background bars
  x = tile(1, ['Single-read', 'estimate'], ['observed minus', 'control background']);
  R(x + 30, dy + 112, 40, 72, { fill: MUT }); R(x + 100, dy + 160, 40, 24, { fill: '#d07a7e' });
  SEG(x + 14, dy + 184, x + tw - 14, dy + 184, 1.6, INK);
  // clones: reads grouped by the V(D)J junction at their 5' end (coloured block); a shared junction marks one
  // clone (three reads of clone 1, two of clone 2), and the clone-weighted rate counts each clone once
  x = tile(2, ['Clones from the', 'V(D)J junction'], ['reads grouped by', 'junction; each clone', 'counted once']);
  const JN = [CASE, CASE, CASE, '#b8742a', '#b8742a'], TICKS = [[0.35, 0.7], [0.35, 0.7], [0.35], [0.55], [0.55, 0.85]];
  JN.forEach((jc, k) => {
    const yy = dy + 104 + k * 16 + (k >= 3 ? 4 : 0);
    R(x + 14, yy, 30, 10, { fill: jc, radius: 3 });
    R(x + 46, yy + 2, tw - 64, 6, { fill: FAINT, radius: 3 });
    TICKS[k].forEach((f) => R(x + 46 + f * (tw - 64) - 1.5, yy - 1, 3, 12, { fill: MUT }));
  });
  // mutation tree
  x = tile(3, ['Mutation', 'trees'], ['full-length reads', 'of one clone']);
  const rx = x + 24, ry = dy + 146;
  R(rx - 6, ry - 6, 12, 12, { fill: INK }); SEG(rx, ry, rx + 50, ry, 2, INK);
  [[ry - 34, 96], [ry - 12, 132], [ry + 12, 74], [ry + 34, 112]].forEach(([yy, len]) => {
    SEG(rx + 50, ry, rx + 50, yy, 1.6, SOFT); SEG(rx + 50, yy, rx + 50 + Math.min(len, tw - 90), yy, 2, INK);
  });
  // descriptive spectrum and sequence-context summaries
  x = tile(4, ['Spectrum and', 'motif shares'], ['G:C / A:T, Ts / Tv,', 'WRC/GYW, WA/TW']);
  [[0.35, '#2c5d8a'], [0.25, '#7fa7cf'], [0.4, '#d9952b']].reduce((acc, [f, c]) => { R(x + 40, dy + 104 + acc * 80, 50, f * 80, { fill: c }); return acc + f; }, 0);
  [[0.66, '#2c5d8a'], [0.29, '#7fa7cf'], [0.05, '#d9952b']].reduce((acc, [f, c]) => { R(x + 110, dy + 104 + acc * 80, 50, f * 80, { fill: c }); return acc + f; }, 0);
  // reports and interfaces
  x = tile(5, ['Reports, app', 'and command line'], ['HTML report, tables,', 'project files']);
  o.push(Make.image(LOGO, x + tw / 2 - 34, dy + 102, 68, 68));
  return { id: 'crm-fig1', name: 'Figure 1', width: W, height: dy + 280, background: '#ffffff', objects: o,
    alt: 'Figure 1. SHMtrace: control-calibrated somatic hypermutation calling, from inputs and read preparation through calling, quality checks and outputs.' };
}
