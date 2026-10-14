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

  // ───────── A. Where the background comes from ─────────
  panel('A', 'Where the background comes from', 30, 10);
  const ax = 30, aw = 660;
  // A1: read-level test
  const a1 = card(ax, 48, aw, 232);
  T('Read-level test', ax + 16, 58, { fontSize: HEAD, bold: true });
  T('(e.g. LoFreq, Geneious)', ax + 184, 61, { fontSize: SM, color: SOFT });
  const cols = [[0.3, [0, 2, 4], MUT], [0.72, [1, 3, 5], ERR]];
  pile(ax + 20, 100, 240, 7, cols);
  mark(ax + 20 + 0.3 * 240, 196, true); mark(ax + 20 + 0.72 * 240, 196, false);
  [['mutation', 0.3, MUT], ['error', 0.72, ERR]].forEach(([w, f, c]) => { TC(w, ax + 20 + f * 240, 228, { fontSize: SM, color: c }); TC('called', ax + 20 + f * 240, 251, { fontSize: SM, color: c }); });
  TC('✓ correct decision   ✕ wrong decision', ax + 20 + 0.5 * 240 + 330, 246, { fontSize: SM, color: SOFT });
  lines(['Each base is judged only by its own quality', 'score, so errors that recur at one position',
    '(PCR, homopolymers, paralog reads) are', 'called as mutations, and reads from one', 'animal are treated as independent',
    'replicates.'], ax + 290, 104);
  // A2: control-calibrated test
  const a2 = card(ax, 294, aw, 250, { stroke: TEAL, strokeWidth: 2 });
  T('Control-calibrated test (SHMtrace)', ax + 16, 304, { fontSize: HEAD, bold: true, color: TEAL });
  pile(ax + 20, 346, 200, 6, [[0.3, [0, 2, 4], MUT], [0.75, [1, 3, 5], ERR]], 'sample');
  pile(ax + 240, 346, 120, 6, [[0.75, [0, 2, 3, 5], ERR]], 'controls');
  mark(ax + 20 + 0.3 * 200, 458, true); mark(ax + 20 + 0.75 * 200, 458, true);
  [['mutation', 'called', 0.3, MUT], ['error', 'rejected', 0.75, ERR]].forEach(([w, v, f, c]) => { TC(w, ax + 20 + f * 200, 488, { fontSize: SM, color: c }); TC(v, ax + 20 + f * 200, 511, { fontSize: SM, color: c }); });
  lines(['Error rates per substitution,', 'strand and base quality are', 'measured in matched AID-null',
    'libraries from the same run.', 'A call needs support beyond', 'this background (FDR ≤ 0.05),', 'and groups are compared with',
    'one value per animal.'], ax + 378, 346);

  // ───────── B. Inputs and read preparation ─────────
  panel('B', 'Inputs and read preparation', 720, 10);
  const bx = 720, bw = 650;
  const inputs = [
    ['s-nanopore-sequencing-strand-through-pore', 'Nanopore', 'full-length amplicons', INK, 74],
    ['s-sanger-chromatogram-trace-only', 'Sanger', '.ab1 traces + sample sheet', SANG, 74],
    ['s-production-sequencer-novaseq-style', 'Illumina', 'paired reads, merged', ILL, 100],
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
  const prep = card(bx, 294, bw, 250, { stroke: TEAL });
  T('Read preparation', bx + 16, 304, { fontSize: HEAD, bold: true, color: TEAL });
  const prepSteps = [
    ['Length and expected-error gates', 'over-length reads and reads averaging > 0.01 errors per base'],
    ['Competitive alignment', 'minimap2 against target and paralog references; chimeras out'],
    ['Paralog-allele filter', 'reads with ≥ 10 paralog-specific bases removed'],
    ['Matched depth', 'random draw of the same number of usable reads per library'],
  ];
  prepSteps.forEach(([t, s], k) => {
    const y = 344 + k * 48;
    badge(k + 1, bx + 30, y + 16, TEAL);
    T(t, bx + 56, y, { bold: true, fontSize: SM });
    T(s, bx + 56, y + 22, { fontSize: SM, color: MID });
  });
  inCards.forEach((k) => LINK(k, 's', prep, 'n', { color: MID, style: 'elbow', headSize: 0.6 }));

  // ───────── C. Calling and checks ─────────
  const cy = 566;
  panel('C', 'Calling and quality checks', 30, cy);
  const csteps = [
    ['Control error model', ['per substitution, strand and', 'base quality; a margin chosen by', 'leave-one-control-out']],
    ['Test and FDR', ['support against expected errors;', 'FDR ≤ 0.05 set from the controls;', 'each control called vs the others']],
    ['Germline and window', ['sites at ≥ 50% removed; calls', 'restricted to the covered window;', 'indels optional']],
    ['Automatic checks', ['spectrum (C→T/G→A share),', 'positional clustering, 5-mer', 'targeting index, Sanger', 'bookkeeping flags']],
  ];
  const cw = (W - 60 - 3 * 18) / 4, ccards = [];
  csteps.forEach(([t, sub], k) => {
    const x = 30 + k * (cw + 18);
    const c = card(x, cy + 38, cw, 152, { stroke: k === 3 ? SANG : TEAL });
    ccards.push(c);
    badge(k + 1, x + 24, cy + 62, k === 3 ? SANG : TEAL);
    T(t, x + 46, cy + 50, { bold: true });
    lines(sub, x + 14, cy + 90);
  });
  ccards.slice(1).forEach((c, k) => LINK(ccards[k], 'e', c, 'w', { headSize: 0.7 }));
  // (no B-to-C connector: an elbow from B crosses the panel C title)
  // LINK(prep, 's', ccards[0], 'n', { color: MID, style: 'elbow', headSize: 0.6 });

  // ───────── D. Outputs ─────────
  const dy = cy + 214;
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
  // repair indices: stacked spectrum bar
  x = tile(4, ['Repair-route', 'indices'], ['G:C Ts / Tv, A:T,', 'WA/TW, decomposition']);
  [[0.35, '#2c5d8a'], [0.25, '#7fa7cf'], [0.4, '#d9952b']].reduce((acc, [f, c]) => { R(x + 40, dy + 104 + acc * 80, 50, f * 80, { fill: c }); return acc + f; }, 0);
  [[0.66, '#2c5d8a'], [0.29, '#7fa7cf'], [0.05, '#d9952b']].reduce((acc, [f, c]) => { R(x + 110, dy + 104 + acc * 80, 50, f * 80, { fill: c }); return acc + f; }, 0);
  // reports and interfaces
  x = tile(5, ['Reports, app', 'and command line'], ['HTML report, tables,', 'project files']);
  o.push(Make.image(LOGO, x + tw / 2 - 34, dy + 102, 68, 68));
  return { id: 'crm-fig1', name: 'Figure 1', width: W, height: dy + 280, background: '#ffffff', objects: o,
    alt: 'Figure 1. SHMtrace: control-calibrated somatic hypermutation calling, from inputs and read preparation through calling, quality checks and outputs.' };
}
