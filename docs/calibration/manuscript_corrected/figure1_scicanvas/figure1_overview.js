() => {
  // Figure 1 (overview): samples -> sequencing -> SHMtrace -> readouts, read left to right.
  // Cell Press production rules (checked by build.js): 174 mm wide (1,400 units, 1 unit = 0.352 pt),
  // Arial, text 6-8 pt (17.5-22.5 units), strokes 0.5-1.5 pt (>= 1.5 units), no red with green.
  // Colours: cases blue and AID-null controls red, as in Figures 2 and 4; lanes blue (J_H4) and
  // purple (N114); Sanger ochre; SHMtrace teal.
  const W = 1400, H = 616, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8';
  const JH4 = '#2c5d8a', N114 = '#6a4c93', SANG = '#a3621b', TEAL = '#1f6f8b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const F = 'arial', TITLE = 22.5, HEAD = 21, LBL = 19.5, SM = 17.5;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const TR = (text, rx, y, e = {}) => { const t = T(text, rx, y, e); t.x = rx - t.w; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const A = (x1, y1, x2, y2, e = {}) => { const c = Make.connector({ x: x1, y: y1 }, { x: x2, y: y2 }, { color: INK, width: 2.4, headSize: 0.9, ...e }); o.push(c); return c; };
  const L = (x1, y1, x2, y2, e = {}) => A(x1, y1, x2, y2, { head: 'none', ...e });
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const badge = (n, cx, cy, fill) => o.push(Make.badge(n, cx, cy, { fill, labelSize: SM, w: 28, h: 28, x: cx - 14, y: cy - 14 }));
  const stage = (n, title, x, w) => {
    badge(n, x + 14, 26, INK);
    T(title, x + 36, 15, { fontSize: TITLE, bold: true });
    L(x, 50, x + w, 50, { color: OUTL, width: 1.6 });
  };

  const c1 = 30, w1 = 392, c2 = 452, w2 = 220, c3 = 702, w3 = 392, c4 = 1124, w4 = 246;
  const top1 = 64, h1 = 256, top2 = 334, h2 = 256;

  // ───────── 1. Samples ─────────
  stage(1, 'Samples', c1, w1);
  card(c1, top1, w1, h1);
  T('Mouse germinal-centre B cells', c1 + 16, top1 + 12, { fontSize: HEAD, bold: true, color: JH4 });
  I('s-mouse-c57bl-6-black-realistic', c1 + 52, top1 + 78, 70, { flipX: true });
  T('{i|Wt}', c1 + 92, top1 + 56, { bold: true, color: CASE }); T('n = 3', c1 + 122, top1 + 57, { fontSize: SM, color: MID });
  I('s-mouse-c57bl-6-black-realistic', c1 + 52, top1 + 146, 70, { flipX: true });
  T('{i|Aicda}^{−/−}', c1 + 92, top1 + 122, { bold: true, color: CTRL }); T('n = 3', c1 + 92, top1 + 146, { fontSize: SM, color: MID });
  T('AID-null control', c1 + 92, top1 + 166, { fontSize: SM, color: CTRL });
  A(c1 + 200, top1 + 112, c1 + 218, top1 + 112, { color: MID, width: 2 });
  I('s-germinal-centre-b-cells-group', c1 + 262, top1 + 102, 70);
  TC('GC B cells', c1 + 262, top1 + 142, { fontSize: SM, color: MID });
  A(c1 + 302, top1 + 112, c1 + 320, top1 + 112, { color: MID, width: 2 });
  I('s-thermal-cycler-lid-open', c1 + 356, top1 + 100, 68);
  TC('PCR', c1 + 356, top1 + 142, { fontSize: SM, color: MID });
  // amplicons drawn to one scale (0.27 units per bp)
  T('{i|J_{H}4} intron amplicon · 1.23 kb', c1 + 16, top1 + 196, { fontSize: SM, bold: true, color: JH4 });
  R(c1 + 16, top1 + 222, 1230 * 0.27, 14, { fill: JH4, radius: 3 });

  card(c1, top2, w1, h2);
  T('N114 hybridoma cell line', c1 + 16, top2 + 12, { fontSize: HEAD, bold: true, color: N114 });
  I('s-retroviral-plasmid-insert-gfp-puror', c1 + 42, top2 + 72, 46);
  T('pMX-{i|Aicda}', c1 + 72, top2 + 60, { bold: true, color: CASE });
  I('s-retroviral-plasmid-empty-gfp-puror', c1 + 42, top2 + 136, 46);
  T('Empty vector', c1 + 72, top2 + 122, { bold: true, color: CTRL });
  T('AID-null control', c1 + 72, top2 + 146, { fontSize: SM, color: CTRL });
  A(c1 + 200, top2 + 104, c1 + 218, top2 + 104, { color: MID, width: 2 });
  I('s-hybridoma-cell', c1 + 262, top2 + 94, 60);
  TC('3 subclones', c1 + 262, top2 + 130, { fontSize: SM, color: MID });
  TC('per vector', c1 + 262, top2 + 150, { fontSize: SM, color: MID });
  A(c1 + 302, top2 + 104, c1 + 320, top2 + 104, { color: MID, width: 2 });
  I('s-thermal-cycler-lid-open', c1 + 356, top2 + 92, 68);
  TC('PCR', c1 + 356, top2 + 130, { fontSize: SM, color: MID });
  T('V-region amplicon · 0.65 kb', c1 + 16, top2 + 196, { fontSize: SM, bold: true, color: N114 });
  R(c1 + 16, top2 + 222, 655 * 0.27, 14, { fill: N114, radius: 3 });

  // ───────── 2. Sequencing ─────────
  stage(2, 'Sequencing', c2, w2);
  const sY = top1, sH = 184;
  card(c2, sY, w2, sH, { stroke: SANG, dashStyle: 'dashed' });
  I('s-sanger-chromatogram-trace-only', c2 + w2 / 2, sY + 38, 86);
  TC('Sanger benchmark', c2 + w2 / 2, sY + 74, { bold: true, color: SANG });
  TC('cloned amplicons from', c2 + w2 / 2, sY + 100, { fontSize: SM, color: MID });
  TC('the same mice', c2 + w2 / 2, sY + 120, { fontSize: SM, color: MID });
  TC('26–40 clones per mouse', c2 + w2 / 2, sY + 146, { fontSize: SM, color: MID });
  A(c1 + w1 + 4, sY + 112, c2 - 4, sY + 112, { color: SANG, width: 2 });
  const oY = sY + sH + 18, oH = top2 + h2 - oY;
  card(c2, oY, w2, oH);
  const oc = oY + oH / 2 - 108;
  I('s-nanopore-sequencer-benchtop-promethion-style', c2 + w2 / 2, oc + 50, 96);
  TC('Oxford Nanopore', c2 + w2 / 2, oc + 106, { bold: true });
  ['full-length reads,', 'no cloning', 'R10.4.1 flow cells, Kit 14', '500 reads per library'].forEach((s, k) =>
    TC(s, c2 + w2 / 2, oc + 134 + k * 22, { fontSize: SM, color: MID }));
  A(c1 + w1 + 4, top1 + 229, c2 - 4, oY + 40);
  A(c1 + w1 + 4, top2 + 104, c2 - 4, top2 + 104);

  // ───────── 3. SHMtrace ─────────
  stage(3, 'Analysis with SHMtrace', c3, w3);
  o.push(Make.image(LOGO, c3, 62, 50, 50));
  T('SHMtrace', c3 + 60, 64, { fontSize: TITLE, bold: true, color: TEAL });
  T('open source · app and command line', c3 + 60, 92, { fontSize: SM, color: MID });
  card(c3, 124, w3, 64, { stroke: SANG, dashStyle: 'dashed' });
  T('Sanger clones', c3 + 14, 132, { bold: true, color: SANG });
  T('same steps; confident base calls only', c3 + 14, 158, { fontSize: SM, color: MID });
  A(c2 + w2 + 4, 156, c3 - 4, 156, { color: SANG, width: 2 });
  const steps = [
    ['Keep the right reads', ['remove reads from a related gene copy,', 'joined (chimeric) reads and', 'low-quality reads']],
    ['Separate mutations from errors', ['keep a change only if it appears more', 'often than sequencing errors do in the', 'AID-null controls (at most 5% of calls', 'expected to be errors)']],
    ['Count mice, not reads', ['one value per mouse or subclone']],
  ];
  let sy = 202;
  steps.forEach(([title, sub], i) => {
    const h = i === steps.length - 1 ? top2 + h2 - sy : 38 + sub.length * 21 + 10;   // last card ends with the columns
    card(c3, sy, w3, h, { stroke: TEAL });
    badge(i + 1, c3 + 22, sy + 24, TEAL);
    T(title, c3 + 44, sy + 13, { bold: true });
    sub.forEach((s, k) => T(s, c3 + 44, sy + 40 + k * 21, { fontSize: SM, color: MID }));
    if (i < steps.length - 1) A(c3 + w3 / 2, sy + h + 1, c3 + w3 / 2, sy + h + 13, { color: TEAL, width: 2, headSize: 0.7 });
    sy += h + 14;
  });
  A(c2 + w2 + 4, oY + 110, c3 - 4, oY + 110);

  // ───────── 4. Readouts ─────────
  stage(4, 'Readouts', c4, w4);
  A(c3 + w3 + 4, oY + 110, c4 - 4, oY + 110);
  const tileH = 126, gap = 8;
  const tile = (i, title, fig) => {
    const y = top1 + i * (tileH + gap);
    card(c4, y, w4, tileH);
    T(title, c4 + 12, y + 9, { bold: true });
    T(fig, c4 + 12, y + 33, { fontSize: SM, color: SOFT });
    return y + 54;
  };
  // AID dependence: case group high, control group low
  let y = tile(0, 'AID dependence', 'Figs 2, 4');
  [[c4 + 56, [y + 8, y + 16, y + 26], CASEF, CASE, 'case'], [c4 + 142, [y + 56, y + 60, y + 62], CTRLF, CTRL, 'control']].forEach(([cx, ys, fill, edge, lab]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 12, yy, 11, { fill, stroke: edge, strokeWidth: 1.5 }));
    const m = ys.reduce((a, b) => a + b) / 3;
    L(cx - 24, m, cx + 24, m, { color: edge, width: 2.4 });
    T(lab, cx + 32, m - 11, { fontSize: SM, color: edge });
  });
  // mutation spectrum: a 4 x 4 From -> To matrix
  y = tile(1, 'Mutation spectrum', 'Fig. 3');
  const shades = [[0, 3, 2, 4], [3, 0, 4, 1], [1, 2, 0, 1], [4, 1, 2, 0]];
  const ramp = ['#ffffff', '#cddff0', '#9fc0df', '#5f93c4', '#2c5d8a'];
  shades.forEach((row, r) => row.forEach((v, c) => R(c4 + 20 + c * 16, y + 2 + r * 16, 15, 15, { fill: r === c ? '#e8ecf0' : ramp[v] })));
  T('Sanger vs ONT', c4 + 98, y + 8, { fontSize: SM, color: MID });
  T('r = 0.91', c4 + 98, y + 32, { fontSize: SM, bold: true, color: TEAL });
  // full-length map: mutations thin out along the 1.23-kb intron
  y = tile(2, 'Full-length map', 'Fig. 2C');
  const pts = [[0, 40], [20, 8], [40, 20], [60, 6], [80, 18], [100, 16], [120, 34], [140, 32], [160, 48], [180, 54], [204, 60]];
  pts.slice(1).forEach(([px, py], k) => L(c4 + 22 + pts[k][0], y + 0.72 * pts[k][1], c4 + 22 + px, y + 0.72 * py, { color: JH4, width: 2.2 }));
  L(c4 + 22, y + 48, c4 + 226, y + 48, { color: INK, width: 1.6 });
  R(c4 + 22, y + 52, 72, 5, { fill: SANG });
  T('Sanger window', c4 + 100, y + 45, { fontSize: SM, color: SANG });
  // clone trees
  y = tile(3, 'Clone trees', 'Fig. S10');
  const tx = c4 + 28, root = y + 30;
  R(tx - 6, root - 6, 12, 12, { fill: INK });
  L(tx, root, tx + 46, root, { color: INK, width: 2 });
  [[y + 4, 104], [y + 22, 150], [y + 40, 86], [y + 58, 124]].forEach(([yy, len], k) => {
    L(tx + 46, root, tx + 46, yy, { color: SOFT, width: 1.6 });
    L(tx + 46, yy, tx + 46 + len, yy, { color: INK, width: 2 });
    for (let j = 1; j <= 3; j++) R(tx + 46 + (len * j) / 4 - 1.5, yy - 6, 3, 12, { fill: ['#c0392b', '#e67e22', '#2471a3'][(j + k) % 3] });
  });
  return { id: 'fig1-overview', name: 'Figure 1', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure 1. Study overview: two biological systems, nanopore and Sanger sequencing, analysis with SHMtrace, and the readouts shown in later figures.' };
}
