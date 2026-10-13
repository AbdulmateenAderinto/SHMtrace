() => {
  // Figure 1 (overview): samples -> sequencing -> SHMtrace -> readouts, read left to right.
  // 1,400 units = 174 mm, so 1 unit = 0.35 pt; labels 17-19 units (6-6.7 pt), headings 22-24.
  const W = 1400, H = 640, o = [];
  const INK = '#1d2733', MID = '#4f5b68', SOFT = '#8a97a6', LINE = '#c9d2dc';
  const JH4 = '#2f6690', N114 = '#2f7a5a', AID = '#c0533a', CTRL = '#8a97a6', TEAL = '#1f6f8b', SANG = '#6d5ba6';
  const F = 'sans', LBL = 19, SM = 17, HEAD = 21;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const TR = (text, rx, y, e = {}) => { const t = T(text, rx, y, e); t.x = rx - t.w; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const A = (x1, y1, x2, y2, e = {}) => { const c = Make.connector({ x: x1, y: y1 }, { x: x2, y: y2 }, { color: INK, width: 2.4, headSize: 1, ...e }); o.push(c); return c; };
  const L = (x1, y1, x2, y2, e = {}) => A(x1, y1, x2, y2, { head: 'none', ...e });
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const stage = (n, title, x, w) => {
    o.push(Make.badge(n, x + 17, 32, { fill: INK, labelSize: 17 }));
    T(title, x + 42, 20, { fontSize: 24, bold: true });
    L(x, 60, x + w, 60, { color: LINE, width: 1.5 });
  };

  // columns
  const c1 = 30, w1 = 392, c2 = 456, w2 = 236, c3 = 728, w3 = 340, c4 = 1104, w4 = 266;
  const top1 = 76, h1 = 262, top2 = 352, h2 = 262;     // the two sample cards

  // ───────── 1. Samples ─────────
  stage(1, 'Samples', c1, w1);
  R(c1, top1, w1, h1, { fill: '#f0f5f9', radius: 16 });
  T('Mouse germinal-centre B cells', c1 + 18, top1 + 14, { fontSize: HEAD, bold: true, color: JH4 });
  I('s-mouse-c57bl-6-black-realistic', c1 + 62, top1 + 86, 80, { flipX: true });
  T('{i|Wt}', c1 + 108, top1 + 62, { bold: true }); T('n = 3', c1 + 142, top1 + 63, { fontSize: SM, color: MID });
  I('s-mouse-c57bl-6-black-realistic', c1 + 62, top1 + 160, 80, { flipX: true });
  T('{i|Aicda}^{−/−}', c1 + 108, top1 + 134, { bold: true }); T('n = 3', c1 + 186, top1 + 135, { fontSize: SM, color: MID });
  T('AID-null control', c1 + 108, top1 + 160, { fontSize: SM, color: CTRL });
  A(c1 + 236, top1 + 118, c1 + 262, top1 + 118, { color: MID, width: 2 });
  I('s-germinal-centre-b-cells-group', c1 + 322, top1 + 106, 86);
  TC("Peyer's patch", c1 + 322, top1 + 152, { fontSize: SM, color: MID });
  TC('GC B cells', c1 + 322, top1 + 172, { fontSize: SM, color: MID });
  // amplicons drawn to one scale (0.27 units per bp)
  T('{i|J_{H}4} intron · 1.23 kb', c1 + 18, top1 + 206, { fontSize: SM, bold: true, color: JH4 });
  R(c1 + 18, top1 + 232, 1230 * 0.27, 14, { fill: JH4, radius: 4 });

  R(c1, top2, w1, h2, { fill: '#f0f6f2', radius: 16 });
  T('N114 hybridoma cell line', c1 + 18, top2 + 14, { fontSize: HEAD, bold: true, color: N114 });
  I('s-retroviral-plasmid-insert-gfp-puror', c1 + 50, top2 + 78, 52);
  T('pMX-{i|Aicda}', c1 + 86, top2 + 66, { bold: true, color: AID });
  I('s-retroviral-plasmid-empty-gfp-puror', c1 + 50, top2 + 146, 52);
  T('Empty vector', c1 + 86, top2 + 132, { bold: true, color: CTRL });
  T('AID-null control', c1 + 86, top2 + 158, { fontSize: SM, color: CTRL });
  A(c1 + 236, top2 + 112, c1 + 262, top2 + 112, { color: MID, width: 2 });
  I('s-hybridoma-cell', c1 + 322, top2 + 100, 72);
  TC('3 subclones', c1 + 322, top2 + 146, { fontSize: SM, color: MID });
  TC('per vector', c1 + 322, top2 + 166, { fontSize: SM, color: MID });
  T('V region · 0.65 kb', c1 + 18, top2 + 206, { fontSize: SM, bold: true, color: N114 });
  R(c1 + 18, top2 + 232, 655 * 0.27, 14, { fill: N114, radius: 4 });

  // ───────── 2. Sequencing ─────────
  stage(2, 'Sequencing', c2, w2);
  // Sanger beside the J_H4 card (same mice); ONT level with both cards
  R(c2, top1, w2, 196, { fill: '#ffffff', stroke: SANG, strokeWidth: 1.6, radius: 14, dashStyle: 'dashed' });
  I('s-sanger-chromatogram-trace-only', c2 + w2 / 2, top1 + 46, 92);
  TC('Sanger benchmark', c2 + w2 / 2, top1 + 86, { bold: true, color: SANG });
  TC('cloned amplicons from', c2 + w2 / 2, top1 + 114, { fontSize: SM, color: MID });
  TC('the same mice', c2 + w2 / 2, top1 + 134, { fontSize: SM, color: MID });
  TC('26–40 clones per {i|Wt} mouse', c2 + w2 / 2, top1 + 160, { fontSize: SM, color: MID });
  A(c1 + w1 + 4, 179, c2 - 4, 179, { color: SANG, width: 2 });
  const ontY = 300;
  R(c2, ontY, w2, 300, { fill: '#f5f7f9', radius: 14 });
  I('s-nanopore-sequencer-benchtop-promethion-style', c2 + w2 / 2, ontY + 70, 104);
  TC('Oxford Nanopore', c2 + w2 / 2, ontY + 130, { bold: true });
  ['full-length reads,', 'no cloning', 'R10.4.1 · Kit 14 · SUP', '500 reads per library'].forEach((s, k) =>
    TC(s, c2 + w2 / 2, ontY + 160 + k * 22, { fontSize: SM, color: MID }));
  A(c1 + w1 + 4, top1 + 220, c2 - 4, ontY + 40);
  A(c1 + w1 + 4, top2 + 120, c2 - 4, top2 + 120);

  // ───────── 3. SHMtrace ─────────
  stage(3, 'SHMtrace analysis', c3, w3);
  o.push(Make.image(LOGO, c3 + 4, 76, 58, 58));
  T('SHMtrace', c3 + 72, 80, { fontSize: 23, bold: true, color: TEAL });
  T('open source · app and command line', c3 + 72, 110, { fontSize: SM, color: MID });
  const steps = [
    ['Keep target reads', ['remove paralogue, chimeric and', 'low-quality reads']],
    ['Test every substitution', ['against errors measured in the', 'AID-null controls (FDR 0.05)']],
    ['Compare animals, not reads', ['one value per mouse or subclone']],
  ];
  // Sanger traces take the same route, one clone at a time
  R(c3, 148, w3, 62, { fill: '#ffffff', stroke: SANG, strokeWidth: 1.6, radius: 12, dashStyle: 'dashed' });
  T('Sanger mode', c3 + 16, 156, { bold: true, color: SANG });
  T('trimmed traces, substitutions at Q ≥ 20', c3 + 16, 182, { fontSize: SM, color: MID });
  A(c2 + w2 + 4, 179, c3 - 4, 179, { color: SANG, width: 2 });
  let sy = 228;
  steps.forEach(([title, sub], i) => {
    const h = 40 + sub.length * 22 + 14;
    R(c3, sy, w3, h, { fill: '#eef5f8', radius: 12 });
    o.push(Make.badge(i + 1, c3 + 26, sy + 28, { fill: TEAL, labelSize: 15, w: 28, h: 28, x: c3 + 12, y: sy + 14 }));
    T(title, c3 + 50, sy + 15, { bold: true });
    sub.forEach((s, k) => T(s, c3 + 50, sy + 44 + k * 22, { fontSize: SM, color: MID }));
    if (i < steps.length - 1) A(c3 + w3 / 2, sy + h + 2, c3 + w3 / 2, sy + h + 16, { color: TEAL, width: 2, headSize: 0.7 });
    sy += h + 18;
  });
  A(c2 + w2 + 4, ontY + 150, c3 - 4, ontY + 150);

  // ───────── 4. Readouts ─────────
  stage(4, 'Readouts', c4, w4);
  A(c3 + w3 + 4, ontY + 150, c4 - 4, ontY + 150);
  const tileH = 124, gap = 14;
  const tile = (i, title, fig) => {
    const y = 76 + i * (tileH + gap);
    R(c4, y, w4, tileH, { fill: '#ffffff', stroke: LINE, strokeWidth: 1.5, radius: 12 });
    T(title, c4 + 14, y + 10, { fontSize: SM + 1, bold: true });
    TR(fig, c4 + w4 - 14, y + 11, { fontSize: SM, color: SOFT });
    return y + 42;                         // top of the drawing area (~70 units tall)
  };
  // AID dependence: two groups of dots
  let y = tile(0, 'AID dependence', 'Figs 2, 4');
  [[c4 + 70, [y + 8, y + 16, y + 26], JH4, '{i|Wt}'], [c4 + 170, [y + 56, y + 60, y + 62], CTRL, 'control']].forEach(([cx, ys, col, lab]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 12, yy, 11, { fill: col }));
    const m = ys.reduce((a, b) => a + b) / 3;
    L(cx - 24, m, cx + 24, m, { color: INK, width: 2 });
    T(lab, cx + 32, m - 10, { fontSize: SM - 1, color: col });
  });
  // mutation spectrum: a 4 x 4 From -> To matrix
  y = tile(1, 'Mutation spectrum', 'Fig. 3');
  const shades = [[0, 3, 2, 4], [3, 0, 4, 1], [1, 2, 0, 1], [4, 1, 2, 0]];
  const ramp = ['#ffffff', '#dbe8f3', '#a9c6e0', '#6c9cc6', '#2f6690'];
  shades.forEach((row, r) => row.forEach((v, c) => R(c4 + 22 + c * 18, y - 4 + r * 18, 17, 17, { fill: r === c ? '#f1f3f5' : ramp[v] })));
  T('Sanger vs ONT', c4 + 110, y + 4, { fontSize: SM, color: MID });
  T('r = 0.91', c4 + 110, y + 28, { fontSize: SM, bold: true, color: TEAL });
  // full-length map: density falls along the 1.23-kb intron
  y = tile(2, 'Full-length map', 'Fig. 2C');
  const pts = [[0, 40], [20, 8], [40, 20], [60, 6], [80, 18], [100, 16], [120, 34], [140, 32], [160, 48], [180, 54], [210, 60]];
  pts.slice(1).forEach(([px, py], k) => L(c4 + 24 + pts[k][0], y + pts[k][1], c4 + 24 + px, y + py, { color: JH4, width: 2.2 }));
  L(c4 + 24, y + 64, c4 + 236, y + 64, { color: INK, width: 1.4 });
  R(c4 + 24, y + 67, 92, 5, { fill: '#b9add9' });
  T('Sanger window', c4 + 124, y + 62, { fontSize: SM - 2, color: SANG });
  // clone trees
  y = tile(3, 'Clone trees', 'Fig. S10');
  const tx = c4 + 30, root = y + 32;
  R(tx - 6, root - 6, 12, 12, { fill: INK });
  L(tx, root, tx + 50, root, { color: INK, width: 2 });
  [[y + 4, 120], [y + 24, 170], [y + 44, 100], [y + 64, 140]].forEach(([yy, len], k) => {
    L(tx + 50, root, tx + 50, yy, { color: SOFT, width: 1.4 });
    L(tx + 50, yy, tx + 50 + len, yy, { color: INK, width: 2 });
    for (let j = 1; j <= 3; j++) R(tx + 50 + (len * j) / 4 - 1.5, yy - 6, 3, 12, { fill: [AID, '#e67e22', JH4][(j + k) % 3] });
  });
  return { id: 'fig1-overview', name: 'Figure 1 (overview)', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure 1. Study overview: two biological systems, ONT and Sanger sequencing, analysis with SHMtrace, and the readouts shown in later figures.' };
}
