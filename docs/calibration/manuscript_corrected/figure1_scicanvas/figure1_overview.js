() => {
  // Figure 1 (overview): samples -> sequencing -> SHMtrace -> readouts, read left to right.
  // Cell Press production rules (checked by build.js): 174 mm wide (1,400 units, 1 unit = 0.352 pt),
  // Arial, text 6-8 pt (17.5-22.5 units), strokes 0.5-1.5 pt, no red with green. Cases blue and
  // AID-null controls red, as in Figures 2 and 4; J_H4 blue, N114 purple, Sanger ochre, SHMtrace teal.
  // Layout: cards in each column share their neighbours' centre lines, so flow arrows are attached to
  // the cards they join (they follow the cards when edited in SciCanvas) and run straight or at right
  // angles. Decorative lines (axes, means, branches) are thin shapes, not connectors.
  const W = 1400, H = 616, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8';
  const JH4 = '#2c5d8a', N114 = '#6a4c93', SANG = '#a3621b', TEAL = '#1f6f8b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const F = 'arial', TITLE = 22.5, HEAD = 21, LBL = 19.5, SM = 17.5;
  const GFP_TEAL = { '#3fae5a': '#16938f', '#bfe8c3': '#bfe3e1' };   // GFP drawn teal: no red-green pairs

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  // a decorative line: a thin rotated rectangle
  const SEG = (x1, y1, x2, y2, w, fill) => {
    const len = Math.hypot(x2 - x1, y2 - y1), cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    return R(cx - len / 2, cy - w / 2, len, w, { fill, radius: w / 2, rot: (Math.atan2(y2 - y1, x2 - x1) * 180) / Math.PI });
  };
  // a flow arrow attached to two objects (port n, e, s or w)
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: INK, width: 2.4, headSize: 0.9, ...e }); o.push(c); return c; };
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, { colorMap: GFP_TEAL, ...e }); if (GFP_TEAL[ic.color]) ic.color = GFP_TEAL[ic.color]; ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const badge = (n, cx, cy, fill) => { const b = Make.badge(n, cx, cy, { fill, labelSize: SM, w: 34, h: 34, x: cx - 17, y: cy - 17 }); o.push(b); return b; };
  const stage = (n, title, x, w) => { badge(n, x + 17, 30, INK); T(title, x + 42, 18, { fontSize: TITLE, bold: true }); SEG(x, 52, x + w, 52, 2, OUTL); };

  const c1 = 30, w1 = 392, c2 = 452, w2 = 220, c3 = 702, w3 = 392, c4 = 1124, w4 = 246;
  const top1 = 66, h1 = 254, top2 = 334, h2 = 254;             // centre lines 193 and 461

  // ───────── 1. Samples ─────────
  stage(1, 'Samples', c1, w1);
  const k1 = card(c1, top1, w1, h1);
  T('Mouse germinal-centre B cells', c1 + 16, top1 + 12, { fontSize: HEAD, bold: true, color: JH4 });
  I('s-mouse-c57bl-6-black-realistic', c1 + 52, top1 + 76, 70, { flipX: true });
  T('{i|Wt}', c1 + 92, top1 + 54, { bold: true, color: CASE }); T('n = 3', c1 + 122, top1 + 55, { fontSize: SM, color: MID });
  I('s-mouse-c57bl-6-black-realistic', c1 + 52, top1 + 146, 70, { flipX: true });
  T('{i|Aicda}^{−/−}', c1 + 92, top1 + 116, { bold: true, color: CTRL });
  T('n = 3', c1 + 92, top1 + 144, { fontSize: SM, color: MID });
  T('AID-null control', c1 + 92, top1 + 168, { fontSize: SM, color: CTRL });
  const gc = I('s-germinal-centre-b-cells-group', c1 + 262, top1 + 100, 70);
  TC('GC B cells', c1 + 262, top1 + 142, { fontSize: SM, color: MID });
  const pcr1 = I('s-thermal-cycler-lid-open', c1 + 352, top1 + 100, 62);
  TC('PCR', c1 + 352, top1 + 142, { fontSize: SM, color: MID });
  LINK(gc, 'e', pcr1, 'w', { color: MID });
  T('{i|J_{H}4} intron amplicon · 1.23 kb', c1 + 16, top1 + 196, { fontSize: SM, bold: true, color: JH4 });
  R(c1 + 16, top1 + 222, 1230 * 0.27, 14, { fill: JH4, radius: 3 });     // amplicons on one scale

  const k2 = card(c1, top2, w1, h2);
  T('N114 hybridoma cell line', c1 + 16, top2 + 12, { fontSize: HEAD, bold: true, color: N114 });
  I('s-retroviral-plasmid-insert-gfp-puror', c1 + 42, top2 + 70, 46);
  T('pMX-{i|Aicda}', c1 + 72, top2 + 58, { bold: true, color: CASE });
  I('s-retroviral-plasmid-empty-gfp-puror', c1 + 42, top2 + 134, 46);
  T('Empty vector', c1 + 72, top2 + 120, { bold: true, color: CTRL });
  T('AID-null control', c1 + 72, top2 + 146, { fontSize: SM, color: CTRL });
  const hyb = I('s-hybridoma-cell', c1 + 262, top2 + 92, 60);
  TC('3 subclones', c1 + 262, top2 + 128, { fontSize: SM, color: MID });
  TC('per vector', c1 + 262, top2 + 152, { fontSize: SM, color: MID });
  const pcr2 = I('s-thermal-cycler-lid-open', c1 + 352, top2 + 92, 62);
  TC('PCR', c1 + 352, top2 + 128, { fontSize: SM, color: MID });
  LINK(hyb, 'e', pcr2, 'w', { color: MID });
  T('V-region amplicon · 0.65 kb', c1 + 16, top2 + 196, { fontSize: SM, bold: true, color: N114 });
  R(c1 + 16, top2 + 222, 655 * 0.27, 14, { fill: N114, radius: 3 });

  // ───────── 2. Sequencing ─────────
  stage(2, 'Sequencing', c2, w2);
  const sg = card(c2, top1, w2, h1, { stroke: SANG, dashStyle: 'dashed' });
  I('s-sanger-chromatogram-trace-only', c2 + w2 / 2, top1 + 74, 86);
  TC('Sanger benchmark', c2 + w2 / 2, top1 + 116, { bold: true, color: SANG });
  TC('cloned amplicons from', c2 + w2 / 2, top1 + 146, { fontSize: SM, color: MID });
  TC('the same mice;', c2 + w2 / 2, top1 + 170, { fontSize: SM, color: MID });
  TC('26–40 clones per mouse', c2 + w2 / 2, top1 + 194, { fontSize: SM, color: MID });
  const ont = card(c2, top2, w2, h2);
  I('s-nanopore-sequencing-strand-through-pore', c2 + w2 / 2, top2 + 70, 92);
  TC('Oxford Nanopore', c2 + w2 / 2, top2 + 122, { bold: true });
  ['full-length amplicon reads,', 'no cloning', 'R10.4.1 flow cells, Kit 14'].forEach((s, k) =>
    TC(s, c2 + w2 / 2, top2 + 152 + k * 24, { fontSize: SM, color: MID }));
  LINK(k1, 'e', sg, 'w', { color: SANG });
  LINK(k1, 'e', ont, 'w', { style: 'elbow' });
  LINK(k2, 'e', ont, 'w');

  // ───────── 3. SHMtrace ─────────
  stage(3, 'Analysis with SHMtrace', c3, w3);
  o.push(Make.image(LOGO, c3, 66, 46, 46));
  T('SHMtrace', c3 + 58, 66, { fontSize: TITLE, bold: true, color: TEAL });
  T('open source · app and command line', c3 + 58, 96, { fontSize: SM, color: MID });
  const sm = card(c3, 161, w3, 64, { stroke: SANG, dashStyle: 'dashed' });       // centre 193
  T('Sanger mode', c3 + 14, 168, { bold: true, color: SANG });
  T('trimmed traces; substitutions at Q ≥ 20', c3 + 14, 194, { fontSize: SM, color: MID });
  LINK(sg, 'e', sm, 'w', { color: SANG });
  const steps = [
    ['Target-read selection', ['removes over-length, low-quality,', 'chimeric and paralogous reads']],
    ['Control-calibrated variant calling', ['each substitution is tested against the', 'error rate measured in the AID-null',
      'controls (Benjamini–Hochberg FDR ≤ 0.05)']],
    ['Animal-level statistics', ['one value per mouse or subclone;', 'reads are not treated as replicates']],
  ];
  let sy = 239;
  const cards = steps.map(([title, sub], i) => {
    const h = i === steps.length - 1 ? top2 + h2 - sy : 44 + sub.length * 23 + 10;
    const k = card(c3, sy, w3, h, { stroke: TEAL });
    badge(i + 1, c3 + 24, sy + 26, TEAL);
    T(title, c3 + 48, sy + 14, { bold: true });
    sub.forEach((s, j) => T(s, c3 + 48, sy + 42 + j * 23, { fontSize: SM, color: MID }));
    sy += h + 14;
    return k;
  });
  cards.slice(1).forEach((k, i) => LINK(cards[i], 's', k, 'n', { color: TEAL, headSize: 0.7 }));
  LINK(ont, 'e', cards[0], 'w', { style: 'elbow' });

  // ───────── 4. Readouts ─────────
  stage(4, 'Readouts', c4, w4);
  const tileH = 124, gap = 10, tiles = [];
  const tile = (i, title, fig) => {
    const y = top1 + i * (tileH + gap);
    tiles.push(card(c4, y, w4, tileH));
    T(title, c4 + 12, y + 9, { bold: true });
    T(fig, c4 + 12, y + 33, { fontSize: SM, color: SOFT });
    return y + 56;
  };
  // AID-dependent SHM: case group high, control group low
  let y = tile(0, 'AID-dependent SHM', 'Figures 2 and 4');
  [[c4 + 52, [y + 4, y + 12, y + 22], CASEF, CASE, 'case'], [c4 + 126, [y + 48, y + 52, y + 54], CTRLF, CTRL, 'control']].forEach(([cx, ys, fill, edge, lab]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 12, yy, 11, { fill, stroke: edge, strokeWidth: 1.5 }));
    const m = ys.reduce((a, b) => a + b) / 3;
    SEG(cx - 24, m, cx + 24, m, 2.4, edge);
    T(lab, cx + 30, m - 11, { fontSize: SM, color: edge });
  });
  // substitution spectrum: a 4 x 4 From -> To matrix
  y = tile(1, 'Substitution spectrum', 'Figure 3');
  const shades = [[0, 3, 2, 4], [3, 0, 4, 1], [1, 2, 0, 1], [4, 1, 2, 0]];
  const ramp = ['#ffffff', '#cddff0', '#9fc0df', '#5f93c4', '#2c5d8a'];
  shades.forEach((row, r) => row.forEach((v, c) => R(c4 + 20 + c * 15, y - 2 + r * 15, 14, 14, { fill: r === c ? '#e8ecf0' : ramp[v] })));
  T('Sanger vs ONT', c4 + 92, y + 2, { fontSize: SM, color: MID });
  T('Pearson r = 0.93', c4 + 92, y + 26, { fontSize: SM, bold: true, color: TEAL });
  // full-length map: mutations thin out along the 1.23-kb intron
  y = tile(2, 'Full-length map', 'Figure 2C');
  const pts = [[0, 30], [20, 6], [40, 15], [60, 4], [80, 13], [100, 12], [120, 25], [140, 23], [160, 35], [180, 39], [204, 44]];
  pts.slice(1).forEach(([px, py], k) => SEG(c4 + 22 + pts[k][0], y - 6 + 0.85 * pts[k][1], c4 + 22 + px, y - 6 + 0.85 * py, 2.2, JH4));
  SEG(c4 + 22, y + 34, c4 + 226, y + 34, 1.8, INK);
  R(c4 + 22, y + 39, 56, 5, { fill: SANG });
  T('Sanger window', c4 + 84, y + 32, { fontSize: SM, color: SANG });
  // clonal lineage trees
  y = tile(3, 'Clonal lineage trees', 'Supplementary Figure S10');
  const tx = c4 + 28, root = y + 26;
  R(tx - 6, root - 6, 12, 12, { fill: INK });
  SEG(tx, root, tx + 46, root, 2, INK);
  [[y + 4, 104], [y + 19, 150], [y + 34, 86], [y + 49, 124]].forEach(([yy, len], k) => {
    SEG(tx + 46, root, tx + 46, yy, 1.6, SOFT);
    SEG(tx + 46, yy, tx + 46 + len, yy, 2, INK);
    for (let j = 1; j <= 3; j++) R(tx + 46 + (len * j) / 4 - 1.5, yy - 5, 3, 10, { fill: ['#c0392b', '#e67e22', '#2471a3'][(j + k) % 3] });
  });
  LINK(cards[1], 'e', tiles[2], 'w', { style: 'elbow' });
  return { id: 'fig1-overview', name: 'Figure 1', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure 1. Study overview: two biological systems, nanopore and Sanger sequencing, analysis with SHMtrace, and the readouts shown in later figures.' };
}
