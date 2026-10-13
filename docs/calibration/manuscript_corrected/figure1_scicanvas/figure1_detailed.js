() => {
  // Detailed workflow (supplementary): design, preprocessing, SHMtrace calling, Sanger route, reads.
  // Cell Press production rules (checked by build.js): 174 mm wide (1,400 units, 1 unit = 0.352 pt), at most
  // 200 mm tall, Arial, text 6-8 pt (17.5-22.5 units), strokes 0.5-1.5 pt, no red with green. Cases blue and
  // AID-null controls red, as in Figures 2 and 4; J_H4 lanes blue, N114 purple, Sanger ochre, SHMtrace teal.
  // Numbers: amplicon coordinates and windows (Methods 2.7-2.8); delivered and passing reads from the
  // deposit's pipeline_summary_study_libraries.csv and read counts; leave-one-control-out factor 5
  // (project.json "calibration"); paralogue filter contaminant_min_alleles = 10.
  const W = 1400, o = [];
  const INK = '#1d2733', MID = '#4a5562', OUTL = '#b4bec8';
  const JH4 = '#2c5d8a', N114 = '#6a4c93', SANG = '#a3621b', TEAL = '#1f6f8b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const TGT = '#5b8db8', CON = '#d9952b', ERR = '#c0392b', VAR = '#6a4c93', GFPC = MID;
  const F = 'arial', TITLE = 22.5, LBL = 19.5, SM = 17.5, HEAD = 21;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const A = (x1, y1, x2, y2, e = {}) => { const c = Make.connector({ x: x1, y: y1 }, { x: x2, y: y2 }, { color: MID, width: 2, headSize: 0.9, ...e }); o.push(c); return c; };
  const L = (x1, y1, x2, y2, e = {}) => A(x1, y1, x2, y2, { head: 'none', ...e });
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const panel = (letter, title, y) => { T(letter, 30, y, { fontSize: TITLE, bold: true }); T(title, 58, y, { fontSize: TITLE, bold: true }); };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const read = (x, y, w, c, h = 9) => R(x, y - h / 2, w, h, { fill: c, radius: h / 2 });
  const tick = (x, y, ok) => TC(ok ? '✓' : '✕', x, y - 11, { fontSize: 20, bold: true, color: ok ? CASE : ERR });
  const bracket = (x1, x2, y, c, label, sub) => {
    L(x1, y, x2, y, { color: c, width: 2, head: 'bar', tail: 'bar', headSize: 0.5 });
    const t = T(label, x1, y + 5, { fontSize: SM, color: c, bold: true });
    if (sub) T(sub, t.x + t.w + 6, y + 5, { fontSize: SM, color: MID });
  };
  const BP = 0.27;                                          // both amplicons on one scale
  // a row of numbered step cards; returns the x of each card
  const cards = (steps, y, h, colour) => {
    const n = steps.length, gap = 16, cw = (W - 60 - gap * (n - 1)) / n;
    return steps.map(([num, title, rules], i) => {
      const x = 30 + i * (cw + gap);
      card(x, y, cw, h);
      o.push(Make.badge(num, x + 24, y + 24, { fill: colour, labelSize: SM, w: 28, h: 28, x: x + 10, y: y + 10 }));
      T(title, x + 46, y + 13, { bold: true });
      rules.forEach((r, k) => T(r, x + 14, y + h - 26 - (rules.length - 1 - k) * 21, { fontSize: SM, color: k === 0 ? colour : MID, bold: k === 0 }));
      if (i < n - 1) A(x + cw + 2, y + h / 2, x + cw + gap - 2, y + h / 2, { color: INK, headSize: 0.7 });
      return { x, w: cw };
    });
  };

  // ───────── A. Experimental design ─────────
  panel('A', 'Experimental design', 10);
  const l1 = 46, l1h = 226, r1 = l1 + 120;
  card(30, l1, W - 60, l1h);
  T('{i|J_{H}4} intron · mouse germinal-centre B cells ({i|in vivo})', 46, l1 + 12, { fontSize: HEAD, bold: true, color: JH4 });
  I('s-mouse-c57bl-6-black-realistic', 86, r1 - 36, 74, { flipX: true });
  T('{i|Wt}', 132, r1 - 58, { bold: true, color: CASE }); T('n = 3', 132, r1 - 34, { fontSize: SM, color: MID });
  I('s-mouse-c57bl-6-black-realistic', 86, r1 + 38, 74, { flipX: true });
  T('{i|Aicda}^{−/−}', 132, r1 + 14, { bold: true, color: CTRL }); T('n = 3', 132, r1 + 36, { fontSize: SM, color: MID }); T('AID-null control', 132, r1 + 56, { fontSize: SM, color: MID });
  T('4-month-old, C57BL/6', 46, r1 + 84, { fontSize: SM, color: MID });
  A(300, r1, 320, r1);
  I('s-peyer-s-patch', 360, r1 - 8, 82); TC("Peyer's patch", 360, r1 + 36, { fontSize: SM });
  A(424, r1, 444, r1);
  I('s-facs-sorter-nozzle-laser-plates-tubes', 494, r1 - 10, 74);
  TC('sort GC B cells', 494, r1 + 36, { fontSize: SM }); TC('B220^{+}GL7^{+}Fas^{+}', 494, r1 + 58, { fontSize: SM, color: MID });
  A(560, r1, 578, r1);
  I('s-thermal-cycler-lid-open', 616, r1 - 12, 66); TC('DNA, PCR', 616, r1 + 36, { fontSize: SM }); TC('SuperFi II', 616, r1 + 58, { fontSize: SM, color: MID });
  A(656, r1, 672, r1);
  const ax = 684, ay = r1 - 30;
  T('{i|J_{H}4} intron amplicon · 1,230 bp', ax, ay - 42, { bold: true });
  R(ax, ay - 7, 1230 * BP, 14, { fill: JH4, radius: 3 });
  bracket(ax + 89 * BP, ax + 1230 * BP, ay + 24, JH4, '1,141 bp', 'nanopore, positions 90–1,230');
  bracket(ax + 100 * BP, ax + 629 * BP, ay + 66, SANG, '529 bp', 'shared with Sanger, 101–629');
  A(1028, r1 - 26, 1056, r1 - 52); A(1028, r1 + 18, 1056, r1 + 44);
  I('s-nanopore-sequencer-benchtop-promethion-style', 1100, r1 - 58, 62);
  T('Nanopore', 1142, r1 - 80, { bold: true }); T('R10.4.1, Kit 14', 1142, r1 - 56, { fontSize: SM, color: MID });
  I('s-sanger-chromatogram-trace-only', 1100, r1 + 54, 70);
  T('Sanger', 1142, r1 + 32, { bold: true, color: SANG }); T('cloned amplicons,', 1142, r1 + 56, { fontSize: SM, color: MID });
  T('separate PCRs', 1142, r1 + 76, { fontSize: SM, color: MID });

  const l2 = l1 + l1h + 10, l2h = 214, r2 = l2 + 112;
  card(30, l2, W - 60, l2h);
  T('N114 hybridoma · fixed V-region substrate ({i|in vitro})', 46, l2 + 12, { fontSize: HEAD, bold: true, color: N114 });
  I('s-retroviral-plasmid-insert-gfp-puror', 68, r2 - 28, 46);
  T('pMX-{i|Aicda}', 100, r2 - 40, { bold: true, color: CASE });
  I('s-retroviral-plasmid-empty-gfp-puror', 68, r2 + 28, 46);
  T('Empty vector', 100, r2 + 16, { bold: true, color: CTRL });
  T('both carry GFP and Puro^{R}', 46, r2 + 62, { fontSize: SM, color: GFPC });
  A(300, r2, 320, r2);
  I('s-retroviral-transduction-gfp', 360, r2 - 10, 92); TC('transduce', 360, r2 + 40, { fontSize: SM });
  A(424, r2, 444, r2);
  I('s-gfp-histogram-with-high-gate', 494, r2 - 10, 64); TC('sort GFP-high', 494, r2 + 40, { fontSize: SM });
  A(560, r2, 578, r2);
  I('s-96-well-plate-top-view-flat', 616, r2 - 30, 56); I('s-puromycin-selection-flat', 616, r2 + 10, 44);
  TC('subclone', 616, r2 + 40, { fontSize: SM }); TC('puromycin', 616, r2 + 60, { fontSize: SM, color: MID });
  A(656, r2, 672, r2);
  const nx = 684, ny = r2 - 26;
  T('V-region amplicon · 655 bp', nx, ny - 42, { bold: true });
  R(nx, ny - 7, 655 * BP, 14, { fill: N114, radius: 3 });
  R(nx, ny - 7, 20 * BP, 14, { fill: '#c9bcdc', radius: 2 }); R(nx + 635 * BP, ny - 7, 20 * BP, 14, { fill: '#c9bcdc', radius: 2 });
  bracket(nx + 20 * BP, nx + 635 * BP, ny + 24, N114, '615 bp', 'positions 21–635');
  [0, 1, 2].forEach((i) => { E(nx + 8 + i * 20, ny + 74, 15, { fill: CASEF, stroke: CASE, strokeWidth: 1.5 }); E(nx + 92 + i * 20, ny + 74, 15, { fill: CTRLF, stroke: CTRL, strokeWidth: 1.5 }); });
  T('3 subclones per vector', nx + 160, ny + 64, { fontSize: SM, color: MID });
  A(1028, r2, 1056, r2);
  I('s-nanopore-sequencer-benchtop-promethion-style', 1100, r2 - 4, 62);
  T('Nanopore', 1142, r2 - 26, { bold: true }); T('R10.4.1, Kit 14', 1142, r2 - 2, { fontSize: SM, color: MID });
  L(1142, l2 + l2h - 24, 1142 + 200 * BP, l2 + l2h - 24, { color: INK, width: 2 });
  T('200 bp', 1142 + 200 * BP + 8, l2 + l2h - 35, { fontSize: SM, color: MID });

  // ───────── B. Read preprocessing (deposited pipeline) ─────────
  const by = l2 + l2h + 16;
  panel('B', 'Preparing the reads', by);
  const bh = 184, b0 = by + 36;
  const bc = cards([
    [1, 'Keep target reads', ['closer to the target than to', 'a related copy ({i|J_{H}4}: C57BL/6', 'vs Agouti; N114: one reference)']],
    [2, 'Full-length reads only', ['{i|J_{H}4} 900–2,000 bp', 'N114 400–900 bp']],
    [3, 'Few predicted errors', ['≤ 0.01 errors per base', 'expected from the base calls']],
    [4, 'Same depth for every library', ['500 reads per library', 'one random draw, deposited']],
  ], b0, bh, JH4);
  // 1: target / contaminant / chimera
  let x = bc[0].x;
  [['target', [TGT]], ['related copy', [CON]], ['chimera', [TGT, CON]]].forEach(([lab, cols], k) => {
    const y = b0 + 56 + k * 22;
    if (cols.length === 1) read(x + 14, y, 84, cols[0]); else { read(x + 14, y, 44, cols[0]); read(x + 56, y, 42, cols[1]); }
    T(lab, x + 108, y - 10, { fontSize: SM });
    tick(x + 270, y, k === 0);
  });
  // 2: size windows on a 0-2 kb axis
  x = bc[1].x; const x0 = x + 70, k2 = 200 / 2000;
  [['{i|J_{H}4}', 900, 2000, JH4, b0 + 62], ['N114', 400, 900, N114, b0 + 92]].forEach(([lab, a, b, c, y]) => {
    T(lab, x + 14, y - 11, { fontSize: SM, bold: true, color: c });
    R(x0, y - 5, 200, 10, { fill: '#d2d8de', radius: 5 });
    R(x0 + a * k2, y - 5, (b - a) * k2, 10, { fill: c, radius: 5 });
  });
  // 3: expected-error dots
  x = bc[2].x;
  [[6, false], [1, true], [0, true]].forEach(([n, ok], k) => {
    const y = b0 + 56 + k * 22;
    read(x + 14, y, 240, TGT);
    for (let j = 0; j < n; j++) E(x + 34 + j * 36, y, 7, { fill: ERR });
    tick(x + 284, y, ok);
  });
  // 4: down-sampling
  x = bc[3].x;
  for (let j = 0; j < 8; j++) read(x + 18, b0 + 50 + j * 8, 80, TGT, 5);
  A(x + 108, b0 + 78, x + 148, b0 + 78, { color: INK });
  for (let j = 0; j < 3; j++) read(x + 160, b0 + 70 + j * 8, 80, TGT, 5);
  TC('911–2,903', x + 58, b0 + 112, { fontSize: SM, bold: true }); TC('500', x + 200, b0 + 112, { fontSize: SM, bold: true });

  // ───────── C. SHMtrace calling ─────────
  const cy = b0 + bh + 16;
  panel('C', 'Finding mutations with SHMtrace', cy);
  o.push(Make.image(LOGO, 436, cy - 4, 30, 30));
  const ch = 228, c0 = cy + 36;
  const cc = cards([
    [5, 'Line up the reads', ['align to the reference;', 'drop reads that match a', 'related gene copy']],
    [6, 'List every change', ['every base that differs;', 'a change seen in one read', 'is enough to be listed']],
    [7, 'Measure errors', ['in the AID-null controls,', 'per type of change, strand', 'and quality; 5× margin']],
    [8, 'Keep real mutations', ['seen more than errors', 'explain (≤ 5% false);', 'inherited changes dropped']],
    [9, 'Per-animal summary', ['mutations per base read;', 'spectra, clones, trees;', 'each mouse counts once']],
  ], c0, ch, TEAL);
  // 5: reads on the reference, one paralogue-like read crossed out
  x = cc[0].x; let rw = cc[0].w - 40;
  R(x + 14, c0 + 54, rw, 3, { fill: INK });
  for (let j = 0; j < 4; j++) {
    const y = c0 + 72 + j * 16;
    read(x + 14, y, rw, j === 2 ? CON : TGT, 8);
  }
  tick(x + rw + 26, c0 + 104, false);
  // 6: candidates: every mismatch marked
  x = cc[1].x; rw = cc[1].w - 32;
  R(x + 14, c0 + 54, rw, 3, { fill: INK });
  [[0.2, 0.55], [0.4], [0.2, 0.8], [0.62]].forEach((ms, j) => {
    const y = c0 + 72 + j * 16;
    read(x + 14, y, rw, '#cfd8df', 8);
    ms.forEach((m) => E(x + 14 + m * rw, y, 9, { fill: VAR }));
  });
  // 7: control libraries -> error rate table
  x = cc[2].x;
  [0, 1, 2].forEach((k) => { for (let j = 0; j < 4; j++) read(x + 14 + k * 30, c0 + 60 + j * 9, 24, CTRLF, 6); });
  A(x + 108, c0 + 74, x + 134, c0 + 74, { color: INK });
  const grid = [[1, 3, 2], [2, 1, 3], [3, 2, 1]];
  grid.forEach((row, r) => row.forEach((v, c) => R(x + 144 + c * 20, c0 + 52 + r * 16, 19, 15, { fill: ['#eef3f7', '#cfe0ec', '#9fc0d8', '#5f93b9'][v] })));
  T('error rate', x + 136, c0 + 104, { fontSize: SM, color: MID });
  // 8: observed vs expected for a candidate
  x = cc[3].x;
  const bx = x + 24, base = c0 + 118;
  L(bx - 6, base, bx + 196, base, { color: INK, width: 1.6 });
  R(bx + 10, base - 60, 34, 60, { fill: VAR }); TC('observed', bx + 27, base + 4, { fontSize: SM, color: MID });
  R(bx + 100, base - 13, 34, 13, { fill: CTRLF }); TC('errors', bx + 117, base + 4, { fontSize: SM, color: MID });
  TC('keep ✓', bx + 172, base - 52, { fontSize: SM, bold: true, color: CASE });
  // 9: per-animal dot plot
  x = cc[4].x;
  L(x + 26, c0 + 122, x + 216, c0 + 122, { color: INK, width: 1.6 });
  [[x + 76, [c0 + 56, c0 + 64, c0 + 74], CASEF, CASE], [x + 166, [c0 + 104, c0 + 108, c0 + 110], CTRLF, CTRL]].forEach(([cx, ys, fill, edge]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 12, yy, 11, { fill, stroke: edge, strokeWidth: 1.5 }));
    const m = ys.reduce((a, b) => a + b) / 3;
    L(cx - 24, m, cx + 24, m, { color: edge, width: 2.4 });
  });

  // ───────── D. Sanger route ─────────
  const dy = c0 + ch + 16;
  panel('D', 'Sanger clones take the same route', dy);
  const d0 = dy + 36, dh = 92;
  card(30, d0, W - 60, dh, { stroke: SANG, dashStyle: 'dashed' });
  const sanger = [
    ['s-sanger-chromatogram-trace-only', 'Sanger traces', '447 traces (.ab1)'],
    [null, 'Assign clones to mice', 'sample sheet, QC checks'],
    [null, 'Trim poor ends', 'align to {i|J_{H}4}'],
    [null, 'Confident calls only', 'base calls ≥ 99% accurate'],
    [null, 'Calls per mouse', 'bases each clone read'],
  ];
  const sw = (W - 60) / sanger.length;
  sanger.forEach(([icon, title, sub], i) => {
    const cx = 30 + sw * i + sw / 2;
    if (icon) I(icon, cx - 96, d0 + 46, 50);
    TC(title, cx + (icon ? 26 : 0), d0 + 22, { bold: true, color: i === 0 ? SANG : INK });
    TC(sub, cx + (icon ? 26 : 0), d0 + 50, { fontSize: SM, color: MID });
    if (i < sanger.length - 1) A(30 + sw * (i + 1) - 12, d0 + dh / 2, 30 + sw * (i + 1) + 12, d0 + dh / 2, { color: SANG });
  });

  // ───────── E. Reads per library ─────────
  const ey = d0 + dh + 16;
  panel('E', 'Reads per library', ey);
  const lgx = 420, ly = ey + 14;
  R(lgx, ly, 26, 12, { fill: '#d2d8de', radius: 2 }); T('delivered', lgx + 34, ly - 5, { fontSize: SM, color: MID });
  R(lgx + 140, ly, 13, 12, { fill: CASE }); R(lgx + 153, ly, 13, 12, { fill: CTRL }); T('passed steps 1–3', lgx + 174, ly - 5, { fontSize: SM, color: MID });
  L(lgx + 340, ly + 6, lgx + 366, ly + 6, { color: INK, width: 2, dashStyle: 'dashed' }); T('500 analysed', lgx + 374, ly - 5, { fontSize: SM, color: MID });
  const groups = [
    { title: '{i|J_{H}4} (GC B cells)', x: 30, rows: [['{i|Wt} 2', 5000, 2729, CASE], ['{i|Wt} 3', 5000, 1906, CASE], ['{i|Wt} 4', 4558, 911, CASE], ['{i|Aicda}^{−/−} 1', 5000, 1492, CTRL], ['{i|Aicda}^{−/−} 2', 5000, 1163, CTRL], ['{i|Aicda}^{−/−} 3', 5000, 1239, CTRL]] },
    { title: 'N114 (subclones)', x: 720, rows: [['pMX-{i|Aicda} C1', 4483, 2903, CASE], ['pMX-{i|Aicda} D3', 3924, 1062, CASE], ['pMX-{i|Aicda} D5', 4235, 1602, CASE], ['empty C3', 2336, 1050, CTRL], ['empty C4', 3153, 1558, CTRL], ['empty D1', 3939, 1412, CTRL]] },
  ];
  let bottom = 0;
  groups.forEach((g) => {
    const top = ey + 72, bx0 = g.x + 170, bw = 380, s = bw / 5000, rh = 27;
    T(g.title, g.x, ey + 36, { fontSize: HEAD, bold: true });
    g.rows.forEach(([lab, del, pass, c], i) => {
      const y = top + i * rh;
      const t = T(lab, 0, y - 2, { fontSize: SM }); t.x = bx0 - 10 - t.w;
      R(bx0, y, del * s, 18, { fill: '#d2d8de', radius: 2 });
      R(bx0, y, pass * s, 18, { fill: c, radius: 2 });
      T(pass.toLocaleString('en-GB'), bx0 + del * s + 8, y - 0.5, { fontSize: SM, color: MID });
    });
    const b = top + g.rows.length * rh - 6;
    L(bx0 + 500 * s, top - 8, bx0 + 500 * s, b, { color: INK, width: 2, dashStyle: 'dashed' });
    L(bx0, b + 4, bx0 + bw, b + 4, { color: INK, width: 1.6 });
    [0, 1000, 2000, 3000, 4000, 5000].forEach((v) => { L(bx0 + v * s, b + 4, bx0 + v * s, b + 10, { color: INK, width: 1.6 }); TC(v ? `${v / 1000}k` : '0', bx0 + v * s, b + 13, { fontSize: SM, color: MID }); });
    TC('reads', bx0 + bw / 2, b + 36, { fontSize: SM, color: MID });
    bottom = Math.max(bottom, b + 62);
  });
  return { id: 'fig1-detailed', name: 'Supplementary Figure S11', width: W, height: bottom, background: '#ffffff', objects: o,
    alt: 'Detailed workflow: experimental design, read preprocessing, SHMtrace calling, the Sanger route and reads per library.' };
}
