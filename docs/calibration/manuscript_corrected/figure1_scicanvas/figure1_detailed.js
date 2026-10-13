() => {
  // Detailed workflow (supplementary): design, preprocessing, SHMtrace calling, Sanger route, reads.
  // 1,400 units = 174 mm (1 unit = 0.35 pt); labels 17-19 units (6-6.7 pt).
  // Numbers: amplicon coordinates and windows (Methods 2.7-2.8); delivered and passing reads from the
  // deposit's pipeline_summary_study_libraries.csv and read counts; leave-one-control-out factor 5
  // (project.json "calibration"); paralogue filter contaminant_min_alleles = 10.
  const W = 1400, o = [];
  const INK = '#1d2733', MID = '#4f5b68', SOFT = '#8a97a6', LINE = '#c9d2dc';
  const JH4 = '#2f6690', N114 = '#2f7a5a', AID = '#c0533a', CTRL = '#8a97a6', GFPC = '#2f8f47';
  const TGT = '#3a8f85', CON = '#d39b3a', ERR = '#c0533a', VAR = '#7b4fa3', TEAL = '#1f6f8b', SANG = '#6d5ba6';
  const F = 'sans', LBL = 19, SM = 17, HEAD = 21;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const A = (x1, y1, x2, y2, e = {}) => { const c = Make.connector({ x: x1, y: y1 }, { x: x2, y: y2 }, { color: MID, width: 2, headSize: 0.9, ...e }); o.push(c); return c; };
  const L = (x1, y1, x2, y2, e = {}) => A(x1, y1, x2, y2, { head: 'none', ...e });
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const panel = (letter, title, y) => { T(letter, 30, y, { fontSize: 34, bold: true }); T(title, 70, y + 6, { fontSize: 25, bold: true }); };
  const read = (x, y, w, c, h = 9) => R(x, y - h / 2, w, h, { fill: c, radius: h / 2 });
  const tick = (x, y, ok) => TC(ok ? '✓' : '✕', x, y - 12, { fontSize: 20, bold: true, color: ok ? TGT : ERR });
  const bracket = (x1, x2, y, c, label, sub) => {
    L(x1, y, x2, y, { color: c, width: 2, head: 'bar', tail: 'bar', headSize: 0.5 });
    const t = T(label, x1, y + 6, { fontSize: SM, color: c, bold: true });
    if (sub) T(sub, t.x + t.w + 6, y + 6, { fontSize: SM, color: MID });
  };
  const BP = 0.27;                                          // both amplicons on one scale
  // a row of numbered step cards; returns the x of each card
  const cards = (steps, y, h, colour) => {
    const n = steps.length, gap = 18, cw = (W - 60 - gap * (n - 1)) / n;
    return steps.map(([num, title, rules], i) => {
      const x = 30 + i * (cw + gap);
      R(x, y, cw, h, { fill: '#f5f7f9', radius: 14 });
      o.push(Make.badge(num, x + 28, y + 30, { fill: colour, labelSize: 16, w: 30, h: 30, x: x + 13, y: y + 15 }));
      T(title, x + 52, y + 15, { fontSize: HEAD - 1, bold: true });
      rules.forEach((r, k) => T(r, x + 14, y + h - 22 - (rules.length - 1 - k) * 22, { fontSize: SM, color: k === 0 ? colour : MID, bold: k === 0 }));
      if (i < n - 1) A(x + cw + 3, y + h / 2, x + cw + gap - 3, y + h / 2, { color: INK });
      return { x, w: cw };
    });
  };

  // ───────── A. Experimental design ─────────
  panel('A', 'Experimental design', 16);
  const l1 = 66, l1h = 250, r1 = l1 + 132;
  R(30, l1, W - 60, l1h, { fill: '#f0f5f9', radius: 16 });
  T('{i|J_{H}4} intron · germinal-centre B cells ({i|in vivo})', 50, l1 + 12, { fontSize: HEAD, bold: true, color: JH4 });
  I('s-mouse-c57bl-6-black-realistic', 90, r1 - 40, 84, { flipX: true });
  T('{i|Wt}', 140, r1 - 62, { bold: true }); T('n = 3', 140, r1 - 37, { fontSize: SM, color: MID });
  I('s-mouse-c57bl-6-black-realistic', 90, r1 + 40, 84, { flipX: true });
  T('{i|Aicda}^{−/−}', 140, r1 + 16, { bold: true }); T('n = 3 · control', 140, r1 + 41, { fontSize: SM, color: MID });
  T('4-month-old, C57BL/6', 50, r1 + 82, { fontSize: SM, color: MID });
  A(272, r1, 296, r1);
  I('s-peyer-s-patch', 350, r1 - 8, 96); TC("Peyer's patches", 350, r1 + 40, { fontSize: SM });
  A(404, r1, 428, r1);
  I('s-facs-sorter-nozzle-laser-plates-tubes', 482, r1 - 10, 88);
  TC('sort GC B cells', 482, r1 + 40, { fontSize: SM }); TC('B220^{+}GL7^{+}Fas^{+}', 482, r1 + 62, { fontSize: SM, color: MID });
  A(536, r1, 560, r1);
  I('s-pcr-tube-strip', 606, r1 - 8, 70); TC('gDNA, PCR', 606, r1 + 40, { fontSize: SM }); TC('SuperFi II', 606, r1 + 62, { fontSize: SM, color: MID });
  A(652, r1, 676, r1);
  const ax = 690, ay = r1 - 30;
  T('{i|J_{H}4} intron amplicon · 1,230 bp', ax, ay - 46, { bold: true });
  R(ax, ay - 8, 1230 * BP, 16, { fill: JH4, radius: 4 });
  bracket(ax + 89 * BP, ax + 1230 * BP, ay + 26, TGT, '1,141 bp', 'ONT, positions 90–1,230');
  bracket(ax + 100 * BP, ax + 629 * BP, ay + 72, SANG, '529 bp', 'shared with Sanger, 101–629');
  A(1030, r1 - 26, 1060, r1 - 56); A(1030, r1 + 22, 1060, r1 + 52);
  I('s-nanopore-sequencer-benchtop-promethion-style', 1104, r1 - 64, 66);
  T('ONT', 1146, r1 - 88, { bold: true }); T('R10.4.1 · Kit 14 · SUP', 1146, r1 - 62, { fontSize: SM, color: MID });
  I('s-sanger-chromatogram-trace-only', 1104, r1 + 62, 78);
  T('Sanger', 1146, r1 + 40, { bold: true, color: SANG }); T('cloned amplicons,', 1146, r1 + 65, { fontSize: SM, color: MID });
  T('separate PCR', 1146, r1 + 86, { fontSize: SM, color: MID });

  const l2 = l1 + l1h + 14, l2h = 236, r2 = l2 + 124;
  R(30, l2, W - 60, l2h, { fill: '#f0f6f2', radius: 16 });
  T('N114 hybridoma · fixed V-region substrate ({i|in vitro})', 50, l2 + 12, { fontSize: HEAD, bold: true, color: N114 });
  I('s-retroviral-plasmid-insert-gfp-puror', 74, r2 - 30, 50);
  T('pMX-{i|Aicda}', 106, r2 - 42, { bold: true, color: AID });
  I('s-retroviral-plasmid-empty-gfp-puror', 74, r2 + 30, 50);
  T('Empty vector', 106, r2 + 18, { bold: true, color: CTRL });
  T('both: GFP, Puro^{R}', 50, r2 + 68, { fontSize: SM, color: GFPC });
  A(256, r2, 280, r2);
  I('s-retroviral-transduction-gfp', 344, r2 - 10, 112); TC('transduce N114', 344, r2 + 44, { fontSize: SM });
  A(408, r2, 432, r2);
  I('s-gfp-histogram-with-high-gate', 482, r2 - 10, 76); TC('sort GFP-high', 482, r2 + 44, { fontSize: SM });
  A(530, r2, 554, r2);
  I('s-96-well-plate-top-view-flat', 606, r2 - 30, 68); I('s-puromycin-selection-flat', 606, r2 + 16, 54);
  TC('limiting dilution,', 606, r2 + 44, { fontSize: SM }); TC('puromycin, 1 month', 606, r2 + 64, { fontSize: SM, color: MID });
  A(652, r2, 676, r2);
  const nx = 690, ny = r2 - 26;
  T('V-region amplicon · 655 bp', nx, ny - 46, { bold: true });
  R(nx, ny - 8, 655 * BP, 16, { fill: N114, radius: 4 });
  R(nx, ny - 8, 20 * BP, 16, { fill: '#a9c9b8', radius: 2 }); R(nx + 635 * BP, ny - 8, 20 * BP, 16, { fill: '#a9c9b8', radius: 2 });
  bracket(nx + 20 * BP, nx + 635 * BP, ny + 26, N114, '615 bp', 'positions 21–635');
  [0, 1, 2].forEach((i) => { E(nx + 6 + i * 22, ny + 82, 16, { fill: AID }); E(nx + 96 + i * 22, ny + 82, 16, { fill: CTRL }); });
  T('3 subclones per vector', nx + 170, ny + 72, { fontSize: SM, color: MID });
  A(1030, r2, 1060, r2);
  I('s-nanopore-sequencer-benchtop-promethion-style', 1104, r2 - 2, 66);
  T('ONT', 1146, r2 - 26, { bold: true }); T('R10.4.1 · Kit 14 · SUP', 1146, r2 + 0, { fontSize: SM, color: MID });
  L(1146, l2 + l2h - 26, 1146 + 200 * BP, l2 + l2h - 26, { color: INK, width: 2 });
  T('200 bp', 1146 + 200 * BP + 10, l2 + l2h - 37, { fontSize: SM, color: MID });

  // ───────── B. Read preprocessing (deposited pipeline) ─────────
  const by = l2 + l2h + 24;
  panel('B', 'Read preprocessing', by);
  const bh = 214, b0 = by + 52;
  const bc = cards([
    [1, 'Target selection', ['best match to target kept', '{i|J_{H}4}: C57BL/6 vs Agouti;', 'N114: one reference']],
    [2, 'Length window', ['full-length reads kept', '{i|J_{H}4} 900–2,000 bp; N114 400–900']],
    [3, 'Expected errors', ['≤ 0.01 per base (≈ Q20)']],
    [4, 'Matched depth', ['500 reads per library', 'one random draw, deposited']],
  ], b0, bh, TGT);
  // 1: target / contaminant / chimera
  let x = bc[0].x;
  [['target', [TGT]], ['paralogue', [CON]], ['chimera', [TGT, CON]]].forEach(([lab, cols], k) => {
    const y = b0 + 70 + k * 26;
    if (cols.length === 1) read(x + 16, y, 84, cols[0]); else { read(x + 16, y, 44, cols[0]); read(x + 58, y, 42, cols[1]); }
    T(lab, x + 110, y - 11, { fontSize: SM });
    tick(x + 230, y, k === 0);
  });
  // 2: size windows on a 0-2 kb axis
  x = bc[1].x; const x0 = x + 74, k2 = 170 / 2000;
  [['{i|J_{H}4}', 900, 2000, JH4, b0 + 76], ['N114', 400, 900, N114, b0 + 112]].forEach(([lab, a, b, c, y]) => {
    T(lab, x + 16, y - 12, { fontSize: SM, bold: true, color: c });
    R(x0, y - 5, 170, 10, { fill: '#e1e6eb', radius: 5 });
    R(x0 + a * k2, y - 5, (b - a) * k2, 10, { fill: c, radius: 5 });
  });
  // 3: expected-error dots
  x = bc[2].x;
  [[6, false], [1, true], [0, true]].forEach(([n, ok], k) => {
    const y = b0 + 70 + k * 30;
    read(x + 16, y, 200, TGT);
    for (let j = 0; j < n; j++) E(x + 34 + j * 30, y, 7, { fill: ERR });
    tick(x + 250, y, ok);
  });
  // 4: down-sampling
  x = bc[3].x;
  for (let j = 0; j < 9; j++) read(x + 18, b0 + 62 + j * 9, 76, TGT, 6);
  A(x + 104, b0 + 98, x + 144, b0 + 98, { color: INK });
  for (let j = 0; j < 4; j++) read(x + 156, b0 + 84 + j * 9, 76, TGT, 6);
  TC('911–2,903', x + 56, b0 + 144, { fontSize: SM, bold: true }); TC('500', x + 194, b0 + 144, { fontSize: SM, bold: true });

  // ───────── C. SHMtrace calling ─────────
  const cy = b0 + bh + 26;
  panel('C', 'Calling with SHMtrace', cy);
  o.push(Make.image(LOGO, 380, cy - 2, 40, 40));
  const ch = 258, c0 = cy + 52;
  const cc = cards([
    [5, 'Align, filter', ['minimap2 to the reference;', 'reads with ≥ 10 paralogue', 'alleles removed']],
    [6, 'Candidates', ['every substitution', 'one read is enough;', 'window as in A']],
    [7, 'Control background', ['per substitution, strand, Q', '× 5 (leave one control out);', 'controls vs the other two']],
    [8, 'Test and filter', ['FDR ≤ 0.05', 'Benjamini–Hochberg;', 'fixed sites (≥ 50%) removed']],
    [9, 'Per animal', ['calls ÷ (reads × bp)', 'spectrum, clones, trees;', 'animal = unit of analysis']],
  ], c0, ch, TEAL);
  // 5: reads on the reference, one paralogue-like read crossed out
  x = cc[0].x; let rw = cc[0].w - 40;
  R(x + 16, c0 + 62, rw, 3, { fill: INK });
  for (let j = 0; j < 4; j++) {
    const y = c0 + 82 + j * 18;
    read(x + 16, y, rw, j === 2 ? CON : TGT, 8);
  }
  tick(x + rw + 28, c0 + 118, false);
  // 6: candidates: every mismatch marked
  x = cc[1].x; rw = cc[1].w - 32;
  R(x + 16, c0 + 62, rw, 3, { fill: INK });
  [[0.2, 0.55], [0.4], [0.2, 0.8], [0.62]].forEach((ms, j) => {
    const y = c0 + 82 + j * 18;
    read(x + 16, y, rw, '#cfd8df', 8);
    ms.forEach((m) => E(x + 16 + m * rw, y, 9, { fill: VAR }));
  });
  // 7: control libraries -> error rate table
  x = cc[2].x;
  [0, 1, 2].forEach((k) => { for (let j = 0; j < 4; j++) read(x + 16 + k * 30, c0 + 66 + j * 9, 24, CTRL, 6); });
  A(x + 110, c0 + 80, x + 136, c0 + 80, { color: INK });
  const grid = [[1, 3, 2], [2, 1, 3], [3, 2, 1]];
  grid.forEach((row, r) => row.forEach((v, c) => R(x + 146 + c * 20, c0 + 58 + r * 16, 19, 15, { fill: ['#eef3f7', '#cfe0ec', '#9fc0d8', '#5f93b9'][v] })));
  T('error rate', x + 146, c0 + 108, { fontSize: SM - 2, color: MID });
  // 8: observed vs expected for a candidate
  x = cc[3].x;
  const bx = x + 30, base = c0 + 128;
  L(bx - 6, base, bx + 180, base, { color: INK, width: 1.4 });
  R(bx + 10, base - 64, 34, 64, { fill: VAR }); TC('observed', bx + 22, base + 4, { fontSize: SM - 2, color: MID });
  R(bx + 90, base - 14, 34, 14, { fill: CTRL }); TC('expected', bx + 107, base + 4, { fontSize: SM - 2, color: MID });
  TC('call ✓', bx + 165, base - 50, { fontSize: SM, bold: true, color: TGT });
  // 9: per-animal dot plot
  x = cc[4].x;
  L(x + 30, c0 + 130, x + 200, c0 + 130, { color: INK, width: 1.4 });
  [[x + 80, [c0 + 66, c0 + 74, c0 + 84], JH4], [x + 160, [c0 + 112, c0 + 116, c0 + 118], CTRL]].forEach(([cx, ys, col]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 12, yy, 11, { fill: col }));
    const m = ys.reduce((a, b) => a + b) / 3;
    L(cx - 24, m, cx + 24, m, { color: INK, width: 2 });
  });

  // ───────── D. Sanger route ─────────
  const dy = c0 + ch + 26;
  panel('D', 'Sanger clones (same calls, one clone per sequence)', dy);
  const d0 = dy + 50, dh = 110;
  R(30, d0, W - 60, dh, { fill: '#f6f4fb', radius: 14 });
  const sanger = [
    ['s-sanger-chromatogram-trace-only', '.ab1 chromatograms', '447 traces, 8 orders'],
    [null, 'Sample sheet and QC', 'duplicates, labels, lineages'],
    [null, 'Trim to clean core', 'aligned to {i|J_{H}4}'],
    [null, 'Substitutions at Q ≥ 20', 'fixed sites (≥ 50%) removed'],
    [null, 'Calls per mouse', 'bases each clone read'],
  ];
  const sw = (W - 60) / sanger.length;
  sanger.forEach(([icon, title, sub], i) => {
    const cx = 30 + sw * i + sw / 2;
    if (icon) I(icon, cx - 112, d0 + 44, 54);
    TC(title, cx + (icon ? 30 : 0), d0 + 30, { bold: true, color: i === 0 ? SANG : INK });
    TC(sub, cx + (icon ? 30 : 0), d0 + 58, { fontSize: SM, color: MID });
    if (i < sanger.length - 1) A(30 + sw * (i + 1) - 12, d0 + dh / 2, 30 + sw * (i + 1) + 12, d0 + dh / 2, { color: SANG });
  });

  // ───────── E. Reads per library ─────────
  const ey = d0 + dh + 26;
  panel('E', 'Reads per library', ey);
  const lgx = 420, ly = ey + 14;
  R(lgx, ly, 26, 12, { fill: '#dfe4ea', radius: 2 }); T('delivered', lgx + 34, ly - 5, { fontSize: SM, color: MID });
  R(lgx + 140, ly, 13, 12, { fill: JH4 }); R(lgx + 153, ly, 13, 12, { fill: AID }); T('passed steps 1–3', lgx + 174, ly - 5, { fontSize: SM, color: MID });
  L(lgx + 340, ly + 6, lgx + 366, ly + 6, { color: INK, width: 2, dashStyle: 'dashed' }); T('500 analysed', lgx + 374, ly - 5, { fontSize: SM, color: MID });
  const groups = [
    { title: '{i|J_{H}4} (GC B cells)', x: 30, rows: [['{i|Wt} 2', 5000, 2729, JH4], ['{i|Wt} 3', 5000, 1906, JH4], ['{i|Wt} 4', 4558, 911, JH4], ['{i|Aicda}^{−/−} 1', 5000, 1492, CTRL], ['{i|Aicda}^{−/−} 2', 5000, 1163, CTRL], ['{i|Aicda}^{−/−} 3', 5000, 1239, CTRL]] },
    { title: 'N114 (subclones)', x: 720, rows: [['pMX-{i|Aicda} C1', 4483, 2903, AID], ['pMX-{i|Aicda} D3', 3924, 1062, AID], ['pMX-{i|Aicda} D5', 4235, 1602, AID], ['empty C3', 2336, 1050, CTRL], ['empty C4', 3153, 1558, CTRL], ['empty D1', 3939, 1412, CTRL]] },
  ];
  let bottom = 0;
  groups.forEach((g) => {
    const top = ey + 86, bx0 = g.x + 170, bw = 380, s = bw / 5000, rh = 30;
    T(g.title, g.x, ey + 48, { fontSize: HEAD - 1, bold: true });
    g.rows.forEach(([lab, del, pass, c], i) => {
      const y = top + i * rh;
      const t = T(lab, 0, y - 2, { fontSize: SM }); t.x = bx0 - 10 - t.w;
      R(bx0, y, del * s, 19, { fill: '#dfe4ea', radius: 2 });
      R(bx0, y, pass * s, 19, { fill: c, radius: 2 });
      T(pass.toLocaleString('en-GB'), bx0 + del * s + 8, y - 0.5, { fontSize: SM, color: MID });
    });
    const b = top + g.rows.length * rh - 6;
    L(bx0 + 500 * s, top - 8, bx0 + 500 * s, b, { color: INK, width: 2, dashStyle: 'dashed' });
    L(bx0, b + 4, bx0 + bw, b + 4, { color: INK, width: 1.4 });
    [0, 1000, 2000, 3000, 4000, 5000].forEach((v) => { L(bx0 + v * s, b + 4, bx0 + v * s, b + 10, { color: INK, width: 1.4 }); TC(v ? `${v / 1000}k` : '0', bx0 + v * s, b + 13, { fontSize: SM, color: MID }); });
    TC('reads', bx0 + bw / 2, b + 36, { fontSize: SM, color: MID });
    bottom = Math.max(bottom, b + 64);
  });
  return { id: 'fig1-detailed', name: 'Workflow (detailed)', width: W, height: bottom, background: '#ffffff', objects: o,
    alt: 'Detailed workflow: experimental design, read preprocessing, SHMtrace calling, the Sanger route and reads per library.' };
}
