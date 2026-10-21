() => {
  // Figure S1 of the Cell Reports Methods article (a copy of the JIM manuscript's detailed workflow, kept
  // apart so that submission's figure is untouched): the full workflow as SHMtrace runs it, from the delivered reads.
  // A design; B read preparation and C variant calling (SHMtrace's own steps and settings, from the
  // presets jh4_intron / n114_vregion and the recommended mode); D the Sanger route; E what happened to
  // every delivered read (READS, SHMtrace's per-library accounting, injected by build.js).
  // Cell Press production rules (checked by build.js): 174 mm wide (1 unit = 0.352 pt), at most 200 mm
  // tall, Arial, text 6-8 pt, strokes 0.5-1.5 pt, no red with green. Flow arrows are attached to the
  // objects they join; decorative lines are thin shapes.
  const W = 1400, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8';
  const JH4 = '#2c5d8a', N114 = '#6a4c93', SANG = '#a3621b', TEAL = '#1f6f8b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const READ = '#5b8db8', PARA = '#d9952b', ERR = '#c0392b', VAR = '#6a4c93';
  const QCG = '#9aa5b1', ALG = '#ccd3da';                        // removed: length/quality; alignment/paralogue
  const F = 'arial', TITLE = 22.5, HEAD = 21, LBL = 19.5, SM = 17.5;
  const GFP_TEAL = { '#3fae5a': '#16938f', '#3f9e48': '#16938f', '#bfe8c3': '#bfe3e1' };

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const TR = (text, rx, y, e = {}) => { const t = T(text, rx, y, e); t.x = rx - t.w; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const SEG = (x1, y1, x2, y2, w, fill) => {
    const len = Math.hypot(x2 - x1, y2 - y1), cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    return R(cx - len / 2, cy - w / 2, len, w, { fill, radius: w / 2, rot: (Math.atan2(y2 - y1, x2 - x1) * 180) / Math.PI });
  };
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: MID, width: 2.4, headSize: 0.8, ...e }); o.push(c); return c; };
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, { colorMap: GFP_TEAL, ...e }); if (GFP_TEAL[ic.color]) ic.color = GFP_TEAL[ic.color]; ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const badge = (n, cx, cy, fill) => { const b = Make.badge(n, cx, cy, { fill, labelSize: SM, w: 32, h: 32, x: cx - 16, y: cy - 16 }); o.push(b); return b; };
  const panel = (letter, title, y) => { T(letter, 30, y, { fontSize: TITLE, bold: true }); T(title, 58, y, { fontSize: TITLE, bold: true }); };
  const read = (x, y, w, c, h = 8) => R(x, y - h / 2, w, h, { fill: c, radius: h / 2 });
  const mark = (x, y, ok) => TC(ok ? '✓' : '✕', x, y - 11, { fontSize: 20, bold: true, color: ok ? CASE : ERR });
  const bracket = (x1, x2, y, c, label, sub) => {
    SEG(x1, y, x2, y, 2, c); SEG(x1, y - 6, x1, y + 6, 2, c); SEG(x2, y - 6, x2, y + 6, 2, c);
    const t = T(label, x1, y + 6, { fontSize: SM, color: c, bold: true });
    if (sub) T(sub, t.x + t.w + 6, y + 6, { fontSize: SM, color: MID });
  };
  // a row of numbered step cards, linked left to right; returns the cards
  const cards = (steps, y, h, colour) => {
    const n = steps.length, gap = 18, cw = (W - 60 - gap * (n - 1)) / n;
    const out = steps.map(([num, title, rules], i) => {
      const x = 30 + i * (cw + gap), k = card(x, y, cw, h);
      badge(num, x + 24, y + 24, colour);
      T(title, x + 46, y + 13, { bold: true });
      rules.forEach((r, j) => T(r, x + 14, y + h - 26 - (rules.length - 1 - j) * 22, { fontSize: SM, color: MID }));
      k.cx = x; k.cw = cw;
      return k;
    });
    out.slice(1).forEach((k, i) => LINK(out[i], 'e', k, 'w', { color: INK, headSize: 0.7 }));
    return out;
  };
  const BP = 0.27;                                              // both amplicons on one scale

  // ───────── A. Experimental design ─────────
  panel('A', 'Experimental design', 10);
  const l1 = 46, l1h = 236, r1 = l1 + 118;
  card(30, l1, W - 60, l1h);
  T('{i|J_{H}4} intron · mouse germinal-center B cells ({i|in vivo})', 46, l1 + 12, { fontSize: HEAD, bold: true, color: JH4 });
  I('s-mouse-c57bl-6-black-realistic', 84, r1 - 34, 70, { flipX: true });
  T('{i|Wt}', 128, r1 - 56, { bold: true, color: CASE }); T('n = 3', 128, r1 - 32, { fontSize: SM, color: MID });
  I('s-mouse-c57bl-6-black-realistic', 84, r1 + 38, 70, { flipX: true });
  T('{i|Aicda}^{−/−}', 128, r1 + 14, { bold: true, color: CTRL }); T('n = 3', 128, r1 + 38, { fontSize: SM, color: MID });
  T('AID-null control', 128, r1 + 60, { fontSize: SM, color: CTRL });
  T('4-month-old, C57BL/6', 46, r1 + 86, { fontSize: SM, color: MID });
  const br1 = R(284, r1 - 52, 3, 104, { fill: OUTL, radius: 1.5 });
  const pp = I('s-peyer-s-patch', 360, r1 - 8, 80); TC("Peyer's patch", 360, r1 + 36, { fontSize: SM });
  const so = I('s-facs-sorter-nozzle-laser-plates-tubes', 494, r1 - 10, 72);
  TC('sort GC B cells', 494, r1 + 36, { fontSize: SM }); TC('B220^{+}GL7^{+}Fas^{+}', 494, r1 + 58, { fontSize: SM, color: MID });
  const pc = I('s-thermal-cycler-lid-open', 616, r1 - 12, 64);
  TC('DNA, PCR', 616, r1 + 36, { fontSize: SM }); TC('SuperFi II', 616, r1 + 58, { fontSize: SM, color: MID });
  const ax = 690, ay = r1 - 30;
  T('{i|J_{H}4} intron amplicon · 1,230 bp', ax, ay - 42, { bold: true });
  const amp1 = R(ax, ay - 7, 1230 * BP, 14, { fill: JH4, radius: 3 });
  bracket(ax + 89 * BP, ax + 1230 * BP, ay + 24, JH4, '1,141 bp', 'nanopore, positions 90–1,230');
  bracket(ax + 100 * BP, ax + 629 * BP, ay + 66, SANG, '529 bp', 'shared with Sanger, 101–629');
  const on1 = I('s-nanopore-sequencing-strand-through-pore', 1104, r1 - 58, 60);
  T('Nanopore', 1144, r1 - 80, { bold: true }); T('R10.4.1, Kit 14', 1144, r1 - 56, { fontSize: SM, color: MID });
  const sa1 = I('s-sanger-chromatogram-trace-only', 1104, r1 + 50, 68);
  T('Sanger', 1144, r1 + 28, { bold: true, color: SANG }); T('cloned amplicons,', 1144, r1 + 52, { fontSize: SM, color: MID });
  T('separate PCRs', 1144, r1 + 74, { fontSize: SM, color: MID });
  LINK(br1, 'e', pp, 'w'); LINK(pp, 'e', so, 'w'); LINK(so, 'e', pc, 'w'); LINK(pc, 'e', amp1, 'w');
  LINK(amp1, 'e', on1, 'w', { style: 'elbow' }); LINK(amp1, 'e', sa1, 'w', { style: 'elbow', color: SANG });

  const l2 = l1 + l1h + 10, l2h = 206, r2 = l2 + 108;
  card(30, l2, W - 60, l2h);
  T('N114 hybridoma · fixed V-region substrate ({i|in vitro})', 46, l2 + 12, { fontSize: HEAD, bold: true, color: N114 });
  I('s-retroviral-plasmid-insert-gfp-puror', 66, r2 - 28, 44);
  T('pMX-{i|Aicda}', 96, r2 - 40, { bold: true, color: CASE });
  I('s-retroviral-plasmid-empty-gfp-puror', 66, r2 + 28, 44);
  T('Empty vector', 96, r2 + 16, { bold: true, color: CTRL });
  T('both carry GFP and Puro^{R}', 46, r2 + 62, { fontSize: SM, color: MID });
  const br2 = R(284, r2 - 46, 3, 92, { fill: OUTL, radius: 1.5 });
  const td = I('s-retroviral-transduction-gfp', 360, r2 - 10, 88); TC('transduce', 360, r2 + 38, { fontSize: SM });
  const gf = I('s-gfp-histogram-with-high-gate', 494, r2 - 10, 62); TC('sort GFP-high', 494, r2 + 38, { fontSize: SM });
  const pl = I('s-96-well-plate-top-view-flat', 616, r2 - 22, 56);
  TC('subclone,', 616, r2 + 38, { fontSize: SM }); TC('puromycin', 616, r2 + 60, { fontSize: SM, color: MID });
  const nx = 690, ny = r2 - 26;
  T('V-region amplicon · 655 bp', nx, ny - 42, { bold: true });
  const amp2 = R(nx, ny - 7, 655 * BP, 14, { fill: N114, radius: 3 });
  R(nx, ny - 7, 20 * BP, 14, { fill: '#c9bcdc', radius: 2 }); R(nx + 635 * BP, ny - 7, 20 * BP, 14, { fill: '#c9bcdc', radius: 2 });
  bracket(nx + 20 * BP, nx + 635 * BP, ny + 24, N114, '615 bp', 'positions 21–635');
  [0, 1, 2].forEach((i) => { E(nx + 8 + i * 20, ny + 72, 15, { fill: CASEF, stroke: CASE, strokeWidth: 1.5 }); E(nx + 92 + i * 20, ny + 72, 15, { fill: CTRLF, stroke: CTRL, strokeWidth: 1.5 }); });
  T('3 subclones per vector', nx + 160, ny + 62, { fontSize: SM, color: MID });
  const on2 = I('s-nanopore-sequencing-strand-through-pore', 1104, ny, 60);
  T('Nanopore', 1144, ny - 22, { bold: true }); T('R10.4.1, Kit 14', 1144, ny + 2, { fontSize: SM, color: MID });
  SEG(1144, l2 + l2h - 22, 1144 + 200 * BP, l2 + l2h - 22, 2, INK);
  T('200 bp', 1144 + 200 * BP + 8, l2 + l2h - 33, { fontSize: SM, color: MID });
  LINK(br2, 'e', td, 'w'); LINK(td, 'e', gf, 'w'); LINK(gf, 'e', pl, 'w'); LINK(pl, 'e', amp2, 'w'); LINK(amp2, 'e', on2, 'w');

  // ───────── B. Read preparation (SHMtrace) ─────────
  const by = l2 + l2h + 16;
  panel('B', 'Read preparation in SHMtrace (from the delivered reads)', by);
  const bh = 210, b0 = by + 36;
  const bc = cards([
    [1, 'Length and quality gates', ['read length ≤ 2,000 bp ({i|J_{H}4}) or', '≤ 900 bp (N114 V_{H}), no minimum;', 'mean expected error ≤ 0.01 per base']],
    [2, 'Competitive alignment', ['minimap2 against target and', 'paralog references; chimeras', 'removed; aligned ≥ 800 / 400 bp']],
    [3, 'Paralog-allele filter', ['reads carrying ≥ 10 bases specific to', 'the {i|J_{H}4} paralog are removed', '(no paralog reference for N114 V_{H})']],
    [4, 'Depth matching', ['random draw of 500 usable reads per', 'library (fixed seed), from 962–2,905', 'usable reads per library']],
  ], b0, bh, JH4);
  let x = bc[0].cx;
  [[160, 0, true], [250, 0, false], [160, 6, false]].forEach(([len, errs, ok], k) => {
    const y = b0 + 56 + k * 28;
    read(x + 14, y, len, READ);
    for (let j = 0; j < errs; j++) E(x + 26 + j * 24, y, 7, { fill: ERR });
    mark(x + 290, y, ok);
  });
  x = bc[1].cx;
  [['target', [READ]], ['paralog', [PARA]], ['chimera', [READ, PARA]]].forEach(([lab, cols], k) => {
    const y = b0 + 56 + k * 28;
    if (cols.length === 1) read(x + 14, y, 84, cols[0]); else { read(x + 14, y, 44, cols[0]); read(x + 56, y, 42, cols[1]); }
    T(lab, x + 108, y - 10, { fontSize: SM });
    mark(x + 290, y, k === 0);
  });
  x = bc[2].cx;
  [[2, true], [12, false]].forEach(([n, ok], k) => {
    const y = b0 + 62 + k * 34;
    read(x + 14, y, 250, READ);
    for (let j = 0; j < n; j++) R(x + 22 + j * 19, y - 7, 4, 14, { fill: PARA });
    mark(x + 290, y, ok);
  });
  x = bc[3].cx;
  for (let j = 0; j < 8; j++) read(x + 18, b0 + 50 + j * 8, 80, READ, 5);
  for (let j = 0; j < 3; j++) read(x + 172, b0 + 70 + j * 8, 80, READ, 5);
  SEG(x + 110, b0 + 78, x + 152, b0 + 78, 2.4, INK);
  R(x + 150, b0 + 72, 2, 12, { fill: INK });
  TC('usable', x + 58, b0 + 112, { fontSize: SM, bold: true }); TC('500', x + 212, b0 + 112, { fontSize: SM, bold: true });

  // ───────── C. Variant calling (SHMtrace) ─────────
  const cy = b0 + bh + 16;
  panel('C', 'Variant calling in SHMtrace', cy);
  o.push(Make.image(LOGO, 380, cy - 3, 28, 28));
  const ch = 250, c0 = cy + 36;
  const cc = cards([
    [5, 'Candidate substitutions', ['every base that differs from the', 'reference in the analyzed window;', 'one supporting read suffices']],
    [6, 'Control error model', ['error rate per substitution, strand and', 'base quality from the AID-null controls;', 'margin (5× J_{H}4, 10× N114 V_{H}) set by', 'leave-one-control-out']],
    [7, 'Significance and filtering', ['each candidate tested against the', 'control error; FDR ≤ 0.05 set by', 'control decoys; sites in ≥ 50% of', 'reads set aside as germline (heuristic)']],
    [8, 'Animal-level summary', ['unique variants per base sequenced;', 'spectra, clones and read-level trees;', 'the mouse or subclone is the unit']],
  ], c0, ch, TEAL);
  x = cc[0].cx; let rw = cc[0].cw - 40;
  R(x + 14, c0 + 52, rw, 3, { fill: INK });
  [[0.2, 0.55], [0.4], [0.2, 0.8], [0.62]].forEach((ms, j) => {
    const y = c0 + 70 + j * 16;
    read(x + 14, y, rw, '#cfd8df', 8);
    ms.forEach((m) => E(x + 14 + m * rw, y, 9, { fill: VAR }));
  });
  x = cc[1].cx;
  [0, 1, 2].forEach((k) => { for (let j = 0; j < 4; j++) read(x + 14 + k * 30, c0 + 58 + j * 9, 24, CTRLF, 6); });
  SEG(x + 108, c0 + 72, x + 132, c0 + 72, 2.4, INK);
  const grid = [[1, 3, 2], [2, 1, 3], [3, 2, 1]];
  grid.forEach((row, r) => row.forEach((v, c) => R(x + 144 + c * 20, c0 + 50 + r * 16, 19, 15, { fill: ['#eef3f7', '#cfe0ec', '#9fc0d8', '#5f93b9'][v] })));
  T('error rate', x + 212, c0 + 64, { fontSize: SM, color: MID });
  x = cc[2].cx;
  const bx = x + 24, base = c0 + 112;
  SEG(bx - 6, base, bx + 210, base, 1.8, INK);
  R(bx + 10, base - 58, 34, 58, { fill: VAR }); TC('observed', bx + 27, base + 4, { fontSize: SM, color: MID });
  R(bx + 100, base - 13, 34, 13, { fill: CTRLF }); TC('expected', bx + 117, base + 4, { fontSize: SM, color: MID });
  TC('call ✓', bx + 180, base - 50, { fontSize: SM, bold: true, color: CASE });
  x = cc[3].cx;
  SEG(x + 26, c0 + 116, x + 216, c0 + 116, 1.8, INK);
  [[x + 76, [c0 + 54, c0 + 62, c0 + 72], CASEF, CASE], [x + 166, [c0 + 100, c0 + 104, c0 + 106], CTRLF, CTRL]].forEach(([cxx, ys, fill, edge]) => {
    ys.forEach((yy, k) => E(cxx + (k - 1) * 12, yy, 11, { fill, stroke: edge, strokeWidth: 1.5 }));
    const m = ys.reduce((a, b) => a + b) / 3;
    SEG(cxx - 24, m, cxx + 24, m, 2.4, edge);
  });

  // ───────── D. Sanger route ─────────
  const dy = c0 + ch + 16;
  panel('D', 'Sanger clones in SHMtrace (one molecule per clone: counted, not tested)', dy);
  const d0 = dy + 36, dh = 88;
  card(30, d0, W - 60, dh, { stroke: SANG, dashStyle: 'dashed' });
  const sanger = [
    ['Chromatograms', '184 traces (.ab1), these six mice'],
    ['Clone assignment', 'sample sheet with QC flags'],
    ['Quality trimming', 'clean core, aligned to {i|J_{H}4}'],
    ['Base calling', 'every substitution at Q ≥ 20'],
    ['Per-mouse summary', 'mutations / (clones × 529 bp)'],
  ];
  const sw = (W - 60) / sanger.length, anchors = [];
  sanger.forEach(([title, sub], i) => {
    const cxx = 30 + sw * i + sw / 2;
    anchors.push(TC(title, cxx, d0 + 20, { bold: true, color: i === 0 ? SANG : INK }));
    TC(sub, cxx, d0 + 48, { fontSize: SM, color: MID });
  });
  // arrows sit in the gaps between steps, clear of the words
  const pin = (x0) => R(x0, d0 + 31, 1, 2, { fill: '#ffffff' });
  anchors.slice(1).forEach((t, i) => LINK(pin(anchors[i].x + anchors[i].w + 12), 'e', pin(t.x - 14), 'w', { color: SANG, headSize: 0.7 }));

  // ───────── E. What SHMtrace did with every delivered read ─────────
  const ey = d0 + dh + 16;
  panel('E', 'Read accounting in SHMtrace', ey);
  const lgx = 460, ly = ey + 6;
  R(lgx, ly, 14, 12, { fill: CASE }); R(lgx + 14, ly, 14, 12, { fill: CTRL }); T('usable', lgx + 36, ly - 5, { fontSize: SM, color: MID });
  R(lgx + 120, ly, 26, 12, { fill: QCG }); T('length / quality', lgx + 154, ly - 5, { fontSize: SM, color: MID });
  R(lgx + 300, ly, 26, 12, { fill: ALG }); T('alignment / paralog', lgx + 334, ly - 5, { fontSize: SM, color: MID });
  SEG(lgx + 540, ly + 6, lgx + 566, ly + 6, 2, INK); T('500 analyzed', lgx + 574, ly - 5, { fontSize: SM, color: MID });
  const split = (a) => [a.usable, a.too_long + a.low_quality, a.not_target + a.chimeric + a.short_alignment + a.paralogue_alleles + a.minor_haplotype];
  const groups = [
    { title: '{i|J_{H}4} (GC B cells)', x: 30, data: READS.jh4, rows: [['{i|Wt} 2', 'WT2', CASE], ['{i|Wt} 3', 'WT3', CASE], ['{i|Wt} 4', 'WT4', CASE], ['{i|Aicda}^{−/−} 1', 'KO1', CTRL], ['{i|Aicda}^{−/−} 2', 'KO2', CTRL], ['{i|Aicda}^{−/−} 3', 'KO3', CTRL]] },
    { title: 'N114 V_{H} (subclones)', x: 720, data: READS.n114, rows: [['pMX-{i|Aicda} C1', 'C1', CASE], ['pMX-{i|Aicda} D3', 'D3', CASE], ['pMX-{i|Aicda} D5', 'D5', CASE], ['empty C3', 'C3', CTRL], ['empty C4', 'C4', CTRL], ['empty D1', 'D1', CTRL]] },
  ];
  let bottom = 0;
  groups.forEach((g) => {
    const top = ey + 70, bx0 = g.x + 170, bw = 380, s = bw / 5000, rh = 27;
    T(g.title, g.x, ey + 36, { fontSize: HEAD, bold: true });
    g.rows.forEach(([lab, key, c], i) => {
      const y = top + i * rh, [use, qc, al] = split(g.data[key]);
      TR(lab, bx0 - 10, y - 2, { fontSize: SM });
      R(bx0, y, use * s, 18, { fill: c });
      R(bx0 + use * s, y, qc * s, 18, { fill: QCG });
      R(bx0 + (use + qc) * s, y, al * s, 18, { fill: ALG });
      T(use.toLocaleString('en-GB'), bx0 + (use + qc + al) * s + 8, y - 0.5, { fontSize: SM, color: MID });
    });
    const b = top + g.rows.length * rh - 6;
    SEG(bx0 + 500 * s, top - 8, bx0 + 500 * s, b, 2, INK);
    SEG(bx0, b + 4, bx0 + bw, b + 4, 1.8, INK);
    [0, 1000, 2000, 3000, 4000, 5000].forEach((v) => { SEG(bx0 + v * s, b + 4, bx0 + v * s, b + 10, 1.8, INK); TC(v ? `${v / 1000}k` : '0', bx0 + v * s, b + 13, { fontSize: SM, color: MID }); });
    TC('delivered reads', bx0 + bw / 2, b + 36, { fontSize: SM, color: MID });
    bottom = Math.max(bottom, b + 64);
  });
  return { id: 'fig-s11', name: 'Figure S1', width: W, height: bottom, background: '#ffffff', objects: o,
    alt: 'Figure S1. Workflow as run by SHMtrace from the delivered reads: experimental design, read preparation, variant calling, the Sanger route and the fate of every delivered read.' };
}
