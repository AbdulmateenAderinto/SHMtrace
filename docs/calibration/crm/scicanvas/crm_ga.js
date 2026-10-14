() => {
  // Graphical abstract, Cell Reports Methods: square, one panel, read top to bottom; larger type than
  // the figures (Arial 9-11 pt at 120 mm). 1 unit = 0.12 mm when exported at 120 mm square.
  const W = 1000, H = 1000, o = [];
  const INK = '#1d2733', MID = '#4a5562', OUTL = '#b4bec8', TEAL = '#1f6f8b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', CTRL = '#9a3439', CTRLF = '#d07a7e', SANG = '#a3621b';
  const F = 'arial', BIG = 34, MED = 29, SM = 26;
  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: MED, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 2.4, radius: 14, ...e });
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const SEG = (x1, y1, x2, y2, w, fill) => {
    const len = Math.hypot(x2 - x1, y2 - y1), cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    return R(cx - len / 2, cy - w / 2, len, w, { fill, radius: w / 2, rot: (Math.atan2(y2 - y1, x2 - x1) * 180) / Math.PI });
  };
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: INK, width: 3.4, headSize: 1.0, ...e }); o.push(c); return c; };
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };

  // 1. inputs and matched controls
  const top = card(30, 30, 940, 250);
  [['s-nanopore-sequencing-strand-through-pore', 'Nanopore', INK], ['s-sanger-chromatogram-trace-only', 'Sanger', SANG],
   ['s-production-sequencer-novaseq-style', 'Illumina', '#3d5a80']].forEach(([id, lab, c], k) => {
    const cx = 130 + k * 190;
    I(id, cx, 115, k === 2 ? 132 : 96);   // the sequencer is wider than tall: larger box, same visual weight
    TC(lab, cx, 182, { bold: true, color: c, fontSize: SM });
  });
  R(640, 60, 2.4, 190, { fill: OUTL });
  I('s-mouse-c57bl-6-black-realistic', 720, 110, 84, { flipX: true });
  I('s-mouse-c57bl-6-black-realistic', 880, 110, 84, { flipX: true });
  TC('samples', 720, 160, { fontSize: SM, color: CASE, bold: true });
  TC('AID-null', 880, 160, { fontSize: SM, color: CTRL, bold: true });
  TC('controls', 880, 192, { fontSize: SM, color: CTRL, bold: true });
  TC('same sequencing run', 800, 232, { fontSize: SM - 3, color: MID });

  // 2. SHMtrace
  const mid = card(30, 330, 940, 250, { stroke: TEAL, strokeWidth: 3.4 });
  o.push(Make.image(LOGO, 70, 375, 120, 120));
  T('SHMtrace', 220, 372, { fontSize: BIG, bold: true, color: TEAL });
  T('each candidate mutation tested against', 220, 428, { fontSize: MED, color: MID });
  T('the errors of matched AID-null libraries;', 220, 466, { fontSize: MED, color: MID });
  T('groups compared animal by animal', 220, 504, { fontSize: MED, color: MID });

  // 3. outcomes
  const bw = (940 - 2 * 24) / 3, by = 630, bh = 340, outs = [];
  ['Low background', 'Agrees with Sanger', 'Clones and trees'].forEach((t, k) => {
    const x = 30 + k * (bw + 24);
    outs.push(card(x, by, bw, bh));
    TC(t, x + bw / 2, by + 18, { fontSize: SM, bold: true });
  });
  // a. separation
  let x = 30;
  [[x + 90, [by + 110, by + 130, by + 150], CASEF, CASE], [x + 210, [by + 250, by + 256, by + 260], CTRLF, CTRL]].forEach(([cx, ys, f, e]) => {
    ys.forEach((yy, k) => E(cx + (k - 1) * 18, yy, 18, { fill: f, stroke: e, strokeWidth: 2 }));
    const m = ys.reduce((a, b) => a + b) / 3; SEG(cx - 40, m, cx + 40, m, 3.4, e);
  });
  SEG(x + 30, by + 290, x + bw - 30, by + 290, 2.4, INK);
  TC('1–4 calls per control set', x + bw / 2, by + 300, { fontSize: SM - 4, color: MID });
  // b. spectrum agreement
  x = 30 + bw + 24;
  const pts = [[0.1, 0.12], [0.2, 0.17], [0.3, 0.33], [0.42, 0.38], [0.5, 0.55], [0.63, 0.6], [0.75, 0.8], [0.85, 0.9]];
  SEG(x + 40, by + 280, x + bw - 30, by + 280, 2.4, INK); SEG(x + 40, by + 80, x + 40, by + 280, 2.4, INK);
  SEG(x + 40, by + 280, x + bw - 40, by + 90, 2, OUTL);
  pts.forEach(([a, b]) => E(x + 40 + a * (bw - 80), by + 280 - b * 190, 16, { fill: TEAL }));
  TC('nanopore vs Sanger, r = 0.92', x + bw / 2, by + 300, { fontSize: SM - 4, color: MID });
  // c. tree
  x = 30 + 2 * (bw + 24);
  const rx = x + 40, ry = by + 190;
  R(rx - 9, ry - 9, 18, 18, { fill: INK }); SEG(rx, ry, rx + 70, ry, 3, INK);
  [[ry - 80, 140], [ry - 30, 175], [ry + 25, 110], [ry + 75, 155]].forEach(([yy, len], k) => {
    SEG(rx + 70, ry, rx + 70, yy, 2.4, OUTL); SEG(rx + 70, yy, rx + 70 + len, yy, 3, INK);
    for (let j = 1; j <= 3; j++) R(rx + 70 + (len * j) / 4 - 2.5, yy - 9, 5, 18, { fill: ['#c0392b', '#e08a2c', '#2c5d8a'][(j + k) % 3] });
  });
  TC('from full-length reads', x + bw / 2, by + 300, { fontSize: SM - 4, color: MID });
  LINK(top, 's', mid, 'n');
  outs.forEach((c) => LINK(mid, 's', c, 'n', { style: 'elbow' }));
  return { id: 'crm-ga', name: 'Graphical abstract', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Graphical abstract: sample reads and matched AID-null controls from nanopore, Sanger or Illumina sequencing are analyzed by SHMtrace, which suppresses background, agrees with Sanger and adds clones and trees.' };
}
