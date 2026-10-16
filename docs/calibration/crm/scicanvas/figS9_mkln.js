() => {
  // Cell Reports Methods, Figure S9: a second published Sanger study (Barbulescu et al., 2025, Mkln1)
  // re-analysed automatically. (A) frequency per mouse, (B) G:C transitions and transversions,
  // (C) the pooled 12-type spectrum, (D) the published comparisons recovered. The numbers are written
  // into DATA by mkln_data.py from the same inputs as figs/figS9_sanger_mkln.py. Panel D is laid out in
  // measured columns: each value column is right-aligned at its own edge, so a long cell cannot run into
  // the next one. Built with manuscript_corrected/figure1_scicanvas/build.js into crm/figures/FigureS9_sanger_mkln.
  const DATA = /*DATA*/{"A":{"WT":{"sh":[87.64,137.81,81.92],"pub":[72.94,94.79,35.48]},"Mkln1_KO":{"sh":[18.9,8.54,10.88],"pub":[8.86,11.2,6.08]},"Aicda_KO":{"sh":[0.0,2.1,2.1],"pub":[]}},"B":[{"group":"WT","kind":"Ts","sh":[20.26,19.03,26.92],"pub":[24.66,22.56,23.29]},{"group":"WT","kind":"Tv","sh":[19.61,15.04,11.54],"pub":[17.12,16.41,13.7]},{"group":"Mkln1_KO","kind":"Ts","sh":[33.33,16.67,21.43],"pub":[32.35,22.03,29.17]},{"group":"Mkln1_KO","kind":"Tv","sh":[29.63,26.19,30.95],"pub":[26.47,30.51,25.0]}],"C":{"WT":[[10.386,11.387,"A>T"],[7.246,7.66,"A>C"],[16.184,15.735,"A>G"],[6.522,8.696,"T>A"],[13.285,14.493,"T>C"],[6.763,5.176,"T>G"],[2.415,2.277,"C>A"],[10.628,9.524,"C>T"],[2.899,2.899,"C>G"],[12.802,11.594,"G>A"],[5.072,4.348,"G>T"],[5.797,6.211,"G>C"]],"Mkln1_KO":[[9.402,6.306,"A>T"],[5.983,5.405,"A>C"],[11.111,17.117,"A>G"],[8.547,6.306,"T>A"],[6.838,9.91,"T>C"],[3.419,3.604,"T>G"],[5.983,6.306,"C>A"],[13.675,7.207,"C>T"],[6.838,4.505,"C>G"],[12.821,15.315,"G>A"],[7.692,9.009,"G>T"],[7.692,9.009,"G>C"]]},"r":0.82,"D":[["Wt / Mkln1KO","7.8\u00d7","8.0\u00d7",0],["P (one-tailed t)","0.0137","0.0038",1],["G:C Ts P","0.23","0.77",0],["G:C Tv P","0.004","0.008",0],["... S2C bars\u1d43","0.013 / 0.31","\u2013",1],["Clones, Wt / KO","110 / 227","88 / 193",0],["Paralog clones out","\u2013","85",0],["Background explains, Wt","\u2013","2-3%",0],["...                 Mkln1KO","\u2013","12-30%",1],["Analysis","Benchling","15 s",0]],"traces":402,"clones":309}/*END*/;
  const W = 1400, H = 560, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#8a96a3', OUTL = '#b4bec8';
  const SH = '#1f6f8b', PUB = '#a3621b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', KO = '#7a4f9a', KOF = '#bda3d4', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const F = 'arial', TITLE = 22.5, PH = 19.5, SM = 17.5, LW = 2.2, AX = 1.8;
  const COL = { WT: [CASE, CASEF], Mkln1_KO: [KO, KOF], Aicda_KO: [CTRL, CTRLF] };

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: SM, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const TR = (text, rx, y, e = {}) => { const t = T(text, rx, y, { align: 'right', ...e }); t.x = rx - t.w; return t; };
  const TV = (text, cx, cy, e = {}) => { const t = T(text, 0, 0, { align: 'center', rot: -90, ...e }); t.x = cx - t.w / 2; t.y = cy - t.h / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const SQ = (cx, cy, d, e = {}) => R(cx - d / 2, cy - d / 2, d, d, e);
  const PL = (pts, color = INK, width = LW, e = {}) => {
    const p = makePathFromNodes(pts.map(([x, y]) => ({ x, y })), { stroke: color, strokeWidth: width, cap: 'butt', ...e }); o.push(p); return p;
  };
  const L = (x1, y1, x2, y2, color = INK, width = LW, e = {}) => PL([[x1, y1], [x2, y2]], color, width, e);
  const width = (s, e = {}) => Make.text(s, 0, 0, { fontSize: SM, family: F, ...e }).w;
  // a run of text pieces on one line: [text, {italic, sup, bold, color}]
  const PAD = (fs) => 2 * Make.text('x', 0, 0, { fontSize: fs, family: F }).w - Make.text('xx', 0, 0, { fontSize: fs, family: F }).w;
  const kern = (_pieces, _i, e) => PAD(e.fontSize || SM);
  const runW = (pieces, e = {}) => pieces.reduce((a, [s, p = {}], i) => a + width(s, { ...e, italic: !!p.italic, bold: p.bold ?? e.bold }) - kern(pieces, i, e), 0);
  const run = (pieces, x, y, e = {}) => {
    let cx = x;
    pieces.forEach(([s, p = {}], i) => {
      const t = T(s, cx, y + (p.sup ? -6 : 0), { ...e, italic: !!p.italic, bold: p.bold ?? e.bold, color: p.color || e.color || INK });
      cx += t.w - kern(pieces, i, e);
    });
    return cx - x;
  };
  const runC = (pieces, cx, y, e = {}) => run(pieces, cx - runW(pieces, e) / 2, y, e);
  const runR = (pieces, rx, y, e = {}) => run(pieces, rx - runW(pieces, e), y, e);
  // the genotype labels, as runs: Mkln1 and Aicda are gene names, so italic, with a superscript −/−
  const GLAB = { WT: [['Wt']], Mkln1_KO: [['Mkln1', { italic: true }], ['−/−', { sup: true }]],
    Aicda_KO: [['Aicda', { italic: true }], ['−/−', { sup: true }]] };

  // an axis frame; returns the data-to-page maps
  const frame = (x0, y0, x1, y1, xd, yd, yt, ylab) => {
    const X = (v) => x0 + (v - xd[0]) / (xd[1] - xd[0]) * (x1 - x0);
    const Y = (v) => y1 - (v - yd[0]) / (yd[1] - yd[0]) * (y1 - y0);
    L(x0, y0 - 4, x0, y1, INK, AX);
    L(x0, y1, x1, y1, INK, AX);
    yt.forEach((v) => { L(x0 - 6, Y(v), x0, Y(v), INK, AX); TR(String(v), x0 - 9, Y(v) - 10); });
    if (ylab) [].concat(ylab).forEach((s, k, a) => TV(s, x0 - 56 - (a.length - 1 - k) * 22, (y0 + y1) / 2, { color: MID }));
    return [X, Y];
  };
  const head = (letter, x, title) => {
    T(letter, x - 34, 2, { fontSize: TITLE, bold: true });
    T(title, x, 6, { fontSize: PH, bold: true });
  };
  // one group's points: SHMtrace filled circles left, published open squares right, means as bars
  const pair = (X, Y, cx, sh, pub, edge, fill) => {
    const half = (X(1) - X(0)) * 0.17, spread = (X(1) - X(0)) * 0.05;
    const put = (vals, dx, draw, colour) => {
      vals.forEach((v, k) => draw(X(cx) + dx + (vals.length > 1 ? (2 * k / (vals.length - 1) - 1) * spread : 0), Y(v)));
      if (vals.length) R(X(cx) + dx - half * 0.75, Y(vals.reduce((a, b) => a + b, 0) / vals.length) - 1.6, half * 1.5, 3.2, { fill: colour });
    };
    put(sh, -half, (x, y) => E(x, y, 12, { fill, stroke: edge, strokeWidth: 1.6 }), SH);
    put(pub, half, (x, y) => SQ(x, y, 11, { fill: '#ffffff', stroke: edge, strokeWidth: 1.6 }), PUB);
  };

  // ───────── A: frequency per mouse ─────────
  {
    const x0 = 92, x1 = 366, y0 = 56, y1 = 452;
    head('A', x0 - 34, 'Frequency per mouse');
    const [X, Y] = frame(x0, y0, x1, y1, [-0.5, 2.5], [0, 215], [0, 50, 100, 150, 200],
      ['Mutations per 10⁴ bp']);
    ['WT', 'Mkln1_KO', 'Aicda_KO'].forEach((g, i) => {
      pair(X, Y, i, DATA.A[g].sh, DATA.A[g].pub, ...COL[g]);
      runC(GLAB[g], X(i), y1 + 9);
    });
    // key
    E(x0 + 14, y0 + 16, 12, { fill: '#b9c2cc', stroke: INK, strokeWidth: 1.6 });
    T('SHMtrace', x0 + 26, y0 + 6, { color: MID });
    SQ(x0 + 14, y0 + 42, 11, { fill: '#ffffff', stroke: INK, strokeWidth: 1.6 });
    T('Published', x0 + 26, y0 + 32, { color: MID });
  }

  // ───────── B: G:C transitions and transversions ─────────
  {
    const x0 = 466, x1 = 650, y0 = 56, y1 = 452;
    head('B', x0 - 34, 'G:C Ts and Tv');
    const [X, Y] = frame(x0, y0, x1, y1, [-0.5, 3.5], [0, 40], [0, 10, 20, 30, 40],
      ['% of mutations at G:C']);
    DATA.B.forEach((b, i) => {
      pair(X, Y, i, b.sh, b.pub, ...COL[b.group]);
      TC(b.kind, X(i), y1 + 9);
      TC(b.group === 'WT' ? 'Wt' : 'KO', X(i), y1 + 31);
    });
  }

  // ───────── C: pooled substitution spectrum ─────────
  {
    const x0 = 716, x1 = 892, y0 = 56, y1 = 452;
    head('C', x0 - 34, 'Substitution spectrum');
    const [X, Y] = frame(x0, y0, x1, y1, [0, 28], [0, 28], [0, 10, 20], ['SHMtrace (%)']);
    [0, 10, 20].forEach((v) => { L(X(v), y1, X(v), y1 + 6, INK, AX); TC(String(v), X(v), y1 + 9); });
    L(X(0), Y(0), X(28), Y(28), SOFT, AX, { dash: true });
    ['WT', 'Mkln1_KO'].forEach((g) => {
      const [e, f] = COL[g];
      DATA.C[g].forEach(([x, y]) => E(X(x), Y(y), 12, { fill: f, stroke: e, strokeWidth: 1.6 }));
    });
    TC('Published (%)', (x0 + x1) / 2, y1 + 34, { color: MID });
    TR(`r = ${DATA.r.toFixed(2)}`, x1, y0 + 6, { color: MID });
    ['WT', 'Mkln1_KO'].forEach((g, k) => {
      E(X(1.5), Y(26 - k * 4), 12, { fill: COL[g][1], stroke: COL[g][0], strokeWidth: 1.6 });
      run(GLAB[g], X(1.5) + 12, Y(26 - k * 4) - 10, { color: MID });
    });
  }

  // ───────── D: the published comparisons ─────────
  // Each value column is right-aligned at its own edge and the label column is clipped to the gap
  // before it, so no cell can run into its neighbour however long the numbers become.
  {
    const x0 = 990, pubR = 1285, shR = 1400, y0 = 56;
    head('D', x0 - 34, 'Published tests');
    TR('Published', pubR, y0 + 34, { bold: true });
    TR('SHMtrace', shR, y0 + 34, { bold: true });
    L(x0, y0 + 60, shR, y0 + 60, SOFT, 1.6);
    const LAB = { 'Wt / Mkln1KO': [['Wt / '], ['Mkln1', { italic: true }], ['−/−', { sup: true }]],
      '...                 Mkln1KO': [['...        '], ['Mkln1', { italic: true }], ['−/−', { sup: true }]] };
    DATA.D.forEach(([lab, pub, sh, indent], k) => {
      const y = y0 + 76 + k * 33;
      if (LAB[lab]) run(LAB[lab], x0 + indent * 14, y, { color: INK });
      else T(lab, x0 + indent * 14, y, { color: INK });
      TR(pub, pubR, y, { color: PUB });
      TR(sh, shR, y, { color: SH });
    });
    T('ᵃ Ts / Tv; fit the tables only transposed', x0, y0 + 76 + DATA.D.length * 33 + 10, { color: MID });
  }

  return { id: 'crm-figS9', name: 'Figure S9', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure S9. A second published Sanger study re-analysed automatically: frequency per mouse, G:C transitions and transversions, the pooled substitution spectrum, and the published comparisons recovered.' };
}
