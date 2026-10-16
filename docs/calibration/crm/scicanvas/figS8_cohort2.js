() => {
  // Cell Reports Methods, Figure S8: a published Sanger study (Barbulescu et al., 2024, J_H4, cohort 2)
  // re-analysed automatically. (A) frequency per mouse, (B) SHMtrace against published per mouse,
  // (C) mutations per clone as in the published Fig. 2g, (D) G:C transitions and transversions,
  // (E) the pooled 12-type spectrum, (F) the published comparisons recovered. The numbers are written
  // into DATA by cohort2_data.py from the same inputs as figs/figS8_sanger_cohort2.py. The knockout is
  // purple, as in Figure S9 (green beside the AID-null red fails a colour-vision check).
  // Built with manuscript_corrected/figure1_scicanvas/build.js into crm/figures/FigureS8_sanger_cohort2.
  const DATA = /*DATA*/{"A":{"WT":{"sh":[132.33,90.32,90.41,83.18],"pub":[128.1,89.54,92.68,90.78]},"Maea_KO":{"sh":[21.23,20.79,24.95,19.35],"pub":[17.33,17.67,26.0,19.54]},"Aicda_KO":{"sh":[0.0,2.1,0.0,4.73],"pub":[]}},"B":{"points":[[128.099,132.325,"WT"],[89.544,90.317,"WT"],[92.683,90.408,"WT"],[90.782,83.176,"WT"],[17.331,21.225,"Maea_KO"],[17.667,20.794,"Maea_KO"],[25.997,24.953,"Maea_KO"],[19.543,19.354,"Maea_KO"]],"rho":0.88,"ratio":1.03},"C":{"bins":["0","1","2-5","6-10","11-20","21+"],"rows":[{"group":"WT","who":"published","n":96,"pct":[8.0,14.0,34.0,10.0,20.0,13.0]},{"group":"WT","who":"SHMtrace","n":92,"pct":[26.09,8.7,22.83,17.39,13.04,11.96]},{"group":"Maea_KO","who":"published","n":183,"pct":[26.0,30.0,33.0,7.0,4.0,0.0]},{"group":"Maea_KO","who":"SHMtrace","n":179,"pct":[38.55,22.91,30.73,6.15,1.68,0.0]}]},"D":[{"group":"WT","kind":"Ts","sh":[19.48,20.16,17.27,18.18],"pub":[19.41,20.0,18.64,15.83]},{"group":"WT","kind":"Tv","sh":[18.83,13.18,12.73,15.91],"pub":[17.06,13.55,11.86,21.67]},{"group":"Maea_KO","kind":"Ts","sh":[12.5,9.09,9.09,11.63],"pub":[14.29,7.84,9.46,13.46]},{"group":"Maea_KO","kind":"Tv","sh":[20.31,24.24,19.7,18.6],"pub":[18.57,21.57,20.27,17.31]}],"E":{"WT":[[8.881,9.356,"A>T"],[10.124,9.356,"A>C"],[19.716,19.958,"A>G"],[8.171,8.524,"T>A"],[13.499,13.514,"T>C"],[4.973,4.99,"T>G"],[4.618,3.534,"C>A"],[7.993,7.9,"C>T"],[3.552,3.119,"C>G"],[10.657,11.019,"G>A"],[2.842,2.911,"G>T"],[4.973,5.821,"G>C"]],"Maea_KO":[[15.385,14.078,"A>T"],[7.287,7.767,"A>C"],[22.267,22.816,"A>G"],[8.097,9.709,"T>A"],[10.931,10.194,"T>C"],[5.263,4.369,"T>G"],[5.263,6.311,"C>A"],[4.858,4.369,"C>T"],[5.668,5.34,"C>G"],[6.478,6.311,"G>A"],[3.644,4.369,"G>T"],[4.858,4.369,"G>C"]]},"r":0.99,"F":[["Wt / Maea frequency","5.0\u00d7","4.6\u00d7",0],["P (one-tailed t)","0.014","0.0002",1],["G:C Ts P","0.0028","0.0003",0],["G:C Tv P","0.0285","0.025",0],["Mice, Wt / Maea / null","4 / 4 / 3","4 / 4 / 4",0],["Clones, Wt and Maea","284","271",0],["Background explains, Wt","\u2013","1-2%",0],["...                 Maea","\u2013","7-8%",1],["Lineages under 2+ mice","\u2013","14 (49 clones)",0],["Analysis","hand-edited","16 s",0]]}/*END*/;
  const W = 1400, H = 940, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#8a96a3', OUTL = '#b4bec8';
  const SH = '#1f6f8b', PUB = '#a3621b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', KO = '#7a4f9a', KOF = '#bda3d4', CTRL = '#9a3439', CTRLF = '#d07a7e';
  const F = 'arial', TITLE = 22.5, PH = 19.5, SM = 17.5, LW = 2.2, AX = 1.8;
  const COL = { WT: [CASE, CASEF], Maea_KO: [KO, KOF], Aicda_KO: [CTRL, CTRLF] };
  const ROW2 = 472;

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
  const MAEA1 = [['Maea', { italic: true }], ['fl/fl', { sup: true }]], MAEA2 = [['Mb1', { italic: true }], ['cre', { sup: true }]];
  const GLAB = { WT: [['Wt']], Maea_KO: [...MAEA1, [' '], ...MAEA2],
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
  const head = (letter, x, title, top = 0) => {
    T(letter, x - 34, top + 2, { fontSize: TITLE, bold: true });
    T(title, x, top + 6, { fontSize: PH, bold: true });
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
    const x0 = 92, x1 = 366, y0 = 56, y1 = 386;
    head('A', x0 - 34, 'Frequency per mouse');
    const [X, Y] = frame(x0, y0, x1, y1, [-0.5, 2.5], [0, 230], [0, 50, 100, 150, 200], ['Mutations per 10⁴ bp']);
    ['WT', 'Maea_KO', 'Aicda_KO'].forEach((g, i) => {
      pair(X, Y, i, DATA.A[g].sh, DATA.A[g].pub, ...COL[g]);
      if (g === 'Maea_KO') { runC(MAEA1, X(i), y1 + 9); runC(MAEA2, X(i), y1 + 33); }
      else runC(GLAB[g], X(i), y1 + 9);
    });
    E(x1 - 104, y0 + 16, 12, { fill: '#b9c2cc', stroke: INK, strokeWidth: 1.6 });
    T('SHMtrace', x1 - 92, y0 + 6, { color: MID });
    SQ(x1 - 104, y0 + 42, 11, { fill: '#ffffff', stroke: INK, strokeWidth: 1.6 });
    T('Published', x1 - 92, y0 + 32, { color: MID });
  }

  // ───────── B: SHMtrace against published, per mouse (log-log) ─────────
  {
    const x0 = 520, x1 = 776, y0 = 56, y1 = 386, lg = Math.log10, dom = [lg(8), lg(400)];
    head('B', x0 - 34, 'Agreement by mouse');
    const [X, Y] = frame(x0, y0, x1, y1, dom, dom, [], ['SHMtrace (per 10⁴ bp)']);
    [10, 100].forEach((v) => {
      L(x0 - 6, Y(lg(v)), x0, Y(lg(v)), INK, AX); TR(String(v), x0 - 9, Y(lg(v)) - 10);
      L(X(lg(v)), y1, X(lg(v)), y1 + 6, INK, AX); TC(String(v), X(lg(v)), y1 + 9);
    });
    TC('Published (per 10⁴ bp)', (x0 + x1) / 2, y1 + 34, { color: MID });
    L(X(dom[0]), Y(dom[0]), X(dom[1]), Y(dom[1]), SOFT, AX, { dash: true });
    DATA.B.points.forEach(([x, y, g]) => E(X(lg(x)), Y(lg(y)), 13, { fill: COL[g][1], stroke: COL[g][0], strokeWidth: 1.6 }));
    T(`ρ = ${DATA.B.rho.toFixed(2)}, 8 mice`, x0 + 12, y0 + 6, { color: MID });
    T(`ratio ${DATA.B.ratio.toFixed(2)}× (geometric)`, x0 + 12, y0 + 30, { color: MID });
  }

  // ───────── C: mutations per clone, binned as in the published Fig. 2g ─────────
  {
    const x0 = 860, x1 = 1380, y0 = 56;
    head('C', x0 - 34, 'Mutations per clone');
    const SHADES = ['#eef1f4', '#d4e1ee', '#a9c3dd', '#6f99c4', '#3f6fa4', '#1d4377'];
    // the bins as a key along the top
    let kx = x0;
    DATA.C.bins.forEach((b, k) => {
      R(kx, y0 + 40, 14, 14, { fill: SHADES[k], stroke: INK, strokeWidth: 1.5 });
      const t = T(b, kx + 18, y0 + 36, { color: MID });
      kx += 18 + t.w + 10;
    });
    const X = (v) => x0 + v / 100 * (x1 - x0);
    let y = y0 + 74;
    DATA.C.rows.forEach((row, i) => {
      if (i === 2) y += 14;                       // a gap between the two genotypes
      const lab = row.group === 'WT' ? [['Wt']] : [...MAEA1, [' '], ...MAEA2];
      const tail = `, ${row.who} (n = ${row.n})`;
      if (row.group === 'WT') T('Wt' + tail, x0, y);
      else run([...lab, [tail]], x0, y, { color: INK });
      let left = 0;
      row.pct.forEach((v, k) => {
        R(X(left), y + 24, X(left + v) - X(left), 28, { fill: SHADES[k], stroke: INK, strokeWidth: 1.5 });
        if (v >= 9) TC(v.toFixed(0), (X(left) + X(left + v)) / 2, y + 27, { color: k >= 4 ? '#ffffff' : INK });
        left += v;
      });
      y += 62;
    });
    L(x0, y + 4, x1, y + 4, INK, AX);
    [0, 20, 40, 60, 80, 100].forEach((v) => { L(X(v), y + 4, X(v), y + 10, INK, AX); TC(String(v), X(v), y + 13); });
    TC('% of clones', (x0 + x1) / 2, y + 37, { color: MID });
  }

  // ───────── D: G:C transitions and transversions ─────────
  {
    const x0 = 92, x1 = 366, y0 = ROW2 + 56, y1 = ROW2 + 386;
    head('D', x0 - 34, 'G:C Ts and Tv', ROW2);
    const [X, Y] = frame(x0, y0, x1, y1, [-0.5, 3.5], [0, 35], [0, 10, 20, 30], ['% of mutations at G:C']);
    DATA.D.forEach((d, i) => {
      pair(X, Y, i, d.sh, d.pub, ...COL[d.group]);
      TC(d.kind, X(i), y1 + 9);
      if (d.group === 'WT') TC('Wt', X(i), y1 + 31); else runC([['Maea', { italic: true }]], X(i), y1 + 31);
    });
  }

  // ───────── E: pooled substitution spectrum ─────────
  {
    const x0 = 520, x1 = 776, y0 = ROW2 + 56, y1 = ROW2 + 386;
    head('E', x0 - 34, 'Substitution spectrum', ROW2);
    const [X, Y] = frame(x0, y0, x1, y1, [0, 30], [0, 30], [0, 10, 20, 30], ['SHMtrace (%)']);
    [0, 10, 20, 30].forEach((v) => { L(X(v), y1, X(v), y1 + 6, INK, AX); TC(String(v), X(v), y1 + 9); });
    TC('Published (%)', (x0 + x1) / 2, y1 + 34, { color: MID });
    L(X(0), Y(0), X(30), Y(30), SOFT, AX, { dash: true });
    ['WT', 'Maea_KO'].forEach((g) => DATA.E[g].forEach(([x, y]) => E(X(x), Y(y), 12, { fill: COL[g][1], stroke: COL[g][0], strokeWidth: 1.6 })));
    TR(`r = ${DATA.r.toFixed(2)}`, x1, y0 + 6, { color: MID });
    ['WT', 'Maea_KO'].forEach((g, k) => {
      E(X(1.5), Y(28 - k * 4), 12, { fill: COL[g][1], stroke: COL[g][0], strokeWidth: 1.6 });
      run(g === 'WT' ? [['Wt']] : [...MAEA1, [' '], ...MAEA2], X(1.5) + 12, Y(28 - k * 4) - 10, { color: MID });
    });
  }

  // ───────── F: the published comparisons recovered ─────────
  // Each value column is right-aligned at its own edge, so a long cell cannot run into its neighbour.
  {
    const x0 = 860, pubR = 1272, shR = 1400, y0 = ROW2 + 56;
    head('F', x0 - 34, 'Published tests', ROW2);
    TR('Published', pubR, y0 + 4, { bold: true });
    TR('SHMtrace', shR, y0 + 4, { bold: true });
    L(x0, y0 + 30, shR, y0 + 30, SOFT, 1.6);
    DATA.F.forEach(([lab, pub, sh, indent], k) => {
      const y = y0 + 44 + k * 32;
      if (lab.startsWith('...')) run([['...        '], ['Maea', { italic: true }]], x0 + indent * 14, y, { color: INK });
      else T(lab, x0 + indent * 14, y, { color: INK });
      TR(pub, pubR, y, { color: PUB });
      TR(sh, shR, y, { color: SH });
    });
  }

  return { id: 'crm-figS8', name: 'Figure S8', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure S8. A published Sanger study re-analysed automatically: frequency per mouse, agreement by mouse, mutations per clone, G:C transitions and transversions, the pooled substitution spectrum, and the published comparisons recovered.' };
}
