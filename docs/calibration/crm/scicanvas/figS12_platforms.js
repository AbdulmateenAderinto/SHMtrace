() => {
  // Cell Reports Methods, Figure S12: what Sanger clones and nanopore reads can each measure. Top, one unit
  // of data from each platform on one base-pair scale; bottom, what only one of them gives and what both give,
  // with the agreement measured in the same six mice (Figure 5) and the long-read results (Figure 6).
  // Every number is one the article reports. Built with manuscript_corrected/figure1_scicanvas/build.js into
  // crm/figures/FigureS12_platforms. Mutations purple, sequencing errors orange-red, Sanger amber, nanopore teal.
  const W = 1400, H = 742, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', FAINT = '#cfd8df', OUTL = '#b4bec8';
  const TEAL = '#1f6f8b', SANG = '#a3621b', MUT = '#6a4c93', ERR = '#d0632a';
  const F = 'arial', TITLE = 22.5, HEAD = 20.5, SM = 17.5;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: SM, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const I = (id, cx, cy, size) => { const ic = Make.icon(id, 0, 0, size); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const lines = (arr, x, y, e = {}, dy = 22.5) => arr.map((s, k) => T(s, x, y + k * dy, { color: MID, ...e }));

  // ───────── top: one unit of data, same base-pair scale ─────────
  const pw = 560, AMP = 1230, WIN = [100, 628];
  const head = (x0, x1, title, sub, colour) => {
    T(title, x0, 4, { fontSize: TITLE, bold: true, color: colour });
    T(sub, x0, 34, { color: SOFT, italic: true });
    R(x0, 64, x1 - x0, 2.6, { fill: colour });
  };
  // Sanger: one clone is one molecule, read once over the window it reads cleanly
  const sx = 20, sX = (v) => sx + 70 + v / AMP * pw;
  head(sx, sx + 650, 'SANGER: CLONED AMPLICONS', 'one clone is one molecule, read once', SANG);
  I('s-sanger-chromatogram-trace-only', sx + 34, 140, 58);
  const sClones = [[0.18, 0.55], [], [0.3], [0.12, 0.47, 0.8], [0.62], [], [0.25, 0.7]];
  sClones.forEach((dots, r) => {
    const y = 88 + r * 15;
    R(sX(WIN[0]), y, sX(WIN[1]) - sX(WIN[0]), 7, { fill: FAINT, radius: 3.5 });
    dots.forEach((f) => E(sX(WIN[0] + f * (WIN[1] - WIN[0])), y + 3.5, 10, { fill: MUT }));
  });
  E(sX(WIN[0] + 0.9 * (WIN[1] - WIN[0])), 88 + 4 * 15 + 3.5, 10, { fill: ERR });
  lines(['each clone\'s own', 'mutations'], sX(WIN[1]) + 18, 104, { color: SANG });
  lines(['529-bp window that the clones read cleanly;',
    '447 chromatograms from the six mice of Figure 5'], sx, 212);

  // nanopore: every read spans the amplicon; errors outnumber mutations in any one read
  const nx = 730, nX = (v) => nx + 70 + v / AMP * pw;
  head(nx, nx + 650, 'NANOPORE: CLONING-FREE AMPLICONS', 'every read spans the 1,230-bp amplicon', TEAL);
  I('s-nanopore-sequencing-strand-through-pore', nx + 34, 140, 58);
  // the 529-bp window both platforms are compared over, as brackets above and below the reads
  [80, 211].forEach((yy) => R(nX(WIN[0]), yy, nX(WIN[1]) - nX(WIN[0]), 2.4, { fill: SOFT }));
  const errs = [[0.07, 0.33, 0.61, 0.9], [0.21, 0.52, 0.78], [0.12, 0.44, 0.69, 0.95], [0.29, 0.58, 0.84],
    [0.05, 0.38, 0.73], [0.16, 0.49, 0.88], [0.26, 0.64, 0.92], [0.1, 0.41, 0.8]];
  const muts = { 0: [210, 380], 1: [210, 380], 2: [210], 4: [210, 380], 5: [520], 6: [210] };
  errs.forEach((ev, r) => {
    const y = 88 + r * 15;
    R(nX(0), y, pw, 7, { fill: FAINT, radius: 3.5 });
    ev.forEach((f) => E(nX(f * AMP), y + 3.5, 10, { fill: ERR }));
    (muts[r] || []).forEach((p) => E(nX(p), y + 3.5, 10, { fill: MUT }));
  });
  TC('529-bp window', (nX(WIN[0]) + nX(WIN[1])) / 2, 216, { color: SOFT });
  lines(['500 reads per library. A mutation is called where more reads carry it',
    'than the controls explain; those in one read that are not called are estimated together.'], nx, 240);
  E(sx + 8, 271, 10, { fill: MUT }); T('mutation', sx + 20, 261, { color: MUT });
  E(sx + 110, 271, 10, { fill: ERR }); T('sequencing error', sx + 122, 261, { color: ERR });

  // ───────── bottom: what each gives ─────────
  const cy = 306, ch = H - cy - 8, cw = 440, gap = 20, c1 = 20, c2 = c1 + cw + gap, c3 = c2 + cw + gap;
  const card = (x, title, colour, fill, stroke) => {
    R(x, cy, cw, ch, { fill, stroke, strokeWidth: 1.6, radius: 10 });
    T(title, x + 16, cy + 12, { fontSize: HEAD, bold: true, color: colour });
  };
  // an item: a bold line, then its lines of explanation
  const item = (x, y, title, body, colour = INK) => {
    T(title, x + 16, y, { bold: true, color: colour });
    lines(body, x + 16, y + 23);
    return y + 23 + body.length * 22.5 + 14;
  };

  card(c1, 'From the Sanger clones here', SANG, '#fbf4ec', '#e2c8a6');
  let y = cy + 52;
  y = item(c1, y, 'Each molecule\'s own mutations', ['the mutation load of every clone, the share of',
    'clones mutated, and frequency as published:', 'unique mutations over the bases each clone read'], SANG);
  y = item(c1, y, 'Bookkeeping checks', ['copied traces, clones filed under two mice,',
    'mutated clones filed as controls, mislabeled', 'tubes (Figure 5E)'], SANG);
  y = item(c1, y, 'Error reported, not subtracted', ['the AID-null clones\' error rate (about 20 per',
    '10⁵ bases at Q ≥ 20 in both published studies),', 'and how many of each animal\'s calls it explains'], SANG);

  card(c2, 'From both, and they agree', INK, '#ffffff', OUTL);
  y = cy + 52;
  y = item(c2, y, 'Wild type against AID-null', ['Sanger 38-fold; nanopore 141-fold, with controls', 'without calls counted at one event']);
  y = item(c2, y, 'Substitution spectrum', ['r = 0.92 (95% CI 0.74–0.98), slope 0.97']);
  y = item(c2, y, 'Position along the window', ['r = 0.94 in 25-bp bins']);
  y = item(c2, y, 'AID hotspot use', ['G:C at WRC/GYW 58% and 54% (20% of sites);',
    'A:T at WA/TW 70% and 69% (50% of sites)']);
  y = item(c2, y, 'Not comparable: absolute density', ['a clone is read once; nanopore counts unique',
    'mutations among 500 reads. Compare densities', 'within a platform, not between.'], MID);

  card(c3, 'From the nanopore reads here', TEAL, '#eef5f8', '#b3d0db');
  y = cy + 52;
  y = item(c3, y, 'The whole amplicon', ['positions 90–1,230: 229-fold over controls; 69%',
    'of mutations in the first 450 bp, falling 7-fold', 'to the last 250 bp (Figure 6A)'], TEAL);
  y = item(c3, y, 'Clones and mutation trees', ['clones of each mouse, and their mutations',
    'ordered into trees (Figures 6 and S11)'], TEAL);
  y = item(c3, y, 'Calls at an estimated false-discovery rate', ['each candidate tested against the controls,',
    'the threshold set by them as decoys; single-read', 'mutations not called, as one background-', 'subtracted estimate'], TEAL);

  return { id: 'crm-figS12-platforms', name: 'Figure S12', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure S12. What Sanger clones and nanopore reads can each measure, and what both measure alike in the same mice.' };
}
