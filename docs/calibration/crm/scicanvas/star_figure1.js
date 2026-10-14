() => {
  // STAR Protocols, Figure 1: protocol overview. Bench steps (top) and analysis steps (bottom), each a
  // card with an icon, the step range and the timing. Cell Press rules as in crm_figure1.js.
  const W = 1400, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', OUTL = '#b4bec8';
  const TEAL = '#1f6f8b', BENCH = '#2c5d8a', CTRL = '#9a3439';
  const F = 'arial', TITLE = 22.5, HEAD = 21, LBL = 19.5, SM = 17.5;
  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: LBL, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const card = (x, y, w, h, e = {}) => R(x, y, w, h, { fill: '#ffffff', stroke: OUTL, strokeWidth: 1.6, radius: 10, ...e });
  const LINK = (a, pa, b, pb, e = {}) => { const c = Make.connector({ id: a.id, port: pa }, { id: b.id, port: pb }, { color: INK, width: 2.4, headSize: 0.8, ...e }); o.push(c); return c; };
  const I = (id, cx, cy, size, e = {}) => { const ic = Make.icon(id, 0, 0, size, e); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const SEG = (x1, y1, x2, y2, w, fill) => {
    const len = Math.hypot(x2 - x1, y2 - y1), mx = (x1 + x2) / 2, my = (y1 + y2) / 2;
    return R(mx - len / 2, my - w / 2, len, w, { fill, radius: w / 2, rot: (Math.atan2(y2 - y1, x2 - x1) * 180) / Math.PI });
  };
  const badge = (n, cx, cy, fill) => { const b = Make.badge(n, cx, cy, { fill, labelSize: SM, w: 32, h: 32, x: cx - 16, y: cy - 16 }); o.push(b); return b; };

  const row = (title, y, colour, steps, first) => {
    T(title, title.startsWith('Analysis') ? 340 : 30, y, { fontSize: TITLE, bold: true, color: colour });
    const n = steps.length, gap = 36, cw = (W - 60 - gap * (n - 1)) / n, top = y + 40, h = 300;
    const cards = steps.map(([num, name, icon, lines, timing, extra], i) => {
      const x = 30 + i * (cw + gap);
      const c = card(x, top, cw, h, extra && extra.ctrl ? { stroke: CTRL, strokeWidth: 2 } : {});
      badge(first + i, x + 24, top + 24, colour);
      T(name, x + 46, top + 12, { bold: true });
      const cx = x + cw / 2, cy = top + 104;
      if (icon === 'LOGO') o.push(Make.image(LOGO, cx - 42, top + 62, 84, 84));
      else if (icon === 'DOTS') {          // a group comparison: case high, control low
        [[cx - 40, [cy - 30, cy - 20, cy - 8], '#6d9ad0', '#274b7e'], [cx + 40, [cy + 22, cy + 27, cy + 30], '#d07a7e', CTRL]].forEach(([px, ys, f, e]) => {
          ys.forEach((yy, k) => o.push(Make.ellipse(px + (k - 1) * 13 - 6, yy - 6, 12, 12, { fill: f, stroke: e, strokeWidth: 1.5 })));
          const m = ys.reduce((a, b) => a + b) / 3; SEG(px - 26, m, px + 26, m, 2.4, e);
        });
        SEG(cx - 76, cy + 44, cx + 76, cy + 44, 1.6, INK);
      } else if (icon === 'CURVE') {       // power rising with animals per group
        const pts = [[0, 50], [30, 38], [60, 26], [95, 14], [130, 6]];
        pts.slice(1).forEach(([px, py], k) => SEG(cx - 65 + pts[k][0], cy - 10 + pts[k][1], cx - 65 + px, cy - 10 + py, 2.4, TEAL));
        pts.forEach(([px, py]) => o.push(Make.ellipse(cx - 65 + px - 5, cy - 10 + py - 5, 10, 10, { fill: TEAL, stroke: 'none', strokeWidth: 0 })));
        SEG(cx - 72, cy + 48, cx + 72, cy + 48, 1.6, INK); SEG(cx - 72, cy - 34, cx - 72, cy + 48, 1.6, INK);
      } else {
        const ic = I(icon, cx, cy, 92, { colorMap: { '#43a047': TEAL, '#a5d6a7': '#9cc4d4', '#c8e6c9': '#d6e8ef' } });
        if (ic.color === '#43a047') ic.color = TEAL;
      }
      lines.forEach((s, k) => TC(s, x + cw / 2, top + 168 + k * 23, { fontSize: SM, color: MID }));
      TC((num.includes('–') ? 'steps ' : 'step ') + num + ' · ' + timing, x + cw / 2, top + h - 34, { fontSize: SM, bold: true, color: SOFT });
      return c;
    });
    cards.slice(1).forEach((c, i) => LINK(cards[i], 'e', c, 'w', { headSize: 0.7 }));
    return cards;
  };

  const bench = row('Bench', 10, BENCH, [
    ['1–4', 'GC B cells', 's-facs-sorter-nozzle-laser-plates-tubes', ['Peyer\'s patches; sort', 'B220^{+}GL7^{+}Fas^{+},', '~10,000 cells per mouse'], '4–5 h'],
    ['5–6', 'Genomic DNA', 's-microcentrifuge-tube', ['proteinase K lysis;', 'Wt and AID-deficient', 'mice side by side'], '2 h'],
    ['7–10', 'J_{H}4 PCR', 's-thermal-cycler-lid-open', ['1.2-kb amplicon,', 'high-fidelity polymerase;', 'gel check, clean-up'], '3 h'],
    ['11–12', 'Sequencing', 's-nanopore-sequencing-strand-through-pore', ['full-length amplicons,', 'R10.4.1 flow cells;', 'FASTQ per library'], '1–2 days'],
  ], 1);
  const analysis = row('Analysis with SHMtrace', 412, TEAL, [
    ['13–15', 'Project', 's-spreadsheet', ['sample sheet: file,', 'sample, group, control;', 'app or command line'], '15 min'],
    ['16–17', 'Run', 'LOGO', ['read preparation, control-', 'calibrated calling,', 'animal-level tests'], '2–10 min'],
    ['18–23', 'Report', 'DOTS', ['controls, margin, read', 'accounting, warnings;', 'spectra, maps, clones'], '30 min'],
    ['24', 'Plan', 'CURVE', ['power from this', 'analysis as a pilot:', 'animals per group'], '5 min'],
  ], 5);
  LINK(bench[3], 's', analysis[0], 'n', { style: 'elbow', color: MID, headSize: 0.7 });
  return { id: 'star-fig1', name: 'Figure 1', width: W, height: 772, background: '#ffffff', objects: o,
    alt: 'Figure 1. Protocol overview: bench steps and analysis steps with timing.' };
}
