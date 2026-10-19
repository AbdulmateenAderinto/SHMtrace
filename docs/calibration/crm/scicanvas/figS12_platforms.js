() => {
  // Cell Reports Methods, Figure S12: what Sanger clones and nanopore reads can each measure, and how Illumina
  // reads enter the same analysis. Top, one unit of data from each of the first two platforms on one base-pair
  // scale; middle, how an Illumina pair becomes one read and how much of the locus each platform's reads cover;
  // bottom, what only one of Sanger and nanopore gives and what both give, with the agreement measured in the
  // same six mice (Figure 5) and the long-read results (Figure 6).
  // Every number is one the article reports. Built with manuscript_corrected/figure1_scicanvas/build.js into
  // crm/figures/FigureS12_platforms. Mutations purple, sequencing errors orange-red, Sanger amber, nanopore teal.
  const W = 1400, BAND = 392, PB = 232, H = 742 + BAND + PB, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#6f7b88', FAINT = '#cfd8df', OUTL = '#b4bec8';
  const TEAL = '#1f6f8b', PACB = '#9a3439', PACBF = '#e2c2c2', ILL = '#3d5a80', ILLF = '#c9d6e3', SANG = '#a3621b', MUT = '#6a4c93', ERR = '#d0632a', ERRF = '#f0bf9e';
  const F = 'arial', TITLE = 22.5, HEAD = 20.5, SM = 17.5;

  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { fontSize: SM, color: INK, family: F, ...e }); o.push(t); return t; };
  const TC = (text, cx, y, e = {}) => { const t = T(text, cx, y, { align: 'center', ...e }); t.x = cx - t.w / 2; return t; };
  const R = (x, y, w, h, e = {}) => { const r = Make.rect(x, y, w, h, { stroke: 'none', strokeWidth: 0, radius: 0, ...e }); o.push(r); return r; };
  const E = (cx, cy, d, e = {}) => { const r = Make.ellipse(cx - d / 2, cy - d / 2, d, d, { stroke: 'none', strokeWidth: 0, ...e }); o.push(r); return r; };
  const I = (id, cx, cy, size) => { const ic = Make.icon(id, 0, 0, size); ic.x = cx - ic.w / 2; ic.y = cy - ic.h / 2; o.push(ic); return ic; };
  const PAD = (fs) => { const w = (s) => Make.text(s, 0, 0, { fontSize: fs, family: F }).w; return 2 * w('x') - w('xx'); };
  // text with J_H4 written with a subscript H (figS14_engine.js does the same)
  const TJ = (text, x, y, e = {}) => {
    let cx = x;
    const put = (str, dy) => { const t = T(str, cx, y + dy, e); cx += t.w - PAD(t.fontSize); };
    text.split('J_H').forEach((part, k) => { if (k > 0) { put('J', 0); put('H', 6); } if (part) put(part, 0); });
    return cx - x;
  };
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

  // nanopore: eight real reads of wild-type library WT3, re-run over the whole reliably covered amplicon
  // (positions 90-1,230, as Figure 6A; scicanvas/s12_reads.py), with every mutation SHMtrace called in them.
  // A quarter of this library's calls lie beyond the 529-bp Sanger window, so the reads carry mutations there
  // too. Errors are drawn, not read: a read's uncalled mismatches cannot be told from uncalled mutations. They
  // are drawn at the rate the AID-null controls show over the same window once germline sites are removed,
  // 0.135% per base with 15% at Q >= 20, at the expected counts (13 over these reads, 2 at Q >= 20), at random
  // positions clear of the mutations; open, Q < 20; filled, Q >= 20.
  const nx = 730, nX = (v) => nx + 70 + v / AMP * pw;
  head(nx, nx + 650, 'NANOPORE: CLONING-FREE AMPLICONS', 'every read spans the 1,230-bp amplicon', TEAL);
  I('s-nanopore-sequencing-strand-through-pore', nx + 34, 140, 58);
  // the 529-bp window both platforms are compared over, as brackets above and below the reads
  [80, 211].forEach((yy) => R(nX(WIN[0]), yy, nX(WIN[1]) - nX(WIN[0]), 2.4, { fill: SOFT }));
  const READS = [[52,1230,[107,136,166,222,316,320,368,508],[[1006,0]]],[52,1230,[107,138,144,166,275,290,320,323,368,403,422,425,442,446,474],[]],[51,1230,[138,139,160,252,261,479],[[199,0],[309,0],[996,0]]],[51,1230,[138,139,160,252,261,479],[[420,0],[864,0]]],[56,1227,[139,160,165,209,307,363,383,395,422,455,488,490,550,559,565,566,695,739],[[468,0],[623,1],[930,0]]],[56,1230,[139,160,165,166,209,307,363,383,395,422,455,488,490,517,550,559,565,566,695,739,1187],[[185,1]]],[63,1230,[107,146,147,151,157,160,166,167,182,217,273,275,279,285,300,320,330,368,428,431,490,557,696,713],[[1227,0]]],[1,1230,[445,700,867,997],[[97,0],[1095,0]]]];
  READS.forEach(([a, b, mut, err], r) => {
    const y = 88 + r * 15;
    R(nX(a), y, nX(b) - nX(a), 7, { fill: FAINT, radius: 3.5 });
    // errors are drawn smaller than mutations, and the ones below Q 20 are paler: at this locus mutations
    // outnumber errors several-fold, and equal-weight markers would not show that
    mut.forEach((p) => E(nX(p), y + 3.5, 7, { fill: MUT }));
    err.forEach(([p, hq]) => E(nX(p), y + 3.5, 6, hq ? { fill: ERR } : { fill: ERRF }));
  });
  TC('529-bp window', (nX(WIN[0]) + nX(WIN[1])) / 2, 216, { color: SOFT });
  lines(['Eight wild-type reads over the whole amplicon, with every mutation called in them;',
    'errors drawn at the rate the controls show, 1 base in 700, pale below Q 20'], nx, 240);
  E(sx + 8, 271, 10, { fill: MUT }); T('mutation', sx + 20, 261, { color: MUT });
  E(sx + 110, 271, 10, { fill: ERR }); T('sequencing error', sx + 122, 261, { color: ERR });

  // ───────── middle: Illumina ─────────
  // Pair merging as io_utils.merge_pair does it, drawn to scale from the median pair of two public libraries
  // (scicanvas/s12_illumina.py: the first 20,000 pairs of a wild-type J_H4 mouse, PRJNA973112, and of an
  // AID-expressing CH12F3 replicate, PRJNA918596). Inside the overlap each base is read twice: agreeing calls add
  // their qualities (capped at Q60), a disagreement keeps the better call, or becomes N on a tie. Read-through
  // past the fragment (adapter) is dropped, and a pair that does not overlap by 20 bp is dropped whole.
  const iy = 300;
  T('ILLUMINA: PAIRED SHORT READS', 20, iy, { fontSize: TITLE, bold: true, color: ILL });
  T('each pair is merged into one read, which then takes the same read route as nanopore', 20, iy + 30, { color: SOFT, italic: true });
  R(20, iy + 60, 1360, 2.6, { fill: ILL });
  I('s-flow-cell-clusters-illumina', 52, iy + 112, 66);
  const mx = 170, k = 1.2, X = (bp) => mx + bp * k, PALE = '#eef1f4';
  // one pair: R1 from 0, R2 ending at the fragment's end; reads longer than the fragment run into adapter
  const pair = (y, title, frag, rlen, note, nPos) => {
    TJ(title, mx - 50, y, { bold: true });
    const r1 = y + 34, r2 = y + 52, mg = y + 80, a = Math.max(frag - rlen, 0), b = Math.min(rlen, frag);
    R(X(0), r1, X(b) - X(0), 7, { fill: FAINT, radius: 3.5 });
    R(X(a), r2, X(frag) - X(a), 7, { fill: FAINT, radius: 3.5 });
    if (rlen > frag) {
      R(X(frag), r1, (rlen - frag) * k, 7, { fill: PALE, radius: 3.5 });
      R(X(frag - rlen), r2, (rlen - frag) * k, 7, { fill: PALE, radius: 3.5 });
      T('adapter, dropped', X(rlen) + 10, r1 - 6, { color: SOFT });
    }
    T('R1 →', mx - 50, r1 - 6, { color: MID });
    T('← R2', X(Math.max(frag, rlen)) + 10, r2 - 6, { color: MID });
    if (a > 0) R(X(0), mg, X(a) - X(0), 10, { fill: ILLF, radius: 2 });
    R(X(a), mg, X(b) - X(a), 10, { fill: ILL, radius: 2 });
    if (b < frag) R(X(b), mg, X(frag) - X(b), 10, { fill: ILLF, radius: 2 });
    R(X(nPos) - 1.5, mg - 4, 3, 18, { fill: ERR });
    T('merged', mx - 84, mg - 5, { color: ILL, bold: true });
    T(note, X(0), mg + 20, { color: MID });
  };
  pair(iy + 76, 'J_H4 (Ldha study): 143-bp fragment, 2 × 151-bp reads', 143, 151, 'every base read twice; 95% of pairs merge', 70);
  pair(iy + 210, 'CH12F3 VDJ: 269-bp fragment, 2 × 150-bp reads', 269, 150, '31 bp read twice, the rest once; 97% merge', 135);
  // key for the merged bars
  const ky = iy + 352;
  R(20, ky + 4, 26, 10, { fill: ILL, radius: 2 }); T('read twice (Q60)', 54, ky - 2, { color: MID });
  R(204, ky + 4, 26, 10, { fill: ILLF, radius: 2 }); T('read once', 238, ky - 2, { color: MID });
  R(336, ky - 1, 3, 20, { fill: ERR }); T('disagreement: better call kept, N on a tie', 348, ky - 2, { color: MID });
  // how much of the locus one read covers, to scale
  const lx = 760, bx0 = 948, kb = 0.24, ly = iy + 76;
  T('How much of the locus one read covers', lx, ly, { bold: true });
  // A junction-bearing PacBio read runs to several kb, so its bar is clipped and the true median
  // labelled; drawing it to scale would shrink every other platform to a sliver.
  const CLIP = 1300;
  [['PacBio, S\u03bc\u2013S\u03b31', 6126, PACB], ['Nanopore, J_H4', 1141, TEAL],
    ['Sanger, J_H4', 529, SANG], ['Illumina, Ramos VDJ', 421, ILL],
    ['Illumina, CH12F3 VDJ', 229, ILL], ['Illumina, J_H4 (Ldha)', 70, ILL]].forEach(([name, bp, c], r) => {
    const y = ly + 34 + r * 26;
    const shown = Math.min(bp, CLIP);
    TJ(name, lx, y - 4, { color: c });
    R(bx0, y, shown * kb, 14, { fill: c, radius: 2 });
    if (bp > CLIP) {
      const xe = bx0 + shown * kb;
      [0, 5, 10].forEach(d => R(xe + d - 12, y, 2.4, 14, { fill: '#ffffff' }));
      T('\u2026', xe + 4, y - 4, { color: c });
    }
    T(`${bp.toLocaleString('en-US')} bp`, bx0 + shown * kb + (bp > CLIP ? 24 : 8), y - 4,
      { color: MID });
  });
  lines(['Illumina settings: minimap2 sr preset; error gate 0.02 per base;',
    'minimum alignment ≤ 100 bp; then the same control test and FDR.',
    'The Ldha and CH12F3 files carry a flat Q30 (SRA), so their error',
    'model has a single quality level; the Ramos files keep binned ones.'], lx, ly + 192);

  // ───────── PacBio: one read carries a whole switch junction ─────────
  // The only platform here whose read spans two loci at once: a HiFi read runs from the Sμ donor,
  // through the join, into the Sγ1 acceptor. Lengths from crm/out/switch_pacbio/switch_pacbio.json
  // (median junction-bearing read over the 12 samples); microhomology classes are that study's own
  // unique junctions. The donor/acceptor split varies with where the break falls, so the bar is drawn
  // at the median length but its split is schematic.
  const py = 300 + BAND + 8;
  T('PACBIO: ONE READ ACROSS A SWITCH JUNCTION', 20, py, { fontSize: TITLE, bold: true, color: PACB });
  T('the donor, the join and the acceptor on a single HiFi read', 20, py + 30, { color: SOFT, italic: true });
  R(20, py + 60, 1360, 2.6, { fill: PACB });
  const jx = 170, jw = 600, jsplit = 0.44, jy = py + 108;
  R(jx, jy, jw * jsplit, 16, { fill: PACB, radius: 2 });
  R(jx + jw * jsplit, jy, jw * (1 - jsplit), 16, { fill: PACBF, radius: 2 });
  R(jx + jw * jsplit - 1.6, jy - 7, 3.2, 30, { fill: INK });
  TC('join', jx + jw * jsplit, jy - 32, { color: INK, bold: true });
  TC('S\u03bc donor', jx + jw * jsplit / 2, jy + 26, { color: PACB });
  TC('S\u03b31 acceptor', jx + jw * jsplit + jw * (1 - jsplit) / 2, jy + 26, { color: MID });
  T('6,126 bp median', jx + jw + 16, jy - 4, { color: MID });
  // the detail belongs in the legend; two lines here, as in the other bands
  lines(['the join, its microhomology and any inserted bases are read, not inferred;',
    'reads per unique junction run 1.1\u201323.3, so counts track PCR recovery, not switching'], 20, jy + 60);
  // what the junctions look like, from the study's own unique junctions
  const mh = [['blunt', 14.4], ['1\u20133 nt', 38.1], ['\u2265 4 nt', 12.6], ['inserted bases', 34.9]];
  const mhx = 1010, mhbar = 150;
  T('Microhomology at the join', mhx, py + 86, { bold: true });
  mh.forEach(([lab, pct], r) => {
    const y = py + 114 + r * 22;
    TJ(lab, mhx, y - 3, { color: MID });
    R(mhx + 118, y, mhbar * pct / 100, 12, { fill: PACB, radius: 2 });
    T(`${pct}%`, mhx + 118 + mhbar * pct / 100 + 6, y - 3, { color: MID });
  });

  // ───────── bottom: what each gives ─────────
  const cy = 306 + BAND + PB, ch = H - cy - 8, cw = 440, gap = 20, c1 = 20, c2 = c1 + cw + gap, c3 = c2 + cw + gap;
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
    'clones mutated, and frequency as published:', 'unique mutations over clones \u00d7 a fixed length'], SANG);
  y = item(c1, y, 'Bookkeeping checks', ['copied traces, clones filed under two mice,',
    'mutated clones filed as controls, mislabeled', 'tubes (Figure S7)'], SANG);
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
    alt: 'Figure S12. What Sanger clones and nanopore reads can each measure, how Illumina read pairs are merged and analyzed, and what Sanger and nanopore measure alike in the same mice.' };
}
