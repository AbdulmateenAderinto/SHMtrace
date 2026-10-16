() => {
  // Cell Reports Methods, Figure S10: a published Sanger study from another laboratory, re-analysed
  // automatically (Thientosapol et al., PNAS 2018; SAMHD1 and Ig hypermutation, Jolly laboratory).
  // (A) mutations per sequence per arm, (B) the 12-type spectrum, (C) the four classes the paper tests,
  // (D) its claims and what SHMtrace finds. Numbers written into DATA by swhel_data.py from
  // external/out/swhel_samhd1.json and published/thientosapol2018_swhel.json.
  // Panels A-C are SciCanvas 2.0 code plots (matplotlib; re-runnable and editable in the app); D is text.
  // Built with crm/scicanvas/build.js into crm/figures/FigureS10_swhel.
  const HELPERS = /*HELPERS*/"# ruff: noqa: F821  (plt, PAGE and FONT_SIZE come from the SciCanvas plot worker and build.js's preamble)\n# Shared by the CRM code plots (SciCanvas 2.0). Prepended to each plot's code by the builders, so every plot in a\n# .scifig re-runs on its own in SciCanvas. PAGE and the rcParams scaling come from build.js's preamble: sizes below\n# are the printed sizes, in points.\nimport json  # noqa: F401  (each plot reads its numbers with json.loads)\n\nINK, MID, SOFT, FAINT = \"#1d2733\", \"#4a5562\", \"#8a96a3\", \"#cfd8df\"\nSH, PUB = \"#1f6f8b\", \"#a3621b\"           # SHMtrace, published\nCASE, CASEF = \"#274b7e\", \"#6d9ad0\"        # wild type\nKO, KOF = \"#7a4f9a\", \"#bda3d4\"            # knockout\nCTRL, CTRLF = \"#9a3439\", \"#d07a7e\"        # AID-null\nplt.rcParams.update({\"axes.edgecolor\": INK, \"xtick.color\": INK, \"ytick.color\": INK, \"axes.labelcolor\": MID,\n                     \"text.color\": INK})\n\n\ndef area(d):\n    \"\"\"Scatter size for a marker d printed points across.\"\"\"\n    return (d * PAGE) ** 2\n\n\ndef lw(w):\n    \"\"\"A line w printed points wide.\"\"\"\n    return w * PAGE\n\n\ndef plain_log(axis):\n    \"\"\"Log-axis tick labels as plain numbers (10, 100), not 10^n, whose exponents print below 6 pt.\"\"\"\n    axis.set_major_formatter(plt.FuncFormatter(lambda v, _: f\"{v:g}\"))\n    axis.set_minor_formatter(plt.NullFormatter())\n\n\ndef strip(ax, x, vals, face, edge, d=5.0, spread=0.10, marker=\"o\", zorder=3):\n    \"\"\"One group's points, spread sideways so they do not hide each other.\"\"\"\n    vals = list(vals)\n    if not vals:\n        return\n    xs = [x] if len(vals) == 1 else [x - spread + 2 * spread * k / (len(vals) - 1) for k in range(len(vals))]\n    ax.scatter(xs, vals, s=area(d), facecolor=face, edgecolor=edge, linewidth=lw(0.6), marker=marker, zorder=zorder)\n\n\ndef mean_bar(ax, x, vals, color, half=0.16):\n    vals = list(vals)\n    if vals:\n        m = sum(vals) / len(vals)\n        ax.plot([x - half, x + half], [m, m], color=color, linewidth=lw(1.3), solid_capstyle=\"butt\", zorder=4)\n"/*END*/;
  const DATA = /*DATA*/{"A":[{"group":"Samhd1_WT","label":"wt/het","pub":1.79,"sh":1.782,"pub_seqs":565,"sh_seqs":565,"pub_muts":1012,"sh_muts":1007},{"group":"Samhd1_KO","label":"Samhd1","pub":1.4,"sh":1.401,"pub_seqs":354,"sh_seqs":354,"pub_muts":496,"sh_muts":496},{"group":"Samhd1_KO_UGI","label":"Samhd1+UGI","pub":2.48,"sh":2.472,"pub_seqs":218,"sh_seqs":218,"pub_muts":541,"sh_muts":539},{"group":"Msh2_KO","label":"Msh2","pub":1.08,"sh":1.083,"pub_seqs":446,"sh_seqs":445,"pub_muts":480,"sh_muts":482},{"group":"Msh2_Samhd1_KO","label":"Msh2 Samhd1","pub":1.31,"sh":1.315,"pub_seqs":213,"sh_seqs":213,"pub_muts":280,"sh_muts":280}],"B":[[16.798,16.783,"Samhd1_WT","A>T"],[8.597,8.64,"Samhd1_WT","A>C"],[12.352,12.413,"Samhd1_WT","A>G"],[4.348,4.369,"Samhd1_WT","T>A"],[6.324,6.356,"Samhd1_WT","T>C"],[2.767,2.781,"Samhd1_WT","T>G"],[2.668,2.681,"Samhd1_WT","C>A"],[14.032,14.002,"Samhd1_WT","C>T"],[2.569,2.582,"Samhd1_WT","C>G"],[18.775,18.669,"Samhd1_WT","G>A"],[3.458,3.476,"Samhd1_WT","G>T"],[7.312,7.249,"Samhd1_WT","G>C"],[3.024,3.024,"Samhd1_KO","A>T"],[1.411,1.411,"Samhd1_KO","A>C"],[12.097,12.097,"Samhd1_KO","A>G"],[2.218,2.218,"Samhd1_KO","T>A"],[6.452,6.452,"Samhd1_KO","T>C"],[1.613,1.613,"Samhd1_KO","T>G"],[0.403,0.403,"Samhd1_KO","C>A"],[29.032,29.032,"Samhd1_KO","C>T"],[2.016,2.016,"Samhd1_KO","C>G"],[38.508,38.508,"Samhd1_KO","G>A"],[0.806,0.806,"Samhd1_KO","G>T"],[2.419,2.419,"Samhd1_KO","G>C"],[3.697,3.711,"Samhd1_KO_UGI","A>T"],[2.218,2.226,"Samhd1_KO_UGI","A>C"],[11.275,11.317,"Samhd1_KO_UGI","A>G"],[1.109,1.113,"Samhd1_KO_UGI","T>A"],[6.285,6.308,"Samhd1_KO_UGI","T>C"],[1.479,1.484,"Samhd1_KO_UGI","T>G"],[0.185,0.186,"Samhd1_KO_UGI","C>A"],[30.129,30.241,"Samhd1_KO_UGI","C>T"],[0.0,0.0,"Samhd1_KO_UGI","C>G"],[43.438,43.228,"Samhd1_KO_UGI","G>A"],[0.185,0.186,"Samhd1_KO_UGI","G>T"],[0.0,0.0,"Samhd1_KO_UGI","G>C"],[1.042,1.037,"Msh2_KO","A>T"],[0.625,0.622,"Msh2_KO","A>C"],[1.042,1.037,"Msh2_KO","A>G"],[0.417,0.415,"Msh2_KO","T>A"],[0.417,0.415,"Msh2_KO","T>C"],[0.0,0.0,"Msh2_KO","T>G"],[1.458,1.452,"Msh2_KO","C>A"],[43.75,43.776,"Msh2_KO","C>T"],[1.042,1.037,"Msh2_KO","C>G"],[38.125,38.174,"Msh2_KO","G>A"],[3.75,3.734,"Msh2_KO","G>T"],[8.333,8.299,"Msh2_KO","G>C"],[0.357,0.357,"Msh2_Samhd1_KO","A>T"],[0.714,0.714,"Msh2_Samhd1_KO","A>C"],[1.071,1.071,"Msh2_Samhd1_KO","A>G"],[0.357,0.357,"Msh2_Samhd1_KO","T>A"],[0.357,0.357,"Msh2_Samhd1_KO","T>C"],[0.0,0.0,"Msh2_Samhd1_KO","T>G"],[1.786,1.786,"Msh2_Samhd1_KO","C>A"],[43.571,43.571,"Msh2_Samhd1_KO","C>T"],[0.714,0.714,"Msh2_Samhd1_KO","C>G"],[46.071,46.071,"Msh2_Samhd1_KO","G>A"],[1.071,1.071,"Msh2_Samhd1_KO","G>T"],[3.929,3.929,"Msh2_Samhd1_KO","G>C"]],"r":1.0,"C":[{"label":"A:T Tv","wt":[39.08,44.83,33.72,35.21,30.19,38.6,31.08,36.73,30.34,38.98,32.08,32.0],"ko":[6.67,4.26,9.09,9.3,13.04,10.26,15.62,11.54],"pub_wt":32.51,"pub_ko":8.27,"p":9.687037896556261e-11},{"label":"A:T Ts","wt":[25.29,16.09,18.6,21.13,17.92,14.04,24.32,24.49,28.09,20.34,18.87,16.0],"ko":[20.0,27.66,22.73,23.26,23.91,20.51,17.19,36.54],"pub_wt":18.68,"pub_ko":18.55,"p":0.13929817995872312},{"label":"C:G Tv","wt":[12.64,14.94,12.79,15.49,22.64,17.54,20.27,10.2,19.1,13.56,16.98,20.0],"ko":[11.11,8.51,4.55,6.98,6.52,2.56,7.81,7.69],"pub_wt":16.01,"pub_ko":5.65,"p":7.694165705126035e-06},{"label":"C:G Ts","wt":[22.99,24.14,34.88,28.17,29.25,29.82,24.32,28.57,22.47,27.12,32.08,32.0],"ko":[62.22,59.57,63.64,60.47,56.52,66.67,59.38,44.23],"pub_wt":32.81,"pub_ko":67.54,"p":1.2259337589644386e-10}],"D":[["A:T transversions","down ****","down ****"],["A:T transitions","unchanged","unchanged"],["C:G transversions","down ****","down ****"],["C:G transitions","up **","up ****"],["Total mutation","lower, ns","lower ****"]],"foot":"2,804 of the 2,809 mutations the paper counts; r = 1.0000 over the 60 spectrum cells of (B); 1,795 clones, 40 animals, 2 s."}/*END*/;
  const W = 1400, H = 528, o = [];
  const INK = '#1d2733', MID = '#4a5562', SOFT = '#8a96a3', OUTL = '#b4bec8';
  const SH = '#1f6f8b', PUB = '#a3621b';
  const CASE = '#274b7e', CASEF = '#6d9ad0', KO = '#7a4f9a', KOF = '#bda3d4';
  const F = 'arial', TITLE = 22.5, PH = 19.5, SM = 17.5, LW = 2.2, AX = 1.8;
  // one colour per arm, wild type blue through to the double knockout
  const ARM = { Samhd1_WT: [CASE, CASEF], Samhd1_KO: [KO, KOF], Samhd1_KO_UGI: ['#1f6f8b', '#8cc3d6'],
    Msh2_KO: ['#a3621b', '#e0b483'], Msh2_Samhd1_KO: ['#4a5562', '#b9c2cc'] };

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
  const ITAL = (g) => [[g, { italic: true }]];

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

  // a code plot: the shared helpers and this panel's numbers travel with its code
  const CODE = (data, body) => `${HELPERS}\nD = json.loads(${JSON.stringify(JSON.stringify(data))})\n${body}`;

  // ───────── A: mutations per sequence, per arm ─────────
  head('A', 38, 'Mutations per sequence');
  o.push(PLOT({ name: 'A: mutations per sequence', x: 20, y: 40, w: 440, h: 432, code: CODE({ A: DATA.A, arm: ARM }, `
fig, ax = plt.subplots()
rows = D["A"][::-1]                        # wild type at the top
for i, a in enumerate(rows):
    e, f = D["arm"][a["group"]]
    ax.scatter(a["sh"], i + 0.14, s=area(5.5), facecolor=f, edgecolor=e, linewidth=lw(0.7), zorder=3)
    ax.scatter(a["pub"], i - 0.14, s=area(5.0), facecolor="white", edgecolor=e, linewidth=lw(0.8), marker="s", zorder=3)
names = {"wt/het": "wt/het", "Samhd1": "Samhd1", "Samhd1+UGI": "Samhd1 + UGI", "Msh2": "Msh2", "Msh2 Samhd1": "Msh2 Samhd1"}
ax.set_yticks(range(len(rows)), [names[a["label"]] for a in rows])
for lab, a in zip(ax.get_yticklabels(), rows):
    lab.set_color(D["arm"][a["group"]][0])
    if a["label"] != "wt/het":
        lab.set_fontstyle("italic")
ax.set_xlim(0, 3.2); ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_xlabel("Mutations per clone")
ax.tick_params(axis="y", length=0)
ax.scatter([], [], s=area(5.5), facecolor="#b9c2cc", edgecolor=INK, linewidth=lw(0.7), label="SHMtrace")
ax.scatter([], [], s=area(5.0), facecolor="white", edgecolor=INK, linewidth=lw(0.8), marker="s", label="Published")
ax.legend(loc="lower right", frameon=False, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.3)
`) }));

  // ───────── B: the 12-type spectrum, every arm ─────────
  head('B', 494, 'Substitution spectrum');
  o.push(PLOT({ name: 'B: substitution spectrum', x: 476, y: 40, w: 352, h: 432, code: CODE({ B: DATA.B, r: DATA.r, arm: ARM }, `
fig, ax = plt.subplots()
ax.plot([0, 49], [0, 49], color=SOFT, linewidth=lw(0.6), linestyle=(0, (3, 2)), zorder=1)
for x, y, g, t in D["B"]:
    e, f = D["arm"][g]
    ax.scatter(x, y, s=area(4.6), facecolor=f, edgecolor=e, linewidth=lw(0.6), zorder=3)
ax.set_xlim(0, 49); ax.set_ylim(0, 49)
ax.set_xticks([0, 20, 40]); ax.set_yticks([0, 20, 40])
ax.set_xlabel("Published (% of mutations)"); ax.set_ylabel("SHMtrace (% of mutations)")
ax.text(0.04, 0.96, f"r = {D['r']:.4f}\\n12 types × 5 arms", transform=ax.transAxes, va="top", color=MID)
`) }));

  // ───────── C: the four classes the paper tests ─────────
  head('C', 862, 'The four classes tested');
  o.push(PLOT({ name: 'C: the four classes tested', x: 844, y: 40, w: 350, h: 432, code: CODE({ C: DATA.C, arm: ARM }, `
fig, ax = plt.subplots()
(we, wf), (ke, kf) = D["arm"]["Samhd1_WT"], D["arm"]["Samhd1_KO"]
for i, c in enumerate(D["C"]):
    strip(ax, i - 0.2, c["wt"], wf, we, d=4.4, spread=0.09)
    strip(ax, i + 0.2, c["ko"], kf, ke, d=4.4, spread=0.09)
    ax.scatter([i - 0.2], [c["pub_wt"]], s=area(5.2), facecolor="white", edgecolor=we, linewidth=lw(1.0), marker="s", zorder=5)
    ax.scatter([i + 0.2], [c["pub_ko"]], s=area(5.2), facecolor="white", edgecolor=ke, linewidth=lw(1.0), marker="s", zorder=5)
    p = c["p"]
    ax.text(i, 70, "****" if p < 1e-4 else ("*" if p < 0.05 else "ns"), ha="center", va="bottom",
            color=INK if p < 0.05 else MID)
ax.set_xticks(range(4), [c["label"].replace(" ", "\\n") for c in D["C"]])
ax.set_xlim(-0.6, 3.6); ax.set_ylim(0, 104); ax.set_yticks([0, 20, 40, 60])
ax.set_ylabel("% of mutations")
ax.scatter([], [], s=area(4.4), facecolor=wf, edgecolor=we, linewidth=lw(0.6), label="wt/het")
ax.scatter([], [], s=area(4.4), facecolor=kf, edgecolor=ke, linewidth=lw(0.6), label="Samhd1 Δ/Δ")
ax.scatter([], [], s=area(5.2), facecolor="white", edgecolor=INK, linewidth=lw(1.0), marker="s", label="published")
leg = ax.legend(loc="upper left", frameon=False, handletextpad=0.1, borderaxespad=0.0, labelspacing=0.2)
leg.get_texts()[1].set_fontstyle("italic")
`) }));

  // ───────── D: the paper's claims, and what SHMtrace finds ─────────
  {
    const x0 = 1204, pubR = 1300, shR = 1400, y0 = 56;
    head('D', x0 - 34, 'Claims');
    TR('Pub.', pubR, y0 + 34, { bold: true });
    TR('SHMtrace', shR, y0 + 34, { bold: true });
    L(x0, y0 + 60, shR, y0 + 60, SOFT, 1.6);
    DATA.D.forEach(([lab, pub, sh], k) => {
      const y = y0 + 74 + k * 50;
      T(lab, x0, y, { color: INK });
      TR(pub, pubR, y + 24, { color: PUB });
      TR(sh, shR, y + 24, { color: SH });
    });
  }

  // the totals, under the panels
  T(DATA.foot, 92, 478, { color: MID });

  return { id: 'crm-figS10', name: 'Figure S10', width: W, height: H, background: '#ffffff', objects: o,
    alt: 'Figure S10. A published Sanger study from another laboratory re-analysed automatically: mutations per sequence, the substitution spectrum, the four mutation classes the paper tests, and its claims against what SHMtrace finds.' };
}
