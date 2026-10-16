// Build a SciCanvas figure from a builder script, outside the app (the Cell Reports Methods figures).
//
//   node build.js <builder.js> <out-stem>      (SCICANVAS=~/SciCanvas-merge, SciCanvas 2.1, by default)
//
// A copy of manuscript_corrected/figure1_scicanvas/build.js with SciCanvas's code plots added. A builder
// calls PLOT({ code, x, y, w, h, name }) with matplotlib code; this runs it through SciCanvas's own plot worker
// (scripts/plot-worker.py: the same rcParams, palette and helpers the app gives it) at its printed size, after
// the figure's shared plot code (plots_common.py, saved as the document's plot prelude), and places it as a code
// plot: an SVG asset plus the code, so the .scifig reopens in SciCanvas with every plot re-runnable and its labels
// editable. The plot style is in printed points, as in the app.
//
// The builder is a `() => page` function written against SciCanvas's own `Make` API, the same
// form SciCanvas runs from its console, so it can also be pasted into the app. This loads the
// app's source the way SciCanvas's tests do (test/load-app.js), with one change: text is
// measured with Helvetica's glyph widths (SciCanvas draws `sans` as Helvetica) instead of the
// tests' 7-px stub, so centred and right-aligned labels land where the app puts them.
// Writes <out-stem>.scifig (open it in SciCanvas to edit), <out-stem>.svg (SciCanvas's renderer),
// <out-stem>.pdf (vector, fonts embedded; Chrome, headless) and <out-stem>_1000dpi.tif (RGB, LZW),
// and checks the page against Cell Press's figure requirements for production
// (cell.com/information-for-authors/figure-guidelines): 17.4 cm wide, at most 20 cm tall; Arial;
// text 6-8 pt; strokes 0.5-1.5 pt; grey fills between 10% and 80%. Any breach stops the build.
const fs = require('fs'), vm = require('vm'), path = require('path'), os = require('os');
const { execFileSync } = require('child_process');

const SCI = process.env.SCICANVAS || path.join(os.homedir(), 'SciCanvas-merge');
const PYTHON = process.env.PLOT_PYTHON || path.join(__dirname, '..', '..', '..', '..', '.venv', 'bin', 'python');
const CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const PRINT_MM = 174;                                  // full-width figure
const MAX_H_MM = 200;
const TIFF_DPI = 1000;                                 // Cell Press: line art at least 1,000 dpi

// Helvetica (regular) advance widths, per 1000 em, for ASCII 32-126; Helvetica-Bold is ~6% wider
const HELV = [278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556,
  556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667, 611, 778, 722,
  278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469,
  556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278,
  556, 500, 722, 500, 500, 500, 334, 260, 334, 584];
const WIDE = { '–': 556, '—': 1000, '·': 278, '×': 584, '≥': 549, '≤': 549, '−': 584, '′': 250, '→': 1000,
  '✓': 760, '✕': 760, 'µ': 556, 'α': 600, 'β': 570, 'Δ': 612 };
function measure(font, text) {
  const px = parseFloat((/(\d+(?:\.\d+)?)px/.exec(font) || [, 16])[1]);
  const bold = /\bbold\b/.test(font);
  let w = 0;
  for (const ch of String(text)) {
    const c = ch.codePointAt(0);
    w += c >= 32 && c <= 126 ? HELV[c - 32] : (WIDE[ch] || 556);
  }
  return (w / 1000) * px * (bold ? 1.06 : 1);
}

function loadApp() {
  const dir = path.join(SCI, 'src');
  const order = [...fs.readFileSync(path.join(dir, 'index.html'), 'utf8').matchAll(/<script src="([^"]+)"/g)]
    .map((m) => m[1]).filter((f) => !f.startsWith('..'));
  const blob = () => new Proxy(function () {}, {
    get: (t, k) => (k === Symbol.toPrimitive ? () => 0 : k === 'then' ? undefined : k === 'length' ? 0
      : k === Symbol.iterator ? () => [][Symbol.iterator]() : blob()),
    apply: () => blob(), construct: () => blob(), set: () => true, has: () => true,
  });
  const ctx = {
    console, Math, Date, JSON, Promise, Map, Set, WeakMap, Symbol, Array, Object, String, Number, Boolean, RegExp,
    Error, TypeError, parseFloat, parseInt, isFinite, isNaN, encodeURIComponent, decodeURIComponent, escape, unescape,
    btoa: (s) => Buffer.from(s, 'binary').toString('base64'), atob: (s) => Buffer.from(s, 'base64').toString('binary'),
    Uint8Array, Uint16Array, Float32Array, Float64Array, Int32Array, ArrayBuffer, DataView, TextEncoder, TextDecoder,
    Intl, structuredClone, queueMicrotask, setTimeout: () => 0, clearTimeout() {}, setInterval: () => 0,
    clearInterval() {}, requestAnimationFrame: () => 0, cancelAnimationFrame() {},
    localStorage: { getItem: () => null, setItem() {}, removeItem() {} }, sessionStorage: { getItem: () => null, setItem() {} },
    navigator: { platform: 'MacIntel', userAgent: 'node', clipboard: blob() }, location: { search: '', hash: '' },
    matchMedia: () => ({ matches: false, addEventListener() {} }), getComputedStyle: () => ({ getPropertyValue: () => '' }),
    DOMParser: blob(), XMLSerializer: blob(), Image: blob(), FileReader: blob(), Blob: blob(), URL: blob(), fetch: blob(),
    ResizeObserver: blob(), MutationObserver: blob(), IntersectionObserver: blob(), Path2D: blob(), OffscreenCanvas: blob(),
    CustomEvent: blob(), Event: blob(), $3Dmol: blob(), JSZip: blob(), PptxGenJS: blob(), SmilesDrawer: blob(),
    pako: blob(), UTIF: blob(), polygonClipping: require(path.join(SCI, 'node_modules', 'polygon-clipping')),
  };
  const canvas = () => {
    const c2d = { font: '16px sans', fillText() {}, drawImage() {}, getImageData: () => ({ data: [] }) };
    c2d.measureText = (t) => {
      const px = parseFloat((/(\d+(?:\.\d+)?)px/.exec(c2d.font) || [, 16])[1]);
      return { width: measure(c2d.font, t), actualBoundingBoxAscent: px * 0.72, actualBoundingBoxDescent: px * 0.21 };
    };
    return { getContext: () => c2d, width: 0, height: 0, toDataURL: () => '' };
  };
  ctx.document = new Proxy({ createElement: (t) => (t === 'canvas' ? canvas() : blob()), createElementNS: () => blob(),
    addEventListener() {}, getElementById: () => blob(), querySelector: () => blob(), querySelectorAll: () => [],
    body: blob(), documentElement: blob(), fonts: { ready: Promise.resolve(), add() {} } },
  { get: (t, k) => (k in t ? t[k] : blob()) });
  ctx.window = ctx; ctx.globalThis = ctx; ctx.self = ctx; ctx.window.native = blob(); ctx.addEventListener = () => {};
  const code = order.map((f) => `/* ${f} */\n` + fs.readFileSync(path.join(dir, f), 'utf8')).join('\n;\n')
    + '\n;({ Make, makeGroup, pageSvgString, ICON_MAP, state, normalizeColors, check: (p) => checkFigure(p) })';
  return { app: vm.runInNewContext(code, ctx, { filename: 'scicanvas-bundle.js' }), ctx };
}

// Cell Press production requirements, measured at print size. GA=1 checks a graphical abstract:
// the same rules, but type up to 18 pt and strokes up to 2.5 pt, as a one-panel summary reads larger.
const GA = process.env.GA === '1';
function check(page) {
  const pt = (u) => (u * PRINT_MM / page.width) / 25.4 * 72;          // page units -> points
  const bad = [];
  const hMM = page.height * PRINT_MM / page.width;
  if (hMM > MAX_H_MM + 0.05) bad.push(`height ${hMM.toFixed(1)} mm > ${MAX_H_MM} mm`);
  const fonts = new Set(); let fmin = 99, fmax = 0, smin = 99, smax = 0;
  const walk = (list) => list.forEach((o) => {
    if (o.children) walk(o.children);
    if (o.type === 'text' || o.label) {
      const size = pt(o.type === 'text' ? o.fontSize : (o.labelSize || 16));
      fmin = Math.min(fmin, size); fmax = Math.max(fmax, size);
      if (o.type === 'text') fonts.add(o.family);
      if (size < 5.95 || size > (GA ? 18.05 : 8.05)) bad.push(`text ${size.toFixed(2)} pt: ${String(o.text || o.label).slice(0, 30)}`);
    }
    const w = o.type === 'connector' ? o.width : (o.stroke && o.stroke !== 'none' ? o.strokeWidth : 0);
    if (w) {
      smin = Math.min(smin, pt(w)); smax = Math.max(smax, pt(w));
      if (pt(w) < 0.495 || pt(w) > (GA ? 2.505 : 1.505)) bad.push(`stroke ${pt(w).toFixed(2)} pt (${o.type})`);
    }
    // a pure grey fill (r = g = b, not white) must be 10-80% grey
    const hex = /^#([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(o.fill || '');
    if (hex && hex[1] === hex[2] && hex[2] === hex[3] && hex[1].toLowerCase() !== 'ff') {
      const grey = 1 - parseInt(hex[1], 16) / 255;
      if (grey < 0.1 || grey > 0.8) bad.push(`grey fill ${o.fill} is ${(grey * 100).toFixed(0)}%`);
    }
  });
  walk(page.objects);
  // code plots: their SVG's text and strokes, scaled to the printed size
  for (const o of page.objects) {
    if (!o.plot) continue;
    const a = ASSETS[o.iconId], k = (o.w / a.vw) * (PRINT_MM / page.width) / 25.4 * 72;   // svg pt -> printed pt
    const fsz = [...a.svg.matchAll(/font-size:\s*([\d.]+)px/g)].map((m) => +m[1] * k);
    const sw = [...a.svg.matchAll(/stroke-width:\s*([\d.]+)/g)].map((m) => +m[1] * k).filter((v) => v > 0.01);
    for (const v of fsz) { fmin = Math.min(fmin, v); fmax = Math.max(fmax, v); if (v < 5.95 || v > 8.05) bad.push(`plot text ${v.toFixed(2)} pt (${o.name})`); }
    for (const v of sw) { smin = Math.min(smin, v); smax = Math.max(smax, v); if (v < 0.495 || v > 1.505) bad.push(`plot stroke ${v.toFixed(2)} pt (${o.name})`); }
    fonts.add('arial');
  }
  if ([...fonts].some((x) => x !== 'arial')) bad.push(`fonts other than Arial: ${[...fonts].join(', ')}`);
  console.log(`Cell Press check: ${PRINT_MM} x ${hMM.toFixed(1)} mm; text ${fmin.toFixed(1)}-${fmax.toFixed(1)} pt; `
    + `strokes ${smin.toFixed(2)}-${smax.toFixed(2)} pt; fonts ${[...fonts].join(', ')}`);
  if (bad.length) throw new Error('Cell Press requirements not met:\n  ' + bad.join('\n  '));
}

// SciCanvas's own Check tab (journal.js checkFigure and every module that extends it: panel overruns,
// line ends, legends, linked data), run with the Cell Press full-width preset.
async function sciCheck(app, page) {
  app.state.doc = { version: 2, pages: [page], uploads: [], assets: ASSETS, plotStyle: PLOT_STYLE, plotPrelude: PRELUDE,
    journal: { id: 'cell', widthMM: PRINT_MM, commercial: false, guide: false } };
  app.state.pageIndex = 0;
  // images: their pixel width, so the resolution check need not load them (the SHMtrace mark is 1,024 px)
  const each = (l) => l.forEach((o) => { if (o.type === 'image' && !o.nw) o.nw = 1024; if (o.children) each(o.children); });
  each(page.objects);
  const run = () => Promise.race([app.check(page), new Promise((_, no) => setTimeout(() => no(new Error('SciCanvas check timed out')), 20000))]);
  let res = await run();
  // objects 0.5-3 px out of line with a neighbour look like a mistake: apply SciCanvas's own "Line them up" fix
  // and check again
  // (repeated: lining one pair up can bring a neighbour within range of another)
  for (let pass = 0; pass < 4; pass++) {
    const align = res.issues.filter((x) => x.fixLabel === 'Line them up' && typeof x.fix === 'function');
    if (!align.length) break;
    align.forEach((x) => x.fix());
    console.log(`  aligned: ${align.map((x) => x.msg.split(',')[0]).join('; ')}`);
    res = await run();
  }
  const byId = new Map(); const walk = (l) => l.forEach((o) => { byId.set(o.id, o); if (o.children) walk(o.children); }); walk(page.objects);
  const shown = res.issues.filter((x) => x.sev !== 'info');
  for (const x of res.issues) {
    const o = x.id && byId.get(x.id);
    const pt = (e) => (e && e.x != null ? `${Math.round(e.x)},${Math.round(e.y)}` : e && e.id ? `#${(byId.get(e.id) || {}).type || '?'}` : '?');
    const where = !o ? '' : o.type === 'connector' ? ` [line ${pt(o.from)} -> ${pt(o.to)}]` : ` [${o.type} at ${Math.round(o.x)},${Math.round(o.y)}]`;
    console.log(`  SciCanvas check ${x.sev}: ${x.kind}: ${x.msg}${where}`);
    for (const pr of x.pairs || []) {
      const A = byId.get(pr.a), B = byId.get(pr.b), at = (q) => (q ? `${q.type} ${Math.round(q.x)},${Math.round(q.y)} ${Math.round(q.w)}x${Math.round(q.h)}` : '?');
      console.log(`    ${pr.what} off by ${pr.off}: ${at(A)} | ${at(B)}`);
    }
  }
  console.log(`SciCanvas check: ${shown.length} problems, ${res.issues.length - shown.length} notes`);
  return shown.length;
}

// ---------- Code plots ----------
// SciCanvas 2.1 renders a plot at its printed size: the page's width maps to the journal width, so a plot w page
// units wide is drawn w * K / 96 inches wide (K = printed px per page unit) and its style is in printed points.
// The plot records K (plot.k), which is how the app knows it was drawn for this print width.
const K = (PRINT_MM / 25.4 * 96) / 1400;              // 0.47 for these figures
const ASSETS = {};
const PLOT_STYLE = { family: 'Arial', fontSize: 7, tickSize: 6.6, labelSize: 7, titleSize: 7, legendSize: 6.6, minScriptPt: 6,
  axisLineWidth: 0.8, lineWidth: 1.4, markerSize: 4,
  palette: ['#274b7e', '#7a4f9a', '#1f6f8b', '#a3621b', '#9a3439', '#4a5562'] };
// The figure's shared plot code (doc.plotPrelude): run before every plot, in the app as here
const PRELUDE = { python: fs.readFileSync(path.join(__dirname, 'plots_common.py'), 'utf8') };
let _plotSeq = 0;
function runPlot(code, w, h) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'crm-plot-'));
  try {
    fs.writeFileSync(path.join(dir, 'code.py'), code);
    const { palette, family, ...sizes } = PLOT_STYLE;
    fs.writeFileSync(path.join(dir, 'params.json'), JSON.stringify({ width: w * K / 96, height: h * K / 96, ...sizes, palette, family }));
    fs.writeFileSync(path.join(dir, 'prelude.py'), PRELUDE.python);
    fs.writeFileSync(path.join(dir, 'data.json'), '{}');
    const worker = path.join(SCI, 'scripts', 'plot-worker.py');
    execFileSync(PYTHON, ['-c', 'import importlib.util, sys\n'
      + 'spec = importlib.util.spec_from_file_location("plot_worker", sys.argv[1]); m = importlib.util.module_from_spec(spec)\n'
      + 'spec.loader.exec_module(m); m.run(sys.argv[2])', worker, dir], { stdio: ['ignore', 'ignore', 'inherit'], env: { ...process.env, MPLBACKEND: 'Agg' } });
    const err = path.join(dir, 'error.txt');
    if (fs.existsSync(err)) throw new Error('plot failed:\n' + fs.readFileSync(err, 'utf8'));
    return fs.readFileSync(path.join(dir, 'out.svg'), 'utf8');
  } finally { fs.rmSync(dir, { recursive: true, force: true }); }
}
// The asset SciCanvas's normalizeSvg makes of an SVG, for matplotlib's output: the root unwrapped, ids made unique,
// the one CSS rule (line joins and caps) written onto the shapes.
function plotAssetFrom(svgText, normalizeColors) {
  const vb = /viewBox="([^"]+)"/.exec(svgText)[1];
  const [, , vw, vh] = vb.split(/[\s,]+/).map(Number);
  let inner = svgText.replace(/^[\s\S]*?<svg[^>]*>/, '').replace(/<\/svg>\s*$/, '')
    .replace(/<metadata>[\s\S]*?<\/metadata>/, '').replace(/<style[^>]*>[\s\S]*?<\/style>/g, '');
  const prefix = `cp${(_plotSeq++).toString(36)}_`;
  inner = inner.replace(/id="([^"]+)"/g, `id="${prefix}$1"`).replace(/href="#([^"]+)"/g, `href="#${prefix}$1"`)
    .replace(/url\(#([^)]+)\)/g, `url(#${prefix}$1)`);
  inner = `<g stroke-linejoin="round" stroke-linecap="butt">${inner}</g>`;
  return { svg: normalizeColors(inner), vb, vw, vh, raw: svgText, name: 'Python plot', pack: 'codeplot', license: 'own', author: 'You' };
}

async function main() {
  const [builderPath, stem] = process.argv.slice(2);
  const { app, ctx } = loadApp();
  // SciCanvas 2.0 boots asynchronously and loads a fresh document when it does: let that happen first, or it
  // replaces the document (and the code plots' assets) the check and the render use
  for (let i = 0; i < 5; i++) await new Promise((r) => setImmediate(r));
  // the SHMtrace mark, for builders that place it (LOGO); in the app, paste it as an image instead
  const logo = path.join(__dirname, '..', '..', '..', '..', 'packaging', 'icon', 'icon.svg');
  ctx.LOGO = 'data:image/svg+xml;base64,' + fs.readFileSync(logo).toString('base64');
  // SHMtrace's read accounting per library (READS.jh4 / READS.n114), from the analysis outputs
  const figs = path.join(__dirname, '..', '..', 'figures');
  ctx.READS = {};
  for (const [k, run] of [['jh4', 'jh4_shmtrace'], ['n114', 'n114_shmtrace']]) {
    const f = path.join(figs, run, 'table_s1.json');
    if (fs.existsSync(f)) ctx.READS[k] = JSON.parse(fs.readFileSync(f, 'utf8'));
  }
  ctx.PLOT = ({ code, x, y, w, h, name = 'Python plot' }) => {
    const a = plotAssetFrom(runPlot(code, w, h), (m) => app.normalizeColors(m));
    const key = `upload:plot${Object.keys(ASSETS).length}`;
    ASSETS[key] = a;
    return { id: `plot${Object.keys(ASSETS).length}`, type: 'icon', iconId: key, x, y, w, h, rot: 0, name,
      plot: { lang: 'python', code, w, h, k: K, files: [], error: null, textEdits: {} } };
  };
  const builder = vm.runInNewContext(`(${fs.readFileSync(builderPath, 'utf8')})`, ctx);
  const page = builder();
  const missing = page.objects.filter((o) => o.type === 'icon' && !app.ICON_MAP[o.iconId] && !ASSETS[o.iconId]).map((o) => o.iconId);
  if (missing.length) throw new Error(`unknown icons: ${missing.join(', ')}`);
  check(page);
  await sciCheck(app, page);
  const svg = app.pageSvgString(page);
  fs.writeFileSync(`${stem}.svg`, svg);
  fs.writeFileSync(`${stem}.scifig`, JSON.stringify({ version: 2, pages: [page], uploads: [], assets: ASSETS, plotStyle: PLOT_STYLE, plotPrelude: PRELUDE,
    journal: { id: 'cell', widthMM: PRINT_MM, commercial: false, guide: false } }));
  // vector PDF at print size (fonts embedded by Chrome), and a 1,000-dpi TIFF
  const wMM = PRINT_MM, hMM = page.height * PRINT_MM / page.width;
  const html = `${stem}.render.html`;
  fs.writeFileSync(html, `<!doctype html><meta charset="utf-8"><style>@page{size:${wMM}mm ${hMM}mm;margin:0}
html,body{margin:0;background:#fff}svg{display:block;width:${wMM}mm;height:${hMM}mm}</style>${svg}`);
  execFileSync(CHROME, ['--headless=new', '--disable-gpu', '--no-pdf-header-footer',
    `--print-to-pdf=${path.resolve(stem)}.pdf`, `file://${path.resolve(html)}`], { stdio: 'ignore' });
  const scale = (wMM / 25.4 * TIFF_DPI) / page.width;
  fs.writeFileSync(html, `<!doctype html><meta charset="utf-8"><style>html,body{margin:0;background:#fff}
svg{display:block;width:${page.width}px;height:${page.height}px}</style>${svg}`);
  const png = `${path.resolve(stem)}_tmp.png`;
  execFileSync(CHROME, ['--headless=new', '--disable-gpu', '--hide-scrollbars', `--force-device-scale-factor=${scale}`,
    `--window-size=${page.width},${Math.ceil(page.height)}`, `--screenshot=${png}`, `file://${path.resolve(html)}`], { stdio: 'ignore' });
  fs.unlinkSync(html);
  execFileSync('python3', ['-c', `from PIL import Image; im = Image.open(${JSON.stringify(png)}).convert("RGB"); `
    + `im.save(${JSON.stringify(path.resolve(stem) + '_1000dpi.tif')}, dpi=(${TIFF_DPI}, ${TIFF_DPI}), compression="tiff_lzw"); `
    + `im.save(${JSON.stringify(path.resolve(stem) + '_preview.png')}, dpi=(300, 300))`]);
  fs.unlinkSync(png);
  console.log(`${stem}: ${page.objects.length} objects; PDF, ${TIFF_DPI}-dpi TIFF, SVG, scifig`);
}

main().catch((e) => { console.error(e.message); process.exit(1); });

