// Build a SciCanvas figure from a builder script, outside the app.
//
//   node build.js <builder.js> <out-stem>      (SCICANVAS=~/SciCanvas by default)
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

const SCI = process.env.SCICANVAS || path.join(os.homedir(), 'SciCanvas');
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
    + '\n;({ Make, makeGroup, pageSvgString, ICON_MAP, state, check: (p) => checkFigure(p) })';
  return { app: vm.runInNewContext(code, ctx, { filename: 'scicanvas-bundle.js' }), ctx };
}

// Cell Press production requirements, measured at print size
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
      if (size < 5.95 || size > 8.05) bad.push(`text ${size.toFixed(2)} pt: ${String(o.text || o.label).slice(0, 30)}`);
    }
    const w = o.type === 'connector' ? o.width : (o.stroke && o.stroke !== 'none' ? o.strokeWidth : 0);
    if (w) {
      smin = Math.min(smin, pt(w)); smax = Math.max(smax, pt(w));
      if (pt(w) < 0.495 || pt(w) > 1.505) bad.push(`stroke ${pt(w).toFixed(2)} pt (${o.type})`);
    }
    // a pure grey fill (r = g = b, not white) must be 10-80% grey
    const hex = /^#([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(o.fill || '');
    if (hex && hex[1] === hex[2] && hex[2] === hex[3] && hex[1].toLowerCase() !== 'ff') {
      const grey = 1 - parseInt(hex[1], 16) / 255;
      if (grey < 0.1 || grey > 0.8) bad.push(`grey fill ${o.fill} is ${(grey * 100).toFixed(0)}%`);
    }
  });
  walk(page.objects);
  if ([...fonts].some((x) => x !== 'arial')) bad.push(`fonts other than Arial: ${[...fonts].join(', ')}`);
  console.log(`Cell Press check: ${PRINT_MM} x ${hMM.toFixed(1)} mm; text ${fmin.toFixed(1)}-${fmax.toFixed(1)} pt; `
    + `strokes ${smin.toFixed(2)}-${smax.toFixed(2)} pt; fonts ${[...fonts].join(', ')}`);
  if (bad.length) throw new Error('Cell Press requirements not met:\n  ' + bad.join('\n  '));
}

// SciCanvas's own Check tab (journal.js checkFigure and every module that extends it: panel overruns,
// line ends, legends, linked data), run with the Cell Press full-width preset.
async function sciCheck(app, page) {
  app.state.doc = { pages: [page], journal: { id: 'cell', widthMM: PRINT_MM, commercial: false, guide: false } };
  app.state.pageIndex = 0;
  // images: their pixel width, so the resolution check need not load them (the SHMtrace mark is 1,024 px)
  const each = (l) => l.forEach((o) => { if (o.type === 'image' && !o.nw) o.nw = 1024; if (o.children) each(o.children); });
  each(page.objects);
  const run = () => Promise.race([app.check(page), new Promise((_, no) => setTimeout(() => no(new Error('SciCanvas check timed out')), 20000))]);
  let res = await run();
  const byId = new Map(); const walk = (l) => l.forEach((o) => { byId.set(o.id, o); if (o.children) walk(o.children); }); walk(page.objects);
  const shown = res.issues.filter((x) => x.sev !== 'info');
  for (const x of res.issues) {
    const o = x.id && byId.get(x.id);
    const pt = (e) => (e && e.x != null ? `${Math.round(e.x)},${Math.round(e.y)}` : e && e.id ? `#${(byId.get(e.id) || {}).type || '?'}` : '?');
    const where = !o ? '' : o.type === 'connector' ? ` [line ${pt(o.from)} -> ${pt(o.to)}]` : ` [${o.type} at ${Math.round(o.x)},${Math.round(o.y)}]`;
    console.log(`  SciCanvas check ${x.sev}: ${x.kind}: ${x.msg}${where}`);
  }
  console.log(`SciCanvas check: ${shown.length} problems, ${res.issues.length - shown.length} notes`);
  return shown.length;
}

async function main() {
  const [builderPath, stem] = process.argv.slice(2);
  const { app, ctx } = loadApp();
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
  const builder = vm.runInNewContext(`(${fs.readFileSync(builderPath, 'utf8')})`, ctx);
  const page = builder();
  const missing = page.objects.filter((o) => o.type === 'icon' && !app.ICON_MAP[o.iconId]).map((o) => o.iconId);
  if (missing.length) throw new Error(`unknown icons: ${missing.join(', ')}`);
  check(page);
  await sciCheck(app, page);
  const svg = app.pageSvgString(page);
  fs.writeFileSync(`${stem}.svg`, svg);
  fs.writeFileSync(`${stem}.scifig`, JSON.stringify({ version: 2, pages: [page], uploads: [], assets: [] }));
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
