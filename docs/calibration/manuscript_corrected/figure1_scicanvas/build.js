// Build a SciCanvas figure from a builder script, outside the app.
//
//   node build.js <builder.js> <out-stem>      (SCICANVAS=~/SciCanvas by default)
//
// The builder is a `() => page` function written against SciCanvas's own `Make` API, the same
// form SciCanvas runs from its console, so it can also be pasted into the app. This loads the
// app's source the way SciCanvas's tests do (test/load-app.js), with one change: text is
// measured with Helvetica's glyph widths (SciCanvas draws `sans` as Helvetica) instead of the
// tests' 7-px stub, so centred and right-aligned labels land where the app puts them.
// Writes <out-stem>.scifig (open it in SciCanvas to edit), <out-stem>.svg (SciCanvas's renderer)
// and <out-stem>_600dpi.png (Chrome, headless, at 600 dpi for a 174-mm-wide page).
const fs = require('fs'), vm = require('vm'), path = require('path'), os = require('os');
const { execFileSync } = require('child_process');

const SCI = process.env.SCICANVAS || path.join(os.homedir(), 'SciCanvas');
const CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const PRINT_MM = 174;                                  // full-width figure

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
    + '\n;({ Make, makeGroup, pageSvgString, ICON_MAP })';
  return { app: vm.runInNewContext(code, ctx, { filename: 'scicanvas-bundle.js' }), ctx };
}

function main() {
  const [builderPath, stem] = process.argv.slice(2);
  const { app, ctx } = loadApp();
  // the SHMtrace mark, for builders that place it (LOGO); in the app, paste it as an image instead
  const logo = path.join(__dirname, '..', '..', '..', '..', 'packaging', 'icon', 'icon.svg');
  ctx.LOGO = 'data:image/svg+xml;base64,' + fs.readFileSync(logo).toString('base64');
  const builder = vm.runInNewContext(`(${fs.readFileSync(builderPath, 'utf8')})`, ctx);
  const page = builder();
  const missing = page.objects.filter((o) => o.type === 'icon' && !app.ICON_MAP[o.iconId]).map((o) => o.iconId);
  if (missing.length) throw new Error(`unknown icons: ${missing.join(', ')}`);
  const svg = app.pageSvgString(page);
  fs.writeFileSync(`${stem}.svg`, svg);
  fs.writeFileSync(`${stem}.scifig`, JSON.stringify({ version: 2, pages: [page], uploads: [], assets: [] }));
  // 600 dpi at the printed width
  const scale = (PRINT_MM / 25.4 * 600) / page.width;
  const html = `${stem}.render.html`;
  fs.writeFileSync(html, `<!doctype html><meta charset="utf-8"><style>html,body{margin:0;background:#fff}
svg{display:block;width:${page.width}px;height:${page.height}px}</style>${svg}`);
  execFileSync(CHROME, ['--headless=new', '--disable-gpu', '--hide-scrollbars', `--force-device-scale-factor=${scale}`,
    `--window-size=${page.width},${page.height}`, `--screenshot=${path.resolve(stem)}_600dpi.png`,
    `file://${path.resolve(html)}`], { stdio: 'ignore' });
  fs.unlinkSync(html);
  console.log(`${stem}: ${page.objects.length} objects, ${page.width} x ${page.height}; png at ${scale.toFixed(2)}x`);
}

main();
