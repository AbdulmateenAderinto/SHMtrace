const sharp = require('/Users/abdulmateenaderinto/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp');
const path = require('path');
const base = path.join(__dirname, 'codex_version', 'figures', 'FigureS21_platforms');
sharp(base + '.svg', { density: 300 }).resize({ width: 4110 }).png().toFile(base + '.png')
  .then(() => console.log('Updated S21 rendered from its SVG'))
  .catch(e => { console.error(e); process.exit(1); });
