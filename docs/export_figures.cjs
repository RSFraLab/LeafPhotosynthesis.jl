// Optional maintainer helper. The gallery itself needs no Node dependencies.
// Install Playwright as described in docs/README.md, then run:
// node docs/export_figures.cjs
const { chromium } = require('playwright');
const { mkdir, writeFile } = require('node:fs/promises');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');

(async () => {
  const browser = await chromium.launch({headless: true,
    ...(process.env.LEAF_DOCS_CHROMIUM ? {executablePath: process.env.LEAF_DOCS_CHROMIUM} : {})});
  try {
    const page = await browser.newPage({viewport: {width: 1280, height: 960}, deviceScaleFactor: 2});
    await page.goto(pathToFileURL(join(__dirname, 'examples.html')).href);
    await page.waitForFunction(() => document.getElementById('status').textContent.includes('converged'));
    const directory = join(__dirname, 'assets', 'figures');
    await mkdir(directory, {recursive: true});
    const figures = ['light', 'co2', 'temperature', 'conductance', 'gas-exchange', 'electron-use', 'carbon-turnover'];
    for (const id of figures) {
      const image = await page.evaluate(id => chartImage(Chart.getChart(id)), id);
      await writeFile(join(directory, `${id}.png`), Buffer.from(image.split(',')[1], 'base64'));
    }
    console.log(`Exported ${figures.length} Float64 Chart.js figures to docs/assets/figures.`);
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
