// Headless browser smoke test: every scene must load without errors, render a
// non-blank WebGL frame and keep simulating. Screenshots go to tests/e2e/out/.
//   node tests/e2e/smoke.mjs [sceneId ...]
// SITE_ROOT=<dir> serves another directory, e.g. the assembled GitHub Pages site.
import { chromium } from 'playwright';
import { mkdir } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join, resolve } from 'node:path';
import { serve } from './serve.mjs';
import { SCENES } from '../../js/scenes/index.js';
import { EXTRA_VIEWS } from './views.mjs';

const here = dirname(fileURLToPath(import.meta.url));
const root = process.env.SITE_ROOT ? resolve(process.env.SITE_ROOT) : join(here, '..', '..');
const outDir = join(here, 'out');
await mkdir(outDir, { recursive: true });

const only = process.argv.slice(2).filter((a) => !a.startsWith('--'));
const modes = process.argv.includes('--modes') ? ['physical', 'confocal', 'histology'] : ['physical'];
const ids = only.length ? only : [...SCENES.map((s) => s.id), ...(process.argv.includes('--all') ? EXTRA_VIEWS : [])];
const { server, url } = await serve(root);
const preinstalled = '/opt/pw-browsers/chromium-1194/chrome-linux/chrome';
const browser = await chromium.launch({
  executablePath: process.env.CHROMIUM_PATH || (existsSync(preinstalled) ? preinstalled : undefined),
  args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
});
const viewport = { width: 1400, height: 860 };

// Bound a browser call that has no timeout of its own (page.evaluate, close).
function within(ms, what, promise) {
  let timer;
  const limit = new Promise((_, reject) => { timer = setTimeout(() => reject(new Error(`${what} timed out after ${ms / 1000}s`)), ms); });
  return Promise.race([promise, limit]).finally(() => clearTimeout(timer));
}

async function checkScene(page, spec, id, errors) {
  await page.goto(`${url}#${spec}`);
  await page.waitForFunction((sid) => window.__sim?.ready && window.__sim.sceneId === sid, id, { timeout: 90000 });
  const f0 = await within(30000, 'reading the frame count', page.evaluate(() => window.__sim.frames));
  await page.waitForFunction((f) => window.__sim.frames > f + 45, f0, { timeout: 90000 });
  let why = '';
  for (const mode of modes) {
    if (mode !== 'physical') {
      await within(30000, `switching to ${mode}`, page.evaluate((m) => window.__sim.setMode(m), mode));
      const f1 = await within(30000, 'reading the frame count', page.evaluate(() => window.__sim.frames));
      await page.waitForFunction((f) => window.__sim.frames > f + 20, f1, { timeout: 60000 }).catch(() => {});
    }
    const stats = await within(60000, `measuring the ${mode} frame`, page.evaluate(() => {
      const c = document.getElementById('gl');
      const w = 240, h = Math.round((240 * c.height) / c.width);
      const cv = Object.assign(document.createElement('canvas'), { width: w, height: h });
      const g = cv.getContext('2d');
      g.drawImage(c, 0, 0, w, h);
      const d = g.getImageData(0, 0, w, h).data;
      let s = 0, s2 = 0, n = 0;
      for (let i = 0; i < d.length; i += 4) { const v = (d[i] + d[i + 1] + d[i + 2]) / 3; s += v; s2 += v * v; n++; }
      const mean = s / n;
      return { mean, std: Math.sqrt(Math.max(s2 / n - mean * mean, 0)) };
    }));
    const shot = spec.replace(/[?&=]/g, '-');
    await page.screenshot({ path: join(outDir, `${shot}${mode === 'physical' ? '' : '-' + mode}.png`), timeout: 60000 });
    if (stats.std < 2) why ||= `blank frame in ${mode} (std ${stats.std.toFixed(2)})`;
    console.log(`   ${mode.padEnd(9)} mean ${stats.mean.toFixed(1)} std ${stats.std.toFixed(1)}`);
  }
  const simErrors = await within(30000, 'reading scene errors', page.evaluate(() => window.__sim?.errors ?? []));
  const all = [...errors, ...simErrors];
  if (all.length) why ||= all.join(' | ').slice(0, 400);
  return why;
}

let failed = 0;
for (const spec of ids) {
  // spec is a scene id, optionally with a query: "gallery?view=skin".
  // Each scene gets its own browser context (and renderer process), so a scene
  // still inside a long software-GL frame cannot hold up or leak into the next.
  const id = spec.split('?')[0];
  const context = await browser.newContext({ viewport });
  const page = await context.newPage();
  const errors = [];
  page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
  page.on('pageerror', (e) => errors.push(String(e)));
  const t0 = Date.now();
  let why;
  try {
    why = await checkScene(page, spec, id, errors);
  } catch (e) {
    why = e.message.split('\n')[0];
  }
  await within(30000, 'closing the page', context.close()).catch((e) => console.error(`   ${e.message}`));
  console.log(`${why ? 'FAIL' : 'PASS'} ${spec.padEnd(16)} ${((Date.now() - t0) / 1000).toFixed(1)}s ${why}`);
  if (why) failed++;
}
await browser.close();
server.close();
if (failed) { console.error(`${failed} scene(s) failed`); process.exit(1); }
