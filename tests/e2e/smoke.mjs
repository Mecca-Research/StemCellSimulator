// Headless browser smoke test: every scene must load without errors, render a
// non-blank WebGL frame and keep simulating. Screenshots go to tests/e2e/out/.
//   node tests/e2e/smoke.mjs [sceneId ...]
import { chromium } from 'playwright';
import { mkdir } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { serve } from './serve.mjs';
import { SCENES } from '../../js/scenes/index.js';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..', '..');
const outDir = join(here, 'out');
await mkdir(outDir, { recursive: true });

const only = process.argv.slice(2).filter((a) => !a.startsWith('--'));
const modes = process.argv.includes('--modes') ? ['physical', 'confocal', 'histology'] : ['physical'];
const ids = only.length ? only : SCENES.map((s) => s.id);
const { server, url } = await serve(root);
const preinstalled = '/opt/pw-browsers/chromium-1194/chrome-linux/chrome';
const browser = await chromium.launch({
  executablePath: process.env.CHROMIUM_PATH || (existsSync(preinstalled) ? preinstalled : undefined),
  args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
});
const page = await browser.newPage({ viewport: { width: 1400, height: 860 } });
const errors = [];
page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
page.on('pageerror', (e) => errors.push(String(e)));

let failed = 0;
for (const spec of ids) {
  // spec is a scene id, optionally with a query: "gallery?view=skin"
  const id = spec.split('?')[0];
  const shot = spec.replace(/[?&=]/g, '-');
  // a fresh document per scene (a hash-only change would reuse the page)
  await page.goto('about:blank');
  const before = errors.length;
  const t0 = Date.now();
  await page.goto(`${url}#${spec}`);
  let ok = true, why = '';
  try {
    await page.waitForFunction((sid) => window.__sim?.ready && window.__sim.sceneId === sid, id, { timeout: 90000 });
    const f0 = await page.evaluate(() => window.__sim.frames);
    await page.waitForFunction((f) => window.__sim.frames > f + 45, f0, { timeout: 90000 });
  } catch (e) { ok = false; why = e.message.split('\n')[0]; }
  for (const mode of modes) {
    if (mode !== 'physical') {
      await page.evaluate((m) => window.__sim.setMode(m), mode);
      const f1 = await page.evaluate(() => window.__sim.frames);
      await page.waitForFunction((f) => window.__sim.frames > f + 20, f1, { timeout: 60000 }).catch(() => {});
    }
    const stats = await page.evaluate(() => {
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
    });
    await page.screenshot({ path: join(outDir, `${shot}${mode === 'physical' ? '' : '-' + mode}.png`) });
    if (stats.std < 2) { ok = false; why ||= `blank frame in ${mode} (std ${stats.std.toFixed(2)})`; }
    console.log(`   ${mode.padEnd(9)} mean ${stats.mean.toFixed(1)} std ${stats.std.toFixed(1)}`);
  }
  if (modes.length > 1) await page.evaluate(() => window.__sim.setMode('physical'));
  const simErrors = await page.evaluate(() => window.__sim?.errors ?? []);
  const newErrors = [...errors.slice(before), ...simErrors];
  if (newErrors.length) { ok = false; why ||= newErrors.join(' | ').slice(0, 400); }
  console.log(`${ok ? 'PASS' : 'FAIL'} ${spec.padEnd(16)} ${((Date.now() - t0) / 1000).toFixed(1)}s ${why}`);
  if (!ok) failed++;
  await page.evaluate(() => { window.__sim.errors.length = 0; });
}
await browser.close();
server.close();
if (failed) { console.error(`${failed} scene(s) failed`); process.exit(1); }
