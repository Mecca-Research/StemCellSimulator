// Record short animated previews of scenes for the README.
//   node tools/capture_previews.mjs [spec=mode ...]      (default: the README set)
// Frames are captured from the 3D viewport in headless Chromium and assembled
// into animated WebP files in docs/images/ by tools/assemble_previews.py.
import { chromium } from 'playwright';
import { mkdir, rm } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { serve } from '../tests/e2e/serve.mjs';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const DEFAULT = [
  'colony=confocal', 'colony=physical', 'signaling=physical', 'hematopoiesis=physical', 'differentiation=physical',
  'morphogenesis=physical', 'proteins=physical', 'landscape=physical', 'network=physical', 'atlas=physical',
  'kidney=histology', 'gallery?view=mitosis=physical',
];
const specs = process.argv.slice(2).length ? process.argv.slice(2) : DEFAULT;
const FRAMES = Number(process.env.FRAMES ?? 36);
const preinstalled = '/opt/pw-browsers/chromium-1194/chrome-linux/chrome';

const { server, url } = await serve(root);
const browser = await chromium.launch({
  executablePath: process.env.CHROMIUM_PATH || (existsSync(preinstalled) ? preinstalled : undefined),
  args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
});
const page = await browser.newPage({ viewport: { width: 1240, height: 760 } });
for (const spec of specs) {
  const cut = spec.lastIndexOf('=');
  const scene = spec.slice(0, cut), mode = spec.slice(cut + 1);
  const name = `${scene.replace(/[?&=]/g, '-')}-${mode}`;
  const dir = join(root, 'tests', 'e2e', 'out', 'frames', name);
  await rm(dir, { recursive: true, force: true });
  await mkdir(dir, { recursive: true });
  await page.goto('about:blank');
  await page.goto(`${url}#${scene}`);
  await page.waitForFunction((id) => window.__sim?.ready && window.__sim.sceneId === id, scene.split('?')[0], { timeout: 120000 });
  await page.addStyleTag({ content: '#hud, #transport, #loading { display: none !important; }' });
  await page.evaluate((m) => { window.__sim.setMode(m); window.__sim.stage.speed = 3; window.__sim.stage.setFixedPixelRatio(1); }, mode);
  // let the simulation develop before recording
  await page.waitForFunction((f) => window.__sim.frames > f, (await page.evaluate(() => window.__sim.frames)) + 60, { timeout: 180000 });
  const vp = page.locator('#viewport');
  for (let i = 0; i < FRAMES; i++) {
    const f0 = await page.evaluate(() => window.__sim.frames);
    await page.waitForFunction((f) => window.__sim.frames > f + 2, f0, { timeout: 60000 });
    await vp.screenshot({ path: join(dir, `${String(i).padStart(3, '0')}.png`) });
  }
  console.log(`captured ${name}`);
}
await browser.close();
server.close();
execFileSync('python3', [join(root, 'tools', 'assemble_previews.py')], { stdio: 'inherit' });
