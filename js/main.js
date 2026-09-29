// App shell: scene registry, routing, transport controls and lifecycle.
import * as THREE from 'three';
import { Stage } from './engine/stage.js';
import { Panel, setLegend, el } from './engine/ui.js';
import { releaseTree } from './engine/materials.js';
import * as assets from './engine/assets.js';
import { SCENES } from './scenes/index.js';

const $ = (id) => document.getElementById(id);
const debug = (window.__sim = { ready: false, sceneId: null, frames: 0, errors: [], stage: null });
window.addEventListener('error', (e) => debug.errors.push(String(e.message)));
window.addEventListener('unhandledrejection', (e) => debug.errors.push(String(e.reason?.message ?? e.reason)));

let stage;
try {
  stage = new Stage($('gl'));
} catch (err) {
  $('loading').innerHTML = '<span>WebGL 2 is required to run the simulator.</span>';
  throw err;
}
debug.stage = stage;
const panel = new Panel($('scene-controls'));
let active = null; // { def, instance }
let activating = null;

// ------------------------------------------------------------------ sidebar
function buildSidebar() {
  const nav = $('sidebar');
  let group = null;
  for (const s of SCENES) {
    if (s.group !== group) {
      group = s.group;
      nav.append(el('div', { class: 'nav-group', text: group }));
    }
    const b = el('button', { class: 'nav-item', 'data-id': s.id, onclick: () => { location.hash = s.id; document.body.classList.remove('nav-open'); } },
      [el('span', { text: s.title }), el('small', { text: s.summary })]);
    nav.append(b);
  }
}

function setStatus(text) {
  $('sim-time').textContent = text;
}

// --------------------------------------------------------------- lifecycle
async function activate(hash) {
  // "#scene-id?view=x&..." -> scene id + query parameters for the scene
  const [id, qs = ''] = hash.split('?');
  const query = new URLSearchParams(qs);
  const def = SCENES.find((s) => s.id === id) ?? SCENES[0];
  const key = `${def.id}?${query}`;
  if (activating === key) return;
  activating = key;
  $('loading').hidden = false;
  document.querySelectorAll('.nav-item').forEach((b) => b.classList.toggle('active', b.dataset.id === def.id));
  $('scene-title').textContent = def.title;
  $('scene-order').textContent = def.order ?? '';
  document.title = `${def.title} · Stem Cell Simulator`;

  if (active) {
    try { active.instance?.dispose?.(); } catch (e) { console.error(e); }
    releaseTree(stage.root);
    stage.root.clear();
  }
  active = null;
  panel.clear();
  setLegend($('hud-legend'), []);
  stage.setPickables([]);
  stage.clipAxis.set(0, 0, 1);
  $('clip-on').checked = false;
  $('clip-pos').value = 0;
  stage.setClip(false);
  stage.simTime = 0;
  stage.controls.autoRotate = false;
  stage.setBloomScale(1);
  setStatus('t = 0');

  const about = $('scene-about');
  about.replaceChildren();
  try {
    const mod = (await def.load()).default;
    if (activating !== key) return;
    about.append(el('h2', { text: def.group }));
    about.insertAdjacentHTML('beforeend', mod.about ?? '');
    if (mod.paperRef) about.append(el('div', { class: 'paper-ref', html: mod.paperRef }));
    const ctx = {
      THREE, stage, root: stage.root, panel, assets, setStatus, query,
      legend: (items) => setLegend($('hud-legend'), items),
    };
    const instance = await mod.create(ctx);
    if (activating !== key) { instance?.dispose?.(); return; }
    active = { def, instance, mod };
    debug.sceneId = def.id;
    debug.ready = true;
  } catch (err) {
    console.error(err);
    debug.errors.push(`${def.id}: ${err.message}`);
    about.append(el('p', { class: 'note', text: `This scene failed to load: ${err.message}` }));
  } finally {
    if (activating === key) activating = null;
    $('loading').hidden = true;
  }
}

stage.onFrame = (dt) => {
  if (!active?.instance) return;
  stage.simTime += dt;
  try {
    active.instance.update?.(dt, stage.simTime);
  } catch (err) {
    console.error(err);
    debug.errors.push(`${active.def.id} update: ${err.message}`);
    stage.running = false;
  }
  if (stage.frames % 3 === 0) panel.drawCharts();
  debug.frames = stage.frames;
};

// ---------------------------------------------------------------- controls
const playBtn = $('btn-play');
function setRunning(r) {
  stage.running = r;
  playBtn.innerHTML = r ? '&#10074;&#10074;' : '&#9654;';
}
playBtn.addEventListener('click', () => setRunning(!stage.running));
$('btn-reset').addEventListener('click', () => {
  if (active?.instance?.reset) { active.instance.reset(); stage.simTime = 0; }
  else if (active) { active = null; activating = null; activate(location.hash.slice(1)); }
});
$('speed').addEventListener('input', (e) => { stage.speed = parseFloat(e.target.value); });
const syncClip = () => stage.setClip($('clip-on').checked, parseFloat($('clip-pos').value));
$('clip-on').addEventListener('change', syncClip);
$('clip-pos').addEventListener('input', () => { if (!$('clip-on').checked) $('clip-on').checked = true; syncClip(); });
$('btn-shot').addEventListener('click', () => {
  const a = el('a', { href: stage.screenshot(), download: `${active?.def.id ?? 'scene'}.png` });
  a.click();
});

const modeButtons = document.querySelectorAll('#mode-switch button');
function setMode(mode) {
  stage.setMode(mode);
  modeButtons.forEach((b) => b.classList.toggle('active', b.dataset.mode === mode));
  active?.instance?.onMode?.(mode);
}
modeButtons.forEach((b) => b.addEventListener('click', () => setMode(b.dataset.mode)));
debug.setMode = setMode;

$('nav-toggle').addEventListener('click', () => document.body.classList.toggle('nav-open'));
$('panel-toggle').addEventListener('click', () => document.body.classList.toggle('panel-open'));

window.addEventListener('keydown', (e) => {
  if (e.target.matches('input, select, textarea')) return;
  if (e.code === 'Space') { e.preventDefault(); setRunning(!stage.running); }
  else if (e.key === 'r') $('btn-reset').click();
  else if (e.key === '1') setMode('physical');
  else if (e.key === '2') setMode('confocal');
  else if (e.key === '3') setMode('histology');
});

// hover picking -> tooltip
const tip = $('hud-tip');
let lastPick = 0;
$('gl').addEventListener('pointermove', (e) => {
  const now = performance.now();
  if (now - lastPick < 60) return;
  lastPick = now;
  const info = stage.pick(e.clientX, e.clientY);
  if (!info) { tip.hidden = true; return; }
  const r = $('viewport').getBoundingClientRect();
  tip.hidden = false;
  tip.textContent = info;
  tip.style.left = `${e.clientX - r.left}px`;
  tip.style.top = `${e.clientY - r.top}px`;
});
$('gl').addEventListener('pointerleave', () => { tip.hidden = true; });

buildSidebar();
window.addEventListener('hashchange', () => activate(location.hash.slice(1)));
stage.start();
activate(location.hash.slice(1) || SCENES[0].id);
