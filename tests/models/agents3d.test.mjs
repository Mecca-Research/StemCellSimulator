import test from 'node:test';
import assert from 'node:assert/strict';
import { Tissue, Cell, doublingTime, mitosisStage } from '../../js/models/morpho/agents3d.js';

test('two overlapping cells relax to contact distance', () => {
  const t = new Tissue({ noise: 0, growth: false, adhesion: () => 0 });
  t.add(new Cell({ p: [0, 0, 0], R: 5 }));
  t.add(new Cell({ p: [4, 0, 0], R: 5 }));
  for (let i = 0; i < 400; i++) t.step(0.02);
  const d = Math.hypot(...t.cells[0].p.map((v, k) => v - t.cells[1].p[k]));
  assert.ok(Math.abs(d - 10) < 0.05, `distance ${d}`);
});

test('adhesion pulls separated cells into contact', () => {
  const t = new Tissue({ noise: 0, growth: false, adhesion: () => 3 });
  t.add(new Cell({ p: [0, 0, 0], R: 5 }));
  t.add(new Cell({ p: [11.5, 0, 0], R: 5 }));
  for (let i = 0; i < 2000; i++) t.step(0.02);
  const d = Math.hypot(...t.cells[0].p.map((v, k) => v - t.cells[1].p[k]));
  assert.ok(d < 10.2, `distance ${d}`);
});

test('division conserves volume and doubles the population once per cycle', () => {
  const t = new Tissue({ noise: 0, seed: 5 });
  t.add(new Cell({ p: [0, 0, 0], R0: 5, R: 5, cycle: 10 }));
  const v0 = 5 ** 3;
  const ts = [], ns = [];
  for (let i = 0; i < 3000; i++) { t.step(0.02); if (i % 50 === 0) { ts.push(t.t); ns.push(t.cells.length); } }
  // after 60 h with T ~ 10 h (+-10 %) the clone must be close to 2^6
  assert.ok(t.cells.length >= 32 && t.cells.length <= 128, `n = ${t.cells.length}`);
  const dt = doublingTime(ts.slice(20), ns.slice(20));
  assert.ok(dt > 8 && dt < 13, `doubling ${dt}`);
  // a fresh daughter carries half the mother's final volume
  const d = t.cells.find((c) => c.phase === 'G1' && c.phaseT < 0.1);
  if (d) assert.ok(Math.abs(d.R0 ** 3 - v0) / v0 < 0.3);
});

test('monolayer growth keeps cells on the substrate', () => {
  const t = new Tissue({ dims: 2, noise: 0.05, substrate: { y: 5, k: 6 }, seed: 2 });
  t.add(new Cell({ p: [0, 5, 0], R: 6, cycle: 8 }));
  for (let i = 0; i < 1500; i++) t.step(0.02);
  for (const c of t.cells) assert.ok(Math.abs(c.p[1] - 5) < 1.5);
  assert.ok(t.cells.length >= 8);
});

test('mitosis stages are ordered', () => {
  assert.equal(mitosisStage(0.05), 'prophase');
  assert.equal(mitosisStage(0.4), 'metaphase');
  assert.equal(mitosisStage(0.6), 'anaphase');
  assert.equal(mitosisStage(0.95), 'cytokinesis');
});
