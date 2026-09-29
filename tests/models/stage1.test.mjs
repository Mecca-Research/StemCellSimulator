import test from 'node:test';
import assert from 'node:assert/strict';
import { integrate } from '../../js/models/core/ode.js';
import { RNG } from '../../js/models/core/rng.js';
import { paperWnt, canonicalWnt } from '../../js/models/pathways/wnt.js';
import { collier, singleCell } from '../../js/models/pathways/notch.js';

test("paper's Wnt cascade converges to its closed-form steady state", () => {
  const p = paperWnt.defaults;
  const { y } = integrate(paperWnt.rhs(p), [0, 0, 0, 0], 0, 80, 0.01);
  const ss = paperWnt.steadyState(p);
  y.at(-1).forEach((v, i) => assert.ok(Math.abs(v - ss[i]) / ss[i] < 1e-3, `${paperWnt.species[i]} ${v} vs ${ss[i]}`));
});

test("paper's Wnt: W(t) matches the exact solution", () => {
  const p = paperWnt.defaults;
  const { t, y } = integrate(paperWnt.rhs(p), [0, 0, 0, 0], 0, 5, 0.01, { every: 50 });
  t.forEach((ti, i) => assert.ok(Math.abs(y[i][0] - paperWnt.W(ti, p)) < 1e-8));
});

test('canonical Wnt: ligand stabilises nuclear beta-catenin, Axin2 feedback limits it', () => {
  const run = (wntIn, kAx) => integrate(canonicalWnt.rhs({ ...canonicalWnt.defaults, wntIn, kAx }), canonicalWnt.initial(), 0, 200, 0.02).y.at(-1);
  const off = run(0, 0.6), on = run(1, 0.6), noFeedback = run(1, 0);
  assert.ok(on[4] > 3 * off[4], `nuclear bcat on ${on[4]} off ${off[4]}`);
  assert.ok(on[6] > off[6]);
  assert.ok(noFeedback[4] > on[4], 'Axin2 negative feedback lowers beta-catenin');
});

test('Collier lateral inhibition on a hexagonal sheet yields isolated senders', () => {
  // 12 x 12 hexagonal lattice (axial coordinates, periodic)
  const W = 12, n = W * W;
  const idx = (q, r) => ((r + W) % W) * W + ((q + W) % W);
  const nbrs = Array.from({ length: n }, (_, i) => {
    const q = i % W, r = Math.floor(i / W);
    return [[1, 0], [-1, 0], [0, 1], [0, -1], [1, -1], [-1, 1]].map(([dq, dr]) => idx(q + dq, r + dr));
  });
  const s = collier.init(n, new RNG(11), 0.1);
  for (let k = 0; k < 4000; k++) collier.step(s, nbrs, 0.05);
  const frac = collier.senderFraction(s);
  assert.ok(frac > 0.15 && frac < 0.4, `sender fraction ${frac}`);
  assert.ok(collier.senderAdjacency(s, nbrs) < 0.05, 'senders should not touch');
});

test('single-cell Notch: gamma-secretase inhibition (DAPT) abolishes NICD and Hes1', () => {
  const run = (kS3) => integrate(singleCell.rhs({ ...singleCell.defaults, kS3 }), singleCell.initial(), 0, 100, 0.02).y.at(-1);
  const on = run(1.2), dapt = run(0);
  assert.ok(on[4] > 0.1 && on[5] > 0.1);
  assert.ok(dapt[4] < 1e-6 && dapt[5] < 1e-6);
  assert.ok(dapt[6] > on[6], 'without Notch signalling the cell up-regulates its own Delta');
});
