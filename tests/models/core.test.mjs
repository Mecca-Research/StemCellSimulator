import test from 'node:test';
import assert from 'node:assert/strict';
import { RNG } from '../../js/models/core/rng.js';
import { integrate, rk45, gillespie, hill, hillR } from '../../js/models/core/ode.js';
import { eigSym, pca } from '../../js/models/core/linalg.js';
import { SpatialHash } from '../../js/models/core/spatial.js';

test('RNG is deterministic and roughly uniform / normal', () => {
  const a = new RNG(42), b = new RNG(42);
  for (let i = 0; i < 100; i++) assert.equal(a.next(), b.next());
  const r = new RNG(1);
  let s = 0, s2 = 0;
  const n = 20000;
  for (let i = 0; i < n; i++) { const x = r.normal(); s += x; s2 += x * x; }
  assert.ok(Math.abs(s / n) < 0.03);
  assert.ok(Math.abs(s2 / n - 1) < 0.05);
});

test('RK4 reproduces exponential decay to 1e-8', () => {
  const k = 0.7;
  const { y } = integrate((t, u, du) => { du[0] = -k * u[0]; }, [2], 0, 5, 0.01);
  assert.ok(Math.abs(y.at(-1)[0] - 2 * Math.exp(-k * 5)) < 1e-8);
});

test('adaptive RK45 solves the harmonic oscillator', () => {
  const res = rk45((t, u, du) => { du[0] = u[1]; du[1] = -u[0]; }, [1, 0], 0, 2 * Math.PI, { rtol: 1e-9, atol: 1e-12 });
  assert.ok(Math.abs(res.y[0] - 1) < 1e-6);
  assert.ok(Math.abs(res.y[1]) < 1e-6);
});

test('Gillespie birth-death mean matches k/g', () => {
  const rng = new RNG(3);
  const k = 20, g = 0.5;
  let sum = 0, n = 0;
  for (let rep = 0; rep < 40; rep++) {
    const { x } = gillespie([0], [
      { rate: () => k, change: [1] },
      { rate: (x) => g * x[0], change: [-1] },
    ], 30, rng);
    sum += x[0]; n++;
  }
  assert.ok(Math.abs(sum / n - k / g) < 4, `mean ${sum / n}`);
});

test('Hill functions', () => {
  assert.equal(hill(1, 1, 2), 0.5);
  assert.equal(hillR(1, 1, 2), 0.5);
  assert.ok(hill(10, 1, 4) > 0.999);
});

test('Jacobi eigen-decomposition of a symmetric matrix', () => {
  const M = [[4, 1, 0], [1, 3, 1], [0, 1, 2]];
  const { values, vectors } = eigSym(M);
  for (let k = 0; k < 3; k++) {
    const v = vectors[k];
    const Mv = M.map((r) => r.reduce((s, x, j) => s + x * v[j], 0));
    for (let i = 0; i < 3; i++) assert.ok(Math.abs(Mv[i] - values[k] * v[i]) < 1e-9);
  }
  assert.ok(values[0] >= values[1] && values[1] >= values[2]);
});

test('PCA recovers the dominant direction', () => {
  const r = new RNG(9);
  const X = Array.from({ length: 500 }, () => { const t = r.normal() * 5; return [t, 2 * t + r.normal() * 0.1, r.normal() * 0.1]; });
  const { components, explained } = pca(X, 2);
  const c = components[0];
  const cos = Math.abs(c[0] * 1 + c[1] * 2) / Math.sqrt(5);
  assert.ok(cos > 0.999);
  assert.ok(explained[0] > 0.99);
});

test('spatial hash finds all neighbours within one cell size', () => {
  const h = new SpatialHash(2);
  const pts = [[0, 0, 0], [1.5, 0, 0], [5, 5, 5], [-1.9, 0.1, 0]];
  pts.forEach((p, i) => h.insert(i, ...p));
  const found = new Set();
  h.near(0, 0, 0, (i) => found.add(i));
  assert.ok(found.has(0) && found.has(1) && found.has(3));
  assert.ok(!found.has(2));
});

test('spatial hash visits each neighbour exactly once', () => {
  const h = new SpatialHash(14);
  h.insert(0, 0, 0, 0); h.insert(1, 4, 0, 0); h.insert(2, -20, 13, 3);
  const seen = [];
  h.near(0, 0, 0, (i) => seen.push(i));
  assert.deepEqual(seen.sort(), [0, 1]);
});
