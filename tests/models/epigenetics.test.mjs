import test from 'node:test';
import assert from 'node:assert/strict';
import { integrate } from '../../js/models/core/ode.js';
import { marks, circuit, landscape, sampleU } from '../../js/models/pathways/epigenetics.js';

test('epigenetic marks relax to their closed-form steady state', () => {
  const p = { ...marks.defaults, DNMT: 0.8, HAT: 0.3, HDAC: 0.9 };
  const y = integrate(marks.rhs(p), [0, 1], 0, 40, 0.01).y.at(-1);
  const ss = marks.steadyState(p);
  assert.ok(Math.abs(y[0] - ss[0]) < 1e-6 && Math.abs(y[1] - ss[1]) < 1e-6);
});

test('strong self-activation: multipotent + two committed attractors; weak: bistable', () => {
  const early = circuit.attractors(circuit.params(0));
  const late = circuit.attractors(circuit.params(1));
  assert.equal(early.length, 3, JSON.stringify(early));
  const mid = early[1];
  assert.ok(Math.abs(mid[0] - mid[1]) < 0.05, 'central attractor is symmetric (co-expression)');
  assert.equal(late.length, 2, JSON.stringify(late));
  assert.ok(Math.abs(late[0][0] - late[1][1]) < 0.05, 'the two fates are mirror images');
});

test('landscape minima coincide with the deterministic attractors', () => {
  const p = circuit.params(1, 0.5, 0.5);
  const L = landscape(p, { walkers: 300, steps: 700 });
  const att = circuit.attractors(p);
  const saddle = sampleU(L, 1.0, 1.0);
  for (const [x, y] of att) assert.ok(sampleU(L, x, y) < saddle - 1, `U at attractor (${x.toFixed(2)}, ${y.toFixed(2)})`);
});
