import test from 'node:test';
import assert from 'node:assert/strict';
import { filtrateProfile, NEPHRON_SEGMENTS, SEGMENT_BOUNDS } from '../../js/models/organ/nephron.js';

test('filtrate volume falls monotonically and matches textbook segment fractions', () => {
  const f = filtrateProfile(0.5);
  let prev = 1.0001;
  for (let i = 0; i <= 200; i++) { const v = f(i / 200); assert.ok(v <= prev + 1e-12); prev = v; }
  assert.ok(Math.abs(f(0) - 1) < 1e-12);
  assert.ok(Math.abs(f(SEGMENT_BOUNDS[1]) - 0.35) < 1e-9, 'after the proximal tubule 35 % remains');
  assert.ok(Math.abs(f(SEGMENT_BOUNDS[3]) - 0.2) < 1e-9, 'the thick ascending limb reabsorbs no water');
});

test('ADH concentrates urine', () => {
  const dilute = filtrateProfile(0)(1), concentrated = filtrateProfile(1)(1);
  assert.ok(dilute > 10 * concentrated);
  assert.ok(concentrated < 0.01 && dilute > 0.1);
  assert.equal(NEPHRON_SEGMENTS.length, SEGMENT_BOUNDS.length - 1);
});
