import test from 'node:test';
import assert from 'node:assert/strict';
import { LINEAGES, GERM_LAYERS, TISSUES, ORGANS, SYSTEMS, SYSTEM_LINKS, allCellTypes, morphology } from '../../js/models/atlas/atlas.js';

test('cell map is complete and consistent', () => {
  const cells = allCellTypes();
  assert.ok(cells.length >= 150, `cell types: ${cells.length}`);
  const names = new Set(cells.map((c) => c.name));
  assert.equal(names.size, cells.length, 'no duplicate cell types');
  const layers = new Set(GERM_LAYERS.map((g) => g.id));
  for (const l of LINEAGES) assert.ok(layers.has(l.layer), l.id);
  for (const t of TISSUES) for (const c of t.cells) assert.ok(names.has(c), `tissue cell ${c}`);
  const systems = new Set(SYSTEMS.map((s) => s.name));
  for (const o of ORGANS) assert.ok(systems.has(o.system), o.name);
  for (const l of SYSTEM_LINKS) assert.ok(systems.has(l.a) && systems.has(l.b));
});

test('morphology classes', () => {
  assert.equal(morphology('Erythrocyte (red blood cell)'), 'rbc');
  assert.equal(morphology('Motor neuron'), 'star');
  assert.equal(morphology('Fibroblast'), 'spindle');
  assert.equal(morphology('Enterocyte'), 'columnar');
});
