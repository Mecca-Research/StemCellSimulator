import test from 'node:test';
import assert from 'node:assert/strict';
import { allReactions, countsByOrder, buildGraph, embed, crosstalkSpecies, MODULES } from '../../js/models/network/registry.js';

test('every pathway module exports well-formed reactions with orders 1-4', () => {
  for (const [name, m] of Object.entries(MODULES)) assert.ok(Array.isArray(m.REACTIONS) && m.REACTIONS.length > 0, `${name} exports REACTIONS`);
  const rs = allReactions();
  const ids = new Set();
  for (const r of rs) {
    assert.ok([1, 2, 3, 4].includes(r.order), `${r.id} order ${r.order}`);
    assert.ok(r.reactants.length + r.products.length > 0, `${r.id} has species`);
    assert.ok(!ids.has(r.id)); ids.add(r.id);
  }
  const c = countsByOrder(rs);
  for (const k of [1, 2, 3, 4]) assert.ok(c[k] > 0, `order ${k} populated`);
});

test('network graph is bipartite and the embedding is finite', () => {
  const g = buildGraph(allReactions());
  for (const e of g.edges) assert.notEqual(g.nodes[e.a].kind, g.nodes[e.b].kind);
  const { positions, explained } = embed(g, { refine: 30 });
  assert.equal(positions.length, g.nodes.length);
  for (const p of positions) for (const v of p) assert.ok(Number.isFinite(v));
  assert.ok(explained[0] >= explained[1] && explained[1] >= explained[2]);
  assert.ok(crosstalkSpecies(g).length > 0, 'some species are shared between pathways');
});
