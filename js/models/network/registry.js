// Reaction registry: every reaction the simulator models, grouped by the
// paper's four orders, turned into a species-reaction network and embedded in
// 3D with PCA of the adjacency matrix (paper, methodology step 5).
import * as wnt from '../pathways/wnt.js';
import * as notch from '../pathways/notch.js';
import * as hedgehog from '../pathways/hedgehog.js';
import * as mapk from '../pathways/mapk.js';
import * as jakstat from '../pathways/jakstat.js';
import * as hematopoiesis from '../pathways/hematopoiesis.js';
import * as myogenesis from '../pathways/myogenesis.js';
import * as neurogenesis from '../pathways/neurogenesis.js';
import * as hormone from '../pathways/hormone.js';
import * as crosstalk from '../pathways/crosstalk.js';
import * as epigenetics from '../pathways/epigenetics.js';
import * as ecm from '../morpho/ecm.js';
import * as gradient from '../morpho/gradient.js';
import * as protein from '../protein/assembly.js';
import * as nephron from '../organ/nephron.js';
import { powerPCA } from '../core/linalg.js';
import { RNG } from '../core/rng.js';

export const MODULES = { wnt, notch, hedgehog, mapk, jakstat, hematopoiesis, myogenesis, neurogenesis, hormone, crosstalk, epigenetics, ecm, gradient, protein, nephron };

/** Cumulative reaction counts reported in the paper for orders 1-4. */
export const PAPER_CUMULATIVE = { 1: 500, 2: 1051, 3: 2051, 4: 3021 };

export function allReactions() {
  const seen = new Map();
  for (const [module, m] of Object.entries(MODULES)) {
    for (const r of m.REACTIONS ?? []) {
      if (!r || seen.has(r.id)) continue;
      seen.set(r.id, { ...r, module, reactants: r.reactants ?? [], products: r.products ?? [], modifiers: r.modifiers ?? [] });
    }
  }
  return [...seen.values()];
}

export function countsByOrder(reactions) {
  const c = { 1: 0, 2: 0, 3: 0, 4: 0 };
  for (const r of reactions) if (c[r.order] !== undefined) c[r.order]++;
  return c;
}

const norm = (s) => s.trim();

/**
 * Bipartite species-reaction graph.
 * nodes: [{ id, kind: 'species'|'reaction', label, order, pathways:Set }]
 * edges: [{ a, b, role: 'substrate'|'product'|'modifier' }] (indices)
 */
export function buildGraph(reactions) {
  const nodes = [], index = new Map(), edges = [];
  const species = (name, r) => {
    const key = `s:${norm(name)}`;
    if (!index.has(key)) { index.set(key, nodes.length); nodes.push({ id: key, kind: 'species', label: norm(name), order: r.order, pathways: new Set() }); }
    const n = nodes[index.get(key)];
    n.order = Math.min(n.order, r.order);
    n.pathways.add(r.pathway);
    return index.get(key);
  };
  for (const r of reactions) {
    const ri = nodes.length;
    nodes.push({ id: `r:${r.id}`, kind: 'reaction', label: r.label, order: r.order, pathways: new Set([r.pathway]), reaction: r });
    for (const s of r.reactants) edges.push({ a: species(s, r), b: ri, role: 'substrate' });
    for (const s of r.products) edges.push({ a: ri, b: species(s, r), role: 'product' });
    for (const s of r.modifiers) edges.push({ a: species(s, r), b: ri, role: 'modifier' });
  }
  const degree = new Array(nodes.length).fill(0);
  for (const e of edges) { degree[e.a]++; degree[e.b]++; }
  nodes.forEach((n, i) => { n.degree = degree[i]; });
  return { nodes, edges };
}

/**
 * 3D layout: the first three principal components of the (symmetrised,
 * self-looped, two-step diffused) adjacency matrix, refined by a short
 * force-directed relaxation so neighbouring nodes do not overlap.
 */
export function embed(graph, { refine = 250, seed = 4, scale = 60 } = {}) {
  const n = graph.nodes.length;
  const A = Array.from({ length: n }, (_, i) => { const r = new Float64Array(n); r[i] = 1; return r; });
  for (const e of graph.edges) { A[e.a][e.b] = 1; A[e.b][e.a] = 1; }
  // two-step diffusion makes the principal axes follow pathway modules
  const A2 = A.map((row) => {
    const out = new Float64Array(n);
    for (let j = 0; j < n; j++) if (row[j]) { const rj = A[j]; for (let k = 0; k < n; k++) out[k] += rj[k]; }
    let s = 0; for (const v of out) s += v; return out.map((v) => v / (s || 1));
  });
  const { scores, explained } = powerPCA(A2, 3, { iters: 200, seed });
  let max = 1e-9;
  for (const s of scores) for (const v of s) max = Math.max(max, Math.abs(v));
  const pos = scores.map((s) => s.map((v) => (v / max) * scale));
  const pca = pos.map((p) => p.slice());
  // spring refinement (Fruchterman-Reingold style) anchored to the PCA layout
  const rng = new RNG(seed);
  for (const p of pos) for (let k = 0; k < 3; k++) p[k] += (rng.next() - 0.5) * 2;
  const kLen = scale * 0.12;
  const disp = pos.map(() => [0, 0, 0]);
  for (let it = 0; it < refine; it++) {
    const temp = scale * 0.05 * (1 - it / refine) + 0.05;
    for (const d of disp) d.fill(0);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const dx = pos[i][0] - pos[j][0], dy = pos[i][1] - pos[j][1], dz = pos[i][2] - pos[j][2];
      const d2 = dx * dx + dy * dy + dz * dz + 0.01;
      if (d2 > (kLen * 6) ** 2) continue;
      const f = (kLen * kLen) / d2;
      disp[i][0] += dx * f; disp[i][1] += dy * f; disp[i][2] += dz * f;
      disp[j][0] -= dx * f; disp[j][1] -= dy * f; disp[j][2] -= dz * f;
    }
    for (const e of graph.edges) {
      const a = pos[e.a], b = pos[e.b];
      const dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz) + 1e-6;
      const f = (d * d) / kLen / d;
      disp[e.a][0] -= dx * f * 0.5; disp[e.a][1] -= dy * f * 0.5; disp[e.a][2] -= dz * f * 0.5;
      disp[e.b][0] += dx * f * 0.5; disp[e.b][1] += dy * f * 0.5; disp[e.b][2] += dz * f * 0.5;
    }
    for (let i = 0; i < n; i++) {
      const d = disp[i];
      // weak anchor to the PCA coordinates keeps the global structure
      for (let k = 0; k < 3; k++) d[k] += (pca[i][k] - pos[i][k]) * 0.08;
      const len = Math.hypot(d[0], d[1], d[2]) || 1;
      const step = Math.min(len, temp);
      for (let k = 0; k < 3; k++) pos[i][k] += (d[k] / len) * step;
    }
  }
  return { positions: pos, pca, explained };
}

/** Species that appear in more than one pathway: the crosstalk hubs. */
export function crosstalkSpecies(graph) {
  return graph.nodes.filter((n) => n.kind === 'species' && n.pathways.size > 1).sort((a, b) => b.pathways.size - a.pathways.size || b.degree - a.degree);
}
