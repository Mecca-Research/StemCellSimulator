// Reaction network of everything the simulator models, grouped by the paper's
// four orders and embedded in 3D by PCA of the adjacency matrix.
import * as THREE from 'three';
import { cellMaterial, lineMaterial } from '../engine/materials.js';
import { unitSphere } from '../engine/geometry.js';
import { InstancePool } from '../engine/instances.js';
import { allReactions, countsByOrder, buildGraph, embed, crosstalkSpecies, PAPER_CUMULATIVE } from '../models/network/registry.js';

const ORDER_COLORS = { 1: 0x5ee0c1, 2: 0xffb454, 3: 0xff6fae, 4: 0x8ea8ff };
const ORDER_NAMES = { 1: 'Order 1 · stem-cell signalling', 2: 'Order 2 · progenitor differentiation', 3: 'Order 3 · differentiated cells', 4: 'Order 4 · integration & tissue' };

export default {
  about: `
    <p>Every reaction implemented in this simulator's models, as one network:
    <b>species</b> (spheres, size ∝ connections) and <b>reactions</b> (octahedra,
    coloured by the paper's order). Species shared between pathways — the
    <b>crosstalk hubs</b> — carry a halo.</p>
    <p>Layout: the first three <b>principal components of the adjacency matrix</b>
    (the paper's "dimensionality reduction: PCA" step), refined by a short spring
    relaxation. Reveal the orders one after another to watch the network grow from
    core stem-cell signalling to tissue-level integration.</p>
    <p class="note">The paper reports cumulative reaction counts of ~500, ~1,051,
    ~2,051 and ~3,021 for orders 1–4 (generated in batches). The counts here are the
    mechanistic reactions explicitly modelled in this code base; they are shown side
    by side, not conflated.</p>`,
  paperRef: 'Paper: "Methodology" (graph construction, batch combination, PCA, visualisation) and "Order of Chemical Reactions" (orders 1–4).',

  async create(ctx) {
    const { stage, root, panel, legend, setStatus } = ctx;
    const reactions = allReactions();
    const graph = buildGraph(reactions);
    const { positions, explained } = embed(graph, { refine: 220 });
    const hubs = new Set(crosstalkSpecies(graph).map((n) => n.id));
    const params = { upTo: 4, reveal: 1, labels: true, spin: true };

    const group = new THREE.Group();
    root.add(group);
    const speciesPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'molecule', color: 0xffffff }), graph.nodes.length, group);
    const rxnPool = new InstancePool(new THREE.OctahedronGeometry(1, 0), cellMaterial({ role: 'reporter', color: 0xffffff }), graph.nodes.length, group);
    const haloPool = new InstancePool(new THREE.TorusGeometry(1, 0.08, 6, 24), cellMaterial({ role: 'reporter', color: 0xffe066 }), graph.nodes.length, group);

    // edges as one LineSegments with per-vertex colours
    const edgePos = new Float32Array(graph.edges.length * 6);
    const edgeCol = new Float32Array(graph.edges.length * 6);
    const eg = new THREE.BufferGeometry();
    eg.setAttribute('position', new THREE.BufferAttribute(edgePos, 3));
    eg.setAttribute('color', new THREE.BufferAttribute(edgeCol, 3));
    const lines = new THREE.LineSegments(eg, lineMaterial({ vertexColors: true, opacity: 0.55 }));
    group.add(lines);

    const tmp = new THREE.Color();
    const white = new THREE.Color(0xffffff);
    const mRot = new THREE.Matrix4(), mScale = new THREE.Matrix4();
    const nodeVisible = (i) => graph.nodes[i].order <= params.upTo;
    function draw(t) {
      speciesPool.begin(); rxnPool.begin(); haloPool.begin();
      graph.nodes.forEach((n, i) => {
        if (!nodeVisible(i)) return;
        const [x, y, z] = positions[i];
        const col = tmp.set(ORDER_COLORS[n.order]);
        if (n.kind === 'species') {
          const r = 0.9 + Math.sqrt(n.degree) * 0.55;
          speciesPool.put(x, y, z, r, r, r, tmp.copy(col).lerp(white, 0.55));
          if (hubs.has(n.id)) {
            mRot.makeRotationY(t * 0.8 + i).premultiply(mScale.makeScale(r * 1.8, r * 1.8, r * 1.8)).setPosition(x, y, z);
            haloPool.putMatrix(mRot, tmp.set(0xffe066));
          }
        } else {
          rxnPool.put(x, y, z, 1.3, 1.3, 1.3, col);
        }
      });
      speciesPool.end(); rxnPool.end(); haloPool.end();
      let k = 0;
      const c1 = new THREE.Color(), c2 = new THREE.Color();
      for (const e of graph.edges) {
        const vis = nodeVisible(e.a) && nodeVisible(e.b);
        const a = positions[e.a], b = positions[e.b];
        const ord = Math.max(graph.nodes[e.a].order, graph.nodes[e.b].order);
        c1.set(ORDER_COLORS[ord]).multiplyScalar(e.role === 'modifier' ? 0.45 : 0.9);
        c2.copy(c1).multiplyScalar(e.role === 'product' ? 1.2 : 0.6);
        edgePos.set(vis ? a : [0, 0, 0], k * 6); edgePos.set(vis ? b : [0, 0, 0], k * 6 + 3);
        edgeCol.set([c1.r, c1.g, c1.b, c2.r, c2.g, c2.b], k * 6);
        k++;
      }
      eg.attributes.position.needsUpdate = true;
      eg.attributes.color.needsUpdate = true;
      eg.computeBoundingSphere();
    }

    // ---- panel
    const counts = countsByOrder(reactions);
    panel.section('Orders');
    panel.slider({ label: 'Reveal orders up to', min: 1, max: 4, step: 1, value: 4, format: (v) => `order ${v}`, onChange: (v) => { params.upTo = v; draw(0); } });
    panel.toggle({ label: 'Rotate', value: true, onChange: (v) => { params.spin = v; } });
    let cum = 0;
    panel.table(['order', 'modelled here', 'cumulative', 'paper (cumulative)'], [1, 2, 3, 4].map((o) => { cum += counts[o]; return [ORDER_NAMES[o], counts[o], cum, `~${PAPER_CUMULATIVE[o].toLocaleString()}`]; }));
    const species = graph.nodes.filter((n) => n.kind === 'species').length;
    const ro = panel.readouts(['network']);
    ro.set('network', `${reactions.length} reactions · ${species} species · ${graph.edges.length} edges`);
    panel.note(`PCA of the diffused adjacency matrix explains ${(explained.reduce((a, b) => a + b, 0) * 100).toFixed(1)} % of the variance in the first three components (${explained.map((e) => (e * 100).toFixed(1)).join(' / ')} %).`);
    const ch = panel.chart({ title: 'Cumulative reactions by order: modelled vs paper', xLabel: 'order', logY: true, series: [{ name: 'modelled (this code)', color: '#5ee0c1', points: true }, { name: 'paper estimate', color: '#ff6fae', points: true, dash: [4, 3] }] });
    let c2 = 0;
    ch.set([1, 2, 3, 4], [[1, 2, 3, 4].map((o) => (c2 += counts[o])), [1, 2, 3, 4].map((o) => PAPER_CUMULATIVE[o])]);
    panel.section('Crosstalk hubs (species in ≥ 2 pathways)');
    panel.table(['species', 'pathways'], crosstalkSpecies(graph).slice(0, 14).map((n) => [n.label, [...n.pathways].join(', ')]));
    panel.section('Pathways');
    const byPath = {};
    for (const r of reactions) byPath[r.pathway] = (byPath[r.pathway] ?? 0) + 1;
    panel.table(['pathway', 'reactions'], Object.entries(byPath).sort((a, b) => b[1] - a[1]));

    legend([
      ...[1, 2, 3, 4].map((o) => ({ color: `#${new THREE.Color(ORDER_COLORS[o]).getHexString()}`, label: ORDER_NAMES[o] })),
      { color: '#ffe066', label: 'crosstalk hub (shared species)' },
    ]);

    stage.setPickables([
      { get object() { return speciesPool.mesh; }, info: (hit) => nodeInfo('species', hit.instanceId) },
      { get object() { return rxnPool.mesh; }, info: (hit) => nodeInfo('reaction', hit.instanceId) },
    ]);
    function nodeInfo(kind, inst) {
      let k = -1;
      for (let i = 0; i < graph.nodes.length; i++) {
        if (!nodeVisible(i) || graph.nodes[i].kind !== kind) continue;
        if (++k === inst) {
          const n = graph.nodes[i];
          if (kind === 'reaction') {
            const r = n.reaction;
            return `${r.label}  [order ${r.order} · ${r.pathway}]\n${r.reactants.join(' + ') || '∅'} → ${r.products.join(' + ') || '∅'}${r.modifiers.length ? `\nmodifiers: ${r.modifiers.join(', ')}` : ''}`;
          }
          return `${n.label}\n${n.degree} connections · ${[...n.pathways].join(', ')}`;
        }
      }
      return null;
    }

    stage.frame([0, 0, 0], 105, new THREE.Vector3(0.6, 0.5, 1));
    stage.clipAxis.set(0, 0, 1);
    draw(0);
    let t = 0;
    return {
      update(dt) {
        t += dt;
        if (params.spin) group.rotation.y += dt * 0.08;
        if (Math.floor(t * 4) !== Math.floor((t - dt) * 4)) draw(t);
        setStatus(`${reactions.length} modelled reactions · orders ≤ ${params.upTo}`);
      },
    };
  },
};
