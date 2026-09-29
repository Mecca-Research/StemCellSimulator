// Protein foundry: the paper's hemoglobin and collagen pathways as kinetics
// (population view) plus a hero molecule assembling step by step in 3D.
import * as THREE from 'three';
import { cellMaterial, releaseTree } from '../engine/materials.js';
import { Panel } from '../engine/ui.js';
import { unitSphere, foldedChainGeometry, blobGeometry, DynamicTube } from '../engine/geometry.js';
import { InstancePool } from '../engine/instances.js';
import { rk45 } from '../models/core/ode.js';
import { RNG } from '../models/core/rng.js';
import { hemoglobin, collagen } from '../models/protein/assembly.js';
import { fmt } from '../engine/chart.js';

const smooth = (a, b, x) => { const t = Math.min(Math.max((x - a) / (b - a), 0), 1); return t * t * (3 - 2 * t); };

export default {
  about: `
    <p>The paper's protein-construction pathways — <b>hemoglobin</b> (HBA1/HBA2 + HBB →
    α/β globins → αβ dimers → α₂β₂ tetramer → + 4 heme) and <b>type I collagen</b>
    (COL1A1/COL1A2 → pro-α chains → hydroxylation & glycosylation → 2 α1 + 1 α2 triple
    helix → secretion → propeptide cleavage → tropocollagen → fibril) — as mass-action
    kinetics with exact mass balance.</p>
    <p>Left: the molecular population inside a cell, with particle counts set by the
    kinetics. Right: one molecule assembling in 3D. Knobs expose real biology:
    β-globin output (β-thalassaemia), the AHSP chaperone, HRI coupling of translation to
    heme, vitamin C for prolyl-4-hydroxylase (scurvy), and ADAMTS2 (dermatosparaxis).
    Shapes are coarse-grained, not atomic structures.</p>`,
  paperRef: 'Paper: "Constructing Complex Globular Proteins" → "Example: Hemoglobin Synthesis Pathway" and "Collagen" (transcription, translation, post-translational modification, triple helix, secretion, fibril formation).',

  async create(ctx) {
    const { stage, root, panel, legend, setStatus, query } = ctx;
    const rng = new RNG(33);
    const params = { view: query?.get('view') === 'collagen' ? 'collagen' : 'hemoglobin', hero: true };
    let view = null;

    panel.section('Protein');
    panel.select({
      label: 'Pathway', value: params.view,
      options: [{ value: 'hemoglobin', label: 'Hemoglobin (globular, α₂β₂ + 4 heme)' }, { value: 'collagen', label: 'Collagen type I (triple helix → fibril)' }],
      onChange: (v) => { params.view = v; build(); },
    });
    const host = document.createElement('div');
    panel.current.append(host);
    const sub = new Panel(host);

    function build() {
      view?.dispose();
      releaseTree(root);
      root.clear();
      sub.clear();
      view = params.view === 'hemoglobin' ? buildHemoglobin() : buildCollagen();
    }

    // ================================================================ hemoglobin
    function buildHemoglobin() {
      const group = new THREE.Group();
      root.add(group);
      const p = { ...hemoglobin.defaults };
      // start from a developed state (15 time units) so the cell is populated
      let t = 15;
      let y = rk45(hemoglobin.rhs(p), hemoglobin.initial(), 0, t, { rtol: 1e-6, atol: 1e-8, h0: 0.01 }).y;
      const I = hemoglobin.idx;

      // cell interior backdrop
      const cell = new THREE.Group();
      cell.position.set(-70, 0, 0);
      group.add(cell);
      const nucleus = new THREE.Mesh(unitSphere(4), cellMaterial({ role: 'nucleus', color: 0x8fa8ff, noiseScale: 0.08 }));
      nucleus.scale.setScalar(26);
      nucleus.position.set(-36, 0, 0);
      cell.add(nucleus);
      const membrane = new THREE.Mesh(unitSphere(4), cellMaterial({ role: 'membrane', color: 0xffb3c1, opacity: 0.08 }));
      membrane.scale.set(80, 58, 62);
      cell.add(membrane);

      const alphaGeo = foldedChainGeometry(() => rng.next(), { length: 42, radius: 1.8, tube: 0.32 });
      const betaGeo = foldedChainGeometry(() => rng.next(), { length: 44, radius: 1.8, tube: 0.32 });
      const tetGeo = blobGeometry(3.2, 0.22, 3, 3);
      const hemeGeo = new THREE.CylinderGeometry(0.9, 0.9, 0.2, 6);
      const pools = {
        mRNA: new InstancePool(new THREE.CapsuleGeometry(0.35, 6, 3, 6), cellMaterial({ role: 'reporter', color: 0xffe066 }), 256, cell),
        ribo: new InstancePool(blobGeometry(1.4, 0.25, 7, 2), cellMaterial({ role: 'molecule', color: 0xc6b3ff }), 512, cell),
        alpha: new InstancePool(alphaGeo, cellMaterial({ role: 'molecule', color: 0xff5a6e }), 256, cell),
        beta: new InstancePool(betaGeo, cellMaterial({ role: 'molecule', color: 0x5aa8ff }), 256, cell),
        dimer: new InstancePool(blobGeometry(2.4, 0.3, 5, 2), cellMaterial({ role: 'molecule', color: 0xc77dff }), 256, cell),
        apo: new InstancePool(tetGeo, cellMaterial({ role: 'molecule', color: 0x9aa3b5 }), 256, cell),
        hb: new InstancePool(tetGeo, cellMaterial({ role: 'molecule', color: 0xe0303f, emissive: 0.1 }), 512, cell),
        heme: new InstancePool(hemeGeo, cellMaterial({ role: 'molecule', color: 0x8a1c1c }), 512, cell),
        precip: new InstancePool(blobGeometry(1.6, 0.45, 11, 2), cellMaterial({ role: 'solid', color: 0x3d2430 }), 256, cell),
      };
      // stable random positions per pool slot (cytoplasm shell outside the nucleus)
      const slots = Array.from({ length: 600 }, () => {
        let v;
        do { v = [rng.uniform(-70, 70), rng.uniform(-48, 48), rng.uniform(-52, 52)]; }
        while ((v[0] / 76) ** 2 + (v[1] / 54) ** 2 + (v[2] / 58) ** 2 > 1 || Math.hypot(v[0] + 36, v[1], v[2]) < 30);
        return { p: v, q: new THREE.Quaternion().setFromEuler(new THREE.Euler(rng.uniform(0, 6), rng.uniform(0, 6), rng.uniform(0, 6))), w: rng.normal() };
      });
      const m4 = new THREE.Matrix4(), sc = new THREE.Vector3(), pv = new THREE.Vector3(), qq = new THREE.Quaternion();
      function drawPopulation(time) {
        const scale = { mRNA: 1.5, ribo: 4, alpha: 3, beta: 3, dimer: 3, apo: 0.25, hb: 0.25, heme: 1.5, precip: 0.08 };
        const counts = {
          mRNA: y[I.ma] + y[I.mb], ribo: (y[I.ma] + y[I.mb]), alpha: y[I.a] + y[I.aA], beta: y[I.b], dimer: y[I.D],
          apo: y[I.T0] + y[I.T1] + y[I.T2] + y[I.T3], hb: y[I.T4], heme: y[I.heme], precip: y[I.P],
        };
        let slot = 0;
        for (const [k, pool] of Object.entries(pools)) {
          const n = Math.min(Math.round(counts[k] * scale[k]), pool.capacity, 180);
          pool.begin();
          for (let i = 0; i < n; i++) {
            const s = slots[(slot + i * 7) % slots.length];
            const wob = Math.sin(time * 0.7 + s.w * 5 + i) * 1.2;
            pv.set(s.p[0] + wob, s.p[1] + Math.cos(time * 0.5 + i) * 1.2, s.p[2]);
            qq.copy(s.q);
            sc.setScalar(k === 'ribo' ? 1 : 1);
            m4.compose(pv, qq, sc);
            pool.putMatrix(m4);
          }
          pool.end();
          slot += 53;
        }
      }

      // ---- hero hemoglobin assembly (right)
      const hero = new THREE.Group();
      hero.position.set(75, 0, 0);
      group.add(hero);
      const heroScale = 7;
      const chainMat = { a: cellMaterial({ role: 'molecule', color: 0xff5a6e }), b: cellMaterial({ role: 'molecule', color: 0x5aa8ff }) };
      const subunits = [];
      const tetPos = [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]].map((v) => new THREE.Vector3(...v).multiplyScalar(2.1));
      const kinds = ['a', 'b', 'a', 'b'];
      for (let k = 0; k < 4; k++) {
        const g = foldedChainGeometry(() => rng.next(), { length: 46, radius: 2.0, tube: 0.36 });
        const m = new THREE.Mesh(g, chainMat[kinds[k]]);
        m.scale.setScalar(heroScale);
        hero.add(m);
        subunits.push({ mesh: m, geo: g, kind: kinds[k], start: new THREE.Vector3((k - 1.5) * 20, 44, (k % 2) * 14 - 7), rot: new THREE.Euler(rng.uniform(0, 6), rng.uniform(0, 6), 0) });
      }
      const hemes = [], o2s = [];
      for (let k = 0; k < 4; k++) {
        const h = new THREE.Group();
        h.add(new THREE.Mesh(hemeGeo, cellMaterial({ role: 'molecule', color: 0x8a1c1c })));
        const fe = new THREE.Mesh(unitSphere(2), cellMaterial({ role: 'reporter', color: 0xff9f43 }));
        fe.scale.setScalar(0.35);
        h.add(fe);
        h.scale.setScalar(heroScale * 0.9);
        hero.add(h);
        hemes.push(h);
        const o2 = new THREE.Group();
        for (const dx of [-0.35, 0.35]) { const at = new THREE.Mesh(unitSphere(2), cellMaterial({ role: 'reporter', color: 0xff3030 })); at.scale.setScalar(0.42); at.position.x = dx; o2.add(at); }
        o2.scale.setScalar(heroScale);
        hero.add(o2);
        o2s.push(o2);
      }
      const stageNames = ['translation (ribosome makes α and β chains)', 'folding (chaperones)', 'αβ dimer assembly', 'α₂β₂ tetramer assembly', 'heme incorporation (×4)', 'functional hemoglobin binds O₂'];
      function drawHero(time) {
        const T = time % 16;
        const idxStage = T < 2.5 ? 0 : T < 4 ? 1 : T < 6.5 ? 2 : T < 9 ? 3 : T < 12 ? 4 : 5;
        subunits.forEach((s, k) => {
          const g = s.geo;
          const grow = smooth(0, 2.5, T - k * 0.25);
          g.setDrawRange(0, Math.round(grow * g.index.count / 3) * 3);
          const dimerPartner = Math.floor(k / 2); // (0,1) and (2,3)
          const dimerCenter = new THREE.Vector3(dimerPartner ? 22 : -22, 10, 0);
          const inDimer = dimerCenter.clone().add(tetPos[k].clone().multiplyScalar(heroScale * 0.9).setY(tetPos[k].y * heroScale * 0.9 * 0.5));
          const inTet = tetPos[k].clone().multiplyScalar(heroScale * 0.95);
          const a = smooth(4, 6.5, T), b = smooth(6.5, 9, T);
          const pos = s.start.clone().lerp(inDimer, a).lerp(inTet, b);
          s.mesh.position.copy(pos);
          const spin = (1 - smooth(2.5, 4, T)) * 2.5;
          s.mesh.rotation.set(s.rot.x + spin, s.rot.y + spin * 0.6, 0);
        });
        hemes.forEach((h, k) => {
          const f = smooth(9 + k * 0.5, 10.5 + k * 0.5, T);
          const target = tetPos[k].clone().multiplyScalar(heroScale * 1.05);
          const from = target.clone().normalize().multiplyScalar(90);
          h.position.copy(from.lerp(target, f));
          h.visible = T > 8.8;
          h.lookAt(0, 0, 0);
          h.rotateX(Math.PI / 2);
        });
        o2s.forEach((o, k) => {
          const f = smooth(12.5 + k * 0.4, 13.5 + k * 0.4, T);
          const target = tetPos[k].clone().multiplyScalar(heroScale * 1.45);
          o.position.copy(target.clone().normalize().multiplyScalar(100).lerp(target, f));
          o.visible = T > 12.3;
        });
        hero.rotation.y = time * 0.15;
        return stageNames[idxStage];
      }

      // ---- controls
      sub.section('Kinetics');
      sub.slider({ label: 'α-globin gene output (HBA1/HBA2)', min: 0, max: 2, value: 1, onChange: (v) => { p.alphaExpr = v; } });
      sub.slider({ label: 'β-globin gene output (HBB) — lower = β-thalassaemia', min: 0, max: 2, value: 1, onChange: (v) => { p.betaExpr = v; } });
      sub.slider({ label: 'Heme synthesis (ALAS2 → ferrochelatase)', min: 0, max: 100, value: p.hemeSyn, onChange: (v) => { p.hemeSyn = v; } });
      sub.slider({ label: 'AHSP chaperone pool', min: 0, max: 60, value: p.AHSP, onChange: (v) => { p.AHSP = v; } });
      sub.toggle({ label: 'HRI couples globin translation to heme', value: true, onChange: (v) => { p.hri = v ? 1 : 0; } });
      sub.buttons([{ label: 'Restart kinetics', onClick: () => { y = hemoglobin.initial(); t = 0; chart.clear(); } }]);
      sub.equation('α + β → αβ          (k_d)\nαβ + αβ → α₂β₂       (k_t)\nα₂β₂ + 4 heme → Hb   (k_h, stepwise)\nα + AHSP ⇌ α·AHSP ;  α → precipitate\ntranslation × heme/(K + heme)   (HRI)');
      const chart = sub.chart({ title: 'Species (arbitrary units)', xLabel: 't', series: [
        { name: 'free α', color: '#ff5a6e' }, { name: 'free β', color: '#5aa8ff' }, { name: 'apo/partial', color: '#9aa3b5' },
        { name: 'Hb (4 heme)', color: '#ffb454' }, { name: 'α precipitate', color: '#b07a8f', dash: [4, 3] }] });
      const out = sub.readouts(['functional Hb', 'α made / precipitated', 'α : β chains in Hb']);

      legend([
        { color: '#ffe066', label: 'globin mRNA' }, { color: '#c6b3ff', label: 'ribosomes' },
        { color: '#ff5a6e', label: 'α-globin' }, { color: '#5aa8ff', label: 'β-globin' }, { color: '#c77dff', label: 'αβ dimer' },
        { color: '#9aa3b5', label: 'tetramer lacking heme' }, { color: '#e0303f', label: 'hemoglobin (4 heme)' },
        { color: '#8a1c1c', label: 'free heme' }, { color: '#3d2430', label: 'α precipitate (Heinz body-like)' },
      ]);
      stage.frame([0, 0, 0], 135, new THREE.Vector3(0.1, 0.35, 1));
      stage.clipAxis.set(0, 0, 1);
      let heroLabel = '', lastPush = -1;
      return {
        update(dt) {
          const tEnd = t + dt * 2;
          y = rk45(hemoglobin.rhs(p), y, t, tEnd, { rtol: 1e-6, atol: 1e-8, h0: 0.01 }).y;
          t = tEnd;
          drawPopulation(t);
          heroLabel = drawHero(t * 0.5 + 0.0);
          if (t - lastPush > 0.5) {
            lastPush = t;
            chart.push(t, [y[I.a] + y[I.aA], y[I.b], y[I.T0] + y[I.T1] + y[I.T2] + y[I.T3], y[I.T4] / 20, y[I.P] / 20]);
            out.set('functional Hb', y[I.T4]);
            out.set('α made / precipitated', `${fmt(y[I.aSyn])} / ${fmt(y[I.P])}`);
            out.set('α : β chains in Hb', '2 : 2');
          }
          setStatus(`t = ${t.toFixed(1)} · hero: ${heroLabel}`);
        },
        dispose() {},
      };
    }

    // ================================================================ collagen
    function buildCollagen() {
      const group = new THREE.Group();
      root.add(group);
      const p = { ...collagen.defaults };
      let t = 15;
      let y = rk45(collagen.rhs(p), collagen.initial(), 0, t, { rtol: 1e-6, atol: 1e-8, h0: 0.01 }).y;
      const I = collagen.idx;

      // --- population: ER -> secretion -> extracellular fibril
      const pop = new THREE.Group();
      pop.position.set(-80, 0, 0);
      group.add(pop);
      const pm = new THREE.Mesh(new THREE.PlaneGeometry(120, 110), cellMaterial({ role: 'membrane', color: 0xffc3d0, opacity: 0.15 }));
      pm.rotation.y = Math.PI / 2;
      pop.add(pm);
      const rodGeo = new THREE.CylinderGeometry(0.6, 0.6, 30, 8, 1);
      rodGeo.rotateZ(Math.PI / 2);
      const endGeo = blobGeometry(1.8, 0.3, 9, 2);
      const pools = {
        chain: new InstancePool(new THREE.CylinderGeometry(0.3, 0.3, 30, 6).rotateZ(Math.PI / 2), cellMaterial({ role: 'molecule', color: 0xff9ab0 }), 256, pop),
        proc: new InstancePool(rodGeo, cellMaterial({ role: 'molecule', color: 0xffd166 }), 256, pop),
        ends: new InstancePool(endGeo, cellMaterial({ role: 'molecule', color: 0x7bdff2 }), 512, pop),
        tropo: new InstancePool(rodGeo, cellMaterial({ role: 'molecule', color: 0xff7a59 }), 256, pop),
        fibril: new InstancePool(rodGeo, cellMaterial({ role: 'solid', color: 0xf4a9a0, stain: 0xffffff }), 2048, pop),
      };
      const inside = Array.from({ length: 400 }, () => ({ p: [rng.uniform(-55, -6), rng.uniform(-45, 45), rng.uniform(-40, 40)], a: rng.uniform(0, 6) }));
      const outside = Array.from({ length: 400 }, () => ({ p: [rng.uniform(8, 60), rng.uniform(-45, 45), rng.uniform(-40, 40)], a: rng.uniform(0, 6) }));
      const m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), s1 = new THREE.Vector3(1, 1, 1), pv = new THREE.Vector3(), ax = new THREE.Vector3(0, 1, 0);
      function drawPopulation(time) {
        const put = (pool, list, n, extra) => {
          pool.begin();
          for (let i = 0; i < Math.min(n, pool.capacity, list.length); i++) {
            const s = list[i];
            q.setFromAxisAngle(ax, s.a + time * 0.05);
            pv.set(s.p[0], s.p[1] + Math.sin(time + i) * 0.8, s.p[2]);
            m4.compose(pv, q, s1);
            pool.putMatrix(m4);
            extra?.(s, i);
          }
          pool.end();
        };
        const chains = Math.round((y[I.c1] + y[I.c2] + y[I.h1] + y[I.h2]) * 1.2);
        put(pools.chain, inside, chains);
        const pcN = Math.round((y[I.PC]) * 3), pceN = Math.round(y[I.PCe] * 3);
        pools.ends.begin();
        pools.proc.begin();
        for (let i = 0; i < Math.min(pcN + pceN, 150); i++) {
          const s = i < pcN ? inside[(i + 200) % 400] : outside[i % 400];
          q.setFromAxisAngle(ax, s.a); pv.set(...s.p); m4.compose(pv, q, s1); pools.proc.putMatrix(m4);
          const dir = new THREE.Vector3(1, 0, 0).applyQuaternion(q);
          for (const sg of [1, -1]) pools.ends.put(s.p[0] + sg * 16 * dir.x, s.p[1], s.p[2] + sg * 16 * dir.z, 1, 1, 1);
        }
        pools.proc.end(); pools.ends.end();
        put(pools.tropo, outside.slice(150), Math.round(y[I.TC] * 3));
        // fibril: quarter-staggered array of tropocollagen along +x, D = 67 nm (scaled)
        const units = Math.min(Math.round(y[I.F] * 1.5), pools.fibril.capacity);
        pools.fibril.begin();
        const L = 30, gap = 0.54 * (L / 4.46) * 0.6; // 0.54 D gap zone, 4.46 D molecule
        for (let i = 0; i < units; i++) {
          const row = i % 5, col = Math.floor(i / 5);
          const ring = Math.floor(col / 3);
          const ang = (row / 5) * Math.PI * 2 + ring * 0.6;
          const r = 1.3 + ring * 1.3;
          const x = 20 + (col % 3) * (L + gap) + row * (L / 4.46) * 1.0;
          m4.compose(pv.set(x, Math.cos(ang) * r - 20, Math.sin(ang) * r), q.identity(), s1);
          pools.fibril.putMatrix(m4);
        }
        pools.fibril.end();
      }

      // --- hero triple helix
      const hero = new THREE.Group();
      hero.position.set(60, 18, 0);
      hero.scale.setScalar(1.7);
      group.add(hero);
      const Lh = 60, Rh = 1.2, turns = 5;
      const strandMats = [cellMaterial({ role: 'molecule', color: 0xff5a8a }), cellMaterial({ role: 'molecule', color: 0xff5a8a }), cellMaterial({ role: 'molecule', color: 0x5aa8ff })];
      const tubes = strandMats.map((mat) => { const dt = new DynamicTube(160, 8, 0.55); hero.add(new THREE.Mesh(dt.geometry, mat)); return dt; });
      const hyp = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter', color: 0xfff275 }), 200, hero);
      const sugar = new InstancePool(new THREE.CylinderGeometry(0.5, 0.5, 0.2, 6), cellMaterial({ role: 'reporter', color: 0x7bff9e }), 60, hero);
      const propep = [0, 1].map(() => { const m = new THREE.Mesh(blobGeometry(3.4, 0.35, 5, 3), cellMaterial({ role: 'molecule', color: 0x7bdff2 })); hero.add(m); return m; });
      const tmp = new THREE.Vector3();
      const stageNames = ['translation of pro-α1, pro-α1, pro-α2 into the ER', 'prolyl hydroxylation + glycosylation', 'triple helix zips from the C-terminus', 'secretion + propeptide cleavage (ADAMTS2, BMP1)', 'tropocollagen joins a D-banded fibril'];
      const strandPoint = (k, s, T, out) => {
        const x = s * Lh - Lh / 2;
        const sep = new THREE.Vector3(0, (k - 1) * 5, (k === 1 ? 3 : -2));
        const ppII = 0.35 * Math.sin(s * 60 + k);
        const zip = smooth(0, 1, (T - 5) / 3 * 1.2 - (1 - s)); // front moves from C (s=1) to N (s=0)
        const a = s * turns * Math.PI * 2 + (k * 2 * Math.PI) / 3;
        const hx = out.set(x, Math.cos(a) * Rh, Math.sin(a) * Rh);
        const px = tmp.set(x, sep.y + ppII, sep.z);
        return out.copy(px.lerp(hx, zip));
      };
      function drawHero(time) {
        const T = time % 16;
        const hasAsc = p.ascorbate > 0.05;
        const label = T < 3 ? 0 : T < 5 ? 1 : T < 9 ? 2 : T < 12 ? 3 : 4;
        tubes.forEach((tb, k) => {
          tb.update((s, out) => strandPoint(k, s, hasAsc ? T : Math.min(T, 5.2), out));
          tb.reveal(smooth(0, 3, T - k * 0.3));
        });
        hyp.begin(); sugar.begin();
        const ohOn = smooth(3, 4.5, T) * (hasAsc ? 1 : 0);
        if (ohOn > 0) {
          for (let k = 0; k < 3; k++) for (let i = 2; i < 58; i += 3) {
            const s = i / 60; strandPoint(k, s, T, tmp);
            const r = 0.35 * ohOn; hyp.put(tmp.x, tmp.y + 0.6, tmp.z, r, r, r);
            if (k < 2 && i % 21 === 2) sugar.put(tmp.x, tmp.y + 1.2, tmp.z, ohOn, ohOn, ohOn);
          }
        }
        hyp.end(); sugar.end();
        const cut = smooth(9.5, 11.5, T);
        propep[0].position.set(-Lh / 2 - 3 - cut * 25, cut * 10, 0);
        propep[1].position.set(Lh / 2 + 3 + cut * 25, -cut * 10, 0);
        propep.forEach((m) => { m.visible = T > 1 && T < 12.5; });
        const join = smooth(12, 15, T);
        hero.position.set(60 - join * 20, 18 - join * 36, 0);
        hero.rotation.x = time * 0.25 * (1 - join);
        return hasAsc ? stageNames[label] : 'no ascorbate: prolines stay unhydroxylated, the helix cannot zip at 37 °C';
      }

      sub.section('Kinetics');
      sub.slider({ label: 'Ascorbate (vitamin C) for prolyl-4-hydroxylase', min: 0, max: 2, value: p.ascorbate, onChange: (v) => { p.ascorbate = v; } });
      sub.slider({ label: 'COL1A2 expression (α2 supply)', min: 0, max: 2, value: 1, onChange: (v) => { p.col1a2 = v; } });
      sub.slider({ label: 'ADAMTS2 N-proteinase activity', min: 0, max: 2, value: 1, onChange: (v) => { p.ADAMTS2 = v; } });
      sub.buttons([{ label: 'Restart kinetics', onClick: () => { y = collagen.initial(); t = 0; chart.clear(); } }]);
      sub.equation('pro-α → pro-α(OH)       k_OH · asc/(K + asc)\n2 α1(OH) + α2(OH) → procollagen\nprocollagen → secreted → tropocollagen   (ADAMTS2, BMP1)\n2 TC → nucleus ;  TC + fibril → fibril');
      const chart = sub.chart({ title: 'Species (arbitrary units)', xLabel: 't', series: [
        { name: 'chains', color: '#ff9ab0' }, { name: 'procollagen', color: '#ffd166' }, { name: 'tropocollagen', color: '#ff7a59' },
        { name: 'fibril ÷ 10', color: '#5ee0c1' }, { name: 'degraded ÷ 10', color: '#8a98ab', dash: [4, 3] }] });
      legend([
        { color: '#ff5a8a', label: 'pro-α1(I) chains (×2)' }, { color: '#5aa8ff', label: 'pro-α2(I) chain' },
        { color: '#fff275', label: '4-hydroxyproline' }, { color: '#7bff9e', label: 'glycosylated hydroxylysine' },
        { color: '#7bdff2', label: 'N/C propeptides' }, { color: '#ffd166', label: 'procollagen' },
        { color: '#ff7a59', label: 'tropocollagen' }, { color: '#f4a9a0', label: 'fibril (quarter-staggered)' },
      ]);
      stage.frame([0, 0, 0], 130, new THREE.Vector3(0.05, 0.4, 1));
      let lastPush = -1, label = '';
      return {
        update(dt) {
          const tEnd = t + dt * 2;
          y = rk45(collagen.rhs(p), y, t, tEnd, { rtol: 1e-6, atol: 1e-8, h0: 0.01 }).y;
          t = tEnd;
          drawPopulation(t);
          label = drawHero(t * 0.5);
          if (t - lastPush > 0.5) {
            lastPush = t;
            chart.push(t, [y[I.c1] + y[I.c2] + y[I.h1] + y[I.h2], y[I.PC] + y[I.PCe], y[I.TC], y[I.F] / 10, y[I.deg] / 10]);
          }
          setStatus(`t = ${t.toFixed(1)} · hero: ${label}`);
        },
        dispose() {},
      };
    }

    build();
    return {
      update(dt) { view?.update(dt); sub.drawCharts(); },
      dispose() { view?.dispose(); },
    };
  },
};
