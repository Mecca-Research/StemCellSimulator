// The paper's cellular map in 3D: lineage tree from the pluripotent stem cell
// to 150+ specialised cell types, the cells -> tissues -> organs -> systems
// hierarchy, and the organ-system interaction circle.
import * as THREE from 'three';
import { cellMaterial, lineMaterial, releaseTree } from '../engine/materials.js';
import { unitSphere, erythrocyteGeometry, lobedNucleus, mergeGeometries } from '../engine/geometry.js';
import { InstancePool } from '../engine/instances.js';
import { RNG } from '../models/core/rng.js';
import { LINEAGES, GERM_LAYERS, TISSUES, ORGANS, SYSTEMS, SYSTEM_LINKS, allCellTypes, morphology } from '../models/atlas/atlas.js';

function starGeometry() {
  const parts = [new THREE.IcosahedronGeometry(1, 2).toNonIndexed()];
  const dirs = [[1, 0.2, 0], [-0.8, 0.5, 0.3], [0.1, -1, 0.4], [0.2, 0.3, -1], [-0.3, -0.4, 1], [0.6, 0.9, 0.5]];
  for (const d of dirs) {
    const v = new THREE.Vector3(...d).normalize();
    const g = new THREE.CylinderGeometry(0.08, 0.28, 2.2, 6).toNonIndexed();
    g.translate(0, 1.6, 0);
    g.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), v));
    parts.push(g);
  }
  const m = mergeGeometries(parts);
  m.computeVertexNormals();
  return m;
}

function spermGeometry() {
  const head = new THREE.SphereGeometry(1, 16, 12).toNonIndexed();
  head.scale(0.7, 1, 0.45);
  const tail = new THREE.CylinderGeometry(0.06, 0.12, 5, 6).toNonIndexed();
  tail.translate(0, -3.2, 0);
  const g = mergeGeometries([head, tail]);
  g.computeVertexNormals();
  return g;
}

const GLYPHS = {
  round: () => unitSphere(2),
  big: () => unitSphere(3),
  rbc: () => erythrocyteGeometry(2.2, 32),
  disc: () => new THREE.CylinderGeometry(1, 1, 0.3, 16),
  lobed: () => lobedNucleus(3, 0.6, () => 0.5),
  star: starGeometry,
  capsule: () => new THREE.CapsuleGeometry(0.55, 2.4, 4, 10).rotateZ(Math.PI / 2),
  spindle: () => unitSphere(2).scale(1.9, 0.55, 0.55),
  flat: () => new THREE.CylinderGeometry(1.2, 1.2, 0.3, 18),
  columnar: () => new THREE.CapsuleGeometry(0.55, 1.3, 4, 10),
  sperm: spermGeometry,
};

export default {
  about: `
    <p>The paper's <b>cellular map of the human body</b> in 3D. From the pluripotent
    stem cell at the centre, lineages branch through the germ layers to more than 150
    specialised cell types, each drawn with a glyph of its characteristic shape
    (biconcave erythrocytes, star-shaped neurons and osteocytes, spindle fibroblasts,
    columnar epithelia, lobed granulocytes, fibres).</p>
    <p>Other views: the paper's hierarchy <b>cells → tissues → organs → organ systems</b>,
    and the <b>higher-level interactions between organ systems</b> as a circular graph
    with signals flowing along each link. Hover anything for its name and function.</p>`,
  paperRef: 'Paper: "Creating a Cellular Map of the Human Body", "Comprehensive and Detailed List of Specialized Cells", "Next Order of Complex Cellular Structures", "Organizing in a Tree Structure" and "Higher-Level Interactions".',

  async create(ctx) {
    const { stage, root, panel, legend, setStatus, query } = ctx;
    const views = { lineage: 'Lineage tree: stem cell → 150+ cell types', hierarchy: 'Cells → tissues → organs → organ systems', systems: 'Organ-system interactions' };
    const params = { view: views[query?.get('view')] ? query.get('view') : 'lineage', spin: true };
    let group = null, pickables = [], animate = null;

    panel.section('View');
    panel.select({ label: 'Map', value: params.view, options: Object.entries(views).map(([value, label]) => ({ value, label })), onChange: (v) => { params.view = v; build(); } });
    panel.toggle({ label: 'Rotate', value: true, onChange: (v) => { params.spin = v; } });
    const info = document.createElement('div');
    panel.current.append(info);

    const layerColor = Object.fromEntries(GERM_LAYERS.map((g) => [g.id, new THREE.Color(g.color)]));
    const lineageColor = (l, k, n) => layerColor[l.layer].clone().offsetHSL(((k / Math.max(n, 1)) - 0.5) * 0.12, 0, 0);

    function curve(a, b, lift = 0) {
      const mid = a.clone().add(b).multiplyScalar(0.5);
      mid.add(mid.clone().normalize().multiplyScalar(lift));
      return new THREE.QuadraticBezierCurve3(a, mid, b).getPoints(12);
    }
    function linesFrom(segments, color, opacity = 0.5) {
      const pos = [];
      for (const pts of segments) for (let i = 0; i < pts.length - 1; i++) pos.push(pts[i].x, pts[i].y, pts[i].z, pts[i + 1].x, pts[i + 1].y, pts[i + 1].z);
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
      return new THREE.LineSegments(g, lineMaterial({ color, opacity }));
    }

    // ------------------------------------------------------ lineage tree
    function buildLineage() {
      const g = new THREE.Group();
      const cells = allCellTypes();
      const pools = Object.fromEntries(Object.entries(GLYPHS).map(([k, f]) => [k, new InstancePool(f(), cellMaterial({ role: 'molecule', color: 0xffffff, noise: 0.2 }), 256, g)]));
      const hubPool = new InstancePool(unitSphere(3), cellMaterial({ role: 'membrane', color: 0xffffff, opacity: 0.35 }), 64, g);
      const hubNuc = new InstancePool(unitSphere(3), cellMaterial({ role: 'nucleus', color: 0xb8c8ff }), 64, g);
      const segs = [];
      const center = new THREE.Vector3();
      hubPool.begin(); hubNuc.begin();
      hubPool.put(0, 0, 0, 6, 6, 6, new THREE.Color(0xffffff));
      hubNuc.put(0, 0, 0, 3.6, 3.6, 3.6);
      const nodeInfo = [];
      // germ layers on a ring, lineages on a dome, cells on the outer shell
      const L = GERM_LAYERS.length;
      const byLayer = GERM_LAYERS.map((gl) => LINEAGES.filter((l) => l.layer === gl.id));
      const totalCells = cells.length;
      let cellIdx = 0;
      const R1 = 26, R2 = 60, R3 = 100;
      let angle0 = 0;
      GERM_LAYERS.forEach((gl, gi) => {
        const lin = byLayer[gi];
        const nCells = lin.reduce((s, l) => s + l.cells.length, 0);
        const span = (nCells / totalCells) * Math.PI * 2;
        const aMid = angle0 + span / 2;
        const gp = new THREE.Vector3(Math.cos(aMid) * R1, 6, Math.sin(aMid) * R1);
        hubPool.put(gp.x, gp.y, gp.z, 3.6, 3.6, 3.6, layerColor[gl.id]);
        hubNuc.put(gp.x, gp.y, gp.z, 2, 2, 2);
        nodeInfo.push({ kind: 'hub', text: `${gl.name}\n${lin.map((l) => l.name).join('\n')}` });
        segs.push(curve(center, gp, 4));
        let a = angle0;
        lin.forEach((l, li) => {
          const lspan = (l.cells.length / totalCells) * Math.PI * 2;
          const la = a + lspan / 2;
          const lp = new THREE.Vector3(Math.cos(la) * R2, 12, Math.sin(la) * R2);
          const lc = lineageColor(l, li, lin.length);
          hubPool.put(lp.x, lp.y, lp.z, 2.6, 2.6, 2.6, lc);
          hubNuc.put(lp.x, lp.y, lp.z, 1.4, 1.4, 1.4);
          nodeInfo.push({ kind: 'hub', text: `${l.name}\nstem/progenitor: ${l.stem}\n${l.cells.length} cell types` });
          segs.push(curve(gp, lp, 6));
          l.cells.forEach((cell, ci) => {
            const ca = a + ((ci + 0.5) / l.cells.length) * lspan;
            const tier = ci % 3; // stagger into three latitudes so labels do not collide
            const cp = new THREE.Vector3(Math.cos(ca) * (R3 + tier * 9), 8 + tier * 11, Math.sin(ca) * (R3 + tier * 9));
            const shape = morphology(cell.name);
            const s = shape === 'big' ? 2.4 : shape === 'rbc' ? 1.3 : 1.5;
            const m = new THREE.Matrix4().compose(cp, new THREE.Quaternion().setFromEuler(new THREE.Euler(cellIdx * 0.7, ca, 0)), new THREE.Vector3(s, s, s));
            const inst = pools[shape].putMatrix(m, lc);
            pools[shape].mesh.userData.info ??= [];
            pools[shape].mesh.userData.info[inst] = `${cell.name}\n${cell.fn}${cell.note ? `\nnote: ${cell.note}` : ''}\nlineage: ${l.name} (${gl.name})`;
            segs.push(curve(lp, cp, 3));
            cellIdx++;
          });
          a += lspan;
        });
        angle0 += span;
      });
      hubPool.end(); hubNuc.end();
      for (const p of Object.values(pools)) p.end();
      g.add(linesFrom(segs, 0x9fb3cc, 0.35));
      pickables = [
        ...Object.values(pools).map((p) => ({ object: p.mesh, info: (hit) => p.mesh.userData.info?.[hit.instanceId] ?? null })),
        { object: hubPool.mesh, info: (hit) => (hit.instanceId === 0 ? 'Pluripotent stem cell\n(totipotent zygote → pluripotent epiblast / ES / iPS cells)' : nodeInfo[hit.instanceId - 1]?.text ?? null) },
      ];
      legend(GERM_LAYERS.map((gl) => ({ color: `#${new THREE.Color(gl.color).getHexString()}`, label: gl.name })));
      info.replaceChildren();
      const tbl = document.createElement('div');
      tbl.className = 'note';
      tbl.innerHTML = `<b>${cells.length}</b> cell types in <b>${LINEAGES.length}</b> lineages. Glyph shapes: biconcave disc (erythrocyte), star (neurons, osteocytes, glia, podocytes), spindle (fibroblasts, smooth muscle), fibre (muscle), column (epithelia), lobed (granulocytes), large sphere (adipocytes, megakaryocytes, oocytes).`;
      info.append(tbl);
      stage.frame([0, 10, 0], 120, new THREE.Vector3(0.2, 0.9, 1));
      animate = null;
      return g;
    }

    // ------------------------------------------------------ hierarchy rings
    function buildHierarchy() {
      const g = new THREE.Group();
      const cellsList = [...new Set(TISSUES.flatMap((t) => t.cells))];
      const levels = [
        { name: 'Cells', items: cellsList, y: 0, r: 70, color: 0x5ee0c1 },
        { name: 'Tissues', items: TISSUES.map((t) => t.name), y: 32, r: 38, color: 0xffb454 },
        { name: 'Organs', items: ORGANS.map((o) => o.name), y: 64, r: 60, color: 0xff6fae },
        { name: 'Organ systems', items: SYSTEMS.map((s) => s.name), y: 96, r: 80, color: 0x8ea8ff },
      ];
      const pos = {};
      const pools = levels.map((lv, i) => new InstancePool(i === 0 ? unitSphere(2) : i === 1 ? new THREE.BoxGeometry(1.6, 1, 1.6) : i === 2 ? unitSphere(3) : new THREE.IcosahedronGeometry(1, 1), cellMaterial({ role: i === 0 ? 'nucleus' : 'molecule', color: lv.color }), 64, g));
      const infos = levels.map(() => []);
      levels.forEach((lv, li) => {
        pools[li].begin();
        lv.items.forEach((name, k) => {
          const a = (k / lv.items.length) * Math.PI * 2;
          const p = new THREE.Vector3(Math.cos(a) * lv.r, lv.y, Math.sin(a) * lv.r);
          pos[`${li}:${name}`] = p;
          const s = [1.8, 3.6, 4, 4.6][li];
          pools[li].put(p.x, p.y, p.z, s, s, s);
          const extra = li === 1 ? `\n${TISSUES.find((t) => t.name === name).kinds.join(', ')}` : li === 2 ? `\n${ORGANS.find((o) => o.name === name).parts}` : li === 3 ? `\n${SYSTEMS.find((s2) => s2.name === name).organs}` : '';
          infos[li].push(`${lv.name.slice(0, -1)}: ${name}${extra}`);
        });
        pools[li].end();
      });
      const segs = [];
      for (const t of TISSUES) for (const c of t.cells) segs.push(curve(pos[`0:${c}`], pos[`1:${t.name}`], 0));
      for (const o of ORGANS) {
        for (const t of o.tissues) segs.push(curve(pos[`1:${t}`], pos[`2:${o.name}`], 0));
        segs.push(curve(pos[`2:${o.name}`], pos[`3:${o.system}`], 0));
      }
      g.add(linesFrom(segs, 0xc9d6e8, 0.45));
      // level discs
      for (const lv of levels) {
        const ring = new THREE.Mesh(new THREE.RingGeometry(lv.r - 0.4, lv.r + 0.4, 96).rotateX(-Math.PI / 2), new THREE.MeshBasicMaterial({ color: lv.color, transparent: true, opacity: 0.25, side: THREE.DoubleSide }));
        ring.position.y = lv.y;
        g.add(ring);
      }
      pickables = pools.map((p, li) => ({ object: p.mesh, info: (hit) => infos[li][hit.instanceId] }));
      legend(levels.map((lv) => ({ color: `#${new THREE.Color(lv.color).getHexString()}`, label: lv.name })));
      info.replaceChildren();
      const n = document.createElement('div');
      n.className = 'note';
      n.textContent = 'Edges follow the paper: cells group into the four basic tissues, tissues build organs (e.g. skin = epithelial + connective tissue), and organs belong to organ systems.';
      info.append(n);
      stage.frame([0, 48, 0], 110, new THREE.Vector3(0.3, 0.55, 1));
      animate = null;
      return g;
    }

    // ------------------------------------------------------ system circle
    function buildSystems() {
      const g = new THREE.Group();
      const R = 60;
      const pos = {};
      const pool = new InstancePool(new THREE.IcosahedronGeometry(1, 3), cellMaterial({ role: 'molecule', color: 0xffffff }), 16, g);
      pool.begin();
      const rng = new RNG(3);
      SYSTEMS.forEach((s, k) => {
        const a = (k / SYSTEMS.length) * Math.PI * 2;
        const p = new THREE.Vector3(Math.cos(a) * R, 0, Math.sin(a) * R);
        pos[s.name] = p;
        pool.put(p.x, p.y, p.z, 5, 5, 5, new THREE.Color().setHSL(k / SYSTEMS.length, 0.6, 0.6));
      });
      pool.end();
      const arcs = SYSTEM_LINKS.map((l) => {
        const a = pos[l.a], b = pos[l.b];
        const mid = a.clone().add(b).multiplyScalar(0.5).setY(28 + a.distanceTo(b) * 0.2);
        const c = new THREE.QuadraticBezierCurve3(a, mid, b);
        const tube = new THREE.Mesh(new THREE.TubeGeometry(c, 48, 0.6, 8), cellMaterial({ role: 'reporter', color: 0x9fb3cc }));
        g.add(tube);
        return { curve: c, link: l, tube };
      });
      const flow = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter', color: 0xffe066 }), 400, g);
      const particles = arcs.flatMap((arc, i) => Array.from({ length: 16 }, () => ({ arc: i, u: rng.next(), dir: rng.next() < 0.5 ? 1 : -1 })));
      pickables = [
        { object: pool.mesh, info: (hit) => { const s = SYSTEMS[hit.instanceId]; const links = SYSTEM_LINKS.filter((l) => l.a === s.name || l.b === s.name).map((l) => `↔ ${l.a === s.name ? l.b : l.a}`); return `${s.name}\n${s.organs}\n${links.join('\n')}`; } },
        ...arcs.map((a) => ({ object: a.tube, info: () => `${a.link.a} ↔ ${a.link.b}\n${a.link.text}` })),
      ];
      legend(SYSTEMS.map((s, k) => ({ color: `#${new THREE.Color().setHSL(k / SYSTEMS.length, 0.6, 0.6).getHexString()}`, label: s.name })));
      info.replaceChildren();
      const tbl = document.createElement('div');
      tbl.className = 'note';
      tbl.innerHTML = SYSTEM_LINKS.map((l) => `<div>· <b>${l.a.replace(' System', '')} ↔ ${l.b.replace(' System', '')}</b>: ${l.text}</div>`).join('');
      info.append(tbl);
      stage.frame([0, 10, 0], 85, new THREE.Vector3(0, 0.8, 1));
      const tmp = new THREE.Vector3();
      animate = (dt) => {
        flow.begin();
        for (const p of particles) {
          p.u = (p.u + dt * 0.18 * p.dir + 1) % 1;
          arcs[p.arc].curve.getPoint(p.u, tmp);
          flow.put(tmp.x, tmp.y, tmp.z, 0.9, 0.9, 0.9);
        }
        flow.end();
      };
      return g;
    }

    function build() {
      if (group) { root.remove(group); releaseTree(group); }
      group = params.view === 'lineage' ? buildLineage() : params.view === 'hierarchy' ? buildHierarchy() : buildSystems();
      root.add(group);
      stage.setPickables(pickables);
    }
    build();

    return {
      update(dt) {
        if (params.spin && group) group.rotation.y += dt * 0.05;
        animate?.(dt);
        setStatus(`${allCellTypes().length} cell types · ${LINEAGES.length} lineages · ${ORGANS.length} organs · ${SYSTEMS.length} systems`);
      },
    };
  },
};
