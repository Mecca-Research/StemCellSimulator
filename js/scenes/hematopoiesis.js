// Stage 2 - progenitor commitment in a 3D bone-marrow niche.
//
// Every nucleated cell is an agent that lives on the lineage tree of
// js/models/pathways/hematopoiesis.js and integrates its own GATA1-PU.1 toggle
// (Langevin); the CMP -> MEP/GMP decision is read from that toggle. EPO sets the
// JAK2/STAT5 drive (GATA1 input + erythroid survival), SCF sets the MAPK/ERK
// drive (progenitor cycling), GM-CSF biases PU.1. A Gillespie compartment model
// with the same rates runs alongside for the population charts.
import * as THREE from 'three';
import { cellMaterial } from '../engine/materials.js';
import { InstancePool } from '../engine/instances.js';
import { unitSphere, erythrocyteGeometry, lobedNucleus, indentedNucleus, blobGeometry, mergeVertices } from '../engine/geometry.js';
import { fmt } from '../engine/chart.js';
import { RNG } from '../models/core/rng.js';
import * as mapk from '../models/pathways/mapk.js';
import { jakstat, timeCourse, stat5Drive } from '../models/pathways/jakstat.js';
import {
  toggle, LINEAGE, LineageSSA, meanFieldSteadyState, fateDistribution, gmpFate, mepFate, clpFate,
  erkSpeed, apoptosisRate, marrowFraction, SSA_DEFAULTS,
} from '../models/pathways/hematopoiesis.js';

const HOURS_PER_SECOND = 6;      // biological clock: 1 s of animation = 6 h
const TF_HOURS = 6;              // unit of the toggle equations (1/k, TF turnover)
const FLOW_VMAX = 34;            // um/s centre-line speed in the sinusoid (display speed)
const MAX_CELLS = 340;
const K3D = 24;                  // HSC niche slots in the 3D window
const MITOSIS_H = 1.2;
const ENUC_H = 8;                // enucleation (minutes in vivo; slowed for display)
const APOP_H = 4;
const GAIN_G = 0.08, GAIN_P = 0.08;  // signal -> toggle input gains (phenomenological)
const LYMPH_NUC_R_FALLBACK = 4.59;

// ------------------------------------------------------------------ niche geometry (um)
const DOM = { x0: -86, x1: 86, z0: -30, z1: 30, top: 38 };
const VES = { y: 21, z: -1, r: 12, x0: -104, x1: 104 };
const boneY = (x, z) => -30 + 3.2 * Math.sin(0.045 * x + 0.6) * Math.cos(0.06 * z - 0.3)
  + 1.6 * Math.sin(0.12 * x - 0.09 * z + 1.1) + 0.8 * Math.cos(0.21 * x + 0.17 * z);

const TYPE_GROUP = {
  HSC: 'stem', CMP: 'prog', CLP: 'prog', MEP: 'prog', GMP: 'prog', ProEB: 'ery', EB: 'ery', Retic: 'ery', RBC: 'ery',
  MK: 'mk', Neu: 'gran', Eos: 'gran', Baso: 'gran', Mono: 'mono', B: 'lymph', T: 'lymph', NK: 'lymph', Pyr: 'ery',
};

// physical-mode lineage palette (also used by the legend and the charts)
const HEX = {
  HSC: 0x3ee6c4, CMP: 0x86b6ff, CLP: 0x6aa0ff, MEP: 0xc58cff, GMP: 0xffc36b,
  ProEB: 0x6d7fd6, EB: 0xb04a6a, Retic: 0xc9404a, RBC: 0xd8352c,
  MK: 0xe59ad8, PLT: 0xf6b8ec, Neu: 0xffb454, Eos: 0xff7a45, Baso: 0xa77bff, Mono: 0x8fb0d6,
  B: 0x58a8ff, T: 0x4bd0ff, NK: 0x8f9cff,
};
const css = (h) => `#${h.toString(16).padStart(6, '0')}`;
const OPTICAL_HALF = 12;   // half-thickness of the optical section (um)

/** Smooth-shade a displaced, non-indexed geometry: weld coincident vertices, then average normals. */
function smooth(g) {
  for (const k of Object.keys(g.attributes)) if (k !== 'position') g.deleteAttribute(k);
  const m = mergeVertices(g, 1e-4);
  m.computeVertexNormals();
  g.dispose();
  return m;
}

export default {
  about: `
    <p>A slab of <b>bone marrow</b>: trabecular bone (endosteum lined by osteoblasts), a
    <b>sinusoid</b> carrying blood, CXCL12-abundant reticular stromal cells, adipocytes and two
    <b>erythroblastic islands</b> (a central macrophage ringed by maturing erythroblasts; Bessis 1958).
    <b>HSCs</b> sit in niche slots on the endosteum and near the vessel; each daughter stays a stem
    cell with probability <i>a</i>(1 − H/K) (niche-limited self-renewal, Schofield 1978).</p>
    <p>Every cell carries its own <b>GATA1–PU.1 toggle</b> (Huang, Guo, May &amp; Enver 2007),
    integrated with Langevin noise. HSCs sit in the central <i>primed</i> attractor; in CMPs the
    self-activation decays, the attractor vanishes in a pitchfork and the cell falls to the
    GATA1-high (MEP) or PU.1-high (GMP) state. <b>EPO → JAK2/STAT5</b> feeds GATA1 and keeps
    erythroid progenitors alive; <b>GM-CSF</b> feeds PU.1; <b>SCF → c-Kit → Raf/MEK/ERK</b>
    (Huang–Ferrell cascade) sets how fast progenitors cycle.</p>
    <p>Morphology follows maturation: erythroblasts shrink, haemoglobinise and <b>extrude their
    nucleus</b> (the pyrenocyte is eaten by the island macrophage); reticulocytes remodel into
    Evans–Fung biconcave discs and cross into the sinusoid; polyploid <b>megakaryocytes</b> push
    proplatelets through the endothelium and shed platelets into the flow; neutrophil nuclei go
    round → kidney → band → segmented.</p>
    <p class="note">Simplifications: signal → toggle-input, ERK → cycling-speed and STAT5 → survival
    maps are phenomenological; HSC cycling is compressed (dormant HSCs divide every ~5 months);
    motility, blood flow and enucleation run at display speed; the lymphoid TF state is not
    modelled by the toggle; MK fraction and platelet yield are scaled for display.</p>`,
  paperRef: 'Paper: "Stage 2 · Progenitor Cell Differentiation" (MAPK/ERK, JAK/STAT) → worked example "Hematopoietic Stem Cell Differentiation" (HSC → HSC + HSC, HSC → CMP, erythropoiesis CMP → proerythroblast → erythroblast → reticulocyte → erythrocyte) and the lineage graph HSC → CMP/CLP → blood cells.',

  async create(ctx) {
    const { stage, root, panel, assets, legend, setStatus } = ctx;
    const rbcData = await assets.loadJSON('rbc.json').catch(() => null);
    const lymphNucR = (rbcData?.lymphocyte?.lymphocyte_nucleus_diameter_um ?? 2 * LYMPH_NUC_R_FALLBACK) / 2;

    // ============================================================ materials & pools
    const sphere = unitSphere(3), sphereLo = unitSphere(1), sphereMid = unitSphere(2);
    const geoRng = new RNG(4242);
    const gr = () => geoRng.next();
    const G = {
      lobed4: smooth(lobedNucleus(4, 1, gr)), lobed3: smooth(lobedNucleus(3, 1, gr)), lobed2: smooth(lobedNucleus(2, 1.05, gr)),
      indent: smooth(indentedNucleus(1)), mkNuc: smooth(blobGeometry(1, 0.34, 2.3, 4)), mkBody: smooth(blobGeometry(1, 0.06, 5.1, 3)),
      retic: smooth(blobGeometry(1, 0.1, 7.7, 3)), rbc: erythrocyteGeometry(7.82, 26), mac: smooth(blobGeometry(1, 0.2, 3.3, 3)),
      shaft: new THREE.CylinderGeometry(1, 1, 1, 7, 1, true).rotateZ(-Math.PI / 2).translate(0.5, 0, 0),
    };
    const M = {
      blast: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.26, he: [0.3, 0.3], gain: 0.2 }),
      hsc: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.4, emissive: 0.45, rimPower: 1.6, he: [0.34, 0.3], gain: 0.55 }),
      gran: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.3, he: [0.05, 0.5], gain: 0.17 }),
      eryEarly: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.5, he: [0.32, 0.32], gain: 0.22 }),
      eryLate: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.46, he: [0.03, 0.85], rimPower: 1.4, gain: 0.24 }),
      mkBody: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.42, he: [0.1, 0.55], noise: 0.6, noiseScale: 0.9, gain: 0.3 }),
      nucOpen: cellMaterial({ role: 'nucleus', color: 0xa9b8ff, stain: 0x4f86ff, he: [0.8, 0.12], gain: 0.3 }),
      nucDense: cellMaterial({ role: 'nucleus', color: 0x7a6fd8, stain: 0x6a9bff, he: [1.6, 0.18], gain: 0.32 }),
      nucLobed: cellMaterial({ role: 'nucleus', color: 0x8a7ce0, stain: 0x6a9bff, he: [1.45, 0.16], gain: 0.32 }),
      nucMK: cellMaterial({ role: 'nucleus', color: 0x9a8ef0, stain: 0x4f86ff, he: [1.15, 0.14], noiseScale: 0.9, gain: 0.35 }),
      retic: cellMaterial({ role: 'rbc', color: 0xffffff, stain: 0xff6a4a, gain: 0.8 }),
      rbc: cellMaterial({ role: 'rbc', color: 0xffffff, stain: 0xff5533, gain: 0.8 }),
      plt: cellMaterial({ role: 'reporter', he: [0.55, 0.3], gain: 0.7 }),
      granule: cellMaterial({ role: 'molecule', he: [0.1, 0.9], gain: 0.35 }),
      shaft: cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.55, he: [0.1, 0.55], side: THREE.DoubleSide, gain: 0.6 }),
    };
    const cellGroup = new THREE.Group(), flowGroup = new THREE.Group(), envGroup = new THREE.Group();
    root.add(envGroup, cellGroup, flowGroup);

    class Layer {
      constructor(geo, mat, cap, parent, order = 0) {
        this.pool = new InstancePool(geo, mat, cap, parent);
        this.pool.mesh.renderOrder = order;
        this.owners = [];
      }
      get mesh() { return this.pool.mesh; }
      begin() { this.pool.begin(); }
      put(o, x, y, z, sx, sy, sz, c) { const i = this.pool.put(x, y, z, sx, sy, sz, c); this.owners[i] = o; return i; }
      putO(o, p, ax, a, b, c, col) { const i = this.pool.putOriented(p, ax, a, b, c, col); this.owners[i] = o; return i; }
      putM(o, m, col) { const i = this.pool.putMatrix(m, col); this.owners[i] = o; return i; }
      end() { this.pool.end(); this.owners.length = this.pool.n; }
    }
    const L = {
      blast: new Layer(sphere, M.blast, 256, cellGroup, 4),
      hsc: new Layer(sphere, M.hsc, 32, cellGroup, 4),
      gran: new Layer(sphere, M.gran, 256, cellGroup, 4),
      eryEarly: new Layer(sphere, M.eryEarly, 128, cellGroup, 3),
      eryLate: new Layer(sphere, M.eryLate, 128, cellGroup, 3),
      mkBody: new Layer(G.mkBody, M.mkBody, 8, cellGroup, 4),
      nucOpen: new Layer(sphere, M.nucOpen, 256, cellGroup),
      nucDense: new Layer(sphere, M.nucDense, 256, cellGroup),
      nucIndent: new Layer(G.indent, M.nucLobed, 64, cellGroup),
      lobed4: new Layer(G.lobed4, M.nucLobed, 64, cellGroup),
      lobed3: new Layer(G.lobed3, M.nucLobed, 64, cellGroup),
      lobed2: new Layer(G.lobed2, M.nucLobed, 32, cellGroup),
      mkNuc: new Layer(G.mkNuc, M.nucMK, 8, cellGroup),
      retic: new Layer(G.retic, M.retic, 96, cellGroup),
      rbc: new Layer(G.rbc, M.rbc, 256, flowGroup),
      plt: new Layer(sphereLo, M.plt, 256, flowGroup),
      granule: new Layer(sphereLo, M.granule, 256, cellGroup),
      shaft: new Layer(G.shaft, M.shaft, 64, cellGroup, 4),
    };
    const ALL_LAYERS = Object.values(L);

    // ============================================================ static environment
    const rng = new RNG(20240607);
    const envRng = new RNG(99);

    // trabecular bone: a slab whose top surface follows boneY (endosteum)
    const boneGeo = new THREE.BoxGeometry(200, 1, 84, 150, 1, 64);
    {
      const p = boneGeo.attributes.position;
      for (let i = 0; i < p.count; i++) {
        const x = p.getX(i), z = p.getZ(i);
        p.setY(i, p.getY(i) > 0 ? boneY(x, z) : -64);
      }
      boneGeo.computeVertexNormals();
    }
    const boneMat = cellMaterial({ role: 'solid', color: 0xe8dcc2, stain: 0x3f63ff, he: [0.06, 0.95], noise: 0.55, noiseScale: 0.14, side: THREE.DoubleSide, gain: 0.9 });
    const bone = new THREE.Mesh(boneGeo, boneMat);
    envGroup.add(bone);
    // a trabecular strut arching behind the marrow cavity (decorative depth cue)
    const strutCurve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(-104, -36, -46), new THREE.Vector3(-62, -4, -52), new THREE.Vector3(-10, 14, -50),
      new THREE.Vector3(44, 2, -48), new THREE.Vector3(104, -32, -44),
    ]);
    const strutGeo = new THREE.TubeGeometry(strutCurve, 90, 9, 18, false);
    {
      const p = strutGeo.attributes.position, n = strutGeo.attributes.normal;
      for (let i = 0; i < p.count; i++) {
        const x = p.getX(i), y = p.getY(i), z = p.getZ(i);
        const d = 1.6 * Math.sin(0.21 * x + 0.3 * y) * Math.cos(0.17 * z - 0.2 * x);
        p.setXYZ(i, x + n.getX(i) * d, y + n.getY(i) * d, z + n.getZ(i) * d);
      }
      strutGeo.computeVertexNormals();
    }
    const strut = new THREE.Mesh(strutGeo, boneMat);
    envGroup.add(strut);

    // sinusoid endothelium
    const vesselGeo = new THREE.CylinderGeometry(VES.r, VES.r, VES.x1 - VES.x0, 64, 16, true).rotateZ(Math.PI / 2);
    const vesselMat = cellMaterial({ role: 'membrane', color: 0xffc9c9, stain: 0x36f2d2, opacity: 0.09, rimPower: 2.6, he: [0.05, 0.35], gain: 0.7 });
    const vessel = new THREE.Mesh(vesselGeo, vesselMat);
    vessel.position.set(0, VES.y, VES.z);
    vessel.renderOrder = 6;
    envGroup.add(vessel);

    const staticNuc = new InstancePool(sphereMid, M.nucDense, 160, envGroup);
    const staticBody = new InstancePool(sphere, cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.45, he: [0.28, 0.3], gain: 0.2 }), 128, envGroup);
    const staticShaft = new InstancePool(G.shaft, cellMaterial({ role: 'cytoplasm', useInst: 1, opacity: 0.5, he: [0.1, 0.5], side: THREE.DoubleSide, gain: 0.5 }), 128, envGroup);
    staticBody.mesh.renderOrder = 4; staticShaft.mesh.renderOrder = 4;
    const _m4 = new THREE.Matrix4(), _q = new THREE.Quaternion(), _v = new THREE.Vector3(), _s = new THREE.Vector3();
    const _bx = new THREE.Vector3(), _by = new THREE.Vector3(), _bz = new THREE.Vector3();
    const _up = new THREE.Vector3(0, 1, 0);
    const col = new THREE.Color(), col2 = new THREE.Color(), white = new THREE.Color(0xffffff);

    staticNuc.begin(); staticBody.begin(); staticShaft.begin();
    // endothelial nuclei, flattened against the vessel wall and aligned with flow
    for (let i = 0; i < 20; i++) {
      const x = VES.x0 + 8 + (i + envRng.next() * 0.6) * ((VES.x1 - VES.x0 - 16) / 20);
      const th = envRng.uniform(0, Math.PI * 2);
      const ny = Math.cos(th), nz = Math.sin(th);
      _bx.set(1, 0, 0); _by.set(0, ny, nz); _bz.set(0, -nz, ny);
      _m4.makeBasis(_bx, _by, _bz).scale(_s.set(4.2, 0.8, 1.9));
      _m4.setPosition(x, VES.y + ny * (VES.r + 0.4), VES.z + nz * (VES.r + 0.4));
      staticNuc.putMatrix(_m4, white);
    }
    // osteoblasts / bone-lining cells on the endosteum; osteocytes in lacunae inside bone
    const macrophages = [
      { x: -36, y: 0, z: 9, R: 8.5, id: 'M1' },
      { x: 30, y: 0, z: 11, R: 8.5, id: 'M2' },
    ];
    for (const m of macrophages) m.y = boneY(m.x, m.z) + m.R + 6;
    for (let i = 0; i < 40; i++) {
      const x = envRng.uniform(-92, 92), z = envRng.uniform(-36, 36);
      const y = boneY(x, z);
      col.setRGB(0.62, 0.62, 0.95);
      staticBody.put(x, y + 1.4, z, 5.2, 1.9, 3.8, col);
      staticNuc.put(x, y + 1.6, z, 2.1, 1.2, 1.7, white);
    }
    for (let i = 0; i < 50; i++) {
      const x = envRng.uniform(-96, 96), z = envRng.uniform(-40, 40);
      staticNuc.put(x, boneY(x, z) - envRng.uniform(4, 24), z, 2.2, 0.9, 1.4, white);
    }
    // CXCL12-abundant reticular (CAR) stromal cells with long processes
    for (let i = 0; i < 18; i++) {
      const x = envRng.uniform(-80, 80), z = envRng.uniform(-26, 26);
      const y = envRng.uniform(boneY(x, z) + 8, VES.y - 4);
      if (Math.hypot(y - VES.y, z - VES.z) < VES.r + 6) continue;
      if (macrophages.some((m) => Math.hypot(x - m.x, y - m.y, z - m.z) < m.R + 8)) continue;
      col.setRGB(0.85, 0.8, 0.95);
      staticBody.put(x, y, z, 3.6, 2.4, 3.2, col);
      staticNuc.put(x, y, z, 2.4, 1.6, 2.0, white);
      const nProc = 3 + envRng.int(3);
      for (let k = 0; k < nProc; k++) {
        const d = envRng.onSphere([0, 0, 0]);
        if (y + d[1] * 22 > VES.y) d[1] = -Math.abs(d[1]);
        staticShaft.putOriented([x, y, z], d, envRng.uniform(12, 24), 0.5, 0.5, col.setRGB(0.8, 0.76, 0.92));
      }
    }
    // adipocytes: large lipid vacuoles with a flattened peripheral nucleus
    const adipocytes = [
      { x: 68, y: 6, z: -24, R: 15 }, { x: -66, y: 4, z: -22, R: 14 }, { x: 74, y: -16, z: 20, R: 11 },
    ];
    const adipoMat = cellMaterial({ role: 'membrane', color: 0xfff1c4, stain: 0xffd94a, opacity: 0.1, rimPower: 2.8, he: [0.03, 0.3], gain: 0.45 });
    for (const a of adipocytes) {
      const mesh = new THREE.Mesh(new THREE.IcosahedronGeometry(a.R, 4), adipoMat);
      mesh.position.set(a.x, a.y, a.z);
      mesh.renderOrder = 5;
      mesh.userData.info = 'Adipocyte · a single lipid droplet (dissolved away in H&E sections)';
      envGroup.add(mesh);
      a.mesh = mesh;
      _v.set(0.4, 0.3, 0.86).normalize();
      _q.setFromUnitVectors(_up, _v);
      _m4.compose(new THREE.Vector3(a.x + _v.x * (a.R - 0.8), a.y + _v.y * (a.R - 0.8), a.z + _v.z * (a.R - 0.8)), _q, _s.set(3.6, 0.9, 2.2));
      staticNuc.putMatrix(_m4, white);
    }
    // island macrophages
    const macMat = cellMaterial({ role: 'cytoplasm', color: 0xf2cdb4, stain: 0x9dff5a, opacity: 0.38, he: [0.1, 0.5], noise: 0.4, gain: 0.28 });
    for (const m of macrophages) {
      const mesh = new THREE.Mesh(G.mac, macMat);
      mesh.position.set(m.x, m.y, m.z);
      mesh.scale.setScalar(m.R);
      mesh.renderOrder = 4;
      envGroup.add(mesh);
      const nuc = new THREE.Mesh(G.indent, M.nucLobed);
      nuc.position.set(m.x - 2, m.y + 1, m.z - 1);
      nuc.scale.set(4.6, 3.6, 3.8);
      envGroup.add(nuc);
      m.mesh = mesh;
    }
    staticNuc.end(); staticBody.end(); staticShaft.end();

    // HSC niche slots: endosteal (on the bone) and perivascular (next to the sinusoid)
    const slots = [];
    for (let i = 0; slots.length < K3D - 6 && i < 600; i++) {
      const x = rng.uniform(-80, 80), z = rng.uniform(-25, 25);
      if (macrophages.some((m) => Math.hypot(x - m.x, z - m.z) < 16)) continue;
      if (slots.some((s) => Math.hypot(x - s.x, z - s.z) < 12)) continue;
      slots.push({ x, y: boneY(x, z) + 4.6, z, kind: 'endosteal', occ: null });
    }
    for (let i = 0; i < 6; i++) {
      const x = -72 + i * 29 + rng.uniform(-4, 4);
      const th = -Math.PI / 2 + (i % 2 ? -0.7 : 0.7) + rng.uniform(-0.2, 0.2); // below the vessel, front or back
      slots.push({ x, y: VES.y + Math.sin(th) * (VES.r + 4.8), z: VES.z + Math.cos(th) * (VES.r + 4.8), kind: 'perivascular', occ: null });
    }

    // ============================================================ signals
    const params = {
      EPO: 0.5, GM: 1, SCF: 1, sigma: 0.1, a: SSA_DEFAULTS.selfRenewal,
      colorBy: 'lineage', showBlood: true, showStroma: true, section: 'confocal', focusZ: 2,
    };
    const sig = { S5: 0, erk: 0, E1: 0, uG: 0, uP: 0, fate: null, dirty: true, lastChange: -1 };
    const tp = { ...toggle.defaults };
    const fateRng = new RNG(31);

    // ============================================================ agents
    let cells = [], flow = [], nextId = 1, simH = 0, released = { RBC: 0, PLT: 0, WBC: 0 };
    let ssa = null;

    const stageH = (type) => (LINEAGE[type]?.stageDays ?? 0) * 24 * Math.exp(0.25 * rng.normal() - 0.03);
    const divisionsOf = (type) => LINEAGE[type]?.divisions ?? 0;
    const SPEED_TYPES = new Set(['CMP', 'CLP', 'MEP', 'GMP', 'ProEB', 'EB']);
    const ANUCLEATE = new Set(['Retic', 'RBC', 'Pyr']);

    function radiusOf(c) {
      const p = c.stageD > 0 ? Math.min(c.stageT / c.stageD, 1) : 0;
      switch (c.type) {
        case 'EB': return THREE.MathUtils.lerp(LINEAGE.EB.diameter, LINEAGE.EB.diameterEnd, p) / 2;
        case 'MK': return 6 * 2 ** ((4 * Math.min(p / 0.75, 1)) / 3);
        case 'Neu': return THREE.MathUtils.lerp(7.4, 6.4, p);
        case 'B': case 'T': return (c.type === 'T' ? 0.95 : 1) * lymphNucR / 0.9;
        case 'NK': return lymphNucR / 0.82;
        case 'Pyr': return 3.0;
        default: return (LINEAGE[c.type]?.diameter ?? 8) / 2;
      }
    }

    function newCell(type, x, y, z, from = null) {
      const c = {
        id: nextId++, type, x, y, z, fx: 0, fy: 0, fz: 0, R: 4,
        stageT: 0, stageD: stageH(type), divDone: 0, k: divisionsOf(type),
        g: from ? from.g : 1, u: from ? from.u : 1, born: simH, gen: from ? from.gen + 1 : 0,
        mit: -1, phase: 'live', phaseT: 0, ax: rng.onSphere([0, 0, 0]), seed: rng.next(),
        slot: -1, island: -1, egress: false, cycleT: 0, sinceDiv: 99, pp: null, anchor: null,
      };
      if (type === 'HSC') c.cycleT = LINEAGE.HSC.cycleDays * 24 * (0.7 + 0.6 * rng.next());
      c.R = radiusOf(c);
      return c;
    }

    function setType(c, type) {
      c.type = type;
      c.stageT = 0; c.stageD = stageH(type); c.divDone = 0; c.k = divisionsOf(type);
      if (type === 'ProEB' || type === 'EB') c.island = nearestIsland(c);
      if (type === 'MK') assignMKAnchor(c);
      if (type === 'RBC') { c.egress = true; }
    }

    function nearestIsland(c) {
      let best = 0, bd = Infinity;
      macrophages.forEach((m, i) => { const d = Math.hypot(c.x - m.x, c.y - m.y, c.z - m.z); if (d < bd) { bd = d; best = i; } });
      return best;
    }

    function assignMKAnchor(c) {
      // megakaryocytes mature against the sinusoid wall (front / lower side of the vessel)
      const taken = cells.filter((o) => o.type === 'MK' && o !== c && o.anchor).map((o) => o.anchor.x);
      let x = THREE.MathUtils.clamp(c.x, -68, 68);
      for (let tries = 0; tries < 16 && taken.some((t) => Math.abs(t - x) < 34); tries++) x = rng.uniform(-68, 68);
      // below the sinusoid, on its front or back flank
      const th = (rng.next() < 0.8 ? 1 : -1) * rng.uniform(0.5, 0.8);
      c.anchor = { x, ny: -Math.cos(th), nz: Math.sin(th) };
    }

    // stochastic rounding of mean-field counts into a warm-start population
    function seedPopulation() {
      cells = []; flow = []; nextId = 1; simH = 0;
      released = { RBC: 0, PLT: 0, WBC: 0 };
      for (const s of slots) s.occ = null;
      const mf = meanFieldSteadyState({ ...ssa.env, K: K3D });
      const start = toggle.fixedPoints({ ...tp, uG: 0, uP: 0 }).find((f) => f.stable && Math.abs(f.g - f.u) < 1e-3) ?? { g: 1, u: 1 };
      const committedE = { g: 1 + sig.uG, u: 0.06 }, committedM = { g: 0.06, u: 1 + sig.uP };
      const want = (id) => { const v = (mf[id] ?? 0) * (id === 'HSC' ? 1 : marrowFraction(id)); return Math.floor(v) + (rng.next() < v - Math.floor(v) ? 1 : 0); };
      const place = (c) => {
        let x, y, z;
        for (let t = 0; t < 30; t++) {
          x = rng.uniform(DOM.x0 + 6, DOM.x1 - 6); z = rng.uniform(DOM.z0 + 6, DOM.z1 - 6);
          y = rng.uniform(boneY(x, z) + c.R + 1, DOM.top - c.R);
          if (Math.hypot(y - VES.y, z - VES.z) > VES.r + c.R + 1 && !adipocytes.some((a) => Math.hypot(x - a.x, y - a.y, z - a.z) < a.R + c.R)) break;
        }
        c.x = x; c.y = y; c.z = z;
      };
      const nH = Math.min(want('HSC'), slots.length);
      for (let i = 0; i < nH; i++) {
        const s = slots[(i * 7) % slots.length].occ ? slots.find((q) => !q.occ) : slots[(i * 7) % slots.length];
        const c = newCell('HSC', s.x, s.y, s.z);
        c.g = start.g; c.u = start.u; c.slot = slots.indexOf(s); s.occ = c;
        c.cycleT *= rng.next();
        cells.push(c);
      }
      for (const id of ['CMP', 'CLP', 'MEP', 'GMP', 'ProEB', 'EB', 'Retic', 'MK', 'Neu', 'Mono', 'Eos', 'Baso', 'B', 'T', 'NK']) {
        const n = want(id);
        for (let i = 0; i < n; i++) {
          const c = newCell(id, 0, 0, 0);
          c.stageT = rng.next() * c.stageD * (id === 'MK' ? 0.95 : 1);
          c.divDone = Math.min(c.k, Math.floor(c.k * c.stageT / Math.max(c.stageD, 1e-9) + 0.5));
          const grp = TYPE_GROUP[id];
          const st = grp === 'ery' || id === 'MEP' || id === 'MK' ? committedE : grp === 'lymph' ? { g: 0.03, u: 0.35 } : id === 'CMP' ? start : committedM;
          c.g = st.g + 0.05 * rng.normal(); c.u = st.u + 0.05 * rng.normal();
          c.R = radiusOf(c);
          if (grp === 'ery') {
            c.island = rng.int(macrophages.length);
            const m = macrophages[c.island];
            const d = rng.onSphere([0, 0, 0]);
            if (d[1] < -0.3) d[1] = -d[1] * 0.5;
            const rr = m.R + c.R + 1 + (id === 'Retic' ? rng.uniform(4, 14) : rng.uniform(0, 5));
            c.x = m.x + d[0] * rr; c.y = Math.max(m.y + d[1] * rr, boneY(m.x, m.z) + c.R); c.z = m.z + d[2] * rr;
          } else if (id === 'MK') {
            assignMKAnchor(c);
            const R = c.R;
            c.x = c.anchor.x; c.y = VES.y + c.anchor.ny * (VES.r + R + 0.5); c.z = VES.z + c.anchor.nz * (VES.r + R + 0.5);
          } else place(c);
          cells.push(c);
        }
      }
      // background blood in the sinusoid
      for (let i = 0; i < 115; i++) flow.push(newFlow('rbc', rng.uniform(VES.x0, VES.x1), null, true));
      for (let i = 0; i < 24; i++) flow.push(newFlow('plt', rng.uniform(VES.x0, VES.x1), null, true));
    }

    function newFlow(kind, s, agent = null, bg = false, r = null, th = null) {
      const wall = VES.r - (kind === 'plt' ? 1.5 : 4.2);
      const r0 = kind === 'plt' ? wall * Math.sqrt(0.55 + 0.45 * rng.next()) : wall * Math.sqrt(rng.next()) * 0.95;
      return {
        kind, s, agent, bg, r: r ?? r0, r0, th: th ?? rng.uniform(0, Math.PI * 2), spin: rng.uniform(0, 6.3),
        spinRate: rng.uniform(0.6, 1.4), tilt: rng.uniform(-0.6, 0.6), x: 0, y: 0, z: 0,
      };
    }

    // ------------------------------------------------------------ neighbour grid (no per-frame allocation)
    const GRID = { x0: -100, y0: -42, z0: -40, h: 18, nx: 12, ny: 5, nz: 5 };
    const head = new Int32Array(GRID.nx * GRID.ny * GRID.nz);
    let nextIdx = new Int32Array(1024);
    const bucket = (x, y, z) => {
      const ix = Math.min(GRID.nx - 1, Math.max(0, Math.floor((x - GRID.x0) / GRID.h)));
      const iy = Math.min(GRID.ny - 1, Math.max(0, Math.floor((y - GRID.y0) / GRID.h)));
      const iz = Math.min(GRID.nz - 1, Math.max(0, Math.floor((z - GRID.z0) / GRID.h)));
      return (ix * GRID.ny + iy) * GRID.nz + iz;
    };

    const mkList = [];
    function mechanics(dt) {
      const n = cells.length;
      if (nextIdx.length < n) nextIdx = new Int32Array(n * 2);
      head.fill(-1);
      mkList.length = 0;
      for (let i = 0; i < n; i++) {
        const c = cells[i];
        c.fx = c.fy = c.fz = 0;
        c.ox = c.oy = c.oz = 0;
        if (c.type === 'MK') { c.bk = -1; mkList.push(c); continue; }
        if (c.egress || c.type === 'Pyr') { c.bk = -1; continue; } // squeeze between cells (not part of the packing)
        const b = bucket(c.x, c.y, c.z);
        c.bk = b;
        nextIdx[i] = head[b]; head[b] = i;
      }
      const kRep = 2.6;
      for (let i = 0; i < n; i++) {
        const a = cells[i];
        if (a.bk < 0) continue;
        const ix = Math.floor(a.bk / (GRID.ny * GRID.nz)), iy = Math.floor(a.bk / GRID.nz) % GRID.ny, iz = a.bk % GRID.nz;
        for (let dx = -1; dx <= 1; dx++) {
          const jx = ix + dx; if (jx < 0 || jx >= GRID.nx) continue;
          for (let dy = -1; dy <= 1; dy++) {
            const jy = iy + dy; if (jy < 0 || jy >= GRID.ny) continue;
            for (let dz = -1; dz <= 1; dz++) {
              const jz = iz + dz; if (jz < 0 || jz >= GRID.nz) continue;
              for (let j = head[(jx * GRID.ny + jy) * GRID.nz + jz]; j >= 0; j = nextIdx[j]) {
                if (j <= i) continue;
                const b = cells[j];
                const ddx = b.x - a.x, ddy = b.y - a.y, ddz = b.z - a.z;
                const s = a.R + b.R;
                const d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                if (d2 >= s * s) continue;
                const d = Math.sqrt(d2) || 1e-3;
                const f = (kRep * (s - d)) / d;
                a.fx -= f * ddx; a.fy -= f * ddy; a.fz -= f * ddz;
                b.fx += f * ddx; b.fy += f * ddy; b.fz += f * ddz;
              }
            }
          }
        }
      }
      const noise = Math.sqrt(2 * 0.5 * dt);
      for (let i = 0; i < n; i++) {
        const c = cells[i];
        const R = c.R;
        let fx = c.fx, fy = c.fy, fz = c.fz;
        // large obstacles: megakaryocytes, island macrophages, adipocytes
        if (c.type !== 'MK') for (const o of mkList) obstacle(c, o.x, o.y, o.z, o.R, 2.5);
        const eryIsland = TYPE_GROUP[c.type] === 'ery';
        for (const m of macrophages) {
          if (c.type !== 'Pyr') obstacle(c, m.x, m.y, m.z, m.R * 0.95, 2.5);
          // erythroblastic islands are erythroid niches: other lineages are gently kept out of the rosette
          if (!eryIsland && c.type !== 'MK') obstacle(c, m.x, m.y, m.z, m.R + 13, 0.9);
        }
        for (const a of adipocytes) obstacle(c, a.x, a.y, a.z, a.R, 3);
        fx += c.ox; fy += c.oy; fz += c.oz;
        // bone, walls, ceiling
        const by = boneY(c.x, c.z) + R;
        if (c.y < by) fy += 5 * (by - c.y);
        if (c.y > DOM.top - R) fy -= 4 * (c.y - DOM.top + R);
        if (c.x < DOM.x0 + R) fx += 4 * (DOM.x0 + R - c.x);
        if (c.x > DOM.x1 - R) fx -= 4 * (c.x - DOM.x1 + R);
        if (c.z < DOM.z0 + R) fz += 4 * (DOM.z0 + R - c.z);
        if (c.z > DOM.z1 - R) fz -= 4 * (c.z - DOM.z1 + R);
        // sinusoid wall (crossed only by cells that are egressing)
        const vy = c.y - VES.y, vz = c.z - VES.z, vr = Math.hypot(vy, vz) || 1e-3;
        if (!c.egress && vr < VES.r + R) { const pen = VES.r + R - vr; fy += (4 * pen * vy) / vr; fz += (4 * pen * vz) / vr; }
        // chemotaxis / homing
        const tgt = targetOf(c);
        if (tgt) {
          const tx = tgt[0] - c.x, ty = tgt[1] - c.y, tz = tgt[2] - c.z;
          const d = Math.hypot(tx, ty, tz) || 1e-3;
          const sp = Math.min(tgt[3] * d, tgt[4]);
          fx += (sp * tx) / d; fy += (sp * ty) / d; fz += (sp * tz) / d;
        }
        const mob = c.type === 'MK' ? 0.15 : c.type === 'HSC' ? 0.3 : 1;
        c.x += fx * dt + noise * mob * rng.normal();
        c.y += fy * dt + noise * mob * rng.normal();
        c.z += fz * dt + noise * mob * rng.normal();
      }
    }

    function obstacle(c, x, y, z, R, k) {
      const dx = c.x - x, dy = c.y - y, dz = c.z - z;
      const s = R + c.R;
      const d2 = dx * dx + dy * dy + dz * dz;
      if (d2 >= s * s) return;
      const d = Math.sqrt(d2) || 1e-3;
      const f = (k * (s - d)) / d;
      c.ox += f * dx; c.oy += f * dy; c.oz += f * dz;
    }

    const _tgt = [0, 0, 0, 0, 0];
    // returns [x, y, z, gain(1/s), max speed(um/s)] or null
    function targetOf(c) {
      const t = _tgt;
      if (c.type === 'HSC' && c.slot >= 0) { const s = slots[c.slot]; t[0] = s.x; t[1] = s.y; t[2] = s.z; t[3] = 0.8; t[4] = 6; return t; }
      if (c.type === 'Pyr') { const m = macrophages[c.island]; t[0] = m.x; t[1] = m.y; t[2] = m.z; t[3] = 1.2; t[4] = 9; return t; }
      if (c.type === 'MK' && c.anchor) {
        const Rr = VES.r + c.R + 0.3;
        t[0] = c.anchor.x; t[1] = VES.y + c.anchor.ny * Rr; t[2] = VES.z + c.anchor.nz * Rr; t[3] = 0.6; t[4] = 3; return t;
      }
      if (c.egress) {
        // nearest point on the sinusoid axis: crossing the endothelium into the lumen
        t[0] = c.x; t[1] = VES.y; t[2] = VES.z; t[3] = 1; t[4] = 14; return t;
      }
      if ((c.type === 'ProEB' || c.type === 'EB' || c.phase === 'enuc') && c.island >= 0) {
        const m = macrophages[c.island];
        const dx = c.x - m.x, dy = c.y - m.y, dz = c.z - m.z, d = Math.hypot(dx, dy, dz) || 1;
        const ring = m.R + c.R + 0.5;
        t[0] = m.x + (dx / d) * ring; t[1] = m.y + (dy / d) * ring; t[2] = m.z + (dz / d) * ring; t[3] = 0.4; t[4] = 3; return t;
      }
      if (c.type === 'Retic') {
        const m = macrophages[c.island >= 0 ? c.island : 0];
        // drift from the island toward the vessel as the reticulocyte matures
        const p = Math.min(c.stageT / c.stageD, 1);
        t[0] = THREE.MathUtils.lerp(m.x, c.x, 0.9); t[1] = THREE.MathUtils.lerp(m.y + 8, VES.y - VES.r - 5, p);
        t[2] = THREE.MathUtils.lerp(m.z, VES.z + VES.r + 5, p); t[3] = 0.2; t[4] = 4; return t;
      }
      return null;
    }

    // ------------------------------------------------------------ biology
    const _rates = [0, 0];
    function tfStep(c, dtH) {
      if (ANUCLEATE.has(c.type)) return;
      const grp = TYPE_GROUP[c.type];
      const dtu = dtH / TF_HOURS;
      if (grp === 'lymph' || c.type === 'CLP') {
        // lymphoid programme: GATA1 off, low PU.1 (not modelled by the toggle)
        const k = 1 - Math.exp(-dtu);
        c.g += (0.03 - c.g) * k + params.sigma * 0.5 * Math.sqrt(dtu) * rng.normal();
        c.u += (0.35 - c.u) * k + params.sigma * 0.5 * Math.sqrt(dtu) * rng.normal();
        c.g = Math.abs(c.g); c.u = Math.abs(c.u);
        return;
      }
      let a = 0;
      if (c.type === 'HSC') a = 1;
      else if (c.type === 'CMP') a = Math.exp(-(c.stageT / TF_HOURS) / 2);
      const sub = Math.max(1, Math.ceil(dtu / 0.05));
      const h = dtu / sub;
      for (let s = 0; s < sub; s++) toggle.langevinStep(c, h, tp, params.sigma, rng, a, _rates);
    }

    function hscDivide(c) {
      // each daughter stays a stem cell with probability a (1 - H/K): it contacts a random niche slot
      // and keeps stem identity only if that slot is free (Schofield niche; Marciniak-Czochra 2009 fraction a)
      const H = cells.reduce((n, o) => n + (o.type === 'HSC' ? 1 : 0), 0);
      const pStay = params.a * Math.max(0, 1 - H / K3D);
      const daughters = [c, newCell('HSC', c.x + c.ax[0] * 2, c.y + Math.abs(c.ax[1]) * 2, c.z + c.ax[2] * 2, c)];
      const d2 = daughters[1];
      d2.g = c.g + 0.05 * rng.normal(); d2.u = c.u + 0.05 * rng.normal();
      const mySlot = c.slot;
      if (mySlot >= 0) slots[mySlot].occ = null;
      c.slot = -1;
      for (const d of daughters) {
        if (rng.next() < pStay) {
          let s = mySlot >= 0 && !slots[mySlot].occ ? mySlot : -1;
          if (s < 0) {
            let bd = Infinity;
            slots.forEach((q, i) => { if (!q.occ) { const dd = Math.hypot(q.x - d.x, q.z - d.z); if (dd < bd) { bd = dd; s = i; } } });
          }
          if (s >= 0) { d.slot = s; slots[s].occ = d; d.type = 'HSC'; d.cycleT = LINEAGE.HSC.cycleDays * 24 * (0.8 + 0.4 * rng.next()); continue; }
        }
        setType(d, rng.next() < ssa.env.fCLP ? 'CLP' : 'CMP');
        d.stageT = 0;
      }
      cells.push(d2);
      d2.sinceDiv = c.sinceDiv = 0;
    }

    function divide(c) {
      const d = newCell(c.type, c.x + c.ax[0] * c.R * 0.45, c.y + c.ax[1] * c.R * 0.45, c.z + c.ax[2] * c.R * 0.45, c);
      c.x -= c.ax[0] * c.R * 0.45; c.y -= c.ax[1] * c.R * 0.45; c.z -= c.ax[2] * c.R * 0.45;
      d.stageT = c.stageT; d.stageD = c.stageD * (0.85 + 0.3 * rng.next()); d.divDone = c.divDone; d.k = c.k;
      d.g = c.g + 0.03 * rng.normal(); d.u = c.u + 0.03 * rng.normal();
      d.island = c.island; d.gen = c.gen + 1; c.gen++;
      d.sinceDiv = c.sinceDiv = 0;
      d.R = c.R;
      cells.push(d);
    }

    function stageEnd(c) {
      switch (c.type) {
        case 'CMP': setType(c, c.g > c.u ? 'MEP' : 'GMP'); break;
        case 'MEP': setType(c, mepFate(rng)); break;
        case 'GMP': setType(c, gmpFate(c.g, c.u, rng)); break;
        case 'CLP': setType(c, clpFate(rng)); break;
        case 'ProEB': setType(c, 'EB'); break;
        case 'EB': c.phase = 'enuc'; c.phaseT = 0; setEnucAxis(c); break;
        case 'Retic': setType(c, 'RBC'); break;
        case 'MK': c.phase = 'dying'; c.phaseT = 0; break;
        default: c.egress = true; // mature leukocytes leave through the sinusoid
      }
    }

    function setEnucAxis(c) {
      // the nucleus is extruded toward the island macrophage, which engulfs the pyrenocyte
      const m = macrophages[c.island >= 0 ? c.island : 0];
      const dx = c.x - m.x, dy = c.y - m.y, dz = c.z - m.z, d = Math.hypot(dx, dy, dz) || 1;
      c.ax = [-dx / d, -dy / d, -dz / d];
    }

    const died = [];
    const apoRate = {};
    function biology(dtH) {
      for (const id of Object.keys(LINEAGE)) apoRate[id] = apoptosisRate(id, ssa.env);
      const fERK = erkSpeed(sig.erk);
      let n = cells.length;
      for (let i = 0; i < n; i++) {
        const c = cells[i];
        c.sinceDiv += dtH;
        tfStep(c, dtH);
        if (c.phase === 'apop') { c.phaseT += dtH; if (c.phaseT > APOP_H) died.push(c); continue; }
        if (c.phase === 'dying') {
          c.phaseT += dtH;
          if (c.phaseT > 6) died.push(c);
          continue;
        }
        if (c.phase === 'enuc') {
          c.phaseT += dtH;
          if (c.phaseT >= ENUC_H) {
            // pyrenocyte leaves; the cell continues as a reticulocyte
            const off = c.R * 1.1;
            const pyr = newCell('Pyr', c.x + c.ax[0] * off, c.y + c.ax[1] * off, c.z + c.ax[2] * off, c);
            pyr.island = c.island;
            cells.push(pyr);
            c.phase = 'live';
            setType(c, 'Retic');
            c.R = radiusOf(c);
          }
          continue;
        }
        if (c.type === 'Pyr') {
          // phosphatidylserine-exposing pyrenocytes are engulfed by the island macrophage
          const m = macrophages[c.island];
          c.phaseT += dtH;
          if (Math.hypot(c.x - m.x, c.y - m.y, c.z - m.z) < m.R + 2) c.eaten = (c.eaten ?? 0) + dtH;
          if ((c.eaten ?? 0) > 3 || c.phaseT > 16) died.push(c);
          continue;
        }
        if (c.mit >= 0) {
          c.mit += dtH;
          if (c.mit >= MITOSIS_H) {
            c.mit = -1;
            if (c.type === 'HSC') hscDivide(c);
            else { c.divDone++; divide(c); }
          }
          continue;
        }
        if (c.type === 'HSC') {
          c.cycleT -= dtH;
          if (c.cycleT <= 0 && cells.length < MAX_CELLS) { c.mit = 0; c.ax = rng.onSphere(c.ax); if (c.ax[1] < 0) c.ax[1] *= -1; }
          continue;
        }
        if (c.type === 'RBC') continue; // waiting to cross into the sinusoid
        // growth-factor withdrawal apoptosis: EPO for erythroid progenitors (Koury & Bondurant 1990),
        // SCF/c-Kit for c-Kit+ progenitors (same hazards as the compartment model)
        const hz = apoRate[c.type];
        if (hz > 0 && rng.next() < 1 - Math.exp(-(hz / 24) * dtH)) { c.phase = 'apop'; c.phaseT = 0; continue; }
        if (c.egress) continue;
        c.stageT += dtH * (SPEED_TYPES.has(c.type) ? fERK : 1);
        const p = c.stageT / c.stageD;
        if (c.k > 0 && c.divDone < c.k && p >= (c.divDone + 0.5) / c.k) {
          if (cells.length < MAX_CELLS) { c.mit = 0; c.ax = rng.onSphere(c.ax); }
          else c.stageT = Math.min(c.stageT, ((c.divDone + 0.5) / c.k) * c.stageD); // crowded: wait
          continue;
        }
        if (c.type === 'MK') shedPlatelets(c, dtH);
        if (p >= 1) stageEnd(c);
      }
      if (died.length) {
        for (const d of died) { if (d.slot >= 0) slots[d.slot].occ = null; d.dead = true; }
        cells = cells.filter((c) => !c.dead);
        died.length = 0;
      }
      for (const c of cells) {
        if (c.phase === 'apop') c.R = Math.max(radiusOf(c) * (1 - 0.6 * c.phaseT / APOP_H), 1);
        else if (c.phase === 'enuc') c.R = THREE.MathUtils.lerp(LINEAGE.EB.diameterEnd / 2, LINEAGE.Retic.diameter / 2, c.phaseT / ENUC_H);
        else if (c.type === 'Pyr') c.R = 3.0 * (1 - 0.6 * Math.min((c.eaten ?? 0) / 3, 1));
        else c.R = radiusOf(c) * (c.sinceDiv < 2 ? 0.85 + 0.075 * c.sinceDiv : 1);
      }
    }

    // proplatelets: each mature MK pushes 3 beaded extensions into the lumen and sheds their tips
    const NPP = 3;
    function shedPlatelets(c, dtH) {
      const p = c.stageT / c.stageD;
      if (p < 0.75) { c.pp = null; return; }
      if (!c.pp) c.pp = Array.from({ length: NPP }, (_, j) => ({ phase: j / NPP + 0.1 * rng.next(), off: [rng.uniform(-0.45, 0.45), rng.uniform(-0.3, 0.3)] }));
      const yieldN = LINEAGE.MK.plateletYield;
      const rate = yieldN / NPP / (0.25 * c.stageD); // cycles per hour per proplatelet
      for (const pp of c.pp) {
        pp.phase += rate * dtH;
        while (pp.phase >= 1) {
          pp.phase -= 1;
          const tip = proplateletPoint(c, pp, 1, _tip);
          const dy = tip[1] - VES.y, dz = tip[2] - VES.z;
          const r = Math.min(Math.hypot(dy, dz), VES.r - 1.5);
          if (flow.length < 420) flow.push(newFlow('plt', tip[0], null, false, r, Math.atan2(dz, dy)));
          released.PLT++;
        }
      }
    }
    const _tip = [0, 0, 0], _pa = [0, 0, 0];
    // point at fraction s along a proplatelet of MK c: straight into the lumen, bent downstream by flow
    function proplateletPoint(c, pp, s, out) {
      const vy = VES.y - c.y, vz = VES.z - c.z, vr = Math.hypot(vy, vz) || 1;
      let dx = pp.off[0], dy = vy / vr, dz = vz / vr;
      dy += pp.off[1] * (vz / vr); dz -= pp.off[1] * (vy / vr);
      const dn = Math.hypot(dx, dy, dz);
      dx /= dn; dy /= dn; dz /= dn;
      const reach = vr - c.R + VES.r * 0.55;
      const Lp = reach * (0.35 + 0.65 * pp.phase);
      const t = s * Lp;
      out[0] = c.x + dx * c.R * 0.9 + dx * t + 0.5 * (t * t) / reach;
      out[1] = c.y + dy * c.R * 0.9 + dy * t;
      out[2] = c.z + dz * c.R * 0.9 + dz * t;
      return out;
    }

    function egressAndFlow(dt) {
      // cells that reached the lumen join the blood
      for (let i = cells.length - 1; i >= 0; i--) {
        const c = cells[i];
        if (!c.egress) continue;
        const dy = c.y - VES.y, dz = c.z - VES.z, r = Math.hypot(dy, dz);
        if (r < VES.r - 0.3 * c.R) {
          const kind = c.type === 'RBC' ? 'rbc' : 'cell';
          if (flow.length < 440 || kind === 'cell') {
            const e = newFlow(kind, c.x, kind === 'cell' ? c : null, false, Math.min(r, VES.r - 4), Math.atan2(dz, dy));
            e.fresh = true;
            flow.push(e);
          }
          if (kind === 'rbc') released.RBC++; else released.WBC++;
          if (c.slot >= 0) slots[c.slot].occ = null;
          cells.splice(i, 1);
        }
      }
      for (let i = flow.length - 1; i >= 0; i--) {
        const e = flow[i];
        const rr = e.r / VES.r;
        const v = FLOW_VMAX * (1 - rr * rr) * (e.kind === 'cell' ? 0.45 : 1);
        e.s += v * dt;
        // released red cells migrate off the wall (cell-free layer), platelets stay marginated
        e.r += (e.r0 - e.r) * Math.min(1, 0.6 * dt);
        e.spin += e.spinRate * (2 * FLOW_VMAX * rr / VES.r) * dt * 0.6;
        if (e.s > VES.x1 - 2) {
          if (e.bg) { e.s = VES.x0 + 2 + rng.next() * 6; } else { flow[i] = flow[flow.length - 1]; flow.pop(); continue; }
        }
        e.x = e.s; e.y = VES.y + e.r * Math.cos(e.th); e.z = VES.z + e.r * Math.sin(e.th);
        if (e.agent) { e.agent.x = e.x; e.agent.y = e.y; e.agent.z = e.z; }
      }
    }

    // ============================================================ drawing
    const _pos = [0, 0, 0], _ax = [1, 0, 0];
    const reticCol = new THREE.Color();
    const Q = { mode: stage.mode };

    function maturity(c) { return c.stageD > 0 ? Math.min(c.stageT / c.stageD, 1) : 1; }
    function hbOf(c) {
      if (c.type === 'ProEB') return 0.05;
      if (c.type === 'EB') return 0.1 + 0.9 * maturity(c);
      if (c.phase === 'enuc' || c.type === 'Retic' || c.type === 'RBC') return 1;
      return 0;
    }

    // GATA1-mCherry (red) + PU.1-eYFP (yellow-green) fusion reporters, as imaged by Hoppe et al. 2016 (Nature 535:299)
    function reporterColor(c, out) {
      const g = Math.min(c.g / 1.2, 1), u = Math.min(c.u / 1.2, 1);
      return out.setRGB(0.06 + 1.0 * g + 0.72 * u, 0.06 + 0.1 * g + 1.0 * u, 0.1 + 0.28 * g + 0.04 * u);
    }

    function bodyColor(c, out) {
      const mode = Q.mode;
      if (params.colorBy === 'reporter' && !ANUCLEATE.has(c.type)) return reporterColor(c, out);
      const p = maturity(c);
      if (mode === 'histology') {
        switch (TYPE_GROUP[c.type]) {
          case 'ery': return out.setRGB(0.75 + 0.25 * hbOf(c), 0.6 + 0.1 * hbOf(c), 0.95 - 0.25 * hbOf(c));
          case 'gran': return c.type === 'Eos' ? out.setRGB(1, 0.62, 0.45) : c.type === 'Baso' ? out.setRGB(0.72, 0.62, 1) : out.setRGB(1, 0.86, 0.9);
          case 'mk': return out.setRGB(1, 0.8, 0.92);
          case 'lymph': return out.setRGB(0.66, 0.72, 1);
          case 'mono': return out.setRGB(0.82, 0.8, 0.98);
          default: return out.setRGB(0.7, 0.66, 1);
        }
      }
      if (c.type === 'EB' || c.type === 'ProEB' || c.phase === 'enuc') {
        // haemoglobinisation: basophilic blue-grey -> red
        const h = hbOf(c);
        out.setHex(HEX.ProEB);
        return out.lerp(col2.setHex(HEX.RBC), h);
      }
      out.setHex(HEX[c.type] ?? 0xffffff);
      if (c.type === 'Neu' && p < 0.3) out.lerp(col2.setHex(HEX.GMP), 0.5);
      if (mode === 'confocal') out.multiplyScalar(1.15);
      return out;
    }

    const _np = [0, 0, 0], _bp = [0, 0, 0];
    const at = (out, x, y, z) => { out[0] = x; out[1] = y; out[2] = z; return out; };
    function putNucleus(c, layer, x, y, z, sx, sy, sz, ax = null) {
      if (ax) L[layer].putO(c, at(_np, x, y, z), ax, sx, sy, sz, white);
      else L[layer].put(c, x, y, z, sx, sy, sz, white);
    }

    function drawCell(c) {
      const R = c.R;
      const x = c.x, y = c.y, z = c.z;
      const bc = bodyColor(c, col);
      const grp = TYPE_GROUP[c.type];
      const p = maturity(c);
      // mitosis: round up, elongate, pinch into two
      if (c.mit >= 0) {
        const m = c.mit / MITOSIS_H;
        const sep = Math.max(0, (m - 0.45) / 0.55);
        const layer = grp === 'ery' ? (hbOf(c) > 0.5 ? 'eryLate' : 'eryEarly') : grp === 'gran' ? 'gran' : 'blast';
        if (sep <= 0) {
          L[layer].putO(c, at(_bp, x, y, z), c.ax, R * (1 + 0.25 * m), R * (1 - 0.08 * m), R * (1 - 0.08 * m), bc);
          const k = 1 - 0.4 * Math.min(m / 0.45, 1);
          L.nucDense.putO(c, _bp, c.ax, R * 0.18 + R * 0.4 * k * (m > 0.3 ? 0.5 : 1), R * 0.55 * k, R * 0.55 * k, white);
        } else {
          const r2 = R / Math.cbrt(2) * 1.05, d = R * 0.55 * sep;
          for (const sg of [1, -1]) {
            _pos[0] = x + sg * d * c.ax[0]; _pos[1] = y + sg * d * c.ax[1]; _pos[2] = z + sg * d * c.ax[2];
            L[layer].put(c, _pos[0], _pos[1], _pos[2], r2, r2, r2, bc);
            L.nucDense.put(c, _pos[0], _pos[1], _pos[2], r2 * 0.45, r2 * 0.5, r2 * 0.45, white);
          }
        }
        return;
      }
      if (c.phase === 'apop') {
        // apoptosis: shrinkage, condensed and fragmented chromatin (karyorrhexis)
        L.blast.put(c, x, y, z, R, R, R, bc.multiplyScalar(0.6));
        const f = c.phaseT / APOP_H;
        for (let k = 0; k < 3; k++) {
          const a = c.seed * 6.28 + k * 2.1;
          L.nucDense.put(c, x + Math.cos(a) * R * 0.45 * f, y + Math.sin(a * 1.3) * R * 0.3 * f, z + Math.sin(a) * R * 0.45 * f, R * 0.28, R * 0.28, R * 0.28, white);
        }
        return;
      }
      switch (c.type) {
        case 'HSC': case 'CMP': case 'CLP': case 'MEP': case 'GMP': {
          L[c.type === 'HSC' ? 'hsc' : 'blast'].put(c, x, y, z, R, R, R, bc);
          const nr = R * (LINEAGE[c.type].nc ?? 0.7);
          putNucleus(c, 'nucOpen', x + c.ax[0] * R * 0.08, y + c.ax[1] * R * 0.08, z + c.ax[2] * R * 0.08, nr, nr * 0.94, nr);
          return;
        }
        case 'ProEB': case 'EB': {
          const late = hbOf(c) > 0.55;
          L[late ? 'eryLate' : 'eryEarly'].put(c, x, y, z, R, R, R, bc);
          const nc = c.type === 'ProEB' ? LINEAGE.ProEB.nc : THREE.MathUtils.lerp(LINEAGE.EB.nc, LINEAGE.EB.ncEnd, p);
          const nr = R * nc;
          // chromatin condenses and the nucleus becomes eccentric before extrusion
          const ecc = c.type === 'EB' ? R * 0.28 * p : 0;
          putNucleus(c, c.type === 'EB' && p > 0.35 ? 'nucDense' : 'nucOpen', x + c.ax[0] * ecc, y + c.ax[1] * ecc, z + c.ax[2] * ecc, nr, nr, nr);
          return;
        }
        case 'Pyr': {
          // pyrenocyte: condensed nucleus in a thin rim of cytoplasm/membrane
          L.blast.put(c, x, y, z, R * 1.12, R * 1.12, R * 1.12, col.setRGB(0.95, 0.75, 0.8));
          L.nucDense.put(c, x, y, z, R * 0.92, R * 0.92, R * 0.92, white);
          return;
        }
        case 'Retic': {
          // reticulocyte: irregular, then progressively flattened toward the biconcave disc
          const f = 1 - 0.45 * p;
          L.retic.putO(c, at(_bp, x, y, z), c.ax, R * f, R, R * (0.92 + 0.08 * p), reticCol);
          return;
        }
        case 'RBC': {
          drawDisc(L.rbc, c, x, y, z, c.ax, 1, 1);
          return;
        }
        case 'MK': {
          L.mkBody.putO(c, at(_bp, x, y, z), c.ax, R, R * 0.96, R * 1.02, bc);
          const nr = R * 0.5 * (c.phase === 'dying' ? 0.7 : 1);
          putNucleus(c, 'mkNuc', x, y, z, nr * 1.1, nr, nr * 0.95, c.ax);
          if (c.pp && c.phase !== 'dying') drawProplatelets(c, bc);
          return;
        }
        case 'Neu': {
          L.gran.put(c, x, y, z, R, R, R, bc);
          if (p < 0.3) { const nr = R * 0.66; putNucleus(c, 'nucOpen', x, y, z, nr, nr, nr); }
          else if (p < 0.55) { const nr = R * 0.56; putNucleus(c, 'nucDense', x + c.ax[0] * R * 0.2, y + c.ax[1] * R * 0.2, z + c.ax[2] * R * 0.2, nr, nr * 0.9, nr); }
          else if (p < 0.75) putNucleus(c, 'nucIndent', x, y, z, R * 0.5, R * 0.5, R * 0.46, c.ax);
          else if (p < 0.88) putNucleus(c, 'nucIndent', x, y, z, R * 0.66, R * 0.32, R * 0.36, c.ax);
          else putNucleus(c, c.id % 2 ? 'lobed4' : 'lobed3', x, y, z, R * 0.3, R * 0.3, R * 0.3, c.ax);
          return;
        }
        case 'Eos': case 'Baso': {
          L.gran.put(c, x, y, z, R, R, R, bc);
          if (p < 0.45) { const nr = R * 0.6; putNucleus(c, 'nucOpen', x, y, z, nr, nr, nr); }
          else if (c.type === 'Eos') putNucleus(c, 'lobed2', x, y, z, R * 0.36, R * 0.36, R * 0.36, c.ax);
          else putNucleus(c, 'nucIndent', x, y, z, R * 0.45, R * 0.42, R * 0.42, c.ax);
          if (p > 0.25) {
            const gc = c.type === 'Eos' ? col2.setRGB(1, 0.45, 0.2) : col2.setRGB(0.45, 0.3, 0.85);
            const ng = c.type === 'Eos' ? 16 : 12, gr2 = c.type === 'Eos' ? 0.75 : 0.85;
            for (let k = 0; k < ng; k++) {
              const t = c.seed * 17 + k * 2.399963; // golden-angle spiral on a shell
              const zz = 1 - (2 * (k + 0.5)) / ng, rr = Math.sqrt(1 - zz * zz);
              const sr = R * (0.62 + 0.18 * ((k * 7) % 3) / 2);
              L.granule.put(c, x + Math.cos(t) * rr * sr, y + zz * sr, z + Math.sin(t) * rr * sr, gr2, gr2, gr2, gc);
            }
          }
          return;
        }
        case 'Mono': {
          L.blast.put(c, x, y, z, R, R, R, bc);
          if (p < 0.35) { const nr = R * 0.62; putNucleus(c, 'nucOpen', x, y, z, nr, nr, nr); }
          else putNucleus(c, 'nucIndent', x, y, z, R * 0.6, R * 0.52, R * 0.55, c.ax);
          return;
        }
        default: { // lymphocytes: scant cytoplasm around a large round dense nucleus (measured size)
          L.blast.put(c, x, y, z, R, R, R, bc);
          const nr = c.type === 'NK' ? lymphNucR * 0.95 : c.type === 'T' ? lymphNucR * 0.95 : lymphNucR;
          putNucleus(c, 'nucDense', x + c.ax[0] * 0.3, y + c.ax[1] * 0.3, z + c.ax[2] * 0.3, nr, nr * 0.96, nr);
          if (c.type === 'NK') {
            for (let k = 0; k < 5; k++) {
              const t = c.seed * 11 + k * 1.3;
              L.granule.put(c, x + Math.cos(t) * R * 0.85, y + Math.sin(t * 1.7) * R * 0.5, z + Math.sin(t) * R * 0.85, 0.55, 0.55, 0.55, col2.setRGB(0.9, 0.35, 0.6));
            }
          }
        }
      }
    }

    function drawEnucleating(c) {
      // nucleus migrates to the pole, the cell constricts, the pyrenocyte buds off
      const f = Math.min(c.phaseT / ENUC_H, 1);
      const R = c.R;
      const bc = bodyColor(c, col);
      const ax = c.ax;
      const body = R * (1 - 0.1 * f);
      L.eryLate.put(c, c.x - ax[0] * R * 0.15 * f, c.y - ax[1] * R * 0.15 * f, c.z - ax[2] * R * 0.15 * f, body, body * 0.95, body, bc);
      const nr = R * 0.5;
      const d = R * (0.35 + 0.95 * THREE.MathUtils.smoothstep(f, 0.1, 0.95));
      const nx = c.x + ax[0] * d, ny = c.y + ax[1] * d, nz = c.z + ax[2] * d;
      L.nucDense.put(c, nx, ny, nz, nr, nr, nr, white);
      // membrane around the extruding nucleus and the thinning neck
      L.blast.put(c, nx, ny, nz, nr * 1.15, nr * 1.15, nr * 1.15, col2.setRGB(0.95, 0.72, 0.8));
      if (f > 0.35 && f < 0.95) {
        const nk = R * 0.42 * (1 - (f - 0.35) / 0.6);
        L.eryLate.put(c, c.x + ax[0] * d * 0.55, c.y + ax[1] * d * 0.55, c.z + ax[2] * d * 0.55, nk, nk, nk, bc);
      }
    }

    const rbcCol = new THREE.Color();
    function drawDisc(layer, owner, x, y, z, n, s, thick) {
      _v.set(n[0], n[1], n[2]);
      if (_v.lengthSq() < 1e-8) _v.set(0, 1, 0);
      _q.setFromUnitVectors(_up, _v.normalize());
      _m4.compose(_bz.set(x, y, z), _q, _s.set(s, s * thick, s));
      return layer.putM(owner, _m4, rbcCol);
    }

    function drawProplatelets(c, bc) {
      const shaftCol = col2.copy(bc).lerp(white, 0.2);
      for (const pp of c.pp) {
        const segs = 6;
        proplateletPoint(c, pp, 0, _pa);
        for (let k = 1; k <= segs; k++) {
          proplateletPoint(c, pp, k / segs, _tip);
          const dx = _tip[0] - _pa[0], dy = _tip[1] - _pa[1], dz = _tip[2] - _pa[2];
          const len = Math.hypot(dx, dy, dz) || 1e-3;
          _ax[0] = dx / len; _ax[1] = dy / len; _ax[2] = dz / len;
          L.shaft.putO(c, _pa, _ax, len, 0.6, 0.6, shaftCol);
          // beads (platelet-sized swellings) along the shaft; the tip bead is the next platelet
          const br = k === segs ? 1.5 : 1.05;
          L.plt.putO(c, _tip, _ax, br * 1.2, br, br, pltColor(col));
          _pa[0] = _tip[0]; _pa[1] = _tip[1]; _pa[2] = _tip[2];
        }
      }
    }

    function pltColor(out) {
      if (Q.mode === 'histology') return out.setRGB(0.62, 0.45, 0.85);
      if (Q.mode === 'confocal') return out.setRGB(0.55, 0.95, 1);
      return out.setHex(HEX.PLT);
    }

    const _n = [0, 1, 0];
    function drawFlow() {
      for (const e of flow) {
        if (Math.abs(e.z - params.focusZ) > Q.slab) continue;
        if (e.kind === 'cell') { if (params.showBlood) drawCell(e.agent); continue; }
        if (!params.showBlood && e.bg) continue;
        // tumbling in the shear plane spanned by the flow (x) and the radial direction
        const cy = Math.cos(e.th), cz = Math.sin(e.th);
        const cs = Math.cos(e.spin), sn = Math.sin(e.spin);
        _n[0] = sn; _n[1] = cs * cy + e.tilt * -cz; _n[2] = cs * cz + e.tilt * cy;
        if (e.kind === 'rbc') drawDisc(L.rbc, e, e.x, e.y, e.z, _n, 1, 1);
        else {
          _v.set(_n[0], _n[1], _n[2]).normalize();
          _q.setFromUnitVectors(_up, _v);
          _m4.compose(_bz.set(e.x, e.y, e.z), _q, _s.set(1.3, 0.45, 1.3));
          L.plt.putM(e, _m4, pltColor(col));
        }
      }
    }

    function draw() {
      rbcCol.setHex(Q.mode === 'physical' ? HEX.RBC : 0xffffff);
      reticCol.setHex(Q.mode === 'physical' ? 0xc4454f : 0xffffff);
      Q.slab = sectioning() ? OPTICAL_HALF : Infinity;
      for (const l of ALL_LAYERS) l.begin();
      for (const c of cells) {
        if (Math.abs(c.z - params.focusZ) > Q.slab + (c.type === 'MK' ? c.R * 0.6 : 0)) continue;
        if (c.phase === 'enuc') drawEnucleating(c);
        else drawCell(c);
      }
      drawFlow();
      for (const l of ALL_LAYERS) l.end();
    }

    // ============================================================ SSA + signals
    function recomputeSignals() {
      sig.S5 = stat5Drive(params.EPO);
      sig.E1 = mapk.inputFromLigand(params.SCF);
      sig.erk = mapk.nuclearERK(params.SCF);
      sig.uG = GAIN_G * sig.S5;
      sig.uP = (GAIN_P * params.GM) / (params.GM + 1);
      tp.uG = sig.uG; tp.uP = sig.uP;
      sig.fate = fateDistribution({ ...toggle.defaults, uG: sig.uG, uP: sig.uP }, { n: 240, rng: fateRng, sigma: params.sigma });
      ssa?.setEnv({ selfRenewal: params.a, pMEP: sig.fate.MEP, gmp: sig.fate.gmp, erk: sig.erk, stat5: sig.S5 });
      sig.dirty = false;
      updateSignalViews();
    }

    // ============================================================ panel
    panel.section('Signals & niche');
    const markDirty = () => { sig.dirty = true; sig.lastChange = performance.now(); };
    const epoSlider = panel.slider({ label: 'EPO (U/ml) → JAK2/STAT5', min: 0.01, max: 10, value: params.EPO, log: true, onChange: (v) => { params.EPO = v; markDirty(); } });
    panel.slider({ label: 'GM-CSF (ng/ml) → PU.1 bias', min: 0.01, max: 100, value: params.GM, log: true, onChange: (v) => { params.GM = v; markDirty(); } });
    panel.slider({ label: 'SCF / c-Kit ligand (nM) → MAPK/ERK', min: 0.03, max: 10, value: params.SCF, log: true, onChange: (v) => { params.SCF = v; markDirty(); } });
    panel.slider({ label: 'Transcriptional noise σ', min: 0, max: 0.3, value: params.sigma, onChange: (v) => { params.sigma = v; markDirty(); } });
    panel.slider({ label: 'HSC self-renewal probability a', min: 0.3, max: 1, value: params.a, onChange: (v) => { params.a = v; markDirty(); } });
    const sigOut = panel.readouts(['STAT5 drive (EPO)', 'ERK activity (SCF)', 'P(CMP → MEP)', 'GMP → neu / mono', 'HSC niches filled']);
    panel.note('HSC steady state H* = K(1 − 1/2a): below a = ½ the stem-cell pool is exhausted.');

    panel.section('View');
    panel.select({
      label: 'Camera', value: 'overview',
      options: [
        { value: 'overview', label: 'Overview: bone, marrow cords, sinusoid' },
        { value: 'island', label: 'Erythroblastic island (enucleation)' },
        { value: 'mk', label: 'Megakaryocyte shedding platelets into the sinusoid' },
        { value: 'niche', label: 'Endosteal HSC niche' },
      ],
      onChange: (v) => focusCamera(v),
    });
    panel.select({
      label: 'Colour cells by', value: params.colorBy,
      options: [
        { value: 'lineage', label: 'Lineage & maturation (stain colours in H&E)' },
        { value: 'reporter', label: 'GATA1-mCherry / PU.1-eYFP reporters' },
      ],
      onChange: (v) => { params.colorBy = v; updateLegend(); },
    });
    panel.toggle({ label: 'Blood in the sinusoid', value: params.showBlood, onChange: (v) => { params.showBlood = v; } });
    const sectionSel = panel.select({
      label: `Optical section (${2 * OPTICAL_HALF} µm slab of cells around the focal plane)`, value: params.section,
      options: [
        { value: 'confocal', label: 'In confocal mode only' },
        { value: 'always', label: 'In every mode' },
        { value: 'off', label: 'Off (whole slab of marrow)' },
      ],
      onChange: (v) => { params.section = v; applyVisibility(); },
    });
    const focusSlider = panel.slider({ label: 'Focal plane depth z (µm)', min: -28, max: 28, step: 0.5, value: params.focusZ, onChange: (v) => { params.focusZ = v; applyVisibility(); } });
    panel.toggle({ label: 'Stroma, adipocytes & bone lining', value: params.showStroma, onChange: (v) => { params.showStroma = v; applyVisibility(); } });
    panel.buttons([
      { label: 'Restart marrow', primary: true, onClick: () => { reset(); stage.simTime = 0; } },
      { label: 'Anaemia (EPO ×20)', onClick: () => { params.EPO = Math.min(params.EPO * 20, 10); epoSlider.set(params.EPO); markDirty(); } },
    ]);

    panel.section('Fate decision · GATA1–PU.1 toggle');
    panel.equation('dG/dt = a·Gⁿ/(θⁿ+Gⁿ) + b·θⁿ/(θⁿ+Pⁿ) − G + u_EPO + σξ\ndP/dt = a·Pⁿ/(θⁿ+Pⁿ) + b·θⁿ/(θⁿ+Gⁿ) − P + u_GM + σξ');
    const scatter = panel.chart({
      title: 'GATA1 (y) vs PU.1 (x) in every nucleated cell, live · dashed: G = P', xLabel: '', yRange: [0, 2.3], height: 170, capacity: 2000,
      series: [
        { name: 'HSC/CMP', color: '#5ee0c1', points: true, width: 0.01 },
        { name: 'ery/MK', color: '#ff5a6a', points: true, width: 0.01 },
        { name: 'myeloid', color: '#ffb454', points: true, width: 0.01 },
        { name: 'lymph', color: '#58a8ff', points: true, width: 0.01 },
        { name: '', color: '#4a5a70', dash: [3, 3] },
      ],
    });
    const togOut = panel.readouts(['attractors (a = 1)', 'attractors (a = 0)', 'pitchfork a_c']);
    {
      const s1 = toggle.fixedPoints({ ...toggle.defaults }).filter((f) => f.stable).length;
      const s0 = toggle.fixedPoints({ ...toggle.defaults, a1: 0, a2: 0 }).filter((f) => f.stable).length;
      togOut.set('attractors (a = 1)', `${s1} (primed centre + 2 fates)`);
      togOut.set('attractors (a = 0)', `${s0} (bistable switch)`);
      togOut.set('pitchfork a_c', toggle.criticalSelfActivation().toFixed(3));
    }
    panel.note('Huang, Guo, May &amp; Enver 2007 (Dev Biol 305:695): n = 4, θ = 0.5, b = 1. Self-activation a decays in CMPs (a = e<sup>−t/2</sup>); signal inputs as in Chickarmane, Enver &amp; Peterson 2009. The toggle is a model hypothesis: live imaging of PU.1-eYFP / GATA1-mCherry reporters (Hoppe et al. 2016, Nature 535:299) found that early myeloid choice is not initiated by random PU.1:GATA1 ratio changes.');

    panel.section('MAPK/ERK · SCF → c-Kit → Raf → MEK → ERK');
    panel.equation('v = k_cat·E·(S/K_m) / (1 + Σ S_j/K_m)   (competitive Michaelis–Menten)\nn_H = ln 81 / ln(EC90/EC10)');
    const mapkChart = panel.chart({
      title: 'Steady state vs log10 Ras-GTP input E1 (µM); dashed = current SCF', xLabel: '', yRange: [0, 1.05],
      series: [{ name: 'Raf*', color: '#8a98ab' }, { name: 'MEK-PP', color: '#c38bff' }, { name: 'ERK-PP', color: '#5ee0c1' }],
    });
    const hill = mapk.cascadeHill();
    {
      const xs = [], a = [], b = [], c = [];
      for (let i = 0; i <= 90; i++) {
        const lx = -6.5 + (i / 90) * 4.5;
        xs.push(lx);
        const s = mapk.doseResponse([10 ** lx]);
        a.push(s.MKKK[0]); b.push(s.MKK[0]); c.push(s.MAPK[0]);
      }
      mapkChart.set(xs, [a, b, c]);
    }
    const mapkOut = panel.readouts(['n_H Raf* / MEK-PP / ERK-PP', 'E1 from current SCF']);
    mapkOut.set('n_H Raf* / MEK-PP / ERK-PP', `${hill.MKKK.nH.toFixed(2)} / ${hill.MKK.nH.toFixed(2)} / ${hill.MAPK.nH.toFixed(1)}`);
    panel.note('Huang &amp; Ferrell 1996 (PNAS 93:10078) totals and constants. Their mass-action model, which keeps enzyme–substrate complexes, gave n<sub>H</sub> ≈ 1.0 / 1.7 / 4.9; this Michaelis–Menten reduction neglects sequestration and is steeper, with the same ordering.');

    panel.section('JAK2/STAT5 · EPO → EpoR');
    panel.equation('STAT5 →(pEpoR) pSTAT5;  2 pSTAT5 → dimer → nucleus → 2 STAT5 → export\nCIS, SOCS ∝ nuclear pSTAT5 (negative feedback)');
    const jakChart = panel.chart({
      title: 'EPO step at the current dose, 0–120 min: pEpoR, pSTAT5 cytoplasm / nucleus (×4, dashed = CIS/SOCS knocked out), CIS', xLabel: '', yRange: [0, 1],
      series: [
        { name: 'pEpoR', color: '#8a98ab' }, { name: 'cyto', color: '#ffb454' },
        { name: 'nuc', color: '#ff5a6a' }, { name: 'no FB', color: '#ff5a6a', dash: [4, 3] },
        { name: 'CIS', color: '#5ee0c1' },
      ],
    });
    panel.note('After Swameye et al. 2003 (PNAS 100:1028) nucleocytoplasmic cycling, with CIS/SOCS feedback as in Vera et al. 2008 (BMC Syst Biol 2:38). Total STAT5 is conserved.');

    panel.section('Marrow output · Gillespie compartment model');
    const popChart = panel.chart({
      title: `Whole-marrow compartments (log) vs days, K = ${SSA_DEFAULTS.K} niches`, logY: true, xLabel: '', capacity: 500, height: 160,
      series: [
        { name: 'HSC', color: '#5ee0c1' }, { name: 'prog', color: '#86b6ff' }, { name: 'ery', color: '#c05a8a' },
        { name: 'RBC', color: '#ff4a4a' }, { name: 'PLT', color: '#f6b8ec' }, { name: 'WBC', color: '#ffb454' },
      ],
    });
    const popOut = panel.readouts(['HSC (SSA) vs H*', 'M : E ratio (marrow)', 'released into sinusoid (3D)']);

    panel.section('Red-cell geometry vs a real blood smear');
    const rbcChart = panel.chart({
      title: 'Radial profile vs r/R: smear absorbance vs Evans–Fung thickness', xLabel: '', yRange: [0, 1.1], height: 140,
      series: [{ name: 'measured (Chula smear)', color: '#ff6fae', points: true }, { name: 'Evans–Fung model', color: '#5ee0c1' }],
    });
    if (rbcData?.normal) {
      const nd = rbcData.normal;
      const EF = (x) => (x >= 1 ? 0 : Math.sqrt(1 - x * x) * (0.207 + 2.003 * x * x - 1.123 * x ** 4));
      const efMax = Math.max(...nd.r_over_R.map((x) => EF(Math.min(x, 1))));
      rbcChart.set(nd.r_over_R, [nd.absorbance_profile, nd.r_over_R.map((x) => EF(x) / efMax)]);
      panel.table(['quantity', 'value'], [
        ['cells measured', nd.n_cells],
        ['profile correlation', nd.profile_correlation],
        ['central pallor (min / max)', nd.pallor_min_over_max],
        ['lymphocyte nucleus Ø (µm)', rbcData.lymphocyte.lymphocyte_nucleus_diameter_um],
      ]);
      panel.note(`${rbcData.source}: RBC diameter used as a 7.82 µm ruler; the same Evans–Fung profile builds every 3D erythrocyte, and the measured lymphocyte nucleus sets the lymphocyte agents.`);
    } else panel.note('rbc.json could not be loaded; the Evans–Fung geometry is still used.');

    function updateSignalViews() {
      sigOut.set('STAT5 drive (EPO)', sig.S5.toFixed(2));
      sigOut.set('ERK activity (SCF)', sig.erk.toFixed(2));
      sigOut.set('P(CMP → MEP)', sig.fate.MEP.toFixed(2));
      sigOut.set('GMP → neu / mono', `${sig.fate.gmp.Neu.toFixed(2)} / ${sig.fate.gmp.Mono.toFixed(2)}`);
      mapkChart.markers = [{ x: Math.log10(Math.max(sig.E1, 1e-9)), color: '#ffb454' }];
      mapkChart.dirty = true;
      mapkOut.set('E1 from current SCF', `${fmt(sig.E1)} µM`);
      const p = { ...jakstat.defaults, EPO: params.EPO };
      const a = timeCourse(p, 120, { every: 1 }), b = timeCourse({ ...p, fb: 0 }, 120, { every: 1 });
      const cisMax = Math.max(...a.y.map((y) => y[7]), 1e-9);
      jakChart.set(a.t, [
        a.y.map((y) => y[0]), a.y.map((y) => y[2] + 2 * y[3]), a.y.map((y) => 8 * y[4]),
        b.y.map((y) => 8 * y[4]), a.y.map((y) => y[7] / cisMax),
      ]);
    }

    function updateScatter() {
      const pts = [];
      for (const c of cells) if (!ANUCLEATE.has(c.type) && c.phase !== 'apop') pts.push(c);
      for (const e of flow) if (e.kind === 'cell') pts.push(e.agent);
      pts.sort((a, b) => a.u - b.u);
      const xs = [0], cols = [[NaN], [NaN], [NaN], [NaN], [0]];
      const grpIdx = (c) => {
        const g = TYPE_GROUP[c.type];
        if (c.type === 'HSC' || c.type === 'CMP') return 0;
        if (g === 'ery' || g === 'mk' || c.type === 'MEP') return 1;
        if (g === 'lymph' || c.type === 'CLP') return 3;
        return 2;
      };
      for (const c of pts) {
        const k = grpIdx(c);
        xs.push(Math.min(c.u, 2.3));
        for (let s = 0; s < 4; s++) cols[s].push(s === k ? Math.min(c.g, 2.3) : NaN);
        cols[4].push(Math.min(c.u, 2.3));
      }
      xs.push(2.3);
      for (let s = 0; s < 4; s++) cols[s].push(NaN);
      cols[4].push(2.3);
      scatter.set(xs, cols);
    }

    // ============================================================ legend, visibility, picking
    function updateLegend() {
      if (Q.mode === 'histology') {
        legend([
          { color: '#4b3aa0', label: 'haematoxylin · chromatin (dense: lymphocytes, late erythroblasts)' },
          { color: '#9486e0', label: 'basophilic cytoplasm · blasts, proerythroblasts' },
          { color: '#f0306a', label: 'eosin · haemoglobin (reticulocytes, erythrocytes)' },
          { color: '#f07aa8', label: 'eosin · bone matrix, granulocyte cytoplasm, MK' },
          { color: '#f7eef2', label: 'adipocyte · lipid dissolved (empty vacuole)' },
        ]);
        return;
      }
      if (Q.mode === 'confocal' && params.colorBy !== 'reporter') {
        legend([
          { color: '#4f86ff', label: 'DNA (Hoechst) · nuclei' },
          { color: '#3f63ff', label: 'bone collagen (second harmonic)' },
          { color: '#36f2d2', label: 'sinusoidal endothelium' },
          { color: '#ff5533', label: 'red cells (TER-119)' },
          { color: '#8cf2ff', label: 'platelets & proplatelets (CD41)' },
          { color: css(HEX.HSC), label: 'HSC (bright) · membrane rims by lineage as in 3D' },
          { color: '#9dff5a', label: 'island macrophages (F4/80)' },
        ]);
        return;
      }
      if (params.colorBy === 'reporter') {
        legend([
          { color: '#ff3a5a', label: 'GATA1-mCherry high (erythroid / MK)' },
          { color: '#d8ff3a', label: 'PU.1-eYFP high (myeloid)' },
          { color: '#f0c070', label: 'both (primed HSC / CMP)' },
          { color: '#2a2e38', label: 'neither (lymphoid)' },
          { color: css(HEX.RBC), label: 'anucleate red cells (no reporter)' },
        ]);
        return;
      }
      legend([
        { color: css(HEX.HSC), label: 'HSC (endosteal / perivascular niche)' },
        { color: css(HEX.CMP), label: 'CMP · CLP · MEP · GMP progenitors' },
        { color: css(HEX.ProEB), label: 'proerythroblast → erythroblast (haemoglobinising)' },
        { color: css(HEX.RBC), label: 'reticulocyte → erythrocyte (biconcave)' },
        { color: css(HEX.MK), label: 'megakaryocyte · proplatelets · platelets' },
        { color: css(HEX.Neu), label: 'neutrophil · eosinophil · basophil' },
        { color: css(HEX.Mono), label: 'monocyte' },
        { color: css(HEX.B), label: 'B · T-precursor · NK lymphocytes' },
        { color: '#e8dcc2', label: 'bone · endothelium · stroma · adipocytes' },
      ]);
    }

    const OVERVIEW = { c: [3, 0, 0], r: 98, dir: new THREE.Vector3(0.16, 0.5, 1) };
    function focusCamera(v) {
      let c = OVERVIEW.c, r = OVERVIEW.r, dir = OVERVIEW.dir;
      const setFocus = (z) => {
        params.focusZ = THREE.MathUtils.clamp(Math.round(z * 2) / 2, -28, 28);
        focusSlider.set(params.focusZ);
        if (v === 'overview') { params.section = 'confocal'; } else if (params.section === 'confocal') params.section = 'always';
        sectionSel.set(params.section);
        applyVisibility();
      };
      if (v === 'island') {
        const m = macrophages[0];
        c = [m.x, m.y + 3, m.z]; r = 30; dir = new THREE.Vector3(0.15, 0.6, 1);
      } else if (v === 'mk') {
        const mk = cells.filter((o) => o.type === 'MK').sort((a, b) => maturity(b) - maturity(a))[0];
        if (mk) { c = [mk.x + 6, (mk.y + VES.y) / 2, mk.z * 0.7 + VES.z * 0.3]; r = 38; dir = new THREE.Vector3(0.2, 0.35, 1); }
      } else if (v === 'niche') {
        const sl = slots.find((q) => q.occ && q.kind === 'endosteal' && q.z > 0) ?? slots.find((q) => q.occ) ?? slots[0];
        c = [sl.x, sl.y + 6, sl.z]; r = 34; dir = new THREE.Vector3(0.12, 0.5, 1);
      }
      setFocus(v === 'overview' ? 2 : c[2]);
      stage.frame(c, r, dir);
      stage.clipAxis.set(0, 0, 1);
      stage.clipRange = [-32 - c[2], 32 - c[2]];
    }

    function sectioning() { return params.section === 'always' || (params.section === 'confocal' && Q.mode === 'confocal'); }
    function applyVisibility() {
      const sect = sectioning();
      staticBody.mesh.visible = staticShaft.mesh.visible = params.showStroma && !sect;
      for (const a of adipocytes) a.mesh.visible = params.showStroma && !(sect && Math.abs(a.z - params.focusZ) > OPTICAL_HALF + a.R * 0.5);
    }

    const STAGE_NAMES = {
      EB: ['basophilic erythroblast', 'polychromatic erythroblast', 'orthochromatic erythroblast'],
      Neu: ['myeloblast / promyelocyte', 'myelocyte', 'metamyelocyte', 'band neutrophil', 'segmented neutrophil'],
    };
    function stageName(c) {
      const p = maturity(c);
      if (c.type === 'EB') return STAGE_NAMES.EB[Math.min(2, Math.floor(p * 3))];
      if (c.type === 'Neu') return STAGE_NAMES.Neu[p < 0.3 ? 0 : p < 0.55 ? 1 : p < 0.75 ? 2 : p < 0.88 ? 3 : 4];
      if (c.type === 'MK') return `${Math.round(2 ** (1 + 4 * Math.min(p / 0.75, 1)))}N${p >= 0.75 ? ' · shedding platelets' : ' · endomitosis'}`;
      return '';
    }
    function cellInfo(c) {
      if (!c) return null;
      if (c.type === 'Pyr') return 'Pyrenocyte · extruded erythroblast nucleus\nexposes phosphatidylserine and is engulfed by the island macrophage';
      const node = LINEAGE[c.type];
      const ageD = (simH - c.born) / 24;
      const st = stageName(c);
      let s = `${node.name}${c.phase === 'enuc' ? ' · ENUCLEATING' : ''}${c.phase === 'apop' ? ' · apoptotic (growth-factor withdrawal)' : ''}${c.phase === 'dying' ? ' · spent (bare nucleus, phagocytosed)' : ''}${c.mit >= 0 ? ' · dividing' : ''}\n`;
      if (st) s += `${st}\n`;
      if (ANUCLEATE.has(c.type)) s += 'anucleate: no GATA1 / PU.1 transcription\n';
      else s += `GATA1 ${c.g.toFixed(2)} · PU.1 ${c.u.toFixed(2)}${c.type === 'HSC' ? ' (primed)' : ''}\n`;
      s += `age ${ageD.toFixed(1)} d · generation ${c.gen}`;
      if (c.stageD > 0 && !ANUCLEATE.has(c.type)) s += ` · stage ${Math.round(maturity(c) * 100)}%`;
      if (c.type === 'HSC' && c.slot >= 0) s += `\n${slots[c.slot].kind} niche slot`;
      if (c.egress) s += '\n→ crossing into the sinusoid';
      s += `\n${node.markers}`;
      return s;
    }
    const pickCell = (layer) => ({
      get object() { return layer.mesh; },
      info: (hit) => {
        const o = layer.owners[hit.instanceId];
        if (!o) return null;
        if (o.kind === 'rbc') return `Erythrocyte in sinusoidal blood\nEvans–Fung biconcave disc, Ø 7.82 µm\n${o.bg ? 'circulating' : 'released from this marrow'}`;
        if (o.kind === 'plt') return `Platelet · anucleate disc ≈ 2.5 µm\n${o.bg ? 'circulating' : 'just shed from a proplatelet'}`;
        return cellInfo(o);
      },
    });
    stage.setPickables([
      ...['hsc', 'blast', 'gran', 'eryEarly', 'eryLate', 'mkBody', 'retic', 'rbc', 'plt', 'shaft'].map((k) => pickCell(L[k])),
      ...macrophages.map((m) => ({ object: m.mesh, info: () => `Central macrophage of an erythroblastic island (${m.id})\nanchors erythroblasts, engulfs extruded nuclei` })),
      ...adipocytes.map((a) => ({ object: a.mesh, info: () => (params.showStroma ? a.mesh.userData.info : null) })),
      { object: vessel, info: () => `Sinusoid · fenestrated endothelium, Ø ${2 * VES.r} µm\nflow shown at display speed` },
      { object: bone, info: () => 'Trabecular bone · endosteal surface lined by osteoblasts\nosteocytes in lacunae (visible in section)' },
    ]);

    // ============================================================ lifecycle
    let sampleAt = 0, scatterAt = 0, recAt = 0;

    function advance(dtVis) {
      const dtH = dtVis * HOURS_PER_SECOND;
      simH += dtH;
      biology(dtH);
      const sub = Math.max(1, Math.ceil(dtVis / 0.05));
      for (let k = 0; k < sub; k++) mechanics(dtVis / sub);
      egressAndFlow(dtVis);
    }

    function reset() {
      ssa = new LineageSSA({ rng: new RNG(777), env: { selfRenewal: params.a, pMEP: sig.fate?.MEP ?? 0.5, gmp: sig.fate?.gmp ?? SSA_DEFAULTS.gmp, erk: sig.erk, stat5: sig.S5 } });
      seedPopulation();
      // relax the warm start (positions, TF states, a few divisions) before the first frame
      for (let i = 0; i < 90; i++) advance(1 / 15);
      simH = 0;
      for (const c of cells) c.born = -rng.uniform(0, Math.max(c.stageT, 4));
      released = { RBC: 0, PLT: 0, WBC: 0 };
      popChart.clear();
      sampleAt = 0;
      updateScatter();
    }

    recomputeSignals();
    reset();
    applyVisibility();
    updateLegend();
    draw();
    focusCamera('overview');

    return {
      update(dt) {
        if (sig.dirty && performance.now() - sig.lastChange > 120) recomputeSignals();
        if (dt > 0) {
          for (let rem = dt; rem > 1e-6; rem -= 0.08) advance(Math.min(rem, 0.08));
          ssa.step((dt * HOURS_PER_SECOND) / 24, 20000);
          const tD = simH / 24;
          if (tD >= sampleAt) {
            sampleAt = tD + 0.1;
            const g = ssa.groups();
            popChart.push(ssa.t, [g.HSC, g.progenitors, g.erythroid, g.RBC, g.PLT, g.myeloid + g.lymphoid].map((v) => Math.max(v, 0.5)));
            const Hs = meanFieldSteadyState(ssa.env).HSC;
            popOut.set('HSC (SSA) vs H*', `${g.HSC} vs ${Hs.toFixed(0)}`);
            // M:E ratio as on a marrow aspirate: nucleated myeloid vs nucleated erythroid precursors in the marrow
            const my = ssa.count('GMP') + ['Neu', 'Mono', 'Eos', 'Baso'].reduce((a, id) => a + ssa.count(id) * marrowFraction(id), 0);
            const ery = ssa.count('ProEB') + ssa.count('EB');
            popOut.set('M : E ratio (marrow)', ery > 0 ? `${(my / ery).toFixed(1)} : 1` : '∞');
            popOut.set('released into sinusoid (3D)', `${released.RBC} RBC · ${released.PLT} PLT · ${released.WBC} WBC`);
          }
        }
        const now = performance.now();
        if (now - scatterAt > 330) { scatterAt = now; updateScatter(); }
        if (now - recAt > 500) {
          recAt = now;
          const H = cells.reduce((n, c) => n + (c.type === 'HSC' ? 1 : 0), 0);
          sigOut.set('HSC niches filled', `${H} / ${K3D}`);
        }
        draw();
        const d = Math.floor(simH / 24), h = Math.floor(simH % 24);
        const nH = slots.reduce((n, s) => n + (s.occ ? 1 : 0), 0);
        setStatus(`t = ${d} d ${String(h).padStart(2, '0')} h · ${cells.length} marrow cells · ${nH} HSC · ${released.RBC} RBC & ${released.PLT} platelets released`);
      },
      reset() { reset(); },
      onMode(mode) { Q.mode = mode; applyVisibility(); updateLegend(); },
      dispose() {
        for (const g of [G.lobed4, G.lobed3, G.lobed2, G.indent, G.mkNuc, G.mkBody, G.retic, G.rbc, G.mac, G.shaft, sphere, sphereLo, sphereMid]) g.dispose();
      },
    };
  },
};
