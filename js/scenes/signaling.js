// Stage 1 · Notch, Wnt & Hedgehog.
//
// (a) "Inside one stem cell": a cut-away receiver cell touching a Delta-
//     presenting sender cell, with a primary cilium. Three pathway models are
//     integrated live - canonicalWnt (wnt.js), singleCell Notch (notch.js) and
//     the Shh -> PTCH1 -> SMO -> GLI model (hedgehog.js) - and every receptor,
//     ligand and messenger drawn is a *sample of their state*: the number of
//     visible particles of a species is proportional to its concentration
//     (smooth add/remove), bound receptors follow the bound fractions, and
//     particle turnover events follow the model fluxes. The paper's own linear
//     Wnt cascade (paperWnt) runs alongside in the charts.
// (b) "Stem-cell spheroid": ~300 packed cells; Delta-Notch lateral inhibition
//     (Collier et al. 1996) on the 3D contact graph, and two opposing morphogen
//     gradients (Wnt bead / Shh bead) from the paper's reaction-diffusion
//     equation, read out per cell by the canonical Wnt and Hedgehog models.
import * as THREE from 'three';
import { cellMaterial } from '../engine/materials.js';
import { InstancePool } from '../engine/instances.js';
import { unitSphere, blobGeometry, mergeGeometries, mergeVertices } from '../engine/geometry.js';
import { fmt } from '../engine/chart.js';
import { rk4Step, rk4Work } from '../models/core/ode.js';
import { RNG } from '../models/core/rng.js';
import { paperWnt, canonicalWnt } from '../models/pathways/wnt.js';
import { singleCell, collier } from '../models/pathways/notch.js';
import { hedgehog, frenchFlag, NEURAL_TUBE_DOMAINS, NEURAL_TUBE_THRESHOLDS } from '../models/pathways/hedgehog.js';
import { Tissue, Cell } from '../models/morpho/agents3d.js';

const TIME_A = 3;            // model time units per animation second (single cell)
const TIME_B = 6;            // model time units per animation second (spheroid)
const FRAME_SECONDS = 0.04;  // animation seconds per frame of the NPC time-lapse (shape calibration only)
const TURNOVER = 0.1;        // fraction of the model turnover flux drawn as particle events
const GOLDEN = 2.399963229728653;

// ------------------------------------------------------------ single-cell layout (µm)
const R_CELL = 10, R_NUC = 5.5, R_CHROM = 3.75, R_SEND = 7.5, GAP = 1.7, R_EXTRA = 19;
const CONTACT = norm3([-0.72, -0.12, 0.68]);
const CILIUM = norm3([-0.3, 1, -0.3]);
const CILIUM_LEN = 6.5, CILIUM_R = 0.5;
const SENDER_C = CONTACT.map((v) => v * (R_CELL + R_SEND + GAP));

const HEX = {
  wnt: '#ff9a3c', fz: '#b8743c', fzOn: '#ffc46b', dvl: '#fff2a6', bcat: '#ffd84a', bcatP: '#8d6f45', dc: '#a08cff',
  delta: '#ff3b9d', notch: '#39ff88', nicd: '#b8ffd8', next: '#2a9d62',
  shh: '#45c8ff', ptch: '#3a74ff', smo: '#8a5cc0', smoOn: '#f2c6ff', gliA: '#5ef0ff', gliR: '#8878ff',
};
const COL = Object.fromEntries(Object.entries(HEX).map(([k, v]) => [k, new THREE.Color(v)]));

// particle states
const ST = { CYTO: 0, IMPORT: 1, PORE_IN: 2, NUC: 3, EXPORT: 4, PORE_OUT: 5, DOOMED: 6, DYING: 7, EXTRA: 8, CILIUM: 9 };
const ST_TEXT = [
  'diffusing in the cytoplasm', 'heading to a nuclear pore (import)', 'translocating through a nuclear pore',
  'in the nucleus', 'heading to a nuclear pore (export)', 'being exported', 'captured by a destruction complex',
  'being degraded', 'extracellular', 'leaving the primary cilium',
];

// ------------------------------------------------------------ small math helpers
function norm3(a) { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; }
const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const clamp01 = (x) => (x < 0 ? 0 : x > 1 ? 1 : x);
function cross3(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
function basis(n) {
  const t = Math.abs(n[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0];
  const e1 = norm3(cross3(n, t));
  return [e1, cross3(n, e1)];
}
/** n points spread evenly over a spherical cap (or annulus th0..th1) around `axis`. */
function sunflower(axis, th1, n, { th0 = 0, phase = 0 } = {}) {
  const [e1, e2] = basis(axis);
  const out = [];
  for (let k = 0; k < n; k++) {
    const th = Math.sqrt(th0 * th0 + (th1 * th1 - th0 * th0) * ((k + 0.5) / n));
    const ph = k * GOLDEN + phase, s = Math.sin(th), c = Math.cos(th);
    out.push(norm3([0, 1, 2].map((j) => axis[j] * c + s * (Math.cos(ph) * e1[j] + Math.sin(ph) * e2[j]))));
  }
  return out;
}
function fibSphere(n) {
  const out = [];
  for (let i = 0; i < n; i++) {
    const y = 1 - (2 * (i + 0.5)) / n, r = Math.sqrt(1 - y * y), ph = i * GOLDEN;
    out.push([Math.cos(ph) * r, y, Math.sin(ph) * r]);
  }
  return out;
}
const inOctant = (d, m = 0.06) => d[0] > -m && d[1] > -m && d[2] > -m;

// ------------------------------------------------------------ geometry builders (+x = outward)
function alongX(g) { g.rotateZ(-Math.PI / 2); return g; }
function rodX(r, x0, x1, seg = 8) { const g = alongX(new THREE.CylinderGeometry(r, r, x1 - x0, seg, 1)); g.translate((x0 + x1) / 2, 0, 0); return g; }
function ballAt(r, x, y = 0, z = 0, detail = 1) { const g = new THREE.IcosahedronGeometry(r, detail); g.translate(x, y, z); return g; }
function merge(list) {
  return mergeGeometries(list.map((x) => (x.index ? x.toNonIndexed() : x)));
}
/** blobGeometry with shared vertices, so the lumpy surface is smooth-shaded. */
function smoothBlob(r, amp, seed, detail) {
  const g = blobGeometry(r, amp, seed, detail);
  g.deleteAttribute('uv'); g.deleteAttribute('normal');
  const m = mergeVertices(g);
  m.computeVertexNormals();
  return m;
}
function arcTube(radius, a, b, tube) {
  // quarter circle in the plane spanned by unit axes a, b (cut rim of the cut-away)
  class Arc extends THREE.Curve {
    getPoint(t, out = new THREE.Vector3()) {
      const th = t * Math.PI / 2;
      return out.set(...[0, 1, 2].map((j) => radius * (Math.cos(th) * a[j] + Math.sin(th) * b[j])));
    }
  }
  return new THREE.TubeGeometry(new Arc(), 48, tube, 6, false);
}

// ------------------------------------------------------------ particle swarm (structure of arrays)
class Swarm {
  constructor(cap, pool) {
    this.cap = cap; this.pool = pool; this.n = 0; this.acc = 0; this.acc2 = 0;
    const f = () => new Float32Array(cap);
    this.x = f(); this.y = f(); this.z = f(); this.vx = f(); this.vy = f(); this.vz = f(); this.s = f(); this.t = f();
    this.st = new Uint8Array(cap); this.aux = new Int16Array(cap);
  }

  add(x, y, z, st, s0 = 0) {
    if (this.n >= this.cap) return -1;
    const i = this.n++;
    this.x[i] = x; this.y[i] = y; this.z[i] = z;
    this.vx[i] = 0; this.vy[i] = 0; this.vz[i] = 0;
    this.s[i] = s0; this.t[i] = 0; this.st[i] = st; this.aux[i] = -1;
    return i;
  }

  kill(i) {
    const j = --this.n;
    if (i === j) return;
    this.x[i] = this.x[j]; this.y[i] = this.y[j]; this.z[i] = this.z[j];
    this.vx[i] = this.vx[j]; this.vy[i] = this.vy[j]; this.vz[i] = this.vz[j];
    this.s[i] = this.s[j]; this.t[i] = this.t[j]; this.st[i] = this.st[j]; this.aux[i] = this.aux[j];
  }

  pick(state, rng) {
    if (!this.n) return -1;
    const start = rng.int(this.n);
    for (let k = 0; k < this.n; k++) {
      const i = (start + k) % this.n;
      if (this.st[i] === state) return i;
    }
    return -1;
  }

  clear() { this.n = 0; this.acc = 0; this.acc2 = 0; }
}

/** Receptor/ligand slots whose visibility eases toward a target count. */
class SlotSet {
  constructor(n) { this.n = n; this.w = new Float32Array(n); this.on = new Uint8Array(n); }
  setCount(k) { for (let i = 0; i < this.n; i++) this.on[i] = i < k ? 1 : 0; }
  setCountFromEnd(k) { for (let i = 0; i < this.n; i++) this.on[i] = i >= this.n - k ? 1 : 0; }
  ease(dt, rate = 3) { const a = Math.min(1, rate * dt); for (let i = 0; i < this.n; i++) this.w[i] += (this.on[i] - this.w[i]) * a; }
  snap() { for (let i = 0; i < this.n; i++) this.w[i] = this.on[i]; }
}

export default {
  about: `
    <p><b>Inside one stem cell.</b> A cut-away pluripotent stem cell receives three
    first-order signals at once: <b>Wnt</b> (Frizzled/LRP6 → Dishevelled → inhibition of the
    APC/Axin/GSK-3β destruction complex → β-catenin stabilisation → nuclear TCF/LEF targets),
    <b>Notch</b> from the touching sender cell (Delta binding → S2/S3 cleavage → NICD → nucleus
    → Hes1) and <b>Hedgehog</b> at the primary cilium (Shh binds Patched → SMO released →
    GLI activator instead of GLI3 repressor). Every particle is a sample of the running ODEs:
    visible counts are proportional to concentrations, bound receptors follow bound
    fractions and turnover events follow the model fluxes.</p>
    <p><b>Stem-cell spheroid.</b> Lateral inhibition on the real 3D contact graph turns
    near-identical cells into a salt-and-pepper pattern, while opposing Wnt and Shh beads set up
    reaction–diffusion gradients that each cell reads out with the same pathway models.</p>
    <p class="note"><b>Established models:</b> the paper's linear Wnt cascade (exact, with its
    closed-form steady state); canonical Wnt with destruction complex and Axin2 feedback
    (reduced from Lee et al. 2003, PLoS Biol); receptor processing of Notch and Collier et al.
    1996 lateral inhibition (J Theor Biol); catalytic PTCH1 inhibition of SMO (Taipale et al. 2002,
    Nature), GLI activator/repressor promoter logic (Lai, Robertson &amp; Schaffer 2004, Biophys J;
    Saha &amp; Schaffer 2006, Development), Ptch1 feedback &amp; adaptation (Dessaud et al. 2007,
    Nature); diffusion–degradation gradients (the paper's ∂u/∂t = D∇²u + f(u) with f = −ku).
    <b>Phenomenological:</b> particle positions and motions, the drawing scale (particles per
    unit), the drawn fraction of turnover, the nuclear/cytoplasmic split of GLI, the Wnt → Delta
    (Jagged1) crosstalk strength and the Wnt × Notch fate map. All rates are nondimensional
    (time in model units τ).</p>`,
  paperRef: 'Paper: "Stage 1: Basic Stem Cell Reactions" → Notch, Wnt and Hedgehog signalling pathways (~500 reactions); "Mathematical Models – Example: Wnt/β-Catenin Pathway" (dW/dt = k1 − k2W …); "Reaction-diffusion: ∂u/∂t = D<sub>u</sub>∇²u + f(u,v)"; order-4 crosstalk (Wnt–Notch) and morphogen gradients in the spheroid view.',

  async create(ctx) {
    const { stage, root, panel, assets, legend, setStatus } = ctx;
    let npc = null;
    try { npc = await assets.loadJSON('npc_kinetics.json'); } catch { npc = null; }
    const tauFrames = npc?.fit?.tau_frames ?? 13.5;
    const IMPORT_TAU = tauFrames * FRAME_SECONDS; // s of animation for the approach to a pore
    const DWELL = 0.35;                           // s spent translocating through the pore

    const rng = new RNG(20240501);
    const groupA = new THREE.Group(), groupB = new THREE.Group();
    root.add(groupA, groupB);
    const params = { view: ctx.query?.get('view') === 'spheroid' ? 'sph' : 'cell', colorBy: 'fate' };

    // ================================================================ (a) single cell
    const octPlanes = [
      new THREE.Plane(new THREE.Vector3(-1, 0, 0), 0),
      new THREE.Plane(new THREE.Vector3(0, -1, 0), 0),
      new THREE.Plane(new THREE.Vector3(0, 0, -1), 0),
    ];
    const cutAway = (mat) => { mat.clippingPlanes = octPlanes; mat.clipIntersection = true; return mat; };

    // membranes, nucleus, cilium
    const membraneMat = cutAway(cellMaterial({ role: 'membrane', color: 0x8fd3ff, opacity: 0.1, gain: 0.28 }));
    const membrane = new THREE.Mesh(new THREE.SphereGeometry(R_CELL, 96, 64), membraneMat);
    membrane.renderOrder = 4;
    const envMat = cutAway(cellMaterial({ role: 'membrane', color: 0xa9b8ff, stain: 0x8a6cff, opacity: 0.12, useInst: 1, gain: 0.3 }));
    const envelope = new THREE.Mesh(new THREE.SphereGeometry(R_NUC, 72, 48), envMat);
    envelope.renderOrder = 3;
    const chromMat = cellMaterial({ role: 'nucleus', color: 0x98acff, stain: 0x3d6bff, noiseScale: 2.4, noise: 0.9, gain: 0.7 });
    const chromatin = new THREE.Mesh(smoothBlob(R_CHROM, 0.1, 3, 5), chromMat);
    const rimMat = cellMaterial({ role: 'solid', color: 0x9fd8ff, emissive: 0.25 });
    const rimNucMat = cellMaterial({ role: 'solid', color: 0xb9c4ff, emissive: 0.25 });
    const ex = [1, 0, 0], ey = [0, 1, 0], ez = [0, 0, 1];
    const rims = new THREE.Group();
    for (const [a, b] of [[ex, ey], [ey, ez], [ez, ex]]) {
      rims.add(new THREE.Mesh(arcTube(R_CELL, a, b, 0.09), rimMat));
      rims.add(new THREE.Mesh(arcTube(R_NUC, a, b, 0.07), rimNucMat));
    }
    const ciliumMat = cellMaterial({ role: 'membrane', color: 0xc9b6ff, stain: 0xc9b6ff, opacity: 0.22, gain: 0.6 });
    const cilium = new THREE.Mesh(new THREE.CapsuleGeometry(CILIUM_R, CILIUM_LEN, 6, 16), ciliumMat);
    cilium.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), new THREE.Vector3(...CILIUM));
    cilium.position.set(...CILIUM.map((v) => v * (R_CELL + CILIUM_LEN / 2 - 0.1)));
    cilium.renderOrder = 5;
    const basalBody = new THREE.Mesh(rodX(0.32, -0.6, 0.6, 12), cellMaterial({ role: 'solid', color: 0xd9ccff }));
    basalBody.quaternion.setFromUnitVectors(new THREE.Vector3(1, 0, 0), new THREE.Vector3(...CILIUM));
    basalBody.position.set(...CILIUM.map((v) => v * (R_CELL - 0.7)));

    // the neighbouring sender cell (presents Delta)
    const senderMat = cellMaterial({ role: 'membrane', color: 0xf6b3d6, opacity: 0.06, gain: 0.4 });
    const sender = new THREE.Mesh(new THREE.SphereGeometry(R_SEND, 64, 48), senderMat);
    sender.position.set(...SENDER_C);
    sender.renderOrder = 4;
    const senderNuc = new THREE.Mesh(smoothBlob(3.4, 0.1, 7, 4), cellMaterial({ role: 'nucleus', color: 0xd7b6ff, stain: 0x3d6bff, gain: 0.8 }));
    senderNuc.position.set(...SENDER_C.map((v, j) => v + CONTACT[j] * 1.6));
    groupA.add(membrane, envelope, chromatin, rims, cilium, basalBody, sender, senderNuc);

    // nuclear pores (static)
    const poreDirs = fibSphere(84).filter((d) => !inOctant(d, 0.12));
    const poreGeo = new THREE.TorusGeometry(0.4, 0.13, 6, 14);
    poreGeo.rotateY(Math.PI / 2);
    const porePool = new InstancePool(poreGeo, cellMaterial({ role: 'solid', color: 0xdfe6ff, emissive: 0.15 }), poreDirs.length, groupA);
    porePool.begin();
    for (const d of poreDirs) porePool.putOriented(d.map((v) => v * R_NUC), d, 1, 1, 1, 0xdfe6ff);
    porePool.end();

    // target-gene loci on the chromatin surface, facing the cut-away
    const LOCI = [
      { key: 'W', name: 'Wnt / TCF target (Axin2, Lgr5)', dir: norm3([0.85, 0.35, 0.55]), color: COL.bcat },
      { key: 'N', name: 'Notch / CSL target (Hes1)', dir: norm3([0.35, 0.75, 0.6]), color: COL.notch },
      { key: 'H', name: 'GLI target (Ptch1, Gli1, Nkx2.2)', dir: norm3([0.5, 0.45, 0.85]), color: COL.gliA },
    ];
    for (const L of LOCI) L.p = L.dir.map((v) => v * (R_CHROM * 1.1 + 0.1));
    const lociPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'reporter', emissive: 0.6 }), 3, groupA);

    // ---------------------------------------------------------------- receptor geometry
    const molMat = (o = {}) => cellMaterial({ role: 'molecule', emissive: 0.18, gain: 0.5, ...o });
    const fzGeo = merge([rodX(0.34, -0.45, 0.45, 10), ballAt(0.3, 0.72), (() => { const g = rodX(0.1, -0.4, 1.35, 6); g.translate(0, 0, 0.5); return g; })(), ballAt(0.24, 1.42, 0, 0.5)]);
    const notchGeo = merge([rodX(0.11, -0.5, 2.25, 6), ballAt(0.2, 0.35), ballAt(0.15, 1.2), ballAt(0.17, 2.25)]);
    const deltaGeo = merge([rodX(0.1, 0, 1.5, 6), ballAt(0.13, 0.8), ballAt(0.24, 1.55)]);
    const halfBridgeGeo = rodX(0.13, 0, 1, 8);
    const nextGeo = merge([rodX(0.12, -0.55, 0.25, 6), ballAt(0.14, -0.55)]);
    const ptchGeo = merge([rodX(0.42, -0.45, 0.45, 12), ballAt(0.3, 0.62, 0.15, 0), ballAt(0.26, 0.6, -0.2, 0.12)]);
    const smoGeo = merge([rodX(0.3, -0.45, 0.45, 10), ballAt(0.27, 0.72)]);
    const dcGeo = merge([ballAt(0.62, 0, 0, 0, 2), ballAt(0.45, 0.6, 0.32, 0, 2), ballAt(0.46, -0.25, 0.55, 0.38, 2), ballAt(0.38, 0.1, -0.5, 0.42, 2), ballAt(0.34, -0.45, -0.2, -0.45, 2)]);

    const fzPool = new InstancePool(fzGeo, molMat(), 32, groupA);
    const notchPool = new InstancePool(notchGeo, molMat(), 48, groupA);
    const deltaPool = new InstancePool(deltaGeo, molMat(), 64, groupA);
    const bridgePool = new InstancePool(halfBridgeGeo, molMat({ emissive: 0.35 }), 64, groupA);
    const nextPool = new InstancePool(nextGeo, molMat(), 32, groupA);
    const ptchPool = new InstancePool(ptchGeo, molMat(), 64, groupA);
    const smoPool = new InstancePool(smoGeo, molMat(), 32, groupA);
    const dvlPool = new InstancePool(unitSphere(1), molMat({ emissive: 0.4 }), 32, groupA);
    const dcPool = new InstancePool(dcGeo, molMat({ emissive: 0.1, noise: 0.3 }), 16, groupA);

    // ---------------------------------------------------------------- receptor slots
    const slotRng = new RNG(99);
    const shuffle = (arr) => { for (let i = arr.length - 1; i > 0; i--) { const j = slotRng.int(i + 1); [arr[i], arr[j]] = [arr[j], arr[i]]; } return arr; };
    const cosDeg = (d) => Math.cos(THREE.MathUtils.degToRad(d));
    const sphereDirs = fibSphere(640).filter((d) => !inOctant(d));
    const awayContact = (d) => dot3(d, CONTACT) < cosDeg(42);
    const awayCilium = (d) => dot3(d, CILIUM) < cosDeg(40);
    const mkSlot = (d, r = R_CELL) => ({ d, p: d.map((v) => v * r) });

    const fzSlots = shuffle(sphereDirs.filter((d) => awayContact(d) && awayCilium(d))).slice(0, 28).map((d) => mkSlot(d));
    const ptchMemSlots = sunflower(CILIUM, THREE.MathUtils.degToRad(34), 30, { th0: THREE.MathUtils.degToRad(9) })
      .filter((d) => !inOctant(d)).map((d) => mkSlot(d));
    const smoMemSlots = shuffle(sphereDirs.filter((d) => awayContact(d) && dot3(d, CILIUM) < cosDeg(36))).slice(0, 14).map((d) => mkSlot(d));
    const notchOtherSlots = shuffle(sphereDirs.filter((d) => awayContact(d) && awayCilium(d)
      && !fzSlots.some((s) => dot3(s.d, d) > 0.995))).slice(0, 10).map((d) => mkSlot(d));
    // contact zone: Delta:Notch bridges across the gap to the sender
    const bridgeSlots = sunflower(CONTACT, THREE.MathUtils.degToRad(19), 22).map((d) => {
      const p = d.map((v) => v * R_CELL);
      const w = p.map((v, j) => v - SENDER_C[j]);
      const b = dot3(w, CONTACT), disc = b * b - dot3(w, w) + R_SEND * R_SEND;
      const t = -b - Math.sqrt(Math.max(disc, 0));
      const q = p.map((v, j) => v + t * CONTACT[j]);
      return { d, p, q, half: t / 2 };
    });
    const notchFreeSlots = [
      ...sunflower(CONTACT, THREE.MathUtils.degToRad(36), 18, { th0: THREE.MathUtils.degToRad(22), phase: 0.7 }).map((d) => mkSlot(d)),
      ...notchOtherSlots,
    ];
    const nextSlots = sunflower(CONTACT, THREE.MathUtils.degToRad(21), 30, { phase: 1.3 }).map((d) => mkSlot(d));
    const senderAxis = CONTACT.map((v) => -v);
    const deltaSlots = sunflower(senderAxis, THREE.MathUtils.degToRad(52), 44, { phase: 0.4 })
      .map((d) => ({ d, p: d.map((v, j) => SENDER_C[j] + v * R_SEND) }));
    // primary-cilium slots (PTCH1 free, SMO active)
    const [ce1, ce2] = basis(CILIUM);
    const ciliumSlot = (k, n, phase) => {
      const f = 0.12 + 0.8 * (k + 0.5) / n, ph = k * 2.2 + phase;
      const out = [0, 1, 2].map((j) => Math.cos(ph) * ce1[j] + Math.sin(ph) * ce2[j]);
      const p = [0, 1, 2].map((j) => CILIUM[j] * (R_CELL + f * CILIUM_LEN) + out[j] * CILIUM_R);
      return { d: out, p };
    };
    const ptchCilSlots = Array.from({ length: 16 }, (_, k) => ciliumSlot(k, 16, 0));
    const smoCilSlots = Array.from({ length: 14 }, (_, k) => ciliumSlot(k, 14, 1.1));
    const ciliumTip = CILIUM.map((v) => v * (R_CELL + CILIUM_LEN * 0.95));
    const ciliumBase = CILIUM.map((v) => v * (R_CELL - 0.9));

    const S_fzBound = new SlotSet(fzSlots.length);
    const S_notchFree = new SlotSet(notchFreeSlots.length);
    const S_bridge = new SlotSet(bridgeSlots.length);
    const S_next = new SlotSet(nextSlots.length);
    const S_delta = new SlotSet(deltaSlots.length);
    const S_ptchCil = new SlotSet(ptchCilSlots.length);
    const S_ptchMemBound = new SlotSet(ptchMemSlots.length);
    const S_ptchMemFree = new SlotSet(ptchMemSlots.length);
    const S_smoCil = new SlotSet(smoCilSlots.length);
    const S_smoMem = new SlotSet(smoMemSlots.length);
    const allSlotSets = [S_fzBound, S_notchFree, S_bridge, S_next, S_delta, S_ptchCil, S_ptchMemBound, S_ptchMemFree, S_smoCil, S_smoMem];

    // destruction complexes: resting in the cytoplasm, recruited to LRP6 signalosomes by Dvl
    // (Bilic et al. 2007, Science 316:1619)
    const DC_MAX = 12;
    const dcRest = [], dcSig = [];
    for (let j = 0; j < DC_MAX; j++) {
      let p;
      do {
        const d = rng.onSphere();
        const r = THREE.MathUtils.lerp(R_NUC + 1.4, R_CELL - 1.6, rng.next());
        p = d.map((v) => v * r);
      } while (inOctant(p, 0.8) || dot3(norm3(p), CONTACT) > 0.85);
      dcRest.push(p);
      dcSig.push(fzSlots[(j * 5) % fzSlots.length].d.map((v) => v * (R_CELL - 1.15)));
    }
    const dcW = new Float32Array(DC_MAX), dcPulse = new Float32Array(DC_MAX), dcPos = dcRest.map((p) => p.slice());

    // ---------------------------------------------------------------- messengers & ligands
    const sphere2 = unitSphere(2);
    const octa = new THREE.OctahedronGeometry(1, 0);
    const mkPool = (geo, cap, o) => new InstancePool(geo, molMat(o), cap, groupA);
    const SW = {
      bcat: new Swarm(260, mkPool(sphere2, 128, { emissive: 0.3 })),
      nicd: new Swarm(160, mkPool(sphere2, 96, { emissive: 0.3 })),
      gliA: new Swarm(120, mkPool(sphere2, 64, { emissive: 0.3 })),
      gliR: new Swarm(80, mkPool(octa, 48, { emissive: 0.2 })),
      wnt: new Swarm(110, mkPool(sphere2, 128, { emissive: 0.25 })),
      shh: new Swarm(110, mkPool(sphere2, 96, { emissive: 0.25 })),
    };
    const SCALE = { bcat: 40, nicd: 22, gliA: 24, gliR: 24, wnt: 30 };
    const GLI_NUC = 0.65; // drawn nuclear share of GLI (the Hh ODE lumps compartments)

    // ---------------------------------------------------------------- models
    const pW = { ...canonicalWnt.defaults, wntIn: 0.1 };
    const pN = { ...singleCell.defaults, deltaTrans: 1 };
    const pH = { ...hedgehog.defaults, shh: 3 };
    const pP = { ...paperWnt.defaults };
    const WNT_REF = 0.1; // Wnt supply at which the paper cascade uses its published k1
    const fW = canonicalWnt.rhs(pW), fN = singleCell.rhs(pN), fH = hedgehog.rhs(pH), fP = paperWnt.rhs(pP);
    const wW = rk4Work(7), wN = rk4Work(7), wH = rk4Work(6), wP = rk4Work(4);
    let yW, yN, yH, yP, tA = 0, lastSampleA = -1;

    function resetModelsA() {
      yW = Float64Array.from(canonicalWnt.initial());
      yN = Float64Array.from(singleCell.initial());
      yH = Float64Array.from(hedgehog.initial());
      yP = new Float64Array(4);
      tA = 0; lastSampleA = -1;
    }

    function clampNonNeg(y) { for (let i = 0; i < y.length; i++) if (y[i] < 0) y[i] = 0; }
    function integrateA(dtm) {
      const n = Math.max(1, Math.ceil(dtm / 0.02)), h = dtm / n;
      pP.k1 = paperWnt.defaults.k1 * (pW.wntIn / WNT_REF);
      for (let k = 0; k < n; k++) {
        rk4Step(fW, tA, yW, h, wW);
        rk4Step(fN, tA, yN, h, wN);
        rk4Step(fH, tA, yH, h, wH);
        rk4Step(fP, tA, yP, h, wP);
        clampNonNeg(yW); clampNonNeg(yN); clampNonNeg(yH); clampNonNeg(yP);
        tA += h;
      }
    }

    // ---------------------------------------------------------------- particle mechanics
    const P3 = [0, 0, 0], A3 = [0, 0, 0];
    function ou(sw, i, dt, speed, tau) {
      const a = Math.exp(-dt / tau), b = speed * Math.sqrt(1 - a * a);
      sw.vx[i] = a * sw.vx[i] + b * rng.normal();
      sw.vy[i] = a * sw.vy[i] + b * rng.normal();
      sw.vz[i] = a * sw.vz[i] + b * rng.normal();
      sw.x[i] += sw.vx[i] * dt; sw.y[i] += sw.vy[i] * dt; sw.z[i] += sw.vz[i] * dt;
    }
    function shellClamp(sw, i, rmin, rmax, cx = 0, cy = 0, cz = 0) {
      const dx = sw.x[i] - cx, dy = sw.y[i] - cy, dz = sw.z[i] - cz;
      const r = Math.hypot(dx, dy, dz) || 1e-6;
      if (r >= rmin && r <= rmax) return;
      const target = r < rmin ? rmin : rmax, k = target / r;
      sw.x[i] = cx + dx * k; sw.y[i] = cy + dy * k; sw.z[i] = cz + dz * k;
      const vr = (sw.vx[i] * dx + sw.vy[i] * dy + sw.vz[i] * dz) / r;
      if ((r < rmin && vr < 0) || (r > rmax && vr > 0)) {
        sw.vx[i] -= 2 * vr * dx / r; sw.vy[i] -= 2 * vr * dy / r; sw.vz[i] -= 2 * vr * dz / r;
      }
    }
    function avoidOctant(sw, i, m = 0.35) {
      const x = sw.x[i], y = sw.y[i], z = sw.z[i];
      if (!(x > -m && y > -m && z > -m)) return;
      if (x <= y && x <= z) { sw.x[i] = -2 * m - x; sw.vx[i] = -Math.abs(sw.vx[i]); }
      else if (y <= z) { sw.y[i] = -2 * m - y; sw.vy[i] = -Math.abs(sw.vy[i]); }
      else { sw.z[i] = -2 * m - z; sw.vz[i] = -Math.abs(sw.vz[i]); }
    }
    function randomIn(rmin, rmax, avoidOct, out) {
      for (let k = 0; k < 40; k++) {
        rng.onSphere(out);
        const r = Math.cbrt(THREE.MathUtils.lerp(rmin ** 3, rmax ** 3, rng.next()));
        out[0] *= r; out[1] *= r; out[2] *= r;
        if (avoidOct && inOctant(out, 0.4)) continue;
        return out;
      }
      out[0] = -Math.abs(out[0]); return out;
    }
    const randomCyto = (out) => randomIn(R_NUC + 0.7, R_CELL - 0.7, true, out);
    const randomNuc = (out) => randomIn(R_CHROM + 0.4, R_NUC - 0.4, false, out);
    function randomExtra(out) {
      for (let k = 0; k < 40; k++) {
        randomIn(R_CELL + 1.2, R_EXTRA, false, out);
        if (Math.hypot(out[0] - SENDER_C[0], out[1] - SENDER_C[1], out[2] - SENDER_C[2]) > R_SEND + 1) return out;
      }
      return out;
    }
    function nearestPore(x, y, z) {
      const r = Math.hypot(x, y, z) || 1;
      let best = 0, bd = -2;
      for (let k = 0; k < poreDirs.length; k++) {
        const d = poreDirs[k], c = (d[0] * x + d[1] * y + d[2] * z) / r;
        if (c > bd) { bd = c; best = k; }
      }
      return best;
    }

    // spawn locations per species
    const T3 = [0, 0, 0];
    function spawn(key, st, s0 = 0) {
      const sw = SW[key];
      if (key === 'bcat') randomCyto(T3);
      else if (key === 'nicd') {
        // released by gamma-secretase at the cleaved receptor, just inside the membrane
        const pool = S_bridge.on.some((v) => v) ? bridgeSlots : nextSlots;
        const sl = pool[rng.int(Math.min(pool.length, 12))];
        for (let j = 0; j < 3; j++) T3[j] = sl.d[j] * (R_CELL - 0.8) + rng.normal(0, 0.35);
      } else if (key === 'gliA') {
        const i = sw.add(ciliumTip[0], ciliumTip[1], ciliumTip[2], ST.CILIUM, s0);
        return i;
      } else if (key === 'gliR') {
        for (let j = 0; j < 3; j++) T3[j] = ciliumBase[j] * 0.92 + rng.normal(0, 0.8);
      } else randomExtra(T3);
      const i = sw.add(T3[0], T3[1], T3[2], st, s0);
      if (i >= 0 && st === ST.NUC) { randomNuc(T3); sw.x[i] = T3[0]; sw.y[i] = T3[1]; sw.z[i] = T3[2]; }
      if (i >= 0 && key !== 'nicd' && st === ST.CYTO) avoidOctant(sw, i);
      return i;
    }
    function beginImport(sw, i) {
      const k = nearestPore(sw.x[i], sw.y[i], sw.z[i]);
      sw.st[i] = ST.IMPORT; sw.aux[i] = k; sw.t[i] = 0;
    }
    function beginExport(sw, i) {
      const k = nearestPore(sw.x[i], sw.y[i], sw.z[i]);
      sw.st[i] = ST.EXPORT; sw.aux[i] = k; sw.t[i] = 0;
    }
    function doom(sw, i, nDC) {
      if (nDC <= 0) { sw.st[i] = ST.DYING; return; }
      let best = 0, bd = Infinity;
      for (let j = 0; j < nDC; j++) {
        const q = dcPos[j], d = (sw.x[i] - q[0]) ** 2 + (sw.y[i] - q[1]) ** 2 + (sw.z[i] - q[2]) ** 2;
        if (d < bd) { bd = d; best = j; }
      }
      sw.st[i] = ST.DOOMED; sw.aux[i] = best;
    }

    /** Drive a two-compartment swarm toward target cytoplasmic / nuclear counts. */
    function control(key, tCyto, tNuc, dt, removeCyto, nDC) {
      const sw = SW[key];
      tCyto = Math.round(tCyto); tNuc = Math.round(tNuc);
      let cyto = 0, nuc = 0;
      for (let i = 0; i < sw.n; i++) {
        const s = sw.st[i];
        if (s === ST.CYTO || s === ST.EXPORT || s === ST.PORE_OUT || s === ST.CILIUM) cyto++;
        else if (s === ST.NUC || s === ST.IMPORT || s === ST.PORE_IN) nuc++;
      }
      const steps = (d) => Math.min(Math.abs(d), 1 + Math.floor(Math.abs(d) * Math.min(1, dt * 4)));
      let d = tNuc - nuc;
      if (d > 0) {
        for (let k = steps(d); k > 0; k--) {
          let i = sw.pick(ST.CYTO, rng);
          if (i < 0) { i = spawn(key, ST.CYTO); if (i < 0) break; } else cyto--;
          if (sw.st[i] === ST.CILIUM) break;
          beginImport(sw, i);
        }
      } else if (d < 0) {
        for (let k = steps(d); k > 0; k--) {
          const i = sw.pick(ST.NUC, rng);
          if (i < 0) break;
          if (key === 'bcat' && cyto < tCyto) { beginExport(sw, i); cyto++; } else sw.st[i] = ST.DYING;
        }
      }
      d = tCyto - cyto;
      if (d > 0) {
        for (let k = steps(d); k > 0; k--) if (spawn(key, ST.CYTO) < 0) break;
      } else if (d < 0) {
        for (let k = steps(d); k > 0; k--) {
          const i = sw.pick(ST.CYTO, rng);
          if (i < 0) break;
          if (removeCyto === 'doom') doom(sw, i, nDC); else sw.st[i] = ST.DYING;
        }
      }
    }
    function controlExtra(key, target, dt) {
      const sw = SW[key];
      target = Math.round(target);
      let n = 0;
      for (let i = 0; i < sw.n; i++) if (sw.st[i] === ST.EXTRA) n++;
      const d = target - n, steps = Math.min(Math.abs(d), 1 + Math.floor(Math.abs(d) * Math.min(1, dt * 4)));
      if (d > 0) for (let k = 0; k < steps; k++) spawn(key, ST.EXTRA);
      else for (let k = 0; k < steps; k++) { const i = sw.pick(ST.EXTRA, rng); if (i < 0) break; sw.st[i] = ST.DYING; }
    }
    /** Poisson-like turnover: `rate` events per animation second applied by fn(i). */
    function turnover(sw, which, rate, dt, state, fn) {
      sw[which] += rate * dt;
      while (sw[which] >= 1) {
        sw[which] -= 1;
        const i = sw.pick(state, rng);
        if (i >= 0) fn(i);
      }
    }

    function moveSwarm(key, dt, locus) {
      const sw = SW[key];
      const aIm = 1 - Math.exp(-dt / IMPORT_TAU);
      for (let i = 0; i < sw.n; i++) {
        const st = sw.st[i];
        if (st !== ST.DYING && sw.s[i] < 1) sw.s[i] = Math.min(1, sw.s[i] + dt * 3);
        switch (st) {
          case ST.CYTO:
            ou(sw, i, dt, 1.0, 0.6);
            shellClamp(sw, i, R_NUC + 0.55, R_CELL - 0.5);
            avoidOctant(sw, i);
            break;
          case ST.NUC: {
            ou(sw, i, dt, 0.7, 0.6);
            if (locus) {
              const k = 0.8 * dt;
              sw.vx[i] += (locus[0] * 1.25 - sw.x[i]) * k; sw.vy[i] += (locus[1] * 1.25 - sw.y[i]) * k; sw.vz[i] += (locus[2] * 1.25 - sw.z[i]) * k;
            }
            shellClamp(sw, i, R_CHROM + 0.35, R_NUC - 0.35);
            break;
          }
          case ST.EXTRA:
            ou(sw, i, dt, 1.4, 0.8);
            shellClamp(sw, i, R_CELL + 0.9, R_EXTRA);
            shellClamp(sw, i, R_SEND + 0.9, 1e9, SENDER_C[0], SENDER_C[1], SENDER_C[2]);
            break;
          case ST.IMPORT: case ST.EXPORT: {
            const d = poreDirs[sw.aux[i]], r = st === ST.IMPORT ? R_NUC + 0.65 : R_NUC - 0.65;
            const tx = d[0] * r, ty = d[1] * r, tz = d[2] * r;
            sw.x[i] += (tx - sw.x[i]) * aIm; sw.y[i] += (ty - sw.y[i]) * aIm; sw.z[i] += (tz - sw.z[i]) * aIm;
            if (st === ST.IMPORT) shellClamp(sw, i, R_NUC + 0.55, R_CELL - 0.5);
            if ((tx - sw.x[i]) ** 2 + (ty - sw.y[i]) ** 2 + (tz - sw.z[i]) ** 2 < 0.09) { sw.st[i] = st === ST.IMPORT ? ST.PORE_IN : ST.PORE_OUT; sw.t[i] = 0; }
            break;
          }
          case ST.PORE_IN: case ST.PORE_OUT: {
            sw.t[i] += dt;
            const f = Math.min(1, sw.t[i] / DWELL), d = poreDirs[sw.aux[i]];
            const r = st === ST.PORE_IN ? THREE.MathUtils.lerp(R_NUC + 0.65, R_NUC - 0.7, f) : THREE.MathUtils.lerp(R_NUC - 0.65, R_NUC + 0.7, f);
            sw.x[i] = d[0] * r; sw.y[i] = d[1] * r; sw.z[i] = d[2] * r;
            if (f >= 1) { sw.st[i] = st === ST.PORE_IN ? ST.NUC : ST.CYTO; sw.vx[i] = sw.vy[i] = sw.vz[i] = 0; }
            break;
          }
          case ST.CILIUM: {
            sw.t[i] += dt;
            const f = Math.min(1, sw.t[i] / 1.1);
            for (let j = 0; j < 3; j++) T3[j] = THREE.MathUtils.lerp(ciliumTip[j], ciliumBase[j], f);
            sw.x[i] = T3[0]; sw.y[i] = T3[1]; sw.z[i] = T3[2];
            if (f >= 1) sw.st[i] = ST.CYTO;
            break;
          }
          case ST.DOOMED: {
            const q = dcPos[sw.aux[i]] ?? dcPos[0], a = 1 - Math.exp(-dt / 0.35);
            sw.x[i] += (q[0] - sw.x[i]) * a; sw.y[i] += (q[1] - sw.y[i]) * a; sw.z[i] += (q[2] - sw.z[i]) * a;
            if ((q[0] - sw.x[i]) ** 2 + (q[1] - sw.y[i]) ** 2 + (q[2] - sw.z[i]) ** 2 < 0.5) {
              sw.st[i] = ST.DYING;
              dcPulse[sw.aux[i]] = 1;
            }
            break;
          }
          case ST.DYING:
            sw.s[i] -= dt * 3.5;
            if (sw.s[i] <= 0) { sw.kill(i); i--; }
            break;
          default: break;
        }
      }
    }

    const tmpC = new THREE.Color();
    const SIZE_BCAT = [0.5, 0.26, 0.26], SIZE_WNT = [0.5, 0.26, 0.26], SIZE_SMALL = [0.32], SIZE_SHH = [0.36];
    function drawSwarm(key, size, elongated, colorFn) {
      const sw = SW[key], pool = sw.pool;
      pool.begin();
      for (let i = 0; i < sw.n; i++) {
        const s = Math.max(sw.s[i], 0.001);
        const c = colorFn(sw.st[i], tmpC);
        if (elongated) {
          P3[0] = sw.x[i]; P3[1] = sw.y[i]; P3[2] = sw.z[i];
          A3[0] = sw.vx[i] + 1e-3; A3[1] = sw.vy[i]; A3[2] = sw.vz[i];
          pool.putOriented(P3, A3, size[0] * s, size[1] * s, size[2] * s, c);
        } else pool.put(sw.x[i], sw.y[i], sw.z[i], size[0] * s, size[0] * s, size[0] * s, c);
      }
    }

    // ---------------------------------------------------------------- per-frame drawing of receptors
    const up = (p, d, k) => { P3[0] = p[0] + d[0] * k; P3[1] = p[1] + d[1] * k; P3[2] = p[2] + d[2] * k; return P3; };
    const negContact = CONTACT.map((v) => -v);
    const cA = new THREE.Color(), cB = new THREE.Color();
    function drawReceptors(dt, st) {
      // Wnt receptors: every Frizzled/LRP6 is drawn, bound ones light up and carry a Wnt ligand
      fzPool.begin(); dvlPool.begin();
      for (let k = 0; k < fzSlots.length; k++) {
        const sl = fzSlots[k], w = S_fzBound.w[k];
        fzPool.putOriented(sl.p, sl.d, 1, 1, 1, cA.copy(COL.fz).lerp(COL.fzOn, w));
        if (w > 0.05 && st.Dv > 0.03) {
          const r = 0.22 + 0.2 * clamp01(st.Dv / 0.75);
          dvlPool.put(...up(sl.p, sl.d, -0.8), r * w, r * w, r * w, COL.dvl);
        }
      }
      fzPool.end(); dvlPool.end();

      // Notch: free receptors, Delta:Notch bridges, S2-cleaved NEXT stubs
      notchPool.begin();
      for (let k = 0; k < notchFreeSlots.length; k++) {
        const w = S_notchFree.w[k];
        if (w > 0.02) notchPool.putOriented(notchFreeSlots[k].p, notchFreeSlots[k].d, w, w, w, COL.notch);
      }
      notchPool.end();
      bridgePool.begin();
      for (let k = 0; k < bridgeSlots.length; k++) {
        const w = S_bridge.w[k];
        if (w < 0.02) continue;
        const b = bridgeSlots[k], L = (b.half + 0.12) * w;
        bridgePool.putOriented(b.p, CONTACT, L, 1, 1, COL.notch);
        bridgePool.putOriented(b.q, negContact, L, 1, 1, COL.delta);
      }
      bridgePool.end();
      nextPool.begin();
      for (let k = 0; k < nextSlots.length; k++) {
        const w = S_next.w[k];
        if (w > 0.02) nextPool.putOriented(nextSlots[k].p, nextSlots[k].d, w, 1, 1, COL.next);
      }
      nextPool.end();
      deltaPool.begin();
      for (let k = 0; k < deltaSlots.length; k++) {
        const w = S_delta.w[k];
        if (w > 0.02) deltaPool.putOriented(deltaSlots[k].p, deltaSlots[k].d, w, w, w, COL.delta);
      }
      deltaPool.end();

      // Hedgehog: PTCH1 in the cilium (free) and around its base (Shh-bound); SMO enters the cilium when active
      ptchPool.begin();
      for (let k = 0; k < ptchCilSlots.length; k++) {
        const w = S_ptchCil.w[k];
        if (w > 0.02) ptchPool.putOriented(ptchCilSlots[k].p, ptchCilSlots[k].d, 0.8 * w, 0.8 * w, 0.8 * w, COL.ptch);
      }
      for (let k = 0; k < ptchMemSlots.length; k++) {
        const wb = S_ptchMemBound.w[k], wf = S_ptchMemFree.w[k], w = Math.max(wb, wf);
        if (w > 0.02) ptchPool.putOriented(ptchMemSlots[k].p, ptchMemSlots[k].d, w, w, w, cB.copy(COL.ptch).multiplyScalar(wb > wf ? 0.75 : 1));
      }
      ptchPool.end();
      smoPool.begin();
      for (let k = 0; k < smoCilSlots.length; k++) {
        const w = S_smoCil.w[k];
        if (w > 0.02) smoPool.putOriented(smoCilSlots[k].p, smoCilSlots[k].d, 0.8 * w, 0.8 * w, 0.8 * w, COL.smoOn);
      }
      for (let k = 0; k < smoMemSlots.length; k++) {
        const w = S_smoMem.w[k];
        if (w > 0.02) smoPool.putOriented(smoMemSlots[k].p, smoMemSlots[k].d, w, w, w, COL.smo);
      }
      smoPool.end();

      // destruction complexes (count follows Axin2, position follows Dvl recruitment, brightness follows activity)
      dcPool.begin();
      const sig = THREE.MathUtils.smoothstep(st.Dv, 0.08, 0.6);
      const act = clamp01(st.dcAct / 1.6);
      for (let j = 0; j < DC_MAX; j++) {
        dcW[j] += ((j < st.nDC ? 1 : 0) - dcW[j]) * Math.min(1, 2.5 * dt);
        dcPulse[j] = Math.max(0, dcPulse[j] - dt * 3);
        const q = dcPos[j];
        for (let m = 0; m < 3; m++) q[m] = THREE.MathUtils.lerp(dcRest[j][m], dcSig[j][m], sig);
        if (dcW[j] < 0.02) continue;
        const s = dcW[j] * (1 + 0.35 * dcPulse[j]);
        dcPool.putOriented(q, dcSig[j], s, s, s, cA.copy(COL.dc).multiplyScalar(0.45 + 0.55 * act + 0.4 * dcPulse[j]));
      }
      dcPool.end();

      // target-gene loci (size and brightness = target expression)
      lociPool.begin();
      for (const L of LOCI) {
        const a = st.loci[L.key];
        const r = 0.28 + 0.55 * a;
        lociPool.put(L.p[0], L.p[1], L.p[2], r, r, r, cB.copy(L.color).multiplyScalar(0.25 + 0.95 * a));
      }
      lociPool.end();
    }

    // bound ligands are drawn in the ligand pools on top of their receptors
    function drawLigands() {
      const wsw = SW.wnt, ssw = SW.shh;
      drawSwarm('wnt', SIZE_WNT, true, colWnt);
      for (let k = 0; k < fzSlots.length; k++) {
        const w = S_fzBound.w[k];
        if (w > 0.03) wsw.pool.putOriented(up(fzSlots[k].p, fzSlots[k].d, 1.08), fzSlots[k].d, 0.26 * w, 0.5 * w, 0.26 * w, COL.wnt);
      }
      wsw.pool.end();
      drawSwarm('shh', SIZE_SHH, false, colShh);
      for (let k = 0; k < ptchMemSlots.length; k++) {
        const w = S_ptchMemBound.w[k];
        if (w > 0.03) ssw.pool.put(...up(ptchMemSlots[k].p, ptchMemSlots[k].d, 0.98), 0.36 * w, 0.36 * w, 0.36 * w, COL.shh);
      }
      ssw.pool.end();
    }

    // ---------------------------------------------------------------- derived state for drawing + tooltips
    const stA = { Dv: 0, dcAct: 1, nDC: 6, loci: { W: 0, N: 0, H: 0 }, fb: 0 };
    function deriveA() {
      stA.Dv = yW[2];
      stA.dcAct = (1 + yW[5]) / (1 + yW[2] / pW.KDvl);
      stA.nDC = Math.min(DC_MAX, Math.round(3.6 * (1 + yW[5])));
      stA.loci.W = clamp01(yW[6] / 2.4);
      stA.loci.N = clamp01(yN[5] / 2.4);
      stA.loci.H = clamp01(yH[5] / 3.2);
      stA.fb = hedgehog.boundFraction(pH.shh, pH.Kd);
      stA.readout = hedgehog.readout(yH, pH);
    }

    function setSlotTargets() {
      S_fzBound.setCount(Math.round(fzSlots.length * clamp01(yW[1] / pW.FzTot)));
      S_notchFree.setCount(Math.round(8 * yN[0]));
      S_bridge.setCount(Math.round(16 * yN[1]));
      S_next.setCount(Math.min(nextSlots.length, Math.round(2.5 * yN[2])));
      S_delta.setCount(Math.min(deltaSlots.length, Math.round(22 * pN.deltaTrans)));
      const nP = Math.round(7 * yH[0]), nB = Math.round(nP * stA.fb), nF = nP - nB;
      S_ptchCil.setCount(Math.min(nF, ptchCilSlots.length));
      S_ptchMemBound.setCount(Math.min(nB, ptchMemSlots.length));
      S_ptchMemFree.setCountFromEnd(Math.min(Math.max(nF - ptchCilSlots.length, 0), ptchMemSlots.length - Math.min(nB, ptchMemSlots.length)));
      const nAct = Math.round(smoCilSlots.length * clamp01(yH[1]));
      S_smoCil.setCount(nAct);
      S_smoMem.setCount(smoMemSlots.length - nAct);
    }

    function populateA() {
      // fill every swarm at its model target immediately (no empty first frame)
      for (const sw of Object.values(SW)) sw.clear();
      const fill = (key, n, st) => {
        const sw = SW[key];
        for (let k = 0; k < Math.round(n); k++) {
          const i = spawn(key, st, 1);
          if (i < 0) break;
          sw.st[i] = st;
          if (st === ST.NUC) randomNuc(T3);
          else if (st === ST.CYTO) randomCyto(T3);
          else continue;
          sw.x[i] = T3[0]; sw.y[i] = T3[1]; sw.z[i] = T3[2];
        }
      };
      fill('bcat', yW[3] * SCALE.bcat, ST.CYTO); fill('bcat', yW[4] * SCALE.bcat, ST.NUC);
      fill('nicd', yN[3] * SCALE.nicd, ST.CYTO); fill('nicd', yN[4] * SCALE.nicd, ST.NUC);
      fill('gliA', yH[3] * SCALE.gliA * (1 - GLI_NUC), ST.CYTO); fill('gliA', yH[3] * SCALE.gliA * GLI_NUC, ST.NUC);
      fill('gliR', yH[4] * SCALE.gliR * (1 - GLI_NUC), ST.CYTO); fill('gliR', yH[4] * SCALE.gliR * GLI_NUC, ST.NUC);
      fill('wnt', yW[0] * SCALE.wnt, ST.EXTRA);
      fill('shh', shhCount(), ST.EXTRA);
      deriveA();
      setSlotTargets();
      for (const s of allSlotSets) s.snap();
    }
    const shhCount = () => 70 * pH.shh / (pH.shh + 5);

    // per-frame callbacks hoisted out of the frame loop (no closures allocated per frame)
    const colBcat = (s, c) => c.copy(s === ST.DOOMED || s === ST.DYING ? COL.bcatP : COL.bcat);
    const colNicd = (s, c) => c.copy(COL.nicd), colGliA = (s, c) => c.copy(COL.gliA), colGliR = (s, c) => c.copy(COL.gliR);
    const colWnt = (s, c) => c.copy(COL.wnt), colShh = (s, c) => c.copy(COL.shh);
    const doomBcat = (i) => doom(SW.bcat, i, stA.nDC), exportBcat = (i) => beginExport(SW.bcat, i);
    const dieNicd = (i) => { SW.nicd.st[i] = ST.DYING; }, dieGliA = (i) => { SW.gliA.st[i] = ST.DYING; }, dieGliR = (i) => { SW.gliR.st[i] = ST.DYING; };

    function updateCell(dt) {
      if (dt > 0) integrateA(dt * TIME_A);
      deriveA();
      setSlotTargets();
      for (const s of allSlotSets) s.ease(dt);
      if (dt > 0) {
        const perSec = TIME_A * TURNOVER;
        // beta-catenin: destruction-complex degradation and nuclear export (import replaces them)
        const dcRate = (pW.kB + pW.kDC * stA.dcAct) * yW[3] * SCALE.bcat * perSec;
        turnover(SW.bcat, 'acc', dcRate, dt, ST.CYTO, doomBcat);
        turnover(SW.bcat, 'acc2', pW.kout * yW[4] * SCALE.bcat * perSec, dt, ST.NUC, exportBcat);
        // NICD degradation (Fbw7), GLI turnover
        turnover(SW.nicd, 'acc', pN.dNICD * yN[4] * SCALE.nicd * perSec, dt, ST.NUC, dieNicd);
        turnover(SW.gliA, 'acc', pH.dA * yH[3] * GLI_NUC * SCALE.gliA * perSec, dt, ST.NUC, dieGliA);
        turnover(SW.gliR, 'acc', pH.dR * yH[4] * GLI_NUC * SCALE.gliR * perSec, dt, ST.NUC, dieGliR);
        control('bcat', yW[3] * SCALE.bcat, yW[4] * SCALE.bcat, dt, 'doom', stA.nDC);
        control('nicd', yN[3] * SCALE.nicd, yN[4] * SCALE.nicd, dt, 'die', 0);
        control('gliA', yH[3] * SCALE.gliA * (1 - GLI_NUC), yH[3] * SCALE.gliA * GLI_NUC, dt, 'die', 0);
        control('gliR', yH[4] * SCALE.gliR * (1 - GLI_NUC), yH[4] * SCALE.gliR * GLI_NUC, dt, 'die', 0);
        controlExtra('wnt', yW[0] * SCALE.wnt, dt);
        controlExtra('shh', shhCount(), dt);
        moveSwarm('bcat', dt, LOCI[0].p);
        moveSwarm('nicd', dt, LOCI[1].p);
        moveSwarm('gliA', dt, LOCI[2].p);
        moveSwarm('gliR', dt, LOCI[2].p);
        moveSwarm('wnt', dt, null);
        moveSwarm('shh', dt, null);
      }
      drawReceptors(dt, stA);
      drawLigands();
      drawSwarm('bcat', SIZE_BCAT, true, colBcat);
      SW.bcat.pool.end();
      drawSwarm('nicd', SIZE_SMALL, false, colNicd);
      SW.nicd.pool.end();
      drawSwarm('gliA', SIZE_SMALL, false, colGliA);
      SW.gliA.pool.end();
      drawSwarm('gliR', SIZE_SMALL, false, colGliR);
      SW.gliR.pool.end();

      // nucleus glow = weighted target-gene activity
      const { W, N, H } = stA.loci;
      const eu = envMat.uniforms;
      eu.uColor.value.setRGB(0.5, 0.56, 0.85)
        .lerp(COL.bcat, 0.45 * W / (1 + W + N + H))
        .lerp(COL.notch, 0.45 * N / (1 + W + N + H))
        .lerp(COL.gliA, 0.45 * H / (1 + W + N + H));
      eu.uEmissive.value = 0.08 + 0.45 * Math.max(W, N, H);
    }

    // ================================================================ (b) spheroid
    const SPH_N = 300, SPH_R0 = 5;
    const sph = buildSpheroid(SPH_N, SPH_R0);
    const beadR = 8, SRC_K = 10;
    // beads touch the polar facet of the spheroid; the K cells nearest to each bead form
    // its contact patch, held at the bead concentration (Dirichlet source)
    const xCellMax = sph.xMax - SPH_R0, xCellMin = sph.xMin + SPH_R0;
    const wntBeadC = [xCellMax + SPH_R0 + beadR - 1, 0, 0];
    const shhBeadC = [xCellMin - SPH_R0 - beadR + 1, 0, 0];
    const srcW = new Uint8Array(SPH_N), srcS = new Uint8Array(SPH_N);
    const nearest = (C) => sph.pos.map((p, i) => [Math.hypot(p[0] - C[0], p[1] - C[1], p[2] - C[2]), i]).sort((a, b) => a[0] - b[0]).slice(0, SRC_K);
    for (const [, i] of nearest(wntBeadC)) srcW[i] = 1;
    for (const [, i] of nearest(shhBeadC)) srcS[i] = 1;
    const pB = { wnt: 0.1, lamW: 2.2, shh: 20, lamS: 1.5, chi: 0.5, dapt: false, cyclopamine: false, kDecay: 0.25, beads: true };
    const u = new Float64Array(SPH_N), v = new Float64Array(SPH_N), lapTmp = new Float64Array(SPH_N);
    const yWc = new Float64Array(SPH_N * 7), yHc = new Float64Array(SPH_N * 6);
    const yWv = Array.from({ length: SPH_N }, (_, i) => yWc.subarray(i * 7, i * 7 + 7));
    const yHv = Array.from({ length: SPH_N }, (_, i) => yHc.subarray(i * 6, i * 6 + 6));
    const pWc = { ...canonicalWnt.defaults }, pHc = { ...hedgehog.defaults };
    const fWc = canonicalWnt.rhs(pWc), fHc = hedgehog.rhs(pHc);
    const wWc = rk4Work(7), wHc = rk4Work(6);
    const noNbrs = sph.nbrs.map(() => []);
    const readoutB = new Float64Array(SPH_N);
    let notchS, tB = 0, lastSampleB = -1, lastProfileB = -1;
    // Collier et al. (1996) functional forms. With the default feedback (k = h = 2, b = 100)
    // the homogeneous state is only weakly unstable on this contact graph (linearised gain
    // f'(D*)·|g'(N*)|·|λ_min| ≈ 1.2, λ_min ≈ −0.47 = most negative eigenvalue of the
    // neighbour-averaging operator), so the pattern needs ~100 τ to grow; a steeper feedback
    // (k = h = 3, b = 1000) patterns in ~15 τ with the same sender fraction (~24 %).
    const COLLIER_SETS = {
      steep: { a: 0.01, b: 1000, k: 3, h: 3, v: 1 },
      collier: { ...collier.defaults },
    };
    const pColl = { ...COLLIER_SETS.steep };

    function resetSpheroid() {
      const r = new RNG(777);
      notchS = collier.init(SPH_N, r, 0.2);
      u.fill(0); v.fill(0);
      const w0 = canonicalWnt.initial(), h0 = hedgehog.initial();
      for (let i = 0; i < SPH_N; i++) { yWv[i].set(w0); yHv[i].set(h0); }
      tB = 0; lastSampleB = -1; lastProfileB = -1;
      liChart?.clear(); fateChart?.clear();
    }

    function stepSpheroid(dtm) {
      const n = Math.max(1, Math.ceil(dtm / 0.04)), h = dtm / n;
      // the graph Laplacian approximates (<deg>/6)·ℓ²·∇², so the per-contact exchange rate
      // κ = 6 k λ² / <deg> gives a continuum decay length λ (in cell spacings ℓ)
      const DW = (6 * pB.kDecay * pB.lamW * pB.lamW) / sph.meanDeg, DS = (6 * pB.kDecay * pB.lamS * pB.lamS) / sph.meanDeg;
      const hRD = 0.4 / (2 * Math.max(DW, DS) * sph.maxDeg + pB.kDecay);
      const nRD = Math.max(1, Math.ceil(h / hRD)), hr = h / nRD;
      const nb = pB.dapt ? noNbrs : sph.nbrs;
      pHc.kSa = pB.cyclopamine ? 0 : hedgehog.defaults.kSa;
      for (let k = 0; k < n; k++) {
        // morphogens: du/dt = D sum_j (u_j - u_i) - k u, u clamped to c0 at the bead (Dirichlet source)
        for (let r = 0; r < nRD; r++) {
          rdStep(u, DW, hr, srcW, pB.wnt);
          rdStep(v, DS, hr, srcS, pB.shh);
        }
        // lateral inhibition + Wnt -> Delta/Jagged1 crosstalk
        collier.step(notchS, nb, h, pColl);
        if (pB.chi > 0) {
          for (let i = 0; i < SPH_N; i++) {
            const bn = yWc[i * 7 + 4], g = 1 / (1 + pColl.b * notchS.N[i] ** pColl.h);
            notchS.D[i] += h * pColl.v * pB.chi * (bn / (bn + 1)) * g;
          }
        }
        // per-cell pathway read-out
        for (let i = 0; i < SPH_N; i++) {
          pWc.wntIn = u[i];
          rk4Step(fWc, tB, yWv[i], h, wWc);
          pHc.shh = v[i];
          rk4Step(fHc, tB, yHv[i], h, wHc);
        }
        tB += h;
      }
      clampNonNeg(yWc); clampNonNeg(yHc);
    }
    function rdStep(f, D, h, src, c0) {
      const nbrs = sph.nbrs;
      for (let i = 0; i < SPH_N; i++) {
        let s = 0;
        const nb = nbrs[i];
        for (let k = 0; k < nb.length; k++) s += f[nb[k]] - f[i];
        lapTmp[i] = s;
      }
      for (let i = 0; i < SPH_N; i++) f[i] = src[i] ? c0 : Math.max(0, f[i] + h * (D * lapTmp[i] - pB.kDecay * f[i]));
    }

    const WNT_HI = 0.7, SENDER = 0.5;
    const FATES = [
      { name: 'stem / progenitor (Wnt-high, Notch-high)', color: '#5ee0c1' },
      { name: 'Wnt-dependent secretory (Paneth-like)', color: '#ffb454' },
      { name: 'absorptive-like progenitor (Wnt-low, Notch-high)', color: '#4aa3ff' },
      { name: 'secretory (goblet-like; Wnt-low, Delta-high)', color: '#ff3b9d' },
    ];
    const FATE_COL = FATES.map((f) => new THREE.Color(f.color));
    const DOMAIN_COL = NEURAL_TUBE_DOMAINS.map((d) => new THREE.Color(d.color));
    const fateOf = (i) => {
      const wHi = yWc[i * 7 + 4] > WNT_HI, send = notchS.D[i] > SENDER;
      return wHi ? (send ? 1 : 0) : (send ? 3 : 2);
    };

    // rendering
    const bodyMatB = cellMaterial({ role: 'membrane', opacity: 0.1, useInst: 1, gain: 0.3 });
    const nucMatB = cellMaterial({ role: 'nucleus', color: 0xffffff, useInst: 1, gain: 0.55, side: THREE.DoubleSide });
    const bodyPool = new InstancePool(unitSphere(2), bodyMatB, SPH_N, groupB);
    bodyPool.mesh.renderOrder = 3;
    const nucPoolB = new InstancePool(unitSphere(2), nucMatB, SPH_N, groupB);
    // confocal and H&E are sectioning techniques: show the spheroid as a ~1-cell-thick section
    const SLAB = 6;
    const slabPlanes = [new THREE.Plane(new THREE.Vector3(0, 0, -1), SLAB), new THREE.Plane(new THREE.Vector3(0, 0, 1), SLAB)];
    function setSlab(on) { for (const m of [bodyMatB, nucMatB]) { m.clippingPlanes = on ? slabPlanes : null; m.needsUpdate = true; } }
    const beadWnt = new THREE.Mesh(unitSphere(4), cellMaterial({ role: 'reporter', color: 0xff9a3c, emissive: 0.5 }));
    beadWnt.scale.setScalar(beadR); beadWnt.position.set(...wntBeadC);
    const beadShh = new THREE.Mesh(unitSphere(4), cellMaterial({ role: 'reporter', color: 0x45c8ff, emissive: 0.5 }));
    beadShh.scale.setScalar(beadR); beadShh.position.set(...shhBeadC);
    groupB.add(beadWnt, beadShh);
    const haloN = 170;
    const haloPool = new InstancePool(unitSphere(1), molMat({ emissive: 0.5 }), haloN * 2, groupB);
    const halo = Array.from({ length: haloN * 2 }, () => ({ p: [0, 0, 0], age: 0, life: 1 }));
    function respawnHalo(h, k) {
      const bead = k < haloN ? wntBeadC : shhBeadC;
      const lam = (k < haloN ? pB.lamW : pB.lamS) * sph.spacing;
      const d = rng.onSphere();
      const r = beadR + 0.6 + rng.exponential(1 / lam);
      for (let j = 0; j < 3; j++) h.p[j] = bead[j] + d[j] * r;
      h.age = 0; h.life = 1.5 + 2 * rng.next();
    }
    halo.forEach((h, k) => { respawnHalo(h, k); h.age = rng.next() * h.life; });

    const cellCol = new THREE.Color(), nucCol = new THREE.Color();
    // multi-stop ramps with a distinct low-end hue (readable under translucent membranes)
    const RAMP_BCAT = ['#2c3a8c', '#b0507a', '#ff9a3c', '#ffe066'].map((c) => new THREE.Color(c));
    const RAMP_GLI = ['#3a2c6e', '#2f6fb0', '#2fc4d8', '#b8fff4'].map((c) => new THREE.Color(c));
    function ramp(stops, t, out) {
      const x = clamp01(t) * (stops.length - 1), k = Math.min(stops.length - 2, Math.floor(x));
      return out.copy(stops[k]).lerp(stops[k + 1], x - k);
    }
    function colorCell(i, out) {
      switch (params.colorBy) {
        case 'notch': {
          const d = clamp01(notchS.D[i]);
          return out.setRGB(0.22 + 0.78 * d, 1 - 0.75 * d, 0.53 + 0.1 * d);
        }
        case 'bcat': return ramp(RAMP_BCAT, (yWc[i * 7 + 4] - 0.3) / 0.75, out);
        case 'gli': return ramp(RAMP_GLI, readoutB[i] / 0.9, out);
        case 'domains': return out.copy(DOMAIN_COL[frenchFlag(readoutB[i], NEURAL_TUBE_THRESHOLDS)]);
        case 'morph': {
          // each field's log level normalised to its range across the spheroid; hue = which
          // morphogen dominates, running Shh blue -> violet -> magenta -> Wnt orange (stays saturated)
          const w = clamp01((Math.log10(u[i] + 1e-9) - morphRange[0]) / (morphRange[1] - morphRange[0] || 1));
          const s = clamp01((Math.log10(v[i] + 1e-9) - morphRange[2]) / (morphRange[3] - morphRange[2] || 1));
          const f = w / (w + s + 1e-9);
          return out.setHSL((0.56 + 0.52 * f) % 1, 0.9, 0.2 + 0.4 * Math.max(w, s));
        }
        default: return out.copy(FATE_COL[fateOf(i)]);
      }
    }
    const morphRange = [0, 1, 0, 1];
    function drawSpheroid(dt) {
      for (let i = 0; i < SPH_N; i++) readoutB[i] = hedgehog.readout(yHv[i], pHc);
      if (params.colorBy === 'morph') {
        morphRange[0] = morphRange[2] = Infinity; morphRange[1] = morphRange[3] = -Infinity;
        for (let i = 0; i < SPH_N; i++) {
          const lu = Math.log10(u[i] + 1e-9), lv = Math.log10(v[i] + 1e-9);
          if (lu < morphRange[0]) morphRange[0] = lu; if (lu > morphRange[1]) morphRange[1] = lu;
          if (lv < morphRange[2]) morphRange[2] = lv; if (lv > morphRange[3]) morphRange[3] = lv;
        }
      }
      bodyPool.begin(); nucPoolB.begin();
      for (let i = 0; i < SPH_N; i++) {
        const p = sph.pos[i], R = sph.R[i];
        colorCell(i, cellCol);
        bodyPool.put(p[0], p[1], p[2], R * 1.1, R * 1.1, R * 1.1, nucCol.copy(cellCol).lerp(WHITE, 0.05));
        nucPoolB.putOriented(p, sph.axis[i], R * 0.58, R * 0.47, R * 0.47, cellCol);
      }
      bodyPool.end(); nucPoolB.end();
      haloPool.begin();
      if (pB.beads) {
        for (let k = 0; k < halo.length; k++) {
          const h = halo[k];
          h.age += dt;
          if (h.age > h.life) respawnHalo(h, k);
          const f = Math.sin(Math.PI * clamp01(h.age / h.life));
          const r = 0.42 * f;
          const on = k < haloN ? pB.wnt > 1e-4 : pB.shh > 1e-3;
          if (on && r > 0.02) haloPool.put(h.p[0], h.p[1], h.p[2], r, r, r, k < haloN ? COL.wnt : COL.shh);
        }
      }
      haloPool.end();
    }
    const WHITE = new THREE.Color(1, 1, 1);

    // ================================================================ panel
    const sections = { cell: [], sph: [] };
    const addSection = (title, view) => { panel.section(title); sections[view].push(panel.current); return panel.current; };

    panel.section('View');
    panel.select({
      label: 'Show', value: params.view,
      options: [
        { value: 'cell', label: 'Inside one stem cell (cut-away, three pathways)' },
        { value: 'sph', label: 'Stem-cell spheroid (lateral inhibition + gradients)' },
      ],
      onChange: (val) => { params.view = val; layout(); },
    });

    // ---- (a) controls
    addSection('Signals & inhibitors', 'cell');
    panel.slider({ label: 'Wnt3a supply (a.u.)', min: 0.001, max: 1, log: true, value: pW.wntIn, onChange: (val) => { pW.wntIn = val; } });
    panel.slider({ label: 'Delta on the neighbour (trans)', min: 0, max: 2, value: pN.deltaTrans, onChange: (val) => { pN.deltaTrans = val; } });
    panel.slider({ label: 'Shh (units of K_d)', min: 0.01, max: 100, log: true, value: pH.shh, onChange: (val) => { pH.shh = val; doseMarker(); } });
    panel.toggle({ label: 'DAPT (γ-secretase inhibitor: blocks S3 cleavage)', value: false, onChange: (on) => { pN.kS3 = on ? 0 : singleCell.defaults.kS3; } });
    panel.toggle({ label: 'Cyclopamine (SMO inhibitor)', value: false, onChange: (on) => { pH.kSa = on ? 0 : hedgehog.defaults.kSa; computeDoseResponse(); } });
    panel.toggle({ label: 'Axin2 negative feedback (Wnt)', value: true, onChange: (on) => { pW.kAx = on ? canonicalWnt.defaults.kAx : 0; } });
    panel.toggle({ label: 'Ptch1 negative feedback (Hh)', value: true, onChange: (on) => { pH.sP = on ? hedgehog.defaults.sP : 0; computeDoseResponse(); } });
    panel.buttons([{ label: 'Restart from an unstimulated cell', primary: true, onClick: () => { resetA(); } }]);
    panel.note('Particle counts are proportional to the model concentrations (β-catenin 40, NICD 22, GLI 24, Wnt 30 particles per unit); turnover events are drawn at 10 % of the model flux.');

    addSection("Paper's Wnt cascade (linear ODEs)", 'cell');
    panel.equation('dW/dt = k1 − k2·W\ndF/dt = k3·W − k4·F\ndD/dt = k5·F − k6·D\ndB/dt = k7·D − k8·B\nsteady state:\n  W* = k1/k2     F* = k3·W*/k4\n  D* = k5·F*/k6  B* = k7·D*/k8');
    const paperChart = panel.chart({
      title: 'W, F, D, B (k1 ∝ Wnt supply) and analytic B*', xLabel: 'τ',
      series: [{ name: 'W', color: '#ff9a3c' }, { name: 'F', color: '#c07a3e' }, { name: 'D', color: '#fff2a6' }, { name: 'B', color: '#ffd84a' }, { name: 'B*', color: '#8a98ab', dash: [4, 3] }],
    });
    const paperRO = panel.readouts(['W*  (now)', 'F*  (now)', 'D*  (now)', 'B*  (now)']);
    panel.note('The paper\'s cascade is linear: B* grows in proportion to k1 without limit. The mechanistic model below saturates (receptor occupancy, destruction-complex inhibition) and adapts (Axin2).');

    addSection('Canonical Wnt / β-catenin', 'cell');
    panel.equation('Fz*:  kon·W·(Fz_tot − Fz*) − koff·Fz*\nDvl*: kDa·Fz*·(1 − Dvl*) − kDi·Dvl*\nDC  = (1 + Axin2) / (1 + Dvl*/K)\nβc:   vB − (kB + kDC·DC)·βc\n        − kin·βc + kout·βn\nβn:   kin·βc − kout·βn\nAxin2: kAx·βn²/(K² + βn²) − dAx·Axin2');
    const wntChart = panel.chart({
      title: 'β-catenin, Axin2 and TCF target', xLabel: 'τ',
      series: [{ name: 'β-cat cyto', color: '#c9a53a' }, { name: 'β-cat nuc', color: '#ffd84a' }, { name: 'Axin2', color: '#a08cff' }, { name: 'target', color: '#fff2a6', dash: [4, 3] }],
    });
    const wntRO = panel.readouts(['Fz/LRP6 bound', 'Dvl active', 'destruction-complex activity']);

    addSection('Notch (receiver cell)', 'cell');
    panel.equation('Delta(trans) + Notch ⇌ DN\nDN → NEXT        (ADAM10, S2)\nNEXT → NICD      (γ-secretase, S3;\n                  blocked by DAPT)\nNICD_c ⇌ NICD_n → Hes1\nHes1 ⊣ own Delta (lateral inhibition)');
    const notchChart = panel.chart({
      title: 'NICD, Hes1 and the cell\'s own Delta', xLabel: 'τ',
      series: [{ name: 'NICD cyto', color: '#2a9d62' }, { name: 'NICD nuc', color: '#b8ffd8' }, { name: 'Hes1', color: '#39ff88' }, { name: 'own Delta', color: '#ff3b9d', dash: [4, 3] }],
    });
    const notchRO = panel.readouts(['surface Notch', 'Delta:Notch', 'NEXT (membrane stub)']);

    addSection('Hedgehog (primary cilium)', 'cell');
    panel.equation('f_b = Shh/(Shh + Kd)   (quasi-steady)\ndS/dt = kSa·(1 − S)\n        − (kSi + kPtc·(1 − f_b)·P)·S\ndG/dt = sG − (kGA·S + kGR·(1−S) + dG)·G\ndA/dt = kGA·S·G + sA1·φ − dA·A\ndR/dt = kGR·(1 − S)·G − dR·R\ndP/dt = sP0 + sP·φ − dP·P − kInt·f_b·P\nφ = (β + a²) / (1 + a² + r²)\n    a = A/KA,  r = R/KR\nS = SMO*, G = GLI2/3, A = GliA,\nR = GLI3R, P = PTCH1');
    const hhChart = panel.chart({
      title: 'GLI activator / repressor, PTCH1, SMO*', xLabel: 'τ',
      series: [{ name: 'GliA', color: '#5ef0ff' }, { name: 'GliR', color: '#8878ff' }, { name: 'PTCH1', color: '#3a74ff' }, { name: 'SMO*', color: '#f2c6ff' }, { name: 'φ', color: '#ffffff', dash: [4, 3] }],
    });
    const doseChart = panel.chart({
      title: 'Steady-state GLI read-out φ vs log10(Shh/Kd); dotted = French-flag thresholds, marker = current dose', xLabel: 'log', yRange: [0, 1],
      series: [{ name: 'Ptch1 feedback', color: '#5ef0ff' }, { name: 'no feedback', color: '#8a98ab', dash: [4, 3] },
        ...NEURAL_TUBE_THRESHOLDS.map(() => ({ name: '', color: '#ff6fae', dash: [2, 4], width: 1 }))],
    });
    const hhRO = panel.readouts(['PTCH1 bound (f_b)', 'SMO active', 'GLI read-out φ', 'neural-tube domain']);
    panel.note('Shh binding relieves the catalytic inhibition of Smoothened by Patched (Taipale et al. 2002); active SMO moves into the primary cilium while PTCH1 leaves it (Rohatgi, Milenkovic &amp; Scott 2007, Science 317:372). Ptch1 is a GLI target, so signalling desensitises the cell over time (Dessaud et al. 2007).');

    addSection('Nuclear import (measured kinetics)', 'cell');
    const npcChart = panel.chart({
      title: 'Nuclear-envelope enrichment, live-cell data (points) and fit', xLabel: 'frame',
      series: [{ name: 'measured rim/cytoplasm', color: '#ffb454', points: true, width: 1 }, { name: 'E(t) fit', color: '#5ee0c1' }],
    });
    if (npc?.rim_over_cytoplasm) {
      const xs = npc.rim_over_cytoplasm.map((_, i) => i);
      const { E0, E_inf, tau_frames } = npc.fit;
      npcChart.set(xs, [npc.rim_over_cytoplasm, xs.map((t) => E_inf - (E_inf - E0) * Math.exp(-t / tau_frames))]);
      panel.note(`E(t) = E∞ − (E∞ − E0)·e^(−t/τ) with E0 = ${fmt(E0)}, E∞ = ${fmt(E_inf)}, τ = ${fmt(tau_frames)} frames (R² = ${fmt(npc.fit.r2)}; ${npc.source}). Messengers approach their nuclear pore with the same single-exponential shape; the frame interval is not recorded, so τ is used as a relative (shape) calibration: ${fmt(IMPORT_TAU)} s of animation.`);
    } else panel.note('npc_kinetics.json not available: using τ = 13.5 frames.');

    // ---- (b) controls
    addSection('Spheroid', 'sph');
    panel.select({
      label: 'Colour cells by', value: params.colorBy,
      options: [
        { value: 'fate', label: 'Combinatorial fate (Wnt × Notch)' },
        { value: 'notch', label: 'Notch / Delta (lateral inhibition)' },
        { value: 'bcat', label: 'Nuclear β-catenin' },
        { value: 'gli', label: 'GLI read-out φ' },
        { value: 'domains', label: 'Shh French-flag domains' },
        { value: 'morph', label: 'Morphogens (Wnt orange, Shh blue)' },
      ],
      onChange: (val) => { params.colorBy = val; updateLegend(); },
    });
    panel.slider({ label: 'Wnt at the bead c0 (a.u.)', min: 0.001, max: 1, log: true, value: pB.wnt, onChange: (val) => { pB.wnt = val; } });
    panel.slider({ label: 'Wnt decay length λ_W (cells)', min: 0.8, max: 6, value: pB.lamW, onChange: (val) => { pB.lamW = val; } });
    panel.slider({ label: 'Shh at the bead (units of K_d)', min: 0.1, max: 300, log: true, value: pB.shh, onChange: (val) => { pB.shh = val; } });
    panel.slider({ label: 'Shh decay length λ_S (cells)', min: 0.8, max: 6, value: pB.lamS, onChange: (val) => { pB.lamS = val; } });
    panel.slider({ label: 'Wnt → Delta/Jagged1 crosstalk χ', min: 0, max: 2, value: pB.chi, onChange: (val) => { pB.chi = val; } });
    panel.select({
      label: 'Lateral-inhibition feedback', value: 'steep',
      options: [
        { value: 'steep', label: 'Steeper (b = 1000, k = h = 3) · patterns in ~15 τ' },
        { value: 'collier', label: 'Collier 1996 defaults (b = 100, k = h = 2) · ~100 τ' },
      ],
      onChange: (val) => { Object.assign(pColl, COLLIER_SETS[val]); collierEq(); },
    });
    panel.toggle({ label: 'DAPT (no Notch signalling between cells)', value: false, onChange: (on) => { pB.dapt = on; } });
    panel.toggle({ label: 'Cyclopamine (SMO inhibitor)', value: false, onChange: (on) => { pB.cyclopamine = on; } });
    panel.toggle({ label: 'Show beads & ligand clouds', value: true, onChange: (on) => { pB.beads = on; beadWnt.visible = beadShh.visible = on; } });
    panel.buttons([{ label: 'Restart patterning', primary: true, onClick: () => resetSpheroid() }]);

    addSection('Lateral inhibition on the 3D contact graph', 'sph');
    const collierEqEl = panel.html('');
    function collierEq() {
      collierEqEl.innerHTML = `<div class="eq">dN_i/dt = ⟨D⟩_i^k/(a + ⟨D⟩_i^k) − N_i\ndD_i/dt = v·(g_i − D_i)\ng_i = (1 + χ·β_i/(β_i + 1))\n      / (1 + b·N_i^h)\nCollier et al. 1996 form:\n  a = ${pColl.a}, b = ${pColl.b}, k = h = ${pColl.k}, v = ${pColl.v}\n⟨D⟩_i: mean Delta of the 3D contacts\nβ_i: nuclear β-catenin of cell i</div>`;
    }
    collierEq();
    const liChart = panel.chart({
      title: 'Delta-high (sender) fraction and sender–sender contacts', xLabel: 'τ', yRange: [0, 0.6],
      series: [{ name: 'sender fraction', color: '#ff3b9d' }, { name: 'sender–sender contacts', color: '#39ff88' }],
    });
    const liRO = panel.readouts(['cells', 'mean 3D contacts', 'senders', 'sender–sender contacts']);

    addSection('Morphogen gradients', 'sph');
    panel.equation('∂u/∂t = D·∇²u − k·u\n(paper: ∂u/∂t = D_u∇²u + f(u,v),\n here f = −k·u)\ncontact graph:\n  du_i/dt = κ·Σ_j (u_j − u_i) − k·u_i\n  κ = 6·D / (⟨deg⟩·ℓ²)\n  u = c0 on the bead\'s contact patch\ndecay length λ = √(D/k)\n(Crick 1970; Kicheva et al. 2007)');
    const profChart = panel.chart({
      title: 'Profiles along the Shh bead (left) → Wnt bead (right) axis', xLabel: 'x (µm)', yRange: [0, 1.05],
      series: [{ name: 'Wnt/c0', color: '#ff9a3c' }, { name: 'β-cat nuc', color: '#ffd84a', dash: [4, 3] }, { name: 'Shh/c0', color: '#45c8ff' }, { name: 'GLI φ', color: '#5ef0ff', dash: [4, 3] }],
    });
    const gradRO = panel.readouts(['λ_W', 'λ_S', 'cell spacing']);

    addSection('Combinatorial fates (Wnt × Notch)', 'sph');
    const fateChart = panel.chart({
      title: 'Fate fractions', xLabel: 'τ', yRange: [0, 1],
      series: FATES.map((f) => ({ name: f.name.split(' (')[0], color: f.color })),
    });
    panel.note('Phenomenological read-out of the intestinal-crypt logic: Notch keeps progenitors undifferentiated and its loss makes secretory cells (van Es et al. 2005, Nature 435:959; Fre et al. 2005, Nature 435:964); Wnt maintains stem cells and matures Paneth cells (van Es et al. 2005, Nat Cell Biol 7:381). Jagged1 is a β-catenin target (Estrach et al. 2006, Development 133:4427), the basis of the χ crosstalk term.');

    // ---------------------------------------------------------------- helpers bound to the panel
    function computeDoseResponse() {
      const xs = [], a = [], b = [];
      for (let k = 0; k <= 24; k++) xs.push(-2 + (4 * k) / 24);
      const q = { ...pH };
      const opts = { tEnd: 300, h: 0.03 };
      hedgehog.doseResponse({ ...q, sP: hedgehog.defaults.sP }, xs.map((x) => 10 ** x), opts).forEach((r) => a.push(r.readout));
      hedgehog.doseResponse({ ...q, sP: 0 }, xs.map((x) => 10 ** x), opts).forEach((r) => b.push(r.readout));
      doseChart.set(xs, [a, b, ...NEURAL_TUBE_THRESHOLDS.map((t) => xs.map(() => t))]);
      doseMarker();
    }
    function doseMarker() {
      doseChart.markers = [{ x: Math.log10(pH.shh), color: '#ffb454' }];
      doseChart.dirty = true;
    }
    computeDoseResponse();

    function sampleCharts() {
      const ss = paperWnt.steadyState(pP);
      paperChart.push(tA, [yP[0], yP[1], yP[2], yP[3], ss[3]]);
      wntChart.push(tA, [yW[3], yW[4], yW[5], yW[6]]);
      notchChart.push(tA, [yN[3], yN[4], yN[5], yN[6]]);
      hhChart.push(tA, [yH[3], yH[4], yH[0], yH[1], stA.readout]);
      const now = [yP[0], yP[1], yP[2], yP[3]];
      ['W*  (now)', 'F*  (now)', 'D*  (now)', 'B*  (now)'].forEach((k, i) => paperRO.set(k, `${fmt(ss[i])}  (${fmt(now[i])})`));
      wntRO.set('Fz/LRP6 bound', `${(100 * clamp01(yW[1] / pW.FzTot)).toFixed(0)} %`);
      wntRO.set('Dvl active', yW[2]);
      wntRO.set('destruction-complex activity', stA.dcAct);
      notchRO.set('surface Notch', yN[0]);
      notchRO.set('Delta:Notch', yN[1]);
      notchRO.set('NEXT (membrane stub)', yN[2]);
      hhRO.set('PTCH1 bound (f_b)', `${(100 * stA.fb).toFixed(0)} %`);
      hhRO.set('SMO active', `${(100 * yH[1]).toFixed(0)} %`);
      hhRO.set('GLI read-out φ', stA.readout);
      hhRO.set('neural-tube domain', NEURAL_TUBE_DOMAINS[frenchFlag(stA.readout)].name.split(' ·')[0]);
    }

    function sampleSpheroid() {
      const sf = collier.senderFraction(notchS, SENDER), ss = collier.senderAdjacency(notchS, sph.nbrs, SENDER);
      liChart.push(tB, [sf, ss]);
      const counts = [0, 0, 0, 0];
      for (let i = 0; i < SPH_N; i++) counts[fateOf(i)]++;
      fateChart.push(tB, counts.map((c) => c / SPH_N));
      liRO.set('cells', SPH_N);
      liRO.set('mean 3D contacts', sph.meanDeg);
      liRO.set('senders', `${(100 * sf).toFixed(1)} %`);
      liRO.set('sender–sender contacts', `${(100 * ss).toFixed(1)} %`);
      gradRO.set('λ_W', `${fmt(pB.lamW * sph.spacing)} µm`);
      gradRO.set('λ_S', `${fmt(pB.lamS * sph.spacing)} µm`);
      gradRO.set('cell spacing', `${fmt(sph.spacing)} µm`);
      return { sf, ss };
    }
    function profiles() {
      const NB = 16, x0 = sph.xMin, x1 = sph.xMax, sums = Array.from({ length: 5 }, () => new Float64Array(NB));
      for (let i = 0; i < SPH_N; i++) {
        const b = Math.min(NB - 1, Math.max(0, Math.floor(((sph.pos[i][0] - x0) / (x1 - x0)) * NB)));
        sums[0][b] += u[i] / pB.wnt; sums[1][b] += clamp01((yWc[i * 7 + 4] - 0.3) / 0.8);
        sums[2][b] += v[i] / pB.shh; sums[3][b] += readoutB[i]; sums[4][b]++;
      }
      const xs = [], cols = [[], [], [], []];
      for (let b = 0; b < NB; b++) {
        if (!sums[4][b]) continue;
        xs.push(x0 + ((b + 0.5) * (x1 - x0)) / NB);
        for (let k = 0; k < 4; k++) cols[k].push(sums[k][b] / sums[4][b]);
      }
      profChart.set(xs, cols);
    }

    // ---------------------------------------------------------------- legend, layout, picking
    const legendCell = [
      { color: HEX.wnt, label: 'Wnt3a' },
      { color: HEX.fzOn, label: 'Frizzled/LRP6 (bright = bound) · Dvl' },
      { color: HEX.dc, label: 'destruction complex APC·Axin·GSK-3β' },
      { color: HEX.bcat, label: 'β-catenin (brown = doomed)' },
      { color: HEX.delta, label: 'Delta (sender cell)' },
      { color: HEX.notch, label: 'Notch · Delta:Notch · NEXT' },
      { color: HEX.nicd, label: 'NICD' },
      { color: HEX.shh, label: 'Shh' },
      { color: HEX.ptch, label: 'Patched-1' },
      { color: HEX.smoOn, label: 'Smoothened (bright = active)' },
      { color: HEX.gliA, label: 'GLI activator' },
      { color: HEX.gliR, label: 'GLI3 repressor' },
      { color: '#6f8fff', label: 'chromatin · glowing loci = targets' },
    ];
    function updateLegend() {
      if (params.view === 'cell') { legend(legendCell); return; }
      const items = {
        fate: FATES.map((f) => ({ color: f.color, label: f.name })),
        notch: [{ color: '#39ff88', label: 'Notch-high (receiver)' }, { color: '#ff3b9d', label: 'Delta-high (sender)' }],
        bcat: [{ color: '#2c3a8c', label: 'low nuclear β-catenin (≈ 0.3)' }, { color: '#ff9a3c', label: 'intermediate' }, { color: '#ffe066', label: 'high nuclear β-catenin (≥ 1)' }],
        gli: [{ color: '#3a2c6e', label: 'GLI read-out φ ≈ 0' }, { color: '#2fc4d8', label: 'φ ≈ 0.6' }, { color: '#b8fff4', label: 'φ ≈ 0.9' }],
        domains: NEURAL_TUBE_DOMAINS.map((d) => ({ color: d.color, label: d.name })),
        morph: [{ color: '#ff9a3c', label: 'Wnt dominates' }, { color: '#e03cc8', label: 'balanced' }, { color: '#3aa0ff', label: 'Shh dominates' }, { color: '#555555', label: '(log levels, each scaled to its range)' }],
      }[params.colorBy];
      legend([...items, { color: '#ff9a3c', label: 'Wnt3a bead (right)' }, { color: '#45c8ff', label: 'Shh-N bead (left)' }]);
    }

    function layout() {
      const cell = params.view === 'cell';
      groupA.visible = cell; groupB.visible = !cell;
      for (const el of sections.cell) el.style.display = cell ? '' : 'none';
      for (const el of sections.sph) el.style.display = cell ? 'none' : '';
      for (const ch of panel.charts) ch.dirty = true;
      if (cell) stage.frame([-6, 0.5, 4.5], 20.5, new THREE.Vector3(0.95, 0.62, 1.0));
      else stage.frame([0, 0, 0], sph.xMax + beadR * 1.6, new THREE.Vector3(0.12, 0.42, 1));
      stage.clipAxis.set(0, 0, 1);
      updateLegend();
    }

    const fmt2 = (x) => fmt(x);
    stage.setPickables([
      ...[['bcat', 'β-catenin', () => `model: cytoplasm ${fmt2(yW[3])} · nucleus ${fmt2(yW[4])} (1 particle = 1/${SCALE.bcat})`],
        ['nicd', 'NICD (Notch intracellular domain)', () => `model: cytoplasm ${fmt2(yN[3])} · nucleus ${fmt2(yN[4])} · Hes1 ${fmt2(yN[5])}`],
        ['gliA', 'GLI activator (Gli2/Gli1)', () => `model: GliA ${fmt2(yH[3])} · read-out φ ${fmt2(stA.readout)}`],
        ['gliR', 'GLI3 repressor (processed by PKA / GSK-3β / CK1 + β-TrCP)', () => `model: GliR ${fmt2(yH[4])} · SMO active ${(100 * yH[1]).toFixed(0)} %`],
        ['wnt', 'Wnt3a ligand', () => `model: free Wnt ${fmt2(yW[0])} · Fz/LRP6 bound ${(100 * clamp01(yW[1])).toFixed(0)} %`],
        ['shh', 'Sonic hedgehog (Shh-N)', () => `Shh ${fmt2(pH.shh)} K_d · PTCH1 bound ${(100 * stA.fb).toFixed(0)} %`],
      ].map(([key, name, model]) => ({
        get object() { return SW[key].pool.mesh; },
        info: (hit) => {
          if (!groupA.visible) return null;
          const sw = SW[key], i = hit.instanceId;
          if (i === undefined) return null;
          if (i >= sw.n) return `${name}\nbound to its receptor\n${model()}`;
          const extra = sw.st[i] === ST.DOOMED ? '\nphosphorylated by CK1α / GSK-3β → β-TrCP → proteasome' : '';
          return `${name}\n${ST_TEXT[sw.st[i]]}${extra}\n${model()}`;
        },
      })),
      { get object() { return fzPool.mesh; }, info: () => groupA.visible && `Frizzled / LRP6 co-receptor\nbound fraction ${(100 * clamp01(yW[1])).toFixed(0)} % · Dvl active ${fmt2(yW[2])}\nWnt-bound LRP6 recruits Dvl and Axin (signalosome)` },
      { get object() { return dcPool.mesh; }, info: () => groupA.visible && `Destruction complex (APC · Axin · GSK-3β · CK1α)\nactivity (1+Axin2)/(1+Dvl*/K) = ${fmt2(stA.dcAct)}\n${stA.nDC} complexes (∝ 1 + Axin2 = ${fmt2(1 + yW[5])})\n${yW[2] > 0.3 ? 'recruited to LRP6 signalosomes at the membrane' : 'free in the cytoplasm, degrading β-catenin'}` },
      { get object() { return notchPool.mesh; }, info: () => groupA.visible && `Notch receptor (free)\nsurface Notch ${fmt2(yN[0])} · waiting for Delta in trans` },
      { get object() { return bridgePool.mesh; }, info: () => groupA.visible && `Delta:Notch complex across the cell–cell gap\nDN ${fmt2(yN[1])} · S2 cleavage by ADAM10 → NEXT` },
      { get object() { return nextPool.mesh; }, info: () => groupA.visible && `NEXT (S2-cleaved Notch stub)\nNEXT ${fmt2(yN[2])}${pN.kS3 === 0 ? '\nDAPT blocks γ-secretase: NEXT accumulates' : ' · γ-secretase releases NICD'}` },
      { get object() { return deltaPool.mesh; }, info: () => groupA.visible && `Delta-like ligand on the sender cell\nDelta (trans) = ${fmt2(pN.deltaTrans)}` },
      { get object() { return ptchPool.mesh; }, info: () => groupA.visible && `Patched-1 (12-pass transporter-like receptor)\nPTCH1 total ${fmt2(yH[0])} · bound ${(100 * stA.fb).toFixed(0)} %\nfree PTCH1 sits in the cilium and inhibits SMO catalytically` },
      { get object() { return smoPool.mesh; }, info: () => groupA.visible && `Smoothened (GPCR-like)\nactive ${(100 * yH[1]).toFixed(0)} %${pH.kSa === 0 ? ' · cyclopamine-bound' : ''}\nactive SMO accumulates in the primary cilium` },
      { get object() { return lociPool.mesh; }, info: (hit) => groupA.visible && `${LOCI[hit.instanceId]?.name ?? 'target gene'}\nexpression: Wnt ${fmt2(yW[6])} · Hes1 ${fmt2(yN[5])} · GLI target ${fmt2(yH[5])}` },
      { object: chromatin, info: () => groupA.visible && `Nucleus (chromatin)\nTCF target ${fmt2(yW[6])} · Hes1 ${fmt2(yN[5])} · GLI target ${fmt2(yH[5])}\nthe envelope glows with target-gene activity` },
      { object: porePool.mesh, info: () => groupA.visible && `Nuclear pore complex\nimport shaped by measured envelope-enrichment kinetics (τ = ${fmt2(tauFrames)} frames)` },
      { object: cilium, info: () => groupA.visible && `Primary cilium · Hedgehog signalling compartment\nSMO* ${(100 * yH[1]).toFixed(0)} % · free PTCH1 in cilium ${fmt2(yH[0] * (1 - stA.fb))}\nGLI proteins are activated at the tip and travel to the nucleus` },
      { object: sender, info: () => groupA.visible && `Neighbouring sender cell\npresents Delta in trans = ${fmt2(pN.deltaTrans)}` },
      {
        get object() { return bodyPool.mesh; },
        info: (hit) => {
          if (!groupB.visible) return null;
          const i = hit.instanceId;
          if (i === undefined || i >= SPH_N) return null;
          const dom = NEURAL_TUBE_DOMAINS[frenchFlag(readoutB[i])].name;
          return `Cell ${i} · ${sph.nbrs[i].length} contacts · ${FATES[fateOf(i)].name}\n` +
            `Notch ${fmt2(notchS.N[i])} · Delta ${fmt2(notchS.D[i])} → ${notchS.D[i] > SENDER ? 'sender' : 'receiver'}\n` +
            `Wnt ${fmt2(u[i])} → nuclear β-catenin ${fmt2(yWc[i * 7 + 4])}\n` +
            `Shh ${fmt2(v[i])} K_d → GLI φ ${fmt2(readoutB[i])} → ${dom}`;
        },
      },
      { object: beadWnt, info: () => groupB.visible && `Wnt3a-soaked bead\nholds Wnt at c0 = ${fmt2(pB.wnt)} on touching cells; λ_W = ${fmt2(pB.lamW * sph.spacing)} µm` },
      { object: beadShh, info: () => groupB.visible && `Shh-N-soaked bead\nholds Shh at ${fmt2(pB.shh)} K_d on touching cells; λ_S = ${fmt2(pB.lamS * sph.spacing)} µm` },
    ]);

    // ---------------------------------------------------------------- lifecycle
    function resetA() {
      resetModelsA();
      populateA();
      paperChart.clear(); wntChart.clear(); notchChart.clear(); hhChart.clear();
    }
    resetA();
    resetSpheroid();
    stepSpheroid(0.001);
    drawSpheroid(0);
    layout();

    function applyMode(mode) {
      // confocal: the envelope shows a lamin/NPC stain instead of the activity tint;
      // confocal and H&E show the spheroid as a physical/optical section
      envMat.uniforms.uUseInst.value = mode === 'confocal' ? 0 : 1;
      setSlab(mode !== 'physical');
    }
    applyMode(stage.mode);

    let readTimer = 0;
    return {
      update(dt) {
        if (params.view === 'cell') {
          updateCell(dt);
          if (tA - lastSampleA >= 0.25 || lastSampleA < 0) { lastSampleA = tA; sampleCharts(); }
          setStatus(`t = ${tA.toFixed(1)} τ · β-cat ${fmt(yW[4])} · NICD ${fmt(yN[4])} · GLI ${fmt(stA.readout)}`);
        } else {
          if (dt > 0) stepSpheroid(dt * TIME_B);
          drawSpheroid(dt);
          readTimer += dt;
          let stats = null;
          if (tB - lastSampleB >= 0.5 || lastSampleB < 0) { lastSampleB = tB; stats = sampleSpheroid(); }
          if (tB - lastProfileB >= 1 || lastProfileB < 0) { lastProfileB = tB; profiles(); }
          if (stats || readTimer > 0.5) {
            readTimer = 0;
            const sf = stats?.sf ?? collier.senderFraction(notchS, SENDER), ss = stats?.ss ?? collier.senderAdjacency(notchS, sph.nbrs, SENDER);
            setStatus(`t = ${tB.toFixed(1)} τ · ${SPH_N} cells · senders ${(100 * sf).toFixed(0)} % · S–S ${(100 * ss).toFixed(1)} %`);
          }
        }
      },
      reset() {
        if (params.view === 'cell') resetA(); else resetSpheroid();
      },
      onMode: applyMode,
    };
  },
};

/**
 * Packed spheroid: jittered FCC positions inside a ball, relaxed with the
 * centre-based mechanics of agents3d (growth off), then its contact graph.
 */
function buildSpheroid(N, r0) {
  const rng = new RNG(20240917);
  const a = r0 * Math.SQRT2 * 0.97;
  const pts = [];
  for (let i = -8; i <= 8; i++) for (let j = -8; j <= 8; j++) for (let k = -8; k <= 8; k++) {
    if ((i + j + k) & 1) continue;
    pts.push([i * a + rng.normal(0, 0.25), j * a + rng.normal(0, 0.25), k * a + rng.normal(0, 0.25)]);
  }
  pts.sort((p, q) => Math.hypot(...p) - Math.hypot(...q));
  pts.length = N;
  const tissue = new Tissue({ noise: 0.05, growth: false, adhesion: () => 1.0, muRep: 10, rng, hashSize: 2 * r0 * 1.3 });
  for (const p of pts) tissue.add(new Cell({ p, R: r0 * (0.93 + 0.14 * rng.next()) }));
  for (let s = 0; s < 60; s++) tissue.step(0.05);
  tissue.noise = 0;
  tissue.step(0.01);
  const cells = tissue.cells;
  const c = tissue.centroid();
  const index = new Map(cells.map((cell, i) => [cell.id, i]));
  const pos = cells.map((cell) => cell.p.map((v, k) => v - c[k]));
  const nbrs = cells.map((cell) => (tissue.neighbours.get(cell.id) ?? []).map((id) => index.get(id)).filter((j) => j !== undefined));
  let spacing = 0, m = 0, maxDeg = 0, degSum = 0;
  nbrs.forEach((nb, i) => {
    maxDeg = Math.max(maxDeg, nb.length); degSum += nb.length;
    for (const j of nb) { spacing += Math.hypot(pos[i][0] - pos[j][0], pos[i][1] - pos[j][1], pos[i][2] - pos[j][2]); m++; }
  });
  const xs = pos.map((p) => p[0]);
  return {
    pos, nbrs, R: cells.map((cell) => cell.R),
    axis: cells.map(() => rng.onSphere()),
    spacing: spacing / Math.max(m, 1), maxDeg, meanDeg: degSum / N,
    xMin: Math.min(...xs) - r0, xMax: Math.max(...xs) + r0,
  };
}
