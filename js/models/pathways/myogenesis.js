// Myogenesis (paper: "Third order - Myogenesis: activation of MyoD and Myf5,
// then myogenin and MRF4, then muscle-specific genes essential for muscle
// fibre formation").
//
// 1) mrf: the myogenic regulatory factor (MRF) cascade in one nucleus.
//    The hierarchy is established (Bentzinger, Wang & Rudnicki 2012, Cold
//    Spring Harb Perspect Biol 4:a008342; Buckingham & Rigby 2014, Dev Cell
//    28:225):
//      Wnt (+ Shh)  -> Myf5, MyoD                (determination; Munsterberg
//                                                  et al. 1995 Genes Dev 9:2911;
//                                                  Gli on the Myf5 epaxial enhancer,
//                                                  Gustafsson et al. 2002 Genes Dev 16:114)
//      MyoD         -> MyoD                       (positive autoregulation;
//                                                  Thayer et al. 1989 Cell 58:241)
//      MyoD         -> myogenin  -| mitogens      (Id sequestration of E proteins,
//                                                  Benezra et al. 1990 Cell 61:49;
//                                                  MyoD -> p21 cell-cycle exit,
//                                                  Halevy et al. 1995 Science 267:1018)
//      myogenin     -> MRF4, myomaker, MHC        (differentiation / fusion;
//                                                  Millay et al. 2013 Nature 499:301)
//    The KINETIC FORM is phenomenological: normalised concentrations (0..~1),
//    Hill functions and first-order turnover (units: hours), chosen so that
//    the qualitative timing of C2C12 differentiation is reproduced
//    (myogenin within ~12 h of serum withdrawal, MHC-positive myotubes at
//    48-72 h). With v_auto = 0.9 the MyoD autoregulation makes commitment
//    irreversible (the high state survives removal of Wnt/Shh); the fold
//    points of the switch are computed in closed form by mrf.thresholds().
//
// 2) MyoCulture: agent model of myoblasts on a 2D substrate (units um, h).
//    Myoblasts are self-propelled bipolar rods (they glide along their long
//    axis and reverse polarity) with steric repulsion between spherocylinders
//    and nematic alignment torque  dtheta_i/dt = gamma sum_j sin 2(theta_j - theta_i)
//    (Kemkemer et al. 2000 Eur Phys J E 3:101; Duclos et al. 2017 Nat Phys
//    13:58). Myogenin-high (fusion-competent) cells that touch an aligned
//    partner end-to-end or side-by-side fuse with a constant hazard; volume
//    and nuclei are conserved. Myotube radius vs nucleus number, nuclear
//    spreading and fusion rate are phenomenological.
import { hill, hillR } from '../core/ode.js';
import { RNG } from '../core/rng.js';

// ------------------------------------------------------------------ MRF ODE
export const mrf = {
  species: ['Myf5', 'MyoD', 'myogenin', 'MRF4', 'MHC'],
  labels: ['Myf5', 'MyoD', 'myogenin', 'MRF4', 'MHC (myosin heavy chain)'],
  defaults: {
    kF: 0.3,               // Myf5 turnover (1/h)
    kD: 0.35,              // MyoD turnover
    kG: 0.04,              // myogenin
    kR: 0.05,              // MRF4
    kM: 0.035,             // MHC (structural protein, slow)
    Kw: 0.5, Ks: 0.5,      // half-activation of the Wnt / Shh inputs
    aE: 0.3, aF: 0.25,     // weight of environmental signal and Myf5 in MyoD induction
    vAuto: 0.9, KAuto: 0.5, nAuto: 4, // MyoD autoregulation
    KDG: 0.45,             // MyoD needed for myogenin
    Kgf: 0.25,             // growth factor that halves myogenin induction
    KGR: 0.4,              // myogenin -> MRF4
    KM: 0.8, nM: 3,        // (myogenin + MRF4) -> MHC and muscle genes
    KGF5: 0.5,             // Myf5 down-regulation on terminal differentiation
    fuseThreshold: 0.4,    // myogenin level for fusion competence (myomaker/myomixer)
  },

  /** Environmental myogenic signal E in [0, 1): Wnt required, Shh synergises. */
  environment(wnt, shh, p = mrf.defaults) {
    return hill(wnt, p.Kw, 2) * (0.3 + 0.7 * hill(shh, p.Ks, 2));
  },

  /** Input to the MyoD switch from the environment and Myf5. */
  switchInput(E, F, p = mrf.defaults) { return p.aE * E + p.aF * F; },

  /**
   * Right-hand side for one nucleus.
   * @param {object} p parameters
   * @param {{wnt:number, shh:number, gf:number}} env  (env may be mutated between steps)
   * @param {number[]} [k] per-cell rate multipliers (heterogeneity), length 5
   * @param {number} [gGain] per-cell gain of myogenin induction (reserve cells ~0.1)
   */
  rhs(p, env, k = null, gGain = 1) {
    return (t, y, dy) => {
      const [F, D, G, R, M] = y;
      const E = mrf.environment(env.wnt, env.shh, p);
      const s = mrf.switchInput(E, F, p);
      dy[0] = p.kF * (E * hillR(G, p.KGF5, 2) - F);
      dy[1] = p.kD * (s + p.vAuto * hill(D, p.KAuto, p.nAuto) - D);
      dy[2] = p.kG * (gGain * hill(D, p.KDG, 2) * hillR(env.gf, p.Kgf, 2) - G);
      dy[3] = p.kR * (hill(G, p.KGR, 2) - R);
      dy[4] = p.kM * (hill(G + R, p.KM, p.nM) - M);
      if (k) for (let i = 0; i < 5; i++) dy[i] *= k[i];
    };
  },

  /**
   * Fold points of the MyoD switch D = s + v hill(D, K, n)  <=>  s = g(D) = D - v hill(D).
   * Turn-on threshold = local maximum of g (lower fold), turn-off = local minimum
   * (upper fold). A negative turn-off threshold means the committed state is
   * irreversible for every physical input s >= 0.
   */
  thresholds(p = mrf.defaults) {
    const g = (D) => D - p.vAuto * hill(D, p.KAuto, p.nAuto);
    const N = 4000, Dmax = 2;
    let on = null, off = null, prev = g(0), prevD = 0, slopePrev = null;
    const folds = [];
    for (let i = 1; i <= N; i++) {
      const D = (i / N) * Dmax, v = g(D);
      const slope = v - prev;
      if (slopePrev !== null && Math.sign(slope) !== Math.sign(slopePrev)) folds.push({ D: prevD, s: prev });
      slopePrev = slope; prev = v; prevD = D;
    }
    if (folds.length >= 2) { on = folds[0].s; off = folds[1].s; }
    return { on, off, bistable: folds.length >= 2, folds };
  },

  /** MyoD steady state reached from D0 with a constant switch input s (for hysteresis sweeps). */
  relaxMyoD(s, D0, p = mrf.defaults, T = 400, h = 0.05) {
    let D = D0;
    for (let t = 0; t < T; t += h) {
      // RK2 on the scalar switch (fast, accurate enough for a steady state)
      const f = (x) => p.kD * (s + p.vAuto * hill(x, p.KAuto, p.nAuto) - x);
      const k1 = f(D), k2 = f(D + 0.5 * h * k1);
      D += h * k2;
    }
    return D;
  },

  /** Quasi-static up- and down-sweep of the input; returns the measured jump points. */
  hysteresisSweep(p = mrf.defaults, { sMax = 0.6, ds = 0.0025, level = 0.45 } = {}) {
    let D = 0, on = null, off = null;
    const up = [], down = [];
    for (let s = 0; s <= sMax + 1e-12; s += ds) {
      D = mrf.relaxMyoD(s, D, p);
      up.push([s, D]);
      if (on === null && D > level) on = s;
    }
    for (let s = sMax; s >= -1e-12; s -= ds) {
      D = mrf.relaxMyoD(s, D, p);
      down.push([s, D]);
      if (off === null && D < level) off = s;
    }
    return { on, off, up, down };
  },

  /** Stage label from a state vector. */
  stage(y, p = mrf.defaults) {
    if (y[4] > 0.5) return 'MHC+ (differentiated)';
    if (y[2] > p.fuseThreshold) return 'myogenin+ (fusion-competent)';
    if (y[1] > 0.5) return 'MyoD+ (committed myoblast)';
    if (y[0] > 0.2) return 'Myf5+ (specified)';
    return 'uncommitted mesoderm';
  },
};

// ---------------------------------------------------------- geometry helpers
const TAU = Math.PI * 2;

/** Closest points between 2D segments [p1,q1] and [p2,q2] (Ericson, RTCD 5.1.9). */
export function segSeg2D(p1x, p1z, q1x, q1z, p2x, p2z, q2x, q2z, out) {
  const d1x = q1x - p1x, d1z = q1z - p1z, d2x = q2x - p2x, d2z = q2z - p2z;
  const rx = p1x - p2x, rz = p1z - p2z;
  const a = d1x * d1x + d1z * d1z, e = d2x * d2x + d2z * d2z, f = d2x * rx + d2z * rz;
  let s, t;
  if (a <= 1e-12 && e <= 1e-12) { s = t = 0; }
  else if (a <= 1e-12) { s = 0; t = clamp01(f / e); }
  else {
    const c = d1x * rx + d1z * rz;
    if (e <= 1e-12) { t = 0; s = clamp01(-c / a); }
    else {
      const b = d1x * d2x + d1z * d2z, den = a * e - b * b;
      s = den > 1e-12 ? clamp01((b * f - c * e) / den) : 0;
      t = (b * s + f) / e;
      if (t < 0) { t = 0; s = clamp01(-c / a); } else if (t > 1) { t = 1; s = clamp01((b - c) / a); }
    }
  }
  out.s = s; out.t = t;
  out.ax = p1x + d1x * s; out.az = p1z + d1z * s;
  out.bx = p2x + d2x * t; out.bz = p2z + d2z * t;
  out.d = Math.hypot(out.bx - out.ax, out.bz - out.az);
  return out;
}
const clamp01 = (x) => (x < 0 ? 0 : x > 1 ? 1 : x);

// ---------------------------------------------------------- culture model
export const CULTURE_DEFAULTS = {
  dishR: 165,          // um, radius of the culture area
  L: 44, W: 10, H: 3.6, // myoblast spindle: length, width, height (um)
  v0: 18,              // gliding speed along the long axis (um/h)
  reversal: 0.35,      // polarity reversal rate (1/h)
  kRep: 9,             // steric repulsion (1/h)
  align: 0.9,          // nematic alignment gain gamma (1/h)
  Dr: 0.04,            // rotational diffusion (rad^2/h)
  kFuse: 0.3,          // fusion hazard for a competent, aligned, touching pair (1/h)
  fuseCos: Math.cos((32 * Math.PI) / 180),
  tubeFuseRate: 0.15,  // end-to-end myotube-myotube fusion hazard (1/h)
  cycle: 18,           // myoblast doubling time in growth medium (h)
  maxCells: 120,
  grooves: 0,          // strength of contact guidance along x (micro-grooved substrate)
  hetero: 0.18,        // cell-to-cell variability of MRF rates
  // "Reserve cells": a quiescent subpopulation that does not induce myogenin in
  // differentiation medium and stays mononucleated (Yoshida et al. 1998, J Cell
  // Sci 111:769). Fraction and 10x lower myogenin gain are phenomenological.
  reserve: 0.25,
};

/** Myoblast volume: flattened ellipsoid with semi-axes L/2, W/2, H. */
export const myoblastVolume = (c) => (4 / 3) * Math.PI * (c.L / 2) * (c.W / 2) * c.H;

/** Myotube: elliptic cross-section (half-width r, half-height 0.7 r) with tapered hemi-ellipsoid ends of length 2r. */
export const TUBE_FLAT = 0.7;
export function tubeRadius(nuclei, W = CULTURE_DEFAULTS.W) {
  return Math.min((W / 2) * (1 + 0.24 * Math.log(Math.max(nuclei, 1))), 12);
}
export function tubeLength(vol, r) {
  const A = Math.PI * r * r * TUBE_FLAT;
  return vol / A + (4 / 3) * r; // L_body + 2 cap = V/A - (4/3) cap + 2 cap, cap = 2r
}
export const tubeCap = (t) => Math.min(2 * t.r, t.len / 2);

let uid = 1;

export class MyoCulture {
  /**
   * @param {object} o
   * @param {number} [o.n] initial myoblasts
   * @param {'committed'|'naive'} [o.start] C2C12-like committed myoblasts or naive mesoderm
   * @param {number} [o.seed]
   */
  constructor(o = {}) {
    this.p = { ...CULTURE_DEFAULTS, ...(o.culture ?? {}) };
    this.mp = { ...mrf.defaults, ...(o.mrf ?? {}) };
    this.env = { wnt: o.wnt ?? 0, shh: o.shh ?? 0, gf: o.gf ?? 0.05 };
    this.rng = new RNG(o.seed ?? 12345);
    this.cells = [];
    this.tubes = [];
    this.t = 0;
    this.fusions = 0;
    this.divisions = 0;
    this.events = [];   // recent fusion events {x, z, t} for rendering flashes
    this._cp = {};
    this._dy = new Float64Array(5);
    this._k1 = new Float64Array(5); this._k2 = new Float64Array(5); this._tmp = new Float64Array(5);
    const n = o.n ?? 90;
    const start = o.start ?? 'committed';
    const R = this.p.dishR - this.p.L / 2;
    for (let i = 0; i < n; i++) {
      let x, z, tries = 0;
      do {
        const r = R * Math.sqrt(this.rng.next()), a = this.rng.next() * TAU;
        x = r * Math.cos(a); z = r * Math.sin(a);
        tries++;
      } while (tries < 30 && this.cells.some((c) => Math.hypot(c.x - x, c.z - z) < this.p.W * 1.6));
      const s = new Float64Array(5);
      if (start === 'committed') {
        s[0] = 0.25 + 0.1 * this.rng.next();
        s[1] = 0.8 + 0.1 * this.rng.next();
      } else {
        s[0] = 0.02 * this.rng.next();
        s[1] = 0.03 * this.rng.next();
      }
      this.cells.push(this.newCell(x, z, this.rng.next() * Math.PI, s));
    }
    this.nuclei0 = this.totalNuclei();
    this.volume0 = this.totalVolume();
  }

  newCell(x, z, th, s) {
    const h = this.p.hetero;
    const k = new Float64Array(5);
    for (let i = 0; i < 5; i++) k[i] = Math.max(0.3, 1 + h * this.rng.normal());
    return {
      id: uid++, x, z, th, s, k,
      L: this.p.L * (0.85 + 0.3 * this.rng.next()), W: this.p.W * (0.9 + 0.2 * this.rng.next()), H: this.p.H,
      pol: this.rng.next() < 0.5 ? -1 : 1, fx: 0, fz: 0, tq: 0, al: 0, reserve: this.rng.next() < this.p.reserve,
      age: this.rng.next() * this.p.cycle, nucSpin: this.rng.next() * TAU, born: this.t,
    };
  }

  totalNuclei() { return this.cells.length + this.tubes.reduce((s, t) => s + t.nuclei.length, 0); }
  totalVolume() { return this.cells.reduce((s, c) => s + myoblastVolume(c), 0) + this.tubes.reduce((s, t) => s + t.vol, 0); }

  /** Fraction of nuclei inside multinucleated (>= 2 nuclei) myotubes: the standard fusion index. */
  fusionIndex() {
    const inT = this.tubes.reduce((s, t) => s + t.nuclei.length, 0);
    return inT / Math.max(this.totalNuclei(), 1);
  }

  /** Population mean of each MRF species over nuclei (myotube nuclei share the syncytium's state). */
  meanState(out = new Float64Array(5)) {
    out.fill(0);
    let n = 0;
    for (const c of this.cells) { for (let i = 0; i < 5; i++) out[i] += c.s[i]; n++; }
    for (const t of this.tubes) { const w = t.nuclei.length; for (let i = 0; i < 5; i++) out[i] += w * t.s[i]; n += w; }
    for (let i = 0; i < 5; i++) out[i] /= Math.max(n, 1);
    return out;
  }

  /** Nematic order parameter S = |<exp(2 i theta)>| of myoblasts (1 = perfectly aligned). */
  nematicOrder() {
    let cx = 0, sx = 0, n = 0;
    for (const c of this.cells) { cx += Math.cos(2 * c.th); sx += Math.sin(2 * c.th); n++; }
    for (const t of this.tubes) { const w = t.nuclei.length; cx += w * Math.cos(2 * t.th); sx += w * Math.sin(2 * t.th); n += w; }
    return n ? Math.hypot(cx, sx) / n : 0;
  }

  // integrate one nucleus state (midpoint rule, rates <= 0.35/h so dt <= 0.2 h is accurate)
  integrateState(s, k, dt, gGain = 1) {
    const f = mrf.rhs(this.mp, this.env, k, gGain);
    const { _k1: k1, _tmp: tmp } = this;
    f(0, s, k1);
    for (let i = 0; i < 5; i++) tmp[i] = s[i] + 0.5 * dt * k1[i];
    f(0, tmp, k1);
    for (let i = 0; i < 5; i++) s[i] = Math.max(0, s[i] + dt * k1[i]);
  }

  /** Body axis endpoints of a myoblast or myotube (segment of the spherocylinder). */
  segment(b, out) {
    const tube = !!b.nuclei;
    const half = tube ? Math.max(b.len / 2 - b.r * 0.9, 0) : Math.max(b.L / 2 - b.W / 2, 0);
    const cx = Math.cos(b.th) * half, cz = Math.sin(b.th) * half;
    out[0] = b.x - cx; out[1] = b.z - cz; out[2] = b.x + cx; out[3] = b.z + cz;
    return out;
  }

  radius(b) { return b.nuclei ? b.r : b.W / 2; }

  step(dt) {
    const p = this.p, rng = this.rng;
    const bodies = [...this.cells, ...this.tubes];
    for (const b of bodies) { b.fx = 0; b.fz = 0; b.tq = 0; b.al = 0; }
    const segA = [0, 0, 0, 0], segB = [0, 0, 0, 0], cp = this._cp;
    const fusePairs = [];
    const thr = this.mp.fuseThreshold;
    const n = bodies.length;
    const halfExt = bodies.map((b) => (b.nuclei ? b.len / 2 : b.L / 2) + this.radius(b));
    for (let i = 0; i < n; i++) {
      const a = bodies[i];
      this.segment(a, segA);
      const ra = this.radius(a);
      for (let j = i + 1; j < n; j++) {
        const b = bodies[j];
        const dx = b.x - a.x, dz = b.z - a.z;
        const reach = halfExt[i] + halfExt[j] + 8;
        if (dx * dx + dz * dz > reach * reach) continue;
        this.segment(b, segB);
        segSeg2D(segA[0], segA[1], segA[2], segA[3], segB[0], segB[1], segB[2], segB[3], cp);
        const rb = this.radius(b);
        const contact = ra + rb;
        const d = cp.d || 1e-6;
        const cosD = Math.cos(a.th - b.th);
        if (d < contact) {
          // steric repulsion applied at the contact points (produces aligning torques)
          const nx = (cp.bx - cp.ax) / d, nz = (cp.bz - cp.az) / d;
          const F = p.kRep * (contact - d);
          a.fx -= F * nx; a.fz -= F * nz; b.fx += F * nx; b.fz += F * nz;
          a.tq += (cp.ax - a.x) * -F * nz - (cp.az - a.z) * -F * nx;
          b.tq += (cp.bx - b.x) * F * nz - (cp.bz - b.z) * F * nx;
        }
        if (d < contact + 6) {
          // nematic alignment with touching neighbours
          const w = p.align * (1 - Math.max(d - contact, 0) / 6);
          const s2 = Math.sin(2 * (b.th - a.th));
          const la = a.nuclei ? 0.05 : 1, lb = b.nuclei ? 0.05 : 1; // myotubes barely rotate
          a.al += w * s2 * la; b.al -= w * s2 * lb;
        }
        // fusion candidates: touching and aligned
        if (d < contact * 1.2 + 1 && Math.abs(cosD) > p.fuseCos) {
          const aC = a.nuclei ? true : a.s[2] > thr;
          const bC = b.nuclei ? true : b.s[2] > thr;
          if (aC && bC && !(a.nuclei && b.nuclei)) fusePairs.push([a, b]);
          else if (a.nuclei && b.nuclei && (cp.s < 0.05 || cp.s > 0.95) && (cp.t < 0.05 || cp.t > 0.95)) fusePairs.push([a, b, true]);
        }
      }
    }

    // move: self-propulsion + forces, rotation: steric torque + alignment + noise
    const R = p.dishR;
    for (const b of bodies) {
      const tube = !!b.nuclei;
      const len = tube ? b.len : b.L;
      const drag = tube ? 1 + len / 25 : 1;
      let vx = b.fx / drag, vz = b.fz / drag;
      if (!tube) { vx += p.v0 * b.pol * Math.cos(b.th); vz += p.v0 * b.pol * Math.sin(b.th); }
      // steric torque over rotational drag of a rod (~ L^2 / 12) + nematic alignment rate
      const rotDrag = (len * len) / 12 * (tube ? 3 : 1);
      let om = b.tq / rotDrag + b.al;
      if (p.grooves > 0) om += p.grooves * -Math.sin(2 * b.th) * (tube ? 0.1 : 1);
      b.th += om * dt + (tube ? 0.2 : 1) * Math.sqrt(2 * p.Dr * dt) * rng.normal();
      b.x += vx * dt; b.z += vz * dt;
      // dish wall: keep both ends inside
      const half = len / 2;
      for (const sg of [-1, 1]) {
        const ex = b.x + sg * half * Math.cos(b.th), ez = b.z + sg * half * Math.sin(b.th);
        const r = Math.hypot(ex, ez), over = r - (R - this.radius(b));
        if (over > 0) { b.x -= (over * ex) / r * 0.5; b.z -= (over * ez) / r * 0.5; if (!tube) b.pol = -sg * Math.sign(Math.cos(b.th) * ex + Math.sin(b.th) * ez) || -b.pol; }
      }
      if (!tube && rng.next() < p.reversal * dt) b.pol = -b.pol;
    }

    // gene regulation (nuclei in a myotube share its cytoplasm -> one state per syncytium)
    for (const c of this.cells) this.integrateState(c.s, c.k, dt, c.reserve ? 0.1 : 1);
    for (const t of this.tubes) this.integrateState(t.s, null, dt);

    // proliferation in growth medium (cell-cycle exit once myogenin rises)
    if (this.cells.length + this.tubes.length < p.maxCells) {
      const born = [];
      for (const c of this.cells) {
        const rate = (Math.LN2 / p.cycle) * hill(this.env.gf, 0.3, 2) * hillR(c.s[2], 0.15, 4);
        if (rng.next() < rate * dt && this.cells.length + born.length < p.maxCells) born.push(this.divide(c));
      }
      this.cells.push(...born);
    }

    // fusion events (each body fuses at most once per step)
    const used = new Set();
    for (const [a, b, tubeTube] of fusePairs) {
      if (used.has(a) || used.has(b)) continue;
      // secondary fusion of myoblasts with nascent myotubes is more efficient than
      // primary myoblast-myoblast fusion (factor 2: phenomenological)
      const rate = tubeTube ? p.tubeFuseRate : p.kFuse * (a.nuclei || b.nuclei ? 2 : 1);
      if (rng.next() >= 1 - Math.exp(-rate * dt)) continue;
      if (tubeTube && (a.len + b.len) > 1.7 * p.dishR) continue;
      used.add(a); used.add(b);
      this.fuse(a, b);
    }

    // myotube growth to its volume-determined length, nuclear spreading
    for (const t of this.tubes) {
      t.len += (t.lenTarget - t.len) * Math.min(1, dt * 0.6);
      t.r += (tubeRadius(t.nuclei.length, this.p.W) - t.r) * Math.min(1, dt * 0.4);
      t.lenTarget = tubeLength(t.vol, t.r);
      const N = t.nuclei.length;
      t.nuclei.sort((u, v) => u.u - v.u);
      const kSpread = 0.04 + 0.35 * t.s[4];
      for (let k = 0; k < N; k++) {
        const target = N === 1 ? 0 : -0.8 + (1.6 * (k + 0.5)) / N;
        t.nuclei[k].u += (target - t.nuclei[k].u) * Math.min(1, kSpread * dt);
      }
      t.age += dt;
    }
    this.events = this.events.filter((e) => this.t - e.t < 1.5);
    this.t += dt;
  }

  divide(c) {
    const off = c.L / 4;
    const d = this.newCell(c.x + off * Math.cos(c.th), c.z + off * Math.sin(c.th), c.th + 0.1 * this.rng.normal(), Float64Array.from(c.s));
    c.x -= off * Math.cos(c.th); c.z -= off * Math.sin(c.th);
    d.L = c.L; d.W = c.W;
    this.divisions++;
    return d;
  }

  /** Fuse a with b (myoblast/myoblast, myoblast/myotube or myotube/myotube). Conserves volume and nuclei. */
  fuse(a, b) {
    const parts = [a, b];
    const vol = (x) => (x.nuclei ? x.vol : myoblastVolume(x));
    const nuc = (x) => (x.nuclei ? x.nuclei.length : 1);
    const V = vol(a) + vol(b);
    const Na = nuc(a), Nb = nuc(b);
    // new axis: nucleus-weighted nematic average of the two directors
    let ca = Math.cos(a.th), sa = Math.sin(a.th), cb = Math.cos(b.th), sb = Math.sin(b.th);
    if (ca * cb + sa * sb < 0) { cb = -cb; sb = -sb; }
    const wa = Na * (a.nuclei ? 3 : 1), wb = Nb * (b.nuclei ? 3 : 1);
    const th = Math.atan2(wa * sa + wb * sb, wa * ca + wb * cb);
    const ux = Math.cos(th), uz = Math.sin(th);
    // span of both bodies projected on the new axis
    const cx = (a.x * vol(a) + b.x * vol(b)) / V, cz = (a.z * vol(a) + b.z * vol(b)) / V;
    let lo = Infinity, hi = -Infinity;
    const absNuclei = [];
    for (const x of parts) {
      const half = (x.nuclei ? x.len : x.L) / 2;
      for (const sg of [-1, 1]) {
        const ex = x.x + sg * half * Math.cos(x.th) - cx, ez = x.z + sg * half * Math.sin(x.th) - cz;
        const pr = ex * ux + ez * uz;
        lo = Math.min(lo, pr); hi = Math.max(hi, pr);
      }
      if (x.nuclei) {
        const usable = Math.max(x.len / 2 - tubeCap(x) * 0.6, 1);
        for (const nu of x.nuclei) {
          const px = x.x + Math.cos(x.th) * nu.u * usable, pz = x.z + Math.sin(x.th) * nu.u * usable;
          absNuclei.push({ ...nu, pr: (px - cx) * ux + (pz - cz) * uz });
        }
      } else {
        absNuclei.push({ id: x.id, lat: 0.5 * (this.rng.next() - 0.5), h: 0.2 * (this.rng.next() - 0.5), spin: x.nucSpin, pr: (x.x - cx) * ux + (x.z - cz) * uz, born: this.t });
      }
    }
    const mid = (lo + hi) / 2;
    const s = new Float64Array(5);
    for (let i = 0; i < 5; i++) s[i] = (Na * a.s[i] + Nb * b.s[i]) / (Na + Nb);
    const N = Na + Nb;
    const r = Math.max(a.nuclei ? a.r : a.W / 2, b.nuclei ? b.r : b.W / 2);
    const tube = {
      id: a.nuclei ? a.id : b.nuclei ? b.id : uid++,
      x: cx + ux * mid, z: cz + uz * mid, th, r, vol: V, s,
      len: Math.max(hi - lo, 2 * r), lenTarget: 0, nuclei: [], age: a.nuclei ? a.age : b.nuclei ? b.age : 0,
      fx: 0, fz: 0, tq: 0, al: 0,
    };
    tube.lenTarget = tubeLength(V, tubeRadius(N, this.p.W));
    const usable = Math.max(tube.len / 2 - tubeCap(tube) * 0.6, 1);
    for (const nu of absNuclei) {
      const u = Math.max(-0.95, Math.min(0.95, (nu.pr - mid) / usable));
      tube.nuclei.push({ id: nu.id, u, lat: nu.lat, h: nu.h, spin: nu.spin, born: nu.born });
    }
    this.cells = this.cells.filter((c) => c !== a && c !== b);
    this.tubes = this.tubes.filter((t) => t !== a && t !== b);
    this.tubes.push(tube);
    this.fusions++;
    this.events.push({ x: (a.x + b.x) / 2, z: (a.z + b.z) / 2, t: this.t });
    return tube;
  }
}

/** Reaction list for the network view (order 3 = differentiated cells). */
export const REACTIONS = [
  { id: 'myo.myf5', order: 3, pathway: 'Myogenesis', reactants: [], products: ['Myf5'], modifiers: ['Wnt', 'Shh', 'beta-cat:TCF', 'GliA (nuc)'], label: 'Myf5 induction by Wnt/Shh (somite)' },
  { id: 'myo.myod', order: 3, pathway: 'Myogenesis', reactants: [], products: ['MyoD'], modifiers: ['Wnt', 'Myf5'], label: 'MyoD activation (determination)' },
  { id: 'myo.myod.auto', order: 3, pathway: 'Myogenesis', reactants: [], products: ['MyoD'], modifiers: ['MyoD'], label: 'MyoD positive autoregulation (commitment switch)' },
  { id: 'myo.myod.e', order: 3, pathway: 'Myogenesis', reactants: ['MyoD', 'E-protein'], products: ['MyoD:E-protein'], label: 'Heterodimerisation with E12/E47 on E-boxes' },
  { id: 'myo.id', order: 3, pathway: 'Myogenesis', reactants: ['E-protein', 'Id'], products: ['E-protein:Id'], modifiers: ['growth factor (FGF)', 'ERK-PP (nuc)'], label: 'Mitogen-induced Id sequesters E proteins' },
  { id: 'myo.p21', order: 3, pathway: 'Myogenesis', reactants: [], products: ['p21'], modifiers: ['MyoD'], label: 'p21 induction: cell-cycle exit' },
  { id: 'myo.myog', order: 3, pathway: 'Myogenesis', reactants: [], products: ['myogenin'], modifiers: ['MyoD:E-protein', 'MEF2'], label: 'Myogenin transcription (differentiation)' },
  { id: 'myo.mrf4', order: 3, pathway: 'Myogenesis', reactants: [], products: ['MRF4'], modifiers: ['myogenin'], label: 'MRF4 expression' },
  { id: 'myo.fusogen', order: 3, pathway: 'Myogenesis', reactants: [], products: ['myomaker'], modifiers: ['myogenin', 'MyoD'], label: 'Fusogen (myomaker/myomixer) expression' },
  { id: 'myo.mhc', order: 3, pathway: 'Myogenesis', reactants: [], products: ['MHC'], modifiers: ['myogenin', 'MRF4', 'MEF2'], label: 'Muscle-specific genes (MHC, actin, creatine kinase)' },
  { id: 'myo.fuse', order: 3, pathway: 'Myogenesis', reactants: ['myoblast', 'myoblast'], products: ['myotube'], modifiers: ['myomaker'], label: 'Myoblast fusion' },
  { id: 'myo.fuse2', order: 3, pathway: 'Myogenesis', reactants: ['myoblast', 'myotube'], products: ['myotube'], modifiers: ['myomaker'], label: 'Secondary fusion into a myotube' },
  { id: 'myo.sarc', order: 3, pathway: 'Myogenesis', reactants: ['MHC', 'actin'], products: ['sarcomere'], label: 'Sarcomere assembly (striations)' },
];
