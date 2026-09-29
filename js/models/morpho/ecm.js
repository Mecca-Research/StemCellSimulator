// Extracellular matrix and cell migration (paper: "Higher-Level Functional
// Interactions - ECM interactions: cells interact with the ECM to maintain
// structure and function"; "Mechanical forces: cell migration and adhesion,
// Cell Movement = Chemotaxis + Haptotaxis"; collagen pathway: procollagen is
// secreted, propeptides are cleaved -> tropocollagen -> collagen fibrils).
//
// 1) Collagen fibrillogenesis, mass-action (per fibroblast, arbitrary
//    triple-helix units, time in hours):
//       d[P]/dt = sigma - kc [P]                        secretion; N-/C-propeptide
//                                                        cleavage (ADAMTS-2, BMP-1)
//       d[T]/dt = kc [P] - 2 kn [T]^2 - kg [T][F]         tropocollagen
//       d[F]/dt = 2 kn [T]^2 + kg [T][F] - kdeg [F]       fibril: nucleation + elongation
//       d[Pp]/dt = kc [P]                               released propeptide pairs
//    P + T + F is conserved when sigma = kdeg = 0. The nucleation-elongation
//    form is a phenomenological mass-action scheme that reproduces the lag-
//    then-growth kinetics of fibril assembly after C-propeptide cleavage
//    (Kadler, Hojima & Prockop 1987, J Biol Chem 262:15696; review Kadler et al.
//    1996, Biochem J 316:1). Fibrils are laid down parallel to the cell's long
//    axis / direction of motion, as in fibripositors (Canty et al. 2004,
//    J Cell Biol 165:553).
//
// 2) Chemotaxis + haptotaxis, the individual-based form of Anderson & Chaplain
//    1998 (Bull Math Biol 60:857):
//       dx = [ chi / (1 + alpha C) grad C  +  rho grad E ] dt + sqrt(2 Dn dt) xi
//    C = chemoattractant (steady point source, gradient.js), E = ECM (fibril)
//    density. Migrating cells degrade the matrix they cross (MMPs), as in the
//    same model's dE/dt = -mu n E term.
//
// 3) Fibroblast motility: persistent random walk (Ornstein-Uhlenbeck velocity,
//    speed S, persistence time P; Dunn & Brown 1987 J Cell Sci Suppl 8:81;
//    Selmeczi et al. 2005 Biophys J 89:912).
//
// Pure JS (no DOM / three.js), deterministic through core/rng.js.
import { rk4Step, rk4Work } from '../core/ode.js';
import { RNG } from '../core/rng.js';
import { pointSource3D, decayLength } from './gradient.js';

export const collagen = {
  species: ['Procollagen', 'Tropocollagen', 'Collagen fibril', 'Propeptides'],
  defaults: { sigma: 1.0, kc: 1.2, kn: 0.3, kg: 3, kdeg: 0 },
  rhs(p) {
    return (t, y, dy) => {
      const [P, T, F] = y;
      const cleave = p.kc * P, nuc = p.kn * T * T, grow = p.kg * T * F;
      dy[0] = p.sigma - cleave;
      dy[1] = cleave - 2 * nuc - grow;
      dy[2] = 2 * nuc + grow - p.kdeg * F;
      dy[3] = cleave;
    };
  },
  initial() { return [0, 0, 0, 0]; },
  /** Collagen mass in triple-helix units (propeptides excluded). */
  mass(y) { return y[0] + y[1] + y[2]; },
};

/** Anderson-Chaplain chemotactic sensitivity with receptor saturation chi(C) = chi / (1 + alpha C). */
export function chemotacticVelocity(C, gradC, { chi, alpha = 0 }, out = [0, 0, 0]) {
  const s = chi / (1 + alpha * C);
  out[0] = s * gradC[0]; out[1] = s * gradC[1]; out[2] = s * gradC[2];
  return out;
}

/** Normalised steady chemoattractant field of a spherical source (C = 1 at its surface). */
export class PointSourceField {
  constructor({ p = [0, 0, 0], radius = 10, D = 100, k = 0.004 }) {
    Object.assign(this, { p, radius, D, k });
    this.lambda = decayLength(D, k);
    this.c0 = pointSource3D(radius, { Q: 1, D, k });
  }

  value(x) {
    const r = Math.hypot(x[0] - this.p[0], x[1] - this.p[1], x[2] - this.p[2]);
    return pointSource3D(Math.max(r, this.radius), { Q: 1, D: this.D, k: this.k }) / this.c0;
  }

  /** grad C = dC/dr r_hat, with dC/dr = -C (1/r + 1/lambda). */
  gradient(x, out = [0, 0, 0]) {
    const dx = x[0] - this.p[0], dy = x[1] - this.p[1], dz = x[2] - this.p[2];
    const r = Math.hypot(dx, dy, dz) || 1e-9;
    if (r < this.radius) { out[0] = out[1] = out[2] = 0; return out; }
    const C = this.value(x);
    const g = -C * (1 / r + 1 / this.lambda) / r;
    out[0] = g * dx; out[1] = g * dy; out[2] = g * dz;
    return out;
  }
}

/** ECM density on a regular grid, rebuilt from fibre segments (mass splatted trilinearly). */
export class ECMGrid {
  constructor({ size, h = 10, origin = null, massPerUnitDensity = 3 }) {
    this.h = h;
    this.nx = Math.max(2, Math.round(size[0] / h));
    this.ny = Math.max(2, Math.round(size[1] / h));
    this.nz = Math.max(2, Math.round(size[2] / h));
    this.origin = origin ?? [-size[0] / 2, -size[1] / 2, -size[2] / 2];
    const n = this.nx * this.ny * this.nz;
    this.E = new Float64Array(n);
    this.gx = new Float64Array(n); this.gy = new Float64Array(n); this.gz = new Float64Array(n);
    this.massPerUnitDensity = massPerUnitDensity; // fibre mass in one voxel that counts as E = 1
  }

  clear() { this.E.fill(0); }

  _coords(p) {
    return [
      Math.min(Math.max((p[0] - this.origin[0]) / this.h - 0.5, 0), this.nx - 1.0001),
      Math.min(Math.max((p[1] - this.origin[1]) / this.h - 0.5, 0), this.ny - 1.0001),
      Math.min(Math.max((p[2] - this.origin[2]) / this.h - 0.5, 0), this.nz - 1.0001),
    ];
  }

  /** Deposit `mass` at a point (trilinear splat). */
  splat(p, mass) {
    const [fx, fy, fz] = this._coords(p);
    const x0 = Math.floor(fx), y0 = Math.floor(fy), z0 = Math.floor(fz);
    const tx = fx - x0, ty = fy - y0, tz = fz - z0;
    const a = mass / this.massPerUnitDensity;
    const { nx, ny, E } = this;
    for (let c = 0; c < 8; c++) {
      const dx = c & 1, dy = (c >> 1) & 1, dz = (c >> 2) & 1;
      const w = (dx ? tx : 1 - tx) * (dy ? ty : 1 - ty) * (dz ? tz : 1 - tz);
      E[((z0 + dz) * ny + (y0 + dy)) * nx + x0 + dx] += a * w;
    }
  }

  /** Rebuild E from fibres ({ p, dir, len, mass }) and recompute its gradient. */
  rebuild(fibres) {
    this.clear();
    for (const f of fibres) {
      if (f.mass <= 0) continue;
      const m = f.mass / 3;
      for (const s of [-0.33, 0, 0.33]) {
        this.splat([f.p[0] + s * f.len * f.dir[0], f.p[1] + s * f.len * f.dir[1], f.p[2] + s * f.len * f.dir[2]], m);
      }
    }
    this.computeGradient();
  }

  computeGradient() {
    const { nx, ny, nz, E, gx, gy, gz, h } = this;
    for (let z = 0; z < nz; z++) for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
      const i = (z * ny + y) * nx + x;
      const xm = x > 0 ? i - 1 : i, xp = x < nx - 1 ? i + 1 : i;
      const ym = y > 0 ? i - nx : i, yp = y < ny - 1 ? i + nx : i;
      const zm = z > 0 ? i - nx * ny : i, zp = z < nz - 1 ? i + nx * ny : i;
      gx[i] = (E[xp] - E[xm]) / (h * (x > 0 && x < nx - 1 ? 2 : 1));
      gy[i] = (E[yp] - E[ym]) / (h * (y > 0 && y < ny - 1 ? 2 : 1));
      gz[i] = (E[zp] - E[zm]) / (h * (z > 0 && z < nz - 1 ? 2 : 1));
    }
  }

  _interp(arr, p) {
    const [fx, fy, fz] = this._coords(p);
    const x0 = Math.floor(fx), y0 = Math.floor(fy), z0 = Math.floor(fz);
    const tx = fx - x0, ty = fy - y0, tz = fz - z0;
    const { nx, ny } = this;
    let s = 0;
    for (let c = 0; c < 8; c++) {
      const dx = c & 1, dy = (c >> 1) & 1, dz = (c >> 2) & 1;
      const w = (dx ? tx : 1 - tx) * (dy ? ty : 1 - ty) * (dz ? tz : 1 - tz);
      s += w * arr[((z0 + dz) * ny + (y0 + dy)) * nx + x0 + dx];
    }
    return s;
  }

  density(p) { return this._interp(this.E, p); }

  gradient(p, out = [0, 0, 0]) {
    out[0] = this._interp(this.gx, p); out[1] = this._interp(this.gy, p); out[2] = this._interp(this.gz, p);
    return out;
  }

  mean() { let s = 0; for (const v of this.E) s += v; return s / this.E.length; }
}

/**
 * Fibroblasts depositing collagen + migrating cells following chemotaxis and
 * haptotaxis inside a box (time in hours, lengths in micrometres).
 */
export class MigrationSim {
  static defaults = {
    chi: 6e5,        // chemotactic sensitivity (um^2/h per unit C)
    alpha: 60,       // receptor saturation (Anderson & Chaplain's chi(c) = chi/(1 + alpha c))
    rho: 600,        // haptotactic coefficient (um^2/h per unit E)
    Dn: 30,          // random motility of migrating cells (um^2/h)
    vMax: 60,        // maximal crawling speed (um/h)
    mmp: 0.6,        // matrix degradation by migrating cells (1/h within rMMP)
    rMMP: 14,
    fibroSpeed: 18,  // fibroblast RMS speed (um/h)
    fibroPersistence: 1.5, // h
    fibreLength: 22, // um
    quantum: 0.6,    // fibril mass per deposited fibre segment
    maxFibres: 2600,
    deposit: true,
  };

  constructor(o = {}) {
    this.size = o.size ?? [420, 110, 300];
    this.rng = o.rng ?? new RNG(o.seed ?? 11);
    this.params = { ...MigrationSim.defaults, ...(o.params ?? {}) };
    this.collagenParams = { ...collagen.defaults, ...(o.collagen ?? {}) };
    const [Lx, Ly, Lz] = this.size;
    this.source = new PointSourceField({ p: o.sourcePos ?? [Lx * 0.36, 0, 0], radius: o.sourceRadius ?? 14, D: o.D ?? 100, k: o.k ?? 0.004 });
    this.ecm = new ECMGrid({ size: this.size, h: o.h ?? 10 });
    this.fibres = [];
    this.fibro = [];
    this.migr = [];
    this.t = 0;
    this.deposited = 0;   // total fibril mass deposited
    this.degraded = 0;    // total fibril mass degraded
    this.cleaved = 0;     // propeptide pairs released
    this.work = rk4Work(4);
    this._rhs = collagen.rhs(this.collagenParams);
    const r = this.rng;
    for (let i = 0; i < (o.nFibro ?? 36); i++) {
      const p = [r.uniform(-0.3, 0.2) * Lx, r.uniform(-0.35, 0.35) * Ly, r.uniform(-0.42, 0.42) * Lz];
      const d = r.onSphere();
      const s = this.params.fibroSpeed;
      this.fibro.push({ id: i, p, v: [d[0] * s, d[1] * s * 0.4, d[2] * s], y: Float64Array.from([r.uniform(0, 0.5), 0, r.uniform(0, 0.4), 0]), fibres: 0 });
    }
    for (let i = 0; i < (o.nMigr ?? 48); i++) {
      const p = [-Lx * 0.42 + r.uniform(0, 0.08) * Lx, r.uniform(-0.3, 0.3) * Ly, r.uniform(-0.3, 0.3) * Lz];
      this.migr.push({ id: i, p, v: [0, 0, 0], start: p.slice(), arrived: false, trail: [p.slice()] });
    }
    // pre-existing interstitial matrix: sparse random fibres
    for (let i = 0; i < (o.initialFibres ?? 120); i++) {
      const p = [r.uniform(-0.48, 0.48) * Lx, r.uniform(-0.45, 0.45) * Ly, r.uniform(-0.48, 0.48) * Lz];
      const d = r.onSphere();
      this.fibres.push({ p, dir: d, len: this.params.fibreLength * r.uniform(0.6, 1.2), mass: this.params.quantum * r.uniform(0.4, 0.9), born: 0, by: -1 });
    }
    this.ecm.rebuild(this.fibres);
  }

  /** Mean distance of migrating cells to the source surface. */
  meanDistance() {
    let s = 0;
    for (const c of this.migr) s += Math.max(0, Math.hypot(c.p[0] - this.source.p[0], c.p[1] - this.source.p[1], c.p[2] - this.source.p[2]) - this.source.radius);
    return s / (this.migr.length || 1);
  }

  fibreMass() { let s = 0; for (const f of this.fibres) s += f.mass; return s; }

  _confine(p, v) {
    for (let k = 0; k < 3; k++) {
      const half = this.size[k] / 2 - 4;
      if (p[k] > half) { p[k] = half; if (v) v[k] = -Math.abs(v[k]); }
      if (p[k] < -half) { p[k] = -half; if (v) v[k] = Math.abs(v[k]); }
    }
  }

  /** Velocity of a migrating cell from chemotaxis + haptotaxis (no noise). */
  driftVelocity(p, out = [0, 0, 0]) {
    const P = this.params;
    const C = this.source.value(p);
    const gC = this.source.gradient(p);
    const gE = this.ecm.gradient(p);
    const s = P.chi / (1 + P.alpha * C);
    for (let k = 0; k < 3; k++) out[k] = s * gC[k] + P.rho * gE[k];
    const sp = Math.hypot(out[0], out[1], out[2]);
    if (sp > P.vMax) for (let k = 0; k < 3; k++) out[k] *= P.vMax / sp;
    return out;
  }

  step(dt) {
    const P = this.params, r = this.rng;
    let changed = false;
    // ---- fibroblasts: persistent random walk + collagen synthesis / deposition
    const sv = P.fibroSpeed * Math.sqrt(2 / (3 * P.fibroPersistence));
    for (const f of this.fibro) {
      for (let k = 0; k < 3; k++) {
        const flat = k === 1 ? 0.35 : 1; // cells in a thin slab mostly crawl in-plane
        f.v[k] += (-f.v[k] / P.fibroPersistence) * dt + sv * flat * Math.sqrt(dt) * r.normal();
        f.p[k] += f.v[k] * dt;
      }
      this._confine(f.p, f.v);
      const before = f.y[3];
      rk4Step(this._rhs, this.t, f.y, dt, this.work);
      this.cleaved += f.y[3] - before;
      if (P.deposit && f.y[2] >= P.quantum) {
        const sp = Math.hypot(f.v[0], f.v[1], f.v[2]) || 1;
        const dir = [f.v[0] / sp, f.v[1] / sp, f.v[2] / sp];
        f.y[2] -= P.quantum;
        this.deposited += P.quantum;
        const fib = { p: f.p.slice(), dir, len: P.fibreLength * (0.8 + 0.4 * r.next()), mass: P.quantum, born: this.t, by: f.id };
        if (this.fibres.length < P.maxFibres) this.fibres.push(fib);
        else {
          // recycle the lightest (most degraded) fibre
          let j = 0;
          for (let i = 1; i < this.fibres.length; i++) if (this.fibres[i].mass < this.fibres[j].mass) j = i;
          this.fibres[j] = fib;
        }
        f.fibres++;
        changed = true;
      }
    }
    // ---- migrating cells: chemotaxis + haptotaxis + random motility, MMP degradation
    const noise = Math.sqrt(2 * P.Dn * dt);
    const v = [0, 0, 0];
    for (const c of this.migr) {
      const d = Math.hypot(c.p[0] - this.source.p[0], c.p[1] - this.source.p[1], c.p[2] - this.source.p[2]);
      if (d < this.source.radius + 12) { c.arrived = true; c.v[0] = c.v[1] = c.v[2] = 0; continue; }
      this.driftVelocity(c.p, v);
      for (let k = 0; k < 3; k++) {
        const dx = v[k] * dt + noise * r.normal() * (k === 1 ? 0.5 : 1);
        c.p[k] += dx;
        c.v[k] = 0.8 * c.v[k] + 0.2 * (dx / Math.max(dt, 1e-9));
      }
      this._confine(c.p, null);
      if (P.mmp > 0) {
        const r2 = P.rMMP * P.rMMP;
        for (const f of this.fibres) {
          const ex = f.p[0] - c.p[0], ey = f.p[1] - c.p[1], ez = f.p[2] - c.p[2];
          if (ex * ex + ey * ey + ez * ez > r2 + f.len * f.len * 0.25) continue;
          const loss = f.mass * (1 - Math.exp(-P.mmp * dt));
          f.mass -= loss;
          this.degraded += loss;
          changed = true;
        }
      }
    }
    if (changed) {
      this.fibres = this.fibres.filter((f) => f.mass > 0.02 * P.quantum);
      this.ecm.rebuild(this.fibres);
    }
    this.t += dt;
    return changed;
  }

  /** Record positions for trails (call at the display rate). */
  recordTrails(maxLen = 60, minStep = 2) {
    for (const c of this.migr) {
      const last = c.trail[c.trail.length - 1];
      if (Math.hypot(c.p[0] - last[0], c.p[1] - last[1], c.p[2] - last[2]) > minStep) {
        c.trail.push(c.p.slice());
        if (c.trail.length > maxLen) c.trail.shift();
      }
    }
  }
}

export const REACTIONS = [
  { id: 'ecm.secrete', order: 4, pathway: 'ECM', reactants: [], products: ['Procollagen'], modifiers: ['Fibroblast'], label: 'Procollagen secretion' },
  { id: 'ecm.cleave', order: 4, pathway: 'ECM', reactants: ['Procollagen'], products: ['Tropocollagen', 'Propeptides'], modifiers: ['ADAMTS2', 'BMP-1'], label: 'N- and C-propeptide cleavage' },
  { id: 'ecm.nucleate', order: 4, pathway: 'ECM', reactants: ['Tropocollagen', 'Tropocollagen'], products: ['Collagen fibril'], label: 'Fibril nucleation' },
  { id: 'ecm.elongate', order: 4, pathway: 'ECM', reactants: ['Tropocollagen', 'Collagen fibril'], products: ['Collagen fibril'], label: 'Fibril elongation' },
  { id: 'ecm.lox', order: 4, pathway: 'ECM', reactants: ['Collagen fibril'], products: ['Crosslinked collagen fibril'], modifiers: ['Lysyl oxidase'], label: 'Covalent crosslinking' },
  { id: 'ecm.mmp', order: 4, pathway: 'ECM', reactants: ['Collagen fibril'], products: ['Collagen fragments'], modifiers: ['MMP-1'], label: 'Collagenolysis by migrating cells' },
  { id: 'ecm.integrin', order: 4, pathway: 'ECM', reactants: ['Integrin a2b1', 'Collagen fibril'], products: ['Integrin:collagen'], label: 'Integrin adhesion to fibrils (haptotaxis)' },
  { id: 'ecm.chemokine', order: 4, pathway: 'ECM', reactants: ['CXCL12', 'CXCR4'], products: ['CXCL12:CXCR4'], label: 'Chemokine receptor binding (chemotaxis)' },
  { id: 'ecm.move', order: 4, pathway: 'ECM', reactants: ['CXCL12:CXCR4', 'Integrin:collagen'], products: ['Cell movement'], label: 'Cell movement = chemotaxis + haptotaxis' },
];
