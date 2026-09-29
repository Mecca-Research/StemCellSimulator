// Center-based 3D cell mechanics (overdamped, off-lattice), in the spirit of
// Drasdo & Hoehme 2005 and Meineke/Osborne 2017 (Chaste):
//
//   eta * dx_i/dt = sum_j F_ij + F_ext,i + sqrt(2 * eta * kT) * xi(t)
//
// Pair force along the centre line, s = R_i + R_j (rest distance):
//   d <  s            F = mu_rep * s * ln(1 + (d - s) / s)       (repulsion, F < 0)
//   s <= d < s * cut  F = mu_adh_ij * (d - s) * exp(-5 (d - s) / s)  (adhesion)
//
// Cells grow in volume through a G1-S-G2-M cycle and divide along an axis;
// mitosis is resolved into prophase .. cytokinesis so it can be rendered.
import { SpatialHash } from '../core/spatial.js';
import { RNG } from '../core/rng.js';

export const PHASES = ['G1', 'S', 'G2', 'M'];
// Typical mammalian proportions of the cycle (Alberts, Molecular Biology of the Cell)
export const PHASE_FRACTIONS = { G1: 0.45, S: 0.3, G2: 0.19, M: 0.06 };
export const MITOSIS_STAGES = [
  ['prophase', 0.0], ['prometaphase', 0.18], ['metaphase', 0.32], ['anaphase', 0.55],
  ['telophase', 0.72], ['cytokinesis', 0.84],
];

export function mitosisStage(progress) {
  let name = MITOSIS_STAGES[0][0];
  for (const [n, start] of MITOSIS_STAGES) if (progress >= start) name = n;
  return name;
}

let nextId = 1;

export class Cell {
  constructor(o = {}) {
    this.id = o.id ?? nextId++;
    this.p = o.p ? o.p.slice() : [0, 0, 0];
    this.v = [0, 0, 0];
    this.f = [0, 0, 0];
    this.R0 = o.R0 ?? 5;               // radius right after division
    this.R = o.R ?? this.R0 * 1.12;     // current radius
    this.type = o.type ?? 0;
    this.phase = o.phase ?? 'G1';
    this.phaseT = o.phaseT ?? 0;        // time spent in current phase
    this.cycle = o.cycle ?? 16;          // total cycle length
    this.axis = o.axis ? o.axis.slice() : [1, 0, 0];
    this.alive = true;
    this.quiescent = false;
    this.generation = o.generation ?? 0;
    this.parent = o.parent ?? null;
    this.state = o.state ? structuredClone(o.state) : {};
    this.age = 0;
  }

  get phaseLength() { return this.cycle * PHASE_FRACTIONS[this.phase]; }

  /** 0..1 progress through mitosis (only meaningful in M). */
  get mitosis() { return this.phase === 'M' ? Math.min(this.phaseT / this.phaseLength, 1) : 0; }
}

export class Tissue {
  /**
   * @param {object} o
   * @param {number} [o.eta] drag
   * @param {number} [o.muRep] repulsion stiffness
   * @param {(a:number,b:number)=>number} [o.adhesion] adhesion strength by types
   * @param {number} [o.cut] interaction cutoff multiple of rest distance
   * @param {number} [o.noise] thermal/active motility amplitude
   * @param {{y:number, k:number}|null} [o.substrate] harmonic confinement of height
   * @param {number|null} [o.contactInhibition] neighbours above which cells arrest in G1
   */
  constructor(o = {}) {
    this.eta = o.eta ?? 1;
    this.muRep = o.muRep ?? 12;
    this.adhesion = o.adhesion ?? (() => 1.2);
    this.cut = o.cut ?? 1.25;
    this.noise = o.noise ?? 0.3;
    this.substrate = o.substrate ?? null;
    this.contactInhibition = o.contactInhibition ?? null;
    this.bounds = o.bounds ?? null; // {r} spherical container
    this.external = o.external ?? null; // (cell, out[3], t) => void
    this.onDivide = o.onDivide ?? null;  // (mother, a, b) => void
    this.growth = o.growth ?? true;
    // 3 = volumetric growth (R ~ V^1/3); 2 = monolayer of fixed height (lateral R ~ V^1/2)
    this.dims = o.dims ?? 3;
    this.rng = o.rng ?? new RNG(o.seed ?? 7);
    this.cells = [];
    this.hash = new SpatialHash(o.hashSize ?? 14);
    this.t = 0;
    this.neighbours = new Map(); // id -> [ids] within contact distance
    this.divisions = 0;
  }

  add(cell) { this.cells.push(cell); return cell; }

  rebuildHash() {
    this.hash.clear();
    let maxR = 0;
    for (let i = 0; i < this.cells.length; i++) {
      const c = this.cells[i];
      this.hash.insert(i, c.p[0], c.p[1], c.p[2]);
      if (c.R > maxR) maxR = c.R;
    }
    // make sure the bucket size covers the interaction range
    const need = 2 * maxR * this.cut;
    if (need > this.hash.h * 1.001) {
      this.hash = new SpatialHash(need);
      for (let i = 0; i < this.cells.length; i++) this.hash.insert(i, ...this.cells[i].p);
    }
  }

  computeForces() {
    const cells = this.cells;
    this.neighbours.clear();
    for (const c of cells) { c.f[0] = c.f[1] = c.f[2] = 0; }
    for (let i = 0; i < cells.length; i++) {
      const a = cells[i];
      const nb = [];
      this.hash.near(a.p[0], a.p[1], a.p[2], (j) => {
        if (j === i) return;
        const b = cells[j];
        const dx = b.p[0] - a.p[0], dy = b.p[1] - a.p[1], dz = b.p[2] - a.p[2];
        const d = Math.hypot(dx, dy, dz) || 1e-6;
        const s = restDistance(a, b);
        if (d >= s * this.cut) return;
        if (d < s * 1.08) nb.push(b.id);
        let F;
        if (d < s) F = this.muRep * s * Math.log(1 + (d - s) / s);
        else F = this.adhesion(a.type, b.type, a, b) * (d - s) * Math.exp((-5 * (d - s)) / s);
        a.f[0] += (F * dx) / d; a.f[1] += (F * dy) / d; a.f[2] += (F * dz) / d;
      });
      this.neighbours.set(a.id, nb);
    }
  }

  step(dt) {
    this.rebuildHash();
    this.computeForces();
    const ext = [0, 0, 0];
    const sd = Math.sqrt(2 * this.noise * dt) / this.eta;
    for (const c of this.cells) {
      if (this.external) { ext[0] = ext[1] = ext[2] = 0; this.external(c, ext, this.t); c.f[0] += ext[0]; c.f[1] += ext[1]; c.f[2] += ext[2]; }
      if (this.substrate) c.f[1] += -this.substrate.k * (c.p[1] - (this.substrate.y + (this.substrate.lift ? c.R : 0)));
      if (this.bounds) {
        const r = Math.hypot(c.p[0], c.p[1], c.p[2]);
        const over = r + c.R - this.bounds.r;
        if (over > 0) for (let k = 0; k < 3; k++) c.f[k] -= (this.bounds.k ?? 8) * over * c.p[k] / (r || 1);
      }
      for (let k = 0; k < 3; k++) {
        c.v[k] = c.f[k] / this.eta;
        c.p[k] += c.v[k] * dt + sd * this.rng.normal();
      }
    }
    if (this.growth) this.advanceCycle(dt);
    this.t += dt;
  }

  advanceCycle(dt) {
    const born = [];
    for (const c of this.cells) {
      c.age += dt;
      if (c.quiescent) continue;
      if (this.contactInhibition && c.phase === 'G1' && c.phaseT > c.phaseLength * 0.9) {
        const nb = this.neighbours.get(c.id)?.length ?? 0;
        if (nb >= this.contactInhibition) continue; // arrest at the restriction point
      }
      c.phaseT += dt;
      // volume grows linearly from V0 to 2 V0 over interphase
      if (c.phase !== 'M') {
        const interphase = c.cycle * (1 - PHASE_FRACTIONS.M);
        const elapsed = timeInInterphase(c);
        c.R = c.R0 * Math.pow(1 + Math.min(elapsed / interphase, 1), 1 / this.dims);
      }
      if (c.phaseT >= c.phaseLength) {
        const idx = PHASES.indexOf(c.phase);
        if (c.phase === 'M') {
          born.push(...this.divide(c));
        } else {
          c.phase = PHASES[idx + 1];
          c.phaseT = 0;
          if (c.phase === 'M') chooseAxis(c, this);
        }
      }
    }
    if (born.length) this.cells.push(...born);
    this.cells = this.cells.filter((c) => c.alive);
  }

  /** Split a mother into two daughters along its division axis. */
  divide(m) {
    const off = m.R * 0.5;
    const R0 = m.R / Math.pow(2, 1 / this.dims);
    const jitter = () => (this.rng.next() - 0.5) * 0.3;
    const cycle = () => m.cycle * (1 + (this.rng.next() - 0.5) * 0.2);
    const mk = (sgn) => new Cell({
      p: [m.p[0] + sgn * off * m.axis[0], m.p[1] + sgn * off * m.axis[1], m.p[2] + sgn * off * m.axis[2]],
      R0, R: R0, type: m.type, phase: 'G1', phaseT: jitter(), cycle: cycle(),
      generation: m.generation + 1, parent: m.id, state: m.state,
    });
    const a = mk(1), b = mk(-1);
    m.alive = false;
    this.divisions++;
    this.onDivide?.(m, a, b);
    return [a, b];
  }

  count(pred = () => true) { return this.cells.reduce((n, c) => n + (pred(c) ? 1 : 0), 0); }

  centroid() {
    const s = [0, 0, 0];
    for (const c of this.cells) for (let k = 0; k < 3; k++) s[k] += c.p[k];
    const n = this.cells.length || 1;
    return s.map((v) => v / n);
  }
}

function restDistance(a, b) { return a.R + b.R; }

function timeInInterphase(c) {
  let t = c.phaseT;
  for (const ph of PHASES) { if (ph === c.phase) break; t += c.cycle * PHASE_FRACTIONS[ph]; }
  return t;
}

/** Division axis: random in 3D, or in the substrate plane for monolayers (planar division). */
function chooseAxis(c, tissue) {
  const a = tissue.rng.onSphere();
  if (tissue.substrate) a[1] = 0;
  const n = Math.hypot(a[0], a[1], a[2]) || 1;
  c.axis = [a[0] / n, a[1] / n, a[2] / n];
}

/** Exponential-growth doubling time estimate from (t, N) samples (least squares on ln N). */
export function doublingTime(ts, ns) {
  const pts = ts.map((t, i) => [t, Math.log(ns[i])]).filter((p) => Number.isFinite(p[1]));
  const n = pts.length;
  if (n < 2) return NaN;
  const mx = pts.reduce((s, p) => s + p[0], 0) / n, my = pts.reduce((s, p) => s + p[1], 0) / n;
  let sxy = 0, sxx = 0;
  for (const [x, y] of pts) { sxy += (x - mx) * (y - my); sxx += (x - mx) ** 2; }
  const slope = sxy / sxx;
  return Math.LN2 / slope;
}
