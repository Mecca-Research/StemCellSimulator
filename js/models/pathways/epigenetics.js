// Epigenetic modifications and the Waddington landscape.
//
// Paper ("Epigenetic Modifications"):
//   DNA + DNMT -> methylated DNA
//   Histone + HAT -> acetylated histone
//   acetylated histone + HDAC -> histone
// as mass-action fractions:
//   dm/dt = kM * DNMT * (1 - m) - kDM * m          (TET / replication dilution)
//   da/dt = kA * HAT  * (1 - a) - kD * HDAC * a
//
// Fate circuit: two cross-antagonistic, self-activating master regulators
// (Huang, Guo, May & Enver 2007, Dev Biol 305:695):
//   dx/dt = s * x^n/(th^n + x^n) + b * th^n/(th^n + y^n) - k x
//   dy/dt = s * y^n/(th^n + y^n) + b * th^n/(th^n + x^n) - k y
// Strong self-activation s gives a central "multipotent" attractor next to
// the two committed ones; lowering s during development removes it.
//
// Landscape: U(x, y) = -ln P_ss(x, y) estimated from Langevin trajectories
// (Wang, Xu & Wang 2008, PNAS 105:12271).
//
// The coupling of chromatin marks to the circuit is PHENOMENOLOGICAL:
// acetylation raises transcriptional noise (plasticity), DNA methylation
// strengthens repression of the alternative fate (canalisation).
import { RNG } from '../core/rng.js';

export const marks = {
  defaults: { kM: 1.0, kDM: 0.5, kA: 1.0, kD: 1.0, DNMT: 0.5, HAT: 0.5, HDAC: 0.5 },
  rhs(p) {
    return (t, y, dy) => {
      dy[0] = p.kM * p.DNMT * (1 - y[0]) - p.kDM * y[0];
      dy[1] = p.kA * p.HAT * (1 - y[1]) - p.kD * p.HDAC * y[1];
    };
  },
  steadyState(p) {
    return [(p.kM * p.DNMT) / (p.kM * p.DNMT + p.kDM), (p.kA * p.HAT) / (p.kA * p.HAT + p.kD * p.HDAC)];
  },
};

export const circuit = {
  defaults: { sHigh: 1.0, sLow: 0.3, b: 1.0, k: 1.0, th: 0.5, n: 4 },

  /** Circuit parameters at developmental time tau in [0, 1] and chromatin state (m, a). */
  params(tau, m = 0.5, a = 0.5, base = circuit.defaults) {
    return {
      ...base,
      s: base.sHigh + (base.sLow - base.sHigh) * tau,
      b: base.b * (0.7 + 0.6 * m),
      sigma: 0.04 + 0.22 * a,
    };
  },

  drift(p, x, y, out) {
    const hn = p.th ** p.n;
    const xn = Math.max(x, 0) ** p.n, yn = Math.max(y, 0) ** p.n;
    out[0] = (p.s * xn) / (hn + xn) + (p.b * hn) / (hn + yn) - p.k * x;
    out[1] = (p.s * yn) / (hn + yn) + (p.b * hn) / (hn + xn) - p.k * y;
    return out;
  },

  /** Euler-Maruyama step for a walker [x, y]. */
  step(p, w, dt, rng, out = [0, 0]) {
    circuit.drift(p, w[0], w[1], out);
    const sd = p.sigma * Math.sqrt(dt);
    w[0] = Math.max(0, w[0] + out[0] * dt + sd * rng.normal());
    w[1] = Math.max(0, w[1] + out[1] * dt + sd * rng.normal());
    return w;
  },

  /** Deterministic attractors: integrate from a grid of starts and cluster. */
  attractors(p, { grid = 7, max = 2.2, T = 60, dt = 0.02 } = {}) {
    const found = [];
    const d = [0, 0];
    for (let i = 0; i < grid; i++) for (let j = 0; j < grid; j++) {
      // small asymmetric offset: starts exactly on the symmetric diagonal would
      // otherwise stay on it and report the saddle as an attractor
      let x = ((i + 0.5) / grid) * max + 0.0131, y = ((j + 0.5) / grid) * max - 0.0077;
      for (let s = 0; s < T / dt; s++) { circuit.drift(p, x, y, d); x += d[0] * dt; y += d[1] * dt; }
      if (!found.some((f) => Math.hypot(f[0] - x, f[1] - y) < 0.05)) found.push([x, y]);
    }
    return found.sort((u, v) => u[0] - v[0]);
  },
};

/**
 * Quasi-potential U = -ln P on an N x N grid over [0, max]^2 from Langevin
 * walkers started uniformly. Returns { U: Float32Array (row-major, y rows), N, max }.
 */
export function landscape(p, { N = 64, max = 2.2, walkers = 600, steps = 900, burn = 150, dt = 0.02, seed = 3 } = {}) {
  const rng = new RNG(seed);
  const H = new Float64Array(N * N);
  const out = [0, 0];
  for (let w = 0; w < walkers; w++) {
    const s = [rng.uniform(0, max), rng.uniform(0, max)];
    for (let k = 0; k < steps; k++) {
      circuit.step(p, s, dt, rng, out);
      if (k < burn) continue;
      const i = Math.min(N - 1, Math.floor((s[0] / max) * N)), j = Math.min(N - 1, Math.floor((s[1] / max) * N));
      H[j * N + i] += 1;
    }
  }
  // light Gaussian smoothing of the histogram, then U = -ln(P + floor)
  const S = new Float64Array(N * N);
  const kern = [0.25, 0.5, 0.25];
  const tmp = new Float64Array(N * N);
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    let v = 0; for (let d = -1; d <= 1; d++) v += kern[d + 1] * H[j * N + Math.min(N - 1, Math.max(0, i + d))]; tmp[j * N + i] = v;
  }
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    let v = 0; for (let d = -1; d <= 1; d++) v += kern[d + 1] * tmp[Math.min(N - 1, Math.max(0, j + d)) * N + i]; S[j * N + i] = v;
  }
  let total = 0;
  for (const v of S) total += v;
  const U = new Float32Array(N * N);
  const floor = 0.2 / (N * N);
  let lo = Infinity;
  for (let k = 0; k < N * N; k++) { U[k] = -Math.log(S[k] / total + floor); if (U[k] < lo) lo = U[k]; }
  for (let k = 0; k < N * N; k++) U[k] -= lo;
  return { U, N, max };
}

/** Bilinear lookup of U at (x, y). */
export function sampleU(L, x, y) {
  const { U, N, max } = L;
  const fx = Math.min(Math.max((x / max) * N - 0.5, 0), N - 1.001), fy = Math.min(Math.max((y / max) * N - 0.5, 0), N - 1.001);
  const i = Math.floor(fx), j = Math.floor(fy), u = fx - i, v = fy - j;
  const a = U[j * N + i], b = U[j * N + i + 1], c = U[(j + 1) * N + i], d = U[(j + 1) * N + i + 1];
  return (a * (1 - u) + b * u) * (1 - v) + (c * (1 - u) + d * u) * v;
}

export const REACTIONS = [
  { id: 'epi.dnmt', order: 1, pathway: 'Epigenetics', reactants: ['DNA (CpG)'], products: ['Methylated DNA'], modifiers: ['DNMT', 'SAM'], label: 'DNA methylation' },
  { id: 'epi.tet', order: 1, pathway: 'Epigenetics', reactants: ['Methylated DNA'], products: ['DNA (CpG)'], modifiers: ['TET'], label: 'Demethylation / dilution' },
  { id: 'epi.hat', order: 1, pathway: 'Epigenetics', reactants: ['Histone', 'Acetyl-CoA'], products: ['Acetylated histone'], modifiers: ['HAT'], label: 'Histone acetylation' },
  { id: 'epi.hdac', order: 1, pathway: 'Epigenetics', reactants: ['Acetylated histone'], products: ['Histone'], modifiers: ['HDAC'], label: 'Histone deacetylation' },
  { id: 'fate.selfA', order: 2, pathway: 'Fate circuit', reactants: [], products: ['GATA1'], modifiers: ['GATA1'], label: 'Self-activation (fate A)' },
  { id: 'fate.selfB', order: 2, pathway: 'Fate circuit', reactants: [], products: ['PU.1'], modifiers: ['PU.1'], label: 'Self-activation (fate B)' },
  { id: 'fate.AB', order: 2, pathway: 'Fate circuit', reactants: [], products: ['GATA1'], modifiers: ['PU.1'], label: 'Cross-repression B ⊣ A' },
  { id: 'fate.BA', order: 2, pathway: 'Fate circuit', reactants: [], products: ['PU.1'], modifiers: ['GATA1'], label: 'Cross-repression A ⊣ B' },
];
