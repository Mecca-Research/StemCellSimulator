// Morphogen gradients (paper: "Higher-Level Functional Interactions -
// morphogen gradients: gradients of signalling molecules guide tissue
// patterning").
//
// Source-diffusion-degradation (SDD) model
//     dC/dt = D lap C - k C + s(x)
// Steady state in 1D from a source at x = 0:
//     C(x) = C0 exp(-x / lambda),   lambda = sqrt(D / k)
// which, read through concentration thresholds, gives Wolpert's "French flag"
// (Wolpert 1969, J Theor Biol 25:1): boundary i sits at x_i = lambda ln(C0 / T_i).
//
// Measured numbers used as defaults: Dpp in the Drosophila wing disc,
// D = 0.10 um^2/s, k = 2.52e-4 1/s, lambda ~ 20 um (Kicheva et al. 2007,
// Science 315:521, FRAP of Dpp-GFP). Local accumulation time of the SDD model:
// tau(x) = (1 + x / lambda) / (2 k) (Berezhkovskii, Sample & Shvartsman 2010,
// Biophys J 99:L59). Positional error of a threshold read-out of an
// exponential gradient with multiplicative noise CV: sigma_x = lambda * CV
// (e.g. Gregor et al. 2007, Cell 130:153; Bollenbach et al. 2008, Development
// 135:1137).
//
// Pure JS; no DOM / three.js.

/** Kicheva et al. 2007 (Dpp, wing disc): D in um^2/s, k in 1/s. */
export const KICHEVA_DPP = { D: 0.10, k: 2.52e-4 };

export const decayLength = (D, k) => Math.sqrt(D / k);

/** Semi-infinite steady state C0 exp(-x / lambda). */
export function steady1D(x, { C0 = 1, lambda }) {
  return C0 * Math.exp(-Math.max(x, 0) / lambda);
}

/** Finite domain [0, L] with C(0) = C0 and no flux at x = L. */
export function steady1DFinite(x, { C0 = 1, lambda, L }) {
  return (C0 * Math.cosh((L - x) / lambda)) / Math.cosh(L / lambda);
}

/** Amplitude at the source for a flux boundary -D C'(0) = j: C0 = j lambda / D. */
export const fluxToAmplitude = (j, D, k) => (j * decayLength(D, k)) / D;

/**
 * Point source of strength Q (molecules / time) in an infinite 3D medium:
 *   D lap C - k C + Q delta(r) = 0  ->  C(r) = Q exp(-r / lambda) / (4 pi D r)
 * (screened-Poisson / Yukawa Green's function). r is clamped at rMin (the
 * source's own radius) to avoid the singularity.
 */
export function pointSource3D(r, { Q = 1, D, k, rMin = 1e-6 }) {
  const rr = Math.max(r, rMin);
  return (Q * Math.exp(-rr / decayLength(D, k))) / (4 * Math.PI * D * rr);
}

/** dC/dr of the point source: -C (1/r + 1/lambda). */
export function pointSource3DSlope(r, o) {
  const rr = Math.max(r, o.rMin ?? 1e-6);
  return -pointSource3D(rr, o) * (1 / rr + 1 / decayLength(o.D, o.k));
}

/** Plane source emitting flux j into the half-space x > 0. */
export function planeSource3D(x, { j = 1, D, k }) {
  return fluxToAmplitude(j, D, k) * Math.exp(-Math.max(x, 0) / decayLength(D, k));
}

/** Where each threshold is crossed by C0 exp(-x/lambda): x_i = lambda ln(C0 / T_i). */
export function boundaryPositions(C0, lambda, thresholds) {
  return thresholds.map((T) => (T > 0 && T < C0 ? lambda * Math.log(C0 / T) : T >= C0 ? 0 : Infinity));
}

/**
 * French-flag read-out: thresholds sorted high -> low; returns 0 for C above
 * the first threshold (the fate nearest the source), thresholds.length for C
 * below all of them.
 */
export function frenchFlag(C, thresholds) {
  let i = 0;
  while (i < thresholds.length && C < thresholds[i]) i++;
  return i;
}

/**
 * Positional error of a threshold read-out: sigma_x = sigma_C / |dC/dx|.
 * For an exponential gradient and multiplicative noise (constant CV) this is
 * lambda * CV, independent of position. With counting (Poisson) noise of N
 * molecules in the sensing volume, CV = 1 / sqrt(N).
 */
export function positionalError({ lambda, cv = null, C = null, sensingVolume = null }) {
  const CV = cv ?? (C !== null && sensingVolume ? 1 / Math.sqrt(Math.max(C * sensingVolume, 1e-12)) : 0);
  return lambda * CV;
}

/** Local accumulation time of the SDD model (Berezhkovskii et al. 2010). */
export const accumulationTime = (x, { D, k }) => (1 + x / decayLength(D, k)) / (2 * k);

/**
 * Numerical steady state of D C'' - k C = 0 on [0, L] with C(0) = C0 and
 * C'(L) = 0, by second-order finite differences and the Thomas algorithm.
 * Returns { x, C }.
 */
export function solve1D({ n = 200, L, D, k, C0 = 1 }) {
  const h = L / n;
  // unknowns C_1..C_n (C_0 = C0 fixed); ghost C_{n+1} = C_{n-1} (no flux)
  const a = new Float64Array(n), b = new Float64Array(n), c = new Float64Array(n), d = new Float64Array(n);
  const r = D / (h * h);
  for (let i = 0; i < n; i++) {
    a[i] = -r; b[i] = 2 * r + k; c[i] = -r; d[i] = 0;
  }
  d[0] += r * C0; a[0] = 0;
  a[n - 1] = -2 * r; c[n - 1] = 0;
  // Thomas
  for (let i = 1; i < n; i++) {
    const m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    d[i] -= m * d[i - 1];
  }
  const C = new Float64Array(n + 1);
  C[0] = C0;
  C[n] = d[n - 1] / b[n - 1];
  for (let i = n - 2; i >= 0; i--) C[i + 1] = (d[i] - c[i] * C[i + 2]) / b[i];
  const x = Float64Array.from({ length: n + 1 }, (_, i) => i * h);
  return { x, C };
}

/**
 * Numerical 3D source-diffusion-degradation solver on a regular grid with an
 * optional tissue mask (diffusion only between tissue voxels = zero flux at
 * the tissue surface), per-voxel production s and Dirichlet (clamped) voxels
 * for source cells held at a fixed concentration.
 * Explicit Euler with CFL-safe substeps: dt <= h^2 / (6 D).
 */
export class GradientGrid3D {
  constructor({ nx, ny, nz, h = 1, D = 1, k = 0.01, origin = [0, 0, 0], mask = null }) {
    Object.assign(this, { nx, ny, nz, h, D, k, origin });
    const n = nx * ny * nz;
    this.n = n;
    this.C = new Float64Array(n);
    this.s = new Float64Array(n);           // production rate per voxel
    this.fixed = new Uint8Array(n);          // 1 = Dirichlet voxel
    this.fixedValue = new Float64Array(n);
    this.mask = mask ?? new Uint8Array(n).fill(1);
    this.work = new Float64Array(n);
    this.t = 0;
  }

  index(x, y, z) { return (z * this.ny + y) * this.nx + x; }

  /** World position of a voxel centre. */
  centre(i, out = [0, 0, 0]) {
    const x = i % this.nx, y = Math.floor(i / this.nx) % this.ny, z = Math.floor(i / (this.nx * this.ny));
    out[0] = this.origin[0] + (x + 0.5) * this.h;
    out[1] = this.origin[1] + (y + 0.5) * this.h;
    out[2] = this.origin[2] + (z + 0.5) * this.h;
    return out;
  }

  maxDt(safety = 0.9) {
    const dDiff = (this.h * this.h) / (6 * this.D);
    return safety * Math.min(dDiff, 1 / Math.max(this.k, 1e-12));
  }

  step(dt) {
    const { nx, ny, nz, C, s, mask, fixed, fixedValue, work, D, k } = this;
    const r = D / (this.h * this.h);
    const sxy = nx * ny;
    for (let z = 0; z < nz; z++) {
      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          const i = z * sxy + y * nx + x;
          if (!mask[i]) { work[i] = 0; continue; }
          if (fixed[i]) { work[i] = fixedValue[i]; continue; }
          const c = C[i];
          let lap = 0;
          // zero flux towards non-tissue voxels and the grid boundary
          if (x > 0 && mask[i - 1]) lap += C[i - 1] - c;
          if (x < nx - 1 && mask[i + 1]) lap += C[i + 1] - c;
          if (y > 0 && mask[i - nx]) lap += C[i - nx] - c;
          if (y < ny - 1 && mask[i + nx]) lap += C[i + nx] - c;
          if (z > 0 && mask[i - sxy]) lap += C[i - sxy] - c;
          if (z < nz - 1 && mask[i + sxy]) lap += C[i + sxy] - c;
          work[i] = c + dt * (r * lap - k * c + s[i]);
        }
      }
    }
    C.set(work);
    this.t += dt;
  }

  /** Advance by T with CFL-safe substeps (at most maxSteps); returns substeps taken. */
  advance(T, maxSteps = 10000) {
    if (!(T > 0)) return 0;
    const dtMax = this.maxDt();
    const n = Math.min(maxSteps, Math.max(1, Math.ceil(T / dtMax)));
    const dt = Math.min(T / n, dtMax);
    for (let i = 0; i < n; i++) this.step(dt);
    return n;
  }

  /** Trilinear sample at a world position (outside the grid -> nearest edge). */
  sample(p) {
    const { nx, ny, nz, h, origin, C } = this;
    const fx = Math.min(Math.max((p[0] - origin[0]) / h - 0.5, 0), nx - 1.001);
    const fy = Math.min(Math.max((p[1] - origin[1]) / h - 0.5, 0), ny - 1.001);
    const fz = Math.min(Math.max((p[2] - origin[2]) / h - 0.5, 0), nz - 1.001);
    const x0 = Math.floor(fx), y0 = Math.floor(fy), z0 = Math.floor(fz);
    const tx = fx - x0, ty = fy - y0, tz = fz - z0;
    const x1 = Math.min(x0 + 1, nx - 1), y1 = Math.min(y0 + 1, ny - 1), z1 = Math.min(z0 + 1, nz - 1);
    const I = (x, y, z) => C[(z * ny + y) * nx + x];
    const c00 = I(x0, y0, z0) * (1 - tx) + I(x1, y0, z0) * tx;
    const c10 = I(x0, y1, z0) * (1 - tx) + I(x1, y1, z0) * tx;
    const c01 = I(x0, y0, z1) * (1 - tx) + I(x1, y0, z1) * tx;
    const c11 = I(x0, y1, z1) * (1 - tx) + I(x1, y1, z1) * tx;
    return (c00 * (1 - ty) + c10 * ty) * (1 - tz) + (c01 * (1 - ty) + c11 * ty) * tz;
  }
}

/** Neural-tube progenitor domains read out from opposing Shh (ventral) and BMP (dorsal) gradients. */
export const NEURAL_TUBE_DOMAINS = [
  // ventral -> dorsal; markers after Briscoe & Ericson 2001 (Curr Opin Neurobiol 11:43) and
  // Dessaud, McMahon & Briscoe 2008 (Development 135:2489)
  { id: 'FP', label: 'floor plate (FoxA2)', color: '#4a6cff' },
  { id: 'p3', label: 'p3 (Nkx2.2)', color: '#b04dff' },
  { id: 'pMN', label: 'pMN (Olig2)', color: '#ff4fa3' },
  { id: 'p2-p0', label: 'p2-p0 (Nkx6.1 / Dbx)', color: '#ffb454' },
  { id: 'dP', label: 'dorsal progenitors (Pax7)', color: '#5ee0c1' },
  { id: 'RP', label: 'roof plate (Lmx1a)', color: '#e8f0ff' },
];

/**
 * Two-gradient read-out: Shh thresholds (high -> low) carve the ventral
 * domains, a BMP threshold the dorsal ones. Returns an index into
 * NEURAL_TUBE_DOMAINS. Phenomenological threshold logic (Wolpert's French
 * flag applied to two opposing sources).
 */
export function neuralTubeFate(shh, bmp, { shhT = [0.6, 0.3, 0.12], bmpT = [0.55, 0.12] } = {}) {
  if (bmp >= bmpT[0]) return 5;
  const v = frenchFlag(shh, shhT); // 0 FP, 1 p3, 2 pMN, 3 below all Shh thresholds
  if (v < 3) return v;
  return bmp >= bmpT[1] ? 4 : 3;
}

export const REACTIONS = [
  { id: 'grad.shh.sec', order: 4, pathway: 'Morphogen gradient', reactants: [], products: ['Shh'], modifiers: ['Floor plate', 'Notochord'], label: 'Shh secretion by floor plate / notochord' },
  { id: 'grad.shh.deg', order: 4, pathway: 'Morphogen gradient', reactants: ['Shh', 'PTCH1'], products: ['Shh:PTCH1'], label: 'Shh capture by PTCH1 (ligand-dependent degradation)' },
  { id: 'grad.shh.gli', order: 4, pathway: 'Morphogen gradient', reactants: ['GLI3R'], products: ['GliA (nuc)'], modifiers: ['Shh:PTCH1'], label: 'Shh signalling shifts Gli repressor to activator' },
  { id: 'grad.bmp.sec', order: 4, pathway: 'Morphogen gradient', reactants: [], products: ['BMP'], modifiers: ['Roof plate'], label: 'BMP secretion by roof plate' },
  { id: 'grad.bmp.smad', order: 4, pathway: 'Morphogen gradient', reactants: ['BMP', 'BMPR'], products: ['p-Smad1/5/8'], label: 'BMP receptor activation (Smad1/5/8 phosphorylation)' },
  { id: 'grad.nkx22', order: 4, pathway: 'Morphogen gradient', reactants: [], products: ['Nkx2.2'], modifiers: ['GliA (nuc)'], label: 'High Shh: Nkx2.2 (p3)' },
  { id: 'grad.olig2', order: 4, pathway: 'Morphogen gradient', reactants: [], products: ['Olig2'], modifiers: ['GliA (nuc)', 'Nkx2.2'], label: 'Intermediate Shh: Olig2 (pMN), repressed by Nkx2.2' },
  { id: 'grad.pax7', order: 4, pathway: 'Morphogen gradient', reactants: [], products: ['Pax7'], modifiers: ['p-Smad1/5/8', 'GLI3R'], label: 'Dorsal identity: Pax7' },
];
