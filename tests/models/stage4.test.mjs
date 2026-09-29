import test from 'node:test';
import assert from 'node:assert/strict';
import { RNG } from '../../js/models/core/rng.js';
import { integrate } from '../../js/models/core/ode.js';
import {
  laplacian2D, laplacian3D, gridOperator2D, icosphere, meshLaplacian, realSH, sphericalPowerSpectrum,
  schnakenberg, grayScott, giererMeinhardt, analyseModel, turingAnalysis, dispersion, sphereModes,
  RDSolver, variance, dominantWavelength2D, stableDt, substeps,
} from '../../js/models/morpho/reactionDiffusion.js';
import {
  decayLength, steady1D, steady1DFinite, solve1D, GradientGrid3D, pointSource3D, boundaryPositions,
  frenchFlag, positionalError, accumulationTime, neuralTubeFate, KICHEVA_DPP, REACTIONS as GRAD_R,
} from '../../js/models/morpho/gradient.js';
import { collagen, MigrationSim, ECMGrid, PointSourceField, REACTIONS as ECM_R } from '../../js/models/morpho/ecm.js';
import {
  enhancerStates, fateProbabilities, fateIndex, crosstalk, CryptModel, cryptLattice, FATES, REACTIONS as XT_R,
} from '../../js/models/pathways/crosstalk.js';

// ------------------------------------------------------------- Turing analysis

test('Schnakenberg: steady state and Jacobian match the closed forms', () => {
  const p = { a: 0.1, b: 0.9, gamma: 1 };
  const [[u, v]] = schnakenberg.steadyStates(p);
  assert.ok(Math.abs(u - 1) < 1e-12 && Math.abs(v - 0.9) < 1e-12);
  const J = schnakenberg.jacobian(p);
  // fu = (b - a)/(a + b), fv = (a + b)^2, gu = -2b/(a + b), gv = -(a + b)^2
  assert.deepEqual(J.map((r) => r.map((x) => +x.toFixed(12))), [[0.8, 1], [-1.8, -1]]);
});

test('Turing analysis flags diffusion-driven instability for a=0.1, b=0.9, Du=1, Dv=40 and stability for Dv=Du', () => {
  const p = { a: 0.1, b: 0.9, gamma: 1 };
  const an = analyseModel(schnakenberg, p, 1, 40);
  assert.ok(an.homogeneousStable, 'stable without diffusion');
  assert.ok(an.turingUnstable, 'Turing unstable');
  // band = roots of 40 k^4 - 31 k^2 + 1 = 0
  assert.ok(Math.abs(an.band[0] - (31 - Math.sqrt(801)) / 80) < 1e-12);
  assert.ok(Math.abs(an.band[1] - (31 + Math.sqrt(801)) / 80) < 1e-12);
  assert.ok(an.fastest.k2 > an.band[0] && an.fastest.k2 < an.band[1]);
  assert.ok(an.fastest.rate > 0);
  assert.ok(dispersion(an.J, 1, 40, 0).max < 0, 'k = 0 mode decays');
  assert.ok(dispersion(an.J, 1, 40, an.band[1] * 1.2).max < 0, 'short waves decay');
  // the numerically found maximum is a maximum
  assert.ok(an.growth(an.fastest.k2 * 1.05) < an.fastest.rate && an.growth(an.fastest.k2 * 0.95) < an.fastest.rate);

  const eq = analyseModel(schnakenberg, p, 1, 1);
  assert.equal(eq.turingUnstable, false);
  assert.equal(eq.band, null);
  for (let k2 = 0; k2 < 5; k2 += 0.05) assert.ok(eq.growth(k2) < 0);
});

test('Gierer-Meinhardt is Turing unstable for its defaults; Gray-Scott trivial state is linearly stable', () => {
  const gm = analyseModel(giererMeinhardt, giererMeinhardt.defaults);
  assert.ok(gm.turingUnstable);
  const [[a, h]] = giererMeinhardt.steadyStates(giererMeinhardt.defaults);
  const f = new Float64Array(1), g = new Float64Array(1);
  giererMeinhardt.kernel(giererMeinhardt.defaults)([a], [h], f, g, 1);
  assert.ok(Math.abs(f[0]) < 1e-9 && Math.abs(g[0]) < 1e-9, 'steady state is a zero of the kinetics');
  const gs = turingAnalysis(grayScott.jacobian(grayScott.defaults, 1, 0), 1, 0.5);
  assert.ok(gs.homogeneousStable && !gs.turingUnstable);
});

test('sphere modes: accessible wavenumbers are l(l+1)/R^2', () => {
  const an = analyseModel(schnakenberg, { a: 0.1, b: 0.9, gamma: 1 }, 1, 40);
  const R = 20;
  const modes = sphereModes(an, R, 20);
  modes.forEach((m) => assert.ok(Math.abs(m.k2 - (m.l * (m.l + 1)) / (R * R)) < 1e-12));
  assert.ok(modes.some((m) => m.rate > 0));
});

// ------------------------------------------------------------- Laplacians

test('grid Laplacian: linear field ~0 in the interior, x^2 + y^2 -> 4, x^2+y^2+z^2 -> 6', () => {
  const nx = 20, ny = 16, h = 0.5;
  const lin = new Float64Array(nx * ny), quad = new Float64Array(nx * ny), out = new Float64Array(nx * ny);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const X = x * h, Y = y * h;
    lin[y * nx + x] = 2 * X - 3 * Y + 1;
    quad[y * nx + x] = X * X + Y * Y;
  }
  laplacian2D(lin, nx, ny, out, h);
  for (let y = 1; y < ny - 1; y++) for (let x = 1; x < nx - 1; x++) assert.ok(Math.abs(out[y * nx + x]) < 1e-9);
  laplacian2D(quad, nx, ny, out, h);
  for (let y = 1; y < ny - 1; y++) for (let x = 1; x < nx - 1; x++) assert.ok(Math.abs(out[y * nx + x] - 4) < 1e-9);
  // zero-flux boundaries conserve the total (sum of the Laplacian = 0)
  assert.ok(Math.abs(out.reduce((s, v) => s + v, 0)) < 1e-6);

  const n3 = 8, q3 = new Float64Array(n3 ** 3), o3 = new Float64Array(n3 ** 3);
  for (let z = 0; z < n3; z++) for (let y = 0; y < n3; y++) for (let x = 0; x < n3; x++) q3[(z * n3 + y) * n3 + x] = x * x + y * y + z * z;
  laplacian3D(q3, n3, n3, n3, o3, 1);
  for (let z = 1; z < n3 - 1; z++) for (let y = 1; y < n3 - 1; y++) for (let x = 1; x < n3 - 1; x++) assert.ok(Math.abs(o3[(z * n3 + y) * n3 + x] - 6) < 1e-9);
});

test('icosphere + cotangent Laplacian: L Y_lm = -l(l+1)/R^2 Y_lm (Meyer et al. mixed areas)', () => {
  const R = 2;
  const s = icosphere(4, R);
  assert.equal(s.n, 2562);
  const L = meshLaplacian(s.positions, s.indices);
  assert.ok(Math.abs(L.totalArea / (4 * Math.PI * R * R) - 1) < 0.003, 'mixed areas tile the sphere');
  const lmax = 4, Y = new Float64Array((lmax + 1) ** 2);
  for (const [l, m] of [[1, 1], [2, 0], [3, -2], [4, 3]]) {
    const u = new Float64Array(s.n), out = new Float64Array(s.n);
    for (let i = 0; i < s.n; i++) { realSH(s.positions[3 * i], s.positions[3 * i + 1], s.positions[3 * i + 2], lmax, Y); u[i] = Y[l * l + l + m]; }
    L.apply(u, out);
    const ev = -(l * (l + 1)) / (R * R);
    let num = 0, den = 0;
    for (let i = 0; i < s.n; i++) { num += L.area[i] * (out[i] - ev * u[i]) ** 2; den += L.area[i] * (ev * u[i]) ** 2; }
    assert.ok(Math.sqrt(num / den) < 0.02, `l=${l}: relative error ${Math.sqrt(num / den)}`);
  }
  // constants are in the kernel
  const one = new Float64Array(s.n).fill(3), o = new Float64Array(s.n);
  L.apply(one, o);
  assert.ok(o.every((v) => Math.abs(v) < 1e-9));
});

test('real spherical harmonics are orthonormal under the mesh quadrature; power spectrum finds the degree', () => {
  const s = icosphere(4, 1);
  const L = meshLaplacian(s.positions, s.indices);
  const lmax = 5, K = (lmax + 1) ** 2, Y = new Float64Array(K);
  const G = Array.from({ length: K }, () => new Float64Array(K));
  const u = new Float64Array(s.n);
  for (let i = 0; i < s.n; i++) {
    realSH(s.positions[3 * i], s.positions[3 * i + 1], s.positions[3 * i + 2], lmax, Y);
    for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) G[a][b] += L.area[i] * Y[a] * Y[b];
    u[i] = Y[4 * 4 + 4 + 2] + 0.2 * Y[1];
  }
  for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) assert.ok(Math.abs(G[a][b] - (a === b ? 1 : 0)) < 0.01);
  const P = sphericalPowerSpectrum(u, s.positions, L.area, 8);
  assert.ok(Math.abs(P[4] - 1) < 0.02 && Math.abs(P[1] - 0.04) < 0.005 && P[3] < 1e-3);
});

test('uniform graph Laplacian fallback: kernel = constants, exact for x^2 + y^2 on a regular triangular lattice', () => {
  // flat hexagonal patch: 7 x 7 triangular lattice
  const pos = [], idx = [], W = 7;
  for (let j = 0; j < W; j++) for (let i = 0; i < W; i++) pos.push(i + 0.5 * j, (j * Math.sqrt(3)) / 2, 0);
  for (let j = 0; j < W - 1; j++) for (let i = 0; i < W - 1; i++) {
    const a = j * W + i;
    idx.push(a, a + 1, a + W, a + 1, a + W + 1, a + W);
  }
  const P = Float64Array.from(pos), I = Uint32Array.from(idx);
  for (const type of ['uniform', 'cotan']) {
    const L = meshLaplacian(P, I, { type });
    const u = new Float64Array(W * W), out = new Float64Array(W * W);
    for (let k = 0; k < W * W; k++) u[k] = P[3 * k] ** 2 + P[3 * k + 1] ** 2;
    L.apply(u, out);
    const centre = 3 * W + 3;
    assert.ok(Math.abs(out[centre] - 4) < 1e-9, `${type}: ${out[centre]}`);
  }
});

// ---------------------------------------------------- pattern formation

test('grid Schnakenberg: noise grows into a pattern with the predicted wavelength', () => {
  const p = { a: 0.1, b: 0.9, gamma: 2, Du: 1, Dv: 40 };
  const an = analyseModel(schnakenberg, p);
  const N = 64;
  const s = new RDSolver({ op: gridOperator2D({ nx: N, ny: N, h: 1, bc: 'periodic' }), model: schnakenberg, params: p });
  s.seed(new RNG(1), 0.01);
  const v0 = variance(s.u);
  // CFL: explicit step respects dt <= h^2 / (4 Dv)
  assert.ok(s.maxDt() <= 1 / (4 * 40) + 1e-12);
  s.advance(20);
  const v1 = variance(s.u);
  assert.ok(v1 > 1e3 * v0, `variance ${v0} -> ${v1}`);
  const { wavelength } = dominantWavelength2D(s.u, N, N, 1);
  const lo = (2 * Math.PI) / Math.sqrt(an.band[1]), hi = (2 * Math.PI) / Math.sqrt(an.band[0]);
  assert.ok(wavelength > lo && wavelength < hi, 'inside the unstable band');
  assert.ok(Math.abs(wavelength / an.fastest.wavelength - 1) < 0.2, `measured ${wavelength} vs fastest-growing ${an.fastest.wavelength}`);
  // no pattern when Dv = Du
  const eq = new RDSolver({ op: gridOperator2D({ nx: 32, ny: 32, h: 1, bc: 'periodic' }), model: schnakenberg, params: { ...p, Dv: 1 } });
  eq.seed(new RNG(1), 0.01);
  const e0 = variance(eq.u);
  eq.advance(20);
  assert.ok(variance(eq.u) < e0);
});

test('IMEX solver on an icosphere selects the fastest-growing spherical-harmonic degree', () => {
  const p = { a: 0.1, b: 0.9, gamma: 1, Du: 1, Dv: 40 };
  const an = analyseModel(schnakenberg, p);
  const R = Math.sqrt(20 / an.fastest.k2); // puts the fastest mode near l = 4
  const sph = icosphere(3, R);
  const L = meshLaplacian(sph.positions, sph.indices);
  const s = new RDSolver({ op: L, model: schnakenberg, params: p, scheme: 'imex' });
  s.seed(new RNG(4), 0.02);
  const modes = sphereModes(an, R, 10);
  const lBest = modes.reduce((a, b) => (b.rate > a.rate ? b : a)).l;
  for (let k = 0; k < 150; k++) s.step(0.1);
  const P = sphericalPowerSpectrum(s.u, sph.positions, L.area, 10);
  let lDom = 1;
  for (let l = 1; l <= 10; l++) if (P[l] > P[lDom]) lDom = l;
  assert.ok(variance(s.u) > 1e-3);
  assert.ok(Math.abs(lDom - lBest) <= 1, `dominant l ${lDom}, predicted ${lBest}`);
});

test('CFL helpers', () => {
  assert.ok(Math.abs(stableDt({ Dmax: 2, rho: 8, safety: 1 }) - 0.125) < 1e-12);
  const { n, dt } = substeps(1, 0.3);
  assert.equal(n, 4);
  assert.ok(Math.abs(dt - 0.25) < 1e-12);
});

// ------------------------------------------------------------- gradients

test('1D morphogen steady state matches C0 exp(-x / lambda), lambda = sqrt(D/k)', () => {
  const { D, k } = KICHEVA_DPP;
  const lam = decayLength(D, k);
  assert.ok(Math.abs(lam - 19.92) < 0.05, `Kicheva Dpp lambda ${lam}`);
  const L = 12 * lam;
  const { x, C } = solve1D({ n: 1200, L, D, k, C0: 1 });
  for (let i = 0; i < x.length; i += 50) {
    assert.ok(Math.abs(C[i] - steady1DFinite(x[i], { C0: 1, lambda: lam, L })) < 2e-4);
    if (x[i] < 6 * lam) assert.ok(Math.abs(C[i] - steady1D(x[i], { lambda: lam })) < 1e-3);
  }
  // French flag boundaries sit at lambda ln(C0 / T)
  const T = [0.5, 0.2];
  const xb = boundaryPositions(1, lam, T);
  xb.forEach((b, i) => assert.ok(Math.abs(steady1D(b, { lambda: lam }) - T[i]) < 1e-12));
  assert.equal(frenchFlag(0.7, T), 0);
  assert.equal(frenchFlag(0.3, T), 1);
  assert.equal(frenchFlag(0.1, T), 2);
  assert.ok(Math.abs(positionalError({ lambda: lam, cv: 0.1 }) - lam * 0.1) < 1e-12);
  assert.ok(accumulationTime(lam, { D, k }) > accumulationTime(0, { D, k }));
});

test('3D grid solver relaxes to the analytic profile from a source plane', () => {
  const h = 2, D = 1, k = 0.01; // lambda = 10
  const g = new GradientGrid3D({ nx: 60, ny: 3, nz: 3, h, D, k });
  for (let z = 0; z < 3; z++) for (let y = 0; y < 3; y++) { const i = g.index(0, y, z); g.fixed[i] = 1; g.fixedValue[i] = 1; g.C[i] = 1; }
  g.advance(1500);
  const lam = decayLength(D, k), L = 59 * h + h / 2;
  for (let x = 1; x < 50; x += 4) {
    const exact = steady1DFinite(x * h, { C0: 1, lambda: lam, L });
    assert.ok(Math.abs(g.C[g.index(x, 1, 1)] - exact) < 0.01, `x=${x}: ${g.C[g.index(x, 1, 1)]} vs ${exact}`);
  }
});

test('3D point source solves D lap C - k C = 0 away from the source', () => {
  const o = { Q: 1, D: 2, k: 0.05 };
  const c = (x, y, z) => pointSource3D(Math.hypot(x, y, z), o);
  const e = 1e-3;
  for (const r of [2, 5, 9]) {
    const lap = (c(r + e, 0, 0) + c(r - e, 0, 0) + c(r, e, 0) + c(r, -e, 0) + c(r, 0, e) + c(r, 0, -e) - 6 * c(r, 0, 0)) / (e * e);
    assert.ok(Math.abs(o.D * lap - o.k * c(r, 0, 0)) < 1e-4 * c(r, 0, 0) / e, `r=${r}`);
  }
});

test('neural-tube read-out orders domains ventral -> dorsal', () => {
  const lam = 20;
  const fates = [];
  for (let y = 0; y <= 200; y += 10) fates.push(neuralTubeFate(Math.exp(-y / lam), Math.exp(-(200 - y) / 40)));
  for (let i = 1; i < fates.length; i++) assert.ok(fates[i] >= fates[i - 1], fates.join(','));
  assert.equal(fates[0], 0);
  assert.equal(fates.at(-1), 5);
});

// ------------------------------------------------------------- ECM & migration

test('collagen fibrillogenesis conserves collagen mass; propeptides = cleaved procollagen', () => {
  const p = { ...collagen.defaults, sigma: 0 };
  const { y } = integrate(collagen.rhs(p), [1, 0, 0, 0], 0, 30, 0.01, { every: 100 });
  for (const s of y) assert.ok(Math.abs(collagen.mass(s) - 1) < 1e-9);
  const end = y.at(-1);
  assert.ok(end[2] > 0.95, 'nearly all collagen ends up in fibrils');
  assert.ok(Math.abs(end[3] - (1 - end[0])) < 1e-9);
  // with secretion the mass grows exactly as sigma t
  const q = { ...collagen.defaults, sigma: 0.7 };
  const r = integrate(collagen.rhs(q), [0, 0, 0, 0], 0, 10, 0.01).y.at(-1);
  assert.ok(Math.abs(collagen.mass(r) - 7) < 1e-8);
  // nucleation-limited lag: fibril mass starts with zero slope and curvature > 0
  const early = integrate(collagen.rhs(q), [0, 0, 0, 0], 0, 0.5, 0.001, { every: 100 }).y;
  assert.ok(early[1][2] < 0.05 * early[1][1], 'fibrils lag behind tropocollagen');
});

test('chemotaxis moves a cell up the chemoattractant gradient', () => {
  const sim = new MigrationSim({ nFibro: 0, nMigr: 1, initialFibres: 0, rng: new RNG(3), params: { Dn: 0, rho: 0, mmp: 0 } });
  const c = sim.migr[0];
  let prev = Infinity;
  for (let i = 0; i < 200; i++) {
    const d = Math.hypot(c.p[0] - sim.source.p[0], c.p[1] - sim.source.p[1], c.p[2] - sim.source.p[2]);
    assert.ok(d <= prev + 1e-9, 'distance never increases');
    prev = d;
    sim.step(0.05);
  }
  assert.ok(sim.meanDistance() < 0.8 * (Math.hypot(...c.start.map((v, k) => v - sim.source.p[k])) - sim.source.radius));
  // no chemotaxis -> no drift
  const still = new MigrationSim({ nFibro: 0, nMigr: 1, initialFibres: 0, rng: new RNG(3), params: { Dn: 0, rho: 0, chi: 0, mmp: 0 } });
  const p0 = still.migr[0].p.slice();
  for (let i = 0; i < 50; i++) still.step(0.05);
  assert.deepEqual(still.migr[0].p, p0);
});

test('haptotaxis moves a cell up the ECM density gradient', () => {
  const g = new ECMGrid({ size: [200, 40, 40], h: 10 });
  // ramp of fibres: density increases with x
  const fibres = [];
  for (let x = -95; x <= 95; x += 5) for (let n = 0; n < Math.round((x + 100) / 10); n++) fibres.push({ p: [x, 0, 0], dir: [0, 0, 1], len: 5, mass: 1 });
  g.rebuild(fibres);
  const gr = g.gradient([0, 0, 0]);
  assert.ok(gr[0] > 0 && Math.abs(gr[1]) < 1e-9);
  const field = new PointSourceField({ p: [1e4, 0, 0] });
  assert.ok(field.value([0, 0, 0]) < 1e-10, 'far source: negligible chemoattractant');
});

test('fibroblasts lay fibres along their direction of motion', () => {
  const sim = new MigrationSim({ nFibro: 3, nMigr: 0, initialFibres: 0, rng: new RNG(9), params: { fibroPersistence: 1e6, quantum: 0.2 } });
  for (const f of sim.fibro) { f.v = [15, 0, 0]; f.y[2] = 0; }
  for (let i = 0; i < 200; i++) sim.step(0.02);
  assert.ok(sim.fibres.length > 6, `fibres ${sim.fibres.length}`);
  for (const fib of sim.fibres) assert.ok(Math.abs(fib.dir[0]) > 0.99, `fibre dir ${fib.dir}`);
  assert.ok(sim.ecm.mean() > 0);
});

// ------------------------------------------------------------- Wnt-Notch crosstalk

test('thermodynamic enhancer: states sum to 1 and cooperativity favours co-occupancy', () => {
  const a = enhancerStates(2, 3, 1), b = enhancerStates(2, 3, 10);
  for (const s of [a, b]) assert.ok(Math.abs(s.none + s.wnt + s.notch + s.both - 1) < 1e-12);
  assert.ok(b.both > a.both);
});

test('combinatorial fate logic (Wnt x Notch) maps to crypt fates', () => {
  const p = crosstalk.defaults;
  assert.equal(FATES[fateIndex(1.5, 1.5, p)].id, 'stem');
  assert.equal(FATES[fateIndex(1.5, 0.0, p)].id, 'paneth');
  assert.equal(FATES[fateIndex(0.1, 1.5, p)].id, 'absorptive');
  assert.equal(FATES[fateIndex(0.1, 0.0, p)].id, 'differentiated');
  const pr = fateProbabilities(0.8, 0.8, p);
  assert.ok(Math.abs(pr.reduce((s, v) => s + v, 0) - 1) < 1e-12);
});

test('crosstalk: Wnt induces Jagged1 and, through Dishevelled, damps NICD', () => {
  const off = crosstalk.steadyCell(0, 0.6), on = crosstalk.steadyCell(1.5, 0.6);
  assert.ok(on[1] > 2 * off[1], 'beta-catenin stabilised');
  assert.ok(on[6] > 3 * off[6], 'Jagged1 induced');
  assert.ok(on[2] < off[2], 'NICD reduced by Dvl');
  const noX = { ...crosstalk.defaults, xDvlNICD: 0 };
  const a = crosstalk.steadyCell(0, 0.6, noX), b = crosstalk.steadyCell(1.5, 0.6, noX);
  assert.ok(Math.abs(a[2] - b[2]) < 1e-6, 'no crosstalk -> NICD independent of Wnt');
});

function crypt(params = {}) {
  const lat = cryptLattice({ perRing: 18, rings: 16, radius: 22 });
  const wnt = lat.arc.map((s) => 1.5 * Math.exp(-s / 30));
  const m = new CryptModel({ nbrs: lat.nbrs, wnt, rng: new RNG(5), params });
  m.advance(60);
  return { lat, m, f: m.fates() };
}

test('crypt: lateral inhibition + Wnt gradient give stem/Paneth at the base, absorptive/goblet above', () => {
  const { lat, m, f } = crypt();
  const L = lat.length;
  const frac = (lo, hi, ids) => {
    let n = 0, c = 0;
    f.forEach((k, i) => { const a = lat.arc[i] / L; if (a >= lo && a < hi) { n++; if (ids.includes(k)) c++; } });
    return c / n;
  };
  assert.ok(frac(0, 0.25, [0, 1]) > 0.9, 'stem + Paneth at the base');
  assert.ok(frac(0, 0.25, [1]) > 0.15, 'Paneth cells present at the base');
  assert.ok(frac(0.75, 1.01, [2, 3]) > 0.9, 'absorptive + goblet at the top');
  // lateral inhibition: Atoh1-high (secretory) cells are not adjacent to one another
  let ss = 0, sec = 0;
  m.S.forEach((s, i) => { if (s[4] > 0.8) { sec++; for (const j of lat.nbrs[i]) if (m.S[j][4] > 0.8) ss++; } });
  assert.ok(sec / m.n > 0.15 && sec / m.n < 0.45, `secretory fraction ${sec / m.n}`);
  assert.equal(ss, 0);
});

test('crypt perturbations: gamma-secretase inhibition -> all secretory (van Es 2005); NICD gain -> none (Fre 2005)', () => {
  const dapt = crypt({ dapt: 1 }).m.fractions();
  assert.ok(dapt[1] + dapt[3] > 0.99, `DAPT fractions ${dapt}`);
  const nicd = crypt({ nicdOE: 2 }).m.fractions();
  assert.ok(nicd[1] + nicd[3] < 0.01, `NICD fractions ${nicd}`);
});

test('every Stage-4 REACTIONS entry is well formed (order 4, unique ids)', () => {
  const all = [...GRAD_R, ...ECM_R, ...XT_R];
  const ids = new Set();
  for (const r of all) {
    assert.equal(r.order, 4);
    assert.ok(r.id && r.pathway && r.label && Array.isArray(r.reactants) && Array.isArray(r.products));
    assert.ok(!ids.has(r.id), r.id);
    ids.add(r.id);
  }
});
