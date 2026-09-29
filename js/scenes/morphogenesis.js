// Stage 4 - morphogenesis: Turing organoids, morphogen gradients (French
// flag), differential-adhesion sorting, ECM deposition + chemotaxis /
// haptotaxis, and Wnt x Notch crosstalk in the intestinal crypt.
import * as THREE from 'three';
import { cellMaterial, lineMaterial, releaseTree } from '../engine/materials.js';
import { InstancePool } from '../engine/instances.js';
import { unitSphere } from '../engine/geometry.js';
import { fmt } from '../engine/chart.js';
import { RNG } from '../models/core/rng.js';
import {
  icosphere, meshLaplacian, meshNeighbours, sphericalPowerSpectrum, SpectrumAccumulator, schnakenberg, grayScott, giererMeinhardt,
  analyseModel, RDSolver, chooseScheme,
} from '../models/morpho/reactionDiffusion.js';
import {
  GradientGrid3D, decayLength, boundaryPositions, neuralTubeFate, NEURAL_TUBE_DOMAINS, KICHEVA_DPP,
  positionalError, accumulationTime, steady1D,
} from '../models/morpho/gradient.js';
import { MigrationSim, PointSourceField } from '../models/morpho/ecm.js';
import { Tissue, Cell } from '../models/morpho/agents3d.js';
import { CryptModel, cryptLattice, crosstalk, FATES, fateProbabilities } from '../models/pathways/crosstalk.js';

const VIEWS = [
  { value: 'organoid', label: 'Turing organoid · reaction–diffusion buds' },
  { value: 'gradient', label: 'Morphogen gradient & French flag · neural tube' },
  { value: 'sorting', label: 'Differential adhesion · cell sorting' },
  { value: 'ecm', label: 'ECM & migration · chemotaxis + haptotaxis' },
  { value: 'crypt', label: 'Crypt fate logic · Wnt × Notch crosstalk' },
];

// ---------------------------------------------------------------- helpers
const lin = (hex) => { const c = new THREE.Color(hex); return [c.r, c.g, c.b]; };
function makeRamp(stops) {
  const s = stops.map(([t, hex]) => [t, ...lin(hex)]);
  return (t, out, o = 0) => {
    const x = t <= 0 ? 0 : t >= 1 ? 1 : t;
    let k = 1;
    while (k < s.length - 1 && x > s[k][0]) k++;
    const a = s[k - 1], b = s[k];
    const f = (x - a[0]) / (b[0] - a[0] || 1);
    out[o] = a[1] + (b[1] - a[1]) * f;
    out[o + 1] = a[2] + (b[2] - a[2]) * f;
    out[o + 2] = a[3] + (b[3] - a[3]) * f;
    return out;
  };
}
const smooth = (x) => { const t = x < 0 ? 0 : x > 1 ? 1 : x; return t * t * (3 - 2 * t); };
const clamp01 = (x) => (x < 0 ? 0 : x > 1 ? 1 : x);

/** Per-vertex normals of an indexed triangle mesh (area weighted), written in place. */
function computeNormals(pos, idx, nrm) {
  nrm.fill(0);
  for (let f = 0; f < idx.length; f += 3) {
    const a = idx[f] * 3, b = idx[f + 1] * 3, c = idx[f + 2] * 3;
    const ux = pos[b] - pos[a], uy = pos[b + 1] - pos[a + 1], uz = pos[b + 2] - pos[a + 2];
    const vx = pos[c] - pos[a], vy = pos[c + 1] - pos[a + 1], vz = pos[c + 2] - pos[a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    for (const k of [a, b, c]) { nrm[k] += nx; nrm[k + 1] += ny; nrm[k + 2] += nz; }
  }
  for (let i = 0; i < nrm.length; i += 3) {
    const l = Math.hypot(nrm[i], nrm[i + 1], nrm[i + 2]) || 1;
    nrm[i] /= l; nrm[i + 1] /= l; nrm[i + 2] /= l;
  }
}

function boxLines(size, color = 0x3a4a61, opacity = 0.5) {
  const g = new THREE.EdgesGeometry(new THREE.BoxGeometry(...size));
  return new THREE.LineSegments(g, lineMaterial({ color, opacity }));
}

function groundGrid(size, y, divisions = 24) {
  const grid = new THREE.GridHelper(size, divisions, 0x2a3a50, 0x16202c);
  grid.material = lineMaterial({ color: 0x2a3a50, opacity: 0.45 });
  grid.position.y = y;
  return grid;
}

/** Unit cylinder along +x (for InstancePool.putOriented). */
function xCylinder(radial = 6) {
  const g = new THREE.CylinderGeometry(1, 1, 1, radial, 1, false);
  g.rotateZ(-Math.PI / 2);
  return g;
}

// ======================================================================
// (a) Turing organoid
// ======================================================================
function organoidView(env) {
  const { stage, panel, group } = env;
  const sph = icosphere(5, 1);
  const N = sph.n;
  const L1 = meshLaplacian(sph.positions, sph.indices);           // unit sphere
  const sph4 = icosphere(4, 1);                                     // the first 2562 level-5 vertices
  const area4 = meshLaplacian(sph4.positions, sph4.indices).area;  // quadrature weights for the spectrum
  const nbrs = meshNeighbours(N, sph.indices);
  const RV = 46;          // drawn radius (scene units)
  const SHELL = 0.14;     // epithelium thickness / RV

  const PRESETS = {
    schnakenberg: {
      label: 'Schnakenberg · spots → buds (Turing)', model: schnakenberg, turing: true,
      params: { a: 0.1, b: 0.9, gamma: 1, Du: 1, Dv: 40 }, targetL: 7, rate: 5, seed: 'noise', act: 'u',
      units: '1/γ',
    },
    gm: {
      label: 'Gierer–Meinhardt · activator–inhibitor peaks (Turing)', model: giererMeinhardt, turing: true,
      params: { ...giererMeinhardt.defaults }, targetL: 7, rate: 5, seed: 'noise', act: 'u', units: '1/μa',
    },
    gsSpots: {
      label: 'Gray–Scott · self-replicating spots (Pearson)', model: grayScott, turing: false,
      params: { ...grayScott.presets.spots, Du: 2.5e-4, Dv: 1.25e-4 }, R: 1, rate: 900, seed: 'patches', act: 'v', units: 'a.u.',
    },
    gsStripes: {
      label: 'Gray–Scott · stripes / worms (Pearson)', model: grayScott, turing: false,
      params: { ...grayScott.presets.stripes, Du: 2.5e-4, Dv: 1.25e-4 }, R: 1, rate: 900, seed: 'patches', act: 'v', units: 'a.u.', feature: 'ridge maxima',
    },
    gsLaby: {
      label: 'Gray–Scott · labyrinth (Pearson)', model: grayScott, turing: false,
      params: { ...grayScott.presets.labyrinth, Du: 2.5e-4, Dv: 1.25e-4 }, R: 1, rate: 900, seed: 'patches', act: 'v', units: 'a.u.', feature: 'ridge maxima',
    },
  };
  for (const pr of Object.values(PRESETS)) {
    pr.base = { ...pr.params };
    if (pr.turing) {
      // choose R (in diffusion lengths) so the fastest-growing mode sits at l ~ targetL
      const an = analyseModel(pr.model, pr.params);
      pr.R = Math.sqrt((pr.targetL * (pr.targetL + 1)) / an.fastest.k2);
    }
    pr.R0 = pr.R;
  }
  const S = { preset: 'schnakenberg', budHeight: 0.36, turntable: true, sizeFactor: 1 };
  let P = PRESETS[S.preset];
  let solver, an, lStar = 0, lDom = 0, buds = 0, stepDt = 0.1, spectrum = null, lastSpec = -1, specClock = 0;
  const rng = new RNG(20240607);

  // ------------------------------------------------ geometry (outer + apical surfaces)
  const makeGeo = () => {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(N * 3), 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute('normal', new THREE.BufferAttribute(new Float32Array(N * 3), 3).setUsage(THREE.DynamicDrawUsage));
    g.setIndex(new THREE.BufferAttribute(Uint32Array.from(sph.indices), 1));
    g.boundingSphere = new THREE.Sphere(new THREE.Vector3(), RV * 1.8);
    return g;
  };
  const geo = makeGeo(), geoIn = makeGeo();
  const colors = new THREE.BufferAttribute(new Float32Array(N * 3), 3).setUsage(THREE.DynamicDrawUsage);
  geo.setAttribute('color', colors);
  geoIn.setAttribute('color', colors);
  geoIn.setAttribute('normal', geo.attributes.normal); // apical surface is a radial offset: same normals
  const mat = cellMaterial({ role: 'field', vertexColors: true, side: THREE.DoubleSide, noise: 0.38, noiseScale: 0.42, gain: 0.42, he: [0.12, 0.5] });
  const matIn = cellMaterial({ role: 'field', vertexColors: true, side: THREE.DoubleSide, color: 0x8c86a8, noise: 0.1, gain: 0.1, he: [0.04, 0.3] });
  const mesh = new THREE.Mesh(geo, mat);
  const inner = new THREE.Mesh(geoIn, matIn);
  group.add(mesh, inner);
  group.add(groundGrid(420, -RV * 1.75, 28));
  const disp = new Float32Array(N);
  const target = new Float32Array(N), tmpT = new Float32Array(N);
  const act = new Float64Array(N);
  let lo = 0, hi = 1, cmin = 0, cmax = 1;

  const RAMPS = {
    physical: makeRamp([[0, '#3f4f94'], [0.3, '#6a58a8'], [0.55, '#c0508a'], [0.8, '#f58a5c'], [1, '#ffe2a6']]),
    confocal: makeRamp([[0, '#10081a'], [0.35, '#3c1040'], [0.62, '#1a8f4e'], [1, '#79ffae']]),
    histology: makeRamp([[0, '#ffffff'], [0.45, '#dcd2f6'], [1, '#5a5cff']]),
  };

  function buildSolver(reseed = true) {
    const p = { ...P.params };
    const R = P.R0 * S.sizeFactor;
    P.R = R;
    // solve on the unit sphere with D / R^2 (equivalent to a sphere of radius R)
    const eff = { ...p, Du: p.Du / (R * R), Dv: p.Dv / (R * R) };
    const scheme = chooseScheme(L1, P.model, eff);
    if (!solver || reseed || solver.model !== P.model || solver.scheme !== scheme) {
      const old = solver;
      solver = new RDSolver({ op: L1, model: P.model, params: eff, scheme });
      if (!reseed && old && old.model === P.model) { solver.u.set(old.u); solver.v.set(old.v); solver.t = old.t; }
    } else solver.setParams(eff);
    stepDt = solver.maxDt();
    an = analyseModel(P.model, p, p.Du, p.Dv, 0);
    lStar = 0;
    let best = -Infinity;
    for (let l = 1; l <= 60; l++) { const g = an.growth((l * (l + 1)) / (R * R)); if (g > best) { best = g; lStar = l; } }
    if (!(best > 0)) lStar = 0;
    if (reseed) seed();
  }

  function seed() {
    if (P.seed === 'noise') solver.seed(rng, 0.08);
    else {
      // Gray-Scott: trivial state (u, v) = (1, 0) plus a few autocatalyst patches
      solver.u.fill(1); solver.v.fill(0);
      for (let c = 0; c < 30; c++) {
        const i = rng.int(N);
        const x = sph.positions[3 * i], y = sph.positions[3 * i + 1], z = sph.positions[3 * i + 2];
        for (let j = 0; j < N; j++) {
          const d = Math.hypot(sph.positions[3 * j] - x, sph.positions[3 * j + 1] - y, sph.positions[3 * j + 2] - z);
          if (d < 0.09) { solver.u[j] = 0.5; solver.v[j] = 0.25 * (1 + 0.1 * rng.next()); }
        }
      }
      solver.t = 0;
    }
    disp.fill(0);
    lastSpec = -1;
    if (specAcc) specAcc.i = specAcc.n; // abandon a spectrum of the previous pattern
    specChart?.clear();
  }

  function activator() { return P.act === 'u' ? solver.u : solver.v; }

  function updateSurface(relax) {
    const a = activator();
    let mean = 0, mx = -Infinity, mn = Infinity;
    for (let i = 0; i < N; i++) { const x = a[i]; mean += x; if (x > mx) mx = x; if (x < mn) mn = x; }
    mean /= N;
    // smoothed display window (stable colours while the pattern evolves)
    const w = relax >= 1 ? 1 : 0.15;
    lo += (mean - lo) * w; hi += (Math.max(mx, mean + 1e-6) - hi) * w;
    cmin += (mn - cmin) * w; cmax += (Math.max(mx, mn + 1e-6) - cmax) * w;
    const pos = geo.attributes.position.array, posIn = geoIn.attributes.position.array, col = colors.array;
    const ramp = RAMPS[stage.mode] ?? RAMPS.physical;
    const span = Math.max(hi - lo, 1e-9), cspan = Math.max(cmax - cmin, 1e-9);
    const H = S.budHeight;
    // bud shape: dome mapping of the activator (steep flanks, rounded tip) then
    // four passes of neighbour averaging, so buds are smooth, finger-like protrusions
    // and small activator peaks do not become thin spikes
    for (let i = 0; i < N; i++) { act[i] = a[i]; target[i] = Math.sqrt(smooth(((a[i] - lo) / span - 0.05) / 0.6)); }
    for (let pass = 0; pass < 4; pass++) {
      for (let i = 0; i < N; i++) {
        const nb = nbrs[i];
        let s = target[i] * 2;
        for (let k = 0; k < nb.length; k++) s += target[nb[k]];
        tmpT[i] = s / (nb.length + 2);
      }
      target.set(tmpT);
    }
    for (let i = 0; i < N; i++) {
      disp[i] += (H * target[i] - disp[i]) * relax;
      const r = RV * (1 + disp[i]);
      const rin = r - RV * SHELL;
      const x = sph.positions[3 * i], y = sph.positions[3 * i + 1], z = sph.positions[3 * i + 2];
      pos[3 * i] = x * r; pos[3 * i + 1] = y * r; pos[3 * i + 2] = z * r;
      posIn[3 * i] = x * rin; posIn[3 * i + 1] = y * rin; posIn[3 * i + 2] = z * rin;
      ramp((a[i] - cmin) / cspan, col, 3 * i);
    }
    computeNormals(pos, sph.indices, geo.attributes.normal.array);
    geo.attributes.position.needsUpdate = true; geo.attributes.normal.needsUpdate = true;
    geoIn.attributes.position.needsUpdate = true;
    colors.needsUpdate = true;
  }

  function countBuds() {
    const a = act;
    let m = 0, s2 = 0;
    for (let i = 0; i < N; i++) m += a[i];
    m /= N;
    for (let i = 0; i < N; i++) s2 += (a[i] - m) ** 2;
    const sd = Math.sqrt(s2 / N);
    let mx = -Infinity, mn = Infinity;
    for (let i = 0; i < N; i++) { if (a[i] > mx) mx = a[i]; if (a[i] < mn) mn = a[i]; }
    // no buds until the pattern has a real amplitude (seed noise has many tiny maxima)
    if (mx - mn < 0.3 * Math.max(Math.abs(m), 1e-3)) return 0;
    let c = 0;
    for (let i = 0; i < N; i++) {
      if (a[i] < m + 0.6 * sd || a[i] - m < 0.25 * (mx - m)) continue;
      let top = true;
      for (const j of nbrs[i]) if (a[j] > a[i] || (a[j] === a[i] && j < i)) { top = false; break; }
      if (top) c++;
    }
    return c;
  }

  // the spectrum is accumulated over several frames (~650 vertices each) so the live view never stalls
  const specAcc = new SpectrumAccumulator(sph4.positions, area4, 24);
  function finishSpectrum(P24) {
    let mx = 0;
    lDom = 1;
    for (let l = 1; l <= 24; l++) { if (P24[l] > mx) { mx = P24[l]; lDom = l; } }
    spectrum = Array.from(P24, (v) => (mx > 0 ? v / mx : 0));
    return spectrum;
  }

  // ------------------------------------------------ panel
  let dispChart, specChart, ro;
  function drawDispersion() {
    if (!dispChart) return;
    // x = wavenumber k (1 / model length unit); dots = the discrete modes a sphere of radius R can hold,
    // k_l = sqrt(l (l + 1)) / R (Laplace-Beltrami eigenvalues)
    const R = P.R;
    const lmax = Math.max(24, Math.ceil(lStar * 2.2));
    const kmax = Math.sqrt(lmax * (lmax + 1)) / R;
    const pts = new Map();
    for (let l = 0; l <= lmax; l++) pts.set(Math.sqrt(l * (l + 1)) / R, l);
    const xs = [...new Set([...Array.from({ length: 241 }, (_, i) => (kmax * i) / 240), ...pts.keys()])].sort((a, b) => a - b);
    const line = xs.map((k) => an.growth(k * k));
    const dots = xs.map((k, i) => (pts.has(k) ? line[i] : NaN));
    const top = Math.max(...line.filter(Number.isFinite));
    dispChart.yRange = top > 0 ? [-top * 1.6, top * 1.35] : [Math.min(...line) * 0.4, Math.max(Math.abs(top), 1e-3) * 0.5];
    dispChart.set(xs, [line, dots]);
    dispChart.markers = lStar > 0 ? [{ x: Math.sqrt(lStar * (lStar + 1)) / R, color: '#ff6fae' }] : [];
  }

  function updateReadouts() {
    if (!ro) return;
    const band = an.band;
    ro.set('homogeneous state (u*, v*)', `(${fmt(an.steady[0])}, ${fmt(an.steady[1])})`);
    ro.set('Turing unstable', P.turing ? (an.turingUnstable ? 'yes' : 'no') : 'n/a (finite-amplitude)');
    const R = P.R;
    const toL = (k2) => (-1 + Math.sqrt(1 + 4 * k2 * R * R)) / 2;
    ro.set('unstable band (degree l)', band ? `${toL(band[0]).toFixed(1)} – ${toL(band[1]).toFixed(1)}` : '—');
    ro.set('fastest mode: k, wavelength', an.fastest.rate > 0 ? `${fmt(an.fastest.k)}, ${fmt(an.fastest.wavelength)}` : '—');
    ro.set('predicted l* (sphere)', lStar || '—');
    ro.set('organoid radius R', `${fmt(R)} (${fmt(S.sizeFactor)} × median)`);
    ro.set('solver', `${solver.scheme === 'imex' ? 'IMEX (implicit diffusion, PCG)' : 'explicit Euler (CFL)'} · Δt ${fmt(stepDt)}`);
  }

  function buildPanel() {
    panel.section('Turing organoid');
    panel.select({
      label: 'Kinetics / preset', value: S.preset,
      options: Object.entries(PRESETS).map(([value, p]) => ({ value, label: p.label })),
      onChange: (v) => { S.preset = v; P = PRESETS[v]; S.sizeFactor = 1; buildSolver(true); env.rebuildPanel(); env.legend(legendItems()); },
    });
    const p = P.params;
    const reparam = () => { buildSolver(false); drawDispersion(); updateReadouts(); };
    if (P.model === schnakenberg) {
      panel.slider({ label: 'Dv / Du (inhibitor diffuses faster)', min: 1, max: 100, value: p.Dv / p.Du, log: true, onChange: (v) => { p.Dv = v * p.Du; reparam(); } });
      panel.slider({ label: 'a (activator source)', min: 0.02, max: 0.4, value: p.a, onChange: (v) => { p.a = v; reparam(); } });
      panel.slider({ label: 'b (substrate source)', min: 0.3, max: 1.6, value: p.b, onChange: (v) => { p.b = v; reparam(); } });
    } else if (P.model === giererMeinhardt) {
      panel.slider({ label: 'Dh / Da (inhibitor diffuses faster)', min: 1, max: 100, value: p.Dv / p.Du, log: true, onChange: (v) => { p.Dv = v * p.Du; reparam(); } });
      panel.slider({ label: 'μh (inhibitor decay)', min: 1.05, max: 4, value: p.muH, onChange: (v) => { p.muH = v; reparam(); } });
      panel.slider({ label: 'κ (activator saturation)', min: 0, max: 0.3, value: p.kappa, onChange: (v) => { p.kappa = v; reparam(); } });
    } else {
      panel.slider({ label: 'F (feed rate)', min: 0.01, max: 0.08, step: 0.0005, value: p.F, format: (v) => v.toFixed(4), onChange: (v) => { p.F = v; reparam(); } });
      panel.slider({ label: 'k (removal rate)', min: 0.045, max: 0.07, step: 0.0005, value: p.k, format: (v) => v.toFixed(4), onChange: (v) => { p.k = v; reparam(); } });
    }
    panel.slider({
      label: 'Organoid size (× median)', min: 0.5, max: 2, value: S.sizeFactor,
      onChange: (v) => { S.sizeFactor = v; reparam(); },
    });
    panel.slider({ label: 'Bud growth (outward displacement)', min: 0, max: 0.7, value: S.budHeight, onChange: (v) => { S.budHeight = v; } });
    panel.toggle({ label: 'Turntable rotation', value: S.turntable, onChange: (v) => { S.turntable = v; stage.controls.autoRotate = v; } });
    panel.buttons([
      { label: 'Re-seed pattern', primary: true, onClick: () => { seed(); } },
      { label: 'New organoid (measured sizes)', onClick: () => newOrganoid() },
    ]);

    panel.section('Linear stability (Turing analysis)');
    panel.equation('∂u/∂t = Du ∇²u + f(u,v)\n∂v/∂t = Dv ∇²v + g(u,v)');
    if (P.model === schnakenberg) panel.equation('f = γ(a − u + u²v),  g = γ(b − u²v)\nu* = a + b,  v* = b/(a + b)²');
    else if (P.model === giererMeinhardt) panel.equation('f = ρa a²/(h(1 + κa²)) − μa a\ng = ρh a² − μh h');
    else panel.equation('f = −uv² + F(1 − u),  g = uv² − (F + k)v');
    panel.equation('λ(k²) = eig[ J − k² diag(Du, Dv) ]\nsphere:  k² = l(l + 1) / R²');
    dispChart = panel.chart({
      title: 'Dispersion Re λ(k); dots = modes the organoid can hold, k² = l(l+1)/R²', xLabel: 'k',
      series: [{ name: 'Re λ(k)', color: '#5ee0c1' }, { name: 'sphere modes l', color: '#ffb454', points: true, width: 0.01 }],
    });
    specChart = panel.chart({
      title: 'Measured angular power spectrum of the activator', xLabel: 'l', yRange: [0, 1.05],
      series: [{ name: 'measured P_l', color: '#ff6fae', points: true }, { name: 'linear theory (2l+1)e^{2Re λ t}', color: '#8a98ab', dash: [4, 3] }],
    });
    ro = panel.readouts(['homogeneous state (u*, v*)', 'Turing unstable', 'unstable band (degree l)', 'fastest mode: k, wavelength',
      'predicted l* (sphere)', 'measured dominant l', 'buds (activator maxima)', 'organoid radius R', 'solver']);
    panel.note(P.turing
      ? 'Only whole spherical harmonics fit on the organoid: the pattern selects the degree l whose k² = l(l+1)/R² grows fastest. Larger organoids hold higher l, i.e. more buds.'
      : 'Gray–Scott patterns grow from finite-amplitude seeds: the trivial state (1, 0) is linearly stable (all Re λ < 0), so there is no Turing band — spots self-replicate instead (Pearson 1993).');

    if (env.orgMeta) {
      const m = env.orgMeta;
      panel.section('Real organoids (reference)');
      panel.table(['OrgaQuant bright-field image (MIT)', 'value'], [
        ['organoids detected (Hough circles)', m.n_organoids],
        ['radius CV', m.radius_cv],
        ['p90 / median radius', m.radius_p90_over_median],
        ['median radius (px)', m.radius_px_median],
      ]);
      panel.note('“New organoid” draws a size from a log-normal with the measured CV, so the accessible modes and the bud count change with size as in the real culture. Pixel size is unknown, so only relative sizes are used.');
    }
    drawDispersion();
    updateReadouts();
    if (spectrum) pushSpectrum();
  }

  function newOrganoid() {
    const cv = env.orgMeta?.radius_cv ?? 0.33;
    const sigma = Math.sqrt(Math.log(1 + cv * cv));
    S.sizeFactor = Math.min(2, Math.max(0.5, Math.exp(sigma * rng.normal())));
    buildSolver(true);
    env.rebuildPanel();
  }

  function pushSpectrum() {
    const xs = spectrum.map((_, l) => l);
    let theory = null;
    if (P.turing) {
      const t = lastSpec >= 0 ? lastSpec : solver.t, R = P.R;
      const raw = xs.map((l) => (l === 0 ? 0 : (2 * l + 1) * Math.exp(Math.min(2 * an.growth((l * (l + 1)) / (R * R)) * t, 600))));
      const mx = Math.max(...raw);
      theory = raw.map((v) => v / mx);
    }
    specChart.set(xs, [spectrum, theory ?? xs.map(() => NaN)]);
    specChart.markers = lStar > 0 ? [{ x: lStar, color: '#ff6fae' }] : [];
    ro?.set('measured dominant l', lDom);
    ro?.set('buds (activator maxima)', buds);
  }

  // ------------------------------------------------ pre-roll so the first frame already shows a pattern
  buildSolver(true);
  {
    const t0 = performance.now();
    while (solver.t < 11 && performance.now() - t0 < 900) solver.step(stepDt);
    updateSurface(1);
    for (let k = 0; k < 3; k++) updateSurface(0.5);
    buds = countBuds();
    finishSpectrum(sphericalPowerSpectrum(activator().subarray(0, sph4.n), sph4.positions, area4, 24));
  }

  function legendItems() {
    const pal = {
      physical: ['#3f4f94', '#c0508a', '#ffe2a6'],
      confocal: ['#3c1040', '#1a8f4e', '#79ffae'],
      histology: ['#f0bfd8', '#b98fd0', '#7a4fb0'],
    }[stage.mode] ?? ['#3f4f94', '#c0508a', '#ffe2a6'];
    const sp = P.model.species;
    return [
      { color: pal[0], label: `${P.act === 'u' ? sp[0] : sp[1]}: low (inter-bud epithelium)` },
      { color: pal[1], label: 'intermediate' },
      { color: pal[2], label: 'high → bud tip (crypt-like domain)' },
    ];
  }

  return {
    frame() {
      stage.frame([0, 0, 0], RV * 1.75, new THREE.Vector3(0.55, 0.45, 1));
      stage.clipAxis.set(0, 0, 1);
      stage.clipRange = [-RV * 1.6, RV * 1.6];
      stage.controls.autoRotate = S.turntable;
      stage.controls.autoRotateSpeed = 0.6;
    },
    buildPanel,
    legend: legendItems,
    pickables: [{
      object: mesh,
      info: (hit) => {
        const i = hit.face?.a;
        if (i === undefined) return null;
        const sp = P.model.species;
        return `Organoid surface · vertex ${i}\n${sp[0]} ${fmt(solver.u[i])} · ${sp[1]} ${fmt(solver.v[i])}\n` +
          `bud height ${fmt(disp[i] * 100)} % of radius${disp[i] > 0.5 * S.budHeight ? ' (bud)' : ''}`;
      },
    }],
    update(dt) {
      if (dt > 0) {
        const tEnd = solver.t + dt * P.rate;
        const t0 = performance.now();
        // as many substeps as fit in ~7 ms of the frame (the RD never blocks rendering)
        while (solver.t < tEnd && performance.now() - t0 < 7) solver.step(stepDt);
        updateSurface(1 - Math.exp(-dt * 2.2));
        specClock += dt;
        if (!specAcc.done) {
          if (specAcc.step(650)) { finishSpectrum(specAcc.result()); buds = countBuds(); pushSpectrum(); }
        } else if (specClock > 0.5 || lastSpec < 0) {
          specClock = 0;
          lastSpec = solver.t;
          specAcc.begin(activator().subarray(0, sph4.n));
        }
      }
      const pat = P.turing ? `l* = ${lStar} predicted · dominant l = ${lDom}` : `dominant l = ${lDom}`;
      env.setStatus(`t = ${fmt(solver.t)} ${P.units} · ${pat} · ${buds} ${P.feature ?? 'buds'}`);
    },
    reset() { seed(); },
    onMode() { updateSurface(0); env.legend(legendItems()); },
    dispose() { stage.controls.autoRotate = false; },
  };
}

// ======================================================================
// (b) Morphogen gradients & French flag: neural tube
// ======================================================================
function gradientView(env) {
  const { stage, panel, group } = env;
  const AX = 55, AY = 100, AXL = 6, AYL = 78, LZ = 150, h = 5;
  const nx = 26, ny = 44, nz = Math.round(LZ / h);
  const origin = [-nx * h / 2, -ny * h / 2, -LZ / 2];
  const inTissue = (x, y) => (x * x) / (AX * AX) + (y * y) / (AY * AY) <= 1 && (x * x) / (AXL * AXL) + (y * y) / (AYL * AYL) > 1;
  const mask = new Uint8Array(nx * ny * nz);
  const tmp = [0, 0, 0];
  const S = {
    D: KICHEVA_DPP.D, k: KICHEVA_DPP.k, lamBMP: 30, shhT: [0.5, 0.25, 0.1], bmpT: [0.5, 0.08],
    bead: false, colorBy: 'fate', cv: 0.15, hoursPerSecond: 0.35,
  };
  const shh = new GradientGrid3D({ nx, ny, nz, h, D: S.D, k: S.k, origin, mask });
  const bmp = new GradientGrid3D({ nx, ny, nz, h, D: S.D, k: S.D / (S.lamBMP * S.lamBMP), origin, mask });
  const BEAD = [AX * 0.62, 28, 0], BEAD_R = 8;
  for (let i = 0; i < mask.length; i++) {
    shh.centre(i, tmp);
    mask[i] = inTissue(tmp[0], tmp[1]) ? 1 : 0;
  }
  function setSources() {
    for (let i = 0; i < mask.length; i++) {
      shh.centre(i, tmp);
      const [x, y, z] = tmp;
      const fp = mask[i] && y < -AY + 16 && Math.abs(x) < 18;
      const rp = mask[i] && y > AY - 16 && Math.abs(x) < 18;
      const bd = S.bead && Math.hypot(x - BEAD[0], y - BEAD[1], z - BEAD[2]) < BEAD_R;
      if (bd) mask[i] = 1;
      shh.fixed[i] = fp || bd ? 1 : 0; shh.fixedValue[i] = 1;
      bmp.fixed[i] = rp ? 1 : 0; bmp.fixedValue[i] = 1;
      if (shh.fixed[i]) shh.C[i] = 1;
      if (bmp.fixed[i]) bmp.C[i] = 1;
      if (!bd && !inTissue(x, y)) mask[i] = 0;
    }
  }
  setSources();

  // ------------------------------------------------ cells (pseudostratified neuroepithelium: nuclei)
  const rng = new RNG(77);
  const cells = [];
  const sp = 8.2;
  for (let z = -LZ / 2 + sp / 2; z < LZ / 2; z += sp) {
    for (let row = 0, y = -AY; y <= AY; y += sp * 0.866, row++) {
      for (let x = -AX + ((row % 2) * sp) / 2; x <= AX; x += sp) {
        const px = x + (rng.next() - 0.5) * 2.4, py = y + (rng.next() - 0.5) * 2.4, pz = z + (rng.next() - 0.5) * 3;
        if (!inTissue(px, py)) continue;
        // inside the tissue with a margin (nuclei are not cut by the surfaces)
        const e = (px * px) / ((AX - 3) ** 2) + (py * py) / ((AY - 3) ** 2);
        if (e > 1) continue;
        // apico-basal axis ~ normal of the ellipse family
        const ax = [px / (AX * AX), py / (AY * AY), 0];
        const n = Math.hypot(ax[0], ax[1]) || 1;
        cells.push({ p: [px, py, pz], axis: [ax[0] / n, ax[1] / n, 0], fate: 3, shh: 0, bmp: 0, s: 0.85 + 0.3 * rng.next() });
      }
    }
  }
  const nucPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', useInst: 1, noise: 0.5, gain: 0.2 }), cells.length + 8, group);
  const tubeOuter = new THREE.CylinderGeometry(1, 1, LZ, 72, 1, true);
  tubeOuter.rotateX(Math.PI / 2);
  const outerMesh = new THREE.Mesh(tubeOuter, cellMaterial({ role: 'membrane', color: 0x7fa6d9, stain: 0x5a7bd6, opacity: 0.04, rimPower: 3.5, gain: 0.3 }));
  outerMesh.scale.set(AX + 2, AY + 2, 1);
  const lumenMesh = new THREE.Mesh(tubeOuter.clone(), cellMaterial({ role: 'membrane', color: 0xd5c2ff, stain: 0x8a6cff, opacity: 0.08, gain: 0.3 }));
  lumenMesh.scale.set(AXL, AYL, 1);
  outerMesh.renderOrder = lumenMesh.renderOrder = 3;
  const noto = new THREE.Mesh(new THREE.CylinderGeometry(15, 15, LZ + 10, 40, 1), cellMaterial({ role: 'solid', color: 0x8fa3bf, stain: 0x4a6cff, noise: 0.4 }));
  noto.rotation.x = Math.PI / 2;
  noto.position.set(0, -AY - 22, 0);
  const bead = new THREE.Mesh(unitSphere(3), cellMaterial({ role: 'reporter', color: 0x9d7bff, emissive: 0.6 }));
  bead.scale.setScalar(BEAD_R);
  bead.position.set(...BEAD);
  bead.visible = S.bead;
  // ectoderm / surface context: a faint plate above the roof
  group.add(outerMesh, lumenMesh, noto, bead);
  group.add(groundGrid(520, -AY - 45, 26));

  // concentration section (vertex-coloured field) at the posterior end
  const secGeo = new THREE.PlaneGeometry(nx * h, ny * h, nx - 1, ny - 1);
  const secCol = new THREE.BufferAttribute(new Float32Array(secGeo.attributes.position.count * 3), 3);
  secGeo.setAttribute('color', secCol);
  const section = new THREE.Mesh(secGeo, cellMaterial({ role: 'field', vertexColors: true, side: THREE.DoubleSide }));
  section.position.set(0, 0, LZ / 2 + 30);
  section.visible = false;
  group.add(section);

  const domainCol = NEURAL_TUBE_DOMAINS.map((d) => new THREE.Color(d.color));
  const shhRamp = makeRamp([[0, '#0b0f2a'], [0.25, '#2b2f8f'], [0.6, '#8b5cff'], [1, '#e9dcff']]);
  const bmpRamp = makeRamp([[0, '#06140f'], [0.3, '#0f5f4c'], [0.7, '#38d6a8'], [1, '#d8fff1']]);
  const tmpC = new THREE.Color();
  const rgb = [0, 0, 0];

  function sampleCells() {
    for (const c of cells) {
      c.shh = shh.sample(c.p);
      c.bmp = bmp.sample(c.p);
      c.fate = neuralTubeFate(c.shh, c.bmp, { shhT: S.shhT, bmpT: S.bmpT });
    }
  }

  function draw() {
    nucPool.begin();
    for (const c of cells) {
      let col;
      if (S.colorBy === 'fate') col = domainCol[c.fate];
      else if (S.colorBy === 'shh') col = tmpC.fromArray(shhRamp(Math.sqrt(clamp01(c.shh)), rgb));
      else col = tmpC.fromArray(bmpRamp(Math.sqrt(clamp01(c.bmp)), rgb));
      nucPool.putOriented(c.p, c.axis, 3.6 * c.s, 2.3, 2.3, col);
    }
    nucPool.end();
    if (!section.visible) return;
    // section: Shh (violet) + BMP (teal) of the middle A-P slice
    const zi = Math.floor(nz / 2);
    const arr = secCol.array;
    const pos = secGeo.attributes.position;
    for (let v = 0; v < pos.count; v++) {
      const ix = v % nx, iy = ny - 1 - Math.floor(v / nx);
      const i = (zi * ny + iy) * nx + ix;
      if (!mask[i]) { arr[3 * v] = arr[3 * v + 1] = arr[3 * v + 2] = 0.02; continue; }
      shhRamp(Math.sqrt(clamp01(shh.C[i])), rgb);
      const r0 = rgb[0], g0 = rgb[1], b0 = rgb[2];
      bmpRamp(Math.sqrt(clamp01(bmp.C[i])), rgb);
      arr[3 * v] = Math.min(1, r0 + rgb[0]); arr[3 * v + 1] = Math.min(1, g0 + rgb[1]); arr[3 * v + 2] = Math.min(1, b0 + rgb[2]);
    }
    secCol.needsUpdate = true;
  }

  // ------------------------------------------------ profile along the lateral wall
  const wallX = (y) => {
    const xo = AX * Math.sqrt(Math.max(0, 1 - (y * y) / (AY * AY)));
    const xl = Math.abs(y) < AYL ? AXL * Math.sqrt(Math.max(0, 1 - (y * y) / (AYL * AYL))) : 0;
    return (xo + xl) / 2;
  };
  const Y0 = -AY + 16; // edge of the floor-plate source
  function profile() {
    const xs = [], cs = [], cb = [];
    for (let d = 0; d <= 2 * AY - 32; d += 4) {
      const y = Y0 + d;
      xs.push(d);
      cs.push(shh.sample([-wallX(y), y, -LZ / 4]));
      cb.push(bmp.sample([-wallX(y), y, -LZ / 4]));
    }
    return { xs, cs, cb };
  }
  const crossing = (xs, cs, T) => { for (let i = 0; i < cs.length; i++) if (cs[i] < T) return xs[i]; return NaN; };

  let profChart, bChart, ro;
  function buildPanel() {
    panel.section('Morphogen gradients (source–diffusion–degradation)');
    panel.equation('∂C/∂t = D ∇²C − k C   (source cells: C = C₀)\nsteady 1D:  C(x) = C₀ e^(−x/λ),  λ = √(D/k)');
    panel.slider({ label: 'Shh diffusion D (µm²/s)', min: 0.02, max: 1, value: S.D, log: true, onChange: (v) => { S.D = v; shh.D = v; bmp.D = v; bmp.k = v / (S.lamBMP ** 2); refresh(); } });
    panel.slider({ label: 'Shh degradation k (1/s)', min: 2e-5, max: 2e-3, value: S.k, log: true, format: (v) => v.toExponential(1), onChange: (v) => { S.k = v; shh.k = v; refresh(); } });
    panel.slider({ label: 'BMP decay length λ_BMP (µm)', min: 10, max: 80, value: S.lamBMP, onChange: (v) => { S.lamBMP = v; bmp.k = S.D / (v * v); refresh(); } });
    panel.note('Defaults: Dpp in the fly wing disc, D = 0.10 µm²/s and k = 2.5×10⁻⁴ s⁻¹ measured by FRAP (Kicheva et al. 2007) → λ ≈ 20 µm, the same order as the Shh gradient in the neural tube.');

    panel.section('French flag read-out (thresholds)');
    const names = ['floor plate / p3', 'p3 / pMN', 'pMN / p2–p0'];
    S.shhT.forEach((T, i) => panel.slider({
      label: `Shh threshold ${names[i]}`, min: 0.01, max: 0.95, value: T,
      onChange: (v) => { S.shhT[i] = v; S.shhT.sort((a, b) => b - a); refresh(); },
    }));
    panel.slider({ label: 'BMP threshold roof plate', min: 0.1, max: 0.95, value: S.bmpT[0], onChange: (v) => { S.bmpT[0] = v; refresh(); } });
    panel.slider({ label: 'BMP threshold dorsal (Pax7)', min: 0.01, max: 0.5, value: S.bmpT[1], onChange: (v) => { S.bmpT[1] = v; refresh(); } });
    panel.select({
      label: 'Colour nuclei by', value: S.colorBy,
      options: [{ value: 'fate', label: 'Progenitor domain (French flag)' }, { value: 'shh', label: 'Shh concentration' }, { value: 'bmp', label: 'BMP concentration' }],
      onChange: (v) => { S.colorBy = v; env.legend(legendItems()); draw(); },
    });
    panel.toggle({ label: 'Ectopic Shh source (bead in the dorsal wall)', value: S.bead, onChange: (v) => { S.bead = v; bead.visible = v; setSources(); } });
    panel.toggle({ label: 'Concentration section in front (Shh violet, BMP teal)', value: section.visible, onChange: (v) => { section.visible = v; draw(); } });
    panel.buttons([{ label: 'Restart (no morphogen)', primary: true, onClick: () => restart() }]);

    panel.section('Profiles');
    profChart = panel.chart({
      title: 'C / C₀ along the lateral wall vs distance from the floor plate', xLabel: 'µm',
      yRange: [0, 1.02],
      series: [
        { name: 'Shh (3D grid)', color: '#9d7bff' }, { name: 'C₀e^(−x/λ)', color: '#e9dcff', dash: [4, 3] },
        { name: 'BMP (3D grid)', color: '#38d6a8' },
        { name: 'thresholds', color: '#b9a6ff', dash: [2, 4], width: 1 }, { name: '', color: '#b9a6ff', dash: [2, 4], width: 1 }, { name: '', color: '#b9a6ff', dash: [2, 4], width: 1 },
        { name: '', color: '#7fe8c6', dash: [2, 4], width: 1 }, { name: '', color: '#7fe8c6', dash: [2, 4], width: 1 },
      ],
    });
    bChart = panel.chart({
      title: 'Domain boundaries vs time (dashed = steady state λ ln(C₀/T))', xLabel: 'h',
      series: [
        { name: 'p3/pMN', color: '#ff4fa3' }, { name: 'pMN/p2', color: '#ffb454' },
        { name: '', color: '#ff4fa3', dash: [4, 3], width: 1 }, { name: '', color: '#ffb454', dash: [4, 3], width: 1 },
      ],
    });
    ro = panel.readouts(['λ_Shh = √(D/k)', 'predicted boundaries (µm)', 'positional error σx = λ·CV', 'accumulation time at pMN/p2', 'cells: FP / p3 / pMN / p2–p0 / dP / RP']);
    panel.note('Positional error assumes a 15 % cell-to-cell CV of the read-out (σx = λ·CV for an exponential gradient; Bollenbach et al. 2008). Accumulation time τ(x) = (1 + x/λ)/(2k) (Berezhkovskii et al. 2010).');
    refresh();
  }

  function refresh() {
    const lam = decayLength(S.D, S.k);
    if (!ro) return;
    const xb = boundaryPositions(1, lam, S.shhT);
    ro.set('λ_Shh = √(D/k)', `${fmt(lam)} µm`);
    ro.set('predicted boundaries (µm)', xb.map((x) => fmt(x)).join(' / '));
    ro.set('positional error σx = λ·CV', `${fmt(positionalError({ lambda: lam, cv: S.cv }))} µm`);
    ro.set('accumulation time at pMN/p2', `${fmt(accumulationTime(xb[2], { D: S.D, k: S.k }) / 3600)} h`);
    bChart?.clear();
  }

  function restart() {
    shh.C.fill(0); bmp.C.fill(0); shh.t = bmp.t = 0;
    setSources();
    bChart?.clear();
  }

  function legendItems() {
    if (S.colorBy === 'fate') return NEURAL_TUBE_DOMAINS.map((d) => ({ color: d.color, label: d.label }));
    return [
      { color: S.colorBy === 'shh' ? '#8b5cff' : '#38d6a8', label: `${S.colorBy === 'shh' ? 'Shh (floor plate / notochord)' : 'BMP (roof plate)'} concentration` },
      { color: '#8fa3bf', label: 'notochord' },
      { color: '#d5c2ff', label: 'lumen (apical surface)' },
    ];
  }

  let sampleClock = 1, chartClock = 1;
  // start from an established gradient (~3 h, beyond the accumulation time of the pMN boundary);
  // 'Restart' shows it forming from zero
  shh.advance(3 * 3600, 2000);
  bmp.advance(3 * 3600, 2000);
  sampleCells();
  draw();

  return {
    frame() {
      stage.frame([0, -12, 0], 165, new THREE.Vector3(0.8, 0.4, 1));
      stage.clipAxis.set(0, 0, 1);
      stage.clipRange = [-LZ / 2 - 30, LZ / 2];
    },
    buildPanel,
    legend: legendItems,
    pickables: [{
      get object() { return nucPool.mesh; },
      info: (hit) => {
        const c = cells[hit.instanceId];
        if (!c) return null;
        const d = NEURAL_TUBE_DOMAINS[c.fate];
        return `Neuroepithelial progenitor · ${d.label}\nShh ${fmt(c.shh)} · BMP ${fmt(c.bmp)} (× source)\n` +
          `${fmt(c.p[1] - Y0)} µm dorsal of the floor plate`;
      },
    }, { object: noto, info: () => 'Notochord: source of Shh that induces the floor plate' }],
    update(dt) {
      if (dt > 0) {
        const T = dt * S.hoursPerSecond * 3600;
        shh.advance(T, 14);
        bmp.advance(T, 14);
      }
      sampleClock += dt;
      if (sampleClock > 0.1) { sampleClock = 0; sampleCells(); draw(); }
      chartClock += dt;
      if ((chartClock > 0.3 || dt === 0) && profChart) {
        chartClock = 0;
        const lam = decayLength(S.D, S.k);
        const { xs, cs, cb } = profile();
        profChart.set(xs, [cs, xs.map((x) => steady1D(x, { lambda: lam })), cb, ...S.shhT.map((T) => xs.map(() => T)), ...S.bmpT.map((T) => xs.map(() => T))]);
        const tH = shh.t / 3600;
        const xb = boundaryPositions(1, lam, S.shhT);
        if (dt > 0) bChart.push(tH, [crossing(xs, cs, S.shhT[1]), crossing(xs, cs, S.shhT[2]), xb[1], xb[2]]);
        const counts = [0, 0, 0, 0, 0, 0];
        for (const c of cells) counts[c.fate]++;
        ro?.set('cells: FP / p3 / pMN / p2–p0 / dP / RP', counts.join(' / '));
      }
      env.setStatus(`t = ${fmt(shh.t / 3600)} h · λ = ${fmt(decayLength(S.D, S.k))} µm · ${cells.length} progenitors`);
    },
    reset() { restart(); },
    onMode() { draw(); },
  };
}

// ======================================================================
// (c) Differential adhesion sorting (Steinberg)
// ======================================================================
function sortingView(env) {
  const { stage, panel, group } = env;
  const PRESETS = {
    engulf: { label: 'Engulfment: most cohesive type sorts inside', types: 2, J: [[8, 3, 0], [3, 1, 0], [0, 0, 0]] },
    mix: { label: 'Mixing: heterotypic > homotypic (checkerboard)', types: 2, J: [[3, 8, 0], [8, 3, 0], [0, 0, 0]] },
    separate: { label: 'Separation: no heterotypic adhesion', types: 2, J: [[6, 0.3, 0], [0.3, 6, 0], [0, 0, 0]] },
    layers: { label: 'Three types: layered (A core, B middle, C outer)', types: 3, J: [[9, 5, 1], [5, 4, 1.8], [1, 1.8, 1]] },
  };
  const COLORS = [new THREE.Color('#39e0a0'), new THREE.Color('#ff4fa3'), new THREE.Color('#ffc94a')];
  const TYPE_NAMES = ['A (E-cadherin high)', 'B (E-cadherin low)', 'C (weakly adhesive)'];
  const S = { preset: 'engulf', noise: 8, n: 360, cutaway: true };
  let P = PRESETS[S.preset];
  let J = P.J.map((r) => r.slice());
  let tissue, rng, lastSample = -1;
  const R = 5, BOUND = 46;

  // physical: opaque cells coloured by type; confocal / H&E: membranes + nuclei (cell-tracker dyes / section look)
  const solidPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'solid', useInst: 1, noise: 0.35 }), S.n + 8, group);
  const bodyPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'membrane', opacity: 0.3, useInst: 1, gain: 0.8 }), S.n + 8, group);
  const nucPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', color: 0xb9c8ff, stain: 0x5aa0ff, gain: 0.6 }), S.n + 8, group);
  bodyPool.mesh.renderOrder = 3;
  const applyMode = (mode) => {
    const phys = mode === 'physical';
    solidPool.mesh.visible = phys;
    bodyPool.mesh.visible = nucPool.mesh.visible = !phys;
  };
  group.add(groundGrid(300, -BOUND - 8, 20));
  const white = new THREE.Color(0xffffff);

  function build() {
    rng = new RNG(4242);
    tissue = new Tissue({
      noise: S.noise, growth: false, muRep: 5, cut: 1.4, rng, hashSize: 14,
      adhesion: (a, b) => J[a][b], bounds: { r: BOUND, k: 6 },
    });
    const Rag = R * Math.cbrt(S.n / 0.6);
    for (let i = 0; i < S.n; i++) {
      let p;
      do { p = [rng.uniform(-1, 1) * Rag, rng.uniform(-1, 1) * Rag, rng.uniform(-1, 1) * Rag]; } while (Math.hypot(...p) > Rag);
      tissue.add(new Cell({ p, R, type: i % P.types }));
    }
    tissue.step(1e-4); // populate the contact lists
    lastSample = -1;
    hetChart?.clear(); idxChart?.clear();
  }

  function stats() {
    const byId = new Map(tissue.cells.map((c) => [c.id, c]));
    let het = 0, tot = 0;
    const cen = tissue.centroid();
    const rs = [0, 0, 0], ns = [0, 0, 0];
    for (const c of tissue.cells) {
      for (const id of tissue.neighbours.get(c.id) ?? []) { tot++; if (byId.get(id)?.type !== c.type) het++; }
      rs[c.type] += Math.hypot(c.p[0] - cen[0], c.p[1] - cen[1], c.p[2] - cen[2]);
      ns[c.type]++;
    }
    const f = ns.map((n) => n / tissue.cells.length);
    let random = 1;
    for (const x of f) random -= x * x;
    return { het: tot ? het / tot : 0, random, meanR: rs.map((r, k) => (ns[k] ? r / ns[k] : NaN)), contacts: tot / tissue.cells.length };
  }

  /** Steinberg / DAH outcome for two types from the adhesion matrix (interfacial-tension argument). */
  function prediction() {
    if (P.types === 3) return 'layered by cohesion: A inside, C outside';
    // works of adhesion W ~ J: mix if W_AB > (W_AA + W_BB)/2; the more cohesive type is
    // completely engulfed if min(W_AA, W_BB) <= W_AB; partial engulfment below; separation as W_AB -> 0
    const [a, b, ab] = [J[0][0], J[1][1], J[0][1]];
    const inner = a >= b ? 'A' : 'B', outer = a >= b ? 'B' : 'A';
    if (ab > (a + b) / 2) return ab > Math.max(a, b) ? 'intermixing (checkerboard)' : 'intermixing';
    if (ab >= Math.min(a, b)) return `engulfment: ${inner} core, ${outer} shell`;
    if (ab < 0.15 * Math.min(a, b)) return 'separation into two aggregates';
    return `partial engulfment of ${inner} by ${outer}`;
  }

  let hetChart, idxChart, ro;
  function buildPanel() {
    panel.section('Differential adhesion (Steinberg)');
    panel.select({
      label: 'Adhesion matrix preset', value: S.preset,
      options: Object.entries(PRESETS).map(([value, p]) => ({ value, label: p.label })),
      onChange: (v) => { S.preset = v; P = PRESETS[v]; J = P.J.map((r) => r.slice()); build(); env.rebuildPanel(); },
    });
    if (P.types === 2) {
      const setJ = (i, j, v) => { J[i][j] = J[j][i] = v; ro?.set('DAH prediction', prediction()); };
      panel.slider({ label: 'J_AA (A–A adhesion)', min: 0, max: 12, value: J[0][0], onChange: (v) => setJ(0, 0, v) });
      panel.slider({ label: 'J_BB (B–B adhesion)', min: 0, max: 12, value: J[1][1], onChange: (v) => setJ(1, 1, v) });
      panel.slider({ label: 'J_AB (heterotypic)', min: 0, max: 12, value: J[0][1], onChange: (v) => setJ(0, 1, v) });
    } else {
      panel.table(['J', 'A', 'B', 'C'], J.map((r, i) => ['ABC'[i], ...r]));
    }
    panel.slider({ label: 'Motility (effective temperature)', min: 1, max: 16, value: S.noise, onChange: (v) => { S.noise = v; tissue.noise = v; } });
    panel.toggle({ label: 'Cut-away (show aggregate interior)', value: S.cutaway, onChange: (v) => { S.cutaway = v; draw(); } });
    panel.buttons([{ label: 'Re-mix cells', primary: true, onClick: () => build() }]);
    panel.equation('η dxᵢ/dt = Σⱼ F(dᵢⱼ) + √(2ηT) ξ\nF = μ s ln(1 + (d − s)/s)  (d < s)\nF = J(τᵢ,τⱼ) (d − s) e^(−5(d − s)/s)  (adhesion)');
    panel.note('Differential adhesion hypothesis: Steinberg 1963 (Science 141:401), measured with aggregate surface tensions by Foty & Steinberg 2005 (Dev Biol 278:255). Mechanics: the centre-based model of <code>agents3d.js</code> with growth off; cells sit in a soft spherical low-adhesion well. J values are illustrative.');
    panel.section('Sorting');
    hetChart = panel.chart({
      title: 'Heterotypic contact fraction', xLabel: 'time (a.u.)', yRange: [0, 0.8],
      series: [{ name: 'aggregate', color: '#ff6fae' }, { name: 'random mixing', color: '#8a98ab', dash: [4, 3] }],
    });
    idxChart = panel.chart({
      title: 'Mean distance from the centre, by type', xLabel: 'time (a.u.)',
      series: [0, 1, 2].slice(0, P.types).map((k) => ({ name: 'ABC'[k], color: `#${COLORS[k].getHexString()}` })),
    });
    ro = panel.readouts(['cells', 'contacts per cell', 'heterotypic fraction', 'DAH prediction']);
    ro.set('cells', S.n);
    ro.set('DAH prediction', prediction());
  }

  // confocal = one optical section (|z| < 7 um) through the aggregate, as a microscope would see it
  const shown = (c) => (stage.mode === 'confocal' ? Math.abs(c.p[2]) < 7 : !(S.cutaway && c.p[2] > 2));
  function draw() {
    solidPool.begin(); bodyPool.begin(); nucPool.begin();
    for (const c of tissue.cells) {
      if (!shown(c)) continue;
      solidPool.put(c.p[0], c.p[1], c.p[2], c.R * 1.02, c.R * 1.02, c.R * 1.02, COLORS[c.type]);
      bodyPool.put(c.p[0], c.p[1], c.p[2], c.R * 1.08, c.R * 1.08, c.R * 1.08, COLORS[c.type]);
      nucPool.put(c.p[0], c.p[1], c.p[2], c.R * 0.5, c.R * 0.46, c.R * 0.5, white);
    }
    solidPool.end(); bodyPool.end(); nucPool.end();
  }

  build();
  draw();
  applyMode(stage.mode);
  const visibleCells = () => tissue.cells.filter(shown);

  return {
    frame() {
      stage.frame([0, 0, 0], 62, new THREE.Vector3(0.25, 0.35, 1));
      stage.clipAxis.set(0, 0, 1);
      stage.clipRange = [-55, 55];
    },
    buildPanel,
    legend: () => [
      ...[0, 1, 2].slice(0, P.types).map((k) => ({ color: `#${COLORS[k].getHexString()}`, label: `type ${TYPE_NAMES[k]}` })),
      { color: '#b9c8ff', label: 'nuclei' },
    ],
    pickables: [solidPool, bodyPool].map((pool) => ({
      get object() { return pool.mesh; },
      info: (hit) => {
        if (!pool.mesh.visible) return null;
        const c = visibleCells()[hit.instanceId];
        if (!c) return null;
        const ids = tissue.neighbours.get(c.id) ?? [];
        const byId = new Map(tissue.cells.map((x) => [x.id, x]));
        const het = ids.filter((id) => byId.get(id)?.type !== c.type).length;
        return `Cell ${c.id} · type ${TYPE_NAMES[c.type]}\n${ids.length} contacts, ${het} heterotypic`;
      },
    })),
    onMode(mode) { applyMode(mode); },
    update(dt) {
      if (dt > 0) {
        const t0 = performance.now();
        const target = tissue.t + dt * 12;
        while (tissue.t < target && performance.now() - t0 < 8) tissue.step(0.04);
        if (tissue.t - lastSample > 0.4 || lastSample < 0) {
          lastSample = tissue.t;
          const s = stats();
          hetChart.push(tissue.t, [s.het, s.random]);
          idxChart.push(tissue.t, s.meanR.slice(0, P.types));
          ro.set('contacts per cell', s.contacts);
          ro.set('heterotypic fraction', s.het);
        }
      }
      draw();
      env.setStatus(`t = ${fmt(tissue.t)} (a.u.) · ${tissue.cells.length} cells · ${prediction()}`);
    },
    reset() { build(); },
  };
}

// ======================================================================
// (d) ECM deposition and chemotaxis + haptotaxis
// ======================================================================
function ecmView(env) {
  const { stage, panel, group } = env;
  const SIZE = [420, 110, 300];
  const S = { lambda: 160, floor: 'chem', trails: true, hoursPerSecond: 1.2 };
  let sim, lastSample = -1;
  // collagen imaged by second-harmonic generation in the confocal look (cyan-white fibres)
  const fibrePool = new InstancePool(xCylinder(6), cellMaterial({ role: 'solid', stain: 0x9fdcff, noise: 0.35, gain: 2.6 }), 3000, group);
  const fibroPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'cytoplasm', opacity: 0.55, stain: 0xff9f43, gain: 0.7 }), 64, group);
  const migrPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'membrane', opacity: 0.35, useInst: 1, gain: 0.8 }), 96, group);
  const nucPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', color: 0xb9c8ff, gain: 0.6 }), 160, group);
  fibroPool.mesh.renderOrder = 2; migrPool.mesh.renderOrder = 3;
  group.add(boxLines(SIZE, 0x3a4a61, 0.45));
  const source = new THREE.Mesh(unitSphere(3), cellMaterial({ role: 'reporter', color: 0xff5ab4, emissive: 0.7, gain: 0.8 }));
  group.add(source);
  // iso-concentration shells, cut to the tissue slab by two local clipping planes
  const slab = [new THREE.Plane(new THREE.Vector3(0, -1, 0), SIZE[1] / 2), new THREE.Plane(new THREE.Vector3(0, 1, 0), SIZE[1] / 2)];
  const shells = [0.35, 0.12].map((c, i) => {
    const mat = cellMaterial({ role: 'membrane', color: [0xff5ab4, 0xb45bff][i], stain: [0xff5ab4, 0xb45bff][i], opacity: 0.03, rimPower: 4, side: THREE.DoubleSide, gain: 0.25 });
    mat.clippingPlanes = slab;
    const m = new THREE.Mesh(unitSphere(4), mat);
    m.userData.level = c;
    m.renderOrder = 4;
    group.add(m);
    return m;
  });
  // floor heat map
  const FN = 64;
  const floorGeo = new THREE.PlaneGeometry(SIZE[0], SIZE[2], FN - 1, Math.round(FN * SIZE[2] / SIZE[0]) - 1);
  floorGeo.rotateX(-Math.PI / 2);
  const floorCol = new THREE.BufferAttribute(new Float32Array(floorGeo.attributes.position.count * 3), 3);
  floorGeo.setAttribute('color', floorCol);
  const floor = new THREE.Mesh(floorGeo, cellMaterial({ role: 'field', vertexColors: true, noise: 0.0, gain: 0.45, he: [0.02, 0.2] }));
  floor.position.y = -SIZE[1] / 2 - 0.5;
  group.add(floor);
  // trails
  const TRAIL = 60;
  const trailGeo = new THREE.BufferGeometry();
  const trailPos = new Float32Array(96 * TRAIL * 2 * 3), trailCol = new Float32Array(96 * TRAIL * 2 * 3);
  trailGeo.setAttribute('position', new THREE.BufferAttribute(trailPos, 3).setUsage(THREE.DynamicDrawUsage));
  trailGeo.setAttribute('color', new THREE.BufferAttribute(trailCol, 3).setUsage(THREE.DynamicDrawUsage));
  const trails = new THREE.LineSegments(trailGeo, lineMaterial({ vertexColors: true, opacity: 0.55 }));
  trails.frustumCulled = false;
  group.add(trails);

  const chemRamp = makeRamp([[0, '#07070f'], [0.35, '#2a1846'], [0.7, '#8a2f7a'], [1, '#ff6fb8']]);
  const ecmRamp = makeRamp([[0, '#06080c'], [0.4, '#3d2b1e'], [0.75, '#b5774a'], [1, '#ffd9b0']]);
  const tmpC = new THREE.Color();
  const cMigr = new THREE.Color('#4aa3ff'), cArr = new THREE.Color('#39ff88'), cFib = new THREE.Color('#ff9f43');
  const cFibre = new THREE.Color('#f6d7c3'), cFibreNew = new THREE.Color('#ffffff');

  function build() {
    const k = 100 / (S.lambda * S.lambda);
    const keep = sim?.params;
    sim = new MigrationSim({ size: SIZE, rng: new RNG(1234), D: 100, k, params: keep ? { ...keep } : {} });
    source.position.set(...sim.source.p);
    source.scale.setScalar(sim.source.radius);
    placeShells();
    lastSample = -1;
    distChart?.clear(); ecmChart?.clear();
  }
  function placeShells() {
    // radius where the normalised steady field equals each level (bisection)
    for (const m of shells) {
      let a = sim.source.radius, b = 600;
      for (let i = 0; i < 50; i++) { const r = (a + b) / 2; if (sim.source.value([sim.source.p[0] + r, sim.source.p[1], sim.source.p[2]]) > m.userData.level) a = r; else b = r; }
      m.scale.setScalar((a + b) / 2);
      m.position.set(...sim.source.p);
    }
  }

  function drawFloor() {
    const pos = floorGeo.attributes.position;
    const arr = floorCol.array;
    const p = [0, 0, 0];
    const E = sim.ecm;
    for (let v = 0; v < pos.count; v++) {
      p[0] = pos.getX(v); p[1] = 0; p[2] = pos.getZ(v);
      if (S.floor === 'chem') {
        const c = sim.source.value(p);
        chemRamp(clamp01(1 + Math.log10(Math.max(c, 1e-4)) / 3), arr, 3 * v);
      } else {
        let s = 0;
        for (let y = -SIZE[1] / 2 + 5; y < SIZE[1] / 2; y += 10) { p[1] = y; s += E.density(p); }
        ecmRamp(clamp01(s / (SIZE[1] / 10) / 0.8), arr, 3 * v);
      }
    }
    floorCol.needsUpdate = true;
  }

  function draw() {
    fibrePool.begin();
    for (const f of sim.fibres) {
      const r = 0.55 + 0.9 * Math.sqrt(Math.max(f.mass, 0) / sim.params.quantum);
      const age = sim.t - f.born;
      const col = f.by < 0 ? cFibre : tmpC.copy(cFibreNew).lerp(cFibre, clamp01(age / 6));
      fibrePool.putOriented(f.p, f.dir, f.len, r, r, col);
    }
    fibrePool.end();
    fibroPool.begin(); migrPool.begin(); nucPool.begin();
    for (const f of sim.fibro) {
      const sp = Math.hypot(...f.v) || 1;
      const ax = [f.v[0] / sp, f.v[1] / sp, f.v[2] / sp];
      fibroPool.putOriented(f.p, ax, 15, 4.2, 5, cFib);
      nucPool.putOriented(f.p, ax, 5, 2.4, 2.8, null);
    }
    for (const c of sim.migr) {
      const sp = Math.hypot(...c.v);
      const ax = sp > 1e-6 ? [c.v[0] / sp, c.v[1] / sp, c.v[2] / sp] : [1, 0, 0];
      const stretch = 1 + Math.min(sp / 60, 0.6);
      migrPool.putOriented(c.p, ax, 6.5 * stretch, 6.5 / Math.sqrt(stretch), 6.5 / Math.sqrt(stretch), c.arrived ? cArr : cMigr);
      nucPool.putOriented(c.p, ax, 3.2, 2.8, 2.8, null);
    }
    fibroPool.end(); migrPool.end(); nucPool.end();
    // trails
    let s = 0;
    if (S.trails) {
      for (const c of sim.migr) {
        const tr = c.trail;
        for (let i = 1; i < tr.length; i++) {
          const a = tr[i - 1], b = tr[i];
          const f = i / tr.length;
          for (const [q, w] of [[a, f - 1 / tr.length], [b, f]]) {
            trailPos[s * 3] = q[0]; trailPos[s * 3 + 1] = q[1]; trailPos[s * 3 + 2] = q[2];
            trailCol[s * 3] = 0.05 + 0.15 * w; trailCol[s * 3 + 1] = 0.15 + 0.4 * w; trailCol[s * 3 + 2] = 0.3 + 0.6 * w;
            s++;
          }
        }
      }
    }
    trailGeo.setDrawRange(0, s);
    trailGeo.attributes.position.needsUpdate = true;
    trailGeo.attributes.color.needsUpdate = true;
  }

  let distChart, ecmChart, ro;
  function buildPanel() {
    const P = sim.params;
    panel.section('Cell movement = chemotaxis + haptotaxis');
    panel.equation('v = χ/(1 + αC) ∇C  +  ρ ∇E  +  √(2Dn) ξ\n(Anderson & Chaplain 1998)');
    panel.slider({ label: 'χ chemotactic sensitivity (µm²/h)', min: 0, max: 2e6, step: 1e4, value: P.chi, format: (v) => v.toExponential(1), onChange: (v) => { sim.params.chi = v; } });
    panel.slider({ label: 'ρ haptotactic coefficient (µm²/h)', min: 0, max: 3000, step: 10, value: P.rho, onChange: (v) => { sim.params.rho = v; } });
    panel.slider({ label: 'Random motility Dn (µm²/h)', min: 0, max: 400, value: P.Dn, onChange: (v) => { sim.params.Dn = v; } });
    panel.slider({ label: 'Chemokine decay length λ (µm)', min: 40, max: 400, value: S.lambda, onChange: (v) => { S.lambda = v; rebuildSource(); } });
    panel.slider({ label: 'MMP matrix degradation (1/h)', min: 0, max: 3, value: P.mmp, onChange: (v) => { sim.params.mmp = v; } });
    panel.toggle({ label: 'Fibroblasts deposit collagen', value: P.deposit, onChange: (v) => { sim.params.deposit = v; } });
    panel.toggle({ label: 'Show migration tracks', value: S.trails, onChange: (v) => { S.trails = v; } });
    panel.select({
      label: 'Floor shows', value: S.floor,
      options: [{ value: 'chem', label: 'Chemoattractant C (log scale)' }, { value: 'ecm', label: 'ECM fibril density E (column mean)' }],
      onChange: (v) => { S.floor = v; drawFloor(); env.legend(legendItems()); },
    });
    panel.buttons([{ label: 'Restart', primary: true, onClick: () => { build(); drawFloor(); } }]);

    panel.section('Collagen pathway (per fibroblast)');
    panel.equation('Procollagen →(ADAMTS2, BMP-1) Tropocollagen + propeptides\n2 Tropocollagen → Fibril (nucleation)\nTropocollagen + Fibril → Fibril (elongation)');
    panel.equation('dP/dt = σ − k_c P\ndT/dt = k_c P − 2k_n T² − k_g T F\ndF/dt = 2k_n T² + k_g T F');
    distChart = panel.chart({
      title: 'Mean distance of migrating cells to the source', xLabel: 'h',
      series: [{ name: 'distance (µm)', color: '#4aa3ff' }],
    });
    ecmChart = panel.chart({
      title: 'Collagen fibril mass in the matrix', xLabel: 'h',
      series: [{ name: 'in ECM', color: '#ffd9b0' }, { name: 'deposited', color: '#ff9f43', dash: [4, 3] }, { name: 'degraded (MMP)', color: '#ff6fae', dash: [2, 3] }],
    });
    ro = panel.readouts(['fibres', 'mean ECM density E', 'propeptides released', 'cells at source', 'mean speed']);
  }

  function rebuildSource() {
    // new decay length, source still normalised to C = 1 at its surface
    const s = sim.source;
    sim.source = new PointSourceField({ p: s.p, radius: s.radius, D: s.D, k: s.D / (S.lambda * S.lambda) });
    placeShells();
    drawFloor();
  }

  function legendItems() {
    return [
      { color: '#ff5ab4', label: 'chemoattractant source (e.g. CXCL12 niche) + iso-concentration shells' },
      { color: '#4aa3ff', label: 'migrating cells (chemotaxis + haptotaxis)' },
      { color: '#39ff88', label: 'arrived at the source' },
      { color: '#ff9f43', label: 'fibroblasts (persistent random walk)' },
      { color: '#f6d7c3', label: 'collagen fibres (white = freshly deposited)' },
      S.floor === 'chem' ? { color: '#8a2f7a', label: 'floor: chemoattractant C (log)' } : { color: '#b5774a', label: 'floor: ECM density E' },
    ];
  }

  build();
  drawFloor();
  draw();

  return {
    frame() {
      stage.frame([0, -10, 0], 205, new THREE.Vector3(0.12, 0.85, 0.75));
      stage.clipAxis.set(0, 1, 0);
      stage.clipRange = [-SIZE[1] / 2, SIZE[1] / 2];
    },
    buildPanel,
    legend: legendItems,
    pickables: [
      {
        get object() { return migrPool.mesh; },
        info: (hit) => {
          const c = sim.migr[hit.instanceId];
          if (!c) return null;
          const d = Math.hypot(c.p[0] - sim.source.p[0], c.p[1] - sim.source.p[1], c.p[2] - sim.source.p[2]) - sim.source.radius;
          return `Migrating cell ${c.id}${c.arrived ? ' · arrived' : ''}\n${fmt(d)} µm from the source · C = ${fmt(sim.source.value(c.p))}\nlocal ECM density E = ${fmt(sim.ecm.density(c.p))}`;
        },
      },
      {
        get object() { return fibroPool.mesh; },
        info: (hit) => {
          const f = sim.fibro[hit.instanceId];
          if (!f) return null;
          return `Fibroblast ${f.id} · ${f.fibres} fibres laid\nprocollagen ${fmt(f.y[0])} · tropocollagen ${fmt(f.y[1])} · fibril ${fmt(f.y[2])}`;
        },
      },
      {
        get object() { return fibrePool.mesh; },
        info: (hit) => {
          const f = sim.fibres[hit.instanceId];
          if (!f) return null;
          return `Collagen fibre · mass ${fmt(f.mass)}\n${f.by < 0 ? 'pre-existing interstitial matrix' : `laid by fibroblast ${f.by} at t = ${fmt(f.born)} h`}`;
        },
      },
      { object: source, info: () => 'Chemoattractant source: steady C ∝ e^(−r/λ)/r (point source with degradation)' },
    ],
    update(dt) {
      if (dt > 0) {
        const T = dt * S.hoursPerSecond;
        const n = Math.max(1, Math.ceil(T / 0.05));
        for (let i = 0; i < n; i++) sim.step(T / n);
        sim.recordTrails(TRAIL, 4);
        if (sim.t - lastSample > 0.25 || lastSample < 0) {
          lastSample = sim.t;
          const mass = sim.fibreMass();
          distChart.push(sim.t, [sim.meanDistance()]);
          ecmChart.push(sim.t, [mass, sim.deposited, sim.degraded]);
          ro.set('fibres', sim.fibres.length);
          ro.set('mean ECM density E', sim.ecm.mean());
          ro.set('propeptides released', sim.cleaved);
          ro.set('cells at source', `${sim.migr.filter((c) => c.arrived).length} / ${sim.migr.length}`);
          const v = [0, 0, 0];
          let spd = 0;
          for (const c of sim.migr) if (!c.arrived) spd += Math.hypot(...sim.driftVelocity(c.p, v));
          ro.set('mean speed', `${fmt(spd / Math.max(1, sim.migr.filter((c) => !c.arrived).length))} µm/h`);
          if (S.floor === 'ecm') drawFloor();
        }
      }
      draw();
      env.setStatus(`t = ${fmt(sim.t)} h · ${sim.fibres.length} collagen fibres · ${sim.migr.filter((c) => c.arrived).length}/${sim.migr.length} cells at source`);
    },
    reset() { build(); drawFloor(); },
  };
}

// ======================================================================
// (e) Crypt fate logic: Wnt x Notch crosstalk
// ======================================================================
function cryptView(env) {
  const { stage, panel, group } = env;
  const lat = cryptLattice({ perRing: 18, rings: 16, radius: 22 });
  const n = lat.pos.length;
  const S = { lamW: 30, W0: 1.5, colorBy: 'fate', hoursPerSecond: 2.5, cutaway: true };
  const params = { ...crosstalk.defaults };
  let model, lastSample = -1;
  const fateCol = FATES.map((f) => new THREE.Color(f.color));
  const bodyPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'cytoplasm', opacity: 0.42, useInst: 1, gain: 0.3, he: [0.02, 0.3], side: THREE.FrontSide }), n + 8, group);
  const nucPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', color: 0xc9d6ff, gain: 0.22, he: [1.5, 0.2] }), n + 8, group);
  // Paneth granules are strongly eosinophilic; goblet mucin stays pale in H&E
  const granPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter', he: [0.02, 1.3], gain: 0.8 }), n * 2, group);
  const mucPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter', he: [0.06, 0.06], gain: 0.6 }), n, group);
  bodyPool.mesh.renderOrder = 2;
  // stromal Wnt niche around the crypt base: vertex-coloured translucent sheath
  const sheathGeo = new THREE.CylinderGeometry(33, 33, lat.length - lat.radius + 8, 48, 24, true);
  sheathGeo.translate(0, lat.radius + (lat.length - lat.radius) / 2 - 4, 0);
  const capGeo = new THREE.SphereGeometry(33, 48, 12, 0, Math.PI * 2, Math.PI / 2, Math.PI / 2);
  capGeo.translate(0, lat.radius, 0);
  const sheaths = [sheathGeo, capGeo].map((g) => {
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(g.attributes.position.count * 3), 3));
    const m = new THREE.Mesh(g, cellMaterial({ role: 'membrane', vertexColors: true, opacity: 0.1, useInst: 1, gain: 0.35 }));
    m.renderOrder = 4;
    group.add(m);
    return m;
  });
  group.add(groundGrid(360, -14, 24));
  const wntRamp = makeRamp([[0, '#0c1320'], [0.4, '#6b3fa0'], [1, '#ffd166']]);
  const valRamp = makeRamp([[0, '#101a2c'], [0.5, '#3f7fd6'], [1, '#e8f4ff']]);
  const tmpC = new THREE.Color();
  const rgb = [0, 0, 0];

  const wntAt = (s) => S.W0 * Math.exp(-s / S.lamW);
  const drawn = [];
  const white = new THREE.Color(0xffffff);
  function paintSheath() {
    for (const m of sheaths) {
      const pos = m.geometry.attributes.position, col = m.geometry.attributes.color;
      for (let v = 0; v < pos.count; v++) {
        const y = pos.getY(v);
        // arc length from the base pole of the crypt (cap below y = radius)
        const s = y < lat.radius ? lat.radius * Math.acos(clamp01((lat.radius - y) / lat.radius)) : lat.radius * Math.PI / 2 + (y - lat.radius);
        wntRamp(clamp01(wntAt(s) / 1.5), col.array, 3 * v);
      }
      col.needsUpdate = true;
    }
  }

  function build() {
    const wnt = lat.arc.map(wntAt);
    model = new CryptModel({ nbrs: lat.nbrs, wnt, rng: new RNG(5), params });
    lastSample = -1;
    fracChart?.clear();
    paintSheath();
  }
  function setWnt() {
    lat.arc.forEach((s, i) => { model.wnt[i] = wntAt(s); });
    paintSheath();
  }

  function draw() {
    bodyPool.begin(); nucPool.begin(); granPool.begin(); mucPool.begin();
    drawn.length = 0;
    for (let i = 0; i < n; i++) {
      if (S.cutaway && lat.pos[i][2] > 3) continue;
      drawn.push(i);
      const s = model.S[i];
      const f = model.fate(i);
      let col;
      // H&E is untinted: fates are then recognised by morphology (Paneth granules, goblet mucin)
      if (stage.mode === 'histology') col = white;
      else if (S.colorBy === 'fate') col = fateCol[f];
      else {
        const v = S.colorBy === 'bcat' ? s[1] / 1.2 : S.colorBy === 'nicd' ? s[2] / 3 : s[4] / 2;
        col = tmpC.fromArray(valRamp(clamp01(v), rgb));
      }
      const p = lat.pos[i], nrm = lat.normal[i];
      // columnar cell: long axis along the surface normal, apical side facing the lumen
      bodyPool.putOriented(p, nrm, 7.5, 3.9, 3.9, col);
      const nucOut = f === 1 ? 3.2 : 2.2; // Paneth nuclei sit basally, granules apically
      nucPool.putOriented([p[0] + nrm[0] * nucOut, p[1] + nrm[1] * nucOut, p[2] + nrm[2] * nucOut], nrm, 2.6, 1.9, 1.9, null);
      if (f === 1) {
        // Paneth cell: apical lysozyme / defensin granules
        for (const k of [-1, 1]) {
          granPool.put(p[0] - nrm[0] * 3.2 + k * nrm[2] * 1.4, p[1] - nrm[1] * 3.2 + k * 1.1, p[2] - nrm[2] * 3.2 - k * nrm[0] * 1.4, 1.25, 1.25, 1.25, tmpC.set('#ffd6ec'));
        }
      } else if (f === 3) {
        // goblet cell: apical mucus theca
        mucPool.put(p[0] - nrm[0] * 3.2, p[1] - nrm[1] * 3.2, p[2] - nrm[2] * 3.2, 2.4, 2.4, 2.4, tmpC.set('#fff3c4'));
      }
    }
    bodyPool.end(); nucPool.end(); granPool.end(); mucPool.end();
  }

  let fracChart, posChart, ro;
  function buildPanel() {
    panel.section('Wnt niche & Notch lateral inhibition');
    panel.equation('W(s) = W₀ e^(−s/λ_W)   (s = distance from the crypt base)\nDvl* ⊣ GSK-3β ⊣ β-catenin;  Dll1ⱼ → NICDᵢ → Hes1 ⊣ Atoh1 → Dll1ᵢ');
    panel.slider({ label: 'Wnt decay length λ_W (µm)', min: 8, max: 90, value: S.lamW, onChange: (v) => { S.lamW = v; setWnt(); } });
    panel.slider({ label: 'Wnt at the base W₀', min: 0, max: 3, value: S.W0, onChange: (v) => { S.W0 = v; setWnt(); } });
    panel.section('Crosstalk & combinatorial regulation');
    panel.slider({ label: 'Dvl–NICD sequestration (Axelrod 1996)', min: 0, max: 4, value: params.xDvlNICD, onChange: (v) => { params.xDvlNICD = v; } });
    panel.slider({ label: 'Wnt → Jagged1 → Notch (Rodilla 2009)', min: 0, max: 1, value: params.xJag, onChange: (v) => { params.xJag = v; } });
    panel.slider({ label: 'Enhancer cooperativity ω', min: 0.2, max: 10, value: params.omega, log: true, onChange: (v) => { params.omega = v; } });
    panel.toggle({ label: 'γ-secretase inhibitor (DAPT; van Es 2005)', value: params.dapt > 0, onChange: (v) => { params.dapt = v ? 1 : 0; } });
    panel.toggle({ label: 'Constitutive NICD (Fre 2005)', value: params.nicdOE > 0, onChange: (v) => { params.nicdOE = v ? 2 : 0; } });
    panel.select({
      label: 'Colour cells by', value: S.colorBy,
      options: [{ value: 'fate', label: 'Fate (combinatorial read-out)' }, { value: 'bcat', label: 'Nuclear β-catenin' }, { value: 'nicd', label: 'NICD' }, { value: 'atoh', label: 'Atoh1 (secretory)' }],
      onChange: (v) => { S.colorBy = v; env.legend(legendItems()); },
    });
    panel.toggle({ label: 'Cut-away (longitudinal section through the crypt)', value: S.cutaway, onChange: (v) => { S.cutaway = v; } });
    panel.buttons([{ label: 'Restart crypt', primary: true, onClick: () => build() }]);
    panel.equation('Z = 1 + x + y + ω x y,  x = (β-cat/K_B)², y = (NICD/K_N)²\nstem ∝ ωxy · Paneth ∝ x · absorptive ∝ y · goblet/EE ∝ 1');
    fracChart = panel.chart({
      title: 'Fate fractions over time', xLabel: 'h', yRange: [0, 1],
      series: FATES.map((f, i) => ({ name: ['stem', 'Paneth', 'absorptive', 'goblet/EE'][i], color: f.color })),
    });
    posChart = panel.chart({
      title: 'Fate composition along the crypt axis (base → top)', xLabel: 'µm from base', yRange: [0, 1],
      series: FATES.map((f, i) => ({ name: ['stem', 'Paneth', 'absorptive', 'goblet/EE'][i], color: f.color })),
    });
    ro = panel.readouts(['cells', 'secretory (Atoh1+) fraction', 'secretory–secretory contacts', 'Paneth : stem at the base']);
    ro.set('cells', n);
  }

  function sample() {
    const f = model.fates();
    const fr = [0, 0, 0, 0];
    for (const k of f) fr[k]++;
    fracChart.push(model.t, fr.map((c) => c / n));
    const bins = 8, xs = [], cols = [[], [], [], []];
    for (let b = 0; b < bins; b++) {
      const lo = (b / bins) * lat.length, hi = ((b + 1) / bins) * lat.length;
      const c = [0, 0, 0, 0];
      let m = 0;
      for (let i = 0; i < n; i++) if (lat.arc[i] >= lo && lat.arc[i] < hi) { c[f[i]]++; m++; }
      xs.push((lo + hi) / 2);
      for (let k = 0; k < 4; k++) cols[k].push(m ? c[k] / m : NaN);
    }
    posChart.set(xs, cols);
    let sec = 0, ss = 0, st = 0, pa = 0;
    for (let i = 0; i < n; i++) {
      if (model.S[i][4] > 0.8) { sec++; for (const j of lat.nbrs[i]) if (model.S[j][4] > 0.8) ss++; }
      if (lat.arc[i] < lat.length * 0.25) { if (f[i] === 0) st++; if (f[i] === 1) pa++; }
    }
    ro.set('secretory (Atoh1+) fraction', sec / n);
    ro.set('secretory–secretory contacts', ss / 2);
    ro.set('Paneth : stem at the base', `${pa} : ${st}`);
  }

  function legendItems() {
    if (S.colorBy === 'fate') {
      return [...FATES.map((f) => ({ color: f.color, label: f.label })), { color: '#ffd166', label: 'Wnt (stromal / Paneth niche), sheath' }];
    }
    const name = { bcat: 'nuclear β-catenin', nicd: 'NICD', atoh: 'Atoh1' }[S.colorBy];
    return [{ color: '#101a2c', label: `${name}: low` }, { color: '#3f7fd6', label: 'intermediate' }, { color: '#e8f4ff', label: 'high' }];
  }

  function preRoll(hours) {
    // start from a patterned crypt; the charts keep the history
    while (model.t < hours) {
      model.advance(0.5, 0.05);
      lastSample = model.t;
      sample();
    }
  }

  build();
  draw();
  return {
    frame() {
      stage.frame([0, lat.length / 2 - 6, 0], 100, new THREE.Vector3(0.55, 0.55, 1));
      stage.clipAxis.set(0, 0, 1);
      stage.clipRange = [-40, 40];
    },
    buildPanel,
    legend: legendItems,
    pickables: [{
      get object() { return bodyPool.mesh; },
      info: (hit) => {
        const i = drawn[hit.instanceId];
        if (i === undefined) return null;
        const s = model.S[i];
        const pr = fateProbabilities(s[1], s[2], model.p);
        return `Crypt cell ${i} · ${FATES[model.fate(i)].label}\n${fmt(lat.arc[i])} µm from base · Wnt ${fmt(model.wnt[i])}\n` +
          `β-cat ${fmt(s[1])} · NICD ${fmt(s[2])} · Hes1 ${fmt(s[3])} · Atoh1 ${fmt(s[4])}\n` +
          `P(stem, Paneth, abs, goblet) = ${pr.map((v) => v.toFixed(2)).join(', ')}`;
      },
    }],
    update(dt) {
      if (dt > 0) {
        model.advance(dt * S.hoursPerSecond, 0.05);
        if (model.t - lastSample > 0.5 || lastSample < 0) { lastSample = model.t; sample(); }
      }
      draw();
      env.setStatus(`t = ${fmt(model.t)} h · ${n} crypt cells · λ_W = ${fmt(S.lamW)} µm`);
    },
    reset() { build(); },
    afterPanel() { if (lastSample < 0) { preRoll(30); draw(); } },
    onMode(mode) { for (const m of sheaths) m.visible = mode !== 'histology'; },
  };
}

// ======================================================================
const BUILDERS = { organoid: organoidView, gradient: gradientView, sorting: sortingView, ecm: ecmView, crypt: cryptView };

export default {
  about: `
    <p>Fourth-order behaviour: many pathways and many cells acting together to build tissue.
    Five views (selector below), all simulated live in 3D:</p>
    <p><b>Turing organoid</b> — a reaction–diffusion system solved <i>on the curved surface</i> of a
    10 242-vertex organoid (cotangent Laplace–Beltrami operator, Meyer et al. 2003). Activator peaks
    push the epithelium outward into <b>buds</b>. Linear Turing analysis (Murray, <i>Mathematical Biology II</i>)
    predicts which spherical harmonic the organoid selects; the measured angular spectrum is shown next to it.
    Schnakenberg (1979), Gierer–Meinhardt (1972) and Gray–Scott (Pearson 1993) kinetics are the published models;
    mapping activator to bud height is a <i>phenomenological</i> stand-in for epithelial mechanics.</p>
    <p><b>Morphogen gradient</b> — Shh from the floor plate and BMP from the roof plate spread by
    source–diffusion–degradation on a 3D grid through the neural tube (λ = √(D/k), Kicheva et al. 2007 values);
    thresholds give Wolpert's French flag of progenitor domains. Threshold read-out is the classic
    idealisation — real domains also use temporal adaptation and cross-repression.</p>
    <p><b>Differential adhesion</b> — Steinberg's hypothesis in a centre-based 3D cell model: type-specific
    adhesion sorts a random mixture into core and shell. Adhesion values are illustrative, not measured.</p>
    <p><b>ECM & migration</b> — fibroblasts secrete procollagen → propeptide cleavage → tropocollagen → fibrils
    (mass-action) and lay fibres along their path; other cells follow <i>chemotaxis + haptotaxis</i>
    (Anderson & Chaplain 1998). The nucleation–elongation rates are phenomenological.</p>
    <p><b>Crypt fate logic</b> — Wnt × Notch crosstalk (Dishevelled–NICD, Wnt→Jagged1) and a thermodynamic
    combinatorial enhancer turn a Wnt gradient plus lateral inhibition into stem, Paneth, absorptive and goblet
    cells. The fate logic follows Fre et al. 2005 and van Es et al. 2005; the ODE rates are dimensionless caricatures.</p>`,
  paperRef: 'Paper: "Fourth Order: Complex Pathways and Higher-Level Interactions" → integrated signalling networks (Wnt–Notch crosstalk, combinatorial gene regulation), tissue and organ formation, morphogen gradients, ECM interactions; "Morphogenesis" → Turing reaction–diffusion ∂u/∂t = Du∇²u + f(u,v), mechanical forces "Cell Movement = Chemotaxis + Haptotaxis"; collagen pathway (procollagen → tropocollagen → fibrils).',

  async create(ctx) {
    const { stage, root, panel, assets, legend, setStatus } = ctx;
    const orgMeta = await assets.loadJSON('organoids.json').catch(() => null);
    let current = null, currentId = 'organoid', group = null;

    const env = {
      stage, panel, legend, setStatus, orgMeta,
      get group() { return group; },
      rebuildPanel() { buildPanel(); },
    };

    function header() {
      panel.section('View');
      panel.select({ label: 'Morphogenesis model', value: currentId, options: VIEWS, onChange: (v) => activate(v) });
    }
    function buildPanel() {
      panel.clear();
      header();
      current.buildPanel();
    }

    async function activate(id) {
      if (current) {
        try { current.dispose?.(); } catch (e) { console.error(e); }
        releaseTree(group);
        root.remove(group);
      }
      currentId = id;
      group = new THREE.Group();
      root.add(group);
      current = BUILDERS[id](env);
      buildPanel();
      current.afterPanel?.();
      current.onMode?.(stage.mode);
      current.frame();
      if (stage.clipOn) stage.setClip(true, 0);
      stage.setPickables(current.pickables ?? []);
      legend(current.legend());
      stage.simTime = 0;
    }

    const requested = ctx.query?.get('view');
    await activate(VIEWS.some((v) => v.value === requested) ? requested : 'organoid');

    return {
      update(dt, t) { current?.update(dt, t); },
      reset() { current?.reset?.(); },
      onMode(mode) { current?.onMode?.(mode); legend(current.legend()); },
      dispose() { current?.dispose?.(); },
    };
  },
};
