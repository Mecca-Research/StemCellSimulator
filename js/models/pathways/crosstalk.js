// Wnt-Notch crosstalk and combinatorial gene regulation (paper: fourth order,
// "Integrated signalling networks - crosstalk between pathways like Wnt and
// Notch, pathway convergence, combinatorial gene regulation involving
// multiple transcription factors").
//
// Biology encoded (intestinal crypt):
//  - Wnt maintains crypt stem cells and is required for Paneth-cell maturation
//    (van Es et al. 2005, Nat Cell Biol 7:381); Notch keeps progenitors
//    proliferative and absorptive - Notch gain removes secretory cells
//    (Fre et al. 2005, Nature 435:964), Notch/gamma-secretase loss converts
//    crypt cells to goblet cells (van Es et al. 2005, Nature 435:959).
//  - Notch -> Hes1 represses Atoh1 (Math1), the secretory master regulator
//    (Yang et al. 2001, Science 294:2155); Atoh1+ secretory cells present
//    Dll1/Dll4 to their neighbours = lateral inhibition (Pellegrinet et al.
//    2011, Gastroenterology 140:1230), so Paneth cells alternate with Lgr5
//    stem cells at the crypt base (Sato et al. 2011, Nature 469:415).
//  - Crosstalk: Dishevelled binds NICD and dampens Notch signalling (Axelrod
//    et al. 1996, Science 271:1826); beta-catenin/TCF induces Jagged1
//    (Rodilla et al. 2009, PNAS 106:6315); GSK-3beta, inhibited by Wnt,
//    phosphorylates NICD (Foltz et al. 2002, Curr Biol 12:1006 - reported
//    effects on NICD stability differ between studies, so that term is off by
//    default).
//
// Combinatorial regulation: a thermodynamic (statistical-mechanics) promoter
// model (Shea & Ackers 1985, J Mol Biol 181:211; Bintu et al. 2005, Curr Opin
// Genet Dev 15:116) of an enhancer with a TCF/LEF site (beta-catenin) and an
// RBPJ/CSL site (NICD) with cooperativity omega:
//     Z = 1 + x + y + omega x y,   x = (B / K_B)^nB,  y = (N / K_N)^nN
// The four occupancy states map one-to-one onto the crypt fate logic
//     both bound    -> stem / progenitor maintenance     (Wnt high, Notch high)
//     TCF only      -> Paneth / secretory (Wnt-dependent) (Wnt high, Notch low)
//     RBPJ only     -> absorptive enterocyte             (Wnt low,  Notch high)
//     neither       -> differentiated goblet / EE        (Wnt low,  Notch low)
// That mapping is a phenomenological encoding of the experimental fate logic
// above; the per-cell ODEs are Hill-type caricatures with dimensionless rates
// (time unit ~ 1 h), not fitted to data.
import { hill } from '../core/ode.js';

export const FATES = [
  { id: 'stem', label: 'Stem / progenitor (Wnt+ Notch+; Lgr5, Olfm4)', color: '#39ff88' },
  { id: 'paneth', label: 'Paneth / secretory (Wnt+ Notch-; Atoh1, Sox9)', color: '#ff4fa3' },
  { id: 'absorptive', label: 'Absorptive enterocyte (Wnt- Notch+; Hes1)', color: '#4aa3ff' },
  { id: 'differentiated', label: 'Goblet / enteroendocrine (Wnt- Notch-)', color: '#ffc94a' },
];

/** Occupancy probabilities of the four enhancer states. */
export function enhancerStates(x, y, omega = 1) {
  const Z = 1 + x + y + omega * x * y;
  return { none: 1 / Z, wnt: x / Z, notch: y / Z, both: (omega * x * y) / Z };
}

/** Fate propensities [stem, paneth, absorptive, differentiated] from nuclear beta-catenin B and NICD N. */
export function fateProbabilities(B, N, p = crosstalk.defaults, out = [0, 0, 0, 0]) {
  const x = (Math.max(B, 0) / p.KB) ** p.nB, y = (Math.max(N, 0) / p.KN) ** p.nN;
  const s = enhancerStates(x, y, p.omega);
  out[0] = s.both; out[1] = s.wnt; out[2] = s.notch; out[3] = s.none;
  return out;
}

export function fateIndex(B, N, p = crosstalk.defaults) {
  const pr = fateProbabilities(B, N, p);
  let k = 0;
  for (let i = 1; i < 4; i++) if (pr[i] > pr[k]) k = i;
  return k;
}

export const crosstalk = {
  species: ['Dvl*', 'beta-catenin', 'NICD', 'Hes1', 'Atoh1', 'Dll1', 'Jagged1'],
  defaults: {
    // Wnt -> Dvl -| GSK-3beta -| beta-catenin
    kDvl: 2, KW: 0.4, dDvl: 1,
    sB: 1, kB0: 0.25, kBG: 4, KG: 0.25,
    // Notch: ligand from neighbours -> NICD -> Hes1 -| Atoh1 -> Dll1. The
    // steep operating points (small KA) give the loop gain > 2 needed for
    // lateral inhibition on a hexagonal sheet (cf. Collier et al. 1996).
    // turnover ~2 / h: Hes1 protein half-life ~22 min (Hirata et al. 2002, Science 298:840)
    kN: 6, KL: 0.5, dN: 2,
    kH: 4, KH: 1, dH: 2,
    kA: 4, KA: 0.05, dA: 2,
    kDl: 4, KDl: 1, dDl: 2,
    // crosstalk strengths
    xDvlNICD: 0.6,   // Dvl sequesters NICD (Axelrod 1996): NICD production / (1 + x Dvl)
    xJag: 0.2,       // weight of Wnt-induced Jagged1 in the ligand seen by neighbours (Rodilla 2009)
    xGSK: 0,         // optional GSK-3beta effect on NICD turnover (sign disputed; 0 = off)
    kJ: 1, KJ: 0.8, dJ: 1,
    // perturbations
    dapt: 0,         // fraction of gamma-secretase inhibited (DAPT; van Es et al. 2005 Nature)
    nicdOE: 0,       // constitutive NICD production (Villin-NICD gain of function; Fre et al. 2005)
    // enhancer read-out
    KB: 0.6, KN: 0.5, nB: 2, nN: 2, omega: 2,
  },

  /** In-place Euler update of one cell; s = Float64Array(7), W = Wnt input, L = mean ligand from neighbours. */
  stepCell(s, W, L, dt, p = crosstalk.defaults) {
    const [Dv, B, N, H, A, Dl, J] = s;
    const G = 1 / (1 + Dv / p.KG); // active GSK-3beta fraction
    const dDv = p.kDvl * hill(W, p.KW, 2) * (1 - Dv) - p.dDvl * Dv;
    const dB = p.sB - (p.kB0 + p.kBG * G) * B;
    const dN = (p.kN * (1 - p.dapt) * hill(L, p.KL, 2)) / (1 + p.xDvlNICD * Dv) + p.nicdOE - p.dN * (1 + p.xGSK * G) * N;
    const dH = p.kH * hill(N, p.KH, 2) - p.dH * H;
    const dA = p.kA / (1 + (H / p.KA) ** 2) - p.dA * A;
    const dDl = p.kDl * hill(A, p.KDl, 2) - p.dDl * Dl;
    const dJ = p.kJ * hill(B, p.KJ, 2) - p.dJ * J;
    s[0] = Dv + dt * dDv; s[1] = B + dt * dB; s[2] = N + dt * dN; s[3] = H + dt * dH;
    s[4] = A + dt * dA; s[5] = Dl + dt * dDl; s[6] = J + dt * dJ;
  },

  /** Ligand a cell presents to its neighbours. */
  ligand(s, p = crosstalk.defaults) { return s[5] + p.xJag * s[6]; },

  /** Integrate one isolated cell with fixed Wnt and ligand input to steady state. */
  steadyCell(W, L, p = crosstalk.defaults, T = 60, dt = 0.02) {
    const s = new Float64Array([0, 0.3, 0, 0, 0.5, 0.2, 0]);
    for (let t = 0; t < T; t += dt) crosstalk.stepCell(s, W, L, dt, p);
    return s;
  },
};

/**
 * Multicellular crypt: every cell runs the crosstalk ODEs; Notch input is
 * the mean ligand (Dll1 + weighted Jagged1) of its contact neighbours.
 */
export class CryptModel {
  constructor({ nbrs, wnt, rng, params = {}, noise = 0.25 }) {
    this.n = nbrs.length;
    this.nbrs = nbrs;
    this.wnt = Float64Array.from(wnt);
    this.p = { ...crosstalk.defaults, ...params };
    this.S = Array.from({ length: this.n }, () => new Float64Array(7));
    this.L = new Float64Array(this.n);
    this.t = 0;
    for (const s of this.S) {
      s[0] = 0; s[1] = 0.3; s[2] = 0.3 * (1 + noise * (rng.next() - 0.5));
      s[3] = 0.3; s[4] = 0.5 * (1 + noise * (rng.next() * 2 - 1)); s[5] = 0.3 * (1 + noise * (rng.next() * 2 - 1)); s[6] = 0;
    }
  }

  step(dt) {
    const { n, S, L, nbrs, p } = this;
    for (let i = 0; i < n; i++) {
      const nb = nbrs[i];
      let s = 0;
      for (let k = 0; k < nb.length; k++) s += crosstalk.ligand(S[nb[k]], p);
      L[i] = nb.length ? s / nb.length : 0;
    }
    for (let i = 0; i < n; i++) crosstalk.stepCell(S[i], this.wnt[i], L[i], dt, p);
    this.t += dt;
  }

  advance(T, dtMax = 0.05) {
    const m = Math.max(1, Math.ceil(T / dtMax));
    for (let k = 0; k < m; k++) this.step(T / m);
  }

  fate(i) { return fateIndex(this.S[i][1], this.S[i][2], this.p); }

  fates(out = new Int8Array(this.n)) {
    for (let i = 0; i < this.n; i++) out[i] = this.fate(i);
    return out;
  }

  fractions() {
    const f = [0, 0, 0, 0];
    for (let i = 0; i < this.n; i++) f[this.fate(i)]++;
    return f.map((c) => c / this.n);
  }
}

/**
 * Crypt geometry: a cylinder of `rings` hexagonally offset rings of `perRing`
 * cells capped by a hemispherical base. Returns 3D positions (y = crypt axis,
 * base at y = 0), arc length from the base pole and contact neighbours.
 */
export function cryptLattice({ perRing = 18, rings = 16, radius = 22 } = {}) {
  const pos = [], arc = [], normal = [];
  const dz = (2 * Math.PI * radius / perRing) * Math.sqrt(3) / 2; // hex row spacing
  const spacing = 2 * Math.PI * radius / perRing;
  // hemispherical cap: rings of decreasing circumference down to the pole
  const capRows = Math.max(1, Math.round((Math.PI / 2) * radius / dz));
  for (let r = 0; r <= capRows; r++) {
    const phi = (r / capRows) * (Math.PI / 2); // 0 at pole
    const rr = radius * Math.sin(phi);
    const count = r === 0 ? 1 : Math.max(3, Math.round((2 * Math.PI * rr) / spacing));
    for (let k = 0; k < count; k++) {
      const a = (2 * Math.PI * (k + (r % 2) * 0.5)) / count;
      const nx = Math.sin(phi) * Math.cos(a), ny = -Math.cos(phi), nz = Math.sin(phi) * Math.sin(a);
      pos.push([radius * nx, radius + radius * ny, radius * nz]);
      normal.push([nx, ny, nz]);
      arc.push(radius * phi);
    }
  }
  for (let r = 1; r <= rings; r++) {
    const y = radius + r * dz;
    for (let k = 0; k < perRing; k++) {
      const a = (2 * Math.PI * (k + ((capRows + r) % 2) * 0.5)) / perRing;
      pos.push([radius * Math.cos(a), y, radius * Math.sin(a)]);
      normal.push([Math.cos(a), 0, Math.sin(a)]);
      arc.push(radius * Math.PI / 2 + r * dz);
    }
  }
  // neighbours: centres within 1.35 x lattice spacing
  const cut2 = (1.35 * spacing) ** 2;
  const nbrs = pos.map(() => []);
  for (let i = 0; i < pos.length; i++) {
    for (let j = i + 1; j < pos.length; j++) {
      const d2 = (pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2 + (pos[i][2] - pos[j][2]) ** 2;
      if (d2 < cut2) { nbrs[i].push(j); nbrs[j].push(i); }
    }
  }
  return { pos, arc, normal, nbrs, spacing, length: radius * Math.PI / 2 + rings * dz, radius };
}

export const REACTIONS = [
  { id: 'xt.dvl.nicd', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: ['Dishevelled*', 'NICD'], products: ['Dishevelled:NICD'], label: 'Dishevelled binds NICD, damping Notch signalling' },
  { id: 'xt.gsk.nicd', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: ['NICD'], products: ['p-NICD'], modifiers: ['GSK-3b'], label: 'GSK-3beta phosphorylates NICD' },
  { id: 'xt.jag1', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: [], products: ['Jagged1'], modifiers: ['beta-cat:TCF'], label: 'beta-catenin/TCF induces Jagged1 (Wnt activates Notch in neighbours)' },
  { id: 'xt.hes.atoh', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: [], products: ['Atoh1'], modifiers: ['Hes1'], label: 'Hes1 represses Atoh1 (Math1)' },
  { id: 'xt.atoh.dll', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: [], products: ['Delta/Jagged'], modifiers: ['Atoh1'], label: 'Atoh1 induces Dll1 (lateral inhibition of neighbours)' },
  { id: 'xt.enh.both', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: ['beta-cat:TCF', 'NICD:CSL:MAML'], products: ['Stem-cell program'], label: 'Combinatorial enhancer: TCF + RBPJ both bound (Lgr5/Olfm4, stem)' },
  { id: 'xt.enh.wnt', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: ['beta-cat:TCF', 'Atoh1'], products: ['Paneth program'], label: 'TCF bound, RBPJ free: Paneth / secretory (Sox9, defensins)' },
  { id: 'xt.enh.notch', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: ['NICD:CSL:MAML'], products: ['Absorptive program'], modifiers: ['Hes1'], label: 'RBPJ bound, TCF free: absorptive enterocyte' },
  { id: 'xt.enh.none', order: 4, pathway: 'Wnt-Notch crosstalk', reactants: ['Atoh1'], products: ['Goblet/EE program'], label: 'Neither bound: goblet / enteroendocrine differentiation' },
];
