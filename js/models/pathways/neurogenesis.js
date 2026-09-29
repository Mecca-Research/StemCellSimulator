// Neurogenesis (paper: "Stem cell -(Sox2, Pax6)-> Neural progenitor
// -(Neurogenin, NeuroD)-> Neuron"; "DNA(Neurogenin) -> mRNA -> protein ->
// activation of neuronal genes"; morphogenetic processes "neurite outgrowth,
// synaptogenesis, axonal guidance"; "Cell Movement = Chemotaxis + Haptotaxis").
//
// 1) paperNgn: the paper's transcription-translation equations, EXACTLY:
//      d mRNA_Ngn / dt = k_tx - k_dm * mRNA_Ngn
//      d Ngn / dt      = k_tl * mRNA_Ngn - k_dp * Ngn
//    closed-form steady state m* = k_tx/k_dm, P* = k_tx k_tl / (k_dm k_dp) and
//    the analytic time course for constant rates (paperNgn.analytic).
//
// 2) neuralGRN: the same two equations embedded in the progenitor network.
//    Established wiring (Bertrand, Castro & Guillemot 2002, Nat Rev Neurosci
//    3:517; Kageyama et al. 2008, Nat Neurosci 11:1247; Scardigli et al. 2003,
//    Development 130:3269 for Pax6 -> Ngn2; Bylund et al. 2003, Nat Neurosci
//    6:1162 for SoxB1 vs proneural antagonism):
//      Notch -> Hes1 -| Neurogenin transcription;  Hes1 maintains Sox2
//      Sox2 -> Pax6 -> Neurogenin transcription (k_tx of the paper's ODE)
//      Neurogenin -> NeuroD (commitment), Neurogenin -> Delta (lateral inhibition)
//      NeuroD -| Sox2, Pax6
//    The kinetic form (Hill functions, normalised units, hours) is
//    phenomenological; a cell is scored as committed (post-mitotic neuron)
//    once NeuroD exceeds a threshold. Hes1 oscillations are not modelled.
//
// 3) guidance: steady state of the reaction-diffusion equation of the paper
//      dC/dt = D lap C + R(C),  R(C) = -lambda C + sum_s Q_s delta(x - x_s)
//    (first-order uptake/decay plus point secretion) in unbounded 3D:
//      C(x) = sum_s Q_s / (4 pi D r_s) exp(-r_s / l),  l = sqrt(D / lambda)
//      grad C = -sum_s Q_s exp(-r_s/l) (1 + r_s/l) / (4 pi D r_s^3) (x - x_s)
//    Growth cones read the gradient through receptor binding; the
//    signal-to-noise of gradient sensing across a growth cone of width w with
//    N receptors follows Berg & Purcell 1977 (Biophys J 20:193) as applied to
//    axons by Goodhill & Urbach 1999 (J Neurobiol 41:230):
//      p = C/(C+Kd),  dp = Kd w |grad C| / (C+Kd)^2,  SNR = dp sqrt(N / (4 p (1-p)))
//    The steering rule is phenomenological:
//      d <- normalise( persistence d + kappa SNR/(1+SNR) grad C/|grad C|
//                      + eta_hapto h + sigma xi )
//    with branching as a Poisson process and a synapse on contact with a target.
import { hill, hillR } from '../core/ode.js';
import { RNG } from '../core/rng.js';

// ----------------------------------------------------------- paper's ODE
export const paperNgn = {
  species: ['mRNA_Ngn', 'Ngn'],
  labels: ['Neurogenin mRNA', 'Neurogenin protein'],
  // 1/h; short-lived mRNA and protein (Ngn2 is degraded by ubiquitin-mediated
  // proteolysis, Vosper et al. 2007 Biochem J 407:277) - illustrative values
  defaults: { ktx: 1.2, kdm: 1.2, ktl: 2.0, kdp: 1.6 },
  rhs(p) {
    return (t, y, dy) => {
      dy[0] = p.ktx - p.kdm * y[0];
      dy[1] = p.ktl * y[0] - p.kdp * y[1];
    };
  },
  steadyState(p) {
    const m = p.ktx / p.kdm;
    return [m, (p.ktl * m) / p.kdp];
  },
  /** Exact solution for constant rates from (m0, P0) at t = 0. */
  analytic(t, p, m0 = 0, P0 = 0) {
    const [ms, Ps] = paperNgn.steadyState(p);
    const em = Math.exp(-p.kdm * t), ep = Math.exp(-p.kdp * t);
    const m = ms + (m0 - ms) * em;
    let P;
    if (Math.abs(p.kdp - p.kdm) > 1e-9) {
      const A = (p.ktl * (m0 - ms)) / (p.kdp - p.kdm);
      P = Ps + A * em + (P0 - Ps - A) * ep;
    } else {
      P = Ps + (P0 - Ps) * ep + p.ktl * (m0 - ms) * t * ep;
    }
    return [m, P];
  },
};

// ------------------------------------------------ progenitor network (GRN)
export const neuralGRN = {
  species: ['Sox2', 'Pax6', 'mRNA_Ngn', 'Ngn', 'NeuroD', 'Delta', 'Hes1'],
  labels: ['Sox2', 'Pax6', 'Neurogenin mRNA', 'Neurogenin', 'NeuroD', 'Delta (Dll1)', 'Hes1'],
  defaults: {
    kS: 0.25, kX: 0.2,          // Sox2, Pax6 turnover (1/h)
    ...paperNgn.defaults,        // k_tx (max), k_dm, k_tl, k_dp of the paper's ODE
    kN: 0.35, kDl: 0.8, kH: 1.0, // NeuroD, Delta, Hes1 turnover
    KHS: 0.3, KNS: 0.4,          // Hes1 -> Sox2, NeuroD -| Sox2
    KSX: 0.4, KNX: 0.5,          // Sox2 -> Pax6, NeuroD -| Pax6
    KXm: 0.35, KHm: 0.12,        // Pax6 -> Ngn transcription, Hes1 -| Ngn transcription
    KPN: 0.55, nPN: 3,           // Ngn -> NeuroD
    KPD: 0.5,                    // Ngn -> Delta
    KnH: 0.12,                   // Notch activation -> Hes1 (sharp, cf. Collier et al. 1996)
    coupling: 1.0,               // Delta-Notch coupling (0 = gamma-secretase inhibitor, DAPT)
    notchExt: 0.0,               // exogenous Notch activation (e.g. immobilised ligand)
    drive: 1.0,                  // neurogenic drive, multiplies k_tx
    commit: 0.5,                 // NeuroD level scored as commitment
  },

  /** Initial progenitor state (Sox2/Pax6 high, small random Neurogenin). */
  initial(rng) {
    return Float64Array.from([0.8 + 0.1 * rng.next(), 0.75 + 0.1 * rng.next(), 0.1 * rng.next(), 0.1 * rng.next(), 0, 0.3 * rng.next(), 0.4 + 0.1 * rng.next()]);
  },

  /** Effective transcription rate of Neurogenin (the paper's k_transcription). */
  ktxEff(y, p) { return p.drive * p.ktx * hill(y[1], p.KXm, 2) * hillR(y[6], p.KHm, 2); },

  /** dy/dt for one cell given its Notch input (mean neighbour Delta already coupled). */
  deriv(y, notchIn, committed, p, dy) {
    const [S, X, m, P, N, Dl, H] = y;
    dy[0] = p.kS * ((0.2 + 0.8 * hill(H, p.KHS, 2)) * hillR(N, p.KNS, 2) - S);
    dy[1] = p.kX * ((0.3 + 0.7 * hill(S, p.KSX, 2)) * hillR(N, p.KNX, 2) - X);
    dy[2] = neuralGRN.ktxEff(y, p) - p.kdm * m;
    dy[3] = p.ktl * m - p.kdp * P;
    // after commitment the neuronal programme maintains NeuroD and Delta is lost
    dy[4] = p.kN * ((committed ? 1 : hill(P, p.KPN, p.nPN)) - N);
    dy[5] = p.kDl * ((committed ? 0 : hill(P, p.KPD, 2)) - Dl);
    dy[6] = p.kH * (hill(notchIn, p.KnH, 2) - H);
    return dy;
  },

  /**
   * Advance a population one step (midpoint rule, neighbour Delta frozen over dt).
   * @param {Float64Array[]} Y states
   * @param {boolean[]} committed flags (updated in place when NeuroD > commit)
   * @param {number[][]} nbrs neighbour index lists
   * @returns {number[]} indices that committed during this step
   */
  step(Y, committed, nbrs, dt, p, work = neuralGRN.work()) {
    const n = Y.length;
    if (work.nIn.length < n) work.nIn = new Float64Array(n * 2);
    for (let i = 0; i < n; i++) {
      const nb = nbrs[i];
      let avg = 0;
      if (nb && nb.length) { for (const j of nb) avg += Y[j][5]; avg /= nb.length; }
      work.nIn[i] = p.coupling * avg + p.notchExt;
    }
    const fresh = [];
    const { k1, tmp } = work;
    for (let i = 0; i < n; i++) {
      const y = Y[i];
      neuralGRN.deriv(y, work.nIn[i], committed[i], p, k1);
      for (let s = 0; s < 7; s++) tmp[s] = y[s] + 0.5 * dt * k1[s];
      neuralGRN.deriv(tmp, work.nIn[i], committed[i], p, k1);
      for (let s = 0; s < 7; s++) y[s] = Math.max(0, y[s] + dt * k1[s]);
      if (!committed[i] && y[4] > p.commit) { committed[i] = true; fresh.push(i); }
    }
    return fresh;
  },

  work() { return { nIn: new Float64Array(256), k1: new Float64Array(7), tmp: new Float64Array(7) }; },
};

// ------------------------------------------------------- guidance field
export const NM_PER_MOLECULE_UM3 = 1 / 0.6022; // 1 molecule/um^3 = 1.66 nM

export const guidance = {
  defaults: {
    D: 5,          // diffusion coefficient (um^2/s) of a netrin-sized protein in tissue
    lambda: 2,     // first-order removal (1/h): degradation + binding to ECM
    Q: 1500,       // secretion per target cell (molecules/s)
    Kd: 1,         // receptor dissociation constant (nM), e.g. netrin-1/DCC ~ nM
    receptors: 10000,
    width: 10,     // growth-cone width (um)
  },

  /** Decay length l = sqrt(D/lambda) in um (D in um^2/s, lambda in 1/h). */
  length(p) { return Math.sqrt((p.D * 3600) / p.lambda); },

  /** Concentration (nM) at x from sources [{p:[x,y,z], radius, Q?}]. */
  concentration(x, sources, p) {
    const l = guidance.length(p);
    let C = 0;
    for (const s of sources) {
      const r = Math.max(Math.hypot(x[0] - s.p[0], x[1] - s.p[1], x[2] - s.p[2]), s.radius ?? 1e-3);
      C += ((s.Q ?? p.Q) / (4 * Math.PI * p.D * r)) * Math.exp(-r / l);
    }
    return C * NM_PER_MOLECULE_UM3;
  },

  /** Analytic gradient (nM/um) written into out; returns C (nM). */
  gradient(x, sources, p, out = [0, 0, 0]) {
    const l = guidance.length(p);
    out[0] = out[1] = out[2] = 0;
    let C = 0;
    for (const s of sources) {
      const dx = x[0] - s.p[0], dy = x[1] - s.p[1], dz = x[2] - s.p[2];
      const r0 = Math.hypot(dx, dy, dz);
      const a = s.radius ?? 1e-3;
      const r = Math.max(r0, a);
      const Q = s.Q ?? p.Q;
      const e = Math.exp(-r / l);
      C += (Q / (4 * Math.PI * p.D * r)) * e;
      if (r0 >= a && r0 > 0) {
        const k = (-Q * e * (1 + r / l)) / (4 * Math.PI * p.D * r * r * r);
        out[0] += k * dx; out[1] += k * dy; out[2] += k * dz;
      }
    }
    out[0] *= NM_PER_MOLECULE_UM3; out[1] *= NM_PER_MOLECULE_UM3; out[2] *= NM_PER_MOLECULE_UM3;
    return C * NM_PER_MOLECULE_UM3;
  },

  /** Signal-to-noise of gradient detection across the growth cone (Berg-Purcell / Goodhill-Urbach). */
  snr(C, gradMag, p) {
    if (C <= 0) return 0;
    const occ = C / (C + p.Kd);
    const dp = (p.Kd * p.width * gradMag) / ((C + p.Kd) * (C + p.Kd));
    return dp * Math.sqrt(p.receptors / (4 * occ * (1 - occ) + 1e-12));
  },

  /** Radius around an isolated source at which C falls to level (nM) (bisection). */
  isoRadius(level, p, Q = p.Q, rMin = 0.5, rMax = 2000) {
    const src = [{ p: [0, 0, 0], Q, radius: 1e-6 }];
    let lo = rMin, hi = rMax;
    if (guidance.concentration([lo, 0, 0], src, p) < level) return lo;
    if (guidance.concentration([hi, 0, 0], src, p) > level) return hi;
    for (let i = 0; i < 60; i++) {
      const mid = 0.5 * (lo + hi);
      if (guidance.concentration([mid, 0, 0], src, p) > level) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  },
};

// ------------------------------------------------------ neurite outgrowth
export const OUTGROWTH_DEFAULTS = {
  step: 2.0,             // segment length (um)
  axonSpeed: 40,         // um/h (typical axon elongation in culture: ~1 mm/day)
  dendriteSpeed: 9,
  kappa: 0.4,            // chemotactic steering gain
  persistence: 1.0,
  noise: 0.22,           // random turning per segment
  hapto: 0.0,            // haptotactic bias gain (bound ECM gradient)
  haptoDir: [1, 0, 0],   // direction of increasing substrate-bound ligand
  branchRate: 0.12,      // axon branching probability per hour per growth cone
  maxAxonBranches: 5,
  dendrites: 4,
  dendriteLen: [22, 45],
  dendriteBranchRate: 0.12,
  maxAxonLen: 520,
  bound: 330,            // growth cones collapse beyond this radius from the origin
  contact: 1.5,          // extra reach for synapse formation (um)
  maxSegments: Infinity, // rendering budget: growth stops when reached
};

/**
 * Growth cones of many neurons advancing through the guidance field.
 * Segments are appended (never moved) so they can be drawn incrementally.
 */
export class NeuriteGrowth {
  constructor(o = {}) {
    this.p = { ...OUTGROWTH_DEFAULTS, ...(o.outgrowth ?? {}) };
    this.gp = { ...guidance.defaults, ...(o.guidance ?? {}) };
    this.rng = o.rng ?? new RNG(o.seed ?? 99);
    this.tips = [];
    this.segments = [];
    this.synapses = [];
    this.totalLength = 0;
    this.t = 0;
    this._g = [0, 0, 0];
    this.neurons = new Map(); // id -> {axonBranches, dendriteBranches}
  }

  /** Start neurites of a new neuron at the soma surface. */
  addNeuron(id, soma, R, outward) {
    const rng = this.rng;
    const o = norm3([...outward]);
    this.neurons.set(id, { axonBranches: 0, synapses: 0 });
    // axon: leaves along the outward normal, biased by the local gradient
    const g = this._g;
    guidance.gradient(soma, this.sources ?? [], this.gp, g);
    const gn = Math.hypot(g[0], g[1], g[2]);
    const d = gn > 0 ? norm3([o[0] + 0.6 * g[0] / gn, o[1] + 0.6 * g[1] / gn, o[2] + 0.6 * g[2] / gn]) : o;
    this.tips.push(this.tip(id, 0, add3(soma, d, R), d, 1.1, 0, this.p.maxAxonLen));
    // dendrites: short, unguided, spread around the soma away from the axon
    for (let k = 0; k < this.p.dendrites; k++) {
      const r = rng.onSphere();
      const dd = norm3([r[0] - 0.8 * d[0] + 0.3 * o[0], r[1] - 0.8 * d[1] + 0.3 * o[1], r[2] - 0.8 * d[2] + 0.3 * o[2]]);
      const [a, b] = this.p.dendriteLen;
      this.tips.push(this.tip(id, 1, add3(soma, dd, R * 0.95), dd, 0.9, 0, a + (b - a) * rng.next()));
    }
  }

  tip(neuron, kind, pos, dir, radius, order, maxLen) {
    return { neuron, kind, pos: pos.slice(), dir: dir.slice(), radius, order, len: 0, maxLen, acc: 0, active: true, snr: 0, C: 0, born: this.t };
  }

  /**
   * @param {number} dt hours
   * @param {{p:number[], radius:number, Q?:number, id?:any}[]} sources chemoattractant-secreting targets
   * @returns {number} segments added
   */
  step(dt, sources) {
    this.sources = sources;
    const p = this.p, rng = this.rng, g = this._g;
    const before = this.segments.length;
    const born = [];
    for (const tp of this.tips) {
      if (!tp.active) continue;
      if (this.segments.length >= p.maxSegments) { tp.active = false; continue; }
      const speed = tp.kind === 0 ? p.axonSpeed : p.dendriteSpeed;
      tp.acc += speed * dt * (0.75 + 0.5 * rng.next());
      while (tp.acc >= p.step && tp.active) {
        tp.acc -= p.step;
        const d = tp.dir;
        let nx = p.persistence * d[0], ny = p.persistence * d[1], nz = p.persistence * d[2];
        if (tp.kind === 0) {
          const C = guidance.gradient(tp.pos, sources, this.gp, g);
          const gm = Math.hypot(g[0], g[1], g[2]);
          const snr = guidance.snr(C, gm, this.gp);
          tp.snr = snr; tp.C = C;
          if (gm > 0) {
            const w = (p.kappa * snr) / (1 + snr) / gm;
            nx += w * g[0]; ny += w * g[1]; nz += w * g[2];
          }
          if (p.hapto) { nx += p.hapto * p.haptoDir[0]; ny += p.hapto * p.haptoDir[1]; nz += p.hapto * p.haptoDir[2]; }
        }
        const sig = tp.kind === 0 ? p.noise : p.noise * 1.3;
        nx += sig * rng.normal(); ny += sig * rng.normal(); nz += sig * rng.normal();
        const nn = Math.hypot(nx, ny, nz) || 1;
        d[0] = nx / nn; d[1] = ny / nn; d[2] = nz / nn;
        const a = tp.pos.slice();
        tp.pos[0] += d[0] * p.step; tp.pos[1] += d[1] * p.step; tp.pos[2] += d[2] * p.step;
        tp.len += p.step;
        this.totalLength += p.step;
        // taper: axons thin with distance, dendrites taper to their tips
        const r = tp.kind === 0 ? Math.max(0.55, tp.radius * Math.exp(-tp.len / 400))
          : Math.max(0.35, tp.radius * (1 - tp.len / (tp.maxLen * 1.3)));
        this.segments.push({ a, b: tp.pos.slice(), r, kind: tp.kind, neuron: tp.neuron, order: tp.order });
        // synapse on contact with a target
        if (tp.kind === 0) {
          for (const s of sources) {
            const dx = tp.pos[0] - s.p[0], dy = tp.pos[1] - s.p[1], dz = tp.pos[2] - s.p[2];
            const dist = Math.hypot(dx, dy, dz);
            if (dist < s.radius + p.contact) {
              const k = s.radius / dist;
              this.synapses.push({ pos: [s.p[0] + dx * k, s.p[1] + dy * k, s.p[2] + dz * k], target: s.id, neuron: tp.neuron, t: this.t });
              tp.active = false; tp.synapse = true;
              const nr = this.neurons.get(tp.neuron); if (nr) nr.synapses++;
              break;
            }
          }
        }
        if (tp.len >= tp.maxLen || Math.hypot(tp.pos[0], tp.pos[1], tp.pos[2]) > p.bound) tp.active = false;
      }
      // branching (Poisson in time)
      if (tp.active) {
        const nr = this.neurons.get(tp.neuron);
        if (tp.kind === 0 && nr && nr.axonBranches < p.maxAxonBranches && tp.len > 20 && rng.next() < p.branchRate * dt) {
          nr.axonBranches++;
          born.push(this.tip(tp.neuron, 0, tp.pos, turn(tp.dir, 0.6 + 0.4 * rng.next(), rng), tp.radius * 0.75, tp.order + 1, p.maxAxonLen * 0.7));
        } else if (tp.kind === 1 && tp.order < 2 && tp.len > 8 && rng.next() < p.dendriteBranchRate * dt) {
          born.push(this.tip(tp.neuron, 1, tp.pos, turn(tp.dir, 0.5 + 0.4 * rng.next(), rng), tp.radius * 0.8, tp.order + 1, tp.maxLen * 0.6));
          tp.dir = turn(tp.dir, 0.3, rng);
        }
      }
    }
    this.tips.push(...born);
    this.t += dt;
    return this.segments.length - before;
  }

  activeTips() { return this.tips.filter((t) => t.active); }
}

function norm3(v) { const n = Math.hypot(v[0], v[1], v[2]) || 1; v[0] /= n; v[1] /= n; v[2] /= n; return v; }
function add3(p, d, s) { return [p[0] + d[0] * s, p[1] + d[1] * s, p[2] + d[2] * s]; }
/** Rotate unit vector d by angle a (rad) about a random axis perpendicular to it. */
function turn(d, a, rng) {
  const r = rng.onSphere();
  const k = r[0] * d[0] + r[1] * d[1] + r[2] * d[2];
  const u = norm3([r[0] - k * d[0], r[1] - k * d[1], r[2] - k * d[2]]);
  return norm3([d[0] * Math.cos(a) + u[0] * Math.sin(a), d[1] * Math.cos(a) + u[1] * Math.sin(a), d[2] * Math.cos(a) + u[2] * Math.sin(a)]);
}

/** Reaction list for the network view (order 3 = differentiated cells). */
export const REACTIONS = [
  { id: 'neuro.sox2', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['Sox2'], modifiers: ['Hes1'], label: 'Sox2 maintenance by Notch/Hes (progenitor)' },
  { id: 'neuro.pax6', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['Pax6'], modifiers: ['Sox2'], label: 'Pax6 expression in neural progenitors' },
  { id: 'neuro.tx', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['Neurogenin mRNA'], modifiers: ['DNA (Neurogenin)', 'Pax6', 'Hes1'], label: 'Neurogenin transcription (k_tx)' },
  { id: 'neuro.mdeg', order: 3, pathway: 'Neurogenesis', reactants: ['Neurogenin mRNA'], products: [], label: 'mRNA degradation (k_dm)' },
  { id: 'neuro.tl', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['Neurogenin'], modifiers: ['Neurogenin mRNA'], label: 'Translation (k_tl)' },
  { id: 'neuro.pdeg', order: 3, pathway: 'Neurogenesis', reactants: ['Neurogenin'], products: [], label: 'Protein degradation (k_dp)' },
  { id: 'neuro.neurod', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['NeuroD'], modifiers: ['Neurogenin'], label: 'NeuroD activation (commitment)' },
  { id: 'neuro.delta', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['Delta/Jagged'], modifiers: ['Neurogenin'], label: 'Proneural induction of Delta (lateral inhibition)' },
  { id: 'neuro.hes', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['Hes1'], modifiers: ['NICD:CSL:MAML'], label: 'Notch target Hes1 represses proneural genes' },
  { id: 'neuro.soxrep', order: 3, pathway: 'Neurogenesis', reactants: ['Sox2'], products: [], modifiers: ['NeuroD'], label: 'Proneural down-regulation of Sox2' },
  { id: 'neuro.genes', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['neuronal genes (beta-III tubulin)'], modifiers: ['NeuroD'], label: 'Activation of neuronal genes' },
  { id: 'neuro.secrete', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['chemoattractant (netrin-1)'], modifiers: ['target cell'], label: 'Target-derived chemoattractant secretion (Q)' },
  { id: 'neuro.decay', order: 3, pathway: 'Neurogenesis', reactants: ['chemoattractant (netrin-1)'], products: [], label: 'Chemoattractant removal (-lambda C)' },
  { id: 'neuro.bind', order: 3, pathway: 'Neurogenesis', reactants: ['chemoattractant (netrin-1)', 'DCC'], products: ['chemoattractant:DCC'], label: 'Guidance-receptor binding at the growth cone' },
  { id: 'neuro.steer', order: 3, pathway: 'Neurogenesis', reactants: [], products: ['neurite'], modifiers: ['chemoattractant:DCC', 'laminin'], label: 'Neurite outgrowth (chemotaxis + haptotaxis)' },
  { id: 'neuro.synapse', order: 3, pathway: 'Neurogenesis', reactants: ['axon terminal', 'target cell'], products: ['synapse'], label: 'Synaptogenesis on contact' },
];
