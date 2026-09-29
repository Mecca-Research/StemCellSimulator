// Hedgehog signalling (paper: "Stage 1 - Hedgehog Signaling Pathway": Shh binds
// Patched -> Patched inhibition relieved -> Smoothened activation -> GLI
// transcription factors).
//
// Mechanistic single-cell model (nondimensional; concentrations in units of the
// Shh-PTCH1 dissociation constant, time in model units). Structure:
//
//   1. Shh + PTCH1 <-> Shh:PTCH1 is fast, so the bound fraction is in
//      quasi-steady state:   f_b = Shh / (Shh + Kd),   PTCH1_free = (1 - f_b) P.
//      Ligand-bound PTCH1 is internalised and degraded (Incardona et al. 2000,
//      PNAS 97:12044).
//   2. Free PTCH1 inhibits Smoothened *catalytically* - sub-stoichiometric PTCH1
//      suppresses a large excess of SMO (Taipale, Cooper, Maiti & Beachy 2002,
//      Nature 418:892), so it enters as an enzyme-like inactivation rate:
//        dS/dt = kSa (1 - S) - (kSi + kPtc PTCH1_free) S
//      (S = active SMO fraction; cyclopamine binds SMO and blocks its
//      activation, Chen, Taipale, Cooper & Beachy 2002 Genes Dev 16:2743:
//      kSa -> 0).
//   3. Full-length GLI2/3 (G) is either converted into the activator GliA at a
//      SMO-dependent rate, or - when SMO is off - phosphorylated by PKA / GSK-3b
//      / CK1 and partially processed by beta-TrCP into the Gli3 repressor GliR.
//      GliA also includes Gli1, a transcriptional target that is itself an
//      activator (positive feedback; the circuit stays monostable for the
//      default parameters - see tests - unlike the bistable regime explored by
//      Lai et al. 2004).
//   4. GLI target promoters integrate activator and repressor by a
//      Shea-Ackers occupancy function, as in the GLI circuit models of
//      Lai, Robertson & Schaffer 2004 (Biophys J 86:2748) and
//      Saha & Schaffer 2006 (Development 133:889):
//        phi(A, R) = (beta + (A/KA)^n) / (1 + (A/KA)^n + (R/KR)^n)
//   5. Ptch1 is itself a GLI target, so signalling raises the receptor that
//      inhibits it: the negative feedback that desensitises cells over time
//      and makes them read out both Shh concentration and exposure duration
//      (temporal adaptation; Dessaud et al. 2007, Nature 450:717).
//
//   dP/dt = sP0 + sP phi - dP P - kInt f_b P
//   dS/dt = kSa (1 - S) - (kSi + kPtc (1 - f_b) P) S
//   dG/dt = sG - (kGA S + kGR (1 - S) + dG) G
//   dA/dt = kGA S G + sA1 phi - dA A
//   dR/dt = kGR (1 - S) G - dR R
//   dT/dt = kT phi - dT T                  (target, e.g. a Nkx2.2/Olig2 reporter)
//
// Parameter values are illustrative (nondimensional), chosen by a small
// parameter search so that the steady-state dose-response spans ~0.3-30 Kd and
// the Ptch1 feedback produces a visible adaptation overshoot (~30 %); they are
// not fitted to data. The *structure* follows the cited models.
// The ODE lumps nuclear and cytoplasmic GLI (the REACTIONS list keeps the
// translocation step for the network view).

export const hedgehog = {
  species: ['PTCH1', 'SMO*', 'GLI', 'GliA', 'GliR', 'target'],
  labels: ['PTCH1 (total)', 'SMO active fraction', 'GLI2/3 full-length', 'GLI activator', 'GLI3 repressor', 'GLI target'],
  defaults: {
    shh: 0,                     // extracellular Shh, units of Kd
    Kd: 1,                      // Shh:PTCH1 dissociation constant
    sP0: 0.1, sP: 0.4,          // basal and GLI-induced Ptch1 synthesis (sP = 0: no feedback)
    dP: 0.08, kInt: 0.05,       // PTCH1 turnover, ligand-induced internalisation
    kSa: 1.0, kSi: 0.05,        // SMO activation (0 = cyclopamine), basal inactivation
    kPtc: 13,                   // catalytic SMO inhibition by free PTCH1
    sG: 0.5, dG: 0.1,           // full-length GLI2/3 synthesis and turnover
    kGA: 1.0, kGR: 0.26,        // SMO-dependent activation vs repressor processing
    dA: 0.5, dR: 0.45,          // GliA (SPOP/Cul3) and GliR turnover
    sA1: 0.6,                   // Gli1 transcription (activator; positive feedback)
    KA: 0.4, KR: 0.28, nH: 2, beta: 0.05, // promoter occupancy
    kT: 1.0, dT: 0.3,           // target gene
  },

  /** Quasi-steady-state fraction of PTCH1 bound by Shh. */
  boundFraction(shh, Kd = 1) { return shh > 0 ? shh / (shh + Kd) : 0; },

  /** Shea-Ackers promoter occupancy by GLI activator A and repressor R. */
  promoter(A, R, p = hedgehog.defaults) {
    const a = Math.pow(Math.max(A, 0) / p.KA, p.nH), r = Math.pow(Math.max(R, 0) / p.KR, p.nH);
    return (p.beta + a) / (1 + a + r);
  },

  /** Analytic quasi-steady SMO activity for a given PTCH1 level (used by tests). */
  smoSteady(P, p = hedgehog.defaults) {
    const Pf = (1 - hedgehog.boundFraction(p.shh, p.Kd)) * P;
    return p.kSa / (p.kSa + p.kSi + p.kPtc * Pf);
  },

  rhs(p) {
    return (t, y, dy) => {
      const P = y[0], S = y[1], G = y[2], A = y[3], R = y[4], T = y[5];
      const fb = p.shh > 0 ? p.shh / (p.shh + p.Kd) : 0;
      const a = Math.pow(Math.max(A, 0) / p.KA, p.nH), r = Math.pow(Math.max(R, 0) / p.KR, p.nH);
      const phi = (p.beta + a) / (1 + a + r);
      dy[0] = p.sP0 + p.sP * phi - p.dP * P - p.kInt * fb * P;
      dy[1] = p.kSa * (1 - S) - (p.kSi + p.kPtc * (1 - fb) * P) * S;
      dy[2] = p.sG - (p.kGA * S + p.kGR * (1 - S) + p.dG) * G;
      dy[3] = p.kGA * S * G + p.sA1 * phi - p.dA * A;
      dy[4] = p.kGR * (1 - S) * G - p.dR * R;
      dy[5] = p.kT * phi - p.dT * T;
    };
  },

  /** Unstimulated cell (Shh = 0) steady state, a sensible initial condition. */
  initial() { return [1.406, 0.052, 1.255, 0.167, 0.688, 0.104]; },

  /** Transcriptional GLI read-out (promoter occupancy, 0..1) of a state vector. */
  readout(y, p = hedgehog.defaults) { return hedgehog.promoter(y[3], y[4], p); },

  /**
   * Steady state for parameters p by fixed-step RK4 relaxation (the model is
   * monostable for the default parameters; see tests).
   */
  steadyState(p, { y0 = hedgehog.initial(), tEnd = 400, h = 0.025 } = {}) {
    const f = hedgehog.rhs(p);
    const y = Float64Array.from(y0);
    const n = y.length;
    const k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n), k4 = new Float64Array(n), tmp = new Float64Array(n);
    const steps = Math.ceil(tEnd / h);
    for (let s = 0; s < steps; s++) {
      f(0, y, k1);
      for (let i = 0; i < n; i++) tmp[i] = y[i] + 0.5 * h * k1[i];
      f(0, tmp, k2);
      for (let i = 0; i < n; i++) tmp[i] = y[i] + 0.5 * h * k2[i];
      f(0, tmp, k3);
      for (let i = 0; i < n; i++) tmp[i] = y[i] + h * k3[i];
      f(0, tmp, k4);
      for (let i = 0; i < n; i++) y[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    return Array.from(y);
  },

  /**
   * Steady-state dose-response: for each Shh dose returns
   * { shh, P, S, G, A, R, T, readout }.
   */
  doseResponse(p, doses, opts) {
    return doses.map((shh) => {
      const q = { ...p, shh };
      const y = hedgehog.steadyState(q, opts);
      return { shh, P: y[0], S: y[1], G: y[2], A: y[3], R: y[4], T: y[5], readout: hedgehog.promoter(y[3], y[4], q) };
    });
  },
};

// ---------------------------------------------------------------- French flag
// Wolpert's French-flag model (Wolpert 1969, J Theor Biol 25:1): cells read a
// graded signal against fixed thresholds and adopt a discrete domain identity.
// Applied to Shh in the ventral neural tube, GLI activity above successive
// thresholds specifies Pax6/Pax7 (dorsal, Hh-off) < Olig2 (pMN) < Nkx2.2 (p3)
// < FoxA2 (floor plate) (Dessaud, McMahon & Briscoe 2008, Development 135:2489).
// The thresholds below are illustrative, on the GLI promoter read-out (0..1).

export const NEURAL_TUBE_DOMAINS = [
  { name: 'Pax6/Pax7 · dorsal (Hh-off)', color: '#4a78ff' },
  { name: 'Olig2 · pMN (motor-neuron progenitor)', color: '#f4f1ea' },
  { name: 'Nkx2.2 · p3', color: '#ff4a5e' },
  { name: 'FoxA2 · floor plate', color: '#ffc04a' },
];
export const NEURAL_TUBE_THRESHOLDS = [0.2, 0.45, 0.7];

/**
 * Domain index of a signal value: the number of (ascending) thresholds it
 * reaches. thresholds [t1, t2] -> 0 (below t1), 1 (t1..t2), 2 (>= t2).
 */
export function frenchFlag(value, thresholds = NEURAL_TUBE_THRESHOLDS) {
  let k = 0;
  while (k < thresholds.length && value >= thresholds[k]) k++;
  return k;
}

/** Vectorised French-flag read-out; returns an Int8Array of domain indices. */
export function frenchFlagField(values, thresholds = NEURAL_TUBE_THRESHOLDS) {
  const out = new Int8Array(values.length);
  for (let i = 0; i < values.length; i++) out[i] = frenchFlag(values[i], thresholds);
  return out;
}

/** Reaction list for the network view (order 1 = the paper's first stage). */
export const REACTIONS = [
  { id: 'hh.bind', order: 1, pathway: 'Hedgehog', reactants: ['Shh', 'PTCH1'], products: ['Shh:PTCH1'], label: 'Shh binds Patched' },
  { id: 'hh.endo', order: 1, pathway: 'Hedgehog', reactants: ['Shh:PTCH1'], products: [], label: 'Ligand-induced PTCH1 internalisation' },
  { id: 'hh.smo_act', order: 1, pathway: 'Hedgehog', reactants: ['SMO'], products: ['SMO*'], label: 'Smoothened activation (ciliary entry)' },
  { id: 'hh.ptc_smo', order: 1, pathway: 'Hedgehog', reactants: ['SMO*'], products: ['SMO'], modifiers: ['PTCH1'], label: 'Catalytic inhibition of SMO by Patched' },
  { id: 'hh.gli_syn', order: 1, pathway: 'Hedgehog', reactants: [], products: ['GLI (full-length)'], label: 'GLI2/3 synthesis' },
  { id: 'hh.gli_act', order: 1, pathway: 'Hedgehog', reactants: ['GLI (full-length)'], products: ['GliA'], modifiers: ['SMO*', 'SUFU'], label: 'GLI activation (SUFU release)' },
  { id: 'hh.gli_proc', order: 1, pathway: 'Hedgehog', reactants: ['GLI (full-length)'], products: ['GLI3R'], modifiers: ['PKA', 'GSK-3b', 'CK1', 'beta-TrCP'], label: 'GLI3 repressor processing' },
  { id: 'hh.gli_imp', order: 1, pathway: 'Hedgehog', reactants: ['GliA'], products: ['GliA (nuc)'], label: 'GLI nuclear translocation' },
  { id: 'hh.ptch_tx', order: 1, pathway: 'Hedgehog', reactants: [], products: ['PTCH1'], modifiers: ['GliA (nuc)', 'GLI3R'], label: 'Ptch1 induction (negative feedback)' },
  { id: 'hh.gli1_tx', order: 1, pathway: 'Hedgehog', reactants: [], products: ['GliA'], modifiers: ['GliA (nuc)', 'GLI3R'], label: 'Gli1 transcription (positive feedback)' },
  { id: 'hh.target', order: 1, pathway: 'Hedgehog', reactants: [], products: ['Hh target mRNA'], modifiers: ['GliA (nuc)', 'GLI3R'], label: 'Target gene transcription' },
  { id: 'hh.deg_a', order: 1, pathway: 'Hedgehog', reactants: ['GliA'], products: [], modifiers: ['SPOP'], label: 'GliA degradation' },
  { id: 'hh.deg_r', order: 1, pathway: 'Hedgehog', reactants: ['GLI3R'], products: [], label: 'GLI3R turnover' },
];
