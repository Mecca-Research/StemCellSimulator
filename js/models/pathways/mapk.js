// MAPK/ERK signalling (paper: "Stage 2 - Progenitor cell differentiation:
// MAPK/ERK pathway": receptor tyrosine kinase activation by growth factors,
// the kinase cascade Raf -> MEK -> ERK by sequential phosphorylation, ERK
// translocation to the nucleus and activation of transcription factors).
//
// 1) Core cascade: Huang & Ferrell 1996, PNAS 93:10078 "Ultrasensitivity in
//    the mitogen-activated protein kinase cascade", in the Michaelis-Menten
//    (quasi-steady-state) form. Every step is catalysed by an enzyme E acting
//    on substrate S with  v = kcat * E * (S/Km) / (1 + sum_j S_j/Km)  where the
//    denominator lists every substrate that competes for the same enzyme
//    (MAPKK and MAPKK-P both bind active MAPKKK; MAPKK-PP and MAPKK-P both bind
//    the MAPKK phosphatase, etc.):
//
//      MAPKKK  --E1-->  MAPKKK*          MAPKKK* --E2--> MAPKKK
//      MAPKK --MAPKKK*--> MAPKK-P --MAPKKK*--> MAPKK-PP      (distributive)
//      MAPKK-PP --PPase--> MAPKK-P --PPase--> MAPKK
//      MAPK  --MAPKK-PP--> MAPK-P --MAPKK-PP--> MAPK-PP      (distributive)
//      MAPK-PP --PPase--> MAPK-P --PPase--> MAPK
//
//    Conservation: MAPKKK + MAPKKK* = MAPKKK_tot, MAPKK + MAPKK-P + MAPKK-PP =
//    MAPKK_tot, MAPK + MAPK-P + MAPK-PP = MAPK_tot. Totals, kcat and Km are the
//    values Huang & Ferrell used (Xenopus oocyte concentrations; every Km =
//    300 nM, every kcat = 150 s^-1). In their Fig. 2 the effective Hill
//    coefficients were ~1.0 (MAPKKK), ~1.7 (MAPKK) and ~4.9 (MAPK); the
//    Michaelis-Menten form reproduces the same ordering (sequestration by the
//    enzyme-substrate complexes, which HF96 kept, is neglected here).
//
// 2) Upstream receptor (illustrative parameters): growth factor (SCF, the
//    c-Kit ligand in hematopoietic progenitors) binds the receptor tyrosine
//    kinase, ligand-bound monomers dimerise and trans-autophosphorylate; the
//    active dimer recruits Grb2/SOS and sets the Ras-GTP level that is the
//    cascade input E1.
//
// 3) Downstream: doubly-phosphorylated ERK enters the nucleus, is inactivated
//    there by a nuclear dual-specificity phosphatase (MKP-1) and exported;
//    nuclear ERK-PP phosphorylates Elk-1 which induces an immediate-early gene
//    (c-Fos). Rates are illustrative (minutes).
//
// Units: concentrations in uM, time in s.

/** Huang & Ferrell (1996) Table 1 values. */
export const HF96 = {
  MKKKtot: 0.003,   // MAPKKK (Mos/Raf) 3 nM
  MKKtot: 1.2,      // MAPKK (MEK) 1.2 uM
  MAPKtot: 1.2,     // MAPK (ERK) 1.2 uM
  E2: 0.0003,       // MAPKKK inactivating enzyme 0.3 nM
  MKKPase: 0.0003,  // MAPKK phosphatase 0.3 nM
  MAPKPase: 0.12,   // MAPK phosphatase 120 nM
  kcat: 150,        // s^-1, every step
  Km: 0.3,          // uM, every step
};

// --------------------------------------------------------------- steady state

/** One-step cycle  x0 <-> x1 (e.g. MAPKKK / MAPKKK*): returns active x1. */
export function oneStepSteadyState(E, P, T, kcat, Km) {
  if (E <= 0) return 0;
  if (P <= 0) return T;
  let lo = 0, hi = T;
  for (let i = 0; i < 200; i++) {
    const x1 = 0.5 * (lo + hi), x0 = T - x1;
    const f = (kcat * E * x0) / (Km + x0) - (kcat * P * x1) / (Km + x1);
    if (f > 0) lo = x1; else hi = x1;
    if (hi - lo < T * 1e-15) break;
  }
  return 0.5 * (lo + hi);
}

/**
 * Two-step distributive cycle x0 <-> x1 <-> x2 with one kinase E (acting on
 * x0 and x1 competitively) and one phosphatase P (acting on x2 and x1).
 * At steady state both net fluxes vanish, which forces x1/x0 = x2/x1 = r with
 *   r = (kcat E / kcat P) * D_P(r) / D_E(r),
 *   D_E = 1 + (x0 + x1)/Km,  D_P = 1 + (x1 + x2)/Km,  x0 = T / (1 + r + r^2).
 * Solved by bisection on ln r.
 */
export function twoStepSteadyState(E, P, T, kcat, Km, kcatP = kcat) {
  if (E <= 0) return { x0: T, x1: 0, x2: 0 };
  if (P <= 0) return { x0: 0, x1: 0, x2: T };
  const c = (kcat * E) / (kcatP * P);
  const g = (lr) => {
    const r = Math.exp(lr);
    const x0 = T / (1 + r + r * r), x1 = r * x0, x2 = r * x1;
    const DE = 1 + (x0 + x1) / Km, DP = 1 + (x1 + x2) / Km;
    return lr - Math.log(c * DP / DE);
  };
  let lo = -60, hi = 60;
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (lo + hi);
    if (g(m) > 0) hi = m; else lo = m;
    if (hi - lo < 1e-13) break;
  }
  const r = Math.exp(0.5 * (lo + hi));
  const x0 = T / (1 + r + r * r);
  return { x0, x1: r * x0, x2: r * r * x0 };
}

/** Steady state of the HF96 cascade for a given input E1 (uM). */
export function cascadeSteadyState(E1, p = HF96) {
  const MKKKa = oneStepSteadyState(E1, p.E2, p.MKKKtot, p.kcat, p.Km);
  const mkk = twoStepSteadyState(MKKKa, p.MKKPase, p.MKKtot, p.kcat, p.Km);
  const mapk = twoStepSteadyState(mkk.x2, p.MAPKPase, p.MAPKtot, p.kcat, p.Km);
  return { MKKKa, MKK: mkk.x0, MKKP: mkk.x1, MKKPP: mkk.x2, MAPK: mapk.x0, MAPKP: mapk.x1, MAPKPP: mapk.x2 };
}

const TIERS = { MKKK: 'MKKKa', MKK: 'MKKPP', MAPK: 'MAPKPP' };
const SATURATING_E1 = 1; // uM: >> every EC90, the "maximal stimulus" used for normalisation

/** Normalised response (0..1) of one tier ('MKKK' | 'MKK' | 'MAPK') to input E1. */
export function tierResponse(tier, E1, p = HF96) {
  const key = TIERS[tier];
  return cascadeSteadyState(E1, p)[key] / cascadeSteadyState(SATURATING_E1, p)[key];
}

/** Stimulus-response curves on the given input grid (normalised to maximal stimulus). */
export function doseResponse(E1s, p = HF96) {
  const max = cascadeSteadyState(SATURATING_E1, p);
  const out = { E1: E1s.slice(), MKKK: [], MKK: [], MAPK: [] };
  for (const E of E1s) {
    const s = cascadeSteadyState(E, p);
    out.MKKK.push(s.MKKKa / max.MKKKa);
    out.MKK.push(s.MKKPP / max.MKKPP);
    out.MAPK.push(s.MAPKPP / max.MAPKPP);
  }
  return out;
}

/** Input giving a normalised response `level` (bisection in log E1). */
export function effectiveConcentration(f, level, lo = 1e-9, hi = SATURATING_E1) {
  let a = Math.log(lo), b = Math.log(hi);
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (a + b);
    if (f(Math.exp(m)) < level) a = m; else b = m;
    if (b - a < 1e-12) break;
  }
  return Math.exp(0.5 * (a + b));
}

/**
 * Effective Hill coefficient of a monotone stimulus-response curve
 * (Goldbeter & Koshland 1981; the measure used by Huang & Ferrell 1996):
 *   nH = ln(81) / ln(EC90 / EC10)
 * nH = 1 for a Michaelis-Menten (hyperbolic) response.
 */
export function effectiveHill(f, lo, hi) {
  const EC10 = effectiveConcentration(f, 0.1, lo, hi);
  const EC50 = effectiveConcentration(f, 0.5, lo, hi);
  const EC90 = effectiveConcentration(f, 0.9, lo, hi);
  return { EC10, EC50, EC90, nH: Math.log(81) / Math.log(EC90 / EC10) };
}

/** Effective Hill coefficients of all three tiers. */
export function cascadeHill(p = HF96) {
  const out = {};
  for (const tier of Object.keys(TIERS)) out[tier] = effectiveHill((E) => tierResponse(tier, E, p), 1e-9, SATURATING_E1);
  return out;
}

// --------------------------------------------------------------- receptor

/** Illustrative receptor parameters (time in s; ligand in nM). */
export const RTK = {
  Rtot: 1,        // receptor (normalised)
  kon: 0.02,      // nM^-1 s^-1   ligand binding
  koff: 0.01,     // s^-1         (Kd = 0.5 nM)
  kdim: 0.5,      // s^-1 per unit (ligand-bound monomer dimerisation + trans-autophosphorylation)
  kundim: 0.05,   // s^-1
  E1max: 4.5e-5,  // uM Ras-GTP when every receptor is in an active dimer (places the ERK switch near 0.3 nM SCF)
};

/**
 * Steady-state fraction of receptors in active (ligand-bound, dimerised,
 * autophosphorylated) dimers:  R + L <-> LR,  2 LR <-> D,  R + LR + 2D = Rtot.
 */
export function receptorActivation(L, r = RTK) {
  if (L <= 0) return 0;
  const a = (r.kon * L) / r.koff;       // LR = a R
  const b = r.kdim / r.kundim;          // D = b LR^2
  // R (1 + a) + 2 b a^2 R^2 = Rtot   -> positive root of the quadratic
  const A = 2 * b * a * a, B = 1 + a, C = -r.Rtot;
  const R = A > 0 ? (-B + Math.sqrt(B * B - 4 * A * C)) / (2 * A) : r.Rtot / B;
  const LR = a * R, D = b * LR * LR;
  return (2 * D) / r.Rtot;
}

/** Cascade input E1 (Ras-GTP, uM) produced by growth factor L (nM). */
export function inputFromLigand(L, r = RTK) { return r.E1max * receptorActivation(L, r); }

// --------------------------------------------------------------- full dynamic model

export const model = {
  // 0 LR, 1 D (active RTK dimer), 2 Raf*, 3 MEK-P, 4 MEK-PP, 5 ERK-P, 6 ERK-PP (cyto),
  // 7 ERK-PP (nuc), 8 ERK (nuc), 9 Elk-1-P, 10 c-Fos
  species: ['LR', 'RTKdimer', 'MKKKa', 'MKKP', 'MKKPP', 'MAPKP', 'MAPKPP', 'MAPKPPn', 'MAPKn', 'Elk1P', 'Fos'],
  labels: ['SCF:c-Kit', 'c-Kit dimer (pY)', 'Raf*', 'MEK-P', 'MEK-PP', 'ERK-P', 'ERK-PP (cyto)', 'ERK-PP (nuc)', 'ERK (nuc)', 'Elk-1-P', 'c-Fos'],
  defaults: {
    ...HF96, ...RTK,
    L: 1,             // growth factor, nM
    kin: 0.02,        // ERK-PP nuclear import, s^-1 (lumped)
    kexPP: 0.005,     // ERK-PP export
    kmkp: 0.01,       // nuclear MKP-1 inactivation of ERK-PP
    kex: 0.02,        // ERK export
    kElk: 0.05, KElk: 0.05, dElk: 0.01,  // Elk-1 phosphorylation by nuclear ERK-PP (fraction of Elk-1)
    kFos: 0.002, dFos: 0.001,             // c-Fos induction and turnover (normalised)
  },
  initial() { return new Array(11).fill(0); },
  rhs(p) {
    const K = p.Km, k = p.kcat;
    return (t, y, dy) => {
      const [LR, D, Ka, MP, MPP, EP, EPP, EPPn, En, Elk, Fos] = y;
      const R = Math.max(p.Rtot - LR - 2 * D, 0);
      const bind = p.kon * p.L * R - p.koff * LR;
      const dim = p.kdim * LR * LR - p.kundim * D;
      dy[0] = bind - 2 * dim;
      dy[1] = dim;
      const E1 = (p.E1max * 2 * D) / p.Rtot;
      // tier 1
      const K0 = Math.max(p.MKKKtot - Ka, 0);
      const v1 = (k * E1 * K0) / (K + K0), v2 = (k * p.E2 * Ka) / (K + Ka);
      dy[2] = v1 - v2;
      // tier 2 (competitive Michaelis-Menten)
      const M0 = Math.max(p.MKKtot - MP - MPP, 0);
      const dK = 1 + (M0 + MP) / K, dP = 1 + (MP + MPP) / K;
      const v3 = (k * Ka * M0) / K / dK, v4 = (k * Ka * MP) / K / dK;
      const v5 = (k * p.MKKPase * MPP) / K / dP, v6 = (k * p.MKKPase * MP) / K / dP;
      dy[3] = v3 - v4 + v5 - v6;
      dy[4] = v4 - v5;
      // tier 3 (cytoplasmic ERK; nuclear ERK is not accessible to MEK or the cytoplasmic phosphatase)
      const E0 = Math.max(p.MAPKtot - EP - EPP - EPPn - En, 0);
      const eK = 1 + (E0 + EP) / K, eP = 1 + (EP + EPP) / K;
      const v7 = (k * MPP * E0) / K / eK, v8 = (k * MPP * EP) / K / eK;
      const v9 = (k * p.MAPKPase * EPP) / K / eP, v10 = (k * p.MAPKPase * EP) / K / eP;
      const imp = p.kin * EPP - p.kexPP * EPPn;
      dy[5] = v7 - v8 + v9 - v10;
      dy[6] = v8 - v9 - imp;
      dy[7] = imp - p.kmkp * EPPn;
      dy[8] = p.kmkp * EPPn - p.kex * En;
      // transcription-factor activation and immediate-early gene
      dy[9] = (p.kElk * EPPn * (1 - Elk)) / (p.KElk + EPPn) - p.dElk * Elk;
      dy[10] = p.kFos * Elk - p.dFos * Fos;
    };
  },
  /** Conserved totals (receptor and ERK) for checking an integration. */
  totals(y, p) {
    const [LR, D, , , , EP, EPP, EPPn, En] = y;
    return { RTK: LR + 2 * D, ERKphosOrNuclear: EP + EPP + EPPn + En, ERKtot: p.MAPKtot };
  },
};

/**
 * Normalised ERK activity (0..1) at steady state for growth factor L: the
 * cascade steady state for E1(L) relative to maximal stimulation. With linear
 * import, nuclear ERK-PP is proportional to this cytoplasmic level, so it is
 * used as the drive of nuclear targets (proliferation genes) in the scenes.
 */
export function nuclearERK(L, p = model.defaults) {
  const s = cascadeSteadyState(inputFromLigand(L, p), p);
  const smax = cascadeSteadyState(SATURATING_E1, p);
  return s.MAPKPP / smax.MAPKPP;
}

/** Reaction list for the network view (order 2 = the paper's progenitor stage). */
export const REACTIONS = [
  { id: 'mapk.bind', order: 2, pathway: 'MAPK/ERK', reactants: ['SCF', 'c-Kit'], products: ['SCF:c-Kit'], label: 'Growth factor binds receptor tyrosine kinase' },
  { id: 'mapk.dimer', order: 2, pathway: 'MAPK/ERK', reactants: ['SCF:c-Kit', 'SCF:c-Kit'], products: ['c-Kit dimer (pY)'], label: 'RTK dimerisation and trans-autophosphorylation' },
  { id: 'mapk.ras', order: 2, pathway: 'MAPK/ERK', reactants: ['Ras-GDP'], products: ['Ras-GTP'], modifiers: ['c-Kit dimer (pY)', 'Grb2', 'SOS'], label: 'Ras activation (Grb2/SOS)' },
  { id: 'mapk.raf', order: 2, pathway: 'MAPK/ERK', reactants: ['Raf'], products: ['Raf*'], modifiers: ['Ras-GTP'], label: 'MAPKKK activation (E1)' },
  { id: 'mapk.raf.off', order: 2, pathway: 'MAPK/ERK', reactants: ['Raf*'], products: ['Raf'], modifiers: ['Raf phosphatase'], label: 'MAPKKK inactivation (E2)' },
  { id: 'mapk.mek.p1', order: 2, pathway: 'MAPK/ERK', reactants: ['MEK'], products: ['MEK-P'], modifiers: ['Raf*'], label: 'MEK first phosphorylation' },
  { id: 'mapk.mek.p2', order: 2, pathway: 'MAPK/ERK', reactants: ['MEK-P'], products: ['MEK-PP'], modifiers: ['Raf*'], label: 'MEK second phosphorylation' },
  { id: 'mapk.mek.d2', order: 2, pathway: 'MAPK/ERK', reactants: ['MEK-PP'], products: ['MEK-P'], modifiers: ['MEK phosphatase'], label: 'MEK-PP dephosphorylation' },
  { id: 'mapk.mek.d1', order: 2, pathway: 'MAPK/ERK', reactants: ['MEK-P'], products: ['MEK'], modifiers: ['MEK phosphatase'], label: 'MEK-P dephosphorylation' },
  { id: 'mapk.erk.p1', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK'], products: ['ERK-P'], modifiers: ['MEK-PP'], label: 'ERK first phosphorylation' },
  { id: 'mapk.erk.p2', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK-P'], products: ['ERK-PP'], modifiers: ['MEK-PP'], label: 'ERK second phosphorylation' },
  { id: 'mapk.erk.d2', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK-PP'], products: ['ERK-P'], modifiers: ['ERK phosphatase'], label: 'ERK-PP dephosphorylation' },
  { id: 'mapk.erk.d1', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK-P'], products: ['ERK'], modifiers: ['ERK phosphatase'], label: 'ERK-P dephosphorylation' },
  { id: 'mapk.imp', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK-PP'], products: ['ERK-PP (nuc)'], label: 'ERK nuclear translocation' },
  { id: 'mapk.mkp', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK-PP (nuc)'], products: ['ERK (nuc)'], modifiers: ['MKP-1'], label: 'Nuclear ERK inactivation' },
  { id: 'mapk.exp', order: 2, pathway: 'MAPK/ERK', reactants: ['ERK (nuc)'], products: ['ERK'], label: 'ERK nuclear export' },
  { id: 'mapk.elk', order: 2, pathway: 'MAPK/ERK', reactants: ['Elk-1'], products: ['Elk-1-P'], modifiers: ['ERK-PP (nuc)'], label: 'Transcription-factor activation (Elk-1)' },
  { id: 'mapk.fos', order: 2, pathway: 'MAPK/ERK', reactants: [], products: ['c-Fos'], modifiers: ['Elk-1-P'], label: 'Immediate-early gene induction' },
];
