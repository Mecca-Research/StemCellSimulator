// Hematopoiesis (paper: "Stage 2 - Progenitor cell differentiation", worked
// example "Hematopoietic Stem Cell Differentiation": HSC -> HSC + HSC;
// HSC --(GATA1, PU.1)--> CMP; CMP --(EPO)--> proerythroblast --(EPO, GATA1)-->
// erythroblast -> reticulocyte -> erythrocyte; lineage graph HSC -> CMP, CLP;
// CMP -> red cells, platelets, neutrophils, monocytes; CLP -> B, T cells).
//
// (i)  toggle: the GATA1-PU.1 fate switch with mutual inhibition and
//      self-activation (Huang, Guo, May & Enver 2007, Dev Biol 305:695):
//        dx1/dt = a1 x1^n/(thA^n + x1^n) + b1 thB^n/(thB^n + x2^n) - k1 x1 + uG
//        dx2/dt = a2 x2^n/(thA^n + x2^n) + b2 thB^n/(thB^n + x1^n) - k2 x2 + uP
//      x1 = GATA1, x2 = PU.1, n = 4, thA = thB = 0.5, k = 1. With strong
//      self-activation (a = b = 1) the system is tristable: a central
//      "primed" attractor (low-level co-expression, as in multipotent
//      progenitors) flanked by two committed attractors. Commitment is the
//      decay of a, which destroys the central attractor in a pitchfork
//      bifurcation and leaves a bistable switch. External signals enter as
//      additive production terms, as in Chickarmane, Enver & Peterson 2009,
//      PLoS Comput Biol 5:e1000268 (EPO -> GATA1, GM-CSF -> PU.1). Noise:
//      additive Langevin term integrated by Euler-Maruyama.
// (ii) LINEAGE: the lineage tree as data (names, markers, morphology keys,
//      sizes, stage durations, divisions) and TRANSITIONS (the paper's
//      reactions).
// (iii) LineageSSA: exact stochastic simulation (Gillespie 1977) of the
//      compartment model; HSC self-renewal is niche-limited (Schofield 1978,
//      Blood Cells 4:7): each HSC daughter stays a stem cell with probability
//      a (1 - H/K), the "fraction of self-renewal" of Marciniak-Czochra, Stiehl,
//      Ho, Jaeger & Wagner 2009, Stem Cells Dev 18:377.

// ============================================================= (i) toggle

export const toggle = {
  defaults: { a1: 1, a2: 1, b1: 1, b2: 1, thA: 0.5, thB: 0.5, n: 4, k1: 1, k2: 1, uG: 0, uP: 0 },

  /** Deterministic rates into out = [dGATA1, dPU1]. */
  rates(g, u, p, out = [0, 0]) {
    const n = p.n, tA = p.thA ** n, tB = p.thB ** n;
    const gn = Math.max(g, 0) ** n, un = Math.max(u, 0) ** n;
    out[0] = (p.a1 * gn) / (tA + gn) + (p.b1 * tB) / (tB + un) - p.k1 * g + p.uG;
    out[1] = (p.a2 * un) / (tA + un) + (p.b2 * tB) / (tB + gn) - p.k2 * u + p.uP;
    return out;
  },

  /** ODE right-hand side in the core ode.js convention (y = [GATA1, PU.1]). */
  rhs(p) {
    const out = [0, 0];
    return (t, y, dy) => { toggle.rates(y[0], y[1], p, out); dy[0] = out[0]; dy[1] = out[1]; };
  },

  /** Analytic Jacobian [[dg/dg, dg/du], [du/dg, du/du]]. */
  jacobian(g, u, p) {
    const n = p.n, tA = p.thA ** n, tB = p.thB ** n;
    const act = (x) => { x = Math.max(x, 1e-12); const xn = x ** n; return (n * tA * xn / x) / (tA + xn) ** 2; };
    const rep = (x) => { x = Math.max(x, 1e-12); const xn = x ** n; return (-n * tB * xn / x) / (tB + xn) ** 2; };
    return [[p.a1 * act(g) - p.k1, p.b1 * rep(u)], [p.b2 * rep(g), p.a2 * act(u) - p.k2]];
  },

  /**
   * All fixed points in [0, xmax]^2 by Newton iteration from a grid of seeds,
   * classified by the Jacobian (stable node/focus, saddle, unstable).
   */
  fixedPoints(p, { xmax = 4, grid = 14 } = {}) {
    const found = [];
    const r = [0, 0];
    for (let i = 0; i <= grid; i++) for (let j = 0; j <= grid; j++) {
      let g = (xmax * i) / grid + 1e-3, u = (xmax * j) / grid + 1e-3;
      let ok = false;
      for (let it = 0; it < 100; it++) {
        toggle.rates(g, u, p, r);
        const J = toggle.jacobian(g, u, p);
        const det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (Math.abs(det) < 1e-14) break;
        const dg = (J[1][1] * r[0] - J[0][1] * r[1]) / det;
        const du = (-J[1][0] * r[0] + J[0][0] * r[1]) / det;
        const s = Math.min(1, 0.5 / Math.max(Math.abs(dg), Math.abs(du), 1e-12));
        g -= s * dg; u -= s * du;
        if (g < 0 || u < 0 || g > 2 * xmax || u > 2 * xmax) break;
        if (Math.abs(dg) + Math.abs(du) < 1e-12) { ok = true; break; }
      }
      if (!ok) continue;
      toggle.rates(g, u, p, r);
      if (Math.abs(r[0]) + Math.abs(r[1]) > 1e-9) continue;
      if (found.some((f) => Math.abs(f.g - g) + Math.abs(f.u - u) < 1e-6)) continue;
      const J = toggle.jacobian(g, u, p);
      const tr = J[0][0] + J[1][1], det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
      const type = det < 0 ? 'saddle' : tr < 0 ? 'stable' : 'unstable';
      found.push({ g, u, type, stable: type === 'stable' });
    }
    return found.sort((a, b) => a.g - b.g);
  },

  /**
   * Self-activation strength below which the symmetric central (primed)
   * attractor is lost (pitchfork), for a1 = a2 = a and symmetric parameters.
   */
  criticalSelfActivation(p = toggle.defaults) {
    const hasCentre = (a) => toggle.fixedPoints({ ...p, a1: a, a2: a, uG: 0, uP: 0 })
      .some((f) => f.stable && Math.abs(f.g - f.u) < 1e-4);
    let lo = 0, hi = 1.5;
    if (!hasCentre(hi)) return NaN;
    for (let i = 0; i < 40; i++) { const m = 0.5 * (lo + hi); if (hasCentre(m)) hi = m; else lo = m; }
    return 0.5 * (lo + hi);
  },

  /**
   * Euler-Maruyama step of the Langevin equation
   *   dx = f(x) dt + sigma dW   (reflecting at 0)
   * on a state object {g, u}; `a` overrides a1 = a2 (the self-activation that
   * decays during commitment).
   */
  langevinStep(s, dt, p, sigma, rng, a = null, scratch = [0, 0]) {
    const q = a === null ? p : withA(p, a);
    toggle.rates(s.g, s.u, q, scratch);
    const sd = sigma * Math.sqrt(dt);
    s.g = Math.abs(s.g + scratch[0] * dt + sd * rng.normal());
    s.u = Math.abs(s.u + scratch[1] * dt + sd * rng.normal());
    return s;
  },
};

const _pa = { ...toggle.defaults };
function withA(p, a) { Object.assign(_pa, p); _pa.a1 = a; _pa.a2 = a; return _pa; }

/**
 * Commitment protocol of Huang et al. 2007: cells start in the primed state
 * with self-activation a0, a decays as a0 exp(-t/tauA), noise sigma; the fate
 * is read at time T (GATA1 > PU.1 -> MEP / erythroid-megakaryocytic,
 * otherwise GMP / myelomonocytic). Time in units of 1/k (TF turnover).
 */
export function commitCell(s, p, { sigma = 0.1, a0 = 1, tauA = 2, T = 8, dt = 0.05, rng }) {
  const scratch = [0, 0];
  for (let t = 0; t < T; t += dt) toggle.langevinStep(s, dt, p, sigma, rng, a0 * Math.exp(-t / tauA), scratch);
  return s;
}

/** Primed starting state: the central attractor at a = a0 (with the current inputs). */
export function primedState(p, a0 = 1) {
  const fps = toggle.fixedPoints({ ...p, a1: a0, a2: a0 });
  let best = null;
  for (const f of fps) if (f.stable && (!best || Math.abs(f.g - f.u) < Math.abs(best.g - best.u))) best = f;
  return best ? { g: best.g, u: best.u } : { g: 1, u: 1 };
}

/** Fraction of CMPs that commit to the GATA1-high (MEP) fate. */
export function fateProbability(p, { n = 400, rng, ...opts }) {
  const start = primedState({ ...p, uG: 0, uP: 0 }, opts.a0 ?? 1);
  let ery = 0;
  for (let i = 0; i < n; i++) {
    const s = commitCell({ g: start.g, u: start.u }, p, { rng, ...opts });
    if (s.g > s.u) ery++;
  }
  return ery / n;
}

// Downstream decision rules read from the committed toggle state.
// GMP: PU.1 dosage decides macrophage/monocyte (high) vs neutrophil (low)
// (Dahl et al. 2003, Nat Immunol 4:1029); residual GATA1 in GMPs licenses the
// eosinophil/basophil fates (Iwasaki et al. 2006, Genes Dev 20:3010). The
// sigmoidal read-outs are phenomenological.
export const FATE_RULES = {
  monoK: 1.35, monoN: 6, eosBasoMax: 0.25, eosBasoK: 0.3,
  pMK: 0.15,                                   // MEP -> megakaryocyte (display; see LINEAGE.MK.plateletYield)
  clp: { B: 0.65, T: 0.2, NK: 0.15 },          // CLP output in marrow (T = thymus-seeding precursors)
};

export function gmpFate(g, u, rng, R = FATE_RULES) {
  const pEB = R.eosBasoMax * (g * g) / (R.eosBasoK ** 2 + g * g);
  const x = rng.next();
  if (x < pEB) return rng.next() < 0.75 ? 'Eos' : 'Baso';
  const um = u ** R.monoN;
  return rng.next() < um / (R.monoK ** R.monoN + um) ? 'Mono' : 'Neu';
}

export function mepFate(rng, R = FATE_RULES) { return rng.next() < R.pMK ? 'MK' : 'ProEB'; }

export function clpFate(rng, R = FATE_RULES) {
  const x = rng.next();
  return x < R.clp.B ? 'B' : x < R.clp.B + R.clp.T ? 'T' : 'NK';
}

/**
 * Monte-Carlo distribution over myeloid fates for the current signal inputs:
 * CMP commitment (toggle) followed by the GMP read-out; the GMP state is the
 * committed state relaxed with a = 0 for one more turnover time.
 */
export function fateDistribution(p, { n = 300, rng, ...opts }) {
  const start = primedState({ ...p, uG: 0, uP: 0 }, opts.a0 ?? 1);
  const out = { MEP: 0, GMP: 0, Neu: 0, Mono: 0, Eos: 0, Baso: 0 };
  for (let i = 0; i < n; i++) {
    const s = commitCell({ g: start.g, u: start.u }, p, { rng, ...opts });
    if (s.g > s.u) { out.MEP++; continue; }
    out.GMP++;
    out[gmpFate(s.g, s.u, rng)]++;
  }
  for (const k of Object.keys(out)) out[k] /= n;
  const gm = out.GMP || 1;
  out.gmp = { Neu: out.Neu / gm, Mono: out.Mono / gm, Eos: out.Eos / gm, Baso: out.Baso / gm };
  return out;
}

// ============================================================= (ii) lineage

/**
 * Lineage tree. Sizes are typical diameters (um) from standard hematology
 * atlases; stageDays / divisions are the marrow residence and number of
 * divisions of each stage (e.g. one proerythroblast -> 16 reticulocytes in
 * 4 divisions over ~4-5 days; neutrophil maturation + marrow storage ~6-7
 * days, Dancey et al. 1976 J Clin Invest 58:705); lifespanDays for terminal
 * cells after release (erythrocyte ~120 d, platelet 7-10 d, neutrophil < 1 d
 * in blood). The HSC cycle is
 * compressed for display (dormant murine HSCs divide every ~145 d; Wilson et
 * al. 2008, Cell 135:1118).
 */
export const LINEAGE = {
  HSC: { name: 'Hematopoietic stem cell', group: 'stem', markers: 'Lin− CD34+ CD38− CD90+ CD45RA− CD49f+', morph: 'blast', diameter: 8, nc: 0.72, cycleDays: 6, children: ['CMP', 'CLP'] },
  CMP: { name: 'Common myeloid progenitor', group: 'progenitor', markers: 'CD34+ CD38+ CD123lo CD45RA−', morph: 'blast', diameter: 10, nc: 0.7, stageDays: 2, divisions: 2, children: ['MEP', 'GMP'], kitDependent: 1 },
  CLP: { name: 'Common lymphoid progenitor', group: 'progenitor', markers: 'CD34+ CD38+ CD10+ CD7+ IL-7Rα+', morph: 'blast', diameter: 9, nc: 0.75, stageDays: 2, divisions: 2, children: ['B', 'T', 'NK'] },
  MEP: { name: 'Megakaryocyte–erythroid progenitor', group: 'progenitor', markers: 'CD34+ CD38+ CD123− CD45RA−', morph: 'blast', diameter: 10, nc: 0.7, stageDays: 1.5, divisions: 1, children: ['ProEB', 'MK'], kitDependent: 1 },
  GMP: { name: 'Granulocyte–macrophage progenitor', group: 'progenitor', markers: 'CD34+ CD38+ CD123+ CD45RA+', morph: 'blast', diameter: 11, nc: 0.68, stageDays: 2, divisions: 2, children: ['Neu', 'Mono', 'Eos', 'Baso'], kitDependent: 1 },
  ProEB: { name: 'Proerythroblast', group: 'erythroid', markers: 'CD71hi CD117+ GlyA lo', morph: 'erythroblast', diameter: 16, nc: 0.72, stageDays: 1, divisions: 1, children: ['EB'], epoDependent: 1, kitDependent: 1 },
  EB: { name: 'Erythroblast (baso → poly → orthochromatic)', group: 'erythroid', markers: 'CD71+ GlyA+', morph: 'erythroblast', diameter: 12, diameterEnd: 8.5, nc: 0.6, ncEnd: 0.42, stageDays: 3, divisions: 3, children: ['Retic'], epoDependent: 0.25 },
  Retic: { name: 'Reticulocyte', group: 'erythroid', markers: 'GlyA+ CD71lo, residual RNA (anucleate)', morph: 'reticulocyte', diameter: 8.4, stageDays: 1, divisions: 0, children: ['RBC'] },
  RBC: { name: 'Erythrocyte', group: 'erythroid', markers: 'GlyA+ CD71− (anucleate biconcave disc)', morph: 'erythrocyte', diameter: 7.82, lifespanDays: 120 },
  MK: { name: 'Megakaryocyte', group: 'megakaryocytic', markers: 'CD41+ CD42b+ CD61+ (polyploid 16–64N)', morph: 'megakaryocyte', diameter: 40, stageDays: 4, divisions: 0, children: ['PLT'], plateletYield: 45 },
  PLT: { name: 'Platelet', group: 'megakaryocytic', markers: 'CD41+ CD42b+ (anucleate)', morph: 'platelet', diameter: 2.5, lifespanDays: 9 },
  Neu: { name: 'Neutrophil', group: 'myeloid', markers: 'CD15+ CD16+ CD66b+', morph: 'neutrophil', diameter: 13, nc: 0.5, stageDays: 6.5, divisions: 1, lifespanDays: 1 },
  Mono: { name: 'Monocyte', group: 'myeloid', markers: 'CD14+ CD64+ CX3CR1+', morph: 'monocyte', diameter: 16, nc: 0.55, stageDays: 2, divisions: 1, lifespanDays: 3 },
  Eos: { name: 'Eosinophil', group: 'myeloid', markers: 'Siglec-8+ CCR3+', morph: 'eosinophil', diameter: 14, nc: 0.45, stageDays: 4, divisions: 1, lifespanDays: 6 },
  Baso: { name: 'Basophil', group: 'myeloid', markers: 'FcεRIα+ CD123+ CD203c+', morph: 'basophil', diameter: 12, nc: 0.5, stageDays: 3, divisions: 1, lifespanDays: 3 },
  B: { name: 'B lymphocyte', group: 'lymphoid', markers: 'CD19+ CD20+ (IgM+)', morph: 'lymphocyte', diameter: 10, nc: 0.9, stageDays: 5, divisions: 1, lifespanDays: 20 },
  T: { name: 'T-cell precursor (seeds the thymus)', group: 'lymphoid', markers: 'CD34+ CD7+ CCR9+', morph: 'lymphocyte', diameter: 9.5, nc: 0.9, stageDays: 1, divisions: 0, lifespanDays: 1 },
  NK: { name: 'NK cell', group: 'lymphoid', markers: 'CD56+ CD16+ CD3−', morph: 'lymphocyte', diameter: 11, nc: 0.82, stageDays: 3, divisions: 1, lifespanDays: 10 },
};

export const NODE_IDS = Object.keys(LINEAGE);

/** The paper's transitions (and the refinements through MEP/GMP) as data. */
export const TRANSITIONS = [
  { id: 'hsc.renew', from: 'HSC', to: ['HSC', 'HSC'], kind: 'self-renewal', modifiers: ['SCF', 'TPO', 'CXCL12'], label: 'Stem cell self-renewal HSC → HSC + HSC' },
  { id: 'hsc.cmp', from: 'HSC', to: ['CMP'], kind: 'differentiation', modifiers: ['GATA1', 'PU.1'], label: 'HSC → common myeloid progenitor' },
  { id: 'hsc.clp', from: 'HSC', to: ['CLP'], kind: 'differentiation', modifiers: ['IL-7', 'Ikaros'], label: 'HSC → common lymphoid progenitor' },
  { id: 'cmp.mep', from: 'CMP', to: ['MEP'], kind: 'fate choice', modifiers: ['GATA1', 'EPO'], label: 'CMP → MEP (GATA1-high)' },
  { id: 'cmp.gmp', from: 'CMP', to: ['GMP'], kind: 'fate choice', modifiers: ['PU.1', 'GM-CSF'], label: 'CMP → GMP (PU.1-high)' },
  { id: 'mep.proeb', from: 'MEP', to: ['ProEB'], kind: 'differentiation', modifiers: ['EPO', 'GATA1', 'KLF1'], label: 'CMP/MEP → proerythroblast (EPO)' },
  { id: 'proeb.eb', from: 'ProEB', to: ['EB'], kind: 'maturation', modifiers: ['EPO', 'GATA1'], label: 'Proerythroblast → erythroblast (EPO, GATA1)' },
  { id: 'eb.retic', from: 'EB', to: ['Retic', 'Pyrenocyte'], kind: 'enucleation', modifiers: ['Rac GTPases', 'mDia2'], label: 'Erythroblast → reticulocyte (enucleation)' },
  { id: 'retic.rbc', from: 'Retic', to: ['RBC'], kind: 'maturation', modifiers: [], label: 'Reticulocyte → erythrocyte' },
  { id: 'mep.mk', from: 'MEP', to: ['MK'], kind: 'differentiation', modifiers: ['TPO', 'FLI1'], label: 'MEP → megakaryocyte (endomitosis)' },
  { id: 'mk.plt', from: 'MK', to: ['PLT'], kind: 'fragmentation', modifiers: ['TPO'], label: 'Megakaryocyte → platelets (proplatelets)' },
  { id: 'gmp.neu', from: 'GMP', to: ['Neu'], kind: 'differentiation', modifiers: ['G-CSF', 'C/EBPα'], label: 'GMP → neutrophil' },
  { id: 'gmp.mono', from: 'GMP', to: ['Mono'], kind: 'differentiation', modifiers: ['M-CSF', 'PU.1'], label: 'GMP → monocyte' },
  { id: 'gmp.eos', from: 'GMP', to: ['Eos'], kind: 'differentiation', modifiers: ['IL-5', 'GATA1'], label: 'GMP → eosinophil' },
  { id: 'gmp.baso', from: 'GMP', to: ['Baso'], kind: 'differentiation', modifiers: ['IL-3', 'GATA2'], label: 'GMP → basophil' },
  { id: 'clp.b', from: 'CLP', to: ['B'], kind: 'differentiation', modifiers: ['IL-7', 'PAX5'], label: 'CLP → B cell' },
  { id: 'clp.t', from: 'CLP', to: ['T'], kind: 'differentiation', modifiers: ['Notch1'], label: 'CLP → T-cell precursor' },
  { id: 'clp.nk', from: 'CLP', to: ['NK'], kind: 'differentiation', modifiers: ['IL-15'], label: 'CLP → NK cell' },
];

/**
 * Exponential-waiting-time rates that reproduce a deterministic stage of
 * length D with k evenly spaced divisions (at (i + 1/2) D / k): the same
 * amplification A = 2^k and the same mean cell-days Q per entering cell.
 *   m (exit) = A / Q,  p (division) = (A - 1) / Q
 */
export function stageRates(node) {
  const D = node.stageDays, k = node.divisions ?? 0;
  const A = 2 ** k;
  let Q = D;
  if (k > 0) {
    let s = 0.5 + 0.5 * A;
    for (let j = 1; j < k; j++) s += 2 ** j;
    Q = (D / k) * s;
  }
  return { exit: A / Q, div: (A - 1) / Q, A, Q };
}

// ============================================================= (iii) population SSA

/**
 * Environment of the compartment model (all dimensionless drives in 0..1):
 *   selfRenewal  a, niche-limited: P(daughter stays HSC) = a (1 - H/K)
 *   K            number of HSC niches in the modelled marrow
 *   fCLP         lymphoid fraction of HSC output
 *   pMEP         CMP -> MEP probability (from the GATA1-PU.1 toggle)
 *   gmp          {Neu, Mono, Eos, Baso} split of GMP output (toggle read-out)
 *   erk          MAPK/ERK drive (SCF): speeds up progenitor cycling
 *   stat5        JAK/STAT5 drive (EPO): erythroid progenitor survival
 */
export const SSA_DEFAULTS = {
  selfRenewal: 0.8, K: 200, fCLP: 0.2, pMEP: 0.55, pMK: FATE_RULES.pMK,
  gmp: { Neu: 0.7, Mono: 0.24, Eos: 0.045, Baso: 0.015 }, clp: FATE_RULES.clp,
  erk: 1, stat5: 0.5, apoEry: 1.5, apoEryK: 0.15, apoKit: 1,
};

/** Cycle-speed factor from ERK activity (phenomenological, 0.35..1). */
export const erkSpeed = (erk) => 0.35 + 0.65 * Math.min(Math.max(erk, 0), 1);
/** EPO-withdrawal apoptosis rate (per day) of EPO-dependent erythroid stages (phenomenological). */
export const eryApoptosis = (stat5, env = SSA_DEFAULTS) => (env.apoEry * env.apoEryK) / (env.apoEryK + Math.max(stat5, 0));
/**
 * SCF-withdrawal apoptosis rate (per day) of c-Kit+ progenitors (phenomenological): SCF/c-Kit
 * is a survival as well as a proliferation signal (Steel / W mutant mice are anaemic).
 */
export const kitApoptosis = (erk, env = SSA_DEFAULTS) => env.apoKit * (1 - Math.min(Math.max(erk, 0), 1)) ** 2;

/** Total apoptosis hazard (per day) of a node under the current signals. */
export function apoptosisRate(id, env = SSA_DEFAULTS) {
  const node = LINEAGE[id];
  return (node.epoDependent ?? 0) * eryApoptosis(env.stat5, env) + (node.kitDependent ?? 0) * kitApoptosis(env.erk, env);
}

const PROGENITOR_SPEED = new Set(['CMP', 'CLP', 'MEP', 'GMP', 'ProEB', 'EB']);

/** Per-node rates (per day) for the current environment. */
export function nodeRates(id, env) {
  const node = LINEAGE[id];
  if (id === 'HSC') return { div: 1 / node.cycleDays, exit: 0, death: 0 };
  if (node.stageDays === undefined) return { div: 0, exit: 0, death: 1 / node.lifespanDays };
  if (node.lifespanDays !== undefined) {
    // terminal leukocytes: marrow maturation (with divisions) + circulation, one compartment
    const r = stageRates(node);
    const Q = r.Q + r.A * node.lifespanDays;
    return { div: (r.A - 1) / Q, exit: 0, death: r.A / Q };
  }
  const r = stageRates(node);
  const f = PROGENITOR_SPEED.has(id) ? erkSpeed(env.erk) : 1;
  return { div: r.div * f, exit: r.exit * f, death: apoptosisRate(id, env) };
}

function childWeights(id, env) {
  switch (id) {
    case 'HSC': return [['CMP', 1 - env.fCLP], ['CLP', env.fCLP]];
    case 'CMP': return [['MEP', env.pMEP], ['GMP', 1 - env.pMEP]];
    case 'MEP': return [['ProEB', 1 - env.pMK], ['MK', env.pMK]];
    case 'GMP': return Object.entries(env.gmp);
    case 'CLP': return Object.entries(env.clp);
    default: return (LINEAGE[id].children ?? []).map((c) => [c, 1]);
  }
}

/**
 * Mean-field steady state of the compartment model (analytic): the HSC
 * balance 2 a (1 - H/K) = 1 gives H* = K (1 - 1/(2a)); every downstream
 * compartment is linear in its inflow J: N = J / (exit + death - div).
 */
export function meanFieldSteadyState(env = SSA_DEFAULTS) {
  const N = Object.fromEntries(NODE_IDS.map((id) => [id, 0]));
  const H = env.selfRenewal > 0.5 ? env.K * (1 - 1 / (2 * env.selfRenewal)) : 0;
  N.HSC = H;
  const inflow = Object.fromEntries(NODE_IDS.map((id) => [id, 0]));
  const hscOut = H * nodeRates('HSC', env).div; // 2 (1 - a s*) p H with 2 a s* = 1
  for (const [c, w] of childWeights('HSC', env)) inflow[c] += hscOut * w;
  const order = ['CMP', 'CLP', 'MEP', 'GMP', 'ProEB', 'EB', 'Retic', 'MK', 'RBC', 'PLT', 'Neu', 'Mono', 'Eos', 'Baso', 'B', 'T', 'NK'];
  for (const id of order) {
    const r = nodeRates(id, env);
    const net = r.exit + r.death - r.div;
    N[id] = net > 0 ? inflow[id] / net : Infinity;
    const out = r.exit * N[id];
    if (out > 0) {
      const ws = childWeights(id, env);
      const tot = ws.reduce((s, [, w]) => s + w, 0);
      for (const [c, w] of ws) inflow[c] += (out * w / tot) * (c === 'PLT' ? LINEAGE.MK.plateletYield : 1);
    }
  }
  return N;
}

/**
 * Mean number of cells of a node that are still inside the marrow (for
 * leukocytes the compartment also counts circulating cells): the fraction
 * Q / (Q + A * lifespan) of the compartment.
 */
export function marrowFraction(id) {
  const node = LINEAGE[id];
  if (node.stageDays === undefined) return 0;
  if (node.lifespanDays === undefined) return 1;
  const r = stageRates(node);
  return r.Q / (r.Q + r.A * node.lifespanDays);
}

export class LineageSSA {
  /**
   * @param {object} o
   * @param {import('../core/rng.js').RNG} o.rng
   * @param {object} [o.env]  see SSA_DEFAULTS
   * @param {object} [o.init] initial counts by node id (default: rounded mean-field steady state)
   */
  constructor({ rng, env = {}, init = null }) {
    this.rng = rng;
    this.env = { ...SSA_DEFAULTS, ...env };
    this.ids = NODE_IDS;
    this.index = Object.fromEntries(this.ids.map((id, i) => [id, i]));
    this.x = new Float64Array(this.ids.length);
    const start = init ?? meanFieldSteadyState(this.env);
    for (const id of this.ids) this.x[this.index[id]] = Math.max(0, Math.round(Number.isFinite(start[id] ?? 0) ? start[id] ?? 0 : 0));
    this.t = 0;
    this.events = 0;
    this.setEnv({});
  }

  setEnv(patch) {
    Object.assign(this.env, patch);
    const env = this.env;
    // flatten into channel tables: per node [div, exit, death] rates and child CDFs
    this.rate = this.ids.map((id) => nodeRates(id, env));
    this.children = this.ids.map((id) => {
      const ws = childWeights(id, env);
      const tot = ws.reduce((s, [, w]) => s + w, 0) || 1;
      let acc = 0;
      return ws.map(([c, w]) => { acc += w / tot; return [this.index[c], acc]; });
    });
    this.props = new Float64Array(this.ids.length * 3);
  }

  count(id) { return this.x[this.index[id]]; }

  pickChild(i) {
    const ch = this.children[i];
    const r = this.rng.next();
    for (let k = 0; k < ch.length; k++) if (r < ch[k][1]) return ch[k][0];
    return ch[ch.length - 1][0];
  }

  /** Advance by dt days (exact SSA; stops early after maxEvents and keeps the clock consistent). */
  step(dt, maxEvents = 50000) {
    const tEnd = this.t + dt;
    const x = this.x, n = this.ids.length, props = this.props, env = this.env;
    const iH = this.index.HSC, iPLT = this.index.PLT, iMK = this.index.MK;
    let ev = 0;
    while (ev < maxEvents) {
      let a0 = 0;
      for (let i = 0; i < n; i++) {
        const r = this.rate[i], N = x[i];
        const o = i * 3;
        props[o] = r.div * N; props[o + 1] = r.exit * N; props[o + 2] = r.death * N;
        a0 += props[o] + props[o + 1] + props[o + 2];
      }
      if (a0 <= 0) { this.t = tEnd; break; }
      const tau = this.rng.exponential(a0);
      if (this.t + tau > tEnd) { this.t = tEnd; break; }
      this.t += tau;
      let u = this.rng.next() * a0, ch = 0;
      for (; ch < props.length - 1; ch++) { u -= props[ch]; if (u < 0) break; }
      const i = (ch / 3) | 0, kind = ch % 3;
      if (i === iH && kind === 0) {
        // HSC division: each daughter stays a stem cell with prob a (1 - H/K)
        const s = env.selfRenewal * Math.max(0, 1 - x[iH] / env.K);
        x[iH] -= 1;
        for (let d = 0; d < 2; d++) {
          if (this.rng.next() < s) x[iH] += 1;
          else x[this.pickChild(iH)] += 1;
        }
      } else if (kind === 0) {
        x[i] += 1;
      } else if (kind === 1) {
        x[i] -= 1;
        if (i === iMK) x[iPLT] += LINEAGE.MK.plateletYield;
        else x[this.pickChild(i)] += 1;
      } else {
        x[i] -= 1;
      }
      ev++;
    }
    this.events += ev;
    if (ev >= maxEvents) this.t = tEnd; // saturated: accept the truncated interval
    return ev;
  }

  /** Grouped counts for charts. */
  groups() {
    const c = (id) => this.count(id);
    return {
      HSC: c('HSC'),
      progenitors: c('CMP') + c('CLP') + c('MEP') + c('GMP'),
      erythroid: c('ProEB') + c('EB') + c('Retic'),
      RBC: c('RBC'),
      MK: c('MK'),
      PLT: c('PLT'),
      myeloid: c('Neu') + c('Mono') + c('Eos') + c('Baso'),
      lymphoid: c('B') + c('T') + c('NK'),
    };
  }
}

/** Lineage transitions as order-2 reactions for the network view. */
export const REACTIONS = [
  ...TRANSITIONS.map((tr) => ({
    id: `hema.${tr.id}`, order: 2, pathway: 'Hematopoiesis',
    reactants: [LINEAGE[tr.from].name],
    products: tr.to.map((t) => (t === 'Pyrenocyte' ? 'Pyrenocyte (extruded nucleus)' : LINEAGE[t].name)),
    modifiers: tr.modifiers, label: tr.label,
  })),
  { id: 'hema.gata1.auto', order: 2, pathway: 'Hematopoiesis', reactants: [], products: ['GATA1'], modifiers: ['GATA1'], label: 'GATA1 self-activation' },
  { id: 'hema.pu1.auto', order: 2, pathway: 'Hematopoiesis', reactants: [], products: ['PU.1'], modifiers: ['PU.1'], label: 'PU.1 self-activation' },
  { id: 'hema.gata1.rep', order: 2, pathway: 'Hematopoiesis', reactants: [], products: ['PU.1'], modifiers: ['GATA1'], label: 'GATA1 represses PU.1 (mutual inhibition)' },
  { id: 'hema.pu1.rep', order: 2, pathway: 'Hematopoiesis', reactants: [], products: ['GATA1'], modifiers: ['PU.1'], label: 'PU.1 represses GATA1 (mutual inhibition)' },
  { id: 'hema.gata1.deg', order: 2, pathway: 'Hematopoiesis', reactants: ['GATA1'], products: [], label: 'GATA1 turnover' },
  { id: 'hema.pu1.deg', order: 2, pathway: 'Hematopoiesis', reactants: ['PU.1'], products: [], label: 'PU.1 turnover' },
];
