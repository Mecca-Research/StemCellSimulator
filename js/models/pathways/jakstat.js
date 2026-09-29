// JAK2/STAT5 signalling downstream of the erythropoietin receptor (paper:
// "Stage 2 - JAK/STAT pathway": cytokine binding, receptor dimerisation
// activates JAK kinases, STAT phosphorylation and nuclear translocation;
// worked example "CMP --(EPO)--> proerythroblast --(EPO, GATA1)--> ...").
//
// Structure after Swameye, Mueller, Timmer, Sandra & Klingmueller 2003, PNAS
// 100:1028 "Identification of nucleocytoplasmic cycling as a remote sensor in
// cellular signaling by databased modeling":
//   x1 STAT5 --(pEpoR/JAK2)--> x2 pSTAT5 ;  2 x2 --> x3 pSTAT5 dimer (cyto)
//   x3 --import--> nuclear dimer --dephosphorylation--> 2 nuclear STAT5
//   --export--> x1
// In Swameye et al. the nuclear sojourn is a pure delay (tau ~ 6 min); here it
// is resolved into the two physical steps (nuclear dephosphorylation, export)
// whose mean residence time 1/kdp + 1/kexp equals that delay (linear-chain
// representation). Total STAT5 is conserved:
//   STAT5 + pSTAT5 + 2 [pSTAT5 dimer] + 2 [pSTAT5 dimer (nuc)] + [STAT5 (nuc)] = STATtot.
// Negative feedback after Vera, Bachmann, Pfeifer, Becker, Hormiga, Torres
// Darias, Timmer, Klingmueller & Wolkenhauer 2008, BMC Syst Biol 2:38: nuclear
// pSTAT5 induces CIS (competes with STAT5 for the phospho-tyrosine docking
// sites of EpoR) and SOCS (inhibits JAK2, accelerating receptor inactivation).
// Receptor activation is lumped (EPO binding, preformed-dimer rearrangement,
// JAK2 trans-phosphorylation). Units: minutes; STAT5 as a fraction of total;
// EPO in U/ml. Parameters are illustrative, chosen so that nuclear pSTAT5
// peaks ~10-20 min after stimulation as in Swameye et al.
import { rk4Step, rk4Work } from '../core/ode.js';

export const jakstat = {
  // 0 pEpoR/JAK2 (active fraction), 1 STAT5, 2 pSTAT5, 3 pSTAT5 dimer, 4 pSTAT5 dimer (nuc),
  // 5 STAT5 (nuc), 6 CIS mRNA, 7 CIS, 8 SOCS mRNA, 9 SOCS
  species: ['pEpoR', 'STAT5', 'pSTAT5', 'pSTAT5:pSTAT5', 'pSTAT5:pSTAT5 (nuc)', 'STAT5 (nuc)', 'CIS mRNA', 'CIS', 'SOCS mRNA', 'SOCS'],
  defaults: {
    EPO: 1,          // U/ml
    kact: 0.3,       // (U/ml)^-1 min^-1  receptor/JAK2 activation
    kdeact: 0.1,     // min^-1            basal inactivation (SHP-1, internalisation)
    kp: 0.25,        // min^-1            STAT5 phosphorylation per active receptor
    kdim: 4.0,       // min^-1 per unit   pSTAT5 dimerisation
    kimp: 0.1066,    // min^-1            dimer nuclear import (Swameye et al. k3)
    kdp: 0.31,       // min^-1            nuclear dephosphorylation
    kexp: 0.31,      // min^-1            export; 1/kdp + 1/kexp = 6.4 min (Swameye et al. tau)
    STATtot: 1,
    fb: 1,           // feedback strength (0 = CIS/SOCS knocked out)
    kmC: 1.2, KC: 0.05, dmC: 0.2, ktlC: 0.12, dC: 0.04, KiC: 0.5,     // CIS
    kmS: 1.2, KS: 0.05, dmS: 0.26, ktlS: 0.12, dS: 0.05, KiS: 0.75,   // SOCS
  },
  initial(p = jakstat.defaults) { return [0, p.STATtot, 0, 0, 0, 0, 0, 0, 0, 0]; },
  rhs(p) {
    return (t, y, dy) => {
      const [R, x1, x2, x3, n2, n1, mC, C, mS, S] = y;
      const deact = p.kdeact * (1 + S / p.KiS);
      dy[0] = p.kact * p.EPO * (1 - R) - deact * R;
      const vp = (p.kp * R * x1) / (1 + C / p.KiC);
      const vd = p.kdim * x2 * x2;
      const vi = p.kimp * x3;
      const vdp = p.kdp * n2;
      const ve = p.kexp * n1;
      dy[1] = -vp + ve;
      dy[2] = vp - 2 * vd;
      dy[3] = vd - vi;
      dy[4] = vi - vdp;
      dy[5] = 2 * vdp - ve;
      // pSTAT5-dimer target genes (the pSTAT5 tetramer/dimer on GAS elements)
      const act = n2 / (p.KC + n2), actS = n2 / (p.KS + n2);
      dy[6] = p.fb * p.kmC * act - p.dmC * mC;
      dy[7] = p.ktlC * mC - p.dC * C;
      dy[8] = p.fb * p.kmS * actS - p.dmS * mS;
      dy[9] = p.ktlS * mS - p.dS * S;
    };
  },
  /** Total STAT5 (in monomer units) - conserved by the dynamics. */
  totalSTAT(y) { return y[1] + y[2] + 2 * y[3] + 2 * y[4] + y[5]; },
  /** Nuclear pSTAT5 as a fraction of total STAT5 (monomer units). */
  nuclearPSTAT(y, p = jakstat.defaults) { return (2 * y[4]) / p.STATtot; },
};

/**
 * Time course after an EPO step at t = 0 (fixed-step RK4).
 * @returns {{t:number[], y:number[][]}} samples every `every` minutes
 */
export function timeCourse(p, tEnd = 120, { h = 0.05, every = 1 } = {}) {
  const f = jakstat.rhs(p);
  const y = Float64Array.from(jakstat.initial(p));
  const w = rk4Work(y.length);
  const ts = [0], ys = [Array.from(y)];
  const steps = Math.round(tEnd / h), stride = Math.max(1, Math.round(every / h));
  for (let s = 1; s <= steps; s++) {
    rk4Step(f, (s - 1) * h, y, h, w);
    if (s % stride === 0) { ts.push(s * h); ys.push(Array.from(y)); }
  }
  return { t: ts, y: ys };
}

/** Peak and late (tEnd) nuclear pSTAT5 fraction after an EPO step. */
export function responseSummary(p, tEnd = 240) {
  const { t, y } = timeCourse(p, tEnd, { every: 0.5 });
  let peak = 0, tPeak = 0;
  for (let i = 0; i < t.length; i++) {
    const v = jakstat.nuclearPSTAT(y[i], p);
    if (v > peak) { peak = v; tPeak = t[i]; }
  }
  return { peak, tPeak, late: jakstat.nuclearPSTAT(y[y.length - 1], p) };
}

/**
 * EPO dose -> sustained nuclear pSTAT5 (value at tEnd, normalised to the value
 * at a saturating dose). Used by the scenes as the JAK/STAT drive of GATA1
 * and of erythroid progenitor survival.
 */
export function stat5Drive(EPO, p = jakstat.defaults, tEnd = 240) {
  if (EPO <= 0) return 0;
  const v = responseSummary({ ...p, EPO }, tEnd).late;
  const vmax = responseSummary({ ...p, EPO: 100 }, tEnd).late;
  return vmax > 0 ? v / vmax : 0;
}

/** Reaction list for the network view (order 2 = the paper's progenitor stage). */
export const REACTIONS = [
  { id: 'jak.bind', order: 2, pathway: 'JAK/STAT', reactants: ['EPO', 'EpoR dimer'], products: ['EPO:EpoR'], label: 'Cytokine binds its receptor dimer' },
  { id: 'jak.jak2', order: 2, pathway: 'JAK/STAT', reactants: ['JAK2'], products: ['JAK2-P'], modifiers: ['EPO:EpoR'], label: 'JAK2 trans-phosphorylation' },
  { id: 'jak.epor', order: 2, pathway: 'JAK/STAT', reactants: ['EPO:EpoR'], products: ['pEpoR'], modifiers: ['JAK2-P'], label: 'Receptor tyrosine phosphorylation (docking sites)' },
  { id: 'jak.off', order: 2, pathway: 'JAK/STAT', reactants: ['pEpoR'], products: ['EpoR dimer'], modifiers: ['SHP-1', 'SOCS'], label: 'Receptor inactivation' },
  { id: 'jak.stat5', order: 2, pathway: 'JAK/STAT', reactants: ['STAT5'], products: ['pSTAT5'], modifiers: ['pEpoR', 'JAK2-P'], label: 'STAT5 tyrosine phosphorylation' },
  { id: 'jak.dimer', order: 2, pathway: 'JAK/STAT', reactants: ['pSTAT5', 'pSTAT5'], products: ['pSTAT5:pSTAT5'], label: 'pSTAT5 dimerisation' },
  { id: 'jak.imp', order: 2, pathway: 'JAK/STAT', reactants: ['pSTAT5:pSTAT5'], products: ['pSTAT5:pSTAT5 (nuc)'], label: 'Nuclear translocation' },
  { id: 'jak.dephos', order: 2, pathway: 'JAK/STAT', reactants: ['pSTAT5:pSTAT5 (nuc)'], products: ['STAT5 (nuc)', 'STAT5 (nuc)'], modifiers: ['nuclear phosphatase'], label: 'Nuclear dephosphorylation' },
  { id: 'jak.exp', order: 2, pathway: 'JAK/STAT', reactants: ['STAT5 (nuc)'], products: ['STAT5'], label: 'Nuclear export (cycling)' },
  { id: 'jak.cis', order: 2, pathway: 'JAK/STAT', reactants: [], products: ['CIS mRNA'], modifiers: ['pSTAT5:pSTAT5 (nuc)'], label: 'CIS transcription' },
  { id: 'jak.cis.tl', order: 2, pathway: 'JAK/STAT', reactants: [], products: ['CIS'], modifiers: ['CIS mRNA'], label: 'CIS translation' },
  { id: 'jak.cis.fb', order: 2, pathway: 'JAK/STAT', reactants: ['CIS', 'pEpoR'], products: ['CIS:pEpoR'], label: 'CIS masks STAT5 docking sites (feedback)' },
  { id: 'jak.socs', order: 2, pathway: 'JAK/STAT', reactants: [], products: ['SOCS mRNA'], modifiers: ['pSTAT5:pSTAT5 (nuc)'], label: 'SOCS transcription' },
  { id: 'jak.socs.tl', order: 2, pathway: 'JAK/STAT', reactants: [], products: ['SOCS'], modifiers: ['SOCS mRNA'], label: 'SOCS translation' },
  { id: 'jak.socs.fb', order: 2, pathway: 'JAK/STAT', reactants: ['SOCS', 'JAK2-P'], products: ['SOCS:JAK2'], label: 'SOCS inhibits JAK2 (feedback)' },
  { id: 'jak.gata1', order: 2, pathway: 'JAK/STAT', reactants: [], products: ['GATA1'], modifiers: ['pSTAT5:pSTAT5 (nuc)'], label: 'EPO/STAT5 support of the erythroid programme' },
];
