// Wnt / beta-catenin signalling (paper: "Stage 1 - Wnt Signaling Pathway" and
// "Mathematical Models - Example: Wnt/beta-Catenin Pathway").
//
// 1) paperWnt: the paper's linear cascade, implemented exactly as written
//      dW/dt = k1 - k2 W
//      dF/dt = k3 W - k4 F
//      dD/dt = k5 F - k6 D
//      dB/dt = k7 D - k8 B
//    with the closed-form steady state W* = k1/k2, F* = k3 W*/k4, ...
//
// 2) canonicalWnt: a mechanistic extension of the same cascade that keeps the
//    paper's reaction steps (Wnt + Frizzled/LRP5/6 -> Dishevelled activation ->
//    inhibition of GSK-3beta -> beta-catenin stabilisation -> nuclear
//    translocation -> target genes) and adds the destruction complex and the
//    Axin2 negative feedback (Lee et al. 2003, PLoS Biol 1:e10, reduced form).

export const paperWnt = {
  species: ['W', 'F', 'D', 'B'],
  labels: ['Wnt ligand W', 'Frizzled (active) F', 'Dishevelled D', 'beta-catenin B'],
  defaults: { k1: 1.0, k2: 0.5, k3: 0.8, k4: 0.4, k5: 0.6, k6: 0.3, k7: 0.5, k8: 0.25 },
  rhs(p) {
    return (t, y, dy) => {
      dy[0] = p.k1 - p.k2 * y[0];
      dy[1] = p.k3 * y[0] - p.k4 * y[1];
      dy[2] = p.k5 * y[1] - p.k6 * y[2];
      dy[3] = p.k7 * y[2] - p.k8 * y[3];
    };
  },
  steadyState(p) {
    const W = p.k1 / p.k2, F = (p.k3 * W) / p.k4, D = (p.k5 * F) / p.k6, B = (p.k7 * D) / p.k8;
    return [W, F, D, B];
  },
  /** Exact solution of the first equation, used by the tests. */
  W(t, p, W0 = 0) { const Ws = p.k1 / p.k2; return Ws + (W0 - Ws) * Math.exp(-p.k2 * t); },
};

export const canonicalWnt = {
  // index: 0 Wnt, 1 Fz* (ligand-bound receptor), 2 Dvl*, 3 beta-cat cytoplasm,
  //        4 beta-cat nucleus, 5 Axin2, 6 target mRNA (e.g. Axin2/Lgr5 reporter)
  species: ['Wnt', 'FzLRP*', 'Dvl*', 'bcat_c', 'bcat_n', 'Axin2', 'target'],
  labels: ['Wnt', 'Fz/LRP6 bound', 'Dvl active', 'beta-cat (cyto)', 'beta-cat (nuc)', 'Axin2', 'TCF target'],
  defaults: {
    wntIn: 0.0,     // external Wnt supply (niche signal)
    dW: 0.4,        // Wnt clearance
    kon: 1.2, koff: 0.3, FzTot: 1.0,
    kDa: 1.5, kDi: 0.5,             // Dvl activation by Fz*, inactivation
    vB: 0.6,                         // beta-catenin synthesis
    kDC: 2.5,                        // destruction-complex degradation
    KDvl: 0.15,                      // Dvl inhibition of GSK-3beta (complex)
    kB: 0.05,                        // basal degradation
    kin: 0.8, kout: 0.4,             // nuclear import/export
    kAx: 0.6, KAx: 0.5, dAx: 0.3,    // Axin2 induction (negative feedback)
    kT: 1.0, KT: 0.6, dT: 0.4,       // target transcription (TCF/LEF)
    nH: 2,
  },
  rhs(p) {
    return (t, y, dy) => {
      const [W, Fz, Dv, Bc, Bn, Ax, T] = y;
      const free = Math.max(p.FzTot - Fz, 0);
      dy[0] = p.wntIn - p.dW * W - p.kon * W * free + p.koff * Fz;
      dy[1] = p.kon * W * free - p.koff * Fz;
      dy[2] = p.kDa * Fz * (1 - Dv) - p.kDi * Dv;
      // destruction complex activity: inhibited by active Dvl, scaled up by Axin2
      const dc = (1 + Ax) / (1 + Dv / p.KDvl);
      dy[3] = p.vB - (p.kB + p.kDC * dc) * Bc - p.kin * Bc + p.kout * Bn;
      dy[4] = p.kin * Bc - p.kout * Bn;
      const act = Bn ** p.nH / (p.KT ** p.nH + Bn ** p.nH);
      dy[5] = p.kAx * (Bn ** p.nH / (p.KAx ** p.nH + Bn ** p.nH)) - p.dAx * Ax;
      dy[6] = p.kT * act - p.dT * T;
    };
  },
  initial() { return [0, 0, 0, 0.2, 0.1, 0.1, 0.02]; },
};

/** Reaction list for the network view (order 1 = paper's first-order stage). */
export const REACTIONS = [
  { id: 'wnt.bind', order: 1, pathway: 'Wnt', reactants: ['Wnt', 'Frizzled', 'LRP5/6'], products: ['Wnt:Fz:LRP'], label: 'Ligand-receptor interaction' },
  { id: 'wnt.dvl', order: 1, pathway: 'Wnt', reactants: ['Dishevelled'], products: ['Dishevelled*'], modifiers: ['Wnt:Fz:LRP'], label: 'Dishevelled activation' },
  { id: 'wnt.gsk', order: 1, pathway: 'Wnt', reactants: ['GSK-3b'], products: ['GSK-3b (inhibited)'], modifiers: ['Dishevelled*'], label: 'Inhibition of GSK-3beta' },
  { id: 'wnt.dc', order: 1, pathway: 'Wnt', reactants: ['beta-catenin'], products: ['p-beta-catenin'], modifiers: ['GSK-3b', 'APC', 'Axin'], label: 'Destruction-complex phosphorylation' },
  { id: 'wnt.deg', order: 1, pathway: 'Wnt', reactants: ['p-beta-catenin'], products: [], modifiers: ['beta-TrCP'], label: 'Proteasomal degradation' },
  { id: 'wnt.syn', order: 1, pathway: 'Wnt', reactants: [], products: ['beta-catenin'], label: 'beta-catenin synthesis' },
  { id: 'wnt.imp', order: 1, pathway: 'Wnt', reactants: ['beta-catenin'], products: ['beta-catenin (nuc)'], label: 'Nuclear translocation' },
  { id: 'wnt.tcf', order: 1, pathway: 'Wnt', reactants: ['beta-catenin (nuc)', 'TCF/LEF'], products: ['beta-cat:TCF'], label: 'TCF/LEF complex' },
  { id: 'wnt.tx', order: 1, pathway: 'Wnt', reactants: [], products: ['Wnt target mRNA'], modifiers: ['beta-cat:TCF'], label: 'Target gene transcription' },
  { id: 'wnt.axin2', order: 1, pathway: 'Wnt', reactants: [], products: ['Axin'], modifiers: ['beta-cat:TCF'], label: 'Axin2 negative feedback' },
];
