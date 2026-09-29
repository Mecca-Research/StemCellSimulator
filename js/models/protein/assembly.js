// Construction of complex proteins (paper: "Constructing Complex Globular
// Proteins" - hemoglobin and collagen), as mass-action kinetics that follow the
// paper's reaction schemes step by step.
//
// Hemoglobin (paper):  DNA(HBA1/2, HBB) -> mRNA -> alpha/beta polypeptides
//   -> 2 alpha + 2 beta -> tetramer -> + 4 heme -> functional hemoglobin.
// Added biology that makes the stoichiometry matter:
//   * AHSP (alpha-hemoglobin stabilising protein) chaperones free alpha chains
//     (Kihm et al. 2002, Nature 417:758);
//   * unpaired alpha chains precipitate (the lesion of beta-thalassaemia);
//   * HRI (heme-regulated eIF2alpha kinase) throttles globin translation when
//     heme is scarce, coupling the two supplies (Chen 2007, Blood 109:2693).
//   Heme is added at the tetramer level exactly as the paper writes it; in
//   cells most heme is inserted into nascent/monomeric globins.
//
// Collagen type I (paper): DNA(COL1A1, COL1A2) -> pro-alpha1, pro-alpha2 ->
//   hydroxylation + glycosylation -> 2 pro-alpha1 + 1 pro-alpha2 -> procollagen
//   triple helix -> secretion -> propeptide cleavage -> tropocollagen -> fibril.
// Added biology: prolyl-4-hydroxylase needs ascorbate (vitamin C); triple
// helices of under-hydroxylated chains are unstable at 37 C and are degraded
// (scurvy); ADAMTS2 (N-proteinase) and BMP1 (C-proteinase) remove the
// propeptides; fibrils form by nucleation and elongation (Kadler et al. 1996,
// Biochem J 316:1).

export const hemoglobin = {
  species: ['mRNA_a', 'mRNA_b', 'alpha', 'alpha:AHSP', 'beta', 'ab dimer', 'a2b2 (apo)', 'Hb+1 heme', 'Hb+2 heme', 'Hb+3 heme', 'Hb (4 heme)', 'heme', 'alpha precipitate', 'alpha made'],
  idx: { ma: 0, mb: 1, a: 2, aA: 3, b: 4, D: 5, T0: 6, T1: 7, T2: 8, T3: 9, T4: 10, heme: 11, P: 12, aSyn: 13 },
  defaults: {
    txA: 1.0, txB: 1.0, betaExpr: 1.0, alphaExpr: 1.0, dm: 0.1, // transcription, mRNA decay
    tl: 2.0,                        // translation per mRNA
    hri: 1.0, KhemeHRI: 0.5,        // HRI coupling strength, heme half-activation
    kd: 1.0,                        // alpha + beta -> dimer
    kt: 0.5,                        // dimer + dimer -> tetramer
    AHSP: 20.0, kA: 2.0, kAr: 0.2,  // chaperone pool / binding / release
    kprec: 0.05,                    // free alpha precipitation
    kbDeg: 0.02,                    // excess beta turnover (beta4 / HbH is ignored)
    hemeSyn: 50.0, kh: 0.4, dHeme: 0.05,
  },
  rhs(p) {
    const I = hemoglobin.idx;
    return (t, y, dy) => {
      const a = y[I.a], aA = y[I.aA], b = y[I.b], D = y[I.D], heme = y[I.heme];
      const freeAHSP = Math.max(p.AHSP - aA, 0);
      const pair = p.kd * (a + aA) * b;       // AHSP-bound alpha is handed over to beta
      const pairFree = p.kd * a * b, pairAHSP = p.kd * aA * b;
      const tet = p.kt * D * D;
      const hb = [0, 1, 2, 3].map((n) => p.kh * y[I.T0 + n] * heme);
      const tl = p.tl * (1 - p.hri + p.hri * heme / (p.KhemeHRI + heme));
      dy[I.ma] = p.txA * p.alphaExpr - p.dm * y[I.ma];
      dy[I.mb] = p.txB * p.betaExpr - p.dm * y[I.mb];
      dy[I.a] = tl * y[I.ma] - pairFree - p.kA * a * freeAHSP + p.kAr * aA - p.kprec * a;
      dy[I.aA] = p.kA * a * freeAHSP - p.kAr * aA - pairAHSP;
      dy[I.b] = tl * y[I.mb] - pair - p.kbDeg * b;
      dy[I.D] = pair - 2 * tet;
      dy[I.T0] = tet - hb[0];
      dy[I.T1] = hb[0] - hb[1];
      dy[I.T2] = hb[1] - hb[2];
      dy[I.T3] = hb[2] - hb[3];
      dy[I.T4] = hb[3];
      dy[I.heme] = p.hemeSyn - (hb[0] + hb[1] + hb[2] + hb[3]) - p.dHeme * heme;
      dy[I.P] = p.kprec * a;
      dy[I.aSyn] = tl * y[I.ma];
    };
  },
  initial() { return new Array(hemoglobin.species.length).fill(0); },
  /** alpha chains accounted for in every compartment (must equal alpha made) */
  alphaBalance(y) {
    const I = hemoglobin.idx;
    return y[I.a] + y[I.aA] + y[I.D] + 2 * (y[I.T0] + y[I.T1] + y[I.T2] + y[I.T3] + y[I.T4]) + y[I.P];
  },
};

export const collagen = {
  species: ['mRNA COL1A1', 'mRNA COL1A2', 'pro-a1', 'pro-a2', 'pro-a1 (OH, glyc)', 'pro-a2 (OH, glyc)', 'procollagen (ER)', 'procollagen (secreted)', 'tropocollagen', 'fibril nuclei', 'fibril mass', 'degraded', 'a1 made', 'a2 made'],
  idx: { m1: 0, m2: 1, c1: 2, c2: 3, h1: 4, h2: 5, PC: 6, PCe: 7, TC: 8, N: 9, F: 10, deg: 11, s1: 12, s2: 13 },
  defaults: {
    tx1: 2.0, tx2: 1.0, col1a2: 1.0, dm: 0.1, tl: 1.0,
    kOH: 1.5, ascorbate: 1.0, Kasc: 0.2,   // prolyl-4-hydroxylase (vitamin C dependent)
    kTri: 0.8,                             // 2 a1 + 1 a2 -> triple helix (C-propeptide registration)
    kUnstable: 0.6,                        // under-hydroxylated chains misfold and are degraded
    kSec: 0.5,                             // secretion
    kNP: 0.6, ADAMTS2: 1.0,                // N- and C-propeptide removal
    kNuc: 0.02, kEl: 0.6,                  // fibril nucleation / elongation
    dChain: 0.05,
  },
  rhs(p) {
    const I = collagen.idx;
    return (t, y, dy) => {
      const asc = p.ascorbate / (p.Kasc + p.ascorbate);
      const oh1 = p.kOH * asc * y[I.c1], oh2 = p.kOH * asc * y[I.c2];
      const tri = p.kTri * y[I.h1] * y[I.h1] * y[I.h2];
      const bad1 = p.kUnstable * y[I.c1] * (1 - asc), bad2 = p.kUnstable * y[I.c2] * (1 - asc);
      const cleave = p.kNP * p.ADAMTS2 * y[I.PCe];
      const nuc = p.kNuc * y[I.TC] * y[I.TC];
      const el = p.kEl * y[I.TC] * y[I.N] / (1 + y[I.N]);
      dy[I.m1] = p.tx1 - p.dm * y[I.m1];
      dy[I.m2] = p.tx2 * p.col1a2 - p.dm * y[I.m2];
      dy[I.c1] = p.tl * y[I.m1] - oh1 - bad1 - p.dChain * y[I.c1];
      dy[I.c2] = p.tl * y[I.m2] - oh2 - bad2 - p.dChain * y[I.c2];
      dy[I.h1] = oh1 - 2 * tri - p.dChain * y[I.h1];
      dy[I.h2] = oh2 - tri - p.dChain * y[I.h2];
      dy[I.PC] = tri - p.kSec * y[I.PC];
      dy[I.PCe] = p.kSec * y[I.PC] - cleave;
      dy[I.TC] = cleave - 2 * nuc - el;
      dy[I.N] = nuc;
      dy[I.F] = 2 * nuc + el;
      dy[I.deg] = bad1 + bad2 + p.dChain * (y[I.c1] + y[I.c2] + y[I.h1] + y[I.h2]);
      dy[I.s1] = p.tl * y[I.m1];
      dy[I.s2] = p.tl * y[I.m2];
    };
  },
  initial() { return new Array(collagen.species.length).fill(0); },
  /** chains (alpha1 + alpha2) accounted for in every compartment */
  chainBalance(y) {
    const I = collagen.idx;
    const trimers = y[I.PC] + y[I.PCe] + y[I.TC] + y[I.F];
    return y[I.c1] + y[I.c2] + y[I.h1] + y[I.h2] + 3 * trimers + y[I.deg];
  },
};

const hb = (id, r, pr, lbl, mod) => ({ id: `hb.${id}`, order: 3, pathway: 'Hemoglobin synthesis', reactants: r, products: pr, modifiers: mod, label: lbl });
const col = (id, r, pr, lbl, mod) => ({ id: `col.${id}`, order: 3, pathway: 'Collagen synthesis', reactants: r, products: pr, modifiers: mod, label: lbl });

export const REACTIONS = [
  hb('txA', [], ['mRNA (alpha-globin)'], 'Transcription of HBA1/HBA2', ['RNA polymerase II', 'GATA1']),
  hb('txB', [], ['mRNA (beta-globin)'], 'Transcription of HBB', ['RNA polymerase II', 'GATA1']),
  hb('tlA', ['mRNA (alpha-globin)', 'tRNA'], ['alpha-globin'], 'Translation (alpha)', ['Ribosome']),
  hb('tlB', ['mRNA (beta-globin)', 'tRNA'], ['beta-globin'], 'Translation (beta)', ['Ribosome']),
  hb('ahsp', ['alpha-globin', 'AHSP'], ['alpha-globin:AHSP'], 'Chaperone binding of free alpha'),
  hb('dimer', ['alpha-globin', 'beta-globin'], ['alpha-beta dimer'], 'Folding and dimer assembly', ['Chaperones']),
  hb('tetramer', ['alpha-beta dimer', 'alpha-beta dimer'], ['Hemoglobin tetramer'], 'Tetramer assembly'),
  hb('heme', ['Hemoglobin tetramer', 'Heme'], ['Hemoglobin'], 'Heme incorporation (x4)'),
  hb('hemeSyn', ['Succinyl-CoA', 'Glycine'], ['Heme'], 'Heme synthesis', ['ALAS2', 'Ferrochelatase']),
  hb('precip', ['alpha-globin'], ['alpha-globin precipitate'], 'Unpaired alpha precipitation'),
  col('tx1', [], ['mRNA (pro-alpha1)'], 'Transcription of COL1A1', ['RNA polymerase II']),
  col('tx2', [], ['mRNA (pro-alpha2)'], 'Transcription of COL1A2', ['RNA polymerase II']),
  col('tl1', ['mRNA (pro-alpha1)', 'tRNA'], ['pro-alpha1 chain'], 'Translation into the ER', ['Ribosome']),
  col('tl2', ['mRNA (pro-alpha2)', 'tRNA'], ['pro-alpha2 chain'], 'Translation into the ER', ['Ribosome']),
  col('oh', ['pro-alpha1 chain', 'O2', '2-oxoglutarate'], ['Modified pro-alpha1'], 'Prolyl/lysyl hydroxylation', ['Prolyl-4-hydroxylase', 'Ascorbate']),
  col('glyc', ['Modified pro-alpha1', 'UDP-galactose'], ['Glycosylated pro-alpha1'], 'Hydroxylysine glycosylation'),
  col('trimer', ['Glycosylated pro-alpha1', 'Glycosylated pro-alpha1', 'Modified pro-alpha2'], ['Procollagen triple helix'], 'Triple-helix formation', ['HSP47', 'PDI']),
  col('sec', ['Procollagen triple helix'], ['Procollagen (extracellular)'], 'Secretion'),
  col('cleave', ['Procollagen (extracellular)'], ['Tropocollagen', 'Propeptides'], 'Propeptide cleavage', ['ADAMTS2', 'BMP1']),
  col('fibril', ['Tropocollagen'], ['Collagen fibril'], 'Fibril assembly', ['Lysyl oxidase']),
];
