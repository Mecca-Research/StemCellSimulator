// Notch signalling (paper: "Stage 1 - Notch Signaling Pathway": ligand binding
// -> receptor cleavage -> NICD -> nucleus -> transcription).
//
// 1) collier: juxtacrine lateral inhibition on any contact graph
//    (Collier, Monk, Maini & Lewis 1996, J Theor Biol 183:429)
//      dN_i/dt = f(<D>_i) - N_i          f(x) = x^k / (a + x^k)
//      dD_i/dt = v (g(N_i) - D_i)        g(x) = 1 / (1 + b x^h)
//    Neighbouring cells amplify small differences into a salt-and-pepper
//    pattern of Delta-high (sender, e.g. neural/secretory fate) and
//    Notch-high (receiver, progenitor-maintaining) cells.
//
// 2) singleCell: the mechanistic steps inside one receiver cell
//    Delta (trans) + Notch <-> DN  -> S2 cleavage (ADAM10) -> NEXT
//    -> S3 cleavage (gamma-secretase) -> NICD (cyto) -> NICD (nuc)
//    -> Hes1 transcription; Hes1 represses Delta (cis) = lateral inhibition.

export const collier = {
  defaults: { a: 0.01, b: 100, k: 2, h: 2, v: 1 },

  /** Allocate state for n cells (small random perturbation breaks symmetry). */
  init(n, rng, amp = 0.05) {
    const N = new Float64Array(n), D = new Float64Array(n);
    for (let i = 0; i < n; i++) { N[i] = 0.5 + amp * (rng.next() - 0.5); D[i] = 0.5 + amp * (rng.next() - 0.5); }
    return { N, D };
  },

  /**
   * Explicit Euler step (stable for dt <= 0.1 with v = 1).
   * @param {{N:Float64Array,D:Float64Array}} s
   * @param {number[][]} nbrs neighbour index lists
   */
  step(s, nbrs, dt, p = collier.defaults) {
    const n = s.N.length;
    // scratch buffers are cached on the state to avoid per-step allocation
    if (!s._dN || s._dN.length !== n) { s._dN = new Float64Array(n); s._dD = new Float64Array(n); }
    const dN = s._dN, dD = s._dD;
    for (let i = 0; i < n; i++) {
      const nb = nbrs[i];
      let avg = 0;
      if (nb && nb.length) { for (const j of nb) avg += s.D[j]; avg /= nb.length; }
      const xk = avg ** p.k;
      dN[i] = xk / (p.a + xk) - s.N[i];
      dD[i] = p.v * (1 / (1 + p.b * s.N[i] ** p.h) - s.D[i]);
    }
    for (let i = 0; i < n; i++) { s.N[i] += dt * dN[i]; s.D[i] += dt * dD[i]; }
  },

  /** Fraction of cells that are Delta-high (sender fate). */
  senderFraction(s, thr = 0.5) {
    let c = 0;
    for (let i = 0; i < s.D.length; i++) if (s.D[i] > thr) c++;
    return c / s.D.length;
  },

  /** Lateral-inhibition quality: fraction of sender-sender contacts (0 = perfect). */
  senderAdjacency(s, nbrs, thr = 0.5) {
    let ss = 0, total = 0;
    for (let i = 0; i < nbrs.length; i++) {
      if (s.D[i] <= thr) continue;
      for (const j of nbrs[i]) { total++; if (s.D[j] > thr) ss++; }
    }
    return total ? ss / total : 0;
  },
};

export const singleCell = {
  species: ['Notch', 'Delta:Notch', 'NEXT', 'NICD_c', 'NICD_n', 'Hes1', 'Delta_cis'],
  labels: ['Notch (surface)', 'Delta:Notch', 'NEXT (S2-cleaved)', 'NICD cytoplasm', 'NICD nucleus', 'Hes1', 'Delta (own)'],
  defaults: {
    deltaTrans: 1.0,          // Delta presented by the neighbouring sender cell
    sN: 1.0, dN: 0.3,          // receptor delivery and turnover
    kon: 1.5, koff: 0.2,       // trans binding
    kS2: 0.9,                  // ADAM10 cleavage (pulling force exposes S2 site)
    kS3: 1.2,                  // gamma-secretase (0 = DAPT inhibitor)
    kin: 0.6, kout: 0.1,       // nuclear import/export (scaled to measured import kinetics)
    dNICD: 0.35,               // NICD turnover (Fbw7)
    kH: 1.2, KH: 0.3, dH: 0.5, // Hes1 induction by NICD:CSL:MAML
    sD: 1.0, KD: 0.4, dD: 0.4, // cis Delta, repressed by Hes1
  },
  rhs(p) {
    return (t, y, dy) => {
      const [N, DN, NEXT, Nc, Nn, H, Dc] = y;
      dy[0] = p.sN - p.dN * N - p.kon * p.deltaTrans * N + p.koff * DN;
      dy[1] = p.kon * p.deltaTrans * N - p.koff * DN - p.kS2 * DN;
      dy[2] = p.kS2 * DN - p.kS3 * NEXT - 0.05 * NEXT;
      dy[3] = p.kS3 * NEXT - p.kin * Nc + p.kout * Nn - p.dNICD * Nc;
      dy[4] = p.kin * Nc - p.kout * Nn - p.dNICD * Nn;
      dy[5] = p.kH * (Nn * Nn) / (p.KH * p.KH + Nn * Nn) - p.dH * H;
      dy[6] = p.sD / (1 + (H / p.KD) ** 2) - p.dD * Dc;
    };
  },
  initial() { return [3, 0, 0, 0, 0, 0, 2.5]; },
};

export const REACTIONS = [
  { id: 'notch.bind', order: 1, pathway: 'Notch', reactants: ['Delta/Jagged', 'Notch'], products: ['Delta:Notch'], label: 'Ligand binding (trans)' },
  { id: 'notch.s2', order: 1, pathway: 'Notch', reactants: ['Delta:Notch'], products: ['NEXT', 'Notch ECD'], modifiers: ['ADAM10'], label: 'S2 cleavage' },
  { id: 'notch.endo', order: 1, pathway: 'Notch', reactants: ['Notch ECD'], products: [], modifiers: ['Delta/Jagged'], label: 'Trans-endocytosis of ECD' },
  { id: 'notch.s3', order: 1, pathway: 'Notch', reactants: ['NEXT'], products: ['NICD'], modifiers: ['gamma-secretase'], label: 'S3 cleavage releases NICD' },
  { id: 'notch.imp', order: 1, pathway: 'Notch', reactants: ['NICD'], products: ['NICD (nuc)'], label: 'Nuclear translocation' },
  { id: 'notch.csl', order: 1, pathway: 'Notch', reactants: ['NICD (nuc)', 'CSL', 'MAML'], products: ['NICD:CSL:MAML'], label: 'Transcription-factor activation' },
  { id: 'notch.hes', order: 1, pathway: 'Notch', reactants: [], products: ['Hes1'], modifiers: ['NICD:CSL:MAML'], label: 'Hes/Hey transcription' },
  { id: 'notch.li', order: 1, pathway: 'Notch', reactants: [], products: ['Delta/Jagged'], modifiers: ['Hes1'], label: 'Hes1 represses Delta (lateral inhibition)' },
  { id: 'notch.deg', order: 1, pathway: 'Notch', reactants: ['NICD (nuc)'], products: [], modifiers: ['Fbw7'], label: 'NICD turnover' },
];
