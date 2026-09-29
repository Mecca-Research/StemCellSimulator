// Feedback loops in hormonal signalling (paper: "Third order - hormone
// secretion in response to stimuli, binding to receptors on target cells,
// negative feedback inhibition maintaining homeostasis").
//
// Goodwin negative-feedback oscillator (Goodwin 1965, Adv Enzyme Regul 3:425)
// applied to a hypothalamus -> pituitary -> gland axis (as for the
// hypothalamic-pituitary-gonadal axis by Smith 1980, Bull Math Biol 42:57;
// see also Murray, Mathematical Biology I, 3rd ed., Springer 2002):
//
//   dX/dt = s a / (1 + (Z/K)^n) - b X     releasing hormone (hypothalamus), s = stimulus
//   dY/dt = c X - d Y                     tropic hormone (pituitary)
//   dZ/dt = e Y - f Z                     end hormone (target gland), feeds back on X
//
// Linearising about the steady state gives the characteristic polynomial
//   (lambda + b)(lambda + d)(lambda + f) + G b d f = 0,
//   G = n u / (1 + u),  u = (Z*/K)^n   (loop gain = log-sensitivity of the repression)
// and Routh-Hurwitz puts the Hopf bifurcation at
//   G_c = (b + d + f)(bd + bf + df) / (b d f) - 1,
// which is 8 for equal degradation rates: sustained oscillation needs a Hill
// coefficient n > 8 (Griffith 1968, J Theor Biol 20:202).
//
// Receptor occupancy of target cells by Z (mass action, depletion of Z by
// binding neglected):  dO/dt = k_on Z (1 - O) - k_off O,  O* = Z / (Z + K_d).
// Time unit: minutes. Rates are illustrative (dimensionless axis), not fitted
// to a specific endocrine system.

export const goodwin = {
  species: ['X', 'Y', 'Z', 'O'],
  labels: ['releasing hormone X', 'tropic hormone Y', 'end hormone Z', 'receptor occupancy O'],
  defaults: { a: 1, b: 0.1, c: 1, d: 0.1, e: 1, f: 0.1, K: 1, n: 12, stim: 1, kon: 0.5, koff: 1 },

  /**
   * @param {object} p parameters; p.stim may be a number or a function of t
   */
  rhs(p) {
    const stim = typeof p.stim === 'function' ? p.stim : () => p.stim;
    return (t, y, dy) => {
      const [X, Y, Z, O] = y;
      const zk = Math.pow(Math.max(Z, 0) / p.K, p.n);
      dy[0] = (stim(t) * p.a) / (1 + zk) - p.b * X;
      dy[1] = p.c * X - p.d * Y;
      dy[2] = p.e * Y - p.f * Z;
      dy[3] = p.kon * Z * (1 - O) - p.koff * O;
    };
  },

  /** Feedback-receptor occupancy (Hill form of the repression term). */
  feedbackOccupancy(Z, p) { const u = Math.pow(Math.max(Z, 0) / p.K, p.n); return u / (1 + u); },

  /** Unique steady state (bisection on Z; the right-hand side is monotone). */
  steadyState(p) {
    const s = typeof p.stim === 'function' ? p.stim(Infinity) : p.stim;
    const gain = (p.c * p.e) / (p.d * p.f);
    const F = (Z) => Z - (gain * s * p.a) / (p.b * (1 + Math.pow(Z / p.K, p.n)));
    let lo = 0, hi = Math.max((gain * s * p.a) / p.b, 1e-9);
    for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (F(m) > 0) hi = m; else lo = m; }
    const Z = 0.5 * (lo + hi);
    const Y = (p.f * Z) / p.e, X = (p.d * Y) / p.c;
    return [X, Y, Z, (p.kon * Z) / (p.kon * Z + p.koff)];
  },

  /** Linear stability: loop gain G vs the Hopf threshold G_c. */
  stability(p) {
    const [, , Z] = goodwin.steadyState(p);
    const u = Math.pow(Z / p.K, p.n);
    const G = (p.n * u) / (1 + u);
    const { b, d, f } = p;
    const Gc = ((b + d + f) * (b * d + b * f + d * f)) / (b * d * f) - 1;
    // at the Hopf point the imaginary pair is +-i omega with omega^2 = bd + bf + df
    const omega = Math.sqrt(b * d + b * f + d * f);
    return { Z, gain: G, gainCrit: Gc, stable: G < Gc, omegaHopf: omega, periodHopf: (2 * Math.PI) / omega };
  },

  /** Smallest Hill coefficient that destabilises the steady state (bisection on n). */
  criticalN(p, nMax = 64) {
    const unstable = (n) => !goodwin.stability({ ...p, n }).stable;
    if (!unstable(nMax)) return Infinity;
    let lo = 1, hi = nMax;
    for (let i = 0; i < 60; i++) { const m = 0.5 * (lo + hi); if (unstable(m)) hi = m; else lo = m; }
    return 0.5 * (lo + hi);
  },
};

/**
 * Classify a time series as sustained oscillation or damped response from
 * its late-time peaks.
 * @param {number[]} ts
 * @param {number[]} xs
 * @param {{from?: number}} [o] only samples with t >= from are analysed
 * @returns {{peaks:number[][], amplitude:number, relAmplitude:number, decay:number, period:number, sustained:boolean}}
 */
export function oscillationStats(ts, xs, { from = 0 } = {}) {
  const peaks = [], troughs = [];
  let mean = 0, n = 0, max = -Infinity, min = Infinity;
  for (let i = 1; i < xs.length - 1; i++) {
    if (ts[i] < from) continue;
    const x = xs[i];
    mean += x; n++;
    if (x > max) max = x;
    if (x < min) min = x;
    if (x > xs[i - 1] && x >= xs[i + 1]) peaks.push([ts[i], x]);
    if (x < xs[i - 1] && x <= xs[i + 1]) troughs.push([ts[i], x]);
  }
  mean /= Math.max(n, 1);
  const amplitude = n ? (max - min) / 2 : 0;
  const relAmplitude = mean > 0 ? amplitude / mean : 0;
  // peak-height ratio across the window (1 = limit cycle, < 1 = damped spiral)
  let decay = 0, period = NaN;
  if (peaks.length >= 2 && troughs.length >= 2) {
    const first = peaks[0][1] - troughs[0][1];
    const last = peaks.at(-1)[1] - troughs.at(-1)[1];
    decay = first > 0 ? last / first : 0;
    period = (peaks.at(-1)[0] - peaks[0][0]) / (peaks.length - 1);
  }
  const sustained = relAmplitude > 0.05 && decay > 0.9;
  return { peaks, amplitude, relAmplitude, decay, period, sustained };
}

/** Reaction list for the network view (order 3 = differentiated cells and interactions). */
export const REACTIONS = [
  { id: 'hormone.secX', order: 3, pathway: 'Hormonal feedback', reactants: [], products: ['releasing hormone (X)'], modifiers: ['stimulus'], label: 'Hypothalamic secretion in response to a stimulus' },
  { id: 'hormone.fb', order: 3, pathway: 'Hormonal feedback', reactants: ['end hormone (Z)', 'feedback receptor'], products: ['Z:feedback receptor'], label: 'Feedback-receptor binding in the hypothalamus' },
  { id: 'hormone.inh', order: 3, pathway: 'Hormonal feedback', reactants: ['Z:feedback receptor', 'X gene promoter'], products: ['X gene (repressed)'], label: 'Negative feedback inhibition of X secretion (cooperativity n)' },
  { id: 'hormone.clrX', order: 3, pathway: 'Hormonal feedback', reactants: ['releasing hormone (X)'], products: [], label: 'Clearance of X' },
  { id: 'hormone.bindX', order: 3, pathway: 'Hormonal feedback', reactants: ['releasing hormone (X)', 'X receptor'], products: ['X:receptor'], label: 'X binds pituitary receptors' },
  { id: 'hormone.secY', order: 3, pathway: 'Hormonal feedback', reactants: [], products: ['tropic hormone (Y)'], modifiers: ['X:receptor'], label: 'Pituitary secretion of Y' },
  { id: 'hormone.clrY', order: 3, pathway: 'Hormonal feedback', reactants: ['tropic hormone (Y)'], products: [], label: 'Clearance of Y' },
  { id: 'hormone.bindY', order: 3, pathway: 'Hormonal feedback', reactants: ['tropic hormone (Y)', 'Y receptor'], products: ['Y:receptor'], label: 'Y binds target-gland receptors' },
  { id: 'hormone.secZ', order: 3, pathway: 'Hormonal feedback', reactants: [], products: ['end hormone (Z)'], modifiers: ['Y:receptor'], label: 'Target-gland secretion of Z' },
  { id: 'hormone.clrZ', order: 3, pathway: 'Hormonal feedback', reactants: ['end hormone (Z)'], products: [], label: 'Clearance of Z' },
  { id: 'hormone.bindZ', order: 3, pathway: 'Hormonal feedback', reactants: ['end hormone (Z)', 'Z receptor'], products: ['Z:receptor'], label: 'Receptor binding on target cells' },
  { id: 'hormone.resp', order: 3, pathway: 'Hormonal feedback', reactants: [], products: ['target-cell response'], modifiers: ['Z:receptor'], label: 'Target-cell response' },
];
