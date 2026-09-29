// ODE integrators and stochastic simulation used by every pathway model.
//
// A system is f(t, y, dy) that writes dy/dt into dy (Float64Array).

/** One classical Runge-Kutta 4 step, in place. */
export function rk4Step(f, t, y, h, work) {
  const n = y.length;
  const w = work ?? rk4Work(n);
  const { k1, k2, k3, k4, tmp } = w;
  f(t, y, k1);
  for (let i = 0; i < n; i++) tmp[i] = y[i] + 0.5 * h * k1[i];
  f(t + 0.5 * h, tmp, k2);
  for (let i = 0; i < n; i++) tmp[i] = y[i] + 0.5 * h * k2[i];
  f(t + 0.5 * h, tmp, k3);
  for (let i = 0; i < n; i++) tmp[i] = y[i] + h * k3[i];
  f(t + h, tmp, k4);
  for (let i = 0; i < n; i++) y[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  return y;
}

export function rk4Work(n) {
  return { k1: new Float64Array(n), k2: new Float64Array(n), k3: new Float64Array(n), k4: new Float64Array(n), tmp: new Float64Array(n) };
}

/** Integrate with fixed RK4 steps from t0 to t1; returns {t, y} samples every `every` steps. */
export function integrate(f, y0, t0, t1, h, { every = 1, clampNonNegative = false } = {}) {
  const y = Float64Array.from(y0);
  const work = rk4Work(y.length);
  const ts = [t0], ys = [Array.from(y)];
  const steps = Math.ceil((t1 - t0) / h - 1e-9);
  let t = t0;
  for (let s = 1; s <= steps; s++) {
    const hh = Math.min(h, t1 - t);
    rk4Step(f, t, y, hh, work);
    if (clampNonNegative) for (let i = 0; i < y.length; i++) if (y[i] < 0) y[i] = 0;
    t += hh;
    if (s % every === 0 || s === steps) { ts.push(t); ys.push(Array.from(y)); }
  }
  return { t: ts, y: ys };
}

// Dormand-Prince 5(4) coefficients
const A = [
  [],
  [1 / 5],
  [3 / 40, 9 / 40],
  [44 / 45, -56 / 15, 32 / 9],
  [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
  [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
  [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
];
const C = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
const B5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];
const B4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];

/**
 * Adaptive Dormand-Prince RK45 (for stiff-ish kinetic cascades such as the
 * MAPK model). Returns final state and number of accepted steps.
 */
export function rk45(f, y0, t0, t1, { rtol = 1e-6, atol = 1e-9, h0 = (t1 - t0) / 100, hmax = Infinity, maxSteps = 1e6, onStep } = {}) {
  const n = y0.length;
  let y = Float64Array.from(y0);
  const k = Array.from({ length: 7 }, () => new Float64Array(n));
  const tmp = new Float64Array(n), y5 = new Float64Array(n);
  let t = t0, h = Math.min(h0, hmax), accepted = 0;
  for (let step = 0; step < maxSteps && t < t1; step++) {
    if (t + h > t1) h = t1 - t;
    f(t, y, k[0]);
    for (let s = 1; s < 7; s++) {
      for (let i = 0; i < n; i++) {
        let acc = y[i];
        for (let j = 0; j < s; j++) acc += h * A[s][j] * k[j][i];
        tmp[i] = acc;
      }
      f(t + C[s] * h, tmp, k[s]);
    }
    let err = 0;
    for (let i = 0; i < n; i++) {
      let s5 = y[i], s4 = y[i];
      for (let s = 0; s < 7; s++) { s5 += h * B5[s] * k[s][i]; s4 += h * B4[s] * k[s][i]; }
      y5[i] = s5;
      const sc = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(s5));
      err = Math.max(err, Math.abs(s5 - s4) / sc);
    }
    if (err <= 1 || h < 1e-12) {
      t += h;
      y.set(y5);
      accepted++;
      onStep?.(t, y);
    }
    const factor = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * err ** -0.2));
    h = Math.min(h * factor, hmax);
  }
  return { t, y: Array.from(y), steps: accepted };
}

/**
 * Gillespie stochastic simulation (direct method).
 * @param {number[]} x0 initial copy numbers
 * @param {{rate:(x:number[])=>number, change:number[]}[]} reactions  change = stoichiometry vector
 */
export function gillespie(x0, reactions, tEnd, rng, { maxEvents = 1e6, sampleEvery = 0 } = {}) {
  const x = x0.slice();
  let t = 0, events = 0;
  const trace = sampleEvery ? [{ t: 0, x: x.slice() }] : null;
  let nextSample = sampleEvery;
  const props = new Float64Array(reactions.length);
  while (t < tEnd && events < maxEvents) {
    let a0 = 0;
    for (let r = 0; r < reactions.length; r++) { props[r] = Math.max(0, reactions[r].rate(x)); a0 += props[r]; }
    if (a0 <= 0) break;
    const tau = rng.exponential(a0);
    if (t + tau > tEnd) { t = tEnd; break; }
    t += tau;
    let u = rng.next() * a0, r = 0;
    for (; r < reactions.length - 1; r++) { u -= props[r]; if (u < 0) break; }
    const ch = reactions[r].change;
    for (let i = 0; i < ch.length; i++) x[i] += ch[i];
    events++;
    if (trace) while (t >= nextSample && nextSample <= tEnd) { trace.push({ t: nextSample, x: x.slice() }); nextSample += sampleEvery; }
  }
  return { t, x, events, trace };
}

/** Hill activation x^n / (K^n + x^n). */
export const hill = (x, K, n) => { const xn = Math.pow(Math.max(x, 0), n); return xn / (Math.pow(K, n) + xn); };
/** Hill repression K^n / (K^n + x^n). */
export const hillR = (x, K, n) => 1 - hill(x, K, n);
