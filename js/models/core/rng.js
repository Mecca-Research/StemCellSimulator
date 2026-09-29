// Deterministic, seedable random numbers (mulberry32) so every simulation run
// and every test is reproducible.
export class RNG {
  constructor(seed = 1) {
    this.s = seed >>> 0 || 1;
    this._spare = null;
  }

  /** uniform [0, 1) */
  next() {
    let t = (this.s = (this.s + 0x6d2b79f5) >>> 0);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  uniform(a = 0, b = 1) { return a + (b - a) * this.next(); }

  int(n) { return Math.floor(this.next() * n); }

  /** standard normal (Marsaglia polar) */
  normal(mu = 0, sigma = 1) {
    if (this._spare !== null) { const v = this._spare; this._spare = null; return mu + sigma * v; }
    let u, v, s;
    do { u = this.next() * 2 - 1; v = this.next() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt((-2 * Math.log(s)) / s);
    this._spare = v * m;
    return mu + sigma * u * m;
  }

  exponential(rate) { return -Math.log(1 - this.next()) / rate; }

  /** Poisson sample (Knuth for small means, normal approximation above 50). */
  poisson(mean) {
    if (mean <= 0) return 0;
    if (mean > 50) return Math.max(0, Math.round(this.normal(mean, Math.sqrt(mean))));
    const L = Math.exp(-mean);
    let k = 0, p = 1;
    do { k++; p *= this.next(); } while (p > L);
    return k - 1;
  }

  /** Uniform random unit vector in 3D. */
  onSphere(out = [0, 0, 0]) {
    const z = this.uniform(-1, 1), a = this.uniform(0, Math.PI * 2), r = Math.sqrt(1 - z * z);
    out[0] = r * Math.cos(a); out[1] = r * Math.sin(a); out[2] = z;
    return out;
  }

  pick(arr) { return arr[this.int(arr.length)]; }

  /** Choose an index with probability proportional to weights. */
  weighted(weights) {
    let total = 0;
    for (const w of weights) total += w;
    let r = this.next() * total;
    for (let i = 0; i < weights.length; i++) { r -= weights[i]; if (r < 0) return i; }
    return weights.length - 1;
  }
}
