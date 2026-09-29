// Small dense linear algebra: symmetric eigendecomposition (Jacobi) and PCA,
// used for the reaction-network embedding (the paper's "dimensionality
// reduction: PCA" step) and for fitting cell ellipsoids.

/** Jacobi eigenvalue algorithm for a symmetric matrix (array of rows). */
export function eigSym(M, { maxSweeps = 100, tol = 1e-12 } = {}) {
  const n = M.length;
  const a = M.map((r) => Float64Array.from(r));
  const v = Array.from({ length: n }, (_, i) => { const r = new Float64Array(n); r[i] = 1; return r; });
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
    if (off < tol) break;
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) {
          const akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq;
          v[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  const values = a.map((r, i) => r[i]);
  const order = values.map((_, i) => i).sort((i, j) => values[j] - values[i]);
  return {
    values: order.map((i) => values[i]),
    vectors: order.map((i) => v.map((r) => r[i])), // vectors[k] = k-th eigenvector
  };
}

/**
 * Principal component analysis.
 * @param {number[][]} X rows = observations
 * @param {number} k components
 * @returns {{scores:number[][], components:number[][], explained:number[]}}
 */
export function pca(X, k = 3) {
  const n = X.length, d = X[0].length;
  const mean = new Float64Array(d);
  for (const r of X) for (let j = 0; j < d; j++) mean[j] += r[j] / n;
  const C = Array.from({ length: d }, () => new Float64Array(d));
  for (const r of X) {
    for (let i = 0; i < d; i++) {
      const ri = r[i] - mean[i];
      if (ri === 0) continue;
      for (let j = i; j < d; j++) C[i][j] += ri * (r[j] - mean[j]);
    }
  }
  for (let i = 0; i < d; i++) for (let j = i; j < d; j++) { C[i][j] /= Math.max(n - 1, 1); C[j][i] = C[i][j]; }
  const { values, vectors } = eigSym(C);
  const total = values.reduce((s, v) => s + Math.max(v, 0), 0) || 1;
  const comps = vectors.slice(0, k);
  const scores = X.map((r) => comps.map((c) => c.reduce((s, cj, j) => s + cj * (r[j] - mean[j]), 0)));
  return { scores, components: comps, explained: values.slice(0, k).map((v) => Math.max(v, 0) / total) };
}

export const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const len3 = (a) => Math.hypot(a[0], a[1], a[2]);

/**
 * PCA by power iteration with deflation (top-k components only). Suitable for
 * a few hundred features where a full Jacobi decomposition would be slow.
 * @returns {{scores:number[][], components:Float64Array[], explained:number[]}}
 */
export function powerPCA(X, k = 3, { iters = 300, seed = 1 } = {}) {
  const n = X.length, d = X[0].length;
  const mean = new Float64Array(d);
  for (const r of X) for (let j = 0; j < d; j++) mean[j] += r[j] / n;
  const Xc = X.map((r) => Float64Array.from(r, (v, j) => v - mean[j]));
  // covariance C = Xc^T Xc / (n - 1), applied implicitly
  const apply = (v, out) => {
    const t = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = 0; const r = Xc[i]; for (let j = 0; j < d; j++) s += r[j] * v[j]; t[i] = s; }
    out.fill(0);
    for (let i = 0; i < n; i++) { const r = Xc[i], ti = t[i]; if (ti) for (let j = 0; j < d; j++) out[j] += r[j] * ti; }
    for (let j = 0; j < d; j++) out[j] /= Math.max(n - 1, 1);
  };
  let total = 0;
  for (const r of Xc) for (let j = 0; j < d; j++) total += (r[j] * r[j]) / Math.max(n - 1, 1);
  let s = seed;
  const rand = () => ((s = (s * 16807) % 2147483647) / 2147483647) - 0.5;
  const comps = [], vals = [];
  const w = new Float64Array(d);
  for (let c = 0; c < k; c++) {
    let v = Float64Array.from({ length: d }, rand);
    for (let it = 0; it < iters; it++) {
      apply(v, w);
      for (const u of comps) { let p = 0; for (let j = 0; j < d; j++) p += w[j] * u[j]; for (let j = 0; j < d; j++) w[j] -= p * u[j]; }
      let nrm = 0; for (let j = 0; j < d; j++) nrm += w[j] * w[j];
      nrm = Math.sqrt(nrm) || 1;
      for (let j = 0; j < d; j++) v[j] = w[j] / nrm;
    }
    apply(v, w);
    let lam = 0; for (let j = 0; j < d; j++) lam += w[j] * v[j];
    comps.push(v); vals.push(lam);
  }
  const scores = Xc.map((r) => comps.map((u) => { let p = 0; for (let j = 0; j < d; j++) p += r[j] * u[j]; return p; }));
  return { scores, components: comps, explained: vals.map((v) => (total ? v / total : 0)) };
}
