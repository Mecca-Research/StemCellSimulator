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
