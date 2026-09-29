// Reaction-diffusion (Turing) morphogenesis (paper: "Morphogenesis - Turing
// patterns as a reaction-diffusion system":
//     du/dt = Du lap(u) + f(u, v),   dv/dt = Dv lap(v) + g(u, v)).
//
// Contents
//   1. Discrete Laplacians
//      - regular 2D / 3D grids (5- / 7-point stencils, zero-flux or periodic)
//      - triangle meshes: cotangent weights with the mixed Voronoi area of
//        Meyer, Desbrun, Schroder & Barr 2003 ("Discrete differential-geometry
//        operators for triangulated 2-manifolds", Visualization and
//        Mathematics III, Springer) and a uniform graph Laplacian fallback
//      - icosphere generator and real spherical harmonics (the exact
//        eigenfunctions of the Laplace-Beltrami operator on a sphere)
//   2. Kinetics + linear (Turing) stability analysis
//      - Schnakenberg 1979 (J Theor Biol 81:389), analysed as in Murray,
//        Mathematical Biology II (3rd ed. 2003), ch. 2
//      - Gray & Scott 1984 kinetics in the parameter map of Pearson 1993
//        (Science 261:189, "Complex patterns in a simple system")
//      - Gierer & Meinhardt 1972 (Kybernetik 12:30) activator-inhibitor with
//        optional activator saturation (Meinhardt 1982, Models of biological
//        pattern formation)
//   3. Time stepping: explicit Euler with a CFL-safe substep helper, and a
//      semi-implicit (IMEX) scheme - diffusion implicit via Jacobi-
//      preconditioned conjugate gradients, reactions explicit - for fine
//      meshes where a fast-diffusing inhibitor makes explicit steps tiny.
//
// Pure JS, no DOM / three.js. Arrays are Float64Array; everything that takes
// randomness takes an RNG (core/rng.js) so runs are reproducible.

// ============================================================ 1. Laplacians

/**
 * 5-point Laplacian on an nx x ny cell-centred grid (row-major, index = y*nx + x).
 * bc 'neumann' = zero flux (mirror ghost cells), 'periodic' = torus.
 */
export function laplacian2D(u, nx, ny, out, h = 1, bc = 'neumann') {
  const ih2 = 1 / (h * h);
  const per = bc === 'periodic';
  for (let y = 0; y < ny; y++) {
    const ym = y > 0 ? y - 1 : per ? ny - 1 : 0;
    const yp = y < ny - 1 ? y + 1 : per ? 0 : ny - 1;
    const row = y * nx, rm = ym * nx, rp = yp * nx;
    for (let x = 0; x < nx; x++) {
      const xm = x > 0 ? x - 1 : per ? nx - 1 : 0;
      const xp = x < nx - 1 ? x + 1 : per ? 0 : nx - 1;
      const c = u[row + x];
      out[row + x] = (u[row + xm] + u[row + xp] + u[rm + x] + u[rp + x] - 4 * c) * ih2;
    }
  }
  return out;
}

/** 7-point Laplacian on an nx x ny x nz grid (index = (z*ny + y)*nx + x). */
export function laplacian3D(u, nx, ny, nz, out, h = 1, bc = 'neumann') {
  const ih2 = 1 / (h * h);
  const per = bc === 'periodic';
  const sxy = nx * ny;
  for (let z = 0; z < nz; z++) {
    const zm = z > 0 ? z - 1 : per ? nz - 1 : 0;
    const zp = z < nz - 1 ? z + 1 : per ? 0 : nz - 1;
    for (let y = 0; y < ny; y++) {
      const ym = y > 0 ? y - 1 : per ? ny - 1 : 0;
      const yp = y < ny - 1 ? y + 1 : per ? 0 : ny - 1;
      const base = z * sxy + y * nx;
      for (let x = 0; x < nx; x++) {
        const xm = x > 0 ? x - 1 : per ? nx - 1 : 0;
        const xp = x < nx - 1 ? x + 1 : per ? 0 : nx - 1;
        const i = base + x;
        const c = u[i];
        out[i] = (u[base + xm] + u[base + xp]
          + u[z * sxy + ym * nx + x] + u[z * sxy + yp * nx + x]
          + u[zm * sxy + y * nx + x] + u[zp * sxy + y * nx + x] - 6 * c) * ih2;
      }
    }
  }
  return out;
}

/** Grid Laplacian wrapped as an operator { n, apply, rho } (rho = spectral-radius bound). */
export function gridOperator2D({ nx, ny, h = 1, bc = 'neumann' }) {
  return { n: nx * ny, nx, ny, h, bc, rho: 8 / (h * h), apply: (u, out) => laplacian2D(u, nx, ny, out, h, bc) };
}

export function gridOperator3D({ nx, ny, nz, h = 1, bc = 'neumann' }) {
  return { n: nx * ny * nz, nx, ny, nz, h, bc, rho: 12 / (h * h), apply: (u, out) => laplacian3D(u, nx, ny, nz, out, h, bc) };
}

/**
 * Geodesic icosphere: an icosahedron subdivided `level` times (each triangle
 * split into four, midpoints projected to the sphere). level 4 -> 2562
 * vertices, level 5 -> 10242 vertices / 20480 triangles. Shared (indexed)
 * vertices, counter-clockwise faces seen from outside.
 */
export function icosphere(level = 4, radius = 1) {
  const t = (1 + Math.sqrt(5)) / 2;
  let verts = [
    [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
    [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
    [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
  ].map((v) => { const n = Math.hypot(...v); return v.map((c) => c / n); });
  let faces = [
    [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
    [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
    [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
    [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
  ];
  for (let s = 0; s < level; s++) {
    const cache = new Map();
    const mid = (a, b) => {
      const key = a < b ? a * 1e6 + b : b * 1e6 + a;
      let m = cache.get(key);
      if (m === undefined) {
        const va = verts[a], vb = verts[b];
        const p = [va[0] + vb[0], va[1] + vb[1], va[2] + vb[2]];
        const n = Math.hypot(p[0], p[1], p[2]);
        m = verts.length;
        verts.push([p[0] / n, p[1] / n, p[2] / n]);
        cache.set(key, m);
      }
      return m;
    };
    const next = [];
    for (const [a, b, c] of faces) {
      const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
      next.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]);
    }
    faces = next;
  }
  const positions = new Float64Array(verts.length * 3);
  verts.forEach((v, i) => { positions[3 * i] = v[0] * radius; positions[3 * i + 1] = v[1] * radius; positions[3 * i + 2] = v[2] * radius; });
  const indices = new Uint32Array(faces.length * 3);
  faces.forEach((f, i) => { indices[3 * i] = f[0]; indices[3 * i + 1] = f[1]; indices[3 * i + 2] = f[2]; });
  return { positions, indices, n: verts.length, radius, level };
}

/**
 * Laplace-Beltrami operator of a triangle mesh in CSR form:
 *     (L u)_i = (1 / A_i) * sum_j w_ij (u_j - u_i)
 * type 'cotan':   w_ij = (cot alpha_ij + cot beta_ij) / 2, A_i = mixed Voronoi
 *                 area (Meyer et al. 2003, Fig. 4 / section 3.3)
 * type 'uniform': w_ij = 1, A_i = deg_i * <|e|^2>_i / 4, i.e. the graph
 *                 Laplacian scaled so it is exact for quadratics on a regular
 *                 triangular or square lattice (robust fallback for bad meshes)
 */
export function meshLaplacian(positions, indices, { type = 'cotan' } = {}) {
  const n = positions.length / 3;
  const nbr = Array.from({ length: n }, () => new Map());
  const area = new Float64Array(n);
  const P = (i, k) => positions[3 * i + k];
  const addW = (i, j, w) => {
    nbr[i].set(j, (nbr[i].get(j) ?? 0) + w);
    nbr[j].set(i, (nbr[j].get(i) ?? 0) + w);
  };
  let obtuse = 0;
  for (let f = 0; f < indices.length; f += 3) {
    const tri = [indices[f], indices[f + 1], indices[f + 2]];
    if (type === 'uniform') {
      for (let e = 0; e < 3; e++) {
        const a = tri[e], b = tri[(e + 1) % 3];
        if (!nbr[a].has(b)) { nbr[a].set(b, 1); nbr[b].set(a, 1); }
      }
      continue;
    }
    // edge vectors and cotangents of the three corner angles
    const cot = [0, 0, 0], ang = [0, 0, 0];
    const len2 = [0, 0, 0]; // squared length of the edge opposite corner k
    let dbl = 0;
    for (let k = 0; k < 3; k++) {
      const o = tri[k], a = tri[(k + 1) % 3], b = tri[(k + 2) % 3];
      const ux = P(a, 0) - P(o, 0), uy = P(a, 1) - P(o, 1), uz = P(a, 2) - P(o, 2);
      const vx = P(b, 0) - P(o, 0), vy = P(b, 1) - P(o, 1), vz = P(b, 2) - P(o, 2);
      const d = ux * vx + uy * vy + uz * vz;
      const cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx;
      const cr = Math.hypot(cx, cy, cz) || 1e-300;
      cot[k] = d / cr;
      ang[k] = Math.atan2(cr, d);
      dbl = cr;
      const wx = P(b, 0) - P(a, 0), wy = P(b, 1) - P(a, 1), wz = P(b, 2) - P(a, 2);
      len2[k] = wx * wx + wy * wy + wz * wz;
    }
    const triArea = dbl / 2;
    for (let k = 0; k < 3; k++) addW(tri[(k + 1) % 3], tri[(k + 2) % 3], cot[k] / 2);
    // mixed Voronoi area
    const isObtuse = ang.some((a) => a > Math.PI / 2);
    if (isObtuse) obtuse++;
    for (let k = 0; k < 3; k++) {
      const v = tri[k];
      if (!isObtuse) {
        // Voronoi region of v: (|v a|^2 cot(angle at b) + |v b|^2 cot(angle at a)) / 8
        const a = (k + 1) % 3, b = (k + 2) % 3;
        area[v] += (len2[a] * cot[a] + len2[b] * cot[b]) / 8; // len2[a] = |v b|^2 (opposite a), len2[b] = |v a|^2
      } else {
        area[v] += ang[k] > Math.PI / 2 ? triArea / 2 : triArea / 4;
      }
    }
  }
  // CSR assembly
  const rowPtr = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) rowPtr[i + 1] = rowPtr[i] + nbr[i].size;
  const col = new Int32Array(rowPtr[n]);
  const w = new Float64Array(rowPtr[n]);
  const wsum = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let e = rowPtr[i];
    for (const [j, wij] of nbr[i]) { col[e] = j; w[e] = wij; wsum[i] += wij; e++; }
  }
  if (type === 'uniform') {
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let e = rowPtr[i]; e < rowPtr[i + 1]; e++) {
        const j = col[e];
        s += (P(j, 0) - P(i, 0)) ** 2 + (P(j, 1) - P(i, 1)) ** 2 + (P(j, 2) - P(i, 2)) ** 2;
      }
      const deg = rowPtr[i + 1] - rowPtr[i];
      area[i] = (deg * (s / deg)) / 4;
    }
  }
  return new MeshLaplacian({ n, rowPtr, col, w, wsum, area, type, obtuse });
}

export class MeshLaplacian {
  constructor({ n, rowPtr, col, w, wsum, area, type, obtuse = 0 }) {
    Object.assign(this, { n, rowPtr, col, w, wsum, area, type, obtuse });
    this.invA = new Float64Array(n);
    let rho = 0;
    for (let i = 0; i < n; i++) {
      this.invA[i] = 1 / area[i];
      // Gershgorin bound on the spectrum of M^-1 K (K = stiffness): |lambda| <= 2 max_i sum_j |w_ij| / A_i
      let s = 0;
      for (let e = rowPtr[i]; e < rowPtr[i + 1]; e++) s += Math.abs(w[e]);
      rho = Math.max(rho, (2 * s) / area[i]);
    }
    this.rho = rho;
    this.totalArea = area.reduce((s, a) => s + a, 0);
    this._pcg = null;
  }

  /** out = L u */
  apply(u, out) {
    const { n, rowPtr, col, w, wsum, invA } = this;
    for (let i = 0; i < n; i++) {
      let s = -wsum[i] * u[i];
      for (let e = rowPtr[i]; e < rowPtr[i + 1]; e++) s += w[e] * u[col[e]];
      out[i] = s * invA[i];
    }
    return out;
  }

  /** Area-weighted mean of a vertex field. */
  mean(u) {
    let s = 0;
    for (let i = 0; i < this.n; i++) s += this.area[i] * u[i];
    return s / this.totalArea;
  }

  /**
   * Implicit diffusion step: solve (I - tau L) x = b, i.e. the symmetric
   * system (M + tau K) x = M b, by Jacobi-preconditioned conjugate gradients.
   * x holds the initial guess (warm start) and receives the solution.
   * @returns {number} iterations used
   */
  solveImplicit(tau, b, x, { tol = 1e-6, maxIter = 400 } = {}) {
    const { n, rowPtr, col, w, wsum, area } = this;
    if (!this._pcg || this._pcg.r.length !== n) {
      this._pcg = { r: new Float64Array(n), z: new Float64Array(n), p: new Float64Array(n), q: new Float64Array(n), d: new Float64Array(n) };
    }
    const { r, z, p, q, d } = this._pcg;
    const Ax = (v, out) => {
      for (let i = 0; i < n; i++) {
        let s = (area[i] + tau * wsum[i]) * v[i];
        for (let e = rowPtr[i]; e < rowPtr[i + 1]; e++) s -= tau * w[e] * v[col[e]];
        out[i] = s;
      }
    };
    for (let i = 0; i < n; i++) d[i] = 1 / (area[i] + tau * wsum[i]);
    Ax(x, q);
    let bnorm = 0, rz = 0;
    for (let i = 0; i < n; i++) {
      const bi = area[i] * b[i];
      bnorm += bi * bi;
      r[i] = bi - q[i];
      z[i] = d[i] * r[i];
      p[i] = z[i];
      rz += r[i] * z[i];
    }
    bnorm = Math.sqrt(bnorm) || 1;
    let it = 0;
    for (; it < maxIter; it++) {
      let rn = 0;
      for (let i = 0; i < n; i++) rn += r[i] * r[i];
      if (Math.sqrt(rn) <= tol * bnorm) break;
      Ax(p, q);
      let pq = 0;
      for (let i = 0; i < n; i++) pq += p[i] * q[i];
      const alpha = rz / pq;
      let rzNew = 0;
      for (let i = 0; i < n; i++) {
        x[i] += alpha * p[i];
        r[i] -= alpha * q[i];
        z[i] = d[i] * r[i];
        rzNew += r[i] * z[i];
      }
      const beta = rzNew / rz;
      rz = rzNew;
      for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    }
    return it;
  }
}

/** Vertex adjacency lists of an indexed triangle mesh. */
export function meshNeighbours(n, indices) {
  const sets = Array.from({ length: n }, () => new Set());
  for (let f = 0; f < indices.length; f += 3) {
    const a = indices[f], b = indices[f + 1], c = indices[f + 2];
    sets[a].add(b); sets[a].add(c); sets[b].add(a); sets[b].add(c); sets[c].add(a); sets[c].add(b);
  }
  return sets.map((s) => [...s]);
}

const _shCoef = new Map();
/** Recurrence coefficients of the normalised associated Legendre functions (cached per lmax). */
function shCoefficients(lmax) {
  let c = _shCoef.get(lmax);
  if (c) return c;
  const T = (lmax + 1) * (lmax + 2) / 2;
  const A = new Float64Array(T), B = new Float64Array(T), D = new Float64Array(lmax + 1), E = new Float64Array(lmax + 1);
  const idx = (l, m) => l * (l + 1) / 2 + m;
  for (let m = 1; m <= lmax; m++) D[m] = -Math.sqrt((2 * m + 1) / (2 * m));
  for (let m = 0; m < lmax; m++) E[m] = Math.sqrt(2 * m + 3);
  for (let m = 0; m <= lmax; m++) {
    for (let l = m + 2; l <= lmax; l++) {
      A[idx(l, m)] = Math.sqrt((4 * l * l - 1) / (l * l - m * m));
      B[idx(l, m)] = Math.sqrt(((l - 1) * (l - 1) - m * m) / (4 * (l - 1) * (l - 1) - 1));
    }
  }
  c = { A, B, D, E };
  _shCoef.set(lmax, c);
  return c;
}

/**
 * Real, orthonormal spherical harmonics Y_lm(theta, phi) for all l <= lmax at
 * one direction (x, y, z). Output is indexed l*l + l + m (m = -l..l).
 * Fully normalised associated Legendre recurrences (e.g. Press et al.,
 * Numerical Recipes 3rd ed., section 6.7).
 */
export function realSH(x, y, z, lmax, out = new Float64Array((lmax + 1) * (lmax + 1)), work = null) {
  const r = Math.hypot(x, y, z) || 1;
  const ct = z / r, st = Math.sqrt(Math.max(0, 1 - ct * ct));
  const phi = Math.atan2(y, x);
  const Pl = work ?? new Float64Array((lmax + 1) * (lmax + 1));
  const { A, B, D, E } = shCoefficients(lmax);
  const idx = (l, m) => l * (l + 1) / 2 + m; // triangular storage for m >= 0
  Pl[0] = Math.sqrt(1 / (4 * Math.PI));
  for (let m = 1; m <= lmax; m++) Pl[idx(m, m)] = D[m] * st * Pl[idx(m - 1, m - 1)];
  for (let m = 0; m < lmax; m++) Pl[idx(m + 1, m)] = E[m] * ct * Pl[idx(m, m)];
  for (let m = 0; m <= lmax; m++) {
    for (let l = m + 2; l <= lmax; l++) {
      const k = idx(l, m);
      Pl[k] = A[k] * (ct * Pl[idx(l - 1, m)] - B[k] * Pl[idx(l - 2, m)]);
    }
  }
  // cos(m phi), sin(m phi) by the angle-addition recurrence (two trig calls per direction)
  const c1 = Math.cos(phi), s1 = Math.sin(phi);
  let cm = 1, sm = 0;
  for (let l = 0; l <= lmax; l++) out[l * l + l] = Pl[idx(l, 0)];
  for (let m = 1; m <= lmax; m++) {
    const c = cm * c1 - sm * s1;
    sm = sm * c1 + cm * s1;
    cm = c;
    for (let l = m; l <= lmax; l++) {
      const p = Math.SQRT2 * Pl[idx(l, m)];
      out[l * l + l + m] = p * cm;
      out[l * l + l - m] = p * sm;
    }
  }
  return out;
}

/**
 * Angular power spectrum of a vertex field on a sphere-like mesh:
 * P_l = sum_m c_lm^2 with c_lm = sum_i A_i u_i Y_lm(x_i) / R^2 (mean removed).
 */
export function sphericalPowerSpectrum(u, positions, area, lmax) {
  const acc = new SpectrumAccumulator(positions, area, lmax);
  acc.begin(u);
  acc.step(Infinity);
  return acc.result();
}

/**
 * The same projection spread over several calls (e.g. a few hundred vertices
 * per animation frame) so a live view never stalls: begin(u) snapshots the
 * field, step(count) processes vertices and returns true when finished.
 */
export class SpectrumAccumulator {
  constructor(positions, area, lmax) {
    this.positions = positions;
    this.area = area;
    this.lmax = lmax;
    this.n = area.length;
    const K = (lmax + 1) * (lmax + 1);
    this.c = new Float64Array(K);
    this.Y = new Float64Array(K);
    this.work = new Float64Array(K);
    this.u = new Float64Array(this.n);
    this.i = this.n;
    let At = 0;
    for (let i = 0; i < this.n; i++) At += area[i];
    this.At = At;
    this.R2 = At / (4 * Math.PI);
  }

  begin(u) {
    let mean = 0;
    for (let i = 0; i < this.n; i++) { this.u[i] = u[i]; mean += this.area[i] * u[i]; }
    this.mean = mean / this.At;
    this.c.fill(0);
    this.i = 0;
  }

  get done() { return this.i >= this.n; }

  step(count) {
    const { positions: p, area, lmax, c, Y, work, u } = this;
    const K = c.length, end = Math.min(this.n, this.i + count);
    for (let i = this.i; i < end; i++) {
      realSH(p[3 * i], p[3 * i + 1], p[3 * i + 2], lmax, Y, work);
      const wv = (area[i] * (u[i] - this.mean)) / this.R2;
      for (let k = 0; k < K; k++) c[k] += wv * Y[k];
    }
    this.i = end;
    return this.done;
  }

  result() {
    const P = new Float64Array(this.lmax + 1);
    for (let l = 0; l <= this.lmax; l++) for (let m = -l; m <= l; m++) P[l] += this.c[l * l + l + m] ** 2;
    return P;
  }
}

// ============================================ 2. Kinetics + Turing analysis

/**
 * Linear stability of a homogeneous steady state (u*, v*) of
 *     u_t = Du lap u + f,  v_t = Dv lap v + g
 * to perturbations ~ exp(lambda t + i k.x): lambda(k^2) are the eigenvalues of
 *     A(k^2) = J - k^2 diag(Du, Dv),   J = [[fu, fv], [gu, gv]].
 * Turing (diffusion-driven) instability conditions (Murray 2003, vol. II, ch. 2):
 *     tr J < 0, det J > 0                   (stable without diffusion)
 *     Dv fu + Du gv > 0
 *     (Dv fu + Du gv)^2 > 4 Du Dv det J     (some k^2 with det A(k^2) < 0)
 * Unstable band: roots of h(k^2) = Du Dv k^4 - (Dv fu + Du gv) k^2 + det J.
 */
export function dispersion(J, Du, Dv, k2) {
  const [[fu, fv], [gu, gv]] = J;
  const a = fu - Du * k2, d = gv - Dv * k2;
  const tr = a + d, det = a * d - fv * gu;
  const disc = (tr * tr) / 4 - det;
  if (disc >= 0) {
    const s = Math.sqrt(disc);
    return { re: [tr / 2 + s, tr / 2 - s], im: [0, 0], max: tr / 2 + s };
  }
  const s = Math.sqrt(-disc);
  return { re: [tr / 2, tr / 2], im: [s, -s], max: tr / 2 };
}

export function turingAnalysis(J, Du, Dv) {
  const [[fu, fv], [gu, gv]] = J;
  const trJ = fu + gv, detJ = fu * gv - fv * gu;
  const homogeneousStable = trJ < 0 && detJ > 0;
  const s = Dv * fu + Du * gv;
  const discH = s * s - 4 * Du * Dv * detJ;
  const kc2 = s / (2 * Du * Dv); // minimum of h(k^2): critical wavenumber at onset
  let band = null;
  if (s > 0 && discH > 0) band = [(s - Math.sqrt(discH)) / (2 * Du * Dv), (s + Math.sqrt(discH)) / (2 * Du * Dv)];
  const turingUnstable = homogeneousStable && band !== null;
  const growth = (k2) => dispersion(J, Du, Dv, k2).max;
  // fastest-growing mode: dense scan then golden-section refinement
  const kHi = band ? band[1] * 1.5 : Math.max(Math.abs(kc2), 1) * 4;
  let best = 0, bestG = growth(0);
  const N = 400;
  for (let i = 1; i <= N; i++) { const k2 = (kHi * i) / N; const g = growth(k2); if (g > bestG) { bestG = g; best = k2; } }
  let lo = Math.max(0, best - kHi / N), hi = best + kHi / N;
  const gr = (Math.sqrt(5) - 1) / 2;
  for (let it = 0; it < 60; it++) {
    const m1 = hi - gr * (hi - lo), m2 = lo + gr * (hi - lo);
    if (growth(m1) > growth(m2)) hi = m2; else lo = m1;
  }
  const k2f = (lo + hi) / 2;
  const fastest = { k2: k2f, k: Math.sqrt(k2f), rate: growth(k2f), wavelength: k2f > 0 ? (2 * Math.PI) / Math.sqrt(k2f) : Infinity };
  return { trJ, detJ, homogeneousStable, turingUnstable, band, kc2, fastest, growth, J, Du, Dv };
}

/**
 * Accessible modes on a sphere of radius R: Laplace-Beltrami eigenvalues
 * k^2 = l (l + 1) / R^2 with multiplicity 2l + 1 (the finite-domain mode
 * selection of Murray 2003, vol. II, ch. 2, applied to a sphere). Returns
 * growth rates per l.
 */
export function sphereModes(analysis, R, lmax = 30) {
  const out = [];
  for (let l = 0; l <= lmax; l++) {
    const k2 = (l * (l + 1)) / (R * R);
    out.push({ l, k2, k: Math.sqrt(k2), rate: analysis.growth(k2) });
  }
  return out;
}

/** Schnakenberg (1979) "activator-substrate" kinetics, scaled by gamma (Murray's form). */
export const schnakenberg = {
  id: 'schnakenberg',
  label: 'Schnakenberg (activator-substrate)',
  species: ['activator u', 'substrate v'],
  // a = 0.1, b = 0.9 with Dv/Du = 40: a standard diffusion-driven-instability test set
  defaults: { a: 0.1, b: 0.9, gamma: 1, Du: 1, Dv: 40 },
  steadyStates(p) { const u = p.a + p.b; return [[u, p.b / (u * u)]]; },
  /** Jacobian at (u, v) (defaults to the homogeneous steady state). */
  jacobian(p, u = p.a + p.b, v = p.b / ((p.a + p.b) ** 2)) {
    const g = p.gamma ?? 1;
    return [[g * (-1 + 2 * u * v), g * u * u], [g * (-2 * u * v), g * (-u * u)]];
  },
  kernel(p) {
    const { a, b } = p, g = p.gamma ?? 1;
    return (u, v, fu, fv, n) => {
      for (let i = 0; i < n; i++) {
        const uu = u[i], u2v = uu * uu * v[i];
        fu[i] = g * (a - uu + u2v);
        fv[i] = g * (b - u2v);
      }
    };
  },
};

/**
 * Gray-Scott (Gray & Scott 1984; patterns catalogued by Pearson 1993 Science 261:189):
 *   u_t = Du lap u - u v^2 + F (1 - u)
 *   v_t = Dv lap v + u v^2 - (F + k) v
 * Presets sit in the regions of Pearson's (F, k) map (labels as in the
 * xmorphia catalogue of that map by R. Munafo) and were checked on the
 * level-5 icosphere: spots -> ~50 isolated blobs, stripes -> ~17 elongated
 * worms, labyrinth -> a few large connected components.
 */
export const grayScott = {
  id: 'grayScott',
  label: 'Gray-Scott (Pearson 1993)',
  species: ['substrate u', 'autocatalyst v'],
  defaults: { F: 0.0367, k: 0.0649, Du: 1, Dv: 0.5 },
  presets: {
    spots: { F: 0.0367, k: 0.0649, label: 'self-replicating spots (mitosis)' },
    stripes: { F: 0.046, k: 0.063, label: 'stripes / worms' },
    labyrinth: { F: 0.029, k: 0.057, label: 'labyrinth (maze)' },
  },
  steadyStates(p) {
    const out = [[1, 0]];
    const disc = 1 - (4 * (p.F + p.k) ** 2) / p.F;
    if (disc >= 0) {
      for (const sg of [1, -1]) {
        const u = (1 + sg * Math.sqrt(disc)) / 2;
        out.push([u, (p.F + p.k) / u]);
      }
    }
    return out;
  },
  jacobian(p, u = 1, v = 0) {
    return [[-v * v - p.F, -2 * u * v], [v * v, 2 * u * v - (p.F + p.k)]];
  },
  kernel(p) {
    const { F, k } = p;
    return (u, v, fu, fv, n) => {
      for (let i = 0; i < n; i++) {
        const uvv = u[i] * v[i] * v[i];
        fu[i] = -uvv + F * (1 - u[i]);
        fv[i] = uvv - (F + k) * v[i];
      }
    };
  },
};

/**
 * Gierer-Meinhardt activator-inhibitor (1972), with activator saturation kappa:
 *   a_t = Da lap a + rhoA a^2 / (h (1 + kappa a^2)) - muA a + sigmaA
 *   h_t = Dh lap h + rhoH a^2 - muH h
 */
export const giererMeinhardt = {
  id: 'giererMeinhardt',
  label: 'Gierer-Meinhardt (activator-inhibitor)',
  species: ['activator a', 'inhibitor h'],
  defaults: { rhoA: 1, rhoH: 1, muA: 1, muH: 2, sigmaA: 0.0, kappa: 0.0, Du: 1, Dv: 30 },
  steadyStates(p) {
    // eliminate h* = rhoH a^2 / muH: phi(a) = rhoA muH / (rhoH (1 + kappa a^2)) - muA a + sigmaA = 0 (monotone)
    const phi = (a) => (p.rhoA * p.muH) / (p.rhoH * (1 + p.kappa * a * a)) - p.muA * a + p.sigmaA;
    let lo = 1e-9, hi = (p.rhoA * p.muH / p.rhoH + p.sigmaA) / p.muA + 1;
    for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (phi(m) > 0) lo = m; else hi = m; }
    const a = (lo + hi) / 2;
    return [[a, (p.rhoH * a * a) / p.muH]];
  },
  jacobian(p, a, h) {
    if (a === undefined) [a, h] = giererMeinhardt.steadyStates(p)[0];
    const s = 1 + p.kappa * a * a;
    return [
      [(2 * p.rhoA * a) / (h * s * s) - p.muA, -(p.rhoA * a * a) / (h * h * s)],
      [2 * p.rhoH * a, -p.muH],
    ];
  },
  kernel(p) {
    const { rhoA, rhoH, muA, muH, sigmaA, kappa } = p;
    return (a, h, fa, fh, n) => {
      for (let i = 0; i < n; i++) {
        const ai = a[i], hi = Math.max(h[i], 1e-6), a2 = ai * ai;
        fa[i] = (rhoA * a2) / (hi * (1 + kappa * a2)) - muA * ai + sigmaA;
        fh[i] = rhoH * a2 - muH * h[i];
      }
    };
  },
};

export const MODELS = { schnakenberg, grayScott, giererMeinhardt };

/** Turing analysis of a model's first (or given) homogeneous steady state. */
export function analyseModel(model, p, Du = p.Du, Dv = p.Dv, which = 0) {
  const ss = model.steadyStates(p);
  const [u, v] = ss[Math.min(which, ss.length - 1)];
  return { steady: [u, v], ...turingAnalysis(model.jacobian(p, u, v), Du, Dv) };
}

// ======================================================= 3. Time stepping

/**
 * CFL-safe explicit Euler step size for u_t = D L u + R(u):
 * diffusion needs dt * Dmax * rho(L) <= 2, reactions dt * |J| <= 1.
 */
export function stableDt({ Dmax, rho, reactionRate = 0, safety = 0.8 }) {
  const dDiff = Dmax > 0 ? 2 / (Dmax * rho) : Infinity;
  const dReact = reactionRate > 0 ? 1 / reactionRate : Infinity;
  return safety * Math.min(dDiff, dReact);
}

/** Split an interval T into n equal substeps no larger than dtMax. */
export function substeps(T, dtMax) {
  if (!(T > 0)) return { n: 0, dt: 0 };
  const n = Math.max(1, Math.ceil(T / dtMax - 1e-12));
  return { n, dt: T / n };
}

/**
 * Pick the cheaper stable scheme: explicit Euler costs ~2 operator
 * applications per step, IMEX ~12 (two warm-started PCG solves).
 */
export function chooseScheme(op, model, params) {
  if (!op.solveImplicit) return 'explicit';
  const probe = new RDSolver({ op: { n: 1, rho: op.rho, apply() {} }, model, params });
  const dExp = probe.maxDt();
  probe.scheme = 'imex';
  return dExp * 6 < probe.maxDt() ? 'imex' : 'explicit';
}

/**
 * Two-species reaction-diffusion solver on any Laplacian operator
 * ({ n, apply(u, out), rho } - a grid operator or a MeshLaplacian).
 * scheme 'explicit': forward Euler. scheme 'imex': reactions explicit,
 * diffusion backward Euler (requires op.solveImplicit, i.e. a mesh).
 */
export class RDSolver {
  constructor({ op, model, params, scheme = 'explicit' }) {
    this.op = op;
    this.n = op.n;
    this.model = model;
    this.scheme = scheme;
    this.u = new Float64Array(this.n);
    this.v = new Float64Array(this.n);
    this.fu = new Float64Array(this.n);
    this.fv = new Float64Array(this.n);
    this.lu = new Float64Array(this.n);
    this.lv = new Float64Array(this.n);
    this.t = 0;
    this.iters = 0;
    this.tol = 1e-4; // relative residual of the implicit solves (IMEX)
    this.setParams(params);
  }

  setParams(params) {
    this.p = { ...this.model.defaults, ...params };
    this.react = this.model.kernel(this.p);
    const J = this.model.jacobian(this.p, ...this.model.steadyStates(this.p)[0]);
    // reaction stiffness estimate: row-sum norm of J at the steady state (x2 for nonlinear excursions)
    this.reactionRate = 2 * Math.max(Math.abs(J[0][0]) + Math.abs(J[0][1]), Math.abs(J[1][0]) + Math.abs(J[1][1]), 1e-6);
  }

  /** Homogeneous steady state plus uniform noise of relative amplitude amp. */
  seed(rng, amp = 0.01, which = 0) {
    const [us, vs] = this.model.steadyStates(this.p)[which];
    for (let i = 0; i < this.n; i++) {
      this.u[i] = us * (1 + amp * (rng.next() * 2 - 1));
      this.v[i] = vs * (1 + amp * (rng.next() * 2 - 1));
    }
    this.t = 0;
  }

  maxDt(safety = 0.8) {
    // IMEX: diffusion is unconditionally stable, the step only has to resolve the kinetics
    if (this.scheme === 'imex') return (0.75 * safety) / this.reactionRate;
    return stableDt({ Dmax: Math.max(this.p.Du, this.p.Dv), rho: this.op.rho, reactionRate: this.reactionRate, safety });
  }

  step(dt) {
    const { n, u, v, fu, fv, lu, lv, p } = this;
    this.react(u, v, fu, fv, n);
    if (this.scheme === 'imex') {
      for (let i = 0; i < n; i++) { lu[i] = u[i] + dt * fu[i]; lv[i] = v[i] + dt * fv[i]; }
      // warm start from the reaction-updated field
      u.set(lu); v.set(lv);
      this.iters = this.op.solveImplicit(dt * p.Du, lu, u, { tol: this.tol })
        + this.op.solveImplicit(dt * p.Dv, lv, v, { tol: this.tol });
    } else {
      this.op.apply(u, lu);
      this.op.apply(v, lv);
      const Du = p.Du, Dv = p.Dv;
      for (let i = 0; i < n; i++) {
        u[i] += dt * (Du * lu[i] + fu[i]);
        v[i] += dt * (Dv * lv[i] + fv[i]);
      }
    }
    // non-negative concentrations; flush tiny values (denormals are very slow)
    for (let i = 0; i < n; i++) { if (!(u[i] > 1e-30)) u[i] = 0; if (!(v[i] > 1e-30)) v[i] = 0; }
    this.t += dt;
  }

  /** Advance by T with CFL-safe substeps; returns the number of substeps taken. */
  advance(T, { safety = 0.8, maxSteps = Infinity } = {}) {
    const { n, dt } = substeps(T, this.maxDt(safety));
    const m = Math.min(n, maxSteps);
    for (let s = 0; s < m; s++) this.step(dt);
    return m;
  }
}

// =================================================== pattern measurement

/** Spatial variance of a field. */
export function variance(u) {
  let s = 0, s2 = 0;
  for (let i = 0; i < u.length; i++) { s += u[i]; s2 += u[i] * u[i]; }
  const m = s / u.length;
  return s2 / u.length - m * m;
}

/**
 * Dominant wavelength of a 2D periodic field from the peak of its radially
 * averaged power spectrum (separable DFT; fine for the small test grids).
 */
export function dominantWavelength2D(u, nx, ny, h = 1) {
  let mean = 0;
  for (let i = 0; i < u.length; i++) mean += u[i];
  mean /= u.length;
  // DFT along x for each row
  const re1 = new Float64Array(nx * ny), im1 = new Float64Array(nx * ny);
  for (let y = 0; y < ny; y++) {
    for (let kx = 0; kx < nx; kx++) {
      let sr = 0, si = 0;
      for (let x = 0; x < nx; x++) {
        const a = (-2 * Math.PI * kx * x) / nx, val = u[y * nx + x] - mean;
        sr += val * Math.cos(a); si += val * Math.sin(a);
      }
      re1[y * nx + kx] = sr; im1[y * nx + kx] = si;
    }
  }
  const nb = Math.floor(Math.min(nx, ny) / 2);
  const bins = new Float64Array(nb + 1), counts = new Float64Array(nb + 1);
  for (let kx = 0; kx < nx; kx++) {
    for (let ky = 0; ky < ny; ky++) {
      let sr = 0, si = 0;
      for (let y = 0; y < ny; y++) {
        const a = (-2 * Math.PI * ky * y) / ny, c = Math.cos(a), s = Math.sin(a);
        const r = re1[y * nx + kx], i = im1[y * nx + kx];
        sr += r * c - i * s; si += r * s + i * c;
      }
      const fx = kx <= nx / 2 ? kx / nx : (kx - nx) / nx;
      const fy = ky <= ny / 2 ? ky / ny : (ky - ny) / ny;
      const f = Math.hypot(fx, fy) * Math.min(nx, ny); // in units of the fundamental
      const b = Math.round(f);
      if (b >= 1 && b <= nb) { bins[b] += sr * sr + si * si; counts[b]++; }
    }
  }
  let best = 1;
  for (let b = 1; b <= nb; b++) if (bins[b] > bins[best]) best = b;
  // parabolic refinement of the peak
  let bf = best;
  if (best > 1 && best < nb) {
    const y0 = bins[best - 1], y1 = bins[best], y2 = bins[best + 1];
    const den = y0 - 2 * y1 + y2;
    if (den < 0) bf = best + 0.5 * (y0 - y2) / den;
  }
  const L = Math.min(nx, ny) * h;
  return { wavelength: L / bf, k: (2 * Math.PI * bf) / L, spectrum: bins };
}
