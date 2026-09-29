// Parametric and data-driven geometries for cells and molecules.
import * as THREE from 'three';
import { mergeGeometries, mergeVertices } from 'three/addons/utils/BufferGeometryUtils.js';

/** Smooth unit sphere suitable for instanced ellipsoids. */
export function unitSphere(detail = 3) {
  return new THREE.IcosahedronGeometry(1, detail);
}

/**
 * Human erythrocyte: Evans & Fung (1972) biconcave profile
 *   z(r) = +/- (D/2) sqrt(1 - x^2) (c0 + c1 x^2 + c2 x^4),  x = 2r/D
 * with D = 7.82 um, c0 = 0.207, c1 = 2.003, c2 = -1.123.
 */
export function erythrocyteGeometry(D = 7.82, segments = 48) {
  const c0 = 0.207, c1 = 2.003, c2 = -1.123;
  const pts = [];
  const N = 40;
  const z = (x) => (D / 2) * Math.sqrt(Math.max(1 - x * x, 0)) * (c0 + c1 * x * x + c2 * x ** 4);
  // lower surface from the axis to the rim, then the upper surface back to the axis
  for (let i = 0; i <= N; i++) { const x = Math.sin((i / N) * Math.PI / 2); pts.push(new THREE.Vector2((x * D) / 2 + 1e-4, -z(x))); }
  for (let i = N; i >= 0; i--) { const x = Math.sin((i / N) * Math.PI / 2); pts.push(new THREE.Vector2((x * D) / 2 + 1e-4, z(x))); }
  const g = new THREE.LatheGeometry(pts, segments);
  g.computeVertexNormals();
  return g;
}

/** Multi-lobed (segmented) neutrophil nucleus: n spheres joined by thin bridges. */
export function lobedNucleus(n = 4, r = 1, rng = Math.random) {
  const parts = [];
  let prev = new THREE.Vector3();
  const dir = new THREE.Vector3(1, 0, 0);
  for (let i = 0; i < n; i++) {
    const s = new THREE.IcosahedronGeometry(r * (0.75 + 0.2 * rng()), 2);
    const p = prev.clone();
    s.translate(p.x, p.y, p.z);
    parts.push(s);
    dir.applyAxisAngle(new THREE.Vector3(0, 0, 1), (rng() - 0.3) * 1.4).applyAxisAngle(new THREE.Vector3(0, 1, 0), (rng() - 0.5) * 1.2);
    const next = p.clone().add(dir.clone().multiplyScalar(r * 1.25));
    if (i < n - 1) {
      const bridge = new THREE.CylinderGeometry(r * 0.18, r * 0.18, r * 1.25, 8);
      bridge.translate(0, r * 0.625, 0);
      bridge.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize()));
      bridge.translate(p.x, p.y, p.z);
      parts.push(bridge);
    }
    prev = next;
  }
  const g = mergeGeometries(parts.map((x) => x.toNonIndexed()));
  g.center();
  g.computeVertexNormals();
  return g;
}

/** Kidney-shaped (monocyte) nucleus: a sphere with an indentation. */
export function indentedNucleus(r = 1) {
  const g = new THREE.IcosahedronGeometry(r, 4);
  const p = g.attributes.position;
  const v = new THREE.Vector3();
  for (let i = 0; i < p.count; i++) {
    v.fromBufferAttribute(p, i);
    const d = Math.exp(-((v.x - r) ** 2 + v.z * v.z) / (0.35 * r * r));
    v.x -= d * 0.55 * r;
    v.y *= 0.8;
    p.setXYZ(i, v.x, v.y, v.z);
  }
  g.computeVertexNormals();
  return g;
}

/** Organic, slightly lumpy blob (for folded globular proteins / debris). */
export function blobGeometry(r = 1, amp = 0.18, seed = 1, detail = 3) {
  const g = new THREE.IcosahedronGeometry(r, detail);
  const p = g.attributes.position;
  const v = new THREE.Vector3();
  const f = (x, y, z) => Math.sin(x * 3.1 + seed) * Math.cos(y * 2.7 + seed * 1.7) * Math.sin(z * 3.4 + seed * 0.3);
  for (let i = 0; i < p.count; i++) {
    v.fromBufferAttribute(p, i);
    const n = v.clone().normalize();
    const k = 1 + amp * f(n.x, n.y, n.z) + amp * 0.5 * f(n.y * 2, n.z * 2, n.x * 2);
    v.copy(n).multiplyScalar(r * k);
    p.setXYZ(i, v.x, v.y, v.z);
  }
  g.computeVertexNormals();
  return g;
}

/**
 * Compact random-walk polypeptide chain packed inside a radius, as a tube:
 * gives a folded-globin look (alpha-helical segments are drawn as coils).
 */
export function foldedChainGeometry(rng, { length = 60, radius = 2.2, tube = 0.28, helixFrac = 0.7 } = {}) {
  const pts = [];
  const p = new THREE.Vector3();
  const d = new THREE.Vector3(1, 0, 0);
  for (let i = 0; i < length; i++) {
    const helical = (i % 10) < helixFrac * 10;
    const step = helical ? 0.35 : 0.8;
    d.add(new THREE.Vector3(rng() - 0.5, rng() - 0.5, rng() - 0.5).multiplyScalar(helical ? 0.6 : 1.4)).normalize();
    p.add(d.clone().multiplyScalar(step));
    if (p.length() > radius) { p.multiplyScalar(radius / p.length()); d.multiplyScalar(-1); }
    if (helical) {
      const a = i * 1.7;
      pts.push(p.clone().add(new THREE.Vector3(Math.cos(a), Math.sin(a), 0).multiplyScalar(0.25)));
    } else pts.push(p.clone());
  }
  const curve = new THREE.CatmullRomCurve3(pts);
  return new THREE.TubeGeometry(curve, length * 4, tube, 6, false);
}

/**
 * One chain of the collagen triple helix: a left-handed polyproline-II helix
 * wound around the common right-handed super-helix axis (x).
 */
export function collagenStrandCurve(phase, { length = 30, radius = 0.35, turns = 6 } = {}) {
  class Strand extends THREE.Curve {
    getPoint(t, out = new THREE.Vector3()) {
      const a = t * turns * Math.PI * 2 + phase;
      return out.set(t * length - length / 2, Math.cos(a) * radius, Math.sin(a) * radius);
    }
  }
  return new Strand();
}

/** Decode the quantised hiPSC mesh buffer written by tools/build_assets.py. */
export function decodeMesh(buffer, meta, meshInfo) {
  const { vertex_offset: vo, vertex_count: vc, index_offset: io, index_count: ic } = meshInfo;
  const q = new Uint16Array(buffer, vo, vc * 3);
  const idx = new Uint16Array(buffer, io, ic);
  const [ox, oy, oz] = meta.meshes.quant_origin_um;
  const [bx, by, bz] = meta.meshes.quant_box_um;
  const pos = new Float32Array(vc * 3);
  for (let i = 0; i < vc; i++) {
    pos[i * 3] = ox + (q[i * 3] / 65535) * bx;
    pos[i * 3 + 1] = oy + (q[i * 3 + 1] / 65535) * by;
    pos[i * 3 + 2] = oz + (q[i * 3 + 2] / 65535) * bz;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setIndex(new THREE.BufferAttribute(new Uint16Array(idx), 1));
  g.computeVertexNormals();
  return g;
}

export { mergeGeometries, mergeVertices };

/**
 * Tube with fixed topology whose centre-line can be re-shaped every frame
 * without allocating (parallel-transport frames). Used for chains that fold,
 * wind or grow.
 */
export class DynamicTube {
  constructor(segments = 120, radial = 8, radius = 0.3) {
    this.segments = segments;
    this.radial = radial;
    this.radius = radius;
    const g = new THREE.BufferGeometry();
    const nv = (segments + 1) * radial;
    this.pos = new Float32Array(nv * 3);
    this.nrm = new Float32Array(nv * 3);
    g.setAttribute('position', new THREE.BufferAttribute(this.pos, 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute('normal', new THREE.BufferAttribute(this.nrm, 3).setUsage(THREE.DynamicDrawUsage));
    const idx = [];
    for (let i = 0; i < segments; i++) {
      for (let j = 0; j < radial; j++) {
        const a = i * radial + j, b = i * radial + ((j + 1) % radial), c = (i + 1) * radial + j, d = (i + 1) * radial + ((j + 1) % radial);
        idx.push(a, c, b, b, c, d);
      }
    }
    g.setIndex(idx);
    this.geometry = g;
    this.pts = Array.from({ length: segments + 1 }, () => new THREE.Vector3());
    this._t = new THREE.Vector3(); this._n = new THREE.Vector3(); this._b = new THREE.Vector3(); this._prev = new THREE.Vector3(0, 1, 0);
  }

  /** fn(s in [0,1], out Vector3) writes the centre-line point; radiusFn(s) optional. */
  update(fn, radiusFn = null) {
    const { segments, radial, pts, pos, nrm } = this;
    for (let i = 0; i <= segments; i++) fn(i / segments, pts[i]);
    const t = this._t, n = this._n, b = this._b;
    n.copy(this._prev);
    for (let i = 0; i <= segments; i++) {
      const p0 = pts[Math.max(i - 1, 0)], p1 = pts[Math.min(i + 1, segments)];
      t.subVectors(p1, p0);
      if (t.lengthSq() < 1e-12) t.set(1, 0, 0);
      t.normalize();
      // parallel transport: remove the tangential part of the previous normal
      n.addScaledVector(t, -n.dot(t));
      if (n.lengthSq() < 1e-8) n.set(-t.y, t.x, 0).lengthSq() < 1e-8 ? n.set(0, -t.z, t.y) : null;
      n.normalize();
      b.crossVectors(t, n);
      const r = radiusFn ? radiusFn(i / segments) : this.radius;
      for (let j = 0; j < radial; j++) {
        const a = (j / radial) * Math.PI * 2, ca = Math.cos(a), sa = Math.sin(a);
        const k = (i * radial + j) * 3;
        const nx = ca * n.x + sa * b.x, ny = ca * n.y + sa * b.y, nz = ca * n.z + sa * b.z;
        nrm[k] = nx; nrm[k + 1] = ny; nrm[k + 2] = nz;
        pos[k] = pts[i].x + r * nx; pos[k + 1] = pts[i].y + r * ny; pos[k + 2] = pts[i].z + r * nz;
      }
    }
    this.geometry.attributes.position.needsUpdate = true;
    this.geometry.attributes.normal.needsUpdate = true;
    this.geometry.computeBoundingSphere();
  }

  /** Reveal only the first fraction f of the tube (e.g. a chain being translated). */
  reveal(f) {
    const n = Math.round(Math.min(Math.max(f, 0), 1) * this.segments) * this.radial * 6;
    this.geometry.setDrawRange(0, n);
  }
}
