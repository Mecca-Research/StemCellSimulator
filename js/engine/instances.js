// Growable instanced-mesh pool: write transforms and colours per frame, the
// pool resizes itself (doubling) when a simulation outgrows it.
import * as THREE from 'three';

const _m = new THREE.Matrix4();
const _q = new THREE.Quaternion();
const _s = new THREE.Vector3();
const _p = new THREE.Vector3();
const _c = new THREE.Color();
const _up = new THREE.Vector3(1, 0, 0);
const _a = new THREE.Vector3();

export class InstancePool {
  constructor(geometry, material, capacity = 256, parent = null) {
    this.geometry = geometry;
    this.material = material;
    this.parent = parent;
    this.n = 0;
    this.mesh = null;
    this.alloc(capacity);
  }

  alloc(capacity) {
    const old = this.mesh;
    const mesh = new THREE.InstancedMesh(this.geometry, this.material, capacity);
    mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(capacity * 3).fill(1), 3);
    mesh.instanceColor.setUsage(THREE.DynamicDrawUsage);
    mesh.frustumCulled = false;
    mesh.count = 0;
    if (old) {
      mesh.renderOrder = old.renderOrder;
      mesh.visible = old.visible;
      mesh.userData = old.userData;
      old.parent?.remove(old);
      old.dispose();
    }
    this.mesh = mesh;
    this.capacity = capacity;
    this.parent?.add(mesh);
    return mesh;
  }

  begin() { this.n = 0; }

  /** Ensure room for one more instance; returns its index. */
  next() {
    if (this.n >= this.capacity) {
      const keepM = this.mesh.instanceMatrix.array.slice(0, this.n * 16);
      const keepC = this.mesh.instanceColor.array.slice(0, this.n * 3);
      this.alloc(this.capacity * 2);
      this.mesh.instanceMatrix.array.set(keepM);
      this.mesh.instanceColor.array.set(keepC);
    }
    return this.n++;
  }

  /** Axis-aligned ellipsoid/sphere instance. */
  put(x, y, z, sx, sy = sx, sz = sx, color = null) {
    const i = this.next();
    _m.makeScale(sx, sy, sz).setPosition(x, y, z);
    this.mesh.setMatrixAt(i, _m);
    if (color) this.mesh.setColorAt(i, color.isColor ? color : _c.set(color));
    return i;
  }

  /** Instance whose local +x axis points along `axis` (unit vector). */
  putOriented(pos, axis, sAlong, sPerp1, sPerp2 = sPerp1, color = null) {
    const i = this.next();
    _a.set(axis[0], axis[1], axis[2]);
    if (_a.lengthSq() < 1e-12) _a.set(1, 0, 0);
    _q.setFromUnitVectors(_up, _a.normalize());
    _s.set(sAlong, sPerp1, sPerp2);
    _p.set(pos[0], pos[1], pos[2]);
    _m.compose(_p, _q, _s);
    this.mesh.setMatrixAt(i, _m);
    if (color) this.mesh.setColorAt(i, color.isColor ? color : _c.set(color));
    return i;
  }

  /** Full matrix instance. */
  putMatrix(m, color = null) {
    const i = this.next();
    this.mesh.setMatrixAt(i, m);
    if (color) this.mesh.setColorAt(i, color.isColor ? color : _c.set(color));
    return i;
  }

  end() {
    this.mesh.count = this.n;
    this.mesh.instanceMatrix.needsUpdate = true;
    if (this.mesh.instanceColor) this.mesh.instanceColor.needsUpdate = true;
    // keep bounding sphere valid for raycasting
    this.mesh.boundingSphere = null;
  }
}

/** Hue helper: distinct, stable colours for integer ids (golden-angle hues). */
export function idColor(id, s = 0.62, l = 0.58, out = new THREE.Color()) {
  return out.setHSL(((id * 0.61803398875) % 1 + 1) % 1, s, l);
}
