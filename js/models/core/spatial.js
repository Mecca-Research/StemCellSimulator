// Uniform-grid spatial hash for neighbour queries among thousands of cells.
export class SpatialHash {
  constructor(cellSize) {
    this.h = cellSize;
    this.inv = 1 / cellSize;
    this.map = new Map();
  }

  key(ix, iy, iz) {
    // 12 bits per axis (+-2048 buckets); the product stays far below 2^53 so
    // every bucket has an exact, unique numeric key
    return ((ix + 2048) * 4096 + (iy + 2048)) * 4096 + (iz + 2048);
  }

  clear() { this.map.clear(); }

  insert(id, x, y, z) {
    const k = this.key(Math.floor(x * this.inv), Math.floor(y * this.inv), Math.floor(z * this.inv));
    let b = this.map.get(k);
    if (!b) { b = []; this.map.set(k, b); }
    b.push(id);
  }

  /** Call fn(id) for every id in the 27 buckets around (x, y, z). */
  near(x, y, z, fn) {
    const ix = Math.floor(x * this.inv), iy = Math.floor(y * this.inv), iz = Math.floor(z * this.inv);
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
      const b = this.map.get(this.key(ix + dx, iy + dy, iz + dz));
      if (b) for (let i = 0; i < b.length; i++) fn(b[i]);
    }
  }
}
