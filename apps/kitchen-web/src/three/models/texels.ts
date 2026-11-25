/**
 * Small maths, seeded value noise and cellular noise for the procedural models' textures. Free of three.js and the
 * DOM, so a worker can paint with them too (steakTissue.worker.ts).
 */

export type P2 = [number, number];
export const PI = Math.PI;
export const TAU = Math.PI * 2;
export const clamp = (v: number, a: number, b: number) => Math.max(a, Math.min(b, v));
export const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
export const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
export const smooth = (e0: number, e1: number, x: number) => {
  const t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3 - 2 * t);
};

export function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** `n` values of rng(seed) in turn: the same sequence, without a call per value (noise lattices are 64 K long). */
export function lattice(seed: number, n: number): Float32Array {
  const v = new Float32Array(n);
  let s = seed >>> 0;
  for (let i = 0; i < n; i++) {
    s = (s * 1664525 + 1013904223) >>> 0;
    v[i] = s / 4294967296;
  }
  return v;
}

/** Smooth value noise on a 256-periodic lattice, 0..1. */
export function noise2(seed: number): (x: number, y: number) => number {
  const N = 256;
  const v = lattice(seed, N * N);
  return (x, y) => {
    const xi = Math.floor(x);
    const yi = Math.floor(y);
    const fx = x - xi;
    const fy = y - yi;
    const ux = fx * fx * (3 - 2 * fx);
    const uy = fy * fy * (3 - 2 * fy);
    const x0 = xi & 255;
    const y0 = yi & 255;
    const x1 = (x0 + 1) & 255;
    const y1 = (y0 + 1) & 255;
    const a = v[y0 * N + x0];
    const b = v[y0 * N + x1];
    const c = v[y1 * N + x0];
    const d = v[y1 * N + x1];
    return a + (b - a) * ux + (c - a) * uy + (a - b - c + d) * ux * uy;
  };
}
/** Fractal sum of `oct` octaves of `n`, 0..1. */
export function fbm(n: (x: number, y: number) => number, x: number, y: number, oct = 4): number {
  let s = 0;
  let a = 0.5;
  let f = 1;
  let norm = 0;
  for (let i = 0; i < oct; i++) {
    s += a * n(x * f + i * 17.3, y * f - i * 11.1);
    norm += a;
    a *= 0.5;
    f *= 2.03;
  }
  return s / norm;
}

export const adiff = (a: number, b: number) => {
  let d = (a - b) % TAU;
  if (d > PI) d -= TAU;
  if (d < -PI) d += TAU;
  return d;
};

/** Integer hash of a lattice cell to 0..1. */
export function hash2(i: number, j: number, s: number): number {
  let h = (Math.imul(i, 374761393) + Math.imul(j, 668265263) + Math.imul(s, 1103515245)) | 0;
  h = Math.imul(h ^ (h >>> 13), 1274126177);
  h ^= h >>> 16;
  return (h >>> 0) / 4294967296;
}

/**
 * Cellular (Worley) noise over a box: seeds jittered one to a cell of size `c`. `at` writes the nearest seed's
 * distance, the second's, the nearest's id and the second's (0..1). Muscle bundles, fat cells, crumb, blisters.
 */
export class Cells {
  private readonly px: Float32Array;
  private readonly pz: Float32Array;
  private readonly id: Float32Array;
  private readonly nx: number;
  private readonly nz: number;
  constructor(private readonly x0: number, private readonly z0: number, x1: number, z1: number, private readonly c: number, seed: number) {
    this.nx = Math.ceil((x1 - x0) / c) + 3;
    this.nz = Math.ceil((z1 - z0) / c) + 3;
    const n = this.nx * this.nz;
    this.px = new Float32Array(n);
    this.pz = new Float32Array(n);
    this.id = new Float32Array(n);
    for (let j = 0; j < this.nz; j++)
      for (let i = 0; i < this.nx; i++) {
        const k = j * this.nx + i;
        this.px[k] = i + hash2(i, j, seed);
        this.pz[k] = j + hash2(i, j, seed + 7);
        this.id[k] = hash2(i, j, seed + 13);
      }
  }
  at(x: number, z: number, out: Float64Array): void {
    const gx = (x - this.x0) / this.c + 1;
    const gz = (z - this.z0) / this.c + 1;
    const ix = Math.floor(gx);
    const iz = Math.floor(gz);
    const { px, pz, nx, nz } = this;
    const ids = this.id;
    let f1 = 1e9;
    let f2 = 1e9;
    let id = 0;
    let id2 = 0;
    // The 3 × 3 cells round (ix, iz) that lie in the grid, row by row (the order decides ties).
    const j0 = iz > 0 ? iz - 1 : 0;
    const j1 = iz + 1 < nz ? iz + 1 : nz - 1;
    const i0 = ix > 0 ? ix - 1 : 0;
    const i1 = ix + 1 < nx ? ix + 1 : nx - 1;
    for (let jj = j0; jj <= j1; jj++) {
      const row = jj * nx;
      for (let ii = i0; ii <= i1; ii++) {
        const k = row + ii;
        const dx = px[k] - gx;
        const dz = pz[k] - gz;
        const d = dx * dx + dz * dz;
        if (d < f1) {
          f2 = f1;
          id2 = id;
          f1 = d;
          id = ids[k];
        } else if (d < f2) {
          f2 = d;
          id2 = ids[k];
        }
      }
    }
    out[0] = Math.sqrt(f1) * this.c;
    out[1] = Math.sqrt(f2) * this.c;
    out[2] = id;
    out[3] = id2;
  }
}
