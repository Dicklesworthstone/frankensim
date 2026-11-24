/**
 * Navigation and collision for people on a kitchen floor: where a body may stand, how it gets from one station to
 * another, and when it may go so that it never walks into anyone.
 *
 * Floor frame: three.js world x and z, metres (y is up and plays no part here). A heading `yaw` follows people.ts:
 * 0 faces +z, and a body's forward is (sin yaw, cos yaw), its left (cos yaw, -sin yaw).
 *
 * The pieces:
 *   NavGrid    the room's free floor. Static obstacles are rectangles (fixture footprints, and the room's walls as its
 *              bounds). Every cell stores the exact distance from its centre to the nearest obstacle (its clearance),
 *              so a body of radius r may stand where the clearance is at least r: the obstacles inflated by r.
 *   findPath   A* over the 8-connected cells that are clear for radius r (and of any extra, dynamic obstacles), with a
 *              cost that prefers the middle of an aisle, then string-pulled (each corner dropped while the straight
 *              line past it stays clear), then shifted towards the walker's right where the aisle has room, as people
 *              keep to one side so that two can pass.
 *   Crowd      people as footprints (circles in their body frame, so a body is wider across the shoulders than it is
 *              deep) moving along committed space-time tracks sampled every `dt` seconds. A new leg is planned against
 *              everyone else's committed track: its path is searched first, then the departure time and walking speed
 *              are chosen (a velocity-obstacle search over the candidates, in order of preference: go now at the
 *              preferred speed, slower, or wait and yield) so that its sampled motion stays clear of every other
 *              body, of every door or rack that is held open or parked, and of the room, for its whole length and
 *              afterwards while it stands at its goal. Only then is it committed. A leg that finds no such slot is
 *              not taken: the person keeps working where they are and asks again later.
 *   The motion of a leg is not computed here: the caller's `realize` turns a path and a speed into samples (the
 *   hero walks people.ts Routines and samples their root), so what is checked is exactly what will be drawn.
 *
 * Invariants (crowd.test.ts checks them; `Crowd.audit` measures them on the live tracks):
 *   I1 (walls and fixtures) Every committed sample keeps every footprint circle's centre at least its radius plus
 *      STATIC_MARGIN from every static obstacle and the walls (a `high` circle, from the walls and the tall ones). Between samples the root and heading are
 *      interpolated and checked at points no more than STATIC_MARGIN of circle travel apart, each with room of
 *      at least STATIC_MARGIN / 2; clearance is 1-Lipschitz, so no circle ever enters an obstacle.
 *   I2 (people) For any two people and any time up to the later of their committed ends, the distance between any
 *      circle of one and any circle of the other is at least the sum of their radii plus PAIR_MARGIN. Positions
 *      between samples are taken as the straight line between them and the closest approach on it is checked.
 *   I3 (the future) After the last sample of its track a person stands still at its last pose until a new leg is
 *      committed. Every leg is checked against the others over the union of both tracks' spans with that rule, so
 *      I2 holds for all time by induction over commits: nobody is ever committed into anybody's present or future.
 *   I4 (held and parked things) A timed rectangle (a door swung open, a rack parked) is an obstacle to everyone but
 *      its owner for its interval; it is checked like a person (I2) on both sides, when it is committed and when
 *      anything else is.
 *   I5 (goals) A committed leg ends where the caller's motion ends: `planLeg` never alters the motion it is given,
 *      only when it starts.
 *   I6 (gates) A gate (a walk-in's doorway and the room behind it) is an obstacle to everyone at all times except to
 *      the holder of a window on it, and only during that window; windows on one gate never overlap, so whoever is
 *      inside has it to themselves.
 *   A leg may take over a parked thing (`release`: the rack a person starts to push stops being a parked obstacle
 *   when their hands take it, and moves on as part of their footprint) and may leave things behind it (`timed`: the
 *   door it opens, for as long as it stands open; the rack it parks, from then on, for everyone including itself).
 *   Taking over only ever frees space, and what is left behind is checked against everyone's committed future (I3)
 *   before the leg is committed, so I2-I4 keep holding.
 */

export type XZ = readonly [number, number];

/** An axis-aligned rectangle on the floor. */
export interface Rect {
  minX: number;
  minZ: number;
  maxX: number;
  maxZ: number;
}

/**
 * A circle of a footprint, in the body frame: x to the person's left, z forward, metres. A `high` circle is held up
 * at the waist (a hotel pan carried in both hands): it passes over counters, so of the room only the walls and the
 * fixtures that stand taller than that (`Crowd.tall`) stop it; people stop it like any other circle.
 */
export interface Circle {
  x: number;
  z: number;
  r: number;
  high?: boolean;
}
export type Footprint = readonly Circle[];

/** Clearance kept from walls and fixtures beyond a circle's radius (I1). */
export const STATIC_MARGIN = 0.02;
/** Clearance kept between two people beyond the sum of their circles' radii (I2). */
export const PAIR_MARGIN = 0.03;
/** The fastest any sampled body may move, m/s (a brisk walk). */
export const MAX_SPEED = 1.6;

/** Signed distance from (x, z) to a rectangle: negative inside. */
export function rectDistance(r: Rect, x: number, z: number): number {
  const dx = Math.max(r.minX - x, 0, x - r.maxX);
  const dz = Math.max(r.minZ - z, 0, z - r.maxZ);
  if (dx > 0 || dz > 0) return Math.hypot(dx, dz);
  return -Math.min(x - r.minX, r.maxX - x, z - r.minZ, r.maxZ - z);
}

/** A footprint circle's world centre for a body at (x, z) facing yaw. */
function circleWorld(c: Circle, x: number, z: number, yaw: number, out: { x: number; z: number }): void {
  const s = Math.sin(yaw);
  const co = Math.cos(yaw);
  out.x = x + co * c.x + s * c.z;
  out.z = z - s * c.x + co * c.z;
}

/** Bounding radius of a footprint about the body's root. */
export function footprintReach(f: Footprint): number {
  let r = 0;
  for (const c of f) r = Math.max(r, Math.hypot(c.x, c.z) + c.r);
  return r;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The free floor                                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

class Heap {
  private k: number[] = [];
  private v: number[] = [];
  get size(): number {
    return this.k.length;
  }
  push(key: number, val: number): void {
    const k = this.k;
    const v = this.v;
    let i = k.length;
    k.push(key);
    v.push(val);
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (k[p]! <= key) break;
      k[i] = k[p]!;
      v[i] = v[p]!;
      i = p;
    }
    k[i] = key;
    v[i] = val;
  }
  pop(): number {
    const k = this.k;
    const v = this.v;
    const top = v[0]!;
    const lk = k.pop()!;
    const lv = v.pop()!;
    const n = k.length;
    if (n > 0) {
      let i = 0;
      for (;;) {
        const a = 2 * i + 1;
        if (a >= n) break;
        const b = a + 1;
        const c = b < n && k[b]! < k[a]! ? b : a;
        if (k[c]! >= lk) break;
        k[i] = k[c]!;
        v[i] = v[c]!;
        i = c;
      }
      k[i] = lk;
      v[i] = lv;
    }
    return top;
  }
}

/** Extra obstacles a query must also keep clear of: returns the clearance at (x, z) from them (Infinity if none). */
export type DynamicClearance = (x: number, z: number) => number;

export interface PathOptions {
  /** Other obstacles for this search only (people standing still, a parked rack). */
  dynamic?: DynamicClearance;
  /** How far towards the walker's right of the aisle's middle to walk, 0..1 of the room there (default 0.5). */
  keepRight?: number;
}

export class NavGrid {
  readonly cell: number;
  readonly x0: number;
  readonly z0: number;
  readonly nx: number;
  readonly nz: number;
  /** Exact clearance at each cell centre, metres (negative inside an obstacle), capped at 4 m. */
  readonly clr: Float32Array;

  /** `bounds` is the room's inside face of its walls; `rects` the obstacles on the floor. */
  constructor(
    readonly bounds: Rect,
    readonly rects: readonly Rect[],
    cell = 0.04,
  ) {
    this.cell = cell;
    this.x0 = bounds.minX;
    this.z0 = bounds.minZ;
    this.nx = Math.max(2, Math.ceil((bounds.maxX - bounds.minX) / cell));
    this.nz = Math.max(2, Math.ceil((bounds.maxZ - bounds.minZ) / cell));
    this.clr = new Float32Array(this.nx * this.nz);
    for (let j = 0; j < this.nz; j++)
      for (let i = 0; i < this.nx; i++) this.clr[j * this.nx + i] = Math.min(4, this.clearance(this.x0 + (i + 0.5) * cell, this.z0 + (j + 0.5) * cell));
  }

  /** Exact clearance at a point: the distance to the nearest wall or obstacle, negative inside one. */
  clearance(x: number, z: number): number {
    const b = this.bounds;
    let d = Math.min(x - b.minX, b.maxX - x, z - b.minZ, b.maxZ - z);
    for (const r of this.rects) {
      const q = rectDistance(r, x, z);
      if (q < d) d = q;
    }
    return d;
  }

  /** A fast lower bound on the clearance at a point, from the four cells round it (clearance is 1-Lipschitz). */
  clearanceLB(x: number, z: number): number {
    const fx = (x - this.x0) / this.cell - 0.5;
    const fz = (z - this.z0) / this.cell - 0.5;
    const i0 = Math.floor(fx);
    const j0 = Math.floor(fz);
    let best = -Infinity;
    for (let dj = 0; dj <= 1; dj++)
      for (let di = 0; di <= 1; di++) {
        const i = i0 + di;
        const j = j0 + dj;
        if (i < 0 || j < 0 || i >= this.nx || j >= this.nz) continue;
        const cx = this.x0 + (i + 0.5) * this.cell;
        const cz = this.z0 + (j + 0.5) * this.cell;
        const lb = this.clr[j * this.nx + i]! - Math.hypot(x - cx, z - cz);
        if (lb > best) best = lb;
      }
    return best === -Infinity ? this.clearance(x, z) : best;
  }

  private index(x: number, z: number): [number, number] {
    return [Math.floor((x - this.x0) / this.cell), Math.floor((z - this.z0) / this.cell)];
  }
  private centre(i: number, j: number): [number, number] {
    return [this.x0 + (i + 0.5) * this.cell, this.z0 + (j + 0.5) * this.cell];
  }

  /** Whether the straight segment a-b keeps clearance r (and of the dynamic obstacles), checked every quarter cell. */
  segmentClear(ax: number, az: number, bx: number, bz: number, r: number, dyn?: DynamicClearance): boolean {
    const len = Math.hypot(bx - ax, bz - az);
    const n = Math.max(1, Math.ceil(len / (this.cell * 0.25)));
    for (let k = 0; k <= n; k++) {
      const u = k / n;
      const x = ax + (bx - ax) * u;
      const z = az + (bz - az) * u;
      if (this.clearanceLB(x, z) < r && this.clearance(x, z) < r) return false;
      if (dyn && dyn(x, z) < r) return false;
    }
    return true;
  }

  /** The nearest cell to (x, z) clear for radius r, within `reach` metres, or null. */
  private nearestFree(x: number, z: number, r: number, reach: number, dyn?: DynamicClearance): [number, number] | null {
    const [ci, cj] = this.index(x, z);
    const R = Math.ceil(reach / this.cell);
    let best: [number, number] | null = null;
    let bestD = Infinity;
    for (let dj = -R; dj <= R; dj++)
      for (let di = -R; di <= R; di++) {
        const i = ci + di;
        const j = cj + dj;
        if (i < 0 || j < 0 || i >= this.nx || j >= this.nz) continue;
        if (this.clr[j * this.nx + i]! < r) continue;
        const [px, pz] = this.centre(i, j);
        if (dyn && dyn(px, pz) < r) continue;
        const d = Math.hypot(px - x, pz - z);
        if (d < bestD && d <= reach) {
          bestD = d;
          best = [i, j];
        }
      }
    return best;
  }

  /**
   * A walkable path for a body of radius r from `from` to `to`: A* on the grid, string-pulled, and kept to the
   * right where there is room. Both ends must be clear for r themselves (null if not, or if no way exists). The
   * returned polyline starts at `from` and ends at `to`, and every point on it keeps clearance r.
   */
  findPath(from: XZ, to: XZ, r: number, o: PathOptions = {}): XZ[] | null {
    const path = this.route(from, to, r, o.dynamic);
    return path && this.keepRight(path, r, o.keepRight ?? 0.5, o.dynamic);
  }

  /** `findPath` for several keep-right fractions at the cost of one search: a path (or null) for each. */
  findPaths(from: XZ, to: XZ, r: number, fracs: readonly number[], dyn?: DynamicClearance): (XZ[] | null)[] {
    const path = this.route(from, to, r, dyn);
    return fracs.map((f) => path && this.keepRight(path, r, f, dyn));
  }

  /**
   * As `findPaths`, a slice at a time: the search yields after every `chunk` cells it expands (and before keeping to
   * the right), so that a long search across a big kitchen spreads over frames. The result is the same.
   */
  *findPathsSteps(from: XZ, to: XZ, r: number, fracs: readonly number[], dyn?: DynamicClearance, chunk = 1500): Generator<void, (XZ[] | null)[]> {
    const path = yield* this.routeSteps(from, to, r, dyn, chunk);
    const out: (XZ[] | null)[] = [];
    for (const f of fracs) {
      yield;
      out.push(path && this.keepRight(path, r, f, dyn));
    }
    return out;
  }

  /** The string-pulled shortest comfortable path, before any keeping to the right. */
  private route(from: XZ, to: XZ, r: number, dyn?: DynamicClearance): XZ[] | null {
    const it = this.routeSteps(from, to, r, dyn, 0);
    let x = it.next();
    while (!x.done) x = it.next();
    return x.value;
  }

  /** The search behind `route`, yielding after every `chunk` cells expanded (never, for 0). */
  private *routeSteps(from: XZ, to: XZ, r: number, dyn: DynamicClearance | undefined, chunk: number): Generator<void, XZ[] | null> {
    const okAt = (x: number, z: number) => this.clearance(x, z) >= r && (!dyn || dyn(x, z) >= r);
    if (!okAt(from[0], from[1]) || !okAt(to[0], to[1])) return null;
    if (this.segmentClear(from[0], from[1], to[0], to[1], r, dyn)) return [from, to];
    const s = this.nearestFree(from[0], from[1], r, 0.12, dyn);
    const g = this.nearestFree(to[0], to[1], r, 0.12, dyn);
    if (!s || !g) return null;
    const { nx, nz } = this;
    const N = nx * nz;
    const cost = new Float32Array(N).fill(Infinity);
    const from_ = new Int32Array(N).fill(-1);
    const closed = new Uint8Array(N);
    const blocked = new Uint8Array(N);
    const start = s[1] * nx + s[0];
    const goal = g[1] * nx + g[0];
    const [gx, gz] = this.centre(g[0], g[1]);
    // Prefer a body's width of room beyond the minimum (the aisle's middle); hugging an edge costs up to 2.5x.
    const comfy = 0.25;
    const h = (i: number, j: number) => {
      const dx = Math.abs(i - g[0]);
      const dz = Math.abs(j - g[1]);
      return this.cell * (Math.max(dx, dz) + (Math.SQRT2 - 1) * Math.min(dx, dz));
    };
    const heap = new Heap();
    cost[start] = 0;
    heap.push(h(s[0], s[1]), start);
    const dirs = [
      [1, 0, 1],
      [-1, 0, 1],
      [0, 1, 1],
      [0, -1, 1],
      [1, 1, Math.SQRT2],
      [1, -1, Math.SQRT2],
      [-1, 1, Math.SQRT2],
      [-1, -1, Math.SQRT2],
    ] as const;
    // Each cell's room (static, and from the dynamic obstacles), worked out once.
    const room = new Float32Array(N);
    const free = (c: number): boolean => {
      if (blocked[c] === 1) return false;
      if (blocked[c] === 2) return true;
      let q = this.clr[c]!;
      if (q >= r && dyn) {
        const [px, pz] = this.centre(c % nx, (c / nx) | 0);
        q = Math.min(q, dyn(px, pz));
      }
      room[c] = q;
      const ok = q >= r;
      blocked[c] = ok ? 2 : 1;
      return ok;
    };
    let found = false;
    let pops = 0;
    while (heap.size) {
      const c = heap.pop();
      if (closed[c]) continue;
      closed[c] = 1;
      if (c === goal) {
        found = true;
        break;
      }
      if (++pops > 400000) break;
      if (chunk > 0 && pops % chunk === 0) yield;
      const i = c % nx;
      const j = (c / nx) | 0;
      for (const [di, dj, w] of dirs) {
        const ii = i + di;
        const jj = j + dj;
        if (ii < 0 || jj < 0 || ii >= nx || jj >= nz) continue;
        const n = jj * nx + ii;
        if (closed[n] || !free(n)) continue;
        // No corner cutting past a blocked cell.
        if (di !== 0 && dj !== 0 && (!free(j * nx + ii) || !free(jj * nx + i))) continue;
        const spare = room[n]! - r;
        const k = 1 + 1.5 * Math.max(0, (comfy - spare) / comfy);
        const nc = cost[c]! + w * this.cell * k;
        if (nc < cost[n]!) {
          cost[n] = nc;
          from_[n] = c;
          heap.push(nc + h(ii, jj), n);
        }
      }
    }
    if (!found) return null;
    const cells: XZ[] = [];
    for (let c = goal; c !== -1; c = from_[c]!) cells.push(this.centre(c % nx, (c / nx) | 0));
    cells.reverse();
    // The real end points, joined to the grid path where the straight line to it is clear.
    const pts: XZ[] = [from, ...cells, to];
    if (!this.segmentClear(from[0], from[1], cells[0]![0], cells[0]![1], r, dyn)) return null;
    if (!this.segmentClear(to[0], to[1], cells[cells.length - 1]![0], cells[cells.length - 1]![1], r, dyn)) return null;
    void gx;
    void gz;
    return chunk > 0 ? yield* this.pullSteps(pts, r, dyn, Math.max(8, chunk >> 6)) : this.pull(pts, r, dyn);
  }

  /** String pulling: from each kept point, jump to the farthest later point in clear sight. */
  pull(pts: readonly XZ[], r: number, dyn?: DynamicClearance): XZ[] {
    const it = this.pullSteps(pts, r, dyn, 0);
    let x = it.next();
    while (!x.done) x = it.next();
    return x.value;
  }

  /** `pull`, yielding after every `chunk` sight lines it tests (never, for 0). */
  private *pullSteps(pts: readonly XZ[], r: number, dyn: DynamicClearance | undefined, chunk: number): Generator<void, XZ[]> {
    const out: XZ[] = [pts[0]!];
    let i = 0;
    let tests = 0;
    while (i < pts.length - 1) {
      let j = pts.length - 1;
      while (j > i + 1 && !this.segmentClear(pts[i]![0], pts[i]![1], pts[j]![0], pts[j]![1], r, dyn)) {
        j--;
        if (chunk > 0 && ++tests % chunk === 0) yield;
      }
      out.push(pts[j]!);
      i = j;
    }
    return out;
  }

  /**
   * Walk to the right of the aisle's middle: resample the path, move each point towards its right by `frac` of the
   * room there (eased in and out over the first and last 0.6 m, smoothed along the path), drop the points that are
   * not needed, and keep the result only if every part of it keeps clearance r; else the path as it was.
   */
  keepRight(path: readonly XZ[], r: number, frac: number, dyn?: DynamicClearance): XZ[] {
    if (frac <= 0 || path.length < 2) return path.slice();
    const step = 0.2;
    const pts: { x: number; z: number; nx: number; nz: number; s: number }[] = [];
    let s = 0;
    for (let k = 0; k < path.length - 1; k++) {
      const [ax, az] = path[k]!;
      const [bx, bz] = path[k + 1]!;
      const L = Math.hypot(bx - ax, bz - az);
      if (L < 1e-6) continue;
      const n = Math.max(1, Math.round(L / step));
      // The walker's right: forward (dx, dz) turned clockwise seen from above, (-dz, dx).
      const rx = -(bz - az) / L;
      const rz = (bx - ax) / L;
      for (let q = 0; q < n; q++) pts.push({ x: ax + ((bx - ax) * q) / n, z: az + ((bz - az) * q) / n, nx: rx, nz: rz, s: s + (L * q) / n });
      s += L;
    }
    const last = path[path.length - 1]!;
    pts.push({ x: last[0], z: last[1], nx: pts[pts.length - 1]?.nx ?? 0, nz: pts[pts.length - 1]?.nz ?? 0, s });
    const total = s;
    if (total < 1.4) return path.slice();
    const clearAt = (x: number, z: number) => Math.min(this.clearance(x, z), dyn ? dyn(x, z) : Infinity);
    /**
     * How far a point may move along (nx, nz), up to 1 m, keeping clearance r + 1 cm: marched by sphere tracing
     * (clearance is 1-Lipschitz, so every point within the clearance's excess of a checked one is clear too).
     */
    const room = (p: { x: number; z: number }, nx: number, nz: number) => {
      let u = 0.02;
      let last = 0;
      while (u <= 1) {
        const c = clearAt(p.x + nx * u, p.z + nz * u) - (r + 0.01);
        if (c < 0) break;
        last = Math.min(1, u + c);
        u = last + 0.01;
      }
      return last;
    };
    const right = pts.map((p) => room(p, p.nx, p.nz));
    const left = pts.map((p) => room(p, -p.nx, -p.nz));
    const clamp = (v: number, i: number) => Math.max(-left[i]!, Math.min(right[i]!, v));
    const attempt = (f: number): XZ[] | null => {
      const shift = pts.map((p, i) => {
        const ease = Math.min(1, p.s / 0.6, (total - p.s) / 0.6);
        if (ease <= 0) return 0;
        // The middle of the aisle is (right - left) / 2 off the path; go `f` of the way from it to the right edge.
        const mid = (right[i]! - left[i]!) / 2;
        return clamp(mid + f * (right[i]! - mid), i) * ease;
      });
      // Smooth the shifts along the path (a body does not zig-zag), never beyond what each point allows.
      const sm = shift.map((_, i) => {
        let a = 0;
        let w = 0;
        for (let k = -3; k <= 3; k++) {
          const v = shift[i + k];
          if (v === undefined) continue;
          const wk = 4 - Math.abs(k);
          a += v * wk;
          w += wk;
        }
        return clamp(a / w, i);
      });
      const moved: XZ[] = pts.map((p, i) => [p.x + p.nx * sm[i]!, p.z + p.nz * sm[i]!] as const);
      moved[0] = path[0]!;
      moved[moved.length - 1] = last;
      const clear = (q: readonly XZ[]) => {
        for (let k = 0; k < q.length - 1; k++) if (!this.segmentClear(q[k]![0], q[k]![1], q[k + 1]![0], q[k + 1]![1], r, dyn)) return false;
        return true;
      };
      const simple = this.simplify(moved, 0.03);
      if (clear(simple)) return simple;
      return clear(moved) ? moved : null;
    };
    for (const f of [frac, frac * 0.6, frac * 0.3]) {
      const q = attempt(f);
      if (q) return q;
    }
    return path.slice();
  }

  /** Ramer-Douglas-Peucker simplification within `tol` metres. */
  simplify(pts: readonly XZ[], tol: number): XZ[] {
    if (pts.length < 3) return pts.slice();
    const keep = new Uint8Array(pts.length);
    keep[0] = keep[pts.length - 1] = 1;
    const stack: [number, number][] = [[0, pts.length - 1]];
    while (stack.length) {
      const [a, b] = stack.pop()!;
      const [ax, az] = pts[a]!;
      const [bx, bz] = pts[b]!;
      const L = Math.hypot(bx - ax, bz - az) || 1e-9;
      let best = -1;
      let bi = -1;
      for (let i = a + 1; i < b; i++) {
        const [px, pz] = pts[i]!;
        const d = Math.abs((bx - ax) * (az - pz) - (ax - px) * (bz - az)) / L;
        if (d > best) {
          best = d;
          bi = i;
        }
      }
      if (best > tol && bi > 0) {
        keep[bi] = 1;
        stack.push([a, bi], [bi, b]);
      }
    }
    return pts.filter((_, i) => keep[i]);
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* People in space and time                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A leg's motion, sampled every `Crowd.dt` from its start: the body's root, heading and which of its footprints it
 * shows. After the last sample the body stays as it is there. `timed` lists rectangles the leg holds for a while
 * (a door it opens), in samples from the leg's start.
 */
export interface Motion {
  n: number;
  x: Float32Array;
  z: Float32Array;
  yaw: Float32Array;
  shape: Uint8Array;
  timed?: TimedSpec[];
  /** Timed rectangles (by key) this leg takes over at step k from its start: they stop being obstacles there. */
  release?: { key: string; k: number }[];
  /** Windows (steps from the leg's start) in which this leg passes through a gate (by key). */
  gates?: { key: string; k0: number; k1: number }[];
}

/**
 * A rectangle a leg holds over steps k0..k1 from its start (k1 may be Infinity: parked for good). Its owner is the
 * person moving (an obstacle to everyone else), unless `owner` says otherwise (-1: an obstacle to everyone, the
 * mover too, as a rack it leaves parked). A `key` lets a later leg take it over (`Motion.release`).
 */
export interface TimedSpec {
  rect: Rect;
  k0: number;
  k1: number;
  owner?: number;
  key?: string;
}

/** One person's committed track: samples from step `k0` on; before it unknown (past), after its end, standing. */
class Track {
  k0: number;
  x: number[] = [];
  z: number[] = [];
  yaw: number[] = [];
  shape: number[] = [];
  constructor(k0: number, x: number, z: number, yaw: number, shape: number) {
    this.k0 = k0;
    this.x.push(x);
    this.z.push(z);
    this.yaw.push(yaw);
    this.shape.push(shape);
  }
  /** The last committed step. */
  get end(): number {
    return this.k0 + this.x.length - 1;
  }
  index(k: number): number {
    return Math.max(0, Math.min(this.x.length - 1, k - this.k0));
  }
  /** Drop the samples before step k (keeping at least one). */
  prune(k: number): void {
    const cut = Math.min(this.x.length - 1, Math.max(0, k - this.k0));
    if (cut <= 0) return;
    this.x.splice(0, cut);
    this.z.splice(0, cut);
    this.yaw.splice(0, cut);
    this.shape.splice(0, cut);
    this.k0 += cut;
  }
}

interface TimedRect {
  rect: Rect;
  k0: number;
  k1: number;
  owner: number;
  key?: string;
}

interface Gate {
  rect: Rect;
  windows: { k0: number; k1: number; owner: number }[];
}

/** Whether two rectangles come within `m` of each other. */
function rectsNear(a: Rect, b: Rect, m: number): boolean {
  return a.minX - m < b.maxX && b.minX - m < a.maxX && a.minZ - m < b.maxZ && b.minZ - m < a.maxZ;
}

export interface AgentInfo {
  /** Footprints by shape id (the ids a Motion's `shape` uses). */
  shapes: readonly Footprint[];
}

/** Where a candidate leg failed, for the caller's log and for tests. */
export type LegVerdict = { ok: true; delay: number; speed: number; variant: number } | { ok: false; reason: string };

export interface LegCandidate {
  /** A path variant's motion at a speed, or null when it cannot be walked (the caller's realize failed). */
  motion?: Motion | null;
  /** Or a way to make it, called only if the planner gets to this candidate (realizing a motion costs time). */
  realize?: () => Motion | null;
  speed: number;
  variant: number;
}

const _a = { x: 0, z: 0 };
const _b = { x: 0, z: 0 };
const _c = { x: 0, z: 0 };
const _d = { x: 0, z: 0 };

export class Crowd {
  readonly agents: AgentInfo[] = [];
  readonly tracks: Track[] = [];
  private timed: TimedRect[] = [];
  private gates = new Map<string, Gate>();

  /**
   * `dt` is the sampling step in seconds. `tall` is the room with only its walls and the fixtures that stand above
   * waist height (reach-ins, racks, shelving, ovens), for the `high` circles; by default the whole room.
   */
  constructor(
    readonly grid: NavGrid,
    readonly dt = 0.05,
    readonly tall: NavGrid = grid,
  ) {}

  /** Add a person standing at (x, z) facing yaw, showing footprint `shape`, from step k. Returns its id. */
  addAgent(info: AgentInfo, x: number, z: number, yaw: number, shape: number, k = 0): number {
    this.agents.push(info);
    this.tracks.push(new Track(k, x, z, yaw, shape));
    return this.agents.length - 1;
  }

  /** The step a time falls in. */
  step(t: number): number {
    return Math.round(t / this.dt);
  }
  /** A person's last committed step. */
  end(id: number): number {
    return this.tracks[id]!.end;
  }
  /** A person's committed pose at step k (standing still past the end). */
  poseAt(id: number, k: number): { x: number; z: number; yaw: number; shape: number } {
    const t = this.tracks[id]!;
    const i = t.index(k);
    return { x: t.x[i]!, z: t.z[i]!, yaw: t.yaw[i]!, shape: t.shape[i]! };
  }

  /** Forget what lies before step k (nothing is ever planned into the past). */
  prune(k: number): void {
    for (const t of this.tracks) t.prune(k - 40);
    this.timed = this.timed.filter((r) => r.k1 >= k - 40);
    for (const g of this.gates.values()) g.windows = g.windows.filter((w) => w.k1 >= k - 40);
  }

  /**
   * A rectangle that is an obstacle to everyone but `owner` over steps k0..k1 (a parked rack, an open door; owner
   * -1 for everyone; k1 Infinity for good). A `key` lets a leg take it over later (`Motion.release`).
   */
  addTimedRect(rect: Rect, k0: number, k1: number, owner: number, key?: string): void {
    this.timed.push({ rect, k0, k1, owner, key });
  }

  /** The timed rectangle with this key that is held latest (a rack's current parking place), if any. */
  timedRect(key: string): { rect: Rect; k0: number; k1: number; owner: number } | null {
    let best: TimedRect | null = null;
    for (const r of this.timed) if (r.key === key && (!best || r.k1 > best.k1)) best = r;
    return best;
  }

  /** A gate (I6): a rectangle closed to everyone except the holder of a window on it, during the window. */
  addGate(key: string, rect: Rect): void {
    this.gates.set(key, { rect, windows: [] });
  }
  /** A gate's windows from step k on, for display and checks. */
  gateWindows(key: string, k = -Infinity): readonly { k0: number; k1: number; owner: number }[] {
    return (this.gates.get(key)?.windows ?? []).filter((w) => w.k1 >= k);
  }

  /** Footprint clearance from the static floor for a person `id` at a pose: min over circles of clearance - r. */
  staticRoom(id: number, x: number, z: number, yaw: number, shape: number): number {
    let worst = Infinity;
    for (const c of this.agents[id]!.shapes[shape] ?? []) {
      circleWorld(c, x, z, yaw, _a);
      const g = c.high ? this.tall : this.grid;
      const q = g.clearanceLB(_a.x, _a.z);
      const room = (q >= c.r + STATIC_MARGIN ? q : g.clearance(_a.x, _a.z)) - c.r;
      if (room < worst) worst = room;
    }
    return worst;
  }

  /**
   * Dynamic obstacles for a path search by person `id` from step k: everyone else as they will stand at the end of
   * their committed tracks, and the timed rectangles held for longer than `longHold` seconds after k (a rack parked
   * for the shift; a door that shuts again in a moment is left to the timing, `planLeg`). Returns the clearance.
   */
  standingObstacles(id: number, k: number, longHold = 10, ignore: readonly string[] = []): DynamicClearance {
    const circles: { x: number; z: number; r: number }[] = [];
    for (let j = 0; j < this.tracks.length; j++) {
      if (j === id) continue;
      const p = this.poseAt(j, this.tracks[j]!.end);
      for (const c of this.agents[j]!.shapes[p.shape] ?? []) {
        circleWorld(c, p.x, p.z, p.yaw, _a);
        circles.push({ x: _a.x, z: _a.z, r: c.r + PAIR_MARGIN });
      }
    }
    // Gates are closed to a search unless it means to pass one (and then the timing is left to `planLeg`).
    const rects = this.timed.filter((r) => r.owner !== id && r.k1 >= k + longHold / this.dt && !(r.key && ignore.includes(r.key))).map((r) => r.rect);
    for (const [key, g] of this.gates) if (!ignore.includes(key)) rects.push(g.rect);
    return (x, z) => {
      let d = Infinity;
      for (const c of circles) {
        const q = Math.hypot(x - c.x, z - c.z) - c.r;
        if (q < d) d = q;
      }
      for (const r of rects) {
        const q = rectDistance(r, x, z);
        if (q < d) d = q;
      }
      return d;
    };
  }

  /**
   * I1 for a motion: every sample, and enough points between samples that no circle is ever more than
   * STATIC_MARGIN from a checked position of itself, keeps every circle clear of the static floor by STATIC_MARGIN.
   */
  motionClearOfRoom(id: number, m: Motion): boolean {
    const shapes = this.agents[id]!.shapes;
    for (let i = 0; i < m.n; i++) {
      if (this.staticRoom(id, m.x[i]!, m.z[i]!, m.yaw[i]!, m.shape[i]!) < STATIC_MARGIN) return false;
      if (i + 1 >= m.n) continue;
      const dyaw = wrapPi(m.yaw[i + 1]! - m.yaw[i]!);
      const reach = footprintReach(shapes[m.shape[i]!] ?? []);
      // The farthest any circle moves over the step: the root's travel plus the turn swept at the circle's reach.
      const sweep = Math.hypot(m.x[i + 1]! - m.x[i]!, m.z[i + 1]! - m.z[i]!) + Math.abs(dyaw) * reach;
      const sub = Math.max(2, Math.ceil(sweep / STATIC_MARGIN));
      for (let q = 1; q < sub; q++) {
        const u = q / sub;
        if (this.staticRoom(id, m.x[i]! + (m.x[i + 1]! - m.x[i]!) * u, m.z[i]! + (m.z[i + 1]! - m.z[i]!) * u, m.yaw[i]! + dyaw * u, m.shape[i]!) < STATIC_MARGIN / 2) return false;
      }
    }
    return true;
  }

  /**
   * Whether person `id`, standing at its current end pose until step k0 + delay and then following motion `m`
   * (from step k0 + delay, standing at its last sample afterwards), keeps clear of everyone else's committed track
   * and of the timed rectangles, and whether its own timed rectangles keep clear of everyone else (I2-I4).
   */
  legClear(id: number, k0: number, delay: number, m: Motion): boolean {
    const me = this.tracks[id]!;
    const startPose = this.poseAt(id, me.end);
    const kStart = k0 + delay;
    const kMine = kStart + m.n - 1;
    const myShapes = this.agents[id]!.shapes;
    const myReach = Math.max(...myShapes.map(footprintReach));
    const mine = (k: number, out: { x: number; z: number; yaw: number; shape: number }) => {
      if (k < kStart) {
        out.x = startPose.x;
        out.z = startPose.z;
        out.yaw = startPose.yaw;
        out.shape = startPose.shape;
        return;
      }
      const i = Math.min(m.n - 1, k - kStart);
      out.x = m.x[i]!;
      out.z = m.z[i]!;
      out.yaw = m.yaw[i]!;
      out.shape = m.shape[i]!;
    };
    const A0 = { x: 0, z: 0, yaw: 0, shape: 0 };
    const A1 = { x: 0, z: 0, yaw: 0, shape: 0 };
    const B0 = { x: 0, z: 0, yaw: 0, shape: 0 };
    const B1 = { x: 0, z: 0, yaw: 0, shape: 0 };
    for (let j = 0; j < this.tracks.length; j++) {
      if (j === id) continue;
      const other = this.tracks[j]!;
      const theirShapes = this.agents[j]!.shapes;
      const reach = myReach + Math.max(...theirShapes.map(footprintReach)) + PAIR_MARGIN + 2 * MAX_SPEED * this.dt;
      const kEnd = Math.max(kMine, other.end) + 1;
      for (let k = k0; k < kEnd; k++) {
        mine(k, A0);
        mine(k + 1, A1);
        const ia = other.index(k);
        const ib = other.index(k + 1);
        B0.x = other.x[ia]!;
        B0.z = other.z[ia]!;
        B0.yaw = other.yaw[ia]!;
        B0.shape = other.shape[ia]!;
        B1.x = other.x[ib]!;
        B1.z = other.z[ib]!;
        B1.yaw = other.yaw[ib]!;
        B1.shape = other.shape[ib]!;
        if (Math.hypot(A0.x - B0.x, A0.z - B0.z) > reach) continue;
        if (!pairClear(myShapes, A0, A1, theirShapes, B0, B1)) return false;
      }
    }
    // What this leg takes over stops being an obstacle the step its hands take it.
    const released = new Map<string, number>();
    for (const r of m.release ?? []) released.set(r.key, kStart + r.k);
    const myRectClear = (rect: Rect, ka: number, kb: number) => {
      // Past my last sample I stand still: one more step covers the rest of any hold, however long.
      const kc = Math.min(kb, Math.max(ka, kMine + 1));
      for (let k = ka; k <= kc; k++) {
        mine(k, A0);
        for (const c of myShapes[A0.shape] ?? []) {
          circleWorld(c, A0.x, A0.z, A0.yaw, _a);
          if (rectDistance(rect, _a.x, _a.z) < c.r + PAIR_MARGIN) return false;
        }
      }
      return true;
    };
    // Timed rectangles held by others (or by nobody), against my circles.
    for (const tr of this.timed) {
      if (tr.owner === id) continue;
      const rel = tr.key !== undefined ? released.get(tr.key) : undefined;
      const k1 = rel !== undefined ? Math.min(tr.k1, rel - 1) : tr.k1;
      if (k1 < k0) continue;
      if (!myRectClear(tr.rect, Math.max(k0, tr.k0), k1)) return false;
    }
    // My timed rectangles, against everyone else (they stand still past their tracks' ends), against the other
    // rectangles held at the same time, and against me when I leave them to everyone.
    for (const tr of m.timed ?? []) {
      const ka = kStart + tr.k0;
      const kb = kStart + tr.k1;
      for (let j = 0; j < this.tracks.length; j++) {
        if (j === id) continue;
        const other = this.tracks[j]!;
        const kc = Math.min(kb, Math.max(ka, other.end + 1));
        for (let k = ka; k <= kc; k++) {
          const i = other.index(k);
          for (const c of this.agents[j]!.shapes[other.shape[i]!] ?? []) {
            circleWorld(c, other.x[i]!, other.z[i]!, other.yaw[i]!, _a);
            if (rectDistance(tr.rect, _a.x, _a.z) < c.r + PAIR_MARGIN) return false;
          }
        }
      }
      for (const o of this.timed) {
        const rel = o.key !== undefined ? released.get(o.key) : undefined;
        const k1 = rel !== undefined ? Math.min(o.k1, rel - 1) : o.k1;
        if (o.k0 <= kb && ka <= k1 && rectsNear(o.rect, tr.rect, PAIR_MARGIN)) return false;
      }
      if (tr.owner !== undefined && tr.owner !== id && !myRectClear(tr.rect, ka, kb)) return false;
    }
    // Gates (I6): my windows clash with nobody's, and outside them I keep out.
    for (const [key, g] of this.gates) {
      const wins = (m.gates ?? []).filter((w) => w.key === key).map((w) => [kStart + w.k0, kStart + w.k1] as const);
      for (const [a, b] of wins) for (const w of g.windows) if (w.owner !== id && w.k0 <= b && a <= w.k1) return false;
      for (let k = k0; k <= kMine + 1; k++) {
        if (wins.some(([a, b]) => a <= k && k <= b)) continue;
        if (g.windows.some((w) => w.owner === id && w.k0 <= k && k <= w.k1)) continue;
        mine(k, A0);
        for (const c of myShapes[A0.shape] ?? []) {
          circleWorld(c, A0.x, A0.z, A0.yaw, _a);
          if (rectDistance(g.rect, _a.x, _a.z) < c.r + PAIR_MARGIN) return false;
        }
      }
    }
    return true;
  }

  /**
   * Keep person `id` standing where its track ends up to step k. Always safe: by I3 everyone else's legs were
   * checked against this person standing there for ever.
   */
  wait(id: number, k: number): void {
    const t = this.tracks[id]!;
    const i = t.x.length - 1;
    for (let q = t.end + 1; q <= k; q++) {
      t.x.push(t.x[i]!);
      t.z.push(t.z[i]!);
      t.yaw.push(t.yaw[i]!);
      t.shape.push(t.shape[i]!);
    }
  }

  /** Commit a checked leg: person `id` stands until k0 + delay, then follows `m`. */
  commit(id: number, k0: number, delay: number, m: Motion): void {
    const t = this.tracks[id]!;
    const p = this.poseAt(id, t.end);
    // Standing from the track's end up to the departure.
    for (let k = t.end + 1; k < k0 + delay; k++) {
      t.x.push(p.x);
      t.z.push(p.z);
      t.yaw.push(p.yaw);
      t.shape.push(p.shape);
    }
    for (let i = 0; i < m.n; i++) {
      t.x.push(m.x[i]!);
      t.z.push(m.z[i]!);
      t.yaw.push(m.yaw[i]!);
      t.shape.push(m.shape[i]!);
    }
    const kStart = k0 + delay;
    for (const r of m.release ?? [])
      for (const tr of this.timed) if (tr.key === r.key) tr.k1 = Math.min(tr.k1, kStart + r.k - 1);
    for (const tr of m.timed ?? []) this.timed.push({ rect: tr.rect, k0: kStart + tr.k0, k1: kStart + tr.k1, owner: tr.owner ?? id, key: tr.key });
    for (const w of m.gates ?? []) this.gates.get(w.key)?.windows.push({ k0: kStart + w.k0, k1: kStart + w.k1, owner: id });
  }

  /**
   * Choose when a leg starts: the first candidate (in the caller's order of preference: path variant, then speed)
   * and the shortest wait up to `maxDelay` seconds that keeps the leg clear of the room and of everyone, and commit
   * it. Person `id` must have been committed up to k0 - 1 (its track's end is k0 - 1 or later is not allowed).
   */
  planLeg(id: number, k0: number, candidates: readonly LegCandidate[], maxDelay = 8): LegVerdict {
    const it = this.planLegSteps(id, k0, candidates, maxDelay);
    for (;;) {
      const r = it.next();
      if (r.done) return r.value;
    }
  }

  /**
   * `planLeg` in slices: it yields after each candidate it makes and each slot it checks, and returns the verdict.
   * Nothing is committed until the end, and the caller must not commit anything else for `id` meanwhile.
   */
  *planLegSteps(id: number, k0: number, candidates: readonly LegCandidate[], maxDelay = 8): Generator<void, LegVerdict> {
    const t = this.tracks[id]!;
    if (t.end >= k0) return { ok: false, reason: `track already committed to ${t.end} >= ${k0}` };
    const maxD = Math.round(maxDelay / this.dt);
    // Each candidate is made and checked against the room at most once, when first needed.
    const made = new Map<LegCandidate, Motion | null>();
    const motionOf = (c: LegCandidate): Motion | null => {
      if (made.has(c)) return made.get(c)!;
      const m = c.motion !== undefined ? c.motion : (c.realize?.() ?? null);
      const ok = m && this.motionClearOfRoom(id, m) ? m : null;
      made.set(c, ok);
      return ok;
    };
    // Waits in growing steps: a quarter second at first, a second later on.
    const delays: number[] = [];
    for (let d = 0; d <= maxD; d += d < 2 / this.dt ? Math.max(1, Math.round(0.25 / this.dt)) : Math.round(1 / this.dt)) delays.push(d);
    let any = false;
    for (const d of delays)
      for (const c of candidates) {
        const fresh = !made.has(c);
        const m = motionOf(c);
        if (fresh) yield;
        if (!m) continue;
        any = true;
        const clear = this.legClear(id, k0, d, m);
        yield;
        if (clear) {
          this.commit(id, k0, d, m);
          return { ok: true, delay: d * this.dt, speed: c.speed, variant: c.variant };
        }
      }
    return { ok: false, reason: any ? "every slot meets someone" : "no candidate keeps clear of the room" };
  }

  /**
   * Measure the invariants over steps k0..k1 on the committed tracks: the least gap between any two people's
   * circles (beyond their radii), and the least static room of any circle (beyond its radius).
   */
  audit(k0: number, k1: number): { minPairGap: number; minStaticRoom: number; pair: [number, number] | null; minRectGap: number; gateIntrusions: number } {
    let minPairGap = Infinity;
    let minStaticRoom = Infinity;
    let minRectGap = Infinity;
    let gateIntrusions = 0;
    let pair: [number, number] | null = null;
    for (let k = k0; k <= k1; k++) {
      for (let i = 0; i < this.tracks.length; i++) {
        const a = this.poseAt(i, k);
        minStaticRoom = Math.min(minStaticRoom, this.staticRoom(i, a.x, a.z, a.yaw, a.shape));
        // Timed rectangles not this person's own (I4), and gates outside this person's windows (I6).
        for (const c of this.agents[i]!.shapes[a.shape] ?? []) {
          circleWorld(c, a.x, a.z, a.yaw, _a);
          for (const tr of this.timed) {
            if (tr.owner === i || k < tr.k0 || k > tr.k1) continue;
            minRectGap = Math.min(minRectGap, rectDistance(tr.rect, _a.x, _a.z) - c.r);
          }
          for (const g of this.gates.values()) {
            if (g.windows.some((w) => w.owner === i && w.k0 <= k && k <= w.k1)) continue;
            if (rectDistance(g.rect, _a.x, _a.z) < c.r) gateIntrusions++;
          }
        }
        for (let j = i + 1; j < this.tracks.length; j++) {
          const b = this.poseAt(j, k);
          for (const ca of this.agents[i]!.shapes[a.shape] ?? []) {
            circleWorld(ca, a.x, a.z, a.yaw, _a);
            for (const cb of this.agents[j]!.shapes[b.shape] ?? []) {
              circleWorld(cb, b.x, b.z, b.yaw, _b);
              const g = Math.hypot(_a.x - _b.x, _a.z - _b.z) - ca.r - cb.r;
              if (g < minPairGap) {
                minPairGap = g;
                pair = [i, j];
              }
            }
          }
        }
      }
    }
    return { minPairGap, minStaticRoom, pair, minRectGap, gateIntrusions };
  }
}

export function wrapPi(a: number): number {
  let x = (a + Math.PI) % (Math.PI * 2);
  if (x < 0) x += Math.PI * 2;
  return x - Math.PI;
}

/**
 * Whether two bodies keep PAIR_MARGIN between every pair of their circles over one step, each moving in a straight
 * line from its pose at the step's start to its pose at its end: the closest approach of each pair on that interval.
 * A body that changes footprint over the step (stops walking, picks up a pan) is taken to show both footprints the
 * whole step, each circle moving with the body.
 */
function pairClear(
  sa: readonly Footprint[],
  a0: { x: number; z: number; yaw: number; shape: number },
  a1: { x: number; z: number; yaw: number; shape: number },
  sb: readonly Footprint[],
  b0: { x: number; z: number; yaw: number; shape: number },
  b1: { x: number; z: number; yaw: number; shape: number },
): boolean {
  const fa = sa[a0.shape] ?? [];
  const fb = sb[b0.shape] ?? [];
  const fa1 = a1.shape !== a0.shape ? (sa[a1.shape] ?? []) : null;
  const fb1 = b1.shape !== b0.shape ? (sb[b1.shape] ?? []) : null;
  const na = fa.length + (fa1?.length ?? 0);
  const nb = fb.length + (fb1?.length ?? 0);
  for (let i = 0; i < na; i++) {
    const ca = i < fa.length ? fa[i]! : fa1![i - fa.length]!;
    circleWorld(ca, a0.x, a0.z, a0.yaw, _a);
    circleWorld(ca, a1.x, a1.z, a1.yaw, _c);
    for (let j = 0; j < nb; j++) {
      const cb = j < fb.length ? fb[j]! : fb1![j - fb.length]!;
      circleWorld(cb, b0.x, b0.z, b0.yaw, _b);
      circleWorld(cb, b1.x, b1.z, b1.yaw, _d);
      // Relative motion p(u) = p0 + (p1 - p0) u, u in [0, 1].
      const px = _a.x - _b.x;
      const pz = _a.z - _b.z;
      const vx = _c.x - _d.x - px;
      const vz = _c.z - _d.z - pz;
      const vv = vx * vx + vz * vz;
      const u = vv > 1e-12 ? Math.max(0, Math.min(1, -(px * vx + pz * vz) / vv)) : 0;
      const dx = px + vx * u;
      const dz = pz + vz * u;
      const need = ca.r + cb.r + PAIR_MARGIN;
      if (dx * dx + dz * dz < need * need) return false;
    }
  }
  return true;
}
