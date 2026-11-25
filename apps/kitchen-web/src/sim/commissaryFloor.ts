/**
 * Where the blast chiller should stand: the commissary's own floor plan, its cook-chill trips, and every spot on its
 * walls that a chiller could take.
 *
 * The plan is the commissary exemplar's build manifest (every fixture's floor box, pose and stance socket) and its
 * site metadata (the corridors that must stay clear). Floor frame: three.js x and z, metres; the manifest's Blender y
 * becomes -z.
 *
 * The trips are the week's cooked-ahead batches (commissaryWeek), every one of which goes through the blast chiller:
 * each pan is carried from its kettle, skillet, combi or prep table to the chiller one at a time, and the cook walks
 * back; a chilled load goes to the walk-in on a rack, and the rack comes back. Every walk is measured on the floor
 * round the fixtures for a cook with a pan (NavGrid paths in crowd.ts): distance fields over the free floor rank every
 * spot at once, and the routes drawn and quoted are the string-pulled A* paths with the chiller in place.
 *
 * A spot is a place on a wall where the chiller's footprint fits flush, facing the room, clear of every fixture, every
 * corridor in the site metadata and the floor under a hood, with room for a cook in front of its door.
 */
import { NavGrid, type Rect, type XZ } from "../three/models/crowd";
import type { Batch, Dish } from "./commissaryWeek";

export interface ManifestFixture {
  fixture_id: string;
  archetype_id: string;
  local_extents_m: [number, number, number];
  pose: { x_m: number; y_m: number; z_m: number; yaw_rad: number };
  world_aabb_m?: { min_m: [number, number, number]; max_m: [number, number, number] };
  sockets?: { socket_id: string; local_pos_m: [number, number, number] }[];
}
export interface FloorManifest {
  room_shell: { width_m: number; depth_m: number; ceiling_height_m: number };
  fixtures: ManifestFixture[];
}
export interface FloorMeta {
  walk_mesh_corridors?: { corridor_id: string; min_m: [number, number, number]; max_m: [number, number, number] }[];
}

/** A cook carrying a full hotel pan, half its width plus a little (the hero's walkers use the same idea). */
export const CARRY_R = 0.32;
/** Walking pace with a pan or a loaded rack, m/s (named, not fitted). */
export const CARRY_MPS = 0.9;
/** Pans in one blast-chiller load, wheeled to the walk-in on one rack. */
export const CHILL_LOAD = 15;

const SMALL = /KNIFE|BOARD|HOTEL-PAN|GRATER|SCALE|SLICER|PROCESSOR|BOOSTER/;

/** A fixture's frame on the floor: local (x across its front, y back from its front) to floor x, z. */
export function toFloor(p: ManifestFixture["pose"], lx: number, ly: number): XZ {
  const c = Math.cos(p.yaw_rad);
  const s = Math.sin(p.yaw_rad);
  const bx = p.x_m + lx * c - ly * s;
  const by = p.y_m + lx * s + ly * c;
  return [bx, -by];
}

export interface Station {
  id: string;
  /** Where a cook stands to work it. */
  stance: XZ;
  /** The way its front faces (people.ts convention: 0 faces +z), and how far in front of it the stance is, metres. */
  yaw: number;
  out: number;
}

export interface Floor {
  bounds: Rect;
  /** Obstacles to walking, without the chiller. */
  rects: Rect[];
  /** Where a chiller may not stand: corridors, hood footprints, every fixture's footprint. */
  keepOut: Rect[];
  stations: Map<string, Station>;
  /** The chiller as drawn. */
  chiller: { rect: Rect; stance: XZ; yaw: number; w: number; d: number };
  grid: NavGrid;
}

function rectOf(f: ManifestFixture): Rect | null {
  const b = f.world_aabb_m;
  if (!b) return null;
  return { minX: b.min_m[0], maxX: b.max_m[0], minZ: -b.max_m[1], maxZ: -b.min_m[1] };
}

export function buildFloor(m: FloorManifest, meta: FloorMeta | null, cell = 0.08): Floor {
  const W = m.room_shell.width_m;
  const D = m.room_shell.depth_m;
  const bounds: Rect = { minX: 0, maxX: W, minZ: -D, maxZ: 0 };
  const rects: Rect[] = [];
  const keepOut: Rect[] = [];
  const stations = new Map<string, Station>();
  let chiller: Floor["chiller"] | null = null;
  for (const f of m.fixtures) {
    const r = rectOf(f);
    if (!r || SMALL.test(f.archetype_id)) continue;
    const b = f.world_aabb_m!;
    keepOut.push(r);
    const stance = f.sockets?.find((s) => s.socket_id === "stance_front");
    if (stance)
      stations.set(f.fixture_id, { id: f.fixture_id, stance: toFloor(f.pose, stance.local_pos_m[0], stance.local_pos_m[1]), yaw: f.pose.yaw_rad, out: -stance.local_pos_m[1] });
    if (/BLAST-CHILLER/.test(f.archetype_id)) {
      const [w, d] = f.local_extents_m;
      chiller = { rect: r, stance: stations.get(f.fixture_id)!.stance, yaw: f.pose.yaw_rad, w, d };
      continue;
    }
    // Hoods hang above head height: they keep a chiller out but nobody walks into them.
    if (b.min_m[2] < 1.95) rects.push(r);
  }
  for (const c of meta?.walk_mesh_corridors ?? []) keepOut.push({ minX: c.min_m[0], maxX: c.max_m[0], minZ: -c.max_m[1], maxZ: -c.min_m[1] });
  if (!chiller) throw new Error("the plan has no blast chiller");
  return { bounds, rects, keepOut, stations, chiller, grid: new NavGrid(bounds, rects, cell) };
}

/** Cells of the grid within reach of a walker of radius r from a point: the walking distance to each, metres. */
export function distanceField(g: NavGrid, from: XZ, r: number): Float64Array {
  const { nx, nz, cell, clr } = g;
  const dist = new Float64Array(nx * nz).fill(Infinity);
  const free = (i: number, j: number) => i >= 0 && j >= 0 && i < nx && j < nz && clr[j * nx + i]! >= r;
  // Start from the free cells nearest the point (a stance can sit just inside the inflated edge of its fixture).
  const si = Math.floor((from[0] - g.x0) / cell);
  const sj = Math.floor((from[1] - g.z0) / cell);
  const heap = new MinHeap();
  for (let dj = -4; dj <= 4; dj++)
    for (let di = -4; di <= 4; di++) {
      const i = si + di;
      const j = sj + dj;
      if (!free(i, j)) continue;
      const d = Math.hypot(g.x0 + (i + 0.5) * cell - from[0], g.z0 + (j + 0.5) * cell - from[1]);
      const k = j * nx + i;
      if (d < dist[k]!) {
        dist[k] = d;
        heap.push(k, d);
      }
    }
  // 16 neighbours: the knight's moves cut the 8-neighbour error on long diagonals from 8% to under 3%.
  const STEPS: [number, number][] = [
    [1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [1, -1], [-1, 1], [-1, -1],
    [2, 1], [2, -1], [-2, 1], [-2, -1], [1, 2], [1, -2], [-1, 2], [-1, -2],
  ];
  const LEN = STEPS.map(([a, b]) => Math.hypot(a, b) * cell);
  while (heap.size) {
    const k = heap.pop();
    const d = heap.lastKey;
    if (d > dist[k]!) continue;
    const i = k % nx;
    const j = (k - i) / nx;
    for (let s = 0; s < 16; s++) {
      const [a, b] = STEPS[s]!;
      const ni = i + a;
      const nj = j + b;
      if (!free(ni, nj)) continue;
      // A long step passes the cells beside it: both must be free.
      if ((Math.abs(a) === 2 || Math.abs(b) === 2) && !(free(i + Math.sign(a), j + Math.sign(b)) && free(i + (Math.abs(a) === 2 ? Math.sign(a) : 0), j + (Math.abs(b) === 2 ? Math.sign(b) : 0)))) continue;
      if (Math.abs(a) === 1 && Math.abs(b) === 1 && !(free(i + a, j) && free(i, j + b))) continue;
      const nd = d + LEN[s]!;
      const nk = nj * nx + ni;
      if (nd < dist[nk]!) {
        dist[nk] = nd;
        heap.push(nk, nd);
      }
    }
  }
  return dist;
}

/** A field's value at a point: the nearest free cell's, plus the straight step to it. */
export function fieldAt(g: NavGrid, field: Float64Array, p: XZ): number {
  const si = Math.floor((p[0] - g.x0) / g.cell);
  const sj = Math.floor((p[1] - g.z0) / g.cell);
  let best = Infinity;
  for (let dj = -4; dj <= 4; dj++)
    for (let di = -4; di <= 4; di++) {
      const i = si + di;
      const j = sj + dj;
      if (i < 0 || j < 0 || i >= g.nx || j >= g.nz) continue;
      const v = field[j * g.nx + i]!;
      if (v === Infinity) continue;
      const d = v + Math.hypot(g.x0 + (i + 0.5) * g.cell - p[0], g.z0 + (j + 0.5) * g.cell - p[1]);
      if (d < best) best = d;
    }
  return best;
}

class MinHeap {
  private k: number[] = [];
  private v: number[] = [];
  lastKey = 0;
  get size(): number {
    return this.k.length;
  }
  push(item: number, key: number): void {
    const k = this.k;
    const v = this.v;
    let i = k.length;
    k.push(key);
    v.push(item);
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (k[p]! <= key) break;
      k[i] = k[p]!;
      v[i] = v[p]!;
      i = p;
    }
    k[i] = key;
    v[i] = item;
  }
  pop(): number {
    const k = this.k;
    const v = this.v;
    const top = v[0]!;
    this.lastKey = k[0]!;
    const lk = k.pop()!;
    const lv = v.pop()!;
    const n = k.length;
    if (n) {
      let i = 0;
      for (;;) {
        const l = 2 * i + 1;
        if (l >= n) break;
        const r = l + 1;
        const c = r < n && k[r]! < k[l]! ? r : l;
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

/** A place on a wall where the chiller fits, facing the room. */
export interface Spot {
  rect: Rect;
  stance: XZ;
  /** Heading of the chiller's front (people.ts convention: 0 faces +z). */
  yaw: number;
  /** Along the walls, metres from the room's front-left corner, clockwise seen from above: for drawing the strip. */
  along: number;
}

function overlaps(a: Rect, b: Rect, m: number): boolean {
  return a.minX < b.maxX + m && a.maxX > b.minX - m && a.minZ < b.maxZ + m && a.maxZ > b.minZ - m;
}

/** Every spot on the four walls, every `step` metres, where the chiller fits and a cook can stand at its door. */
export function wallSpots(f: Floor, step = 0.1): Spot[] {
  const { w, d } = f.chiller;
  const { minX, maxX, minZ, maxZ } = f.bounds;
  const out: Spot[] = [];
  const W = maxX - minX;
  const D = maxZ - minZ;
  // Four walls; the chiller's back to the wall, its front `d` into the room, `w` along the wall.
  const walls: { len: number; make: (s: number) => { rect: Rect; stance: XZ; yaw: number }; base: number }[] = [
    // Front wall (z = 0), facing -z.
    { len: W, base: 0, make: (s) => ({ rect: { minX: minX + s, maxX: minX + s + w, minZ: maxZ - d, maxZ }, stance: [minX + s + w / 2, maxZ - d - 0.5], yaw: Math.PI }) },
    // Right wall (x = W), facing -x.
    { len: D, base: W, make: (s) => ({ rect: { minX: maxX - d, maxX, minZ: maxZ - s - w, maxZ: maxZ - s }, stance: [maxX - d - 0.5, maxZ - s - w / 2], yaw: -Math.PI / 2 }) },
    // Back wall (z = -D), facing +z.
    { len: W, base: W + D, make: (s) => ({ rect: { minX: maxX - s - w, maxX: maxX - s, minZ, maxZ: minZ + d }, stance: [maxX - s - w / 2, minZ + d + 0.5], yaw: 0 }) },
    // Left wall (x = 0), facing +x.
    { len: D, base: 2 * W + D, make: (s) => ({ rect: { minX, maxX: minX + d, minZ: minZ + s, maxZ: minZ + s + w }, stance: [minX + d + 0.5, minZ + s + w / 2], yaw: Math.PI / 2 }) },
  ];
  for (const wall of walls) {
    for (let s = 0.05; s + w <= wall.len - 0.05 + 1e-9; s += step) {
      const c = wall.make(s);
      if (f.keepOut.some((k) => overlaps(c.rect, k, 0.03))) continue;
      if (f.grid.clearanceLB(c.stance[0], c.stance[1]) < CARRY_R) continue;
      out.push({ ...c, along: wall.base + s + w / 2 });
    }
  }
  return out;
}

/** The fixtures each station's pans come from: a combi dish is split evenly across the three combis, a prep dish across
 *  the six prep tables. */
export const SOURCES: Record<Dish["station"], string[]> = {
  kettle: ["kettle-40gal-1"],
  skillet: ["tilt-skillet-40gal-1"],
  combi: ["combi-rollin-1", "combi-rollin-2", "combi-rollin-3"],
  prep: ["preptable-prep-1", "preptable-prep-2", "preptable-prep-3", "preptable-prep-4", "preptable-prep-5", "preptable-prep-6"],
};
export const WALKIN = "walkin-cooler-14x16";

/** The week's cook-chill traffic: pans to the chiller per source fixture, and chilled loads to the walk-in. */
export interface Traffic {
  pans: Map<string, number>;
  loads: number;
  totalPans: number;
}

export function trafficOf(days: readonly { batches: readonly Batch[] }[], menu: readonly Dish[]): Traffic {
  const pans = new Map<string, number>();
  let loads = 0;
  let total = 0;
  for (const day of days) {
    let chilled = 0;
    for (const b of day.batches) {
      const dish = menu.find((m) => m.key === b.dish);
      if (!dish) throw new Error(`a batch of ${b.dish}, which is not on the menu`);
      const src = SOURCES[dish.station];
      for (const id of src) pans.set(id, (pans.get(id) ?? 0) + b.pans / src.length);
      chilled += b.pans;
      total += b.pans;
    }
    loads += Math.ceil(chilled / CHILL_LOAD - 1e-9);
  }
  return { pans, loads, totalPans: total };
}

/** Distance fields from the walk-in and from every source, for ranking spots. */
export interface Fields {
  walkin: Float64Array;
  sources: Map<string, Float64Array>;
}
export function fieldsOf(f: Floor): Fields {
  const sources = new Map<string, Float64Array>();
  for (const ids of Object.values(SOURCES)) for (const id of ids) sources.set(id, distanceField(f.grid, f.stations.get(id)!.stance, CARRY_R));
  return { walkin: distanceField(f.grid, f.stations.get(WALKIN)!.stance, CARRY_R), sources };
}

/** The week's walking, metres, with the chiller's door at `stance`: every pan there and back, every load there and back. */
export function weekMetres(f: Floor, fields: Fields, t: Traffic, stance: XZ): number {
  let m = 0;
  for (const [id, n] of t.pans) if (n > 0) m += 2 * n * fieldAt(f.grid, fields.sources.get(id)!, stance);
  if (t.loads > 0) m += 2 * t.loads * fieldAt(f.grid, fields.walkin, stance);
  return m;
}

export interface Route {
  from: string;
  to: string;
  path: XZ[];
  metres: number;
  /** Trips a week on this route, each there and back. */
  trips: number;
}

function pathLength(p: readonly XZ[]): number {
  let s = 0;
  for (let i = 1; i < p.length; i++) s += Math.hypot(p[i]![0] - p[i - 1]![0], p[i]![1] - p[i - 1]![1]);
  return s;
}

/** The routes with the chiller standing at a spot (and in the way), found on the floor by A*. */
export function routesFor(f: Floor, t: Traffic, spot: { rect: Rect; stance: XZ }): { routes: Route[]; metres: number } {
  const grid = new NavGrid(f.bounds, [...f.rects, spot.rect], f.grid.cell);
  const routes: Route[] = [];
  const add = (from: string, a: XZ, to: string, b: XZ, trips: number) => {
    const path = grid.findPath(a, b, CARRY_R, { keepRight: 0 }) ?? [a, b];
    routes.push({ from, to, path, metres: pathLength(path), trips });
  };
  for (const [id, n] of t.pans) if (n > 0) add(id, f.stations.get(id)!.stance, "chiller", spot.stance, n);
  add("chiller", spot.stance, WALKIN, f.stations.get(WALKIN)!.stance, t.loads);
  return { routes, metres: routes.reduce((s, r) => s + 2 * r.trips * r.metres, 0) };
}

export interface FloorAnswer {
  spots: Spot[];
  /** Each spot's week of walking, metres (distance-field estimate, for ranking and the heat strip). */
  cost: number[];
  best: number;
  drawn: { routes: Route[]; metres: number };
  atBest: { routes: Route[]; metres: number };
}

export function placeChiller(f: Floor, fields: Fields, t: Traffic, spots: readonly Spot[]): FloorAnswer {
  const cost = spots.map((s) => weekMetres(f, fields, t, s.stance));
  let best = 0;
  cost.forEach((c, i) => {
    if (c < cost[best]!) best = i;
  });
  return { spots: [...spots], cost, best, drawn: routesFor(f, t, f.chiller), atBest: routesFor(f, t, spots[best]!) };
}
