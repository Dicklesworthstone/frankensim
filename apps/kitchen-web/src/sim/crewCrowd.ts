/**
 * How many cooks fit in one kitchen (chapter 07, second figure).
 *
 * One commissary hot line under one hood: two combi ovens, two six-burner ranges, a tilting braising pan and a steam
 * kettle. Opposite it, a prep island; cold storage at the left end, the pack-out at the right. A crew of 1 to 8 cooks
 * each turns out hotel pans in a loop: take the food from a reach-in, prep it at the island, cook it on the line,
 * pack it at the pack-out, and back to the reach-ins. Each cook has a place of their own at the island and on the
 * line; the reach-in doors and the pack-out are shared, first come first served.
 *
 * Every walk is a leg planned by the crowd planner (three/models/crowd.ts) on the cook's own routine
 * (three/models/people.ts), as for the cooks in the hero kitchen. A leg leaves only when its whole path stays clear of
 * everyone else's body for as long as it walks and while it stands at its goal, so a cook who finds the way blocked
 * waits where they are and asks again. Cooks at work stand in the aisle, so a narrow aisle with cooks back to back
 * at the line and the island lets nobody past.
 *
 * The same crew run as ghosts, who walk through each other, sets what the crew would turn out with nobody in anybody's
 * way: their only waits are for a free reach-in door or pack-out place. The difference is what bodies cost.
 *
 * Floor frame: three.js x and z, metres; the hood wall at z = 0, the room towards +z. Yaw as in people.ts: 0 faces +z.
 */
import { Crowd, NavGrid, STATIC_MARGIN, type Footprint, type LegCandidate, type Motion, type Rect, type XZ } from "../three/models/crowd";
import { Routine, planningBody, type Cohort, type RoutineStep, type Staff, type WorkKind } from "../three/models/people";

/** The planner's time step, seconds. */
export const DT = 0.1;
/** Walking pace by cohort, m/s (the hero kitchen's). */
export const COHORT_SPEED: Record<Cohort, number> = { p05: 1.02, p50: 1.08, p95: 1.12 };
/** The crew's cohorts in the order cooks join: an average cook, a short one, a tall one, and round again. */
export const CREW_COHORTS: readonly Cohort[] = ["p50", "p05", "p95"];
/** The most cooks the figure tries. */
export const MAX_CREW = 8;

export type SpotKind = "cold" | "prep" | "hot" | "pack";
/** The loop a cook works, one hotel pan at a time. */
export const LOOP: readonly SpotKind[] = ["cold", "prep", "hot", "pack"];
/** Seconds of work at each kind of place, a pan at a time: take it out, prep it, cook it, pack it. */
export const DWELL: Record<SpotKind, number> = { cold: 8, prep: 40, hot: 45, pack: 8 };

/** A piece of equipment on the floor, by its front: the front's middle (fx, fz), the way it faces, and its size. */
export interface Piece {
  id: string;
  fx: number;
  fz: number;
  /** The way its front faces (people.ts: 0 faces +z). */
  yaw: number;
  w: number;
  d: number;
  h: number;
}
/** A place to work at a piece: `u` metres along its front from the middle (to the right of someone facing it). */
export interface Spot {
  kind: SpotKind;
  piece: number;
  u: number;
  work: WorkKind;
  /** How far into the piece the hands work, and at what height. */
  reach: number;
  surfaceY: number;
}
export interface Kitchen {
  key: string;
  name: string;
  bounds: Rect;
  pieces: Piece[];
  spots: Spot[];
  /** The island's prep places in the order cooks are given them (cook i gets the i-th). */
  prepOrder: number[];
}

/** A piece's floor rectangle. */
export function rectOf(p: Piece): Rect {
  const s = Math.round(Math.sin(p.yaw));
  const c = Math.round(Math.cos(p.yaw));
  if (s === 0) {
    const z0 = c > 0 ? p.fz - p.d : p.fz;
    return { minX: p.fx - p.w / 2, maxX: p.fx + p.w / 2, minZ: z0, maxZ: z0 + p.d };
  }
  const x0 = s > 0 ? p.fx - p.d : p.fx;
  return { minX: x0, maxX: x0 + p.d, minZ: p.fz - p.w / 2, maxZ: p.fz + p.w / 2 };
}

/** Where a body of walking radius `rWalk` stands to work a spot, facing it; and where the hands work. */
export function stanceOf(k: Kitchen, s: Spot, rWalk: number): { x: number; z: number; yaw: number; at: [number, number, number] } {
  const p = k.pieces[s.piece]!;
  const nx = Math.sin(p.yaw);
  const nz = Math.cos(p.yaw);
  // Along the front, to the right of someone facing the piece (who faces -n): right = (-nz, nx)... for yaw 0 that is +x.
  const tx = Math.cos(p.yaw);
  const tz = -Math.sin(p.yaw);
  const ex = p.fx + tx * s.u;
  const ez = p.fz + tz * s.u;
  const back = rWalk + STATIC_MARGIN + 0.016;
  return { x: ex + nx * back, z: ez + nz * back, yaw: p.yaw + Math.PI, at: [ex - nx * s.reach, s.surfaceY, ez - nz * s.reach] };
}

/* ---------------------------------------------------------------- the kitchens ---------------------------------------------------------------- */

const COMBI = { id: "NSF-COMBI-10PAN", w: 0.847, d: 0.776, h: 1.79 };
const RANGE = { id: "NSF-RANGE-GAS-6B", w: 0.914, d: 0.864, h: 1.473 };
const TILT = { id: "NSF-TILT-SKILLET-40GAL", w: 1.219, d: 0.991, h: 1.067 };
const KETTLE = { id: "NSF-KETTLE-40GAL", w: 0.914, d: 0.965, h: 1.067 };
const TABLE96 = { id: "NSF-PREPTABLE-96", w: 2.438, d: 0.762, h: 0.914 };
const TABLE48 = { id: "NSF-PREPTABLE-48", w: 1.219, d: 0.762, h: 0.914 };
const REACHIN = { id: "NSF-REACHIN-2D", w: 1.372, d: 0.813, h: 2.057 };
const CARRIER = { id: "NSF-CARRIER-INSULATED", w: 0.46, d: 0.635, h: 0.625 };
export const HOOD = { id: "NSF-HOOD-TYPE1-CANOPY-288IN", w: 7.315, d: 1.219, h: 0.61 };

/** The line's fronts, the room's walls. */
export const LINE_Z = 1.0;
const ROOM: Rect = { minX: -4.5, maxX: 4.5, minZ: 0, maxZ: 4.4 };
const GAP = 0.05;

/** The hot line under the hood: the same in every layout (the hood does not move). */
function line(): { pieces: Piece[]; spots: Spot[] } {
  const order = [COMBI, COMBI, RANGE, RANGE, TILT, KETTLE];
  const total = order.reduce((a, p) => a + p.w, 0) + GAP * (order.length - 1);
  let x = -total / 2 - 0.02;
  const pieces: Piece[] = [];
  for (const p of order) {
    pieces.push({ id: p.id, fx: x + p.w / 2, fz: LINE_Z, yaw: 0, w: p.w, d: p.d, h: p.h });
    x += p.w + GAP;
  }
  // One cook a piece: two side by side at one range would stand shoulder into shoulder.
  const spots: Spot[] = [
    { kind: "hot", piece: 0, u: 0, work: "place", reach: 0.25, surfaceY: 1.02 },
    { kind: "hot", piece: 1, u: 0, work: "place", reach: 0.25, surfaceY: 1.02 },
    { kind: "hot", piece: 2, u: 0, work: "saute", reach: 0.2, surfaceY: 0.914 },
    { kind: "hot", piece: 3, u: 0, work: "stir", reach: 0.2, surfaceY: 0.914 },
    { kind: "hot", piece: 4, u: 0, work: "stir", reach: 0.35, surfaceY: 1.0 },
    { kind: "hot", piece: 5, u: 0, work: "stir", reach: 0.4, surfaceY: 1.067 },
  ];
  return { pieces, spots };
}

/** Cold storage on the left wall, the pack-out on the right: the same in every layout. */
function ends(first: number): { pieces: Piece[]; spots: Spot[] } {
  const fx = ROOM.minX + REACHIN.d;
  const pieces: Piece[] = [
    { id: REACHIN.id, fx, fz: 1.95, yaw: Math.PI / 2, w: REACHIN.w, d: REACHIN.d, h: REACHIN.h },
    { id: REACHIN.id, fx, fz: 1.95 + REACHIN.w + 0.04, yaw: Math.PI / 2, w: REACHIN.w, d: REACHIN.d, h: REACHIN.h },
    { id: TABLE48.id, fx: ROOM.maxX - TABLE48.d, fz: 2.2, yaw: -Math.PI / 2, w: TABLE48.w, d: TABLE48.d, h: TABLE48.h },
    { id: CARRIER.id, fx: ROOM.maxX - CARRIER.d, fz: 3.25, yaw: -Math.PI / 2, w: CARRIER.w, d: CARRIER.d, h: CARRIER.h },
    { id: CARRIER.id, fx: ROOM.maxX - CARRIER.d, fz: 3.75, yaw: -Math.PI / 2, w: CARRIER.w, d: CARRIER.d, h: CARRIER.h },
  ];
  const spots: Spot[] = [
    { kind: "cold", piece: first, u: -0.343, work: "pick", reach: 0.2, surfaceY: 1.1 },
    { kind: "cold", piece: first, u: 0.343, work: "pick", reach: 0.2, surfaceY: 1.1 },
    { kind: "cold", piece: first + 1, u: -0.343, work: "pick", reach: 0.2, surfaceY: 1.1 },
    { kind: "cold", piece: first + 1, u: 0.343, work: "pick", reach: 0.2, surfaceY: 1.1 },
    { kind: "pack", piece: first + 2, u: -0.3, work: "place", reach: 0.3, surfaceY: TABLE48.h },
    { kind: "pack", piece: first + 2, u: 0.3, work: "place", reach: 0.3, surfaceY: TABLE48.h },
  ];
  return { pieces, spots };
}

/**
 * The island: `aisle` metres of hot aisle from the line's fronts, its tables `split` metres apart in the middle (0:
 * end to end). Four prep places face the line across the hot aisle and four the far wall.
 */
function island(first: number, aisle: number, split: number): { pieces: Piece[]; spots: Spot[] } {
  const fz = LINE_Z + aisle;
  const off = TABLE96.w / 2 + split / 2;
  const pieces: Piece[] = [-1, 1].map((sx) => ({ id: TABLE96.id, fx: sx * off, fz, yaw: Math.PI, w: TABLE96.w, d: TABLE96.d, h: TABLE96.h }));
  const spots: Spot[] = [];
  // The near side (facing the line): two places a table.
  for (const [i, u] of [
    [0, -0.6],
    [0, 0.6],
    [1, -0.6],
    [1, 0.6],
  ] as const)
    spots.push({ kind: "prep", piece: first + i, u, work: "chop", reach: 0.3, surfaceY: TABLE96.h });
  return { pieces, spots };
}

/** The far side of an island table: a piece facing +z over the same floor (only its spots are used). */
function farSide(first: number, tables: readonly Piece[]): { pieces: Piece[]; spots: Spot[] } {
  const pieces: Piece[] = tables.map((t) => ({ ...t, id: "", fz: t.fz + t.d, yaw: 0 }));
  const spots: Spot[] = [];
  for (const [i, u] of [
    [0, -0.6],
    [0, 0.6],
    [1, -0.6],
    [1, 0.6],
  ] as const)
    spots.push({ kind: "prep", piece: first + i, u, work: "chop", reach: 0.3, surfaceY: TABLE96.h });
  return { pieces, spots };
}

/** A kitchen: the line and both ends fixed, the island where the layout puts it. */
function kitchen(key: string, name: string, aisle: number, split: number, prepOrder: number[]): Kitchen {
  const L = line();
  const E = ends(L.pieces.length);
  const I = island(L.pieces.length + E.pieces.length, aisle, split);
  const F = farSide(L.pieces.length + E.pieces.length + I.pieces.length, I.pieces);
  const pieces = [...L.pieces, ...E.pieces, ...I.pieces, ...F.pieces];
  const spots = [...L.spots, ...E.spots, ...I.spots, ...F.spots];
  const prep = spots.flatMap((s, i) => (s.kind === "prep" ? [i] : []));
  return { key, name, bounds: ROOM, pieces, spots, prepOrder: prepOrder.map((i) => prep[i]!) };
}

/** As built: one long island 42 in from the line, its tables end to end. */
export const AS_BUILT = kitchen("built", "As built", 1.067, 0, [1, 2, 0, 3, 5, 6, 4, 7]);
/** Same hood, the same island pulled back 18 in: 60 in from the line, the aisle behind it 18 in narrower. */
export const MOVED = kitchen("moved", "Island pulled back", 1.524, 0, [1, 2, 0, 3, 5, 6, 4, 7]);
export const KITCHENS: Record<string, Kitchen> = { built: AS_BUILT, moved: MOVED };

/** The floor's obstacles: every real piece (the far-side stand-ins are the same tables). */
export function obstacles(k: Kitchen): Rect[] {
  return k.pieces.filter((p) => p.id).map(rectOf);
}

/* ---------------------------------------------------------------- the crew ---------------------------------------------------------------- */

/** A body's footprints (crowd.ts shapes): standing, three discs across the shoulders; walking, one. */
export function footprintsOf(body: Staff): { shapes: Footprint[]; rWalk: number } {
  const sb = body.dims.shoulder_breadth_m;
  const cd = body.dims.chest_depth_m;
  const rr = cd / 2 + 0.02;
  const half = sb / 2;
  const rWalk = half + 0.015;
  const stand: Footprint = [
    { x: -(half - rr), z: 0, r: rr },
    { x: 0, z: 0, r: rr },
    { x: half - rr, z: 0, r: rr },
  ];
  return { shapes: [stand, [{ x: 0, z: 0, r: rWalk }]], rWalk };
}
const SH_STAND = 0;
const SH_WALK = 1;

/** What a leg is for: one of the loop's places, or stepping out of someone's way. */
export type LegKind = SpotKind | "aside";

/** One leg of a cook's day: when they were ready and when they went, the way they walked, and the work at the end. */
export interface Leg {
  /** When they were ready to go, and when they went (s); between, they stood waiting. */
  ready: number;
  go: number;
  from: { x: number; z: number; yaw: number };
  path: [number, number][];
  speed: number;
  yaw: number;
  spot: number;
  kind: LegKind;
  work: WorkKind;
  at: [number, number, number];
  surfaceY: number;
  dwell: number;
  /** When they got there (s), and whether the wait before it was for someone in the way (else for a free place). */
  arrive: number;
  blocked: boolean;
}

export interface CookDay {
  cohort: Cohort;
  legs: Leg[];
  /** Seconds standing while someone was in the way, walking out of someone's way, waiting for a free place, walking,
   *  and at work. */
  waitBody: number;
  aside: number;
  waitPlace: number;
  walk: number;
  work: number;
  pans: number;
}

export interface CrowdRun {
  kitchen: string;
  crew: number;
  ghosts: boolean;
  seconds: number;
  /** Pans packed between `from` and `seconds`, and that rate an hour. */
  from: number;
  pans: number;
  perHour: number;
  cooks: CookDay[];
  /** Where cooks stood while someone was in their way, from `from` on: (x, z, seconds). */
  waits: [number, number, number][];
  /** The least gap between two bodies over the whole run, metres: never below zero (null for ghosts). */
  minGap: number | null;
  /** Whether the crew ever stood for two minutes with nobody able to move. */
  gridlock: boolean;
}

interface Cook {
  id: number;
  cohort: Cohort;
  body: Staff;
  rWalk: number;
  rPath: number;
  speed: number;
  prep: number;
  phase: number;
  spot: number;
  pose: { x: number; z: number; yaw: number };
  readyK: number;
  /** The step they first tried the leg they are waiting to make (-1: not waiting). */
  firstTryK: number;
  blocked: boolean;
  day: CookDay;
}

/** The routine steps of a leg: walk the path, turn to the work, and work. `at` is where the hands work. */
export function legSteps<V>(l: Pick<Leg, "path" | "speed" | "yaw" | "work" | "dwell" | "surfaceY">, at: V): RoutineStep[] {
  const steps: RoutineStep[] = [];
  if (l.path.length) steps.push({ kind: "walk", path: l.path, speed: l.speed });
  steps.push({ kind: "turn", yaw: l.yaw });
  steps.push({ kind: "work", work: l.work, duration: l.dwell, opts: { at: at as never, surfaceY: l.surfaceY } });
  return steps;
}

/** A leg's motion for the planner: the routine's root sampled every DT up to the end of its turn. */
function motionOf(body: Staff, from: { x: number; z: number; yaw: number }, path: [number, number][], speed: number, yaw: number): { m: Motion; walkS: number } {
  const steps: RoutineStep[] = [];
  if (path.length) steps.push({ kind: "walk", path, speed });
  steps.push({ kind: "turn", yaw });
  const r = new Routine(body, from, steps, false);
  const n = Math.ceil(r.duration / DT) + 1;
  const m: Motion = { n, x: new Float32Array(n), z: new Float32Array(n), yaw: new Float32Array(n), shape: new Uint8Array(n) };
  for (let i = 0; i < n; i++) {
    const t = Math.min(i * DT, r.duration);
    const p = r.rootAt(t);
    m.x[i] = p.x;
    m.z[i] = p.z;
    m.yaw[i] = p.yaw;
    m.shape[i] = r.activityAt(t) === "walk" ? SH_WALK : SH_STAND;
  }
  return { m, walkS: r.duration };
}

/** Distance from a point to a polyline. */
function toPolyline(pts: readonly XZ[], x: number, z: number): number {
  let d = Infinity;
  for (let i = 0; i + 1 < pts.length; i++) {
    const [ax, az] = pts[i]!;
    const [bx, bz] = pts[i + 1]!;
    const ex = bx - ax;
    const ez = bz - az;
    const l2 = ex * ex + ez * ez;
    const t = l2 > 0 ? Math.max(0, Math.min(1, ((x - ax) * ex + (z - az) * ez) / l2)) : 0;
    d = Math.min(d, Math.hypot(x - (ax + t * ex), z - (az + t * ez)));
  }
  return pts.length === 1 ? Math.hypot(x - pts[0]![0], z - pts[0]![1]) : d;
}

/** Seconds a cook waits on someone in the way before asking whoever stands idle in their way to step aside. */
const MAKE_WAY_S = 2;
/** A step-aside: the short stand at the place stepped into before trying again. */
const ASIDE_DWELL = 0.6;
/** A run's length and the start-up left out of its count, seconds; the floor grid's cell, metres. */
export const RUN_S = 1200;
export const FROM_S = 240;
export const CELL = 0.08;
/** Seconds between one cook starting their first pan and the next. */
export const START_GAP_S = 7;
/**
 * The mornings a setting is run for: the cooks starting this many seconds apart. One run's loss swings with when the
 * cooks happen to meet (eight cooks as built lose 13% to 23% across starts 4 to 11 s apart); three pooled hold to
 * within a few points. The middle one is the morning the figure plays.
 */
export const MORNINGS: readonly number[] = [5, START_GAP_S, 9];

/**
 * One crew's run: `crew` cooks in kitchen `k` for `seconds`, counting pans packed after `from` seconds. With `ghosts`
 * the cooks walk through each other (the same walks, never delayed, nobody ever in the way).
 */
export function runCrowd(k: Kitchen, crew: number, seconds = RUN_S, from = FROM_S, ghosts = false, cell = CELL, stagger = START_GAP_S): CrowdRun {
  const grid = new NavGrid(k.bounds, obstacles(k), cell);
  const crowd = new Crowd(grid, DT);
  // Who has each place: the cook at it or on the way to it.
  const claimed = new Map<number, number>();
  const cooks: Cook[] = [];
  for (let i = 0; i < crew; i++) {
    const cohort = CREW_COHORTS[i % CREW_COHORTS.length]!;
    const body = planningBody(cohort);
    const { shapes, rWalk } = footprintsOf(body);
    const prep = k.prepOrder[i]!;
    const p = stanceOf(k, k.spots[prep]!, rWalk);
    const id = crowd.addAgent({ shapes }, p.x, p.z, p.yaw, SH_STAND, 0);
    // Everyone starts at their own prep place, a pan in hand, a few seconds apart.
    const startK = Math.round((i * stagger) / DT);
    cooks.push({
      id,
      cohort,
      body,
      rWalk,
      rPath: rWalk + STATIC_MARGIN + 0.012,
      speed: COHORT_SPEED[cohort],
      prep,
      phase: 1,
      spot: prep,
      pose: { x: p.x, z: p.z, yaw: p.yaw },
      readyK: startK + Math.round(DWELL.prep / DT),
      firstTryK: -1,
      blocked: false,
      day: { cohort, legs: [], waitBody: 0, aside: 0, waitPlace: 0, walk: 0, work: DWELL.prep, pans: 0 },
    });
    claimed.set(prep, i);
  }
  const endK = Math.round(seconds / DT);
  const fromK = Math.round(from / DT);
  const waits: [number, number, number][] = [];
  let pans = 0;
  let minGap: number | null = ghosts ? null : Infinity;
  let audited = 0;
  let lastMoveK = 0;
  let gridlock = false;

  /** The ways from a pose to a goal, each at a pace and a little slower, as planner candidates. */
  const candidatesTo = (c: Cook, k0: number, goal: { x: number; z: number; yaw: number }) => {
    const here: XZ = [c.pose.x, c.pose.z];
    const near = Math.hypot(goal.x - here[0], goal.z - here[1]) < 0.02;
    const found = near ? [[here]] : grid.findPaths(here, [goal.x, goal.z], c.rPath, [0.75, 1, 0.4], ghosts ? undefined : crowd.standingObstacles(c.id, k0));
    const cands: (LegCandidate & { path: [number, number][]; walkS?: number })[] = [];
    let variant = 0;
    for (const p of found) {
      variant++;
      if (!p) continue;
      const path = p.slice(1).map(([x, z]) => [x, z] as [number, number]);
      if (cands.some((q) => q.path.length === path.length && q.path.every((pt, i) => Math.abs(pt[0] - path[i]![0]) < 1e-3 && Math.abs(pt[1] - path[i]![1]) < 1e-3))) continue;
      for (const v of path.length ? [c.speed, c.speed * 0.8] : [c.speed]) {
        const cand: LegCandidate & { path: [number, number][]; walkS?: number } = { speed: v, variant, path };
        cand.realize = () => {
          const r = motionOf(c.body, c.pose, path, v, goal.yaw);
          cand.walkS = r.walkS;
          return r.m;
        };
        cands.push(cand);
      }
    }
    return cands;
  };

  /** Commit a planned leg to a cook's day and move them on. */
  const take = (c: Cook, go: number, cand: { path: [number, number][]; speed: number; walkS?: number }, target: number, kind: LegKind, goal: ReturnType<typeof stanceOf>, dwell: number) => {
    const spot = k.spots[target]!;
    const walkS = cand.walkS ?? 0;
    const arrive = go * DT + walkS;
    const waited = (go - c.firstTryK) * DT;
    const blocked = c.blocked || go > c.readyK;
    if (blocked) c.day.waitBody += waited;
    else c.day.waitPlace += waited;
    if (waited > 0 && blocked && go > fromK) waits.push([c.pose.x, c.pose.z, waited]);
    c.day.legs.push({
      ready: c.firstTryK * DT,
      go: go * DT,
      from: { ...c.pose },
      path: cand.path,
      speed: cand.speed,
      yaw: goal.yaw,
      spot: target,
      kind,
      work: kind === "aside" ? "stand" : spot.work,
      at: goal.at,
      surfaceY: spot.surfaceY,
      dwell,
      arrive,
      blocked,
    });
    if (kind === "aside") c.day.aside += walkS;
    else {
      c.day.walk += walkS;
      c.day.work += dwell;
    }
    if (claimed.get(c.spot) === c.id) claimed.delete(c.spot);
    claimed.set(target, c.id);
    c.spot = target;
    c.pose = { x: goal.x, z: goal.z, yaw: goal.yaw };
    // The first step on or after the work is done: never before it.
    c.readyK = Math.ceil((arrive + dwell) / DT - 1e-6);
    lastMoveK = go;
  };

  /**
   * A cook kept waiting on someone: whoever stands idle (done, waiting to go themselves) on the way they would take
   * with nobody there steps into the nearest free place clear of it, if they can get there.
   */
  const makeWay = (c: Cook, k0: number, goal: { x: number; z: number }, target: number): boolean => {
    const way = grid.findPath([c.pose.x, c.pose.z], [goal.x, goal.z], c.rPath);
    if (!way) return false;
    const inWay = (x: number, z: number, r: number) => toPolyline(way, x, z) < r + c.rPath + 0.05;
    const idle = cooks.filter((q) => q !== c && q.readyK <= k0 && inWay(q.pose.x, q.pose.z, q.rWalk)).sort((a, b) => Math.hypot(a.pose.x - c.pose.x, a.pose.z - c.pose.z) - Math.hypot(b.pose.x - c.pose.x, b.pose.z - c.pose.z));
    for (const y of idle) {
      const refuges = k.spots
        .map((s, i) => ({ i, p: stanceOf(k, s, y.rWalk) }))
        .filter((r) => r.i !== target && !claimed.has(r.i) && !inWay(r.p.x, r.p.z, y.rWalk))
        .sort((a, b) => Math.hypot(a.p.x - y.pose.x, a.p.z - y.pose.z) - Math.hypot(b.p.x - y.pose.x, b.p.z - y.pose.z))
        .slice(0, 4);
      for (const r of refuges) {
        const ky = Math.max(crowd.end(y.id) + 1, k0);
        const cands = candidatesTo(y, ky, r.p);
        if (!cands.length) continue;
        const v = crowd.planLeg(y.id, ky, cands, 2);
        if (!v.ok) continue;
        if (y.firstTryK < 0) y.firstTryK = y.readyK;
        y.blocked = true;
        const chosen = cands.find((q) => q.variant === v.variant && q.speed === v.speed)!;
        take(y, ky + Math.round(v.delay / DT), chosen, r.i, "aside", r.p, ASIDE_DWELL);
        // Still on their way where they were going: they try again from here.
        y.firstTryK = y.readyK;
        return true;
      }
    }
    return false;
  };

  for (;;) {
    let next: Cook | null = null;
    for (const q of cooks) if (!next || q.readyK < next.readyK) next = q;
    if (!next || next.readyK >= endK) break;
    const c = next;
    const k0 = c.readyK;
    if (c.firstTryK < 0) c.firstTryK = k0;
    if (!ghosts && k0 - lastMoveK > 120 / DT) gridlock = true;
    // Where next: their own prep place, else the nearest free line place, reach-in door or pack-out place.
    const kind = LOOP[(c.phase + 1) % LOOP.length]!;
    let target = -1;
    if (kind === "prep") target = c.prep;
    else {
      let best = Infinity;
      k.spots.forEach((s, i) => {
        if (s.kind !== kind || claimed.has(i)) return;
        const p = stanceOf(k, s, c.rWalk);
        const d = Math.hypot(p.x - c.pose.x, p.z - c.pose.z);
        if (d < best) {
          best = d;
          target = i;
        }
      });
    }
    const retry = (body: boolean) => {
      c.readyK = k0 + Math.round(1 / DT);
      if (body) c.blocked = true;
      if (!ghosts) crowd.wait(c.id, c.readyK - 1);
    };
    if (target < 0) {
      retry(false);
      continue;
    }
    const goal = stanceOf(k, k.spots[target]!, c.rWalk);
    const cands = candidatesTo(c, k0, goal);
    let go = k0;
    let chosen: (typeof cands)[number] | undefined = cands[0];
    if (ghosts) chosen?.realize!();
    else if (chosen) {
      const v = crowd.planLeg(c.id, k0, cands, 6);
      if (v.ok) {
        go = k0 + Math.round(v.delay / DT);
        chosen = cands.find((q) => q.variant === v.variant && q.speed === v.speed)!;
      } else chosen = undefined;
    }
    if (!chosen) {
      if (k0 - c.firstTryK >= MAKE_WAY_S / DT) makeWay(c, k0, goal, target);
      retry(true);
      continue;
    }
    take(c, go, chosen, target, kind, goal, DWELL[kind]);
    c.phase = (c.phase + 1) % LOOP.length;
    c.firstTryK = -1;
    c.blocked = false;
    if (kind === "pack" && c.readyK <= endK) {
      c.day.pans++;
      if (c.readyK > fromK) pans++;
    }
    // Keep the tracks short, measuring first that nobody ever came within reach of anybody.
    if (!ghosts) {
      let lo = Infinity;
      for (const q of cooks) lo = Math.min(lo, q.readyK, crowd.end(q.id));
      if (lo - audited > 200) {
        minGap = Math.min(minGap!, crowd.audit(audited, lo).minPairGap);
        crowd.prune(lo);
        audited = lo;
      }
    }
  }
  // Whoever is still waiting at the end waited until then.
  for (const c of cooks)
    if (c.firstTryK >= 0 && c.firstTryK < endK) {
      const w = (endK - c.firstTryK) * DT;
      if (c.blocked) {
        c.day.waitBody += w;
        if (endK > fromK) waits.push([c.pose.x, c.pose.z, w]);
      } else c.day.waitPlace += w;
    }
  if (!ghosts) {
    let hi = 0;
    for (const q of cooks) hi = Math.max(hi, crowd.end(q.id));
    if (hi > audited) minGap = Math.min(minGap!, crowd.audit(audited, hi).minPairGap);
  }
  const span = Math.max(1, seconds - from);
  return {
    kitchen: k.key,
    crew,
    ghosts,
    seconds,
    from,
    pans,
    perHour: (pans * 3600) / span,
    cooks: cooks.map((c) => c.day),
    waits,
    minGap: minGap === Infinity ? null : minGap,
    gridlock,
  };
}

/** A setting's mornings pooled: pans an hour over all of them, and the morning the figure plays. */
export interface CrowdSet {
  kitchen: string;
  crew: number;
  ghosts: boolean;
  /** Pans packed in the counted spans of all the mornings, an hour. */
  perHour: number;
  /** Seconds each cook spent, on average, waiting on someone in the way or stepping out of it, an hour. */
  waitPerHour: number;
  /** Where cooks stood while someone was in their way, all the mornings. */
  waits: [number, number, number][];
  minGap: number | null;
  gridlock: boolean;
  /** The morning the figure plays (MORNINGS' middle one). */
  shown: CrowdRun;
}

export function runMornings(k: Kitchen, crew: number, ghosts = false): CrowdSet {
  const runs = MORNINGS.map((st) => runCrowd(k, crew, RUN_S, FROM_S, ghosts, CELL, st));
  const pans = runs.reduce((a, r) => a + r.pans, 0);
  const span = runs.reduce((a, r) => a + (r.seconds - r.from), 0);
  const waited = runs.reduce((a, r) => a + r.cooks.reduce((b, c) => b + c.waitBody + c.aside, 0), 0);
  const gaps = runs.map((r) => r.minGap).filter((g): g is number => g !== null);
  return {
    kitchen: k.key,
    crew,
    ghosts,
    perHour: (pans * 3600) / span,
    waitPerHour: (waited * 3600) / (crew * runs.reduce((a, r) => a + r.seconds, 0)),
    waits: runs.flatMap((r) => r.waits),
    minGap: gaps.length ? Math.min(...gaps) : null,
    gridlock: runs.some((r) => r.gridlock),
    shown: runs[MORNINGS.indexOf(START_GAP_S)]!,
  };
}
