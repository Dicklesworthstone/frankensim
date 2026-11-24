/**
 * The crowd's invariants on a small kitchen: two counter runs either side of a 1.2 m aisle, an island, and people
 * working through stations on both sides. Run with `bun test src/three/models/crowd.test.ts`.
 */
import { describe, expect, test } from "bun:test";
import { Crowd, NavGrid, PAIR_MARGIN, STATIC_MARGIN, wrapPi, type Footprint, type LegCandidate, type Motion, type Rect, type XZ } from "./crowd";

const DT = 0.05;

/** Standing: three circles across the shoulders. Walking: the shoulders' disc with the stride in front and behind. */
const STAND: Footprint = [
  { x: -0.1, z: 0.02, r: 0.13 },
  { x: 0, z: 0.02, r: 0.13 },
  { x: 0.1, z: 0.02, r: 0.13 },
];
const WALK: Footprint = [
  { x: 0, z: 0, r: 0.23 },
  { x: 0, z: 0.15, r: 0.12 },
  { x: 0, z: -0.15, r: 0.12 },
];
const SHAPES = [STAND, WALK];
/** The clearance a walker's path keeps: the shoulders' half-width and the margin (the stride is checked on the motion). */
const R_WALK = 0.23 + STATIC_MARGIN + 0.005;

/** A room 6 m by 8 m: counters along both long walls with a 1.2 m aisle in the middle third, and an island. */
const ROOM: Rect = { minX: 0, minZ: -8, maxX: 6, maxZ: 0 };
const COUNTERS: Rect[] = [
  { minX: 0, minZ: -6, maxX: 2.4, maxZ: -1.5 }, // west run (front face x = 2.4)
  { minX: 3.6, minZ: -6, maxX: 6, maxZ: -1.5 }, // east run (front face x = 3.6): the aisle is 1.2 m
  { minX: 2.0, minZ: -7.6, maxX: 4.0, maxZ: -7.0 }, // a table at the far end
];

/** Stations: where a person stands to work (0.3 m off a counter face), and which way they face. */
const STATIONS: { at: XZ; yaw: number }[] = [
  { at: [2.7, -2.0], yaw: -Math.PI / 2 }, // west run, facing -x
  { at: [2.7, -4.5], yaw: -Math.PI / 2 },
  { at: [3.3, -3.2], yaw: Math.PI / 2 }, // east run, facing +x
  { at: [3.3, -5.6], yaw: Math.PI / 2 },
  { at: [3.0, -6.7], yaw: Math.PI }, // the far table, facing -z
  { at: [1.0, -0.6], yaw: 0 }, // the room's front wall, facing +z
  { at: [5.0, -0.6], yaw: 0 },
];

/**
 * A plain walker: turn on the spot to the path's heading (3 rad/s), walk it at `speed`, turn to face the station,
 * and work there `dwell` seconds. Shape 1 while walking, 0 standing or turning.
 */
function realize(start: { x: number; z: number; yaw: number }, path: readonly XZ[], speed: number, face: number, dwell: number): Motion {
  const xs: number[] = [];
  const zs: number[] = [];
  const ys: number[] = [];
  const sh: number[] = [];
  let x = start.x;
  let z = start.z;
  let yaw = start.yaw;
  const push = (moving: boolean) => {
    xs.push(x);
    zs.push(z);
    ys.push(yaw);
    sh.push(moving ? 1 : 0);
  };
  // Turning on the spot is done standing (feet together), so it shows the standing footprint.
  const turnTo = (goal: number) => {
    let d = wrapPi(goal - yaw);
    while (Math.abs(d) > 1e-6) {
      const st = Math.sign(d) * Math.min(Math.abs(d), 3 * DT);
      yaw += st;
      d -= st;
      push(false);
    }
  };
  for (let k = 1; k < path.length; k++) {
    const [ax, az] = path[k - 1]!;
    const [bx, bz] = path[k]!;
    const L = Math.hypot(bx - ax, bz - az);
    if (L < 1e-6) continue;
    turnTo(Math.atan2(bx - ax, bz - az));
    const n = Math.ceil(L / (speed * DT));
    for (let q = 1; q <= n; q++) {
      x = ax + ((bx - ax) * q) / n;
      z = az + ((bz - az) * q) / n;
      push(true);
    }
  }
  turnTo(face);
  for (let q = 0; q < Math.round(dwell / DT); q++) push(false);
  if (!xs.length) push(false);
  return { n: xs.length, x: Float32Array.from(xs), z: Float32Array.from(zs), yaw: Float32Array.from(ys), shape: Uint8Array.from(sh) };
}

function kitchen(): { grid: NavGrid; crowd: Crowd } {
  const grid = new NavGrid(ROOM, COUNTERS, 0.04);
  return { grid, crowd: new Crowd(grid, DT) };
}

/** Plan a leg the way the hero does: path variants (keep right, hug the right, nearer the middle) at two speeds. */
function plan(crowd: Crowd, grid: NavGrid, id: number, k0: number, goal: (typeof STATIONS)[number], dwell: number) {
  const here = crowd.poseAt(id, crowd.end(id));
  const dyn = crowd.standingObstacles(id, k0);
  const cands: LegCandidate[] = [];
  let variant = 0;
  for (const keep of [0.75, 1, 0.4]) {
    const path = grid.findPath([here.x, here.z], goal.at, R_WALK, { dynamic: dyn, keepRight: keep });
    variant++;
    if (!path) continue;
    for (const v of [1.1, 0.8]) cands.push({ motion: realize(here, path, v, goal.yaw, dwell), speed: v, variant });
  }
  const verdict = crowd.planLeg(id, k0, cands, 10);
  if (process.env.CROWD_DEBUG) console.log("plan", id, k0, cands.length, JSON.stringify(verdict));
  return { verdict, path: cands.length > 0 };
}

describe("NavGrid", () => {
  test("clearance is exact at cell centres and a lower bound between them", () => {
    const { grid } = kitchen();
    for (let k = 0; k < 2000; k++) {
      const x = 0.05 + ((k * 0.6180339887) % 1) * 5.9;
      const z = -7.95 + ((k * 0.7548776662) % 1) * 7.9;
      expect(grid.clearanceLB(x, z)).toBeLessThanOrEqual(grid.clearance(x, z) + 1e-9);
    }
  });

  test("paths keep the walker's radius from every wall and counter, and reach their ends", () => {
    const { grid } = kitchen();
    for (const a of STATIONS)
      for (const b of STATIONS) {
        if (a === b) continue;
        const p = grid.findPath(a.at, b.at, R_WALK);
        expect(p).not.toBeNull();
        expect(p![0]).toEqual(a.at);
        expect(p![p!.length - 1]).toEqual(b.at);
        for (let k = 0; k + 1 < p!.length; k++) {
          const [ax, az] = p![k]!;
          const [bx, bz] = p![k + 1]!;
          const n = Math.ceil(Math.hypot(bx - ax, bz - az) / 0.01);
          for (let q = 0; q <= n; q++) expect(grid.clearance(ax + ((bx - ax) * q) / n, az + ((bz - az) * q) / n)).toBeGreaterThanOrEqual(R_WALK - 1e-6);
        }
      }
  });

  test("no path through a gap narrower than the body", () => {
    const grid = new NavGrid({ minX: 0, minZ: -4, maxX: 4, maxZ: 0 }, [
      { minX: 0, minZ: -2.1, maxX: 1.8, maxZ: -1.9 },
      { minX: 2.2, minZ: -2.1, maxX: 4, maxZ: -1.9 },
    ]);
    // The gap is 0.4 m: a body 0.48 m across cannot pass, one 0.36 m across can.
    expect(grid.findPath([2, -0.5], [2, -3.5], 0.24)).toBeNull();
    expect(grid.findPath([2, -0.5], [2, -3.5], 0.18)).not.toBeNull();
  });
});

describe("Crowd", () => {
  test("two people meeting in the aisle pass each other without waiting, on their own right", () => {
    const { grid, crowd } = kitchen();
    const a = crowd.addAgent({ shapes: SHAPES }, 1.5, -0.7, Math.PI, 0);
    const b = crowd.addAgent({ shapes: SHAPES }, 1.5, -6.6, 0, 0);
    const va = plan(crowd, grid, a, 1, { at: [4.5, -6.6], yaw: Math.PI }, 1);
    // b plans after a is committed, starting at the same moment, through the same aisle the other way.
    const vb = plan(crowd, grid, b, 1, { at: [4.5, -0.7], yaw: 0 }, 1);
    expect(va.verdict.ok).toBe(true);
    expect(vb.verdict.ok).toBe(true);
    expect(vb.verdict.ok && vb.verdict.delay).toBe(0);
    // They meet mid-aisle side by side: find the step where their z are equal and check they are apart in x.
    let met = false;
    for (let k = 1; k < Math.max(crowd.end(a), crowd.end(b)); k++) {
      const pa = crowd.poseAt(a, k);
      const pb = crowd.poseAt(b, k);
      if (Math.abs(pa.z - pb.z) < 0.05 && pa.z < -1.6 && pa.z > -5.9) {
        met = true;
        expect(Math.abs(pa.x - pb.x)).toBeGreaterThan(0.45);
      }
    }
    expect(met).toBe(true);
    const audit = crowd.audit(0, Math.max(crowd.end(a), crowd.end(b)));
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
  });

  test("a shift of five people working stations: no two bodies closer than their radii, nobody in a counter, every leg ends at its station", () => {
    const { grid, crowd } = kitchen();
    const starts = [0, 2, 4, 5, 6];
    const ids = starts.map((s) => crowd.addAgent({ shapes: SHAPES }, STATIONS[s]!.at[0], STATIONS[s]!.at[1], STATIONS[s]!.yaw, 0));
    // Each person holds the station they stand at or are walking to; a free one is chosen for the next task.
    const holds = starts.slice();
    const legs: { id: number; goal: XZ; end: number }[] = [];
    let waits = 0;
    const K = Math.round(240 / DT);
    for (let k = 1; k < K; k++) {
      for (const id of ids) {
        if (crowd.end(id) >= k) continue;
        const free = STATIONS.map((_, s) => s).filter((s) => !holds.includes(s));
        const pick = free[(k * 7 + id * 3) % free.length]!;
        const goal = STATIONS[pick]!;
        const r = plan(crowd, grid, id, k, goal, 3 + (id % 3));
        if (r.verdict.ok) {
          legs.push({ id, goal: goal.at, end: crowd.end(id) });
          holds[id] = pick;
        } else {
          // No way there yet without meeting someone: keep working a second longer and ask again.
          waits++;
          crowd.wait(id, k + Math.round(1 / DT));
        }
      }
    }
    // Work got done, by everyone.
    expect(legs.length).toBeGreaterThan(60);
    for (const id of ids) expect(legs.filter((l) => l.id === id).length).toBeGreaterThan(8);
    expect(waits).toBeLessThan(legs.length);
    // I5: each leg's last sample is its station.
    for (const l of legs) {
      const p = crowd.poseAt(l.id, l.end);
      expect(Math.hypot(p.x - l.goal[0], p.z - l.goal[1])).toBeLessThan(1e-4);
    }
    // I1 and I2 over the whole shift, at every sample.
    const audit = crowd.audit(0, K);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minStaticRoom).toBeGreaterThanOrEqual(STATIC_MARGIN - 1e-6);
  });

  test("a held door keeps everyone else out of its swing while it is open", () => {
    const { grid, crowd } = kitchen();
    const a = crowd.addAgent({ shapes: SHAPES }, 3.0, -1.0, Math.PI, 0);
    // A door swing across the whole aisle from 2 s to 6 s.
    crowd.addTimedRect({ minX: 2.4, minZ: -3.6, maxX: 3.6, maxZ: -3.0 }, Math.round(2 / DT), Math.round(6 / DT), 99);
    const r = plan(crowd, grid, a, 1, { at: [3.0, -6.4], yaw: Math.PI }, 1);
    expect(r.verdict.ok).toBe(true);
    // The walk crosses z = -3.3 only after the door shuts (or before it opens).
    for (let k = Math.round(2 / DT); k <= Math.round(6 / DT); k++) {
      const p = crowd.poseAt(a, k);
      expect(p.z > -3.0 + 0.2 || p.z < -3.6 - 0.2).toBe(true);
    }
  });
});
