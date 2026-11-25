/**
 * The crowd's invariants on a small kitchen: two counter runs either side of a 1.2 m aisle, an island, and people
 * working through stations on both sides. Run with `bun test src/three/models/crowd.test.ts`.
 */
import { describe, expect, test } from "bun:test";
import { Crowd, NavGrid, PAIR_MARGIN, STATIC_MARGIN, wrapPi, type Footprint, type LegCandidate, type Motion, type Rect, type TimedSpec, type XZ } from "./crowd";

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
/** A full hotel pan (0.53 m x 0.325 m) held at the waist in front of the belly: four discs over it, corners inside. */
const PAN: Footprint = [-1, 1].flatMap((sx) => [-1, 1].map((sz) => ({ x: sx * 0.1325, z: 0.33 + sz * 0.081, r: 0.155, high: true })));
/** A speed rack's 0.52 m x 0.66 m footprint, its handle side 0.49 m in front of the pusher: four discs over it. */
const RACK_FACE = 0.49;
const RACK: Footprint = [
  { x: -0.13, z: RACK_FACE + 0.165, r: 0.21 },
  { x: 0.13, z: RACK_FACE + 0.165, r: 0.21 },
  { x: -0.13, z: RACK_FACE + 0.495, r: 0.21 },
  { x: 0.13, z: RACK_FACE + 0.495, r: 0.21 },
];
/** The same pan held turned, to go onto a combi's runners lengthwise: its long side ahead, its near end at the belly. */
const TURNED: Footprint = [-1, 1].flatMap((sx) => [-1, 1].map((sz) => ({ x: sx * 0.081, z: 0.42 + sz * 0.15, r: 0.165, high: true })));
/** Standing, walking; standing and walking with a pan; behind a rack; with a pan held turned. */
const SHAPES = [STAND, WALK, [...STAND, ...PAN], [...WALK, ...PAN], [...STAND, ...RACK], [...WALK, ...RACK], [...STAND, ...TURNED], [...WALK, ...TURNED]];
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
function realize(start: { x: number; z: number; yaw: number }, path: readonly XZ[], speed: number, face: number, dwell: number, held = 0): Motion {
  return realizeLegs(start, [{ path, speed, face, dwell }], held);
}

/**
 * Several walks in one leg (in and out of a walk-in), each ending with a turn and a dwell. `held` 0: hands free;
 * 1: a pan carried (shapes 2 and 3) until the last dwell, where it is set down (shape 0); 2: a rack pushed (shapes 4
 * and 5) the whole way, let go for the last dwell; 3: a pan held turned (shapes 6 and 7), set down at the last.
 */
function realizeLegs(start: { x: number; z: number; yaw: number }, legs: { path: readonly XZ[]; speed: number; face: number; dwell: number }[], held = 0): Motion {
  const xs: number[] = [];
  const zs: number[] = [];
  const ys: number[] = [];
  const sh: number[] = [];
  let x = start.x;
  let z = start.z;
  let yaw = start.yaw;
  let free = false;
  const push = (moving: boolean) => {
    xs.push(x);
    zs.push(z);
    ys.push(yaw);
    sh.push(free || held === 0 ? (moving ? 1 : 0) : held === 1 ? (moving ? 3 : 2) : held === 2 ? (moving ? 5 : 4) : moving ? 7 : 6);
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
  for (const [li, { path, speed, face, dwell }] of legs.entries()) {
    for (let k = 1; k < path.length; k++) {
      const [ax, az] = path[k - 1]!;
      const [bx, bz] = path[k]!;
      const L = Math.hypot(bx - ax, bz - az);
      if (L < 1e-6) continue;
      // A rack goes where it is pushed: straight ahead, no turning on the spot behind it.
      if (held !== 2) turnTo(Math.atan2(bx - ax, bz - az));
      const n = Math.ceil(L / (speed * DT));
      for (let q = 1; q <= n; q++) {
        x = ax + ((bx - ax) * q) / n;
        z = az + ((bz - az) * q) / n;
        push(true);
      }
    }
    turnTo(face);
    if (li === legs.length - 1) free = true;
    for (let q = 0; q < Math.round(dwell / DT); q++) push(false);
  }
  if (!xs.length) push(false);
  return { n: xs.length, x: Float32Array.from(xs), z: Float32Array.from(zs), yaw: Float32Array.from(ys), shape: Uint8Array.from(sh) };
}

/** Whether any of a person's circles comes within `m` of a rectangle at step k. */
function nearRect(crowd: Crowd, id: number, k: number, r: Rect, m: number): boolean {
  const p = crowd.poseAt(id, k);
  for (const c of SHAPES[p.shape]!) {
    const cx = p.x + Math.cos(p.yaw) * c.x + Math.sin(p.yaw) * c.z;
    const cz = p.z - Math.sin(p.yaw) * c.x + Math.cos(p.yaw) * c.z;
    const dx = Math.max(r.minX - cx, 0, cx - r.maxX);
    const dz = Math.max(r.minZ - cz, 0, cz - r.maxZ);
    if (Math.hypot(dx, dz) < c.r + m) return true;
  }
  return false;
}

/** The counters and the table stand at 0.9 m: a carried pan passes over them (the room's tall grid is its walls). */
function kitchen(): { grid: NavGrid; crowd: Crowd } {
  const grid = new NavGrid(ROOM, COUNTERS, 0.04);
  return { grid, crowd: new Crowd(grid, DT, new NavGrid(ROOM, [], 0.04)) };
}

/** Plan a leg the way the hero does: path variants (keep right, hug the right, nearer the middle) at two speeds. */
function plan(crowd: Crowd, grid: NavGrid, id: number, k0: number, goal: (typeof STATIONS)[number], dwell: number, held = 0) {
  const here = crowd.poseAt(id, crowd.end(id));
  const dyn = crowd.standingObstacles(id, k0);
  const cands: LegCandidate[] = [];
  let variant = 0;
  // A pan held out in front needs its own room as well as the body's (held turned, less to the sides, but its reach
  // ahead swings wide where the way turns: kept further off as well).
  const rs = held === 1 ? [R_WALK + 0.06] : held === 3 ? [R_WALK + 0.03, R_WALK + 0.12] : [R_WALK];
  for (const r of rs)
    for (const keep of [0.75, 1, 0.4]) {
      const path = grid.findPath([here.x, here.z], goal.at, r, { dynamic: dyn, keepRight: keep });
      variant++;
      if (!path) continue;
      for (const v of [1.1, 0.8]) cands.push({ motion: realize(here, path, v, goal.yaw, dwell, held), speed: v, variant });
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

  test("every cell holds clearance() at its centre, capped at 4 m, whatever the obstacles", () => {
    // Random rooms with overlapping, touching and wall-hugging rectangles, some far apart (past the 4 m cap).
    let seed = 12345;
    const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
    for (let room = 0; room < 12; room++) {
      const W = 3 + rnd() * 14;
      const D = 3 + rnd() * 12;
      const bounds: Rect = { minX: -rnd() * 2, maxX: W, minZ: -D, maxZ: rnd() * 2 };
      const rects: Rect[] = [];
      for (let n = 0; n < 4 + Math.floor(rnd() * 20); n++) {
        const x = bounds.minX + rnd() * (bounds.maxX - bounds.minX);
        const z = bounds.minZ + rnd() * (bounds.maxZ - bounds.minZ);
        const w = rnd() < 0.2 ? 0 : rnd() * 2.5;
        const d = rnd() * 1.5;
        rects.push({ minX: x - w / 2, maxX: x + w / 2, minZ: z - d / 2, maxZ: z + d / 2 });
      }
      const cell = [0.04, 0.08, 0.13][room % 3]!;
      const grid = new NavGrid(bounds, rects, cell);
      for (let j = 0; j < grid.nz; j++)
        for (let i = 0; i < grid.nx; i++) {
          const want = Math.fround(Math.min(4, grid.clearance(grid.x0 + (i + 0.5) * cell, grid.z0 + (j + 0.5) * cell)));
          if (grid.clr[j * grid.nx + i] !== want) expect([room, i, j, grid.clr[j * grid.nx + i]]).toEqual([room, i, j, want]);
        }
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
    // They meet mid-aisle side by side: find the step where their z are equal and check they are apart in x, each on
    // its own right (a heads -z, its right is +x; b heads +z, its right is -x): a passes on the east side.
    let met = false;
    for (let k = 1; k < Math.max(crowd.end(a), crowd.end(b)); k++) {
      const pa = crowd.poseAt(a, k);
      const pb = crowd.poseAt(b, k);
      if (Math.abs(pa.z - pb.z) < 0.05 && pa.z < -1.6 && pa.z > -5.9) {
        met = true;
        expect(pa.x - pb.x).toBeGreaterThan(0.45);
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

describe("Crowd: doors, racks, pans and walk-ins", () => {
  test("a door opened as part of a leg waits for whoever is in its swing, and keeps everyone out while it is open", () => {
    const { grid, crowd } = kitchen();
    // A works at the west run; B is first committed: up the aisle past A to the front, through the door's swing.
    const a = crowd.addAgent({ shapes: SHAPES }, 2.7, -2.0, -Math.PI / 2, 0);
    const b = crowd.addAgent({ shapes: SHAPES }, 3.0, -4.6, 0, 0);
    expect(plan(crowd, grid, b, 1, { at: [3.0, -0.8], yaw: 0 }, 1).verdict.ok).toBe(true);
    // A swings a door 1.1 m out across the aisle for 3 s.
    const swing: Rect = { minX: 2.4, minZ: -3.7, maxX: 3.5, maxZ: -3.0 };
    const n = Math.round(3 / DT);
    const door: Motion = { n, x: new Float32Array(n).fill(2.7), z: new Float32Array(n).fill(-2.0), yaw: new Float32Array(n).fill(-Math.PI / 2), shape: new Uint8Array(n), timed: [{ rect: swing, k0: 0, k1: n - 1 }] };
    // Opening now would hit B on its way past: the door waits.
    expect(crowd.legClear(a, 1, 0, door)).toBe(false);
    const v = crowd.planLeg(a, 1, [{ motion: door, speed: 0, variant: 1 }], 20);
    expect(v.ok).toBe(true);
    expect(v.ok && v.delay).toBeGreaterThan(0);
    const k0 = 1 + Math.round((v.ok ? v.delay : 0) / DT);
    for (let k = k0; k < k0 + n; k++) expect(nearRect(crowd, b, k, swing, PAIR_MARGIN - 1e-6)).toBe(false);
    // C, planned last, crosses the swing only when the door is shut.
    const c = crowd.addAgent({ shapes: SHAPES }, 4.8, -0.8, -Math.PI / 2, 0);
    expect(plan(crowd, grid, c, k0 + 2, { at: [3.0, -6.4], yaw: Math.PI }, 1).verdict.ok).toBe(true);
    for (let k = k0; k < k0 + n; k++) expect(nearRect(crowd, c, k, swing, PAIR_MARGIN - 1e-6)).toBe(false);
    const audit = crowd.audit(0, Math.max(crowd.end(b), crowd.end(c)) + 10);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minRectGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
  });

  test("a rack held for good is checked in finite time against everyone's standing future", () => {
    const { crowd } = kitchen();
    const a = crowd.addAgent({ shapes: SHAPES }, 3.0, -1.0, 0, 0);
    const b = crowd.addAgent({ shapes: SHAPES }, 4.6, -0.8, 0, 0);
    const n = 5;
    const park = (rect: Rect): Motion => ({ n, x: new Float32Array(n).fill(3.0), z: new Float32Array(n).fill(-1.0), yaw: new Float32Array(n), shape: new Uint8Array(n), timed: [{ rect, k0: n - 1, k1: Infinity, owner: -1, key: "rack" }] });
    // Parked on B for ever: refused (and the check ends).
    expect(crowd.legClear(a, 1, 0, park({ minX: 4.3, minZ: -1.2, maxX: 4.9, maxZ: -0.6 }))).toBe(false);
    // Parked clear of both: taken, and from then on an obstacle to A as well as B.
    expect(crowd.legClear(a, 1, 0, park({ minX: 0.3, minZ: -1.2, maxX: 0.9, maxZ: -0.6 }))).toBe(true);
    crowd.commit(a, 1, 0, park({ minX: 0.3, minZ: -1.2, maxX: 0.9, maxZ: -0.6 }));
    expect(crowd.timedRect("rack")!.k1).toBe(Infinity);
    void b;
  });

  test("a pushed rack is part of its pusher's footprint, parked where it is left, and others keep clear of it", () => {
    const { grid, crowd } = kitchen();
    // The rack parked across the front area, its handle side facing +x; the pusher behind it, facing -x.
    // It stands at the mouth of the aisle; it goes 1.3 m west, out of the way.
    const home: Rect = { minX: 2.7, minZ: -1.2, maxX: 3.36, maxZ: -0.68 };
    crowd.addTimedRect(home, 0, Infinity, -1, "rack");
    const zc = (home.minZ + home.maxZ) / 2;
    const x0 = home.maxX + RACK_FACE;
    const a = crowd.addAgent({ shapes: SHAPES }, x0, zc, -Math.PI / 2, 0);
    const travel = 1.3;
    const x1 = x0 - travel;
    const m = realizeLegs({ x: x0, z: zc, yaw: -Math.PI / 2 }, [{ path: [[x0, zc], [x1, zc]], speed: 0.6, face: -Math.PI / 2, dwell: 1 }], 2);
    const dock: Rect = { minX: home.minX - travel, minZ: home.minZ, maxX: home.maxX - travel, maxZ: home.maxZ };
    const release = m.shape.findIndex((q) => q === 0);
    m.release = [{ key: "rack", k: 0 }];
    m.timed = [{ rect: dock, k0: release, k1: Infinity, owner: -1, key: "rack" } as TimedSpec];
    const kA = 20;
    crowd.wait(a, kA - 1);
    expect(crowd.planLeg(a, kA, [{ motion: m, speed: 0.6, variant: 1 }], 4).ok).toBe(true);
    // The rack's home is free from the moment the hands take it; its new place is held for good.
    expect(crowd.timedRect("rack")!.rect).toEqual(dock);
    // B, coming up the aisle from the moment A takes hold, wants the front corner past the rack's home.
    const b = crowd.addAgent({ shapes: SHAPES }, 3.0, -3.5, 0, 0);
    const vb = plan(crowd, grid, b, 1, { at: [4.8, -0.9], yaw: 0 }, 2);
    expect(vb.verdict.ok).toBe(true);
    // C may not stand where the rack is parked now: no way there; and the walk there on the floor plan alone (which
    // the rack is not on) is refused for the rack (I4), while the same walk stopping short of it is taken.
    const c = crowd.addAgent({ shapes: SHAPES }, 0.5, -0.9, Math.PI / 2, 0);
    expect(plan(crowd, grid, c, 1, { at: [1.75, -0.94], yaw: Math.PI / 2 }, 1).verdict.ok).toBe(false);
    const cStart = { x: 0.5, z: -0.9, yaw: Math.PI / 2 };
    const onDock = grid.findPath([0.5, -0.9], [1.75, -0.94], R_WALK);
    const short = grid.findPath([0.5, -0.9], [0.9, -0.94], R_WALK);
    expect(onDock && short).toBeTruthy();
    expect(crowd.legClear(c, 1, 0, realize(cStart, onDock!, 1, Math.PI / 2, 1))).toBe(false);
    expect(crowd.legClear(c, 1, 0, realize(cStart, short!, 1, Math.PI / 2, 1))).toBe(true);
    const K = Math.max(crowd.end(a), crowd.end(b)) + 20;
    // The rack's own circles move with A and keep clear of B (I2); B never touches the rack, parked or moving (I4).
    for (let k = 0; k <= K; k++) {
      if (k < kA) expect(nearRect(crowd, b, k, home, PAIR_MARGIN - 1e-6)).toBe(false);
      if (k >= kA + release) expect(nearRect(crowd, b, k, dock, PAIR_MARGIN - 1e-6)).toBe(false);
    }
    const audit = crowd.audit(0, K);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minStaticRoom).toBeGreaterThanOrEqual(STATIC_MARGIN - 1e-6);
    expect(audit.minRectGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
  });

  test("a pan carried through a working shift keeps its route clear of everyone, the pan included", () => {
    const { grid, crowd } = kitchen();
    // The cold store at the front wall, the line's pan seat at the far table; the carrier shuttles between them.
    const cold = { at: [1.0, -0.9] as XZ, yaw: 0 };
    const seat = { at: [3.0, -6.4] as XZ, yaw: Math.PI };
    const starts = [2, 3, 6];
    const ids = starts.map((s) => crowd.addAgent({ shapes: SHAPES }, STATIONS[s]!.at[0], STATIONS[s]!.at[1], STATIONS[s]!.yaw, 0));
    const carrier = crowd.addAgent({ shapes: SHAPES }, cold.at[0], cold.at[1], cold.yaw, 0);
    const holds = starts.slice();
    let carries = 0;
    let loaded = true;
    const K = Math.round(180 / DT);
    for (let k = 1; k < K; k++) {
      for (const id of [...ids, carrier]) {
        if (crowd.end(id) >= k) continue;
        if (id === carrier) {
          // To the seat with a pan, back to the store empty-handed.
          const r = plan(crowd, grid, id, k, loaded ? seat : cold, 3, loaded ? 1 : 0);
          if (r.verdict.ok) {
            if (loaded) carries++;
            loaded = !loaded;
          } else crowd.wait(id, k + Math.round(1 / DT));
          continue;
        }
        const i = ids.indexOf(id);
        const free = [0, 1, 2, 3].filter((s) => !holds.includes(s));
        const pick = free[(k * 7 + id * 3) % free.length]!;
        if (plan(crowd, grid, id, k, STATIONS[pick]!, 3 + (id % 3)).verdict.ok) holds[i] = pick;
        else crowd.wait(id, k + Math.round(1 / DT));
      }
    }
    expect(carries).toBeGreaterThan(4);
    const audit = crowd.audit(0, K);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minStaticRoom).toBeGreaterThanOrEqual(STATIC_MARGIN - 1e-6);
  }, 30000);

  test("a walk-in's doorway is one person's at a time, and closed to anyone without a window", () => {
    const { grid, crowd } = kitchen();
    // The walk-in: the floor west of the far table and south of the west run, entered from the aisle's far end.
    const inside: Rect = { minX: 0, minZ: -8, maxX: 1.9, maxZ: -6.05 };
    crowd.addGate("walkin", inside);
    const visit = (id: number, k0: number) => {
      const here = crowd.poseAt(id, crowd.end(id));
      const dyn = crowd.standingObstacles(id, k0, 10, ["walkin"]);
      const shelf: XZ = [0.9, -7.2];
      const p1 = grid.findPath([here.x, here.z], shelf, R_WALK, { dynamic: dyn, keepRight: 0 });
      const p2 = grid.findPath(shelf, [here.x, here.z], R_WALK, { dynamic: dyn, keepRight: 0 });
      expect(p1 && p2).toBeTruthy();
      const m = realizeLegs(here, [{ path: p1!, speed: 1, face: -Math.PI / 2, dwell: 3 }, { path: p2!, speed: 1, face: here.yaw, dwell: 1 }]);
      m.gates = [{ key: "walkin", k0: 0, k1: m.n - 1 }];
      return crowd.planLeg(id, k0, [{ motion: m, speed: 1, variant: 1 }], 30);
    };
    const a = crowd.addAgent({ shapes: SHAPES }, 4.6, -7.3, -Math.PI / 2, 0);
    const b = crowd.addAgent({ shapes: SHAPES }, 5.3, -6.5, -Math.PI / 2, 0);
    const va = visit(a, 1);
    const vb = visit(b, 1);
    expect(va.ok && vb.ok).toBe(true);
    // B goes in only once A is out.
    const wins = crowd.gateWindows("walkin");
    expect(wins.length).toBe(2);
    expect(wins[0]!.k1 < wins[1]!.k0 || wins[1]!.k1 < wins[0]!.k0).toBe(true);
    // Without a window there is no way in.
    const c = crowd.addAgent({ shapes: SHAPES }, 3.3, -3.2, Math.PI / 2, 0);
    expect(plan(crowd, grid, c, 1, { at: [0.9, -7.2], yaw: -Math.PI / 2 }, 1).verdict.ok).toBe(false);
    // The same visit on the floor plan alone, once A and B are out: refused with no window of its own (I6), taken
    // holding one.
    const kFree = Math.max(wins[0]!.k1, wins[1]!.k1) + 1;
    crowd.wait(c, kFree - 1);
    const shelf: XZ = [0.9, -7.2];
    const pIn = grid.findPath([3.3, -3.2], shelf, R_WALK, { keepRight: 0 });
    const pOut = grid.findPath(shelf, [3.3, -3.2], R_WALK, { keepRight: 0 });
    expect(pIn && pOut).toBeTruthy();
    const into = realizeLegs({ x: 3.3, z: -3.2, yaw: Math.PI / 2 }, [{ path: pIn!, speed: 1, face: -Math.PI / 2, dwell: 3 }, { path: pOut!, speed: 1, face: Math.PI / 2, dwell: 1 }]);
    expect(crowd.legClear(c, kFree, 0, into)).toBe(false);
    into.gates = [{ key: "walkin", k0: 0, k1: into.n - 1 }];
    expect(crowd.legClear(c, kFree, 0, into)).toBe(true);
    const audit = crowd.audit(0, Math.max(crowd.end(a), crowd.end(b)) + 10);
    expect(audit.gateIntrusions).toBe(0);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minStaticRoom).toBeGreaterThanOrEqual(STATIC_MARGIN - 1e-6);
  });

  test("a combi's door, opened and shut by hand, holds its swing from the pull to the latch: the loader waits for it to be clear, everyone else keeps out", () => {
    const { grid, crowd } = kitchen();
    // A combi in the east run, its door hinged at z = -4.45 and swinging out 0.65 m into the aisle; the loader opens
    // it, fetches the pan from the front, loads it, and shuts it, all one leg.
    const swing: Rect = { minX: 2.95, minZ: -4.45, maxX: 3.6, maxZ: -3.75 };
    // The loader beside the pull (the side away from the hinge); D at work in front of the oven, in the door's way,
    // for its first 3 s, then off to the front.
    const a = crowd.addAgent({ shapes: SHAPES }, 3.2, -3.45, Math.PI / 2, 0);
    const d = crowd.addAgent({ shapes: SHAPES }, 3.25, -4.1, Math.PI / 2, 0);
    crowd.wait(d, Math.round(3 / DT));
    expect(plan(crowd, grid, d, Math.round(3 / DT) + 1, { at: [4.6, -0.8], yaw: 0 }, 2).verdict.ok).toBe(true);
    const out = grid.findPath([3.2, -3.45], [3.0, -0.95], R_WALK, { keepRight: 0.75 })!;
    const back = grid.findPath([3.0, -0.95], [3.15, -3.5], R_WALK, { keepRight: 0.75 })!;
    expect(out && back).toBeTruthy();
    const m = realizeLegs({ x: 3.2, z: -3.45, yaw: Math.PI / 2 }, [
      { path: [[3.2, -3.45]], speed: 1, face: Math.PI / 2, dwell: 2.5 }, // the pull, the swing to its stop
      { path: out, speed: 1, face: 0, dwell: 3.5 }, // the pan off the counter
      { path: back, speed: 0.9, face: Math.PI / 2, dwell: 4.6 }, // slid in and pushed home
      { path: [[3.15, -3.5], [2.8, -3.4]], speed: 0.8, face: Math.PI, dwell: 2.5 }, // round to the pull, shut
    ]);
    const nOpen = m.n;
    m.timed = [{ rect: swing, k0: 0, k1: nOpen - 1 }];
    // At once it would swing into D: it waits until D has gone.
    expect(crowd.legClear(a, 1, 0, m)).toBe(false);
    const v = crowd.planLeg(a, 1, [{ motion: m, speed: 1, variant: 1 }], 12);
    expect(v.ok).toBe(true);
    const k0 = 1 + Math.round((v.ok ? v.delay : 0) / DT);
    for (let k = k0; k < k0 + nOpen; k++) expect(nearRect(crowd, d, k, swing, PAIR_MARGIN - 1e-6)).toBe(false);
    // C, from the front to the far table through the aisle the leaf narrows, keeps out of its swing while it is
    // open (trying again each second until it has a way).
    const c = crowd.addAgent({ shapes: SHAPES }, 1.0, -0.9, Math.PI, 0);
    let ok = false;
    for (let k = k0 + 2; k < k0 + nOpen + Math.round(30 / DT) && !ok; k += Math.round(1 / DT)) {
      if (crowd.end(c) >= k) continue;
      ok = plan(crowd, grid, c, k, { at: [3.0, -6.4], yaw: Math.PI }, 1).verdict.ok;
      if (!ok) crowd.wait(c, k + Math.round(1 / DT) - 1);
    }
    expect(ok).toBe(true);
    for (let k = k0; k < k0 + nOpen; k++) expect(nearRect(crowd, c, k, swing, PAIR_MARGIN - 1e-6)).toBe(false);
    const K = Math.max(crowd.end(a), crowd.end(c), crowd.end(d)) + 10;
    const audit = crowd.audit(0, K);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minStaticRoom).toBeGreaterThanOrEqual(STATIC_MARGIN - 1e-6);
    expect(audit.minRectGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
  });

  test("a pan held turned reaches further ahead: it keeps clear of everyone, and off a tall oven it walks past", () => {
    // The same kitchen with a tall combi in the east run beside the aisle (the counters pass under a pan; it does not).
    const oven: Rect = { minX: 3.6, minZ: -3.0, maxX: 4.45, maxZ: -2.1 };
    const grid = new NavGrid(ROOM, COUNTERS, 0.04);
    const crowd = new Crowd(grid, DT, new NavGrid(ROOM, [oven], 0.04));
    const cold = { at: [1.0, -0.9] as XZ, yaw: 0 };
    const seat = { at: [3.0, -6.4] as XZ, yaw: Math.PI };
    const starts = [3, 6];
    const ids = starts.map((s) => crowd.addAgent({ shapes: SHAPES }, STATIONS[s]!.at[0], STATIONS[s]!.at[1], STATIONS[s]!.yaw, 0));
    const carrier = crowd.addAgent({ shapes: SHAPES }, cold.at[0], cold.at[1], cold.yaw, 0);
    const holds = starts.slice();
    let carries = 0;
    let loaded = true;
    const K = Math.round(150 / DT);
    for (let k = 1; k < K; k++) {
      for (const id of [...ids, carrier]) {
        if (crowd.end(id) >= k) continue;
        if (id === carrier) {
          const r = plan(crowd, grid, id, k, loaded ? seat : cold, 3, loaded ? 3 : 0);
          if (r.verdict.ok) {
            if (loaded) carries++;
            loaded = !loaded;
          } else crowd.wait(id, k + Math.round(1 / DT));
          continue;
        }
        const i = ids.indexOf(id);
        const free = [0, 1, 2, 3].filter((s) => !holds.includes(s));
        if (!free.length) {
          crowd.wait(id, k + Math.round(2 / DT));
          continue;
        }
        const pick = free[(k * 7 + id * 3) % free.length]!;
        if (plan(crowd, grid, id, k, STATIONS[pick]!, 3 + (id % 3)).verdict.ok) holds[i] = pick;
        else crowd.wait(id, k + Math.round(1 / DT));
      }
    }
    expect(carries).toBeGreaterThan(3);
    // Turned, the pan's far end was a body's depth further ahead than crosswise: it was carried that way.
    let turned = 0;
    for (let k = 0; k < K; k++) if (crowd.poseAt(carrier, k).shape === 7) turned++;
    expect(turned).toBeGreaterThan(Math.round(10 / DT));
    const audit = crowd.audit(0, K);
    expect(audit.minPairGap).toBeGreaterThanOrEqual(PAIR_MARGIN - 1e-6);
    expect(audit.minStaticRoom).toBeGreaterThanOrEqual(STATIC_MARGIN - 1e-6);
  }, 30000);
});
