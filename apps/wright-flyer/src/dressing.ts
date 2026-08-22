// Kitty Hawk scene-dressing mathematics (bead guzez.13). PURE and
// deterministic: gull flock kinematics, the Orville ground-helper
// timeline, the 1903 camp layout, and launch-rail geometry. The
// three.js builders in dressing3d.ts are thin consumers — every
// number that moves comes from here so it can be tested headless.
// Repro: node --test test/dressing.test.ts
//
// World frame (terrainMesh.ts): x east, z south, y up. The launch
// point is the origin for every layout offset. The flight runs +x.

/** Deterministic 32-bit LCG in [0, 1) (no Math.random — replays and
 * tests need the flock to be the same flock every run). */
export function lcg(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** One gull's flight-path parameters (all fixed at spawn). */
export interface GullPath {
  /** Orbit center offset from launch [m]. */
  readonly cx: number;
  readonly cz: number;
  /** Orbit radius [m] and height band [m]. */
  readonly radius: number;
  readonly height: number;
  /** Angular speed [rad/s] (sign = direction) and start phase. */
  readonly omega: number;
  readonly phase: number;
  /** Wing-flap frequency [Hz] and glide duty (fraction of time soaring). */
  readonly flapHz: number;
  readonly glide: number;
  /** Vertical bob amplitude [m]. */
  readonly bob: number;
}

/** Hard cap on flock size (typed refusal above; asserted at cap AND
 * cap+1 by the battery per workspace law). */
export const MAX_GULLS = 24;

/** Height band the flock occupies [m above the launch flat]. */
export const GULL_HEIGHT_MIN = 8;
export const GULL_HEIGHT_MAX = 55;

/** Spawn a deterministic flock: same (n, seed) -> same flock. */
export function gullFleet(n: number, seed: number): readonly GullPath[] {
  if (!Number.isInteger(n) || n < 1 || n > MAX_GULLS) {
    throw new RangeError(`flock size must be an integer in [1, ${MAX_GULLS}], got ${n}`);
  }
  const rand = lcg(seed);
  const out: GullPath[] = [];
  for (let i = 0; i < n; i += 1) {
    const radius = 18 + rand() * 60;
    out.push({
      cx: -40 + rand() * 220,
      cz: -90 + rand() * 180,
      radius,
      height: GULL_HEIGHT_MIN + rand() * (GULL_HEIGHT_MAX - GULL_HEIGHT_MIN),
      omega: (rand() < 0.5 ? -1 : 1) * (0.08 + rand() * 0.14),
      phase: rand() * Math.PI * 2,
      flapHz: 2.2 + rand() * 1.4,
      glide: 0.35 + rand() * 0.4,
      bob: 0.6 + rand() * 1.8,
    });
  }
  return out;
}

/** A gull's frame pose: position (launch-relative), heading, and wing
 * flap angle (0 = level; soaring intervals hold the wings up in a
 * shallow dihedral instead of beating). */
export function gullPose(
  g: GullPath,
  t: number,
): { x: number; y: number; z: number; headingRad: number; flapRad: number } {
  const a = g.phase + g.omega * t;
  const x = g.cx + g.radius * Math.cos(a);
  const z = g.cz + g.radius * Math.sin(a);
  const y = g.height + g.bob * Math.sin(0.31 * t + g.phase * 2);
  // Velocity direction on the circle (d/dt of position).
  const headingRad = Math.atan2(g.omega * Math.cos(a), -g.omega * Math.sin(a));
  // Beat-vs-soar: a slow square-ish cycle decides, the flap is a sine.
  const cycle = Math.sin(0.17 * t + g.phase * 3);
  const soaring = cycle < g.glide * 2 - 1;
  const flapRad = soaring
    ? 0.28
    : 0.65 * Math.sin(2 * Math.PI * g.flapHz * t + g.phase);
  return { x, y, z, headingRad, flapRad };
}

/** Orville's ground-helper timeline. He steadies the right wingtip and
 * runs alongside while the machine is ON THE RAIL (the famous photo
 * pose), lets go as it outruns him or lifts, slows to a stop, and
 * after a beat raises his field glasses to watch. Pure in its inputs:
 * no hidden state — the caller feeds the same numbers on replay and
 * gets the same Orville. */
export interface OrvillePose {
  /** Launch-relative position [m]. */
  readonly x: number;
  readonly z: number;
  /** Facing [rad] about y (0 = +x east, matching the flight). */
  readonly headingRad: number;
  /** Run-cycle phase [rad] (drives leg swing; 0 when standing). */
  readonly gaitRad: number;
  /** True once he has stopped and lifted the field glasses. */
  readonly glassesUp: boolean;
}

/** Top running speed [m/s] — he cannot outrun the machine for long. */
export const ORVILLE_MAX_MPS = 5.2;
/** Lateral offset: just beyond the right wingtip (z north = negative;
 * the right wingtip of an east-flying machine points SOUTH = +z). */
export const ORVILLE_SIDE_OFFSET_M = 7.4;
/** Seconds from letting go to raising the glasses. */
export const GLASSES_DELAY_S = 1.6;

/** How far along the rail Orville can be at time t while chasing a
 * machine at aircraftX: a hand's reach behind the wingtip fitting,
 * bounded by his legs from a standing start. The scene uses this to
 * LATCH his release point the first off-rail frame. */
export function orvilleReachableX(t: number, aircraftX: number): number {
  return Math.min(Math.max(aircraftX - 1.5, 0), ORVILLE_MAX_MPS * t);
}

export function orvillePose(
  t: number,
  onRail: boolean,
  aircraftX: number,
  releaseX: number | null,
  releaseT: number | null,
): OrvillePose {
  // Chasing: clamp his x to what his legs allow from a standing start.
  const reachable = orvilleReachableX(t, aircraftX);
  if (onRail && releaseX === null) {
    return {
      x: reachable,
      z: ORVILLE_SIDE_OFFSET_M,
      headingRad: 0,
      gaitRad: reachable > 0.2 ? 2.6 * Math.PI * t : 0,
      glassesUp: false,
    };
  }
  // Released: coast half a metre past the release point, stand, watch.
  const rx = releaseX ?? reachable;
  const rt = releaseT ?? t;
  const since = Math.max(0, t - rt);
  const coast = Math.min(0.5 * since, 0.6);
  return {
    x: rx + coast,
    z: ORVILLE_SIDE_OFFSET_M,
    headingRad: 0,
    gaitRad: 0,
    glassesUp: since >= GLASSES_DELAY_S,
  };
}

/** One camp prop placement (launch-relative, y from the terrain). */
export interface CampPlacement {
  readonly kind:
    | "hangar"
    | "shack"
    | "campfire"
    | "chair"
    | "barrel"
    | "workbench"
    | "toolchest";
  readonly x: number;
  readonly z: number;
  readonly rotY: number;
}

/** The rail corridor no prop may sit in: the machine travels +x from
 * the launch point with wingtips at ±6.4 m. */
export const RAIL_CLEAR_HALF_WIDTH_M = 8;

/** The 1903 camp: the two wooden buildings (hangar + living shack)
 * north-west of the rail, the campfire circle between them and the
 * rail, tools by the hangar door. Fixed layout — history does not
 * reroll per visit. */
export function campLayout(): readonly CampPlacement[] {
  return [
    { kind: "hangar", x: -26, z: -21, rotY: 0.18 },
    { kind: "shack", x: -40, z: -13, rotY: -0.12 },
    { kind: "campfire", x: -18, z: -11, rotY: 0 },
    { kind: "chair", x: -20.4, z: -9.2, rotY: 2.5 },
    { kind: "chair", x: -15.8, z: -9.4, rotY: -2.3 },
    { kind: "chair", x: -17.6, z: -13.6, rotY: 0.4 },
    { kind: "barrel", x: -29.5, z: -16.5, rotY: 0 },
    { kind: "barrel", x: -28.4, z: -15.2, rotY: 0.9 },
    { kind: "workbench", x: -22.5, z: -19.5, rotY: 0.2 },
    { kind: "toolchest", x: -24.5, z: -17.8, rotY: 1.7 },
  ];
}

/** Launch-rail geometry: the Wrights' 60 ft single-rail track — a
 * 2x4 on edge over half-buried ties. Tie pitch ~1.5 m, first tie AT
 * x=0 and last AT the rail end (both ends supported). */
export function railTies(railLengthM: number): readonly number[] {
  if (!(railLengthM > 0) || !Number.isFinite(railLengthM)) {
    throw new RangeError(`rail length must be positive and finite, got ${railLengthM}`);
  }
  const n = Math.max(2, Math.ceil(railLengthM / 1.5) + 1);
  const out: number[] = [];
  for (let i = 0; i < n; i += 1) {
    out.push((railLengthM * i) / (n - 1));
  }
  return out;
}
