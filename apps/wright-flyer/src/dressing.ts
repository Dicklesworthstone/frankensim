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
  /** When set (machine down, ground contact), he RUNS to the machine
   * after watching — starting GLASSES_DELAY_S after release, clamped
   * to his top speed and a hand's reach short of the wingtip. */
  landedX: number | null = null,
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
  // Released: coast half a metre past the release point, stand, watch…
  const rx = releaseX ?? reachable;
  const rt = releaseT ?? t;
  const since = Math.max(0, t - rt);
  const coast = Math.min(0.5 * since, 0.6);
  if (landedX !== null && since >= GLASSES_DELAY_S) {
    // …then sprint to the machine (the famous run to the landed Flyer).
    const target = landedX - 2.2;
    const runStart = rt + GLASSES_DELAY_S;
    const runT = Math.max(0, t - runStart);
    const x = Math.min(rx + ORVILLE_MAX_MPS * runT, target);
    const moving = x < target - 0.05;
    return {
      x,
      z: ORVILLE_SIDE_OFFSET_M,
      headingRad: moving ? 0 : Math.atan2(ORVILLE_SIDE_OFFSET_M, Math.max(target - x, 0.5)),
      gaitRad: moving ? 2.9 * Math.PI * t : 0,
      glassesUp: false,
    };
  }
  return {
    x: rx + coast,
    z: ORVILLE_SIDE_OFFSET_M,
    headingRad: 0,
    gaitRad: 0,
    glassesUp: since >= GLASSES_DELAY_S,
  };
}

/** Landing-dust burst (T3.5): one grain thrown outward from the
 * touchdown point in MACHINE-LOCAL offsets. A 1.8 s closed burst once
 * `t` passes `t0`; null before it starts. */
export function landingDust(
  i: number,
  t0: number,
  t: number,
): { dx: number; dy: number; dz: number; scale: number; opacity: number } | null {
  const durS = 1.8;
  const u = (t - t0) / durS;
  if (u < 0 || u > 1) {
    return null;
  }
  const ang = hash01(i * 89 + 1) * Math.PI * 2;
  const reach = (1.5 + hash01(i * 97 + 3) * 4.5) * Math.sqrt(u);
  return {
    dx: Math.cos(ang) * reach,
    dy: 0.15 + 0.9 * u * (1 - u) * (1 + hash01(i * 101)),
    dz: Math.sin(ang) * reach * 0.8,
    scale: 0.7 + u * 2.4,
    opacity: (1 - u) * 0.42,
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
    | "toolchest"
    | "crate"
    | "battery"
    | "oilcan"
    | "windpost";
  readonly x: number;
  readonly z: number;
  readonly rotY: number;
}

/** The rail corridor no prop may sit in: the machine travels +x from
 * the launch point with wingtips at ±6.4 m. */
export const RAIL_CLEAR_HALF_WIDTH_M = 8;

/** The 1903 camp: the two wooden buildings (hangar + living shack)
 * north-west of the rail, the campfire circle between them and the
 * rail, tools by the hangar door — plus the working clutter (T2.4):
 * packing crates, the magneto/battery box for the engine, an oil can,
 * and the hand-anemometer post. Fixed layout — history does not
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
    { kind: "crate", x: -21.8, z: -16.9, rotY: 0.3 },
    { kind: "crate", x: -21.1, z: -16.2, rotY: 0.85 },
    { kind: "crate", x: -22.4, z: -16.3, rotY: 0.05 },
    { kind: "battery", x: -23.9, z: -15.5, rotY: -0.4 },
    { kind: "oilcan", x: -22.2, z: -18.8, rotY: 0 },
    { kind: "windpost", x: -11.5, z: -7.0, rotY: 0 },
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

/* ---------- presentation-flight math (gulls v2, particles) ----------
 * Everything here stays PURE and stateless-in-(inputs, t): the scene
 * feeds wall-clock sim time and gets poses; replays and tests see the
 * exact same sky. No THREE.js types leak into this module. */

/** Single-integer avalanche hash in [0, 1) — the per-index constant
 * source for stateless particles (no allocator, no shared RNG state). */
export function hash01(n: number): number {
  let h = Math.imul(n | 0, 2654435761);
  h ^= h >>> 15;
  h = Math.imul(h, 2246822519);
  h ^= h >>> 13;
  return (h >>> 0) / 4294967296;
}

/** A gull's bank/pitch for frame orientation: bank comes from the
 * orbit's centripetal requirement (v·ω / g, capped ~57°), pitch from
 * the climb rate of the bob term (nose down on descent, up on climb,
 * capped ±20°). Roll sign follows ω's sign — verify visually against
 * the turn direction in the scene battery screenshot. */
export function gullAttitude(
  g: GullPath,
  t: number,
): { rollRad: number; pitchRad: number } {
  const speed = Math.abs(g.omega) * g.radius;
  const bank = Math.atan((speed * Math.abs(g.omega)) / 9.81);
  const rollRad = Math.min(1.0, bank) * Math.sign(g.omega);
  const climbMps = g.bob * 0.31 * Math.cos(0.31 * t + g.phase * 2);
  const pitchRad = Math.max(
    -0.35,
    Math.min(0.35, -Math.atan2(climbMps, Math.max(speed, 0.5))),
  );
  return { rollRad, pitchRad };
}

/** Hard cap on low flyby count (asserted at cap AND cap+1). */
export const MAX_FLYBYS = 6;

/** One low pass: a STRAIGHT track through near-launch airspace, looped
 * on a per-bird period so passes recur all day without accumulating
 * state. `dir` is a unit travel vector; the bird lives while its
 * along-track coordinate sits inside ±halfSpanM. */
export interface FlybyPath {
  /** Closest-approach point, launch-relative [m]. */
  readonly cx: number;
  readonly cz: number;
  /** Track altitude above the launch flat [m]. */
  readonly y: number;
  /** Unit travel direction (heading = atan2(-dirZ, dirX) convention
   * matches gullPose: Ry(-heading) faces travel). */
  readonly dirX: number;
  readonly dirZ: number;
  readonly speedMps: number;
  /** Recurrence period [s] and phase offset along the track [m]. */
  readonly periodS: number;
  readonly offsetM: number;
  /** Active half-length of the track [m] (spawn/despawn past camera). */
  readonly halfSpanM: number;
  readonly flapHz: number;
  readonly glide: number;
  readonly phase: number;
}

/** Spawn a deterministic low-flyby fleet: same (n, seed) -> same birds.
 * Tracks thread NEAR the diorama (5–16 m up, passing within ~30 m of
 * the rail) so the flock reads as living around the player, not as
 * distant specks. */
export function flybyFleet(n: number, seed: number): readonly FlybyPath[] {
  if (!Number.isInteger(n) || n < 1 || n > MAX_FLYBYS) {
    throw new RangeError(`flyby count must be an integer in [1, ${MAX_FLYBYS}], got ${n}`);
  }
  const rand = lcg(seed ^ 0x5f3759df);
  const out: FlybyPath[] = [];
  for (let i = 0; i < n; i += 1) {
    const heading = (rand() - 0.5) * 1.2 + (rand() < 0.5 ? 0 : Math.PI);
    out.push({
      cx: -30 + rand() * 60,
      cz: -24 + rand() * 48,
      y: 5 + rand() * 11,
      dirX: Math.cos(heading),
      dirZ: -Math.sin(heading),
      speedMps: 8.5 + rand() * 4.5,
      periodS: 26 + rand() * 30,
      offsetM: rand() * 400,
      halfSpanM: 150 + rand() * 60,
      flapHz: 2.4 + rand() * 1.2,
      glide: 0.25 + rand() * 0.35,
      phase: rand() * Math.PI * 2,
    });
  }
  return out;
}

/** A flyby's pose at time t, or NULL while the bird is outside its
 * active span (the scene hides it — no off-camera mesh churn). */
export function flybyPose(
  p: FlybyPath,
  t: number,
): { x: number; y: number; z: number; headingRad: number; flapRad: number } | null {
  const track = (((p.speedMps * t + p.offsetM) % (2 * p.halfSpanM)) + 2 * p.halfSpanM) %
    (2 * p.halfSpanM);
  const s = track - p.halfSpanM;
  if (Math.abs(s) > p.halfSpanM - 8) {
    return null;
  }
  // Gentle shallow descent across the pass — arriving high, leaving low
  // reads as a bird hunting the flat rather than riding rails.
  const y = p.y - (s / p.halfSpanM) * 2.2;
  const cycle = Math.sin(0.23 * t + p.phase * 3);
  const soaring = cycle < p.glide * 2 - 1;
  const flapRad = soaring ? 0.3 : 0.7 * Math.sin(2 * Math.PI * p.flapHz * t + p.phase);
  return {
    x: p.cx + p.dirX * s,
    y,
    z: p.cz + p.dirZ * s,
    headingRad: Math.atan2(-p.dirZ, p.dirX),
    flapRad,
  };
}

/* ------------------------- stateless particles ---------------------- */

const FRACT = (v: number): number => v - Math.floor(v);

/** Propwash sand: one grain-plume particle behind a pulling machine.
 * Age wraps per-particle (closed loop in t), the plume streams BACK
 * along -x (props push air aft over the rear-mounted engine), spreads
 * laterally, arcs up then settles. `strength` in [0,1] fades the whole
 * plume (0 = engine quiet -> invisible). Launch-relative metres. */
export function propwashPuff(
  i: number,
  t: number,
  machineX: number,
  strength: number,
): { x: number; y: number; z: number; scale: number; opacity: number } {
  const lifeS = 0.9 + hash01(i * 3 + 11) * 0.9;
  const u = FRACT(t / lifeS + hash01(i * 7 + 5));
  const back = 2.5 + u * (9 + hash01(i * 13 + 1) * 8);
  const side = (hash01(i * 17 + 3) - 0.5) * (1.5 + u * 7);
  const lift = 0.2 + 2.4 * u * (1 - u) * (0.5 + hash01(i * 19 + 7));
  const fade = u < 0.08 ? u / 0.08 : 1 - (u - 0.08) / 0.92;
  return {
    x: machineX - back,
    y: lift,
    z: side,
    scale: 0.6 + u * 3.2,
    opacity: Math.max(0, Math.min(1, strength)) * Math.max(0, fade) * 0.5,
  };
}

/** Engine exhaust: one smokelet from the oil-smeared crank area, rising
 * and thinning. Fire-and-forget closed loop, launch-relative offsets
 * from the machine's tail (caller adds machineX). */
export function exhaustPuff(
  i: number,
  t: number,
  rpm01: number,
): { dx: number; dy: number; dz: number; scale: number; opacity: number } {
  const lifeS = 1.4 + hash01(i * 23 + 2) * 1.2;
  const u = FRACT(t / lifeS + hash01(i * 29 + 9));
  const rise = 0.35 + u * (1.6 + 0.8 * rpm01);
  const wander = 0.22 * Math.sin(u * 9 + hash01(i * 31 + 4) * 6.28);
  return {
    dx: -2.9 - u * 1.6,
    dy: rise,
    dz: 0.18 + wander,
    scale: 0.25 + u * 1.5,
    opacity: Math.max(0, Math.min(1, rpm01 * 1.4)) * (1 - u) * 0.42,
  };
}

/** Campfire ember: one spark rising off the flames in FIRE-LOCAL
 * coordinates (caller places at the fire ring). Loops forever. */
export function emberAt(
  i: number,
  t: number,
): { x: number; y: number; z: number; opacity: number } {
  const lifeS = 1.6 + hash01(i * 37 + 6) * 1.6;
  const u = FRACT(t / lifeS + hash01(i * 41 + 8));
  const swirlA = hash01(i * 43 + 10) * Math.PI * 2 + u * 5;
  const r = 0.14 * (1 - u);
  return {
    x: Math.cos(swirlA) * r,
    y: 0.45 + u * 2.3,
    z: Math.sin(swirlA) * r,
    opacity: (1 - u) * (0.55 + 0.45 * hash01(i)),
  };
}

/* -------------------- wind-made-visible math ----------------------- */
/* The December 17 headwind blows FROM the east: the machine flies +x
 * INTO it, so the WIND VELOCITY points -x. Every streamer, flag, and
 * smoke column below advects toward -x, scaled by the scenario's
 * headwind. All closed loops in t — no state, replays identical. */

/** Headwind the scene renders when no scenario value is supplied. */
export const DEFAULT_HEADWIND_MPS = 11;

/** One scrub-vegetation placement (launch-relative; y from terrain). */
export interface ScrubPlacement {
  readonly x: number;
  readonly z: number;
  readonly rotY: number;
  readonly scale: number;
  readonly kind: "tuft" | "bush" | "pine";
}

/** The launch/camp flat stays bare (same law as duneDetail's mask). */
export const SCRUB_FLAT_RADIUS_M = 62;
/** Half-extent of the scrub field around the diorama. */
export const SCRUB_SPAN_M = 470;

/** Deterministic scrub field: ring-rejection sampling outside the
 * launch flat, the rail corridor, and the camp clearing. Same
 * (counts, seed) -> same field, every run. */
export function scrubField(
  counts: { tufts: number; bushes: number; pines: number },
  seed: number,
): readonly ScrubPlacement[] {
  const total = counts.tufts + counts.bushes + counts.pines;
  if (!Number.isInteger(total) || total < 0 || total > 1200) {
    throw new RangeError(`scrub total out of [0, 1200]: ${total}`);
  }
  const rand = lcg(seed ^ 0x9e3779b9);
  const out: ScrubPlacement[] = [];
  const inRailCorridor = (x: number, z: number): boolean =>
    Math.abs(z) < 10 && x > -8 && x < 42;
  const inCampClearing = (x: number, z: number): boolean =>
    x > -48 && x < -8 && z > -27 && z < -3;
  let guard = 0;
  while (out.length < total && guard < total * 40) {
    guard += 1;
    const x = (rand() * 2 - 1) * SCRUB_SPAN_M;
    const z = (rand() * 2 - 1) * SCRUB_SPAN_M;
    if (Math.hypot(x, z) < SCRUB_FLAT_RADIUS_M) {
      continue;
    }
    if (inRailCorridor(x, z) || inCampClearing(x, z)) {
      continue;
    }
    const idx = out.length;
    const kind: ScrubPlacement["kind"] =
      idx < counts.pines ? "pine" : idx < counts.pines + counts.bushes ? "bush" : "tuft";
    out.push({
      x,
      z,
      rotY: rand() * Math.PI * 2,
      scale: kind === "pine" ? 0.8 + rand() * 0.9 : 0.6 + rand() * 0.8,
      kind,
    });
  }
  if (out.length !== total) {
    throw new RangeError(`scrub rejection sampling starved: ${out.length}/${total}`);
  }
  return out;
}

/** Sand streamer ribbon: `segs` points of one low ribbon of sand
 * skittering downwind (-x) from a hash-anchored start. `i` picks the
 * streamer, `seg` the point along it. Amplitude grows with wind. */
export function streamerPoint(
  i: number,
  seg: number,
  segs: number,
  t: number,
  windMps: number,
): { x: number; y: number; z: number } {
  const anchorX = -50 + hash01(i * 53 + 3) * 150;
  const anchorZ = -70 + hash01(i * 59 + 7) * 140;
  const lenM = 3.5 + hash01(i * 61 + 11) * 5.5;
  const s = seg / segs;
  const phase = t * (0.9 + hash01(i * 67 + 13) * 0.7);
  const flutter = Math.sin(phase * 6.1 + s * 9.4 + hash01(i) * 6.28);
  const lift = (0.06 + 0.3 * s) * (0.5 + 0.5 * Math.sin(phase * 3.3 + s * 5.1));
  return {
    x: anchorX - s * lenM * (0.6 + windMps / 14),
    y: 0.05 + lift * (0.4 + windMps / 40) + 0.04 * flutter,
    z: anchorZ + s * lenM * 0.16 * flutter,
  };
}

/** Flag cloth point: `seg`/`segs` along a vertical-hinged strip that
 * streams downwind (-x) from the pole top, with a traveling wave whose
 * speed scales with the wind. FIRE-LOCAL to the pole (caller adds the
 * pole position); y measured DOWN from the truck. */
export function flagPoint(
  seg: number,
  segs: number,
  t: number,
  windMps: number,
): { x: number; y: number; z: number } {
  const s = seg / segs;
  const wave = Math.sin(t * (2.2 + windMps * 0.35) - s * 7.5);
  const droop = (1 - Math.min(1, windMps / 8)) * s * s * 1.1;
  return {
    x: -s * (0.9 + windMps * 0.05),
    y: -droop + 0.06 * wave * s,
    z: 0.10 * s * wave,
  };
}

/** Campfire smoke: one smokelet in FIRE-LOCAL coordinates rising from
 * the flames and bending downwind (-x) with height. Loops forever;
 * `strength` in [0,1] scales the whole column (0 = no fire). */
export function smokePuff(
  i: number,
  t: number,
  windMps: number,
  strength: number,
): { x: number; y: number; z: number; scale: number; opacity: number } {
  const lifeS = 5.5 + hash01(i * 71 + 5) * 3.5;
  const u = FRACT(t / lifeS + hash01(i * 73 + 9));
  const rise = 0.5 + u * (5.5 + hash01(i * 79 + 2) * 2.5);
  const bend = (rise * rise) / 40 * (0.5 + windMps / 12);
  const sway = 0.3 * Math.sin(u * 7 + hash01(i * 83 + 4) * 6.28) * u;
  const fade = u < 0.12 ? u / 0.12 : 1 - (u - 0.12) / 0.88;
  return {
    x: -bend + sway,
    y: rise,
    z: 0.18 * Math.sin(u * 5.3 + hash01(i) * 6.28) * u,
    scale: 0.35 + u * 2.6,
    opacity: Math.max(0, Math.min(1, strength)) * Math.max(0, fade) * 0.34,
  };
}
