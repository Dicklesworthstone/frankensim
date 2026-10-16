// Wright-brother figure mathematics (bead guzez.14). PURE: per-brother
// anthropometry, the running-gait joint model, arm aiming, and the
// prone-pilot/binocular pose constants. figure3d.ts is the thin THREE
// consumer — every joint angle comes from here so the motion is
// headless-testable and identical on replay.
// Repro: node --test test/figure.test.ts
//
// Segment proportions follow the Drillis & Contini (1966) fractions of
// stature — the standard anthropometric table — so the bodies read as
// PEOPLE, not scaled boxes. Statures/builds from the historical record:
// Wilbur ~5'10" (1.78 m) and wiry; Orville ~5'8" (1.73 m) and stockier.

export type Brother = "wilbur" | "orville";

/** Segment lengths [m] + build parameters for one brother. */
export interface FigureSpec {
  readonly heightM: number;
  /** Girth multiplier on limb/torso radii (1 = reference build). */
  readonly build: number;
  readonly headRadiusM: number;
  readonly shoulderWidthM: number;
  readonly hipWidthM: number;
  readonly torsoLenM: number;
  readonly upperArmLenM: number;
  readonly forearmLenM: number;
  readonly thighLenM: number;
  readonly shinLenM: number;
  readonly footLenM: number;
  /** Hip-joint height when standing [m] (thigh+shin+ankle). */
  readonly hipHeightM: number;
  /** Shoulder-joint height when standing [m]. */
  readonly shoulderHeightM: number;
}

/** Drillis-Contini stature fractions (their Fig. 2). */
const F = {
  head: 0.065, // head RADIUS ≈ half of the 0.13H head height
  shoulderW: 0.259,
  hipW: 0.191,
  upperArm: 0.186,
  forearm: 0.146,
  thigh: 0.245,
  shin: 0.246,
  foot: 0.152,
  hipH: 0.53,
  shoulderH: 0.818,
} as const;

export function figureSpec(brother: Brother): FigureSpec {
  const heightM = brother === "wilbur" ? 1.78 : 1.73;
  const build = brother === "wilbur" ? 0.92 : 1.06;
  const hipHeightM = F.hipH * heightM;
  const shoulderHeightM = F.shoulderH * heightM;
  return {
    heightM,
    build,
    headRadiusM: F.head * heightM,
    shoulderWidthM: F.shoulderW * heightM * (0.9 + 0.1 * build),
    hipWidthM: F.hipW * heightM * (0.9 + 0.1 * build),
    torsoLenM: shoulderHeightM - hipHeightM,
    upperArmLenM: F.upperArm * heightM,
    forearmLenM: F.forearm * heightM,
    thighLenM: F.thigh * heightM,
    shinLenM: F.shin * heightM,
    footLenM: F.foot * heightM,
    hipHeightM,
    shoulderHeightM,
  };
}

/** One frame of running/standing joint angles [rad]. Convention:
 * positive hip/shoulder = segment swings FORWARD (+x, the flight
 * direction); knee/elbow are pure flexion, always >= 0. */
export interface GaitPose {
  readonly hipL: number;
  readonly hipR: number;
  readonly kneeL: number;
  readonly kneeR: number;
  readonly shoulderL: number;
  readonly shoulderR: number;
  readonly elbowL: number;
  readonly elbowR: number;
  /** Forward torso lean [rad]. */
  readonly leanRad: number;
  /** Vertical bob above standing height [m] (two per stride). */
  readonly bobM: number;
}

/** Reference top speed the gait saturates at [m/s] (matches Orville's
 * dressing.ts chase cap). */
export const GAIT_MAX_MPS = 5.2;

/** Swing-phase knee flexion cap [rad] (~80 deg at a full run). */
export const KNEE_FLEX_MAX_RAD = 1.4;

/** Joint angles at stride phase `phaseRad` and ground speed. Speed 0
 * is EXACTLY the standing pose (all zeros except a relaxed elbow). */
export function gaitPose(phaseRad: number, speedMps: number): GaitPose {
  if (!Number.isFinite(phaseRad) || !Number.isFinite(speedMps)) {
    throw new RangeError(`gait inputs must be finite, got ${phaseRad}, ${speedMps}`);
  }
  const s = Math.min(1, Math.max(0, speedMps / GAIT_MAX_MPS));
  if (s === 0) {
    return {
      hipL: 0,
      hipR: 0,
      kneeL: 0,
      kneeR: 0,
      shoulderL: 0,
      shoulderR: 0,
      elbowL: 0.1,
      elbowR: 0.1,
      leanRad: 0,
      bobM: 0,
    };
  }
  const hipAmp = 0.62 * s;
  const hipL = hipAmp * Math.sin(phaseRad);
  const hipR = hipAmp * Math.sin(phaseRad + Math.PI);
  // Knee flexes during the leg's SWING (hip moving forward): gate on
  // the swing half-cycle, squared for a soft touchdown.
  const swing = (p: number): number => {
    const g = Math.max(0, Math.sin(p + 0.5));
    return KNEE_FLEX_MAX_RAD * s * g * g;
  };
  const kneeL = swing(phaseRad);
  const kneeR = swing(phaseRad + Math.PI);
  // Arms counter-swing the same-side leg, elbows carried bent.
  const shoulderL = -0.75 * hipL;
  const shoulderR = -0.75 * hipR;
  const elbow = 0.1 + 1.15 * s;
  return {
    hipL,
    hipR,
    kneeL,
    kneeR,
    shoulderL,
    shoulderR,
    elbowL: elbow,
    elbowR: elbow,
    leanRad: 0.22 * s,
    bobM: 0.035 * s * (0.5 + 0.5 * Math.cos(2 * phaseRad)),
  };
}

/** Stride frequency [Hz] from speed and leg length (hip height):
 * stride length ~ 0.8 + 0.45·s legs at speed, so f = v / L. Monotone
 * in speed for fixed legs; 0 at rest. */
export function strideFreqHz(speedMps: number, legLenM: number): number {
  if (!(legLenM > 0)) {
    throw new RangeError(`leg length must be positive, got ${legLenM}`);
  }
  const v = Math.max(0, speedMps);
  if (v === 0) {
    return 0;
  }
  const s = Math.min(1, v / GAIT_MAX_MPS);
  const strideLenM = legLenM * (0.8 + 0.9 * s);
  return v / strideLenM;
}

/** Aim angles for an arm whose rest pose hangs along -y: yaw about y
 * then pitch about z' to point the arm axis at a target in
 * SHOULDER-LOCAL coordinates (x forward, y up, z right). Returns the
 * rotations figure3d applies; pure trigonometry. */
export function armAimAngles(
  dx: number,
  dy: number,
  dz: number,
): { yawRad: number; pitchRad: number } {
  const len = Math.hypot(dx, dy, dz);
  if (!(len > 1e-9)) {
    throw new RangeError("arm aim target coincides with the shoulder");
  }
  // pitch from straight-down (-y): 0 = hanging, pi/2 = horizontal.
  const horiz = Math.hypot(dx, dz);
  const pitchRad = Math.atan2(horiz, -dy);
  // YZX application (figure3d aimLeftArm) composes R = RY(yaw)*RZ(pitch);
  // solving R*(sin p, -cos p, ...) = (dx, dy, dz)/|d| yields
  // yaw = atan2(-dz, dx) — the plain atan2(dz, dx) mirrors the arm
  // across the sagittal plane (102 deg off for a left-side wingtip grip).
  const yawRad = Math.atan2(-dz, dx);
  return { yawRad, pitchRad };
}

/** Wilbur's prone piloting pose (hips on the cradle, head up watching
 * the canard, arms forward to the lever). Angles in the same joint
 * convention as GaitPose. */
export const PRONE_POSE = {
  backArchRad: 0.18,
  // Rz(+θ) tips the +x-facing face toward +y — POSITIVE is head-UP
  // (the -0.55 originally shipped here stared into the wing fabric).
  headPitchRad: 0.55,
  hipFlexRad: 0.08,
  kneeFlexRad: 0.12,
  shoulderForwardRad: 1.35,
  elbowFlexRad: 0.5,
} as const;

/** Orville's field-glasses pose: both hands to the eyes. */
export const BINOCULAR_POSE = {
  shoulderForwardRad: 2.35,
  elbowFlexRad: 1.15,
} as const;
