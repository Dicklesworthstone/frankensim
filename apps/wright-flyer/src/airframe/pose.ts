// E2.2 rig pose core (bead wf-root-guzez.3.2): the PURE mapping from
// physics/control state to the articulated-airframe pose. Headless-testable
// — no three.js here; applyPose (the thin mutator) consumes this.
//
// Laws encoded from the frozen evidence: canard travel ±30° is
// photo-inferred (absent-by-verification) — commands clamp there and the
// clamp is REPORTED; warp limit ±8.5° (verified); the 1903 slaved rudder
// δr = 2.5·δw when coupled (verified); geometry deviating beyond ±25% of
// the reference design flags `schematicPreview` (plan DONE-WHEN) — the
// visual stays honest about being an extrapolated schematic.

export const CANARD_TRAVEL_DEG = 30; // photo-inferred (flyer-reference)
export const WARP_LIMIT_DEG = 8.5; // AIAA 2004-0211 (verified)
export const RUDDER_SLAVING = 2.5; // flown configuration (verified)
export const RUDDER_LIMIT_DEG = WARP_LIMIT_DEG * RUDDER_SLAVING;
export const SCHEMATIC_PREVIEW_FRACTION = 0.25;

/** Reference dims the preview flag measures against (flyer-reference). */
export const REFERENCE_DIMS = {
  span_m: 12.29,
  chord_m: 1.981,
  canard_area_m2: 4.46,
  rudder_area_m2: 1.86,
} as const;

export interface ControlState {
  /** Canard command [deg], positive nose-up. */
  canardDeg: number;
  /** Warp command [deg], positive right-wing-down. */
  warpDeg: number;
  /** Independent rudder command [deg]; ignored when `coupled`. */
  rudderDeg: number;
  /** 1903 slaved-rudder wiring active. */
  coupled: boolean;
  /** Propeller shaft angle [rad] (from the physics rotor state). */
  propAngleRad: number;
}

export interface FlyerPose {
  canardRad: number;
  rudderRad: number;
  /** Tip-section warp rotation [rad] (the 0.6 mode-shape factor is the
   * airframe's; this is the commanded tip twist). */
  warpTipRad: number;
  /** Cradle lateral offset [m] (hip drive, ±0.12 m at full warp). */
  cradleOffsetM: number;
  /** Left/right prop blade angles [rad] — counter-rotating pair. */
  leftPropRad: number;
  rightPropRad: number;
  /** Any command hit a mechanical stop this tick (reported, not silent). */
  clamped: boolean;
  /** Geometry beyond ±25% of the reference design (schematic preview). */
  schematicPreview: boolean;
}

const DEG = Math.PI / 180;

function clampReport(value: number, limit: number): [number, boolean] {
  if (value > limit) return [limit, true];
  if (value < -limit) return [-limit, true];
  return [value, false];
}

/** Compute the pose. `dims` defaults to the reference; passing modified
 * design dims drives the schematic-preview flag. */
export function computePose(
  c: ControlState,
  dims: { span_m: number; chord_m: number; canard_area_m2: number; rudder_area_m2: number } = REFERENCE_DIMS,
): FlyerPose {
  for (const [k, v] of Object.entries(c)) {
    if (typeof v === "number" && !Number.isFinite(v)) {
      throw new RangeError(`control ${k} must be finite, got ${v}`);
    }
  }
  const [canard, c1] = clampReport(c.canardDeg, CANARD_TRAVEL_DEG);
  const [warp, c2] = clampReport(c.warpDeg, WARP_LIMIT_DEG);
  const rudderCmd = c.coupled ? RUDDER_SLAVING * warp : c.rudderDeg;
  const [rudder, c3] = clampReport(rudderCmd, RUDDER_LIMIT_DEG);
  let schematicPreview = false;
  for (const key of Object.keys(REFERENCE_DIMS) as (keyof typeof REFERENCE_DIMS)[]) {
    const rel = Math.abs(dims[key] - REFERENCE_DIMS[key]) / REFERENCE_DIMS[key];
    if (rel > SCHEMATIC_PREVIEW_FRACTION) schematicPreview = true;
  }
  return {
    canardRad: canard * DEG,
    rudderRad: rudder * DEG,
    warpTipRad: warp * DEG,
    cradleOffsetM: (warp / WARP_LIMIT_DEG) * 0.12,
    leftPropRad: c.propAngleRad,
    rightPropRad: -c.propAngleRad, // crossed chain: counter-rotation
    clamped: c1 || c2 || c3,
    schematicPreview,
  };
}
