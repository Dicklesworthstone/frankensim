// Launch-ceremony presentation envelopes (bead guzez.16). PURE and
// stateless-in-(inputs, t): the onboard-camera glance toward Orville,
// the release impulse (FOV punch + shake), and Daniels' flashbulb
// moment. Presentation-plane ONLY — none of these numbers ever touch
// physics; the scene feeds latched times and gets weights. No THREE.js
// types leak into this module. Repro: node --test test/ceremony.test.ts
//
// World frame matches dressing.ts: flight runs WORLD +x from the launch
// origin; Orville chases at the right wingtip (+z side,
// ORVILLE_SIDE_OFFSET_M) while phase='on-rail' and releases at the
// LATCHED (releaseX, releaseT).

/* ------------------------- shared easing ---------------------------- */

/** Smoothstep fade [0,1]: the house ramp curve — flat at both ends, so
 * blends never pop when they start or land. */
function smooth01(u: number): number {
  const c = Math.min(1, Math.max(0, u));
  return c * c * (3 - 2 * c);
}

/** Domain guard in the dressing.ts idiom: ceremony times are seconds
 * measured forward from a latch event; anything negative or
 * non-finite is a caller bug, refused loudly rather than blended. */
function requireAge(name: string, v: number | null): void {
  if (v !== null && (!Number.isFinite(v) || v < 0)) {
    throw new RangeError(`${name} must be null or a finite non-negative time, got ${v}`);
  }
}

/* --------------------- onboard wingtip glance ----------------------- */

/** Nominal on-rail duration [s] the ceremony choreography plans
 * against (the historical Dec-17 rail runs were this order; the real
 * release is latched whenever physics actually lifts). */
export const NOMINAL_RAIL_RUN_S = 3;

/** Length of the glance ramp [s]: the FINAL stretch ON the rail. */
export const RAMP_S = 1.2;
/** Fade-back time [s] after release: the camera swings forward again. */
export const RELEASE_DECAY_S = 0.5;

/**
 * Onboard-camera glance weight toward Orville at the right wingtip.
 *
 * Feed `elapsedOnRailS` (seconds since the rail run began, or null
 * when off the rail) and `sinceReleaseS` (seconds since the latched
 * release, or null before it). While on the rail the weight rides a
 * smoothstep up across the LAST `RAMP_S` before the nominal release;
 * at release it decays back to 0 over `RELEASE_DECAY_S` — the "snap"
 * is the immediate handoff onto that fast decay, so the swing home
 * starts the instant Orville lets go instead of holding. Both inputs
 * null (never boarded) -> 0. Pure; replays identical.
 */
export function glanceBlend(
  elapsedOnRailS: number | null,
  sinceReleaseS: number | null,
): number {
  requireAge("elapsedOnRailS", elapsedOnRailS);
  requireAge("sinceReleaseS", sinceReleaseS);
  let w = 0;
  if (elapsedOnRailS !== null) {
    const rampStart = NOMINAL_RAIL_RUN_S - RAMP_S;
    w = smooth01((elapsedOnRailS - rampStart) / RAMP_S);
  }
  if (sinceReleaseS !== null) {
    // Post-release decay wins: it takes over exactly where the ramp
    // peaked (weight 1 at the latch) and pulls the camera home.
    w = Math.max(w, 1 - smooth01(sinceReleaseS / RELEASE_DECAY_S));
  }
  return w;
}

/** First-form convenience (contract spelling): absolute sim clock plus
 * the phase flag and the seconds-until-release the scene projects.
 * Off-rail -> 0 (the glance SNAPS home at release); on-rail it ramps
 * 0->1 as `releaseImminentT` falls through the final `RAMP_S`. */
export function wingtipGlanceBlend(
  nowS: number,
  onRail: boolean,
  releaseImminentT?: number | null,
): number {
  requireAge("nowS", nowS);
  requireAge("releaseImminentT", releaseImminentT ?? null);
  if (!onRail) return 0;
  if (releaseImminentT === null || releaseImminentT === undefined) return 0;
  return smooth01((RAMP_S - releaseImminentT) / RAMP_S);
}

/* ------------------------- release impulse -------------------------- */

/** Peak FOV punch [deg] at the instant of release. */
export const KICK_FOV_PEAK_DEG = 6;
/** Peak camera-shake amplitude [m] at the instant of release. */
export const KICK_SHAKE_PEAK_M = 0.15;
/** Impulse length [s]: critically-damped exponential settling. */
export const KICK_T_S = 0.8;
/** Exponential rate so the envelope settles through ~1.6% at KICK_T_S
 * (four natural time constants — visually done, mathematically clean). */
const KICK_RATE = 4 / KICK_T_S;

/** One impulse channel's critically-damped envelope [0,1]: exactly 1
 * AT the event, monotone exponential settle, hard 0 at and past
 * KICK_T_S (closed envelope — no asymptotic smear into later frames). */
function kickEnv(sinceReleaseS: number): number {
  if (sinceReleaseS >= KICK_T_S) return 0;
  return Math.exp(-KICK_RATE * sinceReleaseS);
}

export interface ReleaseKick {
  /** FOV punch remaining [deg]; subtract/add around the preset FOV. */
  fovKickDeg: number;
  /** Shake amplitude remaining [m]; the scene turns it into jitter. */
  shakeAmpM: number;
}

/**
 * Release-moment impulse: brief FOV punch plus decaying shake over
 * `KICK_T_S`. Peaks EXACTLY at `sinceReleaseS = 0`; zero before
 * release (`null`). Pure — the scene calls it every frame with the
 * latched age; replays see the identical punch.
 */
export function releaseKick(sinceReleaseS: number | null): ReleaseKick {
  requireAge("sinceReleaseS", sinceReleaseS);
  if (sinceReleaseS === null) return { fovKickDeg: 0, shakeAmpM: 0 };
  const env = kickEnv(sinceReleaseS);
  return { fovKickDeg: KICK_FOV_PEAK_DEG * env, shakeAmpM: KICK_SHAKE_PEAK_M * env };
}

/* ---------------------- Daniels' flashbulb moment ------------------- */

/** Flash attack time [s]: shutter-open rise on the cockpit overlay. */
export const FLASH_ATTACK_S = 0.06;
/** Flash decay time [s] measured AFTER the attack completes. */
export const FLASH_DECAY_S = 0.7;
/** Decay shape constant: normalized exponential hitting EXACTLY 0 at
 * attack + FLASH_DECAY_S (closed envelope like the kick). */
const FLASH_SHAPE = 4;

/**
 * Screen-flash opacity [0,1] for the cockpit overlay — Daniels'
 * one-shot magnesium flash. Smoothstep attack over `FLASH_ATTACK_S`
 * to full white, then a normalized exponential fall reaching exactly
 * 0 at `FLASH_ATTACK_S + FLASH_DECAY_S` and staying there. Fired once
 * per event: the scene passes the age since the flash latch (null =
 * never fired). Pure; replays identical.
 */
export function flashPulse(sinceFlashS: number | null): number {
  requireAge("sinceFlashS", sinceFlashS);
  if (sinceFlashS === null) return 0;
  if (sinceFlashS <= 0) return 0;
  if (sinceFlashS < FLASH_ATTACK_S) {
    return smooth01(sinceFlashS / FLASH_ATTACK_S);
  }
  const sincePeak = sinceFlashS - FLASH_ATTACK_S;
  if (sincePeak >= FLASH_DECAY_S) return 0;
  const u = sincePeak / FLASH_DECAY_S;
  const k = FLASH_SHAPE;
  return (Math.exp(-k * u) - Math.exp(-k)) / (1 - Math.exp(-k));
}
