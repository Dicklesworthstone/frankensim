// WHY IT ROLLS AND YAWS view logic (bead wf-root-guzez.8.6, E7.4b).
// Warp command vs loaded twist, rudder, moment traces, sideslip and
// roll rate, spiral-mode indicator, ADVERSE-YAW decomposition, the
// lateral cue, reversal events, and the A/B warp-rudder linkage
// toggle that preserves the atmosphere realization + input prefix.
//
// Laws (siblings of the porpoises view's):
//   - moment decompositions are PUBLISHED rows summed by the sim
//     plane; the view checks the sum against the published net and
//     refuses a broken split (never patches it);
//   - the adverse-yaw verdict is a SIGN law on the decomposition:
//     with the rudder linkage DECOUPLED (the 1901 configuration),
//     the induced-drag yaw component opposes the warp roll command;
//   - the A/B toggle compares linkage-coupled vs decoupled runs that
//     MUST share realization id and input prefix — mismatches refuse.

export interface LateralRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type LateralResult<T> = { ok: true; value: T } | { ok: false; refusal: LateralRefusal };

function refuse<T>(code: string, message: string, repair: string): LateralResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** One tick's published yaw-moment decomposition [N·m]. */
export interface YawDecomposition {
  readonly tick: number;
  /** warp command [rad] (sign: +right-wing-down). */
  readonly warpCommandRad: number;
  /** loaded (aeroelastic) twist actually achieved [rad]. */
  readonly loadedTwistRad: number;
  /** induced-drag differential yaw moment. */
  readonly inducedDragYawNm: number;
  /** rudder yaw moment. */
  readonly rudderYawNm: number;
  /** profile/other yaw moment. */
  readonly profileYawNm: number;
  /** the sim plane's OWN net yaw moment. */
  readonly netYawNm: number;
}

/** Sum-check tolerance relative to the largest component. */
export const DECOMPOSITION_REL_TOL = 1e-9;

/** Window cap. */
export const MAX_WINDOW = 4_096;

/**
 * Validate a decomposition window: published components must sum to
 * the published net (a broken split is a sim-plane bug the view
 * REFUSES to render, never patches).
 */
export function validateDecomposition(
  rows: readonly YawDecomposition[],
): LateralResult<readonly YawDecomposition[]> {
  if (rows.length === 0 || rows.length > MAX_WINDOW) {
    return refuse(
      "lateral-window-invalid",
      `${rows.length} rows outside [1, ${MAX_WINDOW}]`,
      "window the log",
    );
  }
  for (const r of rows) {
    const sum = r.inducedDragYawNm + r.rudderYawNm + r.profileYawNm;
    const scale = Math.max(
      Math.abs(r.inducedDragYawNm),
      Math.abs(r.rudderYawNm),
      Math.abs(r.profileYawNm),
      1e-9,
    );
    if (Math.abs(sum - r.netYawNm) > DECOMPOSITION_REL_TOL * scale) {
      return refuse(
        "lateral-decomposition-broken",
        `tick ${r.tick}: components sum ${sum} vs net ${r.netYawNm}`,
        "the sim plane owns the split; re-emit it",
      );
    }
  }
  return { ok: true, value: rows };
}

export interface AdverseYawVerdict {
  /** true when the induced-drag yaw OPPOSES the warp command. */
  readonly adverse: boolean;
  /** mean sign product over commanded ticks (the evidence). */
  readonly meanSignProduct: number;
  /** ticks that actually carried a command. */
  readonly commandedTicks: number;
}

/**
 * The adverse-yaw SIGN law over a window: adverse means
 * sign(inducedDragYaw) == −sign(warpCommand) on commanded ticks.
 */
export function adverseYawVerdict(
  rows: readonly YawDecomposition[],
): LateralResult<AdverseYawVerdict> {
  const v = validateDecomposition(rows);
  if (!v.ok) return v;
  let product = 0;
  let commanded = 0;
  for (const r of rows) {
    if (Math.abs(r.warpCommandRad) > 1e-6 && Math.abs(r.inducedDragYawNm) > 1e-9) {
      product += Math.sign(r.warpCommandRad) * Math.sign(r.inducedDragYawNm);
      commanded += 1;
    }
  }
  if (commanded === 0) {
    return refuse(
      "lateral-no-commanded-ticks",
      "no ticks carry both a warp command and an induced-drag moment",
      "the verdict needs a commanded window (non-vacuity)",
    );
  }
  const mean = product / commanded;
  return {
    ok: true,
    value: { adverse: mean < -0.5, meanSignProduct: mean, commandedTicks: commanded },
  };
}

/** Spiral-mode indicator from a published lateral pole (passthrough). */
export interface SpiralIndicator {
  readonly reSigmaPerS: number;
  /** time-to-double for divergent spiral; null when convergent. */
  readonly timeToDoubleS: number | null;
  readonly divergent: boolean;
}

/** Build the spiral indicator (the one display derivation). */
export function spiralIndicator(reSigmaPerS: number): LateralResult<SpiralIndicator> {
  if (!Number.isFinite(reSigmaPerS)) {
    return refuse("lateral-pole-invalid", `sigma ${reSigmaPerS}`, "finite pole");
  }
  return {
    ok: true,
    value: {
      reSigmaPerS,
      timeToDoubleS: reSigmaPerS > 0 ? Math.LN2 / reSigmaPerS : null,
      divergent: reSigmaPerS > 0,
    },
  };
}

/** One A/B run header for the linkage toggle. */
export interface LinkageRunHeader {
  /** warp-rudder linkage coupled (1902+) or decoupled (1901 mode)? */
  readonly linkageCoupled: boolean;
  /** the atmosphere realization id. */
  readonly realizationId: string;
  /** digest of the input prefix both runs must share. */
  readonly inputPrefixDigest: string;
}

/**
 * Admit an A/B linkage pair: SAME realization, SAME input prefix,
 * DIFFERENT linkage — anything else refuses (a comparison that
 * varies two things at once attributes nothing).
 */
export function admitLinkagePair(
  a: LinkageRunHeader,
  b: LinkageRunHeader,
): LateralResult<{ coupled: LinkageRunHeader; decoupled: LinkageRunHeader }> {
  if (a.realizationId === "" || a.realizationId !== b.realizationId) {
    return refuse(
      "linkage-ab-realization-mismatch",
      `'${a.realizationId}' vs '${b.realizationId}'`,
      "the A/B holds the atmosphere realization fixed",
    );
  }
  if (a.inputPrefixDigest === "" || a.inputPrefixDigest !== b.inputPrefixDigest) {
    return refuse(
      "linkage-ab-prefix-mismatch",
      "input prefixes differ",
      "the A/B holds the input prefix fixed",
    );
  }
  if (a.linkageCoupled === b.linkageCoupled) {
    return refuse(
      "linkage-ab-degenerate",
      "both runs have the same linkage",
      "toggle exactly the linkage",
    );
  }
  return {
    ok: true,
    value: a.linkageCoupled ? { coupled: a, decoupled: b } : { coupled: b, decoupled: a },
  };
}
