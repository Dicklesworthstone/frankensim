// Eigenmode view + teaching projection + polar redraw + design-diff
// cards (bead wf-root-guzez.9.2.2, E8.2-ii). Presentation on the
// E8.2-i augmented-linearization engine: the sim plane publishes
// labeled poles (structural-freezing attribution) and eigenvectors;
// this module groups, projects, and diffs — it never re-derives an
// eigenvalue.
//
// Laws:
//   - mode families are the ENGINE's labels (passthrough grouping);
//   - the four-state teaching projection reports its OWN residual
//     ("beyond-4-state content") and labels modes that are mostly
//     not rigid — a teaching view that hid the residual would teach
//     the wrong airplane;
//   - polar redraw is verbatim; design-diffs align on IDENTICAL
//     alpha grids (no interpolation — mismatched grids refuse);
//   - design-diff cards carry named attribution verbatim.

export interface EigenRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type EigenResult<T> = { ok: true; value: T } | { ok: false; refusal: EigenRefusal };

function refuse<T>(code: string, message: string, repair: string): EigenResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** The engine's mode families (mirror of fs-flyer::augmented). */
export type ModeFamily = "rigid" | "actuator" | "rotor" | "pilot";

/** A published labeled pole (engine passthrough). */
export interface PublishedLabeledPole {
  readonly re: number;
  readonly im: number;
  readonly family: ModeFamily;
  readonly attributionShift: number;
}

export interface FamilyGroup {
  readonly family: ModeFamily;
  readonly poles: readonly PublishedLabeledPole[];
}

/** Pole budget per view. */
export const MAX_POLES = 64;

/**
 * Group published poles by their ENGINE-labeled family (order within
 * a family preserved; family order fixed rigid→actuator→rotor→pilot).
 */
export function groupModeFamilies(
  poles: readonly PublishedLabeledPole[],
): EigenResult<readonly FamilyGroup[]> {
  if (poles.length === 0 || poles.length > MAX_POLES) {
    return refuse(
      "pole-count-invalid",
      `${poles.length} poles outside [1, ${MAX_POLES}]`,
      "the augmented engine publishes at most 64 states",
    );
  }
  if (poles.some((p) => !Number.isFinite(p.re) || !Number.isFinite(p.im))) {
    return refuse("pole-invalid", "non-finite pole", "finite poles only");
  }
  const order: ModeFamily[] = ["rigid", "actuator", "rotor", "pilot"];
  const groups = order
    .map((family) => ({ family, poles: poles.filter((p) => p.family === family) }))
    .filter((g) => g.poles.length > 0);
  return { ok: true, value: groups };
}

/** The classic four-state longitudinal basis, frozen order. */
export const FOUR_STATE_LABELS = ["u", "w", "q", "theta"] as const;

/** Below this rigid share a mode is labeled mostly-not-rigid. */
export const RIGID_SHARE_FLOOR = 0.5;

export interface TeachingProjection {
  /** |rigid 4-state part| / |full vector| in [0, 1]. */
  readonly rigidShare: number;
  /** the labeled residual: 1 - rigidShare² energy fraction. */
  readonly beyondFourStateContent: number;
  /** per-label magnitudes of the projected 4-state part. */
  readonly components: readonly { label: string; magnitude: number }[];
  /** the honesty label. */
  readonly caption: string;
}

/**
 * Project a published eigenvector onto the labeled four-state basis.
 * `stateLabels` is the engine's frozen label order; the projection
 * picks the u/w/q/theta rows. The residual is REPORTED, and low
 * rigid share flips the caption — the teaching view never dresses an
 * actuator mode as an airplane mode.
 */
export function teachingProjection(
  stateLabels: readonly string[],
  vectorMagnitudes: Float64Array,
): EigenResult<TeachingProjection> {
  if (stateLabels.length !== vectorMagnitudes.length || stateLabels.length === 0) {
    return refuse(
      "projection-shape-mismatched",
      `${stateLabels.length} labels vs ${vectorMagnitudes.length} magnitudes`,
      "one magnitude per engine state",
    );
  }
  let full2 = 0;
  for (const v of vectorMagnitudes) full2 += v * v;
  if (!(full2 > 0) || !Number.isFinite(full2)) {
    return refuse("projection-vector-degenerate", `norm² ${full2}`, "a nonzero finite mode shape");
  }
  const components: { label: string; magnitude: number }[] = [];
  let rigid2 = 0;
  for (const label of FOUR_STATE_LABELS) {
    const idx = stateLabels.indexOf(label);
    const mag = idx >= 0 ? (vectorMagnitudes[idx] ?? 0) : 0;
    components.push({ label, magnitude: mag });
    rigid2 += mag * mag;
  }
  const rigidShare = Math.sqrt(rigid2 / full2);
  const beyond = 1 - rigid2 / full2;
  return {
    ok: true,
    value: {
      rigidShare,
      beyondFourStateContent: beyond,
      components,
      caption:
        rigidShare >= RIGID_SHARE_FLOOR
          ? `four-state projection (beyond-4-state content ${(beyond * 100).toFixed(1)}%)`
          : `mostly NOT a rigid mode (rigid share ${(rigidShare * 100).toFixed(1)}%)`,
    },
  };
}

/** One published polar point. */
export interface PolarPoint {
  readonly alphaRad: number;
  readonly cl: number;
  readonly cd: number;
}

/** Polar point budget. */
export const MAX_POLAR_POINTS = 512;

/** Validate a polar for redraw (verbatim; alpha strictly ascending). */
export function validatePolar(points: readonly PolarPoint[]): EigenResult<readonly PolarPoint[]> {
  if (points.length < 2 || points.length > MAX_POLAR_POINTS) {
    return refuse(
      "polar-length-invalid",
      `${points.length} points outside [2, ${MAX_POLAR_POINTS}]`,
      "publish the sampled polar as-is",
    );
  }
  for (let i = 1; i < points.length; i += 1) {
    const prev = points[i - 1];
    const cur = points[i];
    if (prev === undefined || cur === undefined || !(cur.alphaRad > prev.alphaRad)) {
      return refuse("polar-unordered", `alpha order breaks at ${i}`, "sort upstream");
    }
  }
  return { ok: true, value: points };
}

export interface PolarDiffPoint {
  readonly alphaRad: number;
  readonly dCl: number;
  readonly dCd: number;
}

/**
 * Design-diff of two polars: IDENTICAL alpha grids required (bitwise
 * — interpolation would invent data between published samples).
 */
export function polarDiff(
  before: readonly PolarPoint[],
  after: readonly PolarPoint[],
): EigenResult<readonly PolarDiffPoint[]> {
  const vb = validatePolar(before);
  if (!vb.ok) return vb;
  const va = validatePolar(after);
  if (!va.ok) return va;
  if (before.length !== after.length) {
    return refuse(
      "polar-grids-mismatched",
      `${before.length} vs ${after.length} points`,
      "re-sample the candidate on the baseline grid upstream",
    );
  }
  const out: PolarDiffPoint[] = [];
  for (let i = 0; i < before.length; i += 1) {
    const b = before[i];
    const a = after[i];
    if (b === undefined || a === undefined || !Object.is(b.alphaRad, a.alphaRad)) {
      return refuse(
        "polar-grids-mismatched",
        `alpha grids differ at ${i}`,
        "re-sample the candidate on the baseline grid upstream",
      );
    }
    out.push({ alphaRad: b.alphaRad, dCl: a.cl - b.cl, dCd: a.cd - b.cd });
  }
  return { ok: true, value: out };
}

/** A design-diff card: named change + verbatim metric deltas. */
export interface DesignDiffCard {
  readonly change: string;
  readonly attribution: string;
  readonly metrics: readonly { metric: string; before: number; after: number }[];
}

/** Card metric budget. */
export const MAX_CARD_METRICS = 32;

/** Validate a card (verbatim display; attribution must be named). */
export function validateCard(card: DesignDiffCard): EigenResult<DesignDiffCard> {
  if (card.change.trim() === "" || card.attribution.trim() === "") {
    return refuse(
      "card-attribution-missing",
      "a design diff without a named change/attribution",
      "name the change and its attributed mechanism",
    );
  }
  if (card.metrics.length === 0 || card.metrics.length > MAX_CARD_METRICS) {
    return refuse(
      "card-metrics-invalid",
      `${card.metrics.length} metrics outside [1, ${MAX_CARD_METRICS}]`,
      "cards summarize, logs enumerate",
    );
  }
  if (card.metrics.some((m) => !Number.isFinite(m.before) || !Number.isFinite(m.after))) {
    return refuse("card-metrics-invalid", "non-finite metric", "finite metrics only");
  }
  return { ok: true, value: card };
}
