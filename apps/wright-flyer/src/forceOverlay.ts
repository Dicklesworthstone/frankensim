// Force overlay + strip loads + strip-chart probes (bead
// wf-root-guzez.8.4, E7.3). THE LAW: overlay values equal sim-plane
// state BIT-FOR-BIT — this module arranges and forwards, it never
// recomputes a force. The net gnomon shows the state's OWN net
// vector even when a doctored state's components do not sum to it
// (that inconsistency is the sim plane's to explain, and the test
// battery uses exactly that to falsify renderer-side recomputation).

export interface OverlayRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type OverlayResult<T> =
  | { ok: true; value: T }
  | { ok: false; refusal: OverlayRefusal };

function refuse<T>(code: string, message: string, repair: string): OverlayResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** Per-strip aerodynamic loads as the SIM PLANE published them. */
export interface StripLoadsState {
  readonly nStrips: number;
  /** application point xyz per strip [m]. */
  readonly positions: Float64Array;
  /** force vector xyz per strip [N]. */
  readonly forces: Float64Array;
  /** prop thrust vector [N] + application point [m]. */
  readonly thrustN: readonly [number, number, number];
  readonly thrustAt: readonly [number, number, number];
  /** weight [N] applied at the CG [m]. */
  readonly weightN: readonly [number, number, number];
  readonly cgAt: readonly [number, number, number];
  /** the sim plane's OWN net force [N] (never recomputed here). */
  readonly netN: readonly [number, number, number];
}

/** Strip budget (1903 twin wings + canard + rudder rows). */
export const MAX_STRIPS = 256;

export interface ForceOverlay {
  /** per-strip arrows: position xyz + vector xyz, VERBATIM. */
  readonly stripPositions: Float64Array;
  readonly stripVectors: Float64Array;
  readonly thrust: { at: readonly [number, number, number]; vec: readonly [number, number, number] };
  readonly weight: { at: readonly [number, number, number]; vec: readonly [number, number, number] };
  /** the net gnomon: the STATE's net, verbatim. */
  readonly net: readonly [number, number, number];
}

/**
 * Arrange the overlay. Every numeric value is a verbatim copy of the
 * state — the battery asserts bit-for-bit equality per strip.
 */
export function buildForceOverlay(state: StripLoadsState): OverlayResult<ForceOverlay> {
  if (state.nStrips < 1 || state.nStrips > MAX_STRIPS) {
    return refuse(
      "strip-count-invalid",
      `${state.nStrips} strips outside [1, ${MAX_STRIPS}]`,
      "the sim plane publishes at most 256 strips",
    );
  }
  if (
    state.positions.length !== 3 * state.nStrips ||
    state.forces.length !== 3 * state.nStrips
  ) {
    return refuse(
      "strip-arrays-mismatched",
      `${state.positions.length}/${state.forces.length} for ${state.nStrips} strips`,
      "3 floats per strip in both arrays",
    );
  }
  // Verbatim copies (same bits; a copy so later state mutation cannot
  // silently retint the overlay).
  return {
    ok: true,
    value: {
      stripPositions: Float64Array.from(state.positions),
      stripVectors: Float64Array.from(state.forces),
      thrust: { at: state.thrustAt, vec: state.thrustN },
      weight: { at: state.cgAt, vec: state.weightN },
      net: state.netN,
    },
  };
}

/**
 * The bit-for-bit audit (the DONE-WHEN oracle, exported so the HUD
 * self-check can run it live): every overlay value must equal its
 * state source EXACTLY. Returns the index of the first divergence,
 * or -1 when faithful.
 */
export function firstOverlayDivergence(overlay: ForceOverlay, state: StripLoadsState): number {
  for (let i = 0; i < state.forces.length; i += 1) {
    const a = overlay.stripVectors[i] ?? Number.NaN;
    const b = state.forces[i] ?? Number.NaN;
    if (!Object.is(a, b)) return i;
  }
  for (let c = 0; c < 3; c += 1) {
    if (!Object.is(overlay.net[c], state.netN[c])) return state.forces.length + c;
    if (!Object.is(overlay.thrust.vec[c], state.thrustN[c])) return state.forces.length + 3 + c;
    if (!Object.is(overlay.weight.vec[c], state.weightN[c])) return state.forces.length + 6 + c;
  }
  return -1;
}

/** Probe chart budgets. */
export const MAX_CHART_SAMPLES = 4_096;

/** One strip-chart probe: a declared reference height + its series. */
export interface ProbeChart {
  readonly label: string;
  /** declared reference height above the certified plane [m]. */
  readonly referenceHeightM: number;
  /** (tick, value) pairs, tick-ordered. */
  readonly ticks: Float64Array;
  readonly values: Float64Array;
}

/**
 * Build a probe chart from the LOGGED series verbatim (the chart
 * matches the log or it refuses — no smoothing, no resampling).
 */
export function buildProbeChart(
  label: string,
  referenceHeightM: number,
  ticks: Float64Array,
  values: Float64Array,
): OverlayResult<ProbeChart> {
  if (!Number.isFinite(referenceHeightM)) {
    return refuse(
      "probe-reference-invalid",
      `reference height ${referenceHeightM}`,
      "declare a finite reference height",
    );
  }
  if (ticks.length !== values.length) {
    return refuse(
      "probe-series-mismatched",
      `${ticks.length} ticks vs ${values.length} values`,
      "one value per tick",
    );
  }
  if (ticks.length === 0 || ticks.length > MAX_CHART_SAMPLES) {
    return refuse(
      "probe-series-length-invalid",
      `${ticks.length} samples outside [1, ${MAX_CHART_SAMPLES}]`,
      "window the log to the chart budget",
    );
  }
  for (let i = 1; i < ticks.length; i += 1) {
    if (!((ticks[i] ?? Number.NaN) > (ticks[i - 1] ?? Number.NaN))) {
      return refuse(
        "probe-series-unordered",
        `tick order breaks at ${i}`,
        "charts render the log as recorded; sort/split upstream",
      );
    }
  }
  return {
    ok: true,
    value: {
      label,
      referenceHeightM,
      ticks: Float64Array.from(ticks),
      values: Float64Array.from(values),
    },
  };
}

/**
 * Chart-vs-log audit (DONE-WHEN oracle): the chart's series must be
 * the logged series bit-for-bit. First divergence index, or -1.
 */
export function firstChartDivergence(
  chart: ProbeChart,
  loggedTicks: Float64Array,
  loggedValues: Float64Array,
): number {
  if (chart.ticks.length !== loggedTicks.length) return 0;
  for (let i = 0; i < loggedTicks.length; i += 1) {
    if (!Object.is(chart.ticks[i], loggedTicks[i])) return i;
    if (!Object.is(chart.values[i], loggedValues[i])) return i;
  }
  return -1;
}
