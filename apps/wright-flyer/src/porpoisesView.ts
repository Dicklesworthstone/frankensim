// WHY IT PORPOISES flagship view logic (bead wf-root-guzez.8.5,
// E7.4). Synchronized: open-loop pole / time-to-double, canard
// command-vs-actual with the delay/phase indicator, hinge moment +
// pilot force, pitch / flight-path / height traces, saturation and
// reversal event markers, live loop-component attribution, and A/B
// with the atmosphere realization held fixed (common-prefix
// semantics from ABComparisonReceiptV1).
//
// Laws:
//   - traces and events are the LOGGED series verbatim (the E7.3
//     passthrough law; charts reuse that module's audits);
//   - time-to-double is the one display DERIVATION allowed here —
//     ln 2 / Re(pole) for unstable poles, absent otherwise — and the
//     battery pins it against closed form;
//   - attribution shares are passthrough; the view only ranks them,
//     and mislabeled shares are caught by per-component oracles;
//   - the A/B lane never annotates divergence before the receipt's
//     divergence tick, and a receipt inconsistent with the traces it
//     claims to describe is REFUSED, not smoothed over.

export interface ViewRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type ViewResult<T> = { ok: true; value: T } | { ok: false; refusal: ViewRefusal };

function refuse<T>(code: string, message: string, repair: string): ViewResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** An open-loop pole the sim plane published (continuous-time). */
export interface PublishedPole {
  readonly reSigmaPerS: number;
  readonly imOmegaRadPerS: number;
}

export interface PoleIndicator {
  readonly pole: PublishedPole;
  /** ln2 / sigma for unstable poles; null when non-growing. */
  readonly timeToDoubleS: number | null;
  /** oscillation period 2*pi/|omega|; null for aperiodic. */
  readonly periodS: number | null;
}

/** Time-to-double: the ONE derivation this view performs. */
export function poleIndicator(pole: PublishedPole): ViewResult<PoleIndicator> {
  if (!Number.isFinite(pole.reSigmaPerS) || !Number.isFinite(pole.imOmegaRadPerS)) {
    return refuse("pole-invalid", `(${pole.reSigmaPerS}, ${pole.imOmegaRadPerS})`, "finite pole");
  }
  return {
    ok: true,
    value: {
      pole,
      timeToDoubleS: pole.reSigmaPerS > 0 ? Math.LN2 / pole.reSigmaPerS : null,
      periodS:
        pole.imOmegaRadPerS !== 0 ? (2 * Math.PI) / Math.abs(pole.imOmegaRadPerS) : null,
    },
  };
}

/** Delay-estimate search half-width [ticks]. */
export const MAX_LAG_TICKS = 120;

/**
 * Delay/phase indicator: the lag (in ticks) maximizing the
 * cross-correlation of command vs actual over ±MAX_LAG_TICKS.
 * Deterministic: first maximum wins on ties (scan order -L..+L).
 */
export function estimateDelayTicks(
  command: Float64Array,
  actual: Float64Array,
): ViewResult<number> {
  const n = command.length;
  if (n !== actual.length || n < 2 * MAX_LAG_TICKS + 4) {
    return refuse(
      "delay-window-invalid",
      `${n} vs ${actual.length} samples (need > ${2 * MAX_LAG_TICKS + 3})`,
      "window the traces before estimating delay",
    );
  }
  let bestLag = -MAX_LAG_TICKS;
  let bestScore = Number.NEGATIVE_INFINITY;
  for (let lag = -MAX_LAG_TICKS; lag <= MAX_LAG_TICKS; lag += 1) {
    let s = 0;
    for (let i = MAX_LAG_TICKS; i < n - MAX_LAG_TICKS; i += 1) {
      s += (command[i] ?? 0) * (actual[i + lag] ?? 0);
    }
    if (s > bestScore) {
      bestScore = s;
      bestLag = lag;
    }
  }
  return { ok: true, value: bestLag };
}

/** A saturation/reversal event as logged (verbatim). */
export interface LoopEvent {
  readonly tick: number;
  readonly kind: "saturation-enter" | "saturation-exit" | "command-reversal";
}

/** Event budget per window. */
export const MAX_EVENTS = 1_024;

/** Validate the logged event list (verbatim display; order enforced). */
export function validateEvents(events: readonly LoopEvent[]): ViewResult<readonly LoopEvent[]> {
  if (events.length > MAX_EVENTS) {
    return refuse(
      "event-count-exceeded",
      `${events.length} events > ${MAX_EVENTS}`,
      "window the log",
    );
  }
  for (let i = 1; i < events.length; i += 1) {
    const prev = events[i - 1];
    const cur = events[i];
    if (prev !== undefined && cur !== undefined && cur.tick < prev.tick) {
      return refuse(
        "events-unordered",
        `tick order breaks at ${i}`,
        "events render as logged; sort upstream",
      );
    }
  }
  return { ok: true, value: events };
}

/** One loop component's published attribution share. */
export interface AttributionShare {
  readonly component: string;
  readonly share: number;
}

export interface AttributionView {
  /** shares ranked descending (stable tie-break by input order). */
  readonly ranked: readonly AttributionShare[];
  readonly dominant: string;
  /** |1 - sum(shares)| — displayed, never hidden. */
  readonly residual: number;
}

/** Residual bound beyond which the attribution is refused as broken. */
export const MAX_ATTRIBUTION_RESIDUAL = 0.05;

/**
 * Rank the published attribution shares. The view never recomputes a
 * share; it displays the residual and refuses only when the shares
 * cannot describe a loop at all.
 */
export function attributionView(
  shares: readonly AttributionShare[],
): ViewResult<AttributionView> {
  if (shares.length === 0) {
    return refuse("attribution-empty", "no shares", "the sim plane publishes the loop split");
  }
  if (shares.some((s) => !Number.isFinite(s.share))) {
    return refuse("attribution-invalid", "non-finite share", "finite shares only");
  }
  const sum = shares.reduce((a, s) => a + s.share, 0);
  const residual = Math.abs(1 - sum);
  if (residual > MAX_ATTRIBUTION_RESIDUAL) {
    return refuse(
      "attribution-residual-exceeded",
      `shares sum to ${sum}`,
      "a split this broken is a sim-plane bug, not a display choice",
    );
  }
  const ranked = shares
    .map((s, i) => ({ s, i }))
    .sort((a, b) => b.s.share - a.s.share || a.i - b.i)
    .map((x) => x.s);
  const first = ranked[0];
  if (first === undefined) {
    return refuse("attribution-empty", "no shares", "unreachable");
  }
  return { ok: true, value: { ranked, dominant: first.component, residual } };
}

/** The A/B receipt fields this view consumes (ABComparisonReceiptV1). */
export interface AbReceiptView {
  /** ticks over which the two runs are bitwise identical. */
  readonly commonPrefixTicks: number;
  /** first divergent tick (== commonPrefixTicks when divergent). */
  readonly divergent: boolean;
}

export interface AbAnnotation {
  /** ticks the view may mark as shared (never past the receipt). */
  readonly sharedUntilTick: number;
  /** tick from which divergence annotations are allowed; null if none. */
  readonly divergenceFromTick: number | null;
}

/**
 * A/B annotation law: divergence markers begin AT the receipt's
 * divergence tick, never earlier; and the receipt must actually
 * describe the traces — the first sample-level difference must sit
 * exactly at the receipt boundary, else REFUSE.
 */
export function abAnnotation(
  receipt: AbReceiptView,
  traceA: Float64Array,
  traceB: Float64Array,
): ViewResult<AbAnnotation> {
  if (traceA.length !== traceB.length || traceA.length === 0) {
    return refuse(
      "ab-traces-mismatched",
      `${traceA.length} vs ${traceB.length}`,
      "compare equal windows",
    );
  }
  let firstDiff = -1;
  for (let i = 0; i < traceA.length; i += 1) {
    if (!Object.is(traceA[i], traceB[i])) {
      firstDiff = i;
      break;
    }
  }
  const boundary = Math.min(receipt.commonPrefixTicks, traceA.length);
  if (firstDiff !== -1 && firstDiff < boundary) {
    return refuse(
      "ab-receipt-inconsistent",
      `traces diverge at ${firstDiff} inside the claimed common prefix ${boundary}`,
      "re-emit the receipt from the actual frozen traces",
    );
  }
  if (!receipt.divergent && firstDiff !== -1) {
    return refuse(
      "ab-receipt-inconsistent",
      `receipt claims identical but traces diverge at ${firstDiff}`,
      "re-emit the receipt from the actual frozen traces",
    );
  }
  return {
    ok: true,
    value: {
      sharedUntilTick: boundary,
      divergenceFromTick: receipt.divergent ? receipt.commonPrefixTicks : null,
    },
  };
}
