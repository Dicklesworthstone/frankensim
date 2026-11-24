/**
 * A fast estimate of the delivery box's two hold times (minutes until the crust goes soggy, and until the food falls
 * below 135 °F), and a measured bound on how far it can be from the full model.
 *
 * The estimate is a table: the full carrier model (`stepCarrierKernel`) is run at every node of a grid over the range
 * the estimate covers (vent 0–100 %, 2–8 kg of food, 15–32 °C outside, food packed at 75–92 °C), and a question is
 * answered by multilinear interpolation between the 16 nodes round it.
 *
 * It refuses where it cannot be trusted, and the full model answers instead:
 * - outside the range, where the table has no nodes;
 * - in a cell whose corners disagree on whether the box goes soggy (or cold) within the horizon at all. There the
 *   time jumps from a few minutes to "not within three hours" over a small change of vent, and interpolating across
 *   that jump would be off by an hour or more.
 *
 * Its error is measured, not assumed: the full model is run again at a fixed set of pseudo-random points spread over
 * the range, and the largest difference found where the table answers is what the page quotes, together with how many
 * points were checked. That is the worst error found, not a proof that nothing is worse.
 */
import { createInitialCarrierState, stepCarrierKernel, type CarrierConfig } from "./physicsKernels";

/** The range the table covers, per input. */
export const FAST_RANGE = {
  vent: [0, 1],
  massKg: [2, 8],
  ambientC: [15, 32],
  initialFoodC: [75, 92],
} as const;

/** Grid nodes per axis: vent, mass, outside, packed temperature. */
const NODES = [11, 7, 8, 4] as const;
/** The longest hold either time is followed to, in minutes; a box still crisp and hot then reads as this. */
export const HOLD_HORIZON_MIN = 180;

export interface HoldTimes {
  soggyMin: number;
  below135Min: number;
}

/** The full model's hold times: the box stepped in 30 s steps until both thresholds are crossed or the horizon ends. */
export function fullModelHoldTimes(cfg: CarrierConfig): HoldTimes {
  const s = createInitialCarrierState(cfg);
  for (let i = 0; i < (HOLD_HORIZON_MIN * 60) / 30; i += 1) {
    stepCarrierKernel(s, cfg, 30.0);
    if (s.minutes_until_soggy !== null && s.minutes_until_below_135f !== null) break;
  }
  return {
    soggyMin: s.minutes_until_soggy ?? HOLD_HORIZON_MIN,
    below135Min: s.minutes_until_below_135f ?? HOLD_HORIZON_MIN,
  };
}

export interface FastEstimateTable {
  /** Worst difference from the full model found at the check points where the table answers, in minutes. */
  boundMin: HoldTimes;
  /** How many check points the table answered (and was compared at). */
  checkedPoints: number;
  /** How many check points fell in cells where it refuses. */
  refusedPoints: number;
  /** Full-model runs it took to build and check the table. */
  modelRuns: number;
}

export type FastEstimate =
  | { inRange: true; times: HoldTimes; boundMin: HoldTimes; checkedPoints: number }
  | { inRange: false; reason: "range" | "edge"; axis?: (typeof AXES)[number] };

/** How many pseudo-random points the error is measured at. */
const CHECK_POINTS = 3000;

const AXES = ["vent", "massKg", "ambientC", "initialFoodC"] as const;

function inputsOf(cfg: CarrierConfig): [number, number, number, number] {
  return [cfg.ventFraction, cfg.foodMassKg ?? 4.5, cfg.ambientC ?? 22.0, cfg.initialFoodC ?? 84.0];
}

function configAt(x: readonly number[]): CarrierConfig {
  return { ventFraction: x[0]!, foodMassKg: x[1]!, ambientC: x[2]!, initialFoodC: x[3]! };
}

function nodeValue(axis: number, i: number): number {
  const [lo, hi] = FAST_RANGE[AXES[axis]!];
  return lo + ((hi - lo) * i) / (NODES[axis]! - 1);
}

/** Builds the table, then measures its error. A few thousand full-model runs: tens of milliseconds. */
export function buildFastEstimate(): { table: FastEstimateTable; estimate: (cfg: CarrierConfig) => FastEstimate } {
  const [n0, n1, n2, n3] = NODES;
  const soggy = new Float64Array(n0 * n1 * n2 * n3);
  const below = new Float64Array(n0 * n1 * n2 * n3);
  const at = (i: number, j: number, k: number, l: number) => ((i * n1 + j) * n2 + k) * n3 + l;
  let runs = 0;
  for (let i = 0; i < n0; i++)
    for (let j = 0; j < n1; j++)
      for (let k = 0; k < n2; k++)
        for (let l = 0; l < n3; l++) {
          const t = fullModelHoldTimes(configAt([nodeValue(0, i), nodeValue(1, j), nodeValue(2, k), nodeValue(3, l)]));
          soggy[at(i, j, k, l)] = t.soggyMin;
          below[at(i, j, k, l)] = t.below135Min;
          runs++;
        }

  /** The 16 corners round x, with their weights; null where the corners disagree on reaching the horizon. */
  const interpolate = (x: readonly number[]): HoldTimes | null => {
    // Cell index and fraction along each axis.
    const idx = [0, 0, 0, 0];
    const frac = [0, 0, 0, 0];
    for (let a = 0; a < 4; a++) {
      const [lo, hi] = FAST_RANGE[AXES[a]!];
      const u = ((x[a]! - lo) / (hi - lo)) * (NODES[a]! - 1);
      const c = Math.min(NODES[a]! - 2, Math.max(0, Math.floor(u)));
      idx[a] = c;
      frac[a] = Math.min(1, Math.max(0, u - c));
    }
    let s = 0;
    let b = 0;
    let soggyAtHorizon = 0;
    let belowAtHorizon = 0;
    for (let corner = 0; corner < 16; corner++) {
      let w = 1;
      const c = [0, 0, 0, 0];
      for (let a = 0; a < 4; a++) {
        const up = (corner >> a) & 1;
        c[a] = idx[a]! + up;
        w *= up ? frac[a]! : 1 - frac[a]!;
      }
      const o = at(c[0]!, c[1]!, c[2]!, c[3]!);
      if (soggy[o]! >= HOLD_HORIZON_MIN) soggyAtHorizon++;
      if (below[o]! >= HOLD_HORIZON_MIN) belowAtHorizon++;
      s += w * soggy[o]!;
      b += w * below[o]!;
    }
    // Every corner of the cell is checked, not only those with weight: the edge is a property of the cell.
    const split = (n: number) => n > 0 && n < 16;
    if (split(soggyAtHorizon) || split(belowAtHorizon)) return null;
    return { soggyMin: s, below135Min: b };
  };

  // The error, measured against fresh full-model runs at fixed pseudo-random points over the whole range.
  let seed = 0x2545f491;
  const rand = () => {
    seed ^= seed << 13;
    seed ^= seed >>> 17;
    seed ^= seed << 5;
    return (seed >>> 0) / 4294967296;
  };
  let worstSoggy = 0;
  let worstBelow = 0;
  let checked = 0;
  let refused = 0;
  for (let n = 0; n < CHECK_POINTS; n++) {
    const x = AXES.map((name) => {
      const [lo, hi] = FAST_RANGE[name];
      return lo + (hi - lo) * rand();
    });
    const est = interpolate(x);
    if (!est) {
      refused++;
      continue;
    }
    const truth = fullModelHoldTimes(configAt(x));
    runs++;
    worstSoggy = Math.max(worstSoggy, Math.abs(est.soggyMin - truth.soggyMin));
    worstBelow = Math.max(worstBelow, Math.abs(est.below135Min - truth.below135Min));
    checked++;
  }
  const boundMin = { soggyMin: worstSoggy, below135Min: worstBelow };

  const estimate = (cfg: CarrierConfig): FastEstimate => {
    const x = inputsOf(cfg);
    for (let a = 0; a < 4; a++) {
      const [lo, hi] = FAST_RANGE[AXES[a]!];
      if (!(x[a]! >= lo && x[a]! <= hi)) return { inRange: false, reason: "range", axis: AXES[a]! };
    }
    const times = interpolate(x);
    if (!times) return { inRange: false, reason: "edge" };
    return { inRange: true, times, boundMin, checkedPoints: checked };
  };

  return { table: { boundMin, checkedPoints: checked, refusedPoints: refused, modelRuns: runs }, estimate };
}
