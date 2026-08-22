// Game-HUD analog instrument mathematics (UI overhaul, root guzez).
// PURE: needle kinematics, tick generation, arc classification, and the
// per-instrument view-model builder the SVG panel consumes. The period
// triad (anemometer, stopwatch, rev counter) reflects the instruments
// the Wrights actually carried; every other dial is stamped MODERN on
// its face so the game HUD never launders a modern aid as history.
// Repro: node --test test/gauges.test.ts

import type { PhaseWord } from "./sim/protocol.ts";

/** Angular layout + value range of one dial. Angles are degrees,
 * CLOCKWISE from 12 o'clock (SVG rotate convention). */
export interface GaugeSpec {
  readonly min: number;
  readonly max: number;
  /** Needle angle at `min`. */
  readonly startDeg: number;
  /** Total clockwise sweep from `min` to `max` (> 0). */
  readonly sweepDeg: number;
  /** Optional danger arc [from, to] in VALUE units (closed interval). */
  readonly redline: readonly [number, number] | null;
}

/** Needle angle for a value: linear in value, PINNED at the end stops
 * (an analog needle cannot leave its dial — out-of-range values sit on
 * the stop, they never wrap or extrapolate). */
export function needleDeg(value: number, spec: GaugeSpec): number {
  if (!(spec.max > spec.min) || !(spec.sweepDeg > 0)) {
    throw new RangeError(`malformed gauge spec: [${spec.min}, ${spec.max}] sweep ${spec.sweepDeg}`);
  }
  const v = Number.isNaN(value) ? spec.min : Math.min(spec.max, Math.max(spec.min, value));
  return spec.startDeg + ((v - spec.min) / (spec.max - spec.min)) * spec.sweepDeg;
}

/** True when the value sits inside the dial's danger arc (closed at
 * both ends: AT the redline is IN the redline). */
export function inRedline(value: number, spec: GaugeSpec): boolean {
  return spec.redline !== null && value >= spec.redline[0] && value <= spec.redline[1];
}

/** Major tick positions: `count` evenly spaced values from min to max
 * inclusive, each with its needle angle and label value. */
export function tickMarks(
  spec: GaugeSpec,
  count: number,
): readonly { value: number; deg: number }[] {
  if (!Number.isInteger(count) || count < 2) {
    throw new RangeError(`tick count must be an integer >= 2, got ${count}`);
  }
  const out: { value: number; deg: number }[] = [];
  for (let i = 0; i < count; i += 1) {
    const value = spec.min + ((spec.max - spec.min) * i) / (count - 1);
    out.push({ value, deg: needleDeg(value, spec) });
  }
  return out;
}

/** Stopwatch face text, mm:ss.t (the Wrights timed 12 to 59 seconds). */
export function clockText(elapsedS: number): string {
  const s = Math.max(0, elapsedS);
  const mm = Math.floor(s / 60);
  const rest = s - mm * 60;
  return `${mm}:${rest < 10 ? "0" : ""}${rest.toFixed(1)}`;
}

/** One dial's frame-ready view model. */
export interface DialView {
  readonly id: string;
  readonly label: string;
  /** "PERIOD" = an instrument the Wrights carried; "MODERN" is stamped
   * on the face (honest-history doctrine). */
  readonly provenance: "PERIOD" | "MODERN";
  readonly spec: GaugeSpec;
  readonly needleDeg: number;
  /** Numeric face text under the pivot (already formatted). */
  readonly reading: string;
  readonly unit: string;
  readonly danger: boolean;
  /** Major tick label multiplier (labels print value/labelDiv). */
  readonly majorTicks: number;
}

/** Linear (lever/warp) indicator view model. */
export interface LeverView {
  readonly id: string;
  readonly label: string;
  /** Pointer position in [-1, 1] (pinned like a needle). */
  readonly fraction: number;
  readonly reading: string;
  readonly atStop: boolean;
}

/** Everything the dial panel needs for one frame. */
export interface HudDialInputs {
  readonly airspeedMps: number;
  readonly engineRpm: number;
  readonly elapsedS: number;
  readonly hM: number;
  readonly thetaRad: number;
  readonly dcRad: number;
  readonly warpRad: number;
}

/** All-zero inputs (attract mode / no sim attached). */
export const IDLE_INPUTS: HudDialInputs = {
  airspeedMps: 0,
  engineRpm: 0,
  elapsedS: 0,
  hM: 0,
  thetaRad: 0,
  dcRad: 0,
  warpRad: 0,
};

const MPH_PER_MPS = 2.23694;
const FT_PER_M = 3.28084;
const DEG_PER_RAD = 180 / Math.PI;
/** Canard travel at the pose stop (±30°) in radians. */
export const CANARD_STOP_RAD = (30 * Math.PI) / 180;
/** Warp travel at the stop (±8.5°) in radians. */
export const WARP_STOP_RAD = (8.5 * Math.PI) / 180;

export const ANEMOMETER_SPEC: GaugeSpec = {
  min: 0,
  max: 45,
  startDeg: -135,
  sweepDeg: 270,
  redline: null,
};
export const REV_COUNTER_SPEC: GaugeSpec = {
  min: 0,
  max: 1600,
  startDeg: -135,
  sweepDeg: 270,
  redline: [1250, 1600],
};
export const ALTIMETER_SPEC: GaugeSpec = {
  min: 0,
  max: 200,
  startDeg: -135,
  sweepDeg: 270,
  redline: null,
};
export const INCLINOMETER_SPEC: GaugeSpec = {
  min: -30,
  max: 30,
  startDeg: -90,
  sweepDeg: 180,
  redline: [15, 30],
};
export const STOPWATCH_SPEC: GaugeSpec = {
  min: 0,
  max: 60,
  startDeg: 0,
  sweepDeg: 360,
  redline: null,
};

/** Build the five dial view models for one frame. */
export function dialSetFrom(input: HudDialInputs): readonly DialView[] {
  const mph = input.airspeedMps * MPH_PER_MPS;
  const ft = input.hM * FT_PER_M;
  const pitchDeg = input.thetaRad * DEG_PER_RAD;
  const secondHand = input.elapsedS % 60;
  return [
    {
      id: "anemometer",
      label: "ANEMOMETER",
      provenance: "PERIOD",
      spec: ANEMOMETER_SPEC,
      needleDeg: needleDeg(mph, ANEMOMETER_SPEC),
      reading: mph.toFixed(0),
      unit: "mph",
      danger: inRedline(mph, ANEMOMETER_SPEC),
      majorTicks: 10,
    },
    {
      id: "revcounter",
      label: "ENGINE",
      provenance: "PERIOD",
      spec: REV_COUNTER_SPEC,
      needleDeg: needleDeg(input.engineRpm, REV_COUNTER_SPEC),
      reading: Math.round(input.engineRpm).toString(),
      unit: "rpm",
      danger: inRedline(input.engineRpm, REV_COUNTER_SPEC),
      majorTicks: 9,
    },
    {
      id: "stopwatch",
      label: "STOPWATCH",
      provenance: "PERIOD",
      spec: STOPWATCH_SPEC,
      needleDeg: needleDeg(secondHand, STOPWATCH_SPEC),
      reading: clockText(input.elapsedS),
      unit: "",
      danger: false,
      majorTicks: 13,
    },
    {
      id: "altimeter",
      label: "ALTITUDE",
      provenance: "MODERN",
      spec: ALTIMETER_SPEC,
      needleDeg: needleDeg(ft, ALTIMETER_SPEC),
      reading: ft.toFixed(0),
      unit: "ft",
      danger: false,
      majorTicks: 9,
    },
    {
      id: "inclinometer",
      label: "PITCH",
      provenance: "MODERN",
      spec: INCLINOMETER_SPEC,
      needleDeg: needleDeg(pitchDeg, INCLINOMETER_SPEC),
      reading: pitchDeg.toFixed(1),
      unit: "°",
      danger: inRedline(Math.abs(pitchDeg), INCLINOMETER_SPEC),
      majorTicks: 7,
    },
  ];
}

/** Build the two control-position indicators for one frame. */
export function leverSetFrom(input: HudDialInputs): readonly LeverView[] {
  const pin = (v: number): number => Math.min(1, Math.max(-1, Number.isNaN(v) ? 0 : v));
  const canard = pin(input.dcRad / CANARD_STOP_RAD);
  const warp = pin(input.warpRad / WARP_STOP_RAD);
  return [
    {
      id: "canard",
      label: "CANARD LEVER",
      fraction: canard,
      reading: `${(input.dcRad * DEG_PER_RAD).toFixed(1)}°`,
      atStop: Math.abs(canard) === 1,
    },
    {
      id: "warp",
      label: "WING WARP",
      fraction: warp,
      reading: `${(input.warpRad * DEG_PER_RAD).toFixed(1)}°`,
      atStop: Math.abs(warp) === 1,
    },
  ];
}

/** Phase banner view: what the big top-center card says and how it is
 * toned. In-flight phases get short game-style directives; terminal
 * phases pass the receipted banner text through UNCHANGED (the honest
 * refusal surface stays exactly as loud as before). */
export function phaseDisplay(
  phase: PhaseWord,
  terminalBanner: string | null,
): { text: string; tone: "info" | "good" | "end" } {
  if (terminalBanner !== null) {
    return { text: terminalBanner, tone: "end" };
  }
  switch (phase) {
    case "on-rail":
      return { text: "RAIL RUN — hold her steady", tone: "info" };
    case "airborne":
      return { text: "AIRBORNE", tone: "good" };
    default:
      // A terminal phase whose banner text has not arrived yet.
      return { text: phase.toUpperCase(), tone: "end" };
  }
}
