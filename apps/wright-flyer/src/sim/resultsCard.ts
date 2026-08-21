// E5.5 results card (bead wf-root-guzez.6.8): PURE KPI computation
// from the sim-plane transcript (a sealed FlightRecording) — the card
// NEVER invents numbers the transcript cannot support, distinguishes
// DOWNRANGE distance from PATH LENGTH, and presents the historical
// Dec-17 flights as DISTRIBUTION CONTEXT lines (both lineages, per
// the WindReference doctrine), never as single "beat this" targets.
// The KPI-vs-recompute twin (battery) recomputes every value through
// an independent pass and refuses a card that disagrees.
// Design-diff/attribution cards are E8.2-ii scope (recorded).

import {
  P_DC_RAD,
  P_H_M,
  P_PHASE,
  P_Q_RAD_S,
  P_THETA_RAD,
  P_U_MPS,
  P_W_MPS,
  P_X_M,
  PAYLOAD_F64S,
  PHASE_CODES,
} from "./protocol.ts";
import type { FlightRecording } from "./replay.ts";

export interface FlightKpis {
  /** Snapshot count the KPIs were computed over (transparency). */
  readonly frames: number;
  /** Terminal phase word (from the sealed recording). */
  readonly terminal: string;
  /** Downrange ground distance at the end [m]. */
  readonly downrangeM: number;
  /** Path length integrated over (x, h) snapshots [m] (>= downrange
   * for any flight that climbs — the §9.2 distinction). */
  readonly pathLengthM: number;
  /** Liftoff: first airborne snapshot (null = never left the rail). */
  readonly liftoff: { tick: number; xM: number } | null;
  /** Airborne duration [s] (transcript ticks at 120 Hz). */
  readonly airborneS: number;
  /** Undulations (pitch-rate sign-flip pairs while airborne). */
  readonly undulations: number;
  /** Ride metrics. */
  readonly maxAbsQRadS: number;
  readonly rmsQRadS: number;
  readonly maxAbsThetaRad: number;
  /** Peak height above ground [m]. */
  readonly maxHM: number;
  /** Control activity: total canard travel [rad] (sum |Δdc|). */
  readonly canardTravelRad: number;
  /** Minimum airspeed while airborne [m/s] (separation-margin proxy
   * at this tier — true separation margins need the strip-regime
   * panel, E7-class; declared). */
  readonly minAirspeedMps: number | null;
}

const PHASE_BY_CODE: readonly string[] = (() => {
  const out: string[] = [];
  for (const [word, code] of Object.entries(PHASE_CODES)) {
    out[code] = word;
  }
  return out;
})();

/** Compute the KPIs from the sealed transcript (pure, deterministic). */
export function computeKpis(rec: FlightRecording): FlightKpis {
  const n = rec.ticks.length;
  const at = (i: number, slot: number): number => rec.frames[i * PAYLOAD_F64S + slot]!;
  let pathLength = 0;
  let liftoff: { tick: number; xM: number } | null = null;
  let airborneFrames = 0;
  let flips = 0;
  let lastSign = 0;
  let maxAbsQ = 0;
  let sumQ2 = 0;
  let maxAbsTheta = 0;
  let maxH = 0;
  let canardTravel = 0;
  let minAirspeed: number | null = null;
  for (let i = 0; i < n; i += 1) {
    if (i > 0) {
      const dx = at(i, P_X_M) - at(i - 1, P_X_M);
      const dh = at(i, P_H_M) - at(i - 1, P_H_M);
      pathLength += Math.hypot(dx, dh);
      canardTravel += Math.abs(at(i, P_DC_RAD) - at(i - 1, P_DC_RAD));
    }
    maxH = Math.max(maxH, at(i, P_H_M));
    const phase = at(i, P_PHASE);
    if (phase === 1) {
      airborneFrames += 1;
      if (liftoff === null) {
        liftoff = { tick: rec.ticks[i]!, xM: at(i, P_X_M) };
      }
      const q = at(i, P_Q_RAD_S);
      maxAbsQ = Math.max(maxAbsQ, Math.abs(q));
      sumQ2 += q * q;
      maxAbsTheta = Math.max(maxAbsTheta, Math.abs(at(i, P_THETA_RAD)));
      const airspeed = Math.hypot(at(i, P_U_MPS), at(i, P_W_MPS));
      minAirspeed = minAirspeed === null ? airspeed : Math.min(minAirspeed, airspeed);
      const sign = q > 1e-3 ? 1 : q < -1e-3 ? -1 : 0;
      if (sign !== 0 && lastSign !== 0 && sign !== lastSign) {
        flips += 1;
      }
      if (sign !== 0) {
        lastSign = sign;
      }
    }
  }
  const lastPhaseWord = PHASE_BY_CODE[at(n - 1, P_PHASE)] ?? rec.terminalPhase;
  return {
    frames: n,
    terminal: lastPhaseWord.startsWith("ended:") ? lastPhaseWord : rec.terminalPhase,
    downrangeM: at(n - 1, P_X_M),
    pathLengthM: pathLength,
    liftoff,
    airborneS: airborneFrames / 120,
    undulations: Math.floor(flips / 2),
    maxAbsQRadS: maxAbsQ,
    rmsQRadS: airborneFrames > 0 ? Math.sqrt(sumQ2 / airborneFrames) : 0,
    maxAbsThetaRad: maxAbsTheta,
    maxHM: maxH,
    canardTravelRad: canardTravel,
    minAirspeedMps: minAirspeed,
  };
}

/** The Dec-17 historical context DISTRIBUTION (both lineages recorded;
 * flight 4 is the only precisely measured one — Orville's diary). */
export const DEC17_CONTEXT = [
  { flight: 1, note: "~12 s, canonical 120 ft (~37 m); diary lineage ~100 ft from track end" },
  { flight: 2, note: "~11-12 s, ~175 ft (~53 m) class (accounts vary)" },
  { flight: 3, note: "~15 s, ~200 ft (~61 m) class (accounts vary)" },
  { flight: 4, note: "59 s, 852 ft (~260 m) over ground — the only precisely measured flight" },
] as const;

/** Render the card lines (formatting only — every number comes from
 * the KPIs object; the battery's hostile twin recomputes them). */
export function cardLines(k: FlightKpis, siteWord: string): string[] {
  const lines = [
    `RESULTS — ${siteWord} (${k.frames} sim snapshots)`,
    `terminal: ${k.terminal}`,
    `downrange ${k.downrangeM.toFixed(1)} m | path length ${k.pathLengthM.toFixed(1)} m`,
    k.liftoff !== null
      ? `liftoff at tick ${k.liftoff.tick} (${k.liftoff.xM.toFixed(1)} m); airborne ${k.airborneS.toFixed(1)} s`
      : "never left the rail",
    `undulations ${k.undulations} | max |q| ${k.maxAbsQRadS.toFixed(2)} rad/s | rms q ${k.rmsQRadS.toFixed(2)} rad/s`,
    `max |theta| ${k.maxAbsThetaRad.toFixed(2)} rad | peak height ${k.maxHM.toFixed(1)} m`,
    `canard travel ${k.canardTravelRad.toFixed(2)} rad${
      k.minAirspeedMps !== null ? ` | min airspeed ${k.minAirspeedMps.toFixed(1)} m/s` : ""
    }`,
    "— Dec 17 1903 context (distributions, not targets) —",
    ...DEC17_CONTEXT.map((c) => `flight ${c.flight}: ${c.note}`),
  ];
  return lines;
}

/**
 * The KPI-vs-recompute gate: a card's KPI object must be bitwise equal
 * to an independent recompute from the same transcript. Returns the
 * first divergent field or null (the hostile twin proves it can fire).
 */
export function kpiRecomputeDivergence(rec: FlightRecording, card: FlightKpis): string | null {
  const fresh = computeKpis(rec);
  for (const key of Object.keys(fresh) as (keyof FlightKpis)[]) {
    const a = JSON.stringify(fresh[key]);
    const b = JSON.stringify(card[key]);
    if (a !== b) {
      return `${key}: recompute ${a} != card ${b}`;
    }
  }
  return null;
}
