// Sim-worker protocol v1 (bead wf-root-guzez.6.3.1, E5.2a): the typed
// seam between the main thread and the worker that drives the REAL
// fs-flyer-wasm engine. Snapshot payloads ride the E0.7 seqlock ring
// (SharedArrayBuffer) when available, else a postMessage fallback —
// EITHER WAY the payload is the frozen v2 layout the native engine
// digests (fs_flyer::simloop::SNAPSHOT_LEN order, exactly).

import type { AppliedControlTrace } from "./appliedControls.ts";

/** Frozen v1 prefix length; retained for stored-recording fallback. */
export const PAYLOAD_F64S_V1 = 12;
/** Frozen v2 payload length (v1 prefix + roll + heading). */
export const PAYLOAD_F64S = 14;
export const P_X_M = 0;
export const P_H_M = 1;
export const P_U_MPS = 2;
export const P_W_MPS = 3;
export const P_Q_RAD_S = 4;
export const P_THETA_RAD = 5;
export const P_DC_RAD = 6;
export const P_WARP_RAD = 7;
export const P_OMEGA_RAD_S = 8;
export const P_GUST_W_MPS = 9;
export const P_ASSIST = 10;
export const P_PHASE = 11;
export const P_PHI_RAD = 12;
export const P_PSI_RAD = 13;

/** Phase words (engine envelope) ↔ payload codes (engine digest layout). */
export const PHASE_CODES = {
  "on-rail": 0,
  airborne: 1,
  "ended:ground-contact": 2,
  "ended:rail-end-without-lift": 3,
  "ended:max-ticks": 4,
  "ended:envelope-exceeded": 5,
  "ended:damage-model-unavailable": 6,
} as const;
export type PhaseWord = keyof typeof PHASE_CODES;

/**
 * Deterministic layout identity for the seqlock header (i32). Derived
 * from the layout descriptor string — any reordering or resize of the
 * payload changes it, so a stale reader refuses before touching floats.
 */
export const PAYLOAD_LAYOUT_V1 =
  "wf-snapshot-v1:x_m,h_m,u_mps,w_mps,q_rad_s,theta_rad,dc_rad,warp_rad,omega_prop_rad_s,gust_w_mps,assist,phase";
export const PAYLOAD_LAYOUT_V2 =
  "wf-snapshot-v2:x_m,h_m,u_mps,w_mps,q_rad_s,theta_rad,dc_rad,warp_rad,omega_prop_rad_s,gust_w_mps,assist,phase,phi_rad,psi_rad";
export function payloadLayoutHash(descriptor: string = PAYLOAD_LAYOUT_V2): number {
  // FNV-1a 32-bit — tiny, deterministic, dependency-free (identity tag,
  // not cryptographic; the engine's blake3 digest carries the truth).
  let h = 0x811c9dc5;
  for (let i = 0; i < descriptor.length; i += 1) {
    h ^= descriptor.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return h | 0;
}

/** Pilot-mode words at the wasm ABI (mirror of engine::MODE_*). */
export const MODE_FIXED = 0;
export const MODE_HISTORICAL = 1;
export const MODE_HUMAN = 2;

/** Explicit synthetic atmosphere for the actual force/propeller solver.
 * The scenario headwind supplies the log-law reference speed. Visualization
 * fields remain separate; this configuration controls physics and replay. */
export interface SpatialWindInit {
  readonly sampling: "cg" | "distributed";
  readonly roughnessLengthM: number;
  readonly displacementHeightM: number;
  readonly referenceHeightM: number;
  /** Spatial strength parameter, not a certified RMS gust speed. */
  readonly sigmaMps: number;
  readonly lengthM: number;
  readonly advectionMps: number;
  /** One positive finite time per spatial mode (1..=256). */
  readonly correlationTimesS: readonly number[];
}

/** Validate serialized/worker input before allocating a WASM mode buffer.
 * Native atmosphere admission remains the authority for the actual solve. */
export function spatialWindError(value: unknown): string | null {
  if (typeof value !== "object" || value === null || Array.isArray(value)) {
    return "spatialWind must be an explicit parameter object";
  }
  const wind = value as Record<string, unknown>;
  if (wind.sampling !== "cg" && wind.sampling !== "distributed") {
    return "spatialWind sampling must be cg or distributed";
  }
  for (const key of [
    "roughnessLengthM", "displacementHeightM", "referenceHeightM",
    "sigmaMps", "lengthM", "advectionMps",
  ]) {
    if (typeof wind[key] !== "number" || !Number.isFinite(wind[key])) {
      return `spatialWind ${key} must be finite`;
    }
  }
  const w = value as SpatialWindInit;
  if (w.roughnessLengthM < 1e-5 || w.roughnessLengthM > 1) {
    return "spatialWind roughnessLengthM must lie in [0.00001, 1]";
  }
  if (w.displacementHeightM < 0 || w.displacementHeightM > 1 ||
      w.referenceHeightM <= w.displacementHeightM + w.roughnessLengthM) {
    return "spatialWind reference height must exceed displacement plus roughness";
  }
  if (w.sigmaMps < 0 || w.lengthM <= 0) {
    return "spatialWind requires nonnegative strength and positive length";
  }
  if (!Array.isArray(w.correlationTimesS) || w.correlationTimesS.length === 0 ||
      w.correlationTimesS.length > 256 ||
      w.correlationTimesS.some((time) => typeof time !== "number" || !Number.isFinite(time) || time <= 0)) {
    return "spatialWind requires 1..256 positive finite correlation times";
  }
  return null;
}

/** Prescribed homogeneous moist air, evaluated by the native material owner.
 * Its derived density replaces rhoKgM3 before trim. This requires spatialWind;
 * it does not introduce thermal evolution or Reynolds-dependent polars. */
export interface MaterialAirInit {
  readonly temperatureK: number;
  readonly pressurePa: number;
  /** Fraction [0, 1], relative to saturation over liquid water. */
  readonly relativeHumidity: number;
}

/** Basic wire admission. Gas-law, saturation and flight-envelope checks remain
 * with the native owners; the browser never computes a second air density. */
export function materialAirError(value: unknown): string | null {
  if (typeof value !== "object" || value === null || Array.isArray(value)) {
    return "materialAir must be an explicit parameter object";
  }
  const air = value as Record<string, unknown>;
  for (const key of ["temperatureK", "pressurePa", "relativeHumidity"]) {
    if (typeof air[key] !== "number" || !Number.isFinite(air[key])) {
      return `materialAir ${key} must be finite`;
    }
  }
  const state = value as MaterialAirInit;
  if (state.temperatureK <= 0 || state.pressurePa <= 0 ||
      state.relativeHumidity < 0 || state.relativeHumidity > 1) {
    return "materialAir requires positive temperature/pressure and humidity in [0, 1]";
  }
  return null;
}

export interface ScenarioInit {
  readonly seed: bigint;
  readonly rhoKgM3: number;
  readonly headwindMps: number;
  readonly mode: number;
  readonly member: number;
  readonly railLengthM: number;
  readonly maxTicks: bigint;
  /** E5.3c: ratified bounded assist (ASSIST_V1, authority 0.3). */
  readonly assist: boolean;
  /** E5.4: registered 1904 catapult tow (CATAPULT_1904_V1). */
  readonly catapult: boolean;
  /** Omitted selects the unchanged legacy scalar-OU flight model. */
  readonly spatialWind?: SpatialWindInit;
  /** When present, the shared gas law supplies density for all flight owners. */
  readonly materialAir?: MaterialAirInit;
}

/** The Dec-17 reference scenario (mirror of dec17_scenario). */
export function dec17Scenario(
  seed: bigint,
  mode: number,
  member = 0,
  assist = false,
): ScenarioInit {
  return {
    seed,
    rhoKgM3: 1.294,
    headwindMps: 11.0,
    mode,
    member,
    railLengthM: 18.3,
    maxTicks: 2400n,
    assist,
    catapult: false,
  };
}

/** The Huffman Prairie 1904-05 preset (mirror of huffman_scenario).
 * HONEST STATUS: the near-calm rail run is outside the coupled aero
 * core's admitted envelope — the engine ends it with a receipted
 * envelope exit within the first rail seconds (tracked guzez.5.7.1);
 * the site loads and the refusal surfaces on the HUD, never silently. */
export function huffmanScenario(seed: bigint, mode: number, member = 0): ScenarioInit {
  return {
    seed,
    rhoKgM3: 1.17,
    headwindMps: 2.0,
    mode,
    member,
    railLengthM: 30.0,
    maxTicks: 2400n,
    assist: false,
    catapult: true,
  };
}

// --- main → worker ----------------------------------------------------------
export interface ReplayInitialization {
  readonly runIntentId: string;
  readonly tick0Digest: string;
  readonly appliedControls?: AppliedControlTrace;
}

export type MainToWorker =
  | {
      readonly kind: "init";
      readonly scenario: ScenarioInit;
      readonly replay?: ReplayInitialization;
      /** Present iff SharedArrayBuffer transport is available. */
      readonly sab?: SharedArrayBuffer;
      readonly slots?: number;
      readonly runEpoch: number;
      /** Client-issued lifecycle generation; receipts echo it verbatim. */
      readonly initGeneration: number;
    }
  | {
      readonly kind: "control";
      /** Active run that this control is allowed to affect. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted the run. */
      readonly initGeneration: number;
      readonly leverForceN: number;
      readonly warpCmdRad: number;
      /** E5.3a: identity of this device sample in the latency ledger. */
      readonly sequence: number;
      /** Device timestamp TRANSLATED into the worker's monotonic clock
       * (main applies the ping-exchange offset before sending). */
      readonly deviceWorkerMs: number;
    }
  | {
      readonly kind: "ping";
      /** Active run requesting this clock sample. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted the run. */
      readonly initGeneration: number;
      readonly nonce: number;
      readonly localSentMs: number;
    }
  | {
      readonly kind: "checkpoint";
      /** Monotonic client-side correlation id for this exact request. */
      readonly requestId: number;
      /** Active run expected by the caller; the worker must not substitute. */
      readonly runIntentId: string;
    }
  | { readonly kind: "pause" }
  | { readonly kind: "resume" };

// --- worker → main ----------------------------------------------------------
export interface RefusalEnvelope {
  readonly code: string;
  readonly message: string;
  readonly ranked_repairs: readonly string[];
}

export type WorkerToMain =
  | {
      readonly kind: "ready";
      readonly runIntentId: string;
      readonly tick0Digest: string;
      readonly trimVMps: number;
      readonly layoutHash: number;
      /** Echo of the client-issued init generation. */
      readonly initGeneration: number;
    }
  | {
      readonly kind: "refusal";
      readonly stage: "init";
      /** Echo of the client-issued init generation. */
      readonly initGeneration: number;
      readonly refusal: RefusalEnvelope;
    }
  | {
      readonly kind: "refusal";
      readonly stage: "step";
      /** Active run that produced this step refusal. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted this run. */
      readonly initGeneration: number;
      readonly refusal: RefusalEnvelope;
    }
  | {
      /** Exact live SimLoop checkpoint bytes, transferable to persistence. */
      readonly kind: "checkpoint";
      readonly requestId: number;
      readonly runIntentId: string;
      readonly bytes: Uint8Array;
    }
  | {
      /** A checkpoint refusal bound to its request and expected run. */
      readonly kind: "checkpoint-refusal";
      readonly requestId: number;
      readonly runIntentId: string;
      readonly refusal: RefusalEnvelope;
    }
  | {
      readonly kind: "terminal";
      /** Active run that produced this terminal receipt. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted this run. */
      readonly initGeneration: number;
      readonly phase: PhaseWord;
      readonly tick: number;
      readonly envelopeRefusalCode?: string;
      readonly digest: string;
      /** Exact terminal state, including when no render frame sampled it. */
      readonly payload: Float64Array;
      readonly appliedControls?: AppliedControlTrace;
    }
  | {
      /** postMessage fallback transport only (no SAB): one snapshot. */
      readonly kind: "snapshot";
      /** Active run that produced this fallback snapshot. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted this run. */
      readonly initGeneration: number;
      readonly tick: number;
      readonly payload: Float64Array;
    }
  | {
      readonly kind: "metrics";
      readonly ticksRun: number;
      readonly reanchors: number;
      readonly maxBacklogObserved: number;
    }
  | {
      /** E5.3a clock-sync reply (estimateClockOffsetMs consumes these). */
      readonly kind: "pong";
      /** Active run that requested this clock sample. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted the run. */
      readonly initGeneration: number;
      readonly nonce: number;
      readonly localSentMs: number;
      readonly remoteMs: number;
    }
  | {
      /** E5.3a ApplyNextEligibleTickAndFlag receipt for one control. */
      readonly kind: "control-ack";
      /** Active run that admitted the control. */
      readonly runIntentId: string;
      /** Client lifecycle generation that admitted the run. */
      readonly initGeneration: number;
      readonly sequence: number;
      readonly appliedTick: number;
      readonly lateByTicks: number;
    };
