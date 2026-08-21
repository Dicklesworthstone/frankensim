// Sim-worker protocol v1 (bead wf-root-guzez.6.3.1, E5.2a): the typed
// seam between the main thread and the worker that drives the REAL
// fs-flyer-wasm engine. Snapshot payloads ride the E0.7 seqlock ring
// (SharedArrayBuffer) when available, else a postMessage fallback —
// EITHER WAY the payload is the frozen 12-float layout the native
// engine digests (fs_flyer::simloop::SNAPSHOT_LEN order, exactly).

/** Frozen 12-float payload layout (mirror of simloop::snapshot_payload). */
export const PAYLOAD_F64S = 12;
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

/** Phase words (engine envelope) ↔ payload codes (engine digest layout). */
export const PHASE_CODES = {
  "on-rail": 0,
  airborne: 1,
  "ended:ground-contact": 2,
  "ended:rail-end-without-lift": 3,
  "ended:max-ticks": 4,
  "ended:envelope-exceeded": 5,
} as const;
export type PhaseWord = keyof typeof PHASE_CODES;

/**
 * Deterministic layout identity for the seqlock header (i32). Derived
 * from the layout descriptor string — any reordering or resize of the
 * payload changes it, so a stale reader refuses before touching floats.
 */
export const PAYLOAD_LAYOUT_V1 =
  "wf-snapshot-v1:x_m,h_m,u_mps,w_mps,q_rad_s,theta_rad,dc_rad,warp_rad,omega_prop_rad_s,gust_w_mps,assist,phase";
export function payloadLayoutHash(descriptor: string = PAYLOAD_LAYOUT_V1): number {
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

export interface ScenarioInit {
  readonly seed: bigint;
  readonly rhoKgM3: number;
  readonly headwindMps: number;
  readonly mode: number;
  readonly member: number;
  readonly railLengthM: number;
  readonly maxTicks: bigint;
}

/** The Dec-17 reference scenario (mirror of dec17_scenario). */
export function dec17Scenario(seed: bigint, mode: number, member = 0): ScenarioInit {
  return {
    seed,
    rhoKgM3: 1.294,
    headwindMps: 11.0,
    mode,
    member,
    railLengthM: 18.3,
    maxTicks: 2400n,
  };
}

// --- main → worker ----------------------------------------------------------
export type MainToWorker =
  | {
      readonly kind: "init";
      readonly scenario: ScenarioInit;
      /** Present iff SharedArrayBuffer transport is available. */
      readonly sab?: SharedArrayBuffer;
      readonly slots?: number;
      readonly runEpoch: number;
    }
  | {
      readonly kind: "control";
      readonly leverForceN: number;
      readonly warpCmdRad: number;
      /** E5.3a: identity of this device sample in the latency ledger. */
      readonly sequence: number;
      /** Device timestamp TRANSLATED into the worker's monotonic clock
       * (main applies the ping-exchange offset before sending). */
      readonly deviceWorkerMs: number;
    }
  | { readonly kind: "ping"; readonly nonce: number; readonly localSentMs: number }
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
    }
  | { readonly kind: "refusal"; readonly stage: "init" | "step"; readonly refusal: RefusalEnvelope }
  | {
      readonly kind: "terminal";
      readonly phase: PhaseWord;
      readonly tick: number;
      readonly envelopeRefusalCode?: string;
      readonly digest: string;
    }
  | {
      /** postMessage fallback transport only (no SAB): one snapshot. */
      readonly kind: "snapshot";
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
      readonly nonce: number;
      readonly localSentMs: number;
      readonly remoteMs: number;
    }
  | {
      /** E5.3a ApplyNextEligibleTickAndFlag receipt for one control. */
      readonly kind: "control-ack";
      readonly sequence: number;
      readonly appliedTick: number;
      readonly lateByTicks: number;
    };
