// Engine-envelope facade (bead wf-root-guzez.6.3.1, E5.2a): PURE
// functions between the wasm engine's JSON envelopes and the typed
// protocol — no DOM, no Worker, no wasm import — so every branch is
// headless-testable in node, and the worker entry stays thin.

import {
  P_ASSIST,
  P_DC_RAD,
  P_GUST_W_MPS,
  P_H_M,
  P_OMEGA_RAD_S,
  P_PHASE,
  P_PHI_RAD,
  P_PSI_RAD,
  P_Q_RAD_S,
  P_THETA_RAD,
  P_U_MPS,
  P_W_MPS,
  P_WARP_RAD,
  P_X_M,
  PAYLOAD_F64S,
  PHASE_CODES,
  type PhaseWord,
  type RefusalEnvelope,
} from "./protocol.ts";

export interface EngineInitOk {
  readonly kind: "ok";
  readonly runIntentId: string;
  readonly tick0Digest: string;
  readonly trimVMps: number;
}
export interface EngineStepOk {
  readonly kind: "ok";
  readonly tick: number;
  readonly phase: PhaseWord;
  readonly ended: boolean;
  readonly envelopeRefusalCode?: string;
  readonly state: {
    readonly xM: number;
    readonly hM: number;
    readonly uMps: number;
    readonly wMps: number;
    readonly qRadS: number;
    readonly thetaRad: number;
    readonly pRadS: number;
    readonly phiRad: number;
    readonly rRadS: number;
    readonly psiRad: number;
    readonly dcRad: number;
    readonly warpRad: number;
    readonly omegaPropRadS: number;
    readonly gustWMps: number;
    readonly assistActive: boolean;
  };
}
export interface EngineRefusal {
  readonly kind: "refusal";
  readonly refusal: RefusalEnvelope;
}
export interface EngineMalformed {
  readonly kind: "malformed";
  readonly detail: string;
}

function malformed(detail: string): EngineMalformed {
  return { kind: "malformed", detail };
}

function parseEnvelope(json: string): { ok?: unknown; refusal?: unknown } | EngineMalformed {
  try {
    const value: unknown = JSON.parse(json);
    if (typeof value !== "object" || value === null) {
      return malformed("envelope is not an object");
    }
    return value as { ok?: unknown; refusal?: unknown };
  } catch (e) {
    return malformed(`envelope is not JSON: ${String(e)}`);
  }
}

function asRefusal(value: unknown): RefusalEnvelope | null {
  if (typeof value !== "object" || value === null) {
    return null;
  }
  const r = value as Record<string, unknown>;
  if (typeof r.code !== "string" || typeof r.message !== "string" || !Array.isArray(r.ranked_repairs)) {
    return null;
  }
  return {
    code: r.code,
    message: r.message,
    ranked_repairs: r.ranked_repairs.filter((x): x is string => typeof x === "string"),
  };
}

/** Parse the init envelope (fail-closed: unknown shape is `malformed`). */
export function parseInitEnvelope(json: string): EngineInitOk | EngineRefusal | EngineMalformed {
  const env = parseEnvelope(json);
  if ("kind" in env) {
    return env;
  }
  if (env.refusal !== undefined) {
    const refusal = asRefusal(env.refusal);
    return refusal ? { kind: "refusal", refusal } : malformed("refusal shape invalid");
  }
  const ok = env.ok as Record<string, unknown> | undefined;
  if (
    ok === undefined ||
    typeof ok.run_intent_id !== "string" ||
    typeof ok.tick0_digest !== "string" ||
    typeof ok.trim_v_mps !== "number"
  ) {
    return malformed("init ok-envelope missing identity fields");
  }
  return {
    kind: "ok",
    runIntentId: ok.run_intent_id,
    tick0Digest: ok.tick0_digest,
    trimVMps: ok.trim_v_mps,
  };
}

const STEP_NUMERIC_KEYS = [
  "tick",
  "x_m",
  "h_m",
  "u_mps",
  "w_mps",
  "q_rad_s",
  "theta_rad",
  "p_rad_s",
  "phi_rad",
  "r_rad_s",
  "psi_rad",
  "dc_rad",
  "warp_rad",
  "omega_prop_rad_s",
  "gust_w_mps",
] as const;

/** Parse the step envelope (fail-closed on any missing/wrong field). */
export function parseStepEnvelope(json: string): EngineStepOk | EngineRefusal | EngineMalformed {
  const env = parseEnvelope(json);
  if ("kind" in env) {
    return env;
  }
  if (env.refusal !== undefined) {
    const refusal = asRefusal(env.refusal);
    return refusal ? { kind: "refusal", refusal } : malformed("refusal shape invalid");
  }
  const ok = env.ok as Record<string, unknown> | undefined;
  if (ok === undefined) {
    return malformed("neither ok nor refusal");
  }
  for (const key of STEP_NUMERIC_KEYS) {
    if (typeof ok[key] !== "number" || !Number.isFinite(ok[key] as number)) {
      return malformed(`step field ${key} missing or non-finite`);
    }
  }
  const phase = ok.phase;
  if (typeof phase !== "string" || !(phase in PHASE_CODES)) {
    return malformed(`unknown phase word ${String(phase)}`);
  }
  if (typeof ok.assist_active !== "boolean") {
    return malformed("assist_active missing");
  }
  const envelopeCode = ok.envelope_refusal_code;
  if (envelopeCode !== undefined && typeof envelopeCode !== "string") {
    return malformed("envelope_refusal_code wrong type");
  }
  const phaseWord = phase as PhaseWord;
  return {
    kind: "ok",
    tick: ok.tick as number,
    phase: phaseWord,
    ended: phaseWord.startsWith("ended:"),
    ...(envelopeCode !== undefined ? { envelopeRefusalCode: envelopeCode } : {}),
    state: {
      xM: ok.x_m as number,
      hM: ok.h_m as number,
      uMps: ok.u_mps as number,
      wMps: ok.w_mps as number,
      qRadS: ok.q_rad_s as number,
      thetaRad: ok.theta_rad as number,
      pRadS: ok.p_rad_s as number,
      phiRad: ok.phi_rad as number,
      rRadS: ok.r_rad_s as number,
      psiRad: ok.psi_rad as number,
      dcRad: ok.dc_rad as number,
      warpRad: ok.warp_rad as number,
      omegaPropRadS: ok.omega_prop_rad_s as number,
      gustWMps: ok.gust_w_mps as number,
      assistActive: ok.assist_active as boolean,
    },
  };
}

/**
 * Assemble the frozen v2 ring payload from a parsed step (the
 * EXACT order the native engine digests — per-field, never a spread).
 */
export function fillPayload(step: EngineStepOk, out: Float64Array): void {
  if (out.length < PAYLOAD_F64S) {
    throw new Error(`payload buffer ${out.length} < ${PAYLOAD_F64S}`);
  }
  out[P_X_M] = step.state.xM;
  out[P_H_M] = step.state.hM;
  out[P_U_MPS] = step.state.uMps;
  out[P_W_MPS] = step.state.wMps;
  out[P_Q_RAD_S] = step.state.qRadS;
  out[P_THETA_RAD] = step.state.thetaRad;
  out[P_DC_RAD] = step.state.dcRad;
  out[P_WARP_RAD] = step.state.warpRad;
  out[P_OMEGA_RAD_S] = step.state.omegaPropRadS;
  out[P_GUST_W_MPS] = step.state.gustWMps;
  out[P_ASSIST] = step.state.assistActive ? 1 : 0;
  out[P_PHASE] = PHASE_CODES[step.phase];
  out[P_PHI_RAD] = step.state.phiRad;
  out[P_PSI_RAD] = step.state.psiRad;
}

/** Parse the digest envelope → 64-hex string or a typed failure. */
export function parseDigestEnvelope(json: string): string | EngineRefusal | EngineMalformed {
  const env = parseEnvelope(json);
  if ("kind" in env) {
    return env;
  }
  if (env.refusal !== undefined) {
    const refusal = asRefusal(env.refusal);
    return refusal ? { kind: "refusal", refusal } : malformed("refusal shape invalid");
  }
  const ok = env.ok as Record<string, unknown> | undefined;
  const digest = ok?.digest;
  if (typeof digest !== "string" || !/^[0-9a-f]{64}$/.test(digest)) {
    return malformed("digest missing or not 64-hex");
  }
  return digest;
}
