// E5.2b snapshot-view core (bead wf-root-guzez.6.3.2): PURE mapping
// from the frozen engine payload to what the presentation
// plane consumes — typed snapshot, render-rate interpolation (sim 120
// Hz, render decoupled), pose/HUD inputs. No three.js, no DOM: every
// branch is headless-tested. The scene applies; this decides.

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
  PAYLOAD_F64S_V1,
  PHASE_CODES,
  type PhaseWord,
} from "./protocol.ts";

export interface SimSnapshot {
  readonly tick: number;
  readonly phase: PhaseWord;
  readonly ended: boolean;
  readonly xM: number;
  readonly hM: number;
  readonly uMps: number;
  readonly wMps: number;
  readonly qRadS: number;
  readonly thetaRad: number;
  readonly phiRad: number;
  readonly psiRad: number;
  readonly dcRad: number;
  readonly warpRad: number;
  readonly omegaPropRadS: number;
  readonly gustWMps: number;
  readonly assistActive: boolean;
}

const CODE_TO_PHASE: readonly PhaseWord[] = (() => {
  const words = Object.keys(PHASE_CODES) as PhaseWord[];
  const byCode: PhaseWord[] = [];
  for (const w of words) {
    byCode[PHASE_CODES[w]] = w;
  }
  return byCode;
})();

/** Decode a payload. Stored v1 records receive an explicit zero-lateral
 * fallback; a live v2 ring is hash-negotiated and never silently falls
 * back. Bad lengths and phase codes fail closed. */
export function decodeSnapshot(tick: number, payload: Float64Array): SimSnapshot {
  if (payload.length !== PAYLOAD_F64S_V1 && payload.length < PAYLOAD_F64S) {
    throw new RangeError(
      `payload ${payload.length} is neither v1 (${PAYLOAD_F64S_V1}) nor v2 (${PAYLOAD_F64S})`,
    );
  }
  const phase = CODE_TO_PHASE[payload[P_PHASE]!];
  if (phase === undefined) {
    throw new RangeError(`unknown phase code ${payload[P_PHASE]}`);
  }
  return {
    tick,
    phase,
    ended: phase.startsWith("ended:"),
    xM: payload[P_X_M]!,
    hM: payload[P_H_M]!,
    uMps: payload[P_U_MPS]!,
    wMps: payload[P_W_MPS]!,
    qRadS: payload[P_Q_RAD_S]!,
    thetaRad: payload[P_THETA_RAD]!,
    phiRad: payload.length >= PAYLOAD_F64S ? payload[P_PHI_RAD]! : 0,
    psiRad: payload.length >= PAYLOAD_F64S ? payload[P_PSI_RAD]! : 0,
    dcRad: payload[P_DC_RAD]!,
    warpRad: payload[P_WARP_RAD]!,
    omegaPropRadS: payload[P_OMEGA_RAD_S]!,
    gustWMps: payload[P_GUST_W_MPS]!,
    assistActive: payload[P_ASSIST]! !== 0,
  };
}

/**
 * Render-rate interpolation between two sim snapshots. Continuous
 * fields lerp; DISCRETE fields (phase, assist) HOLD the older value —
 * and a phase boundary is never crossed early: if the phases differ,
 * all fields hold `a` until alpha reaches 1 (the render plane may not
 * announce a landing the sim hasn't reached yet at its display time).
 */
export function interpolateSnapshots(a: SimSnapshot, b: SimSnapshot, alpha: number): SimSnapshot {
  const t = Math.min(1, Math.max(0, alpha));
  if (a.phase !== b.phase) {
    return t >= 1 ? b : a;
  }
  const lerp = (x: number, y: number): number => x + (y - x) * t;
  return {
    tick: t >= 1 ? b.tick : a.tick,
    phase: a.phase,
    ended: a.ended,
    xM: lerp(a.xM, b.xM),
    hM: lerp(a.hM, b.hM),
    uMps: lerp(a.uMps, b.uMps),
    wMps: lerp(a.wMps, b.wMps),
    qRadS: lerp(a.qRadS, b.qRadS),
    thetaRad: lerp(a.thetaRad, b.thetaRad),
    phiRad: lerp(a.phiRad, b.phiRad),
    psiRad: lerp(a.psiRad, b.psiRad),
    dcRad: lerp(a.dcRad, b.dcRad),
    warpRad: lerp(a.warpRad, b.warpRad),
    omegaPropRadS: lerp(a.omegaPropRadS, b.omegaPropRadS),
    gustWMps: lerp(a.gustWMps, b.gustWMps),
    assistActive: a.assistActive,
  };
}

const RAD_TO_DEG = 180 / Math.PI;
/** Prop chain ratio: engine 23t drives the 8t prop sprocket... inverse —
 * the ENGINE turns 1025 rpm through 23:8 down to the prop; engine rpm =
 * prop omega [rad/s] * 60/(2π) * 23/8. */
const ENGINE_PER_PROP = 23 / 8;

export interface SimDriveState {
  /** Integrated prop shaft angle [rad] (payload carries omega, not angle). */
  propAngleRad: number;
}

/** Advance the integrated prop angle by one render frame. */
export function advanceProp(state: SimDriveState, snap: SimSnapshot, dtS: number): SimDriveState {
  if (!Number.isFinite(dtS) || dtS < 0 || dtS > 1) {
    throw new RangeError(`render dt out of domain: ${dtS}`);
  }
  return { propAngleRad: state.propAngleRad + snap.omegaPropRadS * dtS };
}

/** The rig control-state for computePose (physics → visual, 1:1 units). */
export function controlStateFrom(snap: SimSnapshot, drive: SimDriveState): {
  canardDeg: number;
  warpDeg: number;
  rudderDeg: number;
  coupled: boolean;
  propAngleRad: number;
} {
  return {
    canardDeg: snap.dcRad * RAD_TO_DEG,
    warpDeg: snap.warpRad * RAD_TO_DEG,
    rudderDeg: 0,
    coupled: true, // 1903 slaved-rudder wiring
    propAngleRad: drive.propAngleRad,
  };
}

/** World transform for the airframe group: flight line = scene +x,
 * up = scene +y. Pitch, bank, and heading come from the simulation;
 * this layer never infers attitude from the warp control. */
export function worldTransformFrom(
  snap: SimSnapshot,
  launch: readonly [number, number, number],
): { position: [number, number, number]; pitchRad: number; rollRad: number; headingRad: number } {
  return {
    position: [launch[0] + snap.xM, launch[1] + snap.hM, launch[2]],
    pitchRad: snap.thetaRad,
    rollRad: snap.phiRad,
    headingRad: snap.psiRad,
  };
}

/** HUD inputs from a snapshot (the hudLines formatter consumes these). */
export function hudInputsFrom(snap: SimSnapshot): {
  airspeedMps: number;
  elapsedS: number;
  engineRpm: number;
  phase: PhaseWord;
} {
  return {
    airspeedMps: Math.hypot(snap.uMps, snap.wMps),
    elapsedS: snap.tick / 120,
    engineRpm: snap.omegaPropRadS * (60 / (2 * Math.PI)) * ENGINE_PER_PROP,
    phase: snap.phase,
  };
}

/** Phase banner line (terminal states surface loudly, never silently). */
export function phaseBanner(snap: SimSnapshot, envelopeRefusalCode?: string): string | null {
  switch (snap.phase) {
    case "on-rail":
      return null;
    case "airborne":
      return null;
    case "ended:ground-contact":
      return "LANDED — ground contact";
    case "ended:rail-end-without-lift":
      return "RAN OFF THE RAIL — no lift";
    case "ended:max-ticks":
      return "TIME LIMIT — still flying";
    case "ended:envelope-exceeded":
      return `FLIGHT LEFT THE CERTIFIED ENVELOPE${
        envelopeRefusalCode !== undefined ? ` (${envelopeRefusalCode})` : ""
      }`;
  }
}
