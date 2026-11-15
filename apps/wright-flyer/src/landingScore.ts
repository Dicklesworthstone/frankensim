// Landing quality scoring (B11). PURE: reads the frozen recording
// transcript only — no DOM, no sim imports beyond the protocol layout
// constants and the replay type. Headless-tested.
//
// HONESTY LAW (missions/flights.ts idiom): the verdict describes THIS
// modeled run's touchdown against declared sink-rate bands (no source is
// cited for them); it never claims a historical comparison or a scoreboard.
// Repro: node --test test/landingScore.test.ts

import { P_THETA_RAD, P_U_MPS, P_W_MPS } from "./sim/protocol.ts";
import { type FlightRecording, recordingPayloadWords } from "./sim/replay.ts";

/** Vertical speed [m/s, negative = descending] at the FINAL recorded
 * frame (ground contact is terminal: ended:ground-contact). Null when
 * the transcript is malformed or empty — never an invented number.
 *
 * The frame's u and w are body-axis (FRD) components, so the vertical
 * speed is sin(theta) u - cos(theta) w: the inverse of fs-flyer's
 * `body_air_velocity` (simloop/motion.rs), and the reconstruction its
 * simloop_battery uses for ground vertical velocity. The w slot alone
 * is not the sink rate once the aircraft is pitched. */
export function touchdownVerticalSpeed(rec: FlightRecording): number | null {
  const payloadWords = recordingPayloadWords(rec);
  const n = rec.frames.length / payloadWords;
  if (!Number.isInteger(n) || n < 1) {
    return null;
  }
  const base = (n - 1) * payloadWords;
  const u = rec.frames[base + P_U_MPS];
  const w = rec.frames[base + P_W_MPS];
  const theta = rec.frames[base + P_THETA_RAD];
  if (
    typeof u !== "number" ||
    typeof w !== "number" ||
    typeof theta !== "number" ||
    !Number.isFinite(u) ||
    !Number.isFinite(w) ||
    !Number.isFinite(theta)
  ) {
    return null;
  }
  return Math.sin(theta) * u - Math.cos(theta) * w;
}

export type TouchdownGrade = "buttery" | "firm" | "hard" | "unlogged";

export interface TouchdownReport {
  readonly grade: TouchdownGrade;
  /** One card-ready honest line (no outcome promises). */
  readonly line: string;
}

/** Declared sink-rate bands for a skid-equipped machine on flat sand
 * (no source is cited for them): under 1.2 m/s the skids barely notice;
 * past 2.5 m/s the struts take a real knock. Bands are MODEL judgments
 * about THIS modeled run, labeled as such in the line. */
export function scoreTouchdown(rec: FlightRecording): TouchdownReport | null {
  if (rec.terminalPhase !== "ended:ground-contact") {
    return null;
  }
  const w = touchdownVerticalSpeed(rec);
  if (w === null) {
    return { grade: "unlogged", line: "TOUCHDOWN UNLOGGED — transcript missing sink rate" };
  }
  const sink = Math.abs(w);
  if (sink <= 1.2) {
    return {
      grade: "buttery",
      line: `TOUCHDOWN BUTTERY — sink ${sink.toFixed(2)} m/s; the skids settled like a gull landing (modeled run, not history)`,
    };
  }
  if (sink <= 2.5) {
    return {
      grade: "firm",
      line: `TOUCHDOWN FIRM — sink ${sink.toFixed(2)} m/s; a solid arrival on the skids (modeled run, not history)`,
    };
  }
  return {
    grade: "hard",
    line: `TOUCHDOWN HARD — sink ${sink.toFixed(2)} m/s; the structure took the knock (modeled run, not history)`,
  };
}
