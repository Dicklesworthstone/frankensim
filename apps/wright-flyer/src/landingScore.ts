// Landing quality scoring (B11). PURE: reads the frozen recording
// transcript only — no DOM, no sim imports beyond the protocol layout
// constants and the replay type. Headless-tested.
//
// HONESTY LAW (missions/flights.ts idiom): the verdict describes THIS
// modeled run's touchdown against physics-derived bands; it never
// claims a historical comparison or a scoreboard.
// Repro: node --test test/landingScore.test.ts

import { P_W_MPS } from "./sim/protocol.ts";
import { recordingPayloadWords, type FlightRecording } from "./sim/replay.ts";

/** Vertical speed [m/s, negative = descending] at the FINAL recorded
 * frame (ground contact is terminal: ended:ground-contact). Null when
 * the transcript is malformed or empty — never an invented number. */
export function touchdownVerticalSpeed(rec: FlightRecording): number | null {
  const payloadWords = recordingPayloadWords(rec);
  const n = rec.frames.length / payloadWords;
  if (!Number.isInteger(n) || n < 1) {
    return null;
  }
  const w = rec.frames[(n - 1) * payloadWords + P_W_MPS];
  return typeof w === "number" && Number.isFinite(w) ? w : null;
}

export type TouchdownGrade = "buttery" | "firm" | "hard" | "unlogged";

export interface TouchdownReport {
  readonly grade: TouchdownGrade;
  /** One card-ready honest line (no outcome promises). */
  readonly line: string;
}

/** Physics-derived sink-rate bands for a skid-equipped machine on
 * flat sand: under ~1.2 m/s the skids barely notice; past ~2.5 m/s
 * the struts take a real knock. Bands are MODEL judgments about THIS
 * reconstruction, labeled as such in the line. */
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
