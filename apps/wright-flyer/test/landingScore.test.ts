// Landing score battery (B11): touchdown extraction from the frozen
// transcript layout and the honest grading bands.
// Repro: node --test test/landingScore.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";
import { scoreTouchdown, touchdownVerticalSpeed } from "../src/landingScore.ts";
import {
  P_THETA_RAD,
  P_U_MPS,
  P_W_MPS,
  PAYLOAD_F64S,
  PAYLOAD_F64S_V1,
} from "../src/sim/protocol.ts";
import { type FlightRecording, RECORDING_SCHEMA, RECORDING_SCHEMA_V1 } from "../src/sim/replay.ts";

/** A level frame (theta = 0, u = 0) whose only nonzero slot is the
 * body-axis w [m/s, FRD: positive down], unless `last` sets the final
 * frame's u, w and theta. */
function recWith(
  terminalPhase: string,
  lastWMps: number,
  schema: typeof RECORDING_SCHEMA | typeof RECORDING_SCHEMA_V1 = RECORDING_SCHEMA,
  last?: { uMps: number; wMps: number; thetaRad: number },
): FlightRecording {
  const payloadWords = schema === RECORDING_SCHEMA_V1 ? PAYLOAD_F64S_V1 : PAYLOAD_F64S;
  const frames: number[] = [];
  const push = (u: number, w: number, theta: number): void => {
    for (let i = 0; i < payloadWords; i += 1) {
      frames.push(i === P_U_MPS ? u : i === P_W_MPS ? w : i === P_THETA_RAD ? theta : 0);
    }
  };
  push(0, 3, 0);
  push(0, 1.0, 0);
  if (last === undefined) push(0, lastWMps, 0);
  else push(last.uMps, last.wMps, last.thetaRad);
  return {
    schema,
    scenario: {},
    runIntentId: "test",
    tick0Digest: "0".repeat(64),
    terminalPhase,
    finalDigest: "a".repeat(64),
    ticks: [1, 2, 3],
    frames,
  } as unknown as FlightRecording;
}

test("level touchdown: vertical speed is -w (FRD w is positive down)", () => {
  const r = recWith("ended:ground-contact", 1.8);
  assert.ok(Math.abs(touchdownVerticalSpeed(r)! + 1.8) < 1e-12);
});

test("legacy v1 touchdown uses the 12-word frame stride", () => {
  const r = recWith("ended:ground-contact", 1.8, RECORDING_SCHEMA_V1);
  assert.ok(Math.abs(touchdownVerticalSpeed(r)! + 1.8) < 1e-12);
});

test("pitched touchdown recovers the sink rate, not the w slot", () => {
  // Forward map of fs-flyer simloop/motion.rs body_air_velocity: a
  // horizontal 10 m/s and a 0.5 m/s descent seen at theta = 0.1 rad.
  const [vx, vh, theta] = [10, -0.5, 0.1];
  const [s, c] = [Math.sin(theta), Math.cos(theta)];
  const last = { uMps: c * vx + s * vh, wMps: s * vx - c * vh, thetaRad: theta };
  const r = recWith("ended:ground-contact", 0, RECORDING_SCHEMA, last);
  assert.ok(Math.abs(touchdownVerticalSpeed(r)! - vh) < 1e-12);
  // w alone (~1.49 m/s) would have graded this 0.5 m/s descent "firm".
  assert.ok(last.wMps > 1.2);
  assert.equal(scoreTouchdown(r)!.grade, "buttery");
});

test("grades follow the sink bands and stay honest", () => {
  const buttery = scoreTouchdown(recWith("ended:ground-contact", 1.1))!;
  assert.equal(buttery.grade, "buttery");
  assert.match(buttery.line, /modeled run, not history/);
  const firm = scoreTouchdown(recWith("ended:ground-contact", 2.4))!;
  assert.equal(firm.grade, "firm");
  const hard = scoreTouchdown(recWith("ended:ground-contact", 4.5))!;
  assert.equal(hard.grade, "hard");
});

test("non-landing terminals and empty transcripts refuse", () => {
  assert.equal(scoreTouchdown(recWith("ended:max-ticks", -1)), null);
  const bad = { ...recWith("ended:ground-contact", -1), frames: [] };
  assert.equal(scoreTouchdown(bad as unknown as FlightRecording)!.grade, "unlogged");
});
