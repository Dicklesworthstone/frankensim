// Landing score battery (B11): touchdown extraction from the frozen
// transcript layout and the honest grading bands.
// Repro: node --test test/landingScore.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";
import { PAYLOAD_F64S, PAYLOAD_F64S_V1, P_W_MPS } from "../src/sim/protocol.ts";
import {
  RECORDING_SCHEMA,
  RECORDING_SCHEMA_V1,
  type FlightRecording,
} from "../src/sim/replay.ts";
import { scoreTouchdown, touchdownVerticalSpeed } from "../src/landingScore.ts";

function recWith(
  terminalPhase: string,
  lastWMps: number,
  schema: typeof RECORDING_SCHEMA | typeof RECORDING_SCHEMA_V1 = RECORDING_SCHEMA,
): FlightRecording {
  const payloadWords = schema === RECORDING_SCHEMA_V1 ? PAYLOAD_F64S_V1 : PAYLOAD_F64S;
  const frames: number[] = [];
  const push = (w: number): void => {
    for (let i = 0; i < payloadWords; i += 1) {
      frames.push(i === P_W_MPS ? w : 0);
    }
  };
  push(-3);
  push(-1.0);
  push(lastWMps);
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

test("touchdown vertical speed reads the final frame's w slot", () => {
  const r = recWith("ended:ground-contact", -1.8);
  assert.ok(Math.abs(touchdownVerticalSpeed(r)! + 1.8) < 1e-12);
});

test("legacy v1 touchdown uses the 12-word frame stride", () => {
  const r = recWith("ended:ground-contact", -1.8, RECORDING_SCHEMA_V1);
  assert.ok(Math.abs(touchdownVerticalSpeed(r)! + 1.8) < 1e-12);
});

test("grades follow the sink bands and stay honest", () => {
  const buttery = scoreTouchdown(recWith("ended:ground-contact", -1.1))!;
  assert.equal(buttery.grade, "buttery");
  assert.match(buttery.line, /modeled run, not history/);
  const firm = scoreTouchdown(recWith("ended:ground-contact", -2.4))!;
  assert.equal(firm.grade, "firm");
  const hard = scoreTouchdown(recWith("ended:ground-contact", -4.5))!;
  assert.equal(hard.grade, "hard");
});

test("non-landing terminals and empty transcripts refuse", () => {
  assert.equal(scoreTouchdown(recWith("ended:max-ticks", -1)), null);
  const bad = { ...recWith("ended:ground-contact", -1), frames: [] };
  assert.equal(scoreTouchdown(bad as unknown as FlightRecording)!.grade, "unlogged");
});
