// E5.2c replay battery (bead wf-root-guzez.6.3.3): recorder laws
// (monotone ticks, seal caps), fail-closed parse with hostile twins,
// ghost lookup oracles (before-start null, exact hit, hold-last), and
// the LIVE record→replay identity: the REAL engine runs the same
// Dec-17 scenario twice — digests equal — then a third run with a
// DIFFERENT seed proves the verdict actually discriminates (the
// falsifier: an identity check that cannot fail is not a check).
// Repro: WF_PKG=<pkg>/fs_flyer_wasm.js node --test test/replay.test.ts

import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { test } from "node:test";

import {
  FlightRecorder,
  RECORDING_SCHEMA,
  ghostAt,
  parseRecording,
  recordedToScenario,
  replayVerdict,
  scenarioToRecorded,
  type FlightRecording,
} from "../src/sim/replay.ts";
import { MODE_FIXED, PAYLOAD_F64S, dec17Scenario } from "../src/sim/protocol.ts";
import {
  fillPayload,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "../src/sim/engineFacade.ts";

const jlog = (payload: Record<string, unknown>): void => {
  console.info(JSON.stringify({ suite: "wf-e52c-replay", ...payload }));
};

function payloadAt(tick: number, phaseCode = 1): Float64Array {
  const p = new Float64Array(PAYLOAD_F64S);
  p[0] = tick * 0.1; // x
  p[1] = 2 + tick * 0.01; // h
  p[11] = phaseCode;
  return p;
}

function sealedFixture(): FlightRecording {
  const rec = new FlightRecorder();
  for (const t of [10, 20, 30]) {
    rec.append(t, payloadAt(t));
  }
  return rec.seal({
    scenario: dec17Scenario(1n, MODE_FIXED),
    runIntentId: "intent",
    tick0Digest: "0".repeat(64),
    terminalPhase: "ended:max-ticks",
    finalDigest: "a".repeat(64),
  });
}

test("recorder: monotone ticks enforced; empty seal refuses", () => {
  const rec = new FlightRecorder();
  assert.throws(() => rec.seal({} as never), RangeError, "empty");
  rec.append(5, payloadAt(5));
  assert.throws(() => rec.append(5, payloadAt(5)), RangeError, "equal tick");
  assert.throws(() => rec.append(4, payloadAt(4)), RangeError, "backwards");
  assert.throws(() => rec.append(6, new Float64Array(PAYLOAD_F64S - 1)), RangeError, "short");
  rec.append(6, payloadAt(6));
  assert.equal(rec.frameCount(), 2);
});

test("scenario round-trips through the recorded form (bigints exact)", () => {
  const s = dec17Scenario(18446744073709551615n, MODE_FIXED); // u64::MAX seed
  const round = recordedToScenario(scenarioToRecorded(s));
  assert.deepEqual(round, s);
});

test("parse: fail-closed hostile twins", () => {
  const good = sealedFixture();
  const json = JSON.stringify(good);
  const parsed = parseRecording(json);
  assert.ok(!("error" in parsed), JSON.stringify(parsed));
  const twist = (mutate: (r: Record<string, unknown>) => void): string => {
    const r = JSON.parse(json) as Record<string, unknown>;
    mutate(r);
    return JSON.stringify(r);
  };
  assert.ok("error" in parseRecording("not json {"));
  assert.ok("error" in parseRecording(twist((r) => (r.schema = "v0"))));
  assert.ok("error" in parseRecording(twist((r) => (r.finalDigest = "a".repeat(63)))));
  assert.ok("error" in parseRecording(twist((r) => (r.ticks as number[]).reverse())));
  assert.ok("error" in parseRecording(twist((r) => (r.frames as number[]).push(1))));
  assert.ok("error" in parseRecording(twist((r) => ((r.frames as unknown[])[3] = null))));
  assert.ok("error" in parseRecording(twist((r) => ((r.ticks = []), (r.frames = [])))));
});

test("ghost lookup: before-start null, exact hit, between-hold, hold-last", () => {
  const rec = sealedFixture(); // ticks 10, 20, 30
  assert.equal(ghostAt(rec, 9), null, "before the recording starts");
  assert.equal(ghostAt(rec, 10)?.tick, 10, "exact first");
  assert.equal(ghostAt(rec, 19)?.tick, 10, "hold until the next recorded tick");
  assert.equal(ghostAt(rec, 20)?.tick, 20);
  assert.equal(ghostAt(rec, 25)?.xM, 2, "payload of tick 20 (x = 20*0.1)");
  assert.equal(ghostAt(rec, 30)?.tick, 30);
  assert.equal(ghostAt(rec, 10_000)?.tick, 30, "hold-last past the end, never extrapolate");
});

test("verdict discriminates", () => {
  const rec = sealedFixture();
  assert.equal(replayVerdict(rec, "a".repeat(64)).kind, "identical");
  const bad = replayVerdict(rec, "b".repeat(64));
  assert.equal(bad.kind, "diverged");
  if (bad.kind === "diverged") {
    assert.equal(bad.expectedDigest, "a".repeat(64));
    assert.equal(bad.observedDigest, "b".repeat(64));
  }
});

// ---------------------------------------------------------------------------
// LIVE: record→replay identity against the real engine.
// ---------------------------------------------------------------------------
const pkgPath = process.env.WF_PKG;
if (pkgPath === undefined) {
  jlog({ case: "live-replay", skipped: true, reason: "WF_PKG unset — record→replay identity NOT verified in this run" });
} else {
  const require = createRequire(import.meta.url);
  const wasm = require(pkgPath);

  const runAndRecord = (seed: bigint): FlightRecording => {
    const scenario = { ...dec17Scenario(seed, MODE_FIXED), maxTicks: 30n };
    const init = parseInitEnvelope(
      wasm.flyer_engine_init(
        scenario.seed,
        scenario.rhoKgM3,
        scenario.headwindMps,
        scenario.mode,
        scenario.member,
        scenario.railLengthM,
        scenario.maxTicks,
        scenario.assist,
      ),
    );
    assert.equal(init.kind, "ok");
    if (init.kind !== "ok") throw new Error("unreachable");
    const rec = new FlightRecorder();
    const payload = new Float64Array(PAYLOAD_F64S);
    let terminal = "";
    for (;;) {
      const step = parseStepEnvelope(wasm.flyer_engine_step(false, 0, 0));
      assert.equal(step.kind, "ok");
      if (step.kind !== "ok") throw new Error("unreachable");
      fillPayload(step, payload);
      rec.append(step.tick, payload);
      if (step.ended) {
        terminal = step.phase;
        break;
      }
    }
    const digest = parseDigestEnvelope(wasm.flyer_engine_digest());
    assert.equal(typeof digest, "string");
    return rec.seal({
      scenario,
      runIntentId: init.runIntentId,
      tick0Digest: init.tick0Digest,
      terminalPhase: terminal,
      finalDigest: digest as string,
    });
  };

  test("LIVE: same scenario replays bit-identically; different seed diverges", () => {
    const original = runAndRecord(1903n);
    // Serialize→parse (the storage round trip the ghost path uses).
    const parsed = parseRecording(JSON.stringify(original));
    assert.ok(!("error" in parsed));
    const stored = parsed as FlightRecording;
    // Replay: run the recorded scenario again on the same engine.
    const replayRun = runAndRecord(BigInt(stored.scenario.seed));
    const verdict = replayVerdict(stored, replayRun.finalDigest);
    assert.equal(verdict.kind, "identical", JSON.stringify(verdict));
    assert.deepEqual(replayRun.ticks, stored.ticks);
    assert.deepEqual(replayRun.frames, stored.frames, "full transcript bitwise");
    // The falsifier: a different seed MUST diverge (the check can fail).
    const other = runAndRecord(1904n);
    assert.equal(replayVerdict(stored, other.finalDigest).kind, "diverged");
    // Ghost from the stored recording tracks the replay's states.
    const g = ghostAt(stored, stored.ticks[stored.ticks.length - 1]! + 5);
    assert.equal(g?.tick, stored.ticks[stored.ticks.length - 1]);
    jlog({
      case: "live-replay",
      digest: stored.finalDigest.slice(0, 16),
      frames: stored.ticks.length,
      terminal: stored.terminalPhase,
    });
  });
}
