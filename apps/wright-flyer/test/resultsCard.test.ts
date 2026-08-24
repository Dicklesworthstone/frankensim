// E5.5 results-card battery (bead wf-root-guzez.6.8): per-KPI oracles
// on a hand-built transcript (exact numbers, never totals-only), the
// path-length > downrange discriminator, the KPI-vs-recompute HOSTILE
// TWIN (a tampered card is caught and the message names the field),
// distribution-context wording (no "beat this" targets), and a LIVE
// end-to-end card from a real engine recording.
// Repro: WF_PKG=<pkg>/fs_flyer_wasm.js node --test test/resultsCard.test.ts

import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { test } from "node:test";

import {
  cardLines,
  computeKpis,
  kpiRecomputeDivergence,
  type FlightKpis,
} from "../src/sim/resultsCard.ts";
import {
  FlightRecorder,
  RECORDING_SCHEMA_V1,
  type FlightRecording,
} from "../src/sim/replay.ts";
import {
  MODE_FIXED,
  PAYLOAD_F64S,
  PAYLOAD_F64S_V1,
  dec17Scenario,
} from "../src/sim/protocol.ts";
import {
  fillPayload,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "../src/sim/engineFacade.ts";

function frame(fields: Partial<Record<number, number>>): Float64Array {
  const p = new Float64Array(PAYLOAD_F64S);
  for (const [k, v] of Object.entries(fields)) {
    p[Number(k)] = v as number;
  }
  return p;
}

function fixture(): FlightRecording {
  const rec = new FlightRecorder();
  // Rail: two frames straight along x.
  rec.append(1, frame({ 0: 0, 1: 2.4, 2: 11, 11: 0 }));
  rec.append(2, frame({ 0: 3, 1: 2.4, 2: 12, 11: 0, 6: 0.1 }));
  // Airborne: climb 4 up while going 3 along (path 5 per segment).
  rec.append(3, frame({ 0: 6, 1: 6.4, 2: 13, 3: 1, 4: 0.5, 5: 0.2, 11: 1, 6: 0.2 }));
  rec.append(4, frame({ 0: 9, 1: 10.4, 2: 3, 3: 4, 4: -0.5, 5: -0.1, 11: 1, 6: 0.1 }));
  // Terminal: ground contact straight down 10.4 (path 10.4).
  rec.append(5, frame({ 0: 9, 1: 0, 2: 5, 3: 0, 4: 0.5, 11: 2 }));
  return rec.seal({
    scenario: dec17Scenario(1n, MODE_FIXED),
    runIntentId: "i",
    tick0Digest: "0".repeat(64),
    terminalPhase: "ended:ground-contact",
    finalDigest: "a".repeat(64),
  });
}

function asLegacyV1(current: FlightRecording): FlightRecording {
  const frames: number[] = [];
  for (let frameIndex = 0; frameIndex < current.ticks.length; frameIndex += 1) {
    const start = frameIndex * PAYLOAD_F64S;
    frames.push(...current.frames.slice(start, start + PAYLOAD_F64S_V1));
  }
  return { ...current, schema: RECORDING_SCHEMA_V1, frames };
}

test("per-KPI oracles on the hand-built transcript", () => {
  const k = computeKpis(fixture());
  assert.equal(k.frames, 5);
  assert.equal(k.terminal, "ended:ground-contact");
  assert.equal(k.downrangeM, 9);
  // Path: 3 + 5 + 5 + 10.4 = 23.4 — and it EXCEEDS downrange (the
  // §9.2 discriminator lives).
  assert.ok(Math.abs(k.pathLengthM - 23.4) < 1e-12);
  assert.ok(k.pathLengthM > k.downrangeM);
  assert.deepEqual(k.liftoff, { tick: 3, xM: 6 });
  assert.ok(Math.abs(k.airborneS - 2 / 120) < 1e-15);
  // q: +0.5 then −0.5 airborne = one flip = 0 full undulations; the
  // terminal frame's q is NOT airborne and must not count.
  assert.equal(k.undulations, 0);
  assert.equal(k.maxAbsQRadS, 0.5);
  assert.ok(Math.abs(k.rmsQRadS - 0.5) < 1e-12);
  assert.equal(k.maxAbsThetaRad, 0.2);
  assert.equal(k.maxHM, 10.4);
  // Canard travel: |0.1-0| + |0.2-0.1| + |0.1-0.2| + |0-0.1| = 0.4.
  assert.ok(Math.abs(k.canardTravelRad - 0.4) < 1e-12);
  // Min airspeed while airborne: hypot(3,4)=5 (frame 4).
  assert.equal(k.minAirspeedMps, 5);
});

test("legacy v1 recordings use their 12-word stride for every KPI", () => {
  const current = fixture();
  assert.deepEqual(computeKpis(asLegacyV1(current)), computeKpis(current));
});

test("KPI-vs-recompute hostile twin fires and names the field", () => {
  const rec = fixture();
  const honest = computeKpis(rec);
  assert.equal(kpiRecomputeDivergence(rec, honest), null, "honest card passes");
  const tampered: FlightKpis = { ...honest, downrangeM: honest.downrangeM + 5 };
  const verdict = kpiRecomputeDivergence(rec, tampered);
  assert.ok(verdict !== null, "tampering MUST be caught");
  assert.match(verdict, /^downrangeM:/, "the divergent field is named");
});

test("card wording: distribution context, never a target", () => {
  const lines = cardLines(computeKpis(fixture()), "Kill Devil Hills");
  const text = lines.join("\n");
  assert.match(text, /distributions, not targets/);
  assert.match(text, /852 ft .*only precisely measured/);
  assert.doesNotMatch(text, /beat/i);
  assert.match(text, /downrange 9\.0 m \| path length 23\.4 m/);
});

// LIVE: a real short recording through the card.
const pkgPath = process.env.WF_PKG;
if (pkgPath === undefined) {
  console.info(
    JSON.stringify({ suite: "wf-e55-card", case: "live", skipped: true, reason: "WF_PKG unset" }),
  );
} else {
  const require = createRequire(import.meta.url);
  const wasm = require(pkgPath);
  test("LIVE: card from a real engine run passes the recompute gate", () => {
    const scenario = { ...dec17Scenario(1903n, MODE_FIXED), maxTicks: 30n };
    const init = parseInitEnvelope(
      wasm.flyer_engine_init(
        scenario.seed,
        scenario.rhoKgM3,
        scenario.headwindMps,
        scenario.mode,
        scenario.member,
        scenario.railLengthM,
        scenario.maxTicks,
        false,
        false,
      ),
    );
    assert.equal(init.kind, "ok");
    if (init.kind !== "ok") return;
    const rec = new FlightRecorder();
    const payload = new Float64Array(PAYLOAD_F64S);
    let terminal = "";
    for (;;) {
      const step = parseStepEnvelope(wasm.flyer_engine_step(false, 0, 0));
      assert.equal(step.kind, "ok");
      if (step.kind !== "ok") return;
      fillPayload(step, payload);
      rec.append(step.tick, payload);
      if (step.ended) {
        terminal = step.phase;
        break;
      }
    }
    const digest = parseDigestEnvelope(wasm.flyer_engine_digest());
    const sealed = rec.seal({
      scenario,
      runIntentId: init.runIntentId,
      tick0Digest: init.tick0Digest,
      terminalPhase: terminal,
      finalDigest: digest as string,
    });
    const kpis = computeKpis(sealed);
    assert.equal(kpiRecomputeDivergence(sealed, kpis), null);
    assert.equal(kpis.terminal, "ended:max-ticks");
    assert.equal(kpis.liftoff, null, "30 rail ticks: never airborne");
    assert.ok(kpis.downrangeM > 0);
    const lines = cardLines(kpis, "Kill Devil Hills");
    assert.match(lines.join("\n"), /never left the rail/);
    console.info(
      JSON.stringify({ suite: "wf-e55-card", case: "live", downrange: kpis.downrangeM }),
    );
  });
}
