// E5.2a facade battery (bead wf-root-guzez.6.3.1): the pure seam
// between the wasm engine's JSON envelopes and the worker protocol.
// Shape tests use contract-exact synthetic envelopes PLUS hostile
// twins (missing field, non-finite, unknown phase, malformed refusal);
// the LIVE section drives the REAL nodejs wasm pkg when WF_PKG points
// at it (fs_flyer_wasm.js) and proves parser-vs-engine agreement and
// the payload round trip on actual output — loudly reported either way.
// Repro: WF_PKG=<pkg>/fs_flyer_wasm.js node --test test/simFacade.test.ts

import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { test } from "node:test";

import {
  type EngineInitializationApi,
  fillPayload,
  initializeScenario,
  parseCheckpointEnvelope,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "../src/sim/engineFacade.ts";
import {
  dec17Scenario,
  type MainToWorker,
  MODE_FIXED,
  MODE_HUMAN,
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
  PAYLOAD_LAYOUT_V1,
  PAYLOAD_LAYOUT_V2,
  PHASE_CODES,
  payloadLayoutHash,
  type WorkerToMain,
  type SpatialWindInit,
} from "../src/sim/protocol.ts";
import { SimClient } from "../src/sim/simClient.ts";
import { RECORDING_SCHEMA, RECORDING_SCHEMA_HUMAN } from "../src/sim/replay.ts";
import { SeqlockWriter } from "../src/transport/seqlock.ts";

const jlog = (payload: Record<string, unknown>): void => {
  console.info(JSON.stringify({ suite: "wf-e52a-facade", ...payload }));
};

test("spatial initialization routes the complete field and refuses a legacy-only bundle", () => {
  const scenario = dec17Scenario(1903n, MODE_FIXED);
  const wind: SpatialWindInit = {
    sampling: "distributed", roughnessLengthM: 0.001, displacementHeightM: 0,
    referenceHeightM: 2, sigmaMps: 0.3, lengthM: 15, advectionMps: 11,
    correlationTimesS: [0.2, 0.5, 1.2],
  };
  let legacyCalls = 0;
  let spatialCalls = 0;
  let observed: unknown[] = [];
  const engine: EngineInitializationApi = {
    flyer_engine_init() { legacyCalls += 1; return "legacy"; },
    flyer_engine_init_spatial_wind(...args) {
      spatialCalls += 1;
      observed = args;
      return "spatial";
    },
  };
  assert.equal(initializeScenario(engine, scenario), "legacy");
  for (const sampling of ["cg", "distributed"] as const) {
    assert.equal(initializeScenario(engine, { ...scenario, spatialWind: { ...wind, sampling } }), "spatial");
    assert.deepEqual(observed.slice(0, 9), [
      scenario.seed, scenario.rhoKgM3, scenario.headwindMps, scenario.mode,
      scenario.member, scenario.railLengthM, scenario.maxTicks, scenario.assist, scenario.catapult,
    ]);
    assert.deepEqual(observed.slice(9, 16), [sampling === "distributed", 0.001, 0, 2, 0.3, 15, 11]);
    assert.deepEqual(observed[16], new Float64Array(wind.correlationTimesS));
  }
  const unavailable = parseInitEnvelope(initializeScenario(
    { flyer_engine_init: engine.flyer_engine_init }, { ...scenario, spatialWind: wind },
  ));
  assert.equal(unavailable.kind, "refusal");
  if (unavailable.kind === "refusal") assert.equal(unavailable.refusal.code, "spatial-wind-unavailable");
  const invalid = parseInitEnvelope(initializeScenario(engine, {
    ...scenario, spatialWind: { ...wind, correlationTimesS: [0] },
  }));
  assert.equal(invalid.kind, "refusal");
  if (invalid.kind === "refusal") assert.equal(invalid.refusal.code, "spatial-wind-invalid");
  assert.equal(legacyCalls, 1, "spatial requests never reach the legacy fallback");
  assert.equal(spatialCalls, 2, "invalid fields never reach WASM");
});

test("material atmosphere reaches its own native initializer without a scalar-density fallback", () => {
  const wind: SpatialWindInit = {
    sampling: "distributed", roughnessLengthM: 0.003, displacementHeightM: 0,
    referenceHeightM: 2.4, sigmaMps: 0.3, lengthM: 20, advectionMps: 11,
    correlationTimesS: [1.8, 2.1],
  };
  const air = { temperatureK: 293.15, pressurePa: 98000, relativeHumidity: 0.6 };
  const scenario = { ...dec17Scenario(1903n, MODE_FIXED), spatialWind: wind, materialAir: air };
  let calls = 0;
  let received: unknown[] = [];
  const engine: EngineInitializationApi = {
    flyer_engine_init() { throw new Error("material atmosphere cannot use legacy density"); },
    flyer_engine_init_spatial_wind() { throw new Error("material atmosphere cannot lose its gas inputs"); },
    flyer_engine_init_material_air(...args) { calls += 1; received = args; return "material"; },
  };
  for (const sampling of ["cg", "distributed"] as const) {
    assert.equal(initializeScenario(engine, { ...scenario, spatialWind: { ...wind, sampling } }), "material");
    assert.deepEqual(received.slice(0, 11), [
      scenario.seed, air.temperatureK, air.pressurePa, air.relativeHumidity,
      scenario.headwindMps, scenario.mode, scenario.member, scenario.railLengthM,
      scenario.maxTicks, scenario.assist, scenario.catapult,
    ]);
    assert.deepEqual(received.slice(11, 18), [sampling === "distributed", 0.003, 0, 2.4, 0.3, 20, 11]);
    assert.deepEqual(received[18], Float64Array.from(wind.correlationTimesS));
  }
  const unavailable = parseInitEnvelope(initializeScenario({
    flyer_engine_init: engine.flyer_engine_init,
    flyer_engine_init_spatial_wind: engine.flyer_engine_init_spatial_wind,
  }, scenario));
  assert.equal(unavailable.kind, "refusal");
  if (unavailable.kind === "refusal") assert.equal(unavailable.refusal.code, "material-atmosphere-unavailable");
  for (const invalid of [null, {}, { ...air, temperatureK: NaN },
    { ...air, pressurePa: 0 }, { ...air, relativeHumidity: 1.1 }]) {
    const result = parseInitEnvelope(initializeScenario(engine, { ...scenario, materialAir: invalid as never }));
    assert.equal(result.kind, "refusal");
    if (result.kind === "refusal") assert.equal(result.refusal.code, "material-atmosphere-invalid");
  }
  assert.equal(parseInitEnvelope(initializeScenario(engine, {
    ...dec17Scenario(1903n, MODE_FIXED), materialAir: air,
  })).kind, "refusal");
  assert.equal(calls, 2, "malformed or unavailable material runs cannot initialize physics");
});

const STEP_OK =
  '{"ok":{"tick":42,"phase":"airborne","x_m":10.5,"h_m":3.25,"u_mps":13.5,"w_mps":0.5,' +
  '"q_rad_s":0.1,"theta_rad":0.08,"dc_rad":0.12,"warp_rad":0.01,"omega_prop_rad_s":50.5,' +
  '"p_rad_s":0.03,"phi_rad":0.04,"r_rad_s":-0.02,"psi_rad":-0.05,' +
  '"gust_w_mps":0.02,"assist_active":false}}';

test("phase codes mirror the native payload codes exactly", () => {
  assert.deepEqual(PHASE_CODES, {
    "on-rail": 0,
    airborne: 1,
    "ended:ground-contact": 2,
    "ended:rail-end-without-lift": 3,
    "ended:max-ticks": 4,
    "ended:envelope-exceeded": 5,
    "ended:damage-model-unavailable": 6,
  });
});

test("payload layout hash is pinned (layout identity, not crypto)", () => {
  const h = payloadLayoutHash();
  assert.equal(h, payloadLayoutHash(), "deterministic");
  assert.equal(h, payloadLayoutHash(PAYLOAD_LAYOUT_V2));
  assert.notEqual(h, payloadLayoutHash(PAYLOAD_LAYOUT_V1), "v2 refuses a stale v1 ring");
  jlog({ case: "layout-hash", value: h });
});

test("step envelope parses per-field and fills the frozen payload order", () => {
  const step = parseStepEnvelope(STEP_OK);
  assert.equal(step.kind, "ok");
  if (step.kind !== "ok") return;
  assert.equal(step.tick, 42);
  assert.equal(step.phase, "airborne");
  assert.equal(step.ended, false);
  const out = new Float64Array(PAYLOAD_F64S);
  fillPayload(step, out);
  // Per-slot oracle — NEVER a totals/spread check (sum tests are blind
  // to permutation).
  assert.equal(out[P_X_M], 10.5);
  assert.equal(out[P_H_M], 3.25);
  assert.equal(out[P_U_MPS], 13.5);
  assert.equal(out[P_W_MPS], 0.5);
  assert.equal(out[P_Q_RAD_S], 0.1);
  assert.equal(out[P_THETA_RAD], 0.08);
  assert.equal(out[P_DC_RAD], 0.12);
  assert.equal(out[P_WARP_RAD], 0.01);
  assert.equal(out[P_OMEGA_RAD_S], 50.5);
  assert.equal(out[P_GUST_W_MPS], 0.02);
  assert.equal(out[P_ASSIST], 0);
  assert.equal(out[P_PHASE], 1);
  assert.equal(out[P_PHI_RAD], 0.04);
  assert.equal(out[P_PSI_RAD], -0.05);
});

test("terminal phases parse as ended with the envelope receipt carried", () => {
  const json = STEP_OK.replace('"phase":"airborne"', '"phase":"ended:envelope-exceeded"').replace(
    ',"assist_active":false',
    ',"assist_active":true,"envelope_refusal_code":"PropAirframeCouplingDidNotConverge"',
  );
  const step = parseStepEnvelope(json);
  assert.equal(step.kind, "ok");
  if (step.kind !== "ok") return;
  assert.equal(step.ended, true);
  assert.equal(step.envelopeRefusalCode, "PropAirframeCouplingDidNotConverge");
  const out = new Float64Array(PAYLOAD_F64S);
  fillPayload(step, out);
  assert.equal(out[P_PHASE], 5);
  assert.equal(out[P_ASSIST], 1);
});

test("hostile twins refuse fail-closed (malformed, never a guess)", () => {
  // Missing field.
  assert.equal(parseStepEnvelope(STEP_OK.replace('"h_m":3.25,', "")).kind, "malformed");
  // Non-finite number survives JSON as null → wrong type.
  assert.equal(parseStepEnvelope(STEP_OK.replace("3.25", "null")).kind, "malformed");
  // Unknown phase word.
  assert.equal(
    parseStepEnvelope(STEP_OK.replace('"airborne"', '"ended:new-thing"')).kind,
    "malformed",
  );
  // Not JSON at all.
  assert.equal(parseStepEnvelope("<html>proxy error</html>").kind, "malformed");
  // Refusal with a broken shape is malformed, NOT silently accepted.
  assert.equal(parseStepEnvelope('{"refusal":{"code":7}}').kind, "malformed");
  // Well-formed refusal parses as refusal.
  const refusal = parseStepEnvelope(
    '{"refusal":{"code":"run-ended","message":"m","ranked_repairs":["r"]}}',
  );
  assert.equal(refusal.kind, "refusal");
  if (refusal.kind === "refusal") {
    assert.equal(refusal.refusal.code, "run-ended");
  }
  // Init envelope twins.
  assert.equal(parseInitEnvelope('{"ok":{"run_intent_id":"x"}}').kind, "malformed");
  // Digest: 64-hex or refuse (63 chars = cap-1 twin).
  assert.equal(parseDigestEnvelope(`{"ok":{"digest":"${"a".repeat(64)}"}}`), "a".repeat(64));
  const short = parseDigestEnvelope(`{"ok":{"digest":"${"a".repeat(63)}"}}`);
  assert.equal(typeof short === "string" ? "string" : short.kind, "malformed");
});

test("checkpoint envelope decodes only bounded lowercase hex bytes", () => {
  const checkpoint = parseCheckpointEnvelope(
    `{"ok":{"run_intent_id":"${"a".repeat(64)}","checkpoint_hex":"0001fe"}}`,
  );
  assert.equal(checkpoint.kind, "ok");
  if (checkpoint.kind === "ok") {
    assert.deepEqual([...checkpoint.bytes], [0, 1, 254]);
    assert.equal(checkpoint.runIntentId, "a".repeat(64));
  }
  assert.equal(parseCheckpointEnvelope('{"ok":{"checkpoint_hex":"0"}}').kind, "malformed");
  assert.equal(
    parseCheckpointEnvelope(`{"ok":{"run_intent_id":"${"A".repeat(64)}","checkpoint_hex":"00"}}`)
      .kind,
    "malformed",
  );
  assert.equal(
    parseCheckpointEnvelope(`{"ok":{"run_intent_id":"${"a".repeat(64)}","checkpoint_hex":"00FF"}}`)
      .kind,
    "malformed",
  );
  for (const length of [32, 63, 65]) {
    assert.equal(parseCheckpointEnvelope(
      `{"ok":{"run_intent_id":"${"a".repeat(length)}","checkpoint_hex":"00"}}`,
    ).kind, "malformed", `the native content hash is 64 hex characters, not ${length}`);
  }
  const refusal = parseCheckpointEnvelope(
    '{"refusal":{"code":"checkpoint-after-terminal","message":"m","ranked_repairs":["r"]}}',
  );
  assert.equal(refusal.kind, "refusal");
});

class FakeWorker {
  readonly sent: MainToWorker[] = [];
  private readonly listeners = new Map<string, Array<(event: Event) => void>>();

  addEventListener(type: string, listener: (event: Event) => void): void {
    this.listeners.set(type, [...(this.listeners.get(type) ?? []), listener]);
  }

  postMessage(message: MainToWorker): void {
    this.sent.push(message);
  }

  terminate(): void {}

  emit(message: WorkerToMain): void {
    for (const listener of this.listeners.get("message") ?? []) {
      listener({ data: message } as MessageEvent<WorkerToMain>);
    }
  }
}

test("checkpoint replies bind the request and active run across reinit races", () => {
  (globalThis as { crossOriginIsolated?: boolean }).crossOriginIsolated = false;
  const worker = new FakeWorker();
  const checkpoints: Array<{ requestId: number; runIntentId: string; bytes: Uint8Array }> = [];
  const refusals: string[] = [];
  const readyRuns: string[] = [];
  const terminals: string[] = [];
  const controlAcks: number[] = [];
  const client = new SimClient(
    {
      onReady(info): void {
        readyRuns.push(info.runIntentId);
      },
      onRefusal(_stage, refusal): void {
        refusals.push(refusal.code);
      },
      onTerminal(info): void {
        terminals.push(info.digest);
      },
      onControlAck(ack): void {
        controlAcks.push(ack.sequence);
      },
      onCheckpoint(checkpoint): void {
        checkpoints.push(checkpoint);
      },
    },
    () => worker as unknown as Worker,
  );
  const runA = "a".repeat(64);
  const runB = "b".repeat(64);
  const scenarioA = dec17Scenario(1903n, MODE_FIXED);
  const scenarioB = dec17Scenario(1904n, MODE_FIXED);

  client.start(scenarioA);
  const initA = worker.sent.at(-1);
  assert.equal(initA?.kind, "init");
  if (initA?.kind === "init") {
    assert.equal(initA.initGeneration, 1);
  }
  worker.emit({
    kind: "ready",
    runIntentId: runA,
    tick0Digest: "c".repeat(64),
    trimVMps: 10,
    layoutHash: 1,
    initGeneration: 1,
  });
  worker.emit({
    kind: "snapshot",
    runIntentId: runA,
    initGeneration: 1,
    tick: 7,
    payload: new Float64Array(PAYLOAD_F64S),
  });
  assert.equal(client.latestTick(), 7, "accepted A snapshot reaches the render cache");
  assert.equal(client.sample(0)?.tick, 7);
  assert.equal(client.requestCheckpoint(), true);
  const r1 = worker.sent.at(-1);
  assert.deepEqual(r1, { kind: "checkpoint", requestId: 1, runIntentId: runA });
  client.sendPing();
  const pingA = worker.sent.at(-1);
  assert.equal(pingA?.kind, "ping");
  if (pingA?.kind === "ping") {
    assert.equal(pingA.runIntentId, runA);
    assert.equal(pingA.initGeneration, 1);
  }

  // Starting B revokes A's ready capability; no B checkpoint can be exposed
  // until B's new ready receipt arrives.
  client.start(scenarioB);
  assert.equal(client.latestTick(), 0, "B start must clear A's cached latest tick");
  assert.equal(client.sample(0), null, "B cannot render A while awaiting B's first snapshot");
  const initB = worker.sent.at(-1);
  assert.equal(initB?.kind, "init");
  if (initB?.kind === "init") {
    assert.equal(initB.initGeneration, 2);
  }
  const postsBeforeBReadyControl = worker.sent.length;
  client.sendControl(0, 0, 17, 0);
  assert.equal(
    worker.sent.length,
    postsBeforeBReadyControl,
    "start(B) before B ready must not post a control for either lifecycle",
  );
  // Queued A receipts must not re-establish an A capability after B starts.
  worker.emit({
    kind: "ready",
    runIntentId: runA,
    tick0Digest: "c".repeat(64),
    trimVMps: 10,
    layoutHash: 1,
    initGeneration: 1,
  });
  worker.emit({
    kind: "refusal",
    stage: "init",
    initGeneration: 1,
    refusal: { code: "scenario-invalid", message: "stale", ranked_repairs: [] },
  });
  worker.emit({
    kind: "refusal",
    stage: "step",
    runIntentId: runA,
    initGeneration: 1,
    refusal: { code: "run-ended", message: "stale", ranked_repairs: [] },
  });
  assert.deepEqual(readyRuns, [runA]);
  assert.deepEqual(refusals, []);
  worker.emit({
    kind: "pong",
    runIntentId: runA,
    initGeneration: 1,
    nonce: 1,
    localSentMs: 0,
    remoteMs: 1_000_000,
  });
  worker.emit({
    kind: "control-ack",
    runIntentId: runA,
    initGeneration: 1,
    sequence: 17,
    appliedTick: 8,
    lateByTicks: 0,
  });
  assert.deepEqual(controlAcks, [], "queued A control acknowledgement must not notify B");
  assert.equal(client.requestCheckpoint(), false);
  worker.emit({
    kind: "snapshot",
    runIntentId: runA,
    initGeneration: 1,
    tick: 99,
    payload: new Float64Array(PAYLOAD_F64S),
  });
  worker.emit({
    kind: "terminal",
    runIntentId: runA,
    initGeneration: 1,
    phase: "ended:max-ticks",
    tick: 99,
    digest: "e".repeat(64),
    payload: new Float64Array(PAYLOAD_F64S),
  });
  assert.equal(client.latestTick(), 0, "queued A snapshot must not enter B view state");
  assert.deepEqual(terminals, [], "queued A terminal must not notify or seal B");
  assert.equal(client.takeRecording(), null, "queued A terminal must not create a B recording");
  worker.emit({ kind: "checkpoint", requestId: 1, runIntentId: runA, bytes: new Uint8Array([1]) });
  assert.deepEqual(checkpoints, []);

  worker.emit({
    kind: "ready",
    runIntentId: runB,
    tick0Digest: "d".repeat(64),
    trimVMps: 10,
    layoutHash: 1,
    initGeneration: 2,
  });
  client.sendControl(0, 0, 18, 0);
  const controlB = worker.sent.at(-1);
  assert.deepEqual(controlB, {
    kind: "control",
    runIntentId: runB,
    initGeneration: 2,
    leverForceN: 0,
    warpCmdRad: 0,
    sequence: 18,
    deviceWorkerMs: 0,
  });
  assert.equal(client.requestCheckpoint(), true);
  const r2 = worker.sent.at(-1);
  assert.deepEqual(r2, { kind: "checkpoint", requestId: 2, runIntentId: runB });
  assert.deepEqual(readyRuns, [runA, runB]);
  worker.emit({
    kind: "control-ack",
    runIntentId: runB,
    initGeneration: 2,
    sequence: 18,
    appliedTick: 9,
    lateByTicks: 1,
  });
  assert.deepEqual(controlAcks, [18], "the admitted B acknowledgement remains observable");

  // Out-of-order r1 and a forged r2 run identity are both stale authority.
  worker.emit({
    kind: "checkpoint-refusal",
    requestId: 1,
    runIntentId: runA,
    refusal: { code: "checkpoint-run-mismatch", message: "stale", ranked_repairs: [] },
  });
  worker.emit({ kind: "checkpoint", requestId: 2, runIntentId: runA, bytes: new Uint8Array([2]) });
  assert.deepEqual(checkpoints, []);
  assert.deepEqual(refusals, []);

  worker.emit({ kind: "checkpoint", requestId: 2, runIntentId: runB, bytes: new Uint8Array([3]) });
  assert.deepEqual(
    checkpoints.map(({ requestId, runIntentId, bytes }) => [requestId, runIntentId, [...bytes]]),
    [[2, runB, [3]]],
  );
  assert.equal(client.requestCheckpoint(), true);
  worker.emit({ kind: "checkpoint-refusal", requestId: 3, runIntentId: runB,
    refusal: { code: "checkpoint-unavailable", message: "save refused", ranked_repairs: [] } });
  worker.emit({ kind: "snapshot", runIntentId: runB, initGeneration: 2,
    tick: 10, payload: new Float64Array(PAYLOAD_F64S) });
  assert.equal(client.latestTick(), 10, "a refused checkpoint must not invalidate continuing flight");
  assert.equal(client.requestCheckpoint(), true);
  client.dispose();
});

test("terminal seals human replay without render sampling in fallback and SAB modes", () => {
  const runIntentId = "a".repeat(64);
  const tick0Digest = "b".repeat(64);
  for (const isolated of [false, true]) {
    (globalThis as { crossOriginIsolated?: boolean }).crossOriginIsolated = isolated;
    const worker = new FakeWorker();
    const refusals: string[] = [];
    let terminals = 0;
    const client = new SimClient({
      onReady: () => {}, onRefusal: (_stage, r) => refusals.push(r.code),
      onTerminal: () => { terminals += 1; },
    }, () => worker as unknown as Worker);
    client.start(dec17Scenario(1903n, MODE_HUMAN));
    const init = worker.sent[0];
    assert.equal(init?.kind, "init");
    if (init?.kind !== "init") throw new Error("init missing");
    worker.emit({ kind: "ready", runIntentId, tick0Digest, trimVMps: 10,
      layoutHash: payloadLayoutHash(), initGeneration: 1 });
    client.bindAnchor(tick0Digest);
    const payload = new Float64Array(PAYLOAD_F64S); payload[P_PHASE] = 4; payload[P_X_M] = 12;
    if (isolated) {
      assert.ok(init.sab);
      new SeqlockWriter(init.sab, { slots: 4, payloadF64s: PAYLOAD_F64S }, 1,
        payloadLayoutHash(), Number.parseInt(tick0Digest.slice(0, 8), 16) | 0)
        .publish(3, (out) => out.set(payload));
    }
    const terminal = { kind: "terminal", runIntentId, initGeneration: 1,
      phase: "ended:max-ticks", tick: 3, digest: "c".repeat(64), payload,
      appliedControls: { endTick: 3, changes: [{ tick: 1, leverForceN: "-0", warpCmdRad: 0.01 }] },
    } as const;
    worker.emit(terminal);
    worker.emit(terminal);
    const recording = client.takeRecording();
    assert.ok(recording);
    assert.equal(recording.schema, RECORDING_SCHEMA_HUMAN);
    assert.deepEqual(recording.ticks, [3]);
    assert.equal(client.sample(0)?.xM, 12);
    assert.equal(terminals, 1);
    assert.deepEqual(refusals, []);
    assert.equal(client.startReplay(recording), true);
    const replay = worker.sent.at(-1);
    assert.equal(replay?.kind, "init");
    if (replay?.kind !== "init") throw new Error("replay missing");
    assert.deepEqual(replay.replay, { runIntentId, tick0Digest, appliedControls: terminal.appliedControls });
    worker.emit({ kind: "ready", runIntentId, tick0Digest, trimVMps: 10,
      layoutHash: payloadLayoutHash(), initGeneration: 2 });
    const posted = worker.sent.length;
    client.sendControl(220, 0.148, 3, 0);
    assert.equal(worker.sent.length, posted, "live input cannot contaminate replay");
    worker.emit(terminal);
    assert.equal(client.takeRecording(), null, "old terminal cannot seal the replay lifecycle");
    assert.equal(client.startReplay({ ...recording, schema: RECORDING_SCHEMA, appliedControls: undefined }), false);
    assert.equal(worker.sent.length, posted, "legacy human trace cannot be guessed");
    client.dispose();
  }
  (globalThis as { crossOriginIsolated?: boolean }).crossOriginIsolated = false;
});

// --------------------------------------------------------------------------
// LIVE section: the real wasm engine (nodejs target) when WF_PKG is set.
// --------------------------------------------------------------------------
const pkgPath = process.env.WF_PKG;
if (pkgPath === undefined) {
  // Loud, structured, and impossible to mistake for coverage.
  jlog({
    case: "live-engine",
    skipped: true,
    reason: "WF_PKG unset — parser-vs-engine agreement NOT verified in this run",
  });
} else {
  const require = createRequire(import.meta.url);
  const wasm = require(pkgPath);

  test("LIVE: init envelope from the real engine parses with identity fields", () => {
    const init = parseInitEnvelope(
      wasm.flyer_engine_init(1903n, 1.294, 11.0, MODE_FIXED, 0, 18.3, 24n, false, false),
    );
    assert.equal(init.kind, "ok");
    if (init.kind !== "ok") return;
    assert.match(init.tick0Digest, /^[0-9a-f]{64}$/);
    assert.ok(init.trimVMps > 5.0);
    jlog({ case: "live-init", tick0: init.tick0Digest.slice(0, 16) });
  });

  test("LIVE: every step parses; payload round-trips; terminal in-band", () => {
    let last: ReturnType<typeof parseStepEnvelope> | null = null;
    const out = new Float64Array(PAYLOAD_F64S);
    for (let i = 0; i < 24; i += 1) {
      const step = parseStepEnvelope(wasm.flyer_engine_step(false, 0, 0));
      assert.equal(step.kind, "ok", JSON.stringify(step));
      if (step.kind !== "ok") return;
      fillPayload(step, out);
      assert.equal(out[P_PHASE], PHASE_CODES[step.phase]);
      assert.equal(out[P_U_MPS], step.state.uMps, "bitwise field carry");
      last = step;
    }
    assert.ok(last !== null && last.kind === "ok" && last.ended, "max-ticks terminal parsed");
    const digest = parseDigestEnvelope(wasm.flyer_engine_digest());
    assert.equal(typeof digest, "string", "real digest parses");
    jlog({
      case: "live-lifecycle",
      digest: typeof digest === "string" ? digest.slice(0, 16) : "?",
    });
  });

  test("LIVE: real refusal envelopes parse as refusals with stable codes", () => {
    const past = parseStepEnvelope(wasm.flyer_engine_step(false, 0, 0));
    assert.equal(past.kind, "refusal");
    if (past.kind === "refusal") {
      assert.equal(past.refusal.code, "run-ended");
    }
    const bad = parseInitEnvelope(
      wasm.flyer_engine_init(1n, 1.294, 11.0, 3, 0, 18.3, 24n, false, false),
    );
    assert.equal(bad.kind, "refusal");
    if (bad.kind === "refusal") {
      assert.equal(bad.refusal.code, "mode-invalid");
    }
  });
}
