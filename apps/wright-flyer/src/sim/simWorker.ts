// Sim worker entry (bead wf-root-guzez.6.3.1, E5.2a): drives the REAL
// fs-flyer-wasm engine at the 120 Hz bounded-catch-up schedule (E0.8)
// and publishes the frozen v2 snapshot into the E0.7 seqlock
// ring (SharedArrayBuffer) — postMessage fallback when SAB is
// unavailable. THIN by design: every parse/assemble branch lives in
// engineFacade.ts (headless-tested); this file is transport glue.
//
// wasm pkg: built by `npm run wasm` (wasm-pack --target web) into
// src/wasm-pkg/ (gitignored, derived artifact — the Rust crate is the
// source of truth).

import { TickScheduler } from "../transport/schedule.ts";
import { SeqlockWriter } from "../transport/seqlock.ts";
import {
  type EngineInitializationApi,
  fillPayload,
  initializeScenario,
  parseCheckpointEnvelope,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "./engineFacade.ts";
import { AppliedControls, appliedControlTraceError } from "./appliedControls.ts";
import {
  type MainToWorker,
  MODE_HUMAN,
  PAYLOAD_F64S,
  payloadLayoutHash,
  type RefusalEnvelope,
  type WorkerToMain,
} from "./protocol.ts";

const TICK_MS = 1000 / 120;

interface WasmEngine extends EngineInitializationApi {
  flyer_engine_step(hasInput: boolean, leverN: number, warpRad: number): string;
  flyer_engine_digest(): string;
  flyer_engine_checkpoint(): string;
}

let engine: WasmEngine | null = null;
let writer: SeqlockWriter | null = null;
let scheduler: TickScheduler | null = null;
let mode = 0;
let running = false;
let ended = false;
let activeRunIntentId: string | null = null;
let activeInitGeneration: number | null = null;
let initGeneration = 0;
let pumpGeneration = 0;
// E5.3a: Human-mode hold law (ApplyNextEligibleTickAndFlag) + the sim
// epoch for requested-tick targeting (worker monotonic clock).
let controls = new AppliedControls();
let simEpochMs = 0;

function post(msg: WorkerToMain, transfer?: Transferable[]): void {
  (self as unknown as Worker).postMessage(msg, transfer ?? []);
}

function jlog(stage: string, payload: Record<string, unknown>): void {
  console.info(JSON.stringify({ suite: "wf-sim-worker", stage, ...payload }));
}

async function loadEngine(): Promise<WasmEngine> {
  const pkg = await import("../wasm-pkg/fs_flyer_wasm.js");
  await pkg.default();
  return pkg as unknown as WasmEngine;
}

function refuseStep(runIntentId: string, initGeneration: number, refusal: RefusalEnvelope): false {
  // Native errors may follow partial advancement. Never resume or checkpoint a
  // run whose accepted input/snapshot sequence we can no longer establish.
  invalidateRun();
  post({ kind: "refusal", stage: "step", runIntentId, initGeneration, refusal });
  return false;
}

function runTick(tick: number): boolean {
  const runIntentId = activeRunIntentId;
  const currentInitGeneration = activeInitGeneration;
  if (engine === null || ended || runIntentId === null || currentInitGeneration === null) {
    return false; // wasm still loading / run over: sim time waits
  }
  const nativeEngine = engine;
  let step: ReturnType<AppliedControls["step"]>;
  try {
    step = mode === MODE_HUMAN
      ? controls.step(tick, (value) => nativeEngine.flyer_engine_step(true, value.leverForceN, value.warpCmdRad))
      : parseStepEnvelope(nativeEngine.flyer_engine_step(false, 0, 0));
  } catch (error) {
    return refuseStep(runIntentId, currentInitGeneration, {
      code: "native-step-failed", message: String(error),
      ranked_repairs: ["restart the scenario with a matching wasm package"],
    });
  }
  if (step.kind === "waiting") return false;
  if (step.kind === "refusal") {
    return refuseStep(runIntentId, currentInitGeneration, step.refusal);
  }
  if (step.kind === "malformed" || step.tick !== tick) {
    return refuseStep(runIntentId, currentInitGeneration, {
      code: "step-envelope-malformed",
      message: step.kind === "malformed" ? step.detail : "native tick differs from scheduled tick",
      ranked_repairs: ["restart the scenario with a matching wasm package"],
    });
  }
  if (writer !== null) {
    writer.publish(step.tick, (payload) => fillPayload(step, payload));
  } else if (!step.ended) {
    const copy = new Float64Array(PAYLOAD_F64S);
    fillPayload(step, copy);
    post(
      {
        kind: "snapshot",
        runIntentId,
        initGeneration: currentInitGeneration,
        tick: step.tick,
        payload: copy,
      },
      [copy.buffer],
    );
  }
  if (step.ended) {
    ended = true;
    running = false;
    const digest = parseDigestEnvelope(engine.flyer_engine_digest());
    if (typeof digest !== "string") {
      return refuseStep(runIntentId, currentInitGeneration, {
        code: "terminal-digest-unavailable", message: "native terminal digest could not be read",
        ranked_repairs: ["restart the scenario with a matching wasm package"],
      });
    }
    const payload = new Float64Array(PAYLOAD_F64S);
    fillPayload(step, payload);
    post({
      kind: "terminal",
      runIntentId,
      initGeneration: currentInitGeneration,
      phase: step.phase,
      tick: step.tick,
      ...(step.envelopeRefusalCode !== undefined
        ? { envelopeRefusalCode: step.envelopeRefusalCode }
        : {}),
      digest,
      payload,
      ...(mode === MODE_HUMAN ? { appliedControls: controls.seal() } : {}),
    }, [payload.buffer]);
    jlog("terminal", { phase: step.phase, tick: step.tick });
  }
  return true;
}

function pump(generation = pumpGeneration): void {
  if (generation !== pumpGeneration || !running || scheduler === null) {
    return;
  }
  scheduler.pump(performance.now(), runTick, () => performance.now());
  setTimeout(() => pump(generation), TICK_MS / 2);
}

function invalidateRun(): void {
  pumpGeneration += 1;
  running = false;
  ended = false;
  writer = null;
  scheduler = null;
  activeRunIntentId = null;
  activeInitGeneration = null;
}

function postInitRefusal(initGeneration: number, refusal: RefusalEnvelope): void {
  post({ kind: "refusal", stage: "init", initGeneration, refusal });
}

async function handleInit(msg: Extract<MainToWorker, { kind: "init" }>): Promise<void> {
  const generation = ++initGeneration;
  // Stop publication before loading or admitting B. A failed B must not leave
  // A's scheduler, ring writer, or checkpoint authority live.
  invalidateRun();
  if (msg.replay !== undefined) {
    const { runIntentId, tick0Digest, appliedControls } = msg.replay;
    const traceError = appliedControls === undefined ? null : appliedControlTraceError(appliedControls);
    if (![runIntentId, tick0Digest].every((id) => typeof id === "string" && /^[0-9a-f]{64}$/.test(id)) ||
        (msg.scenario.mode === MODE_HUMAN) !== (appliedControls !== undefined) || traceError !== null ||
        (appliedControls !== undefined && BigInt(appliedControls.endTick) > msg.scenario.maxTicks)) {
      postInitRefusal(msg.initGeneration, {
        code: "replay-input-invalid", message: traceError ?? "replay identity, mode or input extent does not match the scenario",
        ranked_repairs: ["use a complete recording from this engine and scenario"],
      });
      return;
    }
  }
  engine = engine ?? (await loadEngine());
  if (generation !== initGeneration) {
    return;
  }
  mode = msg.scenario.mode;
  ended = false;
  controls = new AppliedControls(msg.replay?.appliedControls);
  const initJson = initializeScenario(engine, msg.scenario);
  const init = parseInitEnvelope(initJson);
  if (init.kind === "refusal") {
    postInitRefusal(msg.initGeneration, init.refusal);
    return;
  }
  if (generation !== initGeneration) {
    return;
  }
  if (init.kind === "malformed") {
    postInitRefusal(msg.initGeneration, {
      code: "envelope-malformed",
      message: init.detail,
      ranked_repairs: ["rebuild the wasm pkg (npm run wasm) to match the app protocol"],
    });
    return;
  }
  if (msg.replay !== undefined && (msg.replay.runIntentId !== init.runIntentId || msg.replay.tick0Digest !== init.tick0Digest)) {
    postInitRefusal(msg.initGeneration, {
      code: "replay-identity-mismatch",
      message: "the initialized native scenario or tick-zero state differs from the recording",
      ranked_repairs: ["use the original scenario and matching native engine build"],
    });
    return;
  }
  const layoutHash = payloadLayoutHash();
  if (msg.sab !== undefined && msg.slots !== undefined) {
    const anchorPrefix = Number.parseInt(init.tick0Digest.slice(0, 8), 16) | 0;
    writer = new SeqlockWriter(
      msg.sab,
      { slots: msg.slots, payloadF64s: PAYLOAD_F64S },
      msg.runEpoch,
      layoutHash,
      anchorPrefix,
    );
  } else {
    writer = null;
  }
  simEpochMs = performance.now();
  scheduler = new TickScheduler(TICK_MS, simEpochMs);
  activeRunIntentId = init.runIntentId;
  activeInitGeneration = msg.initGeneration;
  running = true;
  post({
    kind: "ready",
    runIntentId: init.runIntentId,
    tick0Digest: init.tick0Digest,
    trimVMps: init.trimVMps,
    layoutHash,
    initGeneration: msg.initGeneration,
  });
  jlog("ready", { runIntentId: init.runIntentId, sab: writer !== null });
  pump();
}

function checkpoint(msg: Extract<MainToWorker, { kind: "checkpoint" }>): void {
  if (engine === null || activeRunIntentId === null) {
    post({
      kind: "checkpoint-refusal",
      requestId: msg.requestId,
      runIntentId: msg.runIntentId,
      refusal: {
        code: "engine-not-initialized",
        message: "call init before checkpoint",
        ranked_repairs: ["start a scenario before requesting its checkpoint"],
      },
    });
    return;
  }
  if (msg.runIntentId !== activeRunIntentId) {
    post({
      kind: "checkpoint-refusal",
      requestId: msg.requestId,
      runIntentId: msg.runIntentId,
      refusal: {
        code: "checkpoint-run-mismatch",
        message: "checkpoint request does not name the active run",
        ranked_repairs: ["wait for the current run ready receipt before requesting a checkpoint"],
      },
    });
    return;
  }
  const result = parseCheckpointEnvelope(engine.flyer_engine_checkpoint());
  if (result.kind === "refusal") {
    post({
      kind: "checkpoint-refusal",
      requestId: msg.requestId,
      runIntentId: msg.runIntentId,
      refusal: result.refusal,
    });
    return;
  }
  if (result.kind === "malformed") {
    post({
      kind: "checkpoint-refusal",
      requestId: msg.requestId,
      runIntentId: msg.runIntentId,
      refusal: {
        code: "checkpoint-envelope-malformed",
        message: result.detail,
        ranked_repairs: ["rebuild the wasm pkg (npm run wasm) to match the app protocol"],
      },
    });
    return;
  }
  if (result.runIntentId !== activeRunIntentId) {
    post({
      kind: "checkpoint-refusal",
      requestId: msg.requestId,
      runIntentId: msg.runIntentId,
      refusal: {
        code: "checkpoint-run-mismatch",
        message: "engine checkpoint identity does not match the active run",
        ranked_repairs: ["restart the scenario; do not persist an identity-mismatched checkpoint"],
      },
    });
    return;
  }
  post(
    {
      kind: "checkpoint",
      requestId: msg.requestId,
      runIntentId: result.runIntentId,
      bytes: result.bytes,
    },
    [result.bytes.buffer],
  );
}

self.addEventListener("message", (event: MessageEvent<MainToWorker>) => {
  const msg = event.data;
  switch (msg.kind) {
    case "init":
      void handleInit(msg);
      break;
    case "control": {
      // Requested tick from the device time (already translated into
      // this worker's clock) under the minimum-lead rule; the ack is
      // the ApplyNextEligibleTickAndFlag receipt.
      const runIntentId = activeRunIntentId;
      const currentInitGeneration = activeInitGeneration;
      if (runIntentId === null || currentInitGeneration === null || !running || ended || mode !== MODE_HUMAN) {
        break;
      }
      if (msg.runIntentId !== runIntentId || msg.initGeneration !== currentInitGeneration) {
        break;
      }
      const currentTick = scheduler?.currentTick() ?? 0;
      const requested = Math.ceil((msg.deviceWorkerMs - simEpochMs) / TICK_MS) + 1;
      const receipt = controls.admit(
        msg.sequence,
        { leverForceN: msg.leverForceN, warpCmdRad: msg.warpCmdRad },
        requested,
        currentTick,
      );
      if (receipt === null) break;
      post({
        kind: "control-ack",
        runIntentId,
        initGeneration: currentInitGeneration,
        sequence: msg.sequence,
        ...receipt,
      });
      break;
    }
    case "ping":
      if (msg.runIntentId !== activeRunIntentId || msg.initGeneration !== activeInitGeneration) {
        break;
      }
      post({
        kind: "pong",
        runIntentId: msg.runIntentId,
        initGeneration: msg.initGeneration,
        nonce: msg.nonce,
        localSentMs: msg.localSentMs,
        remoteMs: performance.now(),
      });
      break;
    case "checkpoint":
      checkpoint(msg);
      break;
    case "pause":
      running = false;
      break;
    case "resume":
      if (!ended && scheduler !== null) {
        running = true;
        pump();
      }
      break;
  }
});
