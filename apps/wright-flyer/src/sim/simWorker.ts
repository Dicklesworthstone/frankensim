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

import {
  fillPayload,
  parseCheckpointEnvelope,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "./engineFacade.ts";
import {
  MODE_HUMAN,
  PAYLOAD_F64S,
  payloadLayoutHash,
  type MainToWorker,
  type RefusalEnvelope,
  type WorkerToMain,
} from "./protocol.ts";
import { SeqlockWriter } from "../transport/seqlock.ts";
import { TickScheduler } from "../transport/schedule.ts";
import { ControlHold } from "./humanControls.ts";

const TICK_MS = 1000 / 120;

interface WasmEngine {
  flyer_engine_init(
    seed: bigint,
    rho: number,
    headwind: number,
    mode: number,
    member: number,
    railM: number,
    maxTicks: bigint,
    assist: boolean,
    catapult: boolean,
  ): string;
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
let hold = new ControlHold();
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

function runTick(tick: number): boolean {
  const runIntentId = activeRunIntentId;
  const currentInitGeneration = activeInitGeneration;
  if (engine === null || ended || runIntentId === null || currentInitGeneration === null) {
    return false; // wasm still loading / run over: sim time waits
  }
  // E5.3a Human mode: the ControlHold supplies the zero-order-held
  // control every tick; before the FIRST admitted control the sim
  // waits (no step, no refusal spam — the run starts at first touch).
  let hasInput = false;
  let lever = 0;
  let warp = 0;
  if (mode === MODE_HUMAN) {
    const held = hold.valueAt(tick);
    if (held === null) {
      // H-2: report "did not step" so TickScheduler does not consume
      // the schedule slot — scheduler ticks stay equal to engine ticks.
      return false;
    }
    hasInput = true;
    lever = held.leverForceN;
    warp = held.warpCmdRad;
  }
  const json = engine.flyer_engine_step(hasInput, lever, warp);
  const step = parseStepEnvelope(json);
  if (step.kind === "refusal") {
    // Refusal consumes the tick as before (typed receipt posted); only
    // the pre-first-touch / loading waits return false above.
    post({
      kind: "refusal",
      stage: "step",
      runIntentId,
      initGeneration: currentInitGeneration,
      refusal: step.refusal,
    });
    return true;
  }
  if (step.kind === "malformed") {
    jlog("malformed-step", { detail: step.detail, tick });
    return true;
  }
  if (writer !== null) {
    writer.publish(step.tick, (payload) => fillPayload(step, payload));
  } else {
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
    post({
      kind: "terminal",
      runIntentId,
      initGeneration: currentInitGeneration,
      phase: step.phase,
      tick: step.tick,
      ...(step.envelopeRefusalCode !== undefined
        ? { envelopeRefusalCode: step.envelopeRefusalCode }
        : {}),
      digest: typeof digest === "string" ? digest : "unavailable",
    });
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
  engine = engine ?? (await loadEngine());
  if (generation !== initGeneration) {
    return;
  }
  mode = msg.scenario.mode;
  ended = false;
  hold = new ControlHold();
  const initJson = engine.flyer_engine_init(
    msg.scenario.seed,
    msg.scenario.rhoKgM3,
    msg.scenario.headwindMps,
    msg.scenario.mode,
    msg.scenario.member,
    msg.scenario.railLengthM,
    msg.scenario.maxTicks,
    msg.scenario.assist,
    msg.scenario.catapult,
  );
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
    post({ kind: "checkpoint-refusal", requestId: msg.requestId, runIntentId: msg.runIntentId, refusal: result.refusal });
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
    { kind: "checkpoint", requestId: msg.requestId, runIntentId: result.runIntentId, bytes: result.bytes },
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
      const currentTick = scheduler?.currentTick() ?? 0;
      const requested = Math.ceil((msg.deviceWorkerMs - simEpochMs) / TICK_MS) + 1;
      const receipt = hold.admit(
        msg.sequence,
        { leverForceN: msg.leverForceN, warpCmdRad: msg.warpCmdRad },
        requested,
        currentTick,
      );
      post({ kind: "control-ack", sequence: msg.sequence, ...receipt });
      break;
    }
    case "ping":
      post({
        kind: "pong",
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
