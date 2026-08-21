// Sim worker entry (bead wf-root-guzez.6.3.1, E5.2a): drives the REAL
// fs-flyer-wasm engine at the 120 Hz bounded-catch-up schedule (E0.8)
// and publishes the frozen 12-float snapshot into the E0.7 seqlock
// ring (SharedArrayBuffer) — postMessage fallback when SAB is
// unavailable. THIN by design: every parse/assemble branch lives in
// engineFacade.ts (headless-tested); this file is transport glue.
//
// wasm pkg: built by `npm run wasm` (wasm-pack --target web) into
// src/wasm-pkg/ (gitignored, derived artifact — the Rust crate is the
// source of truth).

import {
  fillPayload,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "./engineFacade.ts";
import {
  MODE_HUMAN,
  PAYLOAD_F64S,
  payloadLayoutHash,
  type MainToWorker,
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
  ): string;
  flyer_engine_step(hasInput: boolean, leverN: number, warpRad: number): string;
  flyer_engine_digest(): string;
}

let engine: WasmEngine | null = null;
let writer: SeqlockWriter | null = null;
let scheduler: TickScheduler | null = null;
let mode = 0;
let running = false;
let ended = false;
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

function runTick(tick: number): void {
  if (engine === null || ended) {
    return;
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
      return;
    }
    hasInput = true;
    lever = held.leverForceN;
    warp = held.warpCmdRad;
  }
  const json = engine.flyer_engine_step(hasInput, lever, warp);
  const step = parseStepEnvelope(json);
  if (step.kind === "refusal") {
    post({ kind: "refusal", stage: "step", refusal: step.refusal });
    return;
  }
  if (step.kind === "malformed") {
    jlog("malformed-step", { detail: step.detail, tick });
    return;
  }
  if (writer !== null) {
    writer.publish(step.tick, (payload) => fillPayload(step, payload));
  } else {
    const copy = new Float64Array(PAYLOAD_F64S);
    fillPayload(step, copy);
    post({ kind: "snapshot", tick: step.tick, payload: copy }, [copy.buffer]);
  }
  if (step.ended) {
    ended = true;
    running = false;
    const digest = parseDigestEnvelope(engine.flyer_engine_digest());
    post({
      kind: "terminal",
      phase: step.phase,
      tick: step.tick,
      ...(step.envelopeRefusalCode !== undefined
        ? { envelopeRefusalCode: step.envelopeRefusalCode }
        : {}),
      digest: typeof digest === "string" ? digest : "unavailable",
    });
    jlog("terminal", { phase: step.phase, tick: step.tick });
  }
}

function pump(): void {
  if (!running || scheduler === null) {
    return;
  }
  scheduler.pump(performance.now(), runTick, () => performance.now());
  setTimeout(pump, TICK_MS / 2);
}

async function handleInit(msg: Extract<MainToWorker, { kind: "init" }>): Promise<void> {
  engine = engine ?? (await loadEngine());
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
  );
  const init = parseInitEnvelope(initJson);
  if (init.kind === "refusal") {
    post({ kind: "refusal", stage: "init", refusal: init.refusal });
    return;
  }
  if (init.kind === "malformed") {
    post({
      kind: "refusal",
      stage: "init",
      refusal: {
        code: "envelope-malformed",
        message: init.detail,
        ranked_repairs: ["rebuild the wasm pkg (npm run wasm) to match the app protocol"],
      },
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
  running = true;
  post({
    kind: "ready",
    runIntentId: init.runIntentId,
    tick0Digest: init.tick0Digest,
    trimVMps: init.trimVMps,
    layoutHash,
  });
  jlog("ready", { runIntentId: init.runIntentId, sab: writer !== null });
  pump();
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
