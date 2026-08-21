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
let control = { leverForceN: 0, warpCmdRad: 0, fresh: false };
let payloadScratch = new Float64Array(PAYLOAD_F64S);

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
  const hasInput = mode === MODE_HUMAN && control.fresh;
  const json = engine.flyer_engine_step(hasInput, control.leverForceN, control.warpCmdRad);
  control = { ...control, fresh: mode === MODE_HUMAN ? control.fresh : false };
  const step = parseStepEnvelope(json);
  if (step.kind === "refusal") {
    post({ kind: "refusal", stage: "step", refusal: step.refusal });
    // control-input-missing in Human mode: hold the schedule, wait for
    // input (ApplyNextEligibleTickAndFlag semantics land in E5.3a).
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
  const initJson = engine.flyer_engine_init(
    msg.scenario.seed,
    msg.scenario.rhoKgM3,
    msg.scenario.headwindMps,
    msg.scenario.mode,
    msg.scenario.member,
    msg.scenario.railLengthM,
    msg.scenario.maxTicks,
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
  scheduler = new TickScheduler(TICK_MS, performance.now());
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
    case "control":
      control = { leverForceN: msg.leverForceN, warpCmdRad: msg.warpCmdRad, fresh: true };
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
