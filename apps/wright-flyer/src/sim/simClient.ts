// E5.2b main-thread sim client (bead wf-root-guzez.6.3.2): spawns the
// sim worker, negotiates the transport (SharedArrayBuffer seqlock ring
// when cross-origin isolation grants it, postMessage fallback
// otherwise — capability.ts decides), and hands the render loop an
// interpolated snapshot per frame. Decision logic lives in
// snapshotView.ts (pure, tested); this file is browser glue.

import {
  PAYLOAD_F64S,
  payloadLayoutHash,
  type MainToWorker,
  type RefusalEnvelope,
  type ScenarioInit,
  type WorkerToMain,
} from "./protocol.ts";
import { SeqlockReader, seqlockBytes } from "../transport/seqlock.ts";
import { decodeSnapshot, interpolateSnapshots, type SimSnapshot } from "./snapshotView.ts";
import { FlightRecorder, type FlightRecording } from "./replay.ts";
import { estimateClockOffsetMs, type ClockSyncSample } from "../transport/inputClock.ts";

const RING_SLOTS = 4;
const SIM_TICK_S = 1 / 120;

export interface SimClientEvents {
  onReady(info: { runIntentId: string; tick0Digest: string; trimVMps: number }): void;
  /** `stage` includes "worker" for boot/load failures of the worker
   * itself — without this surface, a dead worker silently degrades the
   * app to the scripted attract loop and nobody can tell. */
  onRefusal(stage: "init" | "step" | "checkpoint" | "worker", refusal: RefusalEnvelope): void;
  onTerminal(info: {
    phase: string;
    tick: number;
    digest: string;
    envelopeRefusalCode?: string;
  }): void;
  /** E5.3a: ApplyNextEligibleTickAndFlag receipt for one control. */
  onControlAck?(ack: { sequence: number; appliedTick: number; lateByTicks: number }): void;
  /** Exact live SimLoop checkpoint bytes, bound to request and run identity. */
  onCheckpoint?(checkpoint: { requestId: number; runIntentId: string; bytes: Uint8Array }): void;
}

export class SimClient {
  private readonly worker: Worker;
  private sab: SharedArrayBuffer | null = null;
  private reader: SeqlockReader | null = null;
  private readonly scratch = new Float64Array(PAYLOAD_F64S);
  private prev: SimSnapshot | null = null;
  private latest: SimSnapshot | null = null;
  private latestArrivalMs = 0;
  private terminalCode: string | undefined;
  // E5.2c: every run records itself; the sealed recording drives the
  // replay ghost and the digest-identity verdict.
  private recorder = new FlightRecorder();
  private scenario: ScenarioInit | null = null;
  private ready: { runIntentId: string; tick0Digest: string } | null = null;
  private recording: FlightRecording | null = null;
  // E5.3a clock sync (main→worker monotonic offset).
  private readonly syncSamples: ClockSyncSample[] = [];
  private clockOffsetMs = 0;
  private pingNonce = 0;
  private initGeneration = 0;
  private checkpointRequestId = 0;
  private readonly pendingCheckpoints = new Map<number, string>();
  constructor(events: SimClientEvents, workerFactory?: () => Worker) {
    this.worker =
      workerFactory !== undefined
        ? workerFactory()
        : new Worker(new URL("./simWorker.ts", import.meta.url), { type: "module" });
    // Surface worker boot/load failures LOUDLY: a dead worker otherwise
    // degrades the app to the scripted attract loop with no diagnostic
    // (observed live as a silent fallback — fresh-eyes fix).
    this.worker.addEventListener("error", (event) => {
      const message =
        event.message !== undefined && event.message !== ""
          ? event.message
          : "sim worker failed to load or threw during boot";
      events.onRefusal("worker", {
        code: "worker-boot-failed",
        message,
        ranked_repairs: [
          "reload the page (a stale optimized-dependency reload can kill the first boot)",
          "check the dev server console for the underlying module error",
        ],
      });
    });
    this.worker.addEventListener("messageerror", () => {
      events.onRefusal("worker", {
        code: "worker-message-deserialize-failed",
        message: "a worker message could not be deserialized",
        ranked_repairs: ["reload the page"],
      });
    });
    this.worker.addEventListener("message", (event: MessageEvent<WorkerToMain>) => {
      const msg = event.data;
      switch (msg.kind) {
        case "ready":
          if (msg.initGeneration !== this.initGeneration) {
            break;
          }
          this.ready = { runIntentId: msg.runIntentId, tick0Digest: msg.tick0Digest };
          events.onReady(msg);
          break;
        case "refusal":
          if (msg.stage !== "init" || msg.initGeneration === this.initGeneration) {
            events.onRefusal(msg.stage, msg.refusal);
          }
          break;
        case "terminal":
          if (
            this.ready?.runIntentId !== msg.runIntentId ||
            msg.initGeneration !== this.initGeneration
          ) {
            break;
          }
          this.terminalCode = msg.envelopeRefusalCode;
          if (this.scenario !== null && this.ready !== null && this.recorder.frameCount() > 0) {
            this.recording = this.recorder.seal({
              scenario: this.scenario,
              runIntentId: this.ready.runIntentId,
              tick0Digest: this.ready.tick0Digest,
              terminalPhase: msg.phase,
              finalDigest: msg.digest,
            });
          }
          events.onTerminal(msg);
          break;
        case "snapshot": {
          // postMessage fallback transport.
          if (
            this.ready?.runIntentId !== msg.runIntentId ||
            msg.initGeneration !== this.initGeneration
          ) {
            break;
          }
          this.recorder.append(msg.tick, msg.payload);
          this.push(decodeSnapshot(msg.tick, msg.payload));
          break;
        }
        case "checkpoint":
          if (
            this.pendingCheckpoints.get(msg.requestId) === msg.runIntentId &&
            this.ready?.runIntentId === msg.runIntentId
          ) {
            this.pendingCheckpoints.delete(msg.requestId);
            events.onCheckpoint?.(msg);
          }
          break;
        case "checkpoint-refusal":
          if (
            this.pendingCheckpoints.get(msg.requestId) === msg.runIntentId &&
            this.ready?.runIntentId === msg.runIntentId
          ) {
            this.pendingCheckpoints.delete(msg.requestId);
            events.onRefusal("checkpoint", msg.refusal);
          }
          break;
        case "metrics":
          console.info(JSON.stringify({ suite: "wf-sim-client", stage: "metrics", ...msg }));
          break;
        case "pong": {
          // Clock sync: worker ≈ main + offset (min-RTT midpoint).
          this.syncSamples.push({
            localSentMs: msg.localSentMs,
            remoteMs: msg.remoteMs,
            localReceivedMs: performance.now(),
          });
          this.clockOffsetMs = estimateClockOffsetMs(this.syncSamples);
          break;
        }
        case "control-ack":
          events.onControlAck?.(msg);
          break;
      }
    });
  }

  /** Start a run. Uses the SAB ring iff the environment grants SAB. */
  start(scenario: ScenarioInit, runEpoch = 1): void {
    this.initGeneration += 1;
    const initGeneration = this.initGeneration;
    this.scenario = scenario;
    this.ready = null;
    this.pendingCheckpoints.clear();
    this.recorder = new FlightRecorder();
    this.recording = null;
    let init: Extract<MainToWorker, { kind: "init" }>;
    if (typeof SharedArrayBuffer !== "undefined" && crossOriginIsolated) {
      const sab = new SharedArrayBuffer(
        seqlockBytes({ slots: RING_SLOTS, payloadF64s: PAYLOAD_F64S }),
      );
      this.sab = sab;
      this.reader = new SeqlockReader(
        sab,
        { slots: RING_SLOTS, payloadF64s: PAYLOAD_F64S },
        { runEpoch, layoutHash: payloadLayoutHash(), anchorPrefix: 0 },
      );
      init = { kind: "init", scenario, sab, slots: RING_SLOTS, runEpoch, initGeneration };
    } else {
      this.reader = null;
      init = { kind: "init", scenario, runEpoch, initGeneration };
    }
    this.worker.postMessage(init);
  }

  /** Rebind the reader's anchor once tick0 is known (SAB path): the
   * header identity check includes the anchor prefix the worker stamps
   * from tick0, so the reader is recreated with the real value. */
  bindAnchor(tick0Digest: string, runEpoch = 1): void {
    if (this.sab !== null) {
      const anchorPrefix = Number.parseInt(tick0Digest.slice(0, 8), 16) | 0;
      this.reader = new SeqlockReader(
        this.sab,
        { slots: RING_SLOTS, payloadF64s: PAYLOAD_F64S },
        { runEpoch, layoutHash: payloadLayoutHash(), anchorPrefix },
      );
    }
  }

  /** Send a control sample. `deviceMs` is the MAIN-clock device event
   * time; the ping-derived offset translates it into the worker clock. */
  sendControl(leverForceN: number, warpCmdRad: number, sequence: number, deviceMs: number): void {
    this.worker.postMessage({
      kind: "control",
      leverForceN,
      warpCmdRad,
      sequence,
      deviceWorkerMs: deviceMs + this.clockOffsetMs,
    } satisfies MainToWorker);
  }

  /** One clock-sync ping (call a few times; min-RTT sample wins). */
  sendPing(): void {
    this.pingNonce += 1;
    this.worker.postMessage({
      kind: "ping",
      nonce: this.pingNonce,
      localSentMs: performance.now(),
    } satisfies MainToWorker);
  }

  /** Request the current live SimLoop checkpoint from the worker. */
  requestCheckpoint(): boolean {
    const runIntentId = this.ready?.runIntentId;
    if (runIntentId === undefined) {
      return false;
    }
    this.checkpointRequestId += 1;
    const requestId = this.checkpointRequestId;
    this.pendingCheckpoints.set(requestId, runIntentId);
    this.worker.postMessage({ kind: "checkpoint", requestId, runIntentId } satisfies MainToWorker);
    return true;
  }

  private push(snap: SimSnapshot): void {
    if (this.latest === null || snap.tick > this.latest.tick) {
      this.prev = this.latest;
      this.latest = snap;
      this.latestArrivalMs = performance.now();
    }
  }

  /**
   * The frame sample: pulls the newest ring snapshot (SAB path), then
   * interpolates between the two most recent sim states at the render
   * clock. Returns null before the first snapshot.
   */
  sample(nowMs: number): SimSnapshot | null {
    if (this.reader !== null) {
      const result = this.reader.read(this.scratch);
      if (typeof result === "number") {
        if (this.latest === null || result > this.latest.tick) {
          this.recorder.append(result, this.scratch);
          this.push(decodeSnapshot(result, this.scratch));
        }
      }
    }
    if (this.latest === null) {
      return null;
    }
    if (this.prev === null || this.latest.ended) {
      return this.latest;
    }
    const dtTicks = Math.max(1, this.latest.tick - this.prev.tick);
    const alpha = (nowMs - this.latestArrivalMs) / 1000 / (dtTicks * SIM_TICK_S);
    return interpolateSnapshots(this.prev, this.latest, Math.min(1, Math.max(0, alpha)));
  }

  envelopeRefusalCode(): string | undefined {
    return this.terminalCode;
  }

  /** Newest sim tick seen by this client (0 before the first snapshot). */
  latestTick(): number {
    return this.latest?.tick ?? 0;
  }

  /** The sealed recording of the finished run (null until terminal).
   * NOTE (declared): the SAB path records at the RENDER sample rate —
   * frames the reader never leased are not in the transcript; the
   * DIGEST identity is still exact because it is the engine's chained
   * digest, not a transcript hash. */
  takeRecording(): FlightRecording | null {
    return this.recording;
  }

  dispose(): void {
    this.worker.terminate();
  }
}
