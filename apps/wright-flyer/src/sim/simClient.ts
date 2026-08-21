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

const RING_SLOTS = 4;
const SIM_TICK_S = 1 / 120;

export interface SimClientEvents {
  onReady(info: { runIntentId: string; tick0Digest: string; trimVMps: number }): void;
  onRefusal(stage: "init" | "step", refusal: RefusalEnvelope): void;
  onTerminal(info: {
    phase: string;
    tick: number;
    digest: string;
    envelopeRefusalCode?: string;
  }): void;
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

  constructor(events: SimClientEvents, workerFactory?: () => Worker) {
    this.worker =
      workerFactory !== undefined
        ? workerFactory()
        : new Worker(new URL("./simWorker.ts", import.meta.url), { type: "module" });
    this.worker.addEventListener("message", (event: MessageEvent<WorkerToMain>) => {
      const msg = event.data;
      switch (msg.kind) {
        case "ready":
          events.onReady(msg);
          break;
        case "refusal":
          events.onRefusal(msg.stage, msg.refusal);
          break;
        case "terminal":
          this.terminalCode = msg.envelopeRefusalCode;
          events.onTerminal(msg);
          break;
        case "snapshot": {
          // postMessage fallback transport.
          this.push(decodeSnapshot(msg.tick, msg.payload));
          break;
        }
        case "metrics":
          console.info(JSON.stringify({ suite: "wf-sim-client", stage: "metrics", ...msg }));
          break;
      }
    });
  }

  /** Start a run. Uses the SAB ring iff the environment grants SAB. */
  start(scenario: ScenarioInit, runEpoch = 1): void {
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
      init = { kind: "init", scenario, sab, slots: RING_SLOTS, runEpoch };
    } else {
      this.reader = null;
      init = { kind: "init", scenario, runEpoch };
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

  sendControl(leverForceN: number, warpCmdRad: number): void {
    this.worker.postMessage({ kind: "control", leverForceN, warpCmdRad } satisfies MainToWorker);
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

  dispose(): void {
    this.worker.terminate();
  }
}
