// Client for `fs-live-mix` (public/live/mix-worklet.js): several live engine Workers, each a
// lane whose audio is stamped with the context frame it plays at, summed in lock-step with the
// AudioContext clock. See public/live/lane.js for the Worker side.

import { LIVE_BASE_URL, LIVE_SAMPLE_RATE, LiveRefusal, loadWorkletModule } from "./liveAudio";

export interface LaneStats {
  queuedFrames: number;
  /** Frames dropped because they arrived after they played (whole messages and partly late heads). */
  lateFrames: number;
  /**
   * Frames a streaming lane skipped (it fell behind and stepped ahead: its skippedFrames). A lane
   * resting between strokes is not a gap. No frame is both, so late + gap is the audio lost.
   */
  gapFrames: number;
}

export interface LiveMix {
  readonly context: AudioContext;
  readonly node: AudioWorkletNode;
  /** Join a lane: the mixer starts sending it clock ticks (attach once its Worker is ready). */
  attach(id: string, port: MessagePort): void;
  dropLane(id: string): void;
  /**
   * Discard queued audio at context frames >= `from` (default: now) of lane `id`, or of every
   * lane; audio stamped with an epoch below `epoch` is refused from then on (see lane resync).
   */
  flush(opts?: { id?: string; from?: number; epoch?: number }): void;
  onStats(listener: (frame: number, lanes: Record<string, LaneStats>) => void): () => void;
  dispose(): void;
}

export async function createLiveMix(context: AudioContext, destination: AudioNode, channels = 2): Promise<LiveMix> {
  if (context.sampleRate !== LIVE_SAMPLE_RATE) {
    throw new LiveRefusal("sample_rate", `live engines render ${LIVE_SAMPLE_RATE} Hz only; this AudioContext runs at ${context.sampleRate} Hz.`);
  }
  await loadWorkletModule(context, `${LIVE_BASE_URL}/mix-worklet.js`).catch((e: unknown) => {
    throw new Error(`the live mixer could not be loaded (if the site was just updated, reload the page): ${(e as Error)?.message ?? e}`);
  });
  const node = new AudioWorkletNode(context, "fs-live-mix", { numberOfInputs: 0, numberOfOutputs: 1, outputChannelCount: [channels] });
  node.connect(destination);
  const listeners = new Set<(frame: number, lanes: Record<string, LaneStats>) => void>();
  node.port.onmessage = (e: MessageEvent<{ type: string; frame: number; lanes: Record<string, LaneStats> }>) => {
    if (e.data?.type === "stats") for (const l of listeners) l(e.data.frame, e.data.lanes);
  };
  return {
    context,
    node,
    attach(id: string, port: MessagePort): void {
      node.port.postMessage({ type: "lane", id, port }, [port]);
    },
    dropLane(id: string): void {
      node.port.postMessage({ type: "drop", id });
    },
    flush(opts = {}): void {
      node.port.postMessage({ type: "flush", ...opts });
    },
    onStats(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    dispose(): void {
      listeners.clear();
      node.port.postMessage({ type: "dispose" }); // the worklet closes every lane port and stops processing
      node.port.onmessage = null;
      node.disconnect();
    },
  };
}

/** The error a lane's `ready` rejects with when the lane is disposed before its Worker is ready. */
export class LaneDisposed extends Error {
  constructor(label: string) {
    super(`${label} was stopped before it was ready`);
    this.name = "LaneDisposed";
  }
}

/**
 * The Worker side of a lane: the Worker, the channel whose port2 joins the mix once the Worker is
 * ready, and a `ready` that settles exactly once. Disposing it terminates the Worker, drops the mix
 * lane, closes the channel and rejects a `ready` still pending (with LaneDisposed), so nothing waits
 * forever on a Worker that is gone and a `ready` that was already in flight is ignored.
 */
class LaneWorker<R> {
  readonly worker: Worker;
  readonly ready: Promise<R>;
  /**
   * Called once when the Worker fails after it was ready (a wasm trap mid-song: every later turn would
   * throw again, about every 10 ms, while the lane stayed silent). The lane is disposed first.
   */
  onFailure: ((e: Error) => void) | null = null;
  private readonly channel = new MessageChannel();
  private settled = false;
  private attached = false;
  private disposed = false;
  private rejectReady: (e: Error) => void = () => {};

  /** `name`: the Worker's name (devtools); `label`: how errors name it ("the snare engine"). */
  constructor(url: string, name: string, private readonly label: string, private readonly mix: LiveMix, private readonly id: string, onMessage: (m: { type: string; [k: string]: unknown }) => void) {
    try { this.worker = new Worker(url, { type: "module", name }); }
    catch (e) {
      this.channel.port1.close();
      this.channel.port2.close();
      throw e;
    }
    this.ready = new Promise<R>((resolve, reject) => {
      this.rejectReady = reject;
      this.worker.onmessage = (e: MessageEvent) => {
        const m = e.data;
        if (this.disposed || !m) return;
        if (m.type === "ready") {
          if (this.settled) return;
          try { mix.attach(id, this.channel.port2); }
          catch (e) { this.fail(e instanceof Error ? e : new Error(String(e))); return; }
          this.settled = true;
          this.attached = true;
          resolve(m as R);
        } else if (m.type === "error") {
          this.fail(new Error(String(m.message)));
        } else onMessage(m);
      };
      // a module Worker that fails to load (a 404 after the site was updated, say) fires a bare Event
      this.worker.onerror = (e) => {
        e.preventDefault();
        this.fail(new Error(e.message || (this.settled ? `${label} failed` : `${label} could not be loaded (if the site was just updated, reload the page)`)));
      };
    });
    this.ready.catch(() => {}); // a caller that awaits it sees the rejection; nobody else should
  }

  /** If init cannot be transferred, construction leaves no caller able to dispose this lane. */
  initialize(options: Record<string, unknown>): void {
    const port = this.channel.port1;
    try { this.worker.postMessage({ type: "init", ...options, port }, [port]); }
    catch (e) { this.dispose(); throw e; }
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.worker.onmessage = null;
    this.worker.onerror = null;
    this.worker.terminate();
    this.channel.port1.close();
    // only a lane that joined the mix leaves it: a new lane may already be joining under the same id
    if (this.attached) this.mix.dropLane(this.id);
    else this.channel.port2.close();
    if (!this.settled) { this.settled = true; this.rejectReady(new LaneDisposed(this.label)); }
  }

  private fail(e: Error): void {
    if (this.disposed) return;
    const report = this.settled ? this.onFailure : null;
    this.onFailure = null;
    if (!this.settled) { this.settled = true; this.rejectReady(e); }
    this.dispose();
    report?.(e);
  }
}

// ------------------------------------------------------------------ kit piece lanes
export interface KitReady {
  info: { piece: string; gridScale: number; radius: number; surfaceNodes: number; h: number; bodies: { name: string; nodes: number }[];
    gas?: { temperatureK: number; pressurePa: number; densityKgM3: number; soundSpeedMps: number } };
  /** Exact requested atmosphere, retained through every grid trial and rebuild. */
  gasRecord?: string;
  xy: Float32Array;
  h: number;
  radius: number;
  scale: number;
  realtime: number;
  /** Factor on the page's gain (set for the 5x grid) that keeps a stroke as loud on `scale`; 1 on 5x. */
  levelFix?: number;
  /** "none": no level table for this wasm, so the gain is unscaled. */
  levels?: "ok" | "none";
}
export interface KitField {
  frame: number;
  /** The lane epoch it was computed in (see KitLane.resync). */
  epoch?: number;
  v: Float32Array;
  kinetic: number;
  force: number;
}
/** The surface every `stepMs` for the first milliseconds after a stroke, `count` arrays back to back. */
export interface KitBurst {
  frame: number;
  stepMs: number;
  count: number;
  v: Float32Array;
}
export interface LaneWorkerStats {
  realtime: number;
  renderedFrames: number;
  skippedFrames: number;
  lateEvents: number;
  asleep?: boolean;
  voices?: number;
  peakPa?: number;
  dynamicModes?: number;
  noteOns?: number;
  refusedNoteOns?: number;
  forcedRetirements?: number;
  /** Piano: the dynamic string mode cap the load governor grants now, and the one asked for. */
  modeCap?: number;
  fullModeCap?: number;
  courseModeBudget?: number;
  governorDown?: number;
  governorUp?: number;
}

/** A stroke of a fixed score: absolute context frame, point on the surface [m], stick speed [m/s]. */
export interface KitScoreEvent { frame: number; x: number; y: number; speed: number }
/** A key or pedal event of a fixed score (one operation per event; silent lifts a key without a strike). */
export interface PianoScoreEvent { frame: number; on?: number; off?: number; silent?: number; speed?: number; sustain?: number; unaCorda?: boolean; sostenuto?: boolean }

/** One kit piece computed in its own Worker, playing into a lane of the mix. */
export class KitLane {
  readonly piece: string;
  readonly ready: Promise<KitReady>;
  onField: ((f: KitField) => void) | null = null;
  onBurst: ((b: KitBurst) => void) | null = null;
  onStats: ((s: LaneWorkerStats) => void) | null = null;
  private readonly lane: LaneWorker<KitReady>;
  private readonly worker: Worker;

  /**
   * `chunk`: most frames per audio message (default 128; a score lane far ahead can send larger ones, see lane.js).
   * `benchMs`: how much of one stroke each grid is timed over (default 112 ms; see kit/worker.js).
   * `gasRecord`: one canonical gas-state record retained through grid trials and rebuilding; an
   * explicit record requires an engine with material-atmosphere support.
   */
  constructor(mix: LiveMix, piece: string, opts: { gain: number; gasRecord?: string; aheadFrames?: number; chunk?: number; scales?: number[]; minRealtime?: number; benchMs?: number; fieldHz?: number; fieldFrames?: number; bursts?: boolean }) {
    this.piece = piece;
    this.lane = new LaneWorker<KitReady>(`${LIVE_BASE_URL}/kit/worker.js`, `kit-${piece}`, `the ${piece} engine`, mix, piece, (m) => {
      if (m.type === "field") this.onField?.(m as unknown as KitField);
      else if (m.type === "burst") this.onBurst?.(m as unknown as KitBurst);
      else if (m.type === "stats") this.onStats?.(m as unknown as LaneWorkerStats);
    });
    this.worker = this.lane.worker;
    this.ready = this.lane.ready;
    this.lane.initialize({ piece, ...opts });
  }
  strike(x: number, y: number, speed: number, frame?: number): void {
    this.worker.postMessage({ type: "strike", x, y, speed, frame });
  }
  sequence(startFrame: number, loopFrames: number, events: { offset: number; x: number; y: number; speed: number }[]): void {
    this.worker.postMessage({ type: "sequence", startFrame, loopFrames, events });
  }
  /** Replace the fixed score (played once, sorted once in the Worker). */
  score(events: KitScoreEvent[]): void {
    this.worker.postMessage({ type: "score", events });
  }
  stopScore(): void {
    this.worker.postMessage({ type: "stopScore" });
  }
  /** After `mix.flush({ id: piece, epoch })`: render for the clock as it is now. */
  resync(epoch: number): void {
    this.worker.postMessage({ type: "resync", epoch });
  }
  stop(): void {
    this.worker.postMessage({ type: "stop" });
  }
  /** Called once if the engine fails after it was ready; the lane has already left the mix (LaneWorker.onFailure). */
  set onFailure(f: ((e: Error) => void) | null) {
    this.lane.onFailure = f;
  }
  /** Terminate the Worker and leave the mix; a `ready` still pending rejects with LaneDisposed. */
  dispose(): void {
    this.onField = this.onBurst = null;
    this.onStats = null;
    this.lane.dispose();
  }
}

export interface PianoReady {
  info: Record<string, unknown>;
  boardCoordinates: number;
  /** The startup probe's speed (a sustained, damped workload) at the chosen board and mode cap. */
  realtime: number;
  /** The dynamic string mode cap the piano starts with (chosen with the board for that headroom). */
  modeCap?: number;
  /** Present when snapshots were requested: the board's quadrature points and mode shapes. */
  board?: { xy: Float32Array; shapes: Float32Array; modeHz: Float64Array };
  keyGeometry?: Float64Array;
}

/** The live piano computed in its own Worker, playing into a lane of the mix. */
export class PianoLane {
  readonly ready: Promise<PianoReady>;
  onStats: ((s: LaneWorkerStats) => void) | null = null;
  /** engine.snapshot() for the state at context frame `frame` (only with `snapshotHz`), and its lane epoch. */
  onSnapshot: ((frame: number, data: Float32Array, epoch: number) => void) | null = null;
  private readonly lane: LaneWorker<PianoReady>;
  private readonly worker: Worker;

  /** `chunk`: most frames per audio message (default 128; a score lane far ahead can send larger ones, see lane.js). */
  constructor(mix: LiveMix, opts: { aheadFrames?: number; chunk?: number; options?: Record<string, number>; gain?: number; minRealtime?: number; snapshotHz?: number; governor?: boolean } = {}) {
    this.lane = new LaneWorker<PianoReady>(`${LIVE_BASE_URL}/piano/lane-worker.js`, "piano-lane", "the piano engine", mix, "piano", (m) => {
      if (m.type === "snapshot") this.onSnapshot?.(m.frame as number, m.data as Float32Array, (m.epoch as number) ?? 0);
      else if (m.type === "stats") this.onStats?.(m as unknown as LaneWorkerStats);
    });
    this.worker = this.lane.worker;
    this.ready = this.lane.ready;
    this.lane.initialize(opts);
  }
  sequence(startFrame: number, loopFrames: number, events: { offset: number; on?: number; off?: number; speed?: number }[]): void {
    this.worker.postMessage({ type: "sequence", startFrame, loopFrames, events });
  }
  /** Replace the fixed score (played once, sorted once in the Worker). */
  score(events: PianoScoreEvent[]): void {
    this.worker.postMessage({ type: "score", events });
  }
  /** Drop the rest of the score; every key and the pedal are released. */
  stopScore(): void {
    this.worker.postMessage({ type: "stopScore" });
  }
  /** After `mix.flush({ id: "piano", epoch })`: render for the clock as it is now. */
  resync(epoch: number): void {
    this.worker.postMessage({ type: "resync", epoch });
  }
  stop(): void {
    this.worker.postMessage({ type: "stop" });
  }
  /** Called once if the engine fails after it was ready; the lane has already left the mix (LaneWorker.onFailure). */
  set onFailure(f: ((e: Error) => void) | null) {
    this.lane.onFailure = f;
  }
  /** Terminate the Worker and leave the mix; a `ready` still pending rejects with LaneDisposed. */
  dispose(): void {
    this.onSnapshot = null;
    this.onStats = null;
    this.lane.dispose();
  }
}
