// Client for `fs-live-mix` (public/live/mix-worklet.js): several live engine Workers, each a
// lane whose audio is stamped with the context frame it plays at, summed in lock-step with the
// AudioContext clock. See public/live/lane.js for the Worker side.

import { LIVE_BASE_URL, LIVE_SAMPLE_RATE, LiveRefusal } from "./liveAudio";

export interface LaneStats {
  queuedFrames: number;
  lateFrames: number;
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

const loaded = new WeakMap<BaseAudioContext, Promise<void>>();

export async function createLiveMix(context: AudioContext, destination: AudioNode, channels = 2): Promise<LiveMix> {
  if (context.sampleRate !== LIVE_SAMPLE_RATE) {
    throw new LiveRefusal("sample_rate", `live engines render ${LIVE_SAMPLE_RATE} Hz only; this AudioContext runs at ${context.sampleRate} Hz.`);
  }
  let loading = loaded.get(context);
  if (!loading) {
    loading = context.audioWorklet.addModule(`${LIVE_BASE_URL}/mix-worklet.js`);
    loaded.set(context, loading);
  }
  await loading;
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
      node.port.onmessage = null;
      node.disconnect();
    },
  };
}

// ------------------------------------------------------------------ kit piece lanes
export interface KitReady {
  info: { piece: string; gridScale: number; radius: number; surfaceNodes: number; h: number; bodies: { name: string; nodes: number }[] };
  xy: Float32Array;
  h: number;
  radius: number;
  scale: number;
  realtime: number;
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
}

/** A stroke of a fixed score: absolute context frame, point on the surface [m], stick speed [m/s]. */
export interface KitScoreEvent { frame: number; x: number; y: number; speed: number }
/** A key event of a fixed score (exactly one of on / off / sustain). */
export interface PianoScoreEvent { frame: number; on?: number; off?: number; speed?: number; sustain?: number }

/** One kit piece computed in its own Worker, playing into a lane of the mix. */
export class KitLane {
  readonly piece: string;
  readonly ready: Promise<KitReady>;
  onField: ((f: KitField) => void) | null = null;
  onBurst: ((b: KitBurst) => void) | null = null;
  onStats: ((s: LaneWorkerStats) => void) | null = null;
  private worker: Worker;
  private mix: LiveMix;

  constructor(mix: LiveMix, piece: string, opts: { gain: number; aheadFrames?: number; scales?: number[]; minRealtime?: number; fieldHz?: number; fieldFrames?: number; bursts?: boolean }) {
    this.piece = piece;
    this.mix = mix;
    this.worker = new Worker(`${LIVE_BASE_URL}/kit/worker.js`, { type: "module", name: `kit-${piece}` });
    const ch = new MessageChannel();
    const port = ch.port1;
    this.ready = new Promise<KitReady>((resolve, reject) => {
      this.worker.onmessage = (e: MessageEvent) => {
        const m = e.data;
        if (m.type === "ready") { mix.attach(piece, ch.port2); resolve(m as KitReady); }
        else if (m.type === "field") this.onField?.(m as KitField);
        else if (m.type === "burst") this.onBurst?.(m as KitBurst);
        else if (m.type === "stats") this.onStats?.(m as LaneWorkerStats);
        else if (m.type === "error") reject(new Error(m.message));
      };
      this.worker.onerror = (e) => reject(new Error(e.message || "kit worker failed"));
    });
    this.worker.postMessage({ type: "init", piece, port, ...opts }, [port]);
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
  dispose(): void {
    this.worker.terminate();
    this.mix.dropLane(this.piece);
  }
}

export interface PianoReady {
  info: Record<string, unknown>;
  boardCoordinates: number;
  realtime: number;
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
  private worker: Worker;
  private mix: LiveMix;

  constructor(mix: LiveMix, opts: { aheadFrames?: number; options?: Record<string, number>; gain?: number; minRealtime?: number; snapshotHz?: number } = {}) {
    this.mix = mix;
    this.worker = new Worker(`${LIVE_BASE_URL}/piano/lane-worker.js`, { type: "module", name: "piano-lane" });
    const ch = new MessageChannel();
    const port = ch.port1;
    this.ready = new Promise((resolve, reject) => {
      this.worker.onmessage = (e: MessageEvent) => {
        const m = e.data;
        if (m.type === "ready") { mix.attach("piano", ch.port2); resolve(m as PianoReady); }
        else if (m.type === "snapshot") this.onSnapshot?.(m.frame as number, m.data as Float32Array, (m.epoch as number) ?? 0);
        else if (m.type === "stats") this.onStats?.(m as LaneWorkerStats);
        else if (m.type === "error") reject(new Error(m.message));
      };
      this.worker.onerror = (e) => reject(new Error(e.message || "piano worker failed"));
    });
    this.worker.postMessage({ type: "init", port, ...opts }, [port]);
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
  dispose(): void {
    this.worker.terminate();
    this.mix.dropLane("piano");
  }
}
