// Shared audio plumbing for the live FrankenSim engines (drum, piano, ...).
//
// One AudioWorklet (`/live/sink-worklet.js`, processor name `fs-live-sink`)
// plays Float32 chunks exactly as an engine Worker produced them: no gain, no
// limiter, no fade, no resampling. Silence on underrun, and the underrun is
// counted. Engines render at 48 kHz, so the AudioContext must run at 48 kHz.
//
// Protocol on any producer port (see sink-worklet.js for the authoritative text):
//   producer -> sink  {type:'audio', channels: Float32Array[]}   (transfer the buffers)
//   producer -> sink  {type:'flush', epoch?: number}
//   sink -> producer  {type:'level', epoch, queuedFrames, receivedFrames, playedFrames,
//                      underruns, underrunFrames, currentFrame, currentTime, sampleRate}

/** The only sample rate the live engines render. */
export const LIVE_SAMPLE_RATE = 48_000;

/** Absolute URL prefix of the static live files (served from `public/live`). */
export const LIVE_BASE_URL = `${import.meta.env.BASE_URL}live`;

/** A typed refusal from an engine: `code` is stable and machine-readable. */
export class LiveRefusal extends Error {
  readonly code: string;
  /** The command that was refused, when the refusal answers one. */
  readonly command: string | undefined;
  constructor(code: string, message: string, command?: string) {
    super(message);
    this.name = 'LiveRefusal';
    this.code = code;
    this.command = command;
  }
}

/** Periodic report from the sink (every `levelEveryQuanta` render quanta). */
export interface LiveSinkLevel {
  type: 'level';
  /** Increments on every flush; counters below restart with it. */
  epoch: number;
  /** Frames waiting in the FIFO right now. */
  queuedFrames: number;
  /** Frames ever received (cumulative over the node's life). */
  receivedFrames: number;
  /** Frames played since the last flush. */
  playedFrames: number;
  /** Render quanta padded with silence since the last flush. */
  underruns: number;
  /** Silent frames emitted because the FIFO was empty, since the last flush. */
  underrunFrames: number;
  currentFrame: number;
  currentTime: number;
  sampleRate: number;
}

export interface LiveSinkOptions {
  /** Defaults to `/live`. */
  baseUrl?: string;
  /** Output channels of the sink node (default 2). Mono chunks are copied to every channel. */
  channels?: number;
  /** Where the sink plays (default `context.destination`). */
  destination?: AudioNode;
  /** Level report cadence in 128-frame quanta (default 2, about 5.3 ms). */
  levelEveryQuanta?: number;
}

export interface LiveSink {
  readonly context: AudioContext;
  readonly node: AudioWorkletNode;
  /**
   * A fresh MessagePort wired straight into the sink. Transfer it to an engine
   * Worker in its `init` message; audio then never touches the main thread.
   */
  createProducerPort(): MessagePort;
  /** Drop everything queued (exact silence follows). */
  flush(): void;
  /** Subscribe to level reports on the main thread. Returns an unsubscribe function. */
  onLevel(listener: (level: LiveSinkLevel) => void): () => void;
  /** `baseLatency + outputLatency` of the context, in milliseconds (0 where unreported). */
  deviceLatencyMs(): number;
  dispose(): void;
}

/**
 * Create the 48 kHz context the live engines need. Call from a user gesture.
 * A context at any other rate is refused by the engines rather than resampled.
 */
export function createLiveAudioContext(): AudioContext {
  return new AudioContext({ sampleRate: LIVE_SAMPLE_RATE, latencyHint: 'interactive' });
}

const loadedSinkModules = new WeakMap<BaseAudioContext, Promise<void>>();

/** Install the sink worklet (once per context) and create a sink node. */
export async function createLiveSink(context: AudioContext, options: LiveSinkOptions = {}): Promise<LiveSink> {
  if (context.sampleRate !== LIVE_SAMPLE_RATE) {
    throw new LiveRefusal(
      'sample_rate',
      `live engines render ${LIVE_SAMPLE_RATE} Hz only; this AudioContext runs at ${context.sampleRate} Hz. ` +
        'Create it with createLiveAudioContext(); no resampling is applied.',
    );
  }
  const baseUrl = options.baseUrl ?? LIVE_BASE_URL;
  let loading = loadedSinkModules.get(context);
  if (!loading) {
    loading = context.audioWorklet.addModule(`${baseUrl}/sink-worklet.js`);
    loadedSinkModules.set(context, loading);
  }
  await loading;

  const channels = options.channels ?? 2;
  const node = new AudioWorkletNode(context, 'fs-live-sink', {
    numberOfInputs: 0,
    numberOfOutputs: 1,
    outputChannelCount: [channels],
  });
  node.connect(options.destination ?? context.destination);
  if (options.levelEveryQuanta !== undefined) {
    node.port.postMessage({ type: 'config', levelEveryQuanta: options.levelEveryQuanta });
  }

  const listeners = new Set<(level: LiveSinkLevel) => void>();
  node.port.onmessage = (event: MessageEvent<LiveSinkLevel>) => {
    if (event.data && event.data.type === 'level') {
      for (const listener of listeners) listener(event.data);
    }
  };

  return {
    context,
    node,
    createProducerPort(): MessagePort {
      const channel = new MessageChannel();
      node.port.postMessage({ type: 'connect', port: channel.port2 }, [channel.port2]);
      return channel.port1;
    },
    flush(): void {
      node.port.postMessage({ type: 'flush' });
    },
    onLevel(listener): () => void {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    deviceLatencyMs(): number {
      const base = Number.isFinite(context.baseLatency) ? context.baseLatency : 0;
      const output = Number.isFinite(context.outputLatency) ? context.outputLatency : 0;
      return (base + output) * 1000;
    },
    dispose(): void {
      listeners.clear();
      node.port.onmessage = null;
      node.disconnect();
    },
  };
}
