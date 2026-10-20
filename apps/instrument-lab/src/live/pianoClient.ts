// Typed client for the live FrankenSim piano engine (the REDUCED REAL-TIME
// model, explicitly not the offline reference). The engine is WebAssembly in a
// module Worker (`/live/piano/worker.js`); audio flows Worker -> AudioWorklet
// sink without touching the main thread. Everything here is SI.
//
//   const context = createLiveAudioContext();              // inside a user gesture
//   const piano = await LivePianoClient.create(context);
//   piano.onSnapshot((s) => drawStrings(s.keys));
//   piano.noteOn(60, 2.5);                                  // MIDI key, hammer velocity [m/s]
//   piano.sustain(1); piano.noteOff(60);
//
// Strikes are POST-ESCAPEMENT hammer velocities in m/s. Audio samples are
// pressure [Pa] at the declared Rayleigh receiver divided by `info.fullScalePa`
// (a fixed, declared full scale: no normalization, no gain riding).

import { LIVE_BASE_URL, LiveRefusal, createLiveSink, type LiveSink } from './liveAudio';

/** MIDI velocity 127 maps to this hammer speed; linear and UNCALIBRATED (grand_piano/midi.rs). */
export const MIDI_MAX_VELOCITY_M_S = 4.5;

/** The reference's uncalibrated linear MIDI-velocity -> hammer-speed map. */
export function midiVelocityToMS(velocity: number): number {
  return (MIDI_MAX_VELOCITY_M_S * Math.max(0, Math.min(127, velocity))) / 127;
}

/** Exactly what the live model contains; produced by the engine itself. */
export interface PianoInfo {
  model: string;
  sampleRate: number;
  /** Mechanical steps per second (one per audio sample; the offline default is four). */
  mechanicalRateHz: number;
  substeps: number;
  /** Keys of the scale (88) and speaking strings attached to the bridge (234). */
  keys: number;
  speakingStrings: number;
  /** Partials held in the cold tables (all keys, all strings, up to the table ceiling). */
  tablePartials: number;
  /** Soundboard coordinates: `boardExactLowModes` exact FE modes plus bridge-response Ritz vectors. */
  boardCoordinates: number;
  boardExactLowModes: number;
  /** Lowest and highest board coordinate frequency [Hz]. */
  boardCoordinateHz: [number, number];
  /** Where the board came from and how it was reduced (the reference modules' own statement). */
  boardProvenance: string;
  /** Exact shared material gas consumed by the Rayleigh receiver. */
  acousticGas?: { temperatureK: number; pressurePa: number; densityKgM3: number; soundSpeedMS: number; gamma: number };
  /** Hammer felt law (engines built before it was reported omit it). */
  hammerLaw?: string;
  /** String partials above this frequency are not retained [Hz]. */
  partialBandHz: number;
  maxPartialsPerString: number;
  /** Dynamic string modes granted to one sounding key (all its strings and duplex segments). */
  courseModeBudget: number;
  maxSoundingKeys: number;
  /** Dynamic string modes simulated at once; beyond it the quietest key is retired (counted). */
  maxDynamicModes: number;
  /** Quadrature points of the estimated viscous damper pad per string. */
  damperPoints: number;
  duplex: boolean;
  /** Every non-sounding string stays attached to the bridge as a quasi-static chord. */
  staticAttachment: boolean;
  /** Rayleigh half-space receiver position in the board frame [m]. */
  microphoneM: [number, number, number];
  microphoneTravelSamples: number;
  microphoneFilterDelaySamples: number;
  /** Shortest acoustic travel time plus the receiver's declared filter delay [ms]. */
  microphoneLatencyMs: number;
  microphoneMultipliesPerSample: number;
  /** Pressure that maps to sample value 1.0 [Pa]. Fixed. */
  fullScalePa: number;
  midiMaxVelocityMS: number;
  /** String display stations per key in a snapshot. */
  stations: number;
  /** Relative string energy below which a sounding (undamped / damped) string returns to the static image [J]. */
  retireRingingJ: number;
  retireDampedJ: number;
  /** Silent strings: partials up to this frequency move (0: rigid chords) [Hz], and how many there are now. */
  silentBankHz?: number;
  silentBankModes?: number;
}

/** Soundboard drawing data: FE quadrature points and bare mode shapes. */
export interface PianoBoard {
  /** Interleaved x, y [m] of every surface quadrature point. */
  xy: Float32Array;
  /** Mode shapes at those points, row-major [point][mode], 1/sqrt(kg). */
  shapes: Float32Array;
  /** Bare-board mode frequencies [Hz]. */
  modeHz: Float64Array;
}

export interface PianoKeySnapshot {
  midi: number;
  held: boolean;
  /** The key's dampers are on its strings (key up and pedal not fully down). */
  damped: boolean;
  /** Hammer crown position [m]; 0 is the undeflected string plane, rest is -0.020. */
  hammerPositionM: number;
  /** Total felt force on the strings during the last step [N]. */
  hammerForceN: number;
  /** Energy of the key's dynamic strings relative to their chords [J]. */
  relativeEnergyJ: number;
  /** Partials of the displayed (first) string, 0 when only duplex segments still ring. */
  partials: number;
  /** Dynamic speaking strings of the key. */
  strings: number;
  /** Bridge displacement under the key [m]. */
  bridgeDisplacementM: number;
  /** Transverse displacement [m] of the first string at x/L = (j + 0.5)/stations. */
  stringShapeM: Float32Array;
}

/** Energy books [J]. `residualJ` is what the books do not explain. */
export interface PianoLedger {
  inputJ: number;
  storedJ: number;
  dissipatedJ: number;
  residualJ: number;
  stringLossJ: number;
  boardLossJ: number;
  damperLossJ: number;
  /** Felt (crush and relaxation) plus shank loss. */
  hammerLossJ: number;
  catchLossJ: number;
  /** Relative string energy discarded when strings returned to the static image. */
  retiredJ: number;
  /** Numerical work of the silent-string bank's explicit exchange with the board (signed). */
  couplingJ: number;
}

export interface PianoSnapshot {
  sustain: number;
  ledger: PianoLedger;
  dynamicModes: number;
  dynamicStrings: number;
  peakPa: number;
  forcedRetirements: number;
  /** Bare soundboard modal coordinates [m sqrt(kg)]; w(x, y) = sum(shapes * coordinates). */
  boardCoordinates: Float32Array;
  keys: PianoKeySnapshot[];
  renderedFrames: number;
  /** Audio already rendered ahead of the playhead when the snapshot was taken [ms]. */
  aheadMs: number;
  raw: Float32Array;
}

export interface PianoQuality {
  courseModeBudget: number;
  partialBandHz: number;
  maxDynamicModes: number;
}

export interface PianoStats {
  type: 'stats';
  /** Audio time rendered per wall time spent rendering, over the last 250 ms. */
  realtimeFactor: number | null;
  /** Smoothed render time / audio time (1 = exactly real time). */
  load: number;
  worstChunkMs: number;
  queuedMs: number;
  targetMs: number;
  underruns: number;
  underrunFrames: number;
  /** Sounding keys, dynamic string modes and dynamic speaking strings. */
  voices: number;
  dynamicModes: number;
  dynamicStrings: number;
  peakPa: number;
  /** Samples at or above full scale (the browser clips them; nothing is limited here). */
  clippedSamples: number;
  /** Keys retired early because the mode budget was full. */
  forcedRetirements: number;
  retirements: number;
  /** Strikes withdrawn because the contact solve refused at this clock. */
  contactFailures: number;
  rebuilds: number;
  noteOns: number;
  refusedNoteOns: number;
  renderedFrames: number;
  playedFrames: number;
  /** What the load governor currently grants (equals the requested quality unless it stepped down). */
  quality: PianoQuality;
  governorSteps: number;
}

export interface PianoBenchRow {
  keys: number;
  dynamicStrings: number;
  dynamicModes: number;
  realtimeFactor: number;
  worstChunkMs: number;
  chunkMs: number;
  peakFullScale: number;
  peakPa: number;
  error: string | null;
}

export interface LivePianoOptions {
  /** Defaults to LIVE_BASE_URL (`live/<hash>/` in a build). */
  baseUrl?: string;
  destination?: AudioNode;
  /**
   * 0/undefined = the complete embedded board (32 coordinates); 1..16 = that many exact low FE modes only;
   * 48/64/96/128 = that larger reduction of the same FE slice, fetched from `live/piano/boards/`;
   * 'auto' = the largest board this device computes `minRealtime` (default 1.4) times faster than real time.
   */
  boardCoordinates?: number | 'auto';
  minRealtime?: number;
  courseModeBudget?: number;
  maxDynamicModes?: number;
  /** Pressure mapped to 1.0 [Pa] (default 10). */
  fullScalePa?: number;
  /** T[K] P[Pa] RH, or mixture T P mole|mass N SPECIES FRACTION…; resolved in Rust. */
  gas?: string;
  /** Render-ahead target [ms] (default 18). */
  targetMs?: number;
  /** Snapshot cadence [Hz] (default 30; 0 disables). */
  snapshotHz?: number;
  /** Lower the mode budget automatically when rendering exceeds 80 % of real time (default true). */
  governor?: boolean;
}

interface ReadyMessage {
  info: PianoInfo;
  board: PianoBoard;
  keyGeometry: Float64Array;
  initMs: number;
}

interface Pending {
  resolve: (value: unknown) => void;
  reject: (error: Error) => void;
}

const KEY_HEADER = 9;

/** Decode the flat Float32Array snapshot posted by the worker. */
export function parsePianoSnapshot(data: Float32Array, renderedFrames = 0, aheadMs = 0): PianoSnapshot {
  const keyCount = data[0] ?? 0;
  const stations = data[1] ?? 0;
  const boardModes = data[2] ?? 0;
  // The bank adds one signed exchange value after the original 18 header values.
  const header = data.length - boardModes - keyCount * (KEY_HEADER + stations);
  if (header !== 18 && header !== 19) throw new RangeError('unsupported piano snapshot layout; rebuild the piano WASM and client together');
  const at = (index: number): number => data[index] ?? 0;
  const keys: PianoKeySnapshot[] = [];
  let offset = header + boardModes;
  for (let k = 0; k < keyCount; k++) {
    keys.push({
      midi: at(offset),
      held: at(offset + 1) > 0,
      damped: at(offset + 2) > 0,
      hammerPositionM: at(offset + 3),
      hammerForceN: at(offset + 4),
      relativeEnergyJ: at(offset + 5),
      partials: at(offset + 6),
      strings: at(offset + 7),
      bridgeDisplacementM: at(offset + 8),
      stringShapeM: data.subarray(offset + KEY_HEADER, offset + KEY_HEADER + stations),
    });
    offset += KEY_HEADER + stations;
  }
  return {
    sustain: at(3),
    ledger: {
      inputJ: at(4),
      storedJ: at(5),
      dissipatedJ: at(6),
      residualJ: at(7),
      stringLossJ: at(8),
      boardLossJ: at(9),
      damperLossJ: at(10),
      hammerLossJ: at(11),
      catchLossJ: at(12),
      retiredJ: at(13),
      couplingJ: header === 19 ? at(18) : 0,
    },
    dynamicModes: at(14),
    dynamicStrings: at(15),
    peakPa: at(16),
    forcedRetirements: at(17),
    boardCoordinates: data.subarray(header, header + boardModes),
    keys,
    renderedFrames,
    aheadMs,
    raw: data,
  };
}

export class LivePianoClient {
  readonly info: PianoInfo;
  readonly board: PianoBoard;
  /** Per key: midi, speaking length [m], strike station (fraction of length), unison strings, first partial [Hz]. */
  readonly keyGeometry: Float64Array;
  /** Time the engine took to assemble the instrument in the worker [ms]. */
  readonly initMs: number;

  private readonly worker: Worker;
  private readonly sink: LiveSink;
  private readonly pending = new Map<number, Pending>();
  private readonly snapshotListeners = new Set<(snapshot: PianoSnapshot) => void>();
  private readonly statsListeners = new Set<(stats: PianoStats) => void>();
  private readonly errorListeners = new Set<(error: LiveRefusal) => void>();
  private nextId = 1;
  private lastStats: PianoStats | null = null;
  private disposed = false;

  private constructor(worker: Worker, sink: LiveSink, ready: ReadyMessage) {
    this.worker = worker;
    this.sink = sink;
    this.info = ready.info;
    this.board = ready.board;
    this.keyGeometry = ready.keyGeometry;
    this.initMs = ready.initMs;
    worker.onmessage = (event: MessageEvent) => this.receive(event.data);
    worker.onerror = (event) => this.workerFailed(event);
  }

  /**
   * Load the WebAssembly engine, assemble the Model D board and the 88-key
   * string tables (about half a second) and start the silent stream.
   * `audioContext` must run at 48 kHz (see `createLiveAudioContext`).
   */
  static async create(audioContext: AudioContext, options: LivePianoOptions = {}): Promise<LivePianoClient> {
    const baseUrl = options.baseUrl ?? LIVE_BASE_URL;
    const sink = await createLiveSink(audioContext, options.destination ? { baseUrl, destination: options.destination } : { baseUrl });
    const worker = new Worker(`${baseUrl}/piano/worker.js`, { type: 'module', name: 'fs-live-piano' });
    const port = sink.createProducerPort();
    let rejectStartup!: (error: LiveRefusal) => void;
    const failed = new Promise<never>((_resolve, reject) => { rejectStartup = reject; });
    const first = new Promise<ReadyMessage>((resolve) => {
      worker.onmessage = (event: MessageEvent) => {
        const message = event.data;
        if (message?.type === 'ready') {
          resolve({ info: message.info, board: message.board, keyGeometry: message.keyGeometry ?? new Float64Array(0), initMs: message.initMs });
        }
        else if (message?.type === 'error') rejectStartup(new LiveRefusal(message.code ?? 'worker_error', message.message, message.command));
      };
      worker.onerror = (event) => {
        event.preventDefault();
        rejectStartup(new LiveRefusal('worker_error', event.message || 'the piano worker failed to load'));
      };
    });
    const { baseUrl: _baseUrl, destination: _destination, ...engineOptions } = options;
    worker.postMessage({ type: 'init', sampleRate: audioContext.sampleRate, port, options: engineOptions }, [port]);
    let ready: ReadyMessage;
    try {
      ready = await Promise.race([first, failed]);
      // Ready settles before the final device resume; a worker fault must still
      // abort startup promptly while that second asynchronous step is pending.
      if (audioContext.state === 'suspended') await Promise.race([audioContext.resume(), failed]);
    } catch (error) {
      worker.terminate();
      sink.dispose();
      throw error;
    }
    return new LivePianoClient(worker, sink, ready);
  }

  private workerFailed(event: ErrorEvent): void {
    event.preventDefault();
    if (this.disposed) return;
    const refusal = new LiveRefusal('worker_failed', event.message || 'the piano worker stopped');
    const listeners = [...this.errorListeners];
    for (const pending of this.pending.values()) pending.reject(refusal);
    this.pending.clear();
    this.dispose();
    // A render fault can recur on every timer turn; stop the failed worker now.
    this.worker.terminate();
    for (const listener of listeners) listener(refusal);
  }

  private receive(message: { type?: string; [key: string]: unknown }): void {
    switch (message?.type) {
      case 'snapshot':
        if (this.snapshotListeners.size > 0) {
          const snapshot = parsePianoSnapshot(message.data as Float32Array, message.renderedFrames as number, message.aheadMs as number);
          for (const listener of this.snapshotListeners) listener(snapshot);
        }
        break;
      case 'stats':
        this.lastStats = message as unknown as PianoStats;
        for (const listener of this.statsListeners) listener(this.lastStats);
        break;
      case 'reply':
      case 'bench': {
        const pending = this.pending.get(message.id as number);
        if (pending) {
          this.pending.delete(message.id as number);
          pending.resolve(message.type === 'bench' ? message.results : message.value);
        }
        break;
      }
      case 'error': {
        const refusal = new LiveRefusal((message.code as string) ?? 'worker_error', message.message as string, message.command as string | undefined);
        const pending = message.id !== undefined ? this.pending.get(message.id as number) : undefined;
        if (pending) {
          this.pending.delete(message.id as number);
          pending.reject(refusal);
        }
        for (const listener of this.errorListeners) listener(refusal);
        break;
      }
      default:
        break;
    }
  }

  private send(name: string, ...args: unknown[]): void {
    if (!this.disposed) this.worker.postMessage({ type: 'cmd', name, args });
  }

  private request<T>(name: string, ...args: unknown[]): Promise<T> {
    if (this.disposed) return Promise.reject(new LiveRefusal('disposed', 'the piano client is disposed', name));
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject });
      this.worker.postMessage({ type: 'cmd', name, args, id });
    });
  }

  /** Strike a key (MIDI 21..108) with a post-escapement hammer velocity in m/s, (0, 8]. */
  noteOn(midi: number, speedMS: number): void {
    this.send('noteOn', midi, speedMS);
  }

  /** Strike with a MIDI velocity through the reference's uncalibrated linear map (127 -> 4.5 m/s). */
  noteOnMidi(midi: number, velocity: number): void {
    if (velocity > 0) this.send('noteOn', midi, midiVelocityToMS(velocity));
    else this.send('noteOff', midi);
  }

  /** Release the key: its dampers fall back unless the pedal holds them. */
  noteOff(midi: number): void {
    this.send('noteOff', midi);
  }

  /** Sustain pedal travel 0..1; damper drag scales with (1 - travel)^2. */
  sustain(value: number): void {
    this.send('sustain', Math.max(0, Math.min(1, Number(value))));
  }

  /** Sostenuto pedal: pressing it latches the dampers of the keys held down then. */
  sostenuto(on: boolean): void {
    this.send('sostenuto', on);
  }

  /** Una corda pedal: strikes from now on miss the last string of each course. */
  unaCorda(on: boolean): void {
    this.send('unaCorda', on);
  }

  /** Hold a key down without a strike: its dampers lift and its strings ring along. Release with noteOff. */
  silentKeyDown(midi: number): void {
    this.send('silentKeyDown', midi);
  }

  allNotesOff(): void {
    this.send('allNotesOff');
  }

  /** Limits for keys struck afterwards. Also resets what the governor grants. */
  setQuality(quality: Partial<PianoQuality>): void {
    this.send('setQuality', quality);
  }

  setTargetMs(ms: number): void {
    this.send('setTargetMs', ms);
  }

  setSnapshotHz(hz: number): void {
    this.send('setSnapshotHz', hz);
  }

  setGovernor(on: boolean): void {
    this.send('setGovernor', on);
  }

  /** Stiff-string partial frequencies [Hz] the live model retains for a key. */
  partials(midi: number): Promise<Float64Array> {
    return this.request<Float64Array>('partials', midi);
  }

  /** Offline speed test on a separate engine instance: real-time factor versus sounding keys. */
  bench(options: { keys?: number[]; counts?: number[]; seconds?: number; speedMS?: number; sustain?: number } = {}): Promise<PianoBenchRow[]> {
    return this.request<PianoBenchRow[]>('bench', options);
  }

  onSnapshot(listener: (snapshot: PianoSnapshot) => void): () => void {
    this.snapshotListeners.add(listener);
    return () => this.snapshotListeners.delete(listener);
  }

  onStats(listener: (stats: PianoStats) => void): () => void {
    this.statsListeners.add(listener);
    return () => this.statsListeners.delete(listener);
  }

  onError(listener: (error: LiveRefusal) => void): () => void {
    this.errorListeners.add(listener);
    return () => this.errorListeners.delete(listener);
  }

  /** Latest stats message, if any. */
  get stats(): PianoStats | null {
    return this.lastStats;
  }

  /**
   * Estimated strike-to-ear latency [ms]: audio queued ahead of the playhead,
   * the receiver's acoustic travel and filter delay, and the output device.
   * The hammer's own flight to the string (2 mm / speed) comes on top.
   */
  latencyEstimateMs(): { queuedMs: number; receiverMs: number; outputMs: number; totalMs: number } {
    const queuedMs = this.lastStats ? this.lastStats.queuedMs : 0;
    const receiverMs = this.info.microphoneLatencyMs;
    const outputMs = this.sink.deviceLatencyMs();
    return { queuedMs, receiverMs, outputMs, totalMs: queuedMs + receiverMs + outputMs };
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.worker.onmessage = null;
    this.worker.onerror = null;
    try {
      this.worker.postMessage({ type: 'cmd', name: 'dispose' });
    } catch {
      // the worker is already gone
    }
    setTimeout(() => this.worker.terminate(), 200);
    this.sink.flush();
    this.sink.dispose();
    for (const pending of this.pending.values()) pending.reject(new LiveRefusal('disposed', 'the piano client is disposed'));
    this.pending.clear();
    this.snapshotListeners.clear();
    this.statsListeners.clear();
    this.errorListeners.clear();
  }
}
