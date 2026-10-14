// Typed client for the live FrankenSim drum engine (the reduced real-time
// model). The engine itself is WebAssembly in a module Worker
// (`/live/drum/worker.js`); audio flows Worker -> AudioWorklet sink without
// touching the main thread. Everything here is SI.
//
//   const context = createLiveAudioContext();              // inside a user gesture
//   const drum = await LiveDrumClient.create(context);
//   drum.onSnapshot((s) => drawHead(s.batterDisplacementM));
//   await drum.strike(0.06, 0.01, 2.0);                     // x [m], y [m], tip speed [m/s]
//   await drum.configure({ batter_tension_n_m: 4200 });     // rebuilds from geometry, restarts silent

import { LIVE_BASE_URL, LiveRefusal, createLiveSink, type LiveSink } from './liveAudio';

/** Every declared input of the live drum. Flat, SI, no pitch or gain parameter. */
export interface DrumParams {
  /** Clear vibrating head radius [m] (0.02..1). */
  radius_m: number;
  /** Shell depth = head separation [m] (0.01..2); sets the sealed air volume. */
  depth_m: number;
  batter_thickness_m: number;
  resonant_thickness_m: number;
  /** Installed membrane tension [N/m]. */
  batter_tension_n_m: number;
  resonant_tension_n_m: number;
  /** Viscous modal damping ratio, in [0, 1). */
  batter_damping_ratio: number;
  resonant_damping_ratio: number;
  /** Equilibrated tension variation: xx, yy, xy [N/m], then c, d, e, f [N/m^2]. Zeros = uniform. */
  batter_tension_variation: [number, number, number, number, number, number, number];
  resonant_tension_variation: [number, number, number, number, number, number, number];
  /** Film material (both heads). */
  young_pa: number;
  poisson: number;
  density_kg_m3: number;
  /** FE mesh of each head: 1..32 radial intervals, 8..128 vertices per ring. */
  radial_intervals: number;
  azimuths: number;
  /** Retained head eigen-window [Hz]. Modes outside it are absent, not approximated. */
  band_low_hz: number;
  band_high_hz: number;
  /** Mechanical steps per 48 kHz frame: 16 (768 kHz, the offline rate), 8, 4, 2 or 1. */
  substeps: 1 | 2 | 4 | 8 | 16;
  gas_density_kg_m3: number;
  gas_sound_speed_m_s: number;
  /** Microphone above the batter head: head-centred x, y and height z >= 0.05 [m]. */
  mic_batter_m: [number, number, number];
  /** Microphone below the resonant head: x, y and distance below it z >= 0.05 [m]. */
  mic_resonant_m: [number, number, number];
  /** Pressure that maps to sample value 1.0 [Pa]. Fixed; never adaptive. */
  full_scale_pa: number;
  /** Hunt-Crossley contact loss chi [s/m]; 0 = the reference's elastic Hertz tip. */
  contact_loss_s_m: number;
  /** Separation at which the departing stick is caught [m]. */
  catch_gap_m: number;
}

/** Write-only convenience accepted by `configure`: sets both heads' damping ratio. */
export type DrumParamPatch = Partial<DrumParams> & { damping_ratio?: number };

/**
 * Mesh / eigen-window presets. Only those four numbers differ; everything
 * else is the same declared drum. Measured in headless Chrome 155 on an Apple
 * M4 with the hardware-FMA build (auto-calibrated rate, 20 Hz roll scenario).
 */
export const DRUM_PRESETS = {
  /** The offline reference's own head cards: 5 x 32 mesh, 80-500 Hz, 6 + 3 modes. Calibrated to 384 kHz (2.3x). */
  reference: { radial_intervals: 5, azimuths: 32, band_low_hz: 80, band_high_hz: 500 },
  /** 10 x 48 mesh, 60-800 Hz: 17 + 10 modes. Calibrated to 192 kHz (1.8x); 110 ms rebuilds. */
  wide: { radial_intervals: 10, azimuths: 48, band_low_hz: 60, band_high_hz: 800 },
  /** 12 x 64 mesh, 60-1000 Hz: 30 + 17 modes. Calibrated to 96 kHz (2.0x); 0.4 s rebuilds. */
  wider: { radial_intervals: 12, azimuths: 64, band_low_hz: 60, band_high_hz: 1000 },
} as const satisfies Record<string, DrumParamPatch>;

export interface DrumHeadInfo {
  name: 'batter' | 'resonant';
  modes: number;
  /** In-vacuo eigenfrequencies of this head alone [Hz], ascending. */
  frequencies_hz: number[];
  film_mass_kg: number;
  modal_areas_m2_per_sqrt_kg: number[];
  damping_per_s: number[];
  max_eigen_residual: number;
}

export interface DrumInfo {
  model: 'reduced real-time model';
  engine: string;
  /** Which fused multiply-add is executing (results are bit-identical on every path). */
  fma: string;
  sample_rate_hz: number;
  substeps: number;
  mechanical_rate_hz: number;
  params: DrumParams;
  heads: [DrumHeadInfo, DrumHeadInfo];
  /** Small-signal frequencies of both heads joined by the sealed-air spring [Hz]. */
  coupled_frequencies_hz: number[];
  mesh: { nodes: number; triangles: number; radial_intervals: number; azimuths: number };
  stick: {
    full_mass_kg: number;
    effective_mass_kg: number;
    hertz_stiffness_n_m15: number;
    hertz_exponent: number;
    contact_loss_s_m: number;
    launch_gap_m: number;
    catch_gap_m: number;
  };
  cavity: { volume_m3: number; bulk_modulus_pa: number };
  radiation: {
    model: string;
    full_scale_pa: number;
    /** Anti-alias (decimator) latency per microphone, in 48 kHz frames. */
    filter_delay_frames: [number, number];
    /** Acoustic travel time (nearest, farthest head point) per microphone, in frames. Physical. */
    travel_delay_frames: [[number, number], [number, number]];
    multiply_adds_per_frame: number;
  };
  snapshot_layout: { header: number; fields: string[] };
  builds: number;
  /** Text the site may show verbatim. */
  provenance: string;
  no_claim: string;
}

export interface DrumCalibrationResult {
  substeps: number;
  measured: boolean;
  /** Audio seconds per wall second for the scenario; >= 1 is real time. */
  realtimeFactor?: number;
  predictedRealtimeFactor?: number;
  wallMs?: number;
  nsPerStep?: number;
  mechanicalSteps?: number;
  stepsWithStick?: number;
  stepsInContact?: number;
  strikes?: number;
  peakPa?: number;
  balanceResidualJ?: number;
  refused?: string;
}

export interface DrumCalibration {
  results: DrumCalibrationResult[];
  /** Highest substeps whose measured factor met `minFactor` (lowest candidate if none did). */
  chosen: number | null;
  minFactor: number;
  /** False when no candidate met `minFactor`: expect underruns. */
  sustained: boolean;
  /** False for `bench()` (measure only). */
  applied?: boolean;
  scenario?: { frames: number; strikeIntervalFrames: number; x: number; y: number; speed: number };
}

export interface DrumCalibrationOptions {
  /** Default [16, 8, 4, 2, 1]. */
  candidates?: number[];
  /** Required real-time factor (default 1.6). */
  minFactor?: number;
  /** Audio frames per measurement (default 9600 = 0.2 s). */
  frames?: number;
  /** A stroke every this many frames (default 2400 = a 20 Hz roll; 0 = one stroke). */
  strikeIntervalFrames?: number;
  x?: number;
  y?: number;
  speed?: number;
}

export interface DrumReady {
  info: DrumInfo;
  /** Head mesh shared by both heads. nodes = [x0, y0, x1, y1, ...] in metres. */
  mesh: { nodes: Float32Array; triangles: Uint32Array };
  /**
   * Mode shapes [batter, resonant], node-major: shape[node * modes + mode] is the
   * transverse displacement per unit modal coordinate [1/sqrt(kg)], positive downward.
   */
  modeShapes: [Float32Array, Float32Array];
  /** Wall time of the build or rebuild inside the worker [ms]. */
  buildMs: number;
  /** True when the heads were re-meshed and their eigenproblems re-solved. */
  eigenSolved: boolean;
  /** True for every rebuild: the stream restarted at rest and silent. */
  rebuilt: boolean;
  calibration?: DrumCalibration | null;
  engineVariant?: 'relaxed-simd' | 'baseline';
}

export interface DrumLedger {
  /** Energy in the instrument now (heads, air spring, contact, stick while present) [J]. */
  storedJ: number;
  /** Cumulative modal-damping + contact dissipation [J]. */
  dissipatedJ: number;
  /** Cumulative stick kinetic energy supplied by the player at launches [J]. */
  suppliedJ: number;
  /** Cumulative stick kinetic energy taken back by the player at catches [J]. */
  recoveredJ: number;
  /** stored - (supplied - recovered - dissipated): the uncorrected residual [J]. */
  balanceResidualJ: number;
  /** Largest single-step residual reported by the mechanics owner [J]. */
  maxStepResidualJ: number;
}

export interface DrumSnapshot {
  /** Engine audio time since the last (re)build [s]. This is ahead of what is audible by `aheadMs`. */
  timeS: number;
  mechanicalSteps: number;
  ledger: DrumLedger;
  stick: {
    /** The stick is currently a body of the mechanical system. */
    present: boolean;
    inContact: boolean;
    /** A stroke is waiting to launch. */
    pending: boolean;
    xM: number;
    yM: number;
    /** Tip displacement, positive downward, 0 = undeformed head plane [m]. */
    positionM: number;
    velocityMS: number;
    penetrationM: number;
    /** Elastic Hertz force at the accepted endpoint [N]. */
    contactForceN: number;
  };
  /** Sealed-air overpressure [Pa]. */
  cavityPressurePa: number;
  strikesLaunched: number;
  /** Peak |pressure| at [batter mic, resonant mic] since the previous snapshot [Pa]. */
  peakPa: [number, number];
  /** Frames that exceeded full scale since the last (re)build. */
  clippedFrames: number;
  batterModes: number;
  resonantModes: number;
  /** Modal displacements [m sqrt(kg)], batter then resonant. */
  q: Float64Array;
  /** Modal velocities [m sqrt(kg)/s], batter then resonant. */
  v: Float64Array;
  /** Transverse displacement of every mesh node [m], positive downward. */
  batterDisplacementM: Float32Array;
  resonantDisplacementM: Float32Array;
  /** Frames rendered since the last (re)build. */
  renderedFrames: number;
  /** How far this state leads the sink output (queued audio) [ms]; add device latency for the ear. */
  aheadMs: number;
}

export interface DrumStats {
  /** Audio seconds rendered per wall second over the last ~250 ms; null while idle. */
  realtimeFactor: number | null;
  worstChunkMs: number;
  queuedMs: number;
  targetMs: number;
  /** Times the worker raised the queue target (by 4 ms, up to 40 ms) after underruns. */
  targetRaises: number;
  underruns: number;
  underrunFrames: number;
  renderedFrames: number;
  playedFrames: number;
  substeps: number | null;
  mechanicalRateHz: number | null;
  monitor: DrumMonitor;
  engineVariant: 'relaxed-simd' | 'baseline' | null;
  halted: { code: string; message: string } | null;
  /** queue + device (base + output latency) + anti-alias filter [ms]. Excludes acoustic travel time. */
  latencyEstimateMs: number;
}

/** Which microphone signal goes to the speakers. `sum` is a unit-gain sum of both. */
export type DrumMonitor = 'batter' | 'resonant' | 'split' | 'sum';

export interface LiveDrumOptions {
  /** Initial parameters overlaid on the reference defaults. */
  params?: DrumParamPatch;
  /** Defaults to `/live`. */
  baseUrl?: string;
  /** Pick the highest mechanical rate sustaining `calibration.minFactor` at start (default true). */
  autoCalibrate?: boolean;
  calibration?: DrumCalibrationOptions;
  /** Audio kept queued ahead of the sink [ms] (default 12; lower = less latency, more underrun risk). */
  targetMs?: number;
  /** Raise the queue target 4 ms at a time (max 40 ms) when underruns occur (default true). */
  adaptiveTarget?: boolean;
  /** Snapshot rate [Hz] (default 30). */
  snapshotHz?: number;
  monitor?: DrumMonitor;
  /** 'baseline' forces the software-FMA build (for comparison). Default 'auto'. */
  engine?: 'auto' | 'baseline';
  /** Deliver snapshots when their audio is heard rather than when it is rendered (default false). */
  syncSnapshotsToAudio?: boolean;
  destination?: AudioNode;
}

type Pending = { resolve: (value: unknown) => void; reject: (reason: LiveRefusal) => void; settleOn: string };

const HEADER = 24;

function decodeSnapshot(
  state: Float64Array, displacements: Float32Array, renderedFrames: number, aheadMs: number,
): DrumSnapshot {
  const at = (index: number): number => state[index] ?? Number.NaN;
  const batterModes = at(22);
  const resonantModes = at(23);
  const modes = batterModes + resonantModes;
  const nodes = displacements.length / 2;
  return {
    timeS: at(0),
    mechanicalSteps: at(1),
    ledger: {
      storedJ: at(2),
      dissipatedJ: at(3),
      suppliedJ: at(4),
      recoveredJ: at(5),
      balanceResidualJ: at(6),
      maxStepResidualJ: at(7),
    },
    stick: {
      present: at(8) !== 0,
      inContact: at(9) !== 0,
      pending: at(10) !== 0,
      xM: at(11),
      yM: at(12),
      positionM: at(13),
      velocityMS: at(14),
      penetrationM: at(15),
      contactForceN: at(16),
    },
    cavityPressurePa: at(17),
    strikesLaunched: at(18),
    peakPa: [at(19), at(20)],
    clippedFrames: at(21),
    batterModes,
    resonantModes,
    q: state.subarray(HEADER, HEADER + modes),
    v: state.subarray(HEADER + modes, HEADER + 2 * modes),
    batterDisplacementM: displacements.subarray(0, nodes),
    resonantDisplacementM: displacements.subarray(nodes, 2 * nodes),
    renderedFrames,
    aheadMs,
  };
}

export class LiveDrumClient {
  /** Latest build description (updated by every rebuild). */
  ready: DrumReady;

  private readonly worker: Worker;
  private readonly sink: LiveSink;
  private readonly pending = new Map<number, Pending>();
  private readonly snapshotListeners = new Set<(snapshot: DrumSnapshot) => void>();
  private readonly statsListeners = new Set<(stats: DrumStats) => void>();
  private readonly readyListeners = new Set<(ready: DrumReady) => void>();
  private readonly errorListeners = new Set<(error: LiveRefusal) => void>();
  private readonly calibrationListeners = new Set<(calibration: DrumCalibration) => void>();
  private readonly syncSnapshots: boolean;
  private nextId = 1;
  private lastStats: DrumStats | null = null;
  private disposed = false;

  private constructor(worker: Worker, sink: LiveSink, ready: DrumReady, syncSnapshots: boolean) {
    this.worker = worker;
    this.sink = sink;
    this.ready = ready;
    this.syncSnapshots = syncSnapshots;
    worker.onmessage = (event: MessageEvent) => this.receive(event.data);
  }

  /**
   * Load the WebAssembly engine, build the drum from geometry, optionally
   * calibrate the mechanical rate to this device, and start the (silent) stream.
   * `audioContext` must run at 48 kHz (see `createLiveAudioContext`).
   */
  static async create(audioContext: AudioContext, options: LiveDrumOptions = {}): Promise<LiveDrumClient> {
    const baseUrl = options.baseUrl ?? LIVE_BASE_URL;
    const sink = await createLiveSink(audioContext, options.destination ? { baseUrl, destination: options.destination } : { baseUrl });
    const worker = new Worker(`${baseUrl}/drum/worker.js`, { type: 'module', name: 'fs-live-drum' });
    const port = sink.createProducerPort();
    const first = new Promise<DrumReady>((resolve, reject) => {
      worker.onmessage = (event: MessageEvent) => {
        const message = event.data;
        if (message?.type === 'ready') resolve(message as DrumReady);
        else if (message?.type === 'error') reject(new LiveRefusal(message.code, message.message, message.command));
      };
      worker.onerror = (event) => reject(new LiveRefusal('worker_error', event.message || 'the drum worker failed to load'));
    });
    worker.postMessage(
      {
        type: 'init',
        sampleRate: audioContext.sampleRate,
        port,
        params: options.params ?? {},
        autoCalibrate: options.autoCalibrate ?? true,
        calibration: options.calibration,
        targetMs: options.targetMs,
        adaptiveTarget: options.adaptiveTarget,
        snapshotHz: options.snapshotHz,
        monitor: options.monitor,
        engine: options.engine ?? 'auto',
      },
      [port],
    );
    let ready: DrumReady;
    try {
      ready = await first;
    } catch (error) {
      worker.terminate();
      sink.dispose();
      throw error;
    }
    if (audioContext.state === 'suspended') await audioContext.resume();
    return new LiveDrumClient(worker, sink, ready, options.syncSnapshotsToAudio ?? false);
  }

  /** Static description of the built instrument. */
  get info(): DrumInfo {
    return this.ready.info;
  }

  /**
   * Rebuild the instrument from geometry. Any subset of parameters; the rest
   * keep their current values. The old instrument stops and the new one starts
   * at rest and silent (nothing is cross-faded). Rejects with a `LiveRefusal`
   * and leaves the playing instrument untouched when the request is refused.
   */
  configure(params: DrumParamPatch): Promise<DrumReady> {
    return this.request<DrumReady>('configure', { params }, 'ready');
  }

  /**
   * Launch the stick at head-centred (x, y) [m] with a downward tip speed [m/s]
   * in (0, 20]. Resolves 'queued' (launches at the next mechanical step) or
   * 'deferred' (the previous stroke's tip is still in contact; this one launches
   * the moment it separates). Rejects with `strike_outside_head` / `strike_speed`.
   */
  strike(xM: number, yM: number, speedMS: number): Promise<'queued' | 'deferred'> {
    return this.request<{ outcome: 'queued' | 'deferred' }>('strike', { x: xM, y: yM, speed: speedMS }, 'strike').then(
      (reply) => reply.outcome,
    );
  }

  /** Change the mechanical rate (48 kHz x substeps). Restarts silent; the eigen-solve is reused. */
  setSubsteps(substeps: 1 | 2 | 4 | 8 | 16): Promise<DrumReady> {
    return this.request<DrumReady>('setSubsteps', { substeps }, 'ready');
  }

  /**
   * Measure this device and switch to the highest mechanical rate that sustains
   * `minFactor` (default 1.6x real time). Blocks audio while measuring.
   */
  calibrate(options: DrumCalibrationOptions = {}): Promise<DrumCalibration> {
    return this.request<DrumCalibration>('calibrate', options, 'calibration');
  }

  /** Measure every candidate rate without changing anything. Blocks audio while measuring. */
  bench(options: DrumCalibrationOptions = {}): Promise<DrumCalibration> {
    return this.request<DrumCalibration>('bench', { minFactor: 0, ...options }, 'calibration');
  }

  /** Same instrument, at rest and silent. Also clears a halted voice. */
  reset(): Promise<void> {
    return this.request<void>('reset', {}, 'done');
  }

  setMonitor(mode: DrumMonitor): Promise<void> {
    return this.request<void>('setMonitor', { mode }, 'done');
  }

  /** Audio kept queued ahead of the sink, 3..500 ms. */
  setTargetMs(ms: number): Promise<void> {
    return this.request<void>('setTargetMs', { ms }, 'done');
  }

  setSnapshotHz(hz: number): Promise<void> {
    return this.request<void>('setSnapshotHz', { hz }, 'done');
  }

  onSnapshot(listener: (snapshot: DrumSnapshot) => void): () => void {
    this.snapshotListeners.add(listener);
    return () => this.snapshotListeners.delete(listener);
  }

  onStats(listener: (stats: DrumStats) => void): () => void {
    this.statsListeners.add(listener);
    return () => this.statsListeners.delete(listener);
  }

  /** Fires after every successful rebuild (configure, setSubsteps, calibrate). */
  onReady(listener: (ready: DrumReady) => void): () => void {
    this.readyListeners.add(listener);
    return () => this.readyListeners.delete(listener);
  }

  /** Refusals that are not tied to a pending call (e.g. a halted voice). */
  onError(listener: (error: LiveRefusal) => void): () => void {
    this.errorListeners.add(listener);
    return () => this.errorListeners.delete(listener);
  }

  onCalibration(listener: (calibration: DrumCalibration) => void): () => void {
    this.calibrationListeners.add(listener);
    return () => this.calibrationListeners.delete(listener);
  }

  /**
   * Stroke-to-ear estimate [ms]: queued audio + device latency + the engine's
   * anti-alias filter. Acoustic travel from head to microphone (about 1 ms at
   * the default 0.35 m) is physics and is not counted as latency.
   */
  latencyEstimateMs(): number {
    const queued = this.lastStats?.queuedMs ?? 0;
    const filter = this.antiAliasMs();
    return queued + this.sink.deviceLatencyMs() + filter;
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    for (const waiter of this.pending.values()) waiter.reject(new LiveRefusal('disposed', 'the drum client was disposed'));
    this.pending.clear();
    this.worker.postMessage({ type: 'cmd', name: 'dispose' });
    // The worker closes itself after flushing; terminate as a backstop.
    setTimeout(() => this.worker.terminate(), 250);
    this.sink.flush();
    this.sink.dispose();
    this.snapshotListeners.clear();
    this.statsListeners.clear();
    this.readyListeners.clear();
    this.errorListeners.clear();
    this.calibrationListeners.clear();
  }

  private antiAliasMs(): number {
    const info = this.ready.info;
    return ((info.radiation.filter_delay_frames[0] ?? 0) * 1000) / info.sample_rate_hz;
  }

  private request<T>(name: string, args: object, settleOn: string): Promise<T> {
    if (this.disposed) return Promise.reject(new LiveRefusal('disposed', 'the drum client was disposed', name));
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject, settleOn });
      this.worker.postMessage({ type: 'cmd', id, name, args });
    });
  }

  private settle(id: unknown, type: string, value: unknown): void {
    if (typeof id !== 'number') return;
    const waiter = this.pending.get(id);
    if (waiter && waiter.settleOn === type) {
      this.pending.delete(id);
      waiter.resolve(value);
    }
  }

  private receive(message: { type?: string; id?: number } & Record<string, unknown>): void {
    switch (message.type) {
      case 'snapshot': {
        const snapshot = decodeSnapshot(
          message.state as Float64Array, message.displacements as Float32Array,
          message.renderedFrames as number, message.aheadMs as number,
        );
        if (this.syncSnapshots) {
          const delay = snapshot.aheadMs + this.sink.deviceLatencyMs();
          setTimeout(() => this.emitSnapshot(snapshot), delay);
        } else {
          this.emitSnapshot(snapshot);
        }
        break;
      }
      case 'stats': {
        const stats = message as unknown as Omit<DrumStats, 'latencyEstimateMs'>;
        const filter = this.antiAliasMs();
        this.lastStats = { ...stats, latencyEstimateMs: stats.queuedMs + this.sink.deviceLatencyMs() + filter };
        for (const listener of this.statsListeners) listener(this.lastStats);
        break;
      }
      case 'ready': {
        this.ready = message as unknown as DrumReady;
        for (const listener of this.readyListeners) listener(this.ready);
        this.settle(message.id, 'ready', this.ready);
        break;
      }
      case 'calibration': {
        const calibration = message as unknown as DrumCalibration;
        for (const listener of this.calibrationListeners) listener(calibration);
        this.settle(message.id, 'calibration', calibration);
        break;
      }
      case 'strike':
        this.settle(message.id, 'strike', message);
        break;
      case 'done':
        this.settle(message.id, 'done', undefined);
        break;
      case 'error': {
        const refusal = new LiveRefusal(
          String(message.code ?? 'worker_error'), String(message.message ?? 'unknown error'),
          typeof message.command === 'string' ? message.command : undefined,
        );
        const waiter = typeof message.id === 'number' ? this.pending.get(message.id) : undefined;
        if (waiter && typeof message.id === 'number') {
          this.pending.delete(message.id);
          waiter.reject(refusal);
        } else {
          for (const listener of this.errorListeners) listener(refusal);
        }
        break;
      }
      default:
        break;
    }
  }

  private emitSnapshot(snapshot: DrumSnapshot): void {
    if (this.disposed) return;
    for (const listener of this.snapshotListeners) listener(snapshot);
  }
}
