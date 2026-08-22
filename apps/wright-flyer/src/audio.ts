// Flight audio (Track 4, plan §5.6/E9.1): a small WebAudio graph driven
// ONLY by presentation-plane state (prop speed, airspeed, phase) — it
// never touches the sim, the digests, or the replay identity.
//
// HONESTY: this is a synthesized soundscape, not a historical recording.
// The engine tone is a two-oscillator approximation of the 1903 inline
// four's firing cadence; the mix mapping below is PURE and tested
// headless (test/audio.test.ts). AudioContext starts only after a user
// gesture (browser policy); M toggles mute.

/** Engine firing frequency [Hz] from prop omega: engine = prop·23/8,
 * and a 4-stroke inline-4 fires twice per revolution. */
export function engineFreqHz(propOmegaRadS: number): number {
  const engineRps = (propOmegaRadS * (23 / 8)) / (2 * Math.PI);
  return Math.max(8, engineRps * 2);
}

/** Normalized prop speed in [0,1] against the 1025 rpm trim. */
export function rpm01FromOmega(propOmegaRadS: number): number {
  const trim = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  return Math.min(1, Math.max(0, propOmegaRadS / trim));
}

export interface MixLevels {
  /** Engine tone gain [0, 0.5]. */
  readonly engine: number;
  /** Wind bed gain [0, 0.4]. */
  readonly wind: number;
  /** Rail rumble gain [0, 0.3] — only while rolling the rail. */
  readonly rumble: number;
}

/** PURE mix law: gains from the live presentation state. Deterministic
 * in its inputs; clamped so a runaway snapshot can never blast the
 * listener. */
export function mixLevels(
  rpm01: number,
  airspeedMps: number,
  onRail: boolean,
  groundSpeedMps: number,
): MixLevels {
  const r = Math.min(1, Math.max(0, rpm01));
  const v = Math.max(0, airspeedMps);
  const g = Math.max(0, groundSpeedMps);
  return {
    engine: 0.5 * Math.min(1, r * 1.25),
    wind: Math.min(0.4, Math.max(0, (v - 4) / 26) * 0.4),
    rumble: onRail ? Math.min(0.3, (g / 12) * 0.3) : 0,
  };
}

interface AudioCtor {
  // Structural typing keeps the class testable without a real context.
  readonly currentTime: number;
  readonly destination: AudioDestinationNode;
  createOscillator(): OscillatorNode;
  createGain(): GainNode;
  createBiquadFilter(): BiquadFilterNode;
  createBuffer(channels: number, frames: number, rate: number): AudioBuffer;
  createBufferSource(): AudioBufferSourceNode;
}

/** The flight soundscape. All node churn is lazy; `update` only sets
 * AudioParams (no allocations on the frame path). */
/** Master bus level. The original 0.85 was near full scale — with a
 * raw-saw engine it read as a loud buzzer. ~-10 dB is a soundscape. */
const MASTER_LEVEL = 0.3;

export class FlightAudio {
  private ctx: AudioContext | null = null;
  private engineGain: GainNode | null = null;
  private engineOsc: OscillatorNode | null = null;
  private engineSub: OscillatorNode | null = null;
  private windGain: GainNode | null = null;
  private rumbleGain: GainNode | null = null;
  private surfGain: GainNode | null = null;
  private master: GainNode | null = null;
  private muted = false;
  private lastCryS = 0;
  private readonly withOcean: boolean;

  constructor(opts: { withOcean: boolean }) {
    this.withOcean = opts.withOcean;
  }

  get isMuted(): boolean {
    return this.muted;
  }

  /** Must be called from a user-gesture handler. Idempotent. */
  ensureStarted(): void {
    if (this.ctx !== null) {
      if (this.ctx.state === "suspended") {
        void this.ctx.resume();
      }
      return;
    }
    const Ctor = window.AudioContext as unknown as new () => AudioContext | undefined;
    const ctx = new Ctor();
    if (ctx === undefined) {
      return;
    }
    this.ctx = ctx;
    const master = ctx.createGain();
    master.gain.value = MASTER_LEVEL;
    master.connect(ctx.destination);
    this.master = master;
    // Engine: saw at the firing frequency + sub square an octave down,
    // through a warm lowpass.
    const eGain = ctx.createGain();
    eGain.gain.value = 0;
    const lp = ctx.createBiquadFilter();
    lp.type = "lowpass";
    lp.frequency.value = 520; // darker: chug, not buzz
    eGain.connect(lp);
    lp.connect(master);
    const osc = ctx.createOscillator();
    osc.type = "triangle"; // raw saw at speaker level was a buzzer
    osc.frequency.value = 40;
    osc.connect(eGain);
    osc.start();
    const sub = ctx.createOscillator();
    sub.type = "triangle";
    sub.frequency.value = 20;
    const subG = ctx.createGain();
    subG.gain.value = 0.2;
    sub.connect(subG);
    subG.connect(eGain);
    sub.start();
    this.engineGain = eGain;
    this.engineOsc = osc;
    this.engineSub = sub;
    // Shared noise buffer: 2 s of white noise.
    const noise = ctx.createBuffer(1, ctx.sampleRate * 2, ctx.sampleRate);
    const data = noise.getChannelData(0);
    let brown = 0;
    let peak = 1e-6;
    for (let i = 0; i < data.length; i += 1) {
      const white = Math.random() * 2 - 1;
      brown = (brown + 0.02 * white) / 1.02;
      data[i] = brown;
      peak = Math.max(peak, Math.abs(brown));
    }
    // Normalize to a safe peak — the previous fixed x3.2 scale could
    // hard-clip the walk's excursions (audible crunch in the wind bed).
    for (let i = 0; i < data.length; i += 1) {
      data[i] = (data[i]! / peak) * 0.65;
    }
    const mkNoise = (freq: number, type: BiquadFilterType, gain: number): GainNode => {
      const g = ctx.createGain();
      g.gain.value = gain;
      const f = ctx.createBiquadFilter();
      f.type = type;
      f.frequency.value = freq;
      const src = ctx.createBufferSource();
      src.buffer = noise;
      src.loop = true;
      src.connect(f);
      f.connect(g);
      g.connect(master);
      src.start();
      return g;
    };
    this.windGain = mkNoise(650, "bandpass", 0);
    this.rumbleGain = mkNoise(110, "lowpass", 0);
    if (this.withOcean) {
      this.surfGain = mkNoise(280, "lowpass", 0.05);
    }
  }

  /** Per-frame mix update (cheap: AudioParams only). */
  update(state: {
    propOmegaRadS: number;
    airspeedMps: number;
    onRail: boolean;
    groundSpeedMps: number;
    nowS: number;
  }): void {
    if (this.ctx === null || this.engineGain === null || this.engineOsc === null) {
      return;
    }
    const mix = mixLevels(
      rpm01FromOmega(state.propOmegaRadS),
      state.airspeedMps,
      state.onRail,
      state.groundSpeedMps,
    );
    const t = this.ctx.currentTime;
    this.engineGain.gain.setTargetAtTime(this.muted ? 0 : mix.engine * 0.5, t, 0.08);
    this.engineOsc.frequency.setTargetAtTime(engineFreqHz(state.propOmegaRadS), t, 0.06);
    if (this.engineSub !== null) {
      this.engineSub.frequency.setTargetAtTime(engineFreqHz(state.propOmegaRadS) / 2, t, 0.06);
    }
    if (this.windGain !== null) {
      this.windGain.gain.setTargetAtTime(this.muted ? 0 : mix.wind * 0.7, t, 0.15);
    }
    if (this.rumbleGain !== null) {
      this.rumbleGain.gain.setTargetAtTime(this.muted ? 0 : mix.rumble * 0.8, t, 0.1);
    }
    // Occasional gull cry: a two-note descending chirp, at most every
    // ~9 s (presentation flavor — no determinism claim).
    if (!this.muted && state.nowS - this.lastCryS > 28 && Math.random() < 0.002) {
      this.lastCryS = state.nowS;
      this.gullCry();
    }
  }

  /** Daniels' plate: a sharp mechanical click (called on the flash). */
  shutter(): void {
    if (this.ctx === null || this.master === null || this.muted) {
      return;
    }
    const ctx = this.ctx;
    const g = ctx.createGain();
    g.gain.setValueAtTime(0.22, ctx.currentTime);
    g.gain.exponentialRampToValueAtTime(0.001, ctx.currentTime + 0.07);
    const f = ctx.createBiquadFilter();
    f.type = "highpass";
    f.frequency.value = 1800;
    const src = ctx.createBufferSource();
    const buf = ctx.createBuffer(1, ctx.sampleRate * 0.07, ctx.sampleRate);
    const d = buf.getChannelData(0);
    for (let i = 0; i < d.length; i += 1) {
      d[i] = (Math.random() * 2 - 1) * (1 - i / d.length);
    }
    src.buffer = buf;
    src.connect(f);
    f.connect(g);
    g.connect(this.master);
    src.start();
  }

  toggleMute(): boolean {
    this.muted = !this.muted;
    if (this.master !== null && this.ctx !== null) {
      this.master.gain.setTargetAtTime(this.muted ? 0 : MASTER_LEVEL, this.ctx.currentTime, 0.05);
    }
    return this.muted;
  }

  private gullCry(): void {
    if (this.ctx === null || this.master === null) {
      return;
    }
    const ctx = this.ctx;
    const g = ctx.createGain();
    g.gain.value = 0.03;
    g.connect(this.master);
    const osc = ctx.createOscillator();
    osc.type = "sine";
    const t0 = ctx.currentTime;
    osc.frequency.setValueAtTime(1250, t0);
    osc.frequency.exponentialRampToValueAtTime(820, t0 + 0.16);
    osc.frequency.setValueAtTime(1150, t0 + 0.2);
    osc.frequency.exponentialRampToValueAtTime(700, t0 + 0.38);
    g.gain.setValueAtTime(0.03, t0);
    g.gain.setValueAtTime(0, t0 + 0.42);
    osc.connect(g);
    osc.start(t0);
    osc.stop(t0 + 0.45);
  }
}
