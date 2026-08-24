// Flight audio (Track 4, plan §5.6/E9.1): a small WebAudio graph driven
// ONLY by presentation-plane state (prop speed, airspeed, phase) — it
// never touches the sim, the digests, or the replay identity.
//
// HONESTY: this is a synthesized soundscape, not a historical recording.
// The engine tone is a two-oscillator approximation of the 1903 inline
// four's firing cadence; the mix mapping below is PURE and tested
// headless (test/audio.test.ts). AudioContext starts only after a user
// gesture (browser policy); M toggles mute.

/** Labeled sound design disclaimer constant (E9.1). */
export const AUDIO_DISCLAIMER_LABEL =
  "Synthesized sound design (WebAudio/BPF/4-cyl model), not an acoustic simulation claim";

/** Propeller blade passage frequency [Hz]: BPF = blade_count · RPM / 60 = blade_count · (omega / 2π).
 * For the 1903 Flyer with twin 2-blade props, bladeCount = 2.
 */
export function propBpfHz(propOmegaRadS: number, bladeCount = 2): number {
  if (propOmegaRadS <= 0 || bladeCount <= 0) {
    return 0;
  }
  const rps = propOmegaRadS / (2 * Math.PI);
  return rps * bladeCount;
}

/** Twin-propeller partially-coherent BPF frequencies [Hz].
 * The Flyer had two counter-rotating 2-blade propellers driven by crossed/uncrossed
 * chains at a 23:8 engine-to-prop reduction. Slight blade pitch/inflow discrepancies
 * and chain elasticity produce a slight beat frequency (~0.3-0.5 Hz) between the two
 * props, creating characteristic twin-rotor beating.
 */
export function twinPropBpfHz(
  propOmegaRadS: number,
  bladeCount = 2,
  inflowAsymmetry01 = 0.038,
): { readonly bpfLeftHz: number; readonly bpfRightHz: number; readonly beatHz: number } {
  const nominal = propBpfHz(propOmegaRadS, bladeCount);
  if (nominal <= 0) {
    return { bpfLeftHz: 0, bpfRightHz: 0, beatHz: 0 };
  }
  const delta = nominal * Math.min(0.1, Math.max(0, inflowAsymmetry01));
  const left = Math.max(0, nominal - delta * 0.5);
  const right = Math.max(0, nominal + delta * 0.5);
  return {
    bpfLeftHz: left,
    bpfRightHz: right,
    beatHz: right - left,
  };
}

/** Airframe turbulence / wire noise gain [0, 0.35].
 * Scales with dynamic pressure (v^2) above walking speed, plus extra separation
 * turbulence when angle of attack indicates high-lift / near-stall flow.
 */
export function airframeNoiseGain(airspeedMps: number, alphaDeg = 0): number {
  if (airspeedMps <= 0) {
    return 0;
  }
  const u = Math.max(0, airspeedMps - 3.5);
  const base = Math.min(0.25, (u / 24) * (u / 24) * 0.25);
  const separationBoost = Math.max(0, Math.min(0.1, (Math.abs(alphaDeg) - 10) * 0.02));
  return Math.min(0.35, base + separationBoost);
}

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

/** PURE wind-rush law (B12): pressure drag on the ear scales with the
 * SQUARE of airspeed above a walking headwind; clamped well below the
 * engine so it textures rather than dominates. */
export function windRushGain(airspeedMps: number): number {
  const u = Math.max(0, airspeedMps - 4);
  return Math.min(0.32, (u / 26) * (u / 26));
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
  private surfWanted = true;
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
    let ctx: AudioContext;
    try {
      ctx = new AudioContext();
    } catch {
      return; // no WebAudio: the game runs silent, never crashes
    }
    this.ctx = ctx;
    const master = ctx.createGain();
    master.gain.value = MASTER_LEVEL;
    master.connect(ctx.destination);
    this.master = master;
    // Engine: warm triangle at the firing frequency + sub an octave
    // down, through a dark lowpass — a chug, not a buzzer.
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
    if (this.withOcean && this.surfWanted) {
      this.surfGain = mkNoise(280, "lowpass", 0.05);
    }
    // Wire whistling: high bandpass for rigging whistling in the wind.
    this.wireGain = mkNoise(1600, "bandpass", 0);

    // Twin counter-rotating propeller blade passage nodes
    const pGain = ctx.createGain();
    pGain.gain.value = 0;
    const pFilter = ctx.createBiquadFilter();
    pFilter.type = "lowpass";
    pFilter.frequency.value = 180;
    pGain.connect(pFilter);
    pFilter.connect(master);
    
    const propL = ctx.createOscillator();
    propL.type = "sine";
    propL.frequency.value = 12;
    propL.connect(pGain);
    propL.start();

    const propR = ctx.createOscillator();
    propR.type = "sine";
    propR.frequency.value = 12.4;
    propR.connect(pGain);
    propR.start();

    this.propBeatingGain = pGain;
    this.propBeatingOscL = propL;
    this.propBeatingOscR = propR;
  }

  private wireGain: GainNode | null = null;
  private propBeatingOscL: OscillatorNode | null = null;
  private propBeatingOscR: OscillatorNode | null = null;
  private propBeatingGain: GainNode | null = null;

  /** Per-frame mix update (cheap: AudioParams only). Optional fields:
   * gust01 drives structural creaks; surfFacing01 (camera dot toward
   * the Atlantic, [0,1]) crossfades the surf bed by heading. */
  update(state: {
    propOmegaRadS: number;
    airspeedMps: number;
    onRail: boolean;
    groundSpeedMps: number;
    nowS: number;
    gust01?: number;
    surfFacing01?: number;
    alphaDeg?: number;
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
    const engFreq = engineFreqHz(state.propOmegaRadS);
    const twinBpf = twinPropBpfHz(state.propOmegaRadS, 2);

    this.engineGain.gain.setTargetAtTime(this.muted ? 0 : mix.engine * 0.52, t, 0.08);
    this.engineOsc.frequency.setTargetAtTime(engFreq, t, 0.06);
    if (this.engineSub !== null) {
      this.engineSub.frequency.setTargetAtTime(engFreq / 2, t, 0.06);
    }
    if (this.propBeatingGain !== null && this.propBeatingOscL !== null && this.propBeatingOscR !== null) {
      this.propBeatingGain.gain.setTargetAtTime(this.muted ? 0 : mix.engine * 0.28, t, 0.08);
      this.propBeatingOscL.frequency.setTargetAtTime(Math.max(4, twinBpf.bpfLeftHz), t, 0.06);
      this.propBeatingOscR.frequency.setTargetAtTime(Math.max(4, twinBpf.bpfRightHz), t, 0.06);
    }
    if (this.windGain !== null) {
      // Wind bed + square-law rush layered under one clamp.
      const rush = windRushGain(state.airspeedMps);
      this.windGain.gain.setTargetAtTime(
        this.muted ? 0 : Math.min(0.6, mix.wind * 0.75 + rush),
        t,
        0.15,
      );
    }
    if (this.wireGain !== null) {
      const wireLevel = airframeNoiseGain(state.airspeedMps, state.alphaDeg);
      this.wireGain.gain.setTargetAtTime(this.muted ? 0 : wireLevel, t, 0.2);
    }
    if (this.rumbleGain !== null) {
      this.rumbleGain.gain.setTargetAtTime(this.muted ? 0 : mix.rumble * 0.85, t, 0.1);
    }
    if (this.surfGain !== null && this.ctx !== null) {
      // Heading crossfade: the surf swells as the lens faces east.
      const facing = Math.max(0, Math.min(1, state.surfFacing01 ?? 0.5));
      this.surfGain.gain.setTargetAtTime(
        this.surfWanted && !this.muted ? 0.02 + 0.06 * facing : 0,
        t,
        0.5,
      );
    }
    // Structural creaks: the airframe talks when gusts load the rig.
    if (!this.muted && (state.gust01 ?? 0) > 0.5 && Math.random() < (state.gust01 ?? 0) * 0.01) {
      this.creak();
    }
    // Occasional realistic gull cries: randomized between short chirp,
    // melodic double-cry, and long soaring screech.
    if (!this.muted && state.nowS - this.lastCryS > 16 && Math.random() < 0.005) {
      this.lastCryS = state.nowS;
      this.gullCry();
    }
  }

  /** Wire release twang at takeoff launch start. */
  twang(): void {
    if (this.ctx === null || this.master === null || this.muted) {
      return;
    }
    const ctx = this.ctx;
    const osc = ctx.createOscillator();
    const g = ctx.createGain();
    osc.type = "sine";
    const t = ctx.currentTime;
    osc.frequency.setValueAtTime(440, t);
    osc.frequency.exponentialRampToValueAtTime(110, t + 0.18);
    g.gain.setValueAtTime(0.25, t);
    g.gain.exponentialRampToValueAtTime(0.001, t + 0.2);
    osc.connect(g);
    g.connect(this.master);
    osc.start(t);
    osc.stop(t + 0.22);
  }

  /** Landing skid sand crunch foley on touchdown. */
  sandCrunch(): void {
    if (this.ctx === null || this.master === null || this.muted) {
      return;
    }
    const ctx = this.ctx;
    const g = ctx.createGain();
    g.gain.setValueAtTime(0.3, ctx.currentTime);
    g.gain.exponentialRampToValueAtTime(0.001, ctx.currentTime + 0.6);
    const f = ctx.createBiquadFilter();
    f.type = "bandpass";
    f.frequency.value = 850;
    const buf = ctx.createBuffer(1, Math.floor(ctx.sampleRate * 0.6), ctx.sampleRate);
    const d = buf.getChannelData(0);
    for (let i = 0; i < d.length; i += 1) {
      d[i] = (Math.random() * 2 - 1) * Math.pow(1 - i / d.length, 1.5);
    }
    const src = ctx.createBufferSource();
    src.buffer = buf;
    src.connect(f);
    f.connect(g);
    g.connect(this.master);
    src.start();
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
    const buf = ctx.createBuffer(1, Math.floor(ctx.sampleRate * 0.07), ctx.sampleRate);
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

  /** Airframe creak (B12): a dry wooden knock — two detuned triangles
   * through a woody bandpass — for gust-loaded structure. */
  private creak(): void {
    if (this.ctx === null || this.master === null || this.muted) {
      return;
    }
    const ctx = this.ctx;
    const t0 = ctx.currentTime;
    const g = ctx.createGain();
    g.gain.setValueAtTime(0.05 + Math.random() * 0.03, t0);
    g.gain.exponentialRampToValueAtTime(0.001, t0 + 0.32);
    const bp = ctx.createBiquadFilter();
    bp.type = "bandpass";
    bp.frequency.value = 170 + Math.random() * 90;
    bp.Q.value = 2.2;
    g.connect(bp);
    bp.connect(this.master);
    for (const det of [1, 1.047]) {
      const o = ctx.createOscillator();
      o.type = "triangle";
      o.frequency.setValueAtTime(140 * det, t0);
      o.frequency.exponentialRampToValueAtTime(95 * det, t0 + 0.28);
      o.connect(g);
      o.start(t0);
      o.stop(t0 + 0.34);
    }
  }

  /** Site surf gate: Huffman Prairie is landlocked Ohio pasture —
   * the ocean bed plays at Kill Devil Hills only. Call once at scene
   * build; safe before or after ensureStarted. */
  setSurfEnabled(on: boolean): void {
    this.surfWanted = on;
    if (this.surfGain !== null && this.ctx !== null) {
      this.surfGain.gain.setTargetAtTime(
        on ? 0.05 : 0,
        this.ctx.currentTime,
        0.4,
      );
    }
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
    // Spatialize (B12): random azimuth pan + distance attenuation so
    // the flock reads as AROUND the listener, not inside the head.
    const pan = ctx.createStereoPanner();
    pan.pan.value = Math.random() * 1.6 - 0.8;
    g.gain.value = 0.035 * (0.45 + Math.random() * 0.55);
    g.connect(pan);
    pan.connect(this.master);
    const osc = ctx.createOscillator();
    osc.type = "sine";
    const t0 = ctx.currentTime;
    const kind = Math.floor(Math.random() * 3);
    if (kind === 0) {
      // Classic 2-note descending gull call
      osc.frequency.setValueAtTime(1350, t0);
      osc.frequency.exponentialRampToValueAtTime(840, t0 + 0.18);
      osc.frequency.setValueAtTime(1200, t0 + 0.22);
      osc.frequency.exponentialRampToValueAtTime(720, t0 + 0.44);
      g.gain.setValueAtTime(g.gain.value, t0);
      g.gain.setValueAtTime(0, t0 + 0.48);
      osc.connect(g);
      osc.start(t0);
      osc.stop(t0 + 0.5);
    } else if (kind === 1) {
      // High rising-then-falling screech
      osc.frequency.setValueAtTime(900, t0);
      osc.frequency.exponentialRampToValueAtTime(1600, t0 + 0.12);
      osc.frequency.exponentialRampToValueAtTime(750, t0 + 0.35);
      g.gain.setValueAtTime(g.gain.value, t0);
      g.gain.setValueAtTime(0, t0 + 0.38);
      osc.connect(g);
      osc.start(t0);
      osc.stop(t0 + 0.4);
    } else {
      // Gentle warble
      osc.frequency.setValueAtTime(1100, t0);
      osc.frequency.exponentialRampToValueAtTime(950, t0 + 0.15);
      osc.frequency.exponentialRampToValueAtTime(1050, t0 + 0.25);
      osc.frequency.exponentialRampToValueAtTime(680, t0 + 0.5);
      g.gain.setValueAtTime(g.gain.value, t0);
      g.gain.setValueAtTime(0, t0 + 0.52);
      osc.connect(g);
      osc.start(t0);
      osc.stop(t0 + 0.55);
    }
  }
}
