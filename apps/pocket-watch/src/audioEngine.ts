// Real-time 48 kHz virtual microphone: the wasm engine (`PocketWatchSim::render_audio_block`)
// synthesizes the structure-borne + airborne pressure at the probe; this module plays it through an
// AudioWorklet ring buffer (ScriptProcessor fallback) and draws a triggered oscilloscope.

import type { PocketWatchSim } from "./wasm-pkg/fs_pocket_watch";

export const AUDIO_SAMPLE_RATE = 48_000;
const SCOPE_RING = 16_384; // power of two, ~341 ms: always holds at least one full beat (200 ms)
const SCOPE_WINDOW = 1_200; // 25 ms: unlock, impulse and drop of one beat
const TARGET_BUFFERED = 3_072; // ~64 ms of queued audio in the worklet
const MAX_CHUNK = 4_096;

export type MicOutputMode = "off" | "worklet" | "script";

export class VirtualMicrophoneAudio {
  private ctx: AudioContext | null = null;
  private node: AudioWorkletNode | null = null;
  private processor: ScriptProcessorNode | null = null;
  private outputMode: MicOutputMode = "off";
  private active = false;
  private workletFailed = false;
  private readonly scopeRing = new Float32Array(SCOPE_RING);
  private scopeWrite = 0;
  private scopeFill = 0;
  private triggerAt = -1;
  private scopePeak = 0.05;
  // Producer-side estimate of how much audio the worklet still has queued.
  private reportedAvailable = 0;
  private reportedAtMs = 0;
  private producedSinceReport = 0;
  private underruns = 0;
  private scopeCtx: CanvasRenderingContext2D | null = null;
  private scopeCssW = 0;
  private scopeCssH = 0;
  private scopeDpr = 0;
  public lastRmsSplDb = 32.0;

  public isRunning(): boolean {
    return this.active;
  }

  public mode(): MicOutputMode {
    return this.active ? this.outputMode : "off";
  }

  public underrunCount(): number {
    return this.underruns;
  }

  public async toggle(sim: PocketWatchSim): Promise<boolean> {
    if (this.active) {
      this.stop();
      return false;
    }
    await this.start(sim);
    return this.active;
  }

  public async start(sim: PocketWatchSim): Promise<void> {
    if (!this.ctx) {
      const Ctor =
        window.AudioContext ??
        (window as unknown as { webkitAudioContext?: typeof AudioContext }).webkitAudioContext;
      if (!Ctor) return;
      this.ctx = new Ctor({ sampleRate: AUDIO_SAMPLE_RATE, latencyHint: "interactive" });
    }
    if (this.ctx.state === "suspended") {
      await this.ctx.resume();
    }
    if (!this.node && !this.processor) {
      if (!this.workletFailed && this.ctx.audioWorklet) {
        try {
          await this.ctx.audioWorklet.addModule(`${import.meta.env.BASE_URL}mic-worklet.js`);
          this.node = new AudioWorkletNode(this.ctx, "fs-mic-ring", {
            numberOfInputs: 0,
            numberOfOutputs: 1,
            outputChannelCount: [2],
          });
          this.node.port.onmessage = (ev: MessageEvent) => {
            const msg = ev.data as { type?: string; available?: number; underruns?: number };
            if (msg.type === "level") {
              this.reportedAvailable = msg.available ?? 0;
              this.reportedAtMs = performance.now();
              this.producedSinceReport = 0;
              this.underruns = msg.underruns ?? this.underruns;
            }
          };
          this.node.connect(this.ctx.destination);
          this.outputMode = "worklet";
        } catch {
          this.workletFailed = true;
          this.node = null;
        }
      }
      if (!this.node) {
        this.installScriptProcessor(sim);
      }
    }
    this.flushQueue();
    this.active = true;
  }

  /** Fallback for browsers without AudioWorklet: pull blocks on the audio callback. */
  private installScriptProcessor(sim: PocketWatchSim): void {
    if (!this.ctx) return;
    const proc = this.ctx.createScriptProcessor(1024, 0, 2);
    proc.onaudioprocess = (ev: AudioProcessingEvent) => {
      const left = ev.outputBuffer.getChannelData(0);
      const right = ev.outputBuffer.getChannelData(1);
      if (!this.active) {
        left.fill(0);
        right.fill(0);
        return;
      }
      const samples = sim.render_audio_block(left.length);
      left.set(samples.subarray(0, left.length));
      right.set(samples.subarray(0, left.length));
      this.feedScope(samples);
    };
    proc.connect(this.ctx.destination);
    this.processor = proc;
    this.outputMode = "script";
  }

  public stop(): void {
    this.active = false;
    this.flushQueue();
  }

  /** Page hidden: suspend the audio clock so nothing plays stale or underruns in the background. */
  public async suspend(): Promise<void> {
    if (this.ctx && this.ctx.state === "running") {
      try {
        await this.ctx.suspend();
      } catch {
        /* ignore */
      }
    }
  }

  public async resume(): Promise<void> {
    if (this.ctx && this.active && this.ctx.state === "suspended") {
      try {
        await this.ctx.resume();
      } catch {
        /* ignore */
      }
    }
    this.flushQueue();
  }

  private flushQueue(): void {
    this.node?.port.postMessage({ type: "flush" });
    this.reportedAvailable = 0;
    this.reportedAtMs = performance.now();
    this.producedSinceReport = 0;
  }

  /**
   * Call once per animation frame while the simulation runs. Generates audio at the real-time rate:
   * live (worklet) it tops the ring buffer up to the latency target; muted it renders exactly
   * `dtRealS * 48 kHz` samples so the oscilloscope and the engine's loudness meters stay live.
   * The ScriptProcessor fallback pulls on the audio thread instead, so this is a no-op there.
   */
  public pump(sim: PocketWatchSim, dtRealS: number): void {
    if (this.active && this.outputMode === "script") return;
    let frames: number;
    if (this.active && this.node) {
      const elapsed = ((performance.now() - this.reportedAtMs) / 1000) * AUDIO_SAMPLE_RATE;
      const queued = Math.max(0, this.reportedAvailable + this.producedSinceReport - elapsed);
      const lead = dtRealS * AUDIO_SAMPLE_RATE; // cover until the next frame arrives
      frames = Math.round(TARGET_BUFFERED + lead - queued);
    } else {
      frames = Math.round(dtRealS * AUDIO_SAMPLE_RATE);
    }
    if (frames < 64) return;
    frames = Math.min(MAX_CHUNK, frames);
    const samples = sim.render_audio_block(frames);
    this.feedScope(samples);
    if (this.active && this.node) {
      // Transfer a copy: the wasm-bindgen array is already a fresh JS copy, so move it.
      this.node.port.postMessage({ type: "chunk", samples }, [samples.buffer]);
      this.producedSinceReport += samples.length;
    }
  }

  private feedScope(samples: Float32Array): void {
    let sumSq = 0;
    const n = samples.length;
    for (let i = 0; i < n; i++) {
      const s = samples[i] ?? 0;
      sumSq += s * s;
      this.scopeRing[this.scopeWrite] = s;
      this.scopeWrite = (this.scopeWrite + 1) & (SCOPE_RING - 1);
    }
    this.scopeFill = Math.min(SCOPE_RING, this.scopeFill + n);
    if (n > 0) {
      const rms = Math.sqrt(sumSq / n);
      const pa = rms * 0.25;
      const db = 20 * Math.log10(Math.max(pa, 2.0e-5) / 2.0e-5);
      this.lastRmsSplDb = 0.8 * this.lastRmsSplDb + 0.2 * Math.min(110, Math.max(26, db));
    }
  }

  private ensureScopeCanvas(canvas: HTMLCanvasElement): CanvasRenderingContext2D | null {
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const cssW = canvas.clientWidth;
    const cssH = canvas.clientHeight;
    if (cssW < 2 || cssH < 2) return null;
    if (!this.scopeCtx || cssW !== this.scopeCssW || cssH !== this.scopeCssH || dpr !== this.scopeDpr) {
      canvas.width = Math.round(cssW * dpr);
      canvas.height = Math.round(cssH * dpr);
      this.scopeCssW = cssW;
      this.scopeCssH = cssH;
      this.scopeDpr = dpr;
      this.scopeCtx = canvas.getContext("2d");
      this.scopeCtx?.setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    return this.scopeCtx;
  }

  /** Find the newest beat onset that has a full window of samples after it (ring-relative index). */
  private findTrigger(): number {
    const N = SCOPE_RING;
    const span = Math.min(this.scopeFill, 12_000);
    if (span < SCOPE_WINDOW * 2) return -1;
    let peak = 0;
    for (let k = SCOPE_WINDOW; k < span; k++) {
      const v = Math.abs(this.scopeRing[(this.scopeWrite - k + N) & (N - 1)] ?? 0);
      if (v > peak) peak = v;
    }
    if (peak < 1e-4) return -1;
    this.scopePeak = Math.max(peak, this.scopePeak * 0.9);
    const thr = peak * 0.45;
    // Walk backwards from the newest sample that still has a full window after it; the onset of the
    // newest burst is where the signal first rises above threshold after >= 4 ms of quiet.
    let onset = -1;
    let quiet = 0;
    for (let k = SCOPE_WINDOW; k < span; k++) {
      const v = Math.abs(this.scopeRing[(this.scopeWrite - k + N) & (N - 1)] ?? 0);
      if (v >= thr) {
        onset = k;
        quiet = 0;
      } else if (onset >= 0 && ++quiet > 192) {
        break;
      }
    }
    return onset < 0 ? -1 : (this.scopeWrite - onset + N) & (N - 1);
  }

  public renderOscilloscope(canvas: HTMLCanvasElement, micPosMm: [number, number, number]): void {
    const ctx = this.ensureScopeCanvas(canvas);
    if (!ctx) return;
    const W = this.scopeCssW;
    const H = this.scopeCssH;
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#070a10";
    ctx.fillRect(0, 0, W, H);

    // Graticule: 10 divisions horizontally (2.5 ms each), 4 vertically.
    ctx.strokeStyle = "rgba(231, 196, 122, 0.09)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 1; i < 10; i++) {
      const x = Math.round((i / 10) * W) + 0.5;
      ctx.moveTo(x, 0);
      ctx.lineTo(x, H);
    }
    for (let j = 1; j < 4; j++) {
      const y = Math.round((j / 4) * H) + 0.5;
      ctx.moveTo(0, y);
      ctx.lineTo(W, y);
    }
    ctx.stroke();

    const N = SCOPE_RING;
    const trig = this.findTrigger();
    if (trig >= 0) this.triggerAt = trig;
    const pre = Math.round(SCOPE_WINDOW * 0.12);
    const start =
      this.triggerAt >= 0
        ? (this.triggerAt - pre + N) & (N - 1)
        : (this.scopeWrite - SCOPE_WINDOW + N) & (N - 1);
    const gain = (H * 0.42) / Math.max(1e-4, this.scopePeak);
    const mid = H / 2;

    // Min/max envelope per pixel column, like a sampling oscilloscope.
    ctx.strokeStyle = this.active ? "#6ff0c0" : "#e8c77d";
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    const cols = Math.max(2, Math.floor(W));
    for (let x = 0; x < cols; x++) {
      const a = Math.floor((x / cols) * SCOPE_WINDOW);
      const b = Math.max(a + 1, Math.floor(((x + 1) / cols) * SCOPE_WINDOW));
      let lo = Infinity;
      let hi = -Infinity;
      for (let k = a; k < b; k++) {
        const v = this.scopeRing[(start + k) & (N - 1)] ?? 0;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
      const y0 = Math.max(1, Math.min(H - 1, mid - hi * gain));
      const y1 = Math.max(1, Math.min(H - 1, mid - lo * gain));
      if (x === 0) ctx.moveTo(x + 0.5, y0);
      ctx.lineTo(x + 0.5, y0);
      ctx.lineTo(x + 0.5, Math.max(y1, y0 + 0.75));
    }
    ctx.stroke();

    // Trigger marker
    if (this.triggerAt >= 0) {
      const tx = (pre / SCOPE_WINDOW) * W;
      ctx.strokeStyle = "rgba(231, 196, 122, 0.35)";
      ctx.setLineDash([2, 3]);
      ctx.beginPath();
      ctx.moveTo(tx + 0.5, 0);
      ctx.lineTo(tx + 0.5, H);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    ctx.fillStyle = "rgba(226, 222, 212, 0.62)";
    ctx.font = '500 10px "JetBrains Mono", ui-monospace, monospace';
    ctx.textBaseline = "top";
    // The right-hand state tag is drawn first; the left label then drops detail until both fit.
    const tag = this.active ? (this.outputMode === "worklet" ? "LIVE · WORKLET" : "LIVE · FALLBACK") : "MUTED";
    const tagW = ctx.measureText(tag).width;
    const mode = this.triggerAt >= 0 ? "TRIG BEAT" : "FREE RUN";
    const mic = `mic (${micPosMm[0].toFixed(0)}, ${micPosMm[1].toFixed(0)}, ${micPosMm[2].toFixed(0)}) mm`;
    const room = W - 16 - tagW - 12;
    const left =
      [`${mode} · 2.5 ms/div · ${mic}`, `${mode} · 2.5 ms/div`, mode].find((s) => ctx.measureText(s).width <= room) ?? "";
    if (left) ctx.fillText(left, 8, 7);
    ctx.textAlign = "right";
    ctx.fillText(tag, W - 8, 7);
    ctx.textAlign = "left";
  }
}
