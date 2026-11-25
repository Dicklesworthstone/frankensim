// 48 kHz live virtual microphone audio engine + triggered oscilloscope + radix-2 FFT
// for the Siemens SST-600 FrankenSim digital twin.

import type { SteamTurbineSim } from "./wasm-pkg/fs_steam_turbine";
import type { TurbineSnapshot } from "./snapshotLayout";

export const AUDIO_SAMPLE_RATE = 48_000;
const RING_LEN = 16_384;
const FFT_SIZE = 2_048;
const TARGET_BUFFERED = 3_072; // ~64 ms
const MAX_CHUNK = 4_096;

export class TurbineAudioEngine {
  private ctx: AudioContext | null = null;
  private node: AudioWorkletNode | null = null;
  private processor: ScriptProcessorNode | null = null;
  private active = false;
  private workletFailed = false;
  private readonly ring = new Float32Array(RING_LEN);
  private writeIdx = 0;
  private filled = 0;
  private reportedAvailable = 0;
  private reportedAtMs = 0;
  private producedSinceReport = 0;
  private readonly fftRe = new Float64Array(FFT_SIZE);
  private readonly fftIm = new Float64Array(FFT_SIZE);
  private readonly hann = new Float64Array(FFT_SIZE);

  constructor() {
    for (let i = 0; i < FFT_SIZE; i++) {
      this.hann[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (FFT_SIZE - 1)));
    }
  }

  public isRunning(): boolean {
    return this.active;
  }

  public async toggle(sim: SteamTurbineSim): Promise<boolean> {
    if (this.active) {
      this.stop();
      return false;
    }
    await this.start(sim);
    return this.active;
  }

  public async start(sim: SteamTurbineSim): Promise<void> {
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
          this.node = new AudioWorkletNode(this.ctx, "fs-turbine-mic-ring", {
            numberOfInputs: 0,
            numberOfOutputs: 1,
            outputChannelCount: [2],
          });
          this.node.port.onmessage = (ev: MessageEvent) => {
            const msg = ev.data as { type?: string; available?: number };
            if (msg?.type === "level") {
              this.reportedAvailable = msg.available ?? 0;
              this.reportedAtMs = performance.now();
              this.producedSinceReport = 0;
            }
          };
          this.node.connect(this.ctx.destination);
        } catch {
          this.workletFailed = true;
          this.node = null;
        }
      }
      if (!this.node) {
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
          this.pushRing(samples);
        };
        proc.connect(this.ctx.destination);
        this.processor = proc;
      }
    }
    this.flush();
    this.active = true;
  }

  public stop(): void {
    this.active = false;
    this.flush();
  }

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
    this.flush();
  }

  private flush(): void {
    this.node?.port.postMessage({ type: "flush" });
    this.reportedAvailable = 0;
    this.reportedAtMs = performance.now();
    this.producedSinceReport = 0;
  }

  private pushRing(samples: Float32Array): void {
    for (let i = 0; i < samples.length; i++) {
      this.ring[this.writeIdx] = samples[i] ?? 0;
      this.writeIdx = (this.writeIdx + 1) & (RING_LEN - 1);
    }
    this.filled = Math.min(RING_LEN, this.filled + samples.length);
  }

  /** Prefill the visualization ring so Scope & FFT render immediately on first Acoustics tab open. */
  public prefillIfNeeded(sim: SteamTurbineSim): void {
    if (this.filled < FFT_SIZE) {
      const samples = sim.render_audio_block(FFT_SIZE);
      this.pushRing(samples);
    }
  }

  /**
   * Called each animation frame. Produces audio for the AudioWorklet when active, OR
   * synthesizes a modest probe block for live SPL / FFT / 3D acoustic lobe metrology while muted.
   */
  public pump(sim: SteamTurbineSim, simRunning: boolean, needProbe: boolean): void {
    if (!simRunning) return;
    if (this.active && this.node) {
      const now = performance.now();
      const elapsedS = Math.max(0, (now - this.reportedAtMs) * 0.001);
      const estBuffered = Math.max(
        0,
        this.reportedAvailable + this.producedSinceReport - elapsedS * AUDIO_SAMPLE_RATE,
      );
      const need = Math.min(MAX_CHUNK, Math.max(0, Math.round(TARGET_BUFFERED - estBuffered)));
      if (need >= 256) {
        const samples = sim.render_audio_block(need);
        const copy = new Float32Array(samples);
        this.pushRing(copy);
        this.node.port.postMessage({ type: "chunk", samples: copy }, [copy.buffer]);
        this.producedSinceReport += copy.length;
      }
    } else if (!this.active && needProbe) {
      // Keep oscilloscope, FFT, 3D acoustic lobes, and ISO 532-1 SPL/Sone meters live while muted
      const samples = sim.render_audio_block(1024);
      this.pushRing(samples);
    }
  }

  public drawScopeAndFft(
    scopeCanvas: HTMLCanvasElement,
    fftCanvas: HTMLCanvasElement,
    snap: TurbineSnapshot,
  ): void {
    this.drawScope(scopeCanvas);
    this.drawFft(fftCanvas, snap);
  }

  private drawScope(canvas: HTMLCanvasElement): void {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = canvas.clientWidth || 320;
    const h = canvas.clientHeight || 130;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
    }
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    ctx.fillStyle = "#080c12";
    ctx.fillRect(0, 0, w, h);

    // Subtle oscilloscope graticule
    ctx.strokeStyle = "rgba(116, 158, 190, 0.14)";
    ctx.lineWidth = 1;
    for (let i = 1; i < 4; i++) {
      const y = (h * i) / 4;
      ctx.beginPath();
      ctx.moveTo(0, y);
      ctx.lineTo(w, y);
      ctx.stroke();
    }
    for (let i = 1; i < 8; i++) {
      const x = (w * i) / 8;
      ctx.beginPath();
      ctx.moveTo(x, 0);
      ctx.lineTo(x, h);
      ctx.stroke();
    }

    const win = 960; // 20 ms (1 full 50 Hz revolution!)
    if (this.filled < win + 384) return;
    const startSearch = (this.writeIdx - win - 384 + RING_LEN) & (RING_LEN - 1);
    let trig = startSearch;
    for (let i = 1; i < 384; i++) {
      const i0 = (startSearch + i - 1) & (RING_LEN - 1);
      const i1 = (startSearch + i) & (RING_LEN - 1);
      if ((this.ring[i0] ?? 0) <= 0 && (this.ring[i1] ?? 0) > 0) {
        trig = i1;
        break;
      }
    }

    let peak = 0.08;
    for (let i = 0; i < win; i++) {
      const v = Math.abs(this.ring[(trig + i) & (RING_LEN - 1)] ?? 0);
      if (v > peak) peak = v;
    }

    ctx.strokeStyle = "#5ad2f4";
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    for (let i = 0; i < win; i++) {
      const v = this.ring[(trig + i) & (RING_LEN - 1)] ?? 0;
      const x = (i / (win - 1)) * w;
      const y = h * 0.5 - (v / (peak * 1.15)) * (h * 0.42);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    ctx.fillStyle = "rgba(180, 205, 224, 0.72)";
    ctx.font = '10px "JetBrains Mono", monospace';
    ctx.fillText("20 ms (1 rev @ 50 Hz)", 6, 12);
    ctx.fillText(this.active ? "LIVE 48 kHz" : "PROBE (MUTED)", w - 88, 12);
  }

  private drawFft(canvas: HTMLCanvasElement, snap: TurbineSnapshot): void {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = canvas.clientWidth || 320;
    const h = canvas.clientHeight || 130;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
    }
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    ctx.fillStyle = "#080c12";
    ctx.fillRect(0, 0, w, h);

    if (this.filled < FFT_SIZE) return;
    const base = (this.writeIdx - FFT_SIZE + RING_LEN) & (RING_LEN - 1);
    for (let i = 0; i < FFT_SIZE; i++) {
      this.fftRe[i] = (this.ring[(base + i) & (RING_LEN - 1)] ?? 0) * (this.hann[i] ?? 1);
      this.fftIm[i] = 0;
    }
    fftRadix2(this.fftRe, this.fftIm);

    // Log-frequency display from 15 Hz to 6,000 Hz so 23 Hz oil whip, 50 Hz 1x, AND 2.2 kHz BPF are ALL visible!
    const fMin = 15;
    const fMax = 6_000;
    const lnMin = Math.log(fMin);
    const lnSpan = Math.log(fMax) - lnMin;
    const padL = 24;
    const padR = 8;
    const padT = 16;
    const padB = 18;
    const pw = w - padL - padR;
    const ph = h - padT - padB;
    const xOfFreq = (f: number) => padL + ((Math.log(Math.max(fMin, Math.min(fMax, f))) - lnMin) / lnSpan) * pw;

    // Frequency gridlines
    ctx.strokeStyle = "rgba(116, 158, 190, 0.14)";
    ctx.fillStyle = "rgba(154, 177, 198, 0.65)";
    ctx.font = '9px "JetBrains Mono", monospace';
    for (const f of [25, 50, 100, 300, 1000, 3000]) {
      const x = xOfFreq(f);
      ctx.beginPath();
      ctx.moveTo(x, padT);
      ctx.lineTo(x, padT + ph);
      ctx.stroke();
      ctx.fillText(f >= 1000 ? `${f / 1000}k` : `${f}`, x - 7, h - 5);
    }

    // Physical marker lines: 0.46x oil whip, 1x sync, LP BPF
    const fSync = Math.max(5, snap.rpm / 60);
    const markers: Array<{ f: number; label: string; color: string }> = [
      { f: snap.dominantWhirlHz, label: "whirl", color: snap.spectralAbscissaPerS > -2 ? "#ff6464" : "rgba(255,130,100,0.45)" },
      { f: fSync, label: "1×", color: "#5ad2f4" },
      { f: snap.bpfLpHz, label: "LP BPF", color: "#eab663" },
    ];
    for (const m of markers) {
      if (m.f >= fMin && m.f <= fMax) {
        const x = xOfFreq(m.f);
        ctx.strokeStyle = m.color;
        ctx.setLineDash([2, 2]);
        ctx.beginPath();
        ctx.moveTo(x, padT);
        ctx.lineTo(x, padT + ph);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = m.color;
        ctx.fillText(m.label, Math.min(w - 36, x + 2), padT + 9);
      }
    }

    // Draw magnitude curve (dB)
    const half = FFT_SIZE >> 1;
    const binHz = AUDIO_SAMPLE_RATE / FFT_SIZE;
    ctx.strokeStyle = "#7dd3fc";
    ctx.fillStyle = "rgba(92, 164, 196, 0.22)";
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    let started = false;
    let lastX = padL;
    for (let k = 1; k < half; k++) {
      const f = k * binHz;
      if (f < fMin || f > fMax) continue;
      const mag = Math.hypot(this.fftRe[k] ?? 0, this.fftIm[k] ?? 0) / (FFT_SIZE * 0.25);
      const db = Math.max(-75, Math.min(0, 20 * Math.log10(mag + 1e-6)));
      const x = xOfFreq(f);
      const y = padT + ((0 - db) / 75) * ph;
      lastX = x;
      if (!started) {
        ctx.moveTo(x, padT + ph);
        ctx.lineTo(x, y);
        started = true;
      } else {
        ctx.lineTo(x, y);
      }
    }
    if (started) {
      ctx.lineTo(lastX, padT + ph);
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    }
    ctx.fillStyle = "rgba(180, 205, 224, 0.72)";
    ctx.fillText("FFT dB (log Hz)", 6, 11);
  }
}

function fftRadix2(re: Float64Array, im: Float64Array): void {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      const tr = re[i] ?? 0;
      const ti = im[i] ?? 0;
      re[i] = re[j] ?? 0;
      im[i] = im[j] ?? 0;
      re[j] = tr;
      im[j] = ti;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len;
    const wlenRe = Math.cos(ang);
    const wlenIm = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wRe = 1;
      let wIm = 0;
      const half = len >> 1;
      for (let j = 0; j < half; j++) {
        const uRe = re[i + j] ?? 0;
        const uIm = im[i + j] ?? 0;
        const vRe = (re[i + j + half] ?? 0) * wRe - (im[i + j + half] ?? 0) * wIm;
        const vIm = (re[i + j + half] ?? 0) * wIm + (im[i + j + half] ?? 0) * wRe;
        re[i + j] = uRe + vRe;
        im[i + j] = uIm + vIm;
        re[i + j + half] = uRe - vRe;
        im[i + j + half] = uIm - vIm;
        const nextWRe = wRe * wlenRe - wIm * wlenIm;
        wIm = wRe * wlenIm + wIm * wlenRe;
        wRe = nextWRe;
      }
    }
  }
}
