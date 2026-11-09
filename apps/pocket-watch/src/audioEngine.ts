// Real-time 48 kHz WebAudio Virtual Microphone Engine & Horological Oscilloscope
// Driven directly by `PocketWatchSim::render_audio_block` in WASM.

import type { PocketWatchSim } from "./wasm-pkg/fs_pocket_watch";

export class VirtualMicrophoneAudio {
  private ctx: AudioContext | null = null;
  private processor: ScriptProcessorNode | null = null;
  private active = false;
  private readonly scopeRing = new Float32Array(2048);
  private scopeWriteIdx = 0;
  public lastRmsSplDb = 32.0;

  public isRunning(): boolean {
    return this.active;
  }

  public async toggle(sim: PocketWatchSim): Promise<boolean> {
    if (this.active) {
      this.stop();
      return false;
    }
    await this.start(sim);
    return true;
  }

  public async start(sim: PocketWatchSim): Promise<void> {
    if (!this.ctx) {
      this.ctx = new AudioContext({ sampleRate: 48000, latencyHint: "interactive" });
    }
    if (this.ctx.state === "suspended") {
      await this.ctx.resume();
    }
    if (!this.processor) {
      const blockSize = 1024;
      this.processor = this.ctx.createScriptProcessor(blockSize, 0, 2);
      this.processor.onaudioprocess = (ev: AudioProcessingEvent) => {
        if (!this.active) {
          return;
        }
        const left = ev.outputBuffer.getChannelData(0);
        const right = ev.outputBuffer.getChannelData(1);
        const samples = sim.render_audio_block(left.length);
        let sumSq = 0;
        for (let i = 0; i < left.length; i++) {
          const s = samples[i] ?? 0;
          left[i] = s;
          right[i] = s;
          sumSq += s * s;
          this.scopeRing[this.scopeWriteIdx] = s;
          this.scopeWriteIdx = (this.scopeWriteIdx + 1) & (this.scopeRing.length - 1);
        }
        const rms = Math.sqrt(sumSq / left.length);
        // Reference 20 uPa acoustic threshold -> dB SPL
        const pa = rms * 0.25;
        const db = 20 * Math.log10(Math.max(pa, 2.0e-5) / 2.0e-5);
        this.lastRmsSplDb = 0.8 * this.lastRmsSplDb + 0.2 * Math.min(110, Math.max(26, db));
      };
      this.processor.connect(this.ctx.destination);
    }
    this.active = true;
  }

  public stop(): void {
    this.active = false;
  }

  /**
   * When live WebAudio speaker output is muted, still drain a silent/scope block
   * each visual frame so the 3-pulse horological oscilloscope and Timegrapher stay live!
   */
  public pumpSilentScopeBlock(sim: PocketWatchSim): void {
    if (this.active) {
      return;
    }
    const samples = sim.render_audio_block(512);
    let sumSq = 0;
    for (let i = 0; i < samples.length; i++) {
      const s = samples[i] ?? 0;
      sumSq += s * s;
      this.scopeRing[this.scopeWriteIdx] = s;
      this.scopeWriteIdx = (this.scopeWriteIdx + 1) & (this.scopeRing.length - 1);
    }
    const rms = Math.sqrt(sumSq / samples.length);
    const pa = rms * 0.25;
    const db = 20 * Math.log10(Math.max(pa, 2.0e-5) / 2.0e-5);
    this.lastRmsSplDb = 0.8 * this.lastRmsSplDb + 0.2 * Math.min(110, Math.max(26, db));
  }

  public renderOscilloscope(canvas: HTMLCanvasElement, micPosMm: [number, number, number]): void {
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    const W = canvas.width;
    const H = canvas.height;

    ctx.fillStyle = "#06080c";
    ctx.fillRect(0, 0, W, H);

    // Subtle graticule lines
    ctx.strokeStyle = "rgba(217, 177, 95, 0.12)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, H / 2);
    ctx.lineTo(W, H / 2);
    ctx.moveTo(W / 3, 0);
    ctx.lineTo(W / 3, H);
    ctx.moveTo((2 * W) / 3, 0);
    ctx.lineTo((2 * W) / 3, H);
    ctx.stroke();

    // Trigger on peak impulse in the ring buffer so the 3-pulse horological signature is clear
    const N = this.scopeRing.length;
    ctx.strokeStyle = this.active ? "#5ef2b8" : "#f4d284";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let x = 0; x < W; x++) {
      const idx = (this.scopeWriteIdx + Math.floor((x / W) * N)) & (N - 1);
      const val = this.scopeRing[idx] ?? 0;
      const y = H / 2 - val * (H * 1.35);
      if (x === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    ctx.fillStyle = "rgba(236, 232, 223, 0.72)";
    ctx.font = '10px "JetBrains Mono", monospace';
    ctx.fillText(
      `MIC (${micPosMm[0].toFixed(0)}, ${micPosMm[1].toFixed(0)}, ${micPosMm[2].toFixed(0)}) mm · ${this.active ? "LIVE 48 kHz" : "SCOPE ONLY (CLICK TO UNMUTE)"}`,
      7,
      13,
    );
  }
}
