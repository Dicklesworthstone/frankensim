// Two running spectrograms stacked over one time axis: the studio render and the live engines,
// read from AnalyserNodes on each signal before the A/B switch, so both keep scrolling and the
// one being heard is lit. Log frequency, 40 Hz to 16 kHz.

const F_LO = 40;
const F_HI = 16000;
const DB_LO = -112;
const DB_HI = -22;

export type Strip = "studio" | "live";

export interface SpectroOptions {
  /** Pixels (CSS) per second of scroll. */
  pxPerSecond?: number;
  /** Bar and beat positions: returns the musical position (in beats) at a context time, or null. */
  beatAt?: (time: number) => number | null;
  beatsPerBar?: number;
}

export class DualSpectrogram {
  active: Strip = "studio";
  private ctx2d: CanvasRenderingContext2D;
  private bins: { studio: Float32Array<ArrayBuffer>; live: Float32Array<ArrayBuffer> };
  private rows: Int32Array = new Int32Array(0);
  private rowsEnd: Int32Array = new Int32Array(0);
  private lut = magma();
  private raf = 0;
  private lastT = 0;
  private carry = 0;
  private lastBeat: number | null = null;
  private stripH = 0;
  private gap = 0;
  private ro: ResizeObserver;

  constructor(private canvas: HTMLCanvasElement, private analysers: { studio: AnalyserNode; live: AnalyserNode }, private opts: SpectroOptions = {}) {
    this.ctx2d = canvas.getContext("2d", { alpha: false })!;
    this.bins = { studio: new Float32Array(analysers.studio.frequencyBinCount), live: new Float32Array(analysers.live.frequencyBinCount) };
    this.ro = new ResizeObserver(() => this.layout());
    this.ro.observe(canvas);
    this.layout();
  }

  private layout(): void {
    const r = this.canvas.getBoundingClientRect();
    const d = Math.min(2, window.devicePixelRatio || 1);
    const w = Math.max(2, Math.round(r.width * d)), h = Math.max(2, Math.round(r.height * d));
    // a new spectrogram on a canvas an earlier one already sized (Play after Stop) still needs its rows
    if (this.stripH > 0 && this.canvas.width === w && this.canvas.height === h) return;
    this.canvas.width = w;
    this.canvas.height = h;
    this.gap = Math.round(10 * d);
    // 0 on a canvas not laid out (2 px tall): a negative length would throw building the rows
    this.stripH = Math.max(0, Math.floor((h - this.gap) / 2));
    const a = this.analysers.studio;
    const hz = a.context.sampleRate / a.fftSize;
    this.rows = new Int32Array(this.stripH);
    this.rowsEnd = new Int32Array(this.stripH);
    for (let y = 0; y < this.stripH; y++) {
      const f0 = F_LO * Math.pow(F_HI / F_LO, 1 - (y + 1) / this.stripH);
      const f1 = F_LO * Math.pow(F_HI / F_LO, 1 - y / this.stripH);
      this.rows[y] = Math.max(1, Math.floor(f0 / hz));
      this.rowsEnd[y] = Math.max(this.rows[y]! + 1, Math.ceil(f1 / hz));
    }
    this.ctx2d.fillStyle = "#07080b";
    this.ctx2d.fillRect(0, 0, w, h);
  }

  start(): void {
    if (this.raf) return;
    this.lastT = performance.now();
    this.carry = 0;
    const tick = (now: number): void => {
      this.raf = requestAnimationFrame(tick);
      this.step(now);
    };
    this.raf = requestAnimationFrame(tick);
  }

  stop(): void {
    cancelAnimationFrame(this.raf);
    this.raf = 0;
  }

  dispose(): void {
    this.stop();
    this.ro.disconnect();
  }

  private step(now: number): void {
    const c = this.canvas, g = this.ctx2d;
    const d = c.width / Math.max(1, c.getBoundingClientRect().width);
    const pps = (this.opts.pxPerSecond ?? 90) * d;
    // at most one canvas width: after a tab was hidden for minutes, one frame would otherwise build an
    // image tens of thousands of pixels wide
    this.carry = Math.min(c.width, this.carry + ((now - this.lastT) / 1000) * pps);
    this.lastT = now;
    if (this.stripH < 1) { this.carry = 0; return; } // nothing to draw into (a canvas not laid out yet)
    const dx = Math.floor(this.carry);
    if (dx < 1) return;
    this.carry -= dx;
    if (dx < c.width) g.drawImage(c, dx, 0, c.width - dx, c.height, 0, 0, c.width - dx, c.height);
    g.fillStyle = "#07080b";
    g.fillRect(c.width - dx, this.stripH, dx, this.gap);
    for (const s of ["studio", "live"] as const) {
      const arr = this.bins[s];
      this.analysers[s].getFloatFrequencyData(arr);
      const y0 = s === "studio" ? 0 : this.stripH + this.gap;
      const img = g.createImageData(dx, this.stripH);
      const lit = s === this.active;
      for (let y = 0; y < this.stripH; y++) {
        let m = -Infinity;
        const a = this.rows[y]!, b = Math.min(this.rowsEnd[y]!, arr.length);
        for (let k = a; k < b; k++) if (arr[k]! > m) m = arr[k]!;
        let v = (m - DB_LO) / (DB_HI - DB_LO);
        v = v < 0 ? 0 : v > 1 ? 1 : v;
        const i = Math.round(v * 255) * 3;
        let r = this.lut[i]!, gr = this.lut[i + 1]!, bl = this.lut[i + 2]!;
        if (!lit) {
          const l = 0.3 * r + 0.59 * gr + 0.11 * bl;
          r = 0.46 * (0.6 * l + 0.4 * r); gr = 0.46 * (0.6 * l + 0.4 * gr); bl = 0.46 * (0.6 * l + 0.4 * bl);
        }
        for (let x = 0; x < dx; x++) {
          const o = (y * dx + x) * 4;
          img.data[o] = r; img.data[o + 1] = gr; img.data[o + 2] = bl; img.data[o + 3] = 255;
        }
      }
      g.putImageData(img, c.width - dx, y0);
    }
    // bar and beat lines at the moment they pass
    const beat = this.opts.beatAt?.(this.analysers.studio.context.currentTime) ?? null;
    if (beat !== null && this.lastBeat !== null && Math.floor(beat) !== Math.floor(this.lastBeat) && beat > this.lastBeat) {
      const bar = Math.floor(beat) % (this.opts.beatsPerBar ?? 4) === 0;
      g.fillStyle = bar ? "rgba(241,234,219,0.42)" : "rgba(241,234,219,0.12)";
      const w = Math.max(1, Math.round(d * (bar ? 1.5 : 1)));
      g.fillRect(c.width - w, 0, w, this.stripH);
      g.fillRect(c.width - w, this.stripH + this.gap, w, this.stripH);
    }
    this.lastBeat = beat;
  }
}

/** A 256-entry magma-like colour ramp (dark violet, red, orange, pale yellow). */
function magma(): Uint8ClampedArray {
  const stops: [number, number, number, number][] = [
    [0, 7, 8, 11], [0.14, 28, 16, 48], [0.3, 82, 22, 99], [0.46, 150, 38, 102],
    [0.62, 214, 70, 72], [0.78, 246, 133, 47], [0.9, 252, 199, 96], [1, 252, 246, 196],
  ];
  const lut = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    let k = 0;
    while (k < stops.length - 2 && stops[k + 1]![0] < t) k++;
    const [t0, r0, g0, b0] = stops[k]!, [t1, r1, g1, b1] = stops[k + 1]!;
    const u = (t - t0) / (t1 - t0);
    lut[i * 3] = r0 + (r1 - r0) * u;
    lut[i * 3 + 1] = g0 + (g1 - g0) * u;
    lut[i * 3 + 2] = b0 + (b1 - b0) * u;
  }
  return lut;
}
