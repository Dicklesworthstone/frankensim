import { claimPlayback, getContext, registerPlayer, resumeContext } from "./audio";
import { dpr, el, fmt, hidpi } from "./util";

interface MixerMeta {
  sample_rate: number;
  frames: number;
  channels: number;
  scales_pa: number[];
  frequency_hz: (number | null)[];
  mode_index: number[];
  energy_share: number[];
  total_modes: number;
  reconstruction_max_rel_error: number;
  note: string;
}

/**
 * The soundboard mixer: each channel is one loaded-board mode's signed pressure at
 * the microphone, taken from the solver; the last is the exact remainder. All
 * enabled channels are summed by WebAudio, so the full mix IS the rendered note.
 */
export class SoundboardMixer {
  private meta: MixerMeta | null = null;
  private data: Float32Array<ArrayBuffer>[] = [];
  private buffers: AudioBuffer[] = [];
  private enabled: boolean[] = [];
  private gains: GainNode[] = [];
  private sources: AudioBufferSourceNode[] = [];
  private master: GainNode | null = null;
  private playing = false;
  private readonly bars: HTMLCanvasElement;
  private readonly wave: HTMLCanvasElement;
  private readonly readout: HTMLElement;
  private readonly playBtn: HTMLButtonElement;
  private loading: Promise<void> | null = null;
  private hover = -1;
  private list!: HTMLElement;
  private boost = 1;

  constructor(private readonly host: HTMLElement) {
    this.bars = el("canvas", { width: "900", height: "260", "aria-label": "Soundboard modes by frequency; tap a bar to mute or unmute it" });
    this.wave = el("canvas", { width: "900", height: "150", "aria-hidden": "true" });
    this.readout = el("p", { class: "small muted num", "aria-live": "polite" });
    this.playBtn = el("button", { class: "btn primary sm", type: "button" }, "▶ Play the mix");
    const all = el("button", { class: "btn sm", type: "button" }, "All modes");
    const low = el("button", { class: "btn sm", type: "button" }, "Below 300 Hz only");
    const top3 = el("button", { class: "btn sm", type: "button" }, "Three strongest");
    const rest = el("button", { class: "btn sm", type: "button" }, "Remove the three strongest");
    this.playBtn.addEventListener("click", () => void this.toggle());
    all.addEventListener("click", () => this.setAll((i) => true));
    low.addEventListener("click", () => this.setAll((i) => (this.meta!.frequency_hz[i] ?? 1e9) < 300));
    top3.addEventListener("click", () => this.setAll((i) => i < 3));
    rest.addEventListener("click", () => this.setAll((i) => i >= 3));
    const chart = el("div", { class: "chart", style: "margin-top: 12px" }, this.bars);
    const waveBox = el("div", { class: "chart", style: "margin-top: 8px" }, this.wave);
    this.list = el("div", { class: "mode-toggles", role: "group", "aria-label": "Soundboard mode channels" });
    const listBox = el("details", { class: "evidence" }, el("summary", {}, "Mode list (keyboard)"), this.list);
    host.append(el("div", { class: "btn-row", style: "margin-top: 14px" }, this.playBtn, all, low, top3, rest), chart, waveBox, this.readout, listBox);
    this.bars.addEventListener("click", (e) => this.toggleAt(e));
    this.bars.addEventListener("pointermove", (e) => { this.hover = this.indexAt(e); this.drawBars(); });
    this.bars.addEventListener("pointerleave", () => { this.hover = -1; this.drawBars(); });
    new ResizeObserver(() => { this.drawBars(); this.drawWave(); }).observe(host);
    registerPlayer(() => this.stop());
    this.readout.textContent = "Loading the per-mode pressure (4 MB)…";
    void this.load();
  }

  private load(): Promise<void> {
    this.loading ??= (async () => {
      const [meta, bytes] = await Promise.all([
        fetch("/data/soundboard-modes.json").then((r) => r.json() as Promise<MixerMeta>),
        fetch("/data/soundboard-modes.i16").then((r) => r.arrayBuffer()),
      ]);
      this.meta = meta;
      const pcm = new Int16Array(bytes);
      for (let c = 0; c < meta.channels; c++) {
        const ch = new Float32Array(new ArrayBuffer(meta.frames * 4));
        const s = meta.scales_pa[c]! / 32767;
        for (let i = 0; i < meta.frames; i++) ch[i] = pcm[c * meta.frames + i]! * s;
        this.data.push(ch);
      }
      this.enabled = new Array(meta.channels).fill(true);
      this.list.replaceChildren(...meta.frequency_hz.map((f, i) => {
        const box = el("input", { type: "checkbox" });
        box.checked = true;
        box.addEventListener("change", () => { this.enabled[i] = box.checked; this.update(); });
        return el("label", { class: "toggle" }, box, f ? `${f.toFixed(0)} Hz · ${(meta.energy_share[i]! * 100).toFixed(1)}%` : `the remaining ${meta.total_modes - meta.channels + 1} modes`);
      }));
      let peak = 0;
      for (let i = 0; i < meta.frames; i++) {
        let sum = 0;
        for (const ch of this.data) sum += ch[i]!;
        peak = Math.max(peak, Math.abs(sum));
      }
      this.boost = 0.85 / peak;
      this.update();
    })();
    return this.loading;
  }

  private indexAt(e: MouseEvent): number {
    if (!this.meta) return -1;
    const r = this.bars.getBoundingClientRect();
    return Math.min(this.meta.channels - 1, Math.max(0, Math.floor(((e.clientX - r.left) / r.width) * this.meta.channels)));
  }

  private toggleAt(e: MouseEvent): void {
    const i = this.indexAt(e);
    if (i < 0) return;
    this.enabled[i] = !this.enabled[i];
    this.update();
  }

  private setAll(pred: (i: number) => boolean): void {
    if (!this.meta) return;
    this.enabled = this.enabled.map((_, i) => pred(i));
    this.update();
  }

  private update(): void {
    this.list?.querySelectorAll<HTMLInputElement>("input").forEach((b, i) => { b.checked = !!this.enabled[i]; });
    this.gains.forEach((g, i) => g.gain.setTargetAtTime(this.enabled[i] ? 1 : 0, getContext().currentTime, 0.015));
    this.drawBars();
    this.drawWave();
    if (!this.meta) return;
    const on = this.enabled.filter(Boolean).length;
    let sumSq = 0, partsRms = 0;
    const n = this.meta.frames;
    for (let c = 0; c < this.data.length; c++) {
      if (!this.enabled[c]) continue;
      let ss = 0;
      const ch = this.data[c]!;
      for (let i = 0; i < n; i += 4) ss += ch[i]! * ch[i]!;
      partsRms += Math.sqrt(ss / (n / 4));
    }
    for (let i = 0; i < n; i += 4) {
      let s = 0;
      for (let c = 0; c < this.data.length; c++) if (this.enabled[c]) s += this.data[c]![i]!;
      sumSq += s * s;
    }
    const mixRms = Math.sqrt(sumSq / (n / 4));
    this.readout.textContent = `${on} of ${this.meta.channels} channels on (${this.meta.total_modes} loaded modes; the last channel is the sum of the ${this.meta.total_modes - this.meta.channels + 1} weakest). Mix RMS ${fmt(mixRms)} Pa vs ${fmt(partsRms)} Pa if the modes never cancelled: ${on ? Math.round((1 - mixRms / Math.max(partsRms, 1e-30)) * 100) : 0}% of the modal pressure cancels at the microphone. Full-mix reconstruction error ${fmt(this.meta.reconstruction_max_rel_error)}.`;
  }

  private drawBars(): void {
    hidpi(this.bars, 0.26, 150);
    const g = this.bars.getContext("2d");
    if (!g) return;
    const d = dpr();
    const { width: W, height: H } = this.bars;
    g.fillStyle = "#07080b";
    g.fillRect(0, 0, W, H);
    if (!this.meta) return;
    const n = this.meta.channels;
    const bw = W / n;
    const maxShare = Math.max(...this.meta.energy_share);
    g.font = `${9.5 * d}px "JetBrains Mono", monospace`;
    g.textAlign = "center";
    for (let i = 0; i < n; i++) {
      const share = this.meta.energy_share[i]!;
      const h = Math.max(2 * d, Math.sqrt(share / maxShare) * (H - 34 * d));
      const x = i * bw;
      const on = this.enabled[i];
      g.fillStyle = on ? (i === n - 1 ? "rgba(201,194,179,0.55)" : "#d9a94f") : "rgba(201,194,179,0.14)";
      if (i === this.hover) g.fillStyle = on ? "#f0c670" : "rgba(201,194,179,0.3)";
      g.fillRect(x + bw * 0.15, H - 20 * d - h, bw * 0.7, h);
      const f = this.meta.frequency_hz[i];
      const every = Math.max(1, Math.ceil((30 * d) / bw));
      if (i % every === 0 || i === n - 1) {
        g.fillStyle = "rgba(201,194,179,0.7)";
        const label = f ? (f >= 1000 ? (f / 1000).toFixed(1) + "k" : f.toFixed(0)) : "rest";
        const half = g.measureText(label).width / 2 + 2 * d;
        g.fillText(label, Math.min(Math.max(x + bw / 2, half), W - half), H - 8 * d);
      }
    }
    g.textAlign = "left";
    g.fillStyle = "rgba(143,137,125,0.95)";
    g.fillText(W / d < 560 ? "mode frequency [Hz] · height ∝ √energy · tap to mute" : "loaded-board mode reference frequency [Hz] · bar height ∝ √(energy share) · tap to mute", 8 * d, 12 * d);
  }

  private drawWave(): void {
    hidpi(this.wave, 0.15, 90);
    const g = this.wave.getContext("2d");
    if (!g) return;
    const { width: W, height: H } = this.wave;
    g.fillStyle = "#07080b";
    g.fillRect(0, 0, W, H);
    if (!this.meta) return;
    const n = this.meta.frames;
    const cols = Math.floor(W);
    const mid = H / 2;
    const per = n / cols;
    const scale = (this.boost * (H / 2 - 4)) / 0.85;
    for (let x = 0; x < cols; x++) {
      let lo = Infinity, hi = -Infinity;
      for (let i = Math.floor(x * per); i < Math.floor((x + 1) * per); i += 2) {
        let s = 0;
        for (let c = 0; c < this.data.length; c++) if (this.enabled[c]) s += this.data[c]![i]!;
        if (s < lo) lo = s;
        if (s > hi) hi = s;
      }
      g.fillStyle = "rgba(217,169,79,0.85)";
      g.fillRect(x, mid - hi * scale, 1, Math.max(1, (hi - lo) * scale));
    }
  }

  private async toggle(): Promise<void> {
    if (this.playing) return this.stop();
    const ctx = await resumeContext();
    await this.load();
    if (!this.meta) return;
    claimPlayback(() => this.stop());
    if (this.buffers.length === 0) {
      for (const ch of this.data) {
        const b = ctx.createBuffer(1, ch.length, this.meta.sample_rate);
        b.copyToChannel(ch, 0);
        this.buffers.push(b);
      }
    }
    this.master = ctx.createGain();
    this.master.gain.value = this.boost;
    this.master.connect(ctx.destination);
    const at = ctx.currentTime + 0.06;
    this.gains = [];
    this.sources = this.buffers.map((b, i) => {
      const src = ctx.createBufferSource();
      src.buffer = b;
      src.loop = true;
      const g = ctx.createGain();
      g.gain.value = this.enabled[i] ? 1 : 0;
      src.connect(g).connect(this.master!);
      src.start(at);
      this.gains.push(g);
      return src;
    });
    this.playing = true;
    this.playBtn.textContent = "Stop";
  }

  private stop(): void {
    for (const s of this.sources) {
      try { s.stop(); } catch { /* stopped */ }
      s.disconnect();
    }
    this.sources = [];
    this.gains = [];
    this.master?.disconnect();
    this.master = null;
    this.playing = false;
    this.playBtn.textContent = "▶ Play the mix";
  }
}
