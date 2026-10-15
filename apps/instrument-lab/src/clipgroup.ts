import { GroupPlayer, type PlayerState } from "./audio";
import { el, fmt, fmtSci, hidpi, sciNode, asset } from "./util";

export interface ClipEvidence {
  command?: string | null;
  wall_s?: number | null;
  wall_audio_ratio?: number | null;
  rendered_s?: number | null;
  string_modes?: number | null;
  board_modes?: number | null;
  board_band_hz?: number | null;
  peak_pa?: number | null;
  energy_input_j?: number | null;
  energy_closure_j?: number | null;
  worst_substep_defect_j?: number | null;
  felt_loss_j?: number | null;
  damper_loss_j?: number | null;
  clips?: number | null;
  bem_panels?: number | null;
  bem_band_hz?: string | null;
  bem_audit_max?: number | null;
  bem_min_panels_per_wavelength?: number | null;
  head_frequencies_hz?: Record<string, number[]> | null;
  shell_frequencies_hz?: number[] | null;
  mechanical_rate_hz?: number | null;
  stroke?: string | null;
  string_scale?: string | null;
  soundboard?: string | null;
  microphone?: string | null;
  [key: string]: unknown;
}
export interface Clip {
  label: string;
  group: string;
  audio: string;
  spectrogram: string;
  duration_s: number;
  sample_rate: number;
  channels: number;
  peak_pa: number;
  rms_pa: number;
  envelope_pa: { lo: number[]; hi: number[] };
  evidence: ClipEvidence;
}
export interface Group {
  title: string;
  clips: string[];
  playback_gain_db: number;
  spectrogram: { f_lo: number; f_hi: number; rows: number; ref_db: number; span_db: number; max_duration_s: number };
  band: { label: string; lo_hz: number; hi_hz: number; acoustic_lo_hz?: number; acoustic_hi_hz?: number; note: string };
  pitch_tracks?: Record<string, { t: number[]; f: number[] }>;
}
export interface Manifest {
  commit: string;
  groups: Record<string, Group>;
  clips: Record<string, Clip>;
}

const ICON_PLAY = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M4.5 2.5v11l9-5.5z" fill="currentColor"/></svg>';
const ICON_PAUSE = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M4 2.5h3v11H4zM9 2.5h3v11H9z" fill="currentColor"/></svg>';
const ICON_LOAD = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.5a6.5 6.5 0 1 0 6.5 6.5" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"/></svg>';

const images = new Map<string, Promise<HTMLImageElement>>();
function loadImage(src: string): Promise<HTMLImageElement> {
  let p = images.get(src);
  if (!p) {
    p = new Promise((resolve, reject) => {
      const img = new Image();
      img.decoding = "async";
      img.onload = () => resolve(img);
      img.onerror = () => reject(new Error(src));
      img.src = src;
    });
    images.set(src, p);
  }
  return p;
}

export interface ClipGroupOptions {
  /** Offer a labelled playback-only boost for quiet physical tails. */
  listenCloserDb?: number;
  onTime?: (clipId: string, t: number) => void;
  onSelect?: (clipId: string) => void;
}

/** A comparison of offline renders that share one fixed playback gain. */
export class ClipGroupView {
  readonly player: GroupPlayer;
  private readonly canvas: HTMLCanvasElement;
  private readonly wave: HTMLCanvasElement;
  private readonly readout: HTMLElement;
  private readonly timeLabel: HTMLElement;
  private readonly playBtn: HTMLButtonElement;
  private readonly evidence: HTMLElement;
  private readonly meta: HTMLElement;
  private current: string;
  private hover: { x: number; y: number } | null = null;
  private img: HTMLImageElement | null = null;
  private readonly error: HTMLElement;
  private seen = false;
  private t = 0;

  constructor(host: HTMLElement, private readonly manifest: Manifest, private readonly groupId: string, private readonly opts: ClipGroupOptions = {}) {
    const group = manifest.groups[groupId]!;
    const urls: Record<string, string> = {};
    for (const id of group.clips) urls[id] = asset(manifest.clips[id]!.audio);
    this.current = group.clips[0]!;
    this.player = new GroupPlayer(urls, this.current);

    const root = el("div", { class: "clip" });
    if (group.clips.length > 1) {
      const seg = el("div", { class: "seg", role: "group", "aria-label": "Choose a render" });
      for (const id of group.clips) {
        const c = manifest.clips[id]!;
        const b = el("button", { type: "button", "aria-pressed": String(id === this.current) }, c.label);
        b.addEventListener("click", () => this.select(id, seg));
        seg.append(b);
      }
      root.append(seg);
    }
    const view = el("div", { class: "clip-view", tabindex: "0", "aria-label": "Spectrogram: Space plays or pauses, arrow keys seek" });
    this.canvas = el("canvas", { width: "1200", height: "420", role: "img" }) as HTMLCanvasElement;
    this.readout = el("div", { class: "clip-readout" });
    view.append(this.canvas, this.readout);
    this.wave = el("canvas", { width: "1200", height: "90", "aria-hidden": "true" }) as HTMLCanvasElement;
    this.wave.style.borderTop = "1px solid rgba(236,226,205,0.08)";
    view.append(this.wave);
    root.append(view);

    this.playBtn = el("button", { class: "playbtn", type: "button", "aria-label": "Play" }) as HTMLButtonElement;
    this.playBtn.innerHTML = ICON_PLAY;
    this.playBtn.addEventListener("click", () => this.togglePlay());
    this.timeLabel = el("span", { class: "clip-time", "aria-live": "off" }, "0.00 s");
    const controls = el("div", { class: "clip-controls" }, this.playBtn, this.timeLabel);
    const loop = el("label", { class: "toggle" }, el("input", { type: "checkbox" }), "Loop");
    loop.querySelector("input")!.addEventListener("change", (e) => this.player.setLoop((e.target as HTMLInputElement).checked));
    controls.append(loop);
    if (opts.listenCloserDb) {
      const closer = el("label", { class: "toggle", title: "Raises playback level only; the clips are unchanged" }, el("input", { type: "checkbox" }), `Listen closer (+${opts.listenCloserDb} dB playback)`);
      closer.querySelector("input")!.addEventListener("change", (e) => this.player.setExtraGainDb((e.target as HTMLInputElement).checked ? opts.listenCloserDb! : 0));
      controls.append(closer);
    }
    this.meta = el("div", { class: "clip-meta" });
    this.error = el("p", { class: "small", role: "alert", style: "color: var(--ember); display: none" });
    root.append(controls, this.meta, this.error, el("p", { class: "small muted" }, this.group.band.note));

    this.evidence = el("div", { class: "evidence-body" });
    const details = el("details", { class: "evidence" }, el("summary", {}, "Run evidence"), this.evidence);
    root.append(details);
    host.append(root);

    this.player.onState = (s) => this.renderState(s);
    this.player.onTime = (t) => {
      this.t = t;
      this.timeLabel.textContent = `${t.toFixed(2)} s / ${this.clip.duration_s.toFixed(2)} s`;
      this.draw();
      this.opts.onTime?.(this.current, t);
    };
    this.bindPointer(view);
    new ResizeObserver(() => this.layout()).observe(view);
    this.showClip();
  }

  get clip(): Clip {
    return this.manifest.clips[this.current]!;
  }
  get group(): Group {
    return this.manifest.groups[this.groupId]!;
  }

  private select(id: string, seg: HTMLElement): void {
    this.current = id;
    for (const b of seg.querySelectorAll("button")) b.setAttribute("aria-pressed", String(b.textContent === this.manifest.clips[id]!.label));
    this.player.select(id);
    this.showClip();
    this.opts.onSelect?.(id);
  }

  private togglePlay(): void {
    if (this.player.state === "playing" || this.player.state === "loading") {
      this.player.pause();
      return;
    }
    this.error.style.display = "none";
    void this.player.play().catch((err: unknown) => {
      this.error.textContent = `Could not load or decode this clip (${err instanceof Error ? err.message : String(err)}). Check the connection and press play to retry.`;
      this.error.style.display = "block";
    });
  }

  private renderState(s: PlayerState): void {
    this.playBtn.dataset.state = s;
    this.playBtn.innerHTML = s === "playing" ? ICON_PAUSE : s === "loading" ? ICON_LOAD : ICON_PLAY;
    this.playBtn.setAttribute("aria-label", s === "playing" ? "Pause" : "Play");
  }

  private showClip(): void {
    const c = this.clip;
    this.canvas.setAttribute("aria-label", `Spectrogram of “${c.label}”, ${c.duration_s.toFixed(2)} seconds`);
    this.timeLabel.textContent = `${Math.min(this.t, c.duration_s).toFixed(2)} s / ${c.duration_s.toFixed(2)} s`;
    const ev = c.evidence;
    const items: [string, string | Node, string?][] = [["peak", `${fmt(c.peak_pa)} Pa`], ["duration", `${c.duration_s.toFixed(2)} s`]];
    if (ev.wall_s) items.push(["render", `${fmt(ev.wall_s)} s`]);
    if (ev.wall_audio_ratio) items.push(["slower than real time", `${fmt(ev.wall_audio_ratio)}×`]);
    else if (ev.wall_s && c.duration_s) items.push(["slower than real time", `${fmt(ev.wall_s / c.duration_s)}×`]);
    if (ev.energy_closure_j != null) items.push(["energy closure", sciNode(ev.energy_closure_j), " J"]);
    items.push(["playback gain", `${this.group.playback_gain_db.toFixed(1)} dB`]);
    this.meta.replaceChildren(...items.map(([k, v, u]) => el("span", {}, `${k} `, el("b", {}, v, u ?? ""))));
    this.renderEvidence();
    const fetchImage = (): void => void loadImage(asset(c.spectrogram)).then((img) => {
      if (this.clip === c) {
        this.img = img;
        this.draw();
      }
    }).catch(() => undefined);
    if (this.seen) fetchImage();
    else {
      const io = new IntersectionObserver((e) => {
        if (!e.some((x) => x.isIntersecting)) return;
        io.disconnect();
        this.seen = true;
        if (this.clip === c) fetchImage();
        else this.showClip();
      }, { rootMargin: "400px" });
      io.observe(this.canvas);
    }
    this.drawWave();
    this.draw();
  }

  private renderEvidence(): void {
    const ev = this.clip.evidence;
    const rows: [string, string][] = [];
    const add = (k: string, v: unknown, unit = ""): void => {
      if (v === null || v === undefined || v === "") return;
      rows.push([k, typeof v === "number" ? `${Math.abs(v) < 1e-3 && v !== 0 ? fmtSci(v) : fmt(v)}${unit}` : String(v)]);
    };
    add("mechanical rate", ev.mechanical_rate_hz ? ev.mechanical_rate_hz / 1000 : null, " kHz");
    add("string partials", ev.string_modes);
    add("soundboard modes", ev.board_modes);
    add("board band", ev.board_band_hz, " Hz");
    add("energy input", ev.energy_input_j, " J");
    add("energy closure", ev.energy_closure_j, " J");
    add("worst substep defect", ev.worst_substep_defect_j, " J");
    add("felt loss", ev.felt_loss_j, " J");
    add("damper loss", ev.damper_loss_j, " J");
    add("BEM panels", ev.bem_panels);
    add("BEM band", ev.bem_band_hz ? `${ev.bem_band_hz} Hz` : null);
    add("BEM fit audit (max)", ev.bem_audit_max);
    add("min panels / wavelength", ev.bem_min_panels_per_wavelength);
    add("PCM clips", ev.clips);
    if (ev.head_frequencies_hz) for (const [h, f] of Object.entries(ev.head_frequencies_hz)) add(`head ${h === "0" ? "batter" : "resonant"} modes`, f.map((x) => x.toFixed(1)).join(", "), " Hz");
    if (ev.shell_frequencies_hz) add("shell modes", ev.shell_frequencies_hz.filter((f) => f > 0).slice(0, 12).map((x) => x.toFixed(0)).join(", ") + " …", " Hz");
    const kv = el("div", { class: "kv" }, ...rows.map(([k, v]) => el("div", {}, el("span", {}, k + ": "), v)));
    const cmd = ev.command ? el("pre", {}, (this.clip.evidence.command!.startsWith("--") ? "grand_piano " : "percussion ") + ev.command.replace(/\/private\/tmp\/[^ ]*\//g, "")) : null;
    const notes = [ev.string_scale, ev.soundboard, ev.microphone, ev.stroke].filter(Boolean).map((t) => el("p", { class: "small muted" }, String(t)));
    const parts: HTMLElement[] = [kv, ...notes];
    if (cmd) parts.unshift(cmd);
    this.evidence.replaceChildren(...parts);
  }

  private layout(): void {
    hidpi(this.canvas, 0.35);
    hidpi(this.wave, 0.075, 70);
    this.drawWave();
    this.draw();
  }

  private xOfT(t: number, w: number): number {
    return (t / this.group.spectrogram.max_duration_s) * w;
  }

  private yOfF(f: number, h: number): number {
    const { f_lo, f_hi } = this.group.spectrogram;
    return h - (Math.log(f / f_lo) / Math.log(f_hi / f_lo)) * h;
  }

  private draw(): void {
    const g = this.canvas.getContext("2d");
    if (!g) return;
    const { width: w, height: h } = this.canvas;
    const dpr = w / this.canvas.clientWidth || 1;
    g.fillStyle = "#06070a";
    g.fillRect(0, 0, w, h);
    const span = this.clip.duration_s / this.group.spectrogram.max_duration_s;
    if (this.img) {
      g.imageSmoothingEnabled = true;
      g.drawImage(this.img, 0, 0, w * span, h);
    }
    // Shade what the model does not resolve: it is absent, not faked.
    const band = this.group.band;
    const hiF = band.acoustic_hi_hz ?? band.hi_hz;
    const yHi = this.yOfF(hiF, h);
    g.fillStyle = "rgba(8,9,12,0.55)";
    g.fillRect(0, 0, w, Math.max(0, yHi));
    g.strokeStyle = "rgba(217,169,79,0.85)";
    g.setLineDash([6 * dpr, 5 * dpr]);
    g.lineWidth = 1 * dpr;
    g.beginPath();
    g.moveTo(0, yHi);
    g.lineTo(w, yHi);
    g.stroke();
    if (band.lo_hz > 0 && band.hi_hz !== hiF) {
      const y2 = this.yOfF(band.hi_hz, h);
      g.strokeStyle = "rgba(217,169,79,0.4)";
      g.beginPath();
      g.moveTo(0, y2);
      g.lineTo(w, y2);
      g.stroke();
    }
    g.setLineDash([]);
    g.font = `${11 * dpr}px "JetBrains Mono", monospace`;
    g.fillStyle = "rgba(241,234,219,0.86)";
    const narrow = this.canvas.clientWidth < 560;
    const msg = band.acoustic_hi_hz
      ? (narrow ? `radiation fit to ${hiF} Hz` : `radiation fit verified to ${hiF} Hz · above: computed, transfer unverified`)
      : (narrow ? `board modes to ${band.hi_hz} Hz` : `${band.label} to ${band.hi_hz} Hz · strings resolved to 21.6 kHz`);
    g.fillText(msg, 10 * dpr, Math.max(16 * dpr, yHi - 8 * dpr));
    // Frequency axis.
    g.fillStyle = "rgba(201,194,179,0.7)";
    for (const f of [50, 100, 200, 500, 1000, 2000, 5000, 10000]) {
      const y = this.yOfF(f, h);
      if (y < 14 * dpr || y > h - 4 * dpr) continue;
      g.fillRect(w - 6 * dpr, y, 6 * dpr, 1);
      const label = f >= 1000 ? `${f / 1000}k` : String(f);
      g.fillText(label, w - 10 * dpr - g.measureText(label).width, y + 4 * dpr);
    }
    // Playhead.
    const x = this.xOfT(this.t, w);
    g.fillStyle = "rgba(241,234,219,0.95)";
    g.fillRect(x - dpr, 0, 2 * dpr, h);
    if (this.hover) {
      const hx = this.hover.x * dpr;
      g.fillStyle = "rgba(241,234,219,0.35)";
      g.fillRect(hx, 0, dpr, h);
    }
    this.drawWave(x);
  }

  private drawWave(playX?: number): void {
    const g = this.wave.getContext("2d");
    if (!g) return;
    const { width: w, height: h } = this.wave;
    g.fillStyle = "#08090c";
    g.fillRect(0, 0, w, h);
    const env = this.clip.envelope_pa;
    const peak = Math.max(...this.group.clips.map((id) => this.manifest.clips[id]!.peak_pa));
    const span = this.clip.duration_s / this.group.spectrogram.max_duration_s;
    const n = env.hi.length;
    const mid = h / 2;
    const px = playX ?? this.xOfT(this.t, w);
    for (let i = 0; i < n; i++) {
      const x = (i / n) * w * span;
      const top = mid - (env.hi[i]! / peak) * (mid - 2);
      const bot = mid - (env.lo[i]! / peak) * (mid - 2);
      g.fillStyle = x <= px ? "rgba(217,169,79,0.95)" : "rgba(201,194,179,0.45)";
      g.fillRect(x, top, Math.max(1, (w * span) / n + 0.5), Math.max(1, bot - top));
    }
  }

  private bindPointer(view: HTMLElement): void {
    const pos = (e: PointerEvent): { x: number; y: number } => {
      const r = this.canvas.getBoundingClientRect();
      return { x: e.clientX - r.left, y: e.clientY - r.top };
    };
    this.canvas.addEventListener("pointermove", (e) => {
      const p = pos(e);
      this.hover = p;
      const r = this.canvas.getBoundingClientRect();
      const { f_lo, f_hi, max_duration_s } = this.group.spectrogram;
      const t = (p.x / r.width) * max_duration_s;
      const f = f_lo * Math.pow(f_hi / f_lo, 1 - p.y / r.height);
      this.readout.textContent = `${t.toFixed(3)} s · ${f >= 1000 ? (f / 1000).toFixed(2) + " kHz" : f.toFixed(0) + " Hz"}`;
      this.draw();
    });
    this.canvas.addEventListener("pointerleave", () => {
      this.hover = null;
      this.draw();
    });
    this.canvas.addEventListener("click", (e) => {
      const p = pos(e);
      const r = this.canvas.getBoundingClientRect();
      const t = (p.x / r.width) * this.group.spectrogram.max_duration_s;
      this.player.seek(Math.min(t, this.clip.duration_s));
      if (this.player.state !== "playing") void this.player.play(Math.min(t, this.clip.duration_s)).catch(() => undefined);
    });
    view.addEventListener("keydown", (e) => {
      if (e.key === " ") {
        e.preventDefault();
        this.togglePlay();
      } else if (e.key === "ArrowRight" || e.key === "ArrowLeft") {
        e.preventDefault();
        const step = this.clip.duration_s / 20;
        this.player.seek(Math.max(0, Math.min(this.clip.duration_s, this.player.time + (e.key === "ArrowRight" ? step : -step))));
      }
    });
  }
}
