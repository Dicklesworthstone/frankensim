// One AudioContext for the whole page (48 kHz, the rate of every offline render
// and of both live engines), created inside the first user gesture.

let context: AudioContext | null = null;

export function getContext(): AudioContext {
  if (!context) {
    const Ctor: typeof AudioContext =
      window.AudioContext ?? (window as unknown as { webkitAudioContext: typeof AudioContext }).webkitAudioContext;
    context = new Ctor({ sampleRate: 48_000, latencyHint: "interactive" });
  }
  return context;
}

export async function resumeContext(): Promise<AudioContext> {
  const ctx = getContext();
  if (ctx.state !== "running") await ctx.resume();
  return ctx;
}

/** Only one offline clip plays at a time; live engines are not interrupted. */
const claims = new Set<() => void>();
export function claimPlayback(stop: () => void): void {
  for (const other of claims) if (other !== stop) other();
}
export function registerPlayer(stop: () => void): void {
  claims.add(stop);
}

const buffers = new Map<string, Promise<AudioBuffer>>();
export function loadBuffer(url: string): Promise<AudioBuffer> {
  let pending = buffers.get(url);
  if (!pending) {
    pending = fetch(url)
      .then((res) => {
        if (!res.ok) throw new Error(`${url}: HTTP ${res.status}`);
        return res.arrayBuffer();
      })
      .then((bytes) => getContext().decodeAudioData(bytes));
    pending.catch(() => buffers.delete(url));
    buffers.set(url, pending);
  }
  return pending;
}

export type PlayerState = "idle" | "loading" | "playing" | "paused";

/**
 * Plays one clip out of a comparison group. Switching clips keeps the playhead,
 * so an A/B is heard at the same instant of the same gesture. Nothing is
 * level-matched here: every clip of a group was encoded with the same gain.
 */
export class GroupPlayer {
  private source: AudioBufferSourceNode | null = null;
  private gain: GainNode | null = null;
  private startedAt = 0;
  private offset = 0;
  private token = 0;
  private raf = 0;
  state: PlayerState = "idle";
  loop = false;
  extraGainDb = 0;
  current: string;
  duration = 0;
  onState: (state: PlayerState) => void = () => {};
  onTime: (seconds: number, duration: number) => void = () => {};

  constructor(private readonly urls: Record<string, string>, first: string) {
    this.current = first;
    registerPlayer(this.stopHandle);
  }

  private readonly stopHandle = (): void => this.pause();

  get time(): number {
    if (this.state !== "playing") return this.offset;
    const t = this.offset + (getContext().currentTime - this.startedAt);
    return this.loop && this.duration > 0 ? t % this.duration : Math.min(t, this.duration);
  }

  private setState(state: PlayerState): void {
    this.state = state;
    this.onState(state);
  }

  async play(from?: number): Promise<void> {
    const token = ++this.token;
    const ctx = await resumeContext();
    this.setState("loading");
    let buffer: AudioBuffer;
    try {
      buffer = await loadBuffer(this.urls[this.current]!);
    } catch (error) {
      if (token === this.token) this.setState("idle");
      throw error;
    }
    if (token !== this.token) return;
    claimPlayback(this.stopHandle);
    this.teardown();
    this.duration = buffer.duration;
    let start = from ?? this.offset;
    if (start >= buffer.duration - 0.02) start = 0;
    const source = ctx.createBufferSource();
    const gain = ctx.createGain();
    gain.gain.value = 10 ** (this.extraGainDb / 20);
    source.buffer = buffer;
    source.loop = this.loop;
    source.connect(gain).connect(ctx.destination);
    source.onended = () => {
      if (this.source !== source) return;
      this.teardown();
      this.offset = 0;
      this.setState("idle");
      this.onTime(0, this.duration);
    };
    source.start(0, start);
    this.source = source;
    this.gain = gain;
    this.offset = start;
    this.startedAt = ctx.currentTime;
    this.setState("playing");
    const tick = (): void => {
      if (this.state !== "playing") return;
      this.onTime(this.time, this.duration);
      this.raf = requestAnimationFrame(tick);
    };
    this.raf = requestAnimationFrame(tick);
  }

  private teardown(): void {
    cancelAnimationFrame(this.raf);
    if (this.source) {
      this.source.onended = null;
      try {
        this.source.stop();
      } catch {
        /* already stopped */
      }
      this.source.disconnect();
      this.gain?.disconnect();
    }
    this.source = null;
    this.gain = null;
  }

  pause(): void {
    if (this.state === "loading") {
      this.token++;
      this.setState(this.offset > 0 ? "paused" : "idle");
      return;
    }
    if (this.state !== "playing") return;
    this.offset = this.time;
    this.teardown();
    this.setState("paused");
  }

  toggle(): void {
    if (this.state === "playing" || this.state === "loading") this.pause();
    else void this.play().catch(() => undefined);
  }

  seek(seconds: number): void {
    this.offset = Math.max(0, Math.min(seconds, Math.max(0, this.duration - 0.01)));
    if (this.state === "playing") void this.play(this.offset).catch(() => undefined);
    else this.onTime(this.offset, this.duration);
  }

  select(id: string): void {
    if (id === this.current) return;
    const wasPlaying = this.state === "playing" || this.state === "loading";
    const at = this.time;
    this.current = id;
    if (wasPlaying) void this.play(at).catch(() => undefined);
    else void loadBuffer(this.urls[id]!).then((b) => {
      if (this.current === id) {
        this.duration = b.duration;
        this.onTime(Math.min(this.offset, b.duration), b.duration);
      }
    }).catch(() => undefined);
  }

  setExtraGainDb(db: number): void {
    this.extraGainDb = db;
    if (this.gain) this.gain.gain.setTargetAtTime(10 ** (db / 20), getContext().currentTime, 0.02);
  }

  setLoop(loop: boolean): void {
    this.loop = loop;
    if (this.source) this.source.loop = loop;
  }
}
