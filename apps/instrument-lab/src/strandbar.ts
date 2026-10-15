import { claimPlayback, registerPlayer } from "./audio";
import { DualSpectrogram, type Strip } from "./spectro";
import { Switch3D, type SwitchState } from "./switch3d";
import { asset, el } from "./util";
import { createLiveAudioContext } from "./live/liveAudio";
import { createLiveMix, KitLane, PianoLane, type LaneWorkerStats, type LiveMix } from "./live/mix";

interface KitEvent { offset: number; x: number; y: number; speed: number }
interface PianoEvent { offset: number; on?: number; off?: number; speed?: number }
interface StrandbarData {
  title: string;
  credit: string;
  bpm: number;
  beatsPerBar: number;
  bars: number;
  loopFrames: number;
  form: { name: string; chords: { slot: number; name: string }[] }[];
  audio: string;
  piano: PianoEvent[];
  kick: KitEvent[];
  hihat: KitEvent[];
  live: { piano: { gain: number }; kick: { gain: number }; hihat: { gain: number } };
  studio: { part: string; detail: string; computeS: number }[];
  liveParts: { part: string; detail: string }[];
  balance: string;
}

const SR = 48_000;
/** Look-ahead of the live lanes, frames, and how long after Play both versions start (the lanes fill first). */
const AHEAD = 96_000;
/** Most frames per audio message: the lanes run two seconds ahead, so 512 (a quarter of the messages; see public/live/lane.js). */
const CHUNK = 512;
const LEAD_S = 1.75;
/**
 * The kit grids each drum may run on, finest first; the Worker keeps the finest whose timing reaches 1.4x
 * real time (kit/worker.js). Strandbar has no restart on a coarser grid while it plays, so the lists go on
 * to 8x for slow devices: with every render 4x slower, the kick and hat on 6x ran 0.92-0.96x real time and
 * lost 2-4 s of 27, while the kick on 7x and the hat on 8x lost nothing. The hat skips 7x: the Worker times
 * the coarse grids first on every device, a hat's build is mostly its pedal settle (8x 50 ms, 7x 0.11 s),
 * and where 6x is too slow 7x was too (it timed under 1.4 at 4x, so 8x played).
 * The hat also tries 4.5x: this loop on 5x rings a 558 Hz line 27 dB above the local spectrum (the studio
 * stem: 14 dB at 550 Hz); on 4.5x the line is gone (5 dB), as on 4x (see kitStage HAT_SCALES).
 */
const KICK_SCALES = [2, 2.5, 3, 3.5, 4, 5, 6, 7, 8];
const HAT_SCALES = [2.5, 3, 3.5, 4, 4.5, 5, 6, 8];
const IDLE_STATUS = "The live engines start when you press play: a piano and two drum pieces, each in its own Worker.";
type Phase = "idle" | "loading" | "playing";

/** Strandbar on the studio render and on the live engines at once, with a switch between them. */
export async function mountStrandbar(root: HTMLElement): Promise<void> {
  let data: StrandbarData;
  try {
    const res = await fetch(asset("data/strandbar.json"));
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    data = (await res.json()) as StrandbarData;
  } catch (e) {
    root.style.minHeight = "0";
    root.replaceChildren(el("p", { class: "panel small", role: "status" }, `The groove could not be loaded (${(e as Error).message}). Reload to try again.`));
    return;
  }
  const beatFrames = data.loopFrames / (data.bars * data.beatsPerBar);

  // ------------------------------------------------------------------ DOM
  const switchHost = el("div", { class: "sb-switch", role: "switch", tabindex: "0", "aria-checked": "false", "aria-label": "Hear the version computed live in this tab" });
  const nowLabel = el("b", {}, "Studio render");
  const caption = el("div", { class: "sb-caption" }, el("span", { class: "small muted" }, "You are hearing"), nowLabel, el("span", { class: "small muted sb-hint" }, "Tap the switch or drag the lever"));
  const play = el("button", { class: "btn primary sb-play", type: "button" }, "▶ Play Strandbar");
  const barEl = el("span", { class: "sb-bar num" }, "Bar –");
  const chordEl = el("span", { class: "sb-chord" }, "");
  const canvas = el("canvas", { class: "sb-canvas", "aria-hidden": "true" });
  const stripLabel = (s: Strip, text: string): HTMLElement => el("div", { class: `sb-strip ${s}` }, el("span", {}, text), el("em", {}, "hearing"));
  const labelStudio = stripLabel("studio", "Studio render");
  const labelLive = stripLabel("live", "Live in this tab");
  const axis = el("div", { class: "sb-axis", "aria-hidden": "true" });
  for (const half of [0, 1]) for (const [f, t] of [[100, "100"], [1000, "1k"], [10000, "10k"]] as const) {
    const y = 1 - Math.log(f / 40) / Math.log(16000 / 40);
    axis.append(el("span", { style: `top: calc(${half * 50}% + ${(y * 50).toFixed(2)}% - ${half ? 1 : 6}px)` }, t));
  }
  const idleHint = el("div", { class: "sb-idle" }, "Press play: both versions scroll here, the one you hear lit");
  const spectro = el("div", { class: "sb-spectro" }, canvas, labelStudio, labelLive, axis, idleHint);
  const meters: Record<string, { root: HTMLElement; value: HTMLElement; bar: HTMLElement }> = {};
  const lanesEl = el("div", { class: "sb-lanes" });
  for (const [id, name] of [["piano", "Piano"], ["kick", "Kick"], ["hihat", "Hi-hat"]] as const) {
    const value = el("span", { class: "num" }, "—");
    const bar = el("i", {});
    const r = el("div", { class: "sb-lane" }, el("span", { class: "sb-lane-name" }, name), el("span", { class: "sb-meter" }, bar), value);
    meters[id] = { root: r, value, bar };
    lanesEl.append(r);
  }
  const status = el("p", { class: "small muted sb-status", role: "status" }, IDLE_STATUS);
  const evidence = el("details", { class: "evidence" }, el("summary", {}, "What each version computes"),
    el("div", { class: "sb-ev" },
      el("div", {}, el("h4", {}, "Studio render"), el("ul", {}, ...data.studio.map((s) => el("li", {}, el("b", {}, `${s.part}: `), `${s.detail} `, el("span", { class: "muted" }, `(${fmtDur(s.computeS)} of computing)`))))),
      el("div", {}, el("h4", {}, "Live in this tab"), el("ul", {}, ...data.liveParts.map((s) => el("li", {}, el("b", {}, `${s.part}: `), s.detail))))),
    el("p", { class: "small muted" }, data.balance),
    el("p", { class: "small muted" }, data.credit));
  root.replaceChildren(el("div", { class: "panel sb" },
    el("div", { class: "sb-grid" },
      el("div", { class: "sb-left" }, switchHost, caption),
      el("div", { class: "sb-main" },
        el("div", { class: "sb-transport" }, play, el("div", { class: "sb-pos" }, barEl, chordEl)),
        spectro, lanesEl, status)),
    evidence));
  root.style.minHeight = "0";

  // --------------------------------------------------------------- switch
  let heard: SwitchState = "studio";
  let sw: Switch3D | null = null;
  try {
    sw = new Switch3D(switchHost);
  } catch {
    switchHost.classList.add("sb-switch-flat");
    switchHost.textContent = "Studio render ⇄ Live";
  }
  const applySwitch = (s: SwitchState): void => {
    heard = s;
    switchHost.setAttribute("aria-checked", String(s === "live"));
    nowLabel.textContent = s === "studio" ? "Studio render" : "Live in this tab";
    nowLabel.className = s;
    labelStudio.classList.toggle("on", s === "studio");
    labelLive.classList.toggle("on", s === "live");
    if (spec) spec.active = s;
    crossfade();
  };
  if (sw) sw.onChange = applySwitch;
  const flip = (): void => {
    if (sw) sw.toggle();
    else applySwitch(heard === "studio" ? "live" : "studio");
  };
  switchHost.addEventListener("keydown", (e) => {
    if (e.key === " " || e.key === "Enter") { e.preventDefault(); flip(); }
  });
  if (!sw) switchHost.addEventListener("click", flip);

  // ---------------------------------------------------------------- audio
  let phase: Phase = "idle";
  let ctx: AudioContext | null = null;
  let buffer: AudioBuffer | null = null;
  let source: AudioBufferSourceNode | null = null;
  let preGain: GainNode | null = null;
  let liveGain: GainNode | null = null;
  let spec: DualSpectrogram | null = null;
  let mix: LiveMix | null = null;
  let lanes: { piano: PianoLane; kick: KitLane; hihat: KitLane } | null = null;
  let startFrame = 0;
  let posRaf = 0;
  const laneHealth: Record<string, { late: number; gap: number; rt: number; scale?: number; board?: number }> = {};

  function crossfade(): void {
    if (!ctx || !preGain || !liveGain) return;
    const t = ctx.currentTime;
    preGain.gain.setTargetAtTime(heard === "studio" ? 1 : 0, t, 0.006);
    liveGain.gain.setTargetAtTime(heard === "live" ? 1 : 0, t, 0.006);
  }

  function setPhase(p: Phase): void {
    phase = p;
    play.textContent = p === "playing" ? "■ Stop" : p === "loading" ? "Starting…" : "▶ Play Strandbar";
    play.toggleAttribute("disabled", p === "loading");
    root.classList.toggle("sb-playing", p === "playing");
    idleHint.hidden = p === "playing";
  }

  // every start is a session: stop() ends it, and a start still loading finds out at its next await
  let session = 0;
  // the nodes of the session playing now (the two paths and their analysers), disconnected after it stops
  let graph: AudioNode[] = [];
  let unStats: (() => void) | null = null;
  const stop = (): void => {
    const wasLoading = phase === "loading";
    const ended = ++session;
    cancelAnimationFrame(posRaf);
    const nodes = graph;
    graph = [];
    // Both paths fade out now. The live lanes have rendered up to two seconds ahead, so letting the mix
    // ring on played the groove's next notes for as long as it rang and then cut them off.
    if (ctx) {
      const t = ctx.currentTime;
      for (const gn of [preGain, liveGain]) if (gn) { gn.gain.cancelScheduledValues(t); gn.gain.setValueAtTime(gn.gain.value, t); gn.gain.linearRampToValueAtTime(0, t + 0.03); }
      try { source?.stop(t + 0.03); } catch { /* not started */ }
    }
    if (source) nodes.push(source);
    source = null;
    preGain = liveGain = null; // the switch no longer reaches a stopped session
    unStats?.();
    unStats = null;
    const l = lanes;
    lanes = null;
    const m = mix;
    mix = null;
    if (l) for (const lane of [l.piano, l.kick, l.hihat]) { lane.onStats = null; lane.stop(); }
    // the engines and the nodes go once the fade is over (at once while loading: nothing sounds yet)
    const ringMs = wasLoading ? 0 : 150;
    setTimeout(() => {
      if (l) { l.piano.dispose(); l.kick.dispose(); l.hihat.dispose(); }
      m?.dispose();
      for (const n of nodes) n.disconnect();
      // nothing plays: the audio device rests until the next start resumes it (unless that start came first)
      if (session === ended && ctx?.state === "running") void ctx.suspend();
    }, ringMs);
    spec?.stop();
    barEl.textContent = "Bar –";
    chordEl.textContent = "";
    delete chordEl.dataset.name;
    status.textContent = IDLE_STATUS; // not "playing in step" once nothing plays
    setPhase("idle");
  };
  const onClaimed = (): void => { if (phase !== "idle") stop(); };
  registerPlayer(onClaimed);

  async function start(): Promise<void> {
    claimPlayback(onClaimed);
    const my = ++session;
    setPhase("loading");
    try {
      ctx ??= createLiveAudioContext();
      await ctx.resume();
      const c = ctx;
      if (my !== session) return;
      status.textContent = "Loading the studio render…";
      if (!buffer) {
        const res = await fetch(asset(data.audio));
        if (!res.ok) throw new Error(`studio render: HTTP ${res.status}`);
        buffer = await c.decodeAudioData(await res.arrayBuffer());
      }
      if (my !== session) return;
      // the two signal paths: each analysed before its gain, so both spectrograms keep running
      const preIn = c.createGain();
      const liveIn = c.createGain();
      preGain = c.createGain();
      liveGain = c.createGain();
      const aPre = c.createAnalyser(), aLive = c.createAnalyser();
      for (const a of [aPre, aLive]) { a.fftSize = 4096; a.smoothingTimeConstant = 0; a.minDecibels = -130; a.maxDecibels = 0; }
      preIn.connect(aPre); preIn.connect(preGain).connect(c.destination);
      liveIn.connect(aLive); liveIn.connect(liveGain).connect(c.destination);
      graph = [preIn, liveIn, preGain, liveGain, aPre, aLive];
      preGain.gain.value = heard === "studio" ? 1 : 0;
      liveGain.gain.value = heard === "live" ? 1 : 0;
      spec?.dispose();
      spec = new DualSpectrogram(canvas, { studio: aPre, live: aLive }, {
        pxPerSecond: window.matchMedia("(max-width: 720px)").matches ? 60 : 96,
        beatsPerBar: data.beatsPerBar,
        beatAt: (t) => (t * SR >= startFrame ? (t * SR - startFrame) / beatFrames : null),
      });
      spec.active = heard;

      status.textContent = "Starting the live engines and measuring this device…";
      const m = await createLiveMix(c, liveIn);
      if (my !== session) { m.dispose(); return; }
      mix = m;
      for (const k of Object.keys(laneHealth)) delete laneHealth[k];
      // the meters start blank, not with the last session's speeds and dropouts
      for (const mtr of Object.values(meters)) { mtr.value.textContent = "—"; mtr.bar.style.width = "0%"; mtr.root.classList.remove("warn", "bad"); }
      unStats = m.onStats((_f, ls) => {
        for (const [id, s] of Object.entries(ls)) {
          const h = (laneHealth[id] ??= { late: 0, gap: 0, rt: 0 });
          h.late = s.lateFrames;
          h.gap = s.gapFrames;
        }
      });
      // a fixed score needs no low latency: two seconds of look-ahead ride out CPU contention
      // (another tab, a busy machine) that an 85 ms buffer turned into dropouts
      const piano = new PianoLane(m, { aheadFrames: AHEAD, chunk: CHUNK, gain: data.live.piano.gain });
      const kick = new KitLane(m, "kick", { gain: data.live.kick.gain, aheadFrames: AHEAD, chunk: CHUNK, scales: KICK_SCALES, minRealtime: 1.4, fieldHz: 0 });
      const hihat = new KitLane(m, "hihat-closed", { gain: data.live.hihat.gain, aheadFrames: AHEAD, chunk: CHUNK, scales: HAT_SCALES, minRealtime: 1.4, fieldHz: 0 });
      // held from here, so a stop while they start disposes them
      lanes = { piano, kick, hihat };
      // An engine that fails once it is ready (while the others still start, or mid-loop) has left the mix
      // (mix.ts LaneWorker): the status says so from then on rather than play on without it.
      let failure: string | null = null;
      const failed = (part: string) => (e: Error): void => {
        if (my !== session) return;
        failure = `The live ${part} stopped (${e.message}). Stop and play again to restart it.`;
        status.textContent = failure;
      };
      piano.onFailure = failed("piano");
      kick.onFailure = failed("kick");
      hihat.onFailure = failed("hi-hat");
      const meter = (id: string) => (s: LaneWorkerStats): void => {
        const h = (laneHealth[id] ??= { late: 0, gap: 0, rt: 0 });
        h.rt = s.realtime;
        paintMeter(id);
      };
      piano.onStats = meter("piano");
      kick.onStats = meter("kick");
      hihat.onStats = meter("hihat-closed");
      const [pr, kr, hr] = await Promise.all([piano.ready, kick.ready, hihat.ready]);
      if (my !== session) return;
      (laneHealth["piano"] ??= { late: 0, gap: 0, rt: 0 }).board = pr.boardCoordinates;
      (laneHealth["kick"] ??= { late: 0, gap: 0, rt: 0 }).scale = kr.scale;
      (laneHealth["hihat-closed"] ??= { late: 0, gap: 0, rt: 0 }).scale = hr.scale;

      // one clock for everything: the loop and all three engines start on the same frame
      startFrame = Math.round((c.currentTime + LEAD_S) * SR);
      source = c.createBufferSource();
      source.buffer = buffer;
      source.loop = true;
      source.connect(preIn);
      source.start(startFrame / SR);
      piano.sequence(startFrame, data.loopFrames, data.piano);
      kick.sequence(startFrame, data.loopFrames, data.kick);
      hihat.sequence(startFrame, data.loopFrames, data.hihat);
      spec.start();
      setPhase("playing");
      const playingText = `Both versions are playing in step. The drums run on a grid ${kr.scale}× (kick) and ${hr.scale}× (hi-hat) coarser than the studio render, chosen for this device.`;
      status.textContent = failure ?? `Pre-buffering 2.0 s of live piano (a ${pr.boardCoordinates}-coordinate soundboard) and kit (${kr.scale}× kick, ${hr.scale}× hi-hat)…`;
      let announcedPlaying = false;
      const pos = (): void => {
        posRaf = requestAnimationFrame(pos);
        const f = c.currentTime * SR - startFrame;
        if (f < 0) {
          const pct = Math.max(0, Math.min(99, Math.round((1 + f / (LEAD_S * SR)) * 100)));
          barEl.textContent = `Buffering ${pct}%`;
          return;
        }
        if (!announcedPlaying) {
          announcedPlaying = true;
          status.textContent = failure ?? playingText;
        }
        const inLoop = f % data.loopFrames;
        const beat = inLoop / beatFrames;
        const bar = Math.floor(beat / data.beatsPerBar);
        const slot = Math.floor((beat % data.beatsPerBar) * 4);
        const form = data.form[bar];
        if (!form) return;
        barEl.textContent = `Bar ${bar + 1} of ${data.bars}`;
        const chord = [...form.chords].reverse().find((ch) => ch.slot <= slot);
        const name = chord ? chord.name : "";
        if (chordEl.dataset.name !== name) {
          chordEl.dataset.name = name;
          chordEl.replaceChildren(...name.split(/([♭♯])/).filter(Boolean).map((t) => (t === "♭" || t === "♯" ? el("span", { class: "acc" }, t) : document.createTextNode(t))));
        }
      };
      pos();
    } catch (e) {
      if (my !== session) return; // stopped while loading: its lanes were disposed on purpose
      stop();
      status.textContent = `Could not start: ${(e as Error).message}`;
    }
  }

  function paintMeter(id: string): void {
    const key = id === "hihat-closed" ? "hihat" : id;
    const mtr = meters[key];
    const h = laneHealth[id];
    if (!mtr || !h) return;
    const rt = h.rt;
    // the mixer counts late frames (arrived after they played) and gap frames (skipped by a lane that
    // fell behind; a drum resting between strokes is not a gap), each frame once: together, the audio lost
    const lost = (h.late + h.gap) / 48;
    const troubled = h.late + h.gap > 0;
    mtr.bar.style.width = `${Math.min(100, (rt / 6) * 100).toFixed(1)}%`;
    mtr.root.classList.toggle("warn", rt > 0 && rt < 1.15);
    mtr.root.classList.toggle("bad", troubled);
    const grid = h.scale ? ` · grid ${h.scale}×` : h.board ? ` · board ${h.board}` : "";
    mtr.value.textContent = rt > 0 ? `${rt.toFixed(1)}× real time${grid}${troubled ? ` · ${lost < 10 ? lost.toFixed(1) + " ms" : lost < 1000 ? Math.round(lost) + " ms" : (lost / 1000).toFixed(1) + " s"} dropped` : ""}` : "—";
    mtr.root.dataset.late = String(h.late);
    mtr.root.dataset.gap = String(h.gap);
  }

  play.addEventListener("click", () => {
    if (phase === "playing") stop();
    else if (phase === "idle") void start();
  });
  applySwitch("studio");
}

function fmtDur(s: number): string {
  if (s >= 3600) return `${(s / 3600).toFixed(1)} h`;
  if (s >= 90) return `${Math.round(s / 60)} min`;
  return `${Math.round(s)} s`;
}
