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
  const status = el("p", { class: "small muted sb-status", role: "status" }, "The live engines start when you press play: a piano and two drum pieces, each in its own Worker.");
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

  const stop = (): void => {
    cancelAnimationFrame(posRaf);
    try { source?.stop(); } catch { /* not started */ }
    source?.disconnect();
    source = null;
    if (lanes) {
      const l = lanes;
      lanes = null;
      l.piano.stop(); l.kick.stop(); l.hihat.stop();
      // let the last notes ring out through the mix before the engines go
      setTimeout(() => { l.piano.dispose(); l.kick.dispose(); l.hihat.dispose(); }, 1500);
    }
    const m = mix;
    mix = null;
    setTimeout(() => m?.dispose(), 1600);
    spec?.stop();
    barEl.textContent = "Bar –";
    chordEl.textContent = "";
    delete chordEl.dataset.name;
    setPhase("idle");
  };
  registerPlayer(() => { if (phase !== "idle") stop(); });

  async function start(): Promise<void> {
    claimPlayback(stop);
    setPhase("loading");
    try {
      ctx ??= createLiveAudioContext();
      await ctx.resume();
      const c = ctx;
      status.textContent = "Loading the studio render…";
      if (!buffer) {
        const res = await fetch(asset(data.audio));
        if (!res.ok) throw new Error(`studio render: HTTP ${res.status}`);
        buffer = await c.decodeAudioData(await res.arrayBuffer());
      }
      // the two signal paths: each analysed before its gain, so both spectrograms keep running
      const preIn = c.createGain();
      const liveIn = c.createGain();
      preGain = c.createGain();
      liveGain = c.createGain();
      const aPre = c.createAnalyser(), aLive = c.createAnalyser();
      for (const a of [aPre, aLive]) { a.fftSize = 4096; a.smoothingTimeConstant = 0; a.minDecibels = -130; a.maxDecibels = 0; }
      preIn.connect(aPre); preIn.connect(preGain).connect(c.destination);
      liveIn.connect(aLive); liveIn.connect(liveGain).connect(c.destination);
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
      mix = await createLiveMix(c, liveIn);
      const m = mix;
      for (const k of Object.keys(laneHealth)) delete laneHealth[k];
      m.onStats((_f, ls) => {
        for (const [id, s] of Object.entries(ls)) {
          const h = (laneHealth[id] ??= { late: 0, gap: 0, rt: 0 });
          h.late = s.lateFrames;
          h.gap = s.gapFrames;
        }
      });
      // a fixed score needs no low latency: a second of look-ahead rides out CPU contention
      // (another tab, a busy machine) that an 85 ms buffer turned into dropouts
      const piano = new PianoLane(m, { aheadFrames: AHEAD, chunk: CHUNK, gain: data.live.piano.gain });
      const kick = new KitLane(m, "kick", { gain: data.live.kick.gain, aheadFrames: AHEAD, chunk: CHUNK, scales: [2, 2.5, 3, 3.5, 4, 5, 6], minRealtime: 1.4, fieldHz: 0 });
      // the hat also tries 4.5x: this loop on 5x rings a 558 Hz line 27 dB above the local spectrum (the
      // studio stem: 14 dB at 550 Hz); on 4.5x the line is gone (5 dB), as on 4x (see kitStage HAT_SCALES)
      const hihat = new KitLane(m, "hihat-closed", { gain: data.live.hihat.gain, aheadFrames: AHEAD, chunk: CHUNK, scales: [2.5, 3, 3.5, 4, 4.5, 5, 6], minRealtime: 1.4, fieldHz: 0 });
      const meter = (id: string) => (s: LaneWorkerStats): void => {
        const h = (laneHealth[id] ??= { late: 0, gap: 0, rt: 0 });
        h.rt = s.realtime;
        paintMeter(id);
      };
      piano.onStats = meter("piano");
      kick.onStats = meter("kick");
      hihat.onStats = meter("hihat-closed");
      const [pr, kr, hr] = await Promise.all([piano.ready, kick.ready, hihat.ready]);
      (laneHealth["piano"] ??= { late: 0, gap: 0, rt: 0 }).board = pr.boardCoordinates;
      (laneHealth["kick"] ??= { late: 0, gap: 0, rt: 0 }).scale = kr.scale;
      (laneHealth["hihat-closed"] ??= { late: 0, gap: 0, rt: 0 }).scale = hr.scale;
      lanes = { piano, kick, hihat };

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
      status.textContent = `Pre-buffering 2.0 s of live piano (${pr.boardCoordinates} board modes) and kit (${kr.scale}× kick, ${hr.scale}× hi-hat)…`;
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
          status.textContent = playingText;
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
