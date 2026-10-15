// The player piano and drum machine: any MIDI file, punched into a roll that drives FrankenSim's
// live grand piano and drum kit. One AudioContext clock; the piano and every drum piece the file
// uses are computed in their own Workers a second ahead of what you hear (a fixed score needs no
// low latency, and a second of look-ahead rides out a busy machine), and every picture on the page
// is the engines' own state for the moment being heard.

import "./player.css";
import { claimPlayback, registerPlayer } from "../audio";
import { asset, el, noteName } from "../util";
import { createLiveAudioContext } from "../live/liveAudio";
import { createLiveMix, KitLane, PianoLane, type KitReady, type KitScoreEvent, type LaneWorkerStats, type LiveMix, type PianoReady } from "../live/mix";
import { parsePianoSnapshot } from "../live/pianoClient";
import { decodeMidi, MAX_MIDI_BYTES } from "../midi/smf";
import { arrange, gmInstrument, PIECE_NAMES, type Arrangement, type PieceId } from "../midi/arrange";
import { RollView } from "./roll";
import { SoundboardView } from "./soundboard3d";
import { KitView, fmtJ } from "./kit3d";

const SR = 48_000;
/** Look-ahead of every lane [frames], and how long after Play the score starts (cold, warm). */
const AHEAD = 48_000;
/**
 * Most frames per audio message from a lane. The score is computed a second ahead, so a lane can send
 * 512 frames at a time: a quarter of the messages into the mixer's realtime thread, every engine call
 * unchanged (public/live/lane.js). The live kit, which plays what is struck now, keeps 128.
 */
const CHUNK = 512;
const LEAD_COLD_S = 1.5;
const LEAD_WARM_S = 0.85;
/** One fixed playback gain per instrument (see the evidence note). */
const PIANO_GAIN = 0.8;
/** Standard stick speeds [m/s] and the live kit's per-piece gains (src/live/kitStage.ts). */
const STICK: Record<PieceId, number> = { kick: 3.0, snare: 4.5, "snare-off": 4.5, "tom-12": 3.5, "tom-14": 3.5, "hihat-closed": 2.5, "hihat-open": 2.5, crash: 6.0, ride: 2.5 };
// the live kit section's gains (src/live/kitStage.ts), re-derived for the retuned kit (fs-vkplate 96ccf8f1b)
const KIT_GAIN: Record<PieceId, number> = { kick: 0.0319, snare: 0.0129, "snare-off": 0.0115, "tom-12": 0.0276, "tom-14": 0.0267, "hihat-closed": 0.134, "hihat-open": 0.134, crash: 0.0496, ride: 0.182 };
const IDLE_DISPOSE_MS = 90_000;
/**
 * Each piece's grid choice: the grid coarsenings it may use (finest first) and the speed it must
 * reach alone when it measures itself. A coarser grid loses the high end first, and pieces differ
 * in how much they have to lose. Measured on single strokes: a closed hi-hat keeps -18 dB of its
 * energy above 4 kHz at 4x and -91 dB at 6x (a bare click with no shimmer), while a kick has no
 * energy above 4 kHz even at 4x. So cymbals hold on to fine grids and the kick gives way first;
 * a piece that still falls behind while playing is restarted one grid coarser.
 * The crash and ride try 3x and 3.5x as well: their octave balance against the 1x grid improves at
 * every step down from 5x (crash 18.2 / 14.6 / 11.2 / 10.2 dB at 5x / 4x / 3.5x / 3x; ride 19.6 /
 * 12.4 / 10.9 / 7.2), and 5x stays as the fallback for slower devices.
 */
const GRID: Record<PieceId, { scales: number[]; minRt: number }> = {
  kick: { scales: [5, 6, 7, 8], minRt: 1.7 },
  snare: { scales: [4, 5, 6, 7, 8], minRt: 2.2 },
  "snare-off": { scales: [4, 5, 6, 7, 8], minRt: 2.0 },
  "tom-12": { scales: [4, 5, 6, 7], minRt: 1.8 },
  "tom-14": { scales: [4, 5, 6, 7], minRt: 1.8 },
  "hihat-closed": { scales: [4, 5, 6], minRt: 1.5 },
  "hihat-open": { scales: [4, 5, 6], minRt: 1.5 },
  crash: { scales: [3, 3.5, 4, 5, 6], minRt: 1.5 },
  ride: { scales: [3, 3.5, 4, 5, 6], minRt: 1.5 },
};

interface Standard {
  id: string; title: string; kicker: string; note: string; provenance: string; groove: string; grooveName: string;
  tempoBpm: number; bars: number; choruses: number; durationS: number; file: string; pianoNotes: number; drumHits: number;
}
interface Catalog { source: { repo: string; commit: string }; entries: Standard[] }
type Phase = "empty" | "loading" | "ready" | "warming" | "playing" | "paused";
/** late: audio that reached the mixer after its moment; skipped: audio a lane that fell behind never rendered; both are real gaps */
interface Health { rt: number; late: number; gap: number; base: number; skipBase: number; detail: string; slow: number; scale?: number; stats?: LaneWorkerStats; lastLate?: number; lateAt?: number; lostSeen?: number }

const ICON = {
  play: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M4.5 2.6v10.8l8.6-5.4z" fill="currentColor"/></svg>',
  pause: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M4 2.5h3v11H4zM9 2.5h3v11H9z" fill="currentColor"/></svg>',
  rewind: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M3 2.5h1.8v11H3zM13.4 2.6v10.8L5.6 8z" fill="currentColor"/></svg>',
  minus: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M3 7.25h10v1.5H3z" fill="currentColor"/></svg>',
  plus: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M7.25 3h1.5v4.25H13v1.5H8.75V13h-1.5V8.75H3v-1.5h4.25z" fill="currentColor"/></svg>',
  file: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 3v12m0-12-4.5 4.5M12 3l4.5 4.5M4 15v4.5A1.5 1.5 0 0 0 5.5 21h13a1.5 1.5 0 0 0 1.5-1.5V15" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  left: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M10 3 5 8l5 5" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  right: '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="m6 3 5 5-5 5" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>',
};
const svg = (html: string): HTMLElement => { const s = document.createElement("span"); s.className = "pp-ico"; s.innerHTML = html; return s; };

export async function mountPlayer(root: HTMLElement): Promise<void> {
  // ----------------------------------------------------------------- DOM
  const fileInput = el("input", { type: "file", accept: ".mid,.midi,.kar,.rmi,audio/midi,audio/x-midi", class: "sr-only", tabindex: "-1", "aria-hidden": "true" }) as HTMLInputElement;
  const fileTile = el("button", { type: "button", class: "pp-tile pp-tile-file", "aria-describedby": "pp-file-hint" },
    svg(ICON.file), el("b", {}, "Your MIDI file"), el("span", { id: "pp-file-hint" }, "Choose a .mid, or drop one anywhere on the player"));
  const shelf = el("div", { class: "pp-shelf", role: "group", "aria-label": "Rolls" }, fileTile);
  const shelfPrev = el("button", { type: "button", class: "pp-shelf-nav", "aria-label": "Scroll the rolls left" }, svg(ICON.left));
  const shelfNext = el("button", { type: "button", class: "pp-shelf-nav", "aria-label": "Scroll the rolls right" }, svg(ICON.right));
  const shelfCount = el("span", { class: "small muted" }, "");
  const shelfHead = el("div", { class: "pp-shelf-head" }, el("div", {}, el("h3", {}, "Rolls"), shelfCount), el("div", { class: "pp-shelf-navs" }, shelfPrev, shelfNext));

  const rollHost = el("div", { class: "pp-roll-host" });
  const nameplate = el("div", { class: "pp-nameplate", "aria-hidden": "true" }, el("span", {}, "FrankenSim"), el("i", {}), el("span", {}, "Player · Live"));
  const dropVeil = el("div", { class: "pp-drop", hidden: "" }, svg(ICON.file), el("b", {}, "Drop to punch a new roll"), el("span", {}, "Standard MIDI files, format 0 or 1, up to 4 MiB"));
  const face = el("div", { class: "pp-face" }, nameplate, rollHost);

  const eyebrow = el("p", { class: "pp-eyebrow" }, "On the player");
  const titleEl = el("h3", { class: "pp-title" }, "No roll loaded");
  const subEl = el("p", { class: "pp-sub small muted" }, "Pick a roll from the shelf above, or drop in a MIDI file of your own.");
  const chordEl = el("span", { class: "pp-chord", "aria-live": "off" }, "");
  const nextEl = el("span", { class: "pp-next small muted" }, "");
  const clockEl = el("span", { class: "pp-clock num" }, "0:00 / 0:00");
  const barEl = el("span", { class: "pp-barno num small muted" }, "");
  const seek = el("input", { type: "range", class: "live pp-seek", min: "0", max: "1000", step: "1", value: "0", "aria-label": "Position in the roll", disabled: "" }) as HTMLInputElement;
  const playBtn = el("button", { type: "button", class: "btn live pp-play", disabled: "" }, svg(ICON.play), el("span", {}, "Play"));
  const rewindBtn = el("button", { type: "button", class: "btn pp-round", "aria-label": "Back to the start", disabled: "" }, svg(ICON.rewind));
  const tempoDown = el("button", { type: "button", class: "btn pp-round sm", "aria-label": "Slower" }, svg(ICON.minus));
  const tempoUp = el("button", { type: "button", class: "btn pp-round sm", "aria-label": "Faster" }, svg(ICON.plus));
  const tempoVal = el("output", { class: "pp-tempo-val num" }, "100%");
  const tempoBpm = el("span", { class: "pp-tempo-bpm num small muted" }, "");
  const tempo = el("div", { class: "pp-tempo", role: "group", "aria-label": "Tempo" }, el("span", { class: "small muted" }, "Tempo"), tempoDown, el("div", { class: "pp-tempo-read" }, tempoVal, tempoBpm), tempoUp);
  const status = el("p", { class: "pp-status small", role: "status" }, "");
  const lanesEl = el("div", { class: "pp-lanes", "aria-label": "Engines" });
  const keysHint = el("p", { class: "pp-keys small muted" }, el("kbd", {}, "Space"), " play or pause · ", el("kbd", {}, "←"), el("kbd", {}, "→"), " 5 s · ", el("kbd", {}, "["), el("kbd", {}, "]"), " tempo · drag the paper to scrub");
  const side = el("div", { class: "pp-side" },
    el("div", { class: "pp-now" }, eyebrow, titleEl, subEl,
      el("div", { class: "pp-chord-row" }, chordEl, nextEl),
      el("div", { class: "pp-pos" }, clockEl, barEl),
      seek,
      el("div", { class: "pp-transport" }, rewindBtn, playBtn, tempo)),
    status, lanesEl, keysHint);

  // physics: soundboard and kit
  const boardStage = el("div", { class: "pp-stage" });
  const boardHud = el("div", { class: "pp-hud" });
  const stringsCanvas = el("canvas", { class: "pp-strings", "aria-hidden": "true", hidden: "" });
  const boardView = el("article", { class: "pp-view", "aria-label": "The live soundboard" },
    el("header", { class: "pp-view-head" }, el("h4", {}, "Soundboard"), el("span", { class: "pp-view-tag" }, "Model D · 32 coordinates · heights ×", el("b", { class: "pp-board-scale" }, "—"))),
    boardStage, boardHud, stringsCanvas);
  boardStage.append(el("div", { class: "pp-stage-idle" }, "Its motion appears here as the roll plays"));
  const kitStage = el("div", { class: "pp-stage" });
  const kitHud = el("div", { class: "pp-hud" });
  const kitView = el("article", { class: "pp-view", "aria-label": "The live drum kit" },
    el("header", { class: "pp-view-head" }, el("h4", {}, "Drum kit"), el("span", { class: "pp-view-tag" }, "surface velocity, as you hear it")),
    kitStage, kitHud);
  const tabBoard = el("button", { type: "button", "aria-pressed": "true" }, "Soundboard");
  const tabKit = el("button", { type: "button", "aria-pressed": "false" }, "Drum kit");
  const tabs = el("div", { class: "pp-tabs seg live", role: "group", "aria-label": "Show the soundboard or the drum kit" }, tabBoard, tabKit);
  const physics = el("div", { class: "pp-physics", "data-view": "board" }, tabs, boardView, kitView);

  const report = el("div", { class: "pp-report", hidden: "" });
  const evidence = el("details", { class: "evidence pp-evidence" }, el("summary", {}, "How a MIDI file becomes a performance"), evidenceBody());

  const panel = el("div", { class: "panel pp" },
    el("div", { class: "pp-shelf-wrap" }, shelfHead, shelf),
    el("div", { class: "pp-console" }, face, side),
    physics, report, evidence, dropVeil, fileInput);
  const mini = el("div", { class: "pp-mini", hidden: "" });
  const miniPlay = el("button", { type: "button", class: "pp-mini-play", "aria-label": "Pause" }, svg(ICON.pause));
  const miniTitle = el("b", {}, "");
  const miniChord = el("span", { class: "pp-mini-chord" }, "");
  const miniBar = el("i", {});
  mini.append(miniPlay, el("div", { class: "pp-mini-text" }, miniTitle, miniChord), el("div", { class: "pp-mini-progress" }, miniBar));
  root.replaceChildren(panel, mini);
  root.style.minHeight = "0";

  const roll = new RollView(rollHost);
  const board = new SoundboardView(boardStage);
  const kit = new KitView(kitStage);
  if (import.meta.env.DEV) Object.assign(window, { __pp: { board, kit, roll } });

  // ---------------------------------------------------------------- state
  let phase: Phase = "empty";
  let arr: Arrangement | null = null;
  let meta: { kind: "standard" | "file"; id: string; title: string; sub: string; file?: string } | null = null;
  let position = 0;
  let tempoFactor = 1;
  let ctx: AudioContext | null = null;
  let mix: LiveMix | null = null;
  let master: GainNode | null = null;
  let piano: PianoLane | null = null;
  let pianoReady: PianoReady | null = null;
  const kits = new Map<PieceId, { lane: KitLane; ready: KitReady | null }>();
  let epoch = 0;
  let startFrame = 0;
  let startPos = 0;
  let idleTimer = 0;
  let token = 0;
  const snaps: { frame: number; data: Float32Array }[] = [];
  // late: audio that reached the mixer after its moment (a real dropout); base: the count when the score last started
  const health: Record<string, Health> = {};
  let ledgerPrev: { frame: number; input: number; string: number; board: number; damper: number; hammer: number } | null = null;
  let hitCursor = 0;
  let selectedTile: HTMLElement | null = null;

  // --------------------------------------------------------------- library
  // two packs: the jazz editor's standards, and video-game MIDI files as found online
  const packs: { id: string; label: string; blurb: string; entries: Standard[] }[] = [];
  const load = async (url: string): Promise<Catalog | null> => {
    try {
      const res = await fetch(asset(url));
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      return (await res.json()) as Catalog;
    } catch { return null; }
  };
  const [standards, game] = await Promise.all([load("midi/standards.json"), load("midi/game-pack.json")]);
  if (standards) packs.push({ id: "standards", label: "Standards", blurb: "chord charts from the jazz chord editor's library, comped by its own bass and piano generators, with drums", entries: standards.entries });
  if (game) packs.push({ id: "game", label: "Game pack", blurb: (game as unknown as { note: string }).note, entries: game.entries });
  const tiles = new Map<string, HTMLElement>();
  const makeTile = (s: Standard): HTMLElement => {
    const tile = el("button", { type: "button", class: `pp-tile g-${s.groove.split("@")[0]}`, "aria-label": `${s.title}, ${s.kicker}: ${s.grooveName}, ${Math.round(s.tempoBpm)} beats per minute, ${fmtTime(s.durationS)}` },
      el("span", { class: "pp-tile-band" }, s.grooveName),
      el("b", { class: "pp-tile-title" }, s.title),
      el("span", { class: "pp-tile-kicker" }, s.kicker),
      el("span", { class: "pp-tile-meta num" }, `♩ ${Math.round(s.tempoBpm)} · ${fmtTime(s.durationS)}`),
      el("span", { class: "pp-tile-spool", "aria-hidden": "true" }));
    tile.title = s.note;
    tile.addEventListener("click", () => {
      primeAudio();
      void loadStandard(s, tile, true);
    });
    return tile;
  };
  const packButtons = el("div", { class: "seg live pp-packs", role: "group", "aria-label": "Packs of rolls" });
  const showPack = (id: string): void => {
    const pack = packs.find((p) => p.id === id) ?? packs[0];
    if (!pack) return;
    for (const b of packButtons.querySelectorAll("button")) b.setAttribute("aria-pressed", String(b.dataset.pack === pack.id));
    shelfCount.textContent = `${pack.entries.length} rolls · ${pack.blurb}`;
    for (const t of [...shelf.children]) if (t !== fileTile) t.remove();
    for (const s of pack.entries) {
      const tile = tiles.get(s.id) ?? makeTile(s);
      tiles.set(s.id, tile);
      shelf.append(tile);
    }
    shelf.scrollLeft = 0;
    shelfEdges();
  };
  for (const p of packs) {
    const b = el("button", { type: "button", "aria-pressed": "false", "data-pack": p.id }, p.label, el("small", {}, `${p.entries.length} rolls`));
    b.addEventListener("click", () => showPack(p.id));
    packButtons.append(b);
  }
  if (packs.length > 1) shelfHead.firstElementChild!.append(packButtons);
  if (!packs.length) shelfCount.textContent = "The rolls could not be loaded; your own files still play.";
  const scrollShelf = (dir: number): void => shelf.scrollBy({ left: dir * Math.max(240, shelf.clientWidth * 0.8), behavior: "smooth" });
  shelfPrev.addEventListener("click", () => scrollShelf(-1));
  shelfNext.addEventListener("click", () => scrollShelf(1));
  const shelfEdges = (): void => {
    shelfPrev.toggleAttribute("disabled", shelf.scrollLeft < 4);
    shelfNext.toggleAttribute("disabled", shelf.scrollLeft + shelf.clientWidth > shelf.scrollWidth - 4);
  };
  shelf.addEventListener("scroll", shelfEdges, { passive: true });
  new ResizeObserver(shelfEdges).observe(shelf);

  fileTile.addEventListener("click", () => { primeAudio(); fileInput.click(); });
  fileInput.addEventListener("change", () => {
    const f = fileInput.files?.[0];
    fileInput.value = "";
    if (f) void loadFile(f, true);
  });
  let dragDepth = 0;
  panel.addEventListener("dragenter", (e) => { if (hasFiles(e)) { dragDepth++; dropVeil.hidden = false; e.preventDefault(); } });
  panel.addEventListener("dragover", (e) => { if (hasFiles(e)) { e.preventDefault(); e.dataTransfer!.dropEffect = "copy"; } });
  panel.addEventListener("dragleave", () => { dragDepth = Math.max(0, dragDepth - 1); if (!dragDepth) dropVeil.hidden = true; });
  panel.addEventListener("drop", (e) => {
    if (!hasFiles(e)) return;
    e.preventDefault();
    dragDepth = 0;
    dropVeil.hidden = true;
    const f = [...(e.dataTransfer?.files ?? [])].find((x) => /\.(mid|midi|kar|rmi)$/i.test(x.name)) ?? e.dataTransfer?.files[0];
    if (f) void loadFile(f, ctx?.state === "running");
  });

  function primeAudio(): void {
    // inside the click: browsers start audio only from a gesture
    try {
      ctx ??= createLiveAudioContext();
      void ctx.resume();
    } catch { /* reported when Play fails */ }
  }

  // ----------------------------------------------------------------- loading
  async function loadStandard(s: Standard, tile: HTMLElement, autoplay: boolean): Promise<void> {
    const my = ++token;
    setPhase("loading");
    say(`Fetching the ${s.title} roll…`);
    try {
      const res = await fetch(asset(s.file));
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const bytes = new Uint8Array(await res.arrayBuffer());
      if (my !== token) return;
      meta = { kind: "standard", id: s.id, title: s.title, sub: `${s.kicker} · ${s.grooveName} · ♩ ${Math.round(s.tempoBpm)}${s.groove === "game" ? "" : ` · ${s.choruses > 1 ? `${s.choruses} choruses` : "one chorus"}`}` };
      const ok = punch(bytes, `${s.id}.mid`);
      if (!ok || my !== token) return;
      select(tile);
      if (autoplay) await play(0);
    } catch (e) {
      if (my !== token) return;
      setPhase(arr ? "ready" : "empty");
      say(`The roll could not be loaded (${(e as Error).message}). Try again, or pick another.`, "bad");
    }
  }

  async function loadFile(f: File, autoplay: boolean): Promise<void> {
    const my = ++token;
    if (f.size > MAX_MIDI_BYTES) { say(`${f.name} is ${(f.size / 1048576).toFixed(1)} MiB; files of up to 4 MiB are accepted.`, "bad"); return; }
    setPhase("loading");
    say(`Reading ${f.name}…`);
    try {
      const bytes = new Uint8Array(await f.arrayBuffer());
      if (my !== token) return;
      meta = { kind: "file", id: `file:${f.name}`, title: f.name.replace(/\.(mid|midi|kar|rmi)$/i, ""), sub: f.name, file: f.name };
      const ok = punch(bytes, f.name);
      if (!ok || my !== token) return;
      fileTile.querySelector("b")!.textContent = f.name.length > 28 ? `${f.name.slice(0, 26)}…` : f.name;
      fileTile.querySelector("span:not(.pp-ico)")!.textContent = "Your file · choose another or drop one";
      select(fileTile);
      if (autoplay) await play(0);
      else say("Your roll is punched. Press Play: the engines start in this tab.");
    } catch (e) {
      if (my !== token) return;
      setPhase(arr ? "ready" : "empty");
      say(`${f.name} could not be read (${(e as Error).message}).`, "bad");
    }
  }

  /** Decode, salvage and arrange: the roll is punched. */
  function punch(bytes: Uint8Array, name: string): boolean {
    if (phase === "playing") pause();
    const decoded = decodeMidi(bytes);
    if (!decoded.ok) {
      setPhase(arr ? "ready" : "empty");
      say(decoded.message, "bad");
      return false;
    }
    const a = arrange(decoded.midi, name);
    if (!a.notes.length && !a.hits.length) {
      setPhase(arr ? "ready" : "empty");
      say("This file holds no notes the piano or the kit could play.", "bad");
      return false;
    }
    arr = a;
    if (meta?.kind === "file" && a.title && a.title !== meta.title) meta.sub = `${a.title !== meta.title ? `“${a.title}” · ` : ""}${meta.file}`;
    position = 0;
    tempoFactor = 1;
    showTempo();
    roll.setArrangement(a);
    kit.setUsed(a.pieces);
    kit.clear();
    titleEl.textContent = meta?.title ?? a.title;
    eyebrow.textContent = meta?.kind === "file" ? "Your file, on the player" : "On the player";
    subEl.textContent = meta?.sub ?? "";
    seek.disabled = false;
    renderReport(a);
    paintLanes();
    setPhase("ready");
    say(readyLine(a));
    paintPosition(0);
    requestDraw();
    return true;
  }

  function select(tile: HTMLElement): void {
    selectedTile?.classList.remove("on");
    selectedTile?.removeAttribute("aria-current");
    tile.classList.add("on");
    tile.setAttribute("aria-current", "true");
    selectedTile = tile;
    const r = tile.getBoundingClientRect(), s = shelf.getBoundingClientRect();
    if (r.left < s.left || r.right > s.right) shelf.scrollBy({ left: r.left - s.left - 24, behavior: "smooth" });
  }

  // ---------------------------------------------------------------- engines
  async function ensureEngines(a: Arrangement): Promise<boolean> {
    ctx ??= createLiveAudioContext();
    await ctx.resume();
    const c = ctx;
    if (!mix) {
      master = c.createGain();
      master.gain.value = 0;
      master.connect(c.destination);
      mix = await createLiveMix(c, master);
      mix.onStats((_f, lanes) => {
        for (const [id, s] of Object.entries(lanes)) {
          const h = health[id];
          if (!h) continue;
          h.late = s.lateFrames;
          h.gap = s.gapFrames;
        }
      });
    }
    const m = mix;
    let warm = !!piano;
    const started: Promise<unknown>[] = [];
    if (!piano) {
      piano = new PianoLane(m, { aheadFrames: AHEAD, chunk: CHUNK, gain: PIANO_GAIN, snapshotHz: 30, minRealtime: 1.3 });
      piano.onStats = (s) => { const h = health["piano"]; if (h) { h.rt = s.realtime; h.stats = s; } };
      piano.onSnapshot = (frame, data, ep) => { if (ep >= epoch) { snaps.push({ frame, data }); if (snaps.length > 90) snaps.shift(); } };
      health["piano"] = { rt: 0, late: 0, gap: 0, base: 0, skipBase: 0, detail: "", slow: 0 };
      started.push(piano.ready.then((r) => {
        pianoReady = r;
        health["piano"]!.detail = `${r.boardCoordinates} board coordinates`;
        if (r.board) board.setBoard(r.board.xy, r.board.shapes);
        boardStage.querySelector(".pp-stage-idle")?.remove();
      }));
    }
    for (const [p, k] of kits) if (!a.pieces.includes(p)) { k.lane.dispose(); kits.delete(p); delete health[p]; }
    for (const p of a.pieces) {
      if (kits.has(p)) continue;
      warm = false;
      started.push(startPiece(m, p, GRID[p].scales));
    }
    if (started.length) {
      const total = started.length;
      let done = 0;
      const n = a.pieces.length;
      say(`Starting ${!pianoReady ? "the piano" : ""}${!pianoReady && n ? " and " : ""}${n ? `${n} drum piece${n > 1 ? "s" : ""}` : ""}, each in its own Worker, and measuring this device…`);
      started.forEach((p) => void p.then(() => { done++; say(`Starting the engines… ${done} of ${total} ready`); }));
      await Promise.all(started);
    }
    paintLanes();
    return warm;
  }

  /**
   * One drum piece in its Worker. Each piece measures itself and keeps the finest grid that runs
   * GRID[p].minRt times faster than real time alone (see GRID). A piece that still falls behind while playing is
   * restarted one grid coarser (see the stats handler), and the page says so.
   */
  function startPiece(m: LiveMix, p: PieceId, scales: number[]): Promise<void> {
    const lane = new KitLane(m, p, { gain: KIT_GAIN[p], aheadFrames: AHEAD, chunk: CHUNK, minRealtime: GRID[p].minRt, scales, fieldFrames: SR / 40, bursts: false });
    const entry = { lane, ready: null as KitReady | null };
    kits.set(p, entry);
    const h: Health = { rt: 0, late: 0, gap: 0, base: 0, skipBase: 0, detail: "", slow: 0 };
    health[p] = h;
    lane.onStats = (s) => {
      if (kits.get(p) !== entry) return;
      if (!s.asleep) h.rt = s.realtime;
      h.stats = s;
      h.slow = phase === "playing" && !s.asleep && s.realtime > 0 && s.realtime < 1.05 ? h.slow + 1 : 0;
      // audio actually lost (more than 50 ms since the last report) steps the piece down at once:
      // a gap is worse than a duller cymbal
      const lost = lostFrames(h);
      const fresh = lost - (h.lostSeen ?? lost);
      h.lostSeen = lost;
      if ((h.slow >= 3 || (phase === "playing" && fresh > 2400)) && (h.scale ?? 8) < 10) void coarsen(p, h.scale ?? 8);
    };
    lane.onField = (f) => { if ((f.epoch ?? 0) >= epoch && kits.get(p) === entry) kit.field(p, f); };
    return lane.ready.then((r) => {
      entry.ready = r;
      h.scale = r.scale;
      h.detail = `grid ${r.scale}× coarser`;
      kit.setReady(p, r);
    });
  }

  async function coarsen(p: PieceId, scale: number): Promise<void> {
    const old = kits.get(p);
    if (!old || !mix) return;
    const lostBefore = health[p] ? lostFrames(health[p]!) : 0;
    old.lane.dispose(); // drops its mix lane before the new one joins under the same name
    kits.delete(p);
    say(`${PIECE_NAMES[p]}: this device could not compute it in real time, so it restarts on a coarser grid.`, "warn");
    // whole grids only: a 3.5x piece restarts on 4x or 5x, grids with measured levels
    const next = [Math.floor(scale) + 1, Math.floor(scale) + 2, 10].filter((s, i, a) => s <= 10 && a.indexOf(s) === i);
    await startPiece(mix, p, next);
    const h = health[p];
    if (h) h.detail += ` (restarted from ${scale}×${lostBefore > 480 ? `, after ${Math.round(lostBefore / 48)} ms were lost` : ""})`;
    const k = kits.get(p);
    if (!k || !arr || phase !== "playing") return;
    // join the score already playing, from a moment the new Worker can still reach
    k.lane.resync(epoch);
    if (h) { h.base = h.late; h.skipBase = h.stats?.skippedFrames ?? 0; }
    k.lane.score(kitScore(p, scoreTime() + LEAD_WARM_S));
  }

  function disposeEngines(): void {
    piano?.dispose();
    piano = null;
    pianoReady = null;
    for (const k of kits.values()) k.lane.dispose();
    kits.clear();
    mix?.dispose();
    mix = null;
    master?.disconnect();
    master = null;
    for (const k of Object.keys(health)) delete health[k];
    paintLanes();
  }

  // ------------------------------------------------------------- transport
  async function play(from = position): Promise<void> {
    if (!arr || phase === "warming") return;
    clearTimeout(idleTimer);
    claimPlayback(onClaimed);
    const my = token;
    setPhase("warming");
    try {
      const warm = await ensureEngines(arr);
      if (my !== token || !arr) return;
      cue(from >= arr.durationS - 0.05 ? 0 : from, warm ? LEAD_WARM_S : LEAD_COLD_S);
      setPhase("playing");
      say(playingLine());
    } catch (e) {
      setPhase("ready");
      say(`The engines could not start: ${(e as Error).message}`, "bad");
    }
  }

  /** Cut what the lanes rendered ahead and start the score at `from`, `lead` seconds from now. */
  function cue(from: number, lead: number): void {
    if (!ctx || !mix || !master || !arr) return;
    const c = ctx, a = arr;
    epoch++;
    silence();
    const now = c.currentTime;
    startFrame = Math.round((now + lead) * SR);
    startPos = from;
    master.gain.setValueAtTime(0, startFrame / SR - 0.03);
    master.gain.linearRampToValueAtTime(1, startFrame / SR - 0.005);
    const at = (t: number): number => startFrame + ((t - from) / tempoFactor) * SR;
    if (piano) {
      const ev: { frame: number; on?: number; off?: number; speed?: number; sustain?: number }[] = [];
      let pedal = 0;
      for (const p of a.pedal) { if (p.t > from) break; pedal = p.value; }
      if (pedal > 0) ev.push({ frame: at(from), sustain: pedal });
      for (const e of a.events) {
        if (e.t < from) continue;
        if (e.on !== undefined) ev.push({ frame: at(e.t), on: e.on, speed: e.speed });
        else if (e.off !== undefined) ev.push({ frame: at(e.t), off: e.off });
        else if (e.sustain !== undefined) ev.push({ frame: at(e.t), sustain: e.sustain });
      }
      piano.score(ev);
    }
    for (const [p, k] of kits) k.lane.score(kitScore(p, from));
    snaps.length = 0;
    ledgerPrev = null;
    flow = null;
    retiredBase = null;
    kit.clear();
    hitCursor = lowerBound(a.hits.map((h) => h.t), from);
    for (const h of Object.values(health)) { h.base = h.late; h.skipBase = h.stats?.skippedFrames ?? 0; }
  }

  /** A piece's strokes from score time `from` on, on the clock of the score now playing. */
  function kitScore(p: PieceId, from: number): KitScoreEvent[] {
    if (!arr) return [];
    const R = kits.get(p)?.ready?.radius ?? 0.2;
    const at = (t: number): number => startFrame + ((t - startPos) / tempoFactor) * SR;
    return arr.hits.filter((h) => h.piece === p && h.t >= from).map((h) => ({ frame: at(h.t), x: h.u * R, y: h.v * R, speed: STICK[p] * h.g }));
  }

  /** Mute now, drop everything rendered ahead, release every key; the lanes render for the clock as it is. */
  function silence(): void {
    if (!ctx || !mix || !master) return;
    const now = ctx.currentTime;
    master.gain.cancelScheduledValues(now);
    master.gain.setValueAtTime(master.gain.value, now);
    master.gain.linearRampToValueAtTime(0, now + 0.012);
    mix.flush({ epoch, from: Math.round((now + 0.014) * SR) });
    piano?.stopScore();
    piano?.resync(epoch);
    for (const k of kits.values()) { k.lane.stopScore(); k.lane.resync(epoch); }
  }

  function pause(): void {
    if (phase !== "playing") return;
    position = scoreTime();
    epoch++;
    silence();
    kit.clear();
    setPhase("paused");
    say(`Paused at ${fmtTime(position)}. The engines stay warm for a minute and a half.`);
    armIdle();
  }

  function onClaimed(): void {
    // another player on the page starts: stop, and free the CPU for it soon after
    if (phase === "playing") pause();
    clearTimeout(idleTimer);
    idleTimer = window.setTimeout(() => { if (phase !== "playing" && phase !== "warming") disposeEngines(); }, 2000);
  }
  registerPlayer(onClaimed);

  function armIdle(): void {
    clearTimeout(idleTimer);
    idleTimer = window.setTimeout(() => { if (phase !== "playing" && phase !== "warming") { disposeEngines(); say("The engines were stopped to save power; Play starts them again."); } }, IDLE_DISPOSE_MS);
  }

  function seekTo(t: number): void {
    if (!arr) return;
    const x = Math.max(0, Math.min(arr.durationS, t));
    if (phase === "playing") {
      cue(x, LEAD_WARM_S);
      say(playingLine());
    } else {
      position = x;
      paintPosition(x);
      requestDraw();
    }
  }

  function setTempo(f: number): void {
    const next = Math.round(Math.max(0.5, Math.min(1.5, f)) * 100) / 100;
    if (next === tempoFactor) return;
    const here = scoreTime();
    tempoFactor = next;
    showTempo();
    if (phase === "playing") { cue(here, LEAD_WARM_S); say(playingLine()); }
  }
  function showTempo(): void {
    tempoVal.textContent = `${Math.round(tempoFactor * 100)}%`;
    tempoBpm.textContent = arr ? `♩ ${Math.round(arr.tempoBpm * tempoFactor)}` : "";
    tempoDown.toggleAttribute("disabled", tempoFactor <= 0.5);
    tempoUp.toggleAttribute("disabled", tempoFactor >= 1.5);
  }

  function scoreTime(): number {
    if (phase !== "playing" || !ctx) return position;
    const heard = ctx.currentTime * SR - latencyFrames();
    return startPos + (Math.max(0, heard - startFrame) / SR) * tempoFactor;
  }
  function latencyFrames(): number {
    if (!ctx) return 0;
    const l = (Number.isFinite(ctx.outputLatency) ? ctx.outputLatency : 0) + (Number.isFinite(ctx.baseLatency) ? ctx.baseLatency : 0);
    return l * SR;
  }

  playBtn.addEventListener("click", () => {
    if (phase === "playing") pause();
    else { primeAudio(); void play(); }
  });
  miniPlay.addEventListener("click", () => { if (phase === "playing") pause(); else { primeAudio(); void play(); } });
  rewindBtn.addEventListener("click", () => seekTo(0));
  tempoDown.addEventListener("click", () => setTempo(tempoFactor - 0.05));
  tempoUp.addEventListener("click", () => setTempo(tempoFactor + 0.05));
  let seeking = false;
  seek.addEventListener("input", () => {
    if (!arr) return;
    seeking = true;
    const t = (Number(seek.value) / 1000) * arr.durationS;
    paintPosition(t, true);
    if (phase !== "playing") { position = t; requestDraw(); }
  });
  seek.addEventListener("change", () => {
    if (!arr) return;
    seeking = false;
    seekTo((Number(seek.value) / 1000) * arr.durationS);
  });
  roll.scrubFrom = () => scoreTime();
  roll.onScrub = (t, done) => {
    if (phase === "playing") { if (done) seekTo(t); return; }
    position = t;
    paintPosition(t);
    requestDraw();
  };
  window.addEventListener("keydown", (e) => {
    if (!arr || e.metaKey || e.ctrlKey || e.altKey) return;
    const target = e.target as HTMLElement;
    if (target.closest("input, textarea, select, [contenteditable]") || (target.closest("button") && (e.key === " " || e.key === "Enter"))) return;
    const r = panel.getBoundingClientRect();
    if (r.bottom < 80 || r.top > innerHeight - 80) return;
    if (e.key === " ") { e.preventDefault(); if (phase === "playing") pause(); else { primeAudio(); void play(); } }
    else if (e.key === "ArrowRight") { e.preventDefault(); seekTo(scoreTime() + 5); }
    else if (e.key === "ArrowLeft") { e.preventDefault(); seekTo(scoreTime() - 5); }
    else if (e.key === "]") setTempo(tempoFactor + 0.05);
    else if (e.key === "[") setTempo(tempoFactor - 0.05);
  });

  // tabs (phones show one 3D view at a time)
  const pick = (view: "board" | "kit"): void => {
    physics.dataset.view = view;
    tabBoard.setAttribute("aria-pressed", String(view === "board"));
    tabKit.setAttribute("aria-pressed", String(view === "kit"));
  };
  tabBoard.addEventListener("click", () => pick("board"));
  tabKit.addEventListener("click", () => pick("kit"));

  // ----------------------------------------------------------------- phases
  function setPhase(p: Phase): void {
    phase = p;
    const playing = p === "playing";
    playBtn.replaceChildren(svg(playing ? ICON.pause : ICON.play), el("span", {}, playing ? "Pause" : p === "warming" ? "Starting…" : p === "loading" ? "Loading…" : position > 0 && p === "paused" ? "Resume" : "Play"));
    playBtn.toggleAttribute("disabled", !arr || p === "warming" || p === "loading");
    playBtn.setAttribute("aria-label", playing ? "Pause" : "Play");
    rewindBtn.toggleAttribute("disabled", !arr || p === "warming" || p === "loading");
    panel.classList.toggle("pp-playing", playing);
    panel.classList.toggle("pp-busy", p === "warming" || p === "loading");
    roll.scrubbable = !!arr;
    miniPlay.replaceChildren(svg(playing ? ICON.pause : ICON.play));
    miniPlay.setAttribute("aria-label", playing ? "Pause" : "Play");
    updateMini();
  }

  function say(text: string, tone: "" | "warn" | "bad" = ""): void {
    status.textContent = text;
    status.dataset.tone = tone;
  }
  function readyLine(a: Arrangement): string {
    const n = a.pieces.length;
    return `${a.notes.length.toLocaleString("en-US")} piano notes${a.hits.length ? ` and ${a.hits.length.toLocaleString("en-US")} drum strokes on ${n} piece${n > 1 ? "s" : ""}` : ""}, ${fmtTime(a.durationS)}. Press Play.`;
  }
  function playingLine(): string {
    return tempoFactor === 1 ? "Playing. Every note is computed in this tab about a second before you hear it." : `Playing at ${Math.round(tempoFactor * 100)}% of the file's tempo.`;
  }

  // ------------------------------------------------------------------ lanes
  function paintLanes(): void {
    const rows: HTMLElement[] = [];
    const ids = ["piano", ...(arr?.pieces ?? [])];
    for (const id of ids) {
      const h = health[id];
      const name = id === "piano" ? "Piano" : PIECE_NAMES[id as PieceId];
      const lost = h ? lostFrames(h) / 48 : 0;
      const rt = h?.rt ?? 0;
      const tone = !h || rt === 0 ? "idle" : lost > 10 ? "bad" : rt < 1.1 ? "warn" : "ok";
      const meter = el("span", { class: "pp-meter" }, el("i", { style: `width:${Math.min(100, (rt / 4) * 100).toFixed(0)}%` }));
      const cap = id === "piano" && h?.stats?.modeCap && h.stats.fullModeCap && h.stats.modeCap < h.stats.fullModeCap ? h.stats.modeCap : 0;
      const grid = h?.scale ? ` · ${h.scale}× grid` : id === "piano" && pianoReady ? ` · ${pianoReady.boardCoordinates} modes${cap ? ` · strings ${cap}` : ""}` : "";
      const val = !h ? "starts with Play" : rt > 0 ? `${rt.toFixed(1)}×${grid}${lost > 10 ? ` · ${lost < 1000 ? `${Math.round(lost)} ms` : `${(lost / 1000).toFixed(1)} s`} lost` : ""}` : h.stats?.asleep || h.scale ? `resting${grid}` : "measuring…";
      const row = el("div", { class: `pp-lane ${tone}`, "data-id": id, "data-rt": rt.toFixed(2), "data-late": String(h ? lostFrames(h) : 0), "data-gap": String(h ? h.gap : 0), "data-lost-ms": lost.toFixed(1) },
        el("span", { class: "pp-lane-dot", "aria-hidden": "true" }),
        el("span", { class: "pp-lane-name" }, name),
        meter,
        el("span", { class: "pp-lane-val num" }, val));
      if (h?.detail) row.title = `${name}: ${h.detail}${cap ? `; string modes held to ${cap} of ${h?.stats?.fullModeCap} so this device keeps up` : ""}; ${rt > 0 ? `${rt.toFixed(2)}× faster than real time` : "not computing (silent pieces are not stepped)"}`;
      rows.push(row);
    }
    lanesEl.replaceChildren(el("p", { class: "pp-lanes-head" }, el("span", { class: "tag live" }, "Live · WebAssembly"), el("span", { class: "small muted" }, "× real time, one Worker each")), ...rows);
    // a slow machine is said out loud, never hidden: audio that arrived too late is a real gap
    const wall = performance.now();
    for (const id of ids) { const h = health[id]; if (h && lostFrames(h) > (h.lastLate ?? 0)) { h.lastLate = lostFrames(h); h.lateAt = wall; } }
    const late = ids.filter((id) => { const h = health[id]; return h && lostFrames(h) / 48 > 20 && wall - (h.lateAt ?? -1e9) < 3000; });
    if (phase === "playing" && late.length) {
      const who = late.map((id) => (id === "piano" ? "the piano" : PIECE_NAMES[id as PieceId].toLowerCase())).join(", ");
      say(`This device computed ${who} too slowly for a moment, so you heard gaps: they are real, not part of the music. ${late.includes("piano") ? "A slower tempo helps." : "A slow drum piece is restarted on a coarser grid."}`, "warn");
      lateWarned = true;
    } else if (lateWarned && phase === "playing" && status.dataset.tone === "warn" && late.length === 0) {
      lateWarned = false;
      say(playingLine());
    }
  }
  let lateWarned = false;
  const lostFrames = (h: Health): number => Math.max(0, h.late - h.base) + Math.max(0, (h.stats?.skippedFrames ?? 0) - h.skipBase);
  setInterval(() => { if (phase === "playing" || phase === "warming") paintLanes(); }, 500);

  // ---------------------------------------------------------------- report
  function renderReport(a: Arrangement): void {
    const r = a.report;
    const fact = (k: string, v: string, note = ""): HTMLElement => el("div", { class: "pp-fact" }, el("dt", {}, k), el("dd", {}, v, note ? el("small", {}, note) : null));
    const pianoFacts = el("dl", { class: "pp-facts" },
      fact("Piano notes", r.pianoNotesIn.toLocaleString("en-US"), a.notes.length !== r.pianoNotesIn ? `${a.notes.length.toLocaleString("en-US")} struck` : ""),
      fact("Folded into 88 keys", String(r.folded), "by octaves"),
      fact("Doublings merged", String(r.merged), "one key, one strike"),
      fact("Thinned", String(r.thinnedChord + r.thinnedDensity), r.thinnedChord + r.thinnedDensity ? `${r.thinnedChord} in chords over 10, ${r.thinnedDensity} over 16 per ¼ s` : "none needed"),
      fact("Sustain pedal", r.pedalChanges ? `${r.pedalChanges} changes` : "none", r.pedalSkipped ? `${r.pedalSkipped} tiny moves skipped` : ""),
      fact("Keys ringing at once", `up to ${r.peakKeys}`, r.peakKeys > 24 ? "the engine keeps 24; the quietest retire" : "within the engine's 24"),
      fact("Copies left out", String(r.partsLeftOut), r.partsLeftOut ? "parts doubling or echoing another" : "every part plays"),
      fact("Faded out by the file", String(r.silencedByVolume), "its own channel volume"),
      fact("Repeats too fast", String(r.tooFast), "a key repeats about 14 times a second"),
      fact("Pitch bends", r.pitchBends ? String(r.pitchBends) : "none", r.pitchBends ? "not played: a piano key cannot bend" : ""));
    const mapped = Object.entries(r.drumsMapped);
    const skipped = Object.entries(r.drumsSkipped);
    const counts = new Map<PieceId, number>();
    for (const h of a.hits) counts.set(h.piece, (counts.get(h.piece) ?? 0) + 1);
    const maxCount = Math.max(1, ...counts.values());
    const drumBlock = a.hits.length || r.drumNotesIn
      ? el("div", { class: "pp-drums" },
        el("ul", { class: "pp-drum-bars" }, ...[...counts].map(([p, n]) => el("li", {}, el("span", {}, PIECE_NAMES[p]), el("i", { style: `--w:${((n / maxCount) * 100).toFixed(1)}%` }), el("b", { class: "num" }, n.toLocaleString("en-US"))))),
        mapped.length ? el("p", { class: "small muted" }, el("b", {}, "Substituted: "), mapped.map(([k, n]) => `${k} (${n})`).join("; "), ".") : null,
        r.drumsLayered ? el("p", { class: "small muted" }, el("b", {}, "Layered: "), `${r.drumsLayered} drum notes struck together on one piece (a doubled kick) played as one stroke.`) : null,
        skipped.length ? el("p", { class: "small pp-skipped" }, el("b", {}, "Skipped, no matching piece: "), skipped.map(([k, n]) => `${k} ×${n}`).join(", "), ".") : null)
      : el("p", { class: "small muted" }, "No drum part (nothing on General MIDI channel 10).");
    const roleWord: Record<string, string> = { percussion: "drums", bass: "bass", melody: "melody", harmony: "accompaniment", ornament: "ornament", pad: "pad", silent: "silent" };
    const parts = el("table", { class: "pp-parts" },
      el("thead", {}, el("tr", {}, el("th", {}, "Track"), el("th", {}, "Part"), el("th", {}, "Role"), el("th", {}, "Notes"), el("th", {}, "What the player does"))),
      el("tbody", {}, ...a.parts.map((p) => el("tr", { class: p.leftOut ? "out" : "" },
        el("td", { class: "num", "data-label": "Track" }, `${p.track + 1}·${p.channel + 1}`),
        el("td", { "data-label": "Part" }, p.name, p.channel !== 9 && gmInstrument(p.program) && gmInstrument(p.program) !== p.name ? el("small", {}, ` ${gmInstrument(p.program)}`) : null),
        el("td", { "data-label": "Role" }, roleWord[p.role] ?? p.role),
        el("td", { class: "num", "data-label": "Notes" }, p.notesOut !== p.notes ? `${p.notes.toLocaleString("en-US")} → ${p.notesOut.toLocaleString("en-US")}` : p.notes.toLocaleString("en-US")),
        el("td", { class: "pp-part-do" }, `${p.channel === 9 ? "Drum kit" : p.leftOut ? "Not played" : "Grand piano"}: ${p.actions.map((x) => x.replace(/^left out: /, "")).join("; ")}`)))));
    const fileFacts = [
      `${a.parts.length} part${a.parts.length === 1 ? "" : "s"}`,
      `♩ ${a.tempoBpm.toFixed(a.tempoBpm % 1 ? 1 : 0)}${r.tempoChanges ? ` (${r.tempoChanges} tempo changes)` : ""}`,
      a.meter,
      fmtTime(a.durationS),
      r.ignoredControllers ? `${r.ignoredControllers.toLocaleString("en-US")} controller messages (modulation, pan and the like) not used; volume and expression shape each strike` : "",
    ].filter(Boolean).join(" · ");
    const repairs = [r.salvage, ...r.lenient].filter(Boolean) as string[];
    report.replaceChildren(
      el("header", { class: "pp-report-head" }, el("h3", {}, "What the roll holds"), el("p", { class: "small muted" }, fileFacts)),
      el("div", { class: "pp-report-grid" },
        el("div", {}, el("h4", {}, "Grand piano"), pianoFacts),
        el("div", {}, el("h4", {}, "Drum kit"), drumBlock),
        el("div", { class: "pp-parts-wrap" }, el("h4", {}, "Parts"), el("div", { class: "table-scroll" }, parts), repairs.length ? el("p", { class: "small muted" }, el("b", {}, "Read with repairs: "), repairs.join("; "), ".") : null)));
    report.hidden = false;
  }

  // ------------------------------------------------------------- the clock
  let visible = false;
  let lastWall = performance.now();
  let drawQueued = false;
  let lastChord = "";
  function requestDraw(): void {
    if (drawQueued) return;
    drawQueued = true;
    requestAnimationFrame(tick);
  }
  new IntersectionObserver((e) => { visible = e.some((x) => x.isIntersecting); if (visible) requestDraw(); }, { rootMargin: "120px" }).observe(panel);
  const views = { board: false, kit: false };
  new IntersectionObserver((e) => { for (const x of e) views[x.target === boardView ? "board" : "kit"] = x.isIntersecting; }).observe(boardView);
  new IntersectionObserver((e) => { for (const x of e) views[x.target === boardView ? "board" : "kit"] = x.isIntersecting; }).observe(kitView);

  function tick(wall: number): void {
    drawQueued = false;
    const dt = Math.min(0.1, (wall - lastWall) / 1000);
    lastWall = wall;
    const t = scoreTime();
    const a = arr;
    if (a && phase === "playing" && t >= a.durationS) {
      // the roll ran out: the last notes keep ringing in the engines
      position = 0;
      setPhase("ready");
      say(`Finished. ${readyLine(a).replace(/\. Press Play\.$/, "")}. Play to hear it again.`);
      kit.clear();
      armIdle();
    }
    if (visible) roll.draw(t, wall, phase === "playing");
    if (a && phase === "playing") {
      const nowFrame = ctx ? ctx.currentTime * SR - latencyFrames() : 0;
      // strokes being heard now: flash the pad and the stick mark
      while (hitCursor < a.hits.length && a.hits[hitCursor]!.t <= t) {
        const h = a.hits[hitCursor++]!;
        if (t - h.t < 0.25) { roll.flash(h.piece, wall); kit.stroke(h.piece, h.u, h.v, wall); }
      }
      // the soundboard as it is now
      let due: { frame: number; data: Float32Array } | null = null;
      while (snaps.length && snaps[0]!.frame <= nowFrame) due = snaps.shift()!;
      if (due) showSnapshot(due.frame, due.data, wall);
      if (views.kit || physics.dataset.view === "kit") kit.frame(nowFrame, wall, dt, true);
      if (!seeking) paintPosition(t);
      board.pulse(pulseLevel, wall);
    } else {
      kit.frame(Number.POSITIVE_INFINITY, wall, dt, false);
    }
    if (views.board) board.frame(wall, dt, phase === "playing");
    updateKitHud();
    // chord display
    if (a) {
      const k = upperBound(a.chords.map((c) => c.t), t + 0.02) - 1;
      const cur = k >= 0 ? a.chords[k]!.label : "";
      if (cur !== lastChord) {
        lastChord = cur;
        chordEl.replaceChildren(...chordNodes(cur || "—"));
        chordEl.classList.remove("pop");
        void chordEl.offsetWidth;
        chordEl.classList.add("pop");
        let j = k + 1;
        while (j < a.chords.length && a.chords[j]!.label === cur) j++;
        const nx = a.chords[j];
        nextEl.replaceChildren(...(nx ? [document.createTextNode("then "), ...chordNodes(nx.label)] : []));
        miniChord.replaceChildren(...chordNodes(cur));
      }
    }
    if (phase === "playing" || phase === "warming") requestDraw();
  }

  // ------------------------------------------------------- soundboard HUD
  let pulseLevel = 0;
  let flow: { input: number; string: number; board: number; damper: number; hammer: number } | null = null;
  let retiredBase: number | null = null;
  const hudLedger = el("div", { class: "pp-ledger" });
  const hudNums = el("div", { class: "pp-hud-nums" });
  boardHud.append(hudNums, hudLedger);
  function showSnapshot(frame: number, data: Float32Array, wall: number): void {
    const s = parsePianoSnapshot(data);
    board.paint(s.boardCoordinates);
    const L = s.ledger;
    const cur = { frame, input: L.inputJ, string: L.stringLossJ, board: L.boardLossJ, damper: L.damperLossJ, hammer: L.hammerLossJ + L.catchLossJ };
    if (ledgerPrev && frame > ledgerPrev.frame) {
      // power from the books' running totals, smoothed over about a third of a second
      const dt = ((frame - ledgerPrev.frame) / SR);
      const now = {
        input: Math.max(0, (cur.input - ledgerPrev.input) / dt), string: Math.max(0, (cur.string - ledgerPrev.string) / dt),
        board: Math.max(0, (cur.board - ledgerPrev.board) / dt), damper: Math.max(0, (cur.damper - ledgerPrev.damper) / dt), hammer: Math.max(0, (cur.hammer - ledgerPrev.hammer) / dt),
      };
      const k = Math.min(1, dt / 0.33);
      flow = flow ? { input: flow.input + (now.input - flow.input) * k, string: flow.string + (now.string - flow.string) * k, board: flow.board + (now.board - flow.board) * k, damper: flow.damper + (now.damper - flow.damper) * k, hammer: flow.hammer + (now.hammer - flow.hammer) * k } : now;
    }
    ledgerPrev = cur;
    const rates = flow;
    stringsCanvas.hidden = s.keys.length === 0;
    pulseLevel = Math.max(0, Math.min(1, (Math.log10(Math.max(1e-9, L.storedJ)) + 5) / 3));
    const sounding = s.keys.length;
    retiredBase ??= s.forcedRetirements; // the engine counts from its start; this counts from Play
    const retired = s.forcedRetirements - retiredBase;
    hudNums.replaceChildren(
      el("span", {}, el("b", {}, fmtJ(Math.max(0, L.storedJ))), " stored in strings and board"),
      el("span", {}, el("b", {}, `${sounding}`), ` of 24 keys sounding · pedal ${Math.round(s.sustain * 100)}%`, retired > 0 ? ` · ${retired} quiet key${retired > 1 ? "s" : ""} retired early` : ""),
      el("span", {}, el("b", {}, `${s.peakPa.toFixed(2)} Pa`), " peak at the microphone since the engines started"));
    if (rates) {
      const out = rates.string + rates.board + rates.damper + rates.hammer;
      const total = Math.max(1e-12, out);
      const seg = (cls: string, v: number, label: string): HTMLElement => el("i", { class: cls, style: `flex-basis:${((v / total) * 100).toFixed(1)}%`, title: `${label}: ${fmtW(v)}` });
      hudLedger.replaceChildren(
        el("div", { class: "pp-ledger-line small" }, el("span", {}, "hammers in ", el("b", {}, fmtW(rates.input))), el("span", {}, "lost ", el("b", {}, fmtW(out)))),
        el("div", { class: "ledger-bar" }, seg("s", rates.string, "strings"), seg("b", rates.board, "soundboard"), seg("d", rates.damper, "dampers"), seg("r", rates.hammer, "hammer felt")),
        el("div", { class: "pp-ledger-keys" }, el("span", {}, el("i", { class: "s" }), "strings"), el("span", {}, el("i", { class: "b" }), "board"), el("span", {}, el("i", { class: "d" }), "dampers"), el("span", {}, el("i", { class: "r" }), "felt")));
    }
    const scale = board.peakM > 0 ? (board.lift * 1.35) / board.peakM : 0;
    const tag = boardView.querySelector(".pp-board-scale");
    if (tag) tag.textContent = scale > 0 ? fmtExaggeration(scale) : "—";
    drawStrings(s.keys, wall);
  }

  function drawStrings(keys: ReturnType<typeof parsePianoSnapshot>["keys"], wall: number): void {
    const c = stringsCanvas;
    const w = c.clientWidth, h = c.clientHeight;
    if (!w || !h) return;
    const d = Math.min(2, devicePixelRatio || 1);
    if (c.width !== Math.round(w * d)) { c.width = Math.round(w * d); c.height = Math.round(h * d); }
    const g = c.getContext("2d")!;
    g.setTransform(d, 0, 0, d, 0, 0);
    g.clearRect(0, 0, w, h);
    const list = [...keys].sort((x, y) => y.relativeEnergyJ - x.relativeEnergyJ).slice(0, 8).sort((x, y) => x.midi - y.midi);
    if (!list.length) return;
    const lane = h / list.length;
    let peak = 0;
    for (const k of list) for (const y of k.stringShapeM) peak = Math.max(peak, Math.abs(y));
    const sc = peak > 0 ? (lane * 0.42) / peak : 0;
    g.font = '500 9.5px "JetBrains Mono", monospace';
    g.textBaseline = "middle";
    list.forEach((k, i) => {
      const y0 = (i + 0.5) * lane;
      const alpha = Math.max(0.25, Math.min(1, (Math.log10(Math.max(1e-14, k.relativeEnergyJ)) + 9.5) / 7));
      g.fillStyle = `rgba(241, 234, 219, ${0.45 + 0.4 * alpha})`;
      g.textAlign = "right";
      g.fillText(noteName(k.midi), 30, y0);
      g.beginPath();
      const n = k.stringShapeM.length;
      g.moveTo(36, y0);
      for (let j = 0; j < n; j++) g.lineTo(36 + ((j + 0.5) / n) * (w - 44), y0 - (k.stringShapeM[j] ?? 0) * sc);
      g.lineTo(w - 8, y0 - k.bridgeDisplacementM * sc);
      g.strokeStyle = k.hammerForceN > 0.01 ? "#f2a184" : `rgba(99, 210, 192, ${alpha})`;
      g.lineWidth = 1.4;
      g.stroke();
    });
    void wall;
  }

  // one row per piece the file plays: its kinetic energy now (log scale, 1 nJ to 100 mJ) and the
  // largest stick force of its last stroke
  let kitRows: { p: PieceId; bar: HTMLElement; e: HTMLElement; f: HTMLElement }[] = [];
  let kitRowsFor: Arrangement | null = null;
  let lastKitHud = 0;
  function updateKitHud(): void {
    if (!arr || !arr.pieces.length) {
      if (kitRowsFor !== arr) { kitRowsFor = arr; kitRows = []; kitHud.replaceChildren(el("span", { class: "small muted" }, arr ? "This file has no drum part (nothing on General MIDI channel 10)." : "The pieces a file plays light up here; the rest of the kit stays dimmed.")); }
      return;
    }
    if (kitRowsFor !== arr) {
      kitRowsFor = arr;
      kitRows = arr.pieces.map((p) => ({ p, bar: el("i", {}), e: el("span", { class: "num" }, "—"), f: el("span", { class: "num muted" }, "") }));
      kitHud.replaceChildren(
        el("p", { class: "pp-kit-head small muted" }, el("span", {}, "Kinetic energy now"), el("span", {}, "stick force")),
        el("div", { class: "pp-kit-rows" }, ...kitRows.map((r) => el("div", { class: "pp-kit-row" }, el("span", {}, PIECE_NAMES[r.p]), el("span", { class: "pp-kit-meter" }, r.bar), r.e, r.f))));
    }
    const now = performance.now();
    if (now - lastKitHud < 100) return;
    lastKitHud = now;
    for (const r of kitRows) {
      const { kinetic, force } = kit.energyOf(r.p);
      const live = phase === "playing" && kinetic > 1e-9;
      r.bar.style.width = live ? `${Math.max(2, Math.min(100, ((Math.log10(kinetic) + 9) / 8) * 100)).toFixed(1)}%` : "0%";
      r.e.textContent = live ? fmtJ(kinetic) : "at rest";
      r.f.textContent = live && force > 0.5 ? `${force.toFixed(0)} N` : "";
    }
  }

  // ------------------------------------------------------------- position
  function paintPosition(t: number, fromSlider = false): void {
    if (!arr) return;
    const d = arr.durationS;
    const clock = `${fmtTime(Math.min(t, d))} / ${fmtTime(d)}`;
    if (clockEl.textContent !== clock) { clockEl.textContent = clock; seek.setAttribute("aria-valuetext", `${fmtTime(Math.min(t, d))} of ${fmtTime(d)}`); }
    const b = upperBound(arr.bars.map((x) => x.t), t + 1e-6) - 1;
    barEl.textContent = b >= 0 ? `bar ${b + 1}` : "";
    if (!fromSlider) seek.value = String(Math.round((Math.min(t, d) / d) * 1000));
    seek.style.setProperty("--fill", `${((Math.min(t, d) / d) * 100).toFixed(2)}%`);
    miniBar.style.width = `${((Math.min(t, d) / d) * 100).toFixed(2)}%`;
  }

  // ---------------------------------------------------- phone mini player
  let sideVisible = true;
  // the transport card itself (on phones its column is display: contents and never intersects)
  new IntersectionObserver((e) => { sideVisible = e.some((x) => x.isIntersecting); updateMini(); }, { threshold: 0.05 }).observe(side.querySelector(".pp-now")!);
  let panelVisible = false;
  new IntersectionObserver((e) => { panelVisible = e.some((x) => x.isIntersecting); updateMini(); }).observe(panel);
  function updateMini(): void {
    mini.hidden = !(arr && (phase === "playing" || phase === "paused") && !sideVisible && panelVisible);
    miniTitle.textContent = meta?.title ?? "";
  }

  setPhase("empty");
  paintLanes();
  requestDraw();
  showPack(packs[0]?.id ?? "");
  shelfEdges();
  say("Pick a roll to hear it. The first press starts the engines in this tab, which takes a few seconds.");
}

// --------------------------------------------------------------------- helpers
/** A chord symbol with its accidentals in a font that has them (the display serif does not). */
function chordNodes(label: string): Node[] {
  return label.split(/([♭♯])/).filter(Boolean).map((t) => (t === "♭" || t === "♯" ? el("span", { class: "acc" }, t) : document.createTextNode(t)));
}
function hasFiles(e: DragEvent): boolean {
  return !!e.dataTransfer && [...e.dataTransfer.types].includes("Files");
}
function fmtTime(s: number): string {
  const t = Math.max(0, Math.round(s));
  return `${Math.floor(t / 60)}:${String(t % 60).padStart(2, "0")}`;
}
function fmtW(w: number): string {
  return w >= 1 ? `${w.toFixed(2)} W` : w >= 1e-3 ? `${(w * 1e3).toFixed(1)} mW` : w >= 1e-6 ? `${(w * 1e6).toFixed(1)} µW` : "0 W";
}
function fmtExaggeration(x: number): string {
  if (x >= 1e6) return `${(x / 1e6).toFixed(1)} million`;
  if (x >= 1e3) return `${Math.round(x / 1e3).toLocaleString("en-US")},000`;
  return Math.round(x).toLocaleString("en-US");
}
function lowerBound(a: number[], x: number): number {
  let lo = 0, hi = a.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m]! < x) lo = m + 1; else hi = m; }
  return lo;
}
function upperBound(a: number[], x: number): number {
  let lo = 0, hi = a.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m]! <= x) lo = m + 1; else hi = m; }
  return lo;
}

function evidenceBody(): HTMLElement {
  const p = (...c: (string | Node)[]): HTMLElement => el("p", { class: "small" }, ...c);
  return el("div", { class: "pp-ev" },
    el("div", {},
      el("h4", {}, "The mapping"),
      p("Every pitched channel plays the one grand piano; General MIDI channel 10 plays the drum kit. Notes outside A0–C8 are moved by octaves into the keyboard. Notes that start within 30 ms on the same key are one strike, and a key held by two parts is released only when the last of them ends."),
      p("Velocity sets the hammer's speed as it leaves the escapement, log-linear between 0.6 m/s (velocity 1), 2 m/s (64, mezzo-forte) and 5 m/s (127). Controller 64 is the sustain pedal, with half-pedal positions kept. Volume, expression and other controllers are not used."),
      p("Drum notes go to the piece their General MIDI number names: kick, snare (side stick near the rim, claps on the snare), rack and floor toms, crash (china and splash too), ride and its bell; bongos and congas to the toms, timbales to the snare with its wires off. Anything else is skipped and listed. Velocity scales the piece's standard stroke, and each stroke lands a little away from the last, so no two are identical."),
      p("Hi-hat strokes are left out for now, and counted. The live hi-hat rings like a chime: closed, it keeps sounding for seconds in a few sustained partials below 4 kHz, where a real closed hi-hat gives a short broadband chick, and open it rings in the same pure tones for over six seconds, because the solver has no friction between the two cymbals. A groove's sixteenth-note hats would pile up into a wash of chimes."),
      p("The engines have limits, and the roll keeps within them where it can: a chord of more than ten notes keeps its bass, its top and its loudest notes; no more than sixteen strikes in any quarter second. The piano simulates 24 sounding keys at once; beyond that the quietest strings are retired early.")),
    el("div", {},
      el("h4", {}, "What is computed"),
      p("The piano is engine-piano, the grand's own hammer, string, damper and soundboard equations for the keys that are sounding, on 32 soundboard coordinates. Each drum piece is fs-vkplate, the plate solver of the drum kit above, on a grid several times coarser, chosen per piece for this device when the engines start. Only the pieces a file uses are started."),
      p("Levels: one gain per instrument and nothing else. The piano plays at 0.8 of its engine's full scale (10 Pa at the microphone); each drum piece keeps the live kit's own gain, matched to a recorded kit's balance and rescaled for the grid it runs on, so that its standard stroke is as loud on every grid (coarser grids make the crash 5 dB louder and the closed hi-hat quieter; a quiet stroke is raised only while its peak stays under the kit's cap). Measured in these same engines on the first 12 seconds of six rolls, the loudest moments reach 0.44–0.70 of full scale, so nothing clips. There is no reverb, no equaliser and no limiter."),
      p("The roll's library comes from the owner's jazz chord editor: its chord charts, played as choruses through that editor's own bass and comping generators (the same voicings and grooves it plays), with drum parts written for each groove. Uploaded files go through that editor's MIDI salvage pass, which repairs re-struck, orphaned and unterminated notes and says so, and its rules for naming chords and track roles."),
      p("The soundboard view is the engine's modal sum on its own surface points, heights exaggerated by the factor shown. The kit shows each piece's computed surface velocity; each piece's colour scale follows its own motion down as it decays.")));
}
