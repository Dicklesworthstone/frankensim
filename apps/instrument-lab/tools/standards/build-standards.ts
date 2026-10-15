/**
 * build-standards.ts — the "standards" MIDI library for Struck's MIDI player.
 *
 * What it does: every chart in the chord-chart editor's reviewed progression library
 * (PROGRESSION_LIBRARY), plus its Deacon Blues starter chart (STARTER_CHART), is repeated
 * as several choruses inside ONE inserted section of a real studio composition, given the
 * entry's tempo and groove, and compiled through the editor's own pipeline
 * (compileStudioPlaybackPlan -> compilePerformancePlan, comp continuity v2). Bass and piano
 * comping are therefore exactly the notes, gates and velocities the editor plays; nothing
 * here voices a chord. This tool adds only: a one-bar ending that re-strikes the last
 * chord's fullest comp voicing and its root bass, the sustain pedal for ballad-comp, and a
 * drum part generated per groove (seeded by entry id, so rebuilds are byte-identical).
 * Output: SMF format 1, PPQ 960 — conductor (title, tempo, 4/4, a marker per source chord
 * with its exact source text, "Chorus N" markers), Bass on channel 0, Piano on channel 1,
 * Drums on GM channel 10 (index 9) — plus a catalog JSON.
 *
 * Run from the frankensim repo root (bun >= 1.1; the editor repo is only read):
 *
 *   bun apps/instrument-lab/tools/standards/build-standards.ts \
 *       ../jazz_chord_progression_editor_html apps/instrument-lab/public
 *
 * writes <public>/midi/standards/<id>.mid and <public>/midi/standards.json. Verify with
 *
 *   uv run --with mido python apps/instrument-lab/tools/standards/verify_standards.py \
 *       apps/instrument-lab/public
 *
 * Provenance: the editor's modules are imported from the path argument at run time, and
 * its commit (`git -C <path> rev-parse --short HEAD`, 0c3c7aa when this tool was written)
 * is recorded in standards.json as source.commit. Entries the pipeline refuses are listed
 * in standards.json "refused" with the stage and reason, never dropped silently.
 */
import { execFileSync } from "node:child_process";
import { mkdirSync, writeFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { pathToFileURL } from "node:url";

/* ------------------------------------------------------------------ constants */

const PPQ = 960;
const BAR = 4 * PPQ; // every library chart is 4/4; anything else is refused
const MAX_BARS = 128;
const TARGET_SECONDS = 60;
const MAX_CHORUSES = 8;
const SIZE_BUDGET_BYTES = 25 * 1024;

const CH_BASS = 0;
const CH_PIANO = 1;
const CH_DRUMS = 9;

/** Tempo for entries that do not commit one. */
const DEFAULT_TEMPO: Record<string, number> = {
  "ballad-comp@1": 72,
  "medium-swing@1": 152,
  "bossa-nova@1": 132,
  "straight-eighths@1": 104,
  "syncopated-sixteenths@1": 100,
  "uptempo-swing@1": 260,
};

const GROOVE_NAME: Record<string, string> = {
  "ballad-comp@1": "Ballad",
  "medium-swing@1": "Medium swing",
  "bossa-nova@1": "Bossa nova",
  "straight-eighths@1": "Straight eighths",
  "syncopated-sixteenths@1": "Syncopated sixteenths",
  "uptempo-swing@1": "Uptempo swing",
  "block-chords@1": "Block chords",
};

/**
 * Where the drums play an offbeat eighth, as a fraction of the beat. The drums follow each
 * groove's own declared ratio, so they land with the bass and comp the editor generates: its
 * ballad-comp@1 and medium-swing@1 tables play STRAIGHT eighths (1/2), and a 2/3 ride against
 * them put its skip note 66-139 ms away from the piano's offbeats. uptempo-swing@1 declares 3/5.
 */
const DRUM_SWING_OVERRIDE: Record<string, [number, number]> = {};

/** Catalog order: recognisable tunes first, then public-domain works, devices, studies. */
const ORDER = [
  "giant-steps",
  "what-a-fool-believes",
  "hello-its-me",
  "peg",
  "deacon-blues",
  "jazz-blues-f",
  "rhythm-turnaround",
  "bird-blues",
  "minor-blues",
  "gymnopedie",
  "pachelbel",
  "ragtime",
  "tritone-chain",
  "tristan",
  "commendatore",
  "lament-bass",
  "chopin-chromatic",
  "major-third-cycle",
  "two-five-one",
  "one-six-two-five",
  "minor-two-five-one",
  "dorian-vamp",
  "modal-planing",
  "chromatic-mediants",
  "lush-ballad-study",
  "gospel-blues-study",
  "mu-major-study",
  "whole-tone-study",
];

/* GM percussion */
const KICK = 36, SIDE = 37, SNARE = 38, CHH = 42, PHH = 44, OHH = 46, CRASH = 49, RIDE = 51;
const TOM_HI = 48, TOM_MID = 45, TOM_FLOOR = 43;

/* ---------------------------------------------------------------- types (local) */

type Entry = {
  id: string;
  title: string;
  kicker: string;
  note: string;
  provenance: string;
  chartText: string;
  grooveStyleId: string;
  tempoBpm?: number;
};
type Ratio = { numerator: number; denominator: number };
type PlanEvent = {
  eventId: string;
  startTick: number;
  durationTicks: number;
  gateDurationTicks: number;
  midiPitches: readonly number[];
  velocity: number;
};
type Plan = { events: readonly PlanEvent[]; totalTicks: number; tempoBpm: number; midiPpq: number };
type Provenance = { eventId: string; sourceEventId: string; role: string };
type Note = { tick: number; end: number; pitch: number; vel: number; ch: number };
type Hit = { tick: number; note: number; vel: number };
type FillKind = "none" | "small" | "big";

class Refusal extends Error {
  constructor(readonly stage: string, readonly reason: string) {
    super(`${stage}: ${reason}`);
  }
}

/* ------------------------------------------------------------------------ args */

const [repoArg, publicArg] = process.argv.slice(2);
if (repoArg === undefined || publicArg === undefined) {
  console.error("usage: bun build-standards.ts <jazz-repo-path> <instrument-lab-public-dir>");
  process.exit(2);
}
const REPO = resolve(repoArg);
const PUBLIC = resolve(publicArg);
const OUT_DIR = join(PUBLIC, "midi", "standards");

const git = (...args: string[]): string =>
  execFileSync("git", ["-C", REPO, ...args], { encoding: "utf8" }).trim();
const COMMIT = git("rev-parse", "--short", "HEAD");
const DIRTY = git("status", "--porcelain", "--untracked-files=no", "--", "src").length > 0;
if (DIRTY) console.warn(`warning: ${REPO} has uncommitted changes under src/; recording dirty: true`);

const load = async (rel: string): Promise<any> => import(pathToFileURL(join(REPO, rel)).href);
const controllerModule = await load("src/application/studio-controller.ts");
const libraryModule = await load("src/application/studio-progression-library.ts");
const starterModule = await load("src/application/studio-starter-chart.ts");
const studioPlayback = await load("src/application/studio-playback.ts");
const playback = await load("src/playback/index.ts");

const LIBRARY: readonly Entry[] = libraryModule.PROGRESSION_LIBRARY;
const STARTER = starterModule.STARTER_CHART as { title: string; chartText: string };
const STYLES: Record<string, { swingRatio: Ratio }> = playback.PERFORMANCE_STYLES;

/* ------------------------------------------------------------- the real pipeline */

function newComposition(): any {
  const created = controllerModule.createStudioComposition({});
  if (!created.ok) throw new Refusal("composition", created.refusal.code);
  return created.composition;
}

/** The studio's own defaults for a pristine document (Deacon Blues commits neither). */
function studioDefaults(): { tempoBpm: number; grooveStyleId: string } {
  const document = newComposition().readApplicationState().document;
  return {
    tempoBpm: document.tempoBpm,
    grooveStyleId: document.playback.grooveStyleId ?? studioPlayback.STUDIO_PERFORMANCE_STYLE,
  };
}

function chartBars(chartText: string): string[] {
  return chartText.split("|").map((bar) => bar.trim()).filter((bar) => bar.length > 0);
}

/**
 * The library test's insertion path: the opening bar fills the pristine measure, the
 * rest appends at the section end — one section, so comp continuity flows across
 * chorus repeats.
 */
function insertBars(controller: any, bars: readonly string[]): void {
  const snapshot = controller.getSnapshot();
  const sectionId = snapshot.sections[0]?.id;
  const measureId = snapshot.sections[0]?.measures[0]?.id;
  if (sectionId === undefined || measureId === undefined) throw new Refusal("insert", "no pristine measure");
  const steps: [string, object][] = [[`| ${bars[0]} |`, { kind: "measure-start", measureId }]];
  if (bars.length > 1) steps.push([`| ${bars.slice(1).join(" | ")} |`, { kind: "section-end", sectionId }]);
  for (const [text, target] of steps) {
    const preview = controller.previewChartText(text);
    if (preview.status !== "ready") {
      throw new Refusal("parse", `${preview.status} ${preview.issueCodes.join(",")}`);
    }
    const drafted = controller.setQuickEntryDraft(text, target, preview.status, preview.issueCodes);
    if (!drafted.ok) throw new Refusal("draft", drafted.refusal.code);
    const applied = controller.applyQuickEntryPreview();
    if (!applied.ok) throw new Refusal("insert", applied.refusal.code);
  }
}

type Performance = {
  base: Plan;
  plan: Plan;
  roles: Map<string, Provenance>;
  chordText: Map<string, string>;
};

function perform(entry: Entry, tempo: number, bars: readonly string[]): Performance {
  const composition = newComposition();
  const { controller } = composition;
  insertBars(controller, bars);
  // setTempo refuses a same-value change (command.payload_invalid), so skip the no-op.
  const tempoNow = composition.readApplicationState().document.tempoBpm;
  for (const [name, run] of [
    ["tempo", () => (tempoNow === tempo ? { ok: true } : controller.setTempo(tempo))],
    ["groove", () => controller.setPerformanceStyle(entry.grooveStyleId)],
    ["title", () => controller.setTitle(entry.title)],
  ] as const) {
    const result = run();
    if (!result.ok) throw new Refusal(name, result.refusal.code);
  }
  const document = composition.readApplicationState().document;
  const measures = document.sections.reduce((n: number, s: any) => n + s.measures.length, 0);
  if (document.sections.length !== 1 || measures !== bars.length) {
    throw new Refusal("document", `expected 1 section of ${bars.length} bars, got ${document.sections.length}/${measures}`);
  }
  if (document.meter.beatsPerBar !== 4 || document.meter.beatUnit !== 4) {
    throw new Refusal("document", `meter ${document.meter.beatsPerBar}/${document.meter.beatUnit} unsupported`);
  }
  const base = studioPlayback.compileStudioPlaybackPlan(document);
  if (!base.ok) throw new Refusal("playback-plan", `${base.refusal.code}: ${base.refusal.message}`);
  const performed = playback.compilePerformancePlan({
    plan: base.plan,
    styleId: entry.grooveStyleId,
    compContinuityVersion: 2,
  });
  if (!performed.ok) {
    throw new Refusal("performance-plan", `${performed.refusal.code} ${performed.refusal.detail}`);
  }
  const plan: Plan = performed.plan;
  if (plan.midiPpq !== PPQ || plan.tempoBpm !== tempo || plan.totalTicks !== bars.length * BAR) {
    throw new Refusal("performance-plan", `unexpected plan shape ppq=${plan.midiPpq} tempo=${plan.tempoBpm} ticks=${plan.totalTicks}`);
  }
  const roles = new Map<string, Provenance>(
    (performed.eventProvenance as Provenance[]).map((p) => [p.eventId, p]),
  );
  const chordText = new Map<string, string>();
  for (const section of document.sections) {
    for (const measure of section.measures) {
      for (const event of measure.events) chordText.set(event.id, event.chord.sourceText);
    }
  }
  return { base: base.plan, plan, roles, chordText };
}

/* --------------------------------------------------------------------- drums */

function fnv1a(text: string): number {
  let h = 0x811c9dc5;
  for (const byte of new TextEncoder().encode(text)) h = Math.imul(h ^ byte, 0x01000193) >>> 0;
  return h;
}

function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Fill on the last bar of each 4-bar phrase: "big" (snare into toms) at 8-bar phrase ends
 * and chorus ends, "small" (a pickup) at the 4-bar ends between. A 4-bar boundary closer
 * than 4 bars to the chorus end is skipped. Choruses shorter than 8 bars alternate a small
 * and a big fill at their ends, so a 4-bar vamp phrases in 8.
 */
function fillKind(barIndex: number, perChorus: number, choruses: number): FillKind {
  const i = barIndex % perChorus;
  const chorus = Math.floor(barIndex / perChorus);
  const last = perChorus - 1;
  if (i === last) {
    return perChorus >= 8 || chorus === choruses - 1 || chorus % 2 === 1 ? "big" : "small";
  }
  if ((i + 1) % 4 !== 0 || last - i < 4) return "none";
  return (i + 1) % 8 === 0 ? "big" : "small";
}

type BarContext = {
  b: number; // bar start tick
  odd: boolean; // global bar parity
  fill: FillKind;
  sw: number; // offbeat-eighth position inside the beat, ticks
  rand: () => number;
};

const Q = PPQ;
const between = (rand: () => number, lo: number, hi: number): number =>
  lo + Math.floor(rand() * (hi - lo + 1));

/** Hits are bar-relative; the caller adds the bar start. */
type Pattern = (c: BarContext) => Hit[];

const swingKit = (uptempo: boolean): Pattern => (c) => {
  const { sw, fill, rand } = c;
  const hits: Hit[] = [];
  const rideEnd = fill === "big" ? 2 * Q : BAR;
  const ride: [number, number][] = uptempo
    ? [[0, 64], [Q, 74], [Q + sw, 52], [2 * Q, 64], [3 * Q, 74], [3 * Q + sw, 52]]
    : [[0, 68], [Q, 76], [Q + sw, 56], [2 * Q, 68], [3 * Q, 76], [3 * Q + sw, 56]];
  for (const [t, v] of ride) if (t < rideEnd) hits.push({ tick: t, note: RIDE, vel: v });
  hits.push({ tick: Q, note: PHH, vel: 56 }, { tick: 3 * Q, note: PHH, vel: 56 });
  for (const beat of uptempo ? [0, 2] : [0, 1, 2, 3]) {
    hits.push({ tick: beat * Q, note: KICK, vel: between(rand, 25, 35) });
  }
  if (fill === "none") {
    // sparse snare comping: ghosted swung eighths, about one a bar
    const spots = uptempo ? [sw, Q + sw, 2 * Q + sw, 3 * Q + sw] : [sw, Q + sw, 2 * Q, 2 * Q + sw, 3 * Q + sw];
    const r = rand();
    const count = r < 0.1 ? 0 : r < 0.8 ? 1 : 2;
    const pool = [...spots];
    for (let k = 0; k < count; k += 1) {
      const [spot] = pool.splice(Math.floor(rand() * pool.length), 1);
      if (spot !== undefined) hits.push({ tick: spot, note: SNARE, vel: between(rand, 40, 70) });
    }
  } else if (fill === "small") {
    hits.push({ tick: 3 * Q, note: SNARE, vel: 58 }, { tick: 3 * Q + sw, note: SNARE, vel: 70 });
    hits.push({ tick: 3 * Q + sw, note: KICK, vel: 52 });
  } else if (uptempo) {
    // swung eighths: at 260-290 BPM a triplet fill is a blur
    hits.push(
      { tick: 2 * Q, note: SNARE, vel: 68 },
      { tick: 2 * Q + sw, note: SNARE, vel: 60 },
      { tick: 3 * Q, note: TOM_HI, vel: 76 },
      { tick: 3 * Q + sw, note: TOM_FLOOR, vel: 84 },
      { tick: 3 * Q + sw, note: KICK, vel: 58 },
    );
  } else {
    const t3 = Math.round(Q / 3);
    hits.push(
      { tick: 2 * Q, note: SNARE, vel: 66 },
      { tick: 2 * Q + t3, note: SNARE, vel: 52 },
      { tick: 2 * Q + 2 * t3, note: SNARE, vel: 72 },
      { tick: 3 * Q, note: TOM_HI, vel: 78 },
      { tick: 3 * Q + t3, note: TOM_MID, vel: 74 },
      { tick: 3 * Q + 2 * t3, note: TOM_FLOOR, vel: 84 },
      { tick: 3 * Q + 2 * t3, note: KICK, vel: 56 },
    );
  }
  return hits;
};

const PATTERNS: Record<string, Pattern> = {
  "ballad-comp@1": (c) => {
    const { sw, fill, odd } = c;
    const hits: Hit[] = [];
    const rideEnd = fill === "big" ? 3 * Q : BAR;
    for (const [t, v] of [[0, 56], [sw, 46], [Q, 50], [2 * Q, 54], [2 * Q + sw, 46], [3 * Q, 50]] as const) {
      if (t < rideEnd) hits.push({ tick: t, note: RIDE, vel: v });
    }
    hits.push({ tick: Q, note: PHH, vel: 48 }, { tick: 3 * Q, note: PHH, vel: 48 });
    hits.push({ tick: 0, note: KICK, vel: 30 });
    if (fill === "none" && odd) hits.push({ tick: 3 * Q, note: SIDE, vel: 46 });
    if (fill === "small") {
      hits.push({ tick: 3 * Q, note: SNARE, vel: 38 }, { tick: 3 * Q + sw, note: SNARE, vel: 48 });
    } else if (fill === "big") {
      const t3 = Math.round(Q / 3);
      hits.push(
        { tick: 2 * Q + sw, note: SNARE, vel: 40 },
        { tick: 3 * Q, note: TOM_HI, vel: 46 },
        { tick: 3 * Q + t3, note: TOM_MID, vel: 48 },
        { tick: 3 * Q + 2 * t3, note: TOM_FLOOR, vel: 54 },
      );
    }
    return hits;
  },
  "medium-swing@1": swingKit(false),
  "uptempo-swing@1": swingKit(true),
  "bossa-nova@1": (c) => {
    const { fill, odd } = c;
    const E = Q / 2;
    const hits: Hit[] = [];
    for (let k = 0; k < 8; k += 1) hits.push({ tick: k * E, note: CHH, vel: k % 2 === 0 ? 60 : 40 });
    hits.push(
      { tick: 0, note: KICK, vel: 64 },
      { tick: Q + E, note: KICK, vel: 52 },
      { tick: 2 * Q, note: KICK, vel: 60 },
      { tick: 3 * Q + E, note: KICK, vel: 52 },
    );
    const stick: [number, number][] = odd ? [[E, 56], [2 * Q, 60]] : [[0, 60], [Q + E, 56], [3 * Q, 58]];
    const stickEnd = fill === "big" ? 2 * Q : fill === "small" ? 3 * Q : BAR;
    for (const [t, v] of stick) if (t < stickEnd) hits.push({ tick: t, note: SIDE, vel: v });
    if (fill === "small") {
      hits.push({ tick: 3 * Q, note: TOM_HI, vel: 50 }, { tick: 3 * Q + E, note: TOM_MID, vel: 54 });
    } else if (fill === "big") {
      hits.push(
        { tick: 2 * Q, note: SNARE, vel: 52 },
        { tick: 2 * Q + E, note: SNARE, vel: 48 },
        { tick: 3 * Q, note: TOM_HI, vel: 56 },
        { tick: 3 * Q + E, note: TOM_MID, vel: 60 },
      );
    }
    return hits;
  },
  "straight-eighths@1": (c) => {
    const { fill, odd } = c;
    const E = Q / 2, S = Q / 4;
    const hits: Hit[] = [];
    const hatEnd = fill === "big" ? 2 * Q : BAR;
    for (let k = 0; k < 8; k += 1) {
      const t = k * E;
      if (t >= hatEnd) continue;
      if (fill === "small" && t === 3 * Q + E) hits.push({ tick: t, note: OHH, vel: 70 });
      else hits.push({ tick: t, note: CHH, vel: k % 2 === 0 ? 72 : 52 });
    }
    hits.push({ tick: 0, note: KICK, vel: 88 }, { tick: 2 * Q, note: KICK, vel: 80 });
    if (odd && fill !== "big") hits.push({ tick: 2 * Q + E, note: KICK, vel: 66 });
    hits.push({ tick: Q, note: SNARE, vel: 90 });
    if (fill !== "big") hits.push({ tick: 3 * Q, note: SNARE, vel: 94 });
    if (fill === "small") hits.push({ tick: 3 * Q + 3 * S, note: SNARE, vel: 58 });
    if (fill === "big") {
      hits.push(
        { tick: 2 * Q, note: SNARE, vel: 88 },
        { tick: 2 * Q + E, note: SNARE, vel: 80 },
        { tick: 3 * Q, note: TOM_HI, vel: 90 },
        { tick: 3 * Q + S, note: TOM_HI, vel: 78 },
        { tick: 3 * Q + 2 * S, note: TOM_MID, vel: 86 },
        { tick: 3 * Q + 3 * S, note: TOM_FLOOR, vel: 94 },
      );
    }
    return hits;
  },
  "syncopated-sixteenths@1": (c) => {
    const { fill, odd, rand } = c;
    const E = Q / 2, S = Q / 4;
    const hits: Hit[] = [];
    const openHat = odd && fill === "none";
    const hatEnd = fill === "big" ? 3 * Q : fill === "small" ? 3 * Q + S : BAR;
    for (let k = 0; k < 16; k += 1) {
      const t = k * S;
      if (t >= hatEnd) continue;
      if (openHat && t === 3 * Q + E) {
        hits.push({ tick: t, note: OHH, vel: 68 });
        continue;
      }
      if (openHat && t === 3 * Q + 3 * S) continue; // the open hat rings through the "a"
      hits.push({ tick: t, note: CHH, vel: k % 4 === 0 ? 72 : k % 2 === 0 ? 60 : 38 });
    }
    hits.push(
      { tick: 0, note: KICK, vel: 96 },
      { tick: Q + E, note: KICK, vel: 78 },
      { tick: Q + 3 * S, note: KICK, vel: 84 },
      { tick: 2 * Q + E, note: KICK, vel: 86 },
    );
    hits.push({ tick: Q, note: SNARE, vel: 96 }, { tick: 3 * Q, note: SNARE, vel: 98 });
    const ghostEnd = fill === "none" ? BAR : 3 * Q;
    for (const [t, p] of [[3 * S, 0.3], [Q + 3 * S, 0.35], [2 * Q + S, 0.35], [3 * Q + 3 * S, 0.3]] as const) {
      if (t < ghostEnd && rand() < p) hits.push({ tick: t, note: SNARE, vel: between(rand, 25, 35) });
    }
    if (fill === "small") {
      hits.push({ tick: 3 * Q + E, note: SNARE, vel: 72 }, { tick: 3 * Q + 3 * S, note: SNARE, vel: 82 });
    } else if (fill === "big") {
      hits.push(
        { tick: 3 * Q + S, note: SNARE, vel: 70 },
        { tick: 3 * Q + E, note: TOM_HI, vel: 86 },
        { tick: 3 * Q + 3 * S, note: TOM_FLOOR, vel: 94 },
      );
    }
    return hits;
  },
};

/** Crash + kick velocities on chorus downbeats and the ending. */
const ACCENT: Record<string, { crash: number; kick: number }> = {
  "ballad-comp@1": { crash: 70, kick: 50 },
  "medium-swing@1": { crash: 96, kick: 82 },
  "uptempo-swing@1": { crash: 94, kick: 80 },
  "bossa-nova@1": { crash: 82, kick: 72 },
  "straight-eighths@1": { crash: 108, kick: 100 },
  "syncopated-sixteenths@1": { crash: 104, kick: 100 },
};

function drumSwingTicks(styleId: string): number {
  const own = STYLES[styleId]?.swingRatio;
  const [n, d] = DRUM_SWING_OVERRIDE[styleId] ?? (own ? [own.numerator, own.denominator] : [1, 2]);
  return Math.round((PPQ * n) / d);
}

function drums(
  entry: Entry,
  tempo: number,
  perChorus: number,
  choruses: number,
  endingStart: number,
  endTick: number,
): Note[] {
  const pattern = PATTERNS[entry.grooveStyleId];
  const accent = ACCENT[entry.grooveStyleId];
  if (pattern === undefined || accent === undefined) {
    throw new Refusal("drums", `no drum pattern for ${entry.grooveStyleId}`);
  }
  const rand = mulberry32(fnv1a(entry.id));
  const sw = drumSwingTicks(entry.grooveStyleId);
  const grid: Hit[] = [];
  const bars = perChorus * choruses;
  for (let bar = 0; bar < bars; bar += 1) {
    const b = bar * BAR;
    let hits = pattern({ b, odd: bar % 2 === 1, fill: fillKind(bar, perChorus, choruses), sw, rand });
    if (bar % perChorus === 0) {
      // crash + kick on beat 1 of every chorus; the crash takes the ride/hat's place
      hits = hits.filter((h) => !(h.tick === 0 && (h.note === RIDE || h.note === CHH)));
      hits.push({ tick: 0, note: CRASH, vel: accent.crash }, { tick: 0, note: KICK, vel: accent.kick });
    }
    for (const h of hits) grid.push({ ...h, tick: b + h.tick });
  }
  grid.push(
    { tick: endingStart, note: CRASH, vel: Math.min(120, accent.crash + 6) },
    { tick: endingStart, note: KICK, vel: Math.min(120, accent.kick + 4) },
  );
  // one hit per note per grid tick (the louder survives)
  const unique = new Map<string, Hit>();
  for (const h of grid) {
    const key = `${h.tick}:${h.note}`;
    const held = unique.get(key);
    if (held === undefined || h.vel > held.vel) unique.set(key, h);
  }
  // light humanization: +-4 ms, +-6 velocity, drawn in grid order so it is deterministic
  const ticksPerMs = (tempo * PPQ) / 60000;
  const jitter = Math.round(4 * ticksPerMs);
  const ordered = [...unique.values()].sort((a, b) => a.tick - b.tick || a.note - b.note);
  const length = Math.min(120, Math.round(60 * ticksPerMs)); // drums ignore note-off; keep it short
  return ordered.map((h) => {
    const tick = Math.max(0, Math.min(endTick - 1, h.tick + between(rand, -jitter, jitter)));
    const vel = Math.max(1, Math.min(127, h.vel + between(rand, -6, 6)));
    return { tick, end: Math.min(endTick, tick + length), pitch: h.note, vel, ch: CH_DRUMS };
  });
}

/* ------------------------------------------------------------------- SMF writer */

type MidiEvent = { tick: number; rank: number; seq: number; bytes: number[] };

const utf8 = (text: string): number[] => [...new TextEncoder().encode(text)];

function vlq(value: number): number[] {
  if (!Number.isSafeInteger(value) || value < 0 || value > 0x0fffffff) throw new Error(`bad delta ${value}`);
  const out = [value & 0x7f];
  for (let v = value >>> 7; v > 0; v >>>= 7) out.unshift((v & 0x7f) | 0x80);
  return out;
}

const meta = (type: number, payload: number[]): number[] => [0xff, type, ...vlq(payload.length), ...payload];

/**
 * Same channel and pitch never overlap: a later attack cuts the earlier note at its own
 * tick (the note-off sorts first), and a second attack at the same tick is dropped in
 * favour of the louder one. Returns how many notes were cut or dropped.
 */
function resolveOverlaps(notes: Note[]): { notes: Note[]; repaired: number } {
  let repaired = 0;
  const byKey = new Map<string, Note[]>();
  for (const n of notes) {
    const key = `${n.ch}:${n.pitch}`;
    (byKey.get(key) ?? byKey.set(key, []).get(key)!).push({ ...n });
  }
  const out: Note[] = [];
  for (const list of byKey.values()) {
    list.sort((a, b) => a.tick - b.tick || b.vel - a.vel);
    for (let i = 0; i < list.length; i += 1) {
      const n = list[i]!;
      if (i > 0 && list[i - 1]!.tick === n.tick) {
        repaired += 1;
        continue;
      }
      const next = list.slice(i + 1).find((m) => m.tick > n.tick);
      if (next !== undefined && n.end > next.tick) {
        n.end = next.tick;
        repaired += 1;
      }
      out.push(n);
    }
  }
  return { notes: out, repaired };
}

function noteEvents(notes: readonly Note[]): MidiEvent[] {
  const events: MidiEvent[] = [];
  let seq = 0;
  for (const n of notes) {
    if (n.end <= n.tick) throw new Error(`empty note ${JSON.stringify(n)}`);
    events.push({ tick: n.tick, rank: 2, seq: seq++, bytes: [0x90 | n.ch, n.pitch, n.vel] });
    events.push({ tick: n.end, rank: 0, seq: seq++, bytes: [0x80 | n.ch, n.pitch, 64] });
  }
  return events;
}

function trackChunk(name: string, events: MidiEvent[], endTick: number): number[] {
  const sorted = [...events].sort((a, b) => a.tick - b.tick || a.rank - b.rank || a.seq - b.seq);
  const body: number[] = [0, ...meta(0x03, utf8(name))];
  let now = 0;
  for (const e of sorted) {
    if (e.tick > endTick) throw new Error(`event past end in ${name}`);
    body.push(...vlq(e.tick - now), ...e.bytes);
    now = e.tick;
  }
  body.push(...vlq(endTick - now), ...meta(0x2f, []));
  const n = body.length;
  return [0x4d, 0x54, 0x72, 0x6b, (n >>> 24) & 255, (n >>> 16) & 255, (n >>> 8) & 255, n & 255, ...body];
}

/* ----------------------------------------------------------------- one entry */

type Built = {
  catalog: Record<string, unknown>;
  bytes: Uint8Array;
  repaired: number;
  compileMs: number;
};

function build(entry: Entry): Built {
  const tempo = entry.tempoBpm ?? DEFAULT_TEMPO[entry.grooveStyleId];
  if (tempo === undefined) throw new Refusal("tempo", `no default tempo for ${entry.grooveStyleId}`);
  const chorusBars = chartBars(entry.chartText);
  const perChorus = chorusBars.length;
  if (perChorus === 0) throw new Refusal("chart", "empty chart");
  const chorusSeconds = (perChorus * 4 * 60) / tempo;
  let choruses = Math.min(MAX_CHORUSES, Math.max(1, Math.round(TARGET_SECONDS / chorusSeconds)));
  while (choruses > 1 && choruses * perChorus + 1 > MAX_BARS) choruses -= 1;
  if (choruses * perChorus + 1 > MAX_BARS) throw new Refusal("chart", `${perChorus} bars exceeds ${MAX_BARS - 1}`);
  const bars: string[] = [];
  for (let c = 0; c < choruses; c += 1) bars.push(...chorusBars);

  const t0 = performance.now();
  const { base, plan, roles, chordText } = perform(entry, tempo, bars);
  const compileMs = performance.now() - t0;

  const endingStart = plan.totalTicks;
  const endTick = endingStart + BAR;
  const bassNotes: Note[] = [];
  const pianoNotes: Note[] = [];
  const bySource = new Map<string, { bass: PlanEvent[]; comp: PlanEvent[] }>();
  for (const event of plan.events) {
    const origin = roles.get(event.eventId);
    if (origin === undefined) throw new Refusal("provenance", `no role for ${event.eventId}`);
    if (origin.role !== "bass" && origin.role !== "comp") {
      throw new Refusal("provenance", `role ${origin.role} is neither bass nor comp`);
    }
    const group = bySource.get(origin.sourceEventId) ?? { bass: [], comp: [] };
    group[origin.role].push(event);
    bySource.set(origin.sourceEventId, group);
    const target = origin.role === "bass" ? bassNotes : pianoNotes;
    for (const pitch of event.midiPitches) {
      target.push({
        tick: event.startTick,
        end: event.startTick + event.gateDurationTicks,
        pitch,
        vel: event.velocity,
        ch: origin.role === "bass" ? CH_BASS : CH_PIANO,
      });
    }
  }

  // Ending: the last chord's fullest comp voicing (latest on a tie) and its first bass
  // note — the performance compiler makes a chord's first bass note its root — held a bar.
  const sources = [...base.events].sort((a, b) => a.startTick - b.startTick);
  let endingComp: PlanEvent | undefined;
  let endingBass: PlanEvent | undefined;
  for (let i = sources.length - 1; i >= 0 && (endingComp === undefined || endingBass === undefined); i -= 1) {
    const group = bySource.get(sources[i]!.eventId);
    if (group === undefined) continue;
    if (endingComp === undefined && group.comp.length > 0) {
      endingComp = group.comp.reduce((best, e) =>
        e.midiPitches.length > best.midiPitches.length ||
        (e.midiPitches.length === best.midiPitches.length && e.startTick >= best.startTick) ? e : best);
    }
    if (endingBass === undefined && group.bass.length > 0) {
      endingBass = group.bass.reduce((first, e) => (e.startTick < first.startTick ? e : first));
    }
  }
  if (endingComp === undefined || endingBass === undefined) throw new Refusal("ending", "plan has no comp or no bass");
  for (const pitch of endingComp.midiPitches) {
    pianoNotes.push({ tick: endingStart, end: endTick, pitch, vel: endingComp.velocity, ch: CH_PIANO });
  }
  for (const pitch of endingBass.midiPitches) {
    bassNotes.push({ tick: endingStart, end: endTick, pitch, vel: endingBass.velocity, ch: CH_BASS });
  }

  const drumNotes = drums(entry, tempo, perChorus, choruses, endingStart, endTick);

  const bass = resolveOverlaps(bassNotes);
  const piano = resolveOverlaps(pianoNotes);
  const kit = resolveOverlaps(drumNotes);

  // Sustain pedal, ballad only: up 30 ms before each chord change, down 40 ms after.
  const pedal: MidiEvent[] = [];
  if (entry.grooveStyleId === "ballad-comp@1") {
    const ticksPerMs = (tempo * PPQ) / 60000;
    const up = Math.round(30 * ticksPerMs);
    const down = Math.round(40 * ticksPerMs);
    const changes = [...new Set([...sources.map((e) => e.startTick), endingStart])].sort((a, b) => a - b);
    let seq = 0;
    for (const t of changes) {
      if (t > 0) pedal.push({ tick: t - up, rank: 1, seq: seq++, bytes: [0xb0 | CH_PIANO, 64, 0] });
      pedal.push({ tick: t + down, rank: 1, seq: seq++, bytes: [0xb0 | CH_PIANO, 64, 127] });
    }
    pedal.push({ tick: endTick, rank: 1, seq: seq++, bytes: [0xb0 | CH_PIANO, 64, 0] });
  }

  // Conductor: title, tempo, 4/4, "Chorus N" then source-chord markers.
  const microseconds = Math.round(60_000_000 / tempo);
  const conductor: MidiEvent[] = [
    { tick: 0, rank: 0, seq: 0, bytes: meta(0x51, [(microseconds >>> 16) & 255, (microseconds >>> 8) & 255, microseconds & 255]) },
    { tick: 0, rank: 0, seq: 1, bytes: meta(0x58, [4, 2, 24, 8]) },
  ];
  let seq = 2;
  for (let c = 0; c < choruses; c += 1) {
    conductor.push({ tick: c * perChorus * BAR, rank: 1, seq: seq++, bytes: meta(0x06, utf8(`Chorus ${c + 1}`)) });
  }
  for (const event of sources) {
    const text = chordText.get(event.eventId);
    if (text === undefined) throw new Refusal("markers", `no source text for ${event.eventId}`);
    conductor.push({ tick: event.startTick, rank: 2, seq: seq++, bytes: meta(0x06, utf8(text)) });
  }

  const header = [0x4d, 0x54, 0x68, 0x64, 0, 0, 0, 6, 0, 1, 0, 4, (PPQ >>> 8) & 255, PPQ & 255];
  const chunks = [
    header,
    trackChunk(entry.title, conductor, endTick),
    trackChunk("Bass", noteEvents(bass.notes), endTick),
    trackChunk("Piano", [...noteEvents(piano.notes), ...pedal], endTick),
    trackChunk("Drums", noteEvents(kit.notes), endTick),
  ];
  const bytes = Uint8Array.from(chunks.flat());

  const totalBars = choruses * perChorus + 1;
  const durationS = Math.round(((endTick / PPQ) * 60 * 100) / tempo) / 100;
  return {
    bytes,
    repaired: bass.repaired + piano.repaired,
    compileMs,
    catalog: {
      id: entry.id,
      title: entry.title,
      kicker: entry.kicker,
      note: entry.note,
      provenance: entry.provenance,
      groove: entry.grooveStyleId,
      grooveName: GROOVE_NAME[entry.grooveStyleId] ?? entry.grooveStyleId,
      tempoBpm: tempo,
      barsPerChorus: perChorus,
      choruses,
      bars: totalBars,
      durationS,
      file: `midi/standards/${entry.id}.mid`,
      pianoNotes: piano.notes.length,
      bassNotes: bass.notes.length,
      drumHits: kit.notes.length,
      drumNotes: [...new Set(kit.notes.map((n) => n.pitch))].sort((a, b) => a - b),
    },
  };
}

/* ---------------------------------------------------------------------- main */

const defaults = studioDefaults();
const deaconBlues: Entry = {
  id: "deacon-blues",
  title: STARTER.title,
  kicker: "Becker & Fagen, 1977",
  note: "The editor's owner-directed starter chart: the intro's chromatic descent, two chords a bar, each major seventh answered by a minor-seven-sharp-five a semitone down, closing on the Ebmaj7 / E7#9 turnaround.",
  provenance: "owner-directed",
  chartText: STARTER.chartText,
  grooveStyleId: defaults.grooveStyleId,
  tempoBpm: defaults.tempoBpm,
};

const sourceOrder: Entry[] = [...LIBRARY, deaconBlues];
const rank = (id: string): number => {
  const i = ORDER.indexOf(id);
  return i < 0 ? ORDER.length + sourceOrder.findIndex((e) => e.id === id) : i; // unlisted ids keep library order
};
const all = [...sourceOrder].sort((a, b) => rank(a.id) - rank(b.id));

mkdirSync(OUT_DIR, { recursive: true });
const entries: Record<string, unknown>[] = [];
const refused: { id: string; stage: string; reason: string }[] = [];
for (const entry of all) {
  try {
    const built = build(entry);
    writeFileSync(join(OUT_DIR, `${entry.id}.mid`), built.bytes);
    entries.push(built.catalog);
    const c = built.catalog;
    const size = built.bytes.length;
    console.log(
      `${entry.id.padEnd(22)} ${String(c.grooveName).padEnd(22)} ${String(c.tempoBpm).padStart(3)} bpm ` +
        `${String(c.choruses)}x${String(c.barsPerChorus)}+1=${String(c.bars).padStart(3)} bars ` +
        `${String(c.durationS).padStart(6)} s  piano ${String(c.pianoNotes).padStart(4)}  drums ${String(c.drumHits).padStart(4)}  ` +
        `${(size / 1024).toFixed(1)} KB  compile ${built.compileMs.toFixed(0)} ms` +
        (built.repaired > 0 ? `  overlaps-cut ${built.repaired}` : "") +
        (size > SIZE_BUDGET_BYTES ? "  OVER SIZE BUDGET" : ""),
    );
  } catch (error) {
    const stage = error instanceof Refusal ? error.stage : "exception";
    const reason = error instanceof Refusal ? error.reason : String(error);
    refused.push({ id: entry.id, stage, reason });
    console.error(`REFUSED ${entry.id} at ${stage}: ${reason}`);
  }
}

const catalog = {
  source: { repo: "jazz_chord_progression_editor_html", commit: COMMIT, ...(DIRTY ? { dirty: true } : {}) },
  entries,
  refused,
};
writeFileSync(join(PUBLIC, "midi", "standards.json"), `${JSON.stringify(catalog, null, 2)}\n`);
console.log(`\n${entries.length} written, ${refused.length} refused -> ${OUT_DIR}`);
