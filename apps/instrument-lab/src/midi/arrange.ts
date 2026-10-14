// From a decoded MIDI file to a performance on FrankenSim's live piano and drum kit.
//
// Every pitched channel goes to the one grand piano and GM channel 10 to the drum kit. Nothing is
// changed silently: each note folded into the keyboard's range, merged with a doubling, or thinned
// for the engine's budget is counted, and every drum note is either mapped to a kit piece (by its
// General MIDI meaning) or skipped and listed by name.

import { chordTimeline } from "./chords";
import { orchestrate, type Role, type SourcePart } from "./orchestrate";
import { tempoClock, type DecodedMidi, type DecodedTrack, type PairedNote } from "./smf";

export type PieceId = "kick" | "snare" | "snare-off" | "tom-12" | "tom-14" | "hihat-closed" | "hihat-open" | "crash" | "ride";
export const PIECE_NAMES: Record<PieceId, string> = {
  kick: "Kick", snare: "Snare", "snare-off": "Snare, wires off", "tom-12": "Rack tom", "tom-14": "Floor tom",
  "hihat-closed": "Hi-hat", "hihat-open": "Open hi-hat", crash: "Crash", ride: "Ride",
};

/** A key the piano plays: times in seconds of the file at its own tempo. */
export interface PianoNote { key: number; t: number; off: number; vel: number; speed: number }
/** A key event for the engine (strike, release) or a pedal position. */
export interface PianoEvent { t: number; on?: number; off?: number; speed?: number; sustain?: number }
/** A stroke: where on the piece, as fractions of its radius, and a factor on its standard stick speed. */
export interface DrumHit { t: number; piece: PieceId; u: number; v: number; g: number; note: number; vel: number }
export type { Role } from "./orchestrate";
/** A part of the file (one track and channel) and what the player did with it. */
export interface Part { track: number; channel: number; name: string; role: Role; notes: number; program: number | null; actions: string[]; notesOut: number; leftOut: boolean }

export interface Arrangement {
  title: string;
  durationS: number;
  /** End of the last note or stroke (the rest of durationS is ringing). */
  lastEventS: number;
  notes: PianoNote[];
  events: PianoEvent[];
  pedal: { t: number; value: number }[];
  hits: DrumHit[];
  pieces: PieceId[];
  chords: { t: number; label: string }[];
  sections: { t: number; text: string }[];
  /** Bar starts from the file's own meter map (4/4 until a time signature says otherwise). */
  bars: { t: number; n: number }[];
  parts: Part[];
  tempoBpm: number;
  meter: string;
  report: {
    pianoNotesIn: number;
    drumNotesIn: number;
    folded: number;
    merged: number;
    thinnedChord: number;
    thinnedDensity: number;
    pedalChanges: number;
    pedalSkipped: number;
    peakKeys: number;
    tempoChanges: number;
    drumsMapped: Record<string, number>;
    drumsSkipped: Record<string, number>;
    ignoredControllers: number;
    salvage: string | null;
    lenient: string[];
    /** notes the file's own channel volume or expression had faded to silence */
    silencedByVolume: number;
    /** parts left out as copies (doublings or echoes) of another part */
    partsLeftOut: number;
    /** strikes closer together on one key than its action can repeat */
    tooFast: number;
    /** drum notes that were one stroke layered (two kick notes at once) */
    drumsLayered: number;
    pitchBends: number;
  };
}

// ----------------------------------------------------------------- calibrations
/** MIDI velocity to post-escapement hammer speed [m/s]: 1 → 0.6, 64 → 2 (mf), 127 → 5, log-linear between. */
export function hammerSpeed(velocity: number): number {
  const v = Math.max(1, Math.min(127, velocity));
  const [v0, s0, v1, s1] = v <= 64 ? [1, 0.6, 64, 2] : [64, 2, 127, 5];
  return Math.exp(Math.log(s0) + ((v - v0) / (v1 - v0)) * (Math.log(s1) - Math.log(s0)));
}
/** MIDI velocity to a factor on a piece's standard stroke: 100 → 1, 64 → 0.68, 127 → 1.24 (clamped 0.25..1.35). */
export function stickFactor(velocity: number): number {
  return Math.max(0.25, Math.min(1.35, (Math.max(1, velocity) / 100) ** 0.88));
}

/** Notes starting within this of each other on one key are one strike (a doubled part). */
const MERGE_S = 0.03;
/** A chord: onsets within 30 ms. More strikes than this in one chord are thinned (outer and loudest kept). */
const MAX_CHORD = 10;
/** No more than this many strikes in any 250 ms (the least important and quietest go first). */
const MAX_PER_QUARTER_SECOND = 16;
/** A piano key repeats at most about 14 times a second (double escapement). */
const REPEAT_S = 0.07;
/** A drummer's ordinary stroke weights, for drum parts written all at one velocity. */
const DRUM_WEIGHT: Partial<Record<PieceId, number>> = { kick: 92, snare: 96, "snare-off": 90, "tom-12": 90, "tom-14": 90, crash: 88, ride: 80 };

// ------------------------------------------------------------------ drum map
interface DrumRule { piece: PieceId; u: number; v: number; g: number; label?: string }
const DRUM_MAP: Record<number, DrumRule> = {
  35: { piece: "kick", u: 0, v: 0.02, g: 1 },
  36: { piece: "kick", u: 0, v: 0.02, g: 1 },
  37: { piece: "snare", u: 0.86, v: 0.05, g: 0.75, label: "Side stick → snare, near the rim" },
  38: { piece: "snare", u: 0.35, v: 0.1, g: 1 },
  39: { piece: "snare", u: 0.55, v: -0.2, g: 0.8, label: "Hand clap → snare" },
  40: { piece: "snare", u: 0.3, v: 0.05, g: 1.05 },
  41: { piece: "tom-14", u: 0.35, v: 0.1, g: 1 },
  43: { piece: "tom-14", u: 0.3, v: 0.12, g: 1 },
  // closed (42), pedal (44) and open (46) hi-hat: see HAT_RINGS
  45: { piece: "tom-12", u: 0.4, v: 0.1, g: 1 },
  47: { piece: "tom-12", u: 0.35, v: 0.1, g: 1 },
  48: { piece: "tom-12", u: 0.3, v: 0.12, g: 1 },
  50: { piece: "tom-12", u: 0.25, v: 0.1, g: 1.05 },
  49: { piece: "crash", u: 0.85, v: 0, g: 1 },
  57: { piece: "crash", u: 0.8, v: 0.1, g: 1 },
  52: { piece: "crash", u: 0.88, v: -0.1, g: 1, label: "China cymbal → crash" },
  55: { piece: "crash", u: 0.9, v: 0.05, g: 0.8, label: "Splash cymbal → crash" },
  51: { piece: "ride", u: 0.6, v: 0, g: 1 },
  59: { piece: "ride", u: 0.62, v: 0.06, g: 1 },
  53: { piece: "ride", u: 0.12, v: 0, g: 1.1, label: "Ride bell → ride, on the bell" },
  60: { piece: "tom-12", u: 0.7, v: 0, g: 0.8, label: "Bongo → rack tom, off centre" },
  61: { piece: "tom-12", u: 0.75, v: 0.1, g: 0.75, label: "Bongo → rack tom, off centre" },
  62: { piece: "tom-14", u: 0.7, v: 0, g: 0.6, label: "Conga → floor tom" },
  63: { piece: "tom-14", u: 0.55, v: 0.05, g: 0.8, label: "Conga → floor tom" },
  64: { piece: "tom-14", u: 0.5, v: 0.1, g: 0.85, label: "Conga → floor tom" },
  65: { piece: "snare-off", u: 0.4, v: 0.1, g: 1, label: "Timbale → snare with the wires off" },
  66: { piece: "snare-off", u: 0.5, v: 0.1, g: 0.9, label: "Timbale → snare with the wires off" },
};
/**
 * The live hi-hat rings like a chime: struck with the pedal down it keeps sounding for seconds in
 * sustained tonal partials (70 Hz-4 kHz in its spectrogram, not the short broadband chick of a real
 * closed hat), and open it rings in the same few pure tones for over six seconds, because the
 * solver has no friction between the two cymbals. A groove's sixteenth-note hats pile up into a
 * wash of chimes, so hi-hat strokes are left out and listed. Put them back (42 and 44 to
 * "hihat-closed", 46 to "hihat-open") when the solver damps the clamp.
 */
const HAT_RINGS = "Hi-hat, closed, pedal and open (the live hi-hat rings like a chime)";
const LEFT_OUT: Record<number, string> = { 42: HAT_RINGS, 44: HAT_RINGS, 46: HAT_RINGS };
/** "Tempo", "Conductor", "untitled", "tk6" and the like are labels, not names. */
const LABEL = /^(?:\s*(?:conductor|tempo|control|meta|master|global|setup|sequence|track|untitled|copyright|midi|smf|time\s*sig\w*)\b.*|\s*(?:tk|trk|tr|ch)?\s*\d*\s*)$/i;
const GM_PERCUSSION: Record<number, string> = {
  27: "High Q", 28: "Slap", 29: "Scratch push", 30: "Scratch pull", 31: "Sticks", 32: "Square click", 33: "Metronome click", 34: "Metronome bell",
  54: "Tambourine", 56: "Cowbell", 58: "Vibraslap", 67: "High agogô", 68: "Low agogô", 69: "Cabasa", 70: "Maracas", 71: "Short whistle",
  72: "Long whistle", 73: "Short güiro", 74: "Long güiro", 75: "Claves", 76: "High wood block", 77: "Low wood block", 78: "Mute cuíca",
  79: "Open cuíca", 80: "Mute triangle", 81: "Open triangle", 82: "Shaker", 83: "Jingle bell", 84: "Bell tree", 85: "Castanets", 86: "Mute surdo", 87: "Open surdo",
};

// ------------------------------------------------------- track roles (copied)
// Provenance: the M1-ROLE rule order and constants of jazz_chord_progression_editor_html (same
// owner, commit 0c3c7aa), src/export/midi-import-automation.ts `classifyTrack` and
// src/export/midi-import-automation-contract.ts: whole-token names beat statistics; GM channel 10
// is percussion; bass by register at most MIDI 55 when mostly monophonic; melody by a mean key of
// at least 64 when mostly monophonic; otherwise harmony.
const PERCUSSION_TOKENS = ["drum", "drums", "perc", "percussion", "kit"];
const BASS_TOKENS = ["bass", "bs", "contrabass", "upright"];
const MELODY_TOKENS = ["melody", "lead", "vocal", "voice", "sax", "trumpet", "flute", "horn"];
function classify(name: string, notes: readonly PairedNote[], channel: number): Role {
  const attacks = notes.length;
  if (!attacks) return "silent";
  const tokens = name.toLowerCase().split(/[\s\-_0-9]+/).filter(Boolean);
  const has = (list: string[]): boolean => tokens.some((t) => list.includes(t));
  if (channel === 9 || has(PERCUSSION_TOKENS)) return "percussion";
  const onTicks = new Map<number, number>();
  for (const n of notes) onTicks.set(n.onTick, (onTicks.get(n.onTick) ?? 0) + 1);
  let chordAttacks = 0, sumKeys = 0, maxKey = 0;
  for (const n of notes) {
    if ((onTicks.get(n.onTick) ?? 0) > 1) chordAttacks++;
    sumKeys += n.key;
    maxKey = Math.max(maxKey, n.key);
  }
  const monophonic = 4 * chordAttacks <= attacks;
  if (has(BASS_TOKENS)) return "bass";
  if (has(MELODY_TOKENS)) return "melody";
  if (maxKey <= 55 && monophonic) return "bass";
  if (monophonic && 2 * sumKeys >= 2 * 64 * attacks) return "melody";
  return "harmony";
}

// ------------------------------------------------------------------- arrange
export function arrange(midi: DecodedMidi, fileName: string): Arrangement {
  const sec = tempoClock(midi);
  const report: Arrangement["report"] = {
    pianoNotesIn: 0, drumNotesIn: 0, folded: 0, merged: 0, thinnedChord: 0, thinnedDensity: 0, pedalChanges: 0, pedalSkipped: 0,
    peakKeys: 0, tempoChanges: Math.max(0, midi.tempoMap.length - 1), drumsMapped: {}, drumsSkipped: {}, ignoredControllers: 0,
    salvage: midi.salvageNote, lenient: [], silencedByVolume: 0, partsLeftOut: 0, tooFast: 0, drumsLayered: 0, pitchBends: 0,
  };
  const l = midi.lenient;
  if (l.missingEndOfTrack) report.lenient.push(`${l.missingEndOfTrack} track${l.missingEndOfTrack > 1 ? "s" : ""} without an end marker read to the end of the chunk`);
  if (l.truncatedChunks) report.lenient.push(`${l.truncatedChunks} chunk${l.truncatedChunks > 1 ? "s" : ""} longer than the file read as far as the file goes`);
  if (l.openNotesClosed) report.lenient.push(`${l.openNotesClosed} note${l.openNotesClosed > 1 ? "s" : ""} never released, ended with the track`);
  if (l.zeroLengthNotes) report.lenient.push(`${l.zeroLengthNotes} zero-length note${l.zeroLengthNotes > 1 ? "s" : ""} given one tick`);

  // ------------------------------------------------- channel volume and expression
  // MIDI controllers are per channel, whichever track sends them; a note plays at its written
  // velocity times the channel's volume (CC7) and expression (CC11) when it starts. A file that
  // sends neither plays as written.
  const cc = new Map<number, { t: number; c: number; v: number }[]>();
  const pedalByChannel = new Map<number, { t: number; value: number }[]>();
  for (const track of midi.tracks) for (const c of track.controls) {
    if (c.controller === 7 || c.controller === 11) {
      const list = cc.get(c.channel) ?? [];
      list.push({ t: sec(c.tick), c: c.controller, v: c.value });
      cc.set(c.channel, list);
    } else if (c.controller === 64 && c.channel !== 9) {
      const list = pedalByChannel.get(c.channel) ?? [];
      list.push({ t: sec(c.tick), value: c.value / 127 });
      pedalByChannel.set(c.channel, list);
    } else if (c.controller !== 0 && c.controller !== 32 && c.controller < 120) report.ignoredControllers++;
  }
  for (const list of cc.values()) list.sort((a, b) => a.t - b.t);
  const levelAt = (ch: number, t: number): number => {
    const list = cc.get(ch);
    if (!list) return 1;
    let vol = 127, expr = 127;
    for (const e of list) { if (e.t > t + 1e-6) break; if (e.c === 7) vol = e.v; else expr = e.v; }
    return (vol / 127) * (expr / 127);
  };
  report.pitchBends = midi.counters.pitchBends;

  // ---------------------------------------------------------------- parts
  const parts: Part[] = [];
  const sources: SourcePart[] = [];
  const sourceOf: number[] = []; // parts index of each source
  const hits: DrumHit[] = [];
  let jitter = 0x9e3779b9;
  const rand = (): number => { jitter ^= jitter << 13; jitter ^= jitter >>> 17; jitter ^= jitter << 5; return ((jitter >>> 0) / 4294967296) - 0.5; };
  for (const track of midi.tracks) {
    const channels = [...new Set(track.notes.map((n) => n.channel))].sort((a, b) => a - b);
    for (const ch of channels) {
      const notes = track.notes.filter((n) => n.channel === ch);
      const program = programOf(track, ch);
      // a track's name names a part only when it is a name and the track holds one part
      const own = [track.name, track.instrumentName].find((x) => x && !LABEL.test(x) && channels.length === 1);
      const name = own || (ch === 9 ? "Drums" : program !== null ? GM_INSTRUMENTS[program] ?? `Channel ${ch + 1}` : `Channel ${ch + 1}`);
      const fileName = ch === 9 ? "Drums" : name;
      const role = classify(`${track.name ?? ""} ${track.instrumentName ?? ""}`, notes, ch);
      parts.push({ track: track.index, channel: ch, name: fileName, role, notes: notes.length, program, actions: [], notesOut: notes.length, leftOut: false });
      if (ch === 9) {
        drumPart(notes, parts[parts.length - 1]!);
      } else {
        report.pianoNotesIn += notes.length;
        sourceOf.push(parts.length - 1);
        sources.push({ track: track.index, channel: ch, name, program, baseRole: role,
          notes: notes.map((n) => ({ t: sec(n.onTick), off: sec(n.offTick), key: n.key, vel: n.onVelocity, level: levelAt(ch, sec(n.onTick)) })) });
      }
    }
  }

  function drumPart(notes: PairedNote[], part: Part): void {
    report.drumNotesIn += notes.length;
    const written = notes.map((n) => n.onVelocity).sort((a, b) => a - b);
    // all the same velocity (machine velocities): strokes at a drummer's ordinary weights instead
    const flat = written.length > 4 && written[Math.floor(written.length * 0.9)]! - written[Math.floor(written.length * 0.1)]! <= 2;
    let silenced = 0;
    for (const n of notes) {
      const rule = DRUM_MAP[n.key];
      if (!rule) {
        const label = LEFT_OUT[n.key] ?? GM_PERCUSSION[n.key] ?? `Note ${n.key}`;
        report.drumsSkipped[label] = (report.drumsSkipped[label] ?? 0) + 1;
        continue;
      }
      const t = sec(n.onTick);
      const eff = n.onVelocity * levelAt(9, t);
      if (eff < 6) { silenced++; continue; }
      const vel = flat ? (DRUM_WEIGHT[rule.piece] ?? 90) * (eff / n.onVelocity) : eff;
      if (rule.label) report.drumsMapped[rule.label] = (report.drumsMapped[rule.label] ?? 0) + 1;
      // strike near the rule's spot, never twice at exactly the same point
      const r = Math.hypot(rule.u, rule.v), a = Math.atan2(rule.v, rule.u);
      const rr = Math.max(0, Math.min(0.92, r + rand() * 0.07)), aa = a + rand() * (r > 0.05 ? 0.5 : 6);
      hits.push({ t, piece: rule.piece, u: rr * Math.cos(aa), v: rr * Math.sin(aa), g: rule.g * stickFactor(vel) * (1 + rand() * 0.04), note: n.key, vel: Math.round(vel) });
    }
    if (flat) part.actions.push(`every stroke written at velocity ${written[written.length >> 1]}: played at a drummer's ordinary weights`);
    if (silenced) part.actions.push(`${silenced} strokes silenced by the file's own volume`);
  }
  hits.sort((a, b) => a.t - b.t);
  // layered drum notes (two kick notes at once, a kick and a doubled kick) are one stroke, and one
  // piece is not struck again within 45 ms (a fast double stroke)
  {
    const lastAt = new Map<PieceId, DrumHit>();
    const kept: DrumHit[] = [];
    for (const h of hits) {
      const prev = lastAt.get(h.piece);
      if (prev && h.t - prev.t < 0.045) { prev.g = Math.max(prev.g, h.g); report.drumsLayered++; continue; }
      lastAt.set(h.piece, h);
      kept.push(h);
    }
    hits.length = 0;
    hits.push(...kept);
  }
  for (const p of parts) if (p.channel === 9 && !p.actions.length) p.actions.push("played as written");

  // ------------------------------------------------------ one pianist's part
  const orch = orchestrate(sources);
  report.silencedByVolume = orch.silenced;
  orch.plans.forEach((plan, i) => {
    const part = parts[sourceOf[i]!]!;
    part.role = plan.role;
    part.actions = plan.actions;
    part.notesOut = plan.notesOut;
    part.leftOut = plan.leftOut;
    if (plan.leftOut) report.partsLeftOut++;
  });
  const raw: (PianoNote & { prio: number })[] = orch.notes.map((n) => {
    let key = n.key;
    while (key < 21) key += 12;
    while (key > 108) key -= 12;
    if (key !== n.key) report.folded++;
    return { key, t: n.t, off: n.off, vel: n.vel, speed: hammerSpeed(n.vel), prio: n.prio };
  });

  // --------------------------------------------- merge doublings, then thin
  raw.sort((a, b) => a.t - b.t || a.key - b.key);
  const byKey = new Map<number, (PianoNote & { prio: number })[]>();
  for (const n of raw) {
    const list = byKey.get(n.key) ?? [];
    const prev = list[list.length - 1];
    if (prev && n.t - prev.t < MERGE_S) {
      prev.vel = Math.max(prev.vel, n.vel);
      prev.speed = hammerSpeed(prev.vel);
      prev.off = Math.max(prev.off, n.off);
      prev.prio = Math.max(prev.prio, n.prio);
      report.merged++;
      continue;
    }
    // a key cannot be struck again sooner than its action repeats (about 14 times a second)
    if (prev && n.t - prev.t < REPEAT_S) {
      prev.off = Math.max(prev.off, n.off);
      report.tooFast++;
      continue;
    }
    list.push(n);
    byKey.set(n.key, list);
  }
  let notes = [...byKey.values()].flat().sort((a, b) => a.t - b.t || a.key - b.key);
  // chords wider than ten fingers: keep the bass, the top and the most important, loudest notes
  const keep: typeof notes = [];
  for (let i = 0; i < notes.length;) {
    let j = i;
    while (j < notes.length && notes[j]!.t - notes[i]!.t < 0.03) j++;
    const group = notes.slice(i, j);
    if (group.length > MAX_CHORD) {
      const lo = group.reduce((a, b) => (b.key < a.key ? b : a));
      const hi = group.reduce((a, b) => (b.key > a.key ? b : a));
      const rest = group.filter((n) => n !== lo && n !== hi).sort((a, b) => b.prio - a.prio || b.vel - a.vel).slice(0, MAX_CHORD - 2);
      report.thinnedChord += group.length - rest.length - 2;
      keep.push(lo, hi, ...rest);
    } else keep.push(...group);
    i = j;
  }
  notes = keep.sort((a, b) => a.t - b.t || a.key - b.key);
  // density: at most MAX_PER_QUARTER_SECOND strikes in any 250 ms; ornaments and accompaniment
  // yield before the bass and the melody, and within a role the quietest yields
  const dense: typeof notes = [];
  const window: typeof notes = [];
  const weaker = (a: { prio: number; vel: number }, b: { prio: number; vel: number }): boolean => a.prio < b.prio || (a.prio === b.prio && a.vel < b.vel);
  for (const n of notes) {
    while (window.length && n.t - window[0]!.t >= 0.25) window.shift();
    if (window.length >= MAX_PER_QUARTER_SECOND) {
      const weakest = window.reduce((a, b) => (weaker(b, a) ? b : a));
      if (weaker(weakest, n)) {
        const k = dense.lastIndexOf(weakest);
        if (k >= 0) dense.splice(k, 1);
        window.splice(window.indexOf(weakest), 1);
      } else { report.thinnedDensity++; continue; }
      report.thinnedDensity++;
    }
    window.push(n);
    dense.push(n);
  }
  notes = dense;

  // ------------------------------------------------------------- pedal
  // The pedal is down as far as any pitched channel holds it (a two-hand file pedals on one).
  const pedalRaw: { t: number; ch: number; value: number }[] = [];
  for (const [ch, list] of pedalByChannel) for (const e of list) pedalRaw.push({ t: e.t, ch, value: e.value });
  pedalRaw.sort((a, b) => a.t - b.t);
  const level = new Map<number, number>();
  const pedal: { t: number; value: number }[] = [];
  let last = 0;
  for (const e of pedalRaw) {
    level.set(e.ch, e.value);
    const v = Math.max(0, ...level.values());
    const prev = pedal[pedal.length - 1];
    const crossing = (v >= 0.5) !== (last >= 0.5);
    if (!crossing && (Math.abs(v - last) < 0.05 || (prev && e.t - prev.t < 0.012))) { report.pedalSkipped++; continue; }
    pedal.push({ t: e.t, value: Math.round(v * 100) / 100 });
    last = v;
  }
  report.pedalChanges = pedal.length;

  // --------------------------------------------------- engine events
  // one key, one set of strings: a key struck again while still held is re-struck, and released
  // only when the last of its overlapping notes ends
  const events: PianoEvent[] = [];
  const perKey = new Map<number, PianoNote[]>();
  for (const n of notes) { const l2 = perKey.get(n.key) ?? []; l2.push(n); perKey.set(n.key, l2); }
  for (const [key, list] of perKey) {
    let end = -Infinity;
    for (const n of list) {
      if (n.t >= end && end > -Infinity) events.push({ t: end, off: key });
      events.push({ t: n.t, on: key, speed: n.speed });
      end = n.t >= end ? n.off : Math.max(end, n.off);
    }
    if (end > -Infinity) events.push({ t: end, off: key });
  }
  for (const p of pedal) events.push({ t: p.t, sustain: p.value });
  // releases before strikes at the same instant, so a repeated key is struck after it is let go
  events.sort((a, b) => a.t - b.t || rank(a) - rank(b));

  // peak keys sounding: held, or kept ringing by the pedal until it lifts (a lower bound on courses)
  report.peakKeys = peakSounding(notes, pedal);

  // ------------------------------------------------------------ labels
  const markers = midi.tracks.flatMap((t) => t.texts.filter((x) => x.kind === "marker").map((x) => ({ t: sec(x.tick), text: x.text })));
  const chordish = /^[A-G][#b♯♭]?(?:[a-zA-Z0-9#b♯♭+°ø()/,.\-]*)$/;
  const chordMarkers = markers.filter((m) => chordish.test(m.text) && m.text.length <= 16);
  const chords = chordMarkers.length >= 2
    ? chordMarkers.map((m) => ({ t: m.t, label: m.text.replace(/b(?=\d|$|\/|m|maj|sus|dim|aug)/g, "♭").replace(/#/g, "♯") }))
    : chordTimeline(notes.map((n) => ({ t: n.t, key: n.key, off: n.off })));
  const sections = markers.filter((m) => !chordMarkers.includes(m)).slice(0, 64);

  // bar starts, one meter segment at a time (a time signature always opens a new bar, as in M0)
  const lastTick = midi.tracks.reduce((m, t) => Math.max(m, t.notes.reduce((x, n) => Math.max(x, n.offTick), 0)), 0);
  const meters = midi.meterMap.length && midi.meterMap[0]!.tick === 0 ? midi.meterMap : [{ tick: 0, numerator: 4, denominatorPower: 2 }, ...midi.meterMap];
  const bars: { t: number; n: number }[] = [];
  for (let i = 0; i < meters.length && bars.length < 20000; i++) {
    const m = meters[i]!;
    const barTicks = (m.numerator * midi.division * 4) / 2 ** m.denominatorPower;
    const end = meters[i + 1]?.tick ?? lastTick + 1;
    for (let tick = m.tick; tick < end && bars.length < 20000; tick += Math.max(1, barTicks)) bars.push({ t: sec(tick), n: bars.length });
  }

  const lastNote = notes.reduce((m, n) => Math.max(m, n.off), 0);
  const lastHit = hits.length ? hits[hits.length - 1]!.t : 0;
  const lastEventS = Math.max(lastNote, lastHit, pedal.length ? pedal[pedal.length - 1]!.t : 0);
  const usedPieces = [...new Set(hits.map((h) => h.piece))];
  const order: PieceId[] = ["kick", "snare", "snare-off", "hihat-closed", "hihat-open", "tom-12", "tom-14", "ride", "crash"];
  const tempo0 = midi.tempoMap.find((e) => e.tick <= 0) ?? midi.tempoMap[0];
  const meter0 = midi.meterMap.reduce<{ m: (typeof midi.meterMap)[number] | null; span: number }>((best, m, i) => {
    const span = (midi.meterMap[i + 1]?.tick ?? lastTick) - m.tick;
    return span > best.span ? { m, span } : best;
  }, { m: null, span: -1 }).m;
  // a song title is often the name of the conductor track
  const named = midi.tracks.find((t) => t.notes.length === 0 && t.name && !LABEL.test(t.name));
  const title = named?.name ?? fileName.replace(/\.(mid|midi|kar|rmi)$/i, "");
  return {
    title,
    durationS: lastEventS + (usedPieces.some((p) => p === "crash" || p === "ride") ? 2.5 : 1.8),
    lastEventS,
    notes,
    events,
    pedal,
    hits,
    pieces: order.filter((p) => usedPieces.includes(p)),
    chords,
    sections,
    bars,
    parts,
    tempoBpm: tempo0 ? 60e6 / tempo0.microsecondsPerQuarter : 120,
    meter: meter0 ? `${meter0.numerator}/${2 ** meter0.denominatorPower}` : "4/4",
    report,
  };
}

const rank = (e: PianoEvent): number => (e.off !== undefined ? 0 : e.sustain !== undefined ? 1 : 2);

function programOf(track: DecodedTrack, channel: number): number | null {
  const p = track.programs.find((x) => x.channel === channel);
  return p ? p.program : null;
}

function peakSounding(notes: PianoNote[], pedal: { t: number; value: number }[]): number {
  // a note rings until its key is up and the pedal is up (an upper bound on real ringing)
  const ups = pedal.filter((p) => p.value < 0.5).map((p) => p.t);
  const times = pedal.map((p) => p.t);
  const after = (arr: number[], t: number): number => { // index of the first entry > t
    let lo = 0, hi = arr.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m]! <= t) lo = m + 1; else hi = m; }
    return lo;
  };
  const edges: [number, number][] = [];
  for (const n of notes) {
    let end = n.off;
    const k = after(times, n.off) - 1;
    if (k >= 0 && pedal[k]!.value >= 0.5) end = ups[after(ups, n.off)] ?? n.off + 4;
    edges.push([n.t, 1], [Math.min(end, n.t + 6), -1]);
  }
  edges.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  let cur = 0, peak = 0;
  for (const [, d] of edges) { cur += d; peak = Math.max(peak, cur); }
  return peak;
}

/** General MIDI program families (program >> 3). */
const GM_PROGRAMS = ["Piano", "Chromatic percussion", "Organ", "Guitar", "Bass", "Strings", "Ensemble", "Brass", "Reed", "Pipe", "Synth lead", "Synth pad", "Synth effects", "Ethnic", "Percussive", "Sound effects"];
export const gmFamily = (program: number | null): string | null => (program === null ? null : GM_PROGRAMS[program >> 3] ?? null);
/** General MIDI program names (0-127). */
const GM_INSTRUMENTS = ["Acoustic grand", "Bright piano", "Electric grand", "Honky-tonk", "Electric piano 1", "Electric piano 2", "Harpsichord", "Clavinet", "Celesta", "Glockenspiel", "Music box", "Vibraphone", "Marimba", "Xylophone", "Tubular bells", "Dulcimer", "Drawbar organ", "Percussive organ", "Rock organ", "Church organ", "Reed organ", "Accordion", "Harmonica", "Tango accordion", "Nylon guitar", "Steel guitar", "Jazz guitar", "Clean guitar", "Muted guitar", "Overdrive guitar", "Distortion guitar", "Guitar harmonics", "Acoustic bass", "Finger bass", "Pick bass", "Fretless bass", "Slap bass 1", "Slap bass 2", "Synth bass 1", "Synth bass 2", "Violin", "Viola", "Cello", "Contrabass", "Tremolo strings", "Pizzicato strings", "Harp", "Timpani", "Strings", "Slow strings", "Synth strings 1", "Synth strings 2", "Choir aahs", "Voice oohs", "Synth voice", "Orchestra hit", "Trumpet", "Trombone", "Tuba", "Muted trumpet", "French horn", "Brass section", "Synth brass 1", "Synth brass 2", "Soprano sax", "Alto sax", "Tenor sax", "Baritone sax", "Oboe", "English horn", "Bassoon", "Clarinet", "Piccolo", "Flute", "Recorder", "Pan flute", "Blown bottle", "Shakuhachi", "Whistle", "Ocarina", "Square lead", "Saw lead", "Calliope lead", "Chiff lead", "Charang lead", "Voice lead", "Fifths lead", "Bass and lead", "New age pad", "Warm pad", "Polysynth pad", "Choir pad", "Bowed pad", "Metallic pad", "Halo pad", "Sweep pad", "Rain", "Soundtrack", "Crystal", "Atmosphere", "Brightness", "Goblins", "Echoes", "Sci-fi", "Sitar", "Banjo", "Shamisen", "Koto", "Kalimba", "Bagpipe", "Fiddle", "Shanai", "Tinkle bell", "Agogo", "Steel drums", "Woodblock", "Taiko drum", "Melodic tom", "Synth drum", "Reverse cymbal", "Guitar fret noise", "Breath noise", "Seashore", "Bird tweet", "Telephone ring", "Helicopter", "Applause", "Gunshot"];
export const gmInstrument = (program: number | null): string | null => (program === null ? null : GM_INSTRUMENTS[program] ?? null);
