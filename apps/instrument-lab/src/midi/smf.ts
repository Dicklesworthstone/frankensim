// Standard MIDI File reader for the player.
//
// Provenance: follows the M0 import laws of jazz_chord_progression_editor_html (same owner,
// commit 0c3c7aa: src/export/midi-import-contract.ts and docs/M1_MIDI_IMPORT_AUTOMATION_CONTRACT.md)
// and runs every file through that project's salvage pass (./midi-salvage.ts) first. Its decode
// model has the same shape as M0's (header, tempo and meter maps, tracks of paired notes sorted by
// (onTick, channel, key)), extended with what a player needs and an editor discards: the sustain
// pedal and other controllers, program changes and every channel's notes. The M0 byte layer is
// Rust compiled into that project's single embedded wasm; this is a TypeScript reader of the same
// grammar, so the site ships no second decoder module.
//
// A player is more forgiving than an editor: a missing End of Track, a chunk that runs past the
// end of the file, or a note still open at the end are read anyway, and every such leniency is
// counted in `lenient` so the page can say what it did. Format 2 and SMPTE time refuse.

import { salvageMidiBytes } from "./midi-salvage";

/** M0 limits (midi-import-contract.ts). */
export const MAX_MIDI_BYTES = 4_194_304;
const MAX_TRACKS = 64;
const MAX_EVENTS = 524_288;
const MAX_NOTES = 131_072;

export interface TempoEntry { tick: number; microsecondsPerQuarter: number }
export interface MeterEntry { tick: number; numerator: number; denominatorPower: number }
export interface PairedNote { channel: number; key: number; onTick: number; offTick: number; onVelocity: number }
export interface ControlEvent { tick: number; channel: number; controller: number; value: number }
export interface ProgramEvent { tick: number; channel: number; program: number }
export interface TextEvent { tick: number; kind: "marker" | "text" | "lyric" | "cue"; text: string }
export interface DecodedTrack {
  index: number;
  name: string | null;
  instrumentName: string | null;
  notes: PairedNote[];
  controls: ControlEvent[];
  programs: ProgramEvent[];
  texts: TextEvent[];
}
export interface DecodedMidi {
  format: 0 | 1;
  division: number;
  tempoMap: TempoEntry[];
  meterMap: MeterEntry[];
  tracks: DecodedTrack[];
  /** Salvage repairs (restruck notes ended, orphan offs dropped, open notes closed), if any. */
  salvageNote: string | null;
  lenient: { missingEndOfTrack: number; truncatedChunks: number; openNotesClosed: number; zeroLengthNotes: number };
  counters: { bytes: number; events: number; notes: number; sysex: number; ignoredMeta: number; pitchBends: number };
}
export type DecodeResult = { ok: true; midi: DecodedMidi } | { ok: false; message: string };

const strictUtf8 = new TextDecoder("utf-8", { fatal: true });
const latin1 = new TextDecoder("latin1");
/** Meta text is UTF-8 when it decodes as such, otherwise Latin-1 (older files). */
function text(bytes: Uint8Array): string {
  try {
    return strictUtf8.decode(bytes).replace(/\0+$/, "").trim();
  } catch {
    return latin1.decode(bytes).replace(/\0+$/, "").trim();
  }
}

function u32(b: Uint8Array, o: number): number {
  return (((b[o] ?? 0) << 24) | ((b[o + 1] ?? 0) << 16) | ((b[o + 2] ?? 0) << 8) | (b[o + 3] ?? 0)) >>> 0;
}
function tag(b: Uint8Array, o: number): string {
  return String.fromCharCode(b[o] ?? 0, b[o + 1] ?? 0, b[o + 2] ?? 0, b[o + 3] ?? 0);
}
/** M0: a variable-length quantity has at most four bytes. */
function vlq(b: Uint8Array, at: number, end: number): { value: number; next: number } | null {
  let value = 0;
  let p = at;
  for (let i = 0; i < 4; i++) {
    if (p >= end) return null;
    const byte = b[p++]!;
    value = value * 128 + (byte & 0x7f);
    if ((byte & 0x80) === 0) return { value, next: p };
  }
  return null;
}

/** Decode a Standard MIDI File (after the salvage pass). Never throws. */
export function decodeMidi(input: Uint8Array): DecodeResult {
  try {
    return decode(input);
  } catch (e) {
    return { ok: false, message: `This file could not be read (${(e as Error).message}).` };
  }
}

function decode(input: Uint8Array): DecodeResult {
  if (input.byteLength > MAX_MIDI_BYTES) return { ok: false, message: "MIDI files of up to 4 MiB are accepted; this one is larger." };
  // RIFF-wrapped MIDI (.rmi): the SMF is the 'data' chunk
  let bytes = input;
  if (tag(bytes, 0) === "RIFF" && tag(bytes, 8) === "RMID") {
    let o = 12;
    while (o + 8 <= bytes.byteLength && tag(bytes, o) !== "data") o += 8 + ((bytes[o + 4]! | (bytes[o + 5]! << 8) | (bytes[o + 6]! << 16) | (bytes[o + 7]! << 24)) >>> 0);
    if (o + 8 <= bytes.byteLength) bytes = bytes.subarray(o + 8);
  }
  if (bytes.byteLength < 14 || tag(bytes, 0) !== "MThd") return { ok: false, message: "This is not a Standard MIDI File (no MThd header)." };
  let salvageNote: string | null = null;
  const salvaged = salvageMidiBytes(bytes);
  if (salvaged.salvaged) {
    bytes = salvaged.bytes;
    salvageNote = salvaged.report.note.replace(/ — chord names are a best guess\.$/, ".");
  }
  const headerLength = u32(bytes, 4);
  if (headerLength < 6) return { ok: false, message: "The MIDI header is malformed." };
  const format = ((bytes[8] ?? 0) << 8) | (bytes[9] ?? 0);
  const declaredTracks = ((bytes[10] ?? 0) << 8) | (bytes[11] ?? 0);
  const division = ((bytes[12] ?? 0) << 8) | (bytes[13] ?? 0);
  if (format === 2) return { ok: false, message: "Format 2 MIDI files (independent patterns) are not supported." };
  if (format > 2) return { ok: false, message: `Unknown MIDI format ${format}.` };
  if (division & 0x8000) return { ok: false, message: "This file counts time in SMPTE frames, which the player does not read." };
  if (division < 1) return { ok: false, message: "The file declares zero ticks per quarter note." };

  const lenient = { missingEndOfTrack: 0, truncatedChunks: 0, openNotesClosed: 0, zeroLengthNotes: 0 };
  const counters = { bytes: bytes.byteLength, events: 0, notes: 0, sysex: 0, ignoredMeta: 0, pitchBends: 0 };
  const tempoMap: TempoEntry[] = [];
  const meterMap: MeterEntry[] = [];
  const tracks: DecodedTrack[] = [];

  let cursor = 8 + headerLength;
  while (cursor + 8 <= bytes.byteLength && tracks.length < Math.max(declaredTracks, 1) + 8) {
    const t = tag(bytes, cursor);
    let length = u32(bytes, cursor + 4);
    const start = cursor + 8;
    if (start + length > bytes.byteLength) { length = bytes.byteLength - start; lenient.truncatedChunks++; }
    const end = start + length;
    cursor = end;
    if (t !== "MTrk") continue; // alien chunks are tolerated, as in M0
    if (tracks.length >= MAX_TRACKS) return { ok: false, message: `More than ${MAX_TRACKS} tracks.` };
    const track: DecodedTrack = { index: tracks.length, name: null, instrumentName: null, notes: [], controls: [], programs: [], texts: [] };
    const open = new Map<number, { tick: number; velocity: number }>(); // channel*128+key
    let tick = 0;
    let running = -1;
    let p = start;
    let ended = false;
    while (p < end) {
      const d = vlq(bytes, p, end);
      if (!d) break;
      p = d.next;
      tick += d.value;
      let status = bytes[p];
      if (status === undefined) break;
      if (status < 0x80) {
        if (running < 0x80) break; // data byte without a running status: stop reading this track
        status = running;
      } else p++;
      if (++counters.events > MAX_EVENTS) return { ok: false, message: "The file holds more events than the player accepts (524,288)." };
      if (status === 0xff) {
        const type = bytes[p++];
        const len = vlq(bytes, p, end);
        if (type === undefined || !len) break;
        const payload = bytes.subarray(len.next, Math.min(end, len.next + len.value));
        p = len.next + len.value;
        running = -1;
        if (type === 0x2f) { ended = true; break; }
        if (type === 0x51 && payload.length === 3) {
          const us = (payload[0]! << 16) | (payload[1]! << 8) | payload[2]!;
          if (us > 0) tempoMap.push({ tick, microsecondsPerQuarter: us });
        } else if (type === 0x58 && payload.length >= 2) {
          if (payload[0]! > 0 && payload[1]! <= 6) meterMap.push({ tick, numerator: payload[0]!, denominatorPower: payload[1]! });
        } else if (type === 0x03) {
          if (track.name === null) track.name = text(payload);
        } else if (type === 0x04) {
          if (track.instrumentName === null) track.instrumentName = text(payload);
        } else if (type === 0x06 || type === 0x01 || type === 0x05 || type === 0x07) {
          const s = text(payload);
          if (s) track.texts.push({ tick, kind: type === 0x06 ? "marker" : type === 0x05 ? "lyric" : type === 0x07 ? "cue" : "text", text: s.slice(0, 120) });
        } else counters.ignoredMeta++;
        continue;
      }
      if (status === 0xf0 || status === 0xf7) {
        const len = vlq(bytes, p, end);
        if (!len) break;
        p = len.next + len.value;
        running = -1;
        counters.sysex++;
        continue;
      }
      const high = status & 0xf0;
      if (high < 0x80 || high === 0xf0) break;
      running = status;
      const channel = status & 0x0f;
      const a = bytes[p] ?? 0;
      const two = high !== 0xc0 && high !== 0xd0;
      const b = two ? bytes[p + 1] ?? 0 : 0;
      p += two ? 2 : 1;
      if (a > 0x7f || b > 0x7f) break;
      if (high === 0x90 && b > 0) {
        const slot = channel * 128 + a;
        const prior = open.get(slot);
        if (prior) close(track, channel, a, prior, tick, lenient); // salvage normally ended it already
        open.set(slot, { tick, velocity: b });
      } else if (high === 0x80 || high === 0x90) {
        const slot = channel * 128 + a;
        const prior = open.get(slot);
        if (prior) { close(track, channel, a, prior, tick, lenient); open.delete(slot); }
      } else if (high === 0xb0) {
        track.controls.push({ tick, channel, controller: a, value: b });
      } else if (high === 0xc0) {
        track.programs.push({ tick, channel, program: a });
      } else if (high === 0xe0 && (a !== 0 || b !== 64)) {
        counters.pitchBends++; // a piano key cannot bend
      }
    }
    if (!ended) lenient.missingEndOfTrack++;
    for (const [slot, prior] of open) {
      lenient.openNotesClosed++;
      close(track, Math.floor(slot / 128), slot % 128, prior, Math.max(tick, prior.tick + 1), lenient);
    }
    counters.notes += track.notes.length;
    if (counters.notes > MAX_NOTES) return { ok: false, message: "The file holds more notes than the player accepts (131,072)." };
    track.notes.sort((x, y) => x.onTick - y.onTick || x.channel - y.channel || x.key - y.key);
    tracks.push(track);
  }
  if (!tracks.length) return { ok: false, message: "The file has no tracks." };
  tempoMap.sort((x, y) => x.tick - y.tick);
  meterMap.sort((x, y) => x.tick - y.tick);
  return { ok: true, midi: { format: format === 0 ? 0 : 1, division, tempoMap, meterMap, tracks, salvageNote, lenient, counters } };
}

function close(track: DecodedTrack, channel: number, key: number, on: { tick: number; velocity: number }, offTick: number, lenient: DecodedMidi["lenient"]): void {
  let off = offTick;
  if (off <= on.tick) { lenient.zeroLengthNotes++; off = on.tick + 1; }
  track.notes.push({ channel, key, onTick: on.tick, offTick: off, onVelocity: on.velocity });
}

/** Seconds at every tick, from the file's tempo map (120 bpm until the first tempo, as SMF says). */
export function tempoClock(midi: DecodedMidi): (tick: number) => number {
  const ppq = midi.division;
  const segs: { tick: number; sec: number; us: number }[] = [];
  let us = 500_000, sec = 0, last = 0;
  for (const e of midi.tempoMap) {
    sec += ((e.tick - last) * us) / ppq / 1e6;
    last = e.tick;
    us = e.microsecondsPerQuarter;
    segs.push({ tick: e.tick, sec, us });
  }
  return (tick: number): number => {
    let lo = 0, hi = segs.length - 1, k = -1;
    while (lo <= hi) { const mid = (lo + hi) >> 1; if (segs[mid]!.tick <= tick) { k = mid; lo = mid + 1; } else hi = mid - 1; }
    if (k < 0) return (tick * 500_000) / ppq / 1e6;
    const s = segs[k]!;
    return s.sec + ((tick - s.tick) * s.us) / ppq / 1e6;
  };
}
