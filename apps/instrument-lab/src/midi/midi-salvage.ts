// Provenance: copied unchanged from jazz_chord_progression_editor_html (same owner) at commit
// 0c3c7aa, src/export/midi-salvage.ts. The player runs every file through it before decoding.
import {
  MAX_MIDI_SALVAGE_REPAIRS,
  MIDI_SALVAGE_ENGINE_VERSION,
  MIDI_SALVAGE_REPAIR_KINDS,
  type MidiSalvageOutcome,
  type MidiSalvageRepairKind,
  type MidiSalvageReport,
} from "./midi-salvage-contract";

/**
 * The byte-level salvage walker.
 *
 * A deliberately small exact SMF reader: header chunk and alien chunks pass
 * through verbatim; each MTrk is decoded into absolute-tick events, the
 * note stream is repaired with the closed vocabulary in the contract, and
 * the track is re-emitted with explicit status bytes and recomputed deltas.
 * Everything that is not a note on/off — meta, sysex, controllers, tempo,
 * meter — is preserved byte-for-byte in its event payload. Any byte the
 * walker cannot parse aborts the salvage; the strict decoder's original
 * refusal then stands. Determinism is total: same bytes in, same bytes and
 * ledger out.
 */

type TrackEvent = Readonly<{
  tick: number;
  /** Complete event bytes with explicit status, delta excluded. */
  bytes: Uint8Array;
}>;

type RepairTally = { count: number; firstByteOffset: number };

function readVlq(
  bytes: Uint8Array,
  start: number,
  end: number,
): Readonly<{ value: number; next: number }> | null {
  let value = 0;
  let position = start;
  for (let index = 0; index < 4; index += 1) {
    if (position >= end) return null;
    const byte = bytes[position];
    if (byte === undefined) return null;
    position += 1;
    value = value * 128 + (byte & 0x7f);
    if ((byte & 0x80) === 0) return { value, next: position };
  }
  return null;
}

function vlqBytes(value: number): number[] {
  const out = [value & 0x7f];
  let rest = Math.floor(value / 128);
  while (rest > 0) {
    out.unshift((rest & 0x7f) | 0x80);
    rest = Math.floor(rest / 128);
  }
  return out;
}

function readU32(bytes: Uint8Array, offset: number): number {
  return (
    ((bytes[offset] ?? 0) << 24) |
    ((bytes[offset + 1] ?? 0) << 16) |
    ((bytes[offset + 2] ?? 0) << 8) |
    (bytes[offset + 3] ?? 0)
  ) >>> 0;
}

function tagAt(bytes: Uint8Array, offset: number): string {
  return String.fromCharCode(
    bytes[offset] ?? 0,
    bytes[offset + 1] ?? 0,
    bytes[offset + 2] ?? 0,
    bytes[offset + 3] ?? 0,
  );
}

const KIND_ORDER: readonly MidiSalvageRepairKind[] = MIDI_SALVAGE_REPAIR_KINDS;

function composeNote(
  tallies: ReadonlyMap<MidiSalvageRepairKind, RepairTally>,
  total: number,
): string {
  const parts: string[] = [];
  const restruck = tallies.get("restruck-note-ended");
  if (restruck !== undefined) {
    parts.push(
      `${String(restruck.count)} restruck ${restruck.count === 1 ? "note" : "notes"} ended early`,
    );
  }
  const orphan = tallies.get("orphan-off-dropped");
  if (orphan !== undefined) {
    parts.push(
      `${String(orphan.count)} orphan note-${orphan.count === 1 ? "off" : "offs"} dropped`,
    );
  }
  const open = tallies.get("unterminated-note-closed");
  if (open !== undefined) {
    parts.push(
      `${String(open.count)} unterminated ${open.count === 1 ? "note" : "notes"} closed`,
    );
  }
  return `Read after ${String(total)} ${total === 1 ? "repair" : "repairs"}: ${parts.join(", ")} — chord names are a best guess.`;
}

/**
 * Attempts the salvage rewrite. Callers gate on the strict decoder having
 * refused with a content-level code first; this function neither sees nor
 * interprets that refusal, it simply repairs what the walk finds.
 */
export function salvageMidiBytes(bytes: Uint8Array): MidiSalvageOutcome {
  const unreadable = Object.freeze({
    salvaged: false as const,
    reason: "unreadable" as const,
  });
  if (bytes.byteLength < 14 || tagAt(bytes, 0) !== "MThd") return unreadable;
  const headerLength = readU32(bytes, 4);
  const headerEnd = 8 + headerLength;
  if (headerLength < 6 || headerEnd > bytes.byteLength) return unreadable;

  const output: number[] = [...bytes.slice(0, headerEnd)];
  const tallies = new Map<MidiSalvageRepairKind, RepairTally>();
  let totalRepairs = 0;
  let tracksExamined = 0;
  let eventsExamined = 0;

  const repair = (kind: MidiSalvageRepairKind, byteOffset: number): boolean => {
    totalRepairs += 1;
    if (totalRepairs > MAX_MIDI_SALVAGE_REPAIRS) return false;
    const tally = tallies.get(kind);
    if (tally === undefined) {
      tallies.set(kind, { count: 1, firstByteOffset: byteOffset });
    } else {
      tally.count += 1;
    }
    return true;
  };

  let cursor = headerEnd;
  while (cursor < bytes.byteLength) {
    if (cursor + 8 > bytes.byteLength) return unreadable;
    const tag = tagAt(bytes, cursor);
    const chunkLength = readU32(bytes, cursor + 4);
    const dataStart = cursor + 8;
    const dataEnd = dataStart + chunkLength;
    if (dataEnd > bytes.byteLength) return unreadable;
    if (tag !== "MTrk") {
      /* Alien chunks pass through untouched, exactly as the reader tolerates. */
      for (let index = cursor; index < dataEnd; index += 1) {
        output.push(bytes[index] ?? 0);
      }
      cursor = dataEnd;
      continue;
    }

    tracksExamined += 1;
    const events: TrackEvent[] = [];
    /* Open-note slots: channel * 128 + key, true while sounding. */
    const open = new Set<number>();
    let tick = 0;
    let runningStatus = -1;
    let position = dataStart;
    let sawEndOfTrack = false;

    while (position < dataEnd) {
      if (sawEndOfTrack) return unreadable;
      const delta = readVlq(bytes, position, dataEnd);
      if (delta === null) return unreadable;
      position = delta.next;
      tick += delta.value;

      let status = bytes[position];
      if (status === undefined) return unreadable;
      if (status < 0x80) {
        if (runningStatus < 0x80) return unreadable;
        status = runningStatus;
      } else {
        position += 1;
      }
      eventsExamined += 1;

      if (status === 0xff) {
        const type = bytes[position];
        if (type === undefined) return unreadable;
        position += 1;
        const length = readVlq(bytes, position, dataEnd);
        if (length === null) return unreadable;
        const payloadStart = length.next;
        const payloadEnd = payloadStart + length.value;
        if (payloadEnd > dataEnd) return unreadable;
        if (type === 0x2f) {
          /* Close every note still sounding, at the End of Track tick. */
          const openSlots = [...open].sort((left, right) => left - right);
          for (const slot of openSlots) {
            if (!repair("unterminated-note-closed", position - 2)) {
              return Object.freeze({
                salvaged: false as const,
                reason: "repairs-exceeded" as const,
              });
            }
            const channel = Math.floor(slot / 128);
            const key = slot % 128;
            events.push(
              Object.freeze({
                tick,
                bytes: Uint8Array.from([0x80 | channel, key, 0]),
              }),
            );
            open.delete(slot);
          }
          sawEndOfTrack = true;
        }
        events.push(
          Object.freeze({
            tick,
            bytes: Uint8Array.from([
              0xff,
              type,
              ...vlqBytes(length.value),
              ...bytes.slice(payloadStart, payloadEnd),
            ]),
          }),
        );
        position = payloadEnd;
        runningStatus = -1;
        continue;
      }

      if (status === 0xf0 || status === 0xf7) {
        const length = readVlq(bytes, position, dataEnd);
        if (length === null) return unreadable;
        const payloadStart = length.next;
        const payloadEnd = payloadStart + length.value;
        if (payloadEnd > dataEnd) return unreadable;
        events.push(
          Object.freeze({
            tick,
            bytes: Uint8Array.from([
              status,
              ...vlqBytes(length.value),
              ...bytes.slice(payloadStart, payloadEnd),
            ]),
          }),
        );
        position = payloadEnd;
        runningStatus = -1;
        continue;
      }

      const high = status & 0xf0;
      if (high < 0x80 || high === 0xf0) return unreadable;
      runningStatus = status;
      const dataCount = high === 0xc0 || high === 0xd0 ? 1 : 2;
      const firstDataOffset = position;
      const data: number[] = [];
      for (let index = 0; index < dataCount; index += 1) {
        const byte = bytes[position];
        if (byte === undefined || byte >= 0x80 || position >= dataEnd) {
          return unreadable;
        }
        data.push(byte);
        position += 1;
      }

      if (high === 0x80 || high === 0x90) {
        const channel = status & 0x0f;
        const key = data[0] ?? 0;
        const velocity = data[1] ?? 0;
        const slot = channel * 128 + key;
        const isOn = high === 0x90 && velocity > 0;
        if (isOn) {
          if (open.has(slot)) {
            /* The re-strike implicitly ends the prior instance here. */
            if (!repair("restruck-note-ended", firstDataOffset)) {
              return Object.freeze({
                salvaged: false as const,
                reason: "repairs-exceeded" as const,
              });
            }
            events.push(
              Object.freeze({
                tick,
                bytes: Uint8Array.from([0x80 | channel, key, 0]),
              }),
            );
          }
          open.add(slot);
          events.push(
            Object.freeze({
              tick,
              bytes: Uint8Array.from([status, key, velocity]),
            }),
          );
        } else {
          if (!open.has(slot)) {
            /* Orphan off: nothing is sounding there; drop it. */
            if (!repair("orphan-off-dropped", firstDataOffset)) {
              return Object.freeze({
                salvaged: false as const,
                reason: "repairs-exceeded" as const,
              });
            }
            continue;
          }
          open.delete(slot);
          events.push(
            Object.freeze({
              tick,
              bytes: Uint8Array.from([status, key, velocity]),
            }),
          );
        }
        continue;
      }

      /* Other channel voice messages pass through with explicit status. */
      events.push(
        Object.freeze({ tick, bytes: Uint8Array.from([status, ...data]) }),
      );
    }

    if (!sawEndOfTrack) return unreadable;

    /* Re-emit the track: explicit statuses, recomputed deltas. */
    const trackBytes: number[] = [];
    let previousTick = 0;
    for (const event of events) {
      trackBytes.push(...vlqBytes(event.tick - previousTick));
      previousTick = event.tick;
      for (const byte of event.bytes) trackBytes.push(byte);
    }
    output.push(0x4d, 0x54, 0x72, 0x6b);
    output.push(
      (trackBytes.length >>> 24) & 0xff,
      (trackBytes.length >>> 16) & 0xff,
      (trackBytes.length >>> 8) & 0xff,
      trackBytes.length & 0xff,
    );
    /*
     * Byte-by-byte: a spread here overflows the call stack on a large
     * track (each spread element becomes an argument), and readFile's
     * never-throw-on-hostile-bytes contract forbids that escape.
     */
    for (const byte of trackBytes) output.push(byte);
    cursor = dataEnd;
  }

  if (totalRepairs === 0) {
    return Object.freeze({
      salvaged: false as const,
      reason: "nothing-to-repair" as const,
    });
  }

  const repairs = Object.freeze(
    KIND_ORDER.flatMap((kind) => {
      const tally = tallies.get(kind);
      return tally === undefined
        ? []
        : [
            Object.freeze({
              kind,
              count: tally.count,
              firstByteOffset: tally.firstByteOffset,
            }),
          ];
    }),
  );
  const report: MidiSalvageReport = Object.freeze({
    engineVersion: MIDI_SALVAGE_ENGINE_VERSION,
    repairs,
    totalRepairs,
    note: composeNote(tallies, totalRepairs),
    evidence: Object.freeze({
      bytesRead: bytes.byteLength,
      tracksExamined,
      eventsExamined,
      termination: "complete" as const,
    }),
  });
  return Object.freeze({
    salvaged: true as const,
    bytes: Uint8Array.from(output),
    report,
  });
}
