// Provenance: copied from jazz_chord_progression_editor_html (same owner) at commit 0c3c7aa,
// src/export/midi-salvage-contract.ts. Only change: the refusal-code type below is local instead
// of imported from that project's M0 import contract.

/** The strict M0 reader's content-level refusal codes (the only ones salvage can answer). */
export type MidiImportRefusalCode = "smf.note_overlap" | "smf.note_off_unmatched" | "smf.note_on_unterminated" | (string & {});

/**
 * MIDI import salvage: the deterministic repair pass behind the strict M0
 * decoder (jcpe-v2r-midi-salvage-ulus).
 *
 * The strict decoder is total and refuses hostile bytes with structured
 * codes; that law is untouched. What this contract adds is a second chance
 * for files whose STRUCTURE is sound but whose NOTE STREAM carries the
 * defects real exports produce (a re-struck note that never ended, an
 * orphan note-off, a note left sounding at end of track). The salvage pass
 * rewrites the raw bytes with the least-disruptive musical repair for each
 * defect, then feeds the repaired bytes through the SAME strict decoder —
 * so every structural guarantee, limit, and naming law still comes from
 * the one reviewed reader, and a repaired read is stated honestly rather
 * than passed off as the file's own content.
 */

export const MIDI_SALVAGE_ENGINE_VERSION = "midi-salvage@1" as const;

/** The closed repair vocabulary, least disruptive first. */
export const MIDI_SALVAGE_REPAIR_KINDS = Object.freeze([
  /** A Note On for an already-sounding key ends the prior instance there. */
  "restruck-note-ended",
  /** A Note Off with nothing sounding on its key is dropped. */
  "orphan-off-dropped",
  /** A note still sounding at End of Track is closed at that tick. */
  "unterminated-note-closed",
] as const);
export type MidiSalvageRepairKind =
  (typeof MIDI_SALVAGE_REPAIR_KINDS)[number];

/**
 * Refusal-code classification. Only CONTENT-level codes are salvageable:
 * the note-pairing trio, whose defects live entirely inside an otherwise
 * well-formed event stream. Everything else is structural (the bytes
 * cannot be walked safely), environmental (limits), or semantic in a way
 * a byte rewrite must not invent an answer to (zero tempo, malformed
 * meter) — those refusals stand unchanged.
 */
export const MIDI_SALVAGE_CONTENT_REFUSAL_CODES = Object.freeze([
  "smf.note_overlap",
  "smf.note_off_unmatched",
  "smf.note_on_unterminated",
] as const) satisfies readonly MidiImportRefusalCode[];

export function isSalvageableRefusalCode(
  code: MidiImportRefusalCode,
): boolean {
  return (MIDI_SALVAGE_CONTENT_REFUSAL_CODES as readonly string[]).includes(
    code,
  );
}

/** Bound on total repairs; a stream needing more is not trustworthy. */
export const MAX_MIDI_SALVAGE_REPAIRS = 4_096;

export type MidiSalvageRepairRecord = Readonly<{
  kind: MidiSalvageRepairKind;
  count: number;
  /** File byte offset where this kind was first repaired. */
  firstByteOffset: number;
}>;

export type MidiSalvageWorkEvidence = Readonly<{
  bytesRead: number;
  tracksExamined: number;
  eventsExamined: number;
  termination: "complete";
}>;

export type MidiSalvageReport = Readonly<{
  engineVersion: typeof MIDI_SALVAGE_ENGINE_VERSION;
  /** One record per repair kind that fired, in vocabulary order. */
  repairs: readonly MidiSalvageRepairRecord[];
  totalRepairs: number;
  /** The honest one-sentence account for the preview surface. */
  note: string;
  evidence: MidiSalvageWorkEvidence;
}>;

export type MidiSalvageOutcome =
  | Readonly<{ salvaged: true; bytes: Uint8Array; report: MidiSalvageReport }>
  | Readonly<{
      salvaged: false;
      reason:
        /** The walker met bytes it could not parse; the original refusal stands. */
        | "unreadable"
        /** The stream needed more than MAX_MIDI_SALVAGE_REPAIRS repairs. */
        | "repairs-exceeded"
        /** The walk finished but found nothing to repair. */
        | "nothing-to-repair";
    }>;
