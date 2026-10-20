import { test } from "bun:test";
import assert from "node:assert/strict";
import { arrange, pianoEventsFrom, type PianoEvent } from "../src/midi/arrange";
import type { DecodedMidi } from "../src/midi/smf";

test("piano restart restores held keys and pedal without inventing boundary strikes", () => {
  const events: PianoEvent[] = [
    { t: 0, sustain: 0.7 }, { t: 0, on: 60, speed: 1 },
    { t: 1, on: 60, speed: 2 }, { t: 1, on: 64, speed: 1 },
    { t: 2, off: 64 }, { t: 3, off: 60 }, { t: 3, sustain: 0 },
    { t: 3, on: 67, speed: 3 }, { t: 4, off: 67 },
  ];
  const original = structuredClone(events);
  assert.deepEqual(pianoEventsFrom([], 0), []);
  assert.deepEqual(pianoEventsFrom(events, 0), events, "beginning is unchanged");
  assert.deepEqual(pianoEventsFrom(events, 1), [
    { t: 1, sustain: 0.7 }, ...events.slice(2),
  ], "an exact re-strike occurs only once");
  assert.deepEqual(pianoEventsFrom(events, 2), [
    { t: 2, sustain: 0.7 }, events[4], { t: 2, on: 60, speed: 2 }, ...events.slice(5),
  ], "released key is not re-struck; overlapping held key uses latest speed");
  assert.deepEqual(pianoEventsFrom(events, 2.5), [
    { t: 2.5, sustain: 0.7 }, { t: 2.5, on: 60, speed: 2 }, ...events.slice(5),
  ], "a released note ringing under the pedal is not a held key");
  assert.deepEqual(pianoEventsFrom(events, 3), events.slice(5));
  assert.deepEqual(pianoEventsFrom(events, 5), []);
  assert.deepEqual(pianoEventsFrom([
    { t: 0, on: 60, speed: 1 }, { t: 1, sustain: 0.2 }, { t: 2, off: 60 },
  ], 1), [
    { t: 1, sustain: 0.2 }, { t: 1, on: 60, speed: 1 }, { t: 2, off: 60 },
  ], "an exact-boundary half pedal is set before the restored strike");
  assert.deepEqual(events, original, "repeated seeks do not mutate the arrangement");
});

test("MIDI soft pedal uses the existing una-corda model and combines pitched channels", () => {
  const controls = [
    { tick: 0, channel: 0, controller: 67, value: 63 },
    { tick: 480, channel: 0, controller: 67, value: 64 },
    { tick: 720, channel: 1, controller: 67, value: 127 },
    { tick: 960, channel: 0, controller: 67, value: 63 },
    { tick: 1200, channel: 1, controller: 67, value: 63 },
    { tick: 1440, channel: 0, controller: 67, value: 0 },
    { tick: 0, channel: 9, controller: 67, value: 127 },
    { tick: 0, channel: 0, controller: 1, value: 127 },
  ];
  const midi: DecodedMidi = {
    format: 1, division: 480,
    tempoMap: [{ tick: 0, microsecondsPerQuarter: 500_000 }], meterMap: [],
    tracks: [{ index: 0, name: "Piano", instrumentName: null, programs: [], texts: [], controls,
      notes: [{ channel: 0, key: 60, onTick: 480, offTick: 1920, onVelocity: 80 }],
    }],
    salvageNote: null,
    lenient: { missingEndOfTrack: 0, truncatedChunks: 0, openNotesClosed: 0, zeroLengthNotes: 0 },
    counters: { bytes: 0, events: controls.length + 2, notes: 1, sysex: 0, ignoredMeta: 0, pitchBends: 0 },
  };
  const original = structuredClone(controls);
  const a = arrange(midi, "soft-pedal.mid");
  assert.deepEqual(a.events.filter(e => e.unaCorda !== undefined), [
    { t: 0.5, unaCorda: true }, { t: 1.25, unaCorda: false },
  ]);
  assert.equal(a.events[0]?.unaCorda, true, "the pedal is set before a simultaneous strike");
  assert.equal(a.report.softPedalChanges, 2);
  assert.equal(a.report.ignoredControllers, 2, "drum soft pedal and unsupported modulation stay explicitly unused");
  assert.deepEqual(midi.tracks[0]!.controls, original, "arranging does not rewrite the source controls");

  const longHeld = structuredClone(midi);
  longHeld.tracks[0]!.notes = [
    { channel: 0, key: 60, onTick: 0, offTick: 11520, onVelocity: 80 },
    { channel: 0, key: 64, onTick: 960, offTick: 11520, onVelocity: 80 },
    { channel: 0, key: 67, onTick: 7680, offTick: 11520, onVelocity: 80 },
  ];
  assert.equal(arrange(longHeld, "held-chord.mid").report.peakKeys, 3,
    "a six-second ringing estimate must not drop keys still held down");
});

test("MIDI sostenuto is a switch shared by pitched channels and precedes simultaneous strikes", () => {
  const midi: DecodedMidi = {
    format: 1, division: 480,
    tempoMap: [{ tick: 0, microsecondsPerQuarter: 500_000 }], meterMap: [],
    tracks: [{ index: 0, name: "Piano", instrumentName: null, programs: [], texts: [],
      notes: [{ channel: 0, key: 60, onTick: 480, offTick: 1920, onVelocity: 80 }],
      controls: [
        { tick: 0, channel: 0, controller: 66, value: 63 },
        { tick: 480, channel: 0, controller: 66, value: 64 },
        { tick: 720, channel: 1, controller: 66, value: 127 },
        { tick: 960, channel: 0, controller: 66, value: 63 },
        { tick: 1200, channel: 1, controller: 66, value: 63 },
        { tick: 0, channel: 9, controller: 66, value: 127 },
      ],
    }], salvageNote: null,
    lenient: { missingEndOfTrack: 0, truncatedChunks: 0, openNotesClosed: 0, zeroLengthNotes: 0 },
    counters: { bytes: 0, events: 8, notes: 1, sysex: 0, ignoredMeta: 0, pitchBends: 0 },
  };
  const a = arrange(midi, "sostenuto.mid");
  assert.deepEqual(a.events.filter(e => e.sostenuto !== undefined), [
    { t: 0.5, sostenuto: true }, { t: 1.25, sostenuto: false },
  ]);
  assert.equal(a.events[0]?.sostenuto, true);
  assert.equal(a.report.sostenutoChanges, 2);
  assert.equal(a.report.ignoredControllers, 1);
  const handover = structuredClone(midi);
  handover.tracks[0]!.controls = [
    { tick: 480, channel: 0, controller: 66, value: 127 },
    { tick: 960, channel: 0, controller: 66, value: 0 },
    { tick: 960, channel: 1, controller: 66, value: 127 },
    { tick: 1440, channel: 1, controller: 66, value: 0 },
  ];
  assert.deepEqual(arrange(handover, "handover.mid").events.filter(e => e.sostenuto !== undefined), [
    { t: 0.5, sostenuto: true }, { t: 1.5, sostenuto: false },
  ], "a same-frame handover must not release and recapture the shared dampers");
  const tails = structuredClone(midi);
  tails.tracks[0]!.notes = [
    { channel: 0, key: 48, onTick: 0, offTick: 1920, onVelocity: 80 },
    { channel: 0, key: 64, onTick: 2880, offTick: 3840, onVelocity: 80 },
  ];
  tails.tracks[0]!.controls = [
    { tick: 960, channel: 0, controller: 66, value: 127 },
    { tick: 3360, channel: 0, controller: 66, value: 0 },
  ];
  assert.equal(arrange(tails, "captured-bass.mid").report.peakKeys, 2,
    "the overlap estimate includes a released key whose damper remains captured");
});

test("sustain smoothing preserves full travel and release before sostenuto capture", () => {
  const midi: DecodedMidi = {
    format: 1, division: 480,
    tempoMap: [{ tick: 0, microsecondsPerQuarter: 500_000 }], meterMap: [],
    tracks: [{ index: 0, name: "Piano", instrumentName: null, programs: [], texts: [],
      notes: [{ channel: 0, key: 60, onTick: 0, offTick: 1920, onVelocity: 80 }],
      controls: [
        { tick: 0, channel: 0, controller: 64, value: 124 },
        { tick: 1, channel: 0, controller: 64, value: 127 },
        { tick: 2, channel: 0, controller: 66, value: 127 },
        { tick: 100, channel: 0, controller: 64, value: 3 },
        { tick: 101, channel: 0, controller: 64, value: 0 },
      ],
    }], salvageNote: null,
    lenient: { missingEndOfTrack: 0, truncatedChunks: 0, openNotesClosed: 0, zeroLengthNotes: 0 },
    counters: { bytes: 0, events: 7, notes: 1, sysex: 0, ignoredMeta: 0, pitchBends: 0 },
  };
  const a = arrange(midi, "pedal-endpoints.mid");
  assert.deepEqual(a.events.filter(e => e.sustain !== undefined), [
    { t: 0, sustain: 0.98 }, { t: 1 / 960, sustain: 1 },
    { t: 100 / 960, sustain: 0.02 }, { t: 101 / 960, sustain: 0 },
  ], "small or rapid movements must still reach the physical endpoints");
  assert.equal(a.report.pedalChanges, 4);
  assert.deepEqual(pianoEventsFrom(a.events, 0.2).slice(0, 3), [
    { t: 0.2, sustain: 1 }, { t: 0.2, sostenuto: true }, { t: 0.2, sustain: 0 },
  ], "seeking restores the all-damper capture made at full sustain");

  midi.tracks[0]!.controls = [
    { tick: 0, channel: 0, controller: 64, value: 127 },
    { tick: 1, channel: 0, controller: 64, value: 124 },
    { tick: 2, channel: 0, controller: 66, value: 127 },
  ];
  assert.deepEqual(arrange(midi, "pedal-partial.mid").events.filter(e => e.sustain !== undefined), [
    { t: 0, sustain: 1 }, { t: 1 / 960, sustain: 0.98 },
  ], "leaving full travel must not keep all dampers raised at a later capture");

  midi.tracks[0]!.controls = [
    { tick: 0, channel: 0, controller: 64, value: 127 },
    { tick: 480, channel: 0, controller: 64, value: 0 },
    { tick: 480, channel: 1, controller: 64, value: 127 },
    { tick: 960, channel: 1, controller: 64, value: 0 },
  ];
  const handover = arrange(midi, "sustain-handover.mid");
  assert.deepEqual(handover.events.filter(e => e.sustain !== undefined), [
    { t: 0, sustain: 1 }, { t: 1, sustain: 0 },
  ], "simultaneous channel handover must not lift the shared pedal");
  assert.equal(handover.report.pedalChanges + handover.report.pedalSkipped, 4);
  midi.tracks[0]!.controls.reverse();
  assert.deepEqual(arrange(midi, "reordered-handover.mid").events, handover.events,
    "independent channel requests at one instant have no track-order dependency");
});

test("sostenuto restart restores captured dampers without striking released or later notes", () => {
  const events: PianoEvent[] = [
    { t: 0, on: 48, speed: 2 }, { t: 1, sostenuto: true },
    { t: 2, off: 48 }, { t: 3, on: 64, speed: 1 },
    { t: 5, sostenuto: false }, { t: 6, off: 64 },
  ];
  const original = structuredClone(events);
  assert.deepEqual(pianoEventsFrom(events, 4), [
    { t: 4, silent: 48 }, { t: 4, sostenuto: true }, { t: 4, off: 48 },
    { t: 4, on: 64, speed: 1 }, ...events.slice(4),
  ], "the released bass is silently recaptured before the later melody is struck");
  assert.deepEqual(pianoEventsFrom(events, 5), [
    events[4], { t: 5, on: 64, speed: 1 }, events[5],
  ], "an exact-boundary pedal release does not recreate obsolete latches");
  assert.deepEqual(pianoEventsFrom(events, 1), [
    { t: 1, silent: 48 }, events[1], { t: 1, on: 48, speed: 2 }, ...events.slice(2),
  ], "a pedal first pressed at the restart catches the already-held key");
  assert.deepEqual(events, original);
  assert.deepEqual(pianoEventsFrom([{ t: 0, silent: 48 }, { t: 2, off: 48 }], 1), [
    { t: 1, silent: 48 }, { t: 2, off: 48 },
  ], "a held silent key stays silent when a restored score is re-cued");
});

test("sostenuto restart preserves a full-sustain capture without activating 88 silent keys", () => {
  const events: PianoEvent[] = [
    { t: 0, sustain: 1 }, { t: 1, sostenuto: true }, { t: 2, sustain: 0.2 },
    { t: 5, sostenuto: false },
  ];
  assert.deepEqual(pianoEventsFrom(events, 3), [
    { t: 3, sustain: 1 }, { t: 3, sostenuto: true }, { t: 3, sustain: 0 },
    { t: 3, sustain: 0.2 }, events[3],
  ]);
  const subset: PianoEvent[] = [
    { t: 0, on: 48, speed: 2 }, { t: 1, sostenuto: true },
    { t: 2, sustain: 1 }, { t: 5, sostenuto: false },
  ];
  assert.deepEqual(pianoEventsFrom(subset, 3), [
    { t: 3, silent: 48 }, { t: 3, sostenuto: true }, { t: 3, off: 48 },
    { t: 3, sustain: 1 }, { t: 3, on: 48, speed: 2 }, subset[3],
  ], "sustain pressed later must not expand the old captured set");
});

test("piano restart restores soft pedal before held strikes and respects exact-boundary release", () => {
  const events: PianoEvent[] = [
    { t: 0, unaCorda: true }, { t: 0, on: 60, speed: 2 },
    { t: 2, unaCorda: false }, { t: 3, off: 60 },
  ];
  assert.deepEqual(pianoEventsFrom(events, 1), [
    { t: 1, unaCorda: true }, { t: 1, on: 60, speed: 2 }, ...events.slice(2),
  ]);
  assert.deepEqual(pianoEventsFrom(events, 2), [
    events[2], { t: 2, on: 60, speed: 2 }, events[3],
  ]);
  assert.deepEqual(pianoEventsFrom(events, 2.5), [
    { t: 2.5, on: 60, speed: 2 }, events[3],
  ], "the released soft pedal must not be restored on a fresh engine");
});
