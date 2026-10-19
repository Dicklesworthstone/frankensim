import { test } from "bun:test";
import assert from "node:assert/strict";
import { pianoEventsFrom, type PianoEvent } from "../src/midi/arrange";

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
