import { test } from "bun:test";
import assert from "node:assert/strict";
import { arrange, hiHatEventsFrom, hiHatPerformance, type HiHatEvent } from "../src/midi/arrange";
import { HIHAT_CLOSED_FORCE_N, HIHAT_TRANSITION_SECONDS } from "../src/live/hiHat";
import type { DecodedMidi } from "../src/midi/smf";

function midi(notes: [key: number, tick: number, velocity?: number][], controls: DecodedMidi["tracks"][number]["controls"] = []): DecodedMidi {
  return {
    format: 1, division: 480,
    tempoMap: [{ tick: 0, microsecondsPerQuarter: 500_000 }], meterMap: [],
    tracks: [{ index: 0, name: "Drums", instrumentName: null, programs: [], texts: [], controls,
      notes: notes.map(([key, tick, velocity = 90]) => ({ channel: 9, key, onTick: tick, offTick: tick + 48, onVelocity: velocity })),
    }], salvageNote: null,
    lenient: { missingEndOfTrack: 0, truncatedChunks: 0, openNotesClosed: 0, zeroLengthNotes: 0 },
    counters: { bytes: 0, events: 2 * notes.length + controls.length, notes: notes.length, sysex: 0, ignoredMeta: 0, pitchBends: 0 },
  };
}
const pedals = (events: readonly HiHatEvent[]) => events.filter((e) => "pedalForceN" in e);
const strokes = (events: readonly HiHatEvent[]) => events.filter((e) => "g" in e);
const pedal = (t: number, pedalForceN: number) => ({ t, pedalForceN, transitionSeconds: HIHAT_TRANSITION_SECONDS });

test("closed, open and foot notes share one hi-hat; the foot never launches a stick", () => {
  const source = midi([[42, 0], [46, 960], [44, 1920], [38, 1920]]);
  const original = structuredClone(source);
  const a = arrange(source, "shared-hat.mid");
  assert.deepEqual(a.pieces, ["snare", "hihat"]);
  assert.deepEqual(a.hits.filter((h) => h.piece === "hihat").map((h) => h.note), [42, 46, 44]);
  assert.equal(a.hiHat.initialForceN, HIHAT_CLOSED_FORCE_N, "first closed stroke is prepared under load");
  assert.deepEqual(pedals(a.hiHat.events), [pedal(0, 80), pedal(0.88, 0), pedal(2, 80)]);
  assert.deepEqual(strokes(a.hiHat.events).map((e) => e.t), [0, 1]);
  assert.equal(a.hits.find((h) => h.note === 44)?.g, 0);
  assert.equal(a.hits.find((h) => h.note === 38)?.t, 2, "other instruments keep their written onset");
  assert.deepEqual(a.report.drumsMapped, {}, "foot closure no longer reports a substituted stick stroke");
  assert.equal(a.report.drumsLayered, 0);
  assert.deepEqual(source, original, "the source MIDI is not rewritten");
});

test("repeated foot notes reopen the same pair before closing, without fake strikes", () => {
  const a = arrange(midi([[44, 0], [44, 960], [44, 1008]]), "foot.mid");
  assert.equal(a.hiHat.initialForceN, 0, "the first foot closure starts open");
  assert.deepEqual(pedals(a.hiHat.events), [pedal(0, 80), pedal(0.88, 0), pedal(1, 80), pedal(1.05, 80)]);
  assert.equal(strokes(a.hiHat.events).length, 0);
  assert.equal(a.hits.length, 3, "rapid foot commands survive the stick-repeat filter");
  assert.equal(a.report.drumsLayered, 0);
  assert.deepEqual(a.pieces, ["hihat"], "a pedal-only file still starts its physical lane");
  assert.equal(a.lastEventS, 1.05);
  assert.ok(a.durationS > a.lastEventS, "the foot-driven collision has time to ring");
});

test("rapid alternating articulations retain every command and clip preparation after preceding notes", () => {
  const a = arrange(midi([[46, 0], [42, 24], [44, 48], [46, 72]]), "rapid.mid");
  assert.deepEqual(a.hits.map((h) => h.note), [46, 42, 44, 46]);
  assert.deepEqual(pedals(a.hiHat.events), [pedal(0, 0), pedal(0, 80), pedal(0.025, 0), pedal(0.05, 80), pedal(0.075, 0)]);
  assert.deepEqual(strokes(a.hiHat.events).map((e) => e.t), [0, 0.025, 0.075]);
  assert.equal(a.report.drumsLayered, 0, "a new articulation is not a doubled stroke");
  assert.ok(a.hiHat.events.every((e, i, es) => !i || es[i - 1]!.t <= e.t));
  assert.ok(a.hiHat.events.indexOf(strokes(a.hiHat.events)[0]!) < a.hiHat.events.indexOf(pedals(a.hiHat.events)[1]!));
});

test("same-tick closed, foot and open requests preserve deterministic control-before-strike order", () => {
  const a = arrange(midi([[42, 0], [44, 0], [46, 0]]), "same-tick.mid");
  assert.deepEqual(a.hits.map((h) => h.note), [42, 44, 46]);
  assert.deepEqual(a.hiHat.events.map((e) => "pedalForceN" in e ? `pedal:${e.pedalForceN}` : "stick"),
    ["pedal:80", "stick", "pedal:80", "pedal:0", "stick"]);
  assert.ok(a.hiHat.events.every((e) => e.t === 0));
  assert.deepEqual(arrange(midi([[42, 0], [44, 0], [46, 0]]), "same-tick.mid").hiHat, a.hiHat);
  const doubled = arrange(midi([[42, 0], [42, 0], [44, 0], [44, 0]]), "doubling.mid");
  assert.equal(doubled.report.drumsLayered, 1, "duplicate stick notes still merge");
  assert.equal(strokes(doubled.hiHat.events).length, 1);
  assert.equal(doubled.hits.filter((h) => h.note === 44).length, 2, "foot commands do not merge with sticks or each other");
});

test("seeking restores preparatory pedal targets without replaying old strokes", () => {
  const a = arrange(midi([[42, 0], [46, 960], [44, 1920]]), "seek.mid");
  const original = structuredClone(a.hiHat);
  const duringPreparation = hiHatEventsFrom(a.hiHat, 0.9);
  assert.equal(duringPreparation.initialForceN, 0, "the open request preceded its future strike");
  assert.deepEqual(strokes(duringPreparation.events).map((e) => e.t), [1]);
  assert.deepEqual(pedals(duringPreparation.events), [pedal(2, 80)]);
  const onPreparation = hiHatEventsFrom(a.hiHat, 0.88);
  assert.equal(onPreparation.initialForceN, 80);
  assert.deepEqual(onPreparation.events[0], pedal(0.88, 0), "boundary control remains in the score");
  assert.deepEqual(strokes(hiHatEventsFrom(a.hiHat, 1).events).map((e) => e.t), [1], "boundary strike occurs once");
  const onFoot = hiHatEventsFrom(a.hiHat, 2);
  assert.equal(onFoot.initialForceN, 0);
  assert.deepEqual(onFoot.events, [pedal(2, 80)]);
  assert.deepEqual(hiHatEventsFrom(a.hiHat, 3), { initialForceN: 80, events: [] });
  assert.deepEqual(hiHatEventsFrom(a.hiHat, 0), a.hiHat);
  assert.deepEqual(hiHatEventsFrom({ initialForceN: 0, events: [] }, 0), { initialForceN: 0, events: [] });
  assert.deepEqual(a.hiHat, original);
});

test("tempo changes keep the physical pedal ramp and preparation at 120 ms", () => {
  const a = arrange(midi([[42, 0], [46, 960], [44, 1920], [44, 2880]]), "tempo.mid");
  for (const factor of [0.5, 1, 2]) {
    const performance = hiHatPerformance(a.hits, factor);
    const open = pedals(performance.events).find((e) => e.pedalForceN === 0)!;
    assert.ok(Math.abs((1 - open.t) / factor - HIHAT_TRANSITION_SECONDS) < 1e-12);
    assert.ok(pedals(performance.events).every((e) => e.transitionSeconds === HIHAT_TRANSITION_SECONDS));
    assert.deepEqual(strokes(performance.events).map((e) => e.t), [0, 1], "tempo does not rewrite note times");
    assert.equal(hiHatEventsFrom(performance, (open.t + 1) / 2).initialForceN, 0);
  }
  for (const factor of [0, -1, NaN, Infinity]) assert.throws(() => hiHatPerformance(a.hits, factor));
});

test("source volume and unsupported-note filtering happen before hi-hat gestures are authored", () => {
  const a = arrange(midi([[44, 0], [42, 480], [99, 960]], [
    { tick: 0, channel: 9, controller: 7, value: 0 },
    { tick: 480, channel: 9, controller: 7, value: 127 },
  ]), "filtered.mid");
  assert.deepEqual(a.hits.map((h) => h.note), [42]);
  assert.equal(a.hiHat.initialForceN, 80);
  assert.deepEqual(pedals(a.hiHat.events), [pedal(0.5, 80)]);
  assert.equal(strokes(a.hiHat.events).length, 1);
  assert.deepEqual(a.report.drumsSkipped, { "Note 99": 1 });
  assert.equal(a.report.drumNotesIn, 3);
  assert.ok(a.parts[0]!.actions.some((action) => action.includes("1 strokes silenced")));
  const kick = arrange(midi([[35, 0], [36, 0]]), "kick.mid");
  assert.deepEqual(kick.pieces, ["kick"]);
  assert.equal(kick.hits.length, 1);
  assert.equal(kick.report.drumsLayered, 1);
  assert.deepEqual(kick.hiHat, { initialForceN: 0, events: [] });
});
