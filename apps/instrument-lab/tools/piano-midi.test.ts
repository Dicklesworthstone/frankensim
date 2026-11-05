import { test } from 'bun:test';
import assert from 'node:assert/strict';
import { PianoMidi, combinePianoPedals, type PianoMidiPedals } from '../src/live/pianoMidi';
import type { MidiChannelMessage } from '../src/live/midiInput';

const on = (channel: number, key: number, velocity = 90): MidiChannelMessage => ({ kind: 'noteon', channel, key, velocity });
const off = (channel: number, key: number): MidiChannelMessage => ({ kind: 'noteoff', channel, key, velocity: 0 });
const cc = (channel: number, controller: number, value = 0): MidiChannelMessage => ({ kind: 'cc', channel, controller, value });
const up = (): PianoMidiPedals => ({ sustain: 0, sostenuto: false, unaCorda: false });

function rig() {
  const notes: Array<{ key: number; velocity?: number; acquireHold?: boolean }> = [];
  const pedals: PianoMidiPedals[] = [];
  const midi = new PianoMidi({
    noteOn: (key, velocity, acquireHold) => notes.push({ key, velocity, acquireHold }),
    noteOff: key => notes.push({ key }),
    pedals: value => pedals.push(value),
  });
  return { midi, notes, pedals };
}

test('piano MIDI defaults to pitched channels and allows an explicit channel 10 selection', () => {
  const { midi, notes, pedals } = rig();
  midi.receive(on(9, 42));
  midi.receive(cc(9, 64, 127));
  midi.receive(on(0, 20));
  midi.receive(on(0, 109));
  assert.deepEqual(notes, []);
  assert.deepEqual(pedals, []);
  midi.receive(on(0, 21));
  midi.receive(on(15, 108));
  assert.deepEqual(notes, [
    { key: 21, velocity: 90, acquireHold: true }, { key: 108, velocity: 90, acquireHold: true },
  ]);
  midi.selectChannel(9);
  assert.deepEqual(notes.slice(2), [{ key: 21 }, { key: 108 }]);
  midi.receive(on(0, 60));
  midi.receive(on(9, 42));
  assert.deepEqual(notes.at(-1), { key: 42, velocity: 90, acquireHold: true });
  assert.throws(() => midi.selectChannel(16), RangeError);
  assert.throws(() => midi.selectChannel(NaN), RangeError);
  midi.receive(off(9, 42));
  assert.deepEqual(notes.at(-1), { key: 42 });
});

test('repeated note-ons retrigger but one note-off releases each channel owner', () => {
  const { midi, notes } = rig();
  midi.receive(on(0, 60, 40));
  midi.receive(on(0, 60, 96));
  midi.receive(on(1, 60, 70));
  midi.receive(off(2, 60));
  midi.receive(off(0, 64));
  midi.receive(off(0, 60));
  midi.receive(off(0, 60));
  assert.deepEqual(notes, [
    { key: 60, velocity: 40, acquireHold: true },
    { key: 60, velocity: 96, acquireHold: false },
    { key: 60, velocity: 70, acquireHold: false },
  ], 'every velocity is preserved and another channel keeps the key down');
  midi.receive(on(1, 60, 0));
  midi.receive(off(1, 60));
  assert.deepEqual(notes.slice(3), [{ key: 60 }], 'velocity-zero and duplicate releases cannot strand or steal a hold');
});

test('MIDI reset and unmatched releases preserve keys owned by pointer, keyboard or demo', () => {
  const holds = new Map([[60, 1], [64, 1], [67, 1]]);
  const midi = new PianoMidi({
    noteOn: (key, _velocity, acquireHold) => {
      if (acquireHold) holds.set(key, (holds.get(key) ?? 0) + 1);
    },
    noteOff: key => holds.set(key, (holds.get(key) ?? 0) - 1),
    pedals: () => {},
  });
  midi.receive(on(0, 60));
  midi.receive(on(0, 60));
  midi.receive(on(1, 60));
  midi.receive(off(0, 64));
  midi.receive(off(1, 67));
  assert.deepEqual([...holds.values()], [2, 1, 1]);
  midi.reset();
  midi.receive(off(0, 60));
  midi.receive(off(1, 60));
  assert.deepEqual([...holds.values()], [1, 1, 1], 'late releases after silence do not release a new local hold');
  midi.receive(on(0, 64));
  midi.selectChannel(2);
  midi.receive(off(0, 64));
  assert.deepEqual([...holds.values()], [1, 1, 1], 'changing channel releases only MIDI ownership');
});

test('half pedal and the two switched pedals aggregate across selected channels', () => {
  const { midi, pedals } = rig();
  midi.receive(cc(0, 64, 90));
  midi.receive(cc(1, 64, 60));
  midi.receive(cc(0, 66, 63));
  assert.deepEqual(pedals, [{ ...up(), sustain: 90 / 127 }]);
  midi.receive(cc(0, 66, 64));
  midi.receive(cc(1, 66, 127));
  midi.receive(cc(0, 66, 0));
  assert.equal(pedals.length, 2, 'another channel keeps sostenuto down without recapturing keys');
  midi.receive(cc(0, 67, 127));
  midi.receive(cc(1, 67, 64));
  midi.receive(cc(0, 121));
  assert.deepEqual(pedals.at(-1), { sustain: 60 / 127, sostenuto: true, unaCorda: true });
  midi.receive(cc(1, 121));
  assert.deepEqual(pedals.at(-1), up());
});

test('channel mode messages release only their channel and reset controllers without releasing keys', () => {
  const { midi, notes, pedals } = rig();
  midi.receive(on(0, 60));
  midi.receive(on(1, 60));
  midi.receive(on(1, 64));
  midi.receive(cc(0, 64, 127));
  midi.receive(cc(1, 64, 64));
  midi.receive(cc(0, 123));
  assert.equal(notes.length, 3, 'all notes off on channel 1 preserves channel 2 holds');
  assert.equal(pedals.at(-1)?.sustain, 1, 'CC123 respects the held sustain pedal');
  midi.receive(cc(0, 121));
  assert.equal(notes.length, 3, 'reset controllers does not release any held key');
  assert.equal(pedals.at(-1)?.sustain, 64 / 127);
  midi.receive(cc(1, 120));
  assert.deepEqual(notes.slice(3), [{ key: 60 }, { key: 64 }]);
  assert.deepEqual(pedals.at(-1), up(), 'CC120 also releases that channel\'s pedal claims');
  midi.receive(cc(1, 120));
  assert.equal(notes.length, 5, 'repeated all-off messages cannot release local controls');
});

test('source reset and channel selection clear MIDI pedals without disturbing local pedal claims', () => {
  const local = { sustain: 0.35, sostenuto: true, unaCorda: false };
  const { midi, notes, pedals } = rig();
  midi.receive(on(0, 60));
  midi.receive(cc(0, 64, 100));
  midi.receive(cc(0, 67, 127));
  assert.deepEqual(combinePianoPedals(local, pedals.at(-1)!), {
    sustain: 100 / 127, sostenuto: true, unaCorda: true,
  });
  midi.reset();
  assert.deepEqual(combinePianoPedals(local, pedals.at(-1)!), local);
  assert.deepEqual(notes.at(-1), { key: 60 });
  midi.selectChannel(2);
  midi.receive(on(2, 64));
  midi.receive(cc(2, 66, 127));
  const before = notes.length;
  midi.selectChannel(2);
  assert.equal(notes.length, before, 'selecting the same channel preserves performance');
  midi.selectChannel('pitched');
  assert.deepEqual(notes.at(-1), { key: 64 });
  assert.deepEqual(pedals.at(-1), up());
  assert.deepEqual(local, { sustain: 0.35, sostenuto: true, unaCorda: false });
});
