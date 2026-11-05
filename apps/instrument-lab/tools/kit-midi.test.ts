import { test } from "bun:test";
import assert from "node:assert/strict";
import { LiveKitMidi, type LiveKitPieceId } from "../src/live/kitMidi";
import { decodeMidiMessage } from "../src/live/midiInput";

type Action = { kind: "strike"; piece: LiveKitPieceId; u: number; v: number; g: number; frame: number }
  | { kind: "pedal"; forceN: number; transitionSeconds: number; frame: number };

function performance() {
  const actions: Action[] = [], unsupported: number[] = [];
  const midi = new LiveKitMidi({
    strike: (piece, u, v, g, frame) => actions.push({ kind: "strike", piece, u, v, g, frame }),
    pedal: (forceN, transitionSeconds, frame) => actions.push({ kind: "pedal", forceN, transitionSeconds, frame }),
    unsupported: note => unsupported.push(note),
  });
  const send = (data: number[], frame = 1000): void => {
    const message = decodeMidiMessage(data);
    if (message) midi.receive(message, frame);
  };
  return { midi, send, actions, unsupported };
}
const pedal = (forceN: number, frame: number): Action => ({ kind: "pedal", forceN, transitionSeconds: 0.12, frame });
const hat = (frame: number): Action => ({ kind: "strike", piece: "hihat", u: 0.85, v: 0, g: 1, frame });

test("GM hardware notes reach all seven pieces with the existing physical strike locations", () => {
  const cases: [number, LiveKitPieceId, number, number, number][] = [
    [36, "kick", 0, 0.02, 1], [38, "snare", 0.35, 0.1, 1], [37, "snare", 0.86, 0.05, 0.75],
    [41, "tom-14", 0.35, 0.1, 1], [45, "tom-12", 0.4, 0.1, 1],
    [46, "hihat", 0.85, 0, 1], [49, "crash", 0.85, 0, 1], [53, "ride", 0.12, 0, 1.1],
  ];
  for (const [key, piece, u, v, g] of cases) {
    const p = performance();
    p.send([0x99, key, 100], 480);
    assert.deepEqual(p.actions, [
      ...(piece === "hihat" ? [pedal(0, 480)] : []),
      { kind: "strike", piece, u, v, g, frame: 480 },
    ]);
    assert.deepEqual(p.unsupported, []);
  }
});

test("velocity drives stick speed, and rapid repeated hardware hits are neither thinned nor randomized", () => {
  const p = performance();
  for (const velocity of [1, 64, 100, 127, 100, 100]) p.send([0x99, 38, velocity], 500);
  const strikes = p.actions.filter(a => a.kind === "strike");
  assert.equal(strikes.length, 6);
  assert.equal(strikes[0]!.g, 0.25);
  assert.ok(strikes[1]!.g > 0.25 && strikes[1]!.g < 1);
  assert.equal(strikes[2]!.g, 1);
  assert.ok(strikes[3]!.g > 1 && strikes[3]!.g <= 1.35);
  assert.deepEqual(strikes[4], strikes[2]);
  assert.deepEqual(strikes[5], strikes[2]);
});

test("discrete open/closed hits order pedal before strike on one frame; foot notes never launch sticks", () => {
  const p = performance();
  p.send([0x99, 42, 100], 512);
  p.send([0x99, 46, 100], 513);
  p.send([0x99, 44, 100], 514);
  p.send([0x99, 44, 127], 515);
  p.send([0x99, 42, 100], 516);
  assert.deepEqual(p.actions, [pedal(80, 512), hat(512), pedal(0, 513), hat(513),
    pedal(80, 514), pedal(80, 515), pedal(80, 516), hat(516)]);
});

test("CC4 preserves a half-open hat across both stick articulations and foot notifications", () => {
  const p = performance();
  p.send([0xb9, 4, 64], 600);
  p.send([0x99, 42, 100], 601);
  p.send([0x99, 46, 100], 602);
  p.send([0x99, 44, 100], 603);
  p.send([0xb9, 4, 0], 604);
  p.send([0xb9, 4, 127], 605);
  assert.deepEqual(p.actions, [pedal(80 * 64 / 127, 600), hat(601), hat(602), pedal(0, 604), pedal(80, 605)]);
});

test("channel selection defaults to GM10, supports both boundaries/all, and refuses invalid selections", () => {
  const p = performance();
  for (const channel of [-1, 16, 0.5, NaN, Infinity]) assert.throws(() => p.midi.setChannel(channel), RangeError);
  p.send([0x90, 38, 100]);
  p.send([0x99, 38, 100], 9);
  p.midi.setChannel(0);
  p.send([0x99, 38, 100]);
  p.send([0x90, 38, 100], 0);
  p.midi.setChannel(15);
  p.send([0x9f, 38, 100], 15);
  p.midi.setChannel(null);
  p.send([0x90, 38, 100], 16);
  p.send([0x9f, 38, 100], 17);
  assert.deepEqual(p.actions.map(a => a.frame), [9, 0, 15, 16, 17]);
});

test("CC4 authority belongs to its channel and resets without moving or silencing the instrument", () => {
  const p = performance();
  p.midi.setChannel(null);
  p.send([0xb9, 4, 64], 700);
  p.send([0x90, 46, 100], 701);
  p.send([0x99, 42, 100], 702);
  assert.deepEqual(p.actions, [pedal(80 * 64 / 127, 700), pedal(0, 701), hat(701), hat(702)]);
  p.midi.reset(703);
  assert.equal(p.actions.length, 4, "disconnect reset supplies no new physical impulse or audio choke");
  p.send([0x99, 42, 100], 704);
  assert.deepEqual(p.actions.slice(4), [pedal(80, 704), hat(704)]);
  p.send([0xb9, 4, 32], 705);
  p.midi.setChannel(0);
  p.midi.setChannel(9);
  p.send([0x99, 46, 100], 706);
  assert.deepEqual(p.actions.slice(-2), [pedal(0, 706), hat(706)], "device/channel changes cannot retain old pedal authority");
});

for (const controller of [120, 121, 123]) test(`CC${controller} clears only its channel's input state, without a fake strike or choke`, () => {
  const p = performance();
  p.midi.setChannel(null);
  p.send([0xb9, 4, 64]);
  p.send([0xb0, 4, 32]);
  const before = [...p.actions];
  p.send([0xb9, controller, 0]);
  assert.deepEqual(p.actions, before);
  p.send([0x90, 42, 100], 800);
  p.send([0x99, 42, 100], 801);
  assert.deepEqual(p.actions.slice(2), [hat(800), pedal(80, 801), hat(801)]);
});

test("note-offs, velocity-zero, unsupported statuses and malformed MIDI bytes cannot create hits", () => {
  const p = performance();
  for (const data of [[0x99, 42, 0], [0x89, 42, 127], [0x99, 44, 0], [0x89, 44, 64],
    [0xa9, 49, 127], [0xe9, 0, 64], [0xb9, 4, 128], [0x99, 128, 100], [0x99, 42],
    [0x99, 42, -1], [0x99, 42, NaN], [0x99, 42, 1.5], [0xf8], [0xf0, 4, 127]]) p.send(data);
  assert.deepEqual(p.actions, []);
  assert.deepEqual(p.unsupported, []);
});

test("unavailable pieces and unknown notes are reported without substituting another physical instrument", () => {
  const p = performance();
  for (const key of [65, 66, 0, 127]) p.send([0x99, key, 100]);
  assert.deepEqual(p.actions, []);
  assert.deepEqual(p.unsupported, [65, 66, 0, 127]);
  p.send([0x90, 65, 100]);
  p.send([0x99, 65, 0]);
  assert.equal(p.unsupported.length, 4, "inactive channels and releases are not new unsupported strikes");
});

test("invalid frame coordinates refuse before changing input state or sending physical commands", () => {
  const p = performance();
  for (const frame of [-1, NaN, Infinity]) assert.throws(() => p.send([0xb9, 4, 64], frame), RangeError);
  assert.deepEqual(p.actions, []);
  p.send([0x99, 42, 100], 0);
  assert.deepEqual(p.actions, [pedal(80, 0), hat(0)]);
});
