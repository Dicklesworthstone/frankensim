import { test } from "bun:test";
import assert from "node:assert/strict";
import { arrange } from "../src/midi/arrange";
import { decodeMidi, tempoClock } from "../src/midi/smf";

function vlq(n: number): number[] {
  const bytes = [n & 127];
  while ((n = Math.floor(n / 128)) > 0) bytes.unshift(128 | (n & 127));
  return bytes;
}

function score(division: number): Uint8Array {
  const track = (...events: number[][]) => {
    const body = [...events.flat(), 0, 255, 47, 0];
    return [77, 84, 114, 107, 0, 0, 0, body.length, ...body];
  };
  return Uint8Array.from([
    77, 84, 104, 100, 0, 0, 0, 6, 0, 1, 0, 2, division >>> 8, division & 255,
    ...track(
      [...vlq(1000), 255, 81, 3, 15, 66, 64], // tempo becomes 60 bpm
      [...vlq(2500), 255, 88, 4, 3, 2, 24, 8], // open a new 3/4 bar
    ),
    ...track(
      [...vlq(500), 144, 60, 96],
      [...vlq(250), 176, 64, 127],
      [...vlq(5500), 176, 64, 0],
      [...vlq(750), 128, 60, 0],
    ),
  ]);
}

test("frame clocks keep performed notes and pedals in time while tempo and meter set bar lines", () => {
  const result = decodeMidi(score(0xe728)); // -25 fps, 40 ticks per frame
  assert.ok(result.ok);
  const song = arrange(result.midi, "timecode.mid");
  assert.equal(song.tempoBpm, 120); // first authored tempo arrives one second in
  assert.equal(song.notes[0]?.t, 0.5);
  assert.equal(song.notes[0]?.off, 7);
  assert.deepEqual(song.events.filter(e => e.sustain !== undefined).map(e => e.t), [0.75, 6.25]);
  assert.deepEqual(song.bars.map(bar => bar.t), [0, 3, 3.5, 6.5]);
});

test("all four SMPTE frame codes preserve absolute time across tempo changes", () => {
  for (const [code, fps] of [[24, 24], [25, 25], [29, 30000 / 1001], [30, 30]] as const) {
    const result = decodeMidi(score(((256 - code) << 8) | 80));
    assert.ok(result.ok);
    const seconds = tempoClock(result.midi);
    assert.ok(Math.abs(seconds(2_400_000) - 30000 / fps) < 1e-9);
  }
  for (const division of [0xe700, 0xe628]) assert.equal(decodeMidi(score(division)).ok, false);
});

test("a PPQN performance uses the same tempo and meter events on its musical clock", () => {
  const result = decodeMidi(score(1000));
  assert.ok(result.ok);
  const song = arrange(result.midi, "musical.mid");
  assert.equal(song.notes[0]?.t, 0.25);
  assert.equal(song.notes[0]?.off, 6.5);
  assert.deepEqual(song.bars.map(bar => bar.t), [0, 3, 6]);
});
