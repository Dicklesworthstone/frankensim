import { test } from "bun:test";
import assert from "node:assert/strict";
import { orchestrate, type SourcePart } from "../src/midi/orchestrate";

function melody(levels: number[], velocities = levels.map(() => 80)): SourcePart {
  return {
    track: 0, channel: 0, name: "Melody", program: 0, baseRole: "melody",
    notes: levels.map((level, i) => ({
      t: i, off: i + 0.2, key: 60 + i, vel: velocities[i]!, level,
    })),
  };
}

test("uniform written velocities retain the file's channel volume and expression", () => {
  const full = melody([1, 1, 1]);
  const half = melody([0.5, 0.5, 0.5]);
  const quarter = melody([0.25, 0.25, 0.25]);
  const originals = structuredClone([full, half, quarter]);
  assert.deepEqual(orchestrate([full]).notes.map(n => n.vel), [84, 84, 84]);
  assert.deepEqual(orchestrate([half]).notes.map(n => n.vel), [42, 42, 42],
    "a quieter channel must stay quieter after the role's written-velocity calibration");
  assert.deepEqual(orchestrate([quarter]).notes.map(n => n.vel), [21, 21, 21]);
  assert.deepEqual([full, half, quarter], originals, "the source performance stays intact");
});

test("an expression ramp does not get normalized into the velocity ceiling", () => {
  const result = orchestrate([melody([0, 0.25, 0.5, 1])]);
  assert.equal(result.silenced, 1);
  assert.deepEqual(result.notes.map(n => [n.t, n.vel]), [[1, 21], [2, 42], [3, 84]]);
  assert.equal(result.plans[0]!.notesOut, 3);
});

test("nonuniform written dynamics retain their existing controller scaling", () => {
  const result = orchestrate([melody([0.5, 0.5, 0.5], [40, 80, 100])]);
  assert.deepEqual(result.notes.map(n => n.vel), [20, 40, 50]);
});
