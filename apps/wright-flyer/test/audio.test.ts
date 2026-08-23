// Audio mix battery (Track 4): the PURE mapping laws — firing
// frequency from prop omega, the clamped mix, and the mute-safe
// invariants. No AudioContext needed (structural typing in audio.ts).
// Repro: node --test test/audio.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import { engineFreqHz, mixLevels, rpm01FromOmega, windRushGain } from "../src/audio.ts";

test("engine frequency follows the 23:8 chain at two fires per rev", () => {
  // Trim: engine 1025 rpm -> prop 356.5 rpm -> 37.33 rad/s.
  const trimOmega = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  const f = engineFreqHz(trimOmega);
  const expected = 1025 / 60 * 2; // firing Hz at the trim
  assert.ok(Math.abs(f - expected) < 0.01, `trim firing ${f} ~ ${expected}`);
  assert.equal(engineFreqHz(0), 8, "floor at 8 Hz");
});

test("rpm01 normalizes against the trim and clamps", () => {
  const trimOmega = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  assert.ok(Math.abs(rpm01FromOmega(trimOmega) - 1) < 1e-9);
  assert.equal(rpm01FromOmega(0), 0);
  assert.equal(rpm01FromOmega(trimOmega * 5), 1, "clamped high");
  assert.ok(rpm01FromOmega(-3) === 0, "negative clamped to 0");
});

test("mix is clamped, zero at rest, and rumble rides the rail only", () => {
  const idle = mixLevels(0, 0, true, 0);
  assert.equal(idle.engine, 0);
  assert.equal(idle.wind, 0);
  assert.equal(idle.rumble, 0);
  const full = mixLevels(1.4, 40, false, 30); // overdriven inputs
  assert.ok(full.engine <= 0.5 && full.wind <= 0.4 && full.rumble === 0);
  const rolling = mixLevels(0.6, 8, true, 10);
  assert.ok(rolling.rumble > 0 && rolling.rumble <= 0.3);
  assert.ok(rolling.engine > 0 && rolling.wind > 0);
});

test("wind rush obeys the square law and clamps", () => {
  assert.equal(windRushGain(0), 0);
  assert.equal(windRushGain(4), 0, "below the walking-headwind floor");
  const a = windRushGain(17); // the December 17 headwind
  const b = windRushGain(34);
  assert.ok(a > 0 && a < 0.32, "in range at racing speed");
  assert.ok(Math.abs(b - 4 * a) < 1e-9 || b === 0.32, "quadratic until clamped");
  assert.equal(windRushGain(-5), 0, "negative refused to zero");
});
