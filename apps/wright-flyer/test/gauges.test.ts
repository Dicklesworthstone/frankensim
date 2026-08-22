// Gauge-math battery (UI overhaul, root guzez): needle pinning AT the
// stop and one-past-the-stop (the analog cap-and-cap+1), monotone
// value->angle mapping, closed-boundary redline classification, tick
// generation at both ends, per-instrument unit oracles (mph/ft/rpm
// conversions asserted against hand-computed values, not round-trips),
// lever pinning, and the phase display passing terminal banners
// through UNCHANGED.
// Repro: node --test test/gauges.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  ANEMOMETER_SPEC,
  CANARD_STOP_RAD,
  IDLE_INPUTS,
  INCLINOMETER_SPEC,
  REV_COUNTER_SPEC,
  clockText,
  dialSetFrom,
  inRedline,
  leverSetFrom,
  needleDeg,
  phaseDisplay,
  tickMarks,
} from "../src/gauges.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-gauges","case":"${kase}",${payload}}`);
}

test("needle pins at the stop AND one past the stop", () => {
  const spec = ANEMOMETER_SPEC;
  const atMax = needleDeg(spec.max, spec);
  const pastMax = needleDeg(spec.max + 1, spec);
  const farPast = needleDeg(1e9, spec);
  assert.equal(atMax, spec.startDeg + spec.sweepDeg);
  assert.equal(pastMax, atMax, "cap+1 sits ON the stop");
  assert.equal(farPast, atMax);
  const atMin = needleDeg(spec.min, spec);
  assert.equal(atMin, spec.startDeg);
  assert.equal(needleDeg(spec.min - 1, spec), atMin);
  // NaN reads as the rest position, never a poisoned transform.
  assert.equal(needleDeg(Number.NaN, spec), atMin);
  jlog("needle-pins", `"at_max_deg":${atMax}`);
});

test("value->angle is monotone and linear at the midpoint", () => {
  const spec = ANEMOMETER_SPEC;
  let prev = needleDeg(spec.min, spec);
  for (let v = spec.min + 1; v <= spec.max; v += 1) {
    const d = needleDeg(v, spec);
    assert.ok(d > prev, `monotone at ${v}`);
    prev = d;
  }
  const mid = needleDeg((spec.min + spec.max) / 2, spec);
  assert.ok(Math.abs(mid - (spec.startDeg + spec.sweepDeg / 2)) < 1e-12);
  jlog("monotone", `"mid_deg":${mid}`);
});

test("malformed specs refuse", () => {
  assert.throws(() =>
    needleDeg(1, { min: 5, max: 5, startDeg: 0, sweepDeg: 90, redline: null }),
  );
  assert.throws(() =>
    needleDeg(1, { min: 0, max: 1, startDeg: 0, sweepDeg: 0, redline: null }),
  );
  assert.throws(() => tickMarks(ANEMOMETER_SPEC, 1));
  jlog("refusals", `"count":3`);
});

test("redline is closed at BOTH boundaries", () => {
  const spec = REV_COUNTER_SPEC;
  const lo = spec.redline![0];
  assert.equal(inRedline(lo, spec), true, "AT the redline start is IN");
  assert.equal(inRedline(lo - 1e-9, spec), false, "just below is OUT");
  assert.equal(inRedline(spec.max, spec), true);
  assert.equal(inRedline(0, spec), false);
  jlog("redline-boundary", `"from":${lo}`);
});

test("tick marks span min to max inclusive", () => {
  const ticks = tickMarks(ANEMOMETER_SPEC, 10);
  assert.equal(ticks.length, 10);
  assert.equal(ticks[0]!.value, ANEMOMETER_SPEC.min);
  assert.equal(ticks[9]!.value, ANEMOMETER_SPEC.max);
  assert.equal(ticks[0]!.deg, ANEMOMETER_SPEC.startDeg);
  assert.equal(ticks[9]!.deg, ANEMOMETER_SPEC.startDeg + ANEMOMETER_SPEC.sweepDeg);
  jlog("ticks", `"count":${ticks.length}`);
});

test("per-instrument unit oracles (hand-computed, not round-trips)", () => {
  // 12 m/s = 26.84 mph; 30.48 m = 100.0 ft; theta 0.1 rad = 5.73 deg.
  const dials = dialSetFrom({
    airspeedMps: 12,
    engineRpm: 1025,
    elapsedS: 59,
    hM: 30.48,
    thetaRad: 0.1,
    dcRad: 0.1,
    warpRad: 0.05,
  });
  const byId = new Map(dials.map((d) => [d.id, d]));
  assert.equal(byId.get("anemometer")!.reading, "27");
  assert.equal(byId.get("revcounter")!.reading, "1025");
  assert.equal(byId.get("altimeter")!.reading, "100");
  assert.equal(byId.get("inclinometer")!.reading, "5.7");
  assert.equal(byId.get("stopwatch")!.reading, "0:59.0");
  // Provenance stamps: the period triad vs modern overlays.
  assert.equal(byId.get("anemometer")!.provenance, "PERIOD");
  assert.equal(byId.get("revcounter")!.provenance, "PERIOD");
  assert.equal(byId.get("stopwatch")!.provenance, "PERIOD");
  assert.equal(byId.get("altimeter")!.provenance, "MODERN");
  assert.equal(byId.get("inclinometer")!.provenance, "MODERN");
  assert.equal(dials.length, 5);
  jlog("unit-oracles", `"dials":${dials.length}`);
});

test("stopwatch clock text rolls the minute", () => {
  assert.equal(clockText(0), "0:00.0");
  assert.equal(clockText(59.94), "0:59.9");
  assert.equal(clockText(60), "1:00.0");
  assert.equal(clockText(-5), "0:00.0", "negative time pins at zero");
  // Regression (fresh-eyes review): rounding must precede decomposition
  // or 59.96 renders the impossible "0:60.0".
  assert.equal(clockText(59.96), "1:00.0");
  assert.equal(clockText(119.97), "2:00.0");
  jlog("clock", `"minute_roll":"1:00.0"`);
});

test("tick labels divide by labelDiv (tachometer x100 face)", () => {
  const dials = dialSetFrom(IDLE_INPUTS);
  const byId = new Map(dials.map((d) => [d.id, d]));
  assert.equal(byId.get("revcounter")!.labelDiv, 100);
  assert.ok(byId.get("revcounter")!.label.includes("×100"), "face declares its scale");
  for (const id of ["anemometer", "stopwatch", "altimeter", "inclinometer"]) {
    assert.equal(byId.get(id)!.labelDiv, 1, id);
  }
  jlog("label-div", `"revcounter":100`);
});

test("mirrored-redline precondition: pitch dial is symmetric", () => {
  // The panel mirrors the danger arc only when min = -max and the
  // redline sits in the positive half; the pitch spec must satisfy
  // that or nose-down danger loses its face marking silently.
  assert.equal(INCLINOMETER_SPEC.min, -INCLINOMETER_SPEC.max);
  assert.ok(INCLINOMETER_SPEC.redline![0] > 0);
  jlog("mirror-precondition", `"symmetric":true`);
});

test("levers pin at the stop and report AT-STOP", () => {
  const at = leverSetFrom({ ...IDLE_INPUTS, dcRad: CANARD_STOP_RAD });
  assert.equal(at[0]!.fraction, 1);
  assert.equal(at[0]!.atStop, true);
  const past = leverSetFrom({ ...IDLE_INPUTS, dcRad: CANARD_STOP_RAD * 1.5 });
  assert.equal(past[0]!.fraction, 1, "past the stop still reads the stop");
  const neutral = leverSetFrom(IDLE_INPUTS);
  assert.equal(neutral[0]!.fraction, 0);
  assert.equal(neutral[0]!.atStop, false);
  assert.equal(neutral[1]!.id, "warp");
  jlog("levers", `"at_stop":true`);
});

test("engine redline engages above 1250 rpm", () => {
  const hot = dialSetFrom({ ...IDLE_INPUTS, engineRpm: 1300 });
  const cool = dialSetFrom({ ...IDLE_INPUTS, engineRpm: 1025 });
  assert.equal(hot.find((d) => d.id === "revcounter")!.danger, true);
  assert.equal(cool.find((d) => d.id === "revcounter")!.danger, false);
  jlog("engine-redline", `"threshold":1250`);
});

test("pitch redline uses magnitude", () => {
  const specDeg = INCLINOMETER_SPEC.redline![0];
  const noseUp = dialSetFrom({ ...IDLE_INPUTS, thetaRad: (specDeg + 1) / (180 / Math.PI) });
  const noseDown = dialSetFrom({ ...IDLE_INPUTS, thetaRad: -(specDeg + 1) / (180 / Math.PI) });
  assert.equal(noseUp.find((d) => d.id === "inclinometer")!.danger, true);
  assert.equal(noseDown.find((d) => d.id === "inclinometer")!.danger, true);
  jlog("pitch-redline", `"deg":${specDeg}`);
});

test("phase display passes terminal banners through UNCHANGED", () => {
  const banner = "FLIGHT LEFT THE CERTIFIED ENVELOPE (EnvelopeExceeded)";
  const end = phaseDisplay("ended:envelope-exceeded", banner);
  assert.equal(end.text, banner, "receipted refusal text is not rewritten");
  assert.equal(end.tone, "end");
  const rail = phaseDisplay("on-rail", null);
  assert.equal(rail.tone, "info");
  const air = phaseDisplay("airborne", null);
  assert.equal(air.tone, "good");
  jlog("phase-display", `"terminal_text_preserved":true`);
});
