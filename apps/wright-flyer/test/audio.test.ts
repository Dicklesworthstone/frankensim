import assert from "node:assert/strict";
import { test } from "node:test";

import {
  AUDIO_DISCLAIMER_LABEL,
  airframeNoiseGain,
  engineFreqHz,
  mixLevels,
  propBpfHz,
  rpm01FromOmega,
  twinPropBpfHz,
  windRushGain,
} from "../src/audio.ts";

test("engine frequency follows the 23:8 chain at two fires per rev", () => {
  // Trim: engine 1025 rpm -> prop 356.5 rpm -> 37.33 rad/s.
  const trimOmega = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  const f = engineFreqHz(trimOmega);
  const expected = (1025 / 60) * 2; // firing Hz at the trim
  assert.ok(Math.abs(f - expected) < 0.01, `trim firing ${f} ~ ${expected}`);
  assert.equal(engineFreqHz(0), 8, "floor at 8 Hz");
});

test("propBpfHz satisfies exact BPF math = blade_count * RPM / 60", () => {
  // Trim prop omega: 356.52 RPM = 5.942 RPS. 2 blades -> BPF = 11.884 Hz.
  const propOmegaTrim = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  const bpf2 = propBpfHz(propOmegaTrim, 2);
  const expectedBpf2 = (1025 * (8 / 23) * 2) / 60;
  assert.ok(
    Math.abs(bpf2 - expectedBpf2) < 1e-6,
    `2-blade BPF ${bpf2} matches theoretical ${expectedBpf2}`,
  );

  // 3-blade prop variant
  const bpf3 = propBpfHz(propOmegaTrim, 3);
  assert.ok(Math.abs(bpf3 - (3 / 2) * bpf2) < 1e-6, "3-blade BPF scales linearly with blade count");

  // Zero and negative checks (non-negative refusal)
  assert.equal(propBpfHz(0, 2), 0, "zero omega produces zero BPF");
  assert.equal(propBpfHz(-10, 2), 0, "negative omega refused to zero");
  assert.equal(propBpfHz(propOmegaTrim, 0), 0, "zero blade count produces zero BPF");
  assert.equal(propBpfHz(propOmegaTrim, -1), 0, "negative blade count produces zero BPF");
});

test("twinPropBpfHz generates partially-coherent left/right BPF with bounded beat frequency", () => {
  const propOmegaTrim = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;
  const nominal = propBpfHz(propOmegaTrim, 2);
  const twin = twinPropBpfHz(propOmegaTrim, 2, 0.04);

  assert.ok(twin.bpfLeftHz > 0 && twin.bpfRightHz > 0, "both rotor BPFs positive");
  assert.ok(
    Math.abs((twin.bpfLeftHz + twin.bpfRightHz) / 2 - nominal) < 1e-6,
    "average of twin BPFs equals nominal BPF",
  );
  assert.ok(
    Math.abs(twin.beatHz - (twin.bpfRightHz - twin.bpfLeftHz)) < 1e-6,
    "beat frequency equals difference between right and left",
  );
  assert.ok(twin.beatHz > 0.2 && twin.beatHz < 1.0, `beat frequency ${twin.beatHz} Hz in expected acoustic envelope`);

  // At zero speed, twin returns zero frequencies
  const zeroTwin = twinPropBpfHz(0, 2);
  assert.equal(zeroTwin.bpfLeftHz, 0);
  assert.equal(zeroTwin.bpfRightHz, 0);
  assert.equal(zeroTwin.beatHz, 0);
});

test("airframeNoiseGain models airspeed scaling and separation / stall buffet boost", () => {
  // At zero airspeed
  assert.equal(airframeNoiseGain(0, 0), 0);
  assert.equal(airframeNoiseGain(-5, 0), 0, "negative airspeed refused to zero");

  // At clean cruise airspeed (12 m/s, alpha = 4 deg)
  const cruiseGain = airframeNoiseGain(12, 4);
  assert.ok(cruiseGain > 0 && cruiseGain <= 0.25, `cruise noise ${cruiseGain} within base range`);

  // At high alpha (incipient separation / stall buffet, e.g. alpha = 14 deg)
  const highAlphaGain = airframeNoiseGain(12, 14);
  assert.ok(
    highAlphaGain > cruiseGain,
    `high alpha noise ${highAlphaGain} exceeds cruise noise ${cruiseGain}`,
  );
  assert.ok(highAlphaGain <= 0.35, "high alpha noise clamped at maximum envelope");
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

test("scripted flight fixture tracks parameters across flight stages", () => {
  const trimOmega = ((1025 * (8 / 23)) / 60) * 2 * Math.PI;

  interface FrameState {
    timeS: number;
    propOmega: number;
    airspeed: number;
    groundSpeed: number;
    onRail: boolean;
    alphaDeg: number;
  }

  // 5-stage scripted flight profile
  const flightStages: FrameState[] = [
    { timeS: 0.0, propOmega: trimOmega * 0.95, airspeed: 8.0, groundSpeed: 0.0, onRail: true, alphaDeg: 0 },
    { timeS: 2.5, propOmega: trimOmega * 1.0, airspeed: 10.5, groundSpeed: 3.5, onRail: true, alphaDeg: 2 },
    { timeS: 5.0, propOmega: trimOmega * 1.0, airspeed: 14.0, groundSpeed: 6.0, onRail: false, alphaDeg: 5 },
    { timeS: 8.0, propOmega: trimOmega * 0.98, airspeed: 15.0, groundSpeed: 7.0, onRail: false, alphaDeg: 12 },
    { timeS: 12.0, propOmega: trimOmega * 0.7, airspeed: 9.0, groundSpeed: 1.0, onRail: false, alphaDeg: 1 },
  ];

  for (const frame of flightStages) {
    const engF = engineFreqHz(frame.propOmega);
    const bpf = propBpfHz(frame.propOmega, 2);
    const twin = twinPropBpfHz(frame.propOmega, 2);
    const mix = mixLevels(rpm01FromOmega(frame.propOmega), frame.airspeed, frame.onRail, frame.groundSpeed);
    const wireNoise = airframeNoiseGain(frame.airspeed, frame.alphaDeg);

    assert.ok(engF >= 8 && engF <= 40, `engine freq ${engF} valid at t=${frame.timeS}`);
    assert.ok(bpf >= 0 && bpf <= 15, `bpf ${bpf} valid at t=${frame.timeS}`);
    assert.ok(twin.bpfLeftHz <= twin.bpfRightHz, `twin order valid at t=${frame.timeS}`);
    assert.ok(mix.engine >= 0 && mix.engine <= 0.5, `mix.engine valid at t=${frame.timeS}`);
    assert.ok(wireNoise >= 0 && wireNoise <= 0.35, `wireNoise valid at t=${frame.timeS}`);
    if (frame.onRail) {
      assert.ok(mix.rumble >= 0, `rumble active on rail at t=${frame.timeS}`);
    } else {
      assert.equal(mix.rumble, 0, `rumble zero off rail at t=${frame.timeS}`);
    }
  }
});

test("sound design label disclaimer is explicitly defined", () => {
  assert.ok(AUDIO_DISCLAIMER_LABEL.length > 20);
  assert.ok(AUDIO_DISCLAIMER_LABEL.includes("sound design"));
  assert.ok(AUDIO_DISCLAIMER_LABEL.includes("not an acoustic simulation claim"));
});
