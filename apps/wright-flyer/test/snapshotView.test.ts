// E5.2b snapshot-view battery (bead wf-root-guzez.6.3.2): decode
// per-slot oracles, interpolation laws (midpoint math, phase-boundary
// hold — a landing is never announced early, terminal hold), prop
// integration domain caps, pose/HUD/world mappings with exact-unit
// oracles, phase banners for every terminal.
// Repro: node --test test/snapshotView.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  advanceProp,
  controlStateFrom,
  decodeSnapshot,
  hudInputsFrom,
  interpolateSnapshots,
  phaseBanner,
  worldTransformFrom,
  type SimSnapshot,
} from "../src/sim/snapshotView.ts";
import { PAYLOAD_F64S, PAYLOAD_F64S_V1, PHASE_CODES } from "../src/sim/protocol.ts";

function payloadOf(fields: Partial<Record<number, number>>): Float64Array {
  const p = new Float64Array(PAYLOAD_F64S);
  for (const [k, v] of Object.entries(fields)) {
    p[Number(k)] = v as number;
  }
  return p;
}

function snap(overrides: Partial<SimSnapshot>): SimSnapshot {
  return {
    tick: 100,
    phase: "airborne",
    ended: false,
    xM: 10,
    hM: 3,
    uMps: 13,
    wMps: 1,
    qRadS: 0.2,
    thetaRad: 0.1,
    phiRad: 0.04,
    psiRad: -0.03,
    dcRad: 0.15,
    warpRad: 0.02,
    omegaPropRadS: 50,
    gustWMps: 0.05,
    assistActive: false,
    ...overrides,
  };
}

test("decode: per-slot oracle for every field and phase code", () => {
  const p = payloadOf({
    0: 12.5, // x
    1: 4.25, // h
    2: 13.5, // u
    3: -0.5, // w
    4: 0.3, // q
    5: 0.12, // theta
    6: 0.14, // dc
    7: -0.02, // warp
    8: 51.5, // omega
    9: 0.07, // gust
    10: 1, // assist
    11: 1, // airborne
    12: 0.09, // roll
    13: -0.04, // heading
  });
  const s = decodeSnapshot(360, p);
  assert.equal(s.tick, 360);
  assert.equal(s.phase, "airborne");
  assert.equal(s.ended, false);
  assert.equal(s.xM, 12.5);
  assert.equal(s.hM, 4.25);
  assert.equal(s.uMps, 13.5);
  assert.equal(s.wMps, -0.5);
  assert.equal(s.qRadS, 0.3);
  assert.equal(s.thetaRad, 0.12);
  assert.equal(s.dcRad, 0.14);
  assert.equal(s.warpRad, -0.02);
  assert.equal(s.omegaPropRadS, 51.5);
  assert.equal(s.gustWMps, 0.07);
  assert.equal(s.assistActive, true);
  assert.equal(s.phiRad, 0.09);
  assert.equal(s.psiRad, -0.04);
  // Every phase code round-trips; ended:* map to ended=true.
  for (const [word, code] of Object.entries(PHASE_CODES)) {
    const d = decodeSnapshot(1, payloadOf({ 11: code }));
    assert.equal(d.phase, word);
    assert.equal(d.ended, word.startsWith("ended:"));
  }
  // Unknown code and short payload throw (fail-closed).
  assert.throws(() => decodeSnapshot(1, payloadOf({ 11: 6 })), RangeError);
  assert.throws(() => decodeSnapshot(1, new Float64Array(PAYLOAD_F64S - 1)), RangeError);
  const legacy = decodeSnapshot(1, new Float64Array(PAYLOAD_F64S_V1));
  assert.equal(legacy.phiRad, 0, "v1 fallback is explicit zero-lateral");
  assert.equal(legacy.psiRad, 0);
});

test("interpolation: exact midpoint on continuous fields, discrete held", () => {
  const a = snap({
    tick: 100,
    xM: 10,
    hM: 2,
    thetaRad: 0.1,
    phiRad: 0.2,
    psiRad: -0.1,
    assistActive: true,
  });
  const b = snap({
    tick: 101,
    xM: 12,
    hM: 4,
    thetaRad: 0.3,
    phiRad: 0.4,
    psiRad: 0.1,
    assistActive: false,
  });
  const mid = interpolateSnapshots(a, b, 0.5);
  assert.equal(mid.xM, 11);
  assert.equal(mid.hM, 3);
  assert.ok(Math.abs(mid.thetaRad - 0.2) < 1e-15);
  assert.ok(Math.abs(mid.phiRad - 0.3) < 1e-15);
  assert.ok(Math.abs(mid.psiRad) < 1e-15);
  assert.equal(mid.assistActive, true, "discrete holds the older value");
  assert.equal(mid.tick, 100, "tick holds until alpha=1");
  assert.equal(interpolateSnapshots(a, b, 1).xM, 12);
  assert.equal(interpolateSnapshots(a, b, 0).xM, 10);
  // Clamped outside [0,1].
  assert.equal(interpolateSnapshots(a, b, 1.7).xM, 12);
  assert.equal(interpolateSnapshots(a, b, -0.3).xM, 10);
});

test("interpolation NEVER crosses a phase boundary early", () => {
  const flying = snap({ tick: 200, phase: "airborne", hM: 0.4 });
  const landed = snap({
    tick: 201,
    phase: "ended:ground-contact",
    ended: true,
    hM: 0,
  });
  // At any alpha < 1 the render plane still shows the airborne state.
  for (const alpha of [0, 0.25, 0.5, 0.99]) {
    const s = interpolateSnapshots(flying, landed, alpha);
    assert.equal(s.phase, "airborne", `alpha=${alpha}`);
    assert.equal(s.hM, 0.4, "no partial ground blend");
  }
  assert.equal(interpolateSnapshots(flying, landed, 1).phase, "ended:ground-contact");
});

test("prop integration: accumulates omega*dt; render-dt domain capped", () => {
  const s = snap({ omegaPropRadS: 50 });
  let d = { propAngleRad: 0 };
  d = advanceProp(d, s, 1 / 60);
  assert.ok(Math.abs(d.propAngleRad - 50 / 60) < 1e-15);
  d = advanceProp(d, s, 1 / 60);
  assert.ok(Math.abs(d.propAngleRad - 100 / 60) < 1e-15);
  // Domain: dt of exactly 1 s admits; beyond refuses (cap AND cap+ulp).
  assert.doesNotThrow(() => advanceProp(d, s, 1.0));
  assert.throws(() => advanceProp(d, s, 1.0000000000000002), RangeError);
  assert.throws(() => advanceProp(d, s, -0.001), RangeError);
  assert.throws(() => advanceProp(d, s, Number.NaN), RangeError);
});

test("control-state mapping: exact rad->deg, slaved rudder, prop angle", () => {
  const s = snap({ dcRad: Math.PI / 12, warpRad: -Math.PI / 36 });
  const c = controlStateFrom(s, { propAngleRad: 2.5 });
  assert.ok(Math.abs(c.canardDeg - 15) < 1e-12, "π/12 = 15°");
  assert.ok(Math.abs(c.warpDeg - -5) < 1e-12, "−π/36 = −5°");
  assert.equal(c.coupled, true, "1903 slaved-rudder wiring");
  assert.equal(c.rudderDeg, 0);
  assert.equal(c.propAngleRad, 2.5);
});

test("world transform: launch offset + simulated pitch, roll, and heading", () => {
  const s = snap({ xM: 25, hM: 6, thetaRad: 0.2, phiRad: -0.3, psiRad: 0.4 });
  const w = worldTransformFrom(s, [100, 5, -3]);
  assert.deepEqual(w.position, [125, 11, -3]);
  assert.equal(w.pitchRad, 0.2);
  assert.equal(w.rollRad, -0.3);
  assert.equal(w.headingRad, 0.4);
});

test("HUD inputs: airspeed hypot, elapsed from tick, engine rpm via 23:8", () => {
  const s = snap({ tick: 600, uMps: 3, wMps: 4, omegaPropRadS: 2 * Math.PI });
  const h = hudInputsFrom(s);
  assert.equal(h.airspeedMps, 5, "hypot(3,4)");
  assert.equal(h.elapsedS, 5, "600 ticks at 120 Hz");
  // 2π rad/s = 60 prop rpm -> engine 60 * 23/8 = 172.5 rpm.
  assert.ok(Math.abs(h.engineRpm - 172.5) < 1e-9);
  assert.equal(h.phase, "airborne");
});

test("phase banner: silent in flight, loud on every terminal", () => {
  assert.equal(phaseBanner(snap({ phase: "on-rail" })), null);
  assert.equal(phaseBanner(snap({ phase: "airborne" })), null);
  assert.match(String(phaseBanner(snap({ phase: "ended:ground-contact" }))), /LANDED/);
  assert.match(
    String(phaseBanner(snap({ phase: "ended:rail-end-without-lift" }))),
    /RAN OFF THE RAIL/,
  );
  assert.match(String(phaseBanner(snap({ phase: "ended:max-ticks" }))), /TIME LIMIT/);
  const banner = phaseBanner(
    snap({ phase: "ended:envelope-exceeded" }),
    "PropAirframeCouplingDidNotConverge",
  );
  assert.match(String(banner), /ENVELOPE/);
  assert.match(String(banner), /PropAirframeCouplingDidNotConverge/, "receipt code surfaced");
  console.info(JSON.stringify({ suite: "wf-e52b-view", case: "banners", ok: true }));
});
