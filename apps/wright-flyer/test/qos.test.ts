// E5.6 QoS battery (bead wf-root-guzez.6.9): exact dwell-count
// escalation/de-escalation oracles, CHATTER RESISTANCE (oscillation
// around a threshold never flaps the state), one-level-at-a-time law,
// the persistent-Critical typed refusal at cap AND cap+1 (emitted
// exactly once), the tier-immutability hostile twin (a physics knob in
// a profile REFUSES), badge exactness, and frozen atomic profiles.
// Repro: node --test test/qos.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  PROFILES,
  QOS_V1,
  QosGovernor,
  validatePresentationProfile,
  type QosSpec,
} from "../src/qos.ts";

const SPEC: QosSpec = {
  enterConstrainedMs: 22,
  exitConstrainedMs: 15,
  enterCriticalMs: 33,
  exitCriticalMs: 22,
  dwellFrames: 5,
  refusalAfterCriticalFrames: 20,
};

function feed(g: QosGovernor, ms: number, n: number): ReturnType<QosGovernor["sample"]> {
  let last!: ReturnType<QosGovernor["sample"]>;
  for (let i = 0; i < n; i += 1) {
    last = g.sample(ms);
  }
  return last;
}

test("escalation takes EXACTLY dwellFrames; one level at a time", () => {
  const g = new QosGovernor(SPEC);
  // 4 slow frames: still normal (dwell not met).
  assert.equal(feed(g, 30, 4).state, "normal");
  // The 5th escalates — exactly at the dwell count.
  const fifth = g.sample(30);
  assert.equal(fifth.state, "constrained");
  assert.equal(fifth.changed, true);
  assert.equal(fifth.profile.badge, "visual analysis reduced; physics unchanged");
  // 30 ms is below enterCritical (33): stays constrained forever.
  assert.equal(feed(g, 30, 50).state, "constrained");
  // Above 33 for dwell frames: critical (never skipped a level).
  assert.equal(feed(g, 40, 4).state, "constrained");
  assert.equal(g.sample(40).state, "critical");
});

test("chatter resistance: oscillation around a threshold never flaps", () => {
  const g = new QosGovernor(SPEC);
  // Alternate 30 / 10 ms — each fast frame resets the streak, so the
  // state NEVER leaves normal no matter how long this runs.
  for (let i = 0; i < 400; i += 1) {
    const s = g.sample(i % 2 === 0 ? 30 : 10);
    assert.equal(s.state, "normal", `flapped at frame ${i}`);
    assert.equal(s.changed, false);
  }
  // Same resistance inside constrained: hold in the hysteresis band
  // (between exitConstrained 15 and enterConstrained 22) forever.
  feed(g, 30, 5);
  assert.equal(g.current(), "constrained");
  for (let i = 0; i < 400; i += 1) {
    assert.equal(g.sample(18).state, "constrained", "hysteresis band must hold");
  }
});

test("de-escalation needs the dwell below the EXIT threshold", () => {
  const g = new QosGovernor(SPEC);
  feed(g, 40, 10); // → critical
  assert.equal(g.current(), "critical");
  // 20 ms is below exitCritical (22): after dwell → constrained.
  assert.equal(feed(g, 20, 4).state, "critical");
  assert.equal(g.sample(20).state, "constrained");
  // 14 ms below exitConstrained (15): after dwell → normal.
  assert.equal(feed(g, 14, 4).state, "constrained");
  assert.equal(g.sample(14).state, "normal");
});

test("persistent-critical refusal at cap AND cap+1, emitted once", () => {
  const g = new QosGovernor(SPEC);
  feed(g, 40, 10); // normal→constrained (5) →critical (5); streak = 1 at the flip frame
  assert.equal(g.current(), "critical");
  // Feed until one frame BEFORE the budget: no refusal yet (cap−1).
  let refusals = 0;
  let criticalFrames = 1;
  while (criticalFrames < SPEC.refusalAfterCriticalFrames - 1) {
    const s = g.sample(40);
    criticalFrames += 1;
    if (s.refusal !== undefined) {
      refusals += 1;
    }
  }
  assert.equal(refusals, 0, "no refusal before the budget");
  // The budget frame: exactly one typed refusal.
  const at = g.sample(40);
  assert.equal(at.refusal?.code, "performance-budget-missed");
  // And NEVER again while critical persists.
  for (let i = 0; i < 100; i += 1) {
    assert.equal(g.sample(40).refusal, undefined, "emitted once");
  }
  // Recover, degrade again: the refusal re-arms.
  feed(g, 20, 5);
  assert.equal(g.current(), "constrained");
  feed(g, 40, 5); // back to critical
  assert.equal(g.current(), "critical");
  let rearmed = 0;
  for (let i = 0; i < SPEC.refusalAfterCriticalFrames + 5; i += 1) {
    if (g.sample(40).refusal !== undefined) {
      rearmed += 1;
    }
  }
  assert.equal(rearmed, 1, "re-armed after recovery");
});

test("tier immutability: a physics knob in a profile REFUSES", () => {
  // The hostile twin — someone tries to smuggle a solver knob through
  // the presentation profile.
  assert.throws(
    () => validatePresentationProfile({ pixelRatioCap: 1, simTickHz: 60 }),
    /physics-tier immutability/,
  );
  assert.throws(
    () => validatePresentationProfile({ couplingCap: 4 }),
    /physics-tier immutability/,
  );
  // Every shipped profile passes, is frozen, and Normal carries no badge.
  for (const [state, profile] of Object.entries(PROFILES)) {
    assert.doesNotThrow(() =>
      validatePresentationProfile(profile as unknown as Record<string, unknown>),
    );
    assert.ok(Object.isFrozen(profile), `${state} profile must be atomic/frozen`);
  }
  assert.equal(PROFILES.normal.badge, null);
  assert.equal(PROFILES.critical.badge, "visual analysis reduced; physics unchanged");
});

test("spec admission: inverted hysteresis refuses; domain caps on samples", () => {
  assert.throws(
    () => new QosGovernor({ ...SPEC, exitConstrainedMs: 23 }),
    RangeError,
    "exit above enter is not hysteresis",
  );
  const g = new QosGovernor(QOS_V1);
  assert.throws(() => g.sample(Number.NaN), RangeError);
  assert.throws(() => g.sample(-1), RangeError);
});
