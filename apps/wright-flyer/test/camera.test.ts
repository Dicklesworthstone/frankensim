// Camera presentation battery: the easing law is frame-rate
// independent, monotone, endpoint-convergent, and — critically —
// NEVER throws on hostile dt (tab-suspend rAF gaps feed huge values;
// a RangeError here would kill the render loop). Repro:
// node --test test/camera.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  BASE_FOV_DEG,
  easeCameraToward,
  speedFov,
  type CameraState,
} from "../src/camera.ts";

const start: CameraState = { pos: [0, 0, 0], look: [10, 0, 0], fovDeg: BASE_FOV_DEG };
const target = { pos: [20, 10, 6] as [number, number, number], look: [30, 5, 0] as [number, number, number] };

test("easing approaches monotonically and converges at the target", () => {
  let s = start;
  let prevDist = Infinity;
  for (let i = 0; i < 200; i++) {
    s = easeCameraToward(s, target, 1 / 60);
    const d = Math.hypot(
      s.pos[0]! - target.pos[0]!,
      s.pos[1]! - target.pos[1]!,
      s.pos[2]! - target.pos[2]!,
    );
    assert.ok(d <= prevDist + 1e-12, "monotone approach");
    prevDist = d;
  }
  assert.ok(prevDist < 0.01, `converged (dist ${prevDist.toExponential(2)})`);
});
test("hostile dt never throws; huge dt effectively snaps to target", () => {
  const big = easeCameraToward(start, target, 3600); // tab was suspended
  // Clamped dt=1 covers 99.6% of the gap in one frame — a visual snap.
  assert.ok(
    Math.abs(big.pos[0]! - target.pos[0]!) < 0.2,
    `huge dt residue ${Math.abs(big.pos[0]! - target.pos[0]!).toFixed(3)} m`,
  );
  const negative = easeCameraToward(start, target, -5); // clock jitter
  assert.deepEqual(negative.pos, start.pos, "dt<0 holds position");
  const nan = easeCameraToward(start, target, Number.NaN);
  assert.ok(
    Math.abs(nan.pos[0]! - target.pos[0]!) < 0.2,
    "NaN dt recovers instead of crashing",
  );
});

test("frame-rate independence: two 30fps steps equal one 60fps step pair", () => {
  // Exponential law: (1-e^{-k·2dt}) vs composed (1-e^{-k·dt})² — the
  // composed path must land within a few percent of the single step.
  const single = easeCameraToward(start, target, 1 / 15);
  let s = start;
  s = easeCameraToward(s, target, 1 / 30);
  s = easeCameraToward(s, target, 1 / 30);
  const err =
    Math.abs(single.pos[0]! - s.pos[0]!) /
    Math.max(1, Math.abs(target.pos[0]!));
  assert.ok(err < 0.05, `frame-rate gap ${err.toFixed(4)}`);
});

test("speedFov widens with airspeed inside a tight band", () => {
  assert.equal(speedFov(BASE_FOV_DEG, 0), BASE_FOV_DEG);
  assert.equal(speedFov(BASE_FOV_DEG, 9), BASE_FOV_DEG, "no widening below 9 m/s");
  const max = speedFov(BASE_FOV_DEG, 90);
  assert.equal(max, BASE_FOV_DEG + 8.4, "clamped at +8.4 deg");
  assert.ok(speedFov(BASE_FOV_DEG, 20) > BASE_FOV_DEG);
});
