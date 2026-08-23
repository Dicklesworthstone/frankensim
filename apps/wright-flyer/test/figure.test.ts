// Figure-math battery (bead guzez.14): anthropometric identity per
// brother (per-segment oracles against the Drillis-Contini fractions,
// not round-trips), gait antiphase/limits/standing-identity, stride
// cadence monotonicity, arm-aim trigonometry at the axis points, and
// pose-constant sanity.
// Repro: node --test test/figure.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  BINOCULAR_POSE,
  GAIT_MAX_MPS,
  KNEE_FLEX_MAX_RAD,
  PRONE_POSE,
  armAimAngles,
  figureSpec,
  gaitPose,
  strideFreqHz,
} from "../src/figure.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-figure","case":"${kase}",${payload}}`);
}

test("anthropometry: statures, builds, and segment fractions", () => {
  const w = figureSpec("wilbur");
  const o = figureSpec("orville");
  assert.equal(w.heightM, 1.78);
  assert.equal(o.heightM, 1.73);
  assert.ok(w.build < 1 && o.build > 1, "Wilbur wiry, Orville stockier");
  // Per-segment oracles (hand-computed from the D-C fractions).
  assert.ok(Math.abs(w.thighLenM - 0.245 * 1.78) < 1e-12);
  assert.ok(Math.abs(o.shinLenM - 0.246 * 1.73) < 1e-12);
  assert.ok(Math.abs(w.upperArmLenM - 0.186 * 1.78) < 1e-12);
  // Torso closes the chain: hip height + torso = shoulder height.
  assert.ok(Math.abs(w.hipHeightM + w.torsoLenM - w.shoulderHeightM) < 1e-12);
  // Wilbur is taller everywhere despite the slighter build.
  assert.ok(w.thighLenM > o.thighLenM && w.forearmLenM > o.forearmLenM);
  // Orville's shoulders are relatively wider (build term).
  assert.ok(o.shoulderWidthM / o.heightM > w.shoulderWidthM / w.heightM);
  jlog("anthropometry", `"wilbur_m":${w.heightM},"orville_m":${o.heightM}`);
});

test("gait: legs and arms antiphase, knees never negative, caps hold", () => {
  for (let phase = 0; phase < 2 * Math.PI; phase += 0.05) {
    const p = gaitPose(phase, GAIT_MAX_MPS);
    const q = gaitPose(phase + Math.PI, GAIT_MAX_MPS);
    assert.ok(Math.abs(p.hipL - q.hipR) < 1e-12, "left leg mirrors right half a stride later");
    assert.ok(Math.abs(p.kneeL - q.kneeR) < 1e-12);
    assert.ok(p.kneeL >= 0 && p.kneeR >= 0, "no knee hyperextension");
    assert.ok(p.kneeL <= KNEE_FLEX_MAX_RAD + 1e-12, "flexion cap");
    // Same-side arm counter-swings its leg (opposite signs when moving).
    if (Math.abs(p.hipL) > 1e-6) {
      assert.ok(p.shoulderL * p.hipL < 0, "arm opposes leg");
    }
    assert.ok(p.bobM >= 0 && p.bobM <= 0.036, "bob bounded");
    // Bob has TWO peaks per stride: period pi.
    const r = gaitPose(phase + Math.PI, GAIT_MAX_MPS);
    assert.ok(Math.abs(p.bobM - r.bobM) < 1e-12, "bob period is half a stride");
  }
  jlog("gait-antiphase", `"samples":126`);
});

test("gait: zero speed is EXACTLY the standing pose", () => {
  for (const phase of [0, 1.3, 4.4]) {
    const p = gaitPose(phase, 0);
    assert.equal(p.hipL, 0);
    assert.equal(p.kneeR, 0);
    assert.equal(p.leanRad, 0);
    assert.equal(p.bobM, 0);
    assert.equal(p.elbowL, 0.1, "relaxed elbow, not locked");
  }
  assert.throws(() => gaitPose(Number.NaN, 1), RangeError);
  assert.throws(() => gaitPose(0, Number.POSITIVE_INFINITY), RangeError);
  jlog("gait-standing", `"zero_is_identity":true`);
});

test("gait scales monotonically with speed", () => {
  let prevHip = -1;
  let prevLean = -1;
  for (const v of [0.5, 1.5, 3, 4.5, GAIT_MAX_MPS, GAIT_MAX_MPS + 3]) {
    const p = gaitPose(Math.PI / 2, v);
    assert.ok(p.hipL >= prevHip, `hip amplitude monotone at ${v}`);
    assert.ok(p.leanRad >= prevLean, `lean monotone at ${v}`);
    prevHip = p.hipL;
    prevLean = p.leanRad;
  }
  // Above the cap the gait saturates (never cartoonish).
  assert.deepEqual(gaitPose(1, GAIT_MAX_MPS), gaitPose(1, 99));
  jlog("gait-speed", `"saturates_at":${GAIT_MAX_MPS}`);
});

test("stride cadence: zero at rest, monotone in speed", () => {
  const leg = figureSpec("orville").hipHeightM;
  assert.equal(strideFreqHz(0, leg), 0);
  let prev = 0;
  for (const v of [0.5, 1, 2, 3, 4, 5]) {
    const f = strideFreqHz(v, leg);
    assert.ok(f > prev, `cadence rises at ${v} m/s`);
    prev = f;
  }
  assert.ok(prev > 1.5 && prev < 4, `running cadence plausible (${prev.toFixed(2)} Hz)`);
  assert.throws(() => strideFreqHz(1, 0), RangeError);
  jlog("cadence", `"top_hz":${prev.toFixed(3)}`);
});

test("arm aiming: axis points and refusals", () => {
  // Straight down (the rest pose): pitch 0.
  assert.ok(Math.abs(armAimAngles(0, -1, 0).pitchRad) < 1e-12);
  // Straight forward: pitch pi/2, yaw 0.
  const fwd = armAimAngles(1, 0, 0);
  assert.ok(Math.abs(fwd.pitchRad - Math.PI / 2) < 1e-12);
  assert.ok(Math.abs(fwd.yawRad) < 1e-12);
  // Straight out to the RIGHT (+z): yaw -pi/2 under the YZX composition.
  const right = armAimAngles(0, 0, 1);
  assert.ok(Math.abs(right.yawRad + Math.PI / 2) < 1e-12);
  // Straight up: pitch pi.
  assert.ok(Math.abs(armAimAngles(0, 1, 0).pitchRad - Math.PI) < 1e-12);
  assert.throws(() => armAimAngles(0, 0, 0), RangeError);
  jlog("arm-aim", `"axes":4`);
});

test("pose constants stay within joint limits", () => {
  // Rz(+θ) on the +x face tips it toward +y: positive = looking UP
  // (regression: the original -0.55 stared into the wing fabric).
  assert.ok(PRONE_POSE.headPitchRad > 0, "prone head looks UP");
  assert.ok(PRONE_POSE.shoulderForwardRad > 0 && PRONE_POSE.shoulderForwardRad < Math.PI);
  assert.ok(BINOCULAR_POSE.shoulderForwardRad > Math.PI / 2, "hands reach the face");
  assert.ok(BINOCULAR_POSE.elbowFlexRad > 0 && BINOCULAR_POSE.elbowFlexRad < 2.2);
  jlog("poses", `"prone_head":${PRONE_POSE.headPitchRad}`);
});
