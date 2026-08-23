// Launch-ceremony battery (bead guzez.16): null inputs -> zeros,
// RangeError domain refusals on negative/non-finite ages, monotone
// ramps sampled per-item (a mean would hide one bad sample), the kick
// peaking EXACTLY at the latch instant, closed envelopes landing on
// exactly 0, and determinism (same inputs -> same numbers).
// Repro: node --test test/ceremony.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  FLASH_ATTACK_S,
  FLASH_DECAY_S,
  KICK_FOV_PEAK_DEG,
  KICK_SHAKE_PEAK_M,
  KICK_T_S,
  NOMINAL_RAIL_RUN_S,
  RAMP_S,
  RELEASE_DECAY_S,
  flashPulse,
  glanceBlend,
  releaseKick,
  wingtipGlanceBlend,
} from "../src/ceremony.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-ceremony","case":"${kase}",${payload}}`);
}

test("null inputs give zeros everywhere", () => {
  assert.equal(glanceBlend(null, null), 0);
  assert.equal(glanceBlend(1.5, null), 0, "rail age alone still ramps");
  assert.equal(glanceBlend(null, null), 0);
  const k = releaseKick(null);
  assert.equal(k.fovKickDeg, 0);
  assert.equal(k.shakeAmpM, 0);
  assert.equal(flashPulse(null), 0);
  jlog("null-inputs", `"glance":0,"kick":0,"flash":0`);
});

test("negative or non-finite times are refused (RangeError)", () => {
  for (const bad of [-0.001, -1, Number.NaN, Number.POSITIVE_INFINITY]) {
    assert.throws(() => glanceBlend(bad, null), RangeError, `elapsedOnRailS=${bad}`);
    assert.throws(() => glanceBlend(null, bad), RangeError, `sinceReleaseS=${bad}`);
    assert.throws(() => releaseKick(bad), RangeError, `sinceReleaseS=${bad}`);
    assert.throws(() => flashPulse(bad), RangeError, `sinceFlashS=${bad}`);
    assert.throws(
      () => wingtipGlanceBlend(bad, true, 0.5),
      RangeError,
      `nowS=${bad}`,
    );
    assert.throws(
      () => wingtipGlanceBlend(1, true, bad),
      RangeError,
      `releaseImminentT=${bad}`,
    );
  }
  jlog("domain-refusals", `"cases":6`);
});

test("glance ramps monotonically across the FINAL RAMP_S of the rail run", () => {
  const rampStart = NOMINAL_RAIL_RUN_S - RAMP_S;
  // Before the window: flat zero.
  assert.equal(glanceBlend(rampStart - RAMP_S, null), 0);
  assert.equal(glanceBlend(0, null), 0);
  // Through the window: strictly increasing, smoothstep-flat at ends.
  let prev = -1;
  for (let i = 0; i <= 24; i += 1) {
    const t = rampStart + (RAMP_S * i) / 24;
    const w = glanceBlend(t, null);
    assert.ok(w >= prev, `weight must be nondecreasing at t=${t.toFixed(3)}: ${w} < ${prev}`);
    assert.ok(w >= 0 && w <= 1, `weight in [0,1] at t=${t.toFixed(3)}, got ${w}`);
    prev = w;
  }
  assert.equal(prev, 1, "holds at full weight by nominal release");
  // Past the nominal release while still fed rail time: clamped hold.
  assert.equal(glanceBlend(NOMINAL_RAIL_RUN_S + 2, null), 1);
  jlog("glance-ramp", `"ramp_s":${RAMP_S},"peak":${prev}`);
});

test("glance decays from 1 to exactly 0 over RELEASE_DECAY_S after release", () => {
  assert.equal(glanceBlend(null, 0), 1, "decay starts AT full weight");
  let prev = Number.POSITIVE_INFINITY;
  for (let i = 0; i <= 20; i += 1) {
    const since = (RELEASE_DECAY_S * i) / 20;
    const w = glanceBlend(null, since);
    assert.ok(w <= prev, `decay must be nonincreasing at since=${since.toFixed(3)}`);
    assert.ok(w >= 0 && w <= 1, `weight in [0,1] at since=${since.toFixed(3)}, got ${w}`);
    prev = w;
  }
  assert.equal(prev, 0, "fully home by RELEASE_DECAY_S");
  assert.equal(glanceBlend(null, RELEASE_DECAY_S * 3), 0, "stays home after");
  jlog("glance-decay", `"decay_s":${RELEASE_DECAY_S}`);
});

test("wingtipGlanceBlend first form: off-rail snaps to 0, ramps via imminent time", () => {
  assert.equal(wingtipGlanceBlend(10, false, 0.2), 0, "off-rail -> snapped home");
  assert.equal(wingtipGlanceBlend(1, true, null), 0, "no projected release -> no glance");
  assert.equal(wingtipGlanceBlend(1, true, RAMP_S + 1), 0, "far from release -> 0");
  assert.equal(wingtipGlanceBlend(1, true, 0), 1, "at release -> full glance");
  const mid = wingtipGlanceBlend(1, true, RAMP_S / 2);
  assert.ok(mid > 0.4 && mid < 0.7, `mid-ramp near 0.5, got ${mid}`);
  jlog("first-form", `"mid":${mid.toFixed(3)}`);
});

test("releaseKick peaks EXACTLY at the latch and settles to hard zero", () => {
  const peak = releaseKick(0);
  assert.equal(peak.fovKickDeg, KICK_FOV_PEAK_DEG, "FOV punch exact at t=0");
  assert.equal(peak.shakeAmpM, KICK_SHAKE_PEAK_M, "shake amp exact at t=0");
  let prevFov = Number.POSITIVE_INFINITY;
  let prevShake = Number.POSITIVE_INFINITY;
  for (let i = 1; i <= 32; i += 1) {
    const since = (KICK_T_S * i) / 32;
    const k = releaseKick(since);
    assert.ok(k.fovKickDeg <= prevFov, `fov decays monotonically at ${since.toFixed(3)}`);
    assert.ok(k.shakeAmpM <= prevShake, `shake decays monotonically at ${since.toFixed(3)}`);
    assert.ok(k.fovKickDeg >= 0 && k.shakeAmpM >= 0);
    prevFov = k.fovKickDeg;
    prevShake = k.shakeAmpM;
  }
  const done = releaseKick(KICK_T_S);
  assert.equal(done.fovKickDeg, 0, "closed envelope: exactly 0 at KICK_T_S");
  assert.equal(done.shakeAmpM, 0);
  const late = releaseKick(KICK_T_S * 5);
  assert.equal(late.fovKickDeg, 0);
  assert.equal(late.shakeAmpM, 0);
  jlog("kick", `"peak_fov":${KICK_FOV_PEAK_DEG},"peak_shake":${KICK_SHAKE_PEAK_M}`);
});

test("flashPulse: attack rises to 1, decay falls to exactly 0, stays dark", () => {
  assert.equal(flashPulse(0), 0, "dark before the shutter");
  let prev = -1;
  for (let i = 0; i <= 8; i += 1) {
    const t = (FLASH_ATTACK_S * i) / 8;
    const v = flashPulse(t);
    assert.ok(v >= prev, `attack monotone up at t=${t.toFixed(4)}`);
    assert.ok(v >= 0 && v <= 1, `opacity in [0,1] at t=${t.toFixed(4)}, got ${v}`);
    prev = v;
  }
  assert.equal(prev, 1, "full white at end of attack");
  prev = Number.POSITIVE_INFINITY;
  for (let i = 0; i <= 28; i += 1) {
    const t = FLASH_ATTACK_S + (FLASH_DECAY_S * i) / 28;
    const v = flashPulse(t);
    assert.ok(v <= prev, `decay monotone down at t=${t.toFixed(4)}`);
    prev = v;
  }
  assert.equal(prev, 0, "exactly dark at attack+FLASH_DECAY_S");
  assert.equal(flashPulse(FLASH_ATTACK_S + FLASH_DECAY_S + 1), 0, "stays dark");
  jlog("flash", `"attack_s":${FLASH_ATTACK_S},"decay_s":${FLASH_DECAY_S}`);
});

test("every envelope is deterministic: same inputs, same numbers", () => {
  for (const t of [0, 0.03, 0.37, 1.1, 2.9]) {
    assert.deepEqual(glanceBlend(t, t + 0.5), glanceBlend(t, t + 0.5));
    assert.deepEqual(releaseKick(t), releaseKick(t));
    assert.deepEqual(flashPulse(t), flashPulse(t));
  }
  assert.deepEqual(
    glanceBlend(NOMINAL_RAIL_RUN_S - RAMP_S / 2, null),
    glanceBlend(NOMINAL_RAIL_RUN_S - RAMP_S / 2, null),
  );
  jlog("determinism", `"probes":5`);
});
