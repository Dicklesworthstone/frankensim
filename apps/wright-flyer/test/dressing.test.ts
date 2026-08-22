// Scene-dressing battery (bead guzez.13): flock determinism + caps at
// cap AND cap+1, per-gull bounds as PER-ITEM oracles (a sum or mean
// would be blind to one runaway gull), heading tangency against the
// finite-difference velocity, Orville's speed limit sampled along the
// whole chase, camp props out of the rail corridor, and rail ties
// covering both ends.
// Repro: node --test test/dressing.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  GULL_HEIGHT_MAX,
  GULL_HEIGHT_MIN,
  MAX_GULLS,
  ORVILLE_MAX_MPS,
  RAIL_CLEAR_HALF_WIDTH_M,
  SCRUB_FLAT_RADIUS_M,
  campLayout,
  flagPoint,
  gullAttitude,
  gullFleet,
  gullPose,
  landingDust,
  lcg,
  orvillePose,
  railTies,
  scrubField,
  smokePuff,
  streamerPoint,
} from "../src/dressing.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-dressing","case":"${kase}",${payload}}`);
}

test("flock is deterministic and capped at cap AND cap+1", () => {
  const a = gullFleet(MAX_GULLS, 1903);
  const b = gullFleet(MAX_GULLS, 1903);
  assert.deepEqual(a, b, "same seed, same flock");
  const c = gullFleet(8, 1904);
  assert.notDeepEqual(a.slice(0, 8), c, "different seed, different flock");
  assert.equal(a.length, MAX_GULLS, "AT cap admits");
  assert.throws(() => gullFleet(MAX_GULLS + 1, 1903), RangeError, "cap+1 refuses");
  assert.throws(() => gullFleet(0, 1903), RangeError);
  assert.throws(() => gullFleet(2.5, 1903), RangeError);
  jlog("flock-caps", `"max":${MAX_GULLS}`);
});

test("every gull stays in its band (per-item oracle, sampled)", () => {
  const fleet = gullFleet(12, 42);
  for (const [i, g] of fleet.entries()) {
    for (let t = 0; t <= 120; t += 1.7) {
      const p = gullPose(g, t);
      assert.ok(Number.isFinite(p.x) && Number.isFinite(p.z), `gull ${i} finite`);
      assert.ok(
        p.y >= GULL_HEIGHT_MIN - 2 && p.y <= GULL_HEIGHT_MAX + 2,
        `gull ${i} height ${p.y} at t=${t}`,
      );
      assert.ok(Math.abs(p.flapRad) <= 0.66, `gull ${i} flap bounded`);
      const r = Math.hypot(p.x - g.cx, p.z - g.cz);
      assert.ok(Math.abs(r - g.radius) < 1e-9, `gull ${i} on its circle`);
    }
  }
  jlog("gull-bounds", `"gulls":12,"samples":71`);
});

test("gull heading is tangent to its orbit (finite-difference check)", () => {
  const g = gullFleet(1, 7)[0]!;
  for (const t of [0, 3.1, 47.9]) {
    const p0 = gullPose(g, t);
    const p1 = gullPose(g, t + 1e-4);
    // Compare direction vectors, not raw angles (branch cuts): the
    // claimed heading must align with the finite-difference velocity.
    const vx = Math.cos(p0.headingRad);
    const vz = Math.sin(p0.headingRad);
    const dx = (p1.x - p0.x) / 1e-4;
    const dz = (p1.z - p0.z) / 1e-4;
    const mag = Math.hypot(dx, dz);
    assert.ok(mag > 0);
    const dot = (vx * dx + vz * dz) / mag;
    assert.ok(dot > 0.999, `heading aligns with velocity at t=${t} (dot ${dot})`);
  }
  jlog("gull-heading", `"tangent":true`);
});

test("Orville never exceeds his legs and raises the glasses after release", () => {
  // Chase phase: sample the whole run; his x must respect the speed cap
  // from a standing start AND never pass the machine.
  let prevX = 0;
  for (let t = 0.1; t <= 8; t += 0.1) {
    const machineX = 4.0 * t; // the rail run accelerates past him
    const p = orvillePose(t, true, machineX, null, null);
    assert.ok(p.x <= ORVILLE_MAX_MPS * t + 1e-9, `speed cap at t=${t}`);
    assert.ok(p.x <= machineX, `never ahead of the machine at t=${t}`);
    assert.ok(p.x >= prevX - 1e-9, "never runs backward");
    assert.equal(p.glassesUp, false, "no glasses while chasing");
    prevX = p.x;
  }
  // Release: stops within a metre and raises the glasses after the delay.
  const rel = orvillePose(10, false, 80, 12, 10);
  assert.equal(rel.glassesUp, false, "not instantly");
  const later = orvillePose(12.5, false, 120, 12, 10);
  assert.equal(later.glassesUp, true, "glasses up after the beat");
  assert.ok(later.x - 12 <= 0.61, "coasts less than a metre");
  assert.equal(later.gaitRad, 0, "standing still");
  jlog("orville", `"max_mps":${ORVILLE_MAX_MPS}`);
});

test("camp props stay clear of the rail corridor (per-prop oracle)", () => {
  const camp = campLayout();
  for (const p of camp) {
    // The machine travels +x from the launch: anything with x >= -2
    // must clear the wingspan corridor.
    if (p.x >= -2) {
      assert.ok(Math.abs(p.z) > RAIL_CLEAR_HALF_WIDTH_M, `${p.kind} in the corridor`);
    }
    assert.ok(Math.hypot(p.x, p.z) < 60, `${p.kind} within the camp flat`);
  }
  const kinds = new Set(camp.map((p) => p.kind));
  for (const k of ["hangar", "shack", "campfire", "chair", "barrel", "workbench", "toolchest"]) {
    assert.ok(kinds.has(k as never), `camp has a ${k}`);
  }
  // Two buildings, three chairs (the famous camp photos show several).
  assert.equal(camp.filter((p) => p.kind === "chair").length, 3);
  jlog("camp", `"props":${camp.length}`);
});

test("rail ties support both ends and pitch near 1.5 m", () => {
  const ties = railTies(18.3);
  assert.equal(ties[0], 0, "first tie AT the start");
  assert.equal(ties[ties.length - 1], 18.3, "last tie AT the end");
  for (let i = 1; i < ties.length; i += 1) {
    const pitch = ties[i]! - ties[i - 1]!;
    assert.ok(pitch > 1.0 && pitch <= 1.5 + 1e-9, `pitch ${pitch}`);
  }
  assert.throws(() => railTies(0), RangeError);
  assert.throws(() => railTies(Number.NaN), RangeError);
  jlog("rail", `"ties":${ties.length}`);
});

test("lcg is deterministic and in [0,1)", () => {
  const a = lcg(99);
  const b = lcg(99);
  for (let i = 0; i < 1000; i += 1) {
    const va = a();
    assert.equal(va, b());
    assert.ok(va >= 0 && va < 1);
  }
  jlog("lcg", `"draws":1000`);
});

test("scrub field is deterministic and honors the clearing laws", () => {
  const counts = { tufts: 300, bushes: 50, pines: 12 };
  const a = scrubField(counts, 1903);
  const b = scrubField(counts, 1903);
  assert.deepEqual(a, b, "same seed, same field");
  assert.equal(a.length, 362);
  for (const p of a) {
    assert.ok(Math.hypot(p.x, p.z) >= SCRUB_FLAT_RADIUS_M, `outside flat at (${p.x},${p.z})`);
    const inRail = Math.abs(p.z) < 10 && p.x > -8 && p.x < 42;
    assert.ok(!inRail, `out of rail corridor at (${p.x},${p.z})`);
    const inCamp = p.x > -48 && p.x < -8 && p.z > -27 && p.z < -3;
    assert.ok(!inCamp, `out of camp clearing at (${p.x},${p.z})`);
    assert.ok(p.scale > 0.4 && p.scale < 2, `sane scale ${p.scale}`);
  }
  const kinds = new Set(a.map((p) => p.kind));
  assert.deepEqual([...kinds].sort(), ["bush", "pine", "tuft"]);
  assert.throws(() => scrubField({ tufts: 2000, bushes: 0, pines: 0 }, 1), RangeError);
});

test("streamers, flag, and smoke all advect DOWNWIND (-x)", () => {
  for (let i = 0; i < 6; i++) {
    for (let seg = 0; seg <= 8; seg++) {
      const sp = streamerPoint(i, seg, 8, 3.3, 11);
      assert.ok(sp.y >= 0 && sp.y < 1.2, `streamer low at seg ${seg}`);
      const fp = flagPoint(seg, 8, 2.5, 11);
      assert.ok(fp.x <= 0.001, "flag streams toward -x");
      const smoke = smokePuff(i, 9.7, 11, 1);
      assert.ok(smoke.x <= 0.31, "smoke bends downwind");
      assert.ok(smoke.opacity >= 0 && smoke.opacity <= 0.35);
    }
  }
  // Stronger wind stretches the flag farther from the pole.
  const calm = flagPoint(8, 8, 1.0, 2);
  const windy = flagPoint(8, 8, 1.0, 16);
  assert.ok(windy.x < calm.x, "more wind, longer reach");
});

test("Orville runs to the machine after a landing; dust bursts once", () => {
  const watching = orvillePose(10, false, 60, 5, 5);
  assert.equal(watching.glassesUp, true);
  const run = orvillePose(12, false, 60, 5, 5, 60);
  assert.ok(run.gaitRad !== 0, "he is moving toward the machine");
  const arrived = orvillePose(40, false, 60, 5, 5, 60);
  assert.ok(arrived.x >= 57.75, "reached a hand's reach of the wingtip");
  assert.equal(arrived.gaitRad, 0, "stopped when there");
  assert.equal(landingDust(0, 5, 4), null, "no dust before touchdown");
  const d = landingDust(0, 5, 6);
  assert.ok(d !== null && d.opacity > 0 && d.dy >= 0);
  assert.equal(landingDust(0, 5, 9), null, "burst ends");
});
