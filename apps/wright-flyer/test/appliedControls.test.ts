import assert from "node:assert/strict";
import { test } from "node:test";
import { AppliedControls, appliedControlTraceError } from "../src/sim/appliedControls.ts";

// Wire-contract fixture only: no substitute flight physics or digest.
function ok(tick: number): string {
  return JSON.stringify({ ok: {
    tick, phase: "on-rail", assist_active: false,
    x_m: 0, h_m: 0, u_mps: 0, w_mps: 0, q_rad_s: 0, theta_rad: 0,
    p_rad_s: 0, phi_rad: 0, r_rad_s: 0, psi_rad: 0, dc_rad: 0, warp_rad: 0,
    omega_prop_rad_s: 0, gust_w_mps: 0,
  } });
}

test("actual consumed commands survive superseded samples, holds and JSON exactly", () => {
  const live = new AppliedControls();
  assert.equal(live.step(1, () => { throw new Error("waiting must not call the engine"); }).kind, "waiting");
  live.admit(1, { leverForceN: 100, warpCmdRad: 0.1 }, 5, 0);
  live.admit(2, { leverForceN: -0, warpCmdRad: 0.02 }, 8, 0);
  const observed: unknown[] = [];
  for (let tick = 1; tick <= 4; tick += 1) {
    if (tick === 3) live.admit(3, { leverForceN: 0, warpCmdRad: -0 }, -2, 2);
    assert.equal(live.step(tick, (value) => { observed.push(value); return ok(tick); }).kind, "ok");
  }
  const trace = live.seal();
  assert.deepEqual(trace, { endTick: 4, changes: [
    { tick: 1, leverForceN: "-0", warpCmdRad: 0.02 },
    { tick: 3, leverForceN: 0, warpCmdRad: "-0" },
  ] });
  const playback = new AppliedControls(JSON.parse(JSON.stringify(trace)));
  assert.equal(playback.admit(4, { leverForceN: 220, warpCmdRad: 0.148 }, 1, 0), null);
  const replayed: unknown[] = [];
  for (let tick = 1; tick <= 4; tick += 1) {
    assert.equal(playback.step(tick, (value) => { replayed.push(value); return ok(tick); }).kind, "ok");
  }
  assert.deepEqual(replayed, observed, "including negative zero bits, without re-quantization");
  assert.deepEqual(playback.seal(), trace);
  assert.equal(playback.step(5, () => { throw new Error("cannot extend source trace"); }).kind, "malformed");
});

test("refused, malformed and wrong-tick native steps never produce a replayable trace", () => {
  for (const result of ["not-json", '{"refusal":{"code":"x","message":"x","ranked_repairs":[]}}', ok(3)]) {
    const live = new AppliedControls();
    live.admit(1, { leverForceN: 1, warpCmdRad: 0 }, 1, 0);
    assert.equal(live.step(1, () => ok(1)).kind, "ok");
    assert.notEqual(live.step(2, () => result).kind, "ok");
    assert.throws(() => live.seal(), /failed native step/);
    assert.equal(live.step(2, () => { throw new Error("invalid run cannot resume"); }).kind, "malformed");
  }
  const trapped = new AppliedControls();
  trapped.admit(1, { leverForceN: 1, warpCmdRad: 0 }, 1, 0);
  assert.equal(trapped.step(1, () => { throw new Error("native trap after partial advancement"); }).kind, "malformed");
  assert.throws(() => trapped.seal(), /failed native step/);
});

test("trace admission refuses missing tick-one, unordered, excessive and nonfinite inputs", () => {
  const first = { tick: 1, leverForceN: 0, warpCmdRad: 0 };
  for (const trace of [null, {}, { endTick: 0, changes: [] },
    { endTick: 72_001, changes: [first] }, { endTick: 2, changes: [{ ...first, tick: 2 }] },
    { endTick: 2, changes: [first, first] }, { endTick: 1, changes: [{ ...first, warpCmdRad: Infinity }] },
    { endTick: 1, changes: [{ ...first, leverForceN: "0" }] },
  ]) assert.notEqual(appliedControlTraceError(trace), null, JSON.stringify(trace));
});
