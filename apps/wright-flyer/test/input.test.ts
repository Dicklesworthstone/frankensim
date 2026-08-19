// E2.4 input/HUD battery (bead wf-root-guzez.3.4): deterministic slew,
// quantization-grid membership of EVERY emitted command, clamp at full
// travel, recenter decay, key-binding map, HUD triad + flag lines.
// Repro: node --test test/input.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";
import { NEUTRAL, keysFrom, stepCommand } from "../src/input.ts";
import { hudLines } from "../src/hud.ts";
import { computePose } from "../src/airframe/pose.ts";

const K = {
  canardUp: false, canardDown: false, warpLeft: false, warpRight: false, recenter: false,
};

test("slew is deterministic and every command sits on the 1/4096 grid", () => {
  let cmd = NEUTRAL;
  const dt = 1 / 120;
  for (let i = 0; i < 240; i++) {
    cmd = stepCommand(cmd, { ...K, canardUp: true, warpRight: true }, dt);
    assert.ok(Math.abs(cmd.canard * 4096 - Math.round(cmd.canard * 4096)) < 1e-9, "grid");
    assert.ok(Math.abs(cmd.warp * 4096 - Math.round(cmd.warp * 4096)) < 1e-9, "grid");
  }
  // Two seconds at 1.4/s -> clamped at full travel; warp at 1.1/s -> 1.0 too.
  assert.equal(cmd.canard, 1);
  assert.equal(cmd.warp, 1);
  // Determinism: same sequence, same bits.
  let again = NEUTRAL;
  for (let i = 0; i < 240; i++) {
    again = stepCommand(again, { ...K, canardUp: true, warpRight: true }, dt);
  }
  assert.deepEqual(cmd, again);
});

test("opposed keys hold; recenter decays toward neutral", () => {
  let cmd = { ...NEUTRAL, canard: 0.5, warp: -0.5 };
  const held = stepCommand(cmd, { ...K, canardUp: true, canardDown: true }, 1 / 120);
  assert.equal(held.canard, 0.5, "opposed keys cancel");
  for (let i = 0; i < 240; i++) {
    cmd = stepCommand(cmd, { ...K, recenter: true }, 1 / 120);
  }
  assert.ok(Math.abs(cmd.canard) < 0.01 && Math.abs(cmd.warp) < 0.01, "recenter decays");
  assert.throws(() => stepCommand(NEUTRAL, K, Number.NaN), /finite/);
});

test("key bindings map both arrow and WASD, pull = nose up", () => {
  const k = keysFrom(new Set(["ArrowDown", "KeyD", "Space"]));
  assert.ok(k.canardUp && k.warpRight && k.recenter && !k.canardDown && !k.warpLeft);
  const wasd = keysFrom(new Set(["KeyS", "KeyA"]));
  assert.ok(wasd.canardUp && wasd.warpLeft);
});

test("HUD renders the period triad and surfaces the rig flags", () => {
  const pose = computePose({ canardDeg: 31, warpDeg: 0, rudderDeg: 0, coupled: true, propAngleRad: 0 });
  const lines = hudLines({ airspeedMps: 13.86, elapsedS: 12, engineRpm: 1025, camera: "chase", pose });
  assert.match(lines[0]!, /31\.0 mph/);
  assert.match(lines[1]!, /12\.0 s/);
  assert.match(lines[2]!, /1025 rpm/);
  assert.ok(lines.includes("CONTROL AT STOP"), "the clamp flag must surface");
  const clean = hudLines({
    airspeedMps: 0, elapsedS: 0, engineRpm: 0, camera: "arrival",
    pose: computePose({ canardDeg: 0, warpDeg: 0, rudderDeg: 0, coupled: true, propAngleRad: 0 }),
  });
  assert.equal(clean.length, 4, "no flags when clean");
});

test("all five camera presets work against scripted state", async () => {
  const { cameraFor, PRESET_KEYS } = await import("../src/camera.ts");
  const launch: [number, number, number] = [0, 4.1, -625];
  const aircraft: [number, number, number] = [12, 5.3, -625];
  const seen = new Set<string>();
  for (const preset of Object.values(PRESET_KEYS)) {
    const shot = cameraFor(preset, 6.5, launch, aircraft);
    for (const v of [...shot.pos, ...shot.look]) assert.ok(Number.isFinite(v), preset);
    assert.ok(shot.pos[1] > 0, `${preset} camera above ground`);
    seen.add(shot.pos.map((v) => v.toFixed(2)).join(","));
    if (preset !== "daniels" && preset !== "free") {
      assert.ok(Math.hypot(shot.pos[0] - aircraft[0], shot.pos[2] - aircraft[2]) < 25, preset);
    }
  }
  assert.equal(seen.size, 5, "five DISTINCT viewpoints");
  const d1 = cameraFor("daniels", 1, launch, aircraft);
  const d2 = cameraFor("daniels", 9, launch, [50, 8, -600]);
  assert.deepEqual(d1.pos, d2.pos, "the tripod never moves");
});
